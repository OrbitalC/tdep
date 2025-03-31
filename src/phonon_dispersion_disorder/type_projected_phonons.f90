module type_projected_phonons
!!
!! Handles the phonons in the commensurates q-points
!!
use konstanter, only: r8, lo_twopi, lo_freqtol, lo_exitcode_param, lo_sqrectol, lo_rectol, lo_tol, &
                      lo_sqtol
use gottochblandat, only: lo_gauss, lo_linspace, lo_chop, lo_progressbar, lo_progressbar_init, walltime, &
                          lo_trueNtimes, lo_trapezoid_integration, lo_sqnorm
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_phonons_disorder, only: lo_phonons_disorder
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
use type_qpointmesh, only: lo_qpoint
use type_unitcell_disorder_projection, only: lo_unitcell_disorder_projection

use type_forceconstant_disorder_secondorder_dynamicalmatrix, only: eigensolve_projected, dynamicalmatrix_projected
use konstanter, only: lo_frequency_hartree_to_meV, lo_frequency_THz_to_Hartree
use gottochblandat, only: lo_invert3x3matrix

implicit none
private
public :: lo_projected_phonons
public :: lo_projected_phonons_qpoints

type lo_projected_phonons_qpoints
    !> The frequencies
    real(r8), dimension(:), allocatable :: omega
    !> The eigenvectors
    complex(r8), dimension(:, :), allocatable :: egv
    !> The group velocities
    real(r8), dimension(:, :, :), allocatable :: vel
    !> The imaginary self-energy
    real(r8), dimension(:, :), allocatable :: im_se
    !> The imaginary self-energy
    real(r8), dimension(:, :), allocatable :: re_se
    !> The lineshape
    real(r8), dimension(:, :), allocatable :: sf
    !> The q-point
    type(lo_qpoint) :: qpoint
contains
    !> Compute the harmonic phonons at a specific q-point
    procedure :: compute_harmonic_qpoint
    !> Project the phonons on one q-point
    procedure :: project_phonons_qpoint
end type

type lo_projected_phonons
    !> The number of qpoints
    integer :: nq
    !> The number of frequencies on the axe
    integer :: nf
    !> The qpoints
    type(lo_projected_phonons_qpoints), dimension(:), allocatable :: qp
contains
    !> Project the phonons on commensurate q-points
    procedure :: project_phonons
end type

contains

!> Project realspace phonons on commensurate q-points
subroutine project_phonons(ph_pr, ph_rs, fc, proj, ss, fmax, nf, is_sf, mem, mw)
    !> The commensurate phonons
    class(lo_projected_phonons), intent(out) :: ph_pr
    !> The realspace phonons
    type(lo_phonons_disorder), intent(in) :: ph_rs
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The projection
    type(lo_unitcell_disorder_projection), intent(in) :: proj
    !> The supercell structure
    type(lo_crystalstructure), intent(in) :: ss
    !> The max frequency on the energy axis
    real(r8), intent(in) :: fmax
    !> The number of points on the frequency axis
    integer :: nf
    !> Do we use the spectral function in the real space phonons ?
    logical, intent(in) :: is_sf
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The frequency axis
    real(r8), dimension(:), allocatable :: e_axis
    !> For the timer
    real(r8) :: t0
    !> Some integers
    integer :: nb, iq, nq, lq

    nb = proj%na_uc * 3
    call mem%allocate(e_axis, nf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call lo_linspace(0.0_r8, fmax, e_axis)

    ! Get the number of commensurate qpoints
    ph_pr%nq = proj%nqpt
    allocate(ph_pr%qp(ph_pr%nq))

    ! Initialize things for the progressbar
    nq = 0
    do iq=1, ph_pr%nq
        if (mod(iq, mw%n) .eq. mw%r) nq = nq + 1
    end do
    if (mw%talk) call lo_progressbar_init()

    ! And project each commensurate q-points, in parallel
    lq = 0
    t0 = walltime()
    do iq=1, ph_pr%nq
        allocate(ph_pr%qp(iq)%omega(nb))
        allocate(ph_pr%qp(iq)%egv(nb, nb))
        allocate(ph_pr%qp(iq)%sf(nf, nb))
        ph_pr%qp(iq)%omega = 0.0_r8
        ph_pr%qp(iq)%egv = 0.0_r8
        ph_pr%qp(iq)%sf = 0.0_r8

        ph_pr%qp(iq)%qpoint%r = proj%cq(iq)%r

        ! We do it in parallel
        if (mod(iq, mw%n) .ne. mw%r) cycle
        lq = lq + 1

        call ph_pr%qp(iq)%compute_harmonic_qpoint(fc, ss, proj, mem)
        call ph_pr%qp(iq)%project_phonons_qpoint(ph_rs, fc, ss, proj, e_axis, is_sf, mem)
        if (mw%talk .and. lo_trueNtimes(iq, 50, nq)) then
            call lo_progressbar(' ... projecting on commensurate modes', lq, nq, walltime() - t0)
        end if
    end do
    if (mw%talk) call lo_progressbar(' ... projecting on commensurate modes', nq, nq, walltime() - t0)

    ! Now we broadcast ?
    do iq=1, ph_pr%nq
        call mw%allreduce('sum', ph_pr%qp(iq)%omega)
        call mw%allreduce('sum', ph_pr%qp(iq)%egv)
        call mw%allreduce('sum', ph_pr%qp(iq)%sf)
    end do

    nb = 59
    if (mw%talk) then
    open(666, file='sf.dat', status='replace')
        do iq=1, nf
            write(666, '(7(1X,F25.12))') e_axis(iq) * lo_frequency_hartree_to_meV, ph_pr%qp(nb)%sf(iq, :) / lo_frequency_hartree_to_meV
        end do
    close(666)
    end if
end subroutine

!> Project phonons on one qpoint
subroutine compute_harmonic_qpoint(phq, fc, ss, proj, mem, withgroupvel)
    !> The phonons at the qpoints
    class(lo_projected_phonons_qpoints), intent(inout) :: phq
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: ss
    !> The projection
    type(lo_unitcell_disorder_projection), intent(in) :: proj
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> Do we compute the groupvelocities ?
    logical, optional, intent(in) :: withgroupvel

    !> The dynamical matrix
    complex(r8), dimension(:, :), allocatable :: dynmat
    !> The dynamical matrix
    complex(r8), dimension(:, :, :), allocatable :: dynmat_grad
    !> Do we compute the groupvel ?
    logical :: groupvel
    !> The number of modes
    integer :: nb

    if (present(withgroupvel)) then
        groupvel = withgroupvel
    else
        groupvel = .false.
    end if

    nb = proj%na_uc * 3
    call mem%allocate(dynmat, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (groupvel) then
        call mem%allocate(dynmat_grad, [3, nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call dynamicalmatrix_projected(fc, proj, phq%qpoint, dynmat, dynmat_grad, ss=ss)
        call eigensolve_projected(fc, phq%qpoint, dynmat, dynmat_grad, omega=phq%omega, eigenvectors=phq%egv, groupvel=phq%vel, mem=mem)
        call mem%deallocate(dynmat_grad, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    else
        call dynamicalmatrix_projected(fc, proj, phq%qpoint, dynmat, ss=ss)
        call eigensolve_projected(fc, phq%qpoint, dynmat, omega=phq%omega, eigenvectors=phq%egv, mem=mem)
    end if
    call mem%deallocate(dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Project phonons on one qpoint
subroutine project_phonons_qpoint(phq, ph_rs, fc, ss, proj, e_axis, is_sf, mem)
    !> The phonons at the qpoints
    class(lo_projected_phonons_qpoints), intent(inout) :: phq
    !> The realspace phonons
    type(lo_phonons_disorder), intent(in) :: ph_rs
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: ss
    !> The projection
    type(lo_unitcell_disorder_projection), intent(in) :: proj
    !> The frequency axis
    real(r8), dimension(:) :: e_axis
    !> Do we project harmonic or a spectral function ?
    logical, intent(in) :: is_sf
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> Eigenvectors for the projection
    complex(r8), dimension(:, :), allocatable :: egv_uc
    !> The weights
    real(r8), dimension(:), allocatable :: w0
    !> The qpoint
    real(r8), dimension(3) :: qv
    !> A complex little buffer
    complex(r8) :: expiqr, c0
    !> Little buffers
    real(r8) :: qdotr, sigma, invf, f0
    !> And the good ol integer for do loops
    integer :: imode, iat, iat_uc, ii, jj, a, nf, nb, istart, iend

    nb = proj%na_uc * 3

    ! First we prepare some flatten unitcell eigenvectors, with the phase factor already in
    call mem%allocate(egv_uc, [ph_rs%nmodes, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    egv_uc = 0.0_r8
    qv = lo_chop(phq%qpoint%r * lo_twopi, lo_sqtol)
    do iat=1, ss%na
        iat_uc = proj%unitcell_index(iat)

        ! We directly include the Fourier factor
        qdotr = dot_product(ss%rcart(:, iat), qv)
        expiqr = cmplx(cos(qdotr), sin(qdotr), r8)
        do a=1, 3
            ii = 3 * (iat - 1) + a
            jj = 3 * (iat_uc - 1) + a
            egv_uc(ii, :) = egv_uc(ii, :) + phq%egv(jj, :) * expiqr
        end do
    end do

    sigma = 0.08_r8 * lo_frequency_THz_to_Hartree
    ! And now we can do the actual projection
    call mem%allocate(w0, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    f0 = real(proj%na_uc, r8) / real(proj%na_ss, r8)
    nf = size(e_axis, 1)
    invf = real(nf, r8) / e_axis(nf)
    do imode=1, ph_rs%nmodes
        ! Skip acoustic gamma modes
        if (ph_rs%omega(imode) .lt. lo_freqtol) cycle

        ! Compute the projection weights of the mode
        c0 = 0.0_r8
        w0 = 0.0_r8
        do ii=1, nb
            ! Skip acoustic gamma modes
            if (phq%omega(ii) .lt. lo_freqtol) cycle

            c0 = sum(egv_uc(:, ii) * ph_rs%egv(:, imode))
            w0(ii) = abs(c0 * conjg(c0)) * f0
        end do

        ! And project
        istart = max(floor((ph_rs%omega(imode) - 4 * sigma)*invf), 1)
        iend = min(ceiling((ph_rs%omega(imode) + 4 * sigma)*invf), nf)
        if (is_sf) then
            do ii=istart, iend
                phq%sf(ii, :) = phq%sf(ii, :) + w0(:) * ph_rs%sf(ii, imode)
            end do
        else
            do ii=istart, iend
                phq%sf(ii, :) = phq%sf(ii, :) + w0(:) * lo_gauss(e_axis(ii), ph_rs%omega(imode), sigma)
            end do
        end if
    end do
    call mem%deallocate(egv_uc, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(w0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Now we normalize the spectral function, simple trapezoidal rule
    do jj=1, nb
        phq%sf(:, jj) = phq%sf(:, jj) / lo_trapezoid_integration(e_axis, phq%sf(:, jj))
    end do
end subroutine
end module
