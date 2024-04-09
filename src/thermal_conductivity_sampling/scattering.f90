#include "precompilerdefinitions"
module new_scattering
use konstanter, only: r8, lo_freqtol, lo_twopi, lo_exitcode_param, lo_hugeint, lo_pi
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_randomnumbers, only: lo_mersennetwister

use options, only: lo_opts


implicit none
private
public :: lo_scattering_rates
public :: compute_scattering

! Container for isotope scattering rates
type iso
    ! The scattering elements
    real(r8), dimension(:, :), allocatable :: psisq
    ! The scattering rates
    real(r8), dimension(:, :), allocatable :: W
end type
! Container for three phonon scattering rates
type threephonon
    ! The scattering elements
    real(r8), dimension(:, :, :), allocatable :: psisq
    ! The scattering rates
    real(r8), dimension(:, :, :), allocatable :: Wplus, Wminus
    ! The qpoints on this rank
    integer, dimension(:), allocatable :: q2, q3
end type
! Container for four phonon scattering rates
type fourphonon
    ! The scattering elements
    real(r8), dimension(:, :, :, :), allocatable :: psisq
    ! The scattering rates
    real(r8), dimension(:, :, :, :), allocatable :: Wpp, Wpm, Wmm
    ! The qpoints on this rank
    integer, dimension(:), allocatable :: q2, q3, q4
end type

! Container for scattering rates
type lo_scattering_rates
    !> The number of qpoint/mode on this rank
    integer :: nlocal_point
    !> The list of qpoint and modes for this rank
    integer, dimension(:), allocatable :: q1, b1
    !> The iso phonon scattering
    type(iso), dimension(:), allocatable :: iso
    !> The three phonon scattering
    type(threephonon), dimension(:), allocatable :: threephonon
    !> The four phonon scattering
    type(fourphonon), dimension(:), allocatable :: fourphonon
    ! The ratio for maximum likelihood estimation of the scattering strengths
    real(r8) :: mle_ratio3ph, mle_ratio4ph
    !> The number of three and four phonon scattering
    integer :: nqpt3ph, nqpt4ph
end type


contains
subroutine compute_scattering(qp, dr, uc, fct, fcf, opts, mw, mem, sr)
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! The third order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    ! The fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    ! The options
    type(lo_opts) :: opts
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> The scattering rate
    type(lo_scattering_rates), intent(out) :: sr

    ! The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    ! The q-point grid dimension
    integer, dimension(3) :: dims
    !> The random number generator
    type(lo_mersennetwister) :: rng
    ! To initialize the random number generator and timing
    real(r8) :: rseed, t0
    ! Some integers
    integer :: q1, b1, q1f, i, j, k, nlocal_point, ctr, nqpt

    ! We start by initializing the random number generator, we need the same on each mpi rank
    rseed = walltime()
    call mw%bcast(rseed, from=mw%n - 1)
    call rng%init(iseed=0, rseed=rseed)

    if (opts%thirdorder .or. opts%fourthorder) then
        call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        qgridfull1 = (/ (i, i=1, qp%n_full_point) /)
    end if
    if (opts%fourthorder) then
        call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        qgridfull2 = (/ (i, i=1, qp%n_full_point) /)
    end if

    ! grid dimensions
    select type(qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! First we distribute qpoint and modes on mpi ranks
    ctr = 0
    nlocal_point = 0
    do q1 =1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! We skip the acoustic mode at Gamma
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            ctr = ctr + 1

            ! MPI thing
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            nlocal_point = nlocal_point + 1
        end do
    end do

    sr%nlocal_point = nlocal_point
    allocate(sr%q1(nlocal_point))
    allocate(sr%b1(nlocal_point))
    if (opts%isotopescattering) allocate(sr%iso(nlocal_point))
    if (opts%thirdorder) allocate(sr%threephonon(nlocal_point))
    if (opts%fourthorder) allocate(sr%fourphonon(nlocal_point))

    i = 0
    ctr = 0
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            ! We skip the acoustic mode at Gamma
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            ctr = ctr + 1

            ! MPI thing
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            i = i + 1
            sr%q1(i) = q1
            sr%b1(i) = b1
        end do
    end do

    if (opts%isotopescattering) then
        do i=1, sr%nlocal_point
            allocate(sr%iso(i)%psisq(dr%n_mode, qp%n_full_point))
            allocate(sr%iso(i)%W(dr%n_mode, qp%n_full_point))
            sr%iso(i)%psisq = 0.0_r8
            sr%iso(i)%W = 0.0_r8
        end do
    end if

    if (opts%thirdorder) then
        nqpt = minval([opts%nsample3ph, qp%n_full_point])
        sr%nqpt3ph = nqpt
        sr%mle_ratio3ph = real(qp%n_full_point, r8) / real(nqpt, r8)

        do i=1, sr%nlocal_point
            allocate(sr%threephonon(i)%q2(nqpt))
            allocate(sr%threephonon(i)%q3(nqpt))
            allocate(sr%threephonon(i)%psisq(dr%n_mode, dr%n_mode, nqpt))
            allocate(sr%threephonon(i)%Wplus(dr%n_mode, dr%n_mode, nqpt))
            allocate(sr%threephonon(i)%Wminus(dr%n_mode, dr%n_mode, nqpt))
            sr%threephonon(i)%psisq = 0.0_r8
            sr%threephonon(i)%Wplus = 0.0_r8
            sr%threephonon(i)%Wminus = 0.0_r8
            call rng%shuffle_int_array(qgridfull1)
            q1f = qp%ip(sr%q1(i))%full_index
            do j=1, nqpt
                sr%threephonon(i)%q2(j) = qgridfull1(j)
                sr%threephonon(i)%q3(j) = fft_third_grid_index(q1f, qgridfull1(j), dims)
            end do
        end do
    end if

    if (opts%fourthorder) then
        nqpt = minval([opts%nsample4ph, qp%n_full_point**2])
        sr%nqpt4ph = nqpt
        sr%mle_ratio4ph = real(qp%n_full_point**2, r8) / real(nqpt, r8)

        do i=1, sr%nlocal_point
            allocate(sr%fourphonon(i)%q2(nqpt))
            allocate(sr%fourphonon(i)%q3(nqpt))
            allocate(sr%fourphonon(i)%q4(nqpt))
            allocate(sr%fourphonon(i)%psisq(dr%n_mode, dr%n_mode, dr%n_mode, nqpt))
            allocate(sr%fourphonon(i)%Wpp(dr%n_mode, dr%n_mode, dr%n_mode, nqpt))
            allocate(sr%fourphonon(i)%Wpm(dr%n_mode, dr%n_mode, dr%n_mode, nqpt))
            allocate(sr%fourphonon(i)%Wmm(dr%n_mode, dr%n_mode, dr%n_mode, nqpt))
            sr%fourphonon(i)%psisq = 0.0_r8
            sr%fourphonon(i)%Wpp = 0.0_r8
            sr%fourphonon(i)%Wpm = 0.0_r8
            sr%fourphonon(i)%Wmm = 0.0_r8
            call rng%shuffle_int_array(qgridfull1)
            call rng%shuffle_int_array(qgridfull2)
            q1f = qp%ip(sr%q1(i))%full_index
            ctr = 0
            counting4ph: do j=1, qp%n_full_point
                do k=1, qp%n_full_point
                    if (ctr .lt. nqpt) then
                        ctr = ctr + 1
                        sr%fourphonon(i)%q2(ctr) = qgridfull1(j)
                        sr%fourphonon(i)%q3(ctr) = qgridfull2(k)
                        sr%fourphonon(i)%q4(ctr) = fft_fourth_grid_index(q1f, qgridfull1(j), qgridfull2(k), dims)
                    else
                        exit counting4ph
                    end if
                end do
            end do counting4ph
        end do
    end if

    ! Deallocate things
    if (opts%thirdorder .or. opts%fourthorder) then
        call mem%deallocate(qgridfull1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if
    if (opts%fourthorder) then
        call mem%deallocate(qgridfull2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    t0 = walltime()
    do i=1, sr%nlocal_point
        if (opts%isotopescattering) then
            call isotope_scattering(i, sr, qp, dr, uc, mw, mem)
        end if
        if (opts%thirdorder) then
            call threephonon_scattering(i, sr, qp, dr, fct, mw, mem)
        end if
        if (opts%fourthorder) then
            call fourphonon_scattering(i, sr, qp, dr, fcf, mw, mem)
        end if
        if (mw%talk .and. lo_trueNtimes(i, 127, sr%nlocal_point)) then
            call lo_progressbar(' ... computing scattering amplitude', i, sr%nlocal_point, walltime() - t0)
        end if
    end do
end subroutine


subroutine isotope_scattering(il, sr, qp, dr, uc, mw, mem)
    !> The local point
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Isotope prefactor
    real(r8), parameter :: isotope_prefactor = lo_pi / 2.0_r8
    ! prefactor
    real(r8) :: prefactor
    ! Eigenvectors
    complex(r8), dimension(uc%na*3, 2) :: egviso
    ! Integers for do loops
    integer :: q1, b1, q2, b2

    q1 = sr%q1(il)
    b1 = sr%b1(il)
    egviso(:, 1) = dr%iq(q1)%egv(:, b1)
    do q2=1, qp%n_full_point
        prefactor = isotope_prefactor * qp%ap(q2)%integration_weight
        do b2=1, dr%n_mode
            egviso(:, 2) = dr%aq(q2)%egv(:, b2)
            sr%iso(il)%psisq(b2, q2) = isotope_scattering_strength(uc, egviso) * prefactor
        end do
    end do
end subroutine


subroutine threephonon_scattering(il, sr, qp, dr, fct, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Fourth order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: threephonon_prefactor = lo_pi / 4.0_r8
    ! Fourier transform of the matrix elements
    complex(r8), dimension(dr%n_mode**3) :: ptf
    ! Frequency scaled eigenvectors
    complex(r8), dimension(dr%n_mode) :: egv1, egv2, egv3
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**2) :: evp1
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**3) :: evp2
    ! The qpoints and the dimension of the qgrid
    real(r8), dimension(3) :: qv2, qv3
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, plf, prefactor
    !
    complex(r8) :: c0
    ! Integers for do loops
    integer :: i, q1, q2, q3, b1, b2, b3

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)

    do i=1, sr%nqpt4ph
        q2 = sr%threephonon(il)%q2(i)
        q3 = sr%threephonon(il)%q3(i)

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        call pretransform_phi3(fct, qv2, qv3, ptf)
        prefactor = threephonon_prefactor * qp%ap(q2)%integration_weight
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)

                evp1 = 0.0_r8
                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                evp2 = conjg(evp2)
                c0 = dot_product(evp2, ptf)
                sr%threephonon(il)%psisq(b2, b3, i) = abs(c0*conjg(c0)) * prefactor
            end do
        end do
    end do
end subroutine


subroutine fourphonon_scattering(il, sr, qp, dr, fcf, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: fourphonon_prefactor = lo_pi / 8.0_r8
    ! Fourier transform of the matrix elements
    complex(r8), dimension(dr%n_mode**4) :: ptf
    ! Frequency scaled eigenvectors
    complex(r8), dimension(dr%n_mode) :: egv1, egv2, egv3, egv4
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**2) :: evp1
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**3) :: evp2
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**4) :: evp3
    ! The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4
    !
    complex(r8) :: c0
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, prefactor
    ! Integers for do loops
    integer :: i, q1, q2, q3, q4, b1, b2, b3, b4

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)

    do i=1, sr%nqpt4ph
        q2 = sr%fourphonon(il)%q2(i)
        q3 = sr%fourphonon(il)%q3(i)
        q4 = sr%fourphonon(il)%q4(i)
        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        prefactor = fourphonon_prefactor * qp%ap(q2)%integration_weight * qp%ap(q3)%integration_weight
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)

                    evp1 = 0.0_r8
                    evp2 = 0.0_r8
                    evp3 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode,    (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                    call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                    call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                    evp3 = conjg(evp3)
                    c0 = dot_product(evp3, ptf)
                    sr%fourphonon(il)%psisq(b2, b3, b4, i) = abs(c0*conjg(c0)) * prefactor
                end do
            end do
        end do
    end do
end subroutine

!> returns the index on the grid that gives q3=-q1-q2
pure function fft_third_grid_index(i1, i2, dims) result(i3)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q3
    integer :: i3

    integer, dimension(3) :: gi1, gi2, gi3
    integer :: l, k

    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    do l = 1, 3
        gi3(l) = 3 - gi1(l) - gi2(l)
    end do
    do k = 1, 3
    do l = 1, 3
        if (gi3(l) .lt. 1) gi3(l) = gi3(l) + dims(l)
        if (gi3(l) .gt. dims(l)) gi3(l) = gi3(l) - dims(l)
    end do
    end do
    ! convert it back to a singlet
    i3 = triplet_to_singlet(gi3, dims(2), dims(3))
end function

!> returns the index on the grid that gives q4=-q3-q2-q1
pure function fft_fourth_grid_index(i1, i2, i3, dims) result(i4)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> index to q3
    integer, intent(in) :: i3
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q4
    integer :: i4

    integer, dimension(3) :: gi1, gi2, gi3, gi4
    integer :: l, k
    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    gi3 = singlet_to_triplet(i3, dims(2), dims(3))
    do l = 1, 3
         gi4(l) = 4 - gi1(l) - gi2(l) - gi3(l)
   end do
   do k = 1, 3
   do l = 1, 3
       if (gi4(l) .lt. 1) gi4(l) = gi4(l) + dims(l)
       if (gi4(l) .gt. dims(l)) gi4(l) = gi4(l) - dims(l)
   end do
   end do

    ! convert it back to a singlet
    i4 = triplet_to_singlet(gi4, dims(2), dims(3))
end function

!> convert a linear index to a triplet
pure function singlet_to_triplet(l, ny, nz) result(gi)
    !> linear index
    integer, intent(in) :: l
    !> second dimension
    integer, intent(in) :: ny
    !> third dimension
    integer, intent(in) :: nz
    !> grid-index
    integer, dimension(3) :: gi

    integer :: i, j, k

    k = mod(l, nz)
    if (k .eq. 0) k = nz
    j = mod((l - k)/nz, ny) + 1
    i = (l - k - (j - 1)*nz)/(nz*ny) + 1
    gi = [i, j, k]
end function

!> convert a triplet index to a singlet
pure function triplet_to_singlet(gi, ny, nz) result(l)
    !> grid-index
    integer, dimension(3), intent(in) :: gi
    !> second dimension
    integer, intent(in) :: ny
    !> third dimension
    integer, intent(in) :: nz
    !> linear index
    integer :: l

    l = (gi(1) - 1)*ny*nz + (gi(2) - 1)*nz + gi(3)
end function

real(r8) function isotope_scattering_strength(uc, egv)
    type(lo_crystalstructure), intent(in) :: uc
    complex(r8), dimension(:, :), intent(in) :: egv
    !
    integer :: i, j
    real(r8) :: f0, f1
    complex(r8), dimension(3) :: cv0, cv1
    !
    f1 = 0.0_r8
    do i = 1, uc%na
        cv0 = egv((i - 1)*3 + 1:(i*3), 1)
        cv1 = egv((i - 1)*3 + 1:(i*3), 2)
        f0 = 0.0_r8
        do j = 1, 3
            f0 = f0 + abs(conjg(cv0(j))*cv1(j))
        end do
        f0 = f0**2
        f1 = f1 + f0*uc%isotope(i)%disorderparameter
    end do
    isotope_scattering_strength = f1
    !
end function

!> Get the Fourier transform of the third order matrix element
subroutine pretransform_phi3(fct, q2, q3, ptf)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3
    real(r8) :: iqr
    integer :: a1, a2, a3, ia, ib, ic, t, nb

    nb = fct%na*3
    ptf = 0.0_r8
    do a1 = 1, fct%na
    do t = 1, fct%atom(a1)%n
        a2 = fct%atom(a1)%triplet(t)%i2
        a3 = fct%atom(a1)%triplet(t)%i3

        rv2 = fct%atom(a1)%triplet(t)%lv2
        rv3 = fct%atom(a1)%triplet(t)%lv3

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            l = (ia - 1)*nb*nb + (ib - 1)*nb + ic
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%mwm(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
end subroutine

!> Get the Fourier transform of the fourth order matrix element
subroutine pretransform_phi4(fcf, q2, q3, q4, ptf)
    !> third order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3, q4
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l, m

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3, rv4
    real(r8) :: iqr
    integer :: a1, a2, a3, a4, ia, ib, ic, id, q, nb

    nb = fcf%na*3
    ptf = 0.0_r8
    do a1 = 1, fcf%na
    do q = 1, fcf%atom(a1)%n
        a2 = fcf%atom(a1)%quartet(q)%i2
        a3 = fcf%atom(a1)%quartet(q)%i3
        a4 = fcf%atom(a1)%quartet(q)%i4

        rv2 = fcf%atom(a1)%quartet(q)%lv2
        rv3 = fcf%atom(a1)%quartet(q)%lv3
        rv4 = fcf%atom(a1)%quartet(q)%lv4

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3) + dot_product(q4, rv4)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do l = 1, 3
        do k = 1, 3
        do j = 1, 3
        do i = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            id = (a4 - 1)*3 + l
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            m = (ia - 1)*nb*nb*nb + (ib - 1)*nb*nb + (ic - 1)*nb + id
            ptf(m) = ptf(m) + fcf%atom(a1)%quartet(q)%mwm(i, j, k, l)*expiqr
        end do
        end do
        end do
        end do
    end do
    end do
end subroutine

end module
