#include "precompilerdefinitions"
module selfenergy
use konstanter, only: r8, lo_freqtol, lo_huge, lo_hugeint, lo_pi, lo_twopi, lo_exitcode_param
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_gauss, lo_planck
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use hdf5_wrappers, only: lo_hdf5_helper
use type_blas_lapack_wrappers, only: lo_dgemv
use quadratic_programming, only: lo_solve_quadratic_program

implicit none
private
public :: lo_selfenergy

type lo_selfenergy
    !> The weight for each self energy
    real(r8), dimension(:, :, :), allocatable :: im_weight
    !> The frequency of the basis set
    real(r8), dimension(:), allocatable :: omega_n
    !> The max frequency
    real(r8) :: omega_max
    !> The width of each lorentzian in the basis set
    real(r8) :: width
    !> The number of basis elements
    integer :: nbasis
contains
    !> initialize the selfenergy
    procedure :: initialize => initialize_selfenergy
    !> Compute imaginary self energy
    procedure :: compute => compute_selfenergy
    !> Compute the selfenergy for a specific frequency
!   procedure :: predict
    !> destroy
    procedure :: destroy => destroy_selfenergy
    !> Write to hdf5
    procedure :: write_to_hdf5 => write_selfenergy_to_hdf5
end type

contains

!> Initialize the self-energy type
subroutine initialize_selfenergy(ls, qp, dr, nbasis, thirdorder, fourthorder, mw, mem)
    !> The selfenergy
    class(lo_selfenergy), intent(out) :: ls
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> The number of basis function for the selfenergy
    integer, intent(in) :: nbasis
    !> Do we compute third order ?
    logical, intent(in) :: thirdorder
    !> Do we compute fourth order ?
    logical, intent(in) :: fourthorder
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> To help decide on the maximum frequency
    real(r8), dimension(3) :: omega_max
    !> some buffers
    real(r8) :: f0, delta
    !> And some integers for the do loops
    integer :: n

    ! First, let's allocate everything
    ls%nbasis = nbasis
    allocate(ls%im_weight(ls%nbasis, dr%n_mode, qp%n_irr_point))
    allocate(ls%omega_n(ls%nbasis))
    ls%im_weight = 0.0_r8
    ls%omega_n = 0.0_r8

    ! Now, let's decide on a maximum frequency, this depends on the processes that we compute
    omega_max = -lo_huge
    omega_max(1) = dr%omega_max
    if (thirdorder) omega_max(2) = 2 * dr%omega_max
    if (fourthorder) omega_max(3) = 3 * dr%omega_max
    ! We add a bit for extra care, but shouldn't matter
    ls%omega_max = maxval(omega_max) + maxval(dr%default_smearing) * 3.0_r8

    ! Ok, now we can get our basis functions
    delta = ls%omega_max / real(ls%nbasis, r8)
    f0 = 0.0_r8
    do n=1, ls%nbasis
        f0 = f0 + delta
        ls%omega_n(n) = f0
    end do
    ! This allows to have a smooth fitting
    ls%width = 2.0_r8 * delta
    ! Done, we are good to go
end subroutine


!> Compute and fit the self-energy on a basis set of lorentzian
subroutine compute_selfenergy(ls, qp, dr, uc, fct, fcf, temperature, isotope, thirdorder, fourthorder, mw, mem)
    !> The selfenergy
    class(lo_selfenergy), intent(inout) :: ls
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    ! Third order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    ! Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> The temperature
    real(r8), intent(in) :: temperature
    ! Do we compute isotope ?
    logical, intent(in) :: isotope
    !> Do we compute third order ?
    logical, intent(in) :: thirdorder
    !> Do we compute fourth order ?
    logical, intent(in) :: fourthorder
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> Prefactor for the isotope scattering
    real(r8), parameter :: prefactor_iso = lo_pi / 4.0_r8
    !> Prefactor for the threephonon scattering
    real(r8), parameter :: prefactor_3ph = lo_pi / 16.0_r8
    !> Prefactor for the fourphonon scattering
    real(r8), parameter :: prefactor_4ph = lo_pi / 96.0_r8

    !> The matrices for the non-negative least-squares
    real(r8), dimension(:, :), allocatable :: amat, Q, A_inequal, A_equal
    !> The vectors for the non-negative least-squares
    real(r8), dimension(:), allocatable :: ymat, c, B_inequal, B_equal, sol
    ! Eigenvectors for the isotope scattering
    complex(r8), dimension(uc%na*3, 2) :: egviso
    ! Fourier transform of the matrix elements for 3ph
    complex(r8), dimension(dr%n_mode**3) :: ptf3
    ! Fourier transform of the matrix elements for 4ph
    complex(r8), dimension(dr%n_mode**4) :: ptf4
    ! Frequency scaled eigenvectors
    complex(r8), dimension(dr%n_mode) :: egv1, egv2, egv3, egv4
    ! Helper for Fourier transforms
    complex(r8), dimension(dr%n_mode**2) :: evp1
    ! Helper for Fourier transforms
    complex(r8), dimension(dr%n_mode**3) :: evp2
    ! Helper for Fourier transforms
    complex(r8), dimension(dr%n_mode**4) :: evp3
    ! The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4
    ! The q-point grid dimension
    integer, dimension(3) :: dims
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> For the harmonic values
    real(r8) :: om1, om2, om3, om4, n2, n3, n4
    !> For the scattering
    real(r8) :: psisq, sigma, sig1, sig2, sig3, sig4, prefactor, plf1, plf2
    !> Integers for the do loops
    integer :: n, m, q1, b1, q2, b2, q3, b3, q4, b4
    !> Integer for the parallelization
    integer :: ctr
    !> timing
    real(r8) :: t0

    ! grid dimensions
    select type(qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! Let's prepare the non-negative least-squares
    allocate(A_inequal(ls%nbasis, ls%nbasis))
    allocate(amat(ls%nbasis, ls%nbasis))
    allocate(Q(ls%nbasis, ls%nbasis))
    allocate(B_inequal(ls%nbasis))
    allocate(ymat(ls%nbasis))
    allocate(sol(ls%nbasis))
    allocate(c(ls%nbasis))
    allocate(A_equal(1, 1))
    allocate(B_equal(1))

    ! Set the matrices for the non-negative least-squares
    ! The inequal thing is what makes it non-negative
    ! TODO find some sum rules for the imaginary self energy for the equal constraints
    A_inequal = 0.0_r8
    B_inequal = 0.0_r8
    B_equal = 0.0_r8
    do n=1, ls%nbasis
        A_inequal(n, n) = 1
    end do

    ! Prefill the feature matrix for the non-negative least-squares
    do n=1, ls%nbasis
        do m=1, ls%nbasis
            amat(n, m) =  modified_lorentzian(ls%omega_n(n), &
                                              ls%omega_n(m), &
                                              ls%width)
        end do
    end do

    ! Little sanity check, should be useless
    ls%im_weight = 0.0_r8

    ! For the nice progressbar
    t0 = walltime()
    if (mw%talk) call lo_progressbar_init()

    ctr = 0
    do q1=1, qp%n_irr_point
    do b1=1, dr%n_mode
        om1 = dr%iq(q1)%omega(b1)
        if (om1 .lt. lo_freqtol) cycle

        ! Paralellize
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle

        ! Define some thing for the mode we are working on
        egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
        egviso(:, 1) = dr%iq(q1)%egv(:, b1)
        ymat = 0.0_r8
        ! First, we take care of the isotope
        if (isotope) then
            do q2=1, qp%n_full_point
            do b2=1, dr%n_mode
                om2 = dr%aq(q2)%omega(b2)
                if (om2 .lt. lo_freqtol) cycle

                egviso(:, 2) = dr%aq(q2)%egv(:, b2)
                psisq = isotope_scattering_strength(uc, egviso) * prefactor_iso * qp%ap(q2)%integration_weight
                sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)
                do n=1, ls%nbasis
                    ymat(n) = ymat(n) + psisq * lo_gauss(ls%omega_n(n), om2, sigma) * om1 * om2
                end do
            end do
            end do
        end if
        ! Now, let's go for the threephonon
        if (thirdorder) then
            do q2=1, qp%n_full_point
                q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)

                prefactor = prefactor_3ph * qp%ap(q2)%integration_weight
                qv2 = qp%ap(q2)%r
                qv3 = qp%ap(q3)%r
                call pretransform_phi3(fct, qv2, qv3, ptf3)
                do b2=1, dr%n_mode
                    om2 = dr%aq(q2)%omega(b2)
                    if (om2 .lt. lo_freqtol) cycle

                    n2 = lo_planck(temperature, om2)
                    egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
                    sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), &
                                             dr%default_smearing(b2), 1.0_r8)

                    evp1 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                    do b3=1, dr%n_mode
                        om3 = dr%aq(q3)%omega(b3)
                        if (om3 .lt. lo_freqtol) cycle

                        n3 = lo_planck(temperature, om3)
                        egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                        sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), &
                                                 dr%default_smearing(b3), 1.0_r8)
                        sigma = sqrt(sig2**2 + sig3**2)

                        evp2 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                        evp2 = conjg(evp2)
                        c0 = dot_product(evp2, ptf3)
                        psisq = abs(c0*conjg(c0)) * prefactor

                        plf1 = n2 + n3 + 1.0_r8
                        plf2 = n2 - n3
                        do n=1, ls%nbasis
                            ymat(n) = ymat(n) + psisq * plf1 * lo_gauss(ls%omega_n(n), om2 + om3, sigma)
                            ymat(n) = ymat(n) - psisq * plf1 * lo_gauss(ls%omega_n(n),-om2 - om3, sigma)
                            ymat(n) = ymat(n) + psisq * plf2 * lo_gauss(ls%omega_n(n),-om2 + om3, sigma)
                            ymat(n) = ymat(n) - psisq * plf2 * lo_gauss(ls%omega_n(n), om2 - om3, sigma)
                        end do
                    end do
                end do
            end do
        end if
        ! And finally, the fourphonon
        if (fourthorder) then
            do q2=1, qp%n_full_point
            do q3=1, qp%n_full_point
                q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)

                qv2 = qp%ap(q2)%r
                qv3 = qp%ap(q3)%r
                qv4 = qp%ap(q4)%r
                call pretransform_phi4(fcf, qv2, qv3, qv4, ptf4)

                prefactor = prefactor_4ph * qp%ap(q2)%integration_weight * qp%ap(q3)%integration_weight
                do b2=1, dr%n_mode
                    om2 = dr%aq(q2)%omega(b2)
                    if (om2 .lt. lo_freqtol) cycle

                    egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
                    n2 = lo_planck(temperature, om2)
                    sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), &
                                             dr%default_smearing(b2), 1.0_r8)

                    evp1 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, &
                               evp1, dr%n_mode)
                    do b3=1, dr%n_mode
                        om3 = dr%aq(q3)%omega(b3)
                        if (om3 .lt. lo_freqtol) cycle

                        egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                        n3 = lo_planck(temperature, om3)
                        sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), &
                                                 dr%default_smearing(b3), 1.0_r8)
                        evp2 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, &
                                   evp2, dr%n_mode)
                        do b4=1, dr%n_mode
                            om4 = dr%aq(q4)%omega(b4)
                            if (om4 .lt. lo_freqtol) cycle

                            n4 = lo_planck(temperature, om4)
                            egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)
                            sig4 = qp%adaptive_sigma(qp%ap(q4)%radius, dr%aq(q4)%vel(:, b4), &
                                                    dr%default_smearing(b4), 1.0_r8)
                            sigma = sqrt(sig2**2 + sig3**2 + sig4**2)

                            evp3 = 0.0_r8
                            call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, &
                                       evp2, 1, evp3, dr%n_mode)
                            evp3 = conjg(evp3)
                            c0 = dot_product(evp3, ptf4)
                            psisq = abs(c0*conjg(c0)) * prefactor

                            ! For fourphonon, we have four process to take into account, but only two prefactors
                            plf1 = (n2 + 1.0_r8) * (n3 + 1.0_r8) * (n4 + 1.0_r8) - n2 * n3 * n4
                            plf2 = 3.0_r8 * n2 * (n3 + 1.0_r8) * (n4 + 1.0_r8) - (n2 + 1.0_r8) * n3 * n4

                            do n=1, ls%nbasis
                                ymat(n) = ymat(n) + psisq * plf1 * &
                                    lo_gauss(ls%omega_n(n), om2 + om3 + om4, sigma)
                                ymat(n) = ymat(n) - psisq * plf1 * &
                                    lo_gauss(ls%omega_n(n), -om2 - om3 - om4, sigma)
                                ymat(n) = ymat(n) + psisq * plf2 * &
                                    lo_gauss(ls%omega_n(n), -om2 + om3 + om4, sigma)
                                ymat(n) = ymat(n) - psisq * plf2 * &
                                    lo_gauss(ls%omega_n(n), om2 - om3 - om4, sigma)
                            end do
                        end do
                    end do
                end do
            end do
            end do
        end if
        ! Now we can fit the imaginary part of the self-energy
        Q = 0.0_r8
        c = 0.0_r8
        ! First we compute Q=A^T * A
        call dsyrk('U', 'T', ls%nbasis, ls%nbasis, 1.0_r8, amat, ls%nbasis, 0.0_r8, Q, ls%nbasis)
        ! Then we compute c=A^T * y
        call lo_dgemv(transpose(amat), ymat, c)
        ! With this and the inequality constraints, the quadratic programing is equivalent to a non-negative lsq
        call lo_solve_quadratic_program(Q, c, sol, A_equal, A_inequal, &
                                        B_equal, B_inequal, 0, ls%nbasis, 0, 1e-13_r8)
        ! I want my zeros to be zeros
        do n=1, ls%nbasis
            if (sol(n) .lt. 1e-13_r8) sol(n) = 0.0_r8
        end do
        ls%im_weight(:, b1, q1) = sol(:)
        if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', ctr, &
                                         dr%n_mode * qp%n_irr_point, walltime() - t0)
    end do
    end do
    if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', dr%n_mode * qp%n_irr_point, &
                                     dr%n_mode * qp%n_irr_point, walltime() - t0)
    call mw%allreduce('sum', ls%im_weight)

    ! Let's fix the degeneracies
    fixdegen: block
        ! Buffer for the weights
        real(r8), dimension(ls%nbasis) :: buf
        integer :: q1, b1, j

        do q1=1, qp%n_irr_point
            do b1=1, dr%n_mode
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle

                buf = 0.0_r8
                do j=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    buf = buf + ls%im_weight(:, b2, q1)
                end do
                buf = buf / real(dr%iq(q1)%degeneracy(b1), r8)
                do j=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    ls%im_weight(:, b2, q1) = buf
                end do
            end do
        end do
    end block fixdegen
    ! And voila, we have our imaginary self-energy. Now we just have to find what to do with it
end subroutine


!> Write the self energy type to hdf5
subroutine write_selfenergy_to_hdf5(ls, filename)
    !> helper container
    class(lo_selfenergy), intent(in) :: ls
    !> filename
    character(len=*), intent(in) :: filename

    type(lo_hdf5_helper) :: h5
    integer :: nirr, nmode

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    nirr = size(ls%im_weight, 3)
    nmode = size(ls%im_weight, 2)

    ! Store metadata
    call h5%store_attribute(nirr, h5%file_id, 'n_irr_point')
    call h5%store_attribute(nmode, h5%file_id, 'n_mode')
    call h5%store_attribute(ls%nbasis, h5%file_id, 'nbasis')
    call h5%store_attribute(ls%width, h5%file_id, 'width')

    ! Store actual data
    call h5%store_data(ls%im_weight, h5%file_id, 'weight', enhet='unitless')
    call h5%store_data(ls%omega_n, h5%file_id, 'omega_n', enhet='atomic')

    ! Clost the file
    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> Clean destruction of the self energy type
subroutine destroy_selfenergy(ls)
    !> helper container
    class(lo_selfenergy), intent(inout) :: ls

    if (allocated(ls%omega_n)) deallocate(ls%omega_n)
    if (allocated(ls%im_weight)) deallocate(ls%im_weight)
    ls%nbasis = -lo_hugeint
    ls%width = -lo_huge
    ls%omega_max = -lo_huge
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

!> A modified lorentzian function, with the property L(0) = 0
function modified_lorentzian(x, mu, sigma) result (l)
    !> point to evaluate
    real(r8), intent(in) :: x
    !> mean
    real(r8), intent(in) :: mu
    !> Full width half maximum
    real(r8), intent(in) :: sigma
    !> value
    real(r8) :: l
    !> Inverse of pi, we don't want to calculate it everytime
    real(r8), parameter :: invpi=0.318309886183791_r8
    !> Some intermediate value
    real(r8) :: f0, f1, f2

    f0 = x * mu * sigma
    f1 = x**2 - mu**2
    f2 = sigma * x
    l= f0 / (f1**2 + f2**2) * invpi
end function
end module
