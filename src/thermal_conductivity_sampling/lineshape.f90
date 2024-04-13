#include "precompilerdefinitions"
module lineshape
use konstanter, only: r8, lo_freqtol, lo_twopi, lo_exitcode_param, lo_hugeint, lo_pi, lo_opi
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_gauss, lo_planck
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions

use new_scattering, only: lo_scattering_rates

implicit none

private
public :: initialize_lineshape
public :: compute_lineshape
contains


! Little container to simplify the computation of the lineshape
type lineshape
    ! The weight of each mode, dimension qp_irr and n_mode
    real(r8), dimension(:, :, :), allocatable :: weight_sf
    ! The frequencies and linewidth of the basis
    real(r8), dimension(:), allocatable :: omega_n
    ! The max frequency and the linewidth of each lorentzian basis
    real(r8) :: omega_max, gamma_n
    ! The number of basis function
    integer :: nbasis
end type


subroutine initialize_lineshape(qp, dr, ls, nbasis, omega_max, mw, mem)
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The lineshape
    type(lo_lineshape), intent(out) :: ls
    !> The number of basis function for the lineshape
    integer :: nbasis
    !> The max frequency for the basis
    real(r8) :: omega_max
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ls%nbasis = nbasis
    ls%omega_max = omega_max
    allocate(ls%omega_n(nbasis))
    delta = omega_max / real(nbasis, r8)
    ls%gamma_n = sqrt(real(nbasis, r8)) * delta
    allocate(ls%weight_sf(qp%n_irr_point, dr%n_mode, nbasis))

    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            do n=1, ls%nbasis
                ls%weight_sf(q1, b1, n) = lorentzian(ls%omega_n(n), dr%iq(q1)%omega(b1), dr%iq(q1)%linewidth(b1))
            end do
        end do
    end do
end subroutine


subroutine compute_lineshape(qp, dr, sr, ls, mw, mem)
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The scattering matrix elements
    type(lo_scattering_amplitude), intent(in) :: sr
    !> The lineshape
    type(lo_lineshape), intent(out) :: ls
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> Buffer for the imaginary part of the self energy
    real(r8), dimension(:, :, :), allocatable :: weight_imag
    !> Some buffers
    real(r8) :: f0, f1, f2, f3, f4, buf
    !> For the harmonic values
    real(r8) :: om1, om2, om3, n1, n2, n3
    !> Integers to help it make more readable
    integer :: i, q1, q2, q3, b1, b2, b3, n, m, k

    allocate(weight_imag(qp%n_irr_point, dr%n_mode, ls%nbasis))

    ! First we compute the imaginary part of the self energy
    do i=1, sr%nscatter_3ph
        q1 = sr%threephonon(i)%q1
        q2 = sr%threephonon(i)%q2
        q3 = sr%threephonon(i)%q3
        b1 = sr%threephonon(i)%b1
        b2 = sr%threephonon(i)%b2
        b3 = sr%threephonon(i)%b3
        do n=1, ls%nbasis
        do m=1, ls%nbasis
        do k=1, ls%nbasis
            f0 = lorentzian(ls%omega_n(n), ls%omega_n(m) + ls%omega_n(k), ls%gamma_n(m) + ls%gamma_n(k))
            f0 = f0 * (n2 + n3 + 1.0_r8)
            weight_imag(q1, b1, n) = ls%weight_imag(q1, b1, n) + f0 * sr(i)%psisq

            f0 = lorentzian(ls%omega_n(n), ls%omega_n(m) - ls%omega_n(k), ls%gamma_n(m) + ls%gamma_n(k))
            f0 = f0 * (n3 - n2)
            weight_imag(q1, b1, n) = ls%weight_imag(q1, b1, n) + f0 * sr(i)%psisq
        end do
        end do
        end do
    end do
    call mw%allreduce('sum', weight_imag)

    ! And we can get everything inside the spectral function
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            do n=1, ls%nbasis
            f0 = 0.0_r8
            f1 = 0.0_r8
                do m=1, ls%nbasis
                    buf = weight_imag(q1, b1, n) * lorentzian(ls%omega_n(n), ls%omega_n(m), ls%gamma_n)
                    f0 = f0 + buf
                    f1 = f1 + buf * (ls%omega_n(n) - ls%omega_n(m))
                end do
            f2 = 4 * dr%iq(q1)%omega(b1)**2 * f0 / lo_pi
            f3 = ls%omega_n(n)**2 - om1**2 - 2 * om1 * f1
            f4 = 4 * om1**2 * f0**2
            ls%weight_sf(q1, b1, n) = f2 / (f3 + f4)
            end do
        end do
    end do
end subroutine
end module lineshape
