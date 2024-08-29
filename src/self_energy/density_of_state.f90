#include "precompilerdefinitions"
module density_of_state
use konstanter, only: r8, lo_freqtol
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, &
                          lo_linspace, lo_gauss, lo_trapezoid_integration
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos

use selfenergy, only: lo_selfenergy

implicit none
private
public :: compute_density_of_state

contains

subroutine compute_density_of_state(qp, dr, uc, ls, nf, pd, mw, mem)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic properties on this mesh
    type(lo_phonon_dispersions), intent(in) :: dr
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The self-energy
    type(lo_selfenergy), intent(in) :: ls
    !> Number of frequencies on the axis
    integer, intent(in) :: nf
    !> phonon dos
    type(lo_phonon_dos), intent(out) :: pd
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    initdos: block
        pd%n_atom = uc%na
        pd%n_mode = dr%n_mode
        pd%n_dos_point = nf
        pd%dosmin = 0.0_r8
        pd%dosmax = ls%omega_max
        pd%integrationtype = -1
        pd%smearing_prefactor = 1.0_r8
        allocate (pd%omega(pd%n_dos_point))
        allocate (pd%dos(pd%n_dos_point))
        allocate (pd%pdos_site(pd%n_dos_point, pd%n_atom))
        allocate (pd%pdos_mode(pd%n_dos_point, pd%n_mode))
        call lo_linspace(pd%dosmin, pd%dosmax, pd%omega)
        pd%dos = 0.0_r8
        pd%pdos_site = 0.0_r8
        pd%pdos_mode = 0.0_r8
    end block initdos

    evaluate: block
        complex(r8), dimension(3) :: cv0
        real(r8), dimension(:), allocatable :: sf, smeared_sf, buf
        real(r8), dimension(:, :), allocatable :: buf_site, buf_mode
        real(r8) :: om1, siteproj, sigma, f0
        integer :: q1, b1, ctr, iat
        real(r8) :: t0

        allocate (sf(pd%n_dos_point))
        allocate (smeared_sf(pd%n_dos_point))
        allocate (buf(pd%n_dos_point))
        allocate (buf_site(pd%n_dos_point, uc%na))
        allocate (buf_mode(pd%n_dos_point, dr%n_mode))
        buf = 0.0_r8
        buf_site = 0.0_r8
        buf_mode = 0.0_r8

        ! For the nice progressbar
        t0 = walltime()
        if (mw%talk) call lo_progressbar_init()

        ctr = 0
        do q1 = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                om1 = dr%iq(q1)%omega(b1)
                if (om1 .lt. lo_freqtol) cycle

                ! MPI
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                call ls%evaluate_spectralfunction(q1, b1, pd%omega, sf)

                sigma = qp%adaptive_sigma(qp%ip(q1)%radius, dr%iq(q1)%vel(:, b1), dr%default_smearing(b1), 1.0_r8)
                ! Smear the spectral function
                call smear_spectralfunction(pd%omega, sf, sigma, smeared_sf)
                ! Remove the possible contribution at zero frequency that could have come from the smearing
                smeared_sf = max(smeared_sf - smeared_sf(1), 0.0_r8)
                ! Normalize the smeared spectral function
                f0 = lo_trapezoid_integration(pd%omega, smeared_sf)
                ! Add the integration weight
                sf = smeared_sf*qp%ip(q1)%integration_weight/f0
                ! Distribute in every things
                buf = buf + sf
                buf_mode(:, b1) = buf_mode(:, b1) + sf
                do iat = 1, uc%na
                    cv0 = dr%iq(q1)%egv((iat - 1)*3 + 1:iat*3, b1)
                    siteproj = abs(dot_product(cv0, conjg(cv0)))
                    buf_site(:, iat) = buf_site(:, iat) + sf*siteproj
                end do
            end do
            if (mw%talk) call lo_progressbar(' ... computing the dos', q1, qp%n_irr_point, walltime() - t0)
        end do
        call mw%allreduce('sum', buf)
        call mw%allreduce('sum', buf_mode)
        call mw%allreduce('sum', buf_site)
        pd%dos = buf
        pd%pdos_mode = buf_mode
        pd%pdos_site = buf_site

        deallocate (sf)
        deallocate (buf)
        deallocate (buf_site)
        deallocate (buf_mode)
    end block evaluate

contains
    subroutine smear_spectralfunction(eaxis, sf, sigma, smeared_sf)
        !> The energy axis
        real(r8), dimension(:), intent(in) :: eaxis
        !> The raw spectral function
        real(r8), dimension(:), intent(in) :: sf
        !> The smearing width
        real(r8), intent(in) :: sigma
        !> The smeared spectral function
        real(r8), dimension(:), intent(out) :: smeared_sf

        !> some buffers
        real(r8) :: dx
        !> Some integer for loops
        integer :: n, nstep, i, j, a1, a2

        n = size(eaxis, 1)
        dx = eaxis(2) - eaxis(1)
        ! sigma = 1e-4_r8
        nstep = ceiling(6.0_r8*sigma/dx)

        smeared_sf = 0.0_r8
        do i = 1, n
            a1 = max(1, i - nstep)
            a2 = min(n, i + nstep)
            do j = a1, a2
                smeared_sf(i) = smeared_sf(i) + sf(j)*lo_gauss(eaxis(i), eaxis(j), sigma)
            end do
        end do

    end subroutine
end subroutine
end module
