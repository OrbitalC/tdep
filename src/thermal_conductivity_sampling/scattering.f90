#include "precompilerdefinitions"
module scattering
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol, lo_twopi, &
                      lo_exitcode_param, lo_pi
use gottochblandat, only: lo_planck, walltime, lo_gauss, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_lorentz
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
public :: compute_scattering

contains
subroutine compute_scattering(qp, dr, uc, fct, fcf, opts, mw, mem)
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

    !> The random number generator
    type(lo_mersennetwister) :: rng
    ! To initialize the random number generator
    real(r8) :: rseed
    ! Some integers
    integer :: q1, b1, ctr

    ! We start by initializing the random number generator, we need the same on each mpi rank
    rseed = walltime()
    call mw%bcast(rseed, from=mw%n - 1)
    call rng%init(iseed=0, rseed=rseed)

    if (mw%talk) call lo_progressbar_init()
    ctr = 0
    do q1 = 1, qp%n_irr_point
        do b1=1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            if (opts%isotopescattering) then
                call isotope_scattering(q1, b1, qp, dr, uc, opts%temperature, opts%sigma, opts%thres, mw, mem)
            end if
            if (opts%thirdorder) then
                call threephonon_scattering(q1, b1, qp, dr, fct, opts%temperature, opts%nsample3ph, &
                                            opts%sigma, opts%thres, rng, mw, mem)
            end if
            if (opts%fourthorder) then
                call fourphonon_scattering(q1, b1, qp, dr, fcf, opts%temperature, opts%nsample4ph, &
                                           opts%sigma, opts%thres, rng, mw, mem)
            end if
            ctr = ctr + 1
        end do
        if (mw%talk) call lo_progressbar(' ... computing scattering', ctr, dr%n_mode*qp%n_irr_point)
    end do
    if (mw%talk) call lo_progressbar(' ... computing scattering', dr%n_mode*qp%n_irr_point, dr%n_mode*qp%n_irr_point)

    degeneracy: block
        ! Buffers
        real(r8) :: buf_i, buf_p, buf_m, buf_pp, buf_pm, buf_mm
        ! Integer for the loops
        integer :: i, b2
    ! Now let's fix degeneracy
    if (opts%correctionlevel .gt. 2) then
        do q1 = 1, qp%n_irr_point
            do b1=1, dr%n_mode
                buf_i = 0.0_r8
                buf_p = 0.0_r8
                buf_m = 0.0_r8
                if (opts%fourthorder) then
                    buf_pp = 0.0_r8
                    buf_pm = 0.0_r8
                    buf_mm = 0.0_r8
                end if
                do i = 1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(i, b1)
                    buf_i = buf_i + dr%iq(q1)%p_iso(b2)
                    buf_p = buf_p + dr%iq(q1)%p_plus(b2)
                    buf_m = buf_m + dr%iq(q1)%p_minus(b2)
                    if (opts%fourthorder) then
                        buf_pp = buf_pp + dr%iq(q1)%p_plusplus(b2)
                        buf_pm = buf_pm + dr%iq(q1)%p_plusminus(b2)
                        buf_mm = buf_mm + dr%iq(q1)%p_minusminus(b2)
                    end if
                end do
                buf_i = buf_i / real(dr%iq(q1)%degeneracy(b1), r8)
                buf_p = buf_p / real(dr%iq(q1)%degeneracy(b1), r8)
                buf_m = buf_m / real(dr%iq(q1)%degeneracy(b1), r8)
                if (opts%fourthorder) then
                    buf_pp = buf_pp / real(dr%iq(q1)%degeneracy(b1), r8)
                    buf_pm = buf_pm / real(dr%iq(q1)%degeneracy(b1), r8)
                    buf_mm = buf_mm / real(dr%iq(q1)%degeneracy(b1), r8)
                end if
                do i=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(i, b1)
                    dr%iq(q1)%p_iso(b2) = buf_i
                    dr%iq(q1)%p_plus(b2) = buf_p
                    dr%iq(q1)%p_minus(b2) = buf_m
                    if (opts%fourthorder) then
                        dr%iq(q1)%p_plusplus(b2) = buf_pp
                        dr%iq(q1)%p_plusminus(b2) = buf_pm
                        dr%iq(q1)%p_minusminus(b2) = buf_mm
                    end if
                end do
            end do
        end do
    end if
    end block degeneracy
end subroutine

subroutine isotope_scattering(q1, b1, qp, dr, uc, temperature, smearing_prefactor, thres, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: q1, b1
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    ! The temperature
    real(r8), intent(in) :: temperature
    ! The smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    ! The gaussian threshold
    real(r8), intent(in) :: thres
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Isotope prefactor
    real(r8), parameter :: isotope_prefactor = lo_pi / 2.0_r8
    ! Eigenvectors
    complex(r8), dimension(uc%na*3, 2) :: egviso
    ! buffer
    real(r8) :: iso_scatter, scatterstrength, deltafunction, sigma, sig1, sig2
    ! The frequencies and bose-einstein distribution
    real(r8) :: om1, om2, n1, n2
    ! Integers for do loops
    integer :: q2, b2, ctr, n_count

    ctr = 0
    om1 = dr%iq(q1)%omega(b1)
    n1 = lo_planck(temperature, om1)
    iso_scatter = 0.0_r8
    do q2 = 1, qp%n_full_point
        do b2 = 1, dr%n_mode
            ! MPI
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle

            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            ! sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), smearing_prefactor)
            sig1 = dr%iq(q1)%linewidth(b1)
            sig2 = dr%iq(qp%ap(q2)%irreducible_index)%linewidth(b2)
            ! sigma = sqrt(sig1**2 + sig2**2)
            sigma = 2.0 * (sig1 + sig2)

            if (abs(om1 - om2) .lt. thres*sigma) then
                ! deltafunction = lo_gauss(om1, om2, sigma)
                deltafunction = lo_lorentz(om1, om2, sigma)
                egviso(:, 1) = dr%iq(q1)%egv(:, b1)
                egviso(:, 2) = dr%aq(q2)%egv(:, b2)
                n2 = lo_planck(temperature, om2)
                scatterstrength = isotope_scattering_strength(uc, egviso) * om1 * om2 * n1 * (n2 + 1.0_r8)
                iso_scatter = iso_scatter + scatterstrength * isotope_prefactor * deltafunction * &
                                            qp%ap(q2)%integration_weight
            end if

        end do
    end do
    call mw%allreduce('sum', iso_scatter)
    dr%iq(q1)%p_iso(b1) = iso_scatter
end subroutine

subroutine threephonon_scattering(q1, b1, qp, dr, fct, temperature, nsample3ph, &
        smearing_prefactor, thres, rng, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: q1, b1
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Fourth order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    ! The temperature
    real(r8), intent(in) :: temperature
    ! The number of 3ph scattering event actually sampled
    integer, intent(in) :: nsample3ph
    ! The smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    ! The gaussian threshold
    real(r8), intent(in) :: thres
    ! Random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: threephonon_prefactor = lo_pi / 4.0_r8
    ! The scattering buffers
    real(r8) :: bufp, bufm
    ! The number of scattering
    integer :: np_tot, np_sample, nm_tot, nm_sample
    ! Do we compute the scattering ?
    logical :: is_scatter
    ! Fourier transform of the matrix elements
    complex(r8), dimension(dr%n_mode**3) :: ptf
    ! Frequency scaled eigenvectors
    complex(r8), dimension(dr%n_mode) :: egv1, egv2, egv3
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**2) :: evp1
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**3) :: evp2
    ! The full qpoint grid
    integer, dimension(qp%n_full_point) :: qgridfull
    ! The frequencies, qpoints and velocities
    ! The grid dimensions
    integer, dimension(3) :: dims
    ! The qpoints and velocities
    real(r8), dimension(3) :: qv2, qv3, vel1, vel2, vel3
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, n1, n2, n3, psisq, plf, deltafunction, prefactor
    ! gaussian width for the delta approximation
    real(r8) :: sigma, sig1, sig2, sig3
    ! scattering strength before absolute value
    complex(r8) :: c0
    ! Integers for do loops
    integer :: ctr, i, j, q2, q3, b2, b3

    ! grid dimensions
    select type(qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select
    np_tot = 0
    np_sample = 0
    bufp = 0.0_r8
    nm_tot = 0
    nm_sample = 0
    bufm = 0.0_r8

    ! Create the array with the full grid and shuffle them
    qgridfull = (/ (i, i=1, qp%n_full_point) /)
    call rng%shuffle_int_array(qgridfull)

    ! already set some values
    om1 = dr%iq(q1)%omega(b1)
    vel1 = dr%iq(q1)%vel(:, b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    n1 = lo_planck(temperature, om1)
    sig1 = dr%iq(q1)%linewidth(b1)

    ctr = 0
    do i=1, qp%n_full_point
        q2 = qgridfull(i)
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)
        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        prefactor = threephonon_prefactor * qp%ap(q2)%integration_weight
        call pretransform_phi3(fct, qv2, qv3, ptf)
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            n2 = lo_planck(temperature, om2)
            vel2 = dr%aq(q2)%vel(:, b2)
            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = dr%iq(qp%ap(q2)%irreducible_index)%linewidth(b2)

            do b3=1, dr%n_mode
                ! MPI
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                n3 = lo_planck(temperature, om3)
                vel3 = dr%aq(q3)%vel(:, b3)
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                sig3 = dr%iq(qp%ap(q3)%irreducible_index)%linewidth(b3)

                is_scatter = .false.
                if (np_sample .lt. nsample3ph / mw%n / 2) is_scatter = .true.

                sigma = 2.0_r8 * (sig1 + sig2 + sig3)

                np_tot = np_tot + 1
                nm_tot = nm_tot + 1

                if (is_scatter) then
                    evp1 = 0.0_r8
                    evp2 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                    call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                    evp2 = conjg(evp2)
                    c0 = dot_product(evp2, ptf)
                    psisq = abs(c0*conjg(c0))

                    np_sample = np_sample + 1
                    plf = n1 * n2 * (n3 + 1.0_r8) * prefactor
                    deltafunction = lo_lorentz(om1, -om2 + om3, sigma)
                    bufp = bufp + deltafunction * psisq * plf

                    nm_sample = nm_sample + 1
                    plf = n1 * (n2 + 1.0_r8) * (n3 + 1.0_r8) * prefactor
                    deltafunction = lo_lorentz(om1, om2 + om3, sigma)
                    bufm = bufm + deltafunction * psisq * plf
                end if
            end do
        end do
    end do
    call mw%allreduce('sum', np_sample)
    call mw%allreduce('sum', nm_sample)
    call mw%allreduce('sum', np_tot)
    call mw%allreduce('sum', nm_tot)
    call mw%allreduce('sum', bufp)
    call mw%allreduce('sum', bufm)
    bufp = bufp * real(np_tot, r8) / real(np_sample, r8)
    bufm = bufm * real(nm_tot, r8) / real(nm_sample, r8)
    dr%iq(q1)%p_plus(b1) = bufp
    dr%iq(q1)%p_minus(b1) = bufm
end subroutine


subroutine fourphonon_scattering(q1, b1, qp, dr, fcf, temperature, nsample4ph, &
        smearing_prefactor, thres, rng, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: q1, b1
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    ! The temperature
    real(r8), intent(in) :: temperature
    ! The number of 4ph scattering event actually sampled
    integer, intent(in) :: nsample4ph
    ! The smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    ! The gaussian threshold
    real(r8), intent(in) :: thres
    ! Random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Four phonon prefactor
    real(r8), parameter :: fourphonon_prefactor = lo_pi / 8.0_r8
    ! The scattering buffers
    real(r8) :: bufpp, bufpm, bufmm
    ! The number of scattering
    integer :: npp_tot, npp_sample, npm_tot, npm_sample, nmm_tot, nmm_sample
    ! Do we compute the scattering ?
    logical :: is_scatter
    ! Fourier transform of the matrix elements
    complex(r8), dimension(dr%n_mode**4) :: ptf
    ! Eigenvectors
    complex(r8), dimension(dr%n_mode) :: egv1, egv2, egv3, egv4
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**2) :: evp1
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**3) :: evp2
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(dr%n_mode**4) :: evp3
    ! The full qpoint grid
    integer, dimension(qp%n_full_point) :: qgridfull1, qgridfull2
    ! The grid dimensions
    integer, dimension(3) :: dims
    ! The qpoints and velocities
    real(r8), dimension(3) :: qv2, qv3, qv4, vel1, vel2, vel3, vel4
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, n1, n2, n3, n4, psisq, plf, deltafunction, prefactor
    ! gaussian width for the delta approximation
    real(r8) :: sigma, sig1, sig2, sig3, sig4
    ! scattering strength before absolute value
    complex(r8) :: c0
    ! Integers for do loops
    integer :: ctr, i, j, q2, q3, q4, b2, b3, b4

    ! grid dimensions
    select type(qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select
    npp_tot = 0
    npp_sample = 0
    bufpp = 0.0_r8
    npm_tot = 0
    npm_sample = 0
    bufpm = 0.0_r8
    nmm_tot = 0
    nmm_sample = 0
    bufmm = 0.0_r8

    ! Create the array with the full grid and shuffle them
    qgridfull1 = (/ (i, i=1, qp%n_full_point) /)
    qgridfull2 = (/ (i, i=1, qp%n_full_point) /)
    call rng%shuffle_int_array(qgridfull1)
    call rng%shuffle_int_array(qgridfull2)

    ! already set some values
    om1 = dr%iq(q1)%omega(b1)
    vel1 = dr%iq(q1)%vel(:, b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    n1 = lo_planck(temperature, om1)
    sig1 = dr%iq(q1)%linewidth(b1)

    ctr = 0
    do i=1, qp%n_full_point
    q2 = qgridfull1(i)
    qv2 = qp%ap(q2)%r
    do j=1, qp%n_full_point
        q3 = qgridfull2(j)
        qv3 = qp%ap(q3)%r
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
        qv4 = qp%ap(q4)%r

        prefactor = fourphonon_prefactor * qp%ap(q2)%integration_weight * qp%ap(q3)%integration_weight
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            n2 = lo_planck(temperature, om2)
            vel2 = dr%aq(q2)%vel(:, b2)
            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = dr%iq(qp%ap(q2)%irreducible_index)%linewidth(b2)

            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                n3 = lo_planck(temperature, om3)
                vel3 = dr%aq(q3)%vel(:, b3)
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                sig3 = dr%iq(qp%ap(q3)%irreducible_index)%linewidth(b3)

                do b4=1, dr%n_mode
                    ctr = ctr + 1
                    if (mod(ctr, mw%n) .ne. mw%r) cycle

                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle
                    n4 = lo_planck(temperature, om4)
                    vel4 = dr%aq(q4)%vel(:, b4)
                    egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)
                    sig4 = dr%iq(qp%ap(q4)%irreducible_index)%linewidth(b4)

                    is_scatter = .false.
                    if (npp_sample .lt. nsample4ph / mw%n / 3) is_scatter = .true.

                    sigma = 2.0_r8 * (sig1 + sig2 + sig3 + sig4)

                    npp_tot = npp_tot + 1
                    npm_tot = npm_tot + 1
                    nmm_tot = nmm_tot + 1

                    if (is_scatter) then
                        evp1 = 0.0_r8
                        evp2 = 0.0_r8
                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode,    (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                        call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                        call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                        evp3 = conjg(evp3)
                        c0 = dot_product(evp3, ptf)
                        psisq = abs(c0*conjg(c0))

                        npp_sample = npp_sample + 1
                        plf = n1 * n2 * n3 * (n4 + 1.0_r8) * prefactor
                        deltafunction = lo_gauss(om1, -om2 - om3 + om4, sigma)
                        bufpp = bufpp + deltafunction*psisq*plf

                        npm_sample = npm_sample + 1
                        plf = n1 * n2 * (n3 + 1.0_r8) * (n4 + 1.0_r8) * prefactor
                        deltafunction = lo_gauss(om1, -om2 + om3 + om4, sigma)
                        bufpm = bufpm + deltafunction*psisq*plf

                        nmm_sample = nmm_sample + 1
                        plf = n1 * (n2 + 1.0_r8) * (n3 + 1.0_r8) * (n4 + 1.0_r8) * prefactor
                        deltafunction = lo_gauss(om1, om2 + om3 + om4, sigma)
                        bufmm = bufmm + deltafunction*psisq*plf
                    end if
                end do
            end do
        end do
    end do
    end do
    call mw%allreduce('sum', npp_sample)
    call mw%allreduce('sum', npm_sample)
    call mw%allreduce('sum', nmm_sample)
    call mw%allreduce('sum', npp_tot)
    call mw%allreduce('sum', npm_tot)
    call mw%allreduce('sum', nmm_tot)
    call mw%allreduce('sum', bufpp)
    call mw%allreduce('sum', bufpm)
    call mw%allreduce('sum', bufmm)
    bufpp = bufpp * real(npp_tot, r8) / real(npp_sample, r8)
    bufpm = bufpm * real(npm_tot, r8) / real(npm_sample, r8)
    bufmm = bufmm * real(nmm_tot, r8) / real(nmm_sample, r8)
    dr%iq(q1)%p_plusplus(b1) = bufpp
    dr%iq(q1)%p_plusminus(b1) = bufpm
    dr%iq(q1)%p_minusminus(b1) = bufmm
end subroutine


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
