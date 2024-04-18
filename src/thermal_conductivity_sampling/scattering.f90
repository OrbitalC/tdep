#include "precompilerdefinitions"
module new_scattering
use konstanter, only: r8, i8, lo_freqtol, lo_twopi, lo_exitcode_param, lo_hugeint, lo_pi, lo_twopi
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_gauss, lo_planck
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

type lo_psisq
    !> The actual psisq
    real(r8), dimension(:), allocatable :: psisq
    !> The indices for the scattering
    integer, dimension(:), allocatable :: q2, q3, q4, b2, b3, b4
    !> The number of scattering
    integer :: n
end type

! Container for scattering rates
type lo_scattering_rates
    !> The number of qpoint/mode on this rank
    integer :: nlocal_point
    !> The list of qpoint and modes for this rank
    integer, dimension(:), allocatable :: q1, b1
    !> The iso phonon scattering
    type(lo_psisq), dimension(:), allocatable :: iso
    !> The three phonon scattering
    type(lo_psisq), dimension(:), allocatable :: threephonon
    !> The four phonon scattering
    type(lo_psisq), dimension(:), allocatable :: fourphonon

    contains
        !> destroy the scattering amplitues
        procedure :: destroy
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
    integer :: q1, b1, q1f, il, j, k, nlocal_point, ctr, niso, nplus, nminus, npp, npm, nmm
    ! To do some counting
    integer :: bufiso, bufplus, bufminus, bufpp, bufpm, bufmm

    ! grid dimensions
    select type(qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    bufiso = 0
    bufplus = 0
    bufminus = 0
    bufpp = 0
    bufpm = 0
    bufmm = 0

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

    ! We can allocate all we need
    sr%nlocal_point = nlocal_point
    allocate(sr%q1(nlocal_point))
    allocate(sr%b1(nlocal_point))
    if (opts%isotopescattering) allocate(sr%iso(nlocal_point))
    if (opts%thirdorder) allocate(sr%threephonon(nlocal_point))
    if (opts%fourthorder) allocate(sr%fourphonon(nlocal_point))

    ! Let's attribute the q1/b1 indices to the rank
    il = 0
    ctr = 0
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            ! We skip the acoustic mode at Gamma
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            ctr = ctr + 1

            ! MPI thing
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            il = il + 1
            sr%q1(il) = q1
            sr%b1(il) = b1
        end do
    end do

    t0 = walltime()
    if (mw%talk) call lo_progressbar_init()
    do il=1, sr%nlocal_point
        if (opts%isotopescattering) then
            call compute_isotope_scattering(il, sr, qp, dr, uc, opts%temperature, mw, mem)
        end if
        if (opts%thirdorder) then
            call compute_threephonon_scattering(il, sr, qp, dr, fct, dims, opts%temperature, mw, mem)
        end if
        if (opts%fourthorder) then
            call compute_fourphonon_scattering(il, sr, qp, dr, fcf, dims, opts%temperature, mw, mem)
        end if
        if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', il, sr%nlocal_point, walltime() - t0)
    end do
    if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', sr%nlocal_point, sr%nlocal_point, walltime() - t0)


    diagnostic: block
        integer(i8) :: bufiso, buf3ph, buf4ph
        integer(i8) :: totiso, tot3ph, tot4ph
        real(r8) :: ratioiso, ratio3ph, ratio4ph

        bufiso = 0
        buf3ph = 0
        buf4ph = 0
        do il=1, sr%nlocal_point
            if (opts%isotopescattering) bufiso = bufiso + sr%iso(il)%n
            if (opts%thirdorder) buf3ph = buf3ph + sr%threephonon(il)%n
            if (opts%fourthorder) buf4ph = buf4ph + sr%fourphonon(il)%n
        end do
        if (opts%isotopescattering) call mw%allreduce('sum', bufiso)
        if (opts%thirdorder) call mw%allreduce('sum', buf3ph)
        if (opts%fourthorder) call mw%allreduce('sum', buf4ph)

        totiso = qp%n_irr_point * qp%n_full_point * dr%n_mode**2
        tot3ph = qp%n_irr_point * qp%n_full_point * dr%n_mode**3
        tot4ph = qp%n_irr_point * qp%n_full_point**2 * dr%n_mode**4

        ratioiso = 100.0_r8 * real(bufiso, r8) / real(totiso, r8)
        ratio3ph = 100.0_r8 * real(buf3ph, r8) / real(tot3ph, r8)
        ratio4ph = 100.0_r8 * real(buf4ph, r8) / real(tot4ph, r8)

        if (mw%talk) then
            write(*, *) ''
            write(*, *) '     number of isotope scattering matrix elements', bufiso
            write(*, *) '                                            % iso', ratioiso
            write(*, *) ' number of threephonon scattering matrix elements', buf3ph
            write(*, *) '                                            % 3ph', ratio3ph
            write(*, *) '  number of fourphonon scattering matrix elements', buf4ph
            write(*, *) '                                            % 4ph', ratio4ph
        end if
    end block diagnostic
end subroutine


subroutine compute_isotope_scattering(il, sr, qp, dr, uc, temperature, mw, mem)
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
    !> The temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Isotope prefactor
    real(r8), parameter :: isotope_prefactor = lo_pi * 0.5_r8
    ! Eigenvectors
    complex(r8), dimension(uc%na*3, 2) :: egviso
    ! prefactor and phonon buffers
    real(r8) :: om1, om2, sig1, sig2, sigma, psisq, prefactor
    ! Integers for do loops
    integer :: q1, b1, q2, b2, i, niso

    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)

    niso = 0
    do q2=1, qp%n_full_point
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)
            if (abs(om1 - om2) .lt. 4.0_r8 * sigma) niso = niso + 1
        end do
    end do

    sr%iso(il)%n = niso
    allocate(sr%iso(il)%psisq(niso))
    allocate(sr%iso(il)%q2(niso))
    allocate(sr%iso(il)%b2(niso))

    om1 = dr%iq(q1)%omega(b1)
    egviso(:, 1) = dr%iq(q1)%egv(:, b1)

    i = 0
    do q2=1, qp%n_full_point
        prefactor = isotope_prefactor * qp%ap(q2)%integration_weight
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)
            if (abs(om1 - om2) .lt. 4.0_r8 * sigma) then
                i = i + 1

                egviso(:, 2) = dr%aq(q2)%egv(:, b2)
                psisq = isotope_scattering_strength(uc, egviso) * prefactor

                sr%iso(il)%q2(i) = q2
                sr%iso(il)%b2(i) = b2
                sr%iso(il)%psisq(i) = psisq
            end if
        end do
    end do
end subroutine


subroutine compute_threephonon_scattering(il, sr, qp, dr, fct, dims, temperature, mw, mem)
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
    ! The dimension of the grid
    integer, dimension(3), intent(in) :: dims
    ! The temperature
    real(r8), intent(in) :: temperature
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: threephonon_prefactor = lo_pi * 0.25_r8
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
    ! The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sigma
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, n1, n2, n3, plf, psisq, prefactor
    !
    complex(r8) :: c0
    ! Integers for do loops
    integer :: i, iq1, q1, q2, q3, b1, b2, b3, n3ph

    ! Already set some values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    n1 = lo_planck(temperature, om1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    sig1 = qp%adaptive_sigma(qp%ip(q1)%radius, dr%iq(q1)%vel(:, b1), &
                             dr%default_smearing(b1), 1.0_r8)

    n3ph = 0
    do q2=1, qp%n_full_point
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)

        if (q3 .lt. q2) cycle

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle
            sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), &
                                        dr%default_smearing(b2), 1.0_r8)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), &
                                        dr%default_smearing(b3), 1.0_r8)
                sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                if (abs(om1 + om2 - om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 + om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 - om3) .lt. 4.0_r8 * sigma) n3ph = n3ph + 1
            end do
        end do
    end do

    sr%threephonon(il)%n = n3ph
    allocate(sr%threephonon(il)%psisq(n3ph))
    allocate(sr%threephonon(il)%q2(n3ph))
    allocate(sr%threephonon(il)%q3(n3ph))
    allocate(sr%threephonon(il)%b2(n3ph))
    allocate(sr%threephonon(il)%b3(n3ph))

    i = 0
    do q2=1, qp%n_full_point
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)

        if (q3 .lt. q2) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        call pretransform_phi3(fct, qv2, qv3, ptf)
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            n2 = lo_planck(temperature, om2)
            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), &
                                     dr%default_smearing(b2), 1.0_r8)
            prefactor = threephonon_prefactor * qp%ap(q2)%integration_weight

            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                n3 = lo_planck(temperature, om3)
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)

                sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), &
                                         dr%default_smearing(b3), 1.0_r8)
                sigma = sqrt(sig1**2 + sig2**2 + sig3**2)

                ! Do we need to compute the scattering ?
                if (abs(om1 + om2 - om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 + om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 - om3) .lt. 4.0_r8 * sigma) then
                    i = i + 1

                    evp1 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                    evp2 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                    evp2 = conjg(evp2)
                    c0 = dot_product(evp2, ptf)
                    psisq = abs(c0*conjg(c0)) * prefactor

                    sr%threephonon(il)%psisq(i) = psisq
                    sr%threephonon(il)%q2(i) = q2
                    sr%threephonon(il)%q3(i) = q3
                    sr%threephonon(il)%b2(i) = b2
                    sr%threephonon(il)%b3(i) = b3
                end if
            end do
        end do
    end do
end subroutine


subroutine compute_fourphonon_scattering(il, sr, qp, dr, fcf, dims, temperature, mw, mem)
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
    ! The dimension of the grid
    integer, dimension(3), intent(in) :: dims
    ! The temperature
    real(r8), intent(in) :: temperature
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
    !> The complex scattering amplitude
    complex(r8) :: c0
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, n1, n2, n3, n4, psisq, prefactor
    ! The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sig4, sigma, plf
    ! Integers for do loops
    integer :: i, q1, q2, q3, q4, b1, b2, b3, b4
    ! Integers for counting
    integer :: n4ph

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    n1 = lo_planck(temperature, om1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    sig1 = qp%adaptive_sigma(qp%ip(q1)%radius, dr%iq(q1)%vel(:, b1), &
                             dr%default_smearing(b1), 1.0_r8)

    n4ph = 0
    do q2=1, qp%n_full_point
    do q3=q2, qp%n_full_point
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)

        if (q4 .lt. q3) cycle

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), &
                                     dr%default_smearing(b2), 1.0_r8)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), &
                                         dr%default_smearing(b3), 1.0_r8)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    sig4 = qp%adaptive_sigma(qp%ap(q4)%radius, dr%aq(q4)%vel(:, b4), &
                                             dr%default_smearing(b4), 1.0_r8)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)
                    if (abs(om1 + om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 + om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 + om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 + om3 + om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. 4.0_r8 * sigma) n4ph = n4ph + 1
                end do
            end do
        end do
    end do
    end do

    sr%fourphonon(il)%n = n4ph
    allocate(sr%fourphonon(il)%psisq(n4ph))
    allocate(sr%fourphonon(il)%q2(n4ph))
    allocate(sr%fourphonon(il)%q3(n4ph))
    allocate(sr%fourphonon(il)%q4(n4ph))
    allocate(sr%fourphonon(il)%b2(n4ph))
    allocate(sr%fourphonon(il)%b3(n4ph))
    allocate(sr%fourphonon(il)%b4(n4ph))

    i = 0
    do q2=1, qp%n_full_point
    do q3=q2, qp%n_full_point
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
        if (q4 .lt. q3) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        prefactor = fourphonon_prefactor * qp%ap(q2)%integration_weight * qp%ap(q3)%integration_weight
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            n2 = lo_planck(temperature, om2)
            sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), &
                                     dr%default_smearing(b2), 1.0_r8)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                n3 = lo_planck(temperature, om3)
                sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), &
                                            dr%default_smearing(b3), 1.0_r8)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    sig4 = qp%adaptive_sigma(qp%ap(q4)%radius, dr%aq(q4)%vel(:, b4), &
                                             dr%default_smearing(b4), 1.0_r8)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)

                    if (abs(om1 + om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 + om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 + om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 + om3 + om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. 4.0_r8 * sigma) then
                        i = i + 1

                        egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)
                        n4 = lo_planck(temperature, om4)

                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                        evp3 = conjg(evp3)
                        c0 = dot_product(evp3, ptf)
                        psisq = abs(c0*conjg(c0)) * prefactor

                        sr%fourphonon(il)%psisq(i) = psisq
                        sr%fourphonon(il)%q2(i) = q2
                        sr%fourphonon(il)%q3(i) = q3
                        sr%fourphonon(il)%q4(i) = q4
                        sr%fourphonon(il)%b2(i) = b2
                        sr%fourphonon(il)%b3(i) = b3
                        sr%fourphonon(il)%b4(i) = b4
                    end if
                end do
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


subroutine destroy(sr)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr

    integer :: il

    do il=1, sr%nlocal_point
        if (allocated(sr%iso)) then
            deallocate(sr%iso(il)%psisq)
            deallocate(sr%iso(il)%q2)
            deallocate(sr%iso(il)%b2)
        end if
        if (allocated(sr%threephonon)) then
            deallocate(sr%threephonon(il)%psisq)
            deallocate(sr%threephonon(il)%q2)
            deallocate(sr%threephonon(il)%q3)
            deallocate(sr%threephonon(il)%b2)
            deallocate(sr%threephonon(il)%b3)
        end if
        if (allocated(sr%fourphonon)) then
            deallocate(sr%fourphonon(il)%psisq)
            deallocate(sr%fourphonon(il)%q2)
            deallocate(sr%fourphonon(il)%q3)
            deallocate(sr%fourphonon(il)%q4)
            deallocate(sr%fourphonon(il)%b2)
            deallocate(sr%fourphonon(il)%b3)
            deallocate(sr%fourphonon(il)%b4)
        end if
    end do
    if (allocated(sr%iso)) deallocate(sr%iso)
    if (allocated(sr%threephonon)) deallocate(sr%threephonon)
    if (allocated(sr%fourphonon)) deallocate(sr%fourphonon)
    if (allocated(sr%q1)) deallocate(sr%q1)
    if (allocated(sr%b1)) deallocate(sr%b1)
    sr%nlocal_point = -lo_hugeint
end subroutine


end module
