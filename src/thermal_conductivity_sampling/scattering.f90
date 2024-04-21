#include "precompilerdefinitions"
module new_scattering
use konstanter, only: r8, i8, lo_freqtol, lo_twopi, lo_exitcode_param, lo_hugeint, lo_pi, lo_twopi, lo_tol
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

use type_symmetryoperation, only: lo_operate_on_vector
use lo_sorting, only: lo_qsort

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
    !> Let's precompute the Bose-Einstein distribution and adaptive smearing
    real(r8), dimension(:, :), allocatable :: be, sigma_q
    !> The ratio for maximum likelihood estimation
    real(r8) :: mle_4ph

    contains
        !> destroy the scattering amplitues
        procedure :: destroy => sr_destroy
        !> Measure size in memory, in bytes
        procedure :: size_in_mem => sr_size_in_mem
        !> Measure size in memory for isotopes, in bytes
        procedure :: size_in_mem_iso
        !> Measure size in memory for threephonon, in bytes
        procedure :: size_in_mem_3ph
        !> Measure size in memory for fourphonon, in bytes
        procedure :: size_in_mem_4ph
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

    ! Initialize the random number generator
    call rng%init(iseed=mw%r, rseed=walltime())

    bufiso = 0
    bufplus = 0
    bufminus = 0
    bufpp = 0
    bufpm = 0
    bufmm = 0

    ! We can start some precomputation
    allocate(sr%be(qp%n_irr_point, dr%n_mode))
    allocate(sr%sigma_q(qp%n_irr_point, dr%n_mode))
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            sr%be(q1, b1) = lo_planck(opts%temperature, dr%iq(q1)%omega(b1))
            sr%sigma_q(q1, b1) = qp%adaptive_sigma(qp%ip(q1)%radius, dr%iq(q1)%vel(:, b1), &
                                                   dr%default_smearing(b1), opts%sigma)
        end do
    end do

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
            call compute_threephonon_scattering(il, sr, qp, dr, uc, fct, dims, opts%temperature, mw, mem)
        end if
        if (opts%fourthorder) then
            if (opts%nsample4ph .lt. 0) then
                call compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, dims, opts%temperature, rng, mw, mem)
            else
                call compute_fourphonon_sampling(il, sr, qp, dr, uc, fcf, dims, opts%temperature, opts%nsample4ph, rng, mw, mem)
            end if
        end if
        if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', il, sr%nlocal_point, walltime() - t0)
    end do
    if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', sr%nlocal_point, sr%nlocal_point, walltime() - t0)


    diagnostic: block
        integer(i8) :: bufiso, buf3ph, buf4ph
        real(r8) :: ratioiso, ratio3ph, ratio4ph
        real(r8), dimension(mw%n) :: mb, mbi, mb3, mb4

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

        ! We do it in steps otherwise the numbers are too big and it prints garbage
        ratioiso = real(bufiso, r8) / real(qp%n_irr_point, r8)
        ratioiso = ratioiso / real(qp%n_full_point, r8)
        ratioiso = ratioiso / real(dr%n_mode**2, r8)
        ratioiso = ratioiso * 100.0_r8

        ratio3ph = real(buf3ph, r8) / real(qp%n_irr_point, r8)
        ratio3ph = ratio3ph / real(qp%n_full_point)
        ratio3ph = ratio3ph / real(dr%n_mode**3)
        ratio3ph = ratio3ph * 100.0_r8

        ratio4ph = real(buf4ph, r8) / real(qp%n_irr_point, r8)
        ratio4ph = ratio4ph / real(qp%n_full_point, r8)
        ratio4ph = ratio4ph / real(qp%n_full_point, r8)
        ratio4ph = ratio4ph / real(dr%n_mode**4)
        ratio4ph = ratio4ph * 100.0_r8

        ! Compute the size of everything
        mb = 0.0_r8
        mbi = 0.0_r8
        mb3 = 0.0_r8
        mb4 = 0.0_r8
        mb(mw%r+1) = real(sr%size_in_mem(), r8) / 1024_r8**2
        mbi(mw%r+1) = real(sr%size_in_mem_iso(), r8) / 1024_r8**2
        mb3(mw%r+1) = real(sr%size_in_mem_3ph(), r8) / 1024_r8**2
        mb4(mw%r+1) = real(sr%size_in_mem_4ph(), r8) / 1024_r8**2
        call mw%allreduce('sum', mb)
        call mw%allreduce('sum', mbi)
        call mw%allreduce('sum', mb3)
        call mw%allreduce('sum', mb4)

        if (mw%talk) then
            write(*, *) ''
            write(*, *) '     number of isotope scattering matrix elements', bufiso
            write(*, *) '                                            % iso', ratioiso
            write(*, *) ' number of threephonon scattering matrix elements', buf3ph
            write(*, *) '                                            % 3ph', ratio3ph
            write(*, *) '  number of fourphonon scattering matrix elements', buf4ph
            write(*, *) '                                            % 4ph', ratio4ph
            write(*, *) ''
            write(*, "(1X,A30,':',4(1X,A20))") 'Memory usage (MiB)', 'total', 'avg per rank', 'max', 'min'
            write(*, "(1X,A30,':',4(1X,F20.6))") 'isotopes', sum(mbi), sum(mbi) / mw%n, maxval(mbi), minval(mbi)
            write(*, "(1X,A30,':',4(1X,F20.6))") '3ph', sum(mb3), sum(mb3) / mw%n, maxval(mb3), minval(mb3)
            write(*, "(1X,A30,':',4(1X,F20.6))") '4ph', sum(mb4), sum(mb4) / mw%n, maxval(mb4), minval(mb4)
            write(*, "(1X,A30,':',4(1X,F20.6))") 'total', sum(mb), sum(mb) / mw%n, maxval(mb), minval(mb)
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


subroutine compute_threephonon_scattering(il, sr, qp, dr, uc, fct, dims, temperature, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
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
    ! Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2
    ! The qpoints and the dimension of the qgrid
    real(r8), dimension(3) :: qv2, qv3
    ! The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sigma
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, plf, psisq, prefactor
    !
    complex(r8) :: c0
    ! Integers for do loops
    integer :: i, iq1, q1, q2, q3, b1, b2, b3
    !> Is the triplet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Already set some values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    sig1 = sr%sigma_q(q1, b1)

    i = 0
    do q2=1, qp%n_full_point
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)

        if (q3 .lt. q2) cycle

        call triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult)
        if (isred) cycle

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle
            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
                sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                if (abs(om1 + om2 - om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 + om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 - om3) .lt. 4.0_r8 * sigma) i = i + 1
            end do
        end do
    end do

    sr%threephonon(il)%n = i
    allocate(sr%threephonon(il)%psisq(i))
    allocate(sr%threephonon(il)%q2(i))
    allocate(sr%threephonon(il)%q3(i))
    allocate(sr%threephonon(il)%b2(i))
    allocate(sr%threephonon(il)%b3(i))

    i = 0
    do q2=1, qp%n_full_point
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)

        if (q3 .lt. q2) cycle

        call triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        call pretransform_phi3(fct, qv2, qv3, ptf)
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
            prefactor = threephonon_prefactor * qp%ap(q2)%integration_weight * mult

            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)

                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
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
    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine


subroutine compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, dims, temperature, rng, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    ! Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    ! The dimension of the grid
    integer, dimension(3), intent(in) :: dims
    ! The temperature
    real(r8), intent(in) :: temperature
    ! The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: fourphonon_prefactor = lo_pi / 8.0_r8
    ! Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2, evp3
    ! The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    ! The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> The ratio for the maximum likelihood estimation of the scattering rates
    real(r8) :: mle_ratio
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, psisq, prefactor
    ! The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sig4, sigma, plf
    ! Integers for do loops
    integer :: i, q1, q2, q3, q4, b1, b2, b3, b4, qi, qj
    ! Do we compute the scattering amplitude ?
    logical :: is_scatter
    !> Is the quartet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv4, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    sig1 = sr%sigma_q(q1, b1)

    i = 0
   do q2=1, qp%n_full_point
   do q3=q2, qp%n_full_point
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
        if (q4 .lt. q3) cycle

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
        if (isred) cycle
        is_scatter = .false.

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    sig4 = sr%sigma_q(qp%ap(q4)%irreducible_index, b4)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)
                    if (abs(om1 + om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 + om4 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 + om4 - om2) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 - om2 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om4 - om2 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. 4.0_r8 * sigma) i = i + 1
                end do
            end do
        end do
    end do
    end do

    sr%fourphonon(il)%n = i
    allocate(sr%fourphonon(il)%psisq(i))
    allocate(sr%fourphonon(il)%q2(i))
    allocate(sr%fourphonon(il)%q3(i))
    allocate(sr%fourphonon(il)%q4(i))
    allocate(sr%fourphonon(il)%b2(i))
    allocate(sr%fourphonon(il)%b3(i))
    allocate(sr%fourphonon(il)%b4(i))

    i = 0
    do q2=1, qp%n_full_point
    do q3=q2, qp%n_full_point
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
        if (q4 .lt. q3) cycle

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        prefactor = fourphonon_prefactor * qp%ap(q2)%integration_weight * qp%ap(q3)%integration_weight * mult
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    sig4 = sr%sigma_q(qp%ap(q4)%irreducible_index, b4)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)

                    if (abs(om1 + om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 + om4 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 + om4 - om2) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 - om2 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om4 - om2 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. 4.0_r8 * sigma) then
                        i = i + 1

                        egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)

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
    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine


subroutine compute_fourphonon_sampling(il, sr, qp, dr, uc, fcf, dims, temperature, n4ph, rng, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    ! Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    ! The dimension of the grid
    integer, dimension(3), intent(in) :: dims
    ! The temperature
    real(r8), intent(in) :: temperature
    ! Maximum number of process to take into account
    integer, intent(in) :: n4ph
    ! The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: fourphonon_prefactor = lo_pi / 8.0_r8
    ! Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    ! Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2, evp3
    ! The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    ! The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> The ratio for the maximum likelihood estimation of the scattering rates
    real(r8) :: mle_ratio
    ! Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, psisq, prefactor
    ! The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sig4, sigma, plf
    ! Integers for do loops
    integer :: i, j, ip, np, q1, q2, q3, q4, b1, b2, b3, b4, qi, qj
    ! Do we compute the scattering amplitude ?
    logical :: is_scatter
    !> Is the quartet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult
    !> np_prec
    integer :: ntot

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv4, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    sig1 = sr%sigma_q(q1, b1)

    qgridfull1 = [(q2, q2=1, qp%n_full_point)]
    qgridfull2 = [(q3, q3=1, qp%n_full_point)]
    call rng%shuffle_int_array(qgridfull1)
    call rng%shuffle_int_array(qgridfull2)

    ntot = 0
    count_loop: do qi=1, qp%n_full_point
    do qj=1, qp%n_full_point
        q2 = qgridfull1(qi)
        q3 = qgridfull2(qj)
        if (q3 .lt. q2) cycle  ! Permutation symmetry
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
        if (q4 .lt. q3) cycle  ! Permutation symmetry

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
        if (isred) cycle  ! Irreducible quartet

        ntot = ntot + dr%n_mode**4
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle
                        ntot = ntot + 1
                end do
            end do
        end do
    end do
    end do count_loop

    np = min(n4ph, ntot)
    mle_ratio = real(ntot, r8) / real(np, r8)

    sr%fourphonon(il)%n = np
    allocate(sr%fourphonon(il)%psisq(np))
    allocate(sr%fourphonon(il)%q2(np))
    allocate(sr%fourphonon(il)%q3(np))
    allocate(sr%fourphonon(il)%q4(np))
    allocate(sr%fourphonon(il)%b2(np))
    allocate(sr%fourphonon(il)%b3(np))
    allocate(sr%fourphonon(il)%b4(np))
    sr%fourphonon(il)%psisq = 0.0_r8

    ip = 0
    full_loop: do qi=1, qp%n_full_point
    do qj=1, qp%n_full_point
        q2 = qgridfull1(qi)
        q3 = qgridfull2(qj)
        if (q3 .lt. q2) cycle
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
        if (q4 .lt. q3) cycle

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        prefactor = fourphonon_prefactor * qp%ap(q2)%integration_weight * qp%ap(q3)%integration_weight * mult * mle_ratio
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    ip = ip +1
                    if (ip .gt. np) exit full_loop

                    sig4 = sr%sigma_q(qp%ap(q4)%irreducible_index, b4)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)

                    egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)

                    evp3 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                    evp3 = conjg(evp3)
                    c0 = dot_product(evp3, ptf)
                    psisq = abs(c0*conjg(c0)) * prefactor

                    sr%fourphonon(il)%psisq(ip) = psisq
                    sr%fourphonon(il)%q2(ip) = q2
                    sr%fourphonon(il)%q3(ip) = q3
                    sr%fourphonon(il)%q4(ip) = q4
                    sr%fourphonon(il)%b2(ip) = b2
                    sr%fourphonon(il)%b3(ip) = b3
                    sr%fourphonon(il)%b4(ip) = b4
                end do
            end do
        end do
    end do
    end do full_loop
    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
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


subroutine sr_destroy(sr)
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
    if (allocated(sr%be)) deallocate(sr%be)
    if (allocated(sr%sigma_q)) deallocate(sr%sigma_q)
    sr%nlocal_point = -lo_hugeint
end subroutine

! Function to measure the size of the memory
function sr_size_in_mem(sr) result(mem)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr
    !> size in memory, bytes
    integer(i8) :: mem

    integer :: il

    mem = 0
    mem = mem + storage_size(sr)
    if (allocated(sr%q1)) mem = mem + storage_size(sr%q1) * size(sr%q1)
    if (allocated(sr%b1)) mem = mem + storage_size(sr%b1) * size(sr%b1)
    if (allocated(sr%be)) mem = mem + storage_size(sr%be) * size(sr%be)
    if (allocated(sr%sigma_q)) mem = mem + storage_size(sr%sigma_q) * size(sr%sigma_q)
    mem = mem / 8
    if (allocated(sr%iso)) mem = mem + sr%size_in_mem_iso()
    if (allocated(sr%threephonon)) mem = mem + sr%size_in_mem_3ph()
    if (allocated(sr%fourphonon)) mem = mem + sr%size_in_mem_4ph()
end function

function size_in_mem_iso(sr) result(mem)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr
    !> size in memory, bytes
    integer(i8) :: mem

    integer :: il

    mem = 0
    if (allocated(sr%iso)) then
        do il=1, sr%nlocal_point
            mem = mem + storage_size(sr%iso(il)%psisq) * size(sr%iso(il)%psisq)
            mem = mem + storage_size(sr%iso(il)%q2) * size(sr%iso(il)%q2)
            mem = mem + storage_size(sr%iso(il)%b2) * size(sr%iso(il)%b2)
        end do
    end if
    mem = mem / 8
end function

function size_in_mem_3ph(sr) result(mem)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr
    !> size in memory, bytes
    integer(i8) :: mem

    integer :: il

    mem = 0
    if (allocated(sr%threephonon)) then
        do il=1, sr%nlocal_point
            mem = mem + storage_size(sr%threephonon(il)%psisq) * size(sr%threephonon(il)%psisq)
            mem = mem + storage_size(sr%threephonon(il)%q2) * size(sr%threephonon(il)%q2)
            mem = mem + storage_size(sr%threephonon(il)%q3) * size(sr%threephonon(il)%q3)
            mem = mem + storage_size(sr%threephonon(il)%b2) * size(sr%threephonon(il)%b2)
            mem = mem + storage_size(sr%threephonon(il)%b3) * size(sr%threephonon(il)%b3)
        end do
    end if
    mem = mem / 8
end function

function size_in_mem_4ph(sr) result(mem)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr
    !> size in memory, bytes
    integer(i8) :: mem

    integer :: il

    mem = 0
    if (allocated(sr%fourphonon)) then
        do il=1, sr%nlocal_point
            mem = mem + storage_size(sr%fourphonon(il)%psisq) * size(sr%fourphonon(il)%psisq)
            mem = mem + storage_size(sr%fourphonon(il)%q2) * size(sr%fourphonon(il)%q2)
            mem = mem + storage_size(sr%fourphonon(il)%q3) * size(sr%fourphonon(il)%q3)
            mem = mem + storage_size(sr%fourphonon(il)%q4) * size(sr%fourphonon(il)%q4)
            mem = mem + storage_size(sr%fourphonon(il)%b2) * size(sr%fourphonon(il)%b2)
            mem = mem + storage_size(sr%fourphonon(il)%b3) * size(sr%fourphonon(il)%b3)
            mem = mem + storage_size(sr%fourphonon(il)%b4) * size(sr%fourphonon(il)%b4)
        end do
    end if
    mem = mem / 8
end function


subroutine triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2 and q3 are full
    integer, intent(in) :: q1, q2, q3
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> If it's reducible, what is its multiplicity
    real(r8), intent(out) :: mult

    ! To hold the q-point in reduced coordinates
    real(r8), dimension(3) :: qv2, qv3, qv2p, qv3p
    !> To get the index of the new triplet on the fft_grid
    integer, dimension(3) :: gi
    !> The new triplet after the operation
    integer, dimension(2) :: qpp
    !> Integers for the loops
    integer :: j, k

    ! First get the q-points in reduce coordinates
    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
    qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)

    mult = 0.0_r8
    isred = .false.
    ! Let's try all operations that leaves q1 invariant
    do j=1, qp%ip(q1)%n_invariant_operation
        k = qp%ip(q1)%invariant_operation(j)
        select type(qp); type is(lo_fft_mesh)
            ! Rotate q2 and look if it's the on grid
            qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle

            ! Rotate q3 and look if it's the on grid
            qv3p = lo_operate_on_vector(uc%sym%op(k), qv3, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv3p) .eqv. .false.) cycle

            ! If everything is on the grid, get the index of each point
            gi = qp%index_from_coordinate(qv2p)
            qpp(1) = qp%gridind2ind(gi(1), gi(2), gi(3))
            gi = qp%index_from_coordinate(qv3p)
            qpp(2) = qp%gridind2ind(gi(1), gi(2), gi(3))
        end select
        if (minval(qpp) .gt. q2) then
            isred = .true.
        ! Now we have to determine the weight
        ! Two roads are possible
        !      1. Get the ratio of number of red point that can give this reducible point
        !      2. Look at the ratio between total number of operations and the ones
        !         that leaves this irreducible triplet unchanged
        ! The second road doesn't requires me to sum over all other qpoints, so I'll go with this one
        else if (minval(qpp) .eq. q2 .and. maxval(qpp) .eq. q3) then
            mult = mult + 1.0_r8
        end if
    end do
    ! The 2.0_r8 is for time reversal symmetry
    mult = qp%ip(q1)%n_invariant_operation / mult
end subroutine


subroutine quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2, q3 and q4 are full
    integer, intent(in) :: q1, q2, q3, q4
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> If it's reducible, what is its multiplicity
    real(r8), intent(out) :: mult

    ! To hold the q-point in reduced coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4, qv2p, qv3p, qv4p
    !> To get the index of the new triplet on the fft_grid
    integer, dimension(3) :: gi
    !> The new triplet after the operation
    integer, dimension(3) :: qpp
    !> Integers for the loops
    integer :: j, k

    ! First get the reciprocal lattice vectors, in reduce coordinates
    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
    qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)
    qv4 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q4)%r)

    isred = .false.
    mult = 0.0_r8
    ! Let's try all operations that leaves q1 invariant
    do j=1, qp%ip(q1)%n_invariant_operation
        k = qp%ip(q1)%invariant_operation(j)
        select type(qp); type is(lo_fft_mesh)
            ! Transform q2 and check if it's on the grid
            qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle

            ! Transform q3 and check if it's on the grid
            qv3p = lo_operate_on_vector(uc%sym%op(k), qv3, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv3p) .eqv. .false.) cycle

            ! Transform q4 and check if it's on the grid
            qv4p = lo_operate_on_vector(uc%sym%op(k), qv4, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv4p) .eqv. .false.) cycle

            ! If everything is on the grid, get the location of each point
            gi = qp%index_from_coordinate(qv2p)
            qpp(1) = qp%gridind2ind(gi(1), gi(2), gi(3))
            gi = qp%index_from_coordinate(qv3p)
            qpp(2) = qp%gridind2ind(gi(1), gi(2), gi(3))
            gi = qp%index_from_coordinate(qv4p)
            qpp(3) = qp%gridind2ind(gi(1), gi(2), gi(3))
        end select
        ! The sorting allows to include permutation invariance
        call lo_qsort(qpp)
        if (qpp(1) .gt. q2) then
            isred = .true.
        ! For the weights, it's the same idea that for the triplet
        ! I compute the number of operations that let the quartet invariant
        ! And take the ratio between the little group of q1 and this number
        else if (qpp(1) .eq. q2 .and. qpp(2) .eq. q3 .and. qpp(3) .eq. q4) then
            mult = mult + 1.0_r8
        end if
    end do
    mult = qp%ip(q1)%n_invariant_operation / mult
end subroutine

end module
