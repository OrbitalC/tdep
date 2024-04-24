#include "precompilerdefinitions"
module scattering
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

    contains
        !> Generate the scattering amplitudes
        procedure :: generate
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


type lo_montecarlo_grid
    !> The size of the grid
    integer :: npoints
    !> The weight of each point on the monte-carlo grid
    real(r8) :: weight
    !> The dimensions of the grid
    integer, dimension(3) :: mc_dims
    !> The dimensions of the full grid
    integer, dimension(3) :: full_dims
    !> The ratio between the full and mc grids
    real(r8), dimension(3) :: ratio

contains
    procedure :: initialize => initialize_montecarlo_grid
    procedure :: mc_point_to_full
    procedure :: generate_grid
end type

contains
subroutine generate(sr, qp, dr, uc, fct, fcf, opts, mw, mem)
    !> The scattering rate
    class(lo_scattering_rates), intent(out) :: sr
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The third order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> The fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> The options
    type(lo_opts) :: opts
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> The q-point grid dimension
    integer, dimension(3) :: dims
    !> The random number generator
    type(lo_mersennetwister) :: rng
    !> To initialize the random number generator and timing
    real(r8) :: rseed, t0
    !> Some integers
    integer :: q1, b1, q1f, il, j, k, nlocal_point, ctr
    !> The grids for monte-carlo integration
    type(lo_montecarlo_grid) :: mcg3, mcg4

    ! grid dimensions
    select type(qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    ! Initialize the random number generator
    call rng%init(iseed=mw%r, rseed=walltime())

    ! Initialize the monte-carlo grid
    if (opts%thirdorder) then
        call mcg3%initialize(dims, opts%qg3ph)
    end if
   if (opts%fourthorder) then
       call mcg4%initialize(dims, opts%qg4ph)
   end if

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
            call compute_threephonon_scattering(il, sr, qp, dr, uc, fct, mcg3, rng, mw, mem)
        end if
        if (opts%fourthorder) then
            call compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg4, rng, mw, mem)
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


#include "scattering_isotope.f90"
#include "scattering_threephonon.f90"
#include "scattering_fourphonon.f90"

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

subroutine initialize_montecarlo_grid(mcg, full_dims, mc_dims)
    !> The monte carlo grid
    class(lo_montecarlo_grid), intent(out) :: mcg
    !> The dimensions of the full grid
    integer, dimension(3), intent(in) :: full_dims
    !> The dimensions of the monte carlo grid
    integer, dimension(3), intent(in) :: mc_dims

    !> Some integers for the do loop
    integer :: i

    mcg%full_dims = full_dims
    do i=1, 3
        mcg%mc_dims(i) = min(mc_dims(i), mcg%full_dims(i))
        mcg%ratio(i) = real(mcg%full_dims(i), r8) / real(mcg%mc_dims(i), r8)
    end do
    mcg%npoints = mcg%mc_dims(1) * mcg%mc_dims(2) * mcg%mc_dims(3)
    mcg%weight = 1.0_r8 / real(mcg%npoints, r8)
end subroutine

function mc_point_to_full(mcg, imc, rng) result(ifull)
    !> The monte carlo grid
    class(lo_montecarlo_grid), intent(in) :: mcg
    !> The point on the monte-carlo grid
    integer, intent(in) :: imc
    ! The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The point on the full grid
    integer :: ifull
    !> The triplets of point on the monte-carlo and full grids
    integer, dimension(3) :: gi_mc, gi_full

    gi_mc = singlet_to_triplet(imc, mcg%mc_dims(2), mcg%mc_dims(3))
    ! This way of generating the number makes it work even if the ratio is not an integer
    gi_full(1) = ceiling((real(gi_mc(1), r8) - 1.0_r8) * mcg%ratio(1) + rng%rnd_real() * mcg%ratio(1))
    gi_full(2) = ceiling((real(gi_mc(2), r8) - 1.0_r8) * mcg%ratio(2) + rng%rnd_real() * mcg%ratio(2))
    gi_full(3) = ceiling((real(gi_mc(3), r8) - 1.0_r8) * mcg%ratio(3) + rng%rnd_real() * mcg%ratio(3))
    ifull = triplet_to_singlet(gi_full, mcg%full_dims(2), mcg%full_dims(3))
end function

subroutine generate_grid(mcg, qgrid, rng)
    !> The monte carlo grid
    class(lo_montecarlo_grid), intent(in) :: mcg
    ! The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The grid to be generated
    integer, dimension(:), intent(out) :: qgrid

    integer :: qi, qprev, qtest

    ! To improve convergence, we avoid repeating points in the integration grid
    qgrid(1) = mcg%mc_point_to_full(1, rng)
    qprev = qgrid(1)
    qtest = qgrid(1)
    do qi=2, mcg%npoints
        do while(qtest .eq. qprev)
            qtest = mcg%mc_point_to_full(qi, rng)
        end do
        qgrid(qi) = qtest
        qprev = qtest
    end do
end subroutine
end module
