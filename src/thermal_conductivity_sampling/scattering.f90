#include "precompilerdefinitions"
module scattering
use konstanter, only: r8, i8, lo_freqtol, lo_twopi, lo_exitcode_param, lo_hugeint, lo_pi, lo_tol, &
                      lo_phonongroupveltol, lo_tol, lo_frequency_THz_to_Hartree
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
use lo_fftgrid_helper, only: lo_montecarlo_grid, singlet_to_triplet, triplet_to_singlet, &
                             fft_third_grid_index, fft_fourth_grid_index

use options, only: lo_opts

use type_symmetryoperation, only: lo_operate_on_vector
use lo_sorting, only: lo_qsort

implicit none
private
public :: lo_scattering_rates

real(r8), parameter :: isotope_prefactor = lo_pi / 4.0_r8
real(r8), parameter :: threephonon_prefactor = lo_pi / 16.0_r8
real(r8), parameter :: fourphonon_prefactor = lo_pi / 96.0_r8

! Container for scattering rates
type lo_scattering_rates
    !> The number of qpoint/mode on this rank
    integer :: nlocal_point
    !> The list of qpoint and modes for this rank
    integer, dimension(:), allocatable :: q1, b1
    !> Let's precompute the Bose-Einstein distribution
    real(r8), dimension(:, :), allocatable :: be, sigsq
    !> The scattering matrix
    real(r8), dimension(:, :), allocatable :: Xi

    contains
        !> Generate the scattering amplitudes
        procedure :: generate
        !> destroy the scattering amplitues
        procedure :: destroy => sr_destroy
        !> Measure size in memory, in bytes
        procedure :: size_in_mem => sr_size_in_mem
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
    real(r8) :: rseed, t0, sigma
    !> Some integers
    integer :: q1, b1, il, j, nlocal_point, ctr
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
    allocate(sr%sigsq(qp%n_irr_point, dr%n_mode))
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            sr%be(q1, b1) = lo_planck(opts%temperature, dr%iq(q1)%omega(b1))

            sr%sigsq(q1, b1) = qp%smearingparameter(dr%iq(q1)%vel(:, b1), dr%default_smearing(b1), 1.0_r8)**2
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
    allocate(sr%Xi(nlocal_point, qp%n_full_point * dr%n_mode))

    ! Let's set it to zero
    sr%Xi = 0.0_r8

    ! Let's attribute the q1/b1 indices to the ranks
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

    memory_estimate: block
        real(r8), dimension(mw%n) :: mb
        integer :: nrow, ncol

        nrow = qp%n_irr_point * dr%n_mode
        ncol = qp%n_full_point * dr%n_mode

        mb(mw%r+1) = real(sr%size_in_mem(), r8) / 1024_r8**2
        call mw%allreduce('sum', mb)
        if (mw%talk) then
            write(*, *) ''
            write(*, *) 'Lower bound estimate of the memory usage :'
            write(*, "(1X,A31,':',1X,I10,A2,I10)") 'size of the scattering matrix', nrow, 'x', ncol
            write(*, "(1X,A31,':',4(1X,A20))") 'Memory usage (MiB)', 'total', 'avg per rank', 'max', 'min'
            write(*, "(1X,A31,1X,4(1X,F20.6))") '', sum(mb), sum(mb) / mw%n, maxval(mb), minval(mb)
            write(*, *) ''
        end if
    end block memory_estimate

    scatt: block
        !> Buffer to contains the linewidth
        real(r8), dimension(:, :), allocatable :: buf_lw
        !> Buffer for the linewidth of the local point
        real(r8) :: buf, f0, velnorm
        !> Some integers for the loops
        integer :: j, q1, b1, b2

        call mem%allocate(buf_lw, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        buf_lw = 0.0_r8

        t0 = walltime()
        if (mw%talk) call lo_progressbar_init()
        do il=1, sr%nlocal_point
            buf = 0.0_r8
            if (opts%isotopescattering) then
                call compute_isotope_scattering(il, sr, qp, dr, uc, opts%temperature, opts%thres, buf, &
                    opts%integrationtype, opts%sigma, mw, mem)
            end if
            if (opts%thirdorder) then
                call compute_threephonon_scattering(il, sr, qp, dr, uc, fct, mcg3, rng, &
                    opts%thres, buf, opts%integrationtype, opts%sigma, mw, mem)
            end if
            if (opts%fourthorder) then
                call compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg4, rng, &
                    opts%thres, buf, opts%integrationtype, opts%sigma, mw, mem)
            end if
            ! We end with the boundary scattering
            if (opts%mfp_max .gt. 0.0_r8) then
                velnorm = norm2(dr%iq(sr%q1(il))%vel(:, sr%b1(il)))
                if (velnorm .gt. lo_phonongroupveltol) then
                    buf = buf + velnorm / opts%mfp_max
                end if
            end if
            ! Now we can update the linewidth for this mode
            buf_lw(sr%q1(il), sr%b1(il)) = buf

            if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', il, sr%nlocal_point, walltime() - t0)
        end do
        if (mw%talk) call lo_progressbar(' ... computing scattering amplitude', sr%nlocal_point, sr%nlocal_point, walltime() - t0)

        ! Reduce the linewidth
        call mw%allreduce('sum', buf_lw)

        ! Distribute it after fixing the degeneracies
        do q1=1, qp%n_irr_point
            do b1=1, dr%n_mode
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
                ! First we fix the degeneracy
                f0 = 0.0_r8
                do j=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    f0 = f0 + buf_lw(q1, b2)
                end do
                f0 = f0 / real(dr%iq(q1)%degeneracy(b1), r8)
                do j=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    buf_lw(q1, b2) = f0
                end do
                ! Now we can set the linewidth
                dr%iq(q1)%linewidth(b1) = buf_lw(q1, b1)
            end do
        end do
        call mem%deallocate(buf_lw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block scatt
end subroutine

#include "scattering_isotope.f90"
#include "scattering_threephonon.f90"
#include "scattering_fourphonon.f90"

subroutine sr_destroy(sr)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr

    integer :: il

    if (allocated(sr%Xi)) deallocate(sr%Xi)
    if (allocated(sr%q1)) deallocate(sr%q1)
    if (allocated(sr%b1)) deallocate(sr%b1)
    if (allocated(sr%be)) deallocate(sr%be)
    if (allocated(sr%sigsq)) deallocate(sr%sigsq)
    sr%nlocal_point = -lo_hugeint
end subroutine

! Function to measure the size of the memory
function sr_size_in_mem(sr) result(mem)
    !> The scattering amplitudes
    class(lo_scattering_rates), intent(inout) :: sr
    !> size in memory, bytes
    integer(i8) :: mem

    integer :: il

    mem = storage_size(sr)
    if (allocated(sr%q1)) mem = mem + storage_size(sr%q1) * size(sr%q1)
    if (allocated(sr%b1)) mem = mem + storage_size(sr%b1) * size(sr%b1)
    if (allocated(sr%be)) mem = mem + storage_size(sr%be) * size(sr%be)
    if (allocated(sr%Xi)) mem = mem + storage_size(sr%Xi) * size(sr%Xi)
    mem = mem / 8
end function
end module
