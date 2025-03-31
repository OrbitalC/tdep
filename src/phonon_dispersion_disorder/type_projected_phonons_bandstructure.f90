module type_projected_phonons_bandstructure
!!
!! Deal with phonons projected on a unitcell along a path in the BZ
!!
use konstanter, only: r8, lo_hugeint, lo_exitcode_param, lo_exitcode_symmetry, lo_Bohr_to_A, &
                      lo_frequency_Hartree_to_THz, lo_frequency_Hartree_to_icm, lo_frequency_Hartree_to_meV
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar, lo_progressbar_init
use hdf5_wrappers, only: lo_hdf5_helper
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
use type_unitcell_disorder_projection, only: lo_unitcell_disorder_projection
use type_qpointmesh, only: lo_bandstructure
use type_projected_phonons, only: lo_projected_phonons, lo_projected_phonons_qpoints


implicit none
private
public :: lo_projected_phonons_bandstructure

type, extends(lo_bandstructure) :: lo_projected_phonons_bandstructure
    !> The number of phonon modes
    integer :: n_mode = -lo_hugeint
    !> Points that hold dispersion data
    type(lo_projected_phonons_qpoints), dimension(:), allocatable :: p

    !> Energy axis for spectral function
    real(r8), dimension(:), allocatable :: energy_axis
contains
    !> Calculate the band structure
    procedure :: generate
    !> Write to HDF5
    procedure :: write_to_hdf5
end type

contains

!> Calculate the band structure
subroutine generate(bs, fc, uc, ss, proj, ph_com, mw, mem, npts, readpathfromfile)
    !> The band structure
    class(lo_projected_phonons_bandstructure), intent(out) :: bs
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> Crystal structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> Crystal structure
    type(lo_crystalstructure), intent(inout) :: ss
    !> The phonon projected on commensurate q-points
    type(lo_projected_phonons), intent(in) :: ph_com
    !> The unitcell projection
    type(lo_unitcell_disorder_projection), intent(in) :: proj
    !> the MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The number of points per path
    integer, optional, intent(in) :: npts
    !> Do we read the path
    logical, optional, intent(in) :: readpathfromfile

    init: block
        integer :: i
        logical :: readpath

        ! Make sure all the symmetry stuff is there!
        if (uc%info%havewedge .eqv. .false.) then
            call uc%classify('wedge', timereversal=.true.)
        end if
        if (uc%info%decidedtimereversal .eqv. .false.) then
            call lo_stop_gracefully(['time-reversal symmetry needs to be decided for dispersions on a path'], &
                                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        if (uc%sym%timereversal .neqv. .true.) then
            call lo_stop_gracefully(['Conflicting instructions regarding time-reversal symmetry for dispersions on path'], &
                                    lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if

        ! Number of points
        if (present(npts)) then
            bs%n_point_per_path = npts
        else
            bs%n_point_per_path = 100
        end if

        ! read or generate path
        if (present(readpathfromfile)) then
            readpath = readpathfromfile
        else
            readpath = .false.
        end if

        ! Coordinates of the path
        if (readpath) then
            call bs%read_path_from_file(uc, mw, 0)
        else
            call bs%standardpath(uc, mw, 0)
        end if

        bs%n_mode = uc%na * 3
        allocate(bs%p(bs%n_point))
        do i=1, bs%n_point
            allocate(bs%p(i)%omega(bs%n_mode))
            allocate(bs%p(i)%egv(bs%n_mode, bs%n_mode))
            allocate(bs%p(i)%vel(3, bs%n_mode, bs%n_mode))
            bs%p(i)%omega = 0.0_r8
            bs%p(i)%egv = 0.0_r8
            bs%p(i)%vel = 0.0_r8
            bs%p(i)%qpoint%r = bs%q(i)%r
        end do
    end block init

    getdispersions: block
        complex(r8), dimension(:, :, :), allocatable :: Dq
        complex(r8), dimension(:, :), allocatable :: D
        real(r8) :: t0
        integer :: q, lq, nq

        ! Allocate the dynamical matrix and its gradient
        call mem%allocate(D, [bs%n_mode, bs%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(Dq, [3, bs%n_mode, bs%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        D = 0.0_r8
        Dq = 0.0_r8

        ! Count the number of point on this rank
        nq = 0
        do q=1, bs%n_point
            if (mod(q, mw%n) .eq. mw%r) nq = nq + 1
        end do
        if (mw%talk) call lo_progressbar_init()
        t0 = walltime()

        lq = 0
        do q=1, bs%n_point
            if (mod(q, mw%n) .ne. mw%r) cycle
            lq = lq + 1

            ! First we compute the harmonic q-point
            call bs%p(q)%compute_harmonic_qpoint(fc, ss, proj, mem)

            if (mw%talk .and. lo_trueNtimes(q, 127, nq)) then
                call lo_progressbar(' ... calculating frequencies', lq, nq, walltime() - t0)
            end if
        end do
        if (mw%talk) call lo_progressbar(' ... calculating frequencies', nq, nq, walltime() - t0)
        call mem%deallocate(D, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(Dq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Now we distribute
        do q=1, bs%n_point
            call mw%allreduce('sum', bs%p(q)%omega)
            call mw%allreduce('sum', bs%p(q)%egv)
            call mw%allreduce('sum', bs%p(q)%vel)
        end do
    end block getdispersions
end subroutine

subroutine write_to_hdf5(bs, filename, enhet, mem)
    !> The band structure
    class(lo_projected_phonons_bandstructure), intent(in) :: bs
    !> Filename
    character(len=*), intent(in) :: filename
    !> unit
    character(len=3), intent(in) :: enhet
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The HDF5 helper
    type(lo_hdf5_helper) :: h5
    !> The unitfactor
    real(r8) :: unitfactor
    !> The unit character
    character(len=1000) :: unitstr, lblstr
    !> To write the path
    character(len=10), dimension(:), allocatable :: dumlbl
    !> Some buffers for the writing
    real(r8), dimension(:, :), allocatable :: d2
    real(r8), dimension(:, :, :), allocatable :: d3
    !> And some integers
    integer :: i, j, k


    select case (enhet)
    case ('thz')
        unitfactor = lo_frequency_hartree_to_THz
        unitstr = 'THz'
    case ('mev')
        unitfactor = lo_frequency_hartree_to_meV
        unitstr = 'meV'
    case ('icm')
        unitfactor = lo_frequency_hartree_to_icm
        unitstr = 'cm^-^1'
    case default
        call lo_stop_gracefully(['Unknown unit'], lo_exitcode_param, __FILE__, __LINE__)
    end select


    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Let's start with the x-axis
    ! We dump the scalar x-axis
    call h5%store_data(bs%q_axis/lo_bohr_to_A, h5%file_id, 'q_values', enhet='A^-1', dimensions='q-point')
    ! The ticks for the x-axis
    call h5%store_data(bs%q_axis_ticks/lo_bohr_to_A, h5%file_id, 'q_ticks', enhet='A^-1')
    ! The labels for the x-ticks
    allocate(dumlbl(size(bs%q_axis_tick_labels, 1)))
    do i = 1, size(dumlbl, 1)
        k = 0
        dumlbl(i) = '         '
        do j = 1, len(dumlbl(i))
            if (bs%q_axis_tick_labels(i) (j:j) .ne. '|') then
                k = k + 1
                dumlbl(i) (k:k) = trim(bs%q_axis_tick_labels(i) (j:j))
            end if
        end do
        if (trim(dumlbl(i)) .eq. 'Γ') then
            dumlbl(i) = 'G'
        end if
    end do
    lblstr = ''
    do i = 1, size(dumlbl, 1)
        lblstr = trim(lblstr)//' '//trim(dumlbl(i))
    end do
    call h5%store_attribute(trim(adjustl(lblstr)), h5%file_id, 'q_tick_labels')
    deallocate(dumlbl)

    ! Now we dump the frequencies
    if (allocated(bs%p(1)%omega)) then
        call mem%allocate(d2, [bs%n_mode, bs%n_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        d2 = 0.0_r8
        do i = 1, bs%n_point
            d2(:, i) = bs%p(i)%omega*unitfactor
        end do
        call h5%store_data(d2, h5%file_id, 'frequencies', enhet=trim(unitstr), dimensions='q-point,mode')
        call mem%deallocate(d2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)

end subroutine
end module
