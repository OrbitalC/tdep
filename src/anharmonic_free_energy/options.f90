#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_author, lo_version, lo_licence, lo_tiny, lo_status, lo_exitcode_baddim, lo_exitcode_param, &
                      lo_frequency_Hartree_to_THz, lo_frequency_Hartree_to_icm, lo_frequency_Hartree_to_meV
use gottochblandat, only: lo_stop_gracefully
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    real(r8) :: temperature
    integer, dimension(3) :: qgrid, qg3ph, qg4ph
    integer :: integrationtype
    integer :: verbosity
    integer :: nblocks
    logical :: quantum = .false.
    logical :: stochastic = .false.
    logical :: thirdorder = .false.
    logical :: fourthorder = .false.
contains
    procedure :: parse
end type

contains

subroutine parse(opts)
    !> the options
    class(lo_opts), intent(out) :: opts
    !> the helper parser
    type(command_line_interface) :: cli

    logical :: dumlog
    real(r8), dimension(3) :: dumr8v
    integer :: errctr

    ! basic info
    call cli%init(progname='anharmonic_free_energy', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculates the anharmonic free energy.', &
                  examples=["mpirun anharmonic_free_energy"], &
                  epilog=new_line('a')//"...")

    ! Specify some options
    cli_qpoint_grid
    call cli%add(switch='--qpoint_grid3ph', switch_ab='-qg3ph', &
                 help='Dimension of the q-grid for three phonon contribution', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid4ph', switch_ab='-qg4ph', &
                 help='Dimension of the q-grid for four phonon contribution', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature', &
                 help='The temperature at which the thermodynamic properties are computed', &
                 required=.false., act='store', def='300', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--quantum', &
                 help='Use Bose-Einstein occupations to compute the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--stochastic', &
                 help='Add second order cumulant contribution to the free energy with a minus sign for self-consistent sampling', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--thirdorder', &
                 help='Compute third order anharmonic correction to the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder', &
                 help='Compute fourth order anharmonic correction to the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nblocks', &
                 help='Number of blocks used to compute uncertainty', &
                 required=.false., act='store', def='10', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_manpage
    cli_verbose

    ! actually parse it
    errctr = 0
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    ! Should the manpage be generated? In that case, no point to go further than here.
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog, error=lo_status); errctr = errctr + lo_status
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    call cli%get(switch='--verbose', val=dumlog, error=lo_status); errctr = errctr + lo_status; if (dumlog) opts%verbosity = 2

    ! get real options
    call cli%get(switch='--temperature', val=opts%temperature, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--qpoint_grid', val=opts%qgrid, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--quantum', val=opts%quantum, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--stochastic', val=opts%stochastic, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--thirdorder', val=opts%thirdorder, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--fourthorder', val=opts%fourthorder, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--nblocks', val=opts%nblocks, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--qpoint_grid3ph', val=opts%qg3ph, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--qpoint_grid4ph', val=opts%qg4ph, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--nblocks', val=opts%nblocks, error=lo_status); errctr = errctr + lo_status

    ! If we have fourthorder we should also have thirdorder
!   if (opts%fourthorder) opts%thirdorder = .true.

    if (errctr .ne. 0) call lo_stop_gracefully(['Failed parsing the command line options'], lo_exitcode_baddim)

    if (maxval(opts%qg3ph) .gt. 0 .and. .not. opts%thirdorder) then
        write(*, *) 'You have to enable thirdorder to use a three-phonon q-grid, stopping calculation.'
        stop
    end if
    if (minval(opts%qg3ph) .le. 0 .and. opts%thirdorder) then
        opts%qg3ph = opts%qgrid
    end if

    if (maxval(opts%qg4ph) .gt. 0 .and. .not. opts%fourthorder) then
        write(*, *) 'You have to enable fourthorder to use a four-phonon q-grid, stopping calculation.'
        stop
    end if
    if (minval(opts%qg4ph) .le. 0 .and. opts%fourthorder) then
        write(*, *) opts%qg4ph
        opts%qg4ph = opts%qgrid
    end if

end subroutine

end module
