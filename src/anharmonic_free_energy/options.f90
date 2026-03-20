#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_author, lo_version, lo_licence, lo_status, lo_exitcode_baddim
use gottochblandat, only: lo_stop_gracefully
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    real(r8) :: temperature
    integer, dimension(3) :: qgh, qg3ph, qg4ph
    integer :: integrationtype
    integer :: verbosity
    logical :: quantum = .false.
    logical :: nothirdorder = .false.
    logical :: nofourthorder = .false.
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
    call cli%add(switch='--qpoint_harmonic', switch_ab='-qgh', &
                 help='Dimension of the q-grid for harmonic contributions', &
                 nargs='3', required=.false., act='store', def='25 25 25', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid3ph', switch_ab='-qg3ph', &
                 help='Dimension of the q-grid for three phonon contribution', &
                 nargs='3', required=.false., act='store', def='20 20 20', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid4ph', switch_ab='-qg4ph', &
                 help='Dimension of the q-grid for four phonon contribution', &
                 nargs='3', required=.false., act='store', def='20 20 20', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--temperature', &
                 help='The temperature at which the thermodynamic properties are computed', &
                 required=.false., act='store', def='300', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--quantum', &
                 help='Use Bose-Einstein occupations to compute the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nothirdorder', &
                 help='Do NOT compute third order anharmonic correction to the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nofourthorder', &
                 help='Do NOT compute fourth order anharmonic correction to the free energy', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
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
    call cli%get(switch='--quantum', val=opts%quantum, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--nothirdorder', val=opts%nothirdorder, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--nofourthorder', val=opts%nofourthorder, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--qpoint_grid3ph', val=opts%qg3ph, error=lo_status); errctr = errctr + lo_status
    call cli%get(switch='--qpoint_grid4ph', val=opts%qg4ph, error=lo_status); errctr = errctr + lo_status

    if (errctr .ne. 0) call lo_stop_gracefully(['Failed parsing the command line options'], lo_exitcode_baddim)

    if (maxval(opts%qg3ph) .gt. 0 .and. opts%nothirdorder) then
        write(*, *) 'WARNING: You turned off thirdorder but still specified a three-phonon q-grid.'
    end if

    if (maxval(opts%qg4ph) .gt. 0 .and. opts%nofourthorder) then
        write(*, *) 'WARNING: You turned off fourthorder but still specified a four-phonon q-grid.'
    end if

end subroutine

end module
