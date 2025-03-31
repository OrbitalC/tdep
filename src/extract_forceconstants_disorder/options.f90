#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_status, lo_author, lo_version, lo_licence, lo_huge, lo_A_to_Bohr, lo_eV_to_Hartree, &
                      lo_forceconstant_2nd_eVA_to_HartreeBohr
use flap, only: command_line_interface
implicit none

private
public :: lo_opts

type lo_opts
    !> Second order cutoff
    real(r8) :: rc2 = -lo_huge
    !> Third order cutoff
    real(r8) :: rc3 = -lo_huge
    !> The gradient descent prefactor
    real(r8) :: alpha_secondorder = 1e-1_r8
    !> The gradient descent prefactor
    real(r8) :: alpha_pos = 1.0_r8
    !> The gradient descent prefactor for third order
    real(r8) :: alpha_thirdorder = 1e-1_r8
    !> The threshold to stop the gradient descent
    real(r8) :: thresh = 1e-4_r8
    !> The number of steps for second order
    integer :: nsteps = 1000
    !> The number of steps for third order
    integer :: nsteps_thirdorder = 1000
    !> The stride
    integer :: stride = 1
    !> The number of pval for third order
    integer :: p3 = 10
    !> Do we compute average positions ?
    logical :: avg_pos = .false.
    !> Do we minimize the positions ?
    logical :: update_pos = .false.
    !> Do we read previous IFC ?
    logical :: readifc = .false.
    !> Do we prefit IFCs ?
    logical :: prefit = .true.
    !> Do we fit third order IFCs ?
    logical :: thirdorder
    !> how much to talk
    integer :: verbosity
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

    call cli%init(progname='extract_forceconstants_disorder', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Extract force constant for a disordered material with all symmetries broken.', &
                  examples=["extract_forceconstants_disorder -rc2 6"], &
                  epilog=new_line('a')//"...")
    cli_manpage
    cli_verbose
    call cli%add(switch='--secondorder_cutoff', switch_ab='-rc2', &
                 help='Cutoff for the second order force constants.', &
                 required=.false., act='store', def='5.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--thirdorder_cutoff', switch_ab='-rc3', &
                 help='Cutoff for the third order force constants.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--alpha_secondorder', switch_ab='-a2', &
                 help='Prefactor for the gradient descent for second order force constants.', &
                 required=.false., act='store', def='0.1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--alpha_positions', &
                 help='Prefactor for the update of positions in the gradient descent.', &
                 required=.false., act='store', def='1E-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--alpha_thirdorder', switch_ab='-a3', &
                 help='Prefactor for the gradient descent for third order force constants.', &
                 required=.false., act='store', def='1E-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--threshold', &
                 help='Threshold to stop the gradient descent.', &
                 required=.false., act='store', def='5e-4', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nsteps', switch_ab='-ns', &
                 help='Number of steps in the gradient descent for second order force constants.', &
                 required=.false., act='store', def='1000', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nsteps_thirdorder', switch_ab='-ns3', &
                 help='Number of steps in the gradient descent for third order force constants.', &
                 required=.false., act='store', def='1000', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--stride', switch_ab='-s', &
                 help='Use every N configurations instead of all.', &
                 required=.false., act='store', def='1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--pval_thirdorder', switch_ab='-p3', &
                 help='Number of principal values for third order force constants.', &
                 required=.false., act='store', def='10', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--average_positions', &
                 help='Compute the average positions before fitting the force constants.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--update_positions', &
                 help='Update the positions during the gradient descent.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--readifc', &
                 help='Read prefiously fitted IFC as a starting point.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noprefit', &
                 help='Do not nitialize IFC by doing per-atom least-squares fits.', &
                 required=.false., act='store_false', def='.true.', error=lo_status)
    if (lo_status .ne. 0) stop

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop

    ! generate manpage?
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if
    opts%verbosity = 0
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) opts%verbosity = 2

    ! Parse the rest
    call cli%get(switch='--secondorder_cutoff', val=opts%rc2, error=lo_status)
    call cli%get(switch='--thirdorder_cutoff', val=opts%rc3, error=lo_status)
    call cli%get(switch='--alpha_secondorder', val=opts%alpha_secondorder, error=lo_status)
    call cli%get(switch='--alpha_positions', val=opts%alpha_pos, error=lo_status)
    call cli%get(switch='--alpha_thirdorder', val=opts%alpha_thirdorder, error=lo_status)
    call cli%get(switch='--threshold', val=opts%thresh, error=lo_status)
    call cli%get(switch='--nsteps', val=opts%nsteps, error=lo_status)
    call cli%get(switch='--nsteps_thirdorder', val=opts%nsteps_thirdorder, error=lo_status)
    call cli%get(switch='--stride', val=opts%stride, error=lo_status)
    call cli%get(switch='--pval_thirdorder', val=opts%p3, error=lo_status)
    call cli%get(switch='--average_positions', val=opts%avg_pos, error=lo_status)
    call cli%get(switch='--update_positions', val=opts%update_pos, error=lo_status)
    call cli%get(switch='--readifc', val=opts%readifc, error=lo_status)
    call cli%get(switch='--noprefit', val=opts%prefit, error=lo_status)

    ! Some unit conversion
    opts%rc2 = opts%rc2 * lo_A_to_Bohr
    opts%rc3 = opts%rc3 * lo_A_to_Bohr
    opts%thresh = opts%thresh * lo_eV_to_Hartree

    ! Now for some sanity checks
    if (opts%alpha_secondorder .lt. 0.0_r8) then
        write(*, *) 'alpha should be positive'
        stop
    end if

    ! Do we fit the third order ?
    if (opts%rc3 .gt. 0.0_r8) opts%thirdorder = .true.

end subroutine

end module
