#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, lo_m_to_bohr, lo_hugeint
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer, dimension(3) :: qgrid   !< the main q-grid
    integer, dimension(3) :: qg3ph   !< The grid for the threephonon integration
    integer, dimension(3) :: qg4ph   !< The grid for the fourphonon integration
    real(flyt) :: temperature        !< temperature
    real(flyt) :: sigma              !< scaling factor for adaptice gaussian
    integer :: nbasis                !< Number of basis function for the spectral functions
    logical :: readiso               !< read isotope distribution from file
    logical :: thirdorder            !< use fourth order contribution
    logical :: fourthorder           !< use fourth order contribution
    logical :: isotopescattering     !< use isotope scattering

    ! Debugging things
    logical :: timereversal
    logical :: qpsymmetry
    !
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
    !
    logical :: dumlog
    real(flyt) :: f0
    real(flyt), dimension(3) :: dumflytv
    integer :: i

    ! basic info
    call cli%init(progname='self_energy', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Calculates the lattice thermal conductivity from the iterative solution of the &
                              &phonon Boltzmann equation. In addition, cumulative plots and raw data dumps &
                              &of intermediate values are available.', &
                  examples=["mpirun thermal_conductivity --temperature 300                              ", &
                            "mpirun thermal_conductivity --fourthorder --nsample4ph 10000 -qg 30 30 30  "], &
                  epilog=new_line('a')//"...")
    ! real options
    call cli%add(switch='--readiso', &
                 help='Read the isotope distribution from `infile.isotopes`.', &
                 help_markdown='The format is specified [here](../page/files.html#infile.isotopes).', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nothirdorder', &
                 help='Not consider third order contributions to the scattering.',  &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--fourthorder', &
                 help='Consider four-phonon contributions to the scattering.',  &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_qpoint_grid
    call cli%add(switch='--sigma', &
                 help='Global scaling factor for adaptive Gaussian smearing.', &
                 required=.false., act='store', def='1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_readqmesh

    call cli%add(switch='--temperature', &
                 help='Evaluate thermal conductivity at a single temperature.', &
                 required=.false., act='store', def='300', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noisotope', &
                 help='Do not consider isotope scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nbasis', &
                 help='Number of basis function to represent the self-energy on the frequency axis.', &
                 required=.false., act='store', def='100', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid3ph', switch_ab='-qg3ph', &
                 help='Dimension of the grid for the threephonon integration.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--qpoint_grid4ph', switch_ab='-qg4ph', &
                 help='Dimension of the grid for the fourphonon integration.', &
                 nargs='3', required=.false., act='store', def='-1 -1 -1', error=lo_status)
    if (lo_status .ne. 0) stop

    ! hidden
    call cli%add(switch='--notr', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nosym', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop

    cli_manpage
    cli_verbose

    ! actually parse it
    call cli%parse(error=lo_status)
    if (lo_status .ne. 0) stop
    !
    ! Should the manpage be generated? In that case, no point to go further than here.
    !
    dumlog = .false.
    call cli%get(switch='--manpage', val=dumlog)
    if (dumlog) then
        call cli%save_man_page(trim(cli%progname)//'.1')
        call cli%save_usage_to_markdown(trim(cli%progname)//'.md')
        write (*, *) 'Wrote manpage for "'//trim(cli%progname)//'"'
        stop
    end if

    ! store things in the right place

    call cli%get(switch='--temperature', val=opts%temperature)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid)
    call cli%get(switch='--qpoint_grid3ph', val=opts%qg3ph)
    call cli%get(switch='--qpoint_grid4ph', val=opts%qg4ph)
    call cli%get(switch='--nbasis', val=opts%nbasis)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--nothirdorder', val=dumlog)
    opts%thirdorder = .not. dumlog
    call cli%get(switch='--fourthorder', val=opts%fourthorder)
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--notr', val=dumlog)
    opts%timereversal = .not. dumlog
    call cli%get(switch='--nosym', val=dumlog)
    opts%qpsymmetry = .not. dumlog
    call cli%get(switch='--verbose', val=dumlog)
    if (dumlog) then
        opts%verbosity = 2
    else
        opts%verbosity = 0
    end if
    call cli%get(switch='--noisotope', val=dumlog)
    opts%isotopescattering = .not. dumlog

    if (opts%thirdorder) then
        do i=1, 3
            if (opts%qg3ph(i) .lt. 0) opts%qg3ph(i) = opts%qgrid(i)
        end do
    end if
    if (opts%fourthorder) then
        do i=1, 3
            if (opts%qg4ph(i) .lt. 0) opts%qg4ph(i) = opts%qgrid(i)
        end do
    end if

end subroutine

end module
