#include "precompilerdefinitions"
module options
use konstanter, only: flyt, lo_status, lo_author, lo_version, lo_licence, lo_m_to_bohr, lo_hugeint
use flap, only: command_line_interface
implicit none
private
public :: lo_opts

type lo_opts
    integer, dimension(3) :: qgrid   !< the main q-grid
    logical :: readqmesh             !< read q-grid from file
    real(flyt) :: temperature        !< temperature
    real(flyt) :: sigma              !< scaling factor for adaptice gaussian
    real(flyt) :: thres              !< consider Gaussian 0 if x-mu is larger than this number times sigma.
    real(flyt) :: tau_boundary       !< add a constant as boundary scattering
    real(flyt) :: mfp_max            !< add a length as boundary scattering
    real(flyt) :: mixing             !< mixing parameter for self consistent linewidth
    real(flyt) :: scftol             !< tolerance for the self-consistent linewidth
    real(flyt) :: btetol             !< tolerance for the self-consistent linewidth
    integer :: nsample3ph            !< the number of 3ph scattering process to actually compute
    integer :: nsample4ph            !< the number of 4ph scattering process to actually compute
    integer :: niter                 !< Number of iteration for the self consistent linewidths
    integer :: bteniter              !< Number of iteration for the self consistent linewidths
    logical :: readiso               !< read isotope distribution from file
    logical :: thirdorder            !< use fourth order contribution
    logical :: fourthorder           !< use fourth order contribution
    logical :: readlw                !< Do we read the linewidths from a file ?

    integer :: correctionlevel       !< how hard to correct
    integer :: mfppts                !< number of points on mfp-plots
    logical :: dumpgrid              !< print everything on a grid
    !logical :: thinfilm             !< Austins thin film thing

    ! Debugging things
    logical :: timereversal
    logical :: qpsymmetry
    logical :: isotopescattering
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

    ! basic info
    call cli%init(progname='thermal_conductivity_sampling', &
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
    call cli%add(switch='--read_linewidths', &
                 help='Read the linewidths from a infile.grid_thermal_conductivity_sampling.hdf5.',  &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_qpoint_grid
    call cli%add(switch='--sigma', &
                 help='Global scaling factor for adaptive Gaussian smearing.', &
                 required=.false., act='store', def='1.0', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--threshold', &
                 help='Consider a Gaussian distribution to be 0 after this many standard deviations.', &
                 required=.false., act='store', def='4.0', error=lo_status)
    if (lo_status .ne. 0) stop
    cli_readqmesh

    call cli%add(switch='--temperature', &
                 help='Evaluate thermal conductivity at a single temperature.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--max_mfp', &
                 help='Add a limit on the mean free path as an approximation of domain size.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--mixing', &
                 help='Mixing parameter for the self-consistent linewidth.', &
                 required=.false., act='store', def='0.75', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--scftol', &
                 help='Tolerance for the self-consistent linewidth.', &
                 required=.false., act='store', def='1e-2', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--btetol', &
                 help='Tolerance for the iterative BTE solution.', &
                 required=.false., act='store', def='1e-5', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--dumpgrid', &
                 help='Write files with q-vectors, frequencies, eigenvectors and group velocities for a grid.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noisotope', &
                 help='Do not consider isotope scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nsample3ph', &
                 help='The number of 3 phonon scattering to sample to estimate the lifetimes for each mode.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nsample4ph', &
                 help='The number of 4 phonon scattering to sample to estimate the lifetimes for each mode.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--niter', &
                 help='Number of iterations for the self consistent computation of the linewidths.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--bte_niter', &
                 help='Number of iterations for the iterative Boltzmann equation.', &
                 required=.false., act='store', def='200', error=lo_status)
    if (lo_status .ne. 0) stop

    ! hidden
    call cli%add(switch='--tau_boundary', hidden=.true., &
                 help='Add a constant boundary scattering term to the lifetimes.', &
                 required=.false., act='store', def='-1', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--mfppts', hidden=.true., help='', &
                 required=.false., act='store', def='200', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--notr', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--nosym', hidden=.true., help='', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--correctionlevel', hidden=.true., &
                 help='How agressively things are corrected due to broken symmetries.', &
                 required=.false., act='store', def='4', error=lo_status)
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
    call cli%get(switch='--nsample3ph', val=opts%nsample3ph)
    if (opts%nsample3ph .lt. 0) opts%nsample3ph = lo_hugeint
    call cli%get(switch='--nsample4ph', val=opts%nsample4ph)
    if (opts%nsample4ph .lt. 0) opts%nsample4ph = lo_hugeint
    call cli%get(switch='--niter', val=opts%niter)
    call cli%get(switch='--bte_niter', val=opts%bteniter)
    call cli%get(switch='--sigma', val=opts%sigma)
    call cli%get(switch='--threshold', val=opts%thres)
    call cli%get(switch='--tau_boundary', val=opts%tau_boundary)
    call cli%get(switch='--nothirdorder', val=dumlog)
    opts%thirdorder = .not. dumlog
    call cli%get(switch='--fourthorder', val=opts%fourthorder)
    if (opts%tau_boundary .gt. 0.0_flyt) opts%tau_boundary = 1E10_flyt
    call cli%get(switch='--read_linewidths', val=opts%readlw)
    call cli%get(switch='--readqmesh', val=opts%readqmesh)
    call cli%get(switch='--readiso', val=opts%readiso)
    call cli%get(switch='--mfppts', val=opts%mfppts)
    call cli%get(switch='--max_mfp', val=opts%mfp_max)
    call cli%get(switch='--mixing', val=opts%mixing)
    if (opts%mixing .lt. 0.0_flyt .or. opts%mixing .gt. 1.0_flyt) then
        write (*, *) 'Mixing parameter should be between 0.0 and 1.0.'
        stop
    end if
    call cli%get(switch='--scftol', val=opts%scftol)
    call cli%get(switch='--btetol', val=opts%btetol)
    call cli%get(switch='--dumpgrid', val=opts%dumpgrid)
    ! stuff that's not really an option
    call cli%get(switch='--correctionlevel', val=opts%correctionlevel)
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

    ! Get things to atomic units
    opts%mfp_max = opts%mfp_max*lo_m_to_Bohr

end subroutine

end module
