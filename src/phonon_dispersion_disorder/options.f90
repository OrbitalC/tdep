#include "precompilerdefinitions"
module options
use konstanter, only: r8, lo_status, lo_author, lo_version, lo_licence, lo_huge, lo_hugeint
use flap, only: command_line_interface

implicit none

private
public :: lo_opts

type lo_opts
    !> What units are we printing
    character(len=3) :: enhet = '***'
    !> What is the temperature
    real(r8) :: temperature = -lo_huge
    !> Do we project on a unitcell reciprocal space ?
    logical :: projected
    !> Do we compute the structure factor ?
    logical :: structurefactor
    !> Do we compute the third order lineshape ?
    logical :: thirdorder
    !> The size of the qgrid
    integer, dimension(3) :: qgrid
    !> The number of frequencies on the axe
    integer :: nf = -lo_hugeint
    !> The number of qpoint on the band structure
    integer :: nq = -lo_hugeint
    !> Do we use isotope scattering
    logical :: isotopescattering
    !> Do we read the isotope scattering
    logical :: readiso
    !> Will we need the spectral function in real space
    logical :: sf_realspace
    !> Do we read the path from file ?
    logical :: readpathfromfile
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

    call cli%init(progname='phonon_relation_disorder', &
                  authors=lo_author, &
                  version=lo_version, &
                  license=lo_licence, &
                  help='Usage: ', &
                  description='Compute phonons for disordered systems.', &
                  examples=["phonon_relation_disorder --structure_factor"], &
                  epilog=new_line('a')//"...")

    cli_unit
    cli_temperature
    cli_manpage
    cli_verbose
    cli_qpoint_grid
    cli_nq_on_path
    cli_readpath
    call cli%add(switch='--structure_factor', switch_ab='-sf', &
                 help='Compute the structure factor, for an isotropic system without a basis.', &
                 required=.false., act='store_true', def='false', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--unitcell_projection', switch_ab='-up', &
                 help='Project the phonons on the basis of a unitcell', &
                 required=.false., act='store_true', def='false', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--n_energies', switch_ab='-ne', &
                 help='Number of energies for the axis', &
                 required=.false., act='store', def='1200', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--noisotope', &
                 help='Do not consider isotope scattering.', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
    if (lo_status .ne. 0) stop
    call cli%add(switch='--readiso', &
                 help='Read the isotope distribution from `infile.isotopes`.', &
                 help_markdown='The format is specified [here](../page/files.html#infile.isotopes).', &
                 required=.false., act='store_true', def='.false.', error=lo_status)
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

    call cli%get(switch='--readpath', val=opts%readpathfromfile)
    call cli%get(switch='-nq', val=opts%nq)
    call cli%get(switch='--unit', val=opts%enhet)
    call cli%get(switch='--structure_factor', val=opts%structurefactor, error=lo_status)
    call cli%get(switch='--unitcell_projection', val=opts%projected, error=lo_status)
    call cli%get(switch='--n_energies', val=opts%nf, error=lo_status)
    call cli%get(switch='--qpoint_grid', val=opts%qgrid, error=lo_status)
    call cli%get(switch='--noisotope', val=dumlog)
    opts%isotopescattering = .not. dumlog
    call cli%get(switch='--readiso', val=opts%readiso)

    ! Now we check that the input makes sense
    ! We have to chose between structure factor and projection
    if (opts%projected .and. opts%structurefactor) then
        write(*, *) 'You cannot do unitcell projection and structure factor at the same time'
        stop
    end if
    ! A band path only makes sense for a projection
    if (opts%structurefactor .and. opts%readpathfromfile) then
        write(*, *) 'You cannot specify a path with structure factor'
        stop
    end if

    ! Check if we need the spectral function in real space
    opts%sf_realspace = .false.
    opts%isotopescattering = .false.
    if (opts%isotopescattering .or. opts%thirdorder) then
        opts%sf_realspace = .true.
    end if

    ! The default thing to do is the structure factor
    if (opts%projected .eqv. .false. .and. opts%structurefactor .eqv. .false.) opts%structurefactor = .true.
end subroutine
end module
