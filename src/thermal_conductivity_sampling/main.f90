#include "precompilerdefinitions"
program thermal_conductivity_sampling
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol
use gottochblandat, only: walltime, tochar, lo_linspace, lo_logspace, lo_mean
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use dump_data, only: lo_dump_gnuplot_2d_real


implicit none
! Standard from libolle
type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_phonon_dispersions) :: dr
type(lo_phonon_dos) :: pd
type(lo_crystalstructure) :: uc
class(lo_qpoint_mesh), allocatable :: qp
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
! Small stuff
real(r8), dimension(:, :), allocatable :: thermal_cond
real(r8), dimension(:), allocatable :: temperatures
! timers
real(r8) :: timer_init, timer_count, timer_matrixelements, timer_scf
real(r8) :: timer_kappa, timer_qs, timer_cumulative, tt0

! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j
    ! Start MPI and timers
    tt0 = walltime()
    timer_init = tt0
    timer_qs = 0.0_r8
    timer_kappa = 0.0_r8
    timer_scf = 0.0_r8
    timer_cumulative = 0.0_r8
    call mw%init()
    ! Get options
    call opts%parse()
    if (mw%r .ne. 0) opts%verbosity = -100
    ! Init memory tracker
    call mem%init()

    if (mw%talk) write (*, *) '... using ', tochar(mw%n), ' MPI ranks'
    ! There is a bunch of stuff that all ranks need, first the unit cell:
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=opts%timereversal)
    if (mw%talk) write (*, *) '... read unitcell poscar'

    ! Perhaps non-natural isotope distribution
    ! write (*, *) 'FIXME OUTPUT UNITS'
    if (opts%readiso) then
        if (mw%talk) write (*, *) '... reading isotope distribution from file'
        call uc%readisotopefromfile()
        if (mw%talk) then
            do i = 1, uc%na
                do j = 1, uc%isotope(i)%n
                    write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                        j, uc%isotope(i)%conc(j), uc%isotope(i)%mass(j)
                end do
                write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                    trim(uc%atomic_symbol(uc%species(i))), uc%isotope(i)%mean_mass, &
                    uc%isotope(i)%disorderparameter
            end do
        end if
    elseif (opts%verbosity .gt. 0) then
        do i = 1, uc%na
            do j = 1, uc%isotope(i)%n
                write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                    j, uc%isotope(i)%conc(j), uc%isotope(i)%mass(j)
            end do
            write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                trim(uc%atomic_symbol(uc%species(i))), uc%isotope(i)%mean_mass, &
                uc%isotope(i)%disorderparameter
        end do
    end if

    ! Read the force constants
    call fc%readfromfile(uc, 'infile.forceconstant', mem, -1)
    if (mw%talk) write (*, *) '... read second order forceconstant'
    call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (mw%talk) write (*, *) '... read third order forceconstant'

    ! Get q-point mesh
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, &
                           headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity, nosym=.not. opts%qpsymmetry)

    ! Get frequencies in the whole BZ
    if (mw%talk) then
        write (*, *) '... getting the full dispersion relations'
    end if
    call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)
    ! Also the phonon DOS, for diagnostics
    call pd%generate(dr, qp, uc, mw, mem, verbosity=opts%verbosity, &
                     sigma=opts%sigma, n_dos_point=opts%mfppts*2, integrationtype=opts%integrationtype)

    ! Make sure it's stable, no point in going further if it is unstable.
    if (dr%omega_min .lt. -lo_freqtol) then
        write (*, *) ''
        write (*, *) 'FOUND UNSTABLE MODES. WILL STOP NOW.'
        call mpi_barrier(mw%comm, mw%error)
        call mpi_finalize(lo_status)
        stop
    end if

    ! now I have all harmonic things, stop the init timer
    timer_init = walltime() - timer_init
end block initharmonic


scatters: block
    real(r8) :: t0
    t0 = walltime()
    timer_count = walltime()
    ! Make some space to keep intermediate values
    do i = 1, qp%n_irr_point
        allocate (dr%iq(i)%p_plus(dr%n_mode))
        allocate (dr%iq(i)%p_minus(dr%n_mode))
        allocate (dr%iq(i)%p_iso(dr%n_mode))
        allocate (dr%iq(i)%linewidth(dr%n_mode))
    end do
    ! First we compute the scattering strength
    ! For the isotopes
    call compute_scattering_isotopes(qp, dr, uc, mw, mem)
    ! For the third order
    call compute_scattering_threephonon(qp, dr, fct, opts%ratio3ph, mw, mem)

    ! stop counting timer, start matrixelement timer
    timer_count = walltime() - timer_count
    timer_matrixelements = walltime()
end block scatters

kappa: block
    integer :: i

    ! space to store the actual thermal conductivity
    allocate(thermal_cond(10, opts%trangenpts))
    thermal_cond = 0.0_r8

    ! temperature axis
    allocate (temperatures(opts%trangenpts))
    if (opts%logtempaxis) then
        call lo_logspace(opts%trangemin, opts%trangemax, temperatures)
    else
        call lo_linspace(opts%trangemin, opts%trangemax, temperatures)
    end if
end block kappa
end program
