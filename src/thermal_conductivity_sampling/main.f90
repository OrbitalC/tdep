#include "precompilerdefinitions"
program thermal_conductivity_sampling
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol
use gottochblandat, only: walltime, tochar, lo_linspace, lo_logspace, lo_mean
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
use dump_data, only: lo_dump_gnuplot_2d_real
use lo_randomnumbers, only: lo_mersennetwister
!
use options, only: lo_opts
use kappa, only: get_kappa, get_kappa_offdiag
use scattering, only: compute_scattering_isotopes, compute_scattering_threephonon, &
                      compute_scattering_fourphonon

implicit none
! Standard from libolle
type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_forceconstant_fourthorder) :: fcf
type(lo_phonon_dispersions) :: dr
type(lo_phonon_dos) :: pd
type(lo_crystalstructure) :: uc
class(lo_qpoint_mesh), allocatable :: qp
type(lo_mersennetwister) :: rng
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
! Small stuff
real(r8), dimension(:, :), allocatable :: thermal_cond
! timers
real(r8) :: timer_init, timer_scatt, timer_kappa, timer_count, tt0

! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j
    ! Start MPI and timers
    tt0 = walltime()
    timer_init = tt0
    timer_kappa = 0.0_r8
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
    if (opts%fourthorder) then
        call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
        if (mw%talk) write (*, *) '... read fourth order forceconstant'
    end if

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
    integer :: i
    real(r8) :: t0
    t0 = walltime()
    timer_scatt = walltime()

    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Calculating scattering events'
    end if

    ! Make some space to keep intermediate values
    do i = 1, qp%n_irr_point
        allocate (dr%iq(i)%linewidth(dr%n_mode))
        allocate (dr%iq(i)%F0(3, dr%n_mode))
        allocate (dr%iq(i)%Fn(3, dr%n_mode))
        allocate (dr%iq(i)%mfp(3, dr%n_mode))
        allocate (dr%iq(i)%scalar_mfp(dr%n_mode))
        dr%iq(i)%linewidth = 0.0_r8
        dr%iq(i)%F0 = 0.0_r8
        dr%iq(i)%Fn = 0.0_r8
        dr%iq(i)%mfp = 0.0_r8
        dr%iq(i)%scalar_mfp = 0.0_r8
    end do
    do i = 1, qp%n_full_point
        allocate (dr%aq(i)%kappa(3, 3, dr%n_mode))
        dr%aq(i)%kappa = 0.0_r8
    end do
    ! First we compute the scattering strength
    ! We start by initializing the random number generator
    call rng%init(iseed=mw%r, rseed=walltime())
    ! For the isotopes
    if (opts%isotopescattering) then
        call compute_scattering_isotopes(qp, dr, uc, opts%integrationtype, opts%sigma, opts%thres, mw, mem)
    end if
    ! For the third order
    call compute_scattering_threephonon(qp, dr, fct, opts%temperature, opts%ratio3ph, opts%integrationtype, opts%sigma, opts%thres, rng, mw, mem)
    ! For the fourth order
    if (opts%fourthorder) then
        call compute_scattering_fourphonon(qp, dr, fcf, opts%temperature, opts%ratio4ph, opts%integrationtype, opts%sigma, opts%thres, rng, mw, mem)
    end if

    ! stop counting timer, start matrixelement timer
    timer_scatt = walltime() - timer_scatt
end block scatters

kappa: block
    real(r8), dimension(3, 3) :: kappa, kappa_offdiag, kappa_sma, m0
    real(r8) :: t0
    integer :: i

    timer_kappa = walltime()

    ! space to store the actual thermal conductivity
    allocate (thermal_cond(10, 1))
    thermal_cond = 0.0_r8

    ! I might get a silly tiny temperature, then things will break.
    if (opts%temperature .lt. lo_temperaturetol) then
        kappa = 0.0_r8
        kappa_sma = 0.0_r8
        kappa_offdiag = 0.0_r8
        thermal_cond(1, 1) = opts%temperature
        thermal_cond(2:10, 1) = 0.0_r8
    end if

    ! Scattering rates
    t0 = walltime()

    call get_kappa(dr, qp, uc, opts%temperature, kappa_sma)
    call get_kappa_offdiag(dr, qp, uc, fc, opts%temperature, mem, mw, kappa_offdiag)
    kappa = kappa_sma + kappa_offdiag
    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'THERMAL CONDUCTIVITY'
        m0 = kappa_sma*lo_kappa_au_to_SI
        write (*, *) ''
        write (*, "(1X,A52)") 'Decomposition of the thermal conductivity (in W/m/K)'
        write (*, "(1X,A85)") 'Single mode relaxation time approximation (RTA) to Boltzmann transport equation (BTE)'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4),2X,E10.3)") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        m0 = kappa_offdiag*lo_kappa_au_to_SI
        write (*, "(1X,A36)") 'Off diagonal (coherent) contribution'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4),2X,E10.3)") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        m0 = kappa*lo_kappa_au_to_SI
        write (*, "(1X,A26)") 'Total thermal conductivity'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4),2X,E10.3)") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
    end if

    ! Store thermal conductivity tensor
    thermal_cond(1, 1) = opts%temperature
    thermal_cond(2, 1) = kappa(1, 1)*lo_kappa_au_to_SI
    thermal_cond(3, 1) = kappa(2, 2)*lo_kappa_au_to_SI
    thermal_cond(4, 1) = kappa(3, 3)*lo_kappa_au_to_SI
    thermal_cond(5, 1) = kappa(1, 3)*lo_kappa_au_to_SI
    thermal_cond(6, 1) = kappa(2, 3)*lo_kappa_au_to_SI
    thermal_cond(7, 1) = kappa(1, 2)*lo_kappa_au_to_SI
    thermal_cond(8, 1) = kappa(3, 1)*lo_kappa_au_to_SI
    thermal_cond(9, 1) = kappa(3, 2)*lo_kappa_au_to_SI
    thermal_cond(10, 1) = kappa(2, 1)*lo_kappa_au_to_SI

    timer_kappa = walltime() - timer_kappa
end block kappa

finalize_and_write: block
    real(r8) :: t0
    ! Write thermal conductivity to file
    if (mw%talk) call lo_dump_gnuplot_2d_real(thermal_cond, 'outfile.thermal_conductivity_sampling', &
                                              ylabel='\kappa W/mK', xlabel='Temperature (K)')

    ! sum up the total time
    if (mw%talk) tt0 = walltime() - tt0

    ! Print timings
    if (mw%talk) then
        write (*, *) ''
        write (*, '(1X,A21)') 'Suggested citations :'
        write (*, '(1X,A41,A56)') 'Software : ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
        write (*, '(1X,A41,A53)') 'Method : ', 'D. A. Broido et al., Appl Phys Lett 91, 231922 (2007)'
        write (*, '(1X,A41,A43)') 'Iterative Boltzmann transport equation : ', 'M. Omini et al., Phys Rev B 53, 9064 (1996)'
        write (*, '(1X,A41,A49)') 'Algorithm : ', 'A. H. Romero et al., Phys Rev B 91, 214310 (2015)'
        write (*, '(1X,A41,A43)') 'Off diagonal coherent contribution : ', 'L. Isaeva et al., Nat Commun 10 3853 (2019)'

        t0 = timer_init + timer_count + timer_kappa
        write (*, *) ' '
        write (*, *) 'Timings:'
        write (*, "(A,F12.3,A,F7.3,A)") '            initialization:', timer_init, ' s, ', real(timer_init*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '    scattering computation:', timer_scatt, ' s, ', real(timer_scatt*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '                     kappa:', timer_kappa, ' s, ', real(timer_kappa*100/tt0), '%'
        write (*, "(A,F12.3,A)") '                     total:', tt0, ' seconds'
    end if
end block finalize_and_write

! And we are done!
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)

end program
