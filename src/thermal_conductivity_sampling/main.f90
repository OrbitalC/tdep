#include "precompilerdefinitions"
program thermal_conductivity_sampling
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol, &
                      lo_frequency_THz_to_Hartree, lo_huge, lo_Hartree_to_eV
use gottochblandat, only: walltime, tochar, lo_linspace, lo_logspace, lo_mean, lo_does_file_exist
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
use kappa, only: get_kappa, get_kappa_offdiag, iterative_bte
use new_scattering, only: compute_scattering, lo_scattering_rates
use linewidths, only: compute_linewidths !, self_consistent_linewidths

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
! The scattering rates
type(lo_scattering_rates) :: sr
! Small stuff
real(r8), dimension(:, :), allocatable :: thermal_cond
! timers
real(r8) :: timer_init, timer_scatt, timer_kappa, timer_lw, tt0

! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j, q1, b1
    ! Start MPI and timers
    tt0 = walltime()
    timer_init = tt0
    timer_kappa = 0.0_r8
    call mw%init()
    ! Get options
    call opts%parse()
    if (.not. mw%talk) opts%verbosity = -100
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
    elseif (mw%talk .and. opts%verbosity .gt. 0) then
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
                     sigma=opts%sigma, n_dos_point=opts%mfppts*2, integrationtype=2)

    ! Make sure it's stable, no point in going further if it is unstable.
    if (dr%omega_min .lt. -lo_freqtol) then
        write (*, *) ''
        write (*, *) 'FOUND UNSTABLE MODES. WILL STOP NOW.'
        call mpi_barrier(mw%comm, mw%error)
        call mpi_finalize(lo_status)
        stop
    end if


    ! Make some space to keep intermediate values
    do q1 = 1, qp%n_irr_point
        allocate (dr%iq(q1)%linewidth(dr%n_mode))
        allocate (dr%iq(q1)%F0(3, dr%n_mode))
        allocate (dr%iq(q1)%Fn(3, dr%n_mode))
        allocate (dr%iq(q1)%p_plus(dr%n_mode))
        allocate (dr%iq(q1)%p_minus(dr%n_mode))
        allocate (dr%iq(q1)%p_iso(dr%n_mode))
        allocate (dr%iq(q1)%qs(dr%n_mode))
        allocate (dr%iq(q1)%mfp(3, dr%n_mode))
        allocate (dr%iq(q1)%scalar_mfp(dr%n_mode))
        do b1 = 1, dr%n_mode
            if (dr%iq(q1)%linewidth(b1) .lt. lo_freqtol) then
                dr%iq(q1)%linewidth(b1) = 0.0_r8
            else
                dr%iq(q1)%linewidth(b1) = qp%adaptive_sigma(qp%ip(q1)%radius, dr%iq(q1)%vel(:, b1), &
                                                            dr%default_smearing(b1), opts%sigma) / sqrt(3.0_r8)
            end if
        end do
        dr%iq(q1)%F0 = 0.0_r8
        dr%iq(q1)%Fn = 0.0_r8
        dr%iq(q1)%p_plus = 0.0_r8
        dr%iq(q1)%p_minus = 0.0_r8
        dr%iq(q1)%p_iso = 0.0_r8
        dr%iq(q1)%qs = 0.0_r8
        dr%iq(q1)%mfp = 0.0_r8
        dr%iq(q1)%scalar_mfp = 0.0_r8
        if (opts%fourthorder) then
            allocate (dr%iq(q1)%p_plusplus(dr%n_mode))
            allocate (dr%iq(q1)%p_plusminus(dr%n_mode))
            allocate (dr%iq(q1)%p_minusminus(dr%n_mode))
            dr%iq(q1)%p_plusplus = 0.0_r8
            dr%iq(q1)%p_plusminus = 0.0_r8
            dr%iq(q1)%p_minusminus = 0.0_r8
        end if
    end do
    do q1 = 1, qp%n_full_point
        allocate (dr%aq(q1)%kappa(3, 3, dr%n_mode))
        dr%aq(q1)%kappa = 0.0_r8
    end do

    ! now I have all harmonic things, stop the init timer
    timer_init = walltime() - timer_init
end block initharmonic

scatters: block

    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Calculating scattering events'
    end if
    timer_scatt = walltime()
    call compute_scattering(qp, dr, uc, fct, fcf, opts, mw, mem, sr)
    timer_scatt = walltime() - timer_scatt

    if (mw%talk) write(*, "(1X,A,F12.3,A)") '... done in ', timer_scatt, ' s'
    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Calculating linewidths'
    end if
    timer_lw = walltime()
    call compute_linewidths(qp, dr, sr, opts, mw, mem)
!   if (opts%niter .gt. 0) then
!       call self_consistent_linewidths(dr, qp, sr, opts, mw, mem)
!   end if

    timer_lw = walltime() - timer_lw
    if (mw%talk) write(*, "(1X,A,F12.3,A)") '... done in ', timer_lw, ' s'
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
    if (opts%bteniter .gt. 0) then
        if (mw%talk) then
            write(*, *) ''
            write(*, *) '... solving iterative BTE'
            write (*, "(1X,A4,6(1X,A14),2X,A10)") 'iter', &
                       'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   ', 'DeltaF/F'
        end if
        call iterative_bte(sr, dr, qp, uc, opts%temperature, opts%bteniter, opts%btetol, &
                           opts%isotopescattering, opts%thirdorder, opts%fourthorder, mw, mem)
    end if
    call get_kappa(dr, qp, uc, opts%temperature, kappa)
    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'THERMAL CONDUCTIVITY'
        write (*, *) ''
        write (*, "(1X,A52)") 'Decomposition of the thermal conductivity (in W/m/K)'
        m0 = kappa_sma*lo_kappa_au_to_SI
        write (*, "(1X,A85)") 'Single mode relaxation time approximation (RTA) to Boltzmann transport equation (BTE)'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        m0 = (kappa - kappa_sma)*lo_kappa_au_to_SI
        write (*, "(1X,A73)") 'Correction to full solution of the linearized BTE via iterative procedure'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        m0 = kappa_offdiag*lo_kappa_au_to_SI
        write (*, "(1X,A36)") 'Off diagonal (coherent) contribution'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        m0 = (kappa + kappa_offdiag)*lo_kappa_au_to_SI
        write (*, "(1X,A26)") 'Total thermal conductivity'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
    end if
    kappa = kappa + kappa_offdiag

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

    ! sum up the total time
    if (mw%talk) then
        tt0 = walltime() - tt0
        write(*, *) ''
        write(*, *) '... dumping auxiliary data to files'
        call dr%write_to_hdf5(qp, uc, 'outfile.grid_thermal_conductivity_sampling.hdf5', mem, opts%temperature)
        call lo_dump_gnuplot_2d_real(thermal_cond, 'outfile.thermal_conductivity_sampling', &
                                     ylabel='\kappa W/mK', xlabel='Temperature (K)')

    ! Print timings
        write (*, *) ''
        write (*, '(1X,A21)') 'Suggested citations :'
        write (*, '(1X,A41,A56)') 'Software : ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
        write (*, '(1X,A41,A53)') 'Method : ', 'D. A. Broido et al., Appl Phys Lett 91, 231922 (2007)'
        write (*, '(1X,A41,A43)') 'Iterative Boltzmann transport equation : ', 'M. Omini et al., Phys Rev B 53, 9064 (1996)'
        write (*, '(1X,A41,A49)') 'Algorithm : ', 'A. H. Romero et al., Phys Rev B 91, 214310 (2015)'
        write (*, '(1X,A41,A43)') 'Off diagonal coherent contribution : ', 'L. Isaeva et al., Nat Commun 10 3853 (2019)'
        write (*, '(1X,A41,A46)') 'Sampling method for scattering : ', 'Z. Guo et al., npj Comput Matter 10, 31 (2024)'

        t0 = timer_init + timer_lw + timer_scatt + timer_kappa
        write (*, *) ' '
        write (*, *) 'Timings:'
        write (*, "(A,F12.3,A,F7.3,A)") '            initialization:', timer_init, ' s, ', real(timer_init*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '    scattering computation:', timer_scatt, ' s, ', real(timer_scatt*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '                linewidths:', timer_lw, ' s, ', real(timer_lw*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '                     kappa:', timer_kappa, ' s, ', real(timer_kappa*100/tt0), '%'
        write (*, "(A,F12.3,A)") '                     total:', tt0, ' seconds'
    end if
end block finalize_and_write

! And we are done!
call sr%destroy()
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)

end program
