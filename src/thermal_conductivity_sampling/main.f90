#include "precompilerdefinitions"
program thermal_conductivity_sampling
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol, &
                      lo_huge, lo_m_to_Bohr
use gottochblandat, only: walltime, tochar, open_file
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_phonon_dos, only: lo_phonon_dos
!
use options, only: lo_opts
use kappa, only: get_kappa, get_kappa_offdiag, iterative_bte, compute_qs, &
                 symmetrize_kappa
use scattering, only: lo_scattering_rates

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
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
! The scattering rates
type(lo_scattering_rates) :: sr
! timers
real(r8) :: timer_init, timer_scatt, timer_kappa, tt0


! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j, q1
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

    if (mw%talk) then
        write(*, *) 'Recap of the parmeters governing the calculation'
        write(*, '(1X,A40,F20.12)') 'Temperature                             ', opts%temperature
        write(*, '(1X,A40,L3)') 'Thirdorder scattering                   ', opts%thirdorder
        write(*, '(1X,A40,L3)') 'Fourthorder scattering                  ', opts%fourthorder
        write(*, '(1X,A40,L3)') 'Isotope scattering                      ', opts%isotopescattering
        write(*, '(1X,A40,I4,I4,I4)') 'full q-point grid                       ', opts%qgrid
        write(*, '(1X,A40,I4,I4,I4)') 'Monte-Carlo 3rd order q-point grid      ', opts%qg3ph
        write(*, '(1X,A40,I4,I4,I4)') 'Monte-Carlo 4th order q-point grid      ', opts%qg4ph
        write(*, '(1X,A40,I5)') 'Max number of iteration                 ', opts%scfiterations
        write(*, '(1X,A40,E20.12)') 'Max mean free path (in m)               ', opts%mfp_max / lo_m_to_Bohr
        write(*, '(1X,A40,E20.12)') 'Tolerance for the iterative BTE         ', opts%btetol
        select case (opts%integrationtype)
        case(1)
            write(*, '(1X,A40,2X,A)') 'Integration type                        ', 'Gaussian with fixed broadening'
            write(*, '(1X,A40,E20.12)') 'Broadening parameter                    ', opts%sigma
        case(2)
            write(*, '(1X,A40,2X,A)') 'Integration type                        ', 'Adaptive Gaussian'
        end select
        write(*, '(1X,A40,I4)') 'Number of MPI ranks                     ', mw%n
        write(*, *) ''
    end if

    if (mw%talk) write(*, *) 'Initialize calculation'
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

    if (mw%talk) write(*, *) '... generating q-point mesh'
    ! Get q-point mesh
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=opts%timereversal, &
                           headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity, nosym=.not. opts%qpsymmetry)

    ! Get frequencies in the whole BZ
    if (mw%talk) write(*, *) '... generating harmonic properties on the q-point mesh'
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
        allocate (dr%iq(q1)%qs(dr%n_mode))
        allocate (dr%iq(q1)%mfp(3, dr%n_mode))
        allocate (dr%iq(q1)%scalar_mfp(dr%n_mode))
        dr%iq(q1)%linewidth = 0.0_r8
        dr%iq(q1)%F0 = 0.0_r8
        dr%iq(q1)%Fn = 0.0_r8
        dr%iq(q1)%qs = 0.0_r8
        dr%iq(q1)%mfp = 0.0_r8
        dr%iq(q1)%scalar_mfp = 0.0_r8

        allocate(dr%iq(q1)%kappa(3, 3, dr%n_mode))
        dr%iq(q1)%kappa = 0.0_r8
    end do
    do q1 = 1, qp%n_full_point
        allocate (dr%aq(q1)%kappa(3, 3, dr%n_mode))
        dr%aq(q1)%kappa = 0.0_r8
    end do

    ! now I have all harmonic things, stop the init timer
    timer_init = walltime() - timer_init
    if (mw%talk) write(*, "(1X,A,F12.3,A)") '... done in ', timer_init, ' s'
end block initharmonic

scatters: block

    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Calculating scattering events'
    end if
    timer_scatt = walltime()
    call sr%generate(qp, dr, uc, fct, fcf, opts, mw, mem)
    timer_scatt = walltime() - timer_scatt

    if (mw%talk) write(*, "(1X,A,F12.3,A)") '... done in ', timer_scatt, ' s'
end block scatters

kappa: block
    real(r8), dimension(3, 3) :: kappa_bte, kappa_offdiag, kappa_sma, m0
    real(r8) :: t0
    integer :: i, u

    timer_kappa = walltime()

    ! space to store the actual thermal conductivity

    ! I might get a silly tiny temperature, then things will break.
    if (opts%temperature .lt. lo_temperaturetol) then
        kappa_bte = 0.0_r8
        kappa_sma = 0.0_r8
        kappa_offdiag = 0.0_r8
    end if

    ! Scattering rates
    t0 = walltime()

    call compute_qs(dr, qp, opts%temperature)
    call get_kappa(dr, qp, uc, opts%temperature, kappa_sma)
    call get_kappa_offdiag(dr, qp, uc, fc, opts%temperature, mem, mw, kappa_offdiag)
    if (opts%scfiterations .gt. 0) then
        if (mw%talk) then
            write(*, *) ''
            write(*, *) 'Solving iterative BTE'
            write (*, "(1X,A4,6(1X,A14),2X,A10)") 'iter', &
                       'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   ', 'DeltaF/F'
        end if
        call iterative_bte(sr, dr, qp, uc, opts%temperature, opts%scfiterations, opts%btetol, mw, mem)
    end if
    call get_kappa(dr, qp, uc, opts%temperature, kappa_bte)
    if (mw%talk) write(*, *) ''
    if (mw%talk) write(*, *) '... symmetrizing the thermal conductivity tensors'
    call symmetrize_kappa(kappa_bte, uc)
    call symmetrize_kappa(kappa_offdiag, uc)
    call symmetrize_kappa(kappa_sma, uc)
    if (mw%talk) then
        ! First we write in the standard output
        u=open_file('out', 'outfile.thermal_conductivity_sampling')
        write (u, '(A2,A5,15X,A)') '# ', 'Unit:', 'W/m/K'
        write (u, '(A2,A12,8X,E20.12)') '# ', 'Temperature:', opts%temperature

        write (*, *) ''
        write (*, *) 'THERMAL CONDUCTIVITY'
        write (*, *) ''
        write (*, "(1X,A52)") 'Decomposition of the thermal conductivity (in W/m/K)'
        m0 = kappa_sma*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A85)") 'Single mode relaxation time approximation (RTA) to Boltzmann transport equation (BTE)'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A43)") '# Single mode relaxation time approximation'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        m0 = (kappa_bte - kappa_sma)*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A73)") 'Correction to full solution of the linearized BTE via iterative procedure'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A25)") '# Collective contribution'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        m0 = kappa_offdiag*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A36)") 'Off diagonal (coherent) contribution'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A36)") '# Off diagonal coherent contribution'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)

        m0 = (kappa_bte + kappa_offdiag)*lo_kappa_au_to_SI
        ! First in the standard output
        write (*, "(1X,A26)") 'Total thermal conductivity'
        write (*, "(1X,A4,6(1X,A14))") '', 'kxx   ', 'kyy   ', 'kzz   ', 'kxy   ', 'kxz   ', 'kyz   '
        write (*, "(5X,6(1X,F14.4))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
        ! Then in the outfile
        write (u, "(A28)") '# Total thermal conductivity'
        write (u, "(A1,6(1X,A24))") '#', 'kxx', 'kyy', 'kzz', 'kxy', 'kxz', 'kyz'
        write (u, "(1X,6(1X,E24.12))") m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
    end if

    timer_kappa = walltime() - timer_kappa
end block kappa

finalize_and_write: block
    integer :: u
    real(r8) :: t0
    ! Write thermal conductivity to file

    ! sum up the total time
    if (mw%talk) then
        tt0 = walltime() - tt0
        write(*, *) ''
        write(*, *) '... dumping auxiliary data to files'
        call dr%write_to_hdf5(qp, uc, 'outfile.grid_thermal_conductivity_sampling.hdf5', mem, opts%temperature)

        write(*, *) ''
        write(*, '(A61,A)') 'Scattering rates can be found in                             ', 'outfile.grid_thermal_conductivity_sampling.hdf5'
        write(*, '(A61,A)') 'Thermal conductivity tensor can be found in                  ', 'outfile.thermal_conductivity_sampling'

    ! Print timings
        write (*, *) ''
        write (*, '(1X,A21)') 'Suggested citations :'
        write (*, '(1X,A41,A56)') 'Software : ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
        write (*, '(1X,A41,A53)') 'Method : ', 'D. A. Broido et al., Appl Phys Lett 91, 231922 (2007)'
        write (*, '(1X,A41,A43)') 'Iterative Boltzmann transport equation : ', 'M. Omini et al., Phys Rev B 53, 9064 (1996)'
        write (*, '(1X,A41,A49)') 'Algorithm : ', 'A. H. Romero et al., Phys Rev B 91, 214310 (2015)'
        write (*, '(1X,A41,A43)') 'Off diagonal coherent contribution : ', 'L. Isaeva et al., Nat Commun 10 3853 (2019)'
        write (*, '(1X,A41,A46)') 'Sampling method for scattering : ', 'Z. Guo et al., npj Comput Matter 10, 31 (2024)'

        t0 = timer_init + timer_scatt + timer_kappa
        write (*, *) ' '
        write (*, *) 'Timings:'
        write (*, "(A,F12.3,A,F7.3,A)") '            initialization:', timer_init, ' s, ', real(timer_init*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '    scattering computation:', timer_scatt, ' s, ', real(timer_scatt*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '                     kappa:', timer_kappa, ' s, ', real(timer_kappa*100/tt0), '%'
        write (*, "(A,F12.3,A)") '                     total:', tt0, ' seconds'
    end if
end block finalize_and_write

! And we are done!
call sr%destroy()
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)

end program
