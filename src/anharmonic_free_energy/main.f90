program anharmonic_free_energy
!!{!src/anharmonic_free_energy/manual.md!}
use konstanter, only: r8, lo_Hartree_to_eV, lo_kb_Hartree, lo_pressure_HartreeBohr_to_GPa, &
                      lo_tol, lo_status
use gottochblandat, only: open_file, walltime, tochar, lo_does_file_exist, open_file, lo_trace
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure

use options, only: lo_opts
use thirdorder, only: free_energy_thirdorder, elastic_thirdorder
use fourthorder, only: free_energy_fourthorder, free_energy_fourthorder_secondorder
use epot, only: lo_energy_differences
use lo_thermodynamic_helpers, only: lo_thermodynamics, lo_symmetrize_stress, lo_full_to_voigt

implicit none

type(lo_opts) :: opts
type(lo_forceconstant_secondorder) :: fc
type(lo_forceconstant_thirdorder) :: fct
type(lo_forceconstant_fourthorder) :: fcf
type(lo_crystalstructure) :: uc
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
type(lo_timer) :: tmr

type(lo_thermodynamics) :: thermo
real(r8), dimension(3, 5) :: cumulant, cumulant_var
real(r8), dimension(3, 3, 5) :: stress_pot, stress_potvar
real(r8) :: timer_init, timer_total
real(r8) :: U0, U1, energy_unit_factor
logical :: havehighorder, havesim

! Init MPI, timers and options
call mw%init()
timer_total = walltime()
timer_init = walltime()
call opts%parse()
call tmr%start()

init: block
    call mw%init()
    ! only be verbose on the first rank
    if (.not. mw%talk) opts%verbosity = -10
    call tmr%start()

    if (mw%talk) then
        write(*, *) 'ANHARMONIC FREE ENERGY'
        write(*, *) 'Recap of the parameters governing the calculation'
        write(*, *) ''
        write(*, '(1X,A40,F20.12)') 'Temperature                             ', opts%temperature
        write(*, '(1X,A40,I4)') 'Number of blocks                        ', opts%nblocks
        write(*, '(1X,A40,L3)') 'Quantum limit                           ', opts%quantum
        write(*, '(1X,A40,L3)') 'Stochastic sampling                     ', opts%stochastic
        write(*, '(1X,A40,I4,I4,I4)') 'Q-point grid Harmonic                   ', opts%qgrid
        write(*, '(1X,A40,I4,I4,I4)') 'Third order q-point grid                ', opts%qg3ph
        write(*, '(1X,A40,I4,I4,I4)') 'Fourth order q-point grid               ', opts%qg4ph
        write(*, *) ''
    end if

    if (mw%talk) write(*, *) 'INITIALIZATION'
    if (mw%talk) write(*, *) '... reading input files'
    ! Read structure
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=.true.)
    if (mw%talk) write (*, *) '... read unitcell'

    ! Read forceconstants
    call fc%readfromfile(uc, 'infile.forceconstant', mem, verbosity=-1)
    if (mw%talk) write (*, *) '... read second order forceconstant'

    if (opts%thirdorder) then
        call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
        if (mw%talk) write (*, *) '... read third order forceconstant'
    end if
    if (opts%fourthorder) then
        call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
        if (mw%talk) write (*, *) '... read fourth order forceconstant'
    end if
    havehighorder = .false.
    if (opts%thirdorder .or. opts%fourthorder) havehighorder = .true.

    call tmr%tock('reading input')

    timer_init = walltime() - timer_init

    thermo%temperature = opts%temperature
    thermo%stochastic = opts%stochastic
    thermo%thirdorder = opts%thirdorder
    thermo%fourthorder = opts%fourthorder
end block init

latdyn: block
    !> The phonon dispersion relation
    type(lo_phonon_dispersions) :: dr
    !> The qpoint mesh
    class(lo_qpoint_mesh), allocatable :: qp
    !> Some stuffs
    real(r8) :: temperature
    !> Some integers for do loops
    integer :: i

    if (mw%talk) then
        write(*, *) ''
        write(*, *) 'HARMONIC CONTRIBUTION'
    end if

    if (mw%talk) write(*, *) '... generating q-mesh'
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=.true., &
                           headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
    if (mw%talk) write(*, *) '... generating harmonic dispersion'
    call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

    ! Little sanity check
    if (dr%omega_min .lt. 0.0_r8) then
        ! Dump the free energies
        if (mw%talk) then
            write (*, *) 'Found negative eigenvalues. Stopping prematurely since no free energy can be defined.'
        end if
        call mw%destroy()
        stop
    end if

    temperature = opts%temperature

    if (mw%talk) write(*, *) '... computing thermodynamic properties with harmonic dispersion'
    ! We start by computing everything we can from the harmonic phonons
    if (opts%quantum) then
        thermo%harmonic%F(1) = dr%phonon_free_energy(temperature)
        thermo%harmonic%S(1) = dr%phonon_entropy(temperature)
        thermo%harmonic%U(1) = thermo%harmonic%F(1) + temperature * thermo%harmonic%S(1)
        thermo%harmonic%Cv(1) = dr%phonon_cv(temperature)
        call dr%phonon_kinetic_stress(qp, uc, temperature, thermo%harmonic%stress(:, :, 1))
    else
        thermo%harmonic%F(1) = dr%phonon_free_energy_classical(temperature)
        thermo%harmonic%U(1) = 3.0_r8 * lo_kb_Hartree * temperature
        thermo%harmonic%S(1) = (thermo%harmonic%U(1) - thermo%harmonic%F(1)) / temperature
        thermo%harmonic%Cv(1) = 3.0_r8 * lo_kb_Hartree
        thermo%harmonic%stress = 0.0_r8
        do i=1, 3
            thermo%harmonic%stress(i, i, 1) = lo_kb_Hartree * temperature * uc%na / uc%volume
        end do
    end if
    ! Usually not needed here, but always a good idea to clean
    call lo_symmetrize_stress(thermo%harmonic%stress(:, :, 1), uc)

    ! If we have third order IFC, might as well compute elastic things
    if (opts%thirdorder) then
        if (mw%talk) write(*, *) '... computing third order contribution to elastic properties'

        call elastic_thirdorder(uc, fc, fct, qp, dr, opts%temperature, thermo%threephonon%stress(:, :, 1), &
                                thermo%alpha, opts%quantum, mw, mem)
        ! Now we symmetrize
        call lo_symmetrize_stress(thermo%threephonon%stress(:, :, 1), uc)
        call lo_symmetrize_stress(thermo%alpha, uc)
        if (opts%stochastic) thermo%first_order%stress = thermo%threephonon%stress
    end if
    call tmr%tock('harmonic properties')
end block latdyn

calcepot: block
    !> The simulation
    type(lo_mdsim) :: sim
    !> The potentiel energy helper
    type(lo_energy_differences) :: pot
    !> The supercell
    type(lo_crystalstructure) :: ss

    if (mw%talk) write(*, *) ''
    if (lo_does_file_exist('infile.sim.hdf5') .or. lo_does_file_exist('infile.meta')) then
    havesim = .true.

    if (mw%talk) write(*, *) 'REAL-SPACE ENERGY CORRECTIONS'

    if (mw%talk) write(*, *) '... reading simulation files'
    ! Read the supercell
    call ss%readfromfile('infile.ssposcar')
    call ss%classify('supercell', uc)
    ! Setup the ifc with the potential
    if (mw%talk) write(*, *) '... preparing potential energy differences'
    call pot%setup(uc, ss, fc, fct, fcf, mw, opts%verbosity+1)

    ! We fetch the simulation data
    if (lo_does_file_exist('infile.sim.hdf5')) then
        call sim%read_from_hdf5('infile.sim.hdf5', verbosity=-1, stride=-1)
    else
        call sim%read_from_file(verbosity=-1, stride=1, magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
    end if

    ! Let's recap info on the simulation
    if (mw%talk) then
            write (*, *) '... short summary of simulation:'
            write (*, *) '                      number of atoms: ', tochar(sim%na)
            write (*, *) '             number of configurations: ', tochar(sim%nt)
            write (*, *) '                thermostat set to (K): ', tochar(sim%temperature_thermostat)
    end if

    if (abs(opts%temperature - sim%temperature_thermostat) .gt. lo_tol .and. mw%talk) then
        write(*, *)
        write(*, *) 'WARNING: input and simulation temperatures do not match'
        write(*, '(1X,A24,F24.12)') 'Input temperature', opts%temperature
        write(*, '(1X,A24,F24.12)') 'Simulation temperature', sim%temperature_thermostat
        write(*, *) 'Results might be garbage !'
        write(*, *)
    end if

    call pot%compute_realspace_thermo(ss, sim, thermo, opts%nblocks, mw, mem)
    ! We also need to symmetrize the stress tensor
    call lo_symmetrize_stress(thermo%first_order%stress(:, :, 1), uc)
    call lo_symmetrize_stress(thermo%first_order%stress(:, :, 2), uc)
    call lo_symmetrize_stress(thermo%second_order%stress(:, :, 1), uc)
    call lo_symmetrize_stress(thermo%second_order%stress(:, :, 2), uc)

    call tmr%tock('simulation')
    else
        havesim = .false.
        if (mw%talk) then
            write(*, *) 'No simulation found, skipping computation of cumulant correction'
        end if
    end if
end block calcepot


latdyn3ph: block
    !> The phonon dispersion relation
    type(lo_phonon_dispersions) :: dr
    !> The qpoint mesh
    class(lo_qpoint_mesh), allocatable :: qp
    !> The third order contribution
    real(r8) :: fe3, s3, cv3


    if (opts%thirdorder) then
        if (mw%talk) then
            write(*, *) ''
            write(*, *) 'THREE PHONONS CONTRIBUTION'
        end if

        if (mw%talk) write(*, *) '... generating q-mesh'
        call lo_generate_qmesh(qp, uc, opts%qg3ph, 'fft', timereversal=.true., &
                               headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        if (mw%talk) write(*, *) '... generating harmonic dispersion'
        call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

        call free_energy_thirdorder(uc, fct, qp, dr, opts%temperature, fe3, s3, cv3, opts%quantum, mw, mem)
        thermo%threephonon%F(1) = fe3
        thermo%threephonon%S(1) = s3
        thermo%threephonon%U(1) = (fe3 + opts%temperature * s3)
        thermo%threephonon%Cv(1) = cv3
        ! If we are in the stochastic mode, we have to add it to the result
        ! It's a second order correction
        if (opts%stochastic) then
            thermo%second_order%F(1) = thermo%second_order%F(1) +thermo%threephonon%F(1)
            thermo%second_order%S(1) = thermo%second_order%S(1) +thermo%threephonon%S(1)
            thermo%second_order%U(1) = thermo%second_order%U(1) +thermo%threephonon%U(1)
            thermo%second_order%Cv(1) = thermo%second_order%Cv(1) + thermo%threephonon%Cv(1)
        end if
    end if
    call tmr%tock('three-phonon')

end block latdyn3ph

latdyn4ph: block
    !> The phonon dispersion relation
    type(lo_phonon_dispersions) :: dr
    !> The qpoint mesh
    class(lo_qpoint_mesh), allocatable :: qp
    !> The third order contribution
    real(r8) :: fe4, s4, cv4

    if (opts%fourthorder) then
        if (mw%talk) then
            write(*, *) ''
            write(*, *) 'FOUR PHONONS CONTRIBUTION'
        end if

        if (mw%talk) write(*, *) '... generating q-mesh'
        call lo_generate_qmesh(qp, uc, opts%qg4ph, 'fft', timereversal=.true., &
                               headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        if (mw%talk) write(*, *) '... generating harmonic dispersion'
        call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

        call free_energy_fourthorder(uc, fcf, qp, dr, opts%temperature, fe4, s4, cv4, opts%quantum, mw, mem)
        thermo%fourphonon%F = fe4
        thermo%fourphonon%S = s4
        thermo%fourphonon%U = (fe4 + opts%temperature * s4)
        thermo%fourphonon%Cv = cv4
        ! If we are in the stochastic mode, we have to add it to the result
        ! It's a first order correction
        if (opts%stochastic) then
            thermo%first_order%F(1) = thermo%first_order%F(1) + thermo%fourphonon%F(1)
            thermo%first_order%U(1) = thermo%first_order%U(1) + thermo%fourphonon%U(1)
            thermo%first_order%S(1) = thermo%first_order%S(1) + thermo%fourphonon%S(1)
            thermo%first_order%Cv(1) = thermo%first_order%Cv(1) + thermo%fourphonon%Cv(1)
        end if

        ! TODO set fourthorder, second order cumulant its  own qpoint grid density, too expensive otherwise
!       call free_energy_fourthorder_secondorder(uc, fcf, qp, dr, opts%temperature, fe4, s4, cv4, opts%quantum, mw, mem)
    end if
    call tmr%tock('four-phonon')
end block latdyn4ph

summary: block
    !> Unit conversion factor
    real(r8) :: f_unit, e_unit, s_unit, c_unit, p_unit
    !> Buffer for the stress tensor
    real(r8), dimension(3, 3, 2) :: sigma
    !> Some buffer to print the results
    real(r8), dimension(4) :: buf, var
    !> The pressure
    real(r8), dimension(2) :: pressure
    !> To have a pretty logfile
    character(len=1000) :: opfc, opff, opfs
    !> A tolerance to clean-up stress results
    real(r8) :: stol
    !> Some integers
    integer :: i, u

    f_unit = lo_Hartree_to_eV
    e_unit = lo_Hartree_to_eV
    s_unit = 1.0 / lo_kb_Hartree
    c_unit = 1.0 / lo_kb_Hartree
    p_unit = lo_pressure_HartreeBohr_to_GPa

    ! Get the stress tensor
    sigma = (thermo%harmonic%stress + thermo%first_order%stress + thermo%second_order%stress) * p_unit
    ! And the pressure
    pressure(1) = lo_trace(sigma(:, :, 1)) / 3.0_r8
    pressure(2) = lo_trace(sigma(:, :, 2)) / 3.0_r8

    if (mw%talk) then
        u = open_file('out', 'outfile.anharmonic_free_energy')
        write(u, '(A2,A12,8X,E20.12)') '# ', 'Temperature:', opts%temperature

        opfc = '(4(1X,A24))'
        opff = '(4(1X,F24.12))'
        opfs = '(3(1X,F24.12))'
        buf(1) = thermo%harmonic%F(1) * f_unit
        buf(2) = thermo%harmonic%U(1) * e_unit
        buf(3) = thermo%harmonic%S(1) * s_unit
        buf(4) = thermo%harmonic%Cv(1) * c_unit
        var = 0.0_r8
        write(*, *) ''
        write(*, *) 'SUMMARY OF RESULTS'
        write(*, *) ''
        write(*, *) 'Effective harmonic contribution: F = F_harm'
        write(*, opfc) 'Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(*, opff) buf
        write(*, opff) var
        write(u, *) '# Effective harmonic contribution: F = F_harm'
        write(u, opfc) '# Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(u, opff) buf
        write(u, opff) var

        ! Update the values
        buf(1) = buf(1) + thermo%first_order%F(1) * f_unit
        buf(2) = buf(2) + thermo%first_order%U(1) * e_unit
        buf(3) = buf(3) + thermo%first_order%S(1) * s_unit
        buf(4) = buf(4) + thermo%first_order%Cv(1) * c_unit
        ! And the uncertainty
        var(1) = var(1) + thermo%first_order%F(2) * f_unit
        var(2) = var(2) + thermo%first_order%U(2) * e_unit
        var(3) = var(3) + thermo%first_order%S(2) * s_unit
        var(4) = var(4) + thermo%first_order%Cv(2) * c_unit
        write(*, *) ''
        write(*, *) 'Gibbs-Bogoliubov approximation: F = F_harm + <V-V_harm>'
        write(*, opfc) 'Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(*, opff) buf
        write(*, opff) var
        write(u, *) '# Gibbs-Bogoliubov approximation: F = F_harm + <V-V_harm>'
        write(u, opfc) '# Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(u, opff) buf
        write(u, opff) var

        ! Update the values
        buf(1) = buf(1) + thermo%second_order%F(1) * f_unit
        buf(2) = buf(2) + thermo%second_order%U(1) * e_unit
        buf(3) = buf(3) + thermo%second_order%S(1) * s_unit
        buf(4) = buf(4) + thermo%second_order%Cv(1) * c_unit
        ! And the uncertainty
        var(1) = var(1) + thermo%second_order%F(2) * f_unit
        var(2) = var(2) + thermo%second_order%U(2) * e_unit
        var(3) = var(3) + thermo%second_order%S(2) * s_unit
        var(4) = var(4) + thermo%second_order%Cv(2) * c_unit
        write(*, *) ''
        write(*, *) 'Second order cumulant approximation to the free energy'
        write(*, opfc) 'Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(*, opff) buf
        write(*, opff) var
        write(u, *) '# Second order cumulant approximation to the free energy'
        write(u, opfc) '# Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(u, opff) buf
        write(u, opff) var

        write(*, *) ''
        write(*, *) 'Elastic properties'
        write(*, '(2(1X,A24))') 'Pressure [GPa]', 'Uncertainty [GPa]'
        write(*, '(2(1X,F24.12))') -pressure(1), pressure(2)
        write(*, *) 'Stress tensor [GPa]'
        do i=1, 3
            write(*, opfs) sigma(i, :, 1)
        end do
        write(*, *) 'Uncertainty [GPa]'
        do i=1, 3
            write(*, opfs) sigma(i, :, 2)
        end do
        close(u)
    end if

    call tmr%stop()
    if (mw%talk) write(*, *) ''
    call tmr%dump(mw, 'Timings:')
end block summary

! And we are done!
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)
end program
