program anharmonic_free_energy
!!{!src/anharmonic_free_energy/manual.md!}
use konstanter, only: r8, lo_Hartree_to_eV, lo_kb_Hartree, lo_pressure_HartreeBohr_to_GPa, lo_tol, lo_status, lo_freqtol, lo_twopi
use gottochblandat, only: open_file, walltime, lo_linspace, lo_progressbar_init, lo_progressbar, tochar, &
                          lo_does_file_exist, lo_mean, lo_stddev, lo_harmonic_oscillator_internal_energy, lo_chop, &
                          lo_invert_real_matrix, lo_harmonic_oscillator_cv
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh, lo_read_qmesh_from_file, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_phonon_bandstructure_on_path, only: lo_phonon_bandstructure
use type_phonon_dos, only: lo_phonon_dos

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
        write(*, *) 'Recap of the parameters governing the calculation'
        write(*, *) ''
        write (*, '(1X,A40,F20.12)') 'Temperature                             ', opts%temperature
        write (*, '(1X,A40,L3)') 'Classical limit                         ', opts%quantum
        write (*, '(1X,A40,I4,I4,I4)') 'Q-point grid Harmonic                   ', opts%qgrid
        write (*, '(1X,A40,I4,I4,I4)') 'Third order q-point grid                ', opts%qg3ph
        write (*, '(1X,A40,I4,I4,I4)') 'Fourth order q-point grid               ', opts%qg4ph
        write(*, *) ''
    end if

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
    call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=.true., headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
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
        thermo%f0 = dr%phonon_free_energy(temperature)
        thermo%s0 = dr%phonon_entropy(temperature)
        thermo%u0 = thermo%f0 + temperature * thermo%s0
        thermo%cv0 = dr%phonon_cv(temperature)
        call dr%phonon_kinetic_stress(qp, uc, temperature, thermo%stress_kin)
    else
        thermo%f0 = dr%phonon_free_energy_classical(temperature)
        thermo%u0 = 3.0_r8 * lo_kb_Hartree * temperature
        thermo%s0 = (thermo%u0 - thermo%f0) / temperature
        thermo%cv0 = 3.0_r8 * lo_kb_Hartree
        thermo%stress_kin = 0.0_r8
        do i=1, 3
            thermo%stress_kin(i, i) = lo_kb_Hartree * temperature * uc%na / uc%volume
        end do
    end if
    ! Usually not needed here, but always a good idea to clean
    call lo_symmetrize_stress(thermo%stress_kin, uc)
    thermo%stress_kin = lo_chop(thermo%stress_kin, sum(abs(thermo%stress_kin))*1e-6_r8)

    ! If we have third order IFC, might as well compute elastic things
    if (opts%thirdorder) then
        if (mw%talk) write(*, *) '... computing third order contribution to elastic properties'

        call elastic_thirdorder(uc, fc, fct, qp, dr, opts%temperature, thermo%stress_3ph, thermo%alpha, opts%quantum, mw, mem)
        ! Now we symmetrize and store everything
        call lo_symmetrize_stress(thermo%stress_3ph, uc)
        thermo%stress_3ph = lo_chop(thermo%stress_3ph, sum(abs(thermo%stress_3ph))*1e-6_r8)
        call lo_symmetrize_stress(thermo%alpha, uc)
        thermo%alpha = lo_chop(thermo%alpha, sum(abs(thermo%alpha))*1e-6_r8)
    end if
    call tmr%tock('harmonic properties')
end block latdyn

calcepot: block
    type(lo_mdsim) :: sim
    type(lo_energy_differences) :: pot
    type(lo_crystalstructure) :: ss

    real(r8), dimension(:, :, :, :), allocatable :: sdiff
    real(r8), dimension(:, :, :), allocatable :: buf_stress
    real(r8), dimension(:, :), allocatable :: buf
    real(r8), dimension(:, :), allocatable :: ediff
    real(r8), dimension(3, 3) :: s3, sk2, sk3, sk4, skp

    real(r8) :: e2, e3, e4, ep, inverse_kbt, f0
    integer :: i, j, a, b, blocksize
    real(r8) :: stol

    if (mw%talk) write(*, *) ''
    if (lo_does_file_exist('infile.sim.hdf5') .or. lo_does_file_exist('infile.meta')) then
    havesim = .true.

    if (mw%talk) write(*, *) 'Computing energy differences'

    call ss%readfromfile('infile.ssposcar')
    call ss%classify('supercell', uc)
    call pot%setup(uc, ss, fc, fct, fcf, mw, opts%verbosity+1)

    ! We fetch the simulation data
    if (lo_does_file_exist('infile.sim.hdf5')) then
        call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=-1)
    else
        call sim%read_from_file(verbosity=opts%verbosity + 2, stride=1, magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
    end if

    if (abs(opts%temperature - sim%temperature_thermostat) .gt. lo_tol .and. mw%talk) then
        write(*, *)
        write(*, *) 'WARNING: input and simulation temperatures do not match'
        write(*, '(1X,A24,F24.12)') 'Input temperature', opts%temperature
        write(*, '(1X,A24,F24.12)') 'Simulation temperature', sim%temperature_thermostat
        write(*, *) 'Results might be garbage !'
        write(*, *)
    end if

    if (thermo%temperature .gt. 1E-5_r8) then
        inverse_kbt = 1.0_r8/lo_kb_Hartree/thermo%temperature
    else
        inverse_kbt = 0.0_r8
    end if

    if (mw%talk) write(*, *) '... computing energy differences'
    call pot%compute_realspace_thermo(ss, sim, thermo, opts%nblocks, mw)

    stol = sum(abs(thermo%stress_pot)) * 1e-6_r8
    call lo_symmetrize_stress(thermo%stress_pot, uc)
    thermo%stress_pot = lo_chop(thermo%stress_pot, stol)
    call lo_symmetrize_stress(thermo%stress_potvar, uc)
    thermo%stress_potvar = lo_chop(thermo%stress_potvar, stol)
    call lo_symmetrize_stress(thermo%stress_diff, uc)
    thermo%stress_diff = lo_chop(thermo%stress_diff, stol)
    call lo_symmetrize_stress(thermo%stress_diffvar, uc)
    thermo%stress_diffvar = lo_chop(thermo%stress_diffvar, stol)

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
        call lo_generate_qmesh(qp, uc, opts%qg3ph, 'fft', timereversal=.true., headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        if (mw%talk) write(*, *) '... generating harmonic dispersion'
        call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

        call free_energy_thirdorder(uc, fct, qp, dr, opts%temperature, fe3, s3, cv3, opts%quantum, mw, mem)
        thermo%f3 = fe3
        thermo%s3 = s3
        thermo%u3 = (fe3 + opts%temperature * s3)
        thermo%cv3 = cv3
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
        call lo_generate_qmesh(qp, uc, opts%qg4ph, 'fft', timereversal=.true., headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
        if (mw%talk) write(*, *) '... generating harmonic dispersion'
        call dr%generate(qp, fc, uc, mw=mw, mem=mem, verbosity=opts%verbosity)

        call free_energy_fourthorder(uc, fcf, qp, dr, opts%temperature, fe4, s4, cv4, opts%quantum, mw, mem)
        thermo%f4 = fe4
        thermo%s4 = s4
        thermo%u4 = (fe4 + opts%temperature * s4)
        thermo%cv4 = cv4

        call free_energy_fourthorder_secondorder(uc, fcf, qp, dr, opts%temperature, fe4, s4, cv4, opts%quantum, mw, mem)
        thermo%f4 = thermo%f4 + fe4
        thermo%s4 = thermo%s4 + s4
        thermo%u4 = thermo%u4 + (fe4 + opts%temperature * s4)
        thermo%cv4 = thermo%cv4 + cv4
    end if
    call tmr%tock('four-phonon')

end block latdyn4ph

summary: block
    real(r8) :: f0, f1, f2
    character(len=1000) :: opf1, opf2

    call tmr%stop()
    if (mw%talk) write(*, *) ''
    call tmr%dump(mw, 'Timings:')

end block summary

! And we are done!
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)
end program
