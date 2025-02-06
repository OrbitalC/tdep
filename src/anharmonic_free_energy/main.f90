program anharmonic_free_energy
!!{!src/anharmonic_free_energy/manual.md!}
use konstanter, only: r8, lo_Hartree_to_eV, lo_kb_Hartree, lo_pressure_HartreeBohr_to_GPa, lo_tol, lo_status
use gottochblandat, only: open_file, walltime, lo_linspace, lo_progressbar_init, lo_progressbar, tochar, &
                          lo_does_file_exist, lo_mean, lo_stddev
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
use thirdorder, only: free_energy_thirdorder
use fourthorder, only: free_energy_fourthorder, free_energy_fourthorder_secondorder
use epot, only: lo_energy_differences
use thermo, only: lo_thermodynamics

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
logical :: havehighorder

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
    end if

    ! Read structure
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=.true.)
    if (mw%talk) write (*, *) '... read unitcell'

    ! Read forceconstants
    call fc%readfromfile(uc, 'infile.forceconstant', mem, verbosity=-1)
    if (mw%talk) write (*, *) '... read second order forceconstant'

    if (opts%thirdorder) call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
    if (opts%fourthorder) call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
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
            thermo%stress_kin(i, i) = -lo_kb_Hartree * temperature * uc%na / uc%volume
        end do
    end if
    write(*, *) 'Free energy'
    write(*, *) thermo%f0
    write(*, *) 'Entropy'
    write(*, *) thermo%s0
    write(*, *) 'Internal energy'
    write(*, *) thermo%u0
    write(*, *) ''

    write(*, *) 'Potential stress'
    do i=1, 3
        write(*, '(1X,F24.12,1X,F24.12,1X,F24.12)') thermo%stress_pot(i, :) * lo_pressure_HartreeBohr_to_GPa
    end do
    write(*, *) 'Kinetic stress'
    do i=1, 3
        write(*, '(1X,F24.12,1X,F24.12,1X,F24.12)') thermo%stress_kin(i, :) * lo_pressure_HartreeBohr_to_GPa
    end do
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
    real(r8), dimension(3, 3) :: s2, s3, s4, sp

    real(r8) :: e2, e3, e4, ep, inverse_kbt, f0
    integer :: i, j, a, b, blocksize

    if (mw%talk) then
        write(*, *) ''
        write(*, *) 'Computing energy differences'
    end if

    call ss%readfromfile('infile.ssposcar')
    call ss%classify('supercell', uc)
    call pot%setup(uc, ss, fc, fct, fcf, mw, opts%verbosity+1)

    ! We fetch the simulation data
    if (lo_does_file_exist('infile.sim.hdf5')) then
        call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity + 2, stride=-1)
    else
        call sim%read_from_file(verbosity=opts%verbosity + 2, stride=1, magnetic=.false., dielectric=.false., nrand=-1, mw=mw)
    end if

    if (abs(opts%temperature - sim%temperature_thermostat) .gt. lo_tol) then
        write(*, *)
        write(*, *) 'WARNING: input and simulation temperatures do not match'
        write(*, '(1X,A24,F24.12)') 'Input temperature', opts%temperature
        write(*, '(1X,A24,F24.12)') 'Simulation temperature', sim%temperature_thermostat
        write(*, *) 'Results might be garbage !'
        write(*, *)
    end if

    thermo%temperature = sim%temperature_thermostat
    if (thermo%temperature .gt. 1E-5_r8) then
        inverse_kbt = 1.0_r8/lo_kb_Hartree/thermo%temperature
    else
        inverse_kbt = 0.0_r8
    end if

    if (mw%talk) write(*, *) '... computing energy differences'

    ! Calculate the baseline energy
    allocate (ediff(sim%nt, 5))
    allocate (sdiff(sim%nt, 3, 3, 5))
    ediff = 0.0_r8
    sdiff = 0.0_r8

    ! First, we compute the energy differences
    do i = 1, sim%nt
        if (mod(i, mw%n) .ne. mw%r) cycle
        call pot%energies_and_forces(sim%u(:, :, i), e2, e3, e4, ep)
        ! We store the energy differences
        ediff(i, 1) = sim%stat%potential_energy(i)
        ediff(i, 2) = sim%stat%potential_energy(i) - e2
        ediff(i, 3) = sim%stat%potential_energy(i) - e2 - ep
        ediff(i, 4) = sim%stat%potential_energy(i) - e2 - ep - e3
        ediff(i, 5) = sim%stat%potential_energy(i) - e2 - ep - e3 - e4

        ! But also the stress
        sdiff(i, :, :, 1) = sim%stat%stress(:, :, i)
!       sdiff(i, :, :, 2) = sim%stat%stress(:, :, i) - s2
!       sdiff(i, :, :, 2) = sim%stat%stress(:, :, i) - s2 - sp
!       sdiff(i, :, :, 2) = sim%stat%stress(:, :, i) - s2 - sp - s3
!       sdiff(i, :, :, 2) = sim%stat%stress(:, :, i) - s2 - sp - s3 - s4
    end do
    call mw%allreduce('sum', ediff)
    call mw%allreduce('sum', sdiff)

    ! We can now compute all properties, including their errors using block averaging
    cumulant = 0.0_r8
    cumulant_var = 0.0_r8
    blocksize = floor(real(sim%nt, r8) / real(opts%nblocks))
    allocate(buf(3, opts%nblocks))
    allocate(buf_stress(3, 3, opts%nblocks))
    do i=1, 5
        ! Here we compute all properties for each blocks
        do j=1, opts%nblocks
            f0 = lo_mean(ediff(blocksize*j:min(blocksize*j+1, sim%nt), i))
            buf(1, j) = lo_mean(ediff(blocksize*j:min(blocksize*j+1, sim%nt), i))
            buf(2, j) = lo_mean((ediff(blocksize*j:min(blocksize*j+1, sim%nt), i) - f0)**2)
            buf(3, j) = lo_mean((ediff(blocksize*j:min(blocksize*j+1, sim%nt), i) - f0)**3)

            do a=1, 3
            do b=1, 3
                buf_stress(a, b, j) = lo_mean(sdiff(i, a, b, blocksize*j:min(blocksize*j+1, sim%nt)))
            end do
            end do
        end do
        ! First order cumulant
        cumulant(1, i) = lo_mean(buf(1, :))
        cumulant_var(1, i) = lo_mean((buf(1, :) - cumulant(1, i))**2)
        ! Second order cumulant
        cumulant(2, i) = lo_mean(buf(2, :))
        cumulant_var(2, i) = lo_mean((buf(2, :) - cumulant(2, i))**2)
        ! Third order cumulant
        cumulant(3, i) = lo_mean(buf(3, :))
        cumulant_var(3, i) = lo_mean((buf(3, :) - cumulant(3, i))**2)

        do a=1, 3
        do b=1, 3
            stress_pot(a, b, i) = lo_mean(buf_stress(a, b, :))
            stress_potvar(a, b, i) = lo_mean((buf_stress(a, b, :) - thermo%stress_pot(a, b))**2)
        end do
        end do
    end do
    cumulant = cumulant/real(ss%na, r8)
    cumulant(2, :) = 0.5_r8 * cumulant(2, :) * inverse_kbt
    cumulant_var(2, :) = 0.5_r8 * cumulant_var(2, :) * inverse_kbt
    cumulant_var(3, :) = cumulant_var(3, :) * inverse_kbt**2 / 6.0_r8

    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Temperature (K) (from infile.meta): ', sim%temperature_thermostat
        write (*, *) 'Potential energy:'
        write (*, "(1X,A,E21.14,1X,A,F21.14)") '                  input: ', cumulant(1, 1)*lo_Hartree_to_eV, ' upper bound:', cumulant(2, 1)*lo_Hartree_to_eV
        write (*, "(1X,A,E21.14,1X,A,F21.14)") '                  E-fc2: ', cumulant(1, 2)*lo_Hartree_to_eV, ' upper bound:', cumulant(2, 2)*lo_Hartree_to_eV
        write (*, "(1X,A,E21.14,1X,A,F21.14)") '            E-fc2-polar: ', cumulant(1, 3)*lo_Hartree_to_eV, ' upper bound:', cumulant(2, 3)*lo_Hartree_to_eV
        if (opts%thirdorder) then
            write (*, "(1X,A,E21.14,1X,A,F21.14)") '        E-fc2-polar-fc3: ', cumulant(1, 4)*lo_Hartree_to_eV, ' upper bound:', cumulant(2, 4)*lo_Hartree_to_eV
        end if
        if (opts%fourthorder) then
            write (*, "(1X,A,E21.14,1X,A,F21.14)") '    E-fc2-polar-fc3-fc4: ', cumulant(1, 5)*lo_Hartree_to_eV, ' upper bound:', cumulant(2, 5)*lo_Hartree_to_eV
        end if
    end if
    call tmr%tock('cumulants')
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

        if (mw%talk) write(*, *) fe3 * lo_Hartree_to_eV
        if (mw%talk) write(*, *) thermo%u3 * lo_Hartree_to_eV
        if (mw%talk) write(*, *) s3 / lo_kb_Hartree
        if (mw%talk) write(*, *) cv3 / lo_kb_Hartree
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

        if (mw%talk) write(*, *) fe4 * lo_Hartree_to_eV
        if (mw%talk) write(*, *) thermo%u4 * lo_Hartree_to_eV
        if (mw%talk) write(*, *) s4 / lo_kb_Hartree
        if (mw%talk) write(*, *) cv4 / lo_kb_Hartree

        call free_energy_fourthorder_secondorder(uc, fcf, qp, dr, opts%temperature, fe4, s4, cv4, opts%quantum, mw, mem)
        thermo%f4 = thermo%f4 + fe4
        thermo%s4 = thermo%s4 + s4
        thermo%u4 = thermo%u4 + (fe4 + opts%temperature * s4)
        thermo%cv4 = thermo%cv4 + cv4

        if (mw%talk) write(*, *) fe4 * lo_Hartree_to_eV
        if (mw%talk) write(*, *) (fe4 + opts%temperature * s4) * lo_Hartree_to_eV
        if (mw%talk) write(*, *) s4 / lo_kb_Hartree
        if (mw%talk) write(*, *) cv4 / lo_kb_Hartree
    end if
    call tmr%tock('four-phonon')

end block latdyn4ph

summary: block
    real(r8) :: f0, f1, f2
    character(len=1000) :: opf1, opf2

    if (mw%talk) then
        write(*, *) ''
        write(*, *) 'SUMMARY OF RESULTS'
        write(*, *) ''
        opf1 = '(1X,3(A20))'
        opf2 = '(1X,3(F20.12))'
        f0 = (cumulant(1, 2) + thermo%f0) * lo_Hartree_to_eV
        f1 = (cumulant(1, 2) + thermo%f0 + thermo%f3) * lo_Hartree_to_eV
        f2 = (cumulant(1, 2) + thermo%f0 + thermo%f3 + thermo%f4) * lo_Hartree_to_eV
        write(*, *) 'Free energy [eV/at]'
        write(*, opf1) 'zeroth order', 'first order', 'second order'
        write(*, opf2) f0, f1, f2

        f0 = (cumulant(1, 2) + thermo%u0) * lo_Hartree_to_eV
        f1 = (cumulant(1, 2) + thermo%u0 + thermo%u3) * lo_Hartree_to_eV
        f2 = (cumulant(1, 2) + thermo%u0 + thermo%u3 + thermo%u4) * lo_Hartree_to_eV
        write(*, *) 'Internal energy [eV/at]'
        write(*, opf1) 'zeroth order', 'first order', 'second order'
        write(*, opf2) f0, f1, f2

        f0 = (thermo%s0) / lo_kb_Hartree
        f1 = (thermo%s0 + thermo%s3) / lo_kb_Hartree
        f2 = (thermo%s0 + thermo%s3 + thermo%s4) / lo_kb_Hartree
        write(*, *) 'Entropy [eV/at/K]'
        write(*, opf1) 'zeroth order', 'first order', 'second order'
        write(*, opf2) f0, f1, f2

        f0 = (thermo%cv0) / lo_kb_Hartree
        f1 = (thermo%cv0 + thermo%cv3) / lo_kb_Hartree
        f2 = (thermo%cv0 + thermo%cv3 + thermo%cv4) / lo_kb_Hartree
        write(*, *) 'Heat capacity [kb]'
        write(*, opf1) 'zeroth order', 'first order', 'second order'
        write(*, opf2) f0, f1, f2

        write (*, *) ''
        write (*, '(1X,A)') 'SUGGESTED CITATIONS:'
        write (*, '(1X,A41,A)') 'Software: ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'

    end if
    call tmr%stop()
    if (mw%talk) write(*, *) ''
    call tmr%dump(mw, 'Timings:')

end block summary

! And we are done!
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)
end program
