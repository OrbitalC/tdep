program anharmonic_free_energy
!!{!src/anharmonic_free_energy/manual.md!}
use konstanter, only: r8, lo_Hartree_to_eV, lo_kb_Hartree, lo_pressure_HartreeBohr_to_GPa, lo_status
use gottochblandat, only: open_file, walltime, lo_trace
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use lo_timetracker, only: lo_timer
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_generate_qmesh
use type_phonon_dispersions, only: lo_phonon_dispersions

use options, only: lo_opts
use thirdorder, only: free_energy_thirdorder, elastic_thirdorder
use fourthorder, only: free_energy_fourthorder
use lo_thermodynamic_helpers, only: lo_thermodynamics, lo_symmetrize_stress

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
real(r8) :: timer_init, timer_total

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
        write(*, '(1X,A40,L3)') 'Quantum limit                           ', opts%quantum
        write(*, '(1X,A40,I4,I4,I4)') 'Q-point grid Harmonic                   ', opts%qgh
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

    if (.not. opts%nothirdorder) then
        call fct%readfromfile(uc, 'infile.forceconstant_thirdorder')
        if (mw%talk) write (*, *) '... read third order forceconstant'
    end if
    if (.not. opts%nofourthorder) then
        call fcf%readfromfile(uc, 'infile.forceconstant_fourthorder')
        if (mw%talk) write (*, *) '... read fourth order forceconstant'
    end if

    call tmr%tock('reading input')

    timer_init = walltime() - timer_init

    thermo%temperature = opts%temperature
    thermo%thirdorder = .not. opts%nothirdorder
    thermo%fourthorder = .not. opts%nofourthorder
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
    call lo_generate_qmesh(qp, uc, opts%qgh, 'fft', timereversal=.true., &
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
        ! call dr%phonon_kinetic_stress(qp, uc, temperature, thermo%harmonic%stress(:, :, 1))
    else
        thermo%harmonic%F(1) = dr%phonon_free_energy_classical(temperature)
        thermo%harmonic%U(1) = 3.0_r8 * lo_kb_Hartree * temperature
        thermo%harmonic%S(1) = (thermo%harmonic%U(1) - thermo%harmonic%F(1)) / temperature
        thermo%harmonic%Cv(1) = 3.0_r8 * lo_kb_Hartree
        thermo%harmonic%stress = 0.0_r8
        ! do i=1, 3
        !     thermo%harmonic%stress(i, i, 1) = lo_kb_Hartree * temperature * uc%na / uc%volume
        ! end do
    end if
    ! Usually not needed here, but always a good idea to clean
    ! call lo_symmetrize_stress(thermo%harmonic%stress(:, :, 1), uc)

    ! If we have third order IFC, might as well compute elastic things
    ! if (opts%thirdorder) then
    !     if (mw%talk) write(*, *) '... computing third order contribution to elastic properties'

    !     call elastic_thirdorder(uc, fc, fct, qp, dr, opts%temperature, thermo%threephonon%stress(:, :, 1), &
    !                             thermo%alpha, opts%quantum, mw, mem)
    !     ! Now we symmetrize
    !     call lo_symmetrize_stress(thermo%threephonon%stress(:, :, 1), uc)
    !     call lo_symmetrize_stress(thermo%alpha, uc)
    ! end if
    call tmr%tock('harmonic properties')
end block latdyn

! Skipping cumulant/real-space corrections - only computing phonon contributions


latdyn3ph: block
    !> The phonon dispersion relation
    type(lo_phonon_dispersions) :: dr
    !> The qpoint mesh
    class(lo_qpoint_mesh), allocatable :: qp
    !> The third order contribution
    real(r8) :: fe3, s3, cv3


    if (.not. opts%nothirdorder) then
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

    if (.not. opts%nofourthorder) then
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
        thermo%fourphonon%F(1) = fe4
        thermo%fourphonon%S(1) = s4
        thermo%fourphonon%U(1) = (fe4 + opts%temperature * s4)
        thermo%fourphonon%Cv(1) = cv4
    end if
    call tmr%tock('four-phonon')
end block latdyn4ph

summary: block
    !> Unit conversion factor
    real(r8) :: f_unit, e_unit, s_unit, c_unit, p_unit
    !> Buffer for the stress tensor
    real(r8), dimension(3, 3) :: sigma
    !> Some buffer to print the results
    real(r8), dimension(4) :: buf, buf3, buf4
    !> The pressure
    real(r8) :: pressure
    !> To have a pretty logfile
    character(len=1000) :: opfc, opff, opfs
    !> Some integers
    integer :: i, u

    f_unit = lo_Hartree_to_eV
    e_unit = lo_Hartree_to_eV
    s_unit = 1.0 / lo_kb_Hartree
    c_unit = 1.0 / lo_kb_Hartree
    ! p_unit = lo_pressure_HartreeBohr_to_GPa

    ! Get the stress tensor (harmonic + threephonon only)
    ! sigma = (thermo%harmonic%stress(:, :, 1) + thermo%threephonon%stress(:, :, 1)) * p_unit
    ! pressure = lo_trace(sigma) / 3.0_r8

    if (mw%talk) then
        u = open_file('out', 'outfile.anharmonic_thermodynamics')
        write(u, '(A2,A12,8X,E20.12)') '# ', 'Temperature:', opts%temperature

        opfc = '(4(1X,A24))'
        opff = '(4(1X,F24.12))'
        opfs = '(3(1X,F24.12))'

        ! Harmonic contribution
        buf(1) = thermo%harmonic%F(1) * f_unit
        buf(2) = thermo%harmonic%U(1) * e_unit
        buf(3) = thermo%harmonic%S(1) * s_unit
        buf(4) = thermo%harmonic%Cv(1) * c_unit
        write(*, *) ''
        write(*, *) 'SUMMARY OF RESULTS'
        write(*, *) ''
        write(*, *) 'Harmonic contribution'
        write(*, opfc) 'Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(*, opff) buf
        write(u, *) '# Harmonic contribution'
        write(u, opfc) '# Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
        write(u, opff) buf

        ! Fourth order (four-phonon) contribution
        if (.not. opts%nofourthorder) then
            buf4(1) = thermo%fourphonon%F(1) * f_unit
            buf4(2) = thermo%fourphonon%U(1) * e_unit
            buf4(3) = thermo%fourphonon%S(1) * s_unit
            buf4(4) = thermo%fourphonon%Cv(1) * c_unit
            write(*, *) ''
            write(*, *) 'First Cumulant Contribution'
            write(*, opfc) 'Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
            write(*, opff) buf4
            write(u, *) '# First Cumulant Contribution'
            write(u, opfc) '# Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
            write(u, opff) buf4
        end if

        ! Third order (three-phonon) contribution
        if (.not. opts%nothirdorder) then
            buf3(1) = thermo%threephonon%F(1) * f_unit
            buf3(2) = thermo%threephonon%U(1) * e_unit
            buf3(3) = thermo%threephonon%S(1) * s_unit
            buf3(4) = thermo%threephonon%Cv(1) * c_unit
            write(*, *) ''
            write(*, *) 'Second Cumulant Contribution'
            write(*, opfc) 'Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
            write(*, opff) buf3
            write(u, *) '# Second Cumulant Contribution'
            write(u, opfc) '# Free energy [eV/at]', 'Internal energy [eV/at]', 'Entropy [kB]', 'Heat capacity [kB]'
            write(u, opff) buf3
        end if

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
