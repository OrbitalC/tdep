program extract_forceconstants_disorder
!!{!src/extract_forceconstants_disorder/manual.md!}!
use konstanter, only: r8, lo_iou, lo_status, lo_Bohr_to_A, lo_frequency_Hartree_to_meV, &
                      lo_groupvel_Hartreebohr_to_ms, lo_twopi, lo_exitcode_param, lo_tol, &
                      lo_Hartree_to_eV, lo_force_HartreeBohr_to_eVA, lo_forceconstant_2nd_HartreeBohr_to_eVA, &
                      lo_pressure_HartreeBohr_to_GPa
use gottochblandat, only: tochar, walltime, lo_does_file_exist, lo_linspace, lo_chop
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use lo_memtracker, only: lo_mem_helper
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
use type_forceconstant_disorder_thirdorder, only: lo_forceconstant_disorder_thirdorder

use options, only: lo_opts
use fit_helper, only: lo_fit_helper

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: uc
type(lo_forceconstant_disorder_secondorder) :: fc2
type(lo_forceconstant_disorder_thirdorder) :: fc3
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem
real(r8) :: t0
! Specific for this program
type(lo_fit_helper) :: fh

! Recap thing, read the cell and set the parameters in the helpers
init: block
    call mw%init()
    call mem%init()
    t0 = walltime()
    call opts%parse()

    if (mw%talk .eqv. .false.) opts%verbosity = -100
    if (mw%talk) then
        write(lo_iou, *) 'EXTRACT FORCECONSTANTS DISORDER'
        write(lo_iou, *) 'Recap of the parameters governing the calculation'
        write(lo_iou, '(1X,A40,F10.6)') 'Cutoff 2nd order (A)                         ', opts%rc2 * lo_Bohr_to_A
        write(lo_iou, '(1X,A40,F10.6)') 'Conjugate gradient prefactor coefficients    ', opts%alpha_secondorder
        write(lo_iou, '(1X,A40,F10.6)') 'Conjugate gradient prefactor positions       ', opts%alpha_pos
        write(lo_iou, '(1X,A40,F10.6)') 'Conjugate gradient threshold (eV)            ', opts%thresh * lo_Hartree_to_eV
        write(lo_iou, '(1X,A40,I10)')   'Conjugate gradient number of steps           ', opts%nsteps
        write(lo_iou, '(1X,A40,L10)')   'Conjugate gradient update positions          ', opts%update_pos
        write(lo_iou, '(1X,A40,L10)')   'Average positions                            ', opts%avg_pos
        write(lo_iou, '(1X,A40,L10)')   'Read previous IFC                            ', opts%readifc
        write(lo_iou, '(1X,A40,L10)')   'Prefit IFC                                   ', opts%prefit
        if (opts%rc3 .gt. 0) then
        write(lo_iou, *) ''
        write(lo_iou, '(1X,A40,F10.6)') 'Cutoff 3rd order (A)                         ', opts%rc3 * lo_Bohr_to_A
        write(lo_iou, '(1X,A40,F10.6)') 'Conjugate gradient prefactor coefficients    ', opts%alpha_thirdorder
!       write(lo_iou, '(1X,A40,F10.6)') 'Conjugate gradient threshold (eV)            ', opts%thresh * lo_Hartree_to_eV
        write(lo_iou, '(1X,A40,I10)')   'Conjugate gradient number of steps           ', opts%nsteps_thirdorder
!       write(lo_iou, '(1X,A40,L10)')   'Read previous IFC                            ', opts%readifc_thirdorder
        end if
        write(lo_iou, *) ''
    end if

    if (mw%talk) write(lo_iou, *) '... reading cell'
    call uc%readfromfile('infile.ssposcar', verbosity=opts%verbosity)
    if (opts%rc2 .gt. uc%maxcutoff()) then
        opts%rc2 = uc%maxcutoff()
        if (mw%talk) then
            write(lo_iou, *) 'Second order cutoff larger than the maximum available cutoff, reducing it to', opts%rc2 * lo_Bohr_to_A
        end if
    end if
    if (opts%rc3 .gt. uc%maxcutoff()) then
        opts%rc3 = uc%maxcutoff()
        if (mw%talk) then
            write(lo_iou, *) 'Third order cutoff larger than the maximum available cutoff, reducing it to', opts%rc3 * lo_Bohr_to_A
        end if
    end if

    ! Stuffs for the second order
    fh%alpha_secondorder = opts%alpha_secondorder
    fh%alpha_pos = opts%alpha_pos
    fh%nsteps = opts%nsteps
    fh%update_pos = opts%update_pos
    fh%verbosity = opts%verbosity

    ! Stuffs for the third order
    fh%alpha_thirdorder = opts%alpha_thirdorder
    fh%nsteps_thirdorder = opts%nsteps_thirdorder
    fh%thresh = opts%thresh
end block init

! Here we get the simulation and distribute everything on the ranks
get_sim: block
    !> The simulation
    type(lo_mdsim) :: sim
    !> To ease the distribution
    integer, dimension(:), allocatable :: mystep
    !> Integers to get indices and count
    integer :: nsteps_per_rank, istart, istop, istep, k

    ! First we read the files
    if (lo_does_file_exist('infile.sim.hdf5')) then
        call sim%read_from_hdf5('infile.sim.hdf5', verbosity=opts%verbosity+2, stride=opts%stride)
    else
        call sim%read_from_file(verbosity=opts%verbosity+2, stride=opts%stride, mw=mw)
    end if
    ! We need to remove any drift
    call sim%remove_force_and_center_of_mass_drift()

    ! Distribute the simulations steps on the ranks
    call mem%allocate(mystep, sim%nt, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    mystep = 0
    do istep=1, sim%nt
        if (mod(istep, mw%n) .eq. mw%r) mystep(istep) = 1
    end do
    fh%mynsteps = sum(mystep)
    fh%ntotsteps = fh%mynsteps
    call mw%allreduce('sum', fh%ntotsteps)
    allocate(fh%f(3, sim%na, fh%mynsteps))
    allocate(fh%u(3, sim%na, fh%mynsteps))
    allocate(fh%epot(fh%mynsteps))
    fh%f = 0.0_r8
    fh%u = 0.0_r8
    fh%epot = 0.0_r8
    k = 0
    do istep=1, sim%nt
        if (mod(istep, mw%n) .eq. mw%r) then
            k = k + 1
            fh%f(:, :, k) = sim%f(:, :, istep)
            fh%u(:, :, k) = sim%u(:, :, istep)
            fh%epot(k) = sim%stat%potential_energy(istep)
        end if
    end do

    ! Do we average the positions ?
    if (opts%avg_pos) then
        if (mw%talk) write(lo_iou, *) '... Averaging positions'
        call fh%get_average_positions(uc, mem, mw)
    end if
    call mem%deallocate(mystep, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end block get_sim

! Now we can initialize the IFC, maybe we read or start from scratch
init_ifc: block
    !> Do the IFC check out with the input structure ?
    logical :: isok

    if (mw%talk) then
        write(lo_iou, *) ''
        write(lo_iou, *) 'INITIALIZING SECOND ORDER IFC'
    end if

    t0 = walltime()
    if (opts%readifc) then
        ! Here we actually read IFC
        if (mw%talk) write(lo_iou, *) '... reading previous IFC'
        call fc2%read_from_hdf5('infile.forceconstant_secondorder.hdf5', mem)

        ! Check that the input cutoff correspond to the cutoff in the IFC
        if (abs(opts%rc2 - fc2%cutoff) .gt. lo_tol) then
            call lo_stop_gracefully(['Cutoff of the infile forceconstants do not match'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
        ! Check that the structure match with the IFC
        call fc2%check_structure(uc, isok)
        if (.not. isok) then
            call lo_stop_gracefully(['Input structure do not match with input IFC'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
    else
        ! Start from scratch
        if (mw%talk) write(lo_iou, *) '... Computing neighbors'
        call fc2%get_pairs(uc, opts%rc2, mem)
    end if
    ! Now we create the pair for per atom access
    call fc2%create_pairs(mem)
    if (mw%talk) write(lo_iou, *) '... Creating constraints'
    call fc2%create_hermitian_constraints()

    ! Print some info
    if (mw%talk) then
        write(lo_iou, '(A40,2X,I7)') 'Number of pairs: ', fc2%npairs
        write(lo_iou, '(A40,2X,I7)') 'Number of coefficients: ', fc2%npairs * 9
        write(lo_iou, '(A40,2X,F10.5)') 'Average number of neighbors: ', 2.0_r8 * real(fc2%npairs, r8) / real(fc2%na, r8)
        write(lo_iou, '(A40,2X,F10.5)') 'Min. atomic distance (Angs): ', fc2%mindist * lo_Bohr_to_A
    end if
    if (mw%talk) write(lo_iou, *) 'done in '//tochar(walltime() - t0)//' s'
end block init_ifc

! Here we actually do the fit, with the conjugate gradient. We might do a prefit if asked for
fit_secondorder: block
    ! If we didn't read the IFC, we might want to initialize the coefficients
    if (.not. opts%readifc .and. opts%prefit) then
        if (mw%talk) then
            write(lo_iou, *) ''
            write(lo_iou, *) 'PREFIT OF SECOND ORDER IFC'
        end if
        if (mw%talk) write(lo_iou, *) '... Initializing coefficients'
        t0 = walltime()
        call fh%prefit_secondorder(fc2, mem, mw)
        if (mw%talk) write(lo_iou, *) '... Applying Hermitian constraint'
        call fc2%apply_hermitian_constraints(fc2%m, mem)
        if (mw%talk) write(lo_iou, *) 'done in '//tochar(walltime() - t0)//' s'
    end if

    if (fh%nsteps .gt. 0) then
        if (mw%talk) write(lo_iou, *) ''
        if (mw%talk) write(lo_iou, *) 'FITTING SECOND ORDER IFC'
        if (mw%talk) write(lo_iou, *) '... Conjugate gradient for second order force constants'
        t0 = walltime()
        if (fh%nsteps .gt. 0) call fh%optimize_secondorder(fc2, uc, mw, mem)
        if (mw%talk) write(lo_iou, *) 'done in '//tochar(walltime() - t0)//' s'
    end if
end block fit_secondorder

! Now we can take care of the fitting of the third order
fit_thirdorder: block

    if (opts%thirdorder) then
        if (mw%talk) then
            write(lo_iou, *) ''
            write(lo_iou, *) 'FITTING THIRD ORDER IFC'
        end if
        call fc3%get_triplets(uc, opts%rc3, mw)
        call fc3%generate(uc%na, opts%p3)
        if (mw%talk) then
            write(lo_iou, '(A40,2X,I7)') 'Number of triplets: ', fc3%ntriplets
!           write(lo_iou, '(A40,2X,I7)') 'Number of coefficients:', fc3%na * fc3%nval
            write(lo_iou, '(A40,2X,I7)') 'Number of coefficients:', fc3%ntriplets * 27
            write(lo_iou, '(A40,2X,F10.5)') 'Average number of triplets: ', 3.0_r8 * real(fc3%ntriplets, r8) / real(fc3%na, r8)
        end if
        if (mw%talk) write(lo_iou, *) '... Removing harmonic forces'
        call fh%prepare_forces_thirdorer(fc2, mem)

        if (mw%talk) write(lo_iou, *) '... Gradient descent for third order force constants'
        t0 = walltime()
        call fh%optimize_thirdorder(fc3, uc, mw, mem)
        if (mw%talk) write(lo_iou, *) 'done in '//tochar(walltime() - t0)//' s'
    end if
end block fit_thirdorder

! Compute the diagnostic for the fit: anharmonicity measure, std and so on
diagnostic: block
    !> Results for the energy
    real(r8), dimension(3, 2) :: res_e, res_f
    !> Some buffers
    real(r8) :: hermit, rot, huang, r2_e, r2_f, sigma_e, sigma_f, u0, u3
    !> Continued
    real(r8) :: r3_e, r3_f, sigma_e3, sigma_f3
    !> Some integer
    integer :: i

    if (mw%talk) then
        write(lo_iou, *) ''
        write(lo_iou, *) 'DIAGNOSTIC'
    end if
    call fh%compute_diagnostic(fc2, fc3, u0, u3, res_e, res_f, opts%thirdorder, mem, mw)
    call fc2%check_hermitian(fc2%m, hermit)
    call fc2%check_rotational(fc2%m, uc, rot)
    call fc2%check_huang(fc2%m, uc, huang)
    call fc2%get_elastic_constants(uc, mw, 0)
    r2_e = 1.0_r8 - (res_e(2, 2) / res_e(1, 1))**2
    r2_f = 1.0_r8 - (res_f(2, 2) / res_f(1, 1))**2
    r3_e = 1.0_r8 - (res_e(3, 2) / res_e(1, 1))**2
    r3_f = 1.0_r8 - (res_f(3, 2) / res_f(1, 1))**2
    sigma_e = res_e(2, 2) / res_e(1, 1)
    sigma_e3 = res_e(3, 2) / res_e(1, 1)
    sigma_f = res_f(2, 2) / res_f(1, 2)
    sigma_f3 = res_f(3, 2) / res_f(1, 2)
    res_e = res_e * lo_Hartree_to_eV * 1000
    res_f = res_f * lo_force_HartreeBohr_to_eVA
    if (mw%talk) then
        write(lo_iou, '(1X,A)') 'ENERGIES (meV/atom)'
        write(lo_iou, '(21X,A,10X,A,5X,A,5X,A)') 'std', 'std(res)', 'R^2(res)', 'normalized std(res)'
        write(lo_iou, '(1X,A15,1X,F12.6,5X,A,12X,A,12X,A,12X,A)') 'input', res_e(1, 1), '-', '-', '-'
        write(lo_iou, '(1X,A15,4(1X,F12.6))') 'second order', res_e(2, 1), res_e(2, 2), r2_e, sigma_e
        if (opts%thirdorder) then
            write(lo_iou, '(1X,A15,4(1X,F12.6))') 'third order', res_e(3, 1), res_e(3, 2), r3_e, sigma_e3
        end if
        write(lo_iou, *) ''
        write(lo_iou, '(1X,A)') 'FORCES (eV/A)'
        write(lo_iou, '(21X,A,10X,A,5X,A,5X,A)') 'std', 'std(res)', 'R^2(res)', 'normalized std(res)'
        write(lo_iou, '(1X,A15,1X,F12.6,5X,A,12X,A,12X,A)') 'input', res_f(1, 1), '-', '-', '-'
        write(lo_iou, '(1X,A15,4(1X,F12.6),13X,A)') 'second order', res_f(2, 1), res_f(2, 2), r2_f, sigma_f, '<-- anharmonicity measure'
        if (opts%thirdorder) then
            write(lo_iou, '(1X,A15,4(1X,F12.6),13X)') 'third order', res_f(3, 1), res_f(3, 2), r3_f, sigma_f3
        end if

        write(lo_iou, *) ''
        write(lo_iou, '(1X,A40,2X,E15.6)') 'Hermitian sum rule (eV/A^2):', hermit * lo_forceconstant_2nd_HartreeBohr_to_eVA
        write(lo_iou, '(1X,A40,2X,E15.6)') 'Rotational sum rule (eV/A):', rot * lo_force_HartreeBohr_to_eVA
        write(lo_iou, '(1X,A40,2X,E15.6)') 'Huang sum rule (eV):', huang * lo_Hartree_to_eV
        write(lo_iou, '(1X,A40,2X,F15.6)') 'Baseline energy <V - V_2> (eV/at):', u0 * lo_Hartree_to_eV
        if (opts%thirdorder) then
            write(lo_iou, '(1X,A40,2X,F15.6)') 'Baseline energy <V - V_2 - V_3> (eV/at):', u3 * lo_Hartree_to_eV
        end if
        write(lo_iou, *)
        write(lo_iou, *) 'ELASTIC CONSTANTS (GPa)'
        do i = 1, 6
            write (*, "(6(3X,F15.5))") lo_chop(fc2%elastic_constants_voigt(:, i)*lo_pressure_HartreeBohr_to_GPa, lo_tol)
        end do
        write(lo_iou, *) ''
    end if
end block diagnostic

io: block
    if (mw%talk) then
        write(lo_iou, *) ''
        write(lo_iou, *) '... Dumping outfiles'
        call fc2%write_to_hdf5('outfile.forceconstant_secondorder.hdf5')
        write(lo_iou, '(1X,A,A)') 'Second order force constants can be found in  ', 'outfile.forceconstant_secondorder.hdf5'
        if (opts%update_pos .or. opts%avg_pos) then
            call uc%writetofile('outfile.ssposcar', 1)
            write(lo_iou, '(1X,A,A)') 'Structure can be found in                     ', 'outfile.ssposcar'
        end if
        write(lo_iou, *) ''
        write (*, '(1X,A)') 'SUGGESTED CITATIONS:'
        write (*, '(1X,A41,A)') 'Software: ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
    end if
end block io

! And we are done !
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)

end program
