program phonon_relation_disorder
!!{!src/phonon_relation_disorder/manual.md!}!
use konstanter, only: r8, lo_emu_to_amu, lo_iou, lo_status, lo_exitcode_param
use gottochblandat, only: tochar, walltime, lo_linspace
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
! use type_forceconstant_disorder_thirdorder, only: lo_forceconstant_disorder_thirdorder
use type_structurefactor_dispersions, only: lo_structurefactor
use type_phonons_disorder, only: lo_phonons_disorder
use type_unitcell_disorder_projection, only: lo_unitcell_disorder_projection
use type_projected_phonons, only: lo_projected_phonons
use type_qpointmesh, only: lo_qpoint_mesh
use type_projected_phonons_bandstructure, only: lo_projected_phonons_bandstructure

use options, only: lo_opts

implicit none
type(lo_opts) :: opts
type(lo_crystalstructure) :: ss
type(lo_forceconstant_disorder_secondorder) :: fc2
type(lo_mpi_helper) :: mw
type(lo_mem_helper) :: mem

real(r8), dimension(:), allocatable :: e_axis
real(r8) :: t0

init: block
    integer :: i, j
    logical :: isok

    call mw%init()
    call mem%init()
    t0 = walltime()
    call opts%parse()
    if (mw%talk .eqv. .false.) opts%verbosity = -100
    if (mw%talk) then
        write(lo_iou, *) 'PHONON RELATION DISORDER'
        write(lo_iou, *) 'Recap of the parameters governing the calculation'
        if (opts%projected) then
            write(lo_iou, '(1X,A)')  'Computing dispersion projected on a unitcell reciprocal space'
        else if (opts%structurefactor) then
            write(lo_iou, '(1X,A)') 'Computing structure factor'
        end if
        write (*, '(1X,A40,L3)') 'Isotope scattering                      ', opts%isotopescattering
        write (*, '(1X,A40,I4,I4,I4)') 'q-point grid                            ', opts%qgrid
        write (*, *) ''
    end if

    ! We will always need the supercell
    if (mw%talk) write(lo_iou, *) '... Reading supercell'
    call ss%readfromfile('infile.ssposcar', verbosity=opts%verbosity)

    if (mw%talk) write(lo_iou, *) '... Reading force constants'
    call fc2%read_from_hdf5('infile.forceconstant_secondorder.hdf5', mem)
    call fc2%check_structure(ss, isok)
    if (.not. isok) then
        call lo_stop_gracefully(['Input structure do not match with input IFC'], &
                                lo_exitcode_param, __FILE__, __LINE__)
    end if

    ! Isotope things
    ! Perhaps non-natural isotope distribution
    if (opts%isotopescattering) then
        if (opts%readiso) then
            if (mw%talk) write (*, *) '... reading isotope distribution from file'
            call ss%readisotopefromfile()
            if (mw%talk) then
                do i = 1, ss%na
                    do j = 1, ss%isotope(i)%n
                        write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                            j, ss%isotope(i)%conc(j), ss%isotope(i)%mass(j) * lo_emu_to_amu
                    end do
                    write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                        trim(ss%atomic_symbol(ss%species(i))), ss%isotope(i)%mean_mass * lo_emu_to_amu, &
                        ss%isotope(i)%disorderparameter
                end do
            end if
        elseif (mw%talk .and. opts%verbosity .gt. 0) then
            do i = 1, ss%na
                do j = 1, ss%isotope(i)%n
                    write (*, "('    isotope: ',I2,' concentration: ',F8.5,' mass: ',F12.6)") &
                        j, ss%isotope(i)%conc(j), ss%isotope(i)%mass(j) * lo_emu_to_amu
                end do
                write (*, "('    element: ',A2,' mean mass: ',F12.6,' mass disorder parameter',F12.9)") &
                    trim(ss%atomic_symbol(ss%species(i))), ss%isotope(i)%mean_mass * lo_emu_to_amu, &
                    ss%isotope(i)%disorderparameter
            end do
        end if
    end if
end block init

! The first part is if we work with a crystal space
if (opts%projected) then
ucell_proj: block
    !> The commensurate phonons
    type(lo_projected_phonons) :: ph_com
    !> The unitcell projection
    type(lo_unitcell_disorder_projection) :: proj
    !> The unitcell
    type(lo_crystalstructure) :: uc

    if (mw%talk) write(lo_iou, *) '... Reading unitcell'
    call uc%readfromfile('infile.ucposcar', verbosity=opts%verbosity)
    call uc%classify('wedge', timereversal=.true.)

    if (mw%talk) write(lo_iou, *) '... Reading mapping between supercell and unitcell'
    call proj%read_from_file(uc, ss, 'infile.unitcell_projection', mem)
    call proj%commensurate_qpoints()
    if (mw%talk) write(lo_iou, *) 'Number of commensurate points (including equivalent edge):', proj%nqpt

    ! We start by projecting phonons on commensurate modes
    projection: block
        !> The realspace phonons
        type(lo_phonons_disorder) :: ph_rs
        !> The third order force constants, if we need the lineshapes
!       type(lo_forceconstant_disorder_thirdorder) :: fc3
        !> The maximum frequency considered here
        real(r8) :: fmax

        if (mw%talk) then
            write(lo_iou, *) ''
            write(lo_iou, *) '... Computing real space phonons'
        end if
        t0 = walltime()
        call ph_rs%generate(fc2, ss, mem, mw)
        if (mw%talk) write(lo_iou, *) 'done in '//tochar(walltime() - t0)//' s'

        ! Stop if we have imaginary frequencies
        if (minval(ph_rs%omega) .lt. 0.0_r8) then
            call lo_stop_gracefully(['Imaginary modes found, it does not make sense to compute things here'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if

        if (opts%thirdorder) then
            fmax = 2.2_r8 * ph_rs%omega_max
        else
            fmax = 1.2_r8 * ph_rs%omega_max
        end if

        ! Maybe we need the spectral function
        if (opts%sf_realspace) then
            ! First we allocate the frequency axis and the spectral function
            ph_rs%nf = opts%nf
            allocate(ph_rs%e_axis(ph_rs%nf))
            allocate(ph_rs%sf(ph_rs%nmodes, ph_rs%nf))
            call lo_linspace(0.0_r8, fmax, ph_rs%e_axis)

            if (opts%isotopescattering) then
!               if (mw%talk) write(lo_iou, *) '... Computing isotope scattering'
!               call ph_rs%compute_isotope_scattering(ss, mw, mem)
            end if
            if (opts%thirdorder) then
!               if (mw%talk) write(lo_iou, *) '... Reading third order force constants'
!               call fc3%read_from_hdf5('infile.forceconstant_thirdorder.hdf5')
!               if (mw%talk) write(lo_iou, *) '... Computing three phonons scattering'
!               call ph_rs%compute_threephonon_scattering(fc3, mw, mem)
            end if
        end if

        if (mw%talk) then
            write(lo_iou, *) ''
            write(lo_iou, *) '... Projecting commensurate phonons in reciprocal space'
        end if
        call ph_com%project_phonons(ph_rs, fc2, proj, ss, fmax, opts%nf, opts%sf_realspace, mem, mw)
    end block projection

    ! Then we interpolate a band structure
    bandstructure: block
        !> The band structure
        type(lo_projected_phonons_bandstructure) :: bs

        if (mw%talk) then
            write(lo_iou, *) ''
            write(lo_iou, *) '... Interpolating phonons on a path in crystal space'
        end if
        call bs%generate(fc2, uc, ss, proj, ph_com, mw, mem, opts%nq, opts%readpathfromfile)
        if (mw%talk) then
            write(lo_iou, *) '... Dumping to file'
            call bs%write_to_hdf5('outfile.phonons_bandstructure_projected.hdf5', opts%enhet, mem)
        end if
    end block bandstructure

    ! Then we interpolate on a grid
    grid_projected: block
!       !> The qpoint mesh
!       class(lo_qpoint_mesh) :: qp
!       !> The dispersion on the grid
!       type(lo_phonon_disorder_grid) :: dr

!       if (mw%talk) write(lo_iou, *) '... Generating q-point grid'
!       call lo_generate_qmesh(qp, uc, opts%qgrid, 'fft', timereversal=.true., &
!                           headrankonly=.false., mw=mw, mem=mem, verbosity=opts%verbosity)
!       if (mw%talk) write(lo_iou, *) '... Interpolating phonons on a grid in crystal space'
!       call dr%generate()
!       if (mw%talk) write(lo_iou, *) '... Dumping to file'
!       call dr%write_to_hdf5('outfile.phonons_grid_projected.hdf5')
    end block grid_projected

    ! And it's done
    io: block
        if (mw%talk) then
            write(lo_iou, *) ''
            write(lo_iou, '(1X,A,A)') 'Band structure can be found in ', 'outfile.phonons_bandstructure_projected.hdf5'
            write(lo_iou, '(1X,A,A)') 'Grid can be found in           ', 'outfile.phonons_grid_projected.hdf5'
            write(lo_iou, *) ''
        end if
    end block io
end block ucell_proj

! Otherwise we work with the structure factor
else if (opts%structurefactor) then
    sf_realspace: block
        !> The realspace phonons
        type(lo_phonons_disorder) :: ph_rs
        !> The third order force constants, if we need the lineshapes
!       type(lo_forceconstant_disorder_thirdorder) :: fc3

        if (mw%talk) write(lo_iou, *) '... Computing real space phonons'
        call ph_rs%generate(fc2, ss, mem, mw)

        if (opts%thirdorder) then
            if (mw%talk) write(lo_iou, *) '... Reading third order force constants'
!           call fc3%read_from_hdf5('infile.forceconstant_thirdorder.hdf5')
        end if

!       if (mw%talk) write(lo_iou, *) '... Projecting commensurate phonons'
!       call ph_com%project(ph_rs, fc, proj, ss, opts%nf, mem, mw)
    end block sf_realspace

    bandstructurefactor: block
    end block bandstructurefactor

    grid_structurefactor: block
        !> The qpoint mesh
        type(lo_qpoint_mesh) :: qp
    end block grid_structurefactor
end if


if (mw%talk) then
    write(lo_iou, '(1X,A)') 'SUGGESTED CITATIONS:'
    write(lo_iou, '(1X,A41,A)') 'Software: ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
end if

! And we are done!
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)
end program
