#include "precompilerdefinitions"
program self_energy
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
use selfenergy, only: lo_selfenergy

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
!> The selfenergy
type(lo_selfenergy) :: ls
! timers
real(r8) :: timer_init, timer_se, tt0

! Set up all harmonic properties. That involves reading all the input file,
! creating grids, getting the harmonic properties on those grids.
initharmonic: block
    integer :: i, j, q1, b1
    ! Start MPI and timers
    tt0 = walltime()
    timer_init = tt0
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

heavywork: block
    if (mw%talk) then
        write (*, *) ''
        write (*, *) 'Calculating self-energy'
    end if
    timer_se = walltime()
    call ls%initialize(qp, dr, opts%nbasis, opts%thirdorder, opts%fourthorder, mw, mem)
    call ls%compute(qp, dr, uc, fct, fcf, opts%temperature, opts%isotopescattering, &
                    opts%thirdorder, opts%fourthorder, mw, mem)
    timer_se = walltime() - timer_se
    if (mw%talk) write(*, "(1X,A,F12.3,A)") '... done in ', timer_se, ' s'
end block heavywork

finalize_and_write: block
    real(r8) :: t0
    ! Write thermal conductivity to file

    ! sum up the total time
    if (mw%talk) then
        tt0 = walltime() - tt0
        write(*, *) ''
        write(*, *) '... dumping auxiliary data to files'
        call ls%write_to_hdf5('outfile.selfenergy.hdf5')

    ! Print timings
        write (*, *) ''
        write (*, '(1X,A21)') 'Suggested citations :'
        write (*, '(1X,A41,A56)') 'Software : ', 'F. Knoop et al., J. Open Source Softw 9(94), 6150 (2024)'
        write (*, '(1X,A41,A53)') 'Method : ', 'D. A. Broido et al., Appl Phys Lett 91, 231922 (2007)'
        write (*, '(1X,A41,A43)') 'Iterative Boltzmann transport equation : ', 'M. Omini et al., Phys Rev B 53, 9064 (1996)'
        write (*, '(1X,A41,A49)') 'Algorithm : ', 'A. H. Romero et al., Phys Rev B 91, 214310 (2015)'
        write (*, '(1X,A41,A43)') 'Off diagonal coherent contribution : ', 'L. Isaeva et al., Nat Commun 10 3853 (2019)'

        t0 = timer_init + timer_se
        write (*, *) ' '
        write (*, *) 'Timings:'
        write (*, "(A,F12.3,A,F7.3,A)") '            initialization:', timer_init, ' s, ', real(timer_init*100/tt0), '%'
        write (*, "(A,F12.3,A,F7.3,A)") '    self-energy computation:', timer_se, ' s, ', real(timer_se*100/tt0), '%'
        write (*, "(A,F12.3,A)") '                     total:', tt0, ' seconds'
    end if
end block finalize_and_write

! And we are done!
call ls%destroy()
call mpi_barrier(mw%comm, mw%error)
call mpi_finalize(lo_status)

end program