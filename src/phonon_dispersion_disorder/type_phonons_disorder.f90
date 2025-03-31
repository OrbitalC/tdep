module type_phonons_disorder
!!
!! Handles phonons in disordered materials
!!
use konstanter, only: r8, lo_freqtol, lo_huge, lo_exitcode_param, lo_frequency_hartree_to_icm, &
                      lo_frequency_Hartree_to_THz, lo_frequency_Hartree_to_meV, lo_groupvel_Hartreebohr_to_ms, &
                      lo_kb_Hartree
use gottochblandat, only: lo_harmonic_oscillator_free_energy, lo_classical_harmonic_oscillator_free_energy, &
                          lo_harmonic_oscillator_cv, lo_harmonic_oscillator_entropy
use hdf5_wrappers, only: lo_hdf5_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder

use type_forceconstant_disorder_secondorder_dynamicalmatrix, only: dynamicalmatrix, frequencies_eigenvectors_groupvelocities

implicit none
private
public :: lo_phonons_disorder

!> Phonons in a disordered materials
type lo_phonons_disorder
    !> The number of modes
    integer :: nmodes
    !> The number of energy in the frequency axis
    integer :: nf
    !> The frequencies
    real(r8), dimension(:), allocatable :: omega
    !> The eigenvectors
    real(r8), dimension(:, :), allocatable :: egv
    !> The group velocities
    real(r8), dimension(:, :, :), allocatable :: vel
    !> The linewidths
    real(r8), dimension(:), allocatable :: linewidths
    !> The spectral functions
    real(r8), dimension(:, :), allocatable :: sf
    !> The frequency axis
    real(r8), dimension(:), allocatable :: e_axis
    !> Largest frequency
    real(r8) :: omega_max = -lo_huge
    !> Smallest non zero frequency
    real(r8) :: omega_min = lo_huge
    !> The default smearing for adaptive broadening
    real(r8) :: default_smearing
contains
    !> Generate the phonons
    procedure :: generate
    !> Reshape eigenvectors to ease access to atoms
    procedure :: reshape_eigenvectors
    !> Phonon free energy
    procedure :: phonon_free_energy
    !> Phonon heat capacity
    procedure :: phonon_cv
    !> Phonon entropy
    procedure :: phonon_entropy
    !> Write to HDF5
    procedure :: write_to_hdf5
end type

contains
subroutine generate(ph, fc, uc, mem, mw, dovel)
    !> The phonons
    class(lo_phonons_disorder), intent(out) :: ph
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Do we compute eigenvectors ?
    logical, optional, intent(in) :: dovel

    !> The dynamical matrix
    real(r8), dimension(:, :), allocatable :: dynmat
    !> The dynamical matrix gradient
    real(r8), dimension(:, :, :), allocatable :: dynmat_grad
    !> Do we compute the group velocities ?
    logical :: withvel
    !> On which rank do we solve
    integer :: solrnk

    solrnk = mw%n - 1
    ph%nmodes = fc%na * 3

    ! Do we compute the group velocities ?
    if (present(dovel)) then
        withvel = dovel
    else
        withvel = .false.
    end if

    ! We allocate things
    allocate(ph%omega(ph%nmodes))
    allocate(ph%egv(ph%nmodes, ph%nmodes))
    if (withvel) allocate(ph%vel(3, ph%nmodes, ph%nmodes))

    ! And we can solve
    if (mw%r .eq. solrnk) then
        if (withvel) then
            call mem%allocate(dynmat, [ph%nmodes, ph%nmodes], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(dynmat_grad, [3, ph%nmodes, ph%nmodes], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

!           call fc%dynamicalmatrix(uc, dynmat, dynmat_grad)
!           call fc%frequencies_eigenvectors_groupvelocities(dynmat, dynmat_grad, ph%omega, ph%egv, ph%vel, mem)
            call dynamicalmatrix(fc, uc, dynmat, dynmat_grad)
            call frequencies_eigenvectors_groupvelocities(fc, dynmat, dynmat_grad, ph%omega, ph%egv, ph%vel, mem)

            call mem%deallocate(dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(dynmat_grad, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        else
            call mem%allocate(dynmat, [ph%nmodes, ph%nmodes], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

!           call fc%dynamicalmatrix(uc, dynmat)
!           call fc%frequencies_eigenvectors_groupvelocities(dynmat, frequencies=ph%omega, eigenvectors=ph%egv, groupvel=ph%vel, mem=mem)
            call dynamicalmatrix(fc, uc, dynmat, dynmat_grad)
            call frequencies_eigenvectors_groupvelocities(fc, dynmat, dynmat_grad, ph%omega, ph%egv, ph%vel, mem)

            call mem%deallocate(dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
    end if
    ! Now we broadcast on other ranks
    call mw%bcast(ph%omega, solrnk, __FILE__, __LINE__)
    call mw%bcast(ph%egv, solrnk, __FILE__, __LINE__)
    if (withvel) call mw%bcast(ph%vel, solrnk, __FILE__, __LINE__)

    ! We compute the default smearing
    baseline: block
        !> A little buffer
        real(r8) :: f0
        !> For the do loop
        integer :: imode

        ph%omega_min = lo_huge
        ph%omega_max = maxval(ph%omega)

        f0 = -lo_huge
        do imode=1, ph%nmodes-1
            if (abs(ph%omega(imode)) .gt. lo_freqtol) ph%omega_min = min(ph%omega_min, ph%omega(imode))
            if (abs(ph%omega(imode+1) - ph%omega(imode)) .gt. lo_freqtol) then
                f0 = max(f0, abs(ph%omega(imode+1) - ph%omega(imode)))
            end if
        end do
        ph%default_smearing = f0
    end block baseline
end subroutine

!> Reshape eigenvectors to ease access to atoms
subroutine reshape_eigenvectors(ph, egv)
    !> The phonons
    class(lo_phonons_disorder), intent(in) :: ph
    !> The reshaped eigenvectors
    real(r8), dimension(:, :, :), intent(out) :: egv

    !> Some integers for the do loops and indexing
    integer :: iat, a, idx

    egv = 0.0_r8
    do iat=1, ph%nmodes / 3
        do a=1, 3
            idx = 3 * (iat - 1) + a
            egv(a, iat, :) = ph%egv(idx, :)
        end do
    end do
end subroutine

!> Compute the phonon free energy
pure function phonon_free_energy(ph, temperature, classical) result(f)
    !> The phonons
    class(lo_phonons_disorder), intent(in) :: ph
    !> The temperature
    real(r8), intent(in) :: temperature
    !> Classical limit ?
    logical, optional, intent(in) :: classical
    !> The free energy
    real(r8) :: f

    logical :: cl
    integer :: imode

    if (present(classical)) then
        cl = classical
    else
        cl = .false.
    end if

    ! Return a stupid number if we have imaginary modes
    if (ph%omega_min .lt. 0.0_r8) then
        f = 123456789.0_r8
        return
    end if

    f = 0.0_r8
    do imode=1, ph%nmodes
        if (cl) then
            f = f + lo_classical_harmonic_oscillator_free_energy(temperature, ph%omega(imode))
        else
            f = f + lo_harmonic_oscillator_free_energy(temperature, ph%omega(imode))
        end if
    end do
    f = f / real(ph%nmodes, r8)
end function

!> Compute the phonon heat capacity
pure function phonon_cv(ph, temperature, classical) result(cv)
    !> The phonons
    class(lo_phonons_disorder), intent(in) :: ph
    !> The temperature
    real(r8), intent(in) :: temperature
    !> Do we compute it classically ?
    logical, optional, intent(in) :: classical
    !> The heat capacity
    real(r8) :: cv

    logical :: cl
    integer :: imode

    if (present(classical)) then
        cl = classical
    else
        cl = .false.
    end if

    ! Return a stupid number if we have imaginary modes
    if (ph%omega_min .lt. 0.0_r8) then
        cv = 123456789.0_r8
        return
    end if

    if (cl) then
        cv = lo_kb_Hartree
        return
    end if

    cv = 0.0_r8
    do imode=1, ph%nmodes
        cv = cv + lo_harmonic_oscillator_cv(temperature, ph%omega(imode))
    end do
    cv = cv / real(ph%nmodes, r8)
end function

!> Compute the phonon entropy
pure function phonon_entropy(ph, temperature, classical) result(s)
    !> The phonons
    class(lo_phonons_disorder), intent(in) :: ph
    !> The temperature
    real(r8), intent(in) :: temperature
    !> Classical limit ?
    logical, optional, intent(in) :: classical
    !> The entropy
    real(r8) :: s

    logical :: cl
    integer :: imode

    if (present(classical)) then
        cl = classical
    else
        cl = .false.
    end if

    ! Return a stupid number if we have imaginary modes
    if (ph%omega_min .lt. 0.0_r8) then
        s = 123456789.0_r8
        return
    end if

    s = 0.0_r8
    do imode=1, ph%nmodes
        if (cl) then
            s = s + lo_harmonic_oscillator_entropy(temperature, ph%omega(imode))
        else
            s = s + lo_harmonic_oscillator_entropy(temperature, ph%omega(imode))
        end if
    end do
    s = s / real(ph%nmodes, r8)
end function

!> Write the phonons to HDF5
subroutine write_to_hdf5(ph, filename, unit, mem)
    !> The phonons
    class(lo_phonons_disorder), intent(in) :: ph
    !> The filename
    character(len=*), intent(in) :: filename
    !> The unit for the output
    character(len=*), intent(in) :: unit
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> HDF5 helper
    type(lo_hdf5_helper) :: h5
    !> The unit factor
    real(r8) :: unitfactor
    !> the reshaped eigenvectors
    real(r8), dimension(:, :, :), allocatable :: egv

    ! First, we sort out the unit
    select case(trim(unit))
    case('mev')
        unitfactor = lo_frequency_hartree_to_meV
    case('icm')
        unitfactor = lo_frequency_hartree_to_icm
    case('thz')
        unitfactor = lo_frequency_Hartree_to_THz
    case default
        call lo_stop_gracefully(["Unknown unit, try 'thz', 'icm', or 'mev'"], lo_exitcode_param, __FILE__, __LINE__)
    end select

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    call h5%store_attribute(ph%nmodes, h5%file_id, 'nmodes')
    call h5%store_data(ph%omega * unitfactor, h5%file_id, 'frequencies', enhet=trim(unit), dimensions='mode')

    call mem%allocate(egv, [3, ph%nmodes / 3, ph%nmodes], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call ph%reshape_eigenvectors(egv)
    call h5%store_data(egv, h5%file_id, 'eigenvectors', enhet='dimensionless', dimensions='mode,atom,xyz')
    call mem%deallocate(egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (allocated(ph%vel)) then
        call h5%store_data(ph%vel * lo_groupvel_Hartreebohr_to_ms, h5%file_id, 'group_velocities', &
                           enhet=trim(unit), dimensions='mode,mode,xyz')
    end if

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine
end module
