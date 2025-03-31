module type_structurefactor_dispersions
!!
!! Handles the structure factor dispersion
!!
use konstanter, only: r8, lo_huge, lo_hugeint, lo_tol, lo_exitcode_param, lo_freqtol, lo_twopi, &
                      lo_groupvel_Hartreebohr_to_ms, lo_Bohr_to_A, lo_frequency_Thz_to_Hartree, &
                      lo_frequency_Hartree_to_THz, lo_frequency_hartree_to_icm, lo_frequency_hartree_to_meV, &
                      lo_Bohr_to_A
use gottochblandat, only: lo_points_on_sphere, qsort, lo_sqnorm, lo_linspace, lo_gauss, lo_chop, &
                          lo_cross, lo_outerproduct
use geometryfunctions, only: lo_rotation_matrix_from_vector_a_to_b
use hdf5_wrappers, only: lo_hdf5_helper
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_qpointmesh, only: lo_qpoint
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
use type_phonons_disorder, only: lo_phonons_disorder

use type_forceconstant_disorder_secondorder_dynamicalmatrix, only: dynamicalmatrix_structurefactor, eigensolve_structurefactor

implicit none
private
public :: lo_structurefactor
public :: lo_structurefactor_qpoint

!> A q-point in the dispersion relation
type lo_structurefactor_qpoint
    !> The q-distance for this point
    real(r8) :: qdist
    !> The longitudinal frequency
    real(r8) :: omega_longitudinal
    !> The longitudinal velocity
    real(r8) :: vel_longitudinal
    !> The transverse frequency
    real(r8) :: omega_transverse
    !> The transverse velocity
    real(r8) :: vel_transverse
    !> The spectral function for this q-point
    real(r8), dimension(:), allocatable :: chi_longitudinal
    !> The spectral function for this q-point
    real(r8), dimension(:), allocatable :: chi_transverse

contains
    !> Calculate all normal phonon things for one q-point
    procedure :: generate => harmonic_things_for_single_qpoint
    !> Project phonons in real space to the structure factor
    procedure :: project => project_phonons_qpoint

end type

!> Structure factor dispersion on a range of q-point
type lo_structurefactor
    !> Number of points
    integer :: nq=-lo_hugeint
    !> The points
    type(lo_structurefactor_qpoint), dimension(:), allocatable :: q
    !> largest frequency
    real(r8) :: omega_max = -lo_huge
    !> lowest non-zero frequency
    real(r8) :: omega_min = lo_huge
    !> Sensible default smearing
    real(r8) :: default_smearing_longitudinal, default_smearing_transverse
    !> Number of frequencies on the axis
    integer :: nf = -lo_hugeint
    !> The frequency axis
    real(r8), dimension(:), allocatable :: xfreq

contains
    !> Create the full dispersion
    procedure :: generate
    !> Project phonons in real space to the full dispersion
    procedure :: project_phonons
    !> Write to HDF5
    procedure :: write_to_hdf5
!   !> Destroy
!   procedure :: destroy
end type

contains

!> Generate the structure factor for a specific q-point
subroutine harmonic_things_for_single_qpoint(ompoint, fc, uc, qdist, navg, mem)
    !> The point in the dispersion
    class(lo_structurefactor_qpoint), intent(inout) :: ompoint
    !> The forceconstants
    type(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The qpoint
    real(r8), intent(in) :: qdist
    !> The number of point on the sphere to average the structure factor
    integer, intent(in) :: navg
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    init: block
        if (uc%na .ne. fc%na) then
            call lo_stop_gracefully(['Different number of atoms in the structure and force constants'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
        ompoint%omega_longitudinal = 0.0_r8
        ompoint%omega_transverse = 0.0_r8
        ompoint%vel_longitudinal = 0.0_r8
        ompoint%vel_transverse = 0.0_r8
        ompoint%qdist = qdist

        ! If we are at the gamma point, we can already finish here
        if (abs(qdist) .lt. lo_tol) return
    end block init

    ! Here we make a circular average for the q-distance as input
    avg: block
        !> We have to create a qpoints
        type(lo_qpoint) :: qpoint
        !> The dynamical matrix
        complex(r8), dimension(3, 3) :: dynmat
        !> The gradient of the dynamical matrix
        complex(r8), dimension(3, 3, 3) :: dynmat_grad
        !> All qpoints
        real(r8), dimension(:, :), allocatable :: allq
        !> The frequencies for all points
        real(r8), dimension(:), allocatable :: omega_longitudinal, omega_transverse
        !> The group velocities for all points
        real(r8), dimension(:), allocatable :: vel_longitudinal, vel_transverse
        !> The projection of modes on the q-point direction
        real(r8), dimension(3) :: proj
        !> Buffer for frequencies
        real(r8), dimension(3) :: omega
        !> Buffer for eigenvectors
        complex(r8), dimension(3, 3) :: egv
        !> Buffer for group velocities
        real(r8), dimension(3, 3, 3) :: vel
        !> To sort frequencies
        integer, dimension(3) :: ind
        !> Some integer for do loops
        integer :: ii, jj

        call mem%allocate(allq, [3, navg], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(omega_longitudinal, navg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(omega_transverse, navg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(vel_longitudinal, navg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(vel_transverse, navg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        allq = 0.0_r8
        omega_longitudinal = 0.0_r8
        omega_transverse = 0.0_r8
        vel_longitudinal = 0.0_r8
        vel_transverse = 0.0_r8

        ! We generate a bunch a point uniformly distributed on the unit sphere
        call lo_points_on_sphere(allq)
        ! We multiply by the distance to get a uniform distribution of point at the required distance
        allq = allq * qdist
        ! Now we can compute the harmonic properties for all those points
        do ii=1, navg
            qpoint%r = allq(:, ii)

            ! Get eigenvalues and so on for this q-point
!           call fc%dynamicalmatrix_structurefactor(qpoint, dynmat, dynmat_grad, uc)
!           call fc%eigensolve_structurefactor(qpoint, dynmat, dynmat_grad, &
!                                              omega, egv, vel)
            call dynamicalmatrix_structurefactor(fc, qpoint, dynmat, dynmat_grad, uc)
            call eigensolve_structurefactor(fc, qpoint, dynmat, dynmat_grad, &
                                            omega, egv, vel)

            ! Now we need to differentiate between longitudinal and transverse
            ! To do this, first we compute the dot product of eigenvectors and the q-point
            do jj=1, 3
                proj(jj) = abs(dot_product(qpoint%r / qdist, egv(:, jj)))
            end do
            ! Then we sort, so that the higher one is the longitudinal one
            call qsort(proj, ind)
            ! And we distribute
            omega_longitudinal(ii) = omega(ind(3))
            omega_transverse(ii) = 0.5_r8 * (omega(ind(1)) + omega(ind(2)))
            vel_longitudinal(ii) = norm2(vel(:, ind(3), ind(3)))
            vel_transverse(ii) = 0.5_r8 * (norm2(vel(:, ind(1), ind(1))) + norm2(vel(:, ind(2), ind(2))))
        end do

        ! And finally, we average over all q-points
        ompoint%omega_longitudinal = sum(omega_longitudinal) / real(navg, r8)
        ompoint%omega_transverse = sum(omega_transverse) / real(navg, r8)
        ompoint%vel_longitudinal = sum(vel_longitudinal) / real(navg, r8)
        ompoint%vel_transverse = sum(vel_transverse) / real(navg, r8)

        call mem%deallocate(allq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(omega_longitudinal, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(omega_transverse, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(vel_longitudinal, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(vel_transverse, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block avg
end subroutine

!> Generate the structure factor on all points
subroutine generate(sf, qdist, fc, uc, navg, mw, mem)
    !> The structure factor
    class(lo_structurefactor), intent(out) :: sf
    !> The q-point distances
    real(r8), dimension(:), intent(in) :: qdist
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Number of point on the sphere to create average for structure factor
    integer, intent(in) :: navg
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> Some integers
    integer :: iq

    ! Store some things
    sf%nq = size(qdist)
    allocate(sf%q(sf%nq))
    sf%omega_min = lo_huge
    sf%omega_max = -lo_huge

    ! Generate the dispersion on all points
    do iq=1, sf%nq
        sf%q(iq)%omega_longitudinal = 0.0_r8
        sf%q(iq)%omega_transverse = 0.0_r8
        sf%q(iq)%vel_longitudinal = 0.0_r8
        sf%q(iq)%vel_transverse = 0.0_r8
        sf%q(iq)%qdist = 0.0_r8

        ! MPI
        if (mod(iq, mw%n) .ne. mw%r) cycle
        call sf%q(iq)%generate(fc, uc, qdist(iq), navg, mem)
    end do

    ! Now we distribute on all ranks
    do iq=1, sf%nq
        call mw%allreduce('sum', sf%q(iq)%omega_longitudinal)
        call mw%allreduce('sum', sf%q(iq)%omega_transverse)
        call mw%allreduce('sum', sf%q(iq)%vel_longitudinal)
        call mw%allreduce('sum', sf%q(iq)%vel_transverse)
        call mw%allreduce('sum', sf%q(iq)%qdist)
        if (qdist(iq) .gt. lo_tol) then
            sf%omega_min = min(sf%omega_min, minval([sf%q(iq)%omega_longitudinal, sf%q(iq)%omega_transverse]))
            sf%omega_max = max(sf%omega_max, maxval([sf%q(iq)%omega_longitudinal, sf%q(iq)%omega_transverse]))
        end if
    end do

    ! And we can compute a baseline smearing for integrations
    baseline: block
        !> All the frequencies in the dispersion
        real(r8), dimension(:, :), allocatable :: allfreq
        !> A little buffer
        real(r8), dimension(2) :: f0
        !> Integer for direction loop
        integer :: ii, jj

        call mem%allocate(allfreq, [2, sf%nq], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Get all frequencies
        do iq=1, sf%nq
            allfreq(1, iq) = sf%q(iq)%omega_longitudinal
            allfreq(2, iq) = sf%q(iq)%omega_transverse
        end do
        ! Sort them in ascending order
        do ii=1, 2
            call qsort(allfreq(ii, :))
        end do
        ! And get the largest separation as a baseline
        f0 = -lo_huge
        do iq=2, sf%nq-1
            do ii=1, 2
                f0(ii) = max(f0(ii), abs(allfreq(ii, iq+1) - allfreq(ii, iq)))
            end do
        end do
        sf%default_smearing_longitudinal = f0(1)
        sf%default_smearing_transverse = f0(2)
    end block baseline
end subroutine

!> Project phonons on a specific point of the structure factor
subroutine project_phonons_qpoint(ompoint, uc, fc, ph, eig, xfreq, navg, mem)
    !> The point in the dispersion
    class(lo_structurefactor_qpoint), intent(inout) :: ompoint
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The phonons
    type(lo_phonons_disorder), intent(in) :: ph
    !> The reshaped eigenvectors
    real(r8), dimension(:, :, :), intent(in) :: eig
    !> The frequency axis
    real(r8), dimension(:), intent(in) :: xfreq
    !> Number of average direction to take
    integer, intent(in) :: navg
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The complex weight from the projection
    complex(r8) :: cl, ct, expiqr
    !> All q-points for the average
    real(r8), dimension(:, :), allocatable :: allq
    !> Some vectors
    real(r8), dimension(3) :: rij, qv, normqv, normqv2, normqv3
    !> Some buffers
    real(r8) :: wl, wt, sigma, qdotr, invf, f0
    !> Some integers
    integer :: iq, imode, iat, jat, ii, jj, kk, j, nat, nmodes, nf

    ! Test
    real(r8), dimension(3) :: eig1, eig2

    nat = uc%na
    nmodes = nat * 3
    nf = size(ompoint%chi_longitudinal)

    sigma = 0.2_r8 * lo_frequency_Thz_to_Hartree
    invf = real(nf, r8) / xfreq(nf)

    ! We generate a bunch a point uniformly distributed on the unit sphere
    call mem%allocate(allq, [3, navg], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call lo_points_on_sphere(allq)
    ! We multiply by the distance to get a uniform distribution of point at the required distance
    allq = allq * ompoint%qdist
    ! And now we can make an average over all q-points
    do iq=1, navg
        ! Get the qpoint
        qv = lo_chop(allq(:, iq) * lo_twopi, 1e-13_r8)
        qv = allq(:, iq) * lo_twopi
        ! And its normalized version
        normqv = qv / norm2(qv)
        normqv2 = [-normqv(2), normqv(1), 0.0_r8]
        normqv2 = normqv2 / norm2(normqv2)
        normqv3 = lo_cross(normqv, normqv2)
        normqv3 = normqv3 / norm2(normqv3)

        ! We compute the weight of each normal modes for this q-point
        do imode=1, nmodes
            if (ph%omega(imode) .lt. lo_freqtol) cycle
            ! To restrict the number of frequencies on which to compute the projection
            jj = max(floor((ph%omega(imode) - 4.0_r8 * sigma)*invf), 1)
            kk = min(ceiling((ph%omega(imode) + 4.0_r8 * sigma)*invf), nf)
            ! First we compute the weights
            cl = 0.0_r8
            ct = 0.0_r8
            do iat=1, nat
                qdotr = dot_product(qv, uc%rcart(:, iat))
                expiqr = cmplx(cos(qdotr), sin(qdotr), r8)
                cl = cl + dot_product(normqv, eig(:, iat, imode)) * expiqr
                ct = ct + dot_product(normqv2, eig(:, iat, imode)) * expiqr
                ct = ct + dot_product(normqv3, eig(:, iat, imode)) * expiqr
            end do

            wl = abs(cl * conjg(cl)) / real(nat, r8)
            wt = abs(ct * conjg(ct)) / real(nat, r8) * 0.5_r8
            do ii=jj, kk
                f0 = lo_gauss(xfreq(ii), ph%omega(imode), sigma)
                ompoint%chi_longitudinal(ii) = ompoint%chi_longitudinal(ii) + f0 * wl
                ompoint%chi_transverse(ii) = ompoint%chi_transverse(ii) + f0 * wt
            end do
        end do
    end do
    ! And we average over all directions
    ompoint%chi_longitudinal = ompoint%chi_longitudinal / real(navg, r8) / real(nmodes, r8)
    ompoint%chi_transverse = ompoint%chi_transverse / real(navg, r8) / real(nmodes, r8)

    !> We get a cleaner spectral function by setting to 0 what should be zero
    ompoint%chi_longitudinal(1) = 0.0_r8
    ompoint%chi_longitudinal(nf) = 0.0_r8
    ompoint%chi_transverse(1) = 0.0_r8
    ompoint%chi_transverse(nf) = 0.0_r8

    call mem%deallocate(allq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

subroutine project_phonons(sf, uc, fc, ph, nfreq, navg, mw, mem)
    !> The structure factor
    class(lo_structurefactor), intent(inout) :: sf
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The phonons
    type(lo_phonons_disorder), intent(in) :: ph
    !> Number of frequencies on the axis
    integer, intent(in) :: nfreq
    !> Number of average to take
    integer, intent(in) :: navg
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The reshaped eigenvectors
    real(r8), dimension(:, :, :), allocatable :: eig
    !> Some integers
    integer :: iat, a, idx, iq, nat, nmodes

    nat = uc%na
    nmodes = nat * 3

    ! First we reshape the eigenvectors
    call mem%allocate(eig, [3, nat, nmodes], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call ph%reshape_eigenvectors(eig)

    ! We need to allocate some things
    sf%nf = nfreq
    allocate(sf%xfreq(sf%nf))
    call lo_linspace(0.0_r8, ph%omega_max * 1.2_r8, sf%xfreq)

    ! And we project for each q-point
    do iq=1, sf%nq
        allocate(sf%q(iq)%chi_longitudinal(sf%nf))
        allocate(sf%q(iq)%chi_transverse(sf%nf))
        sf%q(iq)%chi_longitudinal = 0.0_r8
        sf%q(iq)%chi_transverse = 0.0_r8

        if (sf%q(iq)%omega_longitudinal .lt. lo_freqtol) cycle
        if (mod(iq, mw%n) .ne. mw%r) cycle

        call sf%q(iq)%project(uc, fc, ph, eig, sf%xfreq, navg, mem)
    end do
    ! Now we can get everything on all ranks
    do iq=1, sf%nq
        call mw%allreduce('sum', sf%q(iq)%chi_longitudinal)
        call mw%allreduce('sum', sf%q(iq)%chi_transverse)
    end do

    call mem%deallocate(eig, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Write the dispersion to HDF5
subroutine write_to_hdf5(sf, filename, unit, mem)
    !> The structure factor
    class(lo_structurefactor), intent(in) :: sf
    !> The name of the file
    character(len=*), intent(in) :: filename
    !> The unit for the output
    character(len=*), intent(in) :: unit
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The HDF5 helper
    type(lo_hdf5_helper) :: h5
    !> Some big buffer
    real(r8), dimension(:, :), allocatable :: cl, ct
    !> Some buffer
    real(r8), dimension(:), allocatable :: f0, f1, f2, f3, f4
    !> The unit factor
    real(r8) :: unitfactor
    !> Some integer
    integer :: iq

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

    call mem%allocate(f0, sf%nq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(f1, sf%nq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(f2, sf%nq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(f3, sf%nq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(f4, sf%nq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    do iq=1, sf%nq
        f0(iq) = sf%q(iq)%omega_longitudinal * unitfactor
        f1(iq) = sf%q(iq)%omega_transverse * unitfactor
        f2(iq) = sf%q(iq)%vel_longitudinal * lo_groupvel_Hartreebohr_to_ms
        f3(iq) = sf%q(iq)%vel_transverse * lo_groupvel_Hartreebohr_to_ms
        f4(iq) = sf%q(iq)%qdist / lo_Bohr_to_A
    end do
    call h5%store_data(f0, h5%file_id, 'longitudinal_frequencies', enhet=trim(unit))
    call h5%store_data(f1, h5%file_id, 'transverse_frequencies', enhet=trim(unit))
    call h5%store_data(f2, h5%file_id, 'longitudinal_group_velocities', enhet='m/s')
    call h5%store_data(f3, h5%file_id, 'longitudinal_group_velocities', enhet='m/s')
    call h5%store_data(f4, h5%file_id, 'reciprocal_distance', enhet='1/angstrom')

    call mem%deallocate(f0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(f1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(f2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(f3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(f4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (sf%nf .gt. 0) then
        call h5%store_data(sf%xfreq * unitfactor, h5%file_id, 'frequency_axis', enhet=trim(unit))
        call mem%allocate(cl, [sf%nq, sf%nf], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ct, [sf%nq, sf%nf], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        do iq=1, sf%nq
            cl(iq, :) = sf%q(iq)%chi_longitudinal / unitfactor
            ct(iq, :) = sf%q(iq)%chi_transverse / unitfactor
        end do
        call h5%store_data(cl, h5%file_id, 'longitudinal_spectralfunction', &
                           dimensions='frequency,qpoint', enhet='1.0'//trim(unit))
        call h5%store_data(ct, h5%file_id, 'transverse_spectralfunction', &
                           dimensions='frequency,qpoint', enhet='1.0'//trim(unit))
        call mem%deallocate(cl, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ct, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine
end module
