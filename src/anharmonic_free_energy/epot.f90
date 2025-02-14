module epot
!! Deal with many kinds of potential energy differences
use konstanter, only: r8, lo_pi, lo_twopi, lo_tol, lo_sqtol, lo_status, lo_Hartree_to_eV, lo_kb_hartree
use gottochblandat, only: tochar, walltime, lo_chop, lo_trueNtimes, lo_progressbar_init, &
                          lo_progressbar, lo_frobnorm, open_file, lo_flattentensor, lo_sqnorm, lo_outerproduct, lo_mean, &
                          lo_points_on_sphere, lo_mean, lo_stddev
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
!use geometryfunctions, only: lo_inscribed_sphere_in_box
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim
use lo_thermodynamic_helpers, only: lo_thermodynamics, lo_symmetrize_stress
implicit none

private
public :: lo_energy_differences

type lo_energy_differences
    ! forceconstants needed for evaluation
    type(lo_forceconstant_secondorder) :: fc2
    type(lo_forceconstant_thirdorder) :: fc3
    type(lo_forceconstant_fourthorder) :: fc4
    real(r8), dimension(:, :, :, :), allocatable :: fcp
    logical :: polar = .false.

contains
    procedure :: setup => setup_potential_energy_differences
    procedure :: energies_and_forces
    procedure :: compute_realspace_thermo
end type

contains

! !> calculate energies and forces for a given configuration
subroutine energies_and_forces(pot, ss, u, e2, e3, e4, ep, s3)
    !> The unitcell
    type(lo_crystalstructure), intent(in) :: ss
    !> container for potential energy differences
    class(lo_energy_differences), intent(in) :: pot
    !> displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> energies
    real(r8), intent(out) :: e2, e3, e4, ep
    !> stress
    real(r8), dimension(3, 3), intent(out) :: s3

    real(r8), dimension(3, 3, 3, 3) :: m4
    real(r8), dimension(3, 3, 3) :: m3
    real(r8), dimension(3, 3) :: m2
    real(r8), dimension(3) :: v0, u1, u2, u3, u4
    integer :: a1, a2, a3, a4, i1, i2, i3, i4, i5, i

    e2 = 0.0_r8
    e3 = 0.0_r8
    e4 = 0.0_r8
    ep = 0.0_r8
    s3 = 0.0_r8

    ! pair term
    do a1 = 1, size(u, 2)
        v0 = 0.0_r8
        do i1 = 1, pot%fc2%atom(a1)%n
            a2 = pot%fc2%atom(a1)%pair(i1)%i2
            m2 = pot%fc2%atom(a1)%pair(i1)%m
            v0 = v0 - matmul(m2, u(:, a2))
        end do
        e2 = e2 - dot_product(u(:, a1), v0)*0.5_r8
    end do

    ! polar term
    if (pot%polar) then
        ep = 0.0_r8
        do a1 = 1, size(u, 2)
            v0 = 0.0_r8
            do a2 = 1, size(u, 2)
                v0 = v0 - matmul(pot%fcp(:, :, a1, a2), u(:, a2))
            end do
            ep = ep - dot_product(u(:, a1), v0)*0.5_r8
        end do
    end if

    ! triplet term
    if (pot%fc3%na .eq. size(u, 2)) then
        do a1 = 1, size(u, 2)
            v0 = 0.0_r8
            do i = 1, pot%fc3%atom(a1)%n
                m3 = pot%fc3%atom(a1)%triplet(i)%m
                a2 = pot%fc3%atom(a1)%triplet(i)%i2
                a3 = pot%fc3%atom(a1)%triplet(i)%i3
                u1 = u(:, a1)
                u2 = u(:, a2)
                u3 = u(:, a3)
                do i1 = 1, 3
                do i2 = 1, 3
                do i3 = 1, 3
                    v0(i1) = v0(i1) - m3(i1, i2, i3)*u2(i2)*u3(i3)
                    do i4 = 1, 3
                        s3(i3, i4) = s3(i3, i4) + m3(i1, i2, i3)*u1(i1)*u2(i2)*pot%fc3%atom(a1)%triplet(i)%rv3(i4) / ss%volume / 6.0_r8
                    end do
                end do
                end do
                end do
            end do
            v0 = v0*0.5_r8
            e3 = e3 - dot_product(v0, u(:, a1))/3.0_r8
        end do
    end if

    ! quartet term
    if (pot%fc4%na .eq. size(u, 2)) then
        do a1 = 1, size(u, 2)
            v0 = 0.0_r8
            do i = 1, pot%fc4%atom(a1)%n
                m4 = pot%fc4%atom(a1)%quartet(i)%m
                a2 = pot%fc4%atom(a1)%quartet(i)%i2
                a3 = pot%fc4%atom(a1)%quartet(i)%i3
                a4 = pot%fc4%atom(a1)%quartet(i)%i4
                u2 = u(:, a2)
                u3 = u(:, a3)
                u4 = u(:, a4)
                do i1 = 1, 3
                do i2 = 1, 3
                do i3 = 1, 3
                do i4 = 1, 3
                    v0(i1) = v0(i1) - m4(i1, i2, i3, i4)*u2(i2)*u3(i3)*u4(i4)
                end do
                end do
                end do
                end do
            end do
            v0 = v0/6.0_r8
            e4 = e4 - dot_product(v0, u(:, a1))/4.0_r8
        end do
    end if
end subroutine

!> Calculate potential energy differences in several ways
subroutine setup_potential_energy_differences(pot, uc, ss, fc2, fc3, fc4, mw, verbosity)
    !> container for potential energy differences
    class(lo_energy_differences), intent(out) :: pot
    !> unitcell
    type(lo_crystalstructure), intent(inout) :: uc
    !> supercell
    type(lo_crystalstructure), intent(inout) :: ss
    !> second order forceconstant
    type(lo_forceconstant_secondorder), intent(in) :: fc2
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fc3
    !> fourth order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fc4
    !> mpi things
    type(lo_mpi_helper), intent(inout) :: mw
    !> how much to talk
    integer, intent(in) :: verbosity

    real(r8) :: timer, t0, t1

    timer = walltime()
    t0 = timer
    t1 = timer

    call fc2%remap(uc, ss, pot%fc2)
    if (verbosity .gt. 0) then
        t1 = walltime()
        write (*, *) '... remapped second order (', tochar(t1 - t0), ')'
        t0 = t1
    end if

    if (fc3%na .gt. 0) then
        call fc3%remap(uc, ss, pot%fc3)
        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) '... remapped third order (', tochar(t1 - t0), ')'
            t0 = t1
        end if
    end if

    if (fc4%na .gt. 0) then
        call fc4%remap(uc, ss, pot%fc4)
        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) '... remapped fourth order (', tochar(t1 - t0), ')'
            t0 = t1
        end if
    end if

    if (fc2%polar) then
        allocate (pot%fcp(3, 3, ss%na, ss%na))
        pot%fcp = 0.0_r8
        call fc2%supercell_longrange_dynamical_matrix_at_gamma(ss, pot%fcp, 1E-15_r8)
        pot%polar = .true.
    else
        pot%polar = .false.
    end if
    if (verbosity .gt. 0) then
        t1 = walltime()
        write (*, *) '... built polar forceconstant (', tochar(t1 - t0), ')'
        t0 = t1
    end if
end subroutine

subroutine compute_realspace_thermo(pot, ss, sim, thermo, nblocks, mw)
    !> container for potential energy differences
    class(lo_energy_differences), intent(in) :: pot
    !> The supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> The simulation
    type(lo_mdsim), intent(in) :: sim
    !> The thermodynamic helper
    type(lo_thermodynamics), intent(inout) :: thermo
    !> Number of blocks
    integer, intent(in) :: nblocks
    !> MPI
    type(lo_mpi_helper), intent(inout) :: mw

    real(r8), dimension(3, 5) :: cumulant, cumulant_var
    real(r8), dimension(3, 3, 5) :: stress_pot, stress_potvar
    real(r8), dimension(:, :, :, :), allocatable :: sdiff
    real(r8), dimension(:, :, :), allocatable :: buf_stress
    real(r8), dimension(:, :), allocatable :: buf
    real(r8), dimension(:, :), allocatable :: ediff
    real(r8), dimension(3, 3) :: s3, sk2, sk3, sk4, skp

    real(r8) :: e2, e3, e4, ep, inverse_kbt, f0
    integer :: i, j, a, b, blocksize, istart, iend

    if (thermo%temperature .gt. 1E-5_r8) then
        inverse_kbt = 1.0_r8/lo_kb_Hartree/thermo%temperature
    else
        inverse_kbt = 0.0_r8
    end if

    allocate(ediff(sim%nt, 5))
    allocate(sdiff(sim%nt, 3, 3, 2))
    ediff = 0.0_r8
    sdiff = 0.0_r8

    if (mw%talk) write(*, *) '... computing energy differences'

    ! Let's compute energies and stress differences
    do i=1, sim%nt
        if (mod(i, mw%n) .ne. mw%r) cycle
        call pot%energies_and_forces(ss, sim%u(:, :, i), e2, e3, e4, ep, s3)
        ! We store the energy differences
        ediff(i, 1) = sim%stat%potential_energy(i)
        ediff(i, 2) = sim%stat%potential_energy(i) - e2
        ediff(i, 3) = sim%stat%potential_energy(i) - e2 - ep
        ediff(i, 4) = sim%stat%potential_energy(i) - e2 - ep - e3
        ediff(i, 5) = sim%stat%potential_energy(i) - e2 - ep - e3 - e4

        ! But also the stress
        sdiff(i, :, :, 1) = sim%stat%stress(:, :, i)
        sdiff(i, :, :, 2) = sim%stat%stress(:, :, i) - s3
    end do
    call mw%allreduce('sum', ediff)
    call mw%allreduce('sum', sdiff)

    if (mw%talk) write(*, *) '... computing free energy correction'
    allocate(buf(3, nblocks))
    allocate(buf_stress(3, 3, nblocks))

    cumulant = 0.0_r8
    cumulant_var = 0.0_r8
    stress_pot = 0.0_r8
    stress_potvar = 0.0_r8
    blocksize = floor((real(sim%nt, r8)) / real(nblocks))
    buf = 0.0_r8
    buf_stress = 0.0_r8

    do i=1, 5
        ! First, let's compute properties for each blocks
        do j=1, nblocks
            istart = blocksize*j
            iend = min(blocksize*(j+1), sim%nt)

            f0 = lo_mean(ediff(istart:iend, i))
            buf(1, j) = lo_mean(ediff(istart:iend, i))
            buf(2, j) = lo_mean((ediff(istart:iend, i) - f0)**2)
            buf(3, j) = lo_mean((ediff(istart:iend, i) - f0)**3)

            do a=1, 3
            do b=1, 3
                buf_stress(a, b, j) = lo_mean(sdiff(istart:iend, a, b, i))
            end do
            end do
        end do
        ! First order cumulant
        cumulant(1, i) = lo_mean(buf(1, :))
        cumulant_var(1, i) = sqrt(lo_mean((buf(1, :) - cumulant(1, i))**2))
        ! Second order cumulant
        cumulant(2, i) = lo_mean(buf(2, :))
        cumulant_var(2, i) = sqrt(lo_mean((buf(2, :) - cumulant(2, i))**2))
        ! Third order cumulant
        cumulant(3, i) = lo_mean(buf(3, :))
        cumulant_var(3, i) = sqrt(lo_mean((buf(3, :) - cumulant(3, i))**2))

        ! Same for the stress
        if (i .lt. 3) then
        do a=1, 3
        do b=1, 3
            stress_pot(a, b, i) = lo_mean(buf_stress(a, b, :))
            stress_potvar(a, b, i) = sqrt(lo_mean((buf_stress(a, b, :) - stress_pot(a, b, i))**2))
        end do
        end do
        end if
    end do
    ! We add the prefactor for the cumulants
    cumulant = cumulant/real(ss%na, r8)
    cumulant_var = cumulant_var/real(ss%na, r8)
    cumulant(2, :) = 0.5_r8 * cumulant(2, :) * inverse_kbt
    cumulant_var(2, :) = 0.5_r8 * cumulant_var(2, :) * inverse_kbt
    cumulant_var(3, :) = cumulant_var(3, :) * inverse_kbt**2 / 6.0_r8
    ! And we can store everything
    thermo%clt_ifc2_1(1) = cumulant(1, 3)
    thermo%clt_ifc2_1(2) = cumulant_var(1, 3)
    thermo%clt_ifc2_2(1) = cumulant(2, 3)
    thermo%clt_ifc2_2(2) = cumulant_var(2, 3)
    thermo%clt_ifc3_1(1) = cumulant(1, 4)
    thermo%clt_ifc3_1(2) = cumulant_var(1, 4)
    thermo%clt_ifc3_2(1) = cumulant(2, 4)
    thermo%clt_ifc3_2(2) = cumulant_var(2, 4)
    thermo%clt_ifc4_1(1) = cumulant(1, 5)
    thermo%clt_ifc4_1(2) = cumulant_var(1, 5)
    thermo%clt_ifc4_2(1) = cumulant(2, 5)
    thermo%clt_ifc4_2(2) = cumulant_var(2, 5)
    ! Stress included
    thermo%stress_pot(:, :) = stress_pot(:, :, 1)
    thermo%stress_potvar(:, :) = stress_potvar(:, :, 1)
    thermo%stress_diff(:, :) = stress_pot(:, :, 2)
    thermo%stress_diffvar(:, :) = stress_potvar(:, :, 2)
end subroutine

end module
