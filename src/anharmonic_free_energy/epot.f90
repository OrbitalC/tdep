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
end type

contains

! !> calculate energies and forces for a given configuration
subroutine energies_and_forces(pot, u, e2, e3, e4, ep)
    !> container for potential energy differences
    class(lo_energy_differences), intent(in) :: pot
    !> displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> energies
    real(r8), intent(out) :: e2, e3, e4, ep

    real(r8), dimension(3, 3, 3, 3) :: m4
    real(r8), dimension(3, 3, 3) :: m3
    real(r8), dimension(3, 3) :: m2
    real(r8), dimension(3) :: v0, u2, u3, u4
    integer :: a1, a2, a3, a4, i1, i2, i3, i4, i

    e2 = 0.0_r8
    e3 = 0.0_r8
    e4 = 0.0_r8
    ep = 0.0_r8

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
                u2 = u(:, a2)
                u3 = u(:, a3)
                do i1 = 1, 3
                do i2 = 1, 3
                do i3 = 1, 3
                    v0(i1) = v0(i1) - m3(i1, i2, i3)*u2(i2)*u3(i3)
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

    if (verbosity .gt. 0) then
        write (*, *) ''
        write (*, *) 'PREPARING POTENTIAL ENERGY DIFFERENCES'
    end if

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

end module
