module epot
!! Deal with many kinds of potential energy differences
use konstanter, only: r8, lo_kb_hartree
use gottochblandat, only: tochar, walltime, lo_mean, lo_stddev, &
                          lo_trueNtimes, lo_progressbar_init, lo_progressbar
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_mdsim, only: lo_mdsim

use konstanter, only: lo_Hartree_to_eV

use lo_thermodynamic_helpers, only: lo_thermodynamics, lo_full_to_voigt_33, lo_voigt_to_full_33
implicit none

private
public :: lo_energy_differences

type lo_energy_differences
    !> forceconstants needed for evaluation
    type(lo_forceconstant_secondorder) :: fc2
    type(lo_forceconstant_thirdorder) :: fc3
    type(lo_forceconstant_fourthorder) :: fc4
    !> Holder for long range part of the force constants
    real(r8), dimension(:, :, :, :), allocatable :: fcp
    !> Is it a polar material ?
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

    !> Fourth order IFC
    real(r8), dimension(3, 3, 3, 3) :: m4
    !> Third order IFC
    real(r8), dimension(3, 3, 3) :: m3
    !> Second order IFC
    real(r8), dimension(3, 3) :: m2
    !> Forces and displacements
    real(r8), dimension(3) :: v0, u1, u2, u3, u4
    !> Some integers
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

!> Compute the thermodynamic corrections from the data
subroutine compute_realspace_thermo(pot, ss, sim, thermo, nblocks, mw, mem)
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
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> The stress difference between true and ifcs
    real(r8), dimension(:, :, :), allocatable :: sdiff
    !> The energy difference between true and ifcs
    real(r8), dimension(:, :), allocatable :: ediff

    ! Allocate stuffs
    call mem%allocate(ediff, [sim%nt, 4], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(sdiff, [sim%nt, 6, 2], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    ediff = 0.0_r8
    sdiff = 0.0_r8

    ! Let's compute energies and stress differences
    energy_computation: block
        !> The stress contributions
        real(r8), dimension(3, 3) :: s3
        !> The individual energy contributions
        real(r8) :: e2, e3, e4, ep
        !> The time
        real(r8) :: t0
        !> Some integers
        integer :: i, j, ntot, a, b

        ! For a nice progressbar, we count the number of configuration on the rank
        ntot = 0
        do i=1, sim%nt
            if (mod(i, mw%n) .eq. mw%r) ntot = ntot + 1
        end do

        t0 = walltime()
        if (mw%talk) call lo_progressbar_init()
        j = 0
        do i=1, sim%nt
            if (mod(i, mw%n) .ne. mw%r) cycle
            j = j + 1
            call pot%energies_and_forces(ss, sim%u(:, :, i), e2, e3, e4, ep, s3)
            ! We store the energy differences
            ediff(i, 1) = sim%stat%potential_energy(i)
            ediff(i, 2) = sim%stat%potential_energy(i) - e2 - ep
            ediff(i, 3) = sim%stat%potential_energy(i) - e2 - ep - e3
            ediff(i, 4) = sim%stat%potential_energy(i) - e2 - ep - e3 - e4

            ! But also the stress differences, in Voigt notation
            sdiff(i, :, 1) = lo_full_to_voigt_33(sim%stat%stress(:, :, i))
            sdiff(i, :, 2) = lo_full_to_voigt_33(sim%stat%stress(:, :, i) - s3)
            if (mw%talk .and. lo_trueNtimes(i, 127, ntot)) then
                call lo_progressbar(' ... computing energy from IFC', j, ntot, walltime() - t0)
            end if
        end do
        if (mw%talk) then
            call lo_progressbar(' ... computing energy from IFC', ntot, ntot, walltime() - t0)
        end if
        ! And reduce everything
        call mw%allreduce('sum', ediff)
        call mw%allreduce('sum', sdiff)
    end block energy_computation

    averaging: block
        !> Buffer for the stress tensor
        real(r8), dimension(6, 2) :: stress_voigt
        !> Buffer for the thermodynamic quantities
        real(r8), dimension(2) :: F, S, Cv, U0
        !> Little buffers for the blocks
        real(r8), dimension(:, :), allocatable :: buf, bufcv, bufstress
        !> Some buffers: energy, inverse kbt blabla
        real(r8) :: inverse_kbt, f0, f1, f2, blocksize, sqnm1
        !> Some integers
        integer :: i, j, a, b, istart, iend

        if (mw%talk) write(*, *) '... computing cumulant corrections'
        !> Let's get the inverse temperature
        if (thermo%temperature .gt. 1E-5_r8) then
            inverse_kbt = 1.0_r8/lo_kb_Hartree/thermo%temperature
        else
            inverse_kbt = 0.0_r8
        end if

        ! We need the length of each block for the block averaging
        blocksize = real(sim%nt, r8) / real(nblocks, r8)
        if (thermo%stochastic) then
            call mem%allocate(buf, [2, nblocks], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufcv, [2, nblocks], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufstress, [6, nblocks], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf = 0.0_r8
            bufcv = 0.0_r8
            bufstress = 0.0_r8
            ! For the stochastic, what we compute depends on the order included
            do j=1, nblocks
                ! Index for the block
                istart = ceiling(blocksize*(j-1) + 1)
                iend = min(ceiling(blocksize*j), sim%nt)
                ! Now we compute the average for each block
                if (thermo%fourthorder) then
                    f0 = lo_mean(ediff(istart:iend, 4))
                    buf(1, j) = f0
                    buf(2, j) = lo_mean((ediff(istart:iend, 3) - f0)**2)
                else if (thermo%thirdorder) then
                    f0 = lo_mean(ediff(istart:iend, 3))
                    buf(1, j) = f0
                    buf(2, j) = lo_mean((ediff(istart:iend, 3) - f0)**2)
                else
                    f0 = lo_mean(ediff(istart:iend, 2))
                    buf(1, j) = f0
                    buf(2, j) = lo_mean((ediff(istart:iend, 2) - f0)**2)
                end if

                ! And the stress
                if (thermo%fourthorder .or. thermo%thirdorder) then
                    do a=1, 6
                        bufstress(a, j) = lo_mean(sdiff(istart:iend, a, 2))
                    end do
                else
                    do a=1, 6
                        bufstress(a, j) = lo_mean(sdiff(istart:iend, a, 1))
                    end do
                end if
            end do
            ! First we get the U0 term
            U0(1) = lo_mean(buf(1, :))
            U0(2) = sqrt(lo_mean((buf(1, :) - U0(1))**2))
            U0 = U0 / sim%na
            ! Then the

            ! And the stress
            do a=1, 6
                stress_voigt(a, 1) = lo_mean(bufstress(a, :))
                stress_voigt(a, 2) = lo_stddev(bufstress(a, :)) * sqnm1
            end do

            thermo%first_order%F = U0
            thermo%first_order%stress(:, :, 1) = thermo%first_order%stress(:, :, 1) +  lo_voigt_to_full_33(stress_voigt(:, 1))
            thermo%first_order%stress(:, :, 2) = thermo%first_order%stress(:, :, 2) + lo_voigt_to_full_33(stress_voigt(:, 2))
        else
            call mem%allocate(buf, [2, nblocks], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufcv, [2, nblocks], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(bufstress, [6, nblocks], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            buf = 0.0_r8
            bufcv = 0.0_r8
            bufstress = 0.0_r8
            do j=1, nblocks
                ! Index for the block
                istart = ceiling(blocksize*(j-1) + 1)
                iend = min(ceiling(blocksize*j), sim%nt)
                ! Now we compute the average for each block
                ! This is to get the U0
                f0 = lo_mean(ediff(istart:iend, 2))
                f1 = lo_mean(ediff(istart:iend, 1))
                buf(1, j) = f0
                ! This is to get the correction to the entropy
                buf(2, j) = lo_mean((ediff(istart:iend, 2) - f0)**2)
                ! This is to get the correction to the heat capacity
                bufcv(1, j) = lo_mean((ediff(istart:iend, 1) - f1) * (ediff(istart:iend, 2) - f0))
                ! And the stress
                do a=1, 6
                    bufstress(a, j) = lo_mean(sdiff(istart:iend, a, 1))
                end do
            end do
            ! Something to normalize the block averaging
            if (nblocks .eq. 1) then
                sqnm1 = 1.0_r8
            else
                sqnm1 = 1.0_r8 / sqrt(real(nblocks, r8) - 1.0_r8)
            end if
            ! First we get the U0 term
            U0(1) = lo_mean(buf(1, :))
            U0(2) = lo_stddev(buf(1, :)) * sqnm1
            U0 = U0 / sim%na

            ! Now we need the entropy correction
            S(1) = lo_mean(buf(2, :))
            S(2) = lo_stddev(buf(2, :)) * sqnm1
            S = 0.5_r8 * lo_kb_Hartree * S * inverse_kbt**2 / sim%na

            ! And the heat capacity correction
            Cv(1) = lo_mean(bufcv(1, :))
            Cv(2) = lo_stddev(bufcv(1, :)) * sqnm1
            Cv = Cv * inverse_kbt**2 * lo_kb_Hartree / sim%na

            ! And the stress
            do a=1, 6
                stress_voigt(a, 1) = lo_mean(bufstress(a, :))
                stress_voigt(a, 2) = lo_stddev(bufstress(a, :)) * sqnm1
            end do

            ! And we store the results
            thermo%first_order%F = U0
            thermo%second_order%F(1) = -thermo%temperature * S(1)
            thermo%second_order%F(2) = thermo%temperature * S(2)
            thermo%first_order%U = U0
            thermo%second_order%S = S
            thermo%second_order%Cv = Cv
            thermo%first_order%stress(:,: , 1) = lo_voigt_to_full_33(stress_voigt(:, 1))
            thermo%first_order%stress(:,: , 2) = lo_voigt_to_full_33(stress_voigt(:, 2))
        end if
    end block averaging
end subroutine

end module
