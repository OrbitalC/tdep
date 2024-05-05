#include "precompilerdefinitions"
module thermal_conductivity
use konstanter, only: r8, lo_tol, lo_sqtol, lo_pi, lo_kb_hartree, lo_freqtol, lo_huge, lo_kappa_au_to_SI, &
                      lo_phonongroupveltol, lo_groupvel_Hartreebohr_to_ms
use gottochblandat, only: lo_sqnorm, lo_planck, lo_outerproduct, lo_chop, lo_gauss, lo_harmonic_oscillator_cv, &
                          lo_trapezoid_integration, lo_linspace, walltime, lo_progressbar_init, lo_progressbar
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_vector, lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_dgelss, lo_gemm
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use hdf5_wrappers, only: lo_hdf5_helper, HID_T
!
use selfenergy, only: lo_selfenergy

implicit none

private
public :: compute_thermal_conductivity
public :: lo_thermalconductivity_helper

!> Little container for all the different approximations to kappa
type lo_thermalconductivity_helper
    !> The Green-Kubo kappa diagonal
    real(r8), dimension(3, 3) :: kappa_gk
    !> The Green-Kubo kappa off-diagonal
    real(r8), dimension(3, 3) :: kappa_gk_od
    !> The RTA kappa diagonal
    real(r8), dimension(3, 3) :: kappa_rta
    !> The RTA kappa off-diagonal
    real(r8), dimension(3, 3) :: kappa_rta_od
    !> The RTA with shift kappa diagonal
    real(r8), dimension(3, 3) :: kappa_srta
    !> The RTA with shift kappa off-diagonal
    real(r8), dimension(3, 3) :: kappa_srta_od
    !> The linewidths in the Green-Kubo approximation
    real(r8), dimension(:, :), allocatable :: lw_gk
    !> The linewidths within perturbation theory
    real(r8), dimension(:, :), allocatable :: lw_pert
    contains
        procedure :: compute_thermal_conductivity
        procedure :: write_to_hdf5 => write_tc_to_hdf5
        procedure :: destroy => destroy_tc
end type

contains

!> Compute kappa including memory effects
subroutine compute_thermal_conductivity(tc, qp, dr, ls, uc, fc, nenergy, temperature, mw, mem)
    !> The thermal conductivity helper
    class(lo_thermalconductivity_helper), intent(out) :: tc
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The selfenergy
    class(lo_selfenergy), intent(in) :: ls
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The second order forceconstants
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> The number of point on the frequency axis
    integer, intent(in) :: nenergy
    !> The temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> The generalized group velocities squared
    real(r8), dimension(:, :, :, :), allocatable :: groupvelsq
    !> Some buffers
    real(r8) :: f0, f1, f2, f3, f4, f5, f6, om1, om2, be, tol, pref, n1, n2, d1, d2, som1, som2, sn1, sn2
    !> Integers for the do loops
    integer :: i, q1, b1, b2
    !> Let's keep track of time
    real(r8) :: t0

    call mem%allocate(groupvelsq, [3, 3, dr%n_mode, dr%n_mode],  persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    allocate(tc%lw_gk(dr%n_mode, qp%n_irr_point))
    allocate(tc%lw_pert(dr%n_mode, qp%n_irr_point))

    pref = lo_pi / (lo_kb_hartree * temperature**2)

    ! For the nice progressbar
    t0 = walltime()
    if (mw%talk) call lo_progressbar_init()

    ! Initialize all the values to 0.0
    tc%lw_gk = 0.0_r8
    tc%lw_pert = 0.0_r8
    tc%kappa_rta = 0.0_r8
    tc%kappa_rta_od = 0.0_r8
    tc%kappa_srta = 0.0_r8
    tc%kappa_srta_od = 0.0_r8
    tc%kappa_gk = 0.0_r8
    tc%kappa_gk_od = 0.0_r8

    do q1=1, qp%n_irr_point
        if (mod(q1, mw%n) .ne. mw%r) cycle

        !> Then we get the generalized eigenvectors
        call get_generalized_eigenvectors(q1, qp, dr, uc, fc, mem, groupvelsq)

        !> Now we can integrate everything
        do b1=1, dr%n_mode
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle

            ! We compute the perturbative lifetime for the first mode
            f1 = ls%evaluate_imag_selfenergy_onepoint(q1, b1, om1)
            d1 = ls%evaluate_real_selfenergy_onepoint(q1, b1, om1)
            n1 = lo_planck(temperature, om1)
            som1 = om1 + d1
            sn1 = lo_planck(temperature, som1)

            do b2=1, dr%n_mode
                om2 = dr%iq(q1)%omega(b2)
                if (om2 .lt. lo_freqtol) cycle

                ! This seems to be a good compromise between accuracy/speed
                tol = 1e-13_r8

                ! We compute the perturbative lifetime for the second mode
                f2 = ls%evaluate_imag_selfenergy_onepoint(q1, b2, om2)
                d2 = ls%evaluate_real_selfenergy_onepoint(q1, b2, om2)
                n2 = lo_planck(temperature, om2)
                som2 = om2 + d2
                sn2 = lo_planck(temperature, som2)

                ! We apply the formula with the Markovian approximation
                f3 = n1 * (n2 + 1.0_r8) + n2 * (n1 + 1.0_r8)
                f3 = f3 * (om1 + om2)**2 / 8.0_r8
                f4 = (f1 + f2) / ((om2 - om1)**2 + (f1 + f2)**2)

                ! We apply the formula with the Markovian approximation, but including the shift
                f5 = sn1 * (sn2 + 1.0_r8) + sn2 * (sn1 + 1.0_r8)
                f5 = f5 * (som1 + som2)**2 / 8.0_r8
                f6 = (f1 + f2) / ((som2 - som1)**2 + (f1 + f2)**2)

                ! We compute the "lifetimes"
                f0 = integrate_spectralfunction(q1, b1, b2, om1, om2, temperature, ls, tol)
                if (b1 .eq. b2) then
                    tc%kappa_gk = tc%kappa_gk + groupvelsq(:, :, b1, b2) * f0 * qp%ip(q1)%integration_weight
                    tc%kappa_rta = tc%kappa_rta + groupvelsq(:, :, b1, b2) * f3 * f4 * qp%ip(q1)%integration_weight
                    tc%kappa_srta = tc%kappa_srta + groupvelsq(:, :, b1, b2) * f5 * f6 * qp%ip(q1)%integration_weight

                    ! While we are at it, we can store the linewidths of the mode
                    ! With full memory effects
                    tc%lw_gk(b1, q1) = 0.5_r8 / pref / f0 *  lo_harmonic_oscillator_cv(temperature, om1)
                    ! Within perturbation theory
                    tc%lw_pert(b1, q1) = f1
                else
                    tc%kappa_gk_od = tc%kappa_gk_od + groupvelsq(:, :, b1, b2) * f0 * qp%ip(q1)%integration_weight
                    tc%kappa_rta_od = tc%kappa_rta_od + groupvelsq(:, :, b1, b2) * f3 * f4 * qp%ip(q1)%integration_weight
                    tc%kappa_srta_od = tc%kappa_srta_od + groupvelsq(:, :, b1, b2) * f3 * f4 * qp%ip(q1)%integration_weight
                end if
            end do
        end do
        if (mw%talk) call lo_progressbar(' ... computing thermal conductivity', q1, &
                                         qp%n_irr_point, walltime() - t0)
    end do
    call mw%allreduce('sum', tc%kappa_gk)
    call mw%allreduce('sum', tc%kappa_gk_od)
    call mw%allreduce('sum', tc%kappa_rta)
    call mw%allreduce('sum', tc%kappa_rta_od)
    call mw%allreduce('sum', tc%kappa_srta)
    call mw%allreduce('sum', tc%kappa_srta_od)
    call mw%allreduce('sum', tc%lw_gk)
    call mw%allreduce('sum', tc%lw_pert)

    ! If we want the linewdith to actually be linewidth we have to take the inverse
    tc%kappa_gk = tc%kappa_gk * pref / uc%volume
    f0 = sum(abs(tc%kappa_gk))
    tc%kappa_gk = lo_chop(tc%kappa_gk, f0*1e-6_r8)
    tc%kappa_gk_od = tc%kappa_gk_od * pref / uc%volume
    f0 = sum(abs(tc%kappa_gk_od))
    tc%kappa_gk_od = lo_chop(tc%kappa_gk_od, f0*1e-6_r8)

    tc%kappa_rta = tc%kappa_rta * pref / uc%volume / lo_pi
    f0 = sum(abs(tc%kappa_rta))
    tc%kappa_rta = lo_chop(tc%kappa_rta, f0*1e-6_r8)
    tc%kappa_rta_od = tc%kappa_rta_od * pref / uc%volume / lo_pi
    f0 = sum(abs(tc%kappa_rta_od))
    tc%kappa_rta_od = lo_chop(tc%kappa_rta_od, f0*1e-6_r8)

    tc%kappa_srta = tc%kappa_srta * pref / uc%volume / lo_pi
    f0 = sum(abs(tc%kappa_srta))
    tc%kappa_srta = lo_chop(tc%kappa_srta, f0*1e-6_r8)
    tc%kappa_srta_od = tc%kappa_srta_od * pref / uc%volume / lo_pi
    f0 = sum(abs(tc%kappa_srta_od))
    tc%kappa_srta_od = lo_chop(tc%kappa_srta_od, f0*1e-6_r8)

    call mem%deallocate(groupvelsq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine


!> Compute the generalized eigenvectors
subroutine get_generalized_eigenvectors(q1, qp, dr, uc, fc, mem, groupvelsq)
    !> For which irreducible q-point are we computing it
    integer, intent(in) :: q1
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> The structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The second orderforce constants
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> The generalized eigenvectors
    real(r8), dimension(:, :, :, :), intent(out) :: groupvelsq

    ! The group velocities
    real(r8), dimension(:, :, :), allocatable :: buf_vel
    !> The derivative of the dynamical matrix
    complex(r8), dimension(:, :, :), allocatable :: buf_grad_dynmat
    !> Some buffers for the eigenvectors and the flattened gradient of the dynamical matrix
    complex(r8), dimension(:, :), allocatable :: kronegv, buf_egv, buf_egw, buf_cm0, buf_cm1, buf_cm2
    !> Other buffer
    complex(r8), dimension(3) :: cv0
    !> Still some buffers
    real(r8), dimension(3) :: v0, v1
    !> Some integers for indices and loop
    integer :: a1, a2, ia, ib, ic, ix, iy, iz, k, iop, ii, b1, b2, bb1, bb2

    !> Let's start by putting the group velocity to zero
    groupvelsq = 0.0_r8

    ! Some buffers
    call mem%allocate(buf_grad_dynmat, [dr%n_mode, dr%n_mode, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_vel, [3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_cm0, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_cm1, [dr%n_mode**2, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_cm2, [dr%n_mode**2, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_egv, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_egw, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(kronegv, [dr%n_mode**2, dr%n_mode**2], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    buf_grad_dynmat = 0.0_r8
    buf_cm0 = 0.0_r8
    buf_cm1 = 0.0_r8
    buf_cm2 = 0.0_r8
    buf_egv = 0.0_r8
    buf_egw = 0.0_r8
    kronegv = 0.0_r8

    ! First we get the dynamical matrix and its derivative
    call fc%dynamicalmatrix(uc, qp%ip(q1), buf_cm0, mem, buf_grad_dynmat, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])

    ! Flatten gradient of dynamical matrix, to be able to use gemm from blas/lapack
    do iz=1, 3
        do a1=1, uc%na
        do a2=1, uc%na
        do ix=1, 3
        do iy=1, 3
            ib = (a1 - 1) * 3 + ix
            ia = (a2 - 1) * 3 + iy
            ic = flattenind(a1, a2, ix, iy, dr%n_mode)
            buf_cm1(ic, iz) = buf_grad_dynmat(ia, ib, iz) / (uc%invsqrtmass(a1) * uc%invsqrtmass(a2))
        end do
        end do
        end do
        end do
    end do

    ! Average over all operations of the little group of q1
    kronegv = 0.0_r8
    do k=1, qp%ip(q1)%n_invariant_operation
        iop = qp%ip(q1)%invariant_operation(k)
        ! Rotate eigenvectors
        call lo_eigenvector_transformation_matrix(buf_cm0, uc%rcart, qp%ip(q1)%r, uc%sym%op(abs(iop)))
        if (iop .lt. 0) then
            call lo_gemm(buf_cm0, dr%iq(q1)%egv, buf_egw)
            buf_egw = conjg(buf_egw)
        else
            call lo_gemm(buf_cm0, dr%iq(q1)%egv, buf_egw)
        end if

        do b1=1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) then
                buf_egw(:, b1) = 0.0_r8
                cycle
            end if
            do a1=1, uc%na
            do ix=1, 3
                ib = (a1 - 1) * 3 + ix
                buf_egw(ib, b1) = buf_egw(ib, b1) * uc%invsqrtmass(a1) / sqrt(dr%iq(q1)%omega(b1) * 2.0_r8)
            end do
            end do
        end do

        do b1=1, dr%n_mode
        do b2=1, dr%n_mode
            buf_egv = buf_egw * conjg(buf_egw(b2, b1))
            do bb1=1, dr%n_mode
            do bb2=1, dr%n_mode
                ia = (b1 - 1) * dr%n_mode + bb1
                ib = (b2 - 1) * dr%n_mode + bb2
                kronegv(ia, ib) = kronegv(ia, ib) + buf_egv(bb2, bb1)
            end do
            end do
        end do
        end do
    end do
    kronegv = kronegv / real(qp%ip(q1)%n_invariant_operation, r8)
    call lo_gemm(kronegv, buf_cm1, buf_cm2)

    do b1=1, dr%n_mode
    do b2=1, dr%n_mode
        ii = (b1 - 1) * dr%n_mode + b2
        cv0 = buf_cm2(ii, :)
        ! remove tiny numbers.
        cv0 = lo_chop(cv0, 1E-10/(lo_groupvel_Hartreebohr_to_ms/1000))
        ! I can take the real part since at the end we sum over
        ! both modes and the imaginary components disappear.
        buf_vel(:, b1, b2) = real(cv0, r8)
    end do
    end do

    ! Now we can finally compute the squared generalized group velocities
    do b1=1, dr%n_mode
    do b2=1, dr%n_mode
        v0 = buf_vel(:, b1, b2)
        do k=1, qp%ip(q1)%n_full_point
            iop = qp%ip(q1)%operation_full_point(k)
            v1 = matmul(uc%sym%op(abs(iop))%m, v0)
            groupvelsq(:, :, b1, b2) = groupvelsq(:, :, b1, b2) + lo_outerproduct(v1, v1)
        end do
    end do
    end do
    groupvelsq = groupvelsq / real(qp%ip(q1)%n_full_point, r8)

    call mem%deallocate(buf_grad_dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_cm0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_cm1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_cm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_egw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(kronegv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

contains
    ! Consistent index flattening? Impossibru to get consistent.
    function flattenind(a1, a2, ix, iy, nb) result(i)
        integer, intent(in) :: a1, a2, ix, iy, nb
        integer :: i

        integer :: ia, ib

        ia = (a1 - 1)*3 + ix
        ib = (a2 - 1)*3 + iy
        i = (ib - 1)*nb + ia
    end function
end subroutine


! Function to do the integration of the product of spectral function, using an adaptive scheme
function integrate_spectralfunction(q1, b1, b2, om1, om2, temperature, ls, tol) result(f0)
    !> The qpoint at which we perform the integration
    integer, intent(in) :: q1
    !> The modes for which we perform the integration
    integer, intent(in) :: b1, b2
    !> The harmonic frequency of the modes
    real(r8), intent(in) :: om1, om2
    !> The temperature
    real(r8), intent(in) :: temperature
    !> The self-energy
    type(lo_selfenergy), intent(in) :: ls
    !> The tolerance for the integration
    real(r8), intent(in) :: tol
    !> The result of the integration
    real(r8) :: f0

    !> Nodes for the integration
    real(r8), dimension(:), allocatable :: node, tmpnode, newnode
    !> Values of integrand at the nodes
    real(r8), dimension(:), allocatable :: values, tmpval, newval
    !> Do we need to split the node ?
    logical, dimension(:), allocatable :: isok, tmpisok
    !> Some buffer for the integral
    real(r8) :: Qt, Qs, fc, c
    !> Number of node that are not converged
    integer :: nnotok
    !> Maximum number of iterations
    integer :: maxiter
    !> Integer for the do loops and indices
    integer :: i, n, ctr, nnode, newnnode


    ! Let's start with a hundred nodes, this should avoid pathological cases
    nnode = 100
    allocate(node(nnode))
    allocate(isok(nnode))
    allocate(values(nnode))
    allocate(newnode(nnode))
    allocate(newval(nnode))

    ! We initialize the nodes
    call lo_linspace(0.0_r8, ls%omega_max, node)
    do n=1, nnode
        values(n) = integrand(node(n), q1, b1, b2, ls, temperature, om1, om2)
    end do
    ! Initialize values
    isok = .false.

    ! To perform the integration, we use an adaptive scheme where we compare the integration
    ! given by current nodes to Simpson's quadrature
    ! If the difference between the two is larger than a tolerance, we subdivise the node in two
    ! TODO it might be possible to reduce the number of evaluation of the integrand
    ! by checking if the nodes is good or not, but in a subtle way
    maxiter = 40
    iterloop: do i=1, maxiter
        nnotok = 0
        newval = 0.0_r8
        newnode = 0.0_r8
        ! We start by checking the nodes to see if we need to add new node and values
        do n=1, nnode-1
            if (isok(n)) cycle  ! DO NOT WORK LIKE THIS
            c = node(n) + 0.5_r8 * (node(n+1) - node(n))
            fc = integrand(c, q1, b1, b2, ls, temperature, om1, om2)
            newnode(n) = c
            newval(n) = fc
            ! Estimation of the integral
            Qt = 0.5_r8 * (node(n+1) - node(n)) * (values(n) + values(n+1))
            ! Estimation of error by Simpson's rule
            Qs = (node(n+1) - node(n)) * (values(n) + values(n+1) + 4.0_r8 * fc) / 6.0_r8
            if(abs(Qt - Qs) .lt. 0.5 * tol) then
                isok(n) = .true.
            else
                isok(n) = .false.
                nnotok = nnotok + 1
            end if
        end do

        ! Now we need to update the nodes
        newnnode = nnode + nnotok
        allocate(tmpnode(newnnode))
        allocate(tmpval(newnnode))
        allocate(tmpisok(newnnode))
        ctr = 0
        do n=1, nnode-1
            ctr = ctr + 1
            tmpnode(ctr) = node(n)
            tmpval(ctr) = values(n)
            tmpisok(ctr) = isok(n)
            if (.not. isok(n)) then
                ctr = ctr + 1
                tmpnode(ctr) = newnode(n)
                tmpval(ctr) = newval(n)
                tmpisok(ctr) = .false.
            end if
        end do
        ! And also the last node
        tmpnode(ctr+1) = node(nnode)

        ! Deallocate old nodes
        deallocate(node)
        deallocate(values)
        deallocate(isok)
        deallocate(newnode)
        deallocate(newval)

        ! Reallocate with new values
        allocate(node(newnnode))
        allocate(values(newnnode))
        allocate(isok(newnnode))
        allocate(newnode(newnnode))
        allocate(newval(newnnode))

        ! And update with the new ones
        node = tmpnode
        values = tmpval
        isok = tmpisok
        nnode = newnnode

        ! If every node are converged, we can exit
        if (nnotok .eq. 0) exit iterloop

        ! And finally, we deallocate the temporary nodes
        deallocate(tmpnode)
        deallocate(tmpval)
        deallocate(tmpisok)
    end do iterloop

    do n=1, nnode
        values(n) = integrand(node(n), q1, b1, b2, ls, temperature, om1, om2)
    end do

    ! Now we have all the nodes to perform the integration
    f0 = 0.0_r8
    do n=1, nnode-1
        f0 = f0 + 0.5_r8 * (node(n+1) - node(n)) * (values(n) + values(n+1))
    end do

   !write(*, *) q1, b1, b2, nnode

    ! And final deallocation
    deallocate(node)
    deallocate(newval)
    deallocate(isok)

    contains

    ! The function to evaluate the integrand
    function integrand(a, q1, b1, b2, ls, temperature, om1, om2) result(res)
        !> The value at which we evaluate
        real(r8), intent(in) :: a
        !> The q-point and mode
        integer, intent(in) :: q1, b1, b2
        !> The self energy
        type(lo_selfenergy), intent(in) :: ls
        !> The temperature at which we evaluate the spectral function
        real(r8), intent(in) :: temperature
        !> The harmonic frequency of the phonons
        real(r8), intent(in) :: om1, om2

        !> The integrand at this point
        real(r8) :: res
        !> The Bose-Einstein distribution
        real(r8) :: be
        !> The spectral functions
        real(r8) :: sf1, sf2

        ! The planck constants, to give the heat-capacity like term
        be = lo_planck(temperature, a)
        ! The two spectral function
        sf1 = ls%evaluate_spectralfunction_onepoint(q1, b1, a)
        if (b1 .eq. b2) then
            sf2 = sf1
        else
            sf2 = ls%evaluate_spectralfunction_onepoint(q1, b2, a)
        end if
        res = a**2 * sf1 * sf2 * be * (be + 1.0_r8)
    end function
end function

!> Write thermal conductivity info to a already open hdf5 file
subroutine write_tc_to_hdf5(tc, input_id)
    !> The thermal conductivity
    class(lo_thermalconductivity_helper), intent(in) :: tc
    !> The id for Hdf5
    integer(HID_T), intent(in):: input_id

    type(lo_hdf5_helper) :: h5
    !> Buffer to contains the kappa tensors
    real(r8), dimension(3, 3) :: m0
    !> The units for kappa
    character(len=5) :: units

    units = "W/m/K"

    h5%file_id = input_id

    m0 = tc%kappa_gk * lo_kappa_au_to_SI
    call h5%store_data(m0, h5%file_id, 'kappa_greenkubo_diagonal', enhet=units)
    m0 = tc%kappa_gk_od * lo_kappa_au_to_SI
    call h5%store_data(m0, h5%file_id, 'kappa_greenkubo_coherent', enhet=units)
    m0 = (tc%kappa_gk + tc%kappa_gk_od) * lo_kappa_au_to_SI
    call h5%store_data(m0, h5%file_id, 'kappa_greenkubo', enhet=units)

    m0 = tc%kappa_rta * lo_kappa_au_to_SI
    call h5%store_data(m0, h5%file_id, 'kappa_rta_diagonal', enhet=units)
    m0 = tc%kappa_rta_od * lo_kappa_au_to_SI
    call h5%store_data(m0, h5%file_id, 'kappa_rta_coherent', enhet=units)
    m0 = (tc%kappa_rta + tc%kappa_rta_od) * lo_kappa_au_to_SI
    call h5%store_data(m0, h5%file_id, 'kappa_rta', enhet=units)

    call h5%store_data(tc%lw_pert, h5%file_id, 'linewidths_perturbative', enhet='Hartree')
    call h5%store_data(tc%lw_gk, h5%file_id, 'linewidths_greenkubo', enhet='Hartree')

    ! We let the privilege to close the file to elsewhere
end subroutine

subroutine destroy_tc(tc)
    !> The thermal conductivity
    class(lo_thermalconductivity_helper), intent(inout) :: tc

    if (allocated(tc%lw_gk)) deallocate(tc%lw_gk)
    if (allocated(tc%lw_pert)) deallocate(tc%lw_pert)
end subroutine
end module
