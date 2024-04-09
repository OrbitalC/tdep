#include "precompilerdefinitions"
module kappa
use konstanter, only: r8, lo_tol, lo_sqtol, lo_pi, lo_kb_hartree, lo_freqtol, lo_huge, lo_kappa_au_to_SI, &
                      lo_phonongroupveltol, lo_groupvel_Hartreebohr_to_ms
use gottochblandat, only: lo_sqnorm, lo_planck, lo_outerproduct, lo_chop
use mpi_wrappers, only: lo_mpi_helper, MPI_SUM, MPI_DOUBLE_PRECISION, MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use type_symmetryoperation, only: lo_operate_on_vector, lo_eigenvector_transformation_matrix
use type_blas_lapack_wrappers, only: lo_dgelss, lo_gemm
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder

use new_scattering, only: lo_scattering_rates

implicit none

private
public :: get_kappa
public :: get_kappa_offdiag
public :: iterative_bte
contains


!> Calculate the thermal conductivity
subroutine get_kappa(dr, qp, uc, temperature, kappa)
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> temperature
    real(r8), intent(in) :: temperature
    !> thermal conductivity tensor
    real(r8), dimension(3, 3), intent(out) :: kappa

    real(r8), dimension(3) :: v0, v1
    real(r8) :: n, f0, omega, omthres, prefactor, velnorm
    integer :: i, j, k, l
    !integer :: iop,ifull

    omthres = dr%omega_min*0.5_r8
    prefactor = 1.0_r8/(uc%volume*lo_kb_hartree*temperature)
    do i = 1, qp%n_full_point
        dr%aq(i)%kappa = 0.0_r8
        k = qp%ap(i)%operation_from_irreducible
        do j = 1, dr%n_mode
            ! Which operation takes this point from the wedge to here
            l = qp%ap(i)%irreducible_index
            ! Skip gamma for acoustic branches
            if (dr%aq(i)%omega(j) .lt. omthres) cycle
            ! Rotate things to this points. Negative is the time reversal thingy, but does not really matter here.
            if (k .gt. 0) then
                v0 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%Fn(:, j), reciprocal=.true.)
                v1 = lo_operate_on_vector(uc%sym%op(k), dr%iq(l)%vel(:, j), reciprocal=.true.)
            else
                v0 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%Fn(:, j), reciprocal=.true.)
                v1 = -lo_operate_on_vector(uc%sym%op(abs(k)), dr%iq(l)%vel(:, j), reciprocal=.true.)
            end if
            ! Get kappa for this q-point
            omega = dr%iq(l)%omega(j)
            n = lo_planck(temperature, omega)
            f0 = omega*(n + 1)*n
            dr%aq(i)%kappa(:, :, j) = prefactor*f0*lo_outerproduct(v1, v0)
        end do
    end do

    ! Sum it ip!
    kappa = 0.0_r8
    do i = 1, qp%n_full_point
        do j = 1, dr%n_mode
            kappa = kappa + dr%aq(i)%kappa(:, :, j) / qp%n_full_point
        end do
    end do
    f0 = sum(abs(kappa))
    kappa = lo_chop(kappa, f0*1E-6_r8)
end subroutine

subroutine get_kappa_offdiag(dr, qp, uc, fc, temperature, mem, mw, kappa_offdiag)
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> second order force constant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> temperature
    real(r8), intent(in) :: temperature
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> thermal conductivity tensor
    real(r8), dimension(3, 3), intent(out) :: kappa_offdiag

    real(r8), dimension(:, :, :), allocatable :: buf_vel
    real(r8), dimension(:, :, :, :), allocatable :: buf_velsq
    real(r8) :: pref, f0, tau, om1, om2, n1, n2, tau1, tau2
    integer :: iq, jmode, kmode

    call mem%allocate(buf_vel, [3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(buf_velsq, [3, 3, dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    kappa_offdiag = 0.0_r8

    do iq = 1, qp%n_irr_point
        if (mod(iq, mw%n) .ne. mw%r) cycle
        buf_vel = 0.0_r8
        buf_velsq = 0.0_r8

        pref = 0.5_r8/(uc%volume*lo_kb_hartree*temperature**2)*qp%ip(iq)%integration_weight

        ! Calculate the off-diagonal group velocity.
        groupvel: block
            complex(r8), dimension(:, :, :), allocatable :: buf_grad_dynmat
            complex(r8), dimension(:, :), allocatable :: kronegv, buf_egv, buf_egw, buf_cm0, buf_cm1, buf_cm2
            complex(r8), dimension(3) :: cv0
            real(r8), dimension(3) :: v0, v1
            integer :: a1, a2, ia, ib, ic, ix, iy, iz, k, iop, i, ii, j, jj

            ! Some buffers
            call mem%allocate(buf_grad_dynmat, [dr%n_mode, dr%n_mode, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
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

            ! Dynamical matrix and derivatives
            call fc%dynamicalmatrix(uc, qp%ip(iq), buf_cm0, mem, buf_grad_dynmat, qdirection=[1.0_r8, 0.0_r8, 0.0_r8])

            ! Flatten gradient of dynamical matrix
            do iz = 1, 3
                do a1 = 1, uc%na !ise%n_atom
                do a2 = 1, uc%na !ise%n_atom
                do ix = 1, 3
                do iy = 1, 3
                    ib = (a1 - 1)*3 + ix
                    ia = (a2 - 1)*3 + iy
                    ic = flattenind(a1, a2, ix, iy, dr%n_mode)
                    buf_cm1(ic, iz) = buf_grad_dynmat(ia, ib, iz)/(uc%invsqrtmass(a1)*uc%invsqrtmass(a2))
                end do
                end do
                end do
                end do
            end do

            ! Average over all operations
            kronegv = 0.0_r8
            do k = 1, qp%ip(iq)%n_invariant_operation
                iop = qp%ip(iq)%invariant_operation(k)
                ! Rotate eigenvectors
                call lo_eigenvector_transformation_matrix(buf_cm0, uc%rcart, qp%ip(iq)%r, uc%sym%op(abs(iop)))
                if (iop .lt. 0) then
                    call lo_gemm(buf_cm0, dr%iq(iq)%egv, buf_egw)
                    buf_egw = conjg(buf_egw)
                else
                    call lo_gemm(buf_cm0, dr%iq(iq)%egv, buf_egw)
                end if

                do i = 1, dr%n_mode
                    if (dr%iq(iq)%omega(i) .gt. lo_freqtol) then
                        do a1 = 1, uc%na
                        do ix = 1, 3
                            ib = (a1 - 1)*3 + ix
                            buf_egw(ib, i) = buf_egw(ib, i)*uc%invsqrtmass(a1)/sqrt(dr%iq(iq)%omega(i)*2.0_r8)
                        end do
                        end do
                    else
                        buf_egw(:, i) = 0.0_r8
                    end if
                end do

                do i = 1, dr%n_mode
                do j = 1, dr%n_mode
                    buf_egv = buf_egw*conjg(buf_egw(j, i))
                    do ii = 1, dr%n_mode
                    do jj = 1, dr%n_mode
                        ia = (i - 1)*dr%n_mode + ii
                        ib = (j - 1)*dr%n_mode + jj
                        kronegv(ia, ib) = kronegv(ia, ib) + buf_egv(jj, ii)
                    end do
                    end do
                end do
                end do
            end do
            kronegv = kronegv/real(qp%ip(iq)%n_invariant_operation, r8)
            ! this means sandwich with eigenvectors,  frequencies,
            ! prefactors and masses are already in there.
            call lo_gemm(kronegv, buf_cm1, buf_cm2)

            ! Keep the group velocities?
            do i = 1, dr%n_mode
            do j = 1, dr%n_mode
                ii = (i - 1)*dr%n_mode + j
                jj = (j - 1)*dr%n_mode + i
                cv0 = buf_cm2(ii, :)
                ! remove tiny numbers.
                cv0 = lo_chop(cv0, 1E-10/(lo_groupvel_Hartreebohr_to_ms/1000))
                ! I can take the real part since at the end we sum over
                ! both modes and the imaginary components disappear.
                buf_vel(:, i, j) = real(cv0, r8)
            end do
            end do

            ! Be very careful with degeneracies. Will give very wrong contribution
            ! to thermal transport if not taken care of properly, I'm pretty sure.
            ! Feels like there could be significant double-counting otherwise.
            ! This seems ok for now, but have to keep an eye out.

            ! Rotate out the group velocities to get the symmetry averaged ones
            buf_velsq = 0.0_r8
            do i = 1, dr%n_mode
            do j = 1, dr%n_mode
                v0 = buf_vel(:, i, j)
                do k = 1, qp%ip(iq)%n_full_point
                    iop = qp%ip(iq)%operation_full_point(k)
                    v1 = matmul(uc%sym%op(abs(iop))%m, v0)
                    buf_velsq(:, :, i, j) = buf_velsq(:, :, i, j) + lo_outerproduct(v1, v1)
                end do
            end do
            end do
            buf_velsq = buf_velsq/real(qp%ip(iq)%n_full_point, r8)

            ! cleanup
            call mem%deallocate(buf_grad_dynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_cm2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_egw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(kronegv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end block groupvel

        do jmode = 1, dr%n_mode
            ! Skip gamma for acoustic branches
            if (dr%iq(iq)%omega(jmode) .lt. lo_freqtol) cycle
            om1 = dr%iq(iq)%omega(jmode)
            do kmode = jmode, dr%n_mode
                ! We only compute the off diagonal contribution
                if (jmode .eq. kmode) cycle
                ! Skip gamma for acoustic branches
                if (dr%iq(iq)%omega(kmode) .lt. lo_freqtol) cycle

                om2 = dr%iq(iq)%omega(kmode)
                n1 = lo_planck(temperature, om1)
                n2 = lo_planck(temperature, om2)
                tau1 = dr%iq(iq)%linewidth(jmode)
                tau2 = dr%iq(iq)%linewidth(kmode)

                f0 = n1*(n2 + 1) + n2*(n1 + 1)
                f0 = f0*(om1 + om2)**2/4.0_r8
                tau = (tau1 + tau2)/((tau1 + tau2)**2 + (om1 - om2)**2)
                kappa_offdiag(:, :) = kappa_offdiag(:, :) + buf_velsq(:, :, jmode, kmode)*tau*f0*pref

                f0 = n1*(n2 + 1) + n2*(n1 + 1)
                f0 = f0*(om1 + om2)**2/4.0_r8
                tau = (tau1 + tau2)/((tau1 + tau2)**2 + (om2 - om1)**2)
                kappa_offdiag(:, :) = kappa_offdiag(:, :) + buf_velsq(:, :, kmode, jmode)*tau*f0*pref
            end do ! k mode
        end do ! j mode
    end do ! i qpt
    call mem%deallocate(buf_vel, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(buf_velsq, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mw%allreduce('sum', kappa_offdiag)

    f0 = sum(abs(kappa_offdiag))
    kappa_offdiag = lo_chop(kappa_offdiag, f0*1E-6_r8)

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


subroutine iterative_bte(sr, dr, qp, uc, temperature, niter, tol, &
                         isotope, threephonon, fourphonon, mw, mem)
    !> integration weights
    type(lo_scattering_rates), intent(inout) :: sr
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Temperature
    real(r8), intent(in) :: temperature
    !> Max number of iterations
    integer, intent(in) :: niter
    !> Tolerance
    real(r8), intent(in) :: tol
    !> What do we compute ?
    logical, intent(in) :: isotope, threephonon, fourphonon
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(:, :, :), allocatable :: Fnb, Fbb
    real(r8), dimension(3, 3) :: kappa
    real(r8), dimension(niter) :: scfcheck
    real(r8) :: mixingparameter
    integer :: iter

    ! set some things and make space
    init: block
        real(r8), dimension(3, 3) :: m0
        mixingparameter = 0.95_r8
        allocate (Fnb(3, dr%n_mode, dr%n_irr_qpoint))
        allocate (Fbb(3, dr%n_mode, dr%n_full_qpoint))
        Fnb = 0.0_r8
        Fbb = 0.0_r8
        scfcheck = 0.0_r8
        ! Get the first kappa-value
        call get_kappa(dr, qp, uc, temperature, kappa)
        m0 = kappa*lo_kappa_au_to_SI
        if (mw%talk) write (*, "(1X,I4,6(1X,F14.4))") 0, m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3)
    end block init

    scfloop: do iter = 1, niter
        ! get the Fn values all across the BZ
        foldout: block
            real(r8), dimension(3) :: v
            integer :: iq, jq, iop, b1

            Fbb = 0.0_r8
            do iq = 1, qp%n_full_point
                if (mod(iq, mw%n) .ne. mw%r) cycle
                iop = qp%ap(iq)%operation_from_irreducible
                jq = qp%ap(iq)%irreducible_index
                do b1 = 1, dr%n_mode
                    if (iop .gt. 0) then
                        v = lo_operate_on_vector(uc%sym%op(iop), dr%iq(jq)%Fn(:, b1), reciprocal=.true.)
                        Fbb(:, b1, iq) = v
                    else
                        v = -lo_operate_on_vector(uc%sym%op(abs(iop)), dr%iq(jq)%Fn(:, b1), reciprocal=.true.)
                        Fbb(:, b1, iq) = v
                    end if
                end do
            end do
            call mw%allreduce('sum', Fbb)
        end block foldout

        ! update F to new values
        updateF: block
            real(r8), dimension(3) :: Fp, Fpp, Fppp, v0
            real(r8) :: iQs, W
            integer :: il, j, b1, b2, b3, b4, q2, q3, q4
            integer :: q1

            Fnb = 0.0_r8
            do il=1, sr%nlocal_point
                q1 = sr%q1(il)
                b1 = sr%b1(il)
                ! prefetch some stuff
                iQS = 1.0_r8 / dr%iq(q1)%qs(b1)
                v0 = 0.0_r8
                if (isotope) then
                    do q2 = 1, qp%n_full_point
                        do b2=1, dr%n_mode
                            Fp = Fbb(:, b2, q2)
                            v0 = v0 + sr%iso(il)%W(b2, q2) * iQs * Fp
                        end do
                    end do
                end if
                if (threephonon) then
                    do j=1, sr%nqpt3ph
                        q2 = sr%threephonon(il)%q2(j)
                        q3 = sr%threephonon(il)%q3(j)
                        do b2=1, dr%n_mode
                        do b3=1, dr%n_mode
                            ! scattering plus
                            Fp = Fbb(:, b2, q2)
                            Fpp = Fbb(:, b3, q3)
                            v0 = v0 + (Fp + Fpp) * sr%threephonon(il)%Wplus(b2, b3, j) * &
                                    iQs * sr%mle_ratio3ph
                            v0 = v0 + (Fp + Fpp) * sr%threephonon(il)%Wminus(b2, b3, j) * &
                                    iQs * sr%mle_ratio3ph * 0.5_r8
                        end do
                        end do
                    end do
                end if
                if (fourphonon) then
                    do j=1, sr%nqpt4ph
                        q2 = sr%fourphonon(il)%q2(j)
                        q3 = sr%fourphonon(il)%q3(j)
                        q4 = sr%fourphonon(il)%q4(j)
                        do b2=1, dr%n_mode
                        do b3=1, dr%n_mode
                        do b4=1, dr%n_mode
                            Fp = Fbb(:, b2, q2)
                            Fpp = Fbb(:, b3, q3)
                            Fppp = Fbb(:, b4, q4)
                            v0 = v0 + (Fp + Fpp + Fppp) * sr%fourphonon(il)%Wpp(b2, b3, b4, j) * &
                                        iQs * sr%mle_ratio4ph * 0.5_r8
                            v0 = v0 + (Fp + Fpp + Fppp) * sr%fourphonon(il)%Wpm(b2, b3, b4, j) * &
                                        iQs * sr%mle_ratio4ph * 0.5_r8
                            v0 = v0 + (Fp + Fpp + Fppp) * sr%fourphonon(il)%Wmm(b2, b3, b4, j) * &
                                        iQs * sr%mle_ratio4ph / 6.0_r8
                        end do
                        end do
                        end do
                    end do
                end if
                Fnb(:, b1, q1) = Fnb(:, b1, q1) - v0
            end do
            call mw%allreduce('sum', Fnb)

            ! make sure degeneracies are satisfied properly and add the previous thing
            do q1=1, qp%n_irr_point
                do b1 = 1, dr%n_mode
                    v0 = 0.0_r8
                    do j = 1, dr%iq(q1)%degeneracy(b1)
                        b2 = dr%iq(q1)%degenmode(j, b1)
                        v0 = v0 + Fnb(:, b2, q1)
                    end do
                    v0 = v0/real(dr%iq(q1)%degeneracy(b1), r8)
                    do j = 1, dr%iq(q1)%degeneracy(b1)
                        b2 = dr%iq(q1)%degenmode(j, b1)
                        Fnb(:, b2, q1) = v0
                    end do
                    ! Add the previous thing
                    Fnb(:, b1, q1) = dr%iq(q1)%F0(:, b1) + Fnb(:, b1, q1)
                end do
            end do
        end block updateF

        ! Add everything together and check convergency
        addandcheck: block
            real(r8), dimension(3, 3) :: m0
            real(r8) :: g0, g1, g2
            integer :: i, j

            g0 = 0.0_r8
            do i = 1, dr%n_irr_qpoint
            do j = 1, dr%n_mode
                g1 = lo_sqnorm(dr%iq(i)%Fn(:, j) - Fnb(:, j, i))
                g2 = lo_sqnorm(dr%iq(i)%Fn(:, j))
                if (g2 .gt. lo_sqtol) then
                    g0 = g0 + g1/g2
                end if
                dr%iq(i)%Fn(:, j) = dr%iq(i)%Fn(:, j)*(1.0_r8 - mixingparameter) + &
                                    mixingparameter*(Fnb(:, j, i))
            end do
            end do
            scfcheck(iter) = g0/qp%n_irr_point/dr%n_mode

            ! Check for convergence. The criterion is that the relative difference between the new
            ! and old Fn is to be one part in 1E-5, for two consecutive iterations.
            if (iter .ge. 3) then
                g0 = sum(scfcheck(iter - 2:iter))
                if (g0 .lt. tol) then
                    exit scfloop
                end if
            end if

            ! We are not converged if we made it here. Get the current kappa, to print to stdout.
            call get_kappa(dr, qp, uc, temperature, kappa)
            m0 = kappa*lo_kappa_au_to_SI
            if (mw%r .eq. 0) write (*, "(1X,I4,6(1X,F14.4),2X,ES10.3)") &
                iter, m0(1, 1), m0(2, 2), m0(3, 3), m0(1, 2), m0(1, 3), m0(2, 3), scfcheck(iter)

            ! If we had too many iterations I want to adjust the mixing a little
            if (iter .gt. 15 .and. mixingparameter .gt. 0.10) then
                mixingparameter = mixingparameter*0.98_r8
            end if
        end block addandcheck
    end do scfloop
    deallocate(Fnb)
    deallocate(Fbb)
end subroutine
end module
