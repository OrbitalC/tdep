module fit_helper
!!
!! A little module to help fit the IFC for a random system
!!
use konstanter, only: r8, lo_huge, lo_exitcode_blaslapack, lo_force_HartreeBohr_to_eVA, lo_Hartree_to_eV, &
                      lo_Bohr_to_A, lo_kb_Hartree, lo_iou
use gottochblandat, only: tochar, lo_trueNtimes, lo_outerproduct, lo_nullspace_coefficient_matrix, walltime
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
use type_forceconstant_disorder_thirdorder, only: lo_forceconstant_disorder_thirdorder
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_crystalstructure, only: lo_crystalstructure
use type_mdsim, only: lo_mdsim
use lo_memtracker, only: lo_mem_helper
use type_blas_lapack_wrappers, only: lo_dgelss
use lo_randomnumbers, only: lo_mersennetwister


implicit none
private
public :: lo_fit_helper

type lo_fit_helper
    !> The total number of steps in the simulation
    integer :: ntotsteps
    !> Number of steps in this rank
    integer :: mynsteps
    !> The forces on this rank
    real(r8), dimension(:, :, :), allocatable :: f, f2
    !> The displacements on this rank
    real(r8), dimension(:, :, :), allocatable :: u
    !> The energies on this rank
    real(r8), dimension(:), allocatable :: epot, epot2
    !> The number of steps for the gradient descent
    integer :: nsteps
    !> The number of steps for the gradient descent at the thirdorder
    integer :: nsteps_thirdorder
    !> The gradient prefactors for the descent at second order
    real(r8) :: alpha_secondorder, alpha_pos
    !> The gradient prefactors for the descent at third order
    real(r8) :: alpha_thirdorder
    !> The stopping threshold for the gradient descent
    real(r8) :: thresh
    !> Do we update positions
    logical :: update_pos
    !> Verbosity
    integer :: verbosity
contains
    !> Prefit the forceconstants by doing fit on individual atoms
    procedure :: prefit_secondorder
    !> Perform the gradient descent to fit the force constants
    procedure :: optimize_secondorder
    !> Compute the gradient of the loss for the second order
    procedure :: compute_gradient_secondorder
    !> Remove harmonic forces to fit third order
    procedure :: prepare_forces_thirdorer
    !> Perform the gradient descent to fit the second order force constants
    procedure :: optimize_thirdorder
    !> Compute the gradient of the loss for the third order
    procedure :: compute_gradient_thirdorder
    !> Compute the average positions
    procedure :: get_average_positions
    !> Compute the diagnositc after the fit
    procedure :: compute_diagnostic
end type

contains
! Prefit of the IFC by LSQ-fitting each atom independently and averaging the coefficients
subroutine prefit_secondorder(fh, fc, mem, mw)
    !> The fit helper
    class(lo_fit_helper), intent(in) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The feature matrix and target vectors
    real(r8), dimension(:, :), allocatable :: amat, yvec
    !> Some integers for indices and do loops
    integer :: iat, jat, i, j, a, b, idx, icoef, n, ii, jj, ic

    ! The target vector won't change size between atom, we can allocate it here
    call mem%allocate(yvec, [fh%mynsteps, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    yvec = 0.0_r8

    ! Here we do the thing for every atom
    fc%m = 0.0_r8
    do iat=1, fc%na
        ! The target vector is simply composed of the forces on iat
        yvec = -transpose(fh%f(:, iat, :))

        call mem%allocate(amat, [fh%mynsteps, fc%atom(iat)%n * 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        amat = 0.0_r8
        do j=1, fc%atom(iat)%n
            ic = fc%atom(iat)%ic(j)
            ii = fc%ij(1, abs(ic))
            jj = fc%ij(2, abs(ic))
            do b=1, 3
                idx = 3 * (j - 1) + b
                amat(:, idx) = (fh%u(b, jj, :) - fh%u(b, ii, :)) * sign(1.0_r8, real(ic, r8))
            end do
        end do

        ! And we fit
        call lo_dgelss(amat, yvec, info=i)
        if (i .ne. 0) then
            call lo_stop_gracefully(['dgelss exit status'//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
        end if
        call mem%deallocate(amat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! And we can distribute. The 0.5 is to average because the coefficient will appear for ij and ji
        do j=1, fc%atom(iat)%n
            ic = fc%atom(iat)%ic(j)
            do b=1, 3
                idx = 3 * (j - 1) + b
                fc%m(b, :, abs(ic)) = fc%m(b, :, abs(ic)) + 0.5_r8 * yvec(idx, :)
            end do
        end do
    end do
    ! And we can average everything together
    call mw%allreduce('sum', fc%m)
    fc%m = fc%m / real(mw%n, r8)

    call mem%deallocate(yvec, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Perform the gradient descent to fit the second order IFC
subroutine optimize_secondorder(fh, fc, uc, mw, mem)
    !> The fit helper
    class(lo_fit_helper), intent(inout) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The gradient of the loss
    real(r8), dimension(:, :, :), allocatable :: grad, gm, hm
    !> The harmonic forces
    real(r8), dimension(:, :), allocatable :: feff, fdiff, uavg, gpos, hpos
    !> The measures for the gradient descent
    real(r8) :: maxgrad, maxavgf, loss, l0, rel_loss, gam, gampos
    !> Integers for do loop
    integer :: istep, itest, iat, a
    !> Do we exit the loop ?
    logical :: converged, speak
    !> For better print
    character(len=1000) :: opf1, opf2

    call mem%allocate(hm, [3, 3, fc%npairs], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(gm, [3, 3, fc%npairs], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(grad, [3, 3, fc%npairs], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(feff, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(fdiff, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if(fh%update_pos) then
        call mem%allocate(uavg, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(gpos, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(hpos, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        uavg = 0.0_r8
    end if

    ! We start by applying the Hermitian constraints, to start in the right subspace
    if (mw%talk) write(lo_iou, *) '... Projecting starting point on Hermitian null space'
    call fc%apply_hermitian_constraints(fc%m, mem)

    ! We need to compute the gradient and loss before the start
    if (mw%talk) write(lo_iou, *) '... Computing gradient of least-squares loss'
    call fh%compute_gradient_secondorder(fc, grad, feff, fdiff, loss, mem, mw)
    call fc%apply_hermitian_constraints(grad, mem)

    ! Initialize the vectors for the Conjugate Gradient descent
    gm = -grad
    hm = -grad
    ! Also for the positions, if needed
    if (fh%update_pos) then
        gpos = -fdiff
        hpos = -fdiff
    end if

    ! For better looking log
    opf1 = '(A7,5X,A19,5X,A27,5X,A11,5X,A15)'
    opf2 = '(I7,2X,E15.6,9X,E15.6,17X,E15.6,1X,E15.6)'
    maxgrad = maxval(abs(grad))
    maxavgf = maxval(abs(fdiff)) * lo_force_HartreeBohr_to_eVA
    if (mw%talk) then
        write(lo_iou, opf1) 'Step', 'Max abs. grad. (eV)', 'Max abs. force diff. (eV/A)', 'Loss (eV/a)', 'Rel. Loss diff.'
        write(lo_iou, '(I7,2X,E15.6,9X,E15.6,17X,E15.6)') 0, maxgrad * lo_Hartree_to_eV, maxavgf, sqrt(loss) * lo_force_HartreeBohr_to_eVA
    end if

    ! And let's do the conjugated gradient descent
    converged = .false.
    gradloop: do istep=1, fh%nsteps
        ! We keep track of the previous loss
        l0 = loss

        ! Update the coefficients
        fc%m = fc%m + fh%alpha_secondorder * hm

        ! If needed, we update the positions
        if (fh%update_pos) then
            ! Update the intermediate vector for the conjugate gradient
            gam = sum((fdiff + gpos) * fdiff) / sum(gpos * gpos)
            gpos = -fdiff
            hpos = gpos + gampos * hpos
            do iat=1, fc%na
                do a=1, 3
                    fh%u(a, iat, :) = fh%u(a, iat, :) + hpos(a, iat) * fh%alpha_pos
                    uavg(a, iat) = uavg(a, iat) + hpos(a, iat) * fh%alpha_pos
                end do
            end do
        end if

        ! Compute the gradient and project on the null space of the constraints
        call fh%compute_gradient_secondorder(fc, grad, feff, fdiff, loss, mem, mw)
        call fc%apply_hermitian_constraints(grad, mem)

        ! And update the CG parameters
        gam = sum((grad + gm) * grad) / sum(gm * gm)
        gm = -grad
        hm = gm + gam * hm

        ! Now we can report how the gradient descent is doing
        maxgrad = maxval(abs(grad))
        maxavgf = maxval(abs(fdiff)) * lo_force_HartreeBohr_to_eVA
        rel_loss = sqrt(abs(loss - l0)) / sqrt(l0)

        if (maxgrad .lt. fh%thresh .or. rel_loss .lt. 1e-12_r8) converged = .true.
        speak = .false.
        if (converged .or. istep .eq. fh%nsteps .or. lo_trueNtimes(istep, 10, fh%nsteps) .or. fh%verbosity .eq. 2) speak = .true.

        ! We print 10 steps in total or if we have to exit the loop
        if (speak .and. mw%talk) then
            write(lo_iou, opf2) istep, maxgrad * lo_Hartree_to_eV, maxavgf, sqrt(loss) * lo_force_HartreeBohr_to_eVA, rel_loss
        end if
        ! And maybe we can exit ?
        if (maxgrad .lt. fh%thresh .or. rel_loss .lt. 1e-12_r8) exit gradloop
    end do gradloop

    ! Report on why we exit
    if (mw%talk) then
        if (maxgrad .lt. fh%thresh) then
            write(lo_iou, *) 'Threshold value reached in '//tochar(istep)//' steps'
        else if (rel_loss .lt. 1e-12_r8) then
            write(lo_iou, *) 'Relative loss difference lower than 1e-12'
        else
            write(lo_iou, *) 'Max number of step reached'
        end if
    end if

    ! We might need to update the structure
    if (fh%update_pos) then
        uc%rcart = uc%rcart - uavg
        do iat=1, fc%na
            uc%r(:, iat) = uc%cartesian_to_fractional(uc%rcart(:, iat))
        end do
    end if

    ! And we apply the hermitian constraints one last time, to be sure
    ! it shouldn't change anything anyway
    call fc%apply_hermitian_constraints(fc%m, mem)

    call mem%deallocate(gm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(hm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(grad, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(feff, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(fdiff, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (fh%update_pos) then
        call mem%deallocate(uavg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(gpos, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(hpos, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if
end subroutine

!> Compute the gradient of the least-squares loss
subroutine compute_gradient_secondorder(fh, fc, grad, feff, fdiff, loss, mem, mw)
    !> The fit helper
    class(lo_fit_helper), intent(in) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The gradient
    real(r8), dimension(:, :, :) :: grad
    !> The harmonic forces
    real(r8), dimension(:, :), intent(inout) :: feff
    !> The mean harmonic forces
    real(r8), dimension(:, :), intent(inout) :: fdiff
    !> The least-squares loss
    real(r8), intent(out) :: loss
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> Some buffer
    real(r8), dimension(3) :: dui, duj, dfi, dfj
    !> Integer for do loops
    integer :: istep, iat, j, jat, idx, a, b

    fdiff = 0.0_r8
    loss = 0.0_r8
    grad = 0.0_r8
    do istep=1, fh%mynsteps
        ! We need to compute the harmonic forces @TODO can we compute everything in one go ?
        call fc%forces(fh%u(:, :, istep), feff)
        fdiff = fdiff + (fh%f(:, :, istep) - feff) / real(fh%ntotsteps, r8)
        loss = loss + sum((fh%f(:, :, istep) - feff)**2) / real(fh%ntotsteps, r8) / real(fc%na, r8) / 3.0_r8
        do j=1, fc%npairs
            iat = fc%ij(1, j)
            jat = fc%ij(2, j)

            ! Acoustic sum ruled displacements
            dui = fh%u(:, iat, istep) - fh%u(:, jat, istep)
            duj = fh%u(:, jat, istep) - fh%u(:, iat, istep)
            ! Real minus harmonic forces
            dfi = fh%f(:, iat, istep) - feff(:, iat)
            dfj = fh%f(:, jat, istep) - feff(:, jat)

            ! Gradient is df x u (with x dyadic product)
            ! First atom, it's the normal ordering
            grad(:, :, j) = grad(:, :, j) + lo_outerproduct(duj, dfi) / real(fh%ntotsteps, r8)
            ! Second atom, it's the transpose
            grad(:, :, j) = grad(:, :, j) + lo_outerproduct(dfj, dui) / real(fh%ntotsteps, r8)
        end do
    end do
    call mw%allreduce('sum', grad)
    call mw%allreduce('sum', fdiff)
    call mw%allreduce('sum', loss)
end subroutine

!> Prepare the forces difference to fit the third order
subroutine prepare_forces_thirdorer(fh, fc, mem)
    !> The fit helper
    class(lo_fit_helper), intent(inout) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The effective forces
    real(r8), dimension(:, :), allocatable :: feff
    !> The potential energy
    real(r8) :: epot
    !> Integers
    integer :: istep

    call mem%allocate(feff, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    allocate(fh%f2(3, fc%na, fh%mynsteps))
    allocate(fh%epot2(fh%mynsteps))
    do istep=1, fh%mynsteps
        call fc%forces(fh%u(:, :, istep), feff)
        epot = fc%potential_energy(fh%u(:, :, istep), feff)
        fh%f2(:, :, istep) = feff
        fh%epot2(istep) = epot
    end do
    call mem%deallocate(feff, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Do the gradient descent to fit the third order
subroutine optimize_thirdorder(fh, fc, uc, mw, mem)
    !> The fit helper
    class(lo_fit_helper), intent(inout) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_thirdorder), intent(inout) :: fc
    !> The structure
    type(lo_crystalstructure), intent(inout) :: uc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The gradient of the loss
!   real(r8), dimension(:, :, :, :), allocatable :: grad, gm, hm
    real(r8), dimension(:, :, :), allocatable :: grad, gm, hm
    !> The harmonic forces
    real(r8), dimension(:, :), allocatable :: feff, fdiff, uavg
    !> The measures for the gradient descent
    real(r8) :: maxgrad, maxavgf, loss, l0, rel_loss, gam
    !> Integers for do loop
    integer :: istep, itest, iat, a
    !> For better print
    character(len=1000) :: opf1, opf2

!   call mem%allocate(hm, [3, 3, 3, fc%ntriplets], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!   call mem%allocate(gm, [3, 3, 3, fc%ntriplets], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
!   call mem%allocate(grad, [3, 3, 3, fc%ntriplets], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(hm, [fc%nval, fc%na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(gm, [fc%nval, fc%na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(grad, [fc%nval, fc%na, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(feff, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(fdiff, [3, fc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! For better looking log
    opf1 = '(A7,5X,A19,5X,A11,5X,A15)'
    opf2 = '(I7,2X,E15.6,9X,E15.6,1X,E15.6)'

    ! Initialize
    init_coefs: block
        !> Random number generator
        type(lo_mersennetwister) :: rng
        !>
        real(r8) :: f0
        !> Some integer
        integer :: il, rnk

        if (mw%talk) write(lo_iou, *) '... Initializing coefficients randomly'

        f0 = sum(abs(fh%f2 * fh%u)) / sum(fh%u**2)

        ! We do it on one rank to be sure that we have the same initial value everywhere
        rnk = mw%n-1
        call rng%init(iseed=mw%r, rseed=walltime())
        if (mw%r .eq. rnk) then
            do il=1, fc%nval
                do iat=1, fc%na
                    do a=1, 3
                        fc%pvec(il, iat, a) = rng%rnd_gaussian_real(0.0_r8, f0)
                    end do
                end do
            end do
        end if
        call mw%bcast(fc%pvec, from=rnk)
    end block init_coefs

    fc%m = 0.0_r8

    ! Compute the gradient and loss before the start
    if (mw%talk) write(lo_iou, *) '... Computing gradient of least-squares loss'
    call fh%compute_gradient_thirdorder(fc, grad, feff, fdiff, loss, mem, mw)

    ! Initialize vectors for the Conjugate Gradient descent
    gm = -grad
    hm = -grad

    maxgrad = maxval(abs(grad))
    if (mw%talk) then
        write(lo_iou, opf1) 'Step', 'Max abs. grad. (eV)', 'Loss (eV/a)', 'Rel. Loss diff.'
        write(lo_iou, '(I7,2X,E15.6,9X,E15.6)') 0, maxgrad * lo_Hartree_to_eV, sqrt(loss) * lo_force_HartreeBohr_to_eVA
    end if
    gradloop: do istep=1, fh%nsteps_thirdorder
        ! We keep track of the previous loss
        l0 = loss

        ! We update the force constants
!       fc%pvec = fc%pvec - fh%alpha_thirdorder * grad
        fc%pvec = fc%pvec + fh%alpha_thirdorder * hm
!       fc%m = fc%m + fh%alpha_thirdorder * hm

        ! And we compute the gradient
        call fh%compute_gradient_thirdorder(fc, grad, feff, fdiff, loss, mem, mw)

        ! We update the Conjugate Gradient parameters
        gam = sum((grad + gm) * grad) / sum(gm * gm)
        gm = -grad
        hm = gm + gam * hm

        maxgrad = maxval(abs(grad))
        rel_loss = sqrt(abs(loss - l0)) / sqrt(l0)
!       if (istep .eq. fh%nsteps_thirdorder .or. lo_trueNtimes(istep, 10, fh%nsteps_thirdorder) .or. &
!           maxgrad .lt. fh%thresh .or. rel_loss .lt. 1e-12) then
            if (mw%talk) then
                write(lo_iou, opf2) istep, maxgrad * lo_Hartree_to_eV, sqrt(loss) * lo_force_HartreeBohr_to_eVA, rel_loss
            end if
!       end if
    end do gradloop

    call mem%deallocate(hm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(gm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(grad, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(feff, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(fdiff, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> The gradient of the IFC for the third order
subroutine compute_gradient_thirdorder(fh, fc, grad, feff, fdiff, loss, mem, mw)
    !> The fit helper
    class(lo_fit_helper), intent(in) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_thirdorder), intent(in) :: fc
    !> The gradient
!   real(r8), dimension(:, :, :, :) :: grad
    real(r8), dimension(:, :, :) :: grad
    !> The harmonic forces
    real(r8), dimension(:, :), intent(inout) :: feff
    !> The mean harmonic forces
    real(r8), dimension(:, :), intent(inout) :: fdiff
    !> The least-squares loss
    real(r8), intent(out) :: loss
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The displacements difference
    real(r8), dimension(3) :: dui, duj, duk
    !> Integer for do loops
    integer :: istep, iat, jat, kat, j, idx, a, b, c

    loss = 0.0_r8
    grad = 0.0_r8
    do istep=1, fh%mynsteps
        call fc%forces(fh%u(:, :, istep), feff)
        loss = loss + sum((fh%f(:, :, istep) - fh%f2(:, :, istep) - feff)**2) / real(fh%ntotsteps, r8) / real(fc%na, r8) / 3.0_r8
        fdiff = fh%f(:, :, istep) - fh%f2(:, :, istep) - feff
!       write(*, *) sum(abs(feff)) * lo_force_HartreeBohr_to_eVA / real(fc%na, r8) / 3.0_r8

        do j=1, fc%ntriplets
            iat = fc%ijk(1, j)
            jat = fc%ijk(2, j)
            kat = fc%ijk(3, j)

            do a=1, 3
            do b=1, 3
            do c=1, 3
                grad(:, iat, a) = grad(:, iat, a) + fdiff(a, iat) * fc%pvec(:, jat, b) * fc%pvec(:, kat, c) * &
                                                    (fh%u(b, jat, istep) - fh%u(b, iat, istep)) * &
                                                    (fh%u(c, kat, istep) - fh%u(c, iat, istep))

                grad(:, jat, b) = grad(:, jat, b) + fdiff(b, jat) * fc%pvec(:, iat, a) * fc%pvec(:, kat, c) * &
                                                    (fh%u(a, iat, istep) - fh%u(a, jat, istep)) * &
                                                    (fh%u(c, kat, istep) - fh%u(c, jat, istep))

                grad(:, kat, c) = grad(:, kat, c) + fdiff(c, kat) * fc%pvec(:, iat, a) * fc%pvec(:, jat, b) * &
                                                    (fh%u(a, iat, istep) - fh%u(a, kat, istep)) * &
                                                    (fh%u(b, jat, istep) - fh%u(b, kat, istep))

!               grad(a, b, c, j) = grad(a, b, c, j) + fdiff(a, iat) * (fh%u(b, jat, istep) - fh%u(b, iat, istep)) * &
!                                                                     (fh%u(c, kat, istep) - fh%u(c, iat, istep))
!               grad(a, b, c, j) = grad(a, b, c, j) + fdiff(b, jat) * (fh%u(a, iat, istep) - fh%u(a, jat, istep)) * &
!                                                                     (fh%u(c, kat, istep) - fh%u(c, jat, istep))
!               grad(a, b, c, j) = grad(a, b, c, j) + fdiff(c, kat) * (fh%u(a, iat, istep) - fh%u(a, kat, istep)) * &
!                                                                     (fh%u(b, jat, istep) - fh%u(b, kat, istep))

!               grad(a, b, c, j) = grad(a, b, c, j) + fdiff(a, iat) * (fh%u(b, jat, istep) - fh%u(b, iat, istep)) * &
!                                                                     (fh%u(c, kat, istep) - fh%u(c, iat, istep))
!               grad(a, b, c, j) = grad(a, b, c, j) + fdiff(a, jat) * (fh%u(b, iat, istep) - fh%u(b, jat, istep)) * &
!                                                                     (fh%u(c, kat, istep) - fh%u(c, jat, istep))
!               grad(a, b, c, j) = grad(a, b, c, j) + fdiff(a, kat) * (fh%u(b, iat, istep) - fh%u(b, kat, istep)) * &
!                                                                     (fh%u(c, jat, istep) - fh%u(c, kat, istep))
            end do
            end do
            end do
        end do
    end do
    grad = grad / real(fh%ntotsteps, r8)
    call mw%allreduce('sum', grad)
    call mw%allreduce('sum', loss)
end subroutine

!> Compute the average positions given the simulation
subroutine get_average_positions(fh, uc, mem, mw)
    !> The fit helper
    class(lo_fit_helper), intent(inout) :: fh
    !> The cell
    type(lo_crystalstructure), intent(inout) :: uc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The positions, in fractional coordinates
    real(r8), dimension(:, :), allocatable :: uavg
    !> Some integers
    integer :: istep, iat

    call mem%allocate(uavg, [3, uc%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    uavg = 0.0_r8

    ! First we compute the average displacements
    do istep=1, fh%mynsteps
        uavg = uavg + fh%u(:, :, istep) / real(fh%ntotsteps, r8)
    end do
    call mw%allreduce('sum', uavg)

    ! Now we update the crystal structure
    uc%rcart = uc%rcart + uavg
    do iat=1, uc%na
        uc%r(:, iat) = uc%cartesian_to_fractional(uc%rcart(:, iat))
    end do

    ! And we update the displacements
    do istep=1, fh%mynsteps
        fh%u(:, :, istep) = fh%u(:, :, istep) - uavg
    end do
    call mem%deallocate(uavg, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Compute several diagnostics to check for the fit
subroutine compute_diagnostic(fh, fc2, fc3, u0, u3, res_e, res_f, thirdorder, mem, mw)
    !> The fit helper
    class(lo_fit_helper), intent(in) :: fh
    !> The force constants
    type(lo_forceconstant_disorder_secondorder), intent(in) :: fc2
    !> The force constants
    type(lo_forceconstant_disorder_thirdorder), intent(in) :: fc3
    !> The mean energy difference
    real(r8), intent(out) :: u0
    !> The mean energy difference
    real(r8), intent(out) :: u3
    !> The results for the energy
    real(r8), dimension(3, 2), intent(out) :: res_e
    !> The results for the forces
    real(r8), dimension(3, 2), intent(out) :: res_f
    !> Do we need to compute everything for third order ?
    logical, intent(in) :: thirdorder
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> Some buffers
    real(r8), dimension(:, :), allocatable :: fharm, f3
    !>
    real(r8) :: eharm, norm, mean_e, emean_harm, e3, emean3
    !> Some integers
    integer :: istep

    res_e = 0.0_r8
    res_f = 0.0_r8
    call mem%allocate(fharm, [3, fc2%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (thirdorder) then
        call mem%allocate(f3, [3, fc2%na], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if

    ! Let's get the normalization already
    norm = 1.0_r8 / real(fh%ntotsteps, r8)

    ! We have to compute simple averages first
    mean_e = 0.0_r8
    emean_harm = 0.0_r8
    emean3 = 0.0_r8
    u0 = 0.0_r8
    u3 = 0.0_r8
    do istep=1, fh%mynsteps
        eharm = fc2%potential_energy(fh%u(:, :, istep))

        mean_e = mean_e + fh%epot(istep) * norm
        emean_harm = emean_harm + eharm * norm
        u0 = u0 + (fh%epot(istep) - eharm) * norm

        if (thirdorder) then
            e3 = fc3%potential_energy(fh%u(:, :, istep))
            emean3 = emean3 + (eharm + e3) * norm
            u3 = u3 + (fh%epot(istep) - eharm - e3) * norm
        end if
    end do
    call mw%allreduce('sum', mean_e)
    call mw%allreduce('sum', emean_harm)
    call mw%allreduce('sum', emean3)
    call mw%allreduce('sum', u0)
    call mw%allreduce('sum', u3)

    res_e = 0.0_r8
    res_f = 0.0_r8
    do istep=1, fh%mynsteps
        call fc2%forces(fh%u(:, :, istep), fharm)
        eharm = fc2%potential_energy(fh%u(:, :, istep), fharm)
        ! Now we start computing
        ! First the energy part
        res_e(1, 1) = res_e(1, 1) + (fh%epot(istep) - mean_e)**2 * norm
        res_e(1, 2) = res_e(1, 2) + (fh%epot(istep) - mean_e)**2 * norm

        res_e(2, 1) = res_e(2, 1) + (eharm - emean_harm)**2 * norm
        res_e(2, 2) = res_e(2, 2) + (fh%epot(istep) - eharm - u0)**2 * norm

        ! Now the forces part
        res_f(1, 1) = res_f(1, 1) + sum(fh%f(:, :, istep)**2) * norm
        res_f(1, 2) = res_f(1, 2) + sum(fh%f(:, :, istep)**2) * norm

        res_f(2, 1) = res_f(2, 1) + sum(fharm**2) * norm
        res_f(2, 2) = res_f(2, 2) + sum((fh%f(:, :, istep) - fharm)**2) * norm

        ! Now the third order
        if (thirdorder) then
            call fc3%forces(fh%u(:, :, istep), f3)
            e3 = fc3%potential_energy(fh%u(:, :, istep), f3)

            ! The energy
            res_e(3, 1) = res_e(3, 1) + (eharm + e3 - emean3)**2 * norm
            res_e(3, 2) = res_e(3, 2) + (fh%epot(istep) - eharm - e3 - u3)**2 * norm

            ! The forces
            res_f(3, 1) = res_f(3, 1) + sum((fharm + f3)**2) * norm
            res_f(3, 2) = res_f(3, 2) + sum((fh%f(:, :, istep) - fharm - f3)**2) * norm
        end if
    end do
    call mw%allreduce('sum', res_e)
    call mw%allreduce('sum', res_f)
    ! Now we need to normaliz things
    res_e = sqrt(res_e / real(fc2%na, r8))
    res_f = sqrt(res_f / real(fc2%na, r8) / 3.0_r8)
    u0 = u0 / real(fc2%na, r8)
    call mem%deallocate(fharm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (thirdorder) then
        u3 = u3 / real(fc3%na, r8)
        call mem%deallocate(f3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end if
end subroutine
end module
