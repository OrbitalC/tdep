! submodule(type_forceconstant_disorder_secondorder) type_forceconstant_disorder_secondorder_dynamicalmatrix
module type_forceconstant_disorder_secondorder_dynamicalmatrix
!!
!! Sort out everything related to dynamical matrix with disordered force constants
!!
use konstanter, only: lo_freqtol, lo_imag, lo_twopi, lo_sqtol, lo_groupvel_Hartreebohr_to_ms
use gottochblandat, only: qsort, lo_negsqrt, lo_real_gram_schmidt, lo_complex_gram_schmidt, &
                          lo_chop, lo_sqnorm
use type_blas_lapack_wrappers, only: lo_dsyev, lo_zheev, lo_dgemm
use lo_sorting, only: lo_qsort

use konstanter, only: r8, lo_exitcode_param, lo_exitcode_baddim, lo_exitcode_blaslapack
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_disorder_secondorder, only: lo_forceconstant_disorder_secondorder
use gottochblandat, only: tochar
use lo_memtracker, only: lo_mem_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use type_qpointmesh, only: lo_qpoint
use type_unitcell_disorder_projection, only: lo_unitcell_disorder_projection

use gottochblandat, only: lo_invert3x3matrix

implicit none
private
public :: dynamicalmatrix
public :: frequencies_eigenvectors_groupvelocities
public :: dynamicalmatrix_structurefactor
public :: eigensolve_structurefactor
public :: dynamicalmatrix_projected
public :: eigensolve_projected
contains

! Compute the dynamical matrix in real space
module subroutine dynamicalmatrix(fc, ss, dynmat, dynmat_grad)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: ss
    !> The dynamical matrix
    real(r8), dimension(:, :), intent(out) :: dynmat
    !> The frequencies
    real(r8), dimension(:, :, :), intent(out), optional :: dynmat_grad

    !> Buffer for the i-j IFC
    real(r8), dimension(3, 3) :: m
    !> Little buffer
    real(r8) :: f0
    !> Some integer for do loops and indices
    integer :: iat, jat, j, ii, jj, a, b

    ! First we fill everything
    dynmat = 0.0_r8
    if (present(dynmat_grad)) dynmat_grad = 0.0_r8
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        m = fc%m(:, :, j) * ss%invsqrtmass(iat) * ss%invsqrtmass(jat)

        ! First we fill the dynamical matrix, and its gradient if needed
        do a=1, 3
        do b=1, 3
            ii = 3 * (iat - 1) + a
            jj = 3 * (jat - 1) + b
            dynmat(ii, jj) = dynmat(ii, jj) + m(a, b)
            dynmat(jj, ii) = dynmat(jj, ii) + m(b, a)
            if (present(dynmat_grad)) then
                dynmat_grad(:, ii, jj) = dynmat_grad(:, ii, jj) + m(a, b) * fc%rij(:, j)
                dynmat_grad(:, jj, ii) = dynmat_grad(:, jj, ii) - m(b, a) * fc%rij(:, j)
            end if

            ! Acoustic sum rule
            ! For the first atom
            ii = 3 * (iat - 1) + a
            jj = 3 * (iat - 1) + b
            dynmat(ii, jj) = dynmat(ii, jj) - m(a, b)
            ! For the second atom, it's the transpose
            ii = 3 * (jat - 1) + a
            jj = 3 * (jat - 1) + b
            dynmat(ii, jj) = dynmat(ii, jj) - m(b, a)
        end do
        end do
    end do

    ! Make sure the dynamical matrix is perfectly symmetric
    ! If everything was right before, this shouldn't change anything
    do ii=1, fc%na*3-1
        do jj=ii+1, fc%na*3
            f0 = 0.5_r8 * (dynmat(ii, jj) + dynmat(jj, ii))
            dynmat(ii, jj) = f0
            dynmat(jj, ii) = f0
            ! And antisymmetric for its gradient
            if (present(dynmat_grad)) then
                do a=1, 3
                    f0 = 0.5_r8 * (dynmat_grad(a, ii, jj) - dynmat_grad(a, jj, ii))
                    dynmat_grad(a, ii, jj) = f0
                    dynmat_grad(a, jj, ii) = -f0
                end do
            end if
        end do
    end do
end subroutine

! Compute the frequencies eigenvectors and group velocities from the dynamical matrix
module subroutine frequencies_eigenvectors_groupvelocities(fc, dynmat, dynmat_grad, frequencies, eigenvectors, groupvel, mem)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The dynamical matrix
    real(r8), dimension(:, :), intent(in) :: dynmat
    !> The gradient of the dynamical matrix
    real(r8), dimension(:, :, :), intent(in), optional :: dynmat_grad
    !> The frequencies
    real(r8), dimension(:) :: frequencies
    !> The eigenvectors
    real(r8), dimension(:, :), intent(out) :: eigenvectors
    !> The group velocities
    real(r8), dimension(:, :, :), intent(out), optional :: groupvel
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The index of the gamma points
    integer, dimension(3) :: ind_gamma
    !> Some do loops integers
    integer :: nb

    nb = fc%na * 3
    frequencies = 0.0_r8
    eigenvectors = 0.0_r8
    ! Do we compute eigenvectors ?
    if (present(groupvel)) then
        if (.not. present(dynmat_grad)) then
            call lo_stop_gracefully(['You need the gradient of the dynamical matrix to get the group velocities'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
        groupvel = 0.0_r8
    end if

    ! First we do the actual diagonalization
    diagonalize: block
        !> Buffer for the result of the diagonalization
        real(r8), dimension(:, :), allocatable :: wdynmat
        !> Buffer for the frequencies
        real(r8), dimension(:), allocatable :: dr
        !> Buffer for the sorting
        integer, dimension(:), allocatable :: ind
        !> For the do loops
        integer :: i

        ! We use dsyev. Would dsyevr be better ?
        call mem%allocate(wdynmat, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        wdynmat = dynmat
        call lo_dsyev(wdynmat, frequencies, 'V', 'U', info=i)
        if (i .ne. 0) then
            call lo_stop_gracefully(['dsyev exit code'//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
        end if

        ! Fix acoustic modes
        call mem%allocate(dr, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ind, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        dr = abs(frequencies)
        call qsort(dr, ind)
        ! We might need the index of the gamma point later
        ind_gamma = ind(1:3)
        do i=1, 3
            frequencies(ind(i)) = 0.0_r8
        end do
        ! We take the square root to get the frequencies
        frequencies = lo_negsqrt(frequencies)
        ! Sort eigenvalues
        call qsort(frequencies, ind)
        ! And eigenvectors
        do i=1, nb
            eigenvectors(:, i) = wdynmat(:, ind(i))
        end do
        call mem%deallocate(wdynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ind, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! Orthogonalize the eigenvectors
        call lo_real_gram_schmidt(eigenvectors)

        ! @TODO low probability of it happening, but would be useful to disentangle degenerate bands
    end block diagonalize

    ! If needed, we get the group velocities
    ! Get the group velocities
    if (present(groupvel)) then
    grad: block
        !> The scaled eigenvectors
        real(r8), dimension(:, :), allocatable :: scaledeig, buf
        !> For the do loops
        integer :: i

        ! Some allocation
        call mem%allocate(buf, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(scaledeig, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        scaledeig = 0.0_r8
        buf = 0.0_r8

        ! We scale the eigenvectors before applying them to the dynamical matrix gradient
        do i=1, nb
            ! We don't want the acoustic modes
            if (frequencies(i) .lt. lo_freqtol) cycle
            scaledeig(:, i) = eigenvectors(:, i) / sqrt(frequencies(i))
        end do

        ! Now we can compute the group velocities, simple application of the formula
        do i=1, 3
            call lo_dgemm(dynmat_grad(i, :, :), scaledeig, buf)
            call lo_dgemm(buf, scaledeig, groupvel(i, :, :))
        end do
        ! We set the group velocities involving the acoustic modes to zero
        do i=1, 3
            groupvel(:, ind_gamma(i), :) = 0.0_r8
            groupvel(:, :, ind_gamma(i)) = 0.0_r8
        end do
        call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(scaledeig, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block grad
    end if
end subroutine

!> Compute the dynamical matrix for the structure factor
module subroutine dynamicalmatrix_structurefactor(fc, qpoint, dynmat, dynmat_grad, ss)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The qpoint
    class(lo_qpoint), intent(in) :: qpoint
    !> The dynamical matrix
    complex(r8), dimension(3, 3), intent(out) :: dynmat
    !> The dynamical matrix gradient
    complex(r8), dimension(3, 3, 3), optional, intent(out) :: dynmat_grad
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: ss

    !> Buffer for the dynamical matrix gradient
    complex(r8), dimension(3, 3) :: gm
    !> The Fourier transform phase
    complex(r8) :: expiqr
    !> Pair ifc
    real(r8), dimension(3, 3) :: m, asr
    !> Rij and q-vectors
    real(r8), dimension(3) :: qv, rij
    !> Dot product of q and r
    real(r8) :: qdotr
    !> Some integers for do loops
    integer :: iat, jat, j, ic, a

    dynmat = 0.0_r8
    if (present(dynmat_grad)) dynmat_grad = 0.0_r8

    ! We don't need to do anything at the gamma point
    if (lo_sqnorm(qpoint%r) .lt. lo_sqtol) return

    ! Else, we can fill everything
    qv = lo_chop(qpoint%r * lo_twopi, lo_sqtol)
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        rij = fc%rij(:, j)

        ! Compute the Fourier phase
        qdotr = dot_product(rij, qv)
        expiqr = cmplx(cos(qdotr), sin(qdotr), r8)

        ! Get the mass scaled IFC tensor for this pair
        m = fc%m(:, :, j) * ss%invsqrtmass(iat) * ss%invsqrtmass(jat) / real(fc%na, r8)

        ! Fill the dynamical matrix, with hermicity
        dynmat = dynmat + (m * expiqr + transpose(m) * conjg(expiqr))
        ! And the acoustic sum rule
        dynmat = dynmat - (m + transpose(m))

        ! And if needed, the gradient of the dynamical matrix
        if (present(dynmat_grad)) then
            do a=1, 3
                gm = lo_imag * rij(a) * m * expiqr
                dynmat_grad(a, :, :) = dynmat_grad(a, :, :) + 0.5_r8 * (gm - conjg(transpose(gm)))
            end do
        end if
    end do

    ! Make sure it's Hermitian
    dynmat = 0.5_r8 * (dynmat + conjg(transpose(dynmat)))
    if (present(dynmat_grad)) then
        ! And here, anti-Hermitian
        do a=1, 3
            dynmat_grad(a, :, :) = 0.5_r8 * (dynmat_grad(a, :, :) - conjg(transpose(dynmat_grad(a, :, :))))
        end do
    end if
end subroutine

!> Get the frequencies, eigenvector and maybe group velocities for the structure factor
module subroutine eigensolve_structurefactor(fc, qpoint, dynmat, dynmat_grad, frequencies, eigenvectors, groupvel)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The qpoint
    class(lo_qpoint), intent(in) :: qpoint
    !> The dynamical matrix
    complex(r8), dimension(3, 3), intent(in) :: dynmat
    !> The dynamical matrix gradient
    complex(r8), dimension(3, 3, 3), optional, intent(in) :: dynmat_grad
    !> The frequencies
    real(r8), dimension(3), intent(out) :: frequencies
    !> The eigenvectors
    complex(r8), dimension(3, 3), intent(out) :: eigenvectors
    !> The group velocities
    real(r8), dimension(3, 3, 3), optional, intent(out) :: groupvel

    ! Buf eigenvectors
    complex(r8), dimension(3, 3) :: bufeigenvec, buf
    ! Buf eigenvalues
    real(r8), dimension(3) :: weigenval
    !> Index to sort
    integer, dimension(3) :: ind
    !> Some integer
    integer :: i

    ! Set everything to zero
    frequencies = 0.0_r8
    eigenvectors = 0.0_r8
    if (present(groupvel)) then
        ! Check that there is no inconsistency
        if (.not. present(dynmat_grad)) then
            call lo_stop_gracefully(['You need the derivative of the dynamical matrix to compute the group velocities'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
        groupvel = 0.0_r8
    end if

    ! If it's the gamma point, we can already stop
    if (lo_sqnorm(qpoint%r) .lt. lo_sqtol) return

    ! Do the actual diagonalization
    bufeigenvec = dynmat
    call lo_zheev(bufeigenvec, weigenval, jobz='V', info=i)
    if (i .ne. 0) call lo_stop_gracefully(['zheev exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)
    ! Orthogonalize eigenvectors
    call lo_complex_gram_schmidt(bufeigenvec)

    ! Sort frequencies, and put inside results array
    weigenval = lo_negsqrt(weigenval)
    call qsort(weigenval, ind)
    do i=1, 3
        frequencies(i) = weigenval(i)
        eigenvectors(:, i) = bufeigenvec(:, ind(i))
    end do

    ! Compute group velocities ?
    if (present(groupvel)) then
        bufeigenvec = 0.0_r8
        buf = 0.0_r8
        ! First we scale the frequencies
        do i=1, 3
            bufeigenvec(:, i) = eigenvectors(:, i) / sqrt(frequencies(i))
        end do
        ! Now we simply apply the formula
        do i=1, 3
            buf = matmul(dynmat_grad(i, :, :), bufeigenvec)
            groupvel(i, :, :) = 0.5_r8 * real(matmul(bufeigenvec, buf), r8)
        end do
    end if
end subroutine

!> Compute the dynamical matrix projected on a unitcell
module subroutine dynamicalmatrix_projected(fc, proj, qpoint, dynmat, dynmat_grad, ss)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The projection
    type(lo_unitcell_disorder_projection), intent(in) :: proj
    !> The qpoint
    class(lo_qpoint), intent(in) :: qpoint
    !> The dynamical matrix
    complex(r8), dimension(:, :), intent(out) :: dynmat
    !> The dynamical matrix gradient
    complex(r8), dimension(:, :, :), optional, intent(out) :: dynmat_grad
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: ss

    !> Pair ifc
    real(r8), dimension(3, 3) :: m, asr
    !> Rij and q-vectors
    real(r8), dimension(3) :: qv, rij
    !> Buffer for the Fourier phase and the gradient of the dynamical matrix
    complex(r8) :: expiqr, c0
    !> Dot product of q and r and a buffer
    real(r8) :: qdotr, f0
    !> Some integer for do loops and index
    integer :: iat, jat, ii, jj, j, a, b, c, iat_uc, jat_uc

    dynmat = 0.0_r8
    if (present(dynmat_grad)) dynmat_grad = 0.0_r8

    ! Prepare some buffer
    f0 = real(proj%na_uc, r8) / real(proj%na_ss, r8)
    qv = lo_chop(qpoint%r * lo_twopi, 1e-13_r8)

    ! And here we go
    do j=1, fc%npairs
        ! Get info on iat, jat and their distance
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        rij = fc%rij(:, j)

        ! Compute the Fourier phase
        qdotr = dot_product(rij, qv)
        expiqr = cmplx(cos(qdotr), sin(qdotr), r8)

        ! Get the atom in the unitcell
        iat_uc = proj%unitcell_index(iat)
        jat_uc = proj%unitcell_index(jat)

        ! And the pair IFC
        m = fc%m(:, :, j) * ss%invsqrtmass(iat) * ss%invsqrtmass(jat) * f0
        do a=1, 3
        do b=1, 3
            ii = 3 * (iat_uc - 1) + a
            jj = 3 * (jat_uc - 1) + b
            dynmat(ii, jj) = dynmat(ii, jj) + m(a, b) * expiqr
            dynmat(jj, ii) = dynmat(jj, ii) + m(b, a) * conjg(expiqr)
            if (present(dynmat_grad)) then
                c0 = lo_imag * rij(c) * m(a, b) * expiqr
                do c=1, 3
                    dynmat_grad(c, ii, jj) = dynmat_grad(c, ii, jj) + c0
                    dynmat_grad(c, jj, ii) = dynmat_grad(c, jj, ii) + conjg(c0)
                end do
            end if

            ! Now for the acoustic sum rule
            ! For the first atom, normal order
            ii = 3 * (iat_uc - 1) + a
            jj = 3 * (iat_uc - 1) + b
            dynmat(ii, jj) = dynmat(ii, jj) - m(a, b)
            ! For the second atom it's the transpose
            ii = 3 * (jat_uc - 1) + a
            jj = 3 * (jat_uc - 1) + b
            dynmat(ii, jj) = dynmat(ii, jj) - m(b, a)
        end do
        end do
    end do

    ! Make sure it's Hermitian
    do ii=1, 3 * proj%na_uc
    do jj=ii, 3 * proj%na_uc
        c0 = 0.5_r8 * (dynmat(ii, jj) + conjg(dynmat(jj, ii)))
        dynmat(ii, jj) = c0
        dynmat(jj, ii) = c0
        if (present(dynmat_grad)) then
            do c=1, 3
                c0 = 0.5_r8 * (dynmat_grad(c, ii, jj) - conjg(dynmat_grad(c, jj, ii)))
                dynmat_grad(c, ii, jj) = c0
                dynmat_grad(c, jj, ii) = -conjg(c0)
            end do
        end if
    end do
    end do
end subroutine

!> Get the phonons projected on a unitcell basis
subroutine eigensolve_projected(fc, qpoint, dynmat, dynmat_grad, omega, eigenvectors, groupvel, mem)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The qpoint
    type(lo_qpoint), optional, intent(in) :: qpoint
    !> The dynamical matrix
    complex(r8), dimension(:, :), intent(in) :: dynmat
    !> The dynamical matrix gradient
    complex(r8), dimension(:, :, :), optional, intent(in) :: dynmat_grad
    !> The frequencies
    real(r8), dimension(:), intent(out) :: omega
    !> The eigenvectors
    complex(r8), dimension(:, :), optional, intent(out) :: eigenvectors
    !> The group velocities
    real(r8), dimension(:, :, :), optional, intent(out) :: groupvel
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> What do we compute ?
    logical :: calceig, calcvel
    !> The group velocities
    complex(r8), dimension(:, :, :), allocatable :: wVelocities
    !> The Eigenvectors
    complex(r8), dimension(:, :), allocatable :: wEigenvector
    !> The eigenvalues
    real(r8), dimension(:), allocatable :: wEigenval
    !> The number of modes
    integer :: nb

    nb = size(dynmat, 1)

    ! Return the eigenvectors ?
    if (present(eigenvectors)) then
        calceig = .true.
    else
        calceig = .false.
    end if

    ! Return group velocities ?
    if (present(groupvel)) then
        calcvel = .true.
        if (.not. present(dynmat_grad)) then
            call lo_stop_gracefully(['Need the gradient of the dynamical matrix to compute group velocities'], &
                                    lo_exitcode_param, __FILE__, __LINE__)
        end if
    else
        calcvel = .false.
    end if

    call mem%allocate(wEigenval, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(wEigenvector, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    diagonalize: block
        !> The dynamical matrix
        complex(r8), dimension(:, :), allocatable :: wDynmat
        !> The eigenvalues
        real(r8), dimension(:), allocatable :: dr
        !> Some indices for the ordering
        integer, dimension(:), allocatable :: ind
        !> Some integers for the do loops
        integer :: i

        call mem%allocate(wDynmat, [nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(dr, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ind, nb, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        wDynmat = 0.0_r8
        dr = 0.0_r8
        ind = 0.0_r8

        wEigenvector = dynmat
        call lo_zheev(wEigenvector, wEigenval, jobz='V', info=i)
        if (i .ne. 0) then
            call lo_stop_gracefully(['zheev exit status '//tochar(i)], &
                                    lo_exitcode_blaslapack, __FILE__, __LINE__)
        end if
        call lo_complex_gram_schmidt(wEigenvector)

        ! Set the acoustic modes to zero if we are at gamma
        if (present(qpoint)) then
            if (lo_sqnorm(qpoint%r) .lt. lo_sqtol) then
                dr = abs(wEigenval)
                call lo_qsort(dr, ind)
                do i=1, 3
                    wEigenval(ind(i)) = 0.0_r8
                end do
            end if
        end if

        ! Then sort
        wDynmat = wEigenvector
        dr = lo_negsqrt(wEigenval)
        call lo_qsort(dr, ind)
        dr = wEigenval
        do i=1, nb
            wEigenval(i) = dr(ind(i))
            wEigenvector(:, i) = wDynmat(:, ind(i))
        end do

        call mem%deallocate(wDynmat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(dr, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ind, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block diagonalize

    calc_groupvel: block
        if (calcvel) then
            call mem%allocate(wVelocities, [3, nb, nb], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
    end block calc_groupvel

    returnvals: block
        real(r8) :: f0
        f0 = 1e-9_r8 / lo_groupvel_Hartreebohr_to_ms
        omega = lo_negsqrt(wEigenval)
        if (present(eigenvectors)) eigenvectors = wEigenvector
        if (present(groupvel)) groupvel = lo_chop(real(wVelocities), f0)

        call mem%deallocate(wEigenval, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(wEigenvector, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        if (present(groupvel)) then
            call mem%deallocate(wVelocities, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        end if
    end block returnvals
end subroutine

! end submodule
end module
