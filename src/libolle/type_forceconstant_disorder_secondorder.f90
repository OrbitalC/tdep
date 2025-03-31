module type_forceconstant_disorder_secondorder
!!
!! Container for the second order forceconstants for disordered systems
!!
use konstanter, only: r8, lo_huge, lo_exitcode_baddim, lo_exitcode_blaslapack, lo_exitcode_param,  &
                      lo_tol, lo_hugeint
use gottochblandat, only: tochar
use hdf5_wrappers, only: lo_hdf5_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint
use type_crystalstructure, only: lo_crystalstructure
use type_distancetable, only: lo_distancetable
use type_blas_lapack_wrappers, only: lo_dposv
use type_unitcell_disorder_projection, only: lo_unitcell_disorder_projection

implicit none

private
public :: lo_forceconstant_disorder_secondorder

! Get some per atom things
type lo_forceconstant_disorder_pair
    !> The number of neighbors
    integer :: n=-lo_hugeint
    !> The index for the interaction in the flatten IFC
    integer, dimension(:), allocatable :: ic
    !> The constraint matrix
    real(r8), dimension(:, :), allocatable :: cstr
    !> The c*c^T matrix for the constraint
    real(r8), dimension(:, :), allocatable :: cct
end type

! A flatten version of the force constants, without symmetries
type lo_forceconstant_disorder_secondorder
    !> The number of atoms
    integer :: na
    !> The total number of coefs
    integer :: npairs
    !> The list of coefficients
    real(r8), dimension(:, :, :), allocatable :: m
    !> The cutoff
    real(r8) :: cutoff=-lo_huge
    !> The min distance in the system
    real(r8) :: mindist=-lo_huge
    !> map coefficients to the atomic pair
    integer, dimension(:, :), allocatable :: ij
    !> All the distance vectors
    real(r8), dimension(:, :), allocatable :: rij
    !> The pairs
    type(lo_forceconstant_disorder_pair), dimension(:), allocatable :: atom
    !> The elastic constants tensor
    real(r8), dimension(3, 3, 3, 3) :: elastic_constants_tensor = lo_huge
    !> The elastic constants in Voigt notation
    real(r8), dimension(6, 6) :: elastic_constants_voigt = lo_huge

contains
    !> Get the neighbors table given a unitcell
    procedure :: get_pairs
    !> Create pairs to allow for per atom access to IFC
    procedure :: create_pairs
    !> Check that a given structure agrees with the IFC
    procedure :: check_structure
    !> Check the Hermitian sum rule
    procedure :: check_hermitian
    !> Check the rotational sum rule
    procedure :: check_rotational
    !> Check the Huang sum rule
    procedure :: check_huang
    !> Apply Hermitian constraints
    procedure :: apply_hermitian_constraints
    !> Create the Hermitian constraints matrix
    procedure :: create_hermitian_constraints
    !> Compute the potential energy
    procedure :: potential_energy
    !> Compute the forces on each atom
    procedure :: forces
!   !> Return a dynamical matrix
!   procedure :: dynamicalmatrix
!   !> Compute the phonons
!   procedure :: frequencies_eigenvectors_groupvelocities
!   !> Return the dynamical matrix for a structure factor
!   procedure :: dynamicalmatrix_structurefactor
!   !> Return the eigensolution of structure factor
!   procedure :: eigensolve_structurefactor
    !> Get the elastic constants
    procedure :: get_elastic_constants
    !> Write the forceconstants to hdf5
    procedure :: write_to_hdf5
    !> And read the forceconstants from hdf5
    procedure :: read_from_hdf5
end type

! Interfaces to type_forceconstant_disorder_secondorder_dynamicalmatrix
interface
! Compute the dynamical matrix in real space
    module subroutine dynamicalmatrix(fc, uc, dynmat, dynmat_grad)
        class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
        type(lo_crystalstructure), intent(in) :: uc
        real(r8), dimension(:, :), intent(out) :: dynmat
        real(r8), dimension(:, :, :), intent(out), optional :: dynmat_grad
    end subroutine
    ! Compute the frequencies eigenvectors and group velocities from the dynamical matrix
    module subroutine frequencies_eigenvectors_groupvelocities(fc, dynmat, dynmat_grad, frequencies, eigenvectors, groupvel, mem)
        class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
        real(r8), dimension(:, :), intent(in) :: dynmat
        real(r8), dimension(:, :, :), intent(in), optional :: dynmat_grad
        real(r8), dimension(:) :: frequencies
        real(r8), dimension(:, :), intent(out) :: eigenvectors
        real(r8), dimension(:, :, :), intent(out), optional :: groupvel
        type(lo_mem_helper), intent(inout) :: mem
    end subroutine
    !> Compute the dynamical matrix for the structure factor
    module subroutine dynamicalmatrix_structurefactor(fc, qpoint, dynmat, dynmat_grad, uc)
        class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
        class(lo_qpoint), intent(in) :: qpoint
        complex(r8), dimension(3, 3), intent(out) :: dynmat
        complex(r8), dimension(3, 3, 3), optional, intent(out) :: dynmat_grad
        type(lo_crystalstructure), intent(in) :: uc
    end subroutine
    !> Get the frequencies, eigenvector and maybe group velocities for the structure factor
    module subroutine eigensolve_structurefactor(fc, qpoint, dynmat, dynmat_grad, frequencies, eigenvectors, groupvel)
        class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
        class(lo_qpoint), intent(in) :: qpoint
        complex(r8), dimension(3, 3), intent(in) :: dynmat
        complex(r8), dimension(3, 3, 3), optional, intent(in) :: dynmat_grad
        real(r8), dimension(3), intent(out) :: frequencies
        complex(r8), dimension(3, 3), intent(out) :: eigenvectors
        real(r8), dimension(3, 3, 3), optional, intent(out) :: groupvel
    end subroutine
    !> Compute the dynamical matrix projected on a unitcell
    module subroutine dynamicalmatrix_projected(fc, proj, qpoint, dynmat, dynmat_grad, uc)
        class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
        type(lo_unitcell_disorder_projection), intent(in) :: proj
        class(lo_qpoint), intent(in) :: qpoint
        complex(r8), dimension(3, 3), intent(out) :: dynmat
        complex(r8), dimension(3, 3, 3), optional, intent(out) :: dynmat_grad
        type(lo_crystalstructure), intent(in) :: uc
    end subroutine

end interface

contains
! Create the neighbor list and tables for the force constants
subroutine get_pairs(fc, uc, cutoff, mem)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The cutoff
    real(r8), intent(in) :: cutoff
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The distance table
    type(lo_distancetable) :: dt
    !> Some integers
    integer :: iat, jat, j, k, ii

    fc%cutoff = cutoff
    fc%na = uc%na
    fc%mindist = lo_huge

    !> First, we have to create a distance table
    call dt%generate(uc%r, uc%latticevectors, fc%cutoff, verbosity=0)

    ! Now we count the number of coefficients
    ii = 0
    do iat=1, fc%na
        do j=1, dt%particle(iat)%n
            jat = dt%particle(iat)%ind(j)
            if (jat .lt. iat) ii = ii + 1
        end do
    end do
    fc%npairs = ii
    allocate(fc%m(3, 3, fc%npairs))
    allocate(fc%ij(2, fc%npairs))
    allocate(fc%rij(3, fc%npairs))
    fc%m = 0.0_r8
    fc%ij = 0
    fc%rij = 0.0_r8

    ! Now we map the coefficients to the interactions
    ii = 0
    do iat=1, fc%na
        do j=1, dt%particle(iat)%n
            jat = dt%particle(iat)%ind(j)
            if (jat .lt. iat) then
                ii = ii + 1
                fc%ij(1, ii) = iat
                fc%ij(2, ii) = jat
                fc%rij(:, ii) = dt%particle(iat)%v(:, j)
                fc%mindist = min(fc%mindist, norm2(fc%rij(:, ii)))
            end if
        end do
    end do
    ! A little sanity check
    if (fc%npairs .ne. ii) then
        call lo_stop_gracefully(['Number of neighbor do not match the distance table'], &
                                lo_exitcode_baddim, __FILE__, __LINE__)
    end if
end subroutine

!> Create pairs to allow for per atom access to IFC
subroutine create_pairs(fc, mem)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> To count the number of pair per atom
    integer, dimension(:), allocatable :: npair
    !> Some integers
    integer :: j, iat, jat

    allocate(fc%atom(fc%na))
    call mem%allocate(npair, fc%na, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    npair = 0
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        npair(iat) = npair(iat) + 1
        npair(jat) = npair(jat) + 1
    end do
    ! Now we can allocate the number of neighbors per pair
    do iat=1, fc%na
        fc%atom(iat)%n = npair(iat)
        allocate(fc%atom(iat)%ic(npair(iat)))
    end do
    ! And we fill
    npair = 0
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        npair(iat) = npair(iat) + 1
        npair(jat) = npair(jat) + 1
        ! Plus sign for direct, minus sign for transpose
        fc%atom(iat)%ic(npair(iat)) = j
        fc%atom(jat)%ic(npair(jat)) = -j
    end do
    call mem%deallocate(npair, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Check that the IFC correspond a given structure
subroutine check_structure(fc, uc, isok)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Does everything checks out ?
    logical, intent(out) :: isok

    !> The distance table
    type(lo_distancetable) :: dt
    !> The vector distance between iat and jat
    real(r8), dimension(3) :: rij
    !> Some integers
    integer :: iat, jat, j, k, n

    !> First, we have to create a distance table
    call dt%generate(uc%r, uc%latticevectors, fc%cutoff, verbosity=0)

    n = 0
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        rij = fc%rij(:, j)
        do k=1, dt%particle(iat)%n
            if (dt%particle(iat)%ind(k) .eq. jat) then
                if (maxval(abs(dt%particle(iat)%v(:, k) - rij)) .lt. lo_tol) n = n + 1
            end if
        end do
    end do
    if (n .eq. fc%npairs) isok = .true.
end subroutine

!> Check the Hermitian sum rule
subroutine check_hermitian(fc, m, sumrule, peratom)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The tensor to which Hermiticity is tested
    real(r8), dimension(:, :, :), intent(in) :: m
    !> The sum of the Hermitian constraint on all atom
    real(r8), intent(out) :: sumrule
    !> The per-atom/xyz Hermitian sum rule, optional
    real(r8), dimension(:, :), optional, intent(out) :: peratom

    integer :: j, iat, jat

    sumrule = 0.0_r8
    if (present(peratom)) peratom = 0.0_r8
    do j=1, fc%npairs
        sumrule = sumrule + (m(1, 2, j) - m(2, 1, j)) / fc%na
        sumrule = sumrule + (m(1, 3, j) - m(3, 1, j)) / fc%na
        sumrule = sumrule + (m(2, 3, j) - m(3, 2, j)) / fc%na
        if (present(peratom)) then
            iat = fc%ij(1, j)
            jat = fc%ij(2, j)
            peratom(1, iat) = peratom(1, iat) + m(1, 2, j) - m(2, 1, j)
            peratom(2, iat) = peratom(2, iat) + m(1, 3, j) - m(3, 1, j)
            peratom(3, iat) = peratom(3, iat) + m(2, 3, j) - m(3, 2, j)
            peratom(1, jat) = peratom(1, jat) - m(1, 2, j) + m(2, 1, j)
            peratom(2, jat) = peratom(2, jat) - m(1, 3, j) + m(3, 1, j)
            peratom(3, jat) = peratom(3, jat) - m(2, 3, j) + m(3, 2, j)
        end if
    end do
end subroutine

!> Check the rotational sum rule
subroutine check_rotational(fc, m, uc, sumrule, peratom)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The tensor to which Hermiticity is tested
    real(r8), dimension(:, :, :), intent(in) :: m
    !> The structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The sum of the Hermitian constraint on all atom
    real(r8), intent(out) :: sumrule
    !> The per-atom/xyz Hermitian sum rule, optional
    real(r8), dimension(:, :), optional, intent(out) :: peratom

    integer :: j, iat, jat, a, b, c

    sumrule = 0.0_r8
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)

        do a=1, 3
        do b=1, 3
        do c=1, 3
            sumrule = sumrule + m(a, b, j) * uc%rcart(c, jat) - m(a, c, j) * uc%rcart(b, jat)
            sumrule = sumrule + m(b, a, j) * uc%rcart(c, iat) - m(b, c, j) * uc%rcart(a, iat)
        end do
        end do
        end do
    end do
end subroutine

!> Check the Huang sum rule
subroutine check_huang(fc, m, uc, sumrule, peratom)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The tensor to which Hermiticity is tested
    real(r8), dimension(:, :, :), intent(in) :: m
    !> The structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The sum of the Hermitian constraint on all atom
    real(r8), intent(out) :: sumrule
    !> The per-atom/xyz Hermitian sum rule, optional
    real(r8), dimension(:, :), optional, intent(out) :: peratom

    !> Buffer for the forceconstants
    real(r8), dimension(3, 3) :: m0
    integer :: j, iat, jat, a, b, c, d

    sumrule = 0.0_r8
    if (present(peratom)) peratom = 0.0_r8
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)

        m0 = fc%m(:, :, j)
        do a=1, 3
        do b=1, 3
        do c=1, 3
        do d=1, 3
            sumrule = sumrule + m0(a, b) * fc%rij(c, j) * fc%rij(d, j) - &
                                m0(c, d) * fc%rij(a, j) * fc%rij(b, j)
            sumrule = sumrule + m0(b, a) * fc%rij(c, j) * fc%rij(d, j) - &
                                m0(d, c) * fc%rij(a, j) * fc%rij(b, j)
        end do
        end do
        end do
        end do
    end do
end subroutine

!> Apply the Hermitian constraints on a IFC tensor
subroutine apply_hermitian_constraints(fc, m, mem)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The tensor on which to apply the constraint
    real(r8), dimension(:, :, :), intent(inout) :: m
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> A buffer for the result
    real(r8), dimension(:, :, :), allocatable :: m0
    !> A flat version of the input
    real(r8), dimension(:), allocatable :: mflat
    !> The c * c^T matrix
    real(r8), dimension(3, 3) :: cct
    !> To solve for the Lagrange multiplier
    real(r8), dimension(3, 1) :: lambda
    !> Some integer
    integer :: iat, i, j, ic, a, b, k, ii, jj

    call mem%allocate(m0, [3, 3, fc%npairs], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    m0 = 0.0_r8
    ! We do it for all atom
    do iat=1, fc%na
        call mem%allocate(mflat, fc%atom(iat)%n*9, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        ! First we flatten the tensor
        mflat = 0.0_r8
        do j=1, fc%atom(iat)%n
            ic = fc%atom(iat)%ic(j)

            do a=1, 3
            do b=1, 3
                if (ic .lt. 0) then
                    ii = 9 * (j - 1) + 3 * (b - 1) + a
                    jj = 9 * (j - 1) + 3 * (a - 1) + b
                else
                    ii = 9 * (j - 1) + 3 * (a - 1) + b
                    jj = 9 * (j - 1) + 3 * (b - 1) + a
                end if
                mflat(ii) = m(a, b, abs(ic))
                mflat(jj) = m(b, a, abs(ic))
            end do
            end do
        end do
        ! Then we solve the Lagrange multiplier
        cct = fc%atom(iat)%cct
        lambda(:, 1) = matmul(fc%atom(iat)%cstr, mflat)
        call lo_dposv(cct, lambda, info=i)
        if (i .ne. 0) call lo_stop_gracefully(['dposv exit status '//tochar(i)], lo_exitcode_blaslapack, __FILE__, __LINE__)

        ! We apply the Lagrange multiplier to enforce the constraints
        mflat = mflat - matmul(transpose(fc%atom(iat)%cstr), lambda(:, 1))

        ! And we distribute
        do j=1, fc%atom(iat)%n
            ic = fc%atom(iat)%ic(j)
            do a=1, 3
            do b=1, 3
                if (ic .lt. 0) then
                    ii = 9 * (j - 1) + 3 * (b - 1) + a
                    jj = 9 * (j - 1) + 3 * (a - 1) + b
                else
                    ii = 9 * (j - 1) + 3 * (a - 1) + b
                    jj = 9 * (j - 1) + 3 * (b - 1) + a
                end if
                m0(a, b, abs(ic)) = m0(a, b, abs(ic)) + 0.25_r8 * (mflat(ii) + mflat(jj))
            end do
            end do
        end do
        call mem%deallocate(mflat, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end do
    m = m0
    call mem%deallocate(m0, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

!> Create the tensor to enforce hermitian constraint
subroutine create_hermitian_constraints(fc)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(inout) :: fc

    !> Some integers for the do loops
    integer :: iat, j, k, a, b, ic, ii, jj

    do iat=1, fc%na
        allocate(fc%atom(iat)%cstr(3, fc%atom(iat)%n * 9))
        allocate(fc%atom(iat)%cct(3, 3))
        fc%atom(iat)%cstr = 0.0_r8
        do j=1, fc%atom(iat)%n
            ic = fc%atom(iat)%ic(j)
            k = 0
            do a=1, 3
            do b=a, 3
                if (a .eq. b) cycle
                k = k + 1
                if (ic .lt. 0) then
                    ii = 9 * (j - 1) + 3 * (b - 1) + a
                    jj = 9 * (j - 1) + 3 * (a - 1) + b
                else
                    ii = 9 * (j - 1) + 3 * (a - 1) + b
                    jj = 9 * (j - 1) + 3 * (b - 1) + a
                end if
                fc%atom(iat)%cstr(k, ii) = fc%atom(iat)%cstr(k, ii) + 1
                fc%atom(iat)%cstr(k, jj) = fc%atom(iat)%cstr(k, jj) - 1
            end do
            end do
        end do
        fc%atom(iat)%cct = matmul(fc%atom(iat)%cstr, transpose(fc%atom(iat)%cstr))
    end do
end subroutine

!> Compute the potential energy given displacements.
function potential_energy(fc, u, f) result(energy)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> The forces
    real(r8), dimension(:, :), optional, intent(in) :: f
    !> The energy
    real(r8) :: energy

    !> The force constant between iat and jat
    real(r8), dimension(3, 3) :: m
    !> Integers for do loop
    integer :: iat, jat, j, idx, a, b

    ! Little sanity check
    if (size(u, 2) .ne. fc%na) then
        write(*, *) 'Displacements and force constants do not seem to match'
        stop
    end if
    if (present(f)) then
        if (size(f, 2) .ne. fc%na) then
            write(*, *) 'Force and force constants do not seem to match'
            stop
        end if
    end if

    energy = 0.0_r8
    if (present(f)) then
        do iat=1, fc%na
            energy = energy - 0.5_r8 * dot_product(f(:, iat), u(:, iat))
        end do
    else
        do j=1, fc%npairs
            iat = fc%ij(1, j)
            jat = fc%ij(2, j)
            m = fc%m(:, :, j)

            do a=1, 3
            do b=1, 3
                energy = energy + 0.5_r8 * m(a, b) * u(a, iat) * (u(b, jat) - u(b, iat))
                energy = energy + 0.5_r8 * m(b, a) * u(b, jat) * (u(a, iat) - u(a, jat))
            end do
            end do
        end do
    end if
end function

!> Compute the forces given displacements
subroutine forces(fc, u, f)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> The displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> The forces
    real(r8), dimension(:, :), intent(out) :: f

    !> Integers for do loop
    integer :: iat, jat, j, idx

    integer :: a, b

    ! Little sanity check
    if (size(u, 2) .ne. fc%na) then
        write(*, *) 'Displacements and force constants do not seem to match'
        stop
    end if
    if (size(f, 2) .ne. fc%na) then
        write(*, *) 'Forces and force constants do not seem to match'
        stop
    end if

    f = 0.0_r8
    do j=1, fc%npairs
        iat = fc%ij(1, j)
        jat = fc%ij(2, j)
        ! Here it's the normal order
        f(:, iat) = f(:, iat) - matmul(fc%m(:, :, j), u(:, jat) - u(:, iat))
        ! Here it's the transpose version
        f(:, jat) = f(:, jat) - matmul(transpose(fc%m(:, :, j)), u(:, iat) - u(:, jat))
    end do
end subroutine

!> Compute the elastic constants
subroutine get_elastic_constants(fc, uc, mw, verbosity)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(inout) :: fc
    !> The structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Talk a lot ?
    integer, intent(in) :: verbosity

    !> Buffer for the elastic constants
    real(r8), dimension(3, 3, 3, 3) :: elc, bracket
    !> Buffer for the elastic constants in Voigt notation
    real(r8), dimension(6, 6) :: elc2
    !> The pair IFC matrix
    real(r8), dimension(3, 3) :: m
    !> The distances
    real(r8), dimension(3) :: rv
    !> Some integers
    integer :: iat, ipair, ic, i, j, al, be, gm, la

    bracket = 0.0_r8
    do iat=1, fc%na
        do ipair=1, fc%atom(iat)%n
            ic = fc%atom(iat)%ic(ipair)
            if (ic .lt. 0) then
                m = transpose(fc%m(:, :, abs(ic)))
                rv = -fc%rij(:, abs(ic))
            else
                m = fc%m(:, :, ic)
                rv = fc%rij(:, ic)
            end if
            do la=1, 3
            do gm=1, 3
            do be=1, 3
            do al=1, 3
                bracket(al, be, gm, la) = bracket(al, be, gm, la) + m(al, be) * rv(gm) * rv(la)
            end do
            end do
            end do
            end do
        end do
    end do
    bracket = -0.5_r8 * bracket / uc%volume

    do la=1, 3
    do gm=1, 3
    do be=1, 3
    do al=1, 3
        elc(al, be, gm, la) = bracket(al, gm, be, la) + bracket(be, gm, al, la) - bracket(al, be, gm, la)
    end do
    end do
    end do
    end do

    elc2 = 0.0_r8
    do la = 1, 3
    do gm = 1, 3
    do be = 1, 3
    do al = 1, 3
        i = contract_elastic_constant_indices(al, be)
        j = contract_elastic_constant_indices(gm, la)
        elc2(i, j) = elc(al, be, gm, la)
    end do
    end do
    end do
    end do

    ! Little symmetrization
    elc2 = 0.5_r8 * (elc2 + transpose(elc2))

    fc%elastic_constants_voigt = elc2
    fc%elastic_constants_tensor = elc
contains
    !> Voigt-notation
    function contract_elastic_constant_indices(mu, nu) result(ind)
        integer, intent(in) :: mu, nu
        integer, dimension(2) :: d
        integer :: ind

        ind = 0
        if (mu < nu) then
            d = [mu, nu]
        else
            d = [nu, mu]
        end if
        if (d(1) .eq. d(2)) ind = d(1)
        if (d(1) .eq. 1 .and. d(2) .eq. 2) ind = 6
        if (d(1) .eq. 1 .and. d(2) .eq. 3) ind = 5
        if (d(1) .eq. 2 .and. d(2) .eq. 3) ind = 4
    end function
end subroutine

!> Write forceconstants to HDF5
subroutine write_to_hdf5(fc, filename)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(in) :: fc
    !> Filename
    character(len=*), intent(in) :: filename

    !> HDF5 helper
    type(lo_hdf5_helper) :: h5

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    ! Get some attribute
    call h5%store_attribute(fc%na, h5%file_id, 'natoms')
    call h5%store_attribute(fc%cutoff, h5%file_id, 'cutoff')
    call h5%store_attribute(fc%mindist, h5%file_id, 'mindist')
    ! Store the coefficients
    call h5%store_data(fc%m, h5%file_id, 'coefficients', enhet='Ha/Bohr^2', dimensions='i,xyz,xyz')
    call h5%store_data(fc%ij, h5%file_id, 'ij', dimensions='i,iat-jat')
    call h5%store_data(fc%rij, h5%file_id, 'rij', enhet='Bohr', dimensions='i,xyz')

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> Read forceconstants from HDF5
subroutine read_from_hdf5(fc, filename, mem)
    !> The force constants
    class(lo_forceconstant_disorder_secondorder), intent(out) :: fc
    !> The name of the file
    character(len=*), intent(in) :: filename
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> HDF5 helper
    type(lo_hdf5_helper) :: h5

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('read', trim(filename))

    ! Read the attribute
    call h5%read_attribute(fc%na, h5%file_id, 'natoms')
    call h5%read_attribute(fc%cutoff, h5%file_id, 'cutoff')
    call h5%read_attribute(fc%mindist, h5%file_id, 'mindist')
    ! Read the data
    call h5%read_data(fc%m, h5%file_id, 'coefficients')
    call h5%read_data(fc%ij, h5%file_id, 'ij')
    call h5%read_data(fc%rij, h5%file_id, 'rij')
    fc%npairs = size(fc%m, 3)

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine
end module
