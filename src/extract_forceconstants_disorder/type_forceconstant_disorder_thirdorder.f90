module type_forceconstant_disorder_thirdorder
!!
!! Container for the third order forceconstants for disordered systems
!!
use konstanter, only: r8, lo_exitcode_param, lo_exitcode_baddim, lo_hugeint, lo_huge
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use hdf5_wrappers, only: lo_hdf5_helper
use type_crystalstructure, only: lo_crystalstructure
use type_distancetable, only: lo_distancetable


implicit none

private
public :: lo_forceconstant_disorder_thirdorder

type lo_forceconstant_disorder_thirdorder
    !> The number of atoms
    integer :: na = -lo_hugeint
    !> The total number of principal values
    integer :: nval = -lo_hugeint
    !> The principal vectors
    real(r8), dimension(:, :, :), allocatable :: pvec
    !> The cutoff
    real(r8) :: cutoff = -lo_huge
    !> The number of triplets
    integer :: ntriplets
    !> The triplets
    integer, dimension(:, :), allocatable :: ijk
    !> The distances between i and j
    real(r8), dimension(:, :), allocatable :: rij
    !> The distances between i and k
    real(r8), dimension(:, :), allocatable :: rik
    !> The distances between j and k
    real(r8), dimension(:, :), allocatable :: rjk
    !> The coefficients
    real(r8), dimension(:, :, :, :), allocatable :: m

contains
    !> Get the triplets given a structure
    procedure :: get_triplets
    !> The the triplets given a structure
    procedure :: generate
    !> Compute the potential energy
    procedure :: potential_energy
    !> Compute the forces
    procedure :: forces
    !> Write the forceconstants to hdf5
    procedure :: write_to_hdf5
    !> And read the forceconstants from hdf5
    procedure :: read_from_hdf5
end type

contains
!> Get all the tripletss for this system
subroutine get_triplets(fc, uc, cutoff, mw)
    !> The force constants
    class(lo_forceconstant_disorder_thirdorder), intent(inout) :: fc
    !> The crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The cutoff
    real(r8), intent(in) :: cutoff
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw

    !> The distance table
    type(lo_distancetable) :: dt
    !> The distances
    real(r8), dimension(3) :: rij, rik, rjk
    !> Some buffer
    real(r8) :: csqr
    !> Some integers
    integer :: iat, jat, kat, j, k, ii, jj
    !> For a sanity test
    logical :: ismapped

    fc%cutoff = cutoff
    fc%na = uc%na
    csqr = fc%cutoff**2

    !> First, we have to create a distance table
    call dt%generate(uc%r, uc%latticevectors, fc%cutoff, verbosity=0)

    ii = 0
    do iat=1, fc%na
        do j=1, dt%particle(iat)%n
        do k=j, dt%particle(iat)%n
            jat = dt%particle(iat)%ind(j)
            kat = dt%particle(iat)%ind(k)

            ! We don't need the self-terms
            if (jat .eq. iat .or. kat .eq. iat) cycle

            rij = dt%particle(iat)%v(:, j)
            rik = dt%particle(iat)%v(:, k)
            rjk = rik - rij
            if (sum(rjk**2) .lt. csqr) ii = ii + 1
        end do
        end do
    end do

    !> We can do some allocation
    fc%ntriplets = ii
    allocate(fc%ijk(3, fc%ntriplets))
    allocate(fc%rij(3, fc%ntriplets))
    allocate(fc%rik(3, fc%ntriplets))
    allocate(fc%rjk(3, fc%ntriplets))
    allocate(fc%m(3, 3, 3, fc%ntriplets))

    ! And we fill the tables
    ii = 0
    do iat=1, fc%na
        do j=1, dt%particle(iat)%n
        do k=j, dt%particle(iat)%n
            jat = dt%particle(iat)%ind(j)
            kat = dt%particle(iat)%ind(k)

            ! We don't need the self-term
            if (jat .eq. iat .or. kat .eq. iat) cycle

            rij = dt%particle(iat)%v(:, j)
            rik = dt%particle(iat)%v(:, k)
            rjk = rik - rij
            if (sum(rjk**2) .lt. csqr) then
                ii = ii + 1
                ! Get the atom involved in this triplet
                fc%ijk(1, ii) = iat
                fc%ijk(2, ii) = jat
                fc%ijk(3, ii) = kat

                ! And the distances
                fc%rij(:, ii) = rij
                fc%rik(:, ii) = rik
                fc%rjk(:, ii) = rjk
            end if
        end do
        end do
    end do
end subroutine

subroutine generate(fc, na, nval, pvec)
    !> The force constants
    class(lo_forceconstant_disorder_thirdorder), intent(inout) :: fc
    !> The number of atoms
    integer, intent(in) :: na
    !> The number of principal values
    integer, intent(in) :: nval
    !> The principal value eigenvectors
    real(r8), dimension(:, :, :), optional, intent(in) :: pvec

    fc%na = na
    fc%nval = nval
    allocate(fc%pvec(fc%nval, fc%na, 3))
    if (present(pvec)) then
        ! Some sanity test
        if (size(pvec, 1) .ne. nval .or. size(pvec, 2) .ne. na) then
            call lo_stop_gracefully(['Principal vectors do not match input'], &
                                     lo_exitcode_baddim, __FILE__, __LINE__)
        end if
        fc%pvec = pvec
    else
        fc%pvec = 0.0_r8
    end if
end subroutine

!> Potential energy
function potential_energy(fc, u, f) result(energy)
    !> The force constants
    class(lo_forceconstant_disorder_thirdorder), intent(in) :: fc
    !> The displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> The forces
    real(r8), dimension(:, :), optional, intent(in) :: f
    !> The energy
    real(r8) :: energy

    !> Buffer for the force constants
    real(r8), dimension(3, 3, 3) :: m
    !> Some integer
    integer :: iat, jat, kat, a, b, c, j

    ! Little sanity check
    if (fc%na .ne. size(u, 2)) then
        write(*, *) 'Displacements and force constants do not seem to match'
        stop
    end if
    if (present(f)) then
        if (fc%na .ne. size(f, 2)) then
            write(*, *) 'Force and force constants do not seem to match'
            stop
        end if
    end if

    energy = 0.0_r8
    if (present(f)) then
        do iat=1, fc%na
            energy = energy - dot_product(f(:, iat), u(:, iat)) / 3.0_r8
        end do
    else
!       do j=1, fc%ntriplets
!           iat = fc%ijk(1, j)
!           jat = fc%ijk(2, j)
!           kat = fc%ijk(3, j)

!           do a=1, 3
!           do b=1, 3
!           do c=1, 3
!               m = sum(fc%pvec(:, iat, a) * fc%pvec(:, jat, b) * fc%pvec(:, kat, c))
!               energy = energy + m * u(a, iat) * (u(b, jat) - u(b, iat)) * (u(c, kat) - u(c, iat)) / 6.0_r8
!               energy = energy + m * u(b, jat) * (u(a, iat) - u(a, jat)) * (u(c, kat) - u(c, jat)) / 6.0_r8
!               energy = energy + m * u(c, kat) * (u(a, iat) - u(a, kat)) * (u(b, jat) - u(b, kat)) / 6.0_r8
!           end do
!           end do
!           end do
!       end do
        do j=1, fc%ntriplets
            iat = fc%ijk(1, j)
            jat = fc%ijk(2, j)
            kat = fc%ijk(3, j)

            m = fc%m(:, :, :, j)
            do a=1, 3
            do b=1, 3
            do c=1, 3
                ! No 1/2, it takes into account the permutations j -> k
                energy = energy + m(a, b, c) * u(a, iat) * (u(b, jat) - u(b, iat)) * (u(c, kat) - u(c, iat)) / 6.0_r8
                energy = energy + m(a, b, c) * u(b, jat) * (u(a, iat) - u(a, jat)) * (u(c, kat) - u(c, jat)) / 6.0_r8
                energy = energy + m(a, b, c) * u(c, kat) * (u(a, iat) - u(a, kat)) * (u(b, jat) - u(b, kat)) / 6.0_r8
            end do
            end do
            end do
        end do
    end if

end function

!> Compute the forces given displacements
subroutine forces(fc, u, f)
    !> The force constants
    class(lo_forceconstant_disorder_thirdorder), intent(in) :: fc
    !> The displacements
    real(r8), dimension(:, :), intent(in) :: u
    !> The forces
    real(r8), dimension(:, :), intent(out) :: f

    !> Some buffer
    real(r8), dimension(3, 3) :: m1, m2, m3
    !> Buffer for the IFC tensor
    real(r8), dimension(3, 3, 3) :: m
    !> Integers for the do loops
    integer :: iat, jat, kat, j, a, b, c

    real(r8) :: f0

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
!   do j=1, fc%ntriplets
!       iat = fc%ijk(1, j)
!       jat = fc%ijk(2, j)
!       kat = fc%ijk(3, j)

!       do a=1, 3
!       do b=1, 3
!       do c=1, 3
!           m = sum(fc%pvec(:, iat, a) * fc%pvec(:, jat, b) * fc%pvec(:, kat, c))
!           ! No 1/2, it takes into account the permutations j -> k
!           f(a, iat) = f(a, iat) - m * (u(b, jat) - u(b, iat)) * (u(c, kat) - u(c, iat))
!           f(b, jat) = f(b, jat) - m * (u(a, iat) - u(a, jat)) * (u(c, kat) - u(c, jat))
!           f(c, kat) = f(c, kat) - m * (u(a, iat) - u(a, kat)) * (u(b, jat) - u(b, kat))
!       end do
!       end do
!       end do
!   end do
    do j=1, fc%ntriplets
        iat = fc%ijk(1, j)
        jat = fc%ijk(2, j)
        kat = fc%ijk(3, j)

        m = fc%m(:, :, :, j)
        do a=1, 3
        do b=1, 3
        do c=1, 3
            ! No 1/2, it takes into account the permutations j -> k
!           f(a, iat) = f(a, iat) - m(a, b, c) * (u(b, jat) - u(b, iat)) * (u(c, kat) - u(c, iat))
!           f(a, jat) = f(a, jat) - m(a, b, c) * (u(b, iat) - u(b, jat)) * (u(c, kat) - u(c, jat))
!           f(a, kat) = f(a, kat) - m(a, b, c) * (u(b, iat) - u(b, kat)) * (u(c, jat) - u(c, kat))

            f0 = sum(fc%pvec(:, iat, a) * fc%pvec(:, jat, b) * fc%pvec(:, kat, c))
            f(a, iat) = f(a, iat) - f0 * (u(b, jat) - u(b, iat)) * (u(c, kat) - u(c, iat))
            f(b, jat) = f(b, jat) - f0 * (u(a, iat) - u(a, jat)) * (u(c, kat) - u(c, jat))
            f(c, kat) = f(c, kat) - f0 * (u(a, iat) - u(c, kat)) * (u(b, jat) - u(b, kat))
        end do
        end do
        end do
    end do
end subroutine

!> Write the third order forceconstants to HDF5
subroutine write_to_hdf5(fc, filename)
    !> The force constants
    class(lo_forceconstant_disorder_thirdorder), intent(in) :: fc
    !> Filename
    character(len=*), intent(in) :: filename

    !> HDF5 helper
    type(lo_hdf5_helper) :: h5

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('write', trim(filename))

    call h5%store_attribute(fc%na, h5%file_id, 'natom')
    call h5%store_attribute(fc%nval, h5%file_id, 'nval')
    ! Store the coefficients
    call h5%store_data(fc%pvec, h5%file_id, 'pvec')

    call h5%close_file()
    call h5%destroy(__FILE__, __LINE__)
end subroutine

!> Read the third order forceconstants from HDF5
subroutine read_from_hdf5(fc, filename)
    !> The force constants
    class(lo_forceconstant_disorder_thirdorder), intent(out) :: fc
    !> The name of the file
    character(len=*), intent(in) :: filename

    !> HDF5 helper
    type(lo_hdf5_helper) :: h5

    call h5%init(__FILE__, __LINE__)
    call h5%open_file('read', trim(filename))

    call h5%read_attribute(fc%na, h5%file_id, 'natom')
    call h5%read_attribute(fc%nval, h5%file_id, 'nval')
    ! Store the coefficients
    call h5%read_data(fc%pvec, h5%file_id, 'pvec')

    call h5%close_file()
end subroutine
end module
