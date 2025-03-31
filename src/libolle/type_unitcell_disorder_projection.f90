module type_unitcell_disorder_projection
!!
!! Container for the mapping between a disordered supercell and a pristine unitcell
!!
use konstanter, only: r8, lo_exitcode_baddim, lo_hugeint, lo_sqtol, lo_sqrectol, &
                      lo_huge, lo_degenvector, lo_rectol
use gottochblandat, only: open_file, lo_invert3x3matrix, lo_chop, lo_determ, lo_sqnorm, &
                          lo_return_unique, lo_clean_fractional_coordinates, lo_reciprocal_basis
use geometryfunctions, only: lo_bounding_sphere_of_box, lo_inscribed_sphere_in_box
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint
use type_crystalstructure, only: lo_crystalstructure

type lo_unitcell_disorder_projection
    !> The number of atoms in the unitcell
    integer :: na_uc=-lo_hugeint
    !> The number of atoms in the supercell
    integer :: na_ss=-lo_hugeint
    !> The number of commensurate q-points
    integer :: nqpt
    !> The list of projection
    integer, dimension(:), allocatable :: unitcell_index
    !> The supercell matrix from unitcell to supercell
    integer, dimension(3, 3) :: supercellmatrix=-lo_hugeint
    !> The unitcell on which to project
    type(lo_crystalstructure) :: uc
    !> The commensurate q-points
    type(lo_qpoint), dimension(:), allocatable :: cq
    !> The supercell lattice vectors
    real(r8), dimension(3, 3) :: supercell_latticevectors
    !> The supercell reciprocal lattice vectors
    real(r8), dimension(3, 3) :: supercell_reciprocal_latticevectors
    !> The q_distance between two adjacent point
    real(r8), dimension(3) :: qdist = -lo_huge

contains
    !> Read from file
    procedure :: read_from_file
    !> Generate the commensurate q-points
    procedure :: commensurate_qpoints
end type

contains

!> Read from file
subroutine read_from_file(proj, uc, ss, filename, mem)
    !> The projection
    class(lo_unitcell_disorder_projection), intent(out) :: proj
    !> The unitcell
    type(lo_crystalstructure), intent(in) :: uc
    !> The supercell
    type(lo_crystalstructure), intent(in) :: ss
    !> The filename
    character(len=*), intent(in) :: filename
    !> The memory helper
    type(lo_mem_helper), intent(inout) :: mem

    real(r8), dimension(3, 3) :: ssm
    !> Some integers
    integer :: i, i_uc, u

    integer, dimension(3) :: ssi

    proj%uc = uc

    u=open_file('in', trim(filename))
        read(u, *) proj%na_uc, proj%na_ss
        if (proj%na_uc .ne. uc%na) then
            call lo_stop_gracefully(['The number of atoms in the unitcell and &
                                     &the infile.unitcell_projection do not match'], &
                                    lo_exitcode_baddim, __FILE__, __LINE__)
        end if
        if (proj%na_ss .ne. ss%na) then
            call lo_stop_gracefully(['The number of atoms in the supercell and &
                                     &the infile.unitcell_projection do not match'], &
                                    lo_exitcode_baddim, __FILE__, __LINE__)
        end if
        allocate(proj%unitcell_index(proj%na_ss))
        ! And read the infile
        do i=1, proj%na_ss
            read(u, *) i_uc
            proj%unitcell_index(i) = i_uc
        end do
    close(u)
    ! Some sanity checks
    if (maxval(proj%unitcell_index) .gt. uc%na .or. minval(proj%unitcell_index) .le. 0) then
        call lo_stop_gracefully(['Trying to map to atoms not on the unitcell'], &
                                    lo_exitcode_baddim, __FILE__, __LINE__)
    end if

    ! Now we get the supercell matrix
    ssm = matmul(uc%inv_latticevectors, ss%latticevectors)
    proj%supercellmatrix = int(lo_chop(anint(ssm), lo_sqtol))
    proj%supercell_latticevectors = ss%latticevectors
    proj%supercell_reciprocal_latticevectors = ss%reciprocal_latticevectors
end subroutine

!> Generate the commensurate q-points
subroutine commensurate_qpoints(proj)
    !> The projection
    class(lo_unitcell_disorder_projection), intent(inout) :: proj

    generate_qpoints: block
        !> Some 3x3 matrix
        real(r8), dimension(3, 3) :: m0, reclat
        !> Some q-points
        real(r8), dimension(3) :: q
        !> Some integers
        integer :: i, iq, irep, nrep, a, b, c

        ! First we get a number of repetition to do
        reclat = proj%supercell_reciprocal_latticevectors
        f0 = lo_bounding_sphere_of_box(proj%uc%reciprocal_latticevectors) * 2
        nrep = 0
        do irep=1, 1000
            do i=1, 3
                m0(:, i) = (2 * real(nrep, r8) + 1) * reclat(:, i)
            end do
            f1 = lo_inscribed_sphere_in_box(m0) - 1e-5_r8
            if (f1 .gt. f0) then
                exit
            else
                nrep = nrep + 1
            end if
        end do

        ! Now we count the number of point in the FBZ of the unitcell
        ! Some located on the edge will be equivalent, but we will keep them for interpolate later on
        iq = 0
        do a=-nrep, nrep
        do b=-nrep, nrep
        do c=-nrep, nrep
            q = matmul(reclat, real([a, b, c], r8))
            ! A q-point is in the FBZ if it's the required shift is 0
            if (lo_sqnorm(proj%uc%bz%gshift(q)) .lt. lo_sqrectol) iq = iq + 1
        end do
        end do
        end do

        ! Now we can allocate the q-points
        proj%nqpt = iq
        allocate(proj%cq(iq))

        ! But we still need to check that we have at least more q-point than what we should have
        if (proj%nqpt .lt. abs(lo_determ(proj%supercellmatrix))) then
            call lo_stop_gracefully(['Could not map all commensurate q-points of the supercell'], &
                                    lo_exitcode_baddim, __FILE__, __LINE__)
        end if

        iq = 0
        do a=-nrep, nrep
        do b=-nrep, nrep
        do c=-nrep, nrep
            q = matmul(reclat, real([a, b, c], r8))
            if (lo_sqnorm(proj%uc%bz%gshift(q)) .lt. lo_sqrectol) then
                iq = iq + 1
                proj%cq(iq)%r = lo_chop(q, lo_rectol)
            end if
        end do
        end do
        end do

        !> We also need the spacing of q-point
        proj%qdist = matmul(reclat, [1.0_r8, 1.0_r8, 1.0_r8])
    end block generate_qpoints
end subroutine
end module
