#include "precompilerdefinitions"
module file_io
use konstanter, only: r8, lo_frequency_Hartree_to_Hz,  lo_exitcode_param, lo_tol
use type_phonon_dispersions, only: lo_phonon_dispersions
use hdf5_wrappers, only: lo_hdf5_helper
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint_mesh


implicit none
private
public :: read_linewidths


contains
subroutine read_linewidths(dr, qp, filename, mw, mem)
    !> The phonon dispersion
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The q-point mesh
    type(lo_qpoint_mesh), intent(in) :: qp
    !> The filename
    character(len=*), intent(in) :: filename
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_hdf5_helper) :: h5
    !> The linewidth
    real(r8), dimension(:, :), allocatable :: lw
    !> The qpoints
    real(r8), dimension(:, :), allocatable :: buf_qp
    !> Some integers
    integer :: readrnk, n_mode, n_full_point, q1, a, qi1, b1
    integer, dimension(2) :: buf_shape


    if (mw%talk) write(*, "(1X,A)") '... reading linewidth'

    readrnk = mw%n - 1
    n_mode = -1
    n_full_point = -1

    if (mw%r .eq. readrnk) then
        call h5%init()
        call h5%open_file('read', trim(filename))

        call h5%read_data(lw, h5%file_id, 'linewidths')
        call h5%read_data(buf_qp, h5%file_id, 'qpoints')

        call h5%close_file()
        call h5%destroy()

        buf_shape = shape(lw)
        n_full_point = buf_shape(2)
        n_mode = buf_shape(1)

        ! Little sanity check to make sure everything agree
        if (n_full_point .ne. qp%n_full_point) then
            call lo_stop_gracefully(['Mismatch in number of q-points'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if
        if (n_mode .ne. dr%n_mode) then
            call lo_stop_gracefully(['Mismatch in number of modes'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
        end if

        do q1=1, qp%n_full_point
            do a=1, 3
                if (abs(buf_qp(a, q1) - qp%ap(q1)%r(a)) .gt. lo_tol) then
                    call lo_stop_gracefully(['Mismatch in value of q-point'], lo_exitcode_param, __FILE__, __LINE__, mw%comm)
                end if
            end do
        end do
    end if

    ! And now we broadcast
    if (mw%r .ne. readrnk) then
        allocate(lw(dr%n_mode, qp%n_full_point))
    end if
    call mw%bcast(lw, from=readrnk)

    do q1 = 1, qp%n_irr_point
        allocate (dr%iq(q1)%linewidth(dr%n_mode))
        dr%iq(q1)%linewidth = 0.0_r8
    end do
    ! Now we can populate the linewidths
    do q1=1, qp%n_full_point
        do b1=1, dr%n_mode
            qi1 = qp%ap(q1)%irreducible_index
            dr%iq(qi1)%linewidth(b1) = lw(b1, q1) / lo_frequency_Hartree_to_Hz
        end do
    end do
end subroutine
end module
