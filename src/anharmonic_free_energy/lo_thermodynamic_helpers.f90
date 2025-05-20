#include "precompilerdefinitions"
module lo_thermodynamic_helpers
use konstanter, only: r8, lo_huge
use gottochblandat, only: lo_chop
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_operate_on_secondorder_tensor

implicit none
private
public :: lo_thermodynamics
public :: lo_symmetrize_stress
public :: lo_full_to_voigt

type lo_thermodynamic_contribution
    ! The different contributions. First dimension is the result, second one the uncertainty
    real(r8), dimension(2) :: F=0.0_r8
    real(r8), dimension(2) :: U=0.0_r8
    real(r8), dimension(2) :: S=0.0_r8
    real(r8), dimension(2) :: Cv=0.0_r8
    real(r8), dimension(3, 3, 2) :: stress=0.0_r8
end type

!> A container for all thermodynamic results
type lo_thermodynamics
    !> The temperature at which everything is evaluated
    real(r8) :: temperature=-lo_huge
    ! The things to hold results
    type(lo_thermodynamic_contribution) :: harmonic
    type(lo_thermodynamic_contribution) :: first_order
    type(lo_thermodynamic_contribution) :: second_order
    type(lo_thermodynamic_contribution) :: threephonon
    type(lo_thermodynamic_contribution) :: fourphonon
    !> The thermal expansion
    real(r8), dimension(3, 3) :: alpha
end type

contains
! Symmetrize 3x3 tensors
subroutine lo_symmetrize_stress(m, uc)
    !> The input matrix
    real(r8), dimension(3, 3), intent(inout) :: m
    !> The unitcell
    type(lo_crystalstructure), intent(in) :: uc

    real(r8), dimension(3, 3) :: tmp
    integer :: iop

    tmp = 0.0_r8
    do iop = 1, uc%sym%n
        tmp = tmp + lo_operate_on_secondorder_tensor(uc%sym%op(iop), m)
    end do
    m = tmp / real(uc%sym%n, r8)
    m = lo_chop(m, sum(abs(m))*1e-6_r8)
end subroutine

! Go from full 3x3 real space to 6 voigt notation
function lo_full_to_voigt(i, j) result(k)
    !> Inputs
    integer, intent(in) :: i
    integer, intent(in) :: j

    integer :: k

    integer, dimension(2) :: d
    k = 0
    if (i < j) then
        d = [i, j]
    else
        d = [j, i]
    end if
    if (d(1) .eq. d(2)) k = d(1)
    if (d(1) .eq. 1 .and. d(2) .eq. 2) k = 6
    if (d(1) .eq. 1 .and. d(2) .eq. 3) k = 5
    if (d(1) .eq. 2 .and. d(2) .eq. 3) k = 4
end function
end module
