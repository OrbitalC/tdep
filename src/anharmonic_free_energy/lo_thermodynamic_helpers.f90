#include "precompilerdefinitions"
module lo_thermodynamic_helpers
use konstanter, only: r8, lo_huge
use type_crystalstructure, only: lo_crystalstructure
use type_symmetryoperation, only: lo_operate_on_secondorder_tensor

implicit none
private
public :: lo_thermodynamics
public :: lo_symmetrize_stress
public :: lo_full_to_voigt

!> A container for all thermodynamic results
type lo_thermodynamics
    !> The temperature at which everything is evaluated
    real(r8) :: temperature=-lo_huge
    !> The volume
    real(r8) :: volume=-lo_huge
    !> The free energies
    real(r8) :: f0, f3, f4
    !> The internal energy
    real(r8) :: u0, u3, u4
    !> The entropy
    real(r8) :: s0, s3, s4
    !> The heat capacity at constant volume
    real(r8) :: cv0, cv3, cv4
    !> The heat capacity at constant pressure
    real(r8) :: cp0, cp3, cp4
    !> The cumulant corrections, first dimension is IFC, second is cumulant order
    real(r8), dimension(3, 4) :: corr_fe=0.0_r8, corr_s=0.0_r8, corr_u=0.0_r8, corr_cv=0.0_r8
    !> The uncertainty for the cumulant corrections
    real(r8), dimension(3, 4) :: corr_fe_var=0.0_r8, corr_s_var=0.0_r8, corr_u_var=0.0_r8, corr_cv_var=0.0_r8
    !> The pressure
    real(r8) :: p=-lo_huge
    !> The stress tensor
    real(r8), dimension(3, 3) :: stress_pot, stress_kin, stress_potvar
    real(r8), dimension(3, 3) :: stress_3ph, stress_diff, stress_diffvar
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
