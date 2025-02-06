#include "precompilerdefinitions"
module thermo
use konstanter, only: r8, lo_huge

implicit none
private
public :: lo_thermodynamics

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
    !> The pressure
    real(r8) :: p=-lo_huge
    !> The stress tensor
    real(r8), dimension(3, 3) :: stress_pot, stress_kin, stress_potvar
    !> The thermal expansion
    real(r8) :: alpha
end type
end module
