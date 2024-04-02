#include "precompilerdefinitions"
module scattering
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol
use mpi_wrappers, only: lo_mpi_helper
use lo_memtracker, only: lo_mem_helper
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_randomnumbers, only: lo_mersennetwister


! Compute the scattering
subroutine compute_scattering_threephonon(qp, dr, fct, ratio3ph, mw, mem)
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Third order forceconstants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    ! The qpoint mesh
    class(lo_qpoint_mesh), allocatable :: qp
    ! The ratio of 3 phonon scattering actually computed
    real(r8), intent(in) :: ratio3ph
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    type(lo_mersennetwister) :: rng
    ! The number of scattering
    integer, dimension(:, :), allocatable :: nplus_tot, nplus_sample, nminus_tot, nminus_sample

    integer :: q2, b2, q3, b3, ctr

    call rng%init(iseed=mw%r, rseed=walltime())

    nplus_tot = 0
    nplus_sample = 0
    nminus_tot = 0
    nminus_sample = 0

    call mem%allocate(nplus_tot, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(nplus_sample, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(scatplus, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(nminus_tot, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(nminus_sample, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(scatminus, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    nplus_tot = 0
    nplus_sample = 0
    scatplus = 0.0_r8
    nminus_tot = 0
    nminus_sample = 0
    scatminus = 0.0_r8

    ctr = 0
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! MPI
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            om1 = dr%iq(q1)%omega
            if (om1 .lt. omthres) cycle
            do q2 = 1, qp%n_full_point
                do b2 = 1, dr%n_mode
                    do b3 = 1, dr%n_mode
                        q3 = fft_third_grid_index(q1, q2, dims)
                        om2 = dr%aq(q2)%omega
                        om3 = dr%aq(q3)%omega
                        vel2 = dr%aq(q2)%vel
                        vel3 = dr%aq(q3)%vel
                        if (om2 .lt. omthres) cycle
                        if (om3 .lt. omthres) cycle

                        select case (integrationtype)
                        case (1)
                            sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
                        case (2)
                            sig1 = qp%adaptive_sigma(qp%ap(q1)%radius, vel1(:, b1), dr%default_smearing, smearing_prefactor)
                            sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, vel2(:, b2), dr%default_smearing, smearing_prefactor)
                            sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, vel3(:, b3), dr%default_smearing, smearing_prefactor)
                            sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                        end select

                        if (abs(om1 + om2 - om3) .lt. thres*sigma) then
                            nplus_tot(q1, b1) = nplus_tot(q1, b1) + 1
                            rnd = rng%rnd_real()
                            if (rnd .lt. ratio3ph) then
                                cycle
                            else
                                nplus_sample(q1, b1) = nplus_sample(q1, b1) + 1
                                c0 = fcf%scattering_amplitude()
                                psisq = abs(c0 * conjg(c0))
                                deltafunction = lo_gauss(om1, -om2 + om3, sigma)
                                scatplus(q1, b1) = scatplus(q1, b1) + deltafunction * psisq * qp%ap(i)%integration_weight
                            end if
                        end if

                        if (abs(om1 - om2 - om3) .lt. thres*sigma) then
                            nminus_tot(q1, b1) = nminus_tot(q1, b1) + 1
                            rnd = rng%rnd_real()
                            if (rnd .lt. ratio3ph) then
                                cycle
                            else
                                nminus_sample(q1, b1) = nminus_sample(q1, b1) + 1
                                c0 = fcf%scattering_amplitude()
                                psisq = abs(c0 * conjg(c0))
                                deltafunction = lo_gauss(om1, om2 + om3, sigma)
                                scatminus(q1, b1) = scatminus(q1, b1) + 0.5_r8 * deltafunction * psisq * qp%ap(i)%integration_weight
                            end if
                        end if
                    end do ! b3
                end do ! b2
            end do ! q2
        end do ! b1
    end do ! q1

    ! Sum every rank
    call mw%allreduce('sum', nplus_tot)
    call mw%allreduce('sum', nplus_sample)
    call mw%allreduce('sum', scatplus)
    call mw%allreduce('sum', nminus_tot)
    call mw%allreduce('sum', nminus_sample)
    call mw%allreduce('sum', scatminus)

    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            dr%ip(q1)%p_plus = scatplus(q1, b1) * real(nplus_tot(q1, b1), r8) / real(nplus_sample(q1, b1), r8)
            dr%ip(q1)%p_minus = scatminus(q1, b1) * real(nminus_tot(q1, b1), r8) / real(nminus_sample(q1, b1), r8)
        end do
    end do
end subroutine


subroutine compute_scattering_isotopes(qp, dr, uc, mw, mem)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> MPI helper
    type(lo_mpi_helper), intent(in) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    integer ::

    call mem%allocate(iso_scatter, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    iso_scatter = 0.0_r8

    ctr = 0
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! MPI
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. omthres) cycle
            do q2 = 1, qp%n_full_point
            do b2 = 1, dr%n_mode
                om2 = dr%aq(q2)%omega(b2)
                if (om1 .lt. omthres) cycle
                select case
                case (1)
                    sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
                case (2)
                    sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2))
                end select

                deltafunction = lo_gauss(om1, om2, sigma)
                egviso(:, 1) = dr%iq(q1)%egv(:, b1)
                egviso(:, 2) = dr%aq(q1)%egv(:, b2)
                scatterstrength = isotope_scattering_strength(uc, egviso)
                iso_scatter(q1, b1) = scatterstrength * deltafunction * qp%ap(q2)%integration_weight
            end do
            end do
        end do
    end do

    call mw%allreduce('sum', iso_scatter)
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            dr%ip(q1)%p_iso = iso_scatter(q1, b1)
        end do
    end do
end subroutine

contains

real(r8) function isotope_scattering_strength(uc, egv)
    type(lo_crystalstructure), intent(in) :: uc
    complex(r8), dimension(:, :), intent(in) :: egv
    !
    integer :: i, j
    real(r8) :: f0, f1
    complex(r8), dimension(3) :: cv0, cv1
    !
    f1 = 0.0_r8
    do i = 1, uc%na
        cv0 = egv((i - 1)*3 + 1:(i*3), 1)
        cv1 = egv((i - 1)*3 + 1:(i*3), 2)
        f0 = 0.0_r8
        do j = 1, 3
            f0 = f0 + abs(conjg(cv0(j))*cv1(j))
        end do
        f0 = f0**2
        f1 = f1 + f0*uc%isotope(i)%disorderparameter
    end do
    isotope_scattering_strength = f1
    !
end function
end module
