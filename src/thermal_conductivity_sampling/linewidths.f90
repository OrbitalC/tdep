#include "precompilerdefinitions"
module linewidths
use konstanter, only: r8, lo_freqtol, lo_huge, lo_phonongroupveltol, lo_pi
use gottochblandat, only: walltime, lo_lorentz, lo_planck, lo_gauss
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
! local module
use options, only: lo_opts
use scattering, only: lo_scattering_rates

implicit none

private
public :: compute_linewidths

contains
subroutine compute_linewidths(qp, dr, sr, opts, mw, mem)
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The options
    type(lo_opts), intent(in) :: opts
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> Some buffer
    real(r8), dimension(:, :), allocatable :: buf_lw
    real(r8) :: maxdif, t0, buf, velnorm, f0
    real(r8) :: sigma, sig1, sig2, sig3, sig4, n1, n2, n3, n4, n2p, n3p, n4p, om1, om2, om3, om4, psisq
    !> Integers for the loops
    integer :: i, j, il, q1, q2, q3, q4, b1, b2, b3, b4, ii, jj, kk

    call mem%allocate(buf_lw, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf_lw = 0.0_r8
    do il=1, sr%nlocal_point
        buf = 0.0_r8
        q1 = sr%q1(il)
        b1 = sr%b1(il)

        om1 = dr%iq(q1)%omega(b1)
        sig1 = sr%sigma_q(q1, b1)
        n1 = sr%be(q1, b1)

        if (opts%isotopescattering) then
            do j = 1, sr%iso(il)%n
                q2 = sr%iso(il)%q2(j)
                b2 = sr%iso(il)%b2(j)
                om2 = dr%aq(q2)%omega(b2)

                n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
                sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)

        !       f0 = sr%iso(il)%psisq(j) * om1 * om2 * n1 * (n2 + 1.0_r8)
                f0 = sr%iso(il)%psisq(j) * om1 * om2
                f0 = f0 * lo_gauss(om1, om2, sigma)
                buf = buf + f0
            end do
        end if
        if (opts%thirdorder) then
            do j=1, sr%threephonon(il)%n
                q2 = sr%threephonon(il)%q2(j)
                q3 = sr%threephonon(il)%q3(j)
                b2 = sr%threephonon(il)%b2(j)
                b3 = sr%threephonon(il)%b3(j)

                psisq = sr%threephonon(il)%psisq(j)
                om2 = dr%aq(q2)%omega(b2)
                om3 = dr%aq(q3)%omega(b3)
                n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)
                n2p = n2 + 1.0_r8
                n3p = n3 + 1.0_r8

                sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
                ! sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                sigma = sqrt(sig2**2 + sig3**2)

                ! We have to take care of the permutation here
        !       buf = buf + psisq * n1 * n2 * n3p * lo_gauss(om1, -om2 + om3, sigma)
        !       buf = buf + psisq * n1 * n3 * n2p * lo_gauss(om1, -om3 + om2, sigma)
                buf = buf + 2.0_r8 * psisq * (n2 + n3 + 1.0_r8) * lo_gauss(om1,  om2 + om3, sigma)
                buf = buf - 2.0_r8 * psisq * (n2 + n3 + 1.0_r8) * lo_gauss(om1, -om2 - om3, sigma)
                ! There is a 2.0 prefactor here from permutation,
                ! But it's cancelled by the 0.5 from the process
        !       buf = buf + psisq * n1 * n2p * n3p * lo_gauss(om1, om2 + om3, sigma)
                buf = buf + 2.0_r8 * psisq * (n2 - n3) * lo_gauss(om1, -om2 + om3, sigma)
                buf = buf - 2.0_r8 * psisq * (n2 - n3) * lo_gauss(om1,  om2 - om3, sigma)
            end do
        end if
        if (opts%fourthorder) then
            do j=1, sr%fourphonon(il)%n
                q2 = sr%fourphonon(il)%q2(j)
                q3 = sr%fourphonon(il)%q3(j)
                q4 = sr%fourphonon(il)%q4(j)
                b2 = sr%fourphonon(il)%b2(j)
                b3 = sr%fourphonon(il)%b3(j)
                b4 = sr%fourphonon(il)%b4(j)

                psisq = sr%fourphonon(il)%psisq(j)
                om2 = dr%aq(q2)%omega(b2)
                om3 = dr%aq(q3)%omega(b3)
                om4 = dr%aq(q4)%omega(b4)
                n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)
                n4 = sr%be(qp%ap(q4)%irreducible_index, b4)
                n2p = n2 + 1.0_r8
                n3p = n3 + 1.0_r8
                n4p = n4 + 1.0_r8

                sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
                sig4 = sr%sigma_q(qp%ap(q4)%irreducible_index, b4)
                sigma = sqrt(sig2**2 + sig3**2 + sig4**2)

                ! This one is invariant with changes of 2<->3 -> factor 2.0_r8, cancels out
                buf = buf + psisq * n1 * n2 * n3 * n4p * lo_gauss(om1, -om2 - om3 + om4, sigma)
                ! The same
                buf = buf + psisq * n1 * n2 * n4 * n3p * lo_gauss(om1, -om2 - om4 + om3, sigma)
                ! The same
                buf = buf + psisq * n1 * n3 * n4 * n2p * lo_gauss(om1, -om3 - om4 + om2, sigma)

                ! This one is invariant with changes of 2<->3 -> factor 2.0_r8, cancels out
                buf = buf + psisq * n1 * n2 * n3p * n4p * lo_gauss(om1, -om2 + om3 + om4, sigma)
                ! The same
                buf = buf + psisq * n1 * n3 * n2p * n4p * lo_gauss(om1, -om3 + om2 + om4, sigma)
                ! The same
                buf = buf + psisq * n1 * n4 * n2p * n3p * lo_gauss(om1, -om4 + om2 + om3, sigma)

                ! Prefactor going away from cancellation from permutation
                buf = buf + psisq * n1 * n2p * n3p * n4p * lo_gauss(om1, om2 + om3 + om4, sigma)
            end do
        end if

        ! Let's add the boundary scattering
        if (opts%mfp_max .gt. 0.0_r8) then
            velnorm = norm2(dr%iq(q1)%vel(:, b1))
            ! buf = buf + n1 * (n1 + 1.0_r8) * velnorm / opts%mfp_max
            buf = buf + velnorm / opts%mfp_max
        end if
        ! buf_lw(q1, b1) = 0.5_r8 * buf / (n1 * (n1 + 1.0_r8))
        ! buf_lw(q1, b1) = buf / (n1 * (n1 + 1.0_r8))
          buf_lw(q1, b1) = buf
    end do
    call mw%allreduce('sum', buf_lw)
    ! Now we distribute the linewidths
    distribute: block
        real(r8) :: buf
        integer :: j, q1, b1, b2
        do q1=1, qp%n_irr_point
            do b1=1, dr%n_mode
                if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
                ! First we fix the degeneracy
                buf = 0.0_r8
                do j=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    buf = buf + buf_lw(q1, b2)
                end do
                buf = buf / real(dr%iq(q1)%degeneracy(b1), r8)
                do j=1, dr%iq(q1)%degeneracy(b1)
                    b2 = dr%iq(q1)%degenmode(j, b1)
                    buf_lw(q1, b2) = buf
                end do
                dr%iq(q1)%linewidth(b1) = buf_lw(q1, b1)
            end do
        end do
    end block distribute
    call mem%deallocate(buf_lw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine
end module
