#include "precompilerdefinitions"
module linewidths
use konstanter, only: r8, lo_freqtol, lo_huge, lo_phonongroupveltol
use gottochblandat, only: walltime, lo_lorentz, lo_planck, lo_gauss
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_qpointmesh, only: lo_qpoint_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
! local module
use options, only: lo_opts
use new_scattering, only: lo_scattering_rates

implicit none

private
public :: compute_linewidths
public :: self_consistent_linewidths

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
    real(r8) :: maxdif, t0, buf, n1, velnorm
    !> Integers for the loops
    integer :: i, j, il, q1, b1, b2

    call mem%allocate(buf_lw, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf_lw = 0.0_r8
    do il=1, sr%nlocal_point
        buf = 0.0_r8
        if (opts%isotopescattering) then
            do j = 1, sr%iso(il)%n
                buf = buf + sr%iso(il)%event(j)%W
            end do
        end if
        if (opts%thirdorder) then
            do j = 1, sr%threephonon(il)%nplus
                buf = buf + sr%threephonon(il)%plus(j)%W
            end do
            do j = 1, sr%threephonon(il)%nminus
                buf = buf + sr%threephonon(il)%minus(j)%W * 0.5_r8
            end do
        end if
        if (opts%fourthorder) then
            do j = 1, sr%fourphonon(il)%npp
                buf = buf + sr%fourphonon(il)%pp(j)%W * 0.5_r8
            end do
            do j = 1, sr%fourphonon(il)%npm
                buf = buf + sr%fourphonon(il)%pm(j)%W * 0.5_r8
            end do
            do j = 1, sr%fourphonon(il)%nmm
                buf = buf + sr%fourphonon(il)%mm(j)%W / 6.0_r8
            end do
        end if

        q1 = sr%q1(il)
        b1 = sr%b1(il)
        n1 = lo_planck(opts%temperature, dr%iq(q1)%omega(b1))
        ! Let's add the boundary scattering
        if (opts%mfp_max .gt. 0.0_r8) then
            velnorm = norm2(dr%iq(q1)%vel(:, b1))
            buf = buf + n1 * (n1 + 1.0_r8) * velnorm / opts%mfp_max
        end if
        buf_lw(q1, b1) = 0.5_r8 * buf / (n1 * (n1 + 1.0_r8))
    end do
    call distribute_linewidths(buf_lw, dr, qp, opts)
    call mem%deallocate(buf_lw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine


subroutine self_consistent_linewidths(dr, qp, sr, opts, mw, mem)
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
    real(r8), dimension(:, :), allocatable :: buf_lw, old_lw
    real(r8) :: maxdif, t0, buf, n1, velnorm
    !> Integers for the loops
    integer :: i, j, il, q1, b1, b2

    call mem%allocate(buf_lw, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(old_lw, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf_lw = 0.0_r8
    old_lw = 0.0_r8

    if (mw%talk) write(*, *) '... self consistent linewidths'
    if (mw%talk) write(*, "(1X,A32,A21,5X, A16)") '', 'max(DeltaGamma/Gamma)', 'elapsed time (s)'
    t0 = walltime()

    ! Let put some values for the linewidths
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            old_lw(q1, b1) = dr%iq(q1)%linewidth(b1)
        end do
    end do

    do i=1, opts%niter
        buf_lw = 0.0_r8
        ! Lets compute the scattering rates
        do il=1, sr%nlocal_point
            buf = 0.0_r8
            if (opts%isotopescattering) then
                call compute_scatteringrate_isotope(il, buf, old_lw, dr, qp, sr, opts%temperature, mw, mem)
            end if
            if (opts%thirdorder) then
                call compute_scatteringrate_threephonon(il, buf, old_lw, dr, qp, sr, opts%temperature, mw, mem)
            end if
            if (opts%fourthorder) then
                call compute_scatteringrate_fourphonon(il, buf, old_lw, dr, qp, sr, opts%temperature, mw, mem)
            end if

            q1 = sr%q1(il)
            b1 = sr%b1(il)
            n1 = lo_planck(opts%temperature, dr%iq(q1)%omega(b1))
            ! Let's add the boundary scattering
            if (opts%mfp_max .gt. 0.0_r8) then
                velnorm = norm2(dr%iq(q1)%vel(:, b1))
                buf = buf + n1 * (n1 + 1.0_r8) * velnorm / opts%mfp_max
            end if
            buf_lw(q1, b1) = 0.5_r8 * buf / (n1 * (n1 + 1.0_r8))
        end do
        call mw%allreduce('sum', buf_lw)

        ! Update linewidth and check convergence
        maxdif = -lo_huge
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

                ! Now we update the value
                buf_lw(q1, b1) = opts%mixing * buf_lw(q1, b1) + (1 - opts%mixing) * old_lw(q1, b1)
                maxdif = maxval([maxdif, abs(buf_lw(q1, b1) - old_lw(q1, b1)) / buf_lw(q1, b1)])
           !    maxdif = 0.0_r8
            end do
        end do
        old_lw = buf_lw

        if (mw%talk) write(*, "(11X,A10,I6,E21.3,5X,E14.3)") 'iteration ', i, maxdif, walltime() - t0
        if (maxdif .lt. opts%scftol) exit
    end do

    call distribute_linewidths(buf_lw, dr, qp, opts)
    call mem%deallocate(buf_lw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(old_lw, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine


subroutine distribute_linewidths(lw, dr, qp, opts)
    !> The linewidhts
    real(r8), dimension(:, :), intent(inout) :: lw
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! The options
    type(lo_opts), intent(in) :: opts

    ! Buffer
    real(r8) :: buf, n1, velnorm
    ! For the do loops
    integer :: q1, b1, b2, j

    ! Now we need to compute the Fn
    do q1=1, qp%n_irr_point
        do b1=1, dr%n_mode
            ! Skip gamma for acoustic branches
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle

            ! First we fix the degeneracy
            buf = 0.0_r8
            do j=1, dr%iq(q1)%degeneracy(b1)
                b2 = dr%iq(q1)%degenmode(j, b1)
                buf = buf + lw(q1, b2)
            end do
            buf = buf / real(dr%iq(q1)%degeneracy(b1), r8)
            do j=1, dr%iq(q1)%degeneracy(b1)
                b2 = dr%iq(q1)%degenmode(j, b1)
                lw(q1, b2) = buf
            end do

            n1 = lo_planck(opts%temperature, dr%iq(q1)%omega(b1))
            dr%iq(q1)%linewidth(b1) = lw(q1, b1)
            dr%iq(q1)%qs(b1) = 2.0_r8 * n1 * (n1 + 1.0_r8) * lw(q1, b1)

            velnorm = norm2(dr%iq(q1)%vel(:, b1))
            if (velnorm .gt. lo_phonongroupveltol) then
                dr%iq(q1)%mfp(:, b1) = dr%iq(q1)%vel(:, b1) * 0.5_r8 / dr%iq(q1)%linewidth(b1)
                dr%iq(q1)%scalar_mfp(b1) = velnorm * 0.5_r8 / dr%iq(q1)%linewidth(b1)
                dr%iq(q1)%F0(:, b1) = dr%iq(q1)%mfp(:, b1) * dr%iq(q1)%omega(b1) / opts%temperature
                dr%iq(q1)%Fn(:, b1) = dr%iq(q1)%F0(:, b1)
            end if
        end do
    end do
end subroutine


subroutine compute_scatteringrate_isotope(il, buf, lw, dr, qp, sr, temperature, mw, mem)
    !> The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The buffer value for the scattering
    real(r8), intent(inout) :: buf
    !> The linewidth
    real(r8), dimension(:, :), intent(in) :: lw
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The temperature
    real(r8) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> buffers
    real(r8) :: buf_iso, om1, n1, sig1, om2, n2, sig2, sigma, f0
    !> Integers
    integer :: i, j, q1, q2, b1, b2, q2i

    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    n1 = lo_planck(temperature, om1)
    sig1 = lw(q1, b1)

    do j=1, sr%iso(il)%n
        q2 = sr%iso(il)%event(j)%q2
        b2 = sr%iso(il)%event(j)%b2
        q2i = qp%ap(q2)%irreducible_index

        om2 = dr%aq(q2)%omega(b2)
        n2 = lo_planck(temperature, om2)
        sig2 = lw(q2i, b2)
        sigma = 2.0_r8 * (sig1 + sig2)

        f0 = sr%iso(il)%event(j)%psisq * lo_lorentz(om1, om2, sigma)
        f0 = f0 * om1 * om2 * n1 * (n2 + 1.0_r8)
        buf = buf + f0
        sr%iso(il)%event(j)%W = f0
    end do
end subroutine

subroutine compute_scatteringrate_threephonon(il, buf, lw, dr, qp, sr, temperature, mw, mem)
    !> The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The buffer value for the scattering
    real(r8), intent(inout) :: buf
    !> The linewidth
    real(r8), dimension(:, :), intent(in) :: lw
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Scattering rate
    type(lo_scattering_rates), intent(inout) :: sr
    ! The temperature
    real(r8), intent(in) :: temperature
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> buffers
    real(r8) :: om1, n1, sig1, om2, n2, sig2, om3, n3, sig3, sigma, plf, f0
    !> Integers
    integer :: j, q1, q2, q3, b1, b2, b3, q2i, q3i

    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    n1 = lo_planck(temperature, om1)

    do j=1, sr%threephonon(il)%nplus
        q2 = sr%threephonon(il)%plus(j)%q2
        b2 = sr%threephonon(il)%plus(j)%b2
        q3 = sr%threephonon(il)%plus(j)%q3
        b3 = sr%threephonon(il)%plus(j)%b3
        q2i = qp%ap(q2)%irreducible_index
        q3i = qp%ap(q3)%irreducible_index

        om2 = dr%aq(q2)%omega(b2)
        om3 = dr%aq(q3)%omega(b3)
        n2 = lo_planck(temperature, om2)
        n3 = lo_planck(temperature, om3)
        sig2 = lw(q2i, b2)
        sig3 = lw(q3i, b3)
        sigma = 2.0_r8 * (sig2 + sig3)

        f0 = sr%threephonon(il)%plus(j)%psisq * lo_lorentz(om1, -om2 + om3, sigma)
        f0 = f0 * n1 * n2 * (n3 + 1.0_r8)
        buf = buf + f0
        sr%threephonon(il)%plus(j)%W = f0
    end do

    do j=1, sr%threephonon(il)%nminus
        q2 = sr%threephonon(il)%minus(j)%q2
        b2 = sr%threephonon(il)%minus(j)%b2
        q3 = sr%threephonon(il)%minus(j)%q3
        b3 = sr%threephonon(il)%minus(j)%b3
        q2i = qp%ap(q2)%irreducible_index
        q3i = qp%ap(q3)%irreducible_index

        om2 = dr%aq(q2)%omega(b2)
        om3 = dr%aq(q3)%omega(b3)
        n2 = lo_planck(temperature, om2)
        n3 = lo_planck(temperature, om3)
        sig2 = lw(q2i, b2)
        sig3 = lw(q3i, b3)
        sigma = 2.0_r8 * (sig2 + sig3)

        f0 = sr%threephonon(il)%minus(j)%psisq * lo_lorentz(om1, om2 + om3, sigma)
        f0 = f0 * n1 * (n2 + 1.0_r8) * (n3 + 1.0_r8)
        buf = buf + f0 * 0.5_r8
        sr%threephonon(il)%minus(j)%W = f0
    end do
end subroutine

subroutine compute_scatteringrate_fourphonon(il, buf, lw, dr, qp, sr, temperature, mw, mem)
    !> The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The buffer value for the scattering
    real(r8), intent(inout) :: buf
    !> The linewidth
    real(r8), dimension(:, :), intent(in) :: lw
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The temperature
    real(r8) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> buffers
    real(r8) :: buf_pp, buf_pm, buf_mm, om1, n1, sig1, om2, n2, sig2, om3, n3, sig3, om4, n4, sig4, sigma, plf, f0
    !> Integers
    integer :: j, q1, q2, q3, q4, b1, b2, b3, b4, q2i, q3i, q4i

    do j=1, sr%fourphonon(il)%npp
        buf = buf + sr%fourphonon(il)%pp(j)%W * 0.5_r8
    end do
    do j=1, sr%fourphonon(il)%npm
        buf = buf + sr%fourphonon(il)%pm(j)%W * 0.5_r8
    end do
    do j=1, sr%fourphonon(il)%nmm
        buf = buf + sr%fourphonon(il)%mm(j)%W / 6.0_r8
    end do

!   q1 = sr%q1(il)
!   b1 = sr%b1(il)
!   om1 = dr%iq(q1)%omega(b1)
!   n1 = lo_planck(temperature, om1)
!   do i=1, sr%fourphonon(il)%n
!       q2 = sr%fourphonon(il)%q2(i)
!       q3 = sr%fourphonon(il)%q3(i)
!       q4 = sr%fourphonon(il)%q4(i)
!       q2i = qp%ap(q2)%irreducible_index
!       q3i = qp%ap(q3)%irreducible_index
!       q4i = qp%ap(q4)%irreducible_index

!       do b2=1, dr%n_mode
!           do b3=1, dr%n_mode
!               do b4=1, dr%n_mode

!                   om2 = dr%aq(q2)%omega(b2)
!                   om3 = dr%aq(q3)%omega(b3)
!                   om4 = dr%aq(q4)%omega(b4)

!                   if (any([om2, om3, om4] .lt. lo_freqtol)) cycle

!   !               n2 = lo_planck(temperature, om2)
!   !               n3 = lo_planck(temperature, om3)
!   !               n4 = lo_planck(temperature, om4)
!   !               sig2 = lw(q2i, b2)
!   !               sig3 = lw(q3i, b3)
!   !               sig4 = lw(q4i, b4)




    !           sig1 = qp%adaptive_sigma(qp%ip(q1)%radius, dr%iq(q1)%vel(:, b1), dr%default_smearing(b1), 1.0_r8)
    !           sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)
    !           sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, dr%aq(q3)%vel(:, b3), dr%default_smearing(b3), 1.0_r8)
    !           sig4 = qp%adaptive_sigma(qp%ap(q4)%radius, dr%aq(q4)%vel(:, b4), dr%default_smearing(b4), 1.0_r8)
    !           sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)



                    !sigma = 2.0_r8 * (sig2 + sig3 + sig4)

    !               plf = n1 * n2 * n3 * (n4 + 1.0_r8) * lo_gauss(om1, -om2 - om3 + om4, sigma)
    !               !plf = n1 * n2 * n3 * (n4 + 1.0_r8) * lo_lorentz(om1, -om2 - om3 + om4, sigma)
    !               f0 = plf * sr%fourphonon(il)%psisq(b2, b3, b4, i)
    !               buf = buf + f0 * sr%mle_ratio4ph * 0.5_r8
    !               sr%fourphonon(il)%Wpp(b2, b3, b4, i) = f0

    !               plf = n1 * n2 * (n3 + 1.0_r8) * (n4 + 1.0_r8) * lo_gauss(om1, -om2 + om3 + om4, sigma)
    !               !plf = n1 * n2 * (n3 + 1.0_r8) * (n4 + 1.0_r8) * lo_lorentz(om1, -om2 + om3 + om4, sigma)
    !               f0 = plf * sr%fourphonon(il)%psisq(b2, b3, b4, i)
    !               buf = buf + f0 * sr%mle_ratio4ph * 0.5_r8
    !               sr%fourphonon(il)%Wpm(b2, b3, b4, i) = f0

    !               plf = n1 * (n2 + 1.0_r8) * (n3 + 1.0_r8) * (n4 + 1.0_r8) * lo_gauss(om1, om2 + om3 + om4, sigma)
    !               !plf = n1 * (n2 + 1.0_r8) * (n3 + 1.0_r8) * (n4 + 1.0_r8) * lo_lorentz(om1, om2 + om3 + om4, sigma)
    !               f0 = plf * sr%fourphonon(il)%psisq(b2, b3, b4, i)
    !               buf = buf + f0 * sr%mle_ratio4ph / 6.0_r8
    !               sr%fourphonon(il)%Wmm(b2, b3, b4, i) = f0

!                   buf = buf + sr%fourphonon(il)%Wmm(b2, b3, b4, i) * 0.5_r8
!                   buf = buf + sr%fourphonon(il)%Wmm(b2, b3, b4, i) * 0.5_r8
!                   buf = buf + sr%fourphonon(il)%Wmm(b2, b3, b4, i) / 6.0_r8
!               end do
!           end do
!       end do
!   end do
end subroutine


subroutine compute_qs(dr, qp, temperature, mfpmax, fourthorder)
    !> dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> q-mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> temperature
    real(r8), intent(in) :: temperature
    ! Maximum mean free path
    real(r8), intent(in) :: mfpmax
    ! Did we do fourthorder ?
    logical, intent(in) :: fourthorder

    ! some buffer
    real(r8) :: qs_boundary, n1, velnorm
    ! Some integers for the do loops
    integer :: q1, b1

    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! Skip gamma for acoustic branches
            if (dr%iq(q1)%omega(b1) .lt. lo_freqtol) cycle
            ! First we get the mfp and F0
            n1 = lo_planck(temperature, dr%iq(q1)%omega(b1))
            velnorm = norm2(dr%iq(q1)%vel(:, b1))
            if (mfpmax .gt. 0.0_r8) then
                qs_boundary = n1 * (n1 + 1.0_r8) * velnorm / mfpmax
            else
                qs_boundary = 0.0_r8
            end if

            dr%iq(q1)%qs(b1) = dr%iq(q1)%p_plus(b1) + &
                               dr%iq(q1)%p_minus(b1) * 0.5_r8 + &
                               dr%iq(q1)%p_iso(b1) + &
                               qs_boundary
            if (fourthorder) then
                dr%iq(q1)%qs(b1) = dr%iq(q1)%qs(b1) + dr%iq(q1)%p_plusplus(b1) * 0.5_r8 + &
                                                      dr%iq(q1)%p_plusminus(b1) * 0.5_r8 + &
                                                      dr%iq(q1)%p_minusminus(b1) / 6.0_r8
            end if
            dr%iq(q1)%linewidth(b1) = 0.5_r8 * dr%iq(q1)%qs(b1) / (n1 * (n1 + 1.0_r8))

            if (velnorm .gt. lo_phonongroupveltol) then
                dr%iq(q1)%mfp(:, b1) = dr%iq(q1)%vel(:, b1) * 0.5_r8 / dr%iq(q1)%linewidth(b1)
                dr%iq(q1)%scalar_mfp(b1) = velnorm * 0.5_r8 / dr%iq(q1)%linewidth(b1)
                dr%iq(q1)%F0(:, b1) = dr%iq(q1)%mfp(:, b1) * dr%iq(q1)%omega(b1) / temperature
                dr%iq(q1)%Fn(:, b1) = dr%iq(q1)%F0(:, b1)
            end if
        end do
    end do
end subroutine

end module
