#include "precompilerdefinitions"
module scattering
use konstanter, only: r8, lo_temperaturetol, lo_status, lo_kappa_au_to_SI, lo_freqtol, lo_twopi, &
                      lo_frequency_THz_to_Hartree, lo_exitcode_param, lo_pi
use gottochblandat, only: lo_planck, walltime, lo_gauss, lo_trueNtimes, lo_progressbar_init, lo_progressbar
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_randomnumbers, only: lo_mersennetwister

implicit none

private
public :: compute_scattering_threephonon
public :: compute_scattering_fourphonon
public :: compute_scattering_isotopes

contains
! Compute the scattering for three phonons
subroutine compute_scattering_threephonon(qp, dr, fct, temperature, nsample3ph, integrationtype, &
                                          smearing_prefactor, thres, rng, mw, mem)
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Third order forceconstants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    ! The temperature
    real(r8), intent(in) :: temperature
    ! The number of 3 phonon scattering actually computed
    integer, intent(in) :: nsample3ph
    ! Integration type
    integer, intent(in) :: integrationtype
    ! The smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    ! The gaussian threshold
    real(r8), intent(in) :: thres
    ! Random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: threephonon_prefactor = lo_pi / 16.0_r8
    ! The number of scattering
    integer, dimension(:, :), allocatable :: nplus_tot, nplus_sample, nminus_tot, nminus_sample
    ! The scattering
    real(r8), dimension(:, :), allocatable :: scatplus, scatminus
    ! Eigenvectors
    complex(r8), dimension(:, :), allocatable :: egv
    ! The full qpoint grid
    integer, dimension(:), allocatable :: qgridfull
    ! The frequencies
    real(r8), dimension(3) :: omega, qvec
    ! The grid dimensions
    integer, dimension(3) :: dims
    ! Some buffers
    real(r8) :: bufplus, bufminus, plf, rnd, omthres, deltafunction, om1, om2, om3, n1, n2, n3, psisq
    real(r8) :: sigma, sig1, sig2, sig3
    real(r8), dimension(3) :: qv2, qv3, qv3frac, vel1, vel2, vel3
    integer :: np_count, nps_count, nm_count, nms_count
    real(r8) :: t0
    complex(r8) :: c0
    logical :: isplus, isminus

    integer :: i, q1, b1, q2, b2, q3, b3, ctr

    t0 = walltime()
    if (mw%talk) call lo_progressbar_init()

    ! grid dimensions
    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

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

    call mem%allocate(egv, [dr%n_mode, 3], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! We start by creating and shuffling an array with the qpoint index
    call mem%allocate(qgridfull, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    qgridfull = (/ (i, i=1, qp%n_full_point) /)

    omthres = dr%omega_min*0.5_r8
    ctr = 0
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! MPI
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. omthres) cycle
            vel1 = dr%iq(q1)%vel(:, b1)
            egv(:, 1) = dr%iq(q1)%egv(:, b1)
            omega(1) = om1
            call rng%shuffle_int_array(qgridfull)
            do i = 1, qp%n_full_point
                q2 = qgridfull(i)
                q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)
                qv2 = -qp%ap(q2)%r*lo_twopi
                qv3 = -qp%ap(q3)%r*lo_twopi
                do b2 = 1, dr%n_mode
                    om2 = dr%aq(q2)%omega(b2)
                    if (om2 .lt. omthres) cycle
                    omega(2) = om2
                    egv(:, 2) = dr%aq(q2)%egv(:, b2)
                    vel2 = dr%aq(q2)%vel(:, b2)
                    n2 = lo_planck(temperature, om2)
                    do b3 = 1, dr%n_mode
                        om3 = dr%aq(q3)%omega(b3)
                        if (om3 .lt. omthres) cycle
                        omega(3) = om3
                        n3 = lo_planck(temperature, om3)
                        egv(:, 3) = dr%aq(q3)%egv(:, b3)

                        isplus = .false.
                        isminus = .false.
                        if (nplus_sample(q1, b1) .lt. nsample3ph) isplus = .true.
                        if (nminus_sample(q1, b1) .lt. nsample3ph) isminus = .true.
                        vel3 = dr%aq(q3)%vel(:, b3)

                        select case (integrationtype)
                        case (1)
                            sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
                        case (2)
                            sig1 = qp%adaptive_sigma(qp%ip(q1)%radius, vel1, dr%default_smearing(b1), smearing_prefactor)
                            sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, vel2, dr%default_smearing(b2), smearing_prefactor)
                            sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, vel3, dr%default_smearing(b3), smearing_prefactor)
                            sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                        end select

                        if (abs(om1 + om2 - om3) .lt. thres*sigma .or. &
                            abs(om1 - om2 - om3) .lt. thres*sigma) then
                            if (isplus .or. isminus) then
                                ! Let's compute the scattering strength only once
                                c0 = fct%scatteringamplitude(omega, egv, qv2, qv3)
                                psisq = abs(c0*conjg(c0))
                            end if
                        else
                            cycle
                        end if

                        if (abs(om1 + om2 - om3) .lt. thres*sigma) then
                            nplus_tot(q1, b1) = nplus_tot(q1, b1) + 1
                            if (isplus) then
                                plf = 2.0_r8 * (n2 - n3) * threephonon_prefactor
                                nplus_sample(q1, b1) = nplus_sample(q1, b1) + 1
                                deltafunction = lo_gauss(om1, -om2 + om3, sigma)
                                scatplus(q1, b1) = scatplus(q1, b1) + deltafunction*psisq*qp%ap(q2)%integration_weight*plf
                            end if
                        end if

                        if (abs(om1 - om2 - om3) .lt. thres*sigma) then
                            nminus_tot(q1, b1) = nminus_tot(q1, b1) + 1
                            if (isminus) then
                                plf = (n2 + n3 + 1) * threephonon_prefactor
                                nminus_sample(q1, b1) = nminus_sample(q1, b1) + 1
                                deltafunction = lo_gauss(om1, om2 + om3, sigma)
                                scatminus(q1, b1) = scatminus(q1, b1) + 0.5_r8*deltafunction*psisq*qp%ap(q2)%integration_weight*plf
                            end if
                        end if
                    end do ! b3
                end do ! b2
            end do ! q2
            if (mw%talk) then
                if (lo_trueNtimes(ctr, 127, qp%n_irr_point*dr%n_mode)) call lo_progressbar(' ... threephonon scattering', ctr, dr%n_mode*qp%n_irr_point)
            end if
        end do ! b1
    end do ! q1
    call mem%deallocate(qgridfull, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (mw%talk) call lo_progressbar(' ... threephonon scattering', dr%n_mode*qp%n_irr_point, dr%n_mode*qp%n_irr_point, walltime() - t0)

    ! Sum every rank
    call mw%allreduce('sum', nplus_tot)
    call mw%allreduce('sum', nplus_sample)
    call mw%allreduce('sum', scatplus)
    call mw%allreduce('sum', nminus_tot)
    call mw%allreduce('sum', nminus_sample)
    call mw%allreduce('sum', scatminus)

    np_count = 0
    nps_count = 0
    nm_count = 0
    nms_count = 0
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .lt. omthres) cycle
            bufplus = 0.0_r8
            bufminus = 0.0_r8
            if (nplus_sample(q1, b1) .ne. 0) then
                bufplus = bufplus + scatplus(q1, b1) * real(nplus_tot(q1, b1), r8) / real(nplus_sample(q1, b1), r8)
                np_count = np_count + nplus_tot(q1, b1)
                nps_count = nps_count + nplus_sample(q1, b1)
            end if
            if (nminus_sample(q1, b1) .ne. 0) then
                bufminus = bufminus + scatminus(q1, b1) * real(nminus_tot(q1, b1), r8) / real(nminus_sample(q1, b1), r8)
                nm_count = nm_count + nminus_tot(q1, b1)
                nms_count = nms_count + nminus_sample(q1, b1)
            end if
            ! Direct application of Mathiessen's rule
            dr%iq(q1)%linewidth(b1) = dr%iq(q1)%linewidth(b1) + bufplus + bufminus
            dr%iq(q1)%p_plus(b1) = bufplus * 2 * n1 * (n1 + 1)
            dr%iq(q1)%p_minus(b1) = bufminus * 2 * n1 * (n1 + 1)
        end do
    end do
    if (mw%talk) then
        write(*, *) 'Number of + scattering event ', np_count, 'sampled ', nps_count, 'ratio ', real(nps_count, r8)/real(np_count, r8)
        write(*, *) 'Number of - scattering event ', nm_count, 'sampled ', nms_count, 'ratio ', real(nms_count, r8)/real(nm_count, r8)
    end if

    ! Deallocate things
    call mem%deallocate(nplus_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(nplus_sample, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(scatplus, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(nminus_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(nminus_sample, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(scatminus, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

! Compute the scattering for isotopic defects
subroutine compute_scattering_isotopes(qp, dr, uc, temperature, integrationtype, smearing_prefactor, thres, mw, mem)
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    ! The temperature
    real(r8), intent(in) :: temperature
    ! Integration type
    integer, intent(in) :: integrationtype
    ! The smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    ! The gaussian threshold
    real(r8), intent(in) :: thres
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Isotope prefactor
    real(r8), parameter :: isotope_prefactor = lo_pi / 4.0_r8

    real(r8), dimension(:, :), allocatable :: iso_scatter
    integer :: q1, b1, q2, b2, ctr
    real(r8) :: n1, om1, om2, sigma, scatterstrength, omthres, deltafunction, t0
    complex(r8), dimension(uc%na*3, 2) :: egviso
    integer :: n_count

    t0 = walltime()
    if (mw%talk) call lo_progressbar_init()

    call mem%allocate(iso_scatter, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    iso_scatter = 0.0_r8

    omthres = dr%omega_min*0.5_r8
    n_count = 0
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
                select case (integrationtype)
                case (1)
                    sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
                case (2)
                    sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), smearing_prefactor)
                end select

                if (abs(om1 - om2) .lt. thres*sigma) then
                    deltafunction = lo_gauss(om1, om2, sigma)
                    egviso(:, 1) = dr%iq(q1)%egv(:, b1)
                    egviso(:, 2) = dr%aq(q2)%egv(:, b2)
                    scatterstrength = isotope_scattering_strength(uc, egviso)*om1*om2
                    iso_scatter(q1, b1) = iso_scatter(q1, b1) + scatterstrength*deltafunction*qp%ap(q2)%integration_weight * isotope_prefactor
                    n_count = n_count + 1
                end if
            end do
            end do
            if (mw%talk) then
                if (lo_trueNtimes(ctr, 127, qp%n_irr_point*dr%n_mode)) call lo_progressbar(' ... isotope scattering', ctr, dr%n_mode*qp%n_irr_point)
            end if
        end do
    end do

    if (mw%talk) call lo_progressbar(' ... isotope scattering', dr%n_mode*qp%n_irr_point, dr%n_mode*qp%n_irr_point, walltime() - t0)

    call mw%allreduce('sum', iso_scatter)
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            dr%iq(q1)%linewidth(b1) = dr%iq(q1)%linewidth(b1) + iso_scatter(q1, b1)
            n1 = lo_planck(temperature, dr%iq(q1)%omega(b1))
            dr%iq(q1)%p_iso(b1) = iso_scatter(q1, b1) * 2 * n1 * (n1 + 1)
        end do
    end do
    call mem%deallocate(iso_scatter, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (mw%talk) then
        write(*, *) 'Number of isotope scattering event ', n_count
    end if

end subroutine

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

!> returns the index on the grid that gives q3=-q1-q2
pure function fft_third_grid_index(i1, i2, dims) result(i3)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q3
    integer :: i3

    integer, dimension(3) :: gi1, gi2, gi3
    integer :: l, k

    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    do l = 1, 3
        gi3(l) = 3 - gi1(l) - gi2(l)
    end do
    do k = 1, 3
    do l = 1, 3
        if (gi3(l) .lt. 1) gi3(l) = gi3(l) + dims(l)
        if (gi3(l) .gt. dims(l)) gi3(l) = gi3(l) - dims(l)
    end do
    end do
    ! convert it back to a singlet
    i3 = triplet_to_singlet(gi3, dims(2), dims(3))

contains
    !> convert a linear index to a triplet
    pure function singlet_to_triplet(l, ny, nz) result(gi)
        !> linear index
        integer, intent(in) :: l
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> grid-index
        integer, dimension(3) :: gi

        integer :: i, j, k

        k = mod(l, nz)
        if (k .eq. 0) k = nz
        j = mod((l - k)/nz, ny) + 1
        i = (l - k - (j - 1)*nz)/(nz*ny) + 1
        gi = [i, j, k]
    end function
    !> convert a triplet index to a singlet
    pure function triplet_to_singlet(gi, ny, nz) result(l)
        !> grid-index
        integer, dimension(3), intent(in) :: gi
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> linear index
        integer :: l

        l = (gi(1) - 1)*ny*nz + (gi(2) - 1)*nz + gi(3)
    end function
end function

! Compute the scattering for three phonons
subroutine compute_scattering_fourphonon(qp, dr, fcf, temperature, nsample4ph, integrationtype, &
                                         smearing_prefactor, thres, rng, mw, mem)
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    ! Third order forceconstants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    ! The temperature
    real(r8), intent(in) :: temperature
    ! The number of 3 phonon scattering actually computed
    integer, intent(in) :: nsample4ph
    ! Integration type
    integer, intent(in) :: integrationtype
    ! The smearing prefactor
    real(r8), intent(in) :: smearing_prefactor
    ! The gaussian threshold
    real(r8), intent(in) :: thres
    ! Random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    ! Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    ! Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    ! Three phonon prefactor
    real(r8), parameter :: fourphonon_prefactor = lo_pi / 96.0_r8
    ! The number of scattering
    integer, dimension(:, :), allocatable :: npp_tot, npp_sample, npm_tot, npm_sample, nmm_tot, nmm_sample
    ! The scattering
    real(r8), dimension(:, :), allocatable :: scatpp, scatpm, scatmm
    ! Eigenvectors
    complex(r8), dimension(:, :), allocatable :: egv
    ! The full qpoint grid
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    ! The frequencies
    real(r8), dimension(4) :: omega
    ! the qvector
    real(r8), dimension(3) :: qvec
    ! The grid dimensions
    integer, dimension(3) :: dims
    ! Some buffers
    real(r8) :: bufpp, bufpm, bufmm, plf, rnd, omthres, deltafunction, om1, om2, om3, om4, n2, n3, n4, psisq
    real(r8) :: sigma, sig1, sig2, sig3, sig4
    real(r8), dimension(3) :: qv2, qv3, qv4, vel1, vel2, vel3, vel4
    integer :: npp_count, npps_count, npm_count, npms_count, nmm_count, nmms_count
    real(r8) :: t0
    complex(r8) :: c0
    logical :: ispp, ispm, ismm

    integer :: i, j, q1, b1, q2, b2, q3, b3, q4, b4, ctr

    t0 = walltime()
    if (mw%talk) call lo_progressbar_init()

    ! grid dimensions
    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    call mem%allocate(npp_tot, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(npp_sample, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(scatpp, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(npm_tot, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(npm_sample, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(scatpm, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(nmm_tot, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(nmm_sample, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(scatmm, [qp%n_irr_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    npp_tot = 0
    npp_sample = 0
    scatpp = 0.0_r8
    npm_tot = 0
    npm_sample = 0
    scatpm = 0.0_r8
    nmm_tot = 0
    nmm_sample = 0
    scatmm = 0.0_r8

    call mem%allocate(egv, [dr%n_mode, 4], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    qgridfull1 = (/ (i, i=1, qp%n_full_point) /)
    qgridfull2 = (/ (i, i=1, qp%n_full_point) /)

    omthres = dr%omega_min*0.5_r8
    ctr = 0
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! MPI
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            om1 = dr%iq(q1)%omega(b1)
            vel1 = dr%iq(q1)%vel(:, b1)
            if (om1 .lt. omthres) cycle
            call rng%shuffle_int_array(qgridfull1)
            call rng%shuffle_int_array(qgridfull2)
            do i = 1, qp%n_full_point
            do j = 1, qp%n_full_point
                q2 = qgridfull1(i)
                q3 = qgridfull1(j)
                do b2 = 1, dr%n_mode
                do b3 = 1, dr%n_mode
                    ispp = .false.
                    ispm = .false.
                    ismm = .false.
                    rnd = rng%rnd_real()
                    if (npp_sample(q1, b1) .lt. nsample4ph) ispp = .true.
                    if (npm_sample(q1, b1) .lt. nsample4ph) ispm = .true.
                    if (nmm_sample(q1, b1) .lt. nsample4ph) ismm = .true.

                    q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, dims)
                    om2 = dr%aq(q2)%omega(b2)
                    om3 = dr%aq(q3)%omega(b3)
                    om4 = dr%aq(q4)%omega(b4)
                    if (minval([om2, om3, om4]) .lt. omthres) cycle
                    vel2 = dr%aq(q2)%vel(:, b2)
                    vel3 = dr%aq(q3)%vel(:, b3)
                    vel4 = dr%aq(q4)%vel(:, b4)

                    n2 = lo_planck(temperature, om2)
                    n3 = lo_planck(temperature, om3)
                    n4 = lo_planck(temperature, om4)
                    omega(1) = om1
                    omega(2) = om2
                    omega(3) = om3
                    omega(4) = om4
                    egv(:, 1) = dr%iq(q1)%egv(:, b1)
                    egv(:, 2) = dr%aq(q2)%egv(:, b2)
                    egv(:, 3) = dr%aq(q3)%egv(:, b3)
                    egv(:, 4) = dr%aq(q4)%egv(:, b4)
                    qv2 = -qp%ap(q2)%r*lo_twopi
                    qv3 = -qp%ap(q3)%r*lo_twopi
                    qv4 = -qp%ap(q4)%r*lo_twopi

                    select case (integrationtype)
                    case (1)
                        sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing_prefactor
                    case (2)
                        sig1 = qp%adaptive_sigma(qp%ip(q1)%radius, vel1, dr%default_smearing(b1), smearing_prefactor)
                        sig2 = qp%adaptive_sigma(qp%ap(q2)%radius, vel2, dr%default_smearing(b2), smearing_prefactor)
                        sig3 = qp%adaptive_sigma(qp%ap(q3)%radius, vel3, dr%default_smearing(b3), smearing_prefactor)
                        sig4 = qp%adaptive_sigma(qp%ap(q4)%radius, vel4, dr%default_smearing(b4), smearing_prefactor)
                        sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)
                    end select

                    if (abs(om1 + om2 + om3 - om4) .lt. thres*sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. thres*sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. thres*sigma) then
                        if (any([ispp, ispm, ismm])) then
                            ! Compute the scattering matrix element
                            c0 = fcf%scatteringamplitude(omega, egv, qv2, qv3, qv4)
                            psisq = abs(c0*conjg(c0))
                        end if
                    else
                        cycle
                    end if

                    if (abs(om1 + om2 + om3 - om4) .lt. thres*sigma) then
                        npp_tot(q1, b1) = npp_tot(q1, b1) + 1
                        if (ispp) then
                            plf = n2 * n3 * (n4 + 1) * fourphonon_prefactor
                            npp_sample(q1, b1) = npp_sample(q1, b1) + 1
                            deltafunction = lo_gauss(om1, -om2 - om3 + om4, sigma)
                            ! scatpp(q1, b1) = scatpp(q1, b1) + 0.5_r8 * deltafunction*psisq*qp%ap(q2)%integration_weight*&
                            !                  qp%ap(q3)%integration_weight*plf
                            scatpp(q1, b1) = scatpp(q1, b1) + 0.5_r8 * deltafunction*psisq*plf*&
                                             qp%ap(q2)%integration_weight
                        end if
                    end if

                    if (abs(om1 + om2 - om3 - om4) .lt. thres*sigma) then
                        npm_tot(q1, b1) = npm_tot(q1, b1) + 1
                        if (ispm) then
                            plf = n2 * (n3 + 1) * (n4 + 1) * fourphonon_prefactor
                            npm_sample(q1, b1) = npm_sample(q1, b1) + 1
                            deltafunction = lo_gauss(om1, -om2 + om3 + om4, sigma)
                            ! scatpm(q1, b1) = scatpm(q1, b1) + 0.5_r8 * deltafunction*psisq*qp%ap(q2)%integration_weight*&
                            !                  qp%ap(q3)%integration_weight*plf
                            scatpm(q1, b1) = scatpm(q1, b1) + 0.5_r8 * deltafunction*psisq*plf*&
                                             qp%ap(q2)%integration_weight
                        end if
                    end if

                    if (abs(om1 - om2 - om3 - om4) .lt. thres*sigma) then
                        nmm_tot(q1, b1) = nmm_tot(q1, b1) + 1
                        if (ismm) then
                            plf = (n2 + 1) * (n3 + 1) * (n4 + 1) * fourphonon_prefactor
                            nmm_sample(q1, b1) = nmm_sample(q1, b1) + 1
                            deltafunction = lo_gauss(om1, om2 + om3 + om4, sigma)
                            ! scatmm(q1, b1) = scatmm(q1, b1) + deltafunction*psisq*qp%ap(q2)%integration_weight*&
                            !                  qp%ap(q3)%integration_weight*plf / 6.0_r8
                            scatmm(q1, b1) = scatmm(q1, b1) + deltafunction*psisq*plf*&
                                             qp%ap(q2)%integration_weight / 6.0_r8
                        end if
                    end if
                end do ! b3
                end do ! b2
            end do !q3
            end do !q2
            if (mw%talk) then
                if (lo_trueNtimes(ctr, 127, qp%n_irr_point*dr%n_mode)) call lo_progressbar(' ... fourphonon scattering', ctr, dr%n_mode*qp%n_irr_point)
            end if
        end do ! b1
    end do ! q1
    call mem%deallocate(qgridfull1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (mw%talk) call lo_progressbar(' ... fourphonon scattering', dr%n_mode*qp%n_irr_point, dr%n_mode*qp%n_irr_point, walltime() - t0)

    call mw%allreduce('sum', npp_tot)
    call mw%allreduce('sum', npp_sample)
    call mw%allreduce('sum', scatpp)
    call mw%allreduce('sum', npm_tot)
    call mw%allreduce('sum', npm_sample)
    call mw%allreduce('sum', scatpm)
    call mw%allreduce('sum', nmm_tot)
    call mw%allreduce('sum', nmm_sample)
    call mw%allreduce('sum', scatmm)

    npp_count = 0
    npps_count = 0
    npm_count = 0
    npms_count = 0
    nmm_count = 0
    nmms_count = 0
    do q1 = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            if (dr%iq(q1)%omega(b1) .lt. omthres) cycle
            bufpp = 0.0_r8
            bufpm = 0.0_r8
            bufmm = 0.0_r8
            if (npp_sample(q1, b1) .ne. 0) then ! shouldn't matter, but we don't want to divide by zero
                bufpp = bufpp + scatpp(q1, b1) * real(npp_tot(q1, b1), r8) / real(npp_sample(q1, b1), r8)
                npp_count = npp_count + npp_tot(q1, b1)
                npps_count = npps_count + npp_sample(q1, b1)
            end if
            if (npm_sample(q1, b1) .ne. 0) then ! shouldn't matter, but we don't want to divide by zero
                bufpm = bufpm + scatpm(q1, b1) * real(npm_tot(q1, b1), r8) / real(npm_sample(q1, b1), r8)
                npm_count = npm_count + npm_tot(q1, b1)
                npms_count = npms_count + npm_sample(q1, b1)
            end if
            if (nmm_sample(q1, b1) .ne. 0) then ! shouldn't matter, but we don't want to divide by zero
                bufmm = bufmm + scatmm(q1, b1) * real(nmm_tot(q1, b1), r8) / real(nmm_sample(q1, b1), r8)
                nmm_count = nmm_count + nmm_tot(q1, b1)
                nmms_count = nmms_count + nmm_sample(q1, b1)
            end if
            ! Direct application of Mathiessen's rule
            dr%iq(q1)%linewidth(b1) = dr%iq(q1)%linewidth(b1) + bufpp + bufpm + bufmm
        end do
    end do
    if (mw%talk) then
        write(*, *) 'Number of ++ scattering event ', npp_count, 'sampled ', npps_count, 'ratio ', real(npps_count, r8)/real(npp_count, r8)
        write(*, *) 'Number of +- scattering event ', npm_count, 'sampled ', npms_count, 'ratio ', real(npms_count, r8)/real(npm_count, r8)
        write(*, *) 'Number of -- scattering event ', nmm_count, 'sampled ', nmms_count, 'ratio ', real(nmms_count, r8)/real(nmm_count, r8)
    end if

    call mem%deallocate(npp_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(npp_sample, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(scatpp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(npm_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(npm_sample, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(scatpm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(nmm_tot, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(nmm_sample, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(scatmm, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

end subroutine

!> returns the index on the grid that gives q4=-q3-q2-q1
pure function fft_fourth_grid_index(i1, i2, i3, dims) result(i4)
    !> index to q1
    integer, intent(in) :: i1
    !> index to q2
    integer, intent(in) :: i2
    !> index to q3
    integer, intent(in) :: i3
    !> dimensions of the grid
    integer, dimension(3), intent(in) :: dims
    !> index to q4
    integer :: i4

    integer, dimension(3) :: gi1, gi2, gi3, gi4
    integer :: l, k
    ! Convert triplet to singlet
    gi1 = singlet_to_triplet(i1, dims(2), dims(3))
    gi2 = singlet_to_triplet(i2, dims(2), dims(3))
    gi3 = singlet_to_triplet(i3, dims(2), dims(3))
    do l = 1, 3
         gi4(l) = 4 - gi1(l) - gi2(l) - gi3(l)
   end do
   do k = 1, 3
   do l = 1, 3
       if (gi4(l) .lt. 1) gi4(l) = gi4(l) + dims(l)
       if (gi4(l) .gt. dims(l)) gi4(l) = gi4(l) - dims(l)
   end do
   end do

    ! convert it back to a singlet
    i4 = triplet_to_singlet(gi4, dims(2), dims(3))

contains
    !> convert a linear index to a triplet
    pure function singlet_to_triplet(l, ny, nz) result(gi)
        !> linear index
        integer, intent(in) :: l
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> grid-index
        integer, dimension(3) :: gi

        integer :: i, j, k

        k = mod(l, nz)
        if (k .eq. 0) k = nz
        j = mod((l - k)/nz, ny) + 1
        i = (l - k - (j - 1)*nz)/(nz*ny) + 1
        gi = [i, j, k]
    end function
    !> convert a triplet index to a singlet
    pure function triplet_to_singlet(gi, ny, nz) result(l)
        !> grid-index
        integer, dimension(3), intent(in) :: gi
        !> second dimension
        integer, intent(in) :: ny
        !> third dimension
        integer, intent(in) :: nz
        !> linear index
        integer :: l

        l = (gi(1) - 1)*ny*nz + (gi(2) - 1)*nz + gi(3)
    end function
end function
end module
