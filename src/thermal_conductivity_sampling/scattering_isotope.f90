

subroutine compute_isotope_scattering(il, sr, qp, dr, uc, temperature, mw, mem)
    !> The local point
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> phonon dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Isotope prefactor
    ! real(r8), parameter :: isotope_prefactor = lo_pi * 0.5_r8
    ! Eigenvectors
    complex(r8), dimension(uc%na*3, 2) :: egviso
    ! prefactor and phonon buffers
    real(r8) :: om1, om2, sig1, sig2, sigma, psisq, prefactor
    ! Integers for do loops
    integer :: q1, b1, q2, b2, i, niso

    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)

    niso = 0
    do q2=1, qp%n_full_point
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)
            if (abs(om1 - om2) .lt. 4.0_r8 * sigma) niso = niso + 1
        end do
    end do

    sr%iso(il)%n = niso
    allocate(sr%iso(il)%psisq(niso))
    allocate(sr%iso(il)%q2(niso))
    allocate(sr%iso(il)%b2(niso))

    om1 = dr%iq(q1)%omega(b1)
    egviso(:, 1) = dr%iq(q1)%egv(:, b1)

    i = 0
    do q2=1, qp%n_full_point
        prefactor = isotope_prefactor * qp%ap(q2)%integration_weight
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), dr%default_smearing(b2), 1.0_r8)
            if (abs(om1 - om2) .lt. 4.0_r8 * sigma) then
                i = i + 1

                egviso(:, 2) = dr%aq(q2)%egv(:, b2)
                psisq = isotope_scattering_strength(uc, egviso) * prefactor

                sr%iso(il)%q2(i) = q2
                sr%iso(il)%b2(i) = b2
                sr%iso(il)%psisq(i) = psisq
            end if
        end do
    end do
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
