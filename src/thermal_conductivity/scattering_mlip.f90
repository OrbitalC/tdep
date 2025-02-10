
subroutine compute_mlip_scattering(il, sr, qp, dr, uc, temperature, sqerr, &
                                   g0, integrationtype, smearing, mw, mem)
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
    !> The variance of the forces
    real(r8), intent(in) :: sqerr
    !> The linewidth for this mode
    real(r8), intent(inout) :: g0
    !> what kind of integration are we doing
    integer, intent(in) :: integrationtype
    !> The smearing width
    real(r8), intent(in) :: smearing
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    ! Eigenvectors
    complex(r8), dimension(uc%na*3, 2) :: egviso
    ! prefactor and phonon buffers
    real(r8) :: om1, om2, sigma, psisq, prefactor, f0
    ! Integers for do loops
    integer :: q1, b1, q2, b2, i, niso

    q1 = sr%my_qpoints(il)
    b1 = sr%my_modes(il)
    om1 = dr%iq(q1)%omega(b1)
    egviso(:, 1) = dr%iq(q1)%egv(:, b1)

    do q2 = 1, qp%n_full_point
        prefactor = isotope_prefactor*qp%ap(q2)%integration_weight
        do b2 = 1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            select case (integrationtype)
            case (1)
                sigma = lo_frequency_THz_to_Hartree*smearing
            case (2)
                sigma = sqrt(sr%sigsq(q1, b1) + &
                             sr%sigsq(qp%ap(q2)%irreducible_index, b2))
            case (6)
                sigma = qp%smearingparameter(dr%aq(q2)%vel(:, b2), &
                                             dr%default_smearing(b2), smearing)
            end select

            i = (q2 - 1)*dr%n_mode + b2

            egviso(:, 2) = dr%aq(q2)%egv(:, b2)

            psisq = sqerr*prefactor

            f0 = psisq*om1*om2*lo_gauss(om1, om2, sigma)
            g0 = g0 - f0
!           sr%Xi(il, i) = sr%Xi(il, i) + f0*om2/om1
        end do
    end do

end subroutine
