
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

    ! g0 = g0 - 2.0 * sqerr / uc%mass(1) / lo_kB_Hartree / temperature
    g0 = g0 - 2.0 * sqerr / uc%mass(1) / om1

end subroutine
