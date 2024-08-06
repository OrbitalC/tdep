

subroutine compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg, rng, thres, g0, mw, mem)
    !> The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    ! The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Fourth order force constants
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> The monte-carlo grid
    type(lo_montecarlo_grid), intent(in) :: mcg
    !> The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The threshold for gaussian integration
    real(r8), intent(in) :: thres
    !> The linewidth for this mode
    real(r8), intent(inout) :: g0
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2, evp3
    !> The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    !> The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4, veldiff
    !> The inverse reciprocal lattice vector
    real(r8), dimension(3, 3) :: gvec
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, psisq, prefactor
    ! The gaussian integration width
    real(r8) :: sigma, pref_sigma
    !> Stuff for the linewidths
    real(r8) :: n2, n3, n4, n2p, n3p, n4p, plf1, plf2, plf3, plf4, plf5, plf6, plf7, n1
    !> Integers for do loops
    integer :: i, q1, q2, q3, q4, b1, b2, b3, b4, qi, qj, i2, i3, i4, q1f
    !> Is the quartet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult
    !> All the prefactors for the scattering
    real(r8) :: f0, f1, f2, f3, f4, f5, f6, f7, f8

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv4, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
  ! pref_sigma = qp%ip(1)%radius * lo_twopi / sqrt(2.0_r8)
    pref_sigma = qp%ip(1)%radius
    q1f = qp%ip(q1)%full_index

    gvec(1, :) = uc%reciprocal_latticevectors(:, 1) / mcg%full_dims(1)
    gvec(2, :) = uc%reciprocal_latticevectors(:, 2) / mcg%full_dims(2)
    gvec(3, :) = uc%reciprocal_latticevectors(:, 3) / mcg%full_dims(3)

    n1 = sr%be(q1, b1)

    ! Prepare the grid for the monte-carlo average
    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull1, rng)
    call mcg%generate_grid(qgridfull2, rng)

    do qi=1, mcg%npoints
    do qj=1, mcg%npoints
        q2 = qgridfull1(qi)
        q3 = qgridfull2(qj)
        if (q3 .lt. q2) cycle
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, mcg%full_dims)
        if (q4 .lt. q3) cycle

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call fcf%pretransform(qv2, qv3, qv4, ptf)

        prefactor = fourphonon_prefactor * mcg%weight**2 * mult
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
            n2p = n2 + 1.0_r8
            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)
                n3p = n3 + 1.0_r8
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    n4 = sr%be(qp%ap(q4)%irreducible_index, b4)
                    n4p = n4 + 1.0_r8

                    ! if (b3 .eq. b4 .and. q3 .eq. q4) cycle
                    ! if (b3 .eq. b4) cycle

                    veldiff = dr%aq(q3)%vel(:, b3) - dr%aq(q4)%vel(:, b4)
                    sigma = max(norm2(veldiff * gvec(1, :)), &
                                norm2(veldiff * gvec(2, :)), &
                                norm2(veldiff * gvec(3, :)))

                    if (sigma .lt. lo_freqtol) cycle

                 !  sigma = norm2(dr%aq(q3)%vel(:, b3) - dr%aq(q4)%vel(:, b4)) * pref_sigma
                 !  sigma = max(0.25_r8 * dr%default_smearing(b3), 0.25_r8 * dr%default_smearing(b4), sigma)
                 !  sigma = min(4.0_r8 * dr%default_smearing(b3), 4.0_r8 * dr%default_smearing(b4), sigma)

                    if (abs(om1 + om2 + om3 + om4) .lt. thres * sigma .or. &
                        abs(om1 - om2 - om4 - om3) .lt. thres * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. thres * sigma .or. &
                        abs(om1 - om2 + om3 + om4) .lt. thres * sigma .or. &
                        abs(om1 + om3 - om2 - om4) .lt. thres * sigma .or. &
                        abs(om1 - om3 + om2 + om4) .lt. thres * sigma .or. &
                        abs(om1 + om4 - om2 - om3) .lt. thres * sigma .or. &
                        abs(om1 - om4 + om2 + om3) .lt. thres) then

                        egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)

                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                        evp3 = conjg(evp3)
                        c0 = dot_product(evp3, ptf)
                        psisq = abs(c0*conjg(c0)) * prefactor

                        ! Flat indices for the scattering matrix
                        i2 = (q2 - 1) * dr%n_mode + b2
                        i3 = (q3 - 1) * dr%n_mode + b3
                        i4 = (q4 - 1) * dr%n_mode + b4

                        ! Prefactors, only the Bose-Einstein distributions
                        plf1 = n2p * n3p * n4p - n2 * n3 * n4
                        plf2 = 3.0_r8 * n2 * n3p * n4p - n2p * n3 * n4
                        plf3 = 3.0_r8 * n3 * n2p * n4p - n3p * n2 * n4
                        plf4 = 3.0_r8 * n4 * n3p * n2p - n4p * n3 * n2

                        ! Prefactors, including the matrix elements and dirac
                        f0 = psisq * plf1 * lo_gauss(om1, om2 + om3 + om4, sigma)
                        f1 = psisq * plf1 * lo_gauss(om1,-om2 - om3 - om4, sigma)
                        f2 = psisq * plf2 * lo_gauss(om1,-om2 + om3 + om4, sigma)
                        f3 = psisq * plf2 * lo_gauss(om1, om2 - om3 - om4, sigma)
                        f4 = psisq * plf3 * lo_gauss(om1,-om3 + om2 + om4, sigma)
                        f5 = psisq * plf3 * lo_gauss(om1, om3 - om2 - om4, sigma)
                        f6 = psisq * plf4 * lo_gauss(om1,-om4 + om3 + om2, sigma)
                        f7 = psisq * plf4 * lo_gauss(om1, om4 - om3 - om2, sigma)

                        ! Add everything to the linewidth, the pref comes from permutations
                        g0 = g0 + 6.0_r8 * (f0 - f1) + 2.0_r8 * (f2 - f3 + f4 - f5 + f6 - f7)

                        ! And then to the scattering matrix
                        if (q1f .ne. q2 .or. b1 .ne. b2) then
                            sr%Xi(il, i2) = sr%Xi(il, i2) + 4.0_r8 * (f0 - f1 + f2 - f3) * om2 / om1
                           !sr%Xi(il, i2) = sr%Xi(il, i2) + (12.0_r8 * (f0 - f1) + 4.0_r8 * (f2 - f3)) * om2 / om1
                           !sr%Xi(il, i2) = sr%Xi(il, i2) + 2.0_r8 * (6.0_r8 * (f0 - f1) + 2.0_r8 * (f2 - f3 + f4 - f5 + f6 - f7)) * om2 / om1
                        end if
                        if (q1f .ne. q3 .or. b1 .ne. b3) then
                            sr%Xi(il, i3) = sr%Xi(il, i3) + 2.0_r8 * (f0 - f1 + f4 - f5) * om3 / om1
                           !sr%Xi(il, i3) = sr%Xi(il, i3) + (12.0_r8 * (f0 - f1) + 4.0_r8 * (f4 - f5)) * om3 / om1
                           !sr%Xi(il, i3) = sr%Xi(il, i3) + 2.0_r8 * (6.0_r8 * (f0 - f1) + 2.0_r8 * (f2 - f3 + f4 - f5 + f6 - f7)) * om3 / om1
                        end if
                        if (q1f .ne. q4 .or. b1 .ne. b4) then
                            sr%Xi(il, i4) = sr%Xi(il, i4) + 2.0_r8 * (f0 - f1 + f6 - f7) * om4 / om1
                           !sr%Xi(il, i4) = sr%Xi(il, i4) + (12.0_r8 * (f0 - f1) + 4.0_r8 * (f6 - f7)) * om4 / om1
                           !sr%Xi(il, i4) = sr%Xi(il, i4) + 2.0_r8 * (6.0_r8 * (f0 - f1) + 2.0_r8 * (f2 - f3 + f4 - f5 + f6 - f7)) * om4 / om1
                        end if

                     ! Those are formulas directly from the FourPhonon paper
                     !  psisq = psisq * 12.0_r8  ! Different way to normalize the matrix elements
                     !  plf1 = n2p * n3p  * n4 / n1 * psisq * lo_gauss(om1, -om2 - om3 + om4, sigma)
                     !  plf2 = n2p * n4p  * n3 / n1 * psisq * lo_gauss(om1, -om2 - om4 + om3, sigma)
                     !  plf3 = n3p * n4p  * n2 / n1 * psisq * lo_gauss(om1, -om3 - om4 + om2, sigma)
                     !  plf4 = n2p * n3 * n4 / n1 * psisq * lo_gauss(om1, -om2 + om3 + om4, sigma)
                     !  plf5 = n3p * n2 * n4 / n1 * psisq * lo_gauss(om1, -om3 + om2 + om4, sigma)
                     !  plf6 = n4p * n2 * n3 / n1 * psisq * lo_gauss(om1, -om4 + om2 + om3, sigma)
                     !  plf7 = n2 * n3 * n4 / n1 * psisq * lo_gauss(om1, om2 + om3 + om4, sigma)

                     !  ! Prefactors (or lack of) comes from permutations
                     !  g0 = g0 + (plf1 + plf2 + plf3) + (plf4 + plf5 + plf6) + plf7

                     !  i2 = (q2 - 1) * dr%n_mode + b2
                     !  i3 = (q3 - 1) * dr%n_mode + b3
                     !  i4 = (q4 - 1) * dr%n_mode + b4

                     !  f0 = 0.5_r8 * (-plf1 - plf2 + plf3 - plf4 + plf5 + plf6) + plf7 / 6.0_r8
                     !  sr%Xi(il, i2) = sr%Xi(il, i2) + f0 * om2 / om1
                     !  f0 = 0.5_r8 * (-plf1 + plf2 - plf3 + plf4 - plf5 + plf6) + plf7 / 6.0_r8
                     !  sr%Xi(il, i3) = sr%Xi(il, i3) + f0 * om3 / om1
                     !  f0 = 0.5_r8 * (plf1 - plf2 - plf3 + plf4 + plf5 - plf6) + plf7 / 6.0_r8
                     !  sr%Xi(il, i4) = sr%Xi(il, i4) + f0 * om4 / om1

                       !sr%Xi(il, i2) = sr%Xi(il, i2) + (0.5_r8 * plf1 + 0.5_r8 * plf4 + plf7 / 6.0_r8) * om2 / om1
                       !sr%Xi(il, i3) = sr%Xi(il, i3) + (0.5_r8 * plf2 + 0.5_r8 * plf5 + plf7 / 6.0_r8) * om3 / om1
                       !sr%Xi(il, i4) = sr%Xi(il, i4) + (0.5_r8 * plf3 + 0.5_r8 * plf6 + plf7 / 6.0_r8) * om4 / om1
                    end if
                end do
            end do
        end do
    end do
    end do
    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

subroutine quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2, q3 and q4 are full
    integer, intent(in) :: q1, q2, q3, q4
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> If it's reducible, what is its multiplicity
    real(r8), intent(out) :: mult

    ! To hold the q-point in reduced coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4, qv2p, qv3p, qv4p
    !> To get the index of the new triplet on the fft_grid
    integer, dimension(3) :: gi
    !> The new triplet after the operation
    integer, dimension(3) :: qpp
    !> Integers for the loops
    integer :: j, k

    ! First get the reciprocal lattice vectors, in reduce coordinates
    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
    qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)
    qv4 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q4)%r)

    isred = .false.
    mult = 0.0_r8
    ! Let's try all operations that leaves q1 invariant
    do j=1, qp%ip(q1)%n_invariant_operation
        k = qp%ip(q1)%invariant_operation(j)
        qpp = -lo_hugeint
        select type(qp); type is(lo_fft_mesh)
            ! Transform q2 and check if it's on the grid
            qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle

            ! Transform q3 and check if it's on the grid
            qv3p = lo_operate_on_vector(uc%sym%op(k), qv3, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv3p) .eqv. .false.) cycle

            ! Transform q4 and check if it's on the grid
            qv4p = lo_operate_on_vector(uc%sym%op(k), qv4, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv4p) .eqv. .false.) cycle

            ! If everything is on the grid, get the location of each point
            gi = qp%index_from_coordinate(qv2p)
            qpp(1) = qp%gridind2ind(gi(1), gi(2), gi(3))
            gi = qp%index_from_coordinate(qv3p)
            qpp(2) = qp%gridind2ind(gi(1), gi(2), gi(3))
            gi = qp%index_from_coordinate(qv4p)
            qpp(3) = qp%gridind2ind(gi(1), gi(2), gi(3))
        end select
        ! The sorting allows to include permutation invariance
        call lo_qsort(qpp)
        if (qpp(1) .gt. q2) then
            isred = .true.
        ! For the weights, it's the same idea that for the triplet
        ! I compute the number of operations that let the quartet invariant
        ! And take the ratio between the little group of q1 and this number
        else if (qpp(1) .eq. q2 .and. qpp(2) .eq. q3 .and. qpp(3) .eq. q4) then
            mult = mult + 1.0_r8
        end if
    end do
    mult = qp%ip(q1)%n_invariant_operation * 1.0_r8 / mult
end subroutine
