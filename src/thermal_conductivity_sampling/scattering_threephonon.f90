

subroutine compute_threephonon_scattering(il, sr, qp, dr, uc, fct, mcg, rng, g0, mw, mem)
    ! The qpoint and mode indices considered here
    integer, intent(in) :: il
    !> The scattering amplitudes
    type(lo_scattering_rates), intent(inout) :: sr
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    ! Harmonic dispersions
    type(lo_phonon_dispersions), intent(inout) :: dr
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> Third order force constants
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> The monte-carlo grid
    type(lo_montecarlo_grid), intent(in) :: mcg
    !> The random number generator
    type(lo_mersennetwister), intent(inout) :: rng
    !> The linewidth for this mode
    real(r8), intent(inout) :: g0
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2
    !> The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull
    !> The qpoints and the dimension of the qgrid
    real(r8), dimension(3) :: qv2, qv3
    !> The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sigma
    !> Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, plf, psisq, prefactor, mle_ratio
    !> The bose-einstein distribution for the modes
    real(r8) :: n2, n3, n2p, n3p
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops
    integer :: i, qi, q1, q2, q3, b1, b2, b3
    !> Is the triplet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Already set some values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1) / sqrt(om1)
    sig1 = sr%sigma_q(q1, b1)

    call mem%allocate(qgridfull, mcg%npoints, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull, rng)

    i = 0
    count_loop: do qi=1, mcg%npoints
        q2 = qgridfull(qi)
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, mcg%full_dims)
        if (q3 .lt. q2) cycle

        call triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult)
        if (isred) cycle

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
               !sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                sigma = sqrt(sig2**2 + sig3**2)
                if (abs(om1 + om2 - om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 + om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 - om3) .lt. 4.0_r8 * sigma) i = i + 1
            end do
        end do
    end do count_loop

    sr%threephonon(il)%n = i
    allocate(sr%threephonon(il)%psisq(i))
    allocate(sr%threephonon(il)%q2(i))
    allocate(sr%threephonon(il)%q3(i))
    allocate(sr%threephonon(il)%b2(i))
    allocate(sr%threephonon(il)%b3(i))

    i = 0
    compute_loop: do qi=1, mcg%npoints
        q2 = qgridfull(qi)
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, mcg%full_dims)
        if (q3 .lt. q2) cycle

        call triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        call pretransform_phi3(fct, qv2, qv3, ptf)
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)
            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
            prefactor = threephonon_prefactor * mcg%weight * mult

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)

                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)
               !sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                sigma = sqrt(sig2**2 + sig3**2)

                ! Do we need to compute the scattering ?
                if (abs(om1 + om2 - om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 + om3) .lt. 4.0_r8 * sigma .or. &
                    abs(om1 - om2 - om3) .lt. 4.0_r8 * sigma) then
                    i = i + 1

                    evp2 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                    evp2 = conjg(evp2)
                    c0 = dot_product(evp2, ptf)
                    psisq = abs(c0*conjg(c0)) * prefactor

                    sr%threephonon(il)%psisq(i) = psisq
                    sr%threephonon(il)%q2(i) = q2
                    sr%threephonon(il)%q3(i) = q3
                    sr%threephonon(il)%b2(i) = b2
                    sr%threephonon(il)%b3(i) = b3

                    n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
                    n3 = sr%be(qp%ap(q3)%irreducible_index, b3)
                    n2p = n2 + 1.0_r8
                    n3p = n3 + 1.0_r8

                    g0 = g0 + 2.0_r8 * psisq * (n2 + n3 + 1.0_r8) * lo_gauss(om1, om2 + om3, sigma)
                    g0 = g0 - 2.0_r8 * psisq * (n2 + n3 + 1.0_r8) * lo_gauss(om1, -om2 - om3, sigma)

                    g0 = g0 + 2.0_r8 * psisq * (n2 - n3) * lo_gauss(om1, -om2 + om3, sigma)
                    g0 = g0 - 2.0_r8 * psisq * (n2 - n3) * lo_gauss(om1,  om2 - om3, sigma)
                end if
            end do
        end do
    end do compute_loop
    ! And we can deallocate everything
    call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(egv3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(qgridfull, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine


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
end function

!> Get the Fourier transform of the third order matrix element
subroutine pretransform_phi3(fct, q2, q3, ptf)
    !> third order forceconstant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3
    real(r8) :: iqr
    integer :: a1, a2, a3, ia, ib, ic, t, nb

    nb = fct%na*3
    ptf = 0.0_r8
    do a1 = 1, fct%na
    do t = 1, fct%atom(a1)%n
        a2 = fct%atom(a1)%triplet(t)%i2
        a3 = fct%atom(a1)%triplet(t)%i3

        rv2 = fct%atom(a1)%triplet(t)%lv2
        rv3 = fct%atom(a1)%triplet(t)%lv3

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do i = 1, 3
        do j = 1, 3
        do k = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            l = (ia - 1)*nb*nb + (ib - 1)*nb + ic
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%mwm(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
end subroutine


subroutine triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2 and q3 are full
    integer, intent(in) :: q1, q2, q3
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> If it's reducible, what is its multiplicity
    real(r8), intent(out) :: mult

    ! To hold the q-point in reduced coordinates
    real(r8), dimension(3) :: qv2, qv3, qv2p, qv3p
    !> To get the index of the new triplet on the fft_grid
    integer, dimension(3) :: gi
    !> The new triplet after the operation
    integer, dimension(2) :: qpp
    !> Integers for the loops
    integer :: j, k

    ! First get the q-points in reduce coordinates
    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
    qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)

    mult = 0.0_r8
    isred = .false.
    ! Let's try all operations that leaves q1 invariant
    do j=1, qp%ip(q1)%n_invariant_operation
        k = qp%ip(q1)%invariant_operation(j)
        qpp = -lo_hugeint
        select type(qp); type is(lo_fft_mesh)
            ! Rotate q2 and look if it's the on grid
            qv2p = lo_operate_on_vector(uc%sym%op(k), qv2, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv2p) .eqv. .false.) cycle

            ! Rotate q3 and look if it's the on grid
            qv3p = lo_operate_on_vector(uc%sym%op(k), qv3, reciprocal=.true., fractional=.true.)
            if (qp%is_point_on_grid(qv3p) .eqv. .false.) cycle

            ! If everything is on the grid, get the index of each point
            gi = qp%index_from_coordinate(qv2p)
            qpp(1) = qp%gridind2ind(gi(1), gi(2), gi(3))
            gi = qp%index_from_coordinate(qv3p)
            qpp(2) = qp%gridind2ind(gi(1), gi(2), gi(3))
        end select
        if (minval(qpp) .gt. q2) then
            isred = .true.
        ! Now we have to determine the weight
        ! Two roads are possible
        !      1. Get the ratio of number of red point that can give this irreducible point
        !      2. Look at the ratio between total number of operations and the ones
        !         that leaves this irreducible triplet unchanged
        ! The second road doesn't requires me to sum over all other qpoints, so I'll go with this one
        else if (minval(qpp) .eq. q2 .and. maxval(qpp) .eq. q3) then
            mult = mult + 1.0_r8
        end if
    end do
    mult = qp%ip(q1)%n_invariant_operation * 1.0_r8 / mult
end subroutine
