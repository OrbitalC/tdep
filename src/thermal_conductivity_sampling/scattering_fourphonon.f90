

subroutine compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg, rng, mw, mem)
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
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> Four phonon prefactor
    ! real(r8), parameter :: fourphonon_prefactor = lo_pi / 8.0_r8
    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3, egv4
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2, evp3
    !> The full qpoint grid, to be shuffled
    integer, dimension(:), allocatable :: qgridfull1, qgridfull2
    !> The qpoints in cartesian coordinates
    real(r8), dimension(3) :: qv2, qv3, qv4
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> The ratio for the maximum likelihood estimation of the scattering rates
    real(r8) :: mle_ratio
    !> Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, psisq, prefactor
    ! The gaussian integration width
    real(r8) :: sig1, sig2, sig3, sig4, sigma
    !> Integers for do loops
    integer :: i, q1, q2, q3, q4, b1, b2, b3, b4, qi, qj
    !> Is the quartet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult

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
    sig1 = sr%sigma_q(q1, b1)

    ! Prepare the grid for the monte-carlo average
    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull1, rng)
    call mcg%generate_grid(qgridfull2, rng)

    i = 0
    count_loop: do qi=1, mcg%npoints
    do qj=1, mcg%npoints
        q2 = qgridfull1(qi)
        q3 = qgridfull2(qj)
        if (q3 .lt. q2) cycle  ! Permutation symmetry
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, mcg%full_dims)
        if (q4 .lt. q3) cycle  ! Permutation symmetry

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, mult)
        if (isred) cycle

        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle
            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)

            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)

                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    sig4 = sr%sigma_q(qp%ap(q4)%irreducible_index, b4)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)
                    if (abs(om1 + om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 + om4 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 + om4 - om2) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 - om2 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om4 - om2 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. 4.0_r8 * sigma) i = i + 1
                end do
            end do
        end do
    end do
    end do count_loop

    sr%fourphonon(il)%n = i
    allocate(sr%fourphonon(il)%psisq(i))
    allocate(sr%fourphonon(il)%q2(i))
    allocate(sr%fourphonon(il)%q3(i))
    allocate(sr%fourphonon(il)%q4(i))
    allocate(sr%fourphonon(il)%b2(i))
    allocate(sr%fourphonon(il)%b3(i))
    allocate(sr%fourphonon(il)%b4(i))
    sr%fourphonon(il)%psisq = 0.0_r8

    i = 0
    full_loop: do qi=1, mcg%npoints
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
        call pretransform_phi4(fcf, qv2, qv3, qv4, ptf)

        prefactor = fourphonon_prefactor * mcg%weight**2
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            sig2 = sr%sigma_q(qp%ap(q2)%irreducible_index, b2)
            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)
                sig3 = sr%sigma_q(qp%ap(q3)%irreducible_index, b3)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4=1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    sig4 = sr%sigma_q(qp%ap(q4)%irreducible_index, b4)
                    sigma = sqrt(sig1**2 + sig2**2 + sig3**2 + sig4**2)
                    if (abs(om1 + om2 + om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 + om4 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 + om4 - om2) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om2 - om3 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om3 - om2 - om4) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 + om4 - om2 - om3) .lt. 4.0_r8 * sigma .or. &
                        abs(om1 - om2 - om3 - om4) .lt. 4.0_r8 * sigma) then
                        i = i + 1

                        egv4 = dr%aq(q4)%egv(:, b4) / sqrt(om4)

                        evp3 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                        evp3 = conjg(evp3)
                        c0 = dot_product(evp3, ptf)
                        psisq = abs(c0*conjg(c0)) * prefactor

                        sr%fourphonon(il)%psisq(i) = psisq
                        sr%fourphonon(il)%q2(i) = q2
                        sr%fourphonon(il)%q3(i) = q3
                        sr%fourphonon(il)%q4(i) = q4
                        sr%fourphonon(il)%b2(i) = b2
                        sr%fourphonon(il)%b3(i) = b3
                        sr%fourphonon(il)%b4(i) = b4

                    end if
                end do
            end do
        end do
    end do
    end do full_loop
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
end function

!> Get the Fourier transform of the fourth order matrix element
subroutine pretransform_phi4(fcf, q2, q3, q4, ptf)
    !> third order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q2, q3, q4
    !> flattened, pretransformed matrix element
    complex(r8), dimension(:), intent(out) :: ptf

    integer :: i, j, k, l, m

    complex(r8) :: expiqr
    real(r8), dimension(3) :: rv2, rv3, rv4
    real(r8) :: iqr
    integer :: a1, a2, a3, a4, ia, ib, ic, id, q, nb

    nb = fcf%na*3
    ptf = 0.0_r8
    do a1 = 1, fcf%na
    do q = 1, fcf%atom(a1)%n
        a2 = fcf%atom(a1)%quartet(q)%i2
        a3 = fcf%atom(a1)%quartet(q)%i3
        a4 = fcf%atom(a1)%quartet(q)%i4

        rv2 = fcf%atom(a1)%quartet(q)%lv2
        rv3 = fcf%atom(a1)%quartet(q)%lv3
        rv4 = fcf%atom(a1)%quartet(q)%lv4

        iqr = dot_product(q2, rv2) + dot_product(q3, rv3) + dot_product(q4, rv4)
        iqr = -iqr*lo_twopi
        expiqr = cmplx(cos(iqr), sin(iqr), r8)
        do l = 1, 3
        do k = 1, 3
        do j = 1, 3
        do i = 1, 3
            ia = (a1 - 1)*3 + i
            ib = (a2 - 1)*3 + j
            ic = (a3 - 1)*3 + k
            id = (a4 - 1)*3 + l
            ! Now for the grand flattening scheme, consistent with the zgeru operations above.
            m = (ia - 1)*nb*nb*nb + (ib - 1)*nb*nb + (ic - 1)*nb + id
            ptf(m) = ptf(m) + fcf%atom(a1)%quartet(q)%mwm(i, j, k, l)*expiqr
        end do
        end do
        end do
        end do
    end do
    end do
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
    mult = qp%ip(q1)%n_invariant_operation / mult
end subroutine
