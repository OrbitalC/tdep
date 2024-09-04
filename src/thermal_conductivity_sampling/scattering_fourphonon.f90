
subroutine compute_fourphonon_scattering(il, sr, qp, dr, uc, fcf, mcg, rng, thres, &
                                         g0, integrationtype, smearing, mctol, mw, mem)
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
    !> what kind of integration are we doing
    integer, intent(in) :: integrationtype
    !> The smearing width
    real(r8), intent(in) :: smearing
    !> The Monte-Carlo integration tolerance
    real(r8), intent(in) :: mctol
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
    real(r8), dimension(3) :: qv2, qv3, qv4
    !> The complex scattering amplitude
    complex(r8) :: c0
    !> Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, om4, psisq, prefactor
    ! The gaussian integration width
    real(r8) :: sigma
    !> Stuff for the linewidths
    real(r8) :: n2, n3, n4, n2p, n3p, n4p, plf1, plf2, plf3, plf4, plf5, plf6, plf7
    !> Integers for do loops
    integer :: q1, q2, q3, q4, b1, b2, b3, b4, qi, qj, i2, i3, i4, i
    !> Is the quartet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    real(r8) :: mult1, mult2, mult3, mult4
    !> All the prefactors for the scattering
    real(r8) :: f0, f1, f2, f3, f4, f5, f6, f7, fall
    !> The reducible triplet corresponding to the currently computed quartet
    integer, dimension(:, :), allocatable :: red_quartet

    real(r8), dimension(:, :), allocatable :: od_terms
!   integer :: n, m
    integer :: m
    real(r8) :: buf
    real(r8), dimension(10) :: buf_iter
    real(r8), dimension(9) :: diff_iter

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv4, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    call mem%allocate(od_terms, [qp%n_full_point, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    od_terms = 0.0_r8

    ! Already set some buffer values for mode (q1, b1)
    q1 = sr%q1(il)
    b1 = sr%b1(il)
    om1 = dr%iq(q1)%omega(b1)
    egv1 = dr%iq(q1)%egv(:, b1)/sqrt(om1)

    ! Prepare the grid for the monte-carlo average
    call mem%allocate(qgridfull1, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(qgridfull2, qp%n_full_point, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull1, rng)
    call mcg%generate_grid(qgridfull2, rng)

    call rng%shuffle_int_array(qgridfull1)
    call rng%shuffle_int_array(qgridfull2)

!   n = 0
    m = 0
    buf = 0.0_r8
    buf_iter = 0.0_r8

    compute_loop: do qi = 1, mcg%npoints
    do qj = 1, mcg%npoints
        if (mctol .gt. 0) then
            q2 = rng%rnd_int(mcg%npoints)
            q3 = rng%rnd_int(mcg%npoints)
        else
            q2 = qgridfull1(qi)
            q3 = qgridfull2(qj)
        end if
        if (q3 .lt. q2) cycle
        q4 = fft_fourth_grid_index(qp%ip(q1)%full_index, q2, q3, mcg%full_dims)
        if (q4 .lt. q3) cycle

        call quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, red_quartet, mw, mem)
        if (isred) cycle

        qv2 = qp%ap(q2)%r
        qv3 = qp%ap(q3)%r
        qv4 = qp%ap(q4)%r
        call fcf%pretransform(qv2, qv3, qv4, ptf)

        ! We can already take care of the multiplicity caused by permutation
        if (q2 .eq. q3 .and. q3 .eq. q4) then
            mult1 = 1.0_r8
            mult2 = 1.0_r8
            mult3 = 1.0_r8
            mult4 = 1.0_r8
        else if ((q2 .ne. q3 .and. q3 .eq. q4) .or. &
                 (q3 .ne. q2 .and. q2 .eq. q4) .or. &
                 (q4 .ne. q2 .and. q2 .eq. q3)) then
            mult1 = 3.0_r8
            mult2 = 1.0_r8
            mult3 = 1.0_r8
            mult4 = 1.0_r8
        else
            mult1 = 6.0_r8
            mult2 = 2.0_r8
            mult3 = 2.0_r8
            mult4 = 2.0_r8
        end if
!       mult1 = 1.0_r8
!       mult2 = 1.0_r8
!       mult3 = 1.0_r8
!       mult4 = 1.0_r8

       !prefactor = fourphonon_prefactor*mcg%weight**2
        prefactor = fourphonon_prefactor
        do b2 = 1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
            n2p = n2 + 1.0_r8
            egv2 = dr%aq(q2)%egv(:, b2)/sqrt(om2)

            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3 = 1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle

                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)
                n3p = n3 + 1.0_r8
                egv3 = dr%aq(q3)%egv(:, b3)/sqrt(om3)

                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                do b4 = 1, dr%n_mode
                    om4 = dr%aq(q4)%omega(b4)
                    if (om4 .lt. lo_freqtol) cycle

                    n4 = sr%be(qp%ap(q4)%irreducible_index, b4)
                    n4p = n4 + 1.0_r8

                    select case (integrationtype)
                    case (1)
                        sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing
                    case (2)
                        sigma = sqrt(sr%sigsq(q1, b1) + &
                                     sr%sigsq(qp%ap(q2)%irreducible_index, b2) + &
                                     sr%sigsq(qp%ap(q3)%irreducible_index, b3) + &
                                     sr%sigsq(qp%ap(q4)%irreducible_index, b4))
                    end select

                    egv4 = dr%aq(q4)%egv(:, b4)/sqrt(om4)

                    evp3 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**3, (1.0_r8, 0.0_r8), egv4, 1, evp2, 1, evp3, dr%n_mode)
                    evp3 = conjg(evp3)
                    c0 = dot_product(evp3, ptf)
                    psisq = abs(c0*conjg(c0))*prefactor

                    ! Prefactors, only the Bose-Einstein distributions
                    plf1 = n2p*n3p*n4p - n2*n3*n4
                    plf2 = 3.0_r8*n2*n3p*n4p - n2p*n3*n4
                    plf3 = 3.0_r8*n3*n2p*n4p - n3p*n2*n4
                    plf4 = 3.0_r8*n4*n3p*n2p - n4p*n3*n2

                    ! Prefactors, including the matrix elements and dirac
                    f0 = mult1*psisq*plf1*lo_gauss(om1, om2 + om3 + om4, sigma)
                    f1 = mult1*psisq*plf1*lo_gauss(om1, -om2 - om3 - om4, sigma)
                    f2 = mult2*psisq*plf2*lo_gauss(om1, -om2 + om3 + om4, sigma)
                    f3 = mult2*psisq*plf2*lo_gauss(om1, om2 - om3 - om4, sigma)
                    f4 = mult3*psisq*plf3*lo_gauss(om1, -om3 + om2 + om4, sigma)
                    f5 = mult3*psisq*plf3*lo_gauss(om1, om3 - om2 - om4, sigma)
                    f6 = mult4*psisq*plf4*lo_gauss(om1, -om4 + om3 + om2, sigma)
                    f7 = mult4*psisq*plf4*lo_gauss(om1, om4 - om3 - om2, sigma)
                    fall = f0 - f1 + f2 - f3 + f4 - f5 + f6 - f7

                    ! Add everything to the linewidth
                    do i = 1, size(red_quartet, 2)
                       !g0 = g0 + fall
                        buf = buf + fall

                        od_terms(red_quartet(1, i), b2) = od_terms(red_quartet(1, i), b2) + 2.0_r8 * fall * om2 / om1
                        od_terms(red_quartet(2, i), b3) = od_terms(red_quartet(2, i), b3) + 2.0_r8 * fall * om3 / om1
                        od_terms(red_quartet(3, i), b4) = od_terms(red_quartet(3, i), b4) + 2.0_r8 * fall * om4 / om1

                       !i2 = (red_quartet(1, i) - 1)*dr%n_mode + b2
                       !sr%Xi(il, i2) = sr%Xi(il, i2) + 2.0_r8*fall*om2/om1

                       !i3 = (red_quartet(2, i) - 1)*dr%n_mode + b3
                       !sr%Xi(il, i3) = sr%Xi(il, i3) + 2.0_r8*fall*om3/om1

                       !i4 = (red_quartet(3, i) - 1)*dr%n_mode + b4
                       !sr%Xi(il, i4) = sr%Xi(il, i4) + 2.0_r8*fall*om4/om1
                    end do
                end do
            end do
        end do
!       n = n + 1
        m = m + size(red_quartet, 2) * (mult1 + mult2 + mult3 + mult4)
!       buf_iter(n) = buf / real(m, r8)
        do i=2, 10
            buf_iter(i-1) = buf_iter(i)
        end do
        buf_iter(10) = buf / real(m, r8)
!       if (n .gt. 11) then
            do i=1, 9
                diff_iter(i) = abs(buf_iter(10) - buf_iter(10-i)) / abs(buf_iter(10))
            end do
            if (maxval(diff_iter) .lt. mctol) then
                exit compute_loop
            end if
!       end if
    end do
    end do compute_loop

    ! And now we add things, with the normalization
    g0 = g0 + buf / real(m, r8)
    od_terms = od_terms / real(m, r8)
    do q2=1, qp%n_full_point
        do b2=1, dr%n_mode
            i2 = (q2 - 1)*dr%n_mode + b2
            sr%Xi(il, i2) = sr%Xi(il, i2) + od_terms(q2, b2)
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
    call mem%deallocate(od_terms, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine

subroutine quartet_is_irreducible(qp, uc, q1, q2, q3, q4, isred, red_quartet, mw, mem)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2, q3 and q4 are full
    integer, intent(in) :: q1, q2, q3, q4
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> The equivalent triplet
    integer, dimension(:, :), allocatable, intent(out) :: red_quartet
    !> Mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem

    !> The new-qpoints and the temporary invariant triplet
    integer, dimension(:, :), allocatable :: newqp, newqp_sort
    !> The number of invariant operation for q1
    integer :: n

    n = qp%ip(q1)%n_invariant_operation

    ! We will need this to keep equivalent point
    call mem%allocate(newqp, [3, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(newqp_sort, [3, n], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    newqp = -lo_hugeint
    newqp_sort = -lo_hugeint

    ! The first part is to generate all qpoint invariant in little star of q1
    get_equivalent: block
        !> To hold the q-point in reduced coordinates
        real(r8), dimension(3) :: qv2, qv3, qv4, qv2p, qv3p, qv4p
        !> To get the index of the new triplet on the fft_grid
        integer, dimension(3) :: gi
        !> The new triplet after the operation
        integer, dimension(3) :: qpp
        !> Integers for the do loops
        integer :: j, k

        ! First get the reciprocal lattice vectors, in reduce coordinates
        qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
        qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)
        qv4 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q4)%r)

        isred = .false.
        ! Let's try all operations that leaves q1 invariant
        do j = 1, n
            k = qp%ip(q1)%invariant_operation(j)
            qpp = -lo_hugeint
            select type (qp); type is (lo_fft_mesh)
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
            ! We need to keep the symmetry equivalent point in the order they came in
            newqp(:, j) = qpp
            ! We also need them sorted, to check for redundancy and reducibility
            call lo_qsort(qpp)
            newqp_sort(:, j) = qpp
            if (qpp(1) .gt. q2 .or. qpp(2) .gt. q3) isred = .true.
        end do

        if (minval(newqp) .lt. 0) then
            do j = 1, size(newqp, 2)
                if (any(newqp(:, j) .lt. 0)) newqp(:, j) = [q2, q3, q4]
                if (any(newqp(:, j) .lt. 0)) newqp_sort(:, j) = [q2, q3, q4]
            end do
        end if
    end block get_equivalent

    ! Then we just get rid of redundant point in the list
    sort_reducible: block
        !> Integers for the loops
        integer :: j, k, ctr, ctr2

        ! Now we have the same problem of permutation as in the third order
        ! So first, we count the quartet equivalent by permutation, a little bit harder
        ctr = 0
        do j = 1, n
            ctr2 = 0
            do k = j, n
                if (k .eq. j) cycle
                if (all(newqp_sort(:, j) .eq. newqp_sort(:, k))) ctr2 = ctr2 + 1
            end do
            if (ctr2 .eq. 0) ctr = ctr + 1
        end do

        ! Now we can create the list of equivalent quartet with permutation removed
        allocate (red_quartet(3, ctr))
        ctr = 0
        do j = 1, n
            ctr2 = 0
            do k = j, n
                if (k .eq. j) cycle
                if (all(newqp_sort(:, j) .eq. newqp_sort(:, k))) ctr2 = ctr2 + 1
            end do
            if (ctr2 .eq. 0) then
                ctr = ctr + 1
                red_quartet(1, ctr) = newqp(1, j)
                red_quartet(2, ctr) = newqp(2, j)
                red_quartet(3, ctr) = newqp(3, j)
            end if
        end do

    end block sort_reducible
    call mem%deallocate(newqp, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(newqp_sort, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
end subroutine
