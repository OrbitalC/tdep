

subroutine compute_threephonon_scattering(il, sr, qp, dr, uc, fct, mcg, rng, thres, &
                                          g0, integrationtype, smearing, mw, mem)
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
    !> The threshold for gaussian integration
    real(r8), intent(in) :: thres
    !> The linewidth for this mode
    real(r8), intent(inout) :: g0
    !> what kind of integration are we doing
    integer, intent(in) :: integrationtype
    !> The smearing width
    real(r8), intent(in) :: smearing
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
    real(r8) :: sigma
    !> Frequencies, bose-einstein occupation and scattering strength
    real(r8) :: om1, om2, om3, plf, psisq, prefactor, fall, f0, f1, f2, f3
    !> The bose-einstein distribution for the modes
    real(r8) :: n2, n3
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops
    integer :: qi, q1, q2, q3, b1, b2, b3, i2, i3, i
    !> Is the triplet irreducible ?
    logical :: isred
    !> If so, what is its multiplicity
    integer :: mult
    !> Also what is the prefactor due to permutation ?
    real(r8) :: permutation_mult
    !> The reducible triplet corresponding to the currently computed triplet
    integer, dimension(:, :), allocatable :: red_triplet

    integer :: aaa

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

    call mem%allocate(qgridfull, mcg%npoints, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mcg%generate_grid(qgridfull, rng)

    prefactor = threephonon_prefactor * mcg%weight

    compute_loop: do qi=1, mcg%npoints
        q2 = qgridfull(qi)
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, mcg%full_dims)
        if (q3 .lt. q2) cycle
        call triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult, mem, red_triplet)
        if (isred) cycle

        ! Let's compute the multiplicity due to permutation now
        if (q2 .eq. q3) then
            permutation_mult = 1.0_r8
        else
            permutation_mult = 2.0_r8
        end if

        ! This get the ifc3 in Fourier space, but not on phonons
        call pretransform_phi3(fct, qp%ap(q2)%r, qp%ap(q3)%r, ptf)
        do b2=1, dr%n_mode
            om2 = dr%aq(q2)%omega(b2)
            if (om2 .lt. lo_freqtol) cycle

            egv2 = dr%aq(q2)%egv(:, b2) / sqrt(om2)

            ! This is the multiplication of eigv of phonons 1 and 2
            evp1 = 0.0_r8
            call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
            do b3=1, dr%n_mode
                om3 = dr%aq(q3)%omega(b3)
                if (om3 .lt. lo_freqtol) cycle
                egv3 = dr%aq(q3)%egv(:, b3) / sqrt(om3)

                select case (integrationtype)
                    case (1)
                        sigma = (1.0_r8*lo_frequency_THz_to_Hartree)*smearing
                    case (2)
                        sigma = sqrt(sr%sigsq(q1, b1) + &
                                     sr%sigsq(qp%ap(q2)%irreducible_index, b2) + &
                                     sr%sigsq(qp%ap(q3)%irreducible_index, b3))
                end select

                ! This is the multiplication of eigv of phonons 1 and 2 and now 3
                evp2 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                evp2 = conjg(evp2)
                ! And with this, we have the IFC3 coming in
                c0 = dot_product(evp2, ptf)
                psisq = abs(c0*conjg(c0)) * prefactor

                ! Let's get the Bose-Einstein distributions
                n2 = sr%be(qp%ap(q2)%irreducible_index, b2)
                n3 = sr%be(qp%ap(q3)%irreducible_index, b3)

                ! The prefactor for the scattering
                f0 = (n2 - n3) * lo_gauss(om1, -om2 + om3, sigma)
                f1 = (n2 - n3) * lo_gauss(om1,  om2 - om3, sigma)
                f2 = (n2 + n3 + 1.0_r8) * lo_gauss(om1, om2 + om3, sigma)
                f3 = (n2 + n3 + 1.0_r8) * lo_gauss(om1, -om2 - om3, sigma)

                ! And this is all scattering reunited
                fall = permutation_mult * psisq * (f0 - f1 + f2 - f3)

                ! And we add everything for each triplet equivalent to the one we are actually computing
                do i=1, mult
                    g0 = g0 + fall

                    i2 = (red_triplet(1, i) - 1) * dr%n_mode + b2
                    sr%Xi(il, i2) = sr%Xi(il, i2) + 2.0_r8 * fall * om2 / om1

                    i3 = (red_triplet(2, i) - 1) * dr%n_mode + b3
                    sr%Xi(il, i3) = sr%Xi(il, i3) + 2.0_r8 * fall * om3 / om1
                end do
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

subroutine triplet_is_irreducible(qp, uc, q1, q2, q3, isred, mult, mem, red_triplet)
    !> The qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> structure
    type(lo_crystalstructure), intent(in) :: uc
    !> The indices of the qpoints - q1 is irreducible, q2 and q3 are full
    integer, intent(in) :: q1, q2, q3
    !> Is the triplet reducible ?
    logical, intent(out) :: isred
    !> If it's reducible, what is its multiplicity
    integer, intent(out) :: mult
    !> Memory helper
    type(lo_mem_helper), intent(inout) :: mem
    !> The equivalent triplet
    integer, dimension(:, :), allocatable, intent(out) :: red_triplet

    !> The new-qpoints and the temporary invariant triplet
    integer, dimension(:, :), allocatable :: newqp, tmp_triplet
    ! To hold the q-point in reduced coordinates
    real(r8), dimension(3) :: qv2, qv3, qv2p, qv3p
    !> To get the index of the new triplet on the fft_grid
    integer, dimension(3) :: gi
    !> The new triplet after the operation
    integer, dimension(2) :: qpp, ind_qsort
    !> Integers for the loops
    integer :: j, k, n, ctr, ctr2

    ! First get the q-points in reduce coordinates
    qv2 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q2)%r)
    qv3 = matmul(uc%inv_reciprocal_latticevectors, qp%ap(q3)%r)

    allocate(newqp(2, qp%ip(q1)%n_invariant_operation))
    newqp = -lo_hugeint

    mult = 0
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
        ! If we got here, it means that the new triplet is on the grid

        ! We need to get the reducible q-triplet before the sorting
        newqp(:, j) = qpp
        call lo_qsort(qpp)
        if (qpp(1) .gt. q2) then
            isred = .true.
        ! Now we have to determine the weight
        ! Two roads are possible
        !      1. Get the ratio of number of red point that can give this irreducible point
        !      2. Look at the ratio between total number of operations and the ones
        !         that leaves this irreducible triplet unchanged
        ! The second road doesn't requires me to sum over all other qpoints, so I'll go with this one
        else if (qpp(1) .eq. q2 .and. qpp(2) .eq. q3) then
            mult = mult + 1
        end if
    end do
    mult = qp%ip(q1)%n_invariant_operation / mult

    ! Now we need the unique triplet after applying symmetries
    call lo_return_unique(newqp, tmp_triplet)
    n = size(tmp_triplet, 2)

    ! We could have cases where triplet (q1, R*q2, R*q3) and (q1, R*q3, R*q2) are equivalent
    ! Since we are using permutation symmetry, we have to remove one of them to avoid double counting
    ! First we count if we have triplet equivalent by permutation, by comparing with all other triplet
    ctr = 0
    do j=1, n
        ctr2 = 0
        do k=j, n
            if (k .eq. j) cycle
            if (tmp_triplet(1, j) .eq. tmp_triplet(2, k)) then
                ctr2 = ctr2 + 1
            end if
        end do
        ! If ctr2 is 0, it means that this guy is unique by permutation
        if (ctr2 .eq. 0) ctr = ctr + 1
    end do

    ! Now we can create the list of equivalent triplet with permutation removed
    allocate(red_triplet(2, ctr))
    ctr = 1
    do j=1, n
        ctr2 = 0
        do k=j, n
            if (k .eq. j) cycle
            if (tmp_triplet(1, j) .eq. tmp_triplet(2, k)) then
                ctr2 = ctr2 + 1
            end if
        end do
        ! If ctr2 is 0, it means that this guy is unique by permutation
        if (ctr2 .eq. 0) then
            red_triplet(1, ctr) = tmp_triplet(1, j)
            red_triplet(2, ctr) = tmp_triplet(2, j)
            ctr = ctr + 1
        end if
    end do
end subroutine
