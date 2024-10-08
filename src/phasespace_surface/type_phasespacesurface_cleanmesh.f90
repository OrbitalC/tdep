
!> slice up the tetrahedrons
subroutine bz_isosurface(qp, fvals, isoval, finr, fint, fingrad, uc, fc, qpoint1, plus, &
                         b1, b2, b3, verbosity, calcintens, fct, psisq, psisq_norm, tot_omega, tot_egv, tot_q, mem)
    !> the kpoint mesh
    type(lo_wedge_mesh), intent(in) :: qp
    !> the function values
    real(flyt), dimension(:), intent(in) :: fvals
    !> the isolevel
    real(flyt), intent(in) :: isoval
    !> the triangles
    integer, dimension(:, :), allocatable, intent(out) :: fint
    !> the points
    real(flyt), dimension(:, :), allocatable, intent(out) :: finr
    !> the gradients
    real(flyt), dimension(:, :), allocatable, intent(out) :: fingrad
    !> to calculate the gradients
    type(lo_crystalstructure), intent(in) :: uc
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    type(lo_qpoint), intent(in) :: qpoint1
    logical, intent(in) :: plus
    integer, intent(in) :: b1, b2, b3
    integer, intent(in), optional :: verbosity
    !> calculate intensities (matrix elements)?
    logical, intent(in), optional :: calcintens
    !> third order force constants (for matrix elements)
    type(lo_forceconstant_thirdorder), intent(in), optional :: fct
    !> the 3-phonon matrix elements
    real(flyt), dimension(:), allocatable, intent(out), optional :: psisq
    !> norm of the 3-phonon matrix elements
    real(flyt), intent(out), optional :: psisq_norm
    !> the frequencies involved (index: point, modelabel)
    real(flyt), dimension(:, :), allocatable, intent(out), optional :: tot_omega
    !> the eigenvectors and q-points involved (index: point, :, modelabel)
    real(flyt), dimension(:, :, :), allocatable, intent(out), optional :: tot_egv, tot_q
    !> memory helper to generate the q-points
    type(lo_mem_helper), intent(inout) :: mem
    !
    integer, dimension(:, :), allocatable :: dumt
    integer, dimension(:), allocatable :: ind !,dumind
    integer, dimension(4) :: teti, sorti
    integer :: i, j, k, l, ii, ll, ntri, nfintri
    real(flyt), dimension(:, :), allocatable :: dumr
    real(flyt), dimension(:, :), allocatable :: dumgrad
    real(flyt), dimension(3, 4) :: cpts
    real(flyt), dimension(4) :: tete
    real(flyt), dimension(3) :: v0, v1, v2, v3, v4
    real(flyt), dimension(3) :: r1, r2, r3, b
    real(flyt) :: e21, e31, e41, e32, e42, e43
    real(flyt) :: f0, f1, f2, f3, t0, tolerance
    logical :: verb
    !for matrix elements (intensities)
    complex(flyt), dimension(fc%na*3, 3) :: egv
    complex(flyt) :: c0
    real(flyt), dimension(3) :: q1, q2, q3, omega
    logical :: calcintensities

    if (present(verbosity)) then
        if (verbosity .gt. 0) then
            verb = .true.
        else
            verb = .false.
        end if
    else
        verb = .false.
    end if
    if (present(calcintens)) then
        calcintensities = calcintens
    else
        calcintensities = .false.
    end if
    if (calcintensities) then
        !assert that third order force constants are included
        if (.not. present(fct)) then
            write (*, *) "you must pass fct if calcintensities=.true."
            stop
        end if
    end if
    ! Fetch the tetrahedron energies
    ! space for possible points
    lo_allocate(dumr(3, qp%n_full_tet*4))
    lo_allocate(dumt(3, qp%n_full_tet*4))
    lo_allocate(dumgrad(3, qp%n_full_tet*4))
    dumr = 0.0_flyt
    dumt = 0
    dumgrad = 0.0_flyt

    tolerance = lo_freqtol
    l = 0
    ll = 0
    if (verb) call lo_progressbar_init()
    do i = 1, qp%n_full_tet
        ! fetch energy at the corners, and coordinates
        do j = 1, 4
            ii = qp%at(i)%full_index(j)
            cpts(:, j) = qp%ap(ii)%r
            tete(j) = fvals(ii) - isoval
        end do
        sorti = sort_four_numbers(tete)
        tete = tete(sorti)
        cpts = cpts(:, sorti)
        ! cycle if outside bounds right away
        if (minval(tete) .gt. 0.0_flyt) cycle
        if (maxval(tete) .lt. 0.0_flyt) cycle
        if (sum(abs(tete)) .lt. lo_freqtol) then
            cycle
        end if

        ! Get the gradient, such that e(k)=e(1)+b.(k-k1) inside the tetrahedron
        ! Just the Lehman & Taut paper.
        f0 = -lo_signed_tetrahedron_volume(cpts)*6
        v0 = cpts(:, 1)
        v1 = cpts(:, 2) - v0
        v2 = cpts(:, 3) - v0
        v3 = cpts(:, 4) - v0
        r1 = lo_cross(v2, v3)/f0
        r2 = lo_cross(v3, v1)/f0
        r3 = lo_cross(v1, v2)/f0
        b = 0.0_flyt
        b = b + (tete(2) - tete(1))*r1
        b = b + (tete(3) - tete(1))*r2
        b = b + (tete(4) - tete(1))*r3
        b = b/norm2(b)

        ! Maybe catch degeneracies
        e21 = tete(2) - tete(1)
        e31 = tete(3) - tete(1)
        e41 = tete(4) - tete(1)
        e32 = tete(3) - tete(2)
        e42 = tete(4) - tete(2)
        e43 = tete(4) - tete(2)

        ! Figure out if I'm going to do something.
        if (tete(1) .lt. 0.0_flyt .and. tete(2) .gt. 0.0_flyt) then
            ! easy case!
            v1 = cpts(:, 1) - tete(1)/(tete(3) - tete(1))*(cpts(:, 3) - cpts(:, 1))
            v2 = cpts(:, 1) - tete(1)/(tete(4) - tete(1))*(cpts(:, 4) - cpts(:, 1))
            v3 = cpts(:, 1) - tete(1)/(tete(2) - tete(1))*(cpts(:, 2) - cpts(:, 1))
            ! store
            ll = ll + 1
            l = l + 1; dumr(:, l) = v1; dumt(1, ll) = l; dumgrad(:, ll) = b
            l = l + 1; dumr(:, l) = v2; dumt(2, ll) = l; dumgrad(:, ll) = b
            l = l + 1; dumr(:, l) = v3; dumt(3, ll) = l; dumgrad(:, ll) = b
        elseif (tete(2) .lt. 0.0_flyt .and. tete(3) .gt. 0.0_flyt) then
            ! annoying case
            if (abs(tete(3) - tete(1)) .gt. tolerance .and. &
                abs(tete(4) - tete(1)) .gt. tolerance .and. &
                abs(tete(4) - tete(2)) .gt. tolerance .and. &
                abs(tete(3) - tete(2)) .gt. tolerance) then

                v1 = cpts(:, 1) - tete(1)/(e31)*(cpts(:, 3) - cpts(:, 1))
                v2 = cpts(:, 1) - tete(1)/(e41)*(cpts(:, 4) - cpts(:, 1))
                v3 = cpts(:, 2) - tete(2)/(e42)*(cpts(:, 4) - cpts(:, 2))
                v4 = cpts(:, 2) - tete(2)/(e32)*(cpts(:, 3) - cpts(:, 2))

                ! Test building two triangles out of this, choose the pair
                ! with as similar area as possible
                f0 = norm2(lo_cross(v2 - v1, v4 - v1)) ! 124
                f1 = norm2(lo_cross(v4 - v2, v3 - v2)) ! 324
                !
                f2 = norm2(lo_cross(v3 - v1, v2 - v1)) ! 132
                f3 = norm2(lo_cross(v2 - v3, v4 - v3)) ! 134
                !
                if (abs(f0 - f1) .lt. abs(f3 - f2)) then
                    ! first case
                    ll = ll + 1
                    l = l + 1; dumr(:, l) = v1; dumt(1, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v2; dumt(2, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v4; dumt(3, ll) = l; dumgrad(:, ll) = b
!                dumind(ll)=i

                    ll = ll + 1
                    l = l + 1; dumr(:, l) = v2; dumt(1, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v3; dumt(2, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v4; dumt(3, ll) = l; dumgrad(:, ll) = b
                else
                    ! second
                    ll = ll + 1
                    l = l + 1; dumr(:, l) = v1; dumt(1, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v2; dumt(2, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v3; dumt(3, ll) = l; dumgrad(:, ll) = b
                    !
                    ll = ll + 1
                    l = l + 1; dumr(:, l) = v1; dumt(1, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v3; dumt(2, ll) = l; dumgrad(:, ll) = b
                    l = l + 1; dumr(:, l) = v4; dumt(3, ll) = l; dumgrad(:, ll) = b
                end if
            end if
            !
        elseif (tete(3) .lt. 0.0_flyt .and. tete(4) .gt. 0.0_flyt) then
            ! easy again!
            v1 = cpts(:, 4) - tete(4)/(tete(2) - tete(4))*(cpts(:, 2) - cpts(:, 4))
            v2 = cpts(:, 4) - tete(4)/(tete(1) - tete(4))*(cpts(:, 1) - cpts(:, 4))
            v3 = cpts(:, 4) - tete(4)/(tete(3) - tete(4))*(cpts(:, 3) - cpts(:, 4))
            ! store
            ll = ll + 1
            l = l + 1; dumr(:, l) = v1; dumt(1, ll) = l; dumgrad(:, ll) = b
            l = l + 1; dumr(:, l) = v2; dumt(2, ll) = l; dumgrad(:, ll) = b
            l = l + 1; dumr(:, l) = v3; dumt(3, ll) = l; dumgrad(:, ll) = b
        end if
        !
        if (verb) then
            if (mod(i, 20) .eq. 0) then
                call lo_progressbar(' ... building triangles', i, qp%n_full_tet)
            end if
        end if
    end do

    if (verb) call lo_progressbar(' ... building triangles', qp%n_full_tet, qp%n_full_tet)

    ! Stop if I found no points
    if (l .eq. 0) return

    ! Add together gradients from triangles.

    ! Remove the redundant points. Do this twice, with two different sets of verlet boxes.
    call remove_redundant_points(dumr(:, 1:l), dumt(:, 1:ll), finr, fint, 6)
    if (size(finr, 2) .lt. 3) then
        ! the reduction might remove all points
        if (allocated(finr)) deallocate (finr)
        if (allocated(fint)) deallocate (fint)
        if (allocated(fingrad)) deallocate (fingrad)
        return
    end if
    deallocate (dumr, dumt)
    call remove_redundant_points(finr, fint, dumr, dumt, 5)
    deallocate (finr, fint)
    if (size(dumr, 2) .lt. 3) then
        ! the reduction might remove all points
        return
    end if

    ! Remove the too tiny triangles.
    ntri = ll
    allocate (ind(ntri))
    ind = 1
    ! find triangles with zero area
    do i = 1, ntri
        ! get the area of this triangle
        v1 = dumr(:, dumt(1, i)) - dumr(:, dumt(3, i))
        v2 = dumr(:, dumt(2, i)) - dumr(:, dumt(3, i))
        if (norm2(lo_cross(v1, v2)) .lt. lo_sqtol) then
            ind(i) = 0
        end if
    end do

    nfintri = sum(ind)
    ! make space for the nice things
    allocate (finr(3, size(dumr, 2)))
    allocate (fingrad(3, nfintri))
    allocate (fint(3, nfintri))
    fingrad = 0.0_flyt
    fint = 0
    finr = 0.0_flyt
    l = 0
    do i = 1, ntri
        if (ind(i) .eq. 1) then
            l = l + 1
            fint(:, l) = dumt(:, i)
            fingrad(:, l) = dumgrad(:, i)
        end if
    end do
    finr = dumr

    ! Fix the gradients
    dumr = 0.0_flyt
    do i = 1, nfintri
        do j = 1, 3
            k = fint(j, i)
            dumr(:, k) = dumr(:, k) + fingrad(:, i)
        end do
    end do
    lo_deallocate(fingrad)
    lo_allocate(fingrad(3, size(dumr, 2)))
    fingrad = dumr

    do i = 1, size(fingrad, 2)
        fingrad(:, i) = fingrad(:, i)/norm2(fingrad(:, i))
    end do

    ! Calculate the intensities? If not, we're done.
    if ( .not. calcintensities) return

    print *, "... ntri = ", ntri, size(fint, 2), ", calculate intensities"
    fixgrad: block
        integer :: n, nb
        type(lo_qpoint) :: qp1, qp2, qp3
        type(lo_phonon_dispersions_qpoint) :: drp1, drp2, drp3
        complex(flyt), dimension(:, :), allocatable :: D
        complex(flyt), dimension(:, :, :), allocatable :: Dq
        n = size(finr, 2)
        ! allocate (fingrad(3, size(dumr, 2)))
        fingrad = 0.0_flyt
        !
        nb = fc%na*3
        lo_allocate(drp1%omega(nb))
        lo_allocate(drp1%egv(nb, nb))
        lo_allocate(drp1%vel(3, nb))
        lo_allocate(drp2%omega(nb))
        lo_allocate(drp2%egv(nb, nb))
        lo_allocate(drp2%vel(3, nb))
        lo_allocate(drp3%omega(nb))
        lo_allocate(drp3%egv(nb, nb))
        lo_allocate(drp3%vel(3, nb))
        lo_allocate(D(nb, nb))
        lo_allocate(Dq(3, nb, nb))
        lo_allocate(psisq(n))
        lo_allocate(tot_omega(n, 3))
        lo_allocate(tot_egv(n, nb, 3))
        lo_allocate(tot_q(n, 3, 3))
        !
        qp1%r = qpoint1%r
        call lo_get_small_group_of_qpoint(qp1, uc)
        call drp1%generate(fc, uc, mem, qp1)
        psisq = 0.0_flyt
        tot_omega = 0.0_flyt
        tot_egv = 0.0_flyt
        tot_q = 0.0_flyt
        !
        if (verb) call lo_progressbar_init()
        do i = 1, n
            qp2%r = finr(:, i)
            qp3%r = qp1%r + qp2%r
            qp3%r = qp3%r - uc%bz%gshift(qp3%r)
            call lo_get_small_group_of_qpoint(qp2, uc)
            call lo_get_small_group_of_qpoint(qp3, uc)
            call drp2%generate(fc, uc, mem, qp2)
            call drp3%generate(fc, uc, mem, qp3)
            if (plus) then
                fingrad(:, i) = drp1%vel(:, b1) + drp2%vel(:, b2) - drp3%vel(:, b3)
            else
                fingrad(:, i) = drp1%vel(:, b1) - drp2%vel(:, b2) - drp3%vel(:, b3)
            end if
            ! Calculate intensities
            if (calcintens) then
                if (verb) then
                    if (mod(i, 20) .eq. 0) then
                        call lo_progressbar(' ... calculating intensities', i, n)
                    end if
                end if
                ! below is like thermal_conductivity/scatteringstrengths.f90
                ! (all this should be the same for plus and minus, right?)
                ! q-vectors with correct sign
                q1 = qp1%r*lo_twopi !is this correct? (really not needed)
                q2 = -qp2%r*lo_twopi
                q3 = -qp3%r*lo_twopi
                ! frequencies, eigenvectors, q-vectors
                omega(1) = drp1%omega(b1)
                omega(2) = drp2%omega(b2)
                omega(3) = drp3%omega(b3)
                egv(:, 1) = drp1%egv(:, b1)
                egv(:, 2) = drp2%egv(:, b2)
                egv(:, 3) = drp3%egv(:, b3)

                ! check for small values in omega
                if (any(omega .lt. lo_freqtol)) then
                    write (*, *) "*** small omega in scattering amplitude, set to 0"
                    c0 = 0.0_flyt
                else
                    ! otherwise compute the scattering amplitude
                    c0 = fct%scatteringamplitude(omega, egv, q2, q3)
                end if

                ! should not happend anymore prayer emoji
                if (c0 /= c0) then  ! <- checks for NaN
                    print *, "omega = ", omega
                    print *, "q2 = ", q2
                    print *, "q3 = ", q3
                    print *, "c0 = ", c0
                    write (*, *) "*** NaN in scattering amplitude, set to lo_huge = ", lo_huge
                    c0 = sqrt(lo_huge)
                end if
                psisq(i) = abs(conjg(c0)*c0)
                ! put everything else together here so it's handy
                tot_omega(i, :) = omega
                tot_egv(i, :, :) = egv
                tot_q(i, :, 1) = q1
                tot_q(i, :, 2) = q2
                tot_q(i, :, 3) = q3
            end if
        end do
        print *, "... found ", n, " points."
        write (*, '(A, E15.5)') " ... max(abs(psisq)) = ", maxval(abs(psisq))
        psisq_norm = norm2(psisq)
        write (*, '(A, E15.5)') " ... psisq_norm      = ", psisq_norm
        if (calcintens) then
            call lo_progressbar(' ... calculating intensities', n, n)
        end if
    end block fixgrad
    !
end subroutine

!> Pretty fast way of sorting four numbers. Stole it from stack overflow.
pure function sort_four_numbers(n) result(ind)
    !> the four numbers
    real(flyt), dimension(4), intent(in) :: n
    !> the output order so that n(ind) is sorted
    integer, dimension(4) :: ind
    !
    integer :: low1, high1, low2, high2, highest, lowest, middle1, middle2
    !
    if (n(1) <= n(2)) then
        low1 = 1
        high1 = 2
    else
        low1 = 2
        high1 = 1
    end if

    if (n(3) <= n(4)) then
        low2 = 3
        high2 = 4
    else
        low2 = 4
        high2 = 3
    end if

    if (n(low1) <= n(low2)) then
        lowest = low1
        middle1 = low2
    else
        lowest = low2
        middle1 = low1
    end if

    if (n(high1) >= n(high2)) then
        highest = high1
        middle2 = high2
    else
        highest = high2
        middle2 = high1
    end if

    if (n(middle1) < n(middle2)) then
        ind = (/lowest, middle1, middle2, highest/)
    else
        ind = (/lowest, middle2, middle1, highest/)
    end if
    !
end function

!> use Verlet boxes to remove redundant points
subroutine remove_redundant_points(r, tri, ur, utri, nb)
    !> original points
    real(flyt), dimension(:, :), intent(in) :: r
    !> original triangles
    integer, dimension(:, :), intent(in) :: tri
    !> the unique points
    real(flyt), dimension(:, :), allocatable, intent(out) :: ur
    !> the new, fancy triangles
    integer, dimension(:, :), allocatable, intent(out) :: utri
    !> number of boxes to use
    integer, intent(in) :: nb
    !
    type(lo_points_in_boxes), dimension(:, :, :), allocatable :: b
    integer(flyt), dimension(:), allocatable :: dumi, dumj, dumk
    integer, dimension(3) :: gi
    integer :: i, j, k, l
    integer :: np, ntri
    !
    real(flyt), dimension(:, :), allocatable :: dum
    real(flyt), dimension(3) :: minv, maxv, v1, v2
    real(flyt) :: t0, t1

    ! Get the boundaries of the box
    t0 = walltime()
    np = size(r, 2)
    ntri = size(tri, 2)

    ! Get the box dimensions
    minv = lo_huge
    maxv = -lo_huge
    do i = 1, np
        do j = 1, 3
            minv(j) = min(r(j, i), minv(j))
            maxv(j) = max(r(j, i), maxv(j))
        end do
    end do

    ! It could be just a single point or something like that
    if (lo_sqnorm(minv - maxv) .lt. lo_sqtol) then
        lo_allocate(ur(3, 1))
        lo_allocate(utri(3, 1))
        ur = 0.0_flyt
        utri = 0
        return
    end if
    ! pad it a little
    minv = minv - abs(maxv - minv)*0.001_flyt
    maxv = maxv + abs(maxv - minv)*0.001_flyt

    ! Sort these into boxes, first reset the counter
    allocate (b(nb, nb, nb))
    do i = 1, nb
    do j = 1, nb
    do k = 1, nb
        b(i, j, k)%n = 0
    end do
    end do
    end do
    ! Count particles per box
    do i = 1, np
        gi = box_from_coordinates(r(:, i), minv, maxv, nb)
        b(gi(1), gi(2), gi(3))%n = b(gi(1), gi(2), gi(3))%n + 1
    end do
    ! Make some space
    l = 0
    do i = 1, nb
    do j = 1, nb
    do k = 1, nb
        if (b(i, j, k)%n .gt. 0) then
            l = l + b(i, j, k)%n
            lo_allocate(b(i, j, k)%ind(b(i, j, k)%n))
        end if
    end do
    end do
    end do
    ! Stuff the boxes
    do i = 1, nb
    do j = 1, nb
    do k = 1, nb
        b(i, j, k)%n = 0
    end do
    end do
    end do
    do i = 1, np
        gi = box_from_coordinates(r(:, i), minv, maxv, nb)
        b(gi(1), gi(2), gi(3))%n = b(gi(1), gi(2), gi(3))%n + 1
        b(gi(1), gi(2), gi(3))%ind(b(gi(1), gi(2), gi(3))%n) = i
    end do

    ! Find the unique
    lo_allocate(dumi(np))
    lo_allocate(dumj(np))
    lo_allocate(dumk(np))
    lo_allocate(dum(3, np))
    !
    dumi = 0
    dumj = 1
    l = 0
    do i = 1, np
        ! Which box is it in?
        gi = box_from_coordinates(r(:, i), minv, maxv, nb)
        ! Compare with points in the same box
        l = 0
        do j = 1, b(gi(1), gi(2), gi(3))%n
            k = b(gi(1), gi(2), gi(3))%ind(j)
            if (k .gt. i) then
                if (lo_sqnorm(r(:, i) - r(:, k)) .lt. lo_sqtol) then
                    l = l + 1
                    dumj(k) = 0
                    dumi(k) = i
                end if
            end if
        end do
        !
    end do

    ! Figure out the mapping
    l = sum(dumj)
    lo_allocate(ur(3, l))
    j = 0
    dumk = 0
    do i = 1, np
        if (dumj(i) .eq. 1) then
            j = j + 1
            ur(:, j) = r(:, i)
            dumk(i) = j
        end if
    end do

    ! And the mapping for the redundant points
    do i = 1, np
        if (dumk(i) .eq. 0) then
            j = i
            do
                j = dumi(j)
                k = dumk(j)
                if (k .ne. 0) exit
            end do
            dumk(i) = k
        end if
    end do

    ! and the triangles
    lo_allocate(utri(3, ntri))
    do i = 1, ntri
        do j = 1, 3
            utri(j, i) = dumk(tri(j, i))
        end do
    end do

    ! and some cleanup
    lo_deallocate(dum)
    lo_deallocate(dumi)
    lo_deallocate(dumj)
    do i = 1, nb
    do j = 1, nb
    do k = 1, nb
        if (b(i, j, k)%n .gt. 0) then
            lo_deallocate(b(i, j, k)%ind)
        end if
    end do
    end do
    end do
    lo_deallocate(b)
contains

    function box_from_coordinates(r, minv, maxv, nb) result(gi)
        real(flyt), dimension(3), intent(in) :: r, minv, maxv
        integer, intent(in) :: nb
        integer, dimension(3) :: gi
        !
        integer :: i
        !
        do i = 1, 3
            gi(i) = floor(nb*(r(i) - minv(i))/(maxv(i) - minv(i))) + 1
        end do
    end function

end subroutine

