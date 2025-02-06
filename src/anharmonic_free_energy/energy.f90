module energy
!! get the anharmonic free energy
use konstanter, only: r8, i8, lo_twopi, lo_freqtol, lo_imag, lo_sqtol, lo_status, lo_kb_Hartree
use gottochblandat, only: lo_stop_gracefully, tochar, walltime, lo_chop, lo_trueNtimes, &
                          lo_progressbar_init, lo_progressbar, lo_planck, lo_trapezoid_integration, open_file, &
                          lo_flattentensor, lo_sqnorm, lo_linspace, lo_mean, lo_clean_fractional_coordinates
use mpi_wrappers, only: lo_mpi_helper, MPI_DOUBLE_PRECISION, MPI_SUM, MPI_DOUBLE_COMPLEX, MPI_IN_PLACE
use lo_memtracker, only: lo_mem_helper
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_secondorder, only: lo_forceconstant_secondorder
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_forceconstant_fourthorder, only: lo_forceconstant_fourthorder
use type_qpointmesh, only: lo_fft_mesh
use type_phonon_dispersions, only: lo_phonon_dispersions

implicit none
private

public :: perturbative_anharmonic_free_energy

contains

!> Calculates the anharmonic contributions to the free energy
subroutine perturbative_anharmonic_free_energy(p, fct, fcf, qp, dr, temperature, free_energy_thirdorder, free_energy_fourthorder, &
                                               fourthorder, quantum, mw, mem, verbosity)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: p
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> fourth order force constant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-point mesh
    type(lo_fft_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8) :: temperature
    !> free energies
    real(r8), intent(out) :: free_energy_thirdorder, free_energy_fourthorder
    !> what to compute
    logical, intent(in) :: quantum, fourthorder
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    complex(r8), dimension(:, :, :), allocatable :: nuvec1, nuvec2
    real(r8), dimension(:, :), allocatable :: en3, en4
    real(r8) :: t0, t1, timer

    ! Start timers
    timer = walltime()
    t0 = timer
    t1 = timer

    ! Set basic things
    init: block
        real(r8) :: f0, f1
        integer :: iq, imode, ctr, iatom, ix, ialpha

        if (verbosity .gt. 0) then
            write (*, *) ''
            write (*, *) 'CALCULATING ANHARMONIC FREE ENERGY'
        end if

        ! Get the scaled eigenvectors
        call mem%allocate(nuvec1, [dr%n_mode, dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(nuvec2, [dr%n_mode, dr%n_mode, qp%n_full_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        nuvec1 = 0.0_r8
        nuvec2 = 0.0_r8

        ctr = 0
        do iq = 1, qp%n_irr_point
        do imode = 1, dr%n_mode
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            if (dr%iq(iq)%omega(imode) .gt. lo_freqtol) then
                f0 = 1.0_r8/sqrt(dr%iq(iq)%omega(imode))
            else
                f0 = 0.0_r8
            end if
            do iatom = 1, p%na
                f1 = p%invsqrtmass(iatom)
                do ix = 1, 3
                    ialpha = (iatom - 1)*3 + ix
                    nuvec1(ialpha, imode, iq) = dr%iq(iq)%egv(ialpha, imode)*f0*f1
                end do
            end do
        end do
        end do
        do iq = 1, qp%n_full_point
        do imode = 1, dr%n_mode
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            if (dr%aq(iq)%omega(imode) .gt. lo_freqtol) then
                f0 = 1.0_r8/sqrt(dr%aq(iq)%omega(imode))
            else
                f0 = 0.0_r8
            end if
            do iatom = 1, p%na
                f1 = p%invsqrtmass(iatom)
                do ix = 1, 3
                    ialpha = (iatom - 1)*3 + ix
                    nuvec2(ialpha, imode, iq) = dr%aq(iq)%egv(ialpha, imode)*f0*f1
                end do
            end do
        end do
        end do
        call mw%allreduce('sum', nuvec1)
        call mw%allreduce('sum', nuvec2)

        ! Space for free energies per mode
        call mem%allocate(en3, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(en4, [dr%n_mode, qp%n_irr_point], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        en3 = 0.0_r8
        en4 = 0.0_r8

        if (verbosity .gt. 0) then
            t1 = walltime()
            write (*, *) '... scaled vectors and made space (', tochar(t1 - t0), 's)'
            t0 = t1
        end if
    end block init

    ! Integrate out the free energy
    third: block
        real(r8), external :: zdotu
        real(r8), parameter :: threephonon_prefactor = -1.0_r8/48.0_r8
        real(r8), parameter :: fourphonon_prefactor = 1.0_r8/32.0_r8
        real(r8), parameter :: onethird = 1.0_r8/3.0_r8
        complex(r8), dimension(:, :), allocatable :: buf_ev1, buf_ev2, buf_ev3
        complex(r8), dimension(:), allocatable :: evp1, evp2, evp3, ptf
        complex(r8) :: c0
        real(r8), dimension(3) :: qv1, qv2, qv3
        real(r8) :: psisq, sigma, om1, om2, om3, n1, n2, n3, f1, f2, prefactor, s1, s2, s3
        integer :: ctr, b1, b2, b3, q1, q2, q3

        ! Space for intermediate products
        call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(ptf, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_ev1, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_ev2, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%allocate(buf_ev3, [dr%n_mode, dr%n_mode], persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        !call mem%allocate(buf_ev4,[dr%n_mode,dr%n_mode],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
        evp1 = 0.0_r8
        evp2 = 0.0_r8
        ptf = 0.0_r8
        buf_ev1 = 0.0_r8
        buf_ev2 = 0.0_r8
        buf_ev3 = 0.0_r8
        !buf_ev4=0.0_r8

        ! Optimizations:
        ! 1) check the prefactors, can maybe do those before in a neater way
        ! 2) only have to go over half the q-vectors I think.
        if (verbosity .gt. 0) call lo_progressbar_init()
        en3 = 0.0_r8
        ctr = 0
        do q1 = 1, qp%n_irr_point
        do q2 = 1, qp%n_full_point
            ! make it parallel
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! locate third q-vector
            q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, qp%griddensity)

            qv1 = qp%ip(q1)%r
            qv2 = qp%ap(q2)%r
            qv3 = qp%ap(q3)%r

            prefactor = threephonon_prefactor*qp%ap(q2)%integration_weight

            ! pre-transform the matrix element
            call pretransform_phi3(fct, qv2, qv3, ptf)
            ! pre-fetch eigenvectors. Maybe a speedup, dunno.
            buf_ev1 = nuvec1(:, :, q1)
            buf_ev2 = nuvec2(:, :, q2)
            buf_ev3 = nuvec2(:, :, q3)

            do b1 = 1, dr%n_mode
                do b2 = 1, dr%n_mode
                    evp1 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), buf_ev2(:, b2), 1, buf_ev1(:, b1), 1, evp1, dr%n_mode)
                    do b3 = 1, dr%n_mode
                        evp2 = 0.0_r8
                        call zgeru(dr%n_mode, dr%n_mode*dr%n_mode, (1.0_r8, 0.0_r8), buf_ev3(:, b3), 1, evp1, 1, evp2, dr%n_mode)
                        evp2 = conjg(evp2)
                        c0 = dot_product(evp2, ptf)
                        !c0=zdotu(dr%n_mode**3,evp2,1,ptf,1)
                        ! And now we have the matrix element.
                        psisq = real(conjg(c0)*c0, r8)
                        ! Get a smearing parameter
                        s1 = dr%default_smearing(b1)
                        s2 = dr%default_smearing(b2)
                        s3 = dr%default_smearing(b3)
                        sigma = sqrt(s1**2 + s2**2 + s3**2)
                        ! Same sigma as the usual case for debugging.
                        !sigma=lo_mean(dr%default_smearing)
                        om1 = dr%iq(q1)%omega(b1)
                        om2 = dr%aq(q2)%omega(b2)
                        om3 = dr%aq(q3)%omega(b3)
                        if (om1 .lt. lo_freqtol) cycle
                        if (om2 .lt. lo_freqtol) cycle
                        if (om3 .lt. lo_freqtol) cycle
                        ! Decide if we take the classical limit for the occupations
                        if (quantum) then
                            n1 = lo_planck(temperature, om1)
                            n2 = lo_planck(temperature, om2)
                            n3 = lo_planck(temperature, om3)

                            ! This is the Wallace expression
                            ! f1=n1*n2+n1+onethird
                            ! f1=3*f1*real(1.0_r8/( om1+om2+om3+lo_imag*sigma ))
                            ! f2=2*n1*n3-n1*n2+n3
                            ! f2=3*f2*real(1.0_r8/( om1+om2-om3+lo_imag*sigma ))
                            ! en3(b1,q1)=en3(b1,q1)+( (f1+f2)*psisq )*prefactor

                            ! Try the Cowley expression instead?
                        else
                            n1 = lo_kb_Hartree*temperature/om1
                            n2 = lo_kb_Hartree*temperature/om2
                            n3 = lo_kb_Hartree*temperature/om3
                        end if
                        f1 = (n1 + 1)*(n2 + n3 + 1) + n2*n3
                        f2 = n1*n2 + n1*n3 - n2*n3 + 1
                        f2 = 3*f2*real(1.0_r8/(om1 + om2 - om3 + lo_imag*sigma))
                        f1 = f1*real(1.0_r8/(om1 + om2 + om3 + lo_imag*sigma))
                        en3(b1, q1) = en3(b1, q1) + ((f1 + f2)*psisq)*prefactor
                    end do
                end do
            end do
        end do
        if (verbosity .gt. 0 .and. q1 .lt. qp%n_irr_point) then
            call lo_progressbar(' ... third order', q1, qp%n_irr_point, walltime() - t0)
        end if
        end do

        ! Some cleanup
        call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

        if (verbosity .gt. 0) then
            t1 = walltime()
            call lo_progressbar(' ... third order', qp%n_irr_point, qp%n_irr_point, t1 - t0)
            t0 = t1
        end if
        call mw%allreduce('sum', en3)

        if (fourthorder) then
            ! Do the fourth order? Might as well.
            call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(evp2, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(evp3, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%allocate(ptf, dr%n_mode**4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            evp1 = 0.0_r8
            evp2 = 0.0_r8
            evp3 = 0.0_r8
            ptf = 0.0_r8

            if (verbosity .gt. 0) call lo_progressbar_init()
            en4 = 0.0_r8
            ctr = 0
            do q1 = 1, qp%n_irr_point
            do q2 = 1, qp%n_full_point
                ! make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle
                qv1 = qp%ip(q1)%r
                qv2 = qp%ap(q2)%r
                prefactor = fourphonon_prefactor*qp%ap(q2)%integration_weight

                ! pre-transform the matrix element
                call pretransform_phi4(fcf, qv1, qv2, ptf)
                ! pre-fetch eigenvectors. Maybe a speedup, dunno.
                buf_ev1 = nuvec1(:, :, q1)
                buf_ev2 = nuvec2(:, :, q2)
                do b1 = 1, dr%n_mode
                do b2 = 1, dr%n_mode
                    om1 = dr%iq(q1)%omega(b1)
                    om2 = dr%aq(q2)%omega(b2)
                    if (om1 .lt. lo_freqtol) cycle
                    if (om2 .lt. lo_freqtol) cycle
                    ! Decide if we take the classical limit for the occupations
                    if (quantum) then
                        n1 = lo_planck(temperature, om1)
                        n2 = lo_planck(temperature, om2)
                    else
                        n1 = lo_kb_Hartree*temperature/om1
                        n2 = lo_kb_Hartree*temperature/om2
                    end if
                    ! Now to get
                    evp1 = 0.0_r8
                    evp2 = 0.0_r8
                    evp3 = 0.0_r8
                    call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), buf_ev1(:, b1), 1, buf_ev1(:, b1), 1, evp1, dr%n_mode)
                    call zgerc(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), buf_ev2(:, b2), 1, buf_ev2(:, b2), 1, evp2, dr%n_mode)
                    call zgeru(dr%n_mode**2, dr%n_mode**2, (1.0_r8, 0.0_r8), evp2, 1, evp1, 1, evp3, dr%n_mode**2)
                    evp3 = conjg(evp3)
                    psisq = real(dot_product(evp3, ptf), r8)
                    f1 = (2*n1 + 1)*(2*n2 + 1)*psisq*prefactor

                    en4(b1, q1) = en4(b1, q1) + f1
                end do
                end do
            end do
            if (verbosity .gt. 0 .and. q1 .lt. qp%n_irr_point) then
                call lo_progressbar(' ... fourth order', q1, qp%n_irr_point, walltime() - t0)
            end if
            end do

            call mem%deallocate(evp1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(evp2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(evp3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(ptf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_ev1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_ev2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
            call mem%deallocate(buf_ev3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

            if (verbosity .gt. 0) then
                t1 = walltime()
                call lo_progressbar(' ... fourth order', qp%n_irr_point, qp%n_irr_point, t1 - t0)
                t0 = t1
            end if
            call mw%allreduce('sum', en4)
        end if
    end block third

    ! add things up
    finalize: block
        integer :: q1

        ! Add it up
        free_energy_thirdorder = 0.0_r8
        free_energy_fourthorder = 0.0_r8
        do q1 = 1, qp%n_irr_point
            free_energy_thirdorder = free_energy_thirdorder + sum(en3(:, q1))*qp%ip(q1)%integration_weight
            if (fourthorder) then
                free_energy_fourthorder = free_energy_fourthorder + sum(en4(:, q1))*qp%ip(q1)%integration_weight
            end if
        end do
        ! And normalize it to be per atom
        free_energy_thirdorder = free_energy_thirdorder/p%na
        free_energy_fourthorder = free_energy_fourthorder/p%na

        ! and cleanup
        call mem%deallocate(nuvec1, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(nuvec2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(en3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
        call mem%deallocate(en4, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    end block finalize

    if (verbosity .gt. 0) then
        t1 = walltime()
        write (*, *) '... got anharmonic free energy (', tochar(t1 - timer), 's)'
    end if

end subroutine

!> pre-transform to get half of the third order matrix element
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
            ptf(l) = ptf(l) + fct%atom(a1)%triplet(t)%m(i, j, k)*expiqr
        end do
        end do
        end do
    end do
    end do
end subroutine

!> pre-transform to get half of the fourth order matrix element
subroutine pretransform_phi4(fcf, q1, q2, ptf)
    !> third order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> q-vectors
    real(r8), dimension(3), intent(in) :: q1, q2
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

        iqr = -dot_product(q1, rv2) + dot_product(q2, rv3) - dot_product(q2, rv4)
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
            ptf(m) = ptf(m) + fcf%atom(a1)%quartet(q)%m(i, j, k, l)*expiqr
        end do
        end do
        end do
        end do
    end do
    end do
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
