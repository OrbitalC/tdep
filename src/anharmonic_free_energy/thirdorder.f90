module thirdorder
!! Compute the third order contribution to the anharmonic free energy
use konstanter, only: r8, lo_freqtol, lo_twopi, lo_imag, lo_kb_Hartree, lo_exitcode_param
use gottochblandat, only: walltime, lo_trueNtimes, lo_progressbar_init, lo_progressbar, lo_planck, &
                          lo_planck_deriv, lo_planck_secondderiv
use mpi_wrappers, only: lo_mpi_helper, lo_stop_gracefully
use lo_memtracker, only: lo_mem_helper
use type_qpointmesh, only: lo_qpoint_mesh, lo_fft_mesh
use type_crystalstructure, only: lo_crystalstructure
use type_forceconstant_thirdorder, only: lo_forceconstant_thirdorder
use type_phonon_dispersions, only: lo_phonon_dispersions
use lo_timetracker, only: lo_timer

implicit none
private

public :: free_energy_thirdorder

contains
!> Calculate the third order free energy
subroutine free_energy_thirdorder(uc, fct, qp, dr, temperature, fe3, s3, cv3, quantum, mw, mem)
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> third order force constant
    type(lo_forceconstant_thirdorder), intent(in) :: fct
    !> q-point mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8) :: temperature
    !> free energies and heat capacity
    real(r8), intent(out) :: fe3, s3, cv3
    !> what to compute
    logical, intent(in) :: quantum
    !> mpi helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem

    !> Frequency scaled eigenvectors
    complex(r8), dimension(:), allocatable :: egv1, egv2, egv3
    !> Helper for Fourier transform of psi3
    complex(r8), dimension(:), allocatable :: ptf, evp1, evp2
    !> Frequencies, bose-einstein occupation and scattering strength and some other buffer
    real(r8) :: sigma, om1, om2, om3, n2, n3, psisq, f0, f1, plf0, plf1, perm, f2, n1, pref, t0, prefactor
    !>
    real(r8) :: sig1, sig2, sig3, mult, dn1, dn2, dn3, ddn1, ddn2, ddn3, df1, df2, ddf1, ddf2, df0, ddf0
    !> The complex threephonon matrix element
    complex(r8) :: c0
    !> Integers for do loops and counting
    integer :: qi, q1, q2, q3, q2p, q3p, b1, b2, b3, i, ctr
    !> The dimension of the q-grid
    integer, dimension(3) :: dims

    ! We start by allocating everything
    call mem%allocate(ptf, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp1, dr%n_mode**2, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(evp2, dr%n_mode**3, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv1, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv2, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv3, dr%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    t0 = walltime()

    select type (qp)
    class is (lo_fft_mesh)
        dims = qp%griddensity
    class default
        call lo_stop_gracefully(['This routine only works with FFT meshes'], lo_exitcode_param, __FILE__, __LINE__)
    end select

    fe3 = 0.0_r8
    cv3 = 0.0_r8

    do q1=1, qp%n_irr_point
    do q2=1, qp%n_full_point
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle
        q3 = fft_third_grid_index(qp%ip(q1)%full_index, q2, dims)
        if (q3 .lt. q2) cycle

        ! The prefactor to take into account what we are skipping
        if (q2 .eq. q3) then
            mult = 1.0_r8
        else
            mult = 2.0_r8
        end if

        prefactor = qp%ip(q1)%integration_weight*qp%ap(q2)%integration_weight*mult/uc%na
        ! pre-transform the matrix element
        call pretransform_phi3(fct, qp%ap(q2)%r, qp%ap(q3)%r, ptf)

        do b1=1, dr%n_mode
            om1 = dr%iq(q1)%omega(b1)
            if (om1 .lt. lo_freqtol) cycle
            egv1 = dr%iq(q1)%egv(:, b1)/sqrt(om1)
            do b2=1, dr%n_mode
                om2 = dr%aq(q2)%omega(b2)
                if (om2 .lt. lo_freqtol) cycle

                egv2 = dr%aq(q2)%egv(:, b2)/sqrt(om2)
                evp1 = 0.0_r8
                call zgeru(dr%n_mode, dr%n_mode, (1.0_r8, 0.0_r8), egv2, 1, egv1, 1, evp1, dr%n_mode)
                do b3=1, dr%n_mode
                    om3 = dr%aq(q3)%omega(b3)
                    if (om3 .lt. lo_freqtol) cycle
                    egv3 = dr%aq(q3)%egv(:, b3)/sqrt(om3)

                    evp2 = 0.0_r8
                    call zgeru(dr%n_mode, dr%n_mode**2, (1.0_r8, 0.0_r8), egv3, 1, evp1, 1, evp2, dr%n_mode)
                    evp2 = conjg(evp2)
                    c0 = dot_product(evp2, ptf)
                    psisq = real(conjg(c0)*c0, r8) * prefactor

                    if (quantum) then
                        n1 = lo_planck(temperature, om1)
                        n2 = lo_planck(temperature, om2)
                        n3 = lo_planck(temperature, om3)

                        dn1 = lo_planck_deriv(temperature, om1)
                        dn2 = lo_planck_deriv(temperature, om1)
                        dn3 = lo_planck_deriv(temperature, om1)

                        ddn1 = lo_planck_secondderiv(temperature, om1)
                        ddn2 = lo_planck_secondderiv(temperature, om1)
                        ddn3 = lo_planck_secondderiv(temperature, om1)

                        f1 = (n1 + 1.0_r8)*(n2 + n3 + 1.0_r8) + n2*n3
                        f2 = n1*n2 + n1*n3 - n2*n3 + 1.0_r8

                        df1 = dn1 * (n2 + n3 + 1.0_r8) + (n1 + 1.0_r8) * (dn2 + dn3) + dn2 * n3 + n2 * dn3
                        df2 = dn1 * n2 + n1 * dn2 + dn1 * n3 + n1 * dn3 - dn2 * n3 - n2 * dn3

                        ddf1 = ddn1 * (n2 + n3 + 1.0_r8) + 2.0_r8 * dn1 * (dn2 + dn3) + (n1 + 1.0_r8) * (ddn2 + ddn3) + &
                               ddn2 * n3 + n2 * ddn3 + 2.0_r8 * dn2 * dn3
                        ddf2 = ddn1 * n2 + n1 * ddn2 + 2.0_r8 * dn1 * dn2 + ddn1 * n3 + n1 * ddn3 + 2.0_r8 * dn1 * dn3 - &
                               ddn2 * n3 - n2 * ddn3 - 2.0_r8 * dn2 * dn3

                        f1 = f1/(om1 + om2 + om3)
                        df1 = df1/(om1 + om2 + om3)
                        ddf1 = ddf1/(om1 + om2 + om3)
                        if (abs(om1 + om2 - om3) .lt. lo_freqtol) then
                            ! Get a smearing parameter
                            sig1 = dr%default_smearing(b1)
                            sig2 = dr%default_smearing(b2)
                            sig3 = dr%default_smearing(b3)
                            sigma = sqrt(sig1**2 + sig2**2 + sig3**2)
                            f2 = 3*f2*real(1.0/(om1 + om2 - om3 + lo_imag*sigma), r8)
                            df2 = 3.0_r8 * df2 * real(1.0/(om1 + om2 - om3 + lo_imag*sigma), r8)
                            ddf2 = 3.0_r8 * ddf2 * real(1.0/(om1 + om2 - om3 + lo_imag*sigma), r8)
                        else
                            f2 = 3.0_r8 * f2 / (om1 + om2 - om3)
                            df2 = 3.0_r8 * df2 / (om1 + om2 - om3)
                            ddf2 = 3.0_r8 * ddf2 / (om1 + om2 - om3)
                        end if
                        f0 = (f1 + f2) / 48.0_r8
                        df0 = (df1 + df2) / 48.0_r8
                        ddf0 = (ddf1 + ddf2) * temperature / 48.0_r8
                    else
                        f0 = (lo_kb_Hartree*temperature)**2 / (om1*om2*om3) / 12.0_r8
                        df0 = lo_kb_Hartree**2 * temperature / (om1*om2*om3) / 6.0_r8
                        ddf0 = lo_kb_Hartree**2 * temperature / (om1*om2*om3) / 6.0_r8
                    end if
                    fe3 = fe3 - f0 * psisq
                    s3 = s3 + df0 * psisq
                    cv3 = cv3 + ddf0 * psisq
                end do
            end do
        end do
    end do
    if (mw%talk .and. lo_trueNtimes(q1, 127, qp%n_irr_point)) then
        call lo_progressbar(' ... third order', q1, qp%n_irr_point, walltime() - t0)
    end if
    end do
    if (mw%talk) then
        call lo_progressbar(' ... third order', qp%n_irr_point, qp%n_irr_point, walltime() - t0)
    end if

    call mw%allreduce('sum', fe3)
    call mw%allreduce('sum', cv3)
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
