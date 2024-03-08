
!> Gaussian integration of the isotope self energy
subroutine isotope_imaginary_selfenergy_gaussian(wp, qp, dr, se, sr, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> scattering rates
    type(lo_listofscatteringrates), intent(in) :: sr
    !> MPI communicator
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: buf, sigma
    real(r8), dimension(3) :: v0
    real(r8) :: prefactor, invf, t0
    integer :: q, b1, b2, i, ii, jj, ctr

    ! Some constants
    t0 = walltime()
    invf = se%n_energy/se%energy_axis(se%n_energy)
    ! Some space
    call mem%allocate(sigma, se%n_mode, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (verbosity .gt. 0) call lo_progressbar_init()
    ! start integrating
    ctr = 0
    se%im_iso = 0.0_r8
    if (sr%atgamma) then
        do q = 1, qp%n_irr_point
        do b1 = 1, dr%n_mode
            ! MPI thing
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! Get the smearing parameter, if adaptive
            select case (se%integrationtype)
            case (1)
                do b2 = 1, dr%n_mode
                    sigma(b2) = dr%default_smearing(b2)*se%smearing_prefactor
                end do
            case (2)
                do b2 = 1, dr%n_mode
                    v0 = dr%iq(q)%vel(:, b2)
                    sigma(b2) = qp%smearingparameter(v0, dr%default_smearing(b2), se%smearing_prefactor)
                end do
            end select
            ! Get the prefactor
            prefactor = isotope_prefactor*qp%ip(q)%integration_weight

            ! sum things up
            do b2 = 1, dr%n_mode
                ! Add it to the self-energy
                ii = max(floor((dr%iq(q)%omega(b2) - 4*sigma(b2))*invf), 1)
                jj = min(ceiling((dr%iq(q)%omega(b2) + 4*sigma(b2))*invf), se%n_energy)
                do i = ii, jj
                    se%im_iso(i, b1) = se%im_iso(i, b1) + prefactor*lo_gauss(se%energy_axis(i), dr%iq(q)%omega(b2), sigma(b2))*sr%psi_iso(b1, b2, q)
                end do
            end do

            ! Report?
            if (verbosity .gt. 0) then
                if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... isotope imaginary selfenergy', ctr, dr%n_mode*qp%n_irr_point)
            end if
        end do
        end do
    else
        do q = 1, qp%n_full_point
        do b1 = 1, dr%n_mode
            ! MPI thing
            ctr = ctr + 1
            if (mod(ctr, mw%n) .ne. mw%r) cycle
            ! Get the smearing parameter, if adaptive
            select case (se%integrationtype)
            case (1)
                do b2 = 1, dr%n_mode
                    sigma(b2) = dr%default_smearing(b2)*se%smearing_prefactor
                end do
            case (2)
                do b2 = 1, dr%n_mode
                    v0 = dr%aq(q)%vel(:, b2)
                    sigma(b2) = qp%smearingparameter(v0, dr%default_smearing(b2), se%smearing_prefactor)
                end do
            end select
            ! Get the prefactor
            prefactor = isotope_prefactor*qp%ap(q)%integration_weight

            ! sum things up
            do b2 = 1, dr%n_mode
                ! Add it to the self-energy
                ii = max(floor((dr%aq(q)%omega(b2) - 4*sigma(b2))*invf), 1)
                jj = min(ceiling((dr%aq(q)%omega(b2) + 4*sigma(b2))*invf), se%n_energy)
                do i = ii, jj
                    se%im_iso(i, b1) = se%im_iso(i, b1) + prefactor*lo_gauss(se%energy_axis(i), dr%aq(q)%omega(b2), sigma(b2))*sr%psi_iso(b1, b2, q)
                end do
            end do

            ! Report?
            if (verbosity .gt. 0) then
                if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... isotope imaginary selfenergy', ctr, dr%n_mode*qp%n_full_point)
            end if
        end do
        end do
    end if
    ! Sum it up
    call mw%allreduce('sum', se%im_iso)

    ! Fix degeneracies
    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8
    do b1 = 1, dr%n_mode
        buf = 0.0_r8
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            buf = buf + se%im_iso(:, b2)
        end do
        buf = buf/real(wp%degeneracy(b1), r8)
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            se%im_iso(:, b2) = buf
        end do
    end do
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%deallocate(sigma, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (verbosity .gt. 0) call lo_progressbar(' ... isotope imaginary selfenergy', dr%n_mode*qp%n_full_point, dr%n_mode*qp%n_full_point, walltime() - t0)
end subroutine

!> threephonon self energy gaussian
subroutine threephonon_imaginary_selfenergy_gaussian(wp, se, sr, qp, dr, temperature, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> scattering rates
    type(lo_listofscatteringrates), intent(inout) :: sr
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: buf
    real(r8) :: psisquare, psisq23, psisq32
    real(r8) :: sigma, om1, om2, om3, invf, s2, s3
    real(r8) :: plf1, plf2, t0, pref
    integer :: q, ctr, i, ii, jj, b1, b2, b3, ilo, ihi

    ! set the unit
    t0 = walltime()

    invf = se%n_energy/se%energy_axis(se%n_energy)

    ctr = 0
    se%im_3ph = 0.0_r8
    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    buf = 0.0_r8

    if (sr%atgamma) then
        ! If we are doing the Gamma-point, we can use symmetry
        qloopirr: do q = 1, qp%n_irr_point
            do b1 = 1, dr%n_mode
                ! Make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                ! Prefactor for this q-point
                pref = threephonon_prefactor*qp%ip(q)%integration_weight
                ! Slightly smarter version that should be a little faster.
                do b2 = 1, dr%n_mode
                do b3 = b2, dr%n_mode
                    ! reset buffer
                    buf = 0.0_r8
                    ilo = lo_hugeint
                    ihi = -lo_hugeint
                    ! Get the smearing parameter, in case it's adaptive
                    select case (se%integrationtype)
                    case (1)
                        !sigma=(dr%default_smearing(b2)+dr%default_smearing(b3))*0.5_r8*se%smearing_prefactor
                        s2 = dr%default_smearing(b2)
                        s3 = dr%default_smearing(b3)
                        sigma = sqrt(s2**2 + s3**2)
                    case (2)
                        s2 = qp%adaptive_sigma(qp%ip(q)%radius, dr%iq(q)%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                        s3 = qp%adaptive_sigma(qp%ip(q)%radius, dr%iq(q)%vel(:, b3), dr%default_smearing(b3), se%smearing_prefactor)
                        sigma = sqrt(s2**2 + s3**2)
                        !v0=dr%iq(q)%vel(:,b2)
                        !v1=dr%iq(q)%vel(:,b3)
                        !sigma=qp%smearingparameter(v0-v1,(dr%default_smearing(b2)+dr%default_smearing(b3))*0.5_r8,se%smearing_prefactor)
                    end select

                    ! fetch frequencies
                    om1 = wp%omega(b1)
                    om2 = dr%iq(q)%omega(b2)
                    om3 = dr%iq(q)%omega(b3)
                    ! Get the S-function
                    plf1 = 1.0_r8 + lo_planck(temperature, om2) + lo_planck(temperature, om3)
                    plf2 = lo_planck(temperature, om3) - lo_planck(temperature, om2)

                    ! delta(bigOM-om2-om3)
                    ii = max(floor((om2 + om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((om2 + om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) + plf1*lo_gauss(se%energy_axis(i), om2 + om3, sigma)
                    end do
                    ! delta(bigOM+om2+om3)
                    ii = max(floor((-om2 - om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((-om2 - om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) - plf1*lo_gauss(se%energy_axis(i), -om2 - om3, sigma)
                    end do
                    ! delta(bigOM-om2+om3)
                    ii = max(floor((om2 - om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((om2 - om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) + plf2*lo_gauss(se%energy_axis(i), om2 - om3, sigma)
                    end do
                    ! delta(bigOM+om2-om3)
                    ii = max(floor((-om2 + om3 - 4*sigma)*invf), 1)
                    jj = min(ceiling((-om2 + om3 + 4*sigma)*invf), se%n_energy)
                    ilo = min(ilo, ii)
                    ihi = max(ihi, jj)
                    do i = ii, jj
                        buf(i) = buf(i) - plf2*lo_gauss(se%energy_axis(i), -om2 + om3, sigma)
                    end do

                    ! Increment the self-energy
                    if (ilo .lt. ihi) then
                        if (b2 .eq. b3) then
                            psisquare = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))*pref
                            se%im_3ph(ilo:ihi, b1) = se%im_3ph(ilo:ihi, b1) + buf(ilo:ihi)*psisquare
                        else
                            psisq23 = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))*pref
                            psisq32 = abs(sr%psi_3ph(b1, b3, b2, q)*conjg(sr%psi_3ph(b1, b3, b2, q)))*pref
                            se%im_3ph(ilo:ihi, b1) = se%im_3ph(ilo:ihi, b1) + buf(ilo:ihi)*(psisq23 + psisq32)
                        end if
                    end if
                end do
                end do

                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(ctr, 127, qp%n_irr_point*dr%n_mode)) call lo_progressbar(' ... threephonon imaginary selfenergy', ctr, dr%n_mode*qp%n_irr_point)
                end if
            end do
        end do qloopirr
    else
        ! If not Gamma we have to do all of it.
        qloopfull: do q = 1, qp%n_full_point
            do b1 = 1, dr%n_mode
                ! Make it parallel
                ctr = ctr + 1
                if (mod(ctr, mw%n) .ne. mw%r) cycle

                ! Prefactor for this q-point
                pref = threephonon_prefactor*qp%ap(q)%integration_weight
                ! Slightly smarter version that should be a little faster.
                do b2 = 1, dr%n_mode
                    !do b3=b2,dr%n_mode
                    do b3 = 1, dr%n_mode
                        ! reset buffer
                        buf = 0.0_r8
                        ilo = lo_hugeint
                        ihi = -lo_hugeint
                        ! Get the smearing parameter, in case it's adaptive
                        select case (se%integrationtype)
                        case (1)
                            ! s2=dr%default_smearing(b2)
                            ! s3=dr%default_smearing(b3)
                            ! sigma=sqrt(s2**2 + s3**2)
                            sigma = 1.0_r8*lo_frequency_THz_to_Hartree
                        case (2)
                            s2 = qp%adaptive_sigma(qp%ap(q)%radius, dr%aq(q)%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                            s3 = qp%adaptive_sigma(qp%ap(q)%radius, sr%vel3(:, b3, q), dr%default_smearing(b3), se%smearing_prefactor)
                            sigma = sqrt(s2**2 + s3**2)
                        end select
                        ! fetch frequencies
                        om1 = wp%omega(b1)
                        om2 = dr%aq(q)%omega(b2)
                        om3 = sr%omega3(b3, q)
                        ! Get the S-function
                        plf1 = 1.0_r8 + lo_planck(temperature, om2) + lo_planck(temperature, om3)
                        plf2 = lo_planck(temperature, om3) - lo_planck(temperature, om2)

                        ! delta(bigOM-om2-om3)
                        ii = max(floor((om2 + om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((om2 + om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) + plf1*lo_gauss(se%energy_axis(i), om2 + om3, sigma)
                        end do
                        ! delta(bigOM+om2+om3)
                        ii = max(floor((-om2 - om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((-om2 - om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) - plf1*lo_gauss(se%energy_axis(i), -om2 - om3, sigma)
                        end do
                        ! delta(bigOM-om2+om3)
                        ii = max(floor((om2 - om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((om2 - om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) + plf2*lo_gauss(se%energy_axis(i), om2 - om3, sigma)
                        end do
                        ! delta(bigOM+om2-om3)
                        ii = max(floor((-om2 + om3 - 4*sigma)*invf), 1)
                        jj = min(ceiling((-om2 + om3 + 4*sigma)*invf), se%n_energy)
                        ilo = min(ilo, ii)
                        ihi = max(ihi, jj)
                        do i = ii, jj
                            buf(i) = buf(i) - plf2*lo_gauss(se%energy_axis(i), -om2 + om3, sigma)
                        end do

                        ! Increment the self-energy
                        if (ilo .lt. ihi) then
                            !if ( b2 .eq. b3 ) then
                            psisquare = abs(sr%psi_3ph(b1, b2, b3, q)*conjg(sr%psi_3ph(b1, b2, b3, q)))*pref
                            se%im_3ph(ilo:ihi, b1) = se%im_3ph(ilo:ihi, b1) + buf(ilo:ihi)*psisquare
                            !else
                            !    psisq23=abs( sr%psi_3ph(b1,b2,b3,q)*conjg(sr%psi_3ph(b1,b2,b3,q)) )*pref
                            !    psisq32=abs( sr%psi_3ph(b1,b3,b2,q)*conjg(sr%psi_3ph(b1,b3,b2,q)) )*pref
                            !    se%im_3ph(ilo:ihi,b1)=se%im_3ph(ilo:ihi,b1)+buf(ilo:ihi)*(psisq23+psisq32)
                            !endif
                        end if
                    end do
                end do

                if (verbosity .gt. 0) then
                    if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... threephonon imaginary selfenergy', ctr, dr%n_mode*qp%n_full_point)
                end if
            end do
        end do qloopfull
    end if

    ! Sum it up
    call mw%allreduce('sum', se%im_3ph)

    ! Fix degeneracies?
    do b1 = 1, dr%n_mode
        buf = 0.0_r8
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            buf = buf + se%im_3ph(:, b2)
        end do
        buf = buf/real(wp%degeneracy(b1), r8)
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            se%im_3ph(:, b2) = buf
        end do
    end do
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    if (verbosity .gt. 0) call lo_progressbar(' ... threephonon imaginary selfenergy', dr%n_mode*qp%n_full_point, dr%n_mode*qp%n_full_point, walltime() - t0)
end subroutine


!> fourphonon self energy gaussian
subroutine fourphonon_imaginary_selfenergy_gaussian(wp, se, fc, fcf, qp, dr, gpoint, qpoint, uc, closedgrid, temperature, mw, mem, verbosity)
    !> harmonic properties at this q-point
    type(lo_phonon_dispersions_qpoint), intent(in) :: wp
    !> self-energy
    type(lo_phonon_selfenergy), intent(inout) :: se
    !> second order forceconstant
    type(lo_forceconstant_secondorder), intent(inout) :: fc
    !> fourth order forceconstant
    type(lo_forceconstant_fourthorder), intent(in) :: fcf
    !> qpoint mesh
    class(lo_qpoint_mesh), intent(in) :: qp
    !> harmonic dispersions
    type(lo_phonon_dispersions), intent(in) :: dr
    !> harmonic properties at Gamma with the right direction
    type(lo_phonon_dispersions_qpoint), intent(in) :: gpoint
    !> q-point in question
    type(lo_qpoint), intent(in) :: qpoint
    !> crystal structure
    type(lo_crystalstructure), intent(in) :: uc
    !> is it a closed grid ?
    logical, intent(in) :: closedgrid
    !> temperature
    real(r8), intent(in) :: temperature
    !> MPI helper
    type(lo_mpi_helper), intent(inout) :: mw
    !> memory tracker
    type(lo_mem_helper), intent(inout) :: mem
    !> talk a lot?
    integer, intent(in) :: verbosity

    real(r8), dimension(:), allocatable :: buf
    complex(r8), dimension(:, :), allocatable :: egv

    complex(r8) :: psisq
    real(r8), dimension(3) :: qv1, qv2, qv3, qv4
    real(r8), dimension(4) :: omega
    real(r8) :: sigma, omdiff, invf, s2, s3, s4, t0, t1
    real(r8) :: plf1, plf2, pref, n2, n3, n4
    real(r8) :: om1, om2, om3, om4
    integer :: q2, b1, b2, b3, b4, ctr, i, ii, jj, ilo, ihi

    ! set the unit
    t0 = walltime()

    invf = se%n_energy/se%energy_axis(se%n_energy)

    if (verbosity .gt. 0) call lo_progressbar_init()
    omega = 0.0_r8
    ctr = 0
    se%im_4ph = 0.0_r8
    call mem%allocate(buf, se%n_energy, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)
    call mem%allocate(egv,[dr%n_mode,4],persistent=.false.,scalable=.false.,file=__FILE__,line=__LINE__)
    buf = 0.0_r8
    egv=0.0_r8

    do q2=1, qp%n_full_point
        pref = fourphonon_imag_prefactor * qp%ap(q2)%integration_weight
        qv1 = qpoint%r * lo_twopi
        qv2 = -qv1
        qv3 = qp%ap(q2)%r * lo_twopi
        qv4 = -qv3
        do b1 = 1, dr%n_mode
        ctr = ctr + 1
        if (mod(ctr, mw%n) .ne. mw%r) cycle

        om1 = wp%omega(b1)
        if (om1 .lt. lo_freqtol) cycle
        do b2 = 1, dr%n_mode
        om2 = wp%omega(b2)
        if (om2 .lt. lo_freqtol) cycle
        do b3 = 1, dr%n_mode
        om3 = dr%aq(q2)%omega(b3)
        if (om3 .lt. lo_freqtol) cycle
        do b4 = b3, dr%n_mode
        om4 = dr%aq(q2)%omega(b4)
        if (om4 .lt. lo_freqtol) cycle

            select case (se%integrationtype)
            case (1)
                sigma = 1.0_r8*lo_frequency_THz_to_Hartree
            case (2)
                !s2 = qp%adaptive_sigma(qpoint%radius, qpoint%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                s2 = qp%adaptive_sigma(qp%ap(q2)%radius, wp%vel(:, b2), dr%default_smearing(b2), se%smearing_prefactor)
                s3 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b3), dr%default_smearing(b3), se%smearing_prefactor)
                s4 = qp%adaptive_sigma(qp%ap(q2)%radius, dr%aq(q2)%vel(:, b4), dr%default_smearing(b4), se%smearing_prefactor)
                sigma = sqrt(s2**2 + s3**2 + s4**2)
            end select

            n2 = lo_planck(temperature, om2)
            n3 = lo_planck(temperature, om3)
            n4 = lo_planck(temperature, om4)

            ilo = lo_hugeint
            ihi = -lo_hugeint
            buf = 0.0_r8

                plf1 = (n2 + 1) * (n3 + 1) * (n4 + 1)
                plf2 = n2 * n3 * n4
                omdiff = om2 +  om3 + om4
                ii = max(floor((omdiff - 4 * sigma)*invf), 1)
                jj = min(ceiling((omdiff + 4 * sigma)*invf), se%n_energy)
                ilo = min(ilo, ii)
                ihi = max(ihi, jj)
                do i = ii, jj
                    buf(i) = buf(i) + (plf1 - plf2) * lo_gauss(se%energy_axis(i), omdiff, sigma)
                end do
                omdiff = -om2 - om3 - om4
                ii = max(floor((omdiff - 4 * sigma)*invf), 1)
                jj = min(ceiling((omdiff + 4 * sigma)*invf), se%n_energy)
                ilo = min(ilo, ii)
                ihi = max(ihi, jj)
                do i = ii, jj
                    buf(i) = buf(i) - (plf1 - plf2) * lo_gauss(se%energy_axis(i), omdiff, sigma)
                end do

                plf1 = n2 * (n3 + 1) * (n4 + 1)
                plf2 = (n2 + 1) * n3 * n4
                omdiff = om2 - om3 - om4
                ii = max(floor((omdiff - 4 * sigma)*invf), 1)
                jj = min(ceiling((omdiff + 4 * sigma)*invf), se%n_energy)
                ilo = min(ilo, ii)
                ihi = max(ihi, jj)
                do i = ii, jj
                    buf(i) = buf(i) - 3 * (plf1 - plf2) * lo_gauss(se%energy_axis(i), omdiff, sigma)
                end do
                omdiff = -om2 + om3 + om4
                ii = max(floor((omdiff - 4 * sigma)*invf), 1)
                jj = min(ceiling((omdiff + 4 * sigma)*invf), se%n_energy)
                ilo = min(ilo, ii)
                ihi = max(ihi, jj)
                do i = ii, jj
                    buf(i) = buf(i) + 3 * (plf1 - plf2) * lo_gauss(se%energy_axis(i), omdiff, sigma)
                end do

            if (ilo .lt. ihi) then
                omega(1) = om1
                omega(2) = om2
                omega(3) = om3
                omega(4) = om4

                egv(:,1)=wp%egv(:,b1)
                egv(:,2)=conjg(wp%egv(:,b2))
                egv(:,3)=dr%aq(q2)%egv(:,b3)
                egv(:,4)=conjg(dr%aq(q2)%egv(:,b4))

                psisq = fcf%scatteringamplitude(omega, egv, qv2, qv3, qv4)
                if (b3 .eq. b4) then
                    se%im_4ph(ilo:ihi, b1) = se%im_4ph(ilo:ihi, b1) + buf(ilo:ihi)*abs(psisq*conjg(psisq))*pref
                else
                    se%im_4ph(ilo:ihi, b1) = se%im_4ph(ilo:ihi, b1) + 2 * buf(ilo:ihi)*abs(psisq*conjg(psisq))*pref
                end if
            end if
        end do
        end do
        end do
        end do
        if (verbosity .gt. 0) then
            if (lo_trueNtimes(ctr, 127, qp%n_full_point*dr%n_mode)) call lo_progressbar(' ... fourphonon imaginary selfenergy', ctr, &
            qp%n_full_point*dr%n_mode)
        end if
    end do
    call mem%deallocate(egv, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    ! Sum it up
    call mw%allreduce('sum', se%im_4ph)

    ! Fix degeneracies ?
    do b1 = 1, dr%n_mode
        buf = 0.0_r8
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            buf = buf + se%im_4ph(:, b2)
        end do
        buf = buf / real(wp%degeneracy(b1), r8)
        do i = 1, wp%degeneracy(b1)
            b2 = wp%degenmode(i, b1)
            se%im_4ph(:, b2) = buf
        end do
    end do
    call mem%deallocate(buf, persistent=.false., scalable=.false., file=__FILE__, line=__LINE__)

    if (verbosity .gt. 0) call lo_progressbar(' ... fourphonon imaginary selfenergy', qp%n_full_point*dr%n_mode,qp%n_full_point*dr%n_mode, walltime() - t0)
end subroutine
