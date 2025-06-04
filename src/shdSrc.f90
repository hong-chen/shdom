module shd_src

    use globals
    use shd_pha
    use shd_util
    use shd_init

    implicit none

    public :: shd_src__computeSrcOnePnt
    public :: shd_src__computeSrcOnePnt_
    public :: shd_src__computeSrcOneCell
    public :: shd_src__computeSrcOneCell_
    public :: shd_src__computeSrcOneDir
    public :: shd_src__computeSrc

    public :: shd_src__integrateSrcOneRay
    public :: shd_src__integrateSrc

    private :: shd_src__computeSolarStokesVect
    private :: shd_src__rotatePolarPlane
    private :: shd_src__findBoundaryRad

contains

subroutine shd_src__computeSrcOneCell(&
    &ICELL, GRIDPTR,&
    &NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,&
    &NPTS, DELTAM, SRCTYPE, SOLARMU,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,&
    &SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
    &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
    &EXTINCT8, SRCEXT8)

    !   Computes the source function times extinction for gridpoints
    ! belonging to cell ICELL in the direction (MU,PHI).  The results
    ! are returned in SRCEXT8 and EXTINCT8.
    ! The spherical harmonic source function series is input in SOURCE.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    integer ICELL, NSTOKES, NSTLEG, NPTS, ML, MM, NLM, NLEG, NUMPHASE
    integer GRIDPTR(8, *), SHPTR(NPTS + 1)
    integer DONETHIS(8), OLDIPTS(8)
    integer*2 IPHASE(NPTS)
    logical DELTAM
    real SOLARMU
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real DIRFLUX(NPTS), SOURCE(NSTOKES, *)
    real YLMDIR(NSTLEG, NLM), YLMSUN(NSTLEG, NLM)
    real SINGSCAT(NSTOKES, NUMPHASE)
    real OEXTINCT8(8), OSRCEXT8(NSTOKES, 8)
    real EXTINCT8(8), SRCEXT8(NSTOKES, 8)
    real(RD_) :: SUNDIRLEG(0:NLEG)
    character(len=1) :: SRCTYPE

    integer IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I
    real SECMU0, F, DA, A, A1, B1

    SECMU0 = 1.0d0/abs(SOLARMU)

    ! Loop over the grid points, computing the source function
    !   at the viewing angle from the spherical harmonics source fun
    do N = 1, 8
        IP = GRIDPTR(N, ICELL)
        I = DONETHIS(N)
        if (I .gt. 0 .and. IP .eq. OLDIPTS(N)) then
            EXTINCT8(N) = OEXTINCT8(I)
            SRCEXT8(:, N) = OSRCEXT8(:, I)
        else if (I .lt. 0) then
            EXTINCT8(N) = EXTINCT8(abs(I))
            SRCEXT8(:, N) = SRCEXT8(:, abs(I))
        else

            OLDIPTS(N) = IP
            IS = SHPTR(IP)
            NS = SHPTR(IP + 1) - IS
            ! Sum over the real generalized spherical harmonic series
            !   of the source function
            SRCEXT8(:, N) = 0.0
            do J = 1, NS
                SRCEXT8(1, N) = SRCEXT8(1, N) + SOURCE(1, IS + J)*YLMDIR(1, J)
            end do
            if (NSTOKES .gt. 1) then
                do J = 1, NS
                    SRCEXT8(2, N) = SRCEXT8(2, N) + SOURCE(2, IS + J)*YLMDIR(2, J)&
                        &+ SOURCE(3, IS + J)*YLMDIR(5, J)
                    SRCEXT8(3, N) = SRCEXT8(3, N) + SOURCE(2, IS + J)*YLMDIR(6, J)&
                        &+ SOURCE(3, IS + J)*YLMDIR(3, J)
                end do
            end if
            if (NSTOKES .eq. 4) then
                do J = 1, NS
                    SRCEXT8(4, N) = SRCEXT8(4, N) + SOURCE(4, IS + J)*YLMDIR(4, J)
                end do
            end if

            ! Special case for solar source and Delta-M
            if (SRCTYPE .ne. 'T' .and. DELTAM) then
                if (NUMPHASE .gt. 0) then
                    IPH = IPHASE(IP)
                else
                    IPH = IP
                end if
                ! First subtract off the truncated single scattering
                DA = ALBEDO(IP)*DIRFLUX(IP)*SECMU0
                J = 1
                do L = 0, ML
                    ME = min(L, MM)
                    MS = -ME
                    A1 = DA*LEGEN(1, L, IPH)
                    B1 = DA*LEGEN(5, L, IPH)
                    if (J .le. NS) then
                        JT = J
                        do M = MS, ME
                            SRCEXT8(1, N) = SRCEXT8(1, N) - A1*YLMDIR(1, J)*YLMSUN(1, J)
                            J = J + 1
                        end do
                        if (NSTOKES .gt. 1) then
                            J = JT
                            do M = MS, ME
                                SRCEXT8(2, N) = SRCEXT8(2, N) - B1*YLMDIR(2, J)*YLMSUN(1, J)
                                SRCEXT8(3, N) = SRCEXT8(3, N) - B1*YLMDIR(6, J)*YLMSUN(1, J)
                                J = J + 1
                            end do
                        end if
                    end if
                end do
                ! Then add in the single scattering contribution for the
                !   original unscaled phase function.
                if (NUMPHASE .gt. 0) then
                    SRCEXT8(:, N) = SRCEXT8(:, N) + DA*SINGSCAT(:, IPH)
                else if (NSTOKES .eq. 1) then
                    F = LEGEN(1, ML + 1, IPH)
                    do L = 0, NLEG
                        if (L .le. ML) then
                            A = DA*(LEGEN(1, L, IPH) + F/(1 - F))
                        else
                            A = DA*LEGEN(1, L, IPH)/(1 - F)
                        end if
                        SRCEXT8(1, N) = SRCEXT8(1, N) + A*SUNDIRLEG(L)
                    end do
                end if
            end if

            SRCEXT8(:, N) = SRCEXT8(:, N)*EXTINCT(IP)
            EXTINCT8(N) = EXTINCT(IP)
        end if
    end do

    return
end subroutine shd_src__computeSrcOneCell

subroutine shd_src__computeSrcOneCell_(&
    &ICELL, GRIDPTR,&
    &ML, MM, NLM, NLEG, NUMPHASE,&
    &NPTS, DELTAM, SRCTYPE, SOLARMU,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,&
    &SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
    &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
    &EXTINCT8, SRCEXT8)

    !   Computes the source function times extinction for gridpoints
    ! belonging to cell ICELL in the direction (MU,PHI).  The results
    ! are returned in SRCEXT8 and EXTINCT8.
    ! The spherical harmonic source function series is input in SOURCE.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    integer ICELL, NPTS, ML, MM, NLM, NLEG, NUMPHASE
    integer GRIDPTR(8, *), SHPTR(NPTS + 1)
    integer DONETHIS(8), OLDIPTS(8)
    integer*2 IPHASE(NPTS)
    logical DELTAM
    real SOLARMU
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(0:NLEG, NPTS)
    real DIRFLUX(NPTS), SOURCE(*)
    real YLMDIR(NLM), YLMSUN(NLM)
    real SINGSCAT(NUMPHASE)
    real OEXTINCT8(8), OSRCEXT8(8)
    real EXTINCT8(8), SRCEXT8(8)
    real(RD_) :: SUNDIRLEG(0:NLEG)
    character(len=1) :: SRCTYPE

    integer IP, J, JT, L, M, MS, ME, IPH, IS, NS, N, I
    real SECMU0, F, DA, A, A1, B1

    SECMU0 = 1.0d0/abs(SOLARMU)

    ! Loop over the grid points, computing the source function
    !   at the viewing angle from the spherical harmonics source fun
    do N = 1, 8
        IP = GRIDPTR(N, ICELL)
        I = DONETHIS(N)
        if (I .gt. 0 .and. IP .eq. OLDIPTS(N)) then
            EXTINCT8(N) = OEXTINCT8(I)
            SRCEXT8(N) = OSRCEXT8(I)
        else if (I .lt. 0) then
            EXTINCT8(N) = EXTINCT8(abs(I))
            SRCEXT8(N) = SRCEXT8(abs(I))
        else

            OLDIPTS(N) = IP
            IS = SHPTR(IP)
            NS = SHPTR(IP + 1) - IS
            ! Sum over the real generalized spherical harmonic series
            !   of the source function
            SRCEXT8(N) = 0.0
            do J = 1, NS
                SRCEXT8(N) = SRCEXT8(N) + SOURCE(IS + J)*YLMDIR(J)
            end do

            ! Special case for solar source and Delta-M
            if (SRCTYPE .ne. 'T' .and. DELTAM) then
                if (NUMPHASE .gt. 0) then
                    IPH = IPHASE(IP)
                else
                    IPH = IP
                end if
                ! First subtract off the truncated single scattering
                DA = ALBEDO(IP)*DIRFLUX(IP)*SECMU0
                J = 1
                do L = 0, ML
                    ME = min(L, MM)
                    MS = -ME
                    A1 = DA*LEGEN(L, IPH)
                    if (J .le. NS) then
                        JT = J
                        do M = MS, ME
                            SRCEXT8(N) = SRCEXT8(N) - A1*YLMDIR(J)*YLMSUN(J)
                            J = J + 1
                        end do
                    end if
                end do
                ! Then add in the single scattering contribution for the
                !   original unscaled phase function.
                if (NUMPHASE .gt. 0) then
                    SRCEXT8(N) = SRCEXT8(N) + DA*SINGSCAT(IPH)
                else
                    F = LEGEN(ML + 1, IPH)
                    do L = 0, NLEG
                        if (L .le. ML) then
                            A = DA*(LEGEN(L, IPH) + F/(1 - F))
                        else
                            A = DA*LEGEN(L, IPH)/(1 - F)
                        end if
                        SRCEXT8(N) = SRCEXT8(N) + A*SUNDIRLEG(L)
                    end do
                end if
            end if

            SRCEXT8(N) = SRCEXT8(N)*EXTINCT(IP)
            EXTINCT8(N) = EXTINCT(IP)
        end if
    end do

    return
end subroutine shd_src__computeSrcOneCell_

subroutine shd_src__computeSrcOneDir(NSTOKES,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NPTS, DELTAM, MU, PHI, SRCTYPE, SOLARMU, SOLARAZ,&
    &ALBEDO, LEGEN, IPHASE, DIRFLUX,&
    &SHPTR, SOURCE, SOURCE1)

    ! Computes the source function (SOURCE1) in the direction (MU,PHI)
    ! for the whole domain (NX,NZ).  The spherical harmonic source function
    ! series is input in SOURCE.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    integer NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer SHPTR(NPTS + 1)
    integer*2 IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ, MU, PHI
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS), DIRFLUX(NPTS)
    real SOURCE(NSTOKES, *), SOURCE1(NSTOKES, NPTS)
    character(len=1) :: SRCTYPE
    integer I, J, JT, K, L, M, N, IPH, IS, NS, ME, MS
    real DA, F, A, A1, B1, SECMU0, COSSCAT
    real, allocatable :: YLMDIR(:, :), YLMSUN(:, :)
    real(RD_), allocatable :: SUNDIRLEG(:)
    real, allocatable :: UNSCLEGEN(:, :), SINGSCAT(:, :)
    logical TMS

    allocate (YLMDIR(NSTLEG, NLM), YLMSUN(NSTLEG, NLM))
    allocate (SUNDIRLEG(0:NLEG))
    allocate (UNSCLEGEN(2, 0:NLEG), SINGSCAT(NSTOKES, NUMPHASE))
    if (abs(MU) .gt. 1.0) stop 'Error [_computeOneSrc]: Bad mu.'

    TMS = DELTAM
    ! TMS = .FALSE.

    ! Precompute Ylm's for output direction and solar direction
    call shd_pha__computeYLM(.false., MU, PHI, ML, MM, NSTLEG, YLMDIR)
    if (SRCTYPE .ne. 'T') then
        call shd_pha__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG,&
            &YLMSUN)
        SECMU0 = 1.0/abs(SOLARMU)
        COSSCAT = SOLARMU*MU + sqrt((1.0d0 - SOLARMU**2)*(1.0d0 - MU**2))&
            &*cos(SOLARAZ - PHI)
        if (NUMPHASE .eq. 0 .and. TMS) then
            ! If non-tabulate phase functions (unpolarized) then
            !   compute the Legendre polynomials for the scattering angle
            !   for the untruncated solar single scattering computation.
            call shd_pha__computeLegendre(dble(COSSCAT), NLEG, SUNDIRLEG)
            do L = 0, NLEG
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)*O4PI_
            end do
        end if
    end if

    ! If solar source and tabulated phase functions then precompute
    !   the single scattered solar Stokes vector to outgoing direction
    if (SRCTYPE .ne. 'T' .and. TMS .and. NUMPHASE .gt. 0) then
        do IPH = 1, NUMPHASE
            ! Unscaled the needed Legendre coefficients, but divide by 1-f
            !   because extinction is still scaled (for TMS method)
            F = LEGEN(1, ML + 1, IPH)
            do L = 0, NLEG
                if (L .le. ML) then
                    UNSCLEGEN(1, L) = LEGEN(1, L, IPH) + F/(1 - F)
                else
                    UNSCLEGEN(1, L) = LEGEN(1, L, IPH)/(1 - F)
                end if
                if (NSTLEG .gt. 1) then
                    if (L .le. ML) then
                        UNSCLEGEN(2, L) = LEGEN(5, L, IPH)
                    else
                        UNSCLEGEN(2, L) = LEGEN(5, L, IPH)/(1 - F)
                    end if
                end if
            end do
            ! Sum the Wigner series (for A1 and B1), rotate the plane of
            !   polarization to the outgoing direction, and return the Stoke
            call shd_src__computeSolarStokesVect(NLEG, UNSCLEGEN, COSSCAT,&
                &SOLARMU, MU, SOLARAZ - PHI, NSTOKES,&
                &SINGSCAT(1, IPH))
            if (SINGSCAT(1, IPH) .le. 0.0) then
                write (6, *) 'Error [_computeOneSrc]: negative source function',&
                    &' for tabulated phase function:', IPH, MU, PHI, SINGSCAT(1, IPH)
                stop
            end if
        end do
    end if

    ! Loop over all the grid points, computing the source function
    !   at the viewing angle from the spherical harmonic source func
    do I = 1, NPTS
        IS = SHPTR(I)
        NS = SHPTR(I + 1) - IS
        ! Sum over the real generalized spherical harmonic series
        !   of the source function
        SOURCE1(:, I) = 0.0
        do J = 1, NS
            SOURCE1(1, I) = SOURCE1(1, I) + SOURCE(1, IS + J)*YLMDIR(1, J)
        end do
        if (NSTOKES .gt. 1) then
            do J = 1, NS
                SOURCE1(2, I) = SOURCE1(2, I) + SOURCE(2, IS + J)*YLMDIR(2, J)&
                    &+ SOURCE(3, IS + J)*YLMDIR(5, J)
                SOURCE1(3, I) = SOURCE1(3, I) + SOURCE(2, IS + J)*YLMDIR(6, J)&
                    &+ SOURCE(3, IS + J)*YLMDIR(3, J)
            end do
        end if
        if (NSTOKES .eq. 4) then
            do J = 1, NS
                SOURCE1(4, I) = SOURCE1(4, I) + SOURCE(4, IS + J)*YLMDIR(4, J)
            end do
        end if

        ! Special case for solar source and Delta-M
        if (SRCTYPE .ne. 'T' .and. TMS) then
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            ! First subtract off the truncated single scattering
            DA = ALBEDO(I)*DIRFLUX(I)*SECMU0
            J = 1
            do L = 0, ML
                ME = min(L, MM)
                MS = -ME
                A1 = DA*LEGEN(1, L, IPH)
                B1 = DA*LEGEN(5, L, IPH)
                if (J .le. NS) then
                    JT = J
                    do M = MS, ME
                        SOURCE1(1, I) = SOURCE1(1, I) - A1*YLMDIR(1, J)*YLMSUN(1, J)
                        J = J + 1
                    end do
                    if (NSTOKES .gt. 1) then
                        J = JT
                        do M = MS, ME
                            SOURCE1(2, I) = SOURCE1(2, I) - B1*YLMDIR(2, J)*YLMSUN(1, J)
                            SOURCE1(3, I) = SOURCE1(3, I) - B1*YLMDIR(6, J)*YLMSUN(1, J)
                            J = J + 1
                        end do
                    end if
                end if
            end do

            ! Then add in the single scattering contribution for the
            !   original unscaled phase matrix.
            if (NUMPHASE .gt. 0) then
                SOURCE1(:, I) = SOURCE1(:, I) + DA*SINGSCAT(:, IPH)
            else if (NSTOKES .eq. 1) then
                ! If we have a phase function for each grid point (which
                !   should be the unpolarized case) then unscale the Legendre
                !   coefficients for L<=ML and sum over the Legendre series
                !   in cosine of the scattering angle.
                F = LEGEN(1, ML + 1, IPH)
                do L = 0, NLEG
                    if (L .le. ML) then
                        A = DA*(LEGEN(1, L, IPH) + F/(1 - F))
                    else
                        A = DA*LEGEN(1, L, IPH)/(1 - F)
                    end if
                    SOURCE1(1, I) = SOURCE1(1, I) + A*SUNDIRLEG(L)
                end do
            end if
        end if
    end do

    deallocate (YLMDIR, YLMSUN, SUNDIRLEG, UNSCLEGEN, SINGSCAT)

    return

end subroutine shd_src__computeSrcOneDir

subroutine shd_src__computeSrc(&
    &NSTOKES, ML, MM, NLM,&
    &NSTLEG, NLEG, NUMPHASE, NPTS, FIXSH,&
    &SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
    &IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
    &FIRST, ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

    ! Computes the source function (SOURCE) in spherical harmonic space
    !   for all the grid points.  The thermal source and/or solar
    !   pseudo-source (in PLANCK or DIRFLUX) is added to the scattering source
    !   (computed from LEGEN and the spherical harmonic expansion in RADIANCE).
    !   Also computes the adaptive spherical harmonic truncation for
    !   the source function, unless FIXSH is true.  The SH truncation is set
    !   so that the higher terms of the source function are smaller in
    !   absolute value than SHACC. The truncation is to the nearest l,
    !   so that a triangular (for M=L) truncation is preserved.  If there
    !   is little or no scattering then the SH truncation can be set to
    !   zero terms.  The source function SH truncation is output in terms
    !   of the pointer SHPTR (note that RSHPTR refers to the RADIANCE array
    !   and OSHPTR to the DELSOURCE array).
    ! To save memory this routine also computes some things for the
    !   series acceleration: the dot products of the difference in source
    !   function between successive iterations and the new difference source
    !   function vector (DELSOURCE).  To avoid using another array and
    !   because of the adaptive truncation pointer stuff, this requires
    !   two extra calculations of the source function.
    ! Computes the sum of squares of the new source function for calculation
    !   of the solution criterion (SOLCRIT, the RMS difference in succesive
    !   source function fields normalized by the RMS of the field).

    integer NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG
    integer NUMPHASE, MAXIV
    integer RSHPTR(NPTS + 1), SHPTR(NPTS + 1), OSHPTR(NPTS + 1)
    integer*2 IPHASE(NPTS)
    logical FIXSH, FIRST, ACCELFLAG
    real SHACC, YLMSUN(NSTLEG, NLM), SOLARMU
    real DELJDOT, DELJOLD, DELJNEW, JNORM
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real PLANCK(NPTS), DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *)
    real SOURCE(NSTOKES, *), DELSOURCE(NSTOKES, *)
    character(len=1) :: SRCTYPE
    integer IS, ISO, IR, I, IPH, J, JS, K, L, LS, M, MS, ME, NS, NR
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)
    real SRCMIN, C, SECMU0, D

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))

    SRCMIN = SHACC
    ! Make source function from scattering integral and real sources:
    !   Solar source: DIRFLUX has F0*exp(-tau)
    !   Thermal source: PLANCK has (1-omega)*B(T)
    SECMU0 = 1.0/abs(SOLARMU)
    C = sqrt(4.0*acos(-1.0))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    ! Compute the dot products of the successive source differences
    if (.not. FIRST) then
        DELJDOT = 0.0
        DELJOLD = 0.0
        DELJNEW = 0.0
        JNORM = 0.0
        do I = 1, NPTS
            ! Compute the temporary source function for this point
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            IR = RSHPTR(I)
            NR = RSHPTR(I + 1) - IR
            if (NR .gt. NLM) then
                write (6, *) 'Error [_computeSrc]: NR>NLM 1', I, IR, RSHPTR(I + 1)
                stop
            end if
            if (NSTOKES .eq. 1) then
                call shd_src__computeSrcOnePnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_src__computeSrcOnePnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Compute the dot product of the new and old source difference,
            !   and the solution criterion (don't use new points for Jnorm).
            if (ACCELFLAG) then
                IS = SHPTR(I)
                ISO = OSHPTR(I)
                NS = min(SHPTR(I + 1) - IS, OSHPTR(I + 1) - ISO)
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJDOT = DELJDOT +&
                            &(SOURCET(K, J) - SOURCE(K, IS + J))*DELSOURCE(K, ISO + J)
                        DELJOLD = DELJOLD + DELSOURCE(K, ISO + J)**2
                        DELJNEW = DELJNEW + (SOURCET(K, J) - SOURCE(K, IS + J))**2
                        JNORM = JNORM + SOURCE(K, IS + J)**2
                    end do
                end do
            else
                IS = SHPTR(I)
                NS = SHPTR(I + 1) - IS
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJNEW = DELJNEW + (SOURCET(K, J) - SOURCE(K, IS + J))**2
                        JNORM = JNORM + SOURCE(K, IS + J)**2
                    end do
                end do
            end if
        end do
    end if

    ! Compute and store the new source function difference vector
    if (.not. FIRST .and. ACCELFLAG) then
        do I = 1, NPTS
            ! Compute the temporary source function for this point
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            IR = RSHPTR(I)
            NR = RSHPTR(I + 1) - IR
            if (NSTOKES .eq. 1) then
                call shd_src__computeSrcOnePnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_src__computeSrcOnePnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Make the source difference vector
            IS = SHPTR(I)
            NS = SHPTR(I + 1) - IS
            OSHPTR(I) = IS
            do J = 1, NS
                DELSOURCE(:, IS + J) = SOURCET(:, J) - SOURCE(:, IS + J)
            end do
        end do
        OSHPTR(NPTS + 1) = SHPTR(NPTS + 1)
    end if

    ! Compute the source function in the temporary array.
    IS = 0
    do I = 1, NPTS
        if (NUMPHASE .gt. 0) then
            IPH = IPHASE(I)
        else
            IPH = I
        end if
        IR = RSHPTR(I)
        NR = RSHPTR(I + 1) - IR
        if (NR .gt. NLM) then
            write (6, *) 'Error [_computeSrc]: NR>NLM 3', I, IR, RSHPTR(I + 1)
            stop
        end if
        if (NSTOKES .eq. 1) then
            call shd_src__computeSrcOnePnt_(NLM, NLEG, LOFJ,&
                &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                &NR, RADIANCE(1, IR + 1), SOURCET)
        else
            call shd_src__computeSrcOnePnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                &NR, RADIANCE(1, IR + 1), SOURCET)
        end if
        if (FIXSH) then
            NS = SHPTR(I + 1) - IS
        else
            ! Then do the adaptive spherical harmonics:
            !   Find how many source function terms to use, and how many
            !   terms wanted for the radiance.
            if (SRCTYPE .eq. 'S') then
                JS = 0
            else
                JS = 1
            end if
            do J = 1, NLM
                do K = 1, NSTOKES
                    if (abs(SOURCET(K, J)) .gt. SRCMIN) JS = J
                end do
            end do
            ! Use all the m's for the last l for the source function.
            if (JS .eq. 0) then
                NS = 0
            else
                LS = LOFJ(JS)
                if (LS .le. MM) then
                    NS = (LS*(LS + 1)) + LS + 1
                else
                    NS = (2*MM + 1)*LS - (MM*(1 + (MM - 1))) + MM + 1
                end if
            end if
            SHPTR(I) = IS
        end if
        if (IS + NS .gt. MAXIV) then
            write (6, *) 'Error [_computeSrc]: MAXIV exceeded.', MAXIV
            write (6, *) '  Out of memory for more spherical harmonic terms.'
            stop
        end if
        ! Transfer the new source function
        do J = 1, NS
            SOURCE(:, IS + J) = SOURCET(:, J)
        end do
        IS = IS + NS
    end do
    SHPTR(NPTS + 1) = IS
    deallocate (LOFJ, SOURCET)
    return
end subroutine shd_src__computeSrc

subroutine shd_src__integrateSrcOneRay(&
    &BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
    &NSTPHASE, NSCATANGLE, PHASETAB,&
    &NX, NY, NZ, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &ML, MM, NLM, NLEG, NUMPHASE,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &DELTAM, SRCTYPE, WAVELEN, SOLARMU, SOLARAZ,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
    &DIRFLUX, SHPTR, SOURCE, YLMSUN,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &MU2, PHI2, X0, Y0, Z0,&
    &XE, YE, ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

    !   Integrates the source function through the extinction field
    ! (EXTINCT) backward from the outgoing direction (MU2,PHI2) to find
    ! radiance (RADIANCE) at the point X0,Y0,Z0.
    ! The transmission and radiance of the ray so far (TRANSMIT, RADIANC
    ! are input and returned after the integration along with the exitti
    ! ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+
    ! 5=-Z,6=+Z).

    integer BCFLAG, IPFLAG, NSTOKES, NSTLEG, NSTPHASE, NSCATANGLE
    integer NX, NY, NZ, NPTS, NCELLS, SIDE
    integer ML, MM, NLM, NLEG, NUMPHASE
    integer MAXNBC, NTOPPTS, NBOTPTS
    integer NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SHPTR(NPTS + 1)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    integer BCPTR(MAXNBC, 2)
    logical DELTAM, VALIDRAD
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, SOLARMU, SOLARAZ
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS), BCRAD(*)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real DIRFLUX(NPTS), SOURCE(NSTOKES, *)
    real YLMSUN(NSTLEG, NLM)
    real PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE)
    real(RD_) :: MU2, PHI2, X0, Y0, Z0, XE, YE, ZE
    real(RD_) :: TRANSMIT, RADIANCE(NSTOKES)
    character(len=1) :: SRCTYPE
    character(len=2) :: SFCTYPE

    integer BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
    integer IOPP, NTAU, IT, I, IPT1, IPT2, J, L
    logical DONE, IPINX, IPINY, OPENBCFACE
    integer JFACE, KFACE, IC, MAXCELLSCROSS, NGRID
    integer OPPFACE(6), OLDIPTS(8), DONETHIS(8)
    integer DONEFACE(8, 7), ONEX(8), ONEY(8)
    real OEXTINCT8(8), OSRCEXT8(NSTOKES, 8)
    real EXTINCT8(8), SRCEXT8(NSTOKES, 8)
    real EXT0, EXT1, EXTN
    real SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RADBND(NSTOKES)
    real XM, YM, F
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XN, YN, ZN, XI, YI, ZI
    real(RD_) :: SO, SOX, SOY, SOZ, EPS
    real(RD_) :: TAUTOL, TAUGRID, S, DELS, TRANSCUT
    real(RD_) :: EXT, TAU, TRANSCELL, ABSCELL
    real(RD_) :: SRC(NSTOKES)
    real(RD_) :: U, V, W, DELX, DELY, DELZ, INVDELX, INVDELY, INVDELZ
    real(RD_) :: COSSCAT
    real, allocatable :: YLMDIR(:, :), SINGSCAT(:, :)
    real(RD_), allocatable :: SUNDIRLEG(:)

    data OPPFACE/2, 1, 4, 3, 6, 5/
    data ONEY/0, 0, -1, -2, 0, 0, -5, -6/, ONEX/0, -1, 0, -3, 0, -5, 0, -7/
    data DONEFACE/0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 3, 0, 5, 0, 7, 2, 0, 4, 0, 6, 0, 8, 0,&
        &0, 0, 1, 2, 0, 0, 5, 6, 3, 4, 0, 0, 5, 6, 0, 0,&
        &0, 0, 0, 0, 1, 2, 3, 4, 5, 6, 7, 8, 0, 0, 0, 0/

    ! TRANSCUT is the transmission to stop the integration at
    TRANSCUT = 5.0e-5
    ! TAUTOL is the maximum optical path for the subgrid intervals
    TAUTOL = 0.2

    EPS = 1.0e-5*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))
    MAXCELLSCROSS = 50*max(NX, NY, NZ)

    ! Calculate the generalized spherical harmonics for this direction
    allocate (YLMDIR(NSTLEG, NLM))
    call shd_pha__computeYLM(.false., sngl(MU2), sngl(PHI2), ML, MM, NSTLEG, YLMDIR)

    if (SRCTYPE .ne. 'T' .and. DELTAM) then
        ! Get the solar single scattering Stokes vector for the outgoin
        !   direction my interpolating in scattering angle in the PHASETA
        !   table and then rotating the Q/U polarization to the outgoing
        COSSCAT = SOLARMU*MU2 + sqrt((1.0 - SOLARMU**2)*(1.0 - MU2**2))&
            &*cos(SOLARAZ - PHI2)
        if (NUMPHASE .gt. 0) then
            allocate (SINGSCAT(NSTOKES, NUMPHASE))
            F = (NSCATANGLE - 1)*(acos(COSSCAT)/DPI_) + 1
            J = min(NSCATANGLE - 1, int(F))
            F = F - J
            do I = 1, NUMPHASE
                SINGSCAT(1:NSTPHASE, I)&
                    &= (1 - F)*PHASETAB(:, I, J) + F*PHASETAB(:, I, J + 1)
                if (NSTOKES .gt. 1) then
                    call shd_src__rotatePolarPlane(NSTOKES, COSSCAT, SOLARMU,&
                        &sngl(MU2), SOLARAZ - sngl(PHI2), SINGSCAT(:, I))
                end if
            end do
        else
            allocate (SUNDIRLEG(0:NLEG))
            call shd_pha__computeLegendre(COSSCAT, NLEG, SUNDIRLEG)
            do L = 0, NLEG
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)/(4*DPI_)
            end do
        end if
    end if

    ! Make the ray direction (opposite to the outgoing direction)
    CX = sqrt(1.0d0 - MU2**2)*cos(PHI2 - DPI_)
    CY = sqrt(1.0d0 - MU2**2)*sin(PHI2 - DPI_)
    CZ = -MU2
    if (abs(CX) .gt. 1.0e-6) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-6) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    if (abs(CZ) .gt. 1.0e-6) then
        CZINV = 1.0d0/CZ
    else
        CZ = 0.0
        CZINV = 1.0e6
    end if

    ! Setup for the ray path direction
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. 0.0) then
        BITZ = 1
    else
        BITZ = 0
    end if
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    XM = 0.5d0*(XGRID(1) + XGRID(NX))
    YM = 0.5d0*(YGRID(1) + YGRID(NY))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    ! ZE = MAX(MIN(Z0,DBLE(ZGRID(NZ))),DBLE(ZGRID(1)))
    call shd_util__findCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
        &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
        &BCFLAG, IPFLAG, XE, YE, ZE, ICELL)
    ! print '(A,5(1X,F8.5),1X,I6)',&
    !     &'_integrateSrcOneRay:',X0,Y0,Z0,MU2,PHI2,ICELL
    IFACE = 0
    NGRID = 0

    ! Loop until reach a Z boundary or transmission is very small
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
        ! Make sure current cell is valid
        if (ICELL .le. 0) then
            write (6, *) 'Error [_integrateSrcOneRay]: ICELL=', ICELL,&
                &MU2, PHI2, XE, YE, ZE
            stop
        end if
        NGRID = NGRID + 1

        ! Decide which of the eight grid points we need the source fun
        do I = 1, 8
            DONETHIS(I) = DONEFACE(I, IFACE + 1)
            if (NX .eq. 1 .and. ONEX(I) .lt. 0) DONETHIS(I) = ONEX(I)
            if (NY .eq. 1 .and. ONEY(I) .lt. 0) DONETHIS(I) = ONEY(I)
            OEXTINCT8(I) = EXTINCT8(I)
            OSRCEXT8(:, I) = SRCEXT8(:, I)
        end do
        ! Compute the source function times extinction in direction (MU2
        if (NSTOKES .eq. 1) then
            call shd_src__computeSrcOneCell_(ICELL, GRIDPTR,&
                &ML, MM, NLM, NLEG, NUMPHASE,&
                &NPTS, DELTAM, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,&
                &SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
                &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
                &EXTINCT8, SRCEXT8)
        else
            call shd_src__computeSrcOneCell(ICELL, GRIDPTR,&
                &NSTOKES, NSTLEG, ML, MM, NLM, NLEG, NUMPHASE,&
                &NPTS, DELTAM, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,&
                &SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
                &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
                &EXTINCT8, SRCEXT8)
        end if

        ! Interpolate the source and extinction to the current point
        IPT1 = GRIDPTR(1, ICELL)
        IPT2 = GRIDPTR(8, ICELL)
        DELX = GRIDPOS(1, IPT2) - GRIDPOS(1, IPT1)
        if (DELX .le. 0.0) then
            INVDELX = 1.0
        else
            INVDELX = 1.0d0/DELX
        end if
        DELY = GRIDPOS(2, IPT2) - GRIDPOS(2, IPT1)
        if (DELY .le. 0.0) then
            INVDELY = 1.0
        else
            INVDELY = 1.0d0/DELY
        end if
        DELZ = GRIDPOS(3, IPT2) - GRIDPOS(3, IPT1)
        INVDELZ = 1.0d0/DELZ
        U = (XE - GRIDPOS(1, IPT1))*INVDELX
        V = (YE - GRIDPOS(2, IPT1))*INVDELY
        W = (ZE - GRIDPOS(3, IPT1))*INVDELZ
        SRCEXT1(:) = (1 - W)*((1 - V)*((1 - U)*SRCEXT8(:, 1) + U*SRCEXT8(:, 2))&
            &+ V*((1 - U)*SRCEXT8(:, 3) + U*SRCEXT8(:, 4)))&
            &+ W*((1 - V)*((1 - U)*SRCEXT8(:, 5) + U*SRCEXT8(:, 6))&
            &+ V*((1 - U)*SRCEXT8(:, 7) + U*SRCEXT8(:, 8)))
        SRCEXT1(1) = max(0.0, SRCEXT1(1))
        EXT1 = (1 - W)*((1 - V)*((1 - U)*EXTINCT8(1) + U*EXTINCT8(2))&
            &+ V*((1 - U)*EXTINCT8(3) + U*EXTINCT8(4)))&
            &+ W*((1 - V)*((1 - U)*EXTINCT8(5) + U*EXTINCT8(6))&
            &+ V*((1 - U)*EXTINCT8(7) + U*EXTINCT8(8)))

        ! This cell is independent pixel if IP mode or open boundary
        !   conditions and ray is leaving domain (i.e. not entering)
        IPINX = btest(int(CELLFLAGS(ICELL)), 0) .and.&
            &.not. (btest(BCFLAG, 0) .and.&
            &((CX .gt. 0 .and. XE .lt. XM) .or. (CX .lt. 0 .and. XE .gt. XM)))
        IPINY = btest(int(CELLFLAGS(ICELL)), 1) .and.&
            &.not. (btest(BCFLAG, 1) .and.&
            &((CY .gt. 0 .and. YE .lt. YM) .or. (CY .lt. 0 .and. YE .gt. YM)))

        ! Find boundaries of the current cell
        ! Find the three possible intersection planes (X,Y,Z)
        !   from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9 - IOCT, ICELL)
        ! Get the distances to the 3 planes and select the closest
        !   (always need to deal with the cell that is wrapped)
        if (IPINX) then
            SOX = 1.0e20
        else
            SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
        end if
        if (IPINY) then
            SOY = 1.0e20
        else
            SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
        end if
        SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Error [_integrateSrcOneRay]: SO<0 ',&
                &MU2, PHI2, XE, YE, ZE, SO, ICELL
            stop
        end if
        XN = XE + SO*CX
        YN = YE + SO*CY
        ZN = ZE + SO*CZ

        ! Find the optical path across the grid cell and figure how
        !   many subgrid intervals to use
        U = (XN - GRIDPOS(1, IPT1))*INVDELX
        V = (YN - GRIDPOS(2, IPT1))*INVDELY
        W = (ZN - GRIDPOS(3, IPT1))*INVDELZ
        EXTN = (1 - W)*((1 - V)*((1 - U)*EXTINCT8(1) + U*EXTINCT8(2))&
            &+ V*((1 - U)*EXTINCT8(3) + U*EXTINCT8(4)))&
            &+ W*((1 - V)*((1 - U)*EXTINCT8(5) + U*EXTINCT8(6))&
            &+ V*((1 - U)*EXTINCT8(7) + U*EXTINCT8(8)))
        TAUGRID = SO*0.5*(EXT1 + EXTN)
        NTAU = max(1, 1 + int(TAUGRID/TAUTOL))
        DELS = SO/NTAU

        ! Loop over the subgrid cells
        do IT = 1, NTAU
            S = IT*DELS
            XI = XE + S*CX
            YI = YE + S*CY
            ZI = ZE + S*CZ
            ! Interpolate extinction and source function along path
            U = (XI - GRIDPOS(1, IPT1))*INVDELX
            V = (YI - GRIDPOS(2, IPT1))*INVDELY
            W = (ZI - GRIDPOS(3, IPT1))*INVDELZ
            EXT0 = (1 - W)*((1 - V)*((1 - U)*EXTINCT8(1) + U*EXTINCT8(2))&
                &+ V*((1 - U)*EXTINCT8(3) + U*EXTINCT8(4)))&
                &+ W*((1 - V)*((1 - U)*EXTINCT8(5) + U*EXTINCT8(6))&
                &+ V*((1 - U)*EXTINCT8(7) + U*EXTINCT8(8)))
            SRCEXT0(:) = (1 - W)*((1 - V)*((1 - U)*SRCEXT8(:, 1) + U*SRCEXT8(:, 2))&
                &+ V*((1 - U)*SRCEXT8(:, 3) + U*SRCEXT8(:, 4)))&
                &+ W*((1 - V)*((1 - U)*SRCEXT8(:, 5) + U*SRCEXT8(:, 6))&
                &+ V*((1 - U)*SRCEXT8(:, 7) + U*SRCEXT8(:, 8)))
            SRCEXT0(1) = max(0.0, SRCEXT0(1))

            ! Compute the subgrid radiance: integration of the source fun
            EXT = 0.5*(EXT0 + EXT1)
            if (EXT .ne. 0.0) then
                TAU = EXT*DELS
                ABSCELL = TAU*(1.0d0 - 0.5d0*TAU*(1.0d0 - FRAC13_*TAU))
                TRANSCELL = 1.0d0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5d0*(SRCEXT0(:) + SRCEXT1(:))&
                    &+ FRAC112_*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*DELS&
                    &*(1.0d0 - 0.05d0*(EXT1 - EXT0)*DELS))/EXT
            else
                ABSCELL = 0.0
                TRANSCELL = 1.0
                SRC(:) = 0.0
            end if

            RADIANCE(:) = RADIANCE(:) + TRANSMIT*SRC(:)*ABSCELL
            TRANSMIT = TRANSMIT*TRANSCELL
            EXT1 = EXT0
            SRCEXT1(:) = SRCEXT0(:)
        ! End of sub grid cell loop
        end do

        ! Get the intersection face number (i.e. neighptr index)
        if (SOX .le. SOZ .and. SOX .le. SOY) then
            IFACE = 2 - BITX
            JFACE = 1
            OPENBCFACE = btest(int(CELLFLAGS(ICELL)), 0) .and. btest(BCFLAG, 0)
        else if (SOY .le. SOZ) then
            IFACE = 4 - BITY
            JFACE = 2
            OPENBCFACE = btest(int(CELLFLAGS(ICELL)), 1) .and. btest(BCFLAG, 1)
        else
            IFACE = 6 - BITZ
            JFACE = 3
            OPENBCFACE = .false.
        end if
        ! Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE, ICELL)
        if (INEXTCELL .lt. 0) then
            call shd_util__findNextCell(XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
        end if
        ! If going to same or larger face then use previous face
        if (NEIGHPTR(IFACE, ICELL) .ge. 0 .and. .not. OPENBCFACE) then
            KFACE = IFACE
            IC = ICELL
        else
            ! If going to smaller face then use next face (more accurate
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
            IFACE = 0
        end if
        ! Get the location coordinate
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XN = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YN = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
            else
                ZN = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
            end if
        end if

        ! If the transmission is greater than zero and not at a
        !   boundary then prepare for next cell
        if (TRANSMIT .lt. TRANSCUT .or. NGRID .gt. MAXCELLSCROSS) then
            VALIDRAD = .true.
        else if (INEXTCELL .eq. 0 .and. IFACE .ge. 5) then
            VALIDRAD = .true.
            call shd_src__findBoundaryRad(NSTOKES, XN, YN,&
                &sngl(MU2), sngl(PHI2),&
                &IC, KFACE, GRIDPTR, GRIDPOS,&
                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                &RADBND)
            RADIANCE(:) = RADIANCE(:) + TRANSMIT*RADBND(:)
        else
            ICELL = INEXTCELL
        end if
        XE = XN
        YE = YN
        ZE = ZN

    end do

    SIDE = IFACE

    if (SRCTYPE .ne. 'T' .and. DELTAM) then
        if (NUMPHASE .gt. 0) then
            deallocate (SINGSCAT)
        else
            deallocate (SUNDIRLEG)
        end if
    end if
    deallocate (YLMDIR)
    return
end subroutine shd_src__integrateSrcOneRay

subroutine shd_src__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
    &NX, NY, NZ, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &MU, PHI, GRIDRAD, EXTINCT, SOURCE,&
    &X0, Y0, Z0, XE, YE, ZE, SIDE,&
    &TRANSMIT, RADIANCE, VALIDRAD)

    ! Integrates the source function (SOURCE) through the extinction
    ! field (EXTINCT) backward in the direction opposite to (MU,PHI)
    ! to find the outgoing radiance (RADIANCE) at the point X0,Y0,Z0.
    ! The radiances at the opposite boundary are input in GRIDRAD.
    ! The transmission and radiance of the ray so far (TRANSMIT, RADIANCE)
    ! are input and returned after the integration along with the exitting
    ! ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
    ! 5=-Z,6=+Z).

    integer BCFLAG, IPFLAG, NSTOKES, NX, NY, NZ, NPTS, NCELLS, SIDE
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    logical VALIDRAD
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real MU, PHI
    real(RD_) :: X0, Y0, Z0, XE, YE, ZE, TRANSMIT, RADIANCE(NSTOKES)
    real EXTINCT(NPTS)
    real SOURCE(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
    integer I1, I2, I3, I4, IOPP, NTAU, IT, K
    logical IPINX, IPINY, OPENBCFACE
    integer GRIDFACE(4, 6), OPPFACE(6), JFACE, KFACE, IC
    real EXT0(1), EXT1(1), EXTN(1)
    real XM, YM
    real SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RAD0(NSTOKES)
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XN, YN, ZN, XI, YI, ZI
    real(RD_) :: SO, SOX, SOY, SOZ, EPS
    real(RD_) :: TAUTOL, TAUGRID, S, DELS
    real(RD_) :: EXT, TAU, TRANSCELL, ABSCELL
    real(RD_) :: SRC(NSTOKES), RAD(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/

    ! TAUTOL is the maximum optical path for the subgrid intervals
    TAUTOL = 0.1
    EPS = 1.0e-4*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    CX = sqrt(1.0 - MU**2)*cos(PHI + DPI_)
    CY = sqrt(1.0 - MU**2)*sin(PHI + DPI_)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-5) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup for the ray path direction
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [_integrateSrc]: Bad MU.'
    end if
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    XM = 0.5*(XGRID(1) + XGRID(NX))
    YM = 0.5*(YGRID(1) + YGRID(NY))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    call shd_util__findCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
        &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
        &BCFLAG, IPFLAG, XE, YE, ZE, ICELL)
    ! print '(A,5(1X,F8.5),1X,I6)', &
    !     &'_integrateSrc:', X0, Y0, Z0, MU, PHI, ICELL
    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
        &XE, YE, ZE, 1, EXTINCT, 1, 1, EXT1)
    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
        &XE, YE, ZE, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT1)
    SRCEXT1(1) = max(0.0, SRCEXT1(1))

    ! Loop until finding a face with known radiances
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
        ! This cell is independent pixel if IP mode or open boundary
        !   conditions and ray is leaving domain (i.e. not entering)
        IPINX = btest(int(CELLFLAGS(ICELL)), 0) .and.&
            &.not. (btest(BCFLAG, 0) .and.&
            &((CX .gt. 0 .and. XE .lt. XM) .or. (CX .lt. 0 .and. XE .gt. XM)))
        IPINY = btest(int(CELLFLAGS(ICELL)), 1) .and.&
            &.not. (btest(BCFLAG, 1) .and.&
            &((CY .gt. 0 .and. YE .lt. YM) .or. (CY .lt. 0 .and. YE .gt. YM)))

        ! Find boundaries of the current cell
        ! Find the three possible intersection planes (X,Y,Z)
        ! from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9 - IOCT, ICELL)
        ! Get the distances to the 3 planes and select the closest
        !   (always need to deal with the cell that is wrapped)
        if (IPINX .or. CX .eq. 0.0) then
            SOX = 1.0e20
        else
            SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
        end if
        if (IPINY .or. CY .eq. 0.0) then
            SOY = 1.0e20
        else
            SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
        end if
        SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Error [_integrateSrc]: SO<0 ',&
                &MU, PHI, X0, Y0, Z0, XE, YE, ZE, SO, SOX, SOY, SOZ, XGRID(1), YGRID(1),&
                &GRIDPOS(1, IOPP), GRIDPOS(2, IOPP), GRIDPOS(3, IOPP)
            stop
        end if
        if (IPINX) then
            XN = XE
        else
            XN = XE + SO*CX
        end if
        if (IPINY) then
            YN = YE
        else
            YN = YE + SO*CY
        end if
        ZN = ZE + SO*CZ

        ! Find the optical path across the grid cell and figure how
        !   many subgrid intervals to use
        call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
            &XN, YN, ZN, 1, EXTINCT, 1, 1, EXTN)
        TAUGRID = SO*0.5*(EXT1(1) + EXTN(1))
        NTAU = max(1, 1 + int(TAUGRID/TAUTOL))
        DELS = SO/NTAU
        ! Loop over the subgrid cells
        do IT = 1, NTAU
            S = IT*DELS
            if (IPINX) then
                XI = XE
            else
                XI = XE + S*CX
            end if
            if (IPINY) then
                YI = YE
            else
                YI = YE + S*CY
            end if
            ZI = ZE + S*CZ
            ! Interpolate extinction and source function along path
            if (IT .eq. NTAU) then
                EXT0(1) = EXTN(1)
            else
                call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
                    &XI, YI, ZI, 1, EXTINCT, 1, 1, EXT0)
            end if
            call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
                &XI, YI, ZI, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT0)
            SRCEXT0(1) = max(0.0, SRCEXT0(1))

            ! Compute the subgrid radiance: integration of the source fun
            EXT = 0.5*(EXT0(1) + EXT1(1))
            if (EXT .ne. 0.0) then
                TAU = EXT*DELS
                ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - FRAC13_*TAU))
                TRANSCELL = 1.0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + FRAC112_*&
                    &(EXT0(1)*SRCEXT1(:) - EXT1(1)*SRCEXT0(:))*DELS&
                    &*(1.0 - 0.05*(EXT1(1) - EXT0(1))*DELS))/EXT
            else
                ABSCELL = 0.0
                TRANSCELL = 1.0
                SRC(:) = 0.0
            end if

            RADIANCE(:) = RADIANCE(:) + TRANSMIT*SRC(:)*ABSCELL
            SRCEXT1(:) = SRCEXT0(:)
            EXT1(1) = EXT0(1)
            TRANSMIT = TRANSMIT*TRANSCELL
        ! End of sub grid cell loop
        end do

        ! Get the intersection face number (i.e. neighptr index)
        if (SOX .le. SOZ .and. SOX .le. SOY) then
            IFACE = 2 - BITX
            JFACE = 1
            OPENBCFACE = btest(int(CELLFLAGS(ICELL)), 0) .and. btest(BCFLAG, 0)
        else if (SOY .le. SOZ) then
            IFACE = 4 - BITY
            JFACE = 2
            OPENBCFACE = btest(int(CELLFLAGS(ICELL)), 1) .and. btest(BCFLAG, 1)
        else
            IFACE = 6 - BITZ
            JFACE = 3
            OPENBCFACE = .false.
        end if
        ! Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE, ICELL)
        if (INEXTCELL .lt. 0) then
            call shd_util__findNextCell(XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
        end if
        ! Get the four grid points on the intersection face for
        !   telling if at the boundary
        !   If going to same or larger face then use previous face
        if (NEIGHPTR(IFACE, ICELL) .ge. 0 .and. .not. OPENBCFACE) then
            KFACE = IFACE
            IC = ICELL
        else
            ! If going to smaller face then use next face (more accurate
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
        end if
        I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
        I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
        I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
        I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
        ! Get the location coordinate
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XN = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YN = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
            else
                ZN = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
            end if
        end if

        ! See if there are valid radiances for the intersection face:
        !   If so, set flag, interpolate radiance, and add in
        !   If not, only prepare for next cell
        if (TRANSMIT .lt. 1.0e-5) then
            VALIDRAD = .true.
        else if (GRIDRAD(1, I1) .ge. 0.0 .and. GRIDRAD(1, I2) .ge. 0.0 .and.&
            &GRIDRAD(1, I3) .ge. 0.0 .and. GRIDRAD(1, I4) .ge. 0.0) then
            VALIDRAD = .true.
            call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS, XI, YI, ZI,&
                &NSTOKES, GRIDRAD, 1, NSTOKES, RAD0)
            RADIANCE(:) = RADIANCE(:) + TRANSMIT*RAD0(:)
        else
            ICELL = INEXTCELL
        end if
        XE = XN
        YE = YN
        ZE = ZN

    end do
    SIDE = IFACE
    return
end subroutine shd_src__integrateSrc

subroutine shd_src__computeSolarStokesVect(&
    &NLEG, SINGLEGEN, COSSCAT,&
    &SOLARMU, MU, DELPHI, NSTOKES,&
    &SCATVECT)

    ! Calculates the Stokes vector for a single scattering of unpolarized
    ! sunlight from the solar direction (SOLARMU) to the outgoing
    ! direction (MU) and with cosine scattering angle (COSSCAT).
    ! The A1 and B1 elements of the phase matrix are calculated from
    ! the Wigner function expansion coefficients (SINGLEGEN).
    ! The plane of polarization is rotated from the scattering plane to
    ! the outgoing meridional plane.  Also divides the output vector by 4*pi.

    integer NLEG, NSTOKES
    real SINGLEGEN(2, 0:NLEG), COSSCAT, SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    integer I, J, K, L
    real(RD_) :: X, FCT, A1, B1
    real(RD_), allocatable :: DMM1(:)
    real(RD_) :: SIN_SCAT, SIN_THETA1, SIN_THETA2
    real(RD_) :: SINPHI, COSPHI, SIN2, COS2, SIN22, COS22

    allocate (DMM1(0:NLEG))

    ! Sum the first Wigner function series, which is actually a Legend
    X = dble(COSSCAT)
    call shd_pha__computeWigner(X, NLEG, 0, 0, DMM1)
    A1 = 0.0d0
    do L = 0, NLEG
        FCT = 2.0d0*L + 1.0d0
        A1 = A1 + FCT*dble(SINGLEGEN(1, L))*DMM1(L)
    end do
    SCATVECT(1) = A1*DO4PI_

    if (NSTOKES .gt. 1) then
        ! If doing polarization, sum the second Wigner function series
        call shd_pha__computeWigner(X, NLEG, 2, 0, DMM1)
        B1 = 0.0d0
        do L = 0, NLEG
            FCT = 2.0d0*L + 1.0d0
            B1 = B1 - FCT*dble(SINGLEGEN(2, L))*DMM1(L)
        end do

        ! Rotate the plane of polarization from the scattering plane
        !   to the outgoing meridional plane.
        SIN_SCAT = dsqrt(max(0.0d0, 1.d0 - COSSCAT**2))
        SIN_THETA1 = dsqrt(1.d0 - SOLARMU**2)
        SIN_THETA2 = dsqrt(1.d0 - MU**2)
        SINPHI = dsin(dble(DELPHI))
        COSPHI = dcos(dble(DELPHI))
        if (SIN_SCAT .eq. 0.0) then
            SIN2 = 0.0d0
            COS2 = -1.0d0
        else
            SIN2 = SIN_THETA1*SINPHI/SIN_SCAT
            COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        end if
        SIN22 = 2.0d0*SIN2*COS2
        COS22 = 1.0d0 - 2.0d0*SIN2**2
        SCATVECT(2) = B1*COS22*DO4PI_
        SCATVECT(3) = B1*SIN22*DO4PI_
    end if
    if (NSTOKES .eq. 4) then
        SCATVECT(4) = 0.0
    end if
    deallocate (DMM1)

    return

end subroutine shd_src__computeSolarStokesVect

subroutine shd_src__findBoundaryRad(&
    &NSTOKES, XB, YB, MU2, PHI2,&
    &ICELL, KFACE, GRIDPTR, GRIDPOS,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &RADBND)

    !   Returns the interpolated Stokes radiance at the boundary (RADBND)
    ! Inputs are the boundary location (XB,YB), ray direction away from
    ! boundary (MU2,PHI2), cell number (ICELL) and face (KFACE) at
    ! the boundary point.

    integer NSTOKES, ICELL, KFACE, MAXNBC, NTOPPTS, NBOTPTS
    integer GRIDPTR(8, *), BCPTR(MAXNBC, 2)
    integer NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
    real MU2, PHI2, RADBND(NSTOKES)
    real(RD_) :: XB, YB
    real GRIDPOS(3, *)
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS), BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE
    character(len=2) :: SFCTYPE

    integer IL, IM, IU, IP, IBC, J
    logical LAMBERTIAN
    real X(4), Y(4), RAD(NSTOKES, 4), U, V
    integer GRIDFACE(4, 6)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'

    ! Loop over the four gridpoints of the cell face on the boundary
    do J = 1, 4
        IP = GRIDPTR(GRIDFACE(J, KFACE), ICELL)
        X(J) = GRIDPOS(1, IP)
        Y(J) = GRIDPOS(2, IP)

        if (MU2 .lt. 0.0) then
            ! Do a binary search to locate the top boundary point
            IL = 1
            IU = NTOPPTS
            do while (IU - IL .gt. 1)
                IM = (IU + IL)/2
                if (IP .ge. BCPTR(IM, 1)) then
                    IL = IM
                else
                    IU = IM
                end if
            end do
            IBC = IL
            if (BCPTR(IBC, 1) .ne. IP) IBC = IU
            if (BCPTR(IBC, 1) .ne. IP)&
                &stop 'Error [_findBoundaryRad]: Not at boundary.'
            RAD(:, J) = BCRAD(:, IBC)

        else

            ! Do a binary search to locate the bottom boundary point
            IL = 1
            IU = NBOTPTS
            do while (IU - IL .gt. 1)
                IM = (IU + IL)/2
                if (IP .ge. BCPTR(IM, 2)) then
                    IL = IM
                else
                    IU = IM
                end if
            end do
            IBC = IL
            if (BCPTR(IBC, 2) .ne. IP) IBC = IU
            if (BCPTR(IBC, 2) .ne. IP)&
                &stop 'Error [_findBoundaryRad]: Not at boundary.'

            if (.not. LAMBERTIAN) then
                call shd_init__makeSfcRadV(NBOTPTS, IBC, IBC, BCPTR(1, 2),&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,&
                    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                    &SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES,&
                    &BCRAD(:, 1 + NTOPPTS))
            end if
            RAD(:, J) = BCRAD(:, NTOPPTS + IBC)
        end if
    end do
    if (X(2) - X(1) .gt. 0.0) then
        U = (XB - X(1))/(X(2) - X(1))
    else
        U = 0.0
    end if
    if (Y(3) - Y(1) .gt. 0.0) then
        V = (YB - Y(1))/(Y(3) - Y(1))
    else
        V = 0.0
    end if
    RADBND(:) = (1 - U)*(1 - V)*RAD(:, 1) + U*(1 - V)*RAD(:, 2)&
        &+ (1 - U)*V*RAD(:, 3) + U*V*RAD(:, 4)
    return
end subroutine shd_src__findBoundaryRad

subroutine shd_src__rotatePolarPlane(NSTOKES, COSSCAT,&
    &SOLARMU, MU, DELPHI, SCATVECT)

    ! Rotates the plane of polarization from the scattering plane
    ! to the outgoing meridional plane for the solar single
    ! scattered Stokes vector in SINGSCAT.

    integer NSTOKES
    real(RD_) :: COSSCAT
    real SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    real B1
    real(RD_) :: SIN_SCAT, SIN_THETA1, SIN_THETA2
    real(RD_) :: SINPHI, COSPHI, SIN2, COS2, SIN22, COS22

    if (NSTOKES .gt. 1) then
        B1 = SCATVECT(2)
        SIN_SCAT = dsqrt(max(0.0d0, 1.d0 - COSSCAT**2))
        SIN_THETA1 = dsqrt(1.d0 - SOLARMU**2)
        SIN_THETA2 = dsqrt(1.d0 - MU**2)
        SINPHI = dsin(dble(DELPHI))
        COSPHI = dcos(dble(DELPHI))
        if (SIN_SCAT .eq. 0.0) then
            SIN2 = 0.0d0
            COS2 = -1.0d0
        else
            SIN2 = SIN_THETA1*SINPHI/SIN_SCAT
            COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        end if
        SIN22 = 2.0d0*SIN2*COS2
        COS22 = 1.0d0 - 2.0d0*SIN2**2
        SCATVECT(2) = B1*COS22
        SCATVECT(3) = B1*SIN22
    end if
    if (NSTOKES .eq. 4) then
        SCATVECT(4) = 0.0
    end if

    return

end subroutine shd_src__rotatePolarPlane

subroutine shd_src__computeSrcOnePnt(&
    &NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    &SOURCET)

    ! Calculates the polarized source function in real generalized
    ! spherical harmonics at one grid point from the radiance in
    ! real generalized spherical harmonics and thermal and/or solar
    ! source information.

    integer, intent(in) :: NSTOKES, NSTLEG, NLEG, NLM, LOFJ(NLM), IPH, NR
    real, intent(in) :: FLUX0, YLMSUN(NSTLEG, NLM), PLANCK, ALBEDO
    real, intent(in) :: LEGEN(NSTLEG, 0:NLEG, *)
    real, intent(in) :: RADIANCE(NSTOKES, NR)
    character(len=1), intent(in) :: SRCTYPE

    real, intent(out) :: SOURCET(NSTOKES, NLM)

    real, parameter :: C = 3.544907703
    integer J

    SOURCET(:, :) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.  Do different parts of the
    ! matrix multiply depending on NSTOKES.
    do J = 1, NR
        SOURCET(1, J) = SOURCET(1, J) + LEGEN(1, LOFJ(J), IPH)*RADIANCE(1, J)
    end do

    if (NSTOKES .gt. 1) then

        do J = 1, NR
            SOURCET(1, J) = SOURCET(1, J) + LEGEN(5, LOFJ(J), IPH)*RADIANCE(2, J)
        end do

        do J = 5, NR
            SOURCET(2, J) = SOURCET(2, J) + LEGEN(5, LOFJ(J), IPH)*RADIANCE(1, J)&
                &+ LEGEN(2, LOFJ(J), IPH)*RADIANCE(2, J)
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(3, LOFJ(J), IPH)*RADIANCE(3, J)
        end do

    end if

    if (NSTOKES .eq. 4) then

        do J = 1, NR
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(6, LOFJ(J), IPH)*RADIANCE(4, J)
            SOURCET(4, J) = SOURCET(4, J) - LEGEN(6, LOFJ(J), IPH)*RADIANCE(3, J)&
                &+ LEGEN(4, LOFJ(J), IPH)*RADIANCE(4, J)
        end do

    end if

    ! Do the final multiplication by the single scattering albedo
    SOURCET(:, :) = ALBEDO*SOURCET(:, :)

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        do J = 1, NLM
            SOURCET(1, J) = SOURCET(1, J)&
                &+ FLUX0*ALBEDO*LEGEN(1, LOFJ(J), IPH)*YLMSUN(1, J)
        end do

        if (NSTOKES .gt. 1) then
            do J = 5, NLM
                SOURCET(2, J) = SOURCET(2, J)&
                    &+ FLUX0*ALBEDO*LEGEN(5, LOFJ(J), IPH)*YLMSUN(1, J)
            end do
        end if

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1, 1) = SOURCET(1, 1) + C*PLANCK
    end if

    return

end subroutine shd_src__computeSrcOnePnt

subroutine shd_src__computeSrcOnePnt_(&
    &NLM, NLEG, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    &SOURCET)

    ! Calculates the scalar source function in spherical harmonics at
    ! one grid point from the radiance in spherical harmonics and
    ! thermal and/or solar source information.

    integer, intent(in) :: NLEG, NLM, LOFJ(NLM), IPH, NR
    character(len=1), intent(in) :: SRCTYPE
    real, intent(in) :: FLUX0, YLMSUN(NLM), PLANCK, ALBEDO
    real, intent(in) :: LEGEN(0:NLEG, *)
    real, intent(in) :: RADIANCE(NR)

    real, intent(out) :: SOURCET(NLM)

    real, parameter :: C = 3.544907703
    integer J

    SOURCET(:) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        do J = 1, NLM
            SOURCET(J) = FLUX0*ALBEDO*LEGEN(LOFJ(J), IPH)*YLMSUN(J)
        end do

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1) = SOURCET(1) + C*PLANCK
    end if

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.
    do J = 1, NR
        SOURCET(J) = SOURCET(J) + ALBEDO*LEGEN(LOFJ(J), IPH)*RADIANCE(J)
    end do

    return

end subroutine shd_src__computeSrcOnePnt_

end module shd_src
