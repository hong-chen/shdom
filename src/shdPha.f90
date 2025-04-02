module shd_pha

    use globals

    implicit none

    public :: shd_pha__computePhase
    public :: shd_pha__computeLegendre
    public :: shd_pha__computeYLM
    public :: shd_pha__computePLM
    public :: shd_pha__computeWigner

    private :: shd_pha__computeYLM_
    private :: shd_pha__computeWignerDM0
    private :: shd_pha__computeWignerDM0DM2PDM2M

    private :: wignerDMM1N0
    private :: wignerDMM0N0

contains

subroutine shd_pha__computePhase(NSCATANGLE, NUMPHASE, NSTPHASE,&
    &NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)

    !   Precomputes the phase matrix elements I-I and I-Q as a function
    ! of scattering angle for solar direct scattering for all the
    ! tabulated phase functions. Output is in PHASETAB.
    ! NSTPHASE=1 for NSTOKES=1 and NSTPHASE=2 for NSTOKES>1.

    integer NSCATANGLE, NUMPHASE, NSTPHASE
    integer NSTOKES, NSTLEG, ML, NLEG
    real LEGEN(NSTLEG, 0:NLEG, NUMPHASE)
    real PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE)
    integer IPH, I, J, L
    double precision PI, OFOURPI, COSSCAT, FCT, F, X, A1, B1
    double precision, allocatable :: UNSCLEGEN(:, :)
    double precision, allocatable :: DMM1(:), DMM2(:)

    allocate (UNSCLEGEN(2, 0:NLEG), DMM1(0:NLEG), DMM2(0:NLEG))
    PI = acos(-1.0d0)
    OFOURPI = 1.0/(4.0*PI)

    do J = 1, NSCATANGLE
        COSSCAT = cos(PI*DFLOAT(J - 1)/(NSCATANGLE - 1))
        X = dble(COSSCAT)
        call shd_pha__computeWigner(X, NLEG, 0, 0, DMM1)
        if (NSTOKES .gt. 1) then
            call shd_pha__computeWigner(X, NLEG, 2, 0, DMM2)
        end if

        ! Unscaled the needed Legendre coefficients, but divide by 1-f
        !   because extinction is still scaled (for TMS method)
        do IPH = 1, NUMPHASE
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

            ! Sum the first Wigner function series, which is actually a Le
            A1 = 0.0d0
            do L = 0, NLEG
                FCT = 2.0d0*L + 1.0d0
                A1 = A1 + FCT*dble(UNSCLEGEN(1, L))*DMM1(L)
            end do
            if (A1 .le. 0.0) then
                write (6, *) 'Error [_computePhase]: negative phase function',&
                    &' for tabulated phase function: ', IPH, J, A1
                stop
            end if
            PHASETAB(1, IPH, J) = SNGL(A1*OFOURPI)

            if (NSTOKES .gt. 1) then
                ! If doing polarization, sum the second Wigner function series
                B1 = 0.0d0
                do L = 0, NLEG
                    FCT = 2.0d0*L + 1.0d0
                    B1 = B1 - FCT*dble(UNSCLEGEN(2, L))*DMM2(L)
                end do
                PHASETAB(2, IPH, J) = SNGL(B1*OFOURPI)
            end if
        end do
    end do

    deallocate (UNSCLEGEN, DMM1, DMM2)
    return

end subroutine shd_pha__computePhase

subroutine shd_pha__computeLegendre(COSSCAT, NLEG, P)

    !   This subroutine computes a set of Legendre polynomials for
    ! a particular scattering angle COSSCAT.  NLEG is the maximum term.
    ! The Legendre functions evaluated at COSSCAT are returned in
    ! P, starting at l=0 and ending with l=NLEG  (NLEG+1 terms).

    integer, intent(in) :: NLEG
    double precision COSSCAT, P(0:NLEG)
    integer L
    double precision X, PL, PL1, PL2

    X = dble(COSSCAT)
    if (X*X .gt. 1.) stop 'Error [_computeLegendre]: <COSSCAT> is greater than 1.'
    ! Use the stable upward recursion on l, starting from P_0
    PL2 = 1.0d0
    P(0) = PL2
    if (NLEG .gt. 1) then
        PL1 = X
        P(1) = X
    end if
    do L = 2, NLEG
        PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        P(L) = PL
        PL2 = PL1
        PL1 = PL
    end do
    return
end subroutine shd_pha__computeLegendre

subroutine shd_pha__computeYLM(TRANSPOSED, MU, PHI, ML, MM, NSTLEG, YR)

    !   This subroutine computes a set of normalized real generalized
    ! spherical harmonic functions (YR) for a particular direction mu,phi.
    ! (Doicu et al. 2013; http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mode.
    ! NSTLEG is the number of matrix elements returned, either 1 for the
    ! 1,1 term or 6 for the 6 unique elements of the 4x4 matrix.
    ! J = 2*(L*(L+1))/2 + M+1  for L<=MM
    ! J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM

    integer ML, MM, NSTLEG
    real MU, PHI
    logical TRANSPOSED
    real YR(NSTLEG, *)
    integer J, M, MABS, L
    double precision X, PI, FCT, P1, P2, P3, COSM, SINM, SIGN
    double precision, allocatable :: DM0(:), DM2P(:), DM2M(:)

    if (NSTLEG .eq. 1) then
        call shd_pha__computeYLM_(MU, PHI, ML, MM, YR)
        return
    end if

    allocate (DM0(0:ML), DM2P(0:ML), DM2M(0:ML))
    X = dble(MU)

    PI = DACOS(-1.d0)
    FCT = 1.d0/DSQRT(2.d0*PI)

    ! sign accounting for transpose matrices
    if (.not. TRANSPOSED) then
        SIGN = 1.d0
    else
        SIGN = -1.d0
    end if

    ! ..................................................................
    !                               M = 0
    M = 0
    call shd_pha__computeWignerDM0DM2PDM2M(X, ML, M, DM0, DM2P, DM2M)
    do L = 0, ML
        if (L .le. MM) then
            J = L*(L + 1) + M + 1
        else
            J = (2*MM + 1)*L - MM**2 + M + 1
        end if
        P1 = FCT*DM0(L)
        P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
        P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

    ! --- Unm = Pnm*cos(m*phi) ---
        YR(1, J) = P1
        if (NSTLEG .eq. 6) then
            YR(2, J) = P2
            YR(3, J) = P2
            YR(4, J) = P1
            YR(5, J) = P3
            YR(6, J) = P3
        end if
    end do

    ! ..................................................................
    !                               M =/ 0
    do MABS = 1, MM
        call shd_pha__computeWignerDM0DM2PDM2M(X, ML, MABS, DM0, DM2P, DM2M)
        COSM = cos(MABS*PHI)
        SINM = sin(MABS*PHI)
        do L = MABS, ML
    ! ..............................................................
    !                               M > 0
            M = MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if
            P1 = FCT*DM0(L)
            P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
            P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

            ! --- Unm = Pnm*cos(m*phi) - D*Pnm*sin(m*phi) ---
            YR(1, J) = P1*COSM - P1*SINM
            if (NSTLEG .eq. 6) then
                YR(2, J) = P2*COSM - P2*SINM
                YR(3, J) = P2*COSM + P2*SINM
                YR(4, J) = P1*COSM + P1*SINM
                YR(5, J) = P3*COSM - SIGN*P3*SINM
                YR(6, J) = P3*COSM + SIGN*P3*SINM
            end if
    ! ..............................................................
    !                               M < 0
            M = -MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if

            ! --- Vnm = Pnm*sin(m*phi) + D*Pnm*cos(m*phi) ---
            YR(1, J) = P1*SINM + P1*COSM
            if (NSTLEG .eq. 6) then
                YR(2, J) = P2*SINM + P2*COSM
                YR(3, J) = P2*SINM - P2*COSM
                YR(4, J) = P1*SINM - P1*COSM
                YR(5, J) = P3*SINM + SIGN*P3*COSM
                YR(6, J) = P3*SINM - SIGN*P3*COSM
            end if
        end do
    end do
    deallocate (DM0, DM2P, DM2M)
    return
end subroutine shd_pha__computeYLM

subroutine shd_pha__computeYLM_(MU, PHI, ML, MM, YR)

    !   This subroutine computes a set of the unpolarized (I-I) element
    ! normalized real generalized spherical harmonic functions (YR) for
    ! a particular direction mu,phi.  (Doicu et al. 2013;
    ! http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mod
    !   J = 2*(L*(L+1))/2 + M+1  for L<=MM
    !   J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM

    integer ML, MM
    real MU, PHI
    real YR(*)
    integer J, M, L
    double precision X, PI, FCT, COSM, SINM
    double precision, allocatable :: DM0(:)

    allocate (DM0(0:ML))
    X = dble(MU)

    PI = DACOS(-1.d0)
    FCT = 1.d0/DSQRT(2.d0*PI)

    do M = 0, MM
        call shd_pha__computeWignerDM0(X, ML, M, DM0)
        DM0(:) = FCT*DM0(:)
        if (M .gt. 0) then
            COSM = cos(M*PHI)
            SINM = sin(M*PHI)
        else
            COSM = 1.0d0
            SINM = 0.0d0
        end if
        do L = M, ML
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if

            ! M>=0:  Unm = Pnm*cos(m*phi) - D*Pnm*sin(m*phi)
            YR(J) = (COSM - SINM)*DM0(L)

            J = J - 2*M
            ! M<0:  Vnm = Pnm*sin(m*phi) + D*Pnm*cos(m*phi)
            YR(J) = (COSM + SINM)*DM0(L)
        end do
    end do
    deallocate (DM0)
    return
end subroutine shd_pha__computeYLM_

subroutine shd_pha__computePLM(TRANSPOSED, MU, ML, MM, PRC)

    integer ML, MM
    real MU
    logical TRANSPOSED
    real PRC(6, *)
    ! --- local variables ---
    integer J, M, MABS, L
    double precision X, PI, FCT, P1, P2, P3, SIGN
    double precision, allocatable :: DM0(:), DM2P(:),&
        &DM2M(:)

    ! --- allocate ---
    allocate (DM0(0:ML), DM2P(0:ML), DM2M(0:ML))

    X = dble(MU)

    ! --- normalization constant ---
    PI = DACOS(-1.d0)
    FCT = 1.d0/DSQRT(2.d0*PI)

    ! --- sign accounting for transpose matrices ---
    if (.not. TRANSPOSED) then
        SIGN = 1.d0
    else
        SIGN = -1.d0
    end if
    ! ..................................................................
    !                            M = 0
    ! ..................................................................
    M = 0
    call shd_pha__computeWignerDM0DM2PDM2M&
        &(X, ML, M,&
        &DM0, DM2P, DM2M)
    do L = 0, ML
        if (L .le. MM) then
            J = L*(L + 1) + M + 1
        else
            J = (2*MM + 1)*L - MM**2 + M + 1
        end if
        P1 = FCT*DM0(L)
        P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
        P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

        ! --- PRCnm = Pnm ---
        PRC(1, J) = P1
        PRC(2, J) = P2
        PRC(3, J) = P2
        PRC(4, J) = P1
        PRC(5, J) = P3
        PRC(6, J) = P3
    end do
    ! ..................................................................
    !                               M =/ 0
    ! ..................................................................
    do MABS = 1, MM
        call shd_pha__computeWignerDM0DM2PDM2M&
            &(X, ML, MABS,&
            &DM0, DM2P, DM2M)
        do L = MABS, ML
    ! ..................................................................
    !                               M > 0
    ! ..................................................................
            M = MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if
            P1 = FCT*DM0(L)
            P2 = -0.5d0*FCT*(DM2P(L) + DM2M(L))
            P3 = -0.5d0*FCT*(DM2P(L) - DM2M(L))

            ! --- PRCnm = Pnm ---
            PRC(1, J) = P1
            PRC(2, J) = P2
            PRC(3, J) = P2
            PRC(4, J) = P1
            PRC(5, J) = P3
            PRC(6, J) = P3
    ! ..................................................................
    !                               M < 0
    ! ..................................................................
            M = -MABS
            if (L .le. MM) then
                J = L*(L + 1) + M + 1
            else
                J = (2*MM + 1)*L - MM**2 + M + 1
            end if

            ! --- PRCnm = D*Pnm ---
            PRC(1, J) = P1
            PRC(2, J) = P2
            PRC(3, J) = -P2
            PRC(4, J) = -P1
            PRC(5, J) = SIGN*P3
            PRC(6, J) = -SIGN*P3
        end do
    end do
    deallocate (DM0, DM2P, DM2M)

    return

end subroutine shd_pha__computePLM

subroutine shd_pha__computeWigner(X, NRANK, M, M1, DMM1)

    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for
    ! |               M >= 0, M1 >= 0,  N = 0,...,NRANK, and
    ! |               -1 < X= COS(BETA) < 1
    ! ----------------------------------------------------------------

    integer NRANK, M, M1
    double precision X, DMM1(0:NRANK)

    integer N0, N
    double precision FACT1, FACT2, wignerDMM1N0

    N0 = max(M, M1)
    DMM1 = 0.d0
    if (N0 .eq. 0) then
        DMM1(0) = 1.d0
        DMM1(1) = X
        do N = 1, NRANK - 1
            FACT1 = dble(2*N + 1)*X/dble(N + 1)
            FACT2 = dble(N)/dble(N + 1)
            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    else
        DMM1(N0) = wignerDMM1N0(X, M, M1)
        do N = N0, NRANK - 1
            FACT1 = dble(N*(N + 1))*X - dble(M*M1)
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M**2))
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M1**2))
            FACT1 = FACT1*dble(2*N + 1)/dble(N)

            FACT2 = DSQRT(dble(N**2 - M**2))&
                &*DSQRT(dble(N**2 - M1**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M1**2))
            FACT2 = FACT2*dble(N + 1)/dble(N)

            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    end if
    return
end subroutine shd_pha__computeWigner

subroutine shd_pha__computeWignerDM0DM2PDM2M(&
    &MIU, NRANK, M, DM0, DM2P, DM2M)

    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for
    ! |               M >= 0; M1 = 0, 2, -2; N = 0,...,NRANK, and
    ! |              -1 <= MIU <= 1
    ! |  For M1 = -2, we use the symmetry relation
    ! |               DM-M1(MIU) = (-1)**(N+M) * DMM1(-MIU)
    ! ----------------------------------------------------------------

    integer NRANK, M
    double precision MIU
    double precision DM0(0:NRANK), DM2P(0:NRANK), DM2M(0:NRANK)
    integer N0, N
    double precision XP, XM, FACT1, FACT2, FACTP, FACTM, wignerDMM1N0

    !for M1 = 0 and M1 = 2
    XP = MIU
    !for M1 = - 2
    XM = -MIU
    ! ................................................................
    ! DM0 for M1 = 0
    ! ................................................................
    DM0 = 0.d0
    if (M .le. NRANK) then
        if (M .eq. 0) then
            DM0(0) = 1.d0
            DM0(1) = XP
            do N = 1, NRANK - 1
                FACT1 = dble(2*N + 1)*XP/dble(N + 1)
                FACT2 = dble(N)/dble(N + 1)
                DM0(N + 1) = FACT1*DM0(N) - FACT2*DM0(N - 1)
            end do
        else
            DM0(M) = wignerDMM1N0(XP, M, 0)
            do N = M, NRANK - 1
                FACT1 = dble(N*(N + 1))*XP
                FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M**2))
                FACT1 = FACT1/dble(N + 1)
                FACT1 = FACT1*dble(2*N + 1)/dble(N)

                FACT2 = DSQRT(dble(N**2 - M**2))*dble(N)
                FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
                FACT2 = FACT2/dble(N + 1)
                FACT2 = FACT2*dble(N + 1)/dble(N)

                DM0(N + 1) = FACT1*DM0(N) - FACT2*DM0(N - 1)
            end do
        end if
        do N = 0, NRANK
            DM0(N) = DSQRT(N + 0.5d0)*DM0(N)
        end do
    end if
    ! ................................................................
    ! DM2P and DM2M for M1 = 2 and M1 = -2
    ! ................................................................
    N0 = max(M, 2)
    DM2P = 0.d0
    DM2M = 0.d0
    if (M .le. NRANK .and. NRANK .ge. 2) then
        DM2P(N0) = wignerDMM1N0(XP, M, 2)
        DM2M(N0) = wignerDMM1N0(XM, M, 2)
        do N = N0, NRANK - 1
            FACTP = dble(N*(N + 1))*XP - dble(2*M)
            FACTM = dble(N*(N + 1))*XM - dble(2*M)

            FACT1 = 1.d0/DSQRT(dble((N + 1)**2 - M**2))
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - 4))
            FACT1 = FACT1*dble(2*N + 1)/dble(N)

            FACT2 = DSQRT(dble(N**2 - M**2))&
                &*DSQRT(dble(N**2 - 4))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - 4))
            FACT2 = FACT2*dble(N + 1)/dble(N)

            DM2P(N + 1) = FACTP*FACT1*DM2P(N) - FACT2*DM2P(N - 1)
            DM2M(N + 1) = FACTM*FACT1*DM2M(N) - FACT2*DM2M(N - 1)
        end do
        do N = 0, NRANK
            DM2M(N) = (-1.d0)**(N + M)*DM2M(N)

            DM2P(N) = DSQRT(N + 0.5d0)*DM2P(N)
            DM2M(N) = DSQRT(N + 0.5d0)*DM2M(N)
        end do
    end if

    return

end subroutine shd_pha__computeWignerDM0DM2PDM2M

subroutine shd_pha__computeWignerDM0(X, NRANK, M, DM0)

    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for      |
    ! |               M >= 0; M1 = 0; N = 0,...,NRANK, and           |
    ! |              -1 <= X <= 1                                    |
    ! ----------------------------------------------------------------

    integer NRANK, M
    double precision X
    double precision DM0(0:NRANK)
    integer N
    double precision wignerDMM0N0

    ! DM0 for M1 = 0
    DM0(:) = 0.d0
    if (M .le. NRANK) then
        if (M .eq. 0) then
            DM0(0) = 1.0d0
            DM0(1) = X
            do N = 1, NRANK - 1
                DM0(N + 1) = ((2*N + 1)*X*DM0(N) - N*DM0(N - 1))/(N + 1)
            end do
        else
            DM0(M) = wignerDMM0N0(X, M)
            do N = M, NRANK - 1
                DM0(N + 1) = ((2*N + 1)*X*DM0(N)&
                    &- DSQRT(dble(N**2 - M**2))*DM0(N - 1))&
                    &/DSQRT(dble((N + 1)**2 - M**2))
            end do
        end if
        do N = 0, NRANK
            DM0(N) = DSQRT(N + 0.5d0)*DM0(N)
        end do
    end if
    return
end subroutine shd_pha__computeWignerDM0

function wignerDMM1N0(X, M, M1) result(DMM1N0)

    ! ---------------------------------------------------------------
    ! | THE ROUTINE COMPUTES THE WIGNER FUNCTIONS FOR               |
    ! |     M >= 0, M1 >= 0, AND N0 = MAX(M,M1)                     |
    ! ---------------------------------------------------------------

    integer M, M1
    double precision X, DMM1N0

    integer P, MAXM, MINM
    double precision PROD, FACT, CMM1

    if (M .eq. M1) then
        FACT = ((1.d0 + X)/2.d0)**M
        DMM1N0 = FACT
    else
        if (M1 .gt. M) then
            CMM1 = 1.d0
        else
            CMM1 = (-1.d0)**(M - M1)
        end if
        MAXM = max(M, M1)
        MINM = min(M, M1)
        PROD = 1.d0
        do P = 1, MAXM - MINM
            FACT = DSQRT(dble(M + M1 + P)/dble(P))
            PROD = PROD*FACT
        end do
        FACT = DSQRT((1.d0 - X)/2.d0)
        DMM1N0 = CMM1*PROD*FACT**(MAXM - MINM)
        FACT = DSQRT((1.d0 + X)/2.d0)
        DMM1N0 = DMM1N0*FACT**(MAXM + MINM)
    end if

    return

end function wignerDMM1N0

function wignerDMM0N0(X, M) result(DMM0N0)

    ! The routine computes the Wigner functions for M >= 0, M1 = 0, a
    ! N0 = M

    integer M
    double precision X, DMM0N0
    integer P
    double precision CM, PROD

    if (M .eq. 0) then
        DMM0N0 = 1.d0
    else
        CM = (-1.d0)**M
        PROD = 1.d0
        do P = 1, M
            PROD = PROD*DSQRT((M + P)/dble(P))
        end do
        DMM0N0 = CM*PROD*(0.5d0*DSQRT((1.d0 - X)*(1.d0 + X)))**M
    end if

    return

end function wignerDMM0N0

end module shd_pha
