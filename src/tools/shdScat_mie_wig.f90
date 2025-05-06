subroutine MIE_ONE(WAVELENGTH, MINDEX, RADIUS,&
    & MAXRANK, EXTINCTION, SCATTER, NRANK, WIGCOEF)
    !   Computes the Mie scattering properties for a single homogeneous
    ! sphere of radius RADIUS.  The six phase matrix elements times the
    ! scattering coefficient is returned as Wigner d-function series
    ! coefficients.
    implicit none
    integer MAXRANK, NRANK
    real WAVELENGTH, RADIUS
    real EXTINCTION, SCATTER, WIGCOEF(6, 0:*)
    complex MINDEX
    integer MAXN
    parameter(MAXN=10000)
    integer NMIE, NQUAD, LASTNQUAD
    integer I, J, L
    double precision X, PI, QEXT, QSCAT, GEOMAREA
    double precision MU(MAXN), WTS(MAXN)
    double precision P1, P2, P3, P4
    double precision D00(0:MAXN), D20(0:MAXN)
    double precision D22P(0:MAXN), D22M(0:MAXN)
    double precision COEF(6, 0:MAXN), F, A2, A3
    double complex MSPHERE
    double complex A(MAXN), B(MAXN)
    save MU, WTS, LASTNQUAD
    data LASTNQUAD/-1/

    PI = dacos(-1.0d0)
    X = 2.0d0*PI*RADIUS/WAVELENGTH
    GEOMAREA = PI*RADIUS**2
    MSPHERE = MINDEX

    ! Compute the An's and Bn's and the cross sections
    NMIE = 0
    call MIECALC(NMIE, X, MSPHERE, A, B)
    call MIECROSS(NMIE, X, A, B, QEXT, QSCAT)
    EXTINCTION = GEOMAREA*QEXT
    SCATTER = GEOMAREA*QSCAT

    ! Get the Gauss-Legendre quadrature abscissas and weights
    NRANK = min(MAXRANK, 2*NMIE)
    NQUAD = (NRANK + 2*NMIE + 2)/2
    if (NQUAD .le. LASTNQUAD .and. NQUAD .ge. nint(0.8*LASTNQUAD)) then
        NQUAD = LASTNQUAD
    else
        NQUAD = min(nint(1.25*NQUAD), MAXRANK)
        if (NQUAD .gt. MAXN) then
            stop 'Error [MIE_ONE]: MAXN exceeded by NQUAD.'
        end if
        call GAUSQUAD(NQUAD, MU, WTS)
        LASTNQUAD = NQUAD
    end if

    ! Calculate the phase function at the quadrature angles and then
    ! integrate over angle to get the Wigner d-function coefficients
    do L = 0, NRANK
        do J = 1, 6
            COEF(J, L) = 0.0d0
        end do
    end do
    do I = 1, NQUAD
        call MIEANGLE(NMIE, A, B, MU(I), P1, P2, P3, P4)

        ! Calculate the Wigner functions D20, D22, D2-2 and D00 for
        !   l=0 to NRANK for MU(I)
        call WIGNERFCT( MU(I), NRANK, 0, 0, D00)
        call WIGNERFCT( MU(I), NRANK, 2, 2, D22P)
        call WIGNERFCT(-MU(I), NRANK, 2, 2, D22M)
        call WIGNERFCT( MU(I), NRANK, 2, 0, D20)
        ! Calculate the series expansion coefficients
        do L = 0, NRANK
            F = (-1.0d0)**L
            COEF(1, L) = COEF(1, L) + (L + 0.5d0)*WTS(I)*P1*D00(L)
            COEF(2, L) = COEF(2, L) + (L + 0.5d0)*WTS(I)*(P1 + P3)*D22P(L)
            COEF(3, L) = COEF(3, L) + (L + 0.5d0)*WTS(I)*(P1 - P3)*F*D22M(L)
            COEF(4, L) = COEF(4, L) + (L + 0.5d0)*WTS(I)*P3*D00(L)
            COEF(5, L) = COEF(5, L) - (L + 0.5d0)*WTS(I)*P2*D20(L)
            COEF(6, L) = COEF(6, L) - (L + 0.5d0)*WTS(I)*P4*D20(L)
        end do
    end do
    do L = 0, NRANK
        A2 = 0.5d0*(COEF(2, L) + COEF(3, L))
        A3 = 0.5d0*(COEF(2, L) - COEF(3, L))
        COEF(2, L) = A2
        COEF(3, L) = A3
    end do
    do L = 0, NRANK
        do J = 1, 6
            WIGCOEF(J, L) = SNGL((WAVELENGTH**2/PI)*COEF(J, L))
        end do
    end do

    return
end subroutine

subroutine MIECALC(NTERMS, X, MN, A, B)
    !   MIECALC calculates the complex Mie coefficients An and Bn
    ! given the dimensionless size parameter X and the complex
    ! index of refraction (Mre,Mim).  The number of terms calculated
    ! is given by NTERMS unless NTERMS <= 0 or in which case the
    ! appropriate number is calculated and returned in NTERMS.
    implicit none
    integer NTERMS
    double precision X
    double complex MN, A(*), B(*)
    integer MAXTERMS, NSTOP, N, NN
    parameter(MAXTERMS=10000)
    double precision PSIN, PSIM, CHIN, CHIM, TMP
    double precision dcos, dsin
    double complex M, Y, D(MAXTERMS + 15), XIN, XIM, CTMP
    double complex DCMPLX

    ! If NTERMS is not specified calculate it
    NSTOP = X + 4.0*X**0.3334 + 2
    if (NTERMS .le. 0) NTERMS = NSTOP
    if (NTERMS .gt. MAXTERMS) then
        stop 'Error [MIECAL]: Mie calculation requires more terms than available.'
    end if

    ! Generate the Dn's by down recurrence  D = d(log(PSI(y)))/dy
    M = DCONJG(MN)
    Y = M*X
    NN = NTERMS + 15
    D(NN) = DCMPLX(0.0d0, 0.0d0)
    do N = NN, 2, -1
        D(N - 1) = N/Y - 1.0d0/(D(N) + N/Y)
    end do

    ! Generate the PSIn's and XIn'S by upward recurrence
    !   and calculate the An's and Bn's from them.
    !   (PSIN = PSI(n), PSIM = PSI(n-1), same for CHI)
    PSIM =  dcos(X)
    PSIN =  dsin(X)
    CHIM = -dsin(X)
    CHIN =  dcos(X)
    do N = 1, NTERMS
        TMP = PSIN
        PSIN = (2*N - 1)/X*PSIN - PSIM
        PSIM = TMP
        TMP = CHIN
        CHIN = (2*N - 1)/X*CHIN - CHIM
        CHIM = TMP
        XIN = DCMPLX(PSIN, -CHIN)
        XIM = DCMPLX(PSIM, -CHIM)
        CTMP = D(N)/M + N/X
        A(N) = (CTMP*PSIN - PSIM)/(CTMP*XIN - XIM)
        CTMP = M*D(N) + N/X
        B(N) = (CTMP*PSIN - PSIM)/(CTMP*XIN - XIM)
    end do

    return
end subroutine

subroutine MIECROSS(NTERMS, X, A, B, QEXT, QSCAT)
    !   MIECROSS calculates the extinction, scattering, and
    ! backscatter efficiencies given the Mie coefficients An and Bn
    ! and the size parameter X.
    implicit none
    integer NTERMS
    double precision X, QEXT, QSCAT
    double complex A(*), B(*)
    integer N
    double precision SUM1, SUM2, DREAL

    SUM1 = 0.0d0
    SUM2 = 0.0d0
    do N = 1, NTERMS
        SUM1 = SUM1 + (2*N + 1)*(DREAL(A(N)) + DREAL(B(N)))
        SUM2 = SUM2 + (2*N + 1)*(DREAL(A(N)*DCONJG(A(N)))&
            & + DREAL(B(N)*DCONJG(B(N))))
    end do
    QEXT  = 2.0d0/X**2*SUM1
    QSCAT = 2.0d0/X**2*SUM2

    return
end subroutine

subroutine MIEANGLE(NTERMS, A, B, MU, P1, P2, P3, P4)
    !   MIEANGLE calculates the scattering matrix elements (P1,P2,P3,P4
    ! for a particular value of MU (cos(theta)) from the Mie coefficien
    ! An's and Bn's.  The matrix elements are for the Stokes vector
    ! (I,Q,U,V) and are calculated from the complex scattering
    ! amplitudes S1 and S2.
    implicit none
    integer NTERMS
    double precision MU, P1, P2, P3, P4
    double complex A(*), B(*)
    integer N
    double precision TMP, PIN, PIM, TAUN, C
    double complex S1, S2

    S1 = DCMPLX(0.0, 0.0)
    S2 = DCMPLX(0.0, 0.0)
    ! Sum up the series using the An's and Bn's
    PIN = 1.0d0
    PIM = 0.0d0
    do N = 1, NTERMS
        TAUN = N*MU*PIN - (N + 1)*PIM
        ! Calculate the scattering functions at +mu and -mu
        !   using the PIn's and the TAUn's.
        C = (2*N + 1)/DFLOAT(N*(N + 1))
        S1 = S1 + C*(A(N)*PIN + B(N)*TAUN)
        S2 = S2 + C*(B(N)*PIN + A(N)*TAUN)
        ! Calculate the angular function PIn by up recurrence
        TMP = PIN
        PIN = ((2*N + 1)*MU*PIN - (N + 1)*PIM)/N
        PIM = TMP
    end do
    ! Calculate the first Stokes parameter scattering matrix eleme
    P1 = 0.5d0*(CDABS(S2)**2 + CDABS(S1)**2)
    P2 = 0.5d0*(CDABS(S2)**2 - CDABS(S1)**2)
    P3 =  DREAL(DCONJG(S1)*S2)
    P4 = -DIMAG(DCONJG(S1)*S2)
    return
end subroutine

subroutine WIGNERFCT(X, NRANK, M, M1, DMM1)
    ! The routine computes the DMM1N vector coefficients for
    ! M >= 0, M1 >= 0,  N = 0,...,NRANK, and -1 < X= COS(BETA) < 1
    implicit none
    integer NRANK, M, M1
    double precision X
    double precision DMM1(0:NRANK)
    integer N0, N
    double precision FACT1, FACT2, DMM1_N0

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
        DMM1(N0) = DMM1_N0(X, M, M1)
        do N = N0, NRANK - 1
            FACT1 = dble(N*(N + 1))*X - dble(M*M1)
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M**2))
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M1**2))
            FACT1 = FACT1*dble(2*N + 1)/dble(N)
            FACT2 = DSQRT(dble(N**2 - M**2))&
                & *DSQRT(dble(N**2 - M1**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M1**2))
            FACT2 = FACT2*dble(N + 1)/dble(N)
            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    end if
    return
end subroutine

double precision function DMM1_N0(X, M, M1)
    ! The routine computes the Wigner functions for
    ! M >= 0, M1 >= 0, AND N0 = MAX(M,M1)
    implicit none
    integer M, M1
    double precision X
    integer P, MAXM, MINM
    double precision PROD, FACT, CMM1, DMM1N0

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
    DMM1_N0 = DMM1N0
    return
end function

subroutine GAUSQUAD(N, XA, WT)
    ! Generates the abscissas (X) and weights (W) for an N point
    ! Gauss-Legendre quadrature.
    implicit none
    integer N
    double precision XA(*), WT(*)
    integer K, I, J, L
    double precision :: X, XP, PL, PL1, PL2, DPL, TINY
    parameter(TINY=3.0d-13)
    double precision, parameter :: PI = 3.141592653589793238462643383279502884197d0

    K = (N + 1)/2
    do 130 J = 1, K
        X = cos(PI*(J - 0.25)/(N + 0.5))
        I = 0
100     continue
        PL1 = 1
        PL = X
        do 120 L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
120     continue
        DPL = N*(X*PL - PL1)/(X*X - 1)
        XP = X
        X = XP - PL/DPL
        I = I + 1
        if (abs(X - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J) = -X
        XA(N - J + 1) = X
        WT(J) = 2.0d0/((1.0d0 - X*X)*DPL*DPL)
        WT(N - J + 1) = WT(J)
130 end do

    return
end subroutine
