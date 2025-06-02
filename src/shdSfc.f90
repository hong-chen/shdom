module shd_sfc

    use globals
    use globals_prp
    use shd_util

    implicit none
    private

    public :: shd_sfc__interpSfc
    public :: shd_sfc__brdfSurface

    ! surface brdf (including ocean brdf )
    !╭────────────────────────────────────────────────╮!
    public :: shd_sfc__brdfFresnel
    public :: shd_sfc__brdfDiner
    public :: shd_sfc__brdfRPV
    public :: shd_sfc__brdfLSRT
    public :: shd_sfc__brdfJiao
    public :: shd_sfc__brdfOcean
    !╰────────────────────────────────────────────────╯!

    ! ocean brdf sub-modules
    !╭────────────────────────────────────────────────╮!
    private :: morcasiwat
    private :: indwat
    private :: sunglint
    private :: Fresnel
    private :: glitalbe
    private :: getbound
    !╰────────────────────────────────────────────────╯!

contains

subroutine shd_sfc__interpSfc(&
    &NBOTPTS, BCPTR, GRIDPOS,&
    &SRCTYPE, WAVENO, WAVELEN, UNITS,&
    &NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
    &SFCGRIDPARMS)

    !   Bilinearly interpolates the surface properties from the regular
    ! surface grid to the boundary points of the adaptive grid.
    ! The locations of the NBOTPTS grid points on the bottom boundary
    ! (pointed to by BCPTR) are in the GRIDPOS array.  The input regular
    ! surface info is specified by the number of grid lines (NXSFC,NYSFC),
    ! the grid spacing (DELXSFC,DELYSFC), the number of parameters at
    ! each grid point (NSFCPAR), and the surface parameters (SFCPARMS).
    ! The output interpolated surface parameters at the NBOTPTS boundary
    ! points are in SFCGRIDPARMS.  The first parameter is input as
    ! temperature, but converted to Planck function for output.

    integer NBOTPTS, BCPTR(*), NXSFC, NYSFC, NSFCPAR
    real WAVENO(2), WAVELEN
    real DELXSFC, DELYSFC
    real GRIDPOS(3, *), SFCPARMS(NSFCPAR, NXSFC+1, NYSFC+1)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IBC, IX, IY, J
    real RX, RY, U, V, TEMP, PLANCK

    ! Loop over all bottom points
    do IBC = 1, NBOTPTS

        I = BCPTR(IBC)

        if ((NPX_ .gt. 1) .and. (NPY_ .gt. 1) .and.&
            &(NXSFC .gt. 1) .and. (NYSFC .gt. 1) .and.&
            &((abs(DELX_-DELXSFC) .gt. 1.0e-6) .or.&
            &(abs(DELY_-DELYSFC) .gt. 1.0e-6))) then
            ! added half dx,dy (DELX_/2.0-DELXSFC/2.0, DELY_/2.0-DELYSFC/2.0) to
            ! find surface point that is closest to grid center
            !╭────────────────────────────────────────────────╮!
            RX = (GRIDPOS(1, I) + DELX_/2.0 - DELXSFC/2.0) / DELXSFC
            RY = (GRIDPOS(2, I) + DELY_/2.0 - DELYSFC/2.0) / DELYSFC
            !╰────────────────────────────────────────────────╯!
        else
            RX = (GRIDPOS(1, I)) / DELXSFC
            RY = (GRIDPOS(2, I)) / DELYSFC
        end if

        IX = max(1, min(NXSFC, int(RX) + 1))
        IY = max(1, min(NYSFC, int(RY) + 1))
        U = max(0.0, min(1.0, RX - (IX - 1)))
        V = max(0.0, min(1.0, RY - (IY - 1)))

        do J = 1, NSFCPAR
            ! bi-linear interpolation
            !╭────────────────────────────────────────────────╮!
            SFCGRIDPARMS(J, IBC) =&
                &  (1.0 - U)*(1.0 - V)*SFCPARMS(J, IX  , IY)&
                &+ (1.0 - U)*       V *SFCPARMS(J, IX  , IY+1)&
                &+        U *(1.0 - V)*SFCPARMS(J, IX+1, IY)&
                &+        U *       V *SFCPARMS(J, IX+1, IY+1)
            !╰────────────────────────────────────────────────╯!

            ! nearest interpolation
            !╭────────────────────────────────────────────────╮!
            !SFCGRIDPARMS(J, IBC) = SFCPARMS(J, IX, IY)
            !╰────────────────────────────────────────────────╯!
        end do

        if (SRCTYPE .eq. 'S') then
            PLANCK = 0.0
        else
            TEMP = SFCGRIDPARMS(1, IBC)
            call shd_util__computePlanck(TEMP, UNITS, WAVENO, WAVELEN, PLANCK)
        end if

        SFCGRIDPARMS(1, IBC) = PLANCK

    end do

    return

end subroutine shd_sfc__interpSfc

subroutine shd_sfc__brdfSurface(SFCTYPE, REFPARMS, WAVELEN,&
    &MU2, PHI2, MU1, PHI1, NSTOKES, REFLECT)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Returns the reflection matrix for the general bidirectional
    ! reflection distribution function of the specified type (SFCTYPE).
    ! The incident direction is (MU1,PHI1), and the outgoing direction
    ! is (MU2,PHI2) (MU is cosine of zenith angle, and PHI is the azimuth
    ! angle in radians).  The incident directions have mu<0 (downward),
    ! while the outgoing directions have mu>0 (upward). The reflection
    ! function is normalized so that for a Lambertian surface (uniform
    ! reflection) the returned value (REFLECT) is simply the albedo.
    ! This routine calls the desired BRDF function, passing the
    ! appropriate parameters.  More BRDF surface types may be added easily
    ! by putting in the appropriate function calls.
    !
    !      Type        Parameters
    ! L  Lambertian    albedo
    ! W  WaveFresnel   Real, Imaginary index of refraction, wind speed
    ! D  Diner et al   a, k, b, zeta, sigma
    ! O  Ocean         Wind speed (m/s), Pigmentation (mg/m^3)
    ! R  RPV-original  rho0, k, Theta
    ! T  LSRT          k_iso, k_geo, k_vol
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer :: NSTOKES
    integer :: K1, K2
    real    :: REFPARMS(*), WAVELEN, MU1, PHI1, MU2, PHI2, REFLECT(4, 4)
    real    :: PHIDIFF
    character(len=1) :: SFCTYPE

    if (SFCTYPE .eq. 'L' .or. SFCTYPE .eq. 'l') then

        ! L or l: Lambertian surface BRDF is constant.
        !   (for testing, as done more efficiently by Lambertian routine
        REFLECT(1:NSTOKES, 1:NSTOKES) = 0.0
        REFLECT(1, 1) = REFPARMS(1)

    else if (SFCTYPE .eq. 'W' .or. SFCTYPE .eq. 'w') then

        ! W: Fresnel reflection with Gaussian distribution of slopes
        call shd_sfc__brdfFresnel(REFPARMS(1), REFPARMS(2),&
            &REFPARMS(3), MU1, MU2, PHI1, PHI2, NSTOKES, REFLECT)

    else if (SFCTYPE .eq. 'D' .or. SFCTYPE .eq. 'd') then

        ! D: Diner et al., depolarizing modified RPV + Fresnel reflection
        !   from microfacets with uniform or Gaussian distribution of s
        PHIDIFF = PHI2 - PHI1
        call shd_sfc__brdfDiner(REFPARMS(1), REFPARMS(2), REFPARMS(3),&
            &REFPARMS(4), REFPARMS(5),&
            &-MU1, MU2, PHIDIFF, NSTOKES, REFLECT)

    else if (SFCTYPE .eq. 'R' .or. SFCTYPE .eq. 'r') then

        ! R: Rahman, Pinty, and Verstraete
        if (NSTOKES .gt. 1) then
            stop 'Error [_brdfSurface]: RPV BRDF is only for unpolarized case'
        end if

        PHIDIFF = PHI1-PHI2-PI_
        call shd_sfc__brdfRPV(REFPARMS(1), REFPARMS(2),&
            &REFPARMS(3), -MU1, MU2, PHIDIFF, REFLECT(1, 1))

    else if (SFCTYPE .eq. 'T' .or. SFCTYPE .eq. 't') then

        ! T: LSRT model
        if (NSTOKES .gt. 1) then
            stop 'Error [_brdfSurface]: LSRT BRDF is only for unpolarized case'
        end if

        PHIDIFF = PHI1-PHI2-PI_
        ! T: LSRT model
        call shd_sfc__brdfLSRT(REFPARMS(1), REFPARMS(2),&
            &REFPARMS(3), -MU1, MU2, PHIDIFF, REFLECT(1, 1))

    else if (SFCTYPE .eq. 'O' .or. SFCTYPE .eq. 'o') then

        ! O: Ocean BRDF from 6S modified by Norm Loeb
        if (NSTOKES .gt. 1) then
            stop 'Error [_brdfSurface]: Ocean BRDF from 6S is only for unpolarized case'
        end if

        PHIDIFF = PHI1-PHI2
        call shd_sfc__brdfOcean(REFPARMS(1), -1., REFPARMS(2), WAVELEN,&
            & -MU1, MU2, PHIDIFF, PHI1, REFLECT(1, 1))

    else
        stop 'Error [_brdfSurface]: Unknown BRDF type.'
    end if

    return
end subroutine shd_sfc__brdfSurface

subroutine shd_sfc__brdfFresnel(MRE, MIM, WINDSPEED,&
    &MUI, MUR, PHII, PHIR, NSTOKES, REFLECT)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !   Computes the polarized Fresnel reflection from a dielectric
    ! interface with a Gaussian distribution of slopes, including shadow
    ! The complex index of refraction of the dielectric is specified wit
    ! MRE,MIM (MIM is negative for an absorbing medium) and the wind spe
    ! is in m/s.  The incident direction is (MUI,PHII), and the outgoing
    ! direction is (MUR,PHIR) (MU is cosine of zenith angle with MU>0 an
    ! PHI is the azimuthal angle in radians).
    !   The code is courtesy of Michael Mishchenko, with some modificati
    ! The mean square surface slope (s^2) or SIGMA2 in Mishchenko and Tr
    ! is obtained from the Cox and Munk formula in the paper.
    !   For all formulas and definitions, see the paper:
    ! M. I. Mishchenko and L. D. Travis, 1997: Satellite retrieval
    ! of aerosol properties over the ocean using polarization as well as
    ! intensity of reflected sunlight.  J. Geophys. Res. 102, 16989-1701
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer NSTOKES
    real MRE, MIM, WINDSPEED, MUI, MUR, PHII, PHIR, REFLECT(4, 4)
    integer I, J
    real*8 SIGMA2, DMUI, DMUR, DCOSI, DSINI, DCOSR, DSINR, DSI, DSR
    real*8 VI1, VI2, VI3, VR1, VR2, VR3, UNIT1, UNIT2, UNIT3
    real*8 FACT1, FACTOR, XI1
    real*8 TI1, TI2, TI3, TR1, TR2, TR3, PI1, PI2, PI3, PR1, PR2, PR3
    real*8 PIKR, PRKI, TIKR, TRKI, E1, E2, E3, E4
    real*8 VP1, VP2, VP3, DMOD, RDZ2, RDZ4, DCOEFF, DEX
    real*8 AF, AF11, AF12, AF21, AF22
    real*8 P, S1, S2, S3, XI, XXI, DCOT, T1, T2
    real*8 SHADOWI, SHADOWR, SHADOW, DERFC
    complex*16 CN1, CN2, CXI2, C1, C2, CRPER, CRPAR
    complex*16 CF11, CF12, CF21, CF22
    complex*16 CI, C21, C22, CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

    CN1 = (1.0d0, 0.0d0)
    CN2 = cmplx(MRE, MIM)
    SIGMA2 = max(0.0005d0, 0.0015d0 + 0.00256d0*WINDSPEED)
    DMUI = DABS(dble(MUI))
    DMUR = dble(MUR)

    ! Cartesian components of the unit vectors of the incident and
    ! scattered beams
    if (DABS(DMUI - 1e0) .lt. 1d-10) DMUI = 0.999999999999d0
    if (DABS(DMUR - 1e0) .lt. 1d-10) DMUR = 0.999999999999d0
    DCOSI = DCOS(dble(PHII))
    DSINI = DSIN(dble(PHII))
    DCOSR = DCOS(dble(PHIR))
    DSINR = DSIN(dble(PHIR))
    DSI = DSQRT(1d0 - DMUI*DMUI)
    DSR = DSQRT(1d0 - DMUR*DMUR)
    VI1 = DSI*DCOSI
    VI2 = DSI*DSINI
    VI3 = -DMUI
    VR1 = DSR*DCOSR
    VR2 = DSR*DSINR
    VR3 = DMUR

    ! Local surface normal for specular reflection
    UNIT1 = VI1 - VR1
    UNIT2 = VI2 - VR2
    UNIT3 = VI3 - VR3
    FACT1 = UNIT1*UNIT1 + UNIT2*UNIT2 + UNIT3*UNIT3
    FACTOR = sqrt(1d0/FACT1)

    ! Fresnel reflection coefficients
    XI1 = FACTOR*(UNIT1*VI1 + UNIT2*VI2 + UNIT3*VI3)
    CXI2 = CDSQRT(1d0 - (1d0 - XI1*XI1)*CN1*CN1/(CN2*CN2))
    C1 = CN1*XI1
    C2 = CN2*CXI2
    CRPER = (C1 - C2)/(C1 + C2)
    C1 = CN2*XI1
    C2 = CN1*CXI2
    CRPAR = (C1 - C2)/(C1 + C2)

    ! Calculation of the amplitude scattering matrix
    TI1 = -DMUI*DCOSI
    TI2 = -DMUI*DSINI
    TI3 = -DSI
    TR1 = DMUR*DCOSR
    TR2 = DMUR*DSINR
    TR3 = -DSR
    PI1 = -DSINI
    PI2 = DCOSI
    PI3 = 0d0
    PR1 = -DSINR
    PR2 = DCOSR
    PR3 = 0d0
    PIKR = PI1*VR1 + PI2*VR2 + PI3*VR3
    PRKI = PR1*VI1 + PR2*VI2 + PR3*VI3
    TIKR = TI1*VR1 + TI2*VR2 + TI3*VR3
    TRKI = TR1*VI1 + TR2*VI2 + TR3*VI3
    E1 = PIKR*PRKI
    E2 = TIKR*TRKI
    E3 = TIKR*PRKI
    E4 = PIKR*TRKI
    CF11 = E1*CRPER + E2*CRPAR
    CF12 = -E3*CRPER + E4*CRPAR
    CF21 = -E4*CRPER + E3*CRPAR
    CF22 = E2*CRPER + E1*CRPAR

    ! Calculation of the Stokes reflection matrix
    VP1 = VI2*VR3 - VI3*VR2
    VP2 = VI3*VR1 - VI1*VR3
    VP3 = VI1*VR2 - VI2*VR1
    DMOD = (VP1*VP1 + VP2*VP2 + VP3*VP3)**2
    RDZ2 = UNIT3*UNIT3
    RDZ4 = RDZ2*RDZ2
    DEX = DEXP(-(UNIT1*UNIT1 + UNIT2*UNIT2)/(2*SIGMA2*RDZ2))
    DCOEFF = FACT1**2*DEX/(4*DMUI*DMUR*DMOD*RDZ4*2*SIGMA2)
    AF = 0.5d0*DCOEFF
    AF11 = CDABS(CF11)**2
    AF12 = CDABS(CF12)**2
    AF21 = CDABS(CF21)**2
    AF22 = CDABS(CF22)**2

    REFLECT(1, 1) = (AF11 + AF12 + AF21 + AF22)*AF
    if (NSTOKES .ge. 2) then
        REFLECT(1, 2) = (AF11 - AF12 + AF21 - AF22)*AF
        REFLECT(2, 1) = (AF11 - AF22 + AF12 - AF21)*AF
        REFLECT(2, 2) = (AF11 - AF12 - AF21 + AF22)*AF
    end if

    CI = (0.0d0, -1.0d0)
    C21 = conjg(CF21)
    C22 = conjg(CF22)
    CTTTP = CF11*conjg(CF12)
    CTTPT = CF11*C21
    CTTPP = CF11*C22
    CTPPT = CF12*C21
    CTPPP = CF12*C22
    CPTPP = CF21*C22

    if (NSTOKES .ge. 3) then
        REFLECT(1, 3) = (-CTTTP - CPTPP)*DCOEFF
        REFLECT(2, 3) = (-CTTTP + CPTPP)*DCOEFF
        REFLECT(3, 1) = (-CTTPT - CTPPP)*DCOEFF
        REFLECT(3, 2) = (-CTTPT + CTPPP)*DCOEFF
        REFLECT(3, 3) = (CTTPP + CTPPT)*DCOEFF
    end if

    if (NSTOKES .ge. 4) then
        REFLECT(1, 4) = -CI*(CTTTP + CPTPP)*DCOEFF
        REFLECT(2, 4) = -CI*(CTTTP - CPTPP)*DCOEFF
        REFLECT(3, 4) = CI*(CTTPP - CTPPT)*DCOEFF
        REFLECT(4, 1) = CI*(CTTPT + CTPPP)*DCOEFF
        REFLECT(4, 2) = CI*(CTTPT - CTPPP)*DCOEFF
        REFLECT(4, 3) = -CI*(CTTPP + CTPPT)*DCOEFF
        REFLECT(4, 4) = (CTTPP - CTPPT)*DCOEFF
    end if

    ! Shadowing
    P = DACOS(-1.0d0)
    S1 = DSQRT(2*SIGMA2/P)
    S3 = 1d0/(DSQRT(2*SIGMA2))
    S2 = S3*S3
    DCOT = DMUI/sqrt(1d0 - DMUI**2)
    T1 = DEXP(-S2*DCOT**2)
    T2 = DERFC(DCOT*S3)
    SHADOWI = 0.5d0*(S1*T1/DCOT - T2)
    DCOT = DMUR/DSQRT(1d0 - DMUR**2)
    T1 = DEXP(-S2*DCOT**2)
    T2 = DERFC(DCOT*S3)
    SHADOWR = 0.5d0*(S1*T1/DCOT - T2)
    SHADOW = 1d0/(1d0 + SHADOWI + SHADOWR)
    REFLECT(1:NSTOKES, 1:NSTOKES) = REFLECT(1:NSTOKES, 1:NSTOKES)*SHADOW

    return

end subroutine shd_sfc__brdfFresnel

subroutine shd_sfc__brdfDiner(A, K, B, ZETA, SIGMA,&
    &MU1, MU2, PHI, NSTOKES, REFLECT)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !   Computes polarized reflection from the Diner et al. (2012) model
    ! which two terms: 1) completely depolarizing volumetric scattering
    ! the modified RPV model, and 2) Fresnel reflection from randomly
    ! oriented microfacets with either a uniform or Gaussian distributio
    ! of slopes.  The three parameters of the modified-RPV model are A,
    ! and B.  The weight given to the Fresnel-microfacet term is ZETA (i
    ! ZETA=0 then only the REFLECT(1,1) element is non-zero with the
    ! modified-RPV model).  SIGMA is the standard deviation of microface
    ! slopes for the Gaussian orientation distribution, but if SIGMA<=0
    ! the uniform orientation distribution is used.  The complex index o
    ! refraction of the Fresnel microfacets is specified in this subrout
    ! The Diner et al. (2012) model leaves out the hotspot factor of the
    ! modified RPV model (and the original RPV), but the hotspot factor
    ! included here.
    !   The cosine of the incident zenith angle is MU1 (MU1>0), the
    ! cosine of the outgoing zenith angle is MU2 (MU2>0), and the relati
    ! azimuthal angle is PHI=PHI2-PHI1 (radians).
    !   For all formulas and definitions, see the paper:
    ! Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
    ! V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
    ! Exploration of a polarized surface bidirectional reflectance model
    ! using the ground-based multiangle spectropolarimetric imager.
    ! Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer NSTOKES
    real A, K, B, ZETA, SIGMA, MU1, MU2, PHI
    real REFLECT(4, 4)
    integer K1, K2
    real SINTH1, SINTH2, COSPHI, COSSCATANG, TAN1, TAN2, CAPG, HOT
    real GAMMA, COSGAMMA, F11, F12, F33, F34
    real COSBETA, H, SINPHI, ALPHA1, ALPHA2
    real COS2ALPHA1, SIN2ALPHA1, COS2ALPHA2, SIN2ALPHA2
    complex SFCINDEX, EPSILON, D, RP, RS

    ! Another free parameter of the Fresnel microfacet model is the
    !   index of refraction specified here
    SFCINDEX = cmplx(1.5, 0.0)

    REFLECT(1:NSTOKES, 1:NSTOKES) = 0.0

    ! The modified RPV model (Martonchik et al., 1998, IEEE-TGARS, 36,
    !   is implemented here.  Equation 6 in Diner et al. leaves out the h
    SINTH1 = sqrt(1.0 - MU1**2)
    SINTH2 = sqrt(1.0 - MU2**2)
    COSPHI = cos(PHI)
    COSSCATANG = -MU1*MU2 + SINTH1*SINTH2*COSPHI
    COSSCATANG = min(1.0, max(-1.0, COSSCATANG))
    TAN1 = SINTH1/MU1
    TAN2 = SINTH2/MU2
    CAPG = sqrt(abs(TAN1**2 + TAN2**2 + 2*TAN1*TAN2*COSPHI))
    HOT = 1 + (1 - A)/(1 + CAPG)
    REFLECT(1, 1) = A*((MU1 + MU2)*MU1*MU2)**(K - 1)*exp(B*COSSCATANG)
    REFLECT(1, 1) = REFLECT(1, 1)*HOT
    if (ZETA .lt. 0.0) return

    ! Calculate incident/reflectance angle relative to microfacet (eq 1
    GAMMA = 0.5d0*acos(-COSSCATANG)
    COSGAMMA = cos(GAMMA)

    ! Calculate the Fresnel reflection coefficients (equivalent to eq 1
    !   [D = n*cos(gamma')]
    EPSILON = SFCINDEX**2
    D = sqrt(EPSILON - 1.0 + COSGAMMA**2)
    RP = (EPSILON*COSGAMMA - D)/(EPSILON*COSGAMMA + D)
    RS = (COSGAMMA - D)/(COSGAMMA + D)

    ! Calculate the Mueller matrix elements for Fresnel reflection (eq
    F11 = 0.5d0*(CABS(RP)**2 + CABS(RS)**2)
    F12 = 0.5d0*(CABS(RP)**2 - CABS(RS)**2)
    F33 = real(RP*conjg(RS))
    F34 = -IMAG(RP*conjg(RS))

    ! Calculate facet tilt angle, beta (eq 16b)
    COSBETA = 0.5d0*(MU1 + MU2)/COSGAMMA

    ! Calculate the factor before the M*P*M matrix including the distri
    !   of microfacet orientations (Gaussian for SIGMA>0) (eq 13 & 14)
    !   Note: Diner's BRDF is multiplied by pi here
    if (SIGMA .gt. 0.0) then
        H = ZETA*exp(-0.5d0*(1/COSBETA**2 - 1)/SIGMA**2)&
            & /(8*SIGMA**2*MU2*MU1*COSBETA**4)
    else
        H = ZETA/(8*MU2*MU1*COSBETA)
    end if

    REFLECT(1, 1) = REFLECT(1, 1) + H*F11

    if (NSTOKES .ge. 2) then
        ! Calculate the polarization rotation angles (eq 20 & 21)
        SINPHI = sin(PHI)
        ALPHA1 = atan(SINTH2*SINPHI/(MU2*SINTH1 + SINTH2*MU1*COSPHI))
        ALPHA2 = atan(SINTH1*SINPHI/(SINTH2*MU1 + MU2*SINTH1*COSPHI))
        COS2ALPHA1 = cos(2*ALPHA1)
        SIN2ALPHA1 = sin(2*ALPHA1)
        COS2ALPHA2 = cos(2*ALPHA2)
        SIN2ALPHA2 = sin(2*ALPHA2)

        ! Calculate the rotated Mueller matrix (M*P*M) for Fresnel reflec
        !   and multiply by the H factor to get the final BRDF
        REFLECT(1, 2) = H*F12*COS2ALPHA1
        REFLECT(2, 1) = H*F12*COS2ALPHA2
        REFLECT(2, 2) = H*(F11*COS2ALPHA1*COS2ALPHA2&
            & + F33*SIN2ALPHA1*SIN2ALPHA2)
    end if

    if (NSTOKES .ge. 3) then
        REFLECT(1, 3) = -H*F12*SIN2ALPHA1
        REFLECT(2, 3) = H*(-F11*SIN2ALPHA1*COS2ALPHA2&
            & + F33*COS2ALPHA1*SIN2ALPHA2)
        REFLECT(3, 1) = -H*F12*SIN2ALPHA2
        REFLECT(3, 2) = H*(-F11*COS2ALPHA1*SIN2ALPHA2&
            & + F33*SIN2ALPHA1*COS2ALPHA2)
        REFLECT(3, 3) = H*(F11*SIN2ALPHA1*SIN2ALPHA2&
            & + F33*COS2ALPHA1*COS2ALPHA2)
    end if

    if (NSTOKES .ge. 4) then
        REFLECT(1, 4) = 0.0
        REFLECT(2, 4) = H*F34*SIN2ALPHA2
        REFLECT(3, 4) = H*F34*COS2ALPHA2
        REFLECT(4, 1) = 0.0
        REFLECT(4, 2) = -H*F34*SIN2ALPHA1
        REFLECT(4, 3) = -H*F34*COS2ALPHA1
        REFLECT(4, 4) = H*F33
    end if

    return

end subroutine shd_sfc__brdfDiner

subroutine shd_sfc__brdfRPV(RHO0, K, THETA, MU1, MU2, PHI, REFLECT)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !   Computes the Rahman, Pinty, Verstraete BRDF.  The incident
    ! and outgoing cosine zenith angles are MU1 and MU2, respectively,
    ! and the relative azimuthal angle is PHI.  In this case the incident
    ! direction is where the radiation is coming from (i.e. opposite of
    ! discrete ordinate), so MU1>0 and the hot spot is MU2=MU1 and PHI=0
    ! The reference is:
    !   Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere
    !   Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable
    !   With NOAA Advanced Very High Resolution Radiometer Data,
    !   J. Geophys. Res., 98, 20791-20801.
    !
    ! Implementation modified by Alexei Lyapustin to prevent negative
    ! values in SQRT function in calculating CAPG, and to prevent
    ! unphysical growth of RPV BRDF at small mu or mu0.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real, intent(in)  :: RHO0, K, THETA, MU1, MU2, PHI
    real, intent(out) :: REFLECT
    real :: mu_min, x1, x2
    real :: M, F, H, COSPHI, SIN1, SIN2, COSG, TAN1, TAN2, CAPG

    ! calculate M
    !╭────────────────────────────────────────────────╮!
    mu_min = 0.03
    x1 = MU1
    if (x1 .lt. mu_min) x1 = mu_min
    x2 = MU2
    if (x2 .lt. mu_min) x2 = mu_min
    M = (x1*x2*(x1 + x2))**(K - 1)
    !╰────────────────────────────────────────────────╯!

    ! calculate F
    !╭────────────────────────────────────────────────╮!
    COSPHI = cos(PHI)
    SIN1 = sqrt(1.0 - x1**2)
    SIN2 = sqrt(1.0 - x2**2)
    COSG = x1*x2 + SIN1*SIN2*COSPHI
    F = (1 - THETA**2)/(1 + 2*THETA*COSG + THETA**2)**1.5
    !╰────────────────────────────────────────────────╯!

    ! calculate H
    !╭────────────────────────────────────────────────╮!
    TAN1 = SIN1/x1
    TAN2 = SIN2/x2
    CAPG = sqrt(abs(TAN1**2 + TAN2**2 - 2*TAN1*TAN2*COSPHI))
    H = 1 + (1 - RHO0)/(1 + CAPG)
    !╰────────────────────────────────────────────────╯!

    REFLECT = RHO0*M*F*H

    return

end subroutine shd_sfc__brdfRPV

subroutine shd_sfc__brdfLSRT(F_ISO, F_GEO, F_VOL, MU1, MU2, PHI, REFLECT)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Computes the RossLi BRDF.
    ! MU1: cosine of zenith angle of the incident light (solar zenith angle)
    ! MU2: cosine of zenith angle of the receiver (viewing zenith angle)
    ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
    !
    ! When MU1>0 and the hot spot is MU2=MU1 and PHI=0
    !
    ! Code reference:
    !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
    ! Formula reference:
    !   MODIS BRDF/Albedo Product ATBD v5.0, equations 37-44
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(4), intent(in)  :: F_ISO, F_GEO, F_VOL, MU1, MU2, PHI
    real(4), intent(out) :: REFLECT
    real(4) :: COS1, COS2, SIN1, SIN2, TAN1, TAN2, COS_PHI, SIN_PHI
    real(4) :: K_GEO, K_VOL
    real(4) :: XI, COS_XI, SIN_XI
    real(4) :: B2R, H2B
    real(4) :: THETA1_, COS1_, SIN1_, TAN1_, SEC1_
    real(4) :: THETA2_, COS2_, SIN2_, TAN2_, SEC2_
    real(4) :: COS_XI_, D_2, COS_T, SIN_T, T, O
    real(4) :: COS_MIN


    ! calculate trigonometry
    !╭────────────────────────────────────────────────╮!
    COS_MIN = 0.03

    COS1 = MU1
    COS2 = MU2

    if (COS1 .lt. COS_MIN) COS1 = COS_MIN
    if (COS2 .lt. COS_MIN) COS2 = COS_MIN

    SIN1 = sqrt(1.0 - COS1**2)
    SIN2 = sqrt(1.0 - COS2**2)
    TAN1 = SIN1/COS1
    TAN2 = SIN2/COS2
    COS_PHI = cos(PHI)
    SIN_PHI = sin(PHI)
    !╰────────────────────────────────────────────────╯!

    ! K_VOL
    !╭────────────────────────────────────────────────╮!
    COS_XI = (COS1*COS2) + (SIN1*SIN2*COS_PHI)
    if (COS_XI .lt. -1.0) COS_XI = -1.0
    if (COS_XI .gt.  1.0) COS_XI =  1.0

    SIN_XI = sqrt(1.0 - COS_XI**2)
    XI = acos(COS_XI)

    K_VOL = ((PI_/2.0 - XI)*COS_XI + SIN_XI)/(COS1+COS2) - (PI_/4.0)
    !╰────────────────────────────────────────────────╯!

    ! K_GEO
    !╭────────────────────────────────────────────────╮!
    B2R = 1.0
    H2B = 2.0

    THETA1_ = atan(B2R*TAN2)
    COS1_ = cos(THETA1_)
    if (COS1_ .lt. -1.0) COS1_ = -1.0
    if (COS1_ .gt.  1.0) COS1_ =  1.0
    SIN1_ = sqrt(1.0 - COS1_**2)
    TAN1_ = SIN1_/COS1_
    SEC1_ = 1.0/COS1_

    THETA2_ = atan(B2R*TAN1)
    COS2_ = cos(THETA2_)
    if (COS2_ .lt. -1.0) COS2_ = -1.0
    if (COS2_ .gt.  1.0) COS2_ =  1.0
    SIN2_ = sqrt(1.0 - COS2_**2)
    TAN2_ = SIN2_/COS2_
    SEC2_ = 1.0/COS2_

    COS_XI_ = (COS1_*COS2_) + (SIN1_*SIN2_*COS_PHI)

    D_2 = TAN1_**2 + TAN2_**2 - 2.0*TAN1_*TAN2_*COS_PHI
    if (D_2 .lt. 0.0) D_2 = 0.0

    COS_T = H2B * (sqrt(D_2 + (TAN1_*TAN2_*SIN_PHI)**2))/(SEC1_ + SEC2_)
    if (COS_T .lt. -1.0) COS_T = -1.0
    if (COS_T .gt.  1.0) COS_T =  1.0
    SIN_T = sqrt(1.0 - COS_T**2)
    T = acos(COS_T)

    O = (1.0/PI_) * (T - SIN_T*COS_T) * (SEC1_ + SEC2_)

    K_GEO = O - SEC1_ - SEC2_ + 0.5d0*(1.0+COS_XI_)*SEC1_*SEC2_
    !╰────────────────────────────────────────────────╯!

    REFLECT = F_ISO + F_VOL*K_VOL + F_GEO*K_GEO

    if (REFLECT .lt. 0.0) REFLECT = 0.0

    return

end subroutine shd_sfc__brdfLSRT


! subroutine shd_sfc__brdfLSRTKZ(F_ISO, F_GEO, F_VOL, F_KZ, MU1, MU2, PHI, REFLECT)

!     !╭────────────────────────────────────────────────────────────────────────────╮!
!     ! Computes the RossLi-KokhanovskyZege BRDF (combines LSRT model above with a fourth kernel for snow from Kokhanovsky-Zege).
!     ! MU1: cosine of zenith angle of the incident light (solar zenith angle)
!     ! MU2: cosine of zenith angle of the receiver (viewing zenith angle)
!     ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
!     !
!     ! When MU1>0 and the hot spot is MU2=MU1 and PHI=0
!     !
!     ! Code reference:
!     !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
!     !
!     ! Formula reference:
!     !    KZ04: Kokhanovsky & Zege (2004), doi:10.1364/ao.43.001589
!     !    KB12: Kokhanovsky & Breon (2012), doi:10.1109/lgrs.2012.2185775
!     !╰────────────────────────────────────────────────────────────────────────────╯!

!     real(4), intent(in)  :: F_ISO, F_GEO, F_VOL, F_KZ, MU1, MU2, PHI
!     real(4), intent(out) :: REFLECT
!     real(4) :: COS1, COS2, SIN1, SIN2, TAN1, TAN2, COS_PHI, SIN_PHI
!     real(4) :: K_GEO, K_VOL, K_KZ
!     real(4) :: XI, COS_XI, SIN_XI
!     real(4) :: B2R, H2B
!     real(4) :: THETA1_, COS1_, SIN1_, TAN1_, SEC1_
!     real(4) :: THETA2_, COS2_, SIN2_, TAN2_, SEC2_
!     real(4) :: COS_XI_, D_2, COS_T, SIN_T, T, O
!     real(4) :: COS_MIN
!     parameter(4) :: a, b, c, b_snow_grain

!     ! calculate trigonometry
!     !╭────────────────────────────────────────────────╮!
!     COS_MIN = 0.03

!     COS1 = MU1
!     COS2 = MU2

!     if (COS1 .lt. COS_MIN) COS1 = COS_MIN
!     if (COS2 .lt. COS_MIN) COS2 = COS_MIN

!     SIN1 = sqrt(1.0 - COS1**2)
!     SIN2 = sqrt(1.0 - COS2**2)
!     TAN1 = SIN1/COS1
!     TAN2 = SIN2/COS2
!     COS_PHI = cos(PHI)
!     SIN_PHI = sin(PHI)
!     !╰────────────────────────────────────────────────╯!

!     ! K_VOL
!     !╭────────────────────────────────────────────────╮!
!     COS_XI = (COS1*COS2) + (SIN1*SIN2*COS_PHI)
!     if (COS_XI .lt. -1.0) COS_XI = -1.0
!     if (COS_XI .gt.  1.0) COS_XI =  1.0

!     SIN_XI = sqrt(1.0 - COS_XI**2)
!     XI = acos(COS_XI)

!     K_VOL = ((PI_/2.0 - XI)*COS_XI + SIN_XI)/(COS1+COS2) - (PI_/4.0)
!     !╰────────────────────────────────────────────────╯!

!     ! K_GEO
!     !╭────────────────────────────────────────────────╮!
!     B2R = 1.0
!     H2B = 2.0

!     THETA1_ = atan(B2R*TAN2)
!     COS1_ = cos(THETA1_)
!     if (COS1_ .lt. -1.0) COS1_ = -1.0
!     if (COS1_ .gt.  1.0) COS1_ =  1.0
!     SIN1_ = sqrt(1.0 - COS1_**2)
!     TAN1_ = SIN1_/COS1_
!     SEC1_ = 1.0/COS1_

!     THETA2_ = atan(B2R*TAN1)
!     COS2_ = cos(THETA2_)
!     if (COS2_ .lt. -1.0) COS2_ = -1.0
!     if (COS2_ .gt.  1.0) COS2_ =  1.0
!     SIN2_ = sqrt(1.0 - COS2_**2)
!     TAN2_ = SIN2_/COS2_
!     SEC2_ = 1.0/COS2_

!     COS_XI_ = (COS1_*COS2_) + (SIN1_*SIN2_*COS_PHI)

!     D_2 = TAN1_**2 + TAN2_**2 - 2.0*TAN1_*TAN2_*COS_PHI
!     if (D_2 .lt. 0.0) D_2 = 0.0

!     COS_T = H2B * (sqrt(D_2 + (TAN1_*TAN2_*SIN_PHI)**2))/(SEC1_ + SEC2_)
!     if (COS_T .lt. -1.0) COS_T = -1.0
!     if (COS_T .gt.  1.0) COS_T =  1.0
!     SIN_T = sqrt(1.0 - COS_T**2)
!     T = acos(COS_T)

!     O = (1.0/PI_) * (T - SIN_T*COS_T) * (SEC1_ + SEC2_)

!     K_GEO = O - SEC1_ - SEC2_ + 0.5d0*(1.0+COS_XI_)*SEC1_*SEC2_
!     !╰────────────────────────────────────────────────╯!

!     ! K_KZ
!     !╭────────────────────────────────────────────────╮!
!     PHI_ = 180.0 - PHI ! because KZ defines RAA with sun at RAA=0, forward scattering denotes 0-90 RAA.

!     ! ensure PHI_ is in [0, 360] range
!     PHI_ = change_range(PHI_, 0.0, 360.0)
!     if (PHI_ .eq. 360.0) PHI_ = 0.0

!     ! set constants for parameterized escape function
!     a = 1.247
!     b = 1.186
!     c = 5.157
!     b_snow_grain = 3.62 ! called b in KZ, 3.62 for fractals, 4.53 for spheres
!     COS_PHI_ = cos(PHI_)

!     ! to calculate R0 (reflectance at 0 absorption by snow medium),
!     ! we need phase function p of scattering angle xi.
!     COS_XI_KZ = (SIN1*SIN2*COS_PHI_) - (COS1*COS2)
!     if (COS_XI_KZ .lt. -1.0) COS_XI_KZ = -1.0
!     if (COS_XI_KZ .gt.  1.0) COS_XI_KZ =  1.0
!     XI_KZ = acos(COS_XI_KZ)
!     XI_KZ_DEG = XI_KZ * 180.0 / PI_ ! convert to degrees to calculate phase function

!     ! calculate phase function p of scattering angle xi
!     P_pha = 11.1 * exp(-0.087 * XI_KZ_DEG) + 1.1 * exp(-0.014 * XI_KZ_DEG)

!     ! calculate R0 (reflectance at 0 absorption by snow medium)
!     R0 = (a + (b * (MU1 + MU2)) + (c * MU1 * MU2) + P_pha) / (4 * (MU1 + MU2))

!     ! calculate f (escape function)
!     ! note that this is a parameterization of the escape function
!     ! but one can run RT simulation to get the full escape function
!     K0_sza  = (3/7) * (1 + (2 * MU1))
!     K0_vza  = (3/7) * (1 + (2 * MU2))
!     f = K0_sza * K0_vza / R0


!     REFLECT = F_ISO + F_VOL*K_VOL + F_GEO*K_GEO + F_KZ*K_KZ

!     if (REFLECT .lt. 0.0) REFLECT = 0.0 ! relfectance cannot be negative

!     return

! end subroutine shd_sfc__brdfLSRTKZ

subroutine shd_sfc__brdfLSRTJ(F_ISO, F_GEO, F_VOL, F_J, MU1, MU2, PHI, ALPHA, REFLECT)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Computes the RossLi-Jiao BRDF (combines LSRT model above with a fourth kernel for snow from Jiao et al.).
    ! MU1: cosine of zenith angle of the incident light (solar zenith angle)
    ! MU2: cosine of zenith angle of the receiver (viewing zenith angle)
    ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
    ! ALPHA: parameter for the empirical function to control the forward scattering. Jiao et al. (2019) suggest ALPHA=0.2-0.3.
    !
    ! When MU1>0 and the hot spot is MU2=MU1 and PHI=0
    !
    ! Code reference:
    !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
    !
    ! Formula reference:
    !    Jiao et al. (2019) https://doi.org/10.1016/j.rse.2018.11.001
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(4), intent(in)  :: F_ISO, F_GEO, F_VOL, F_J, MU1, MU2, PHI, ALPHA
    real(4), intent(out) :: REFLECT
    real(4) :: COS1, COS2, SIN1, SIN2, TAN1, TAN2, COS_PHI, SIN_PHI
    real(4) :: K_GEO, K_VOL, K_J
    real(4) :: XI, COS_XI, SIN_XI
    real(4) :: B2R, H2B
    real(4) :: THETA1_, COS1_, SIN1_, TAN1_, SEC1_
    real(4) :: THETA2_, COS2_, SIN2_, TAN2_, SEC2_
    real(4) :: COS_XI_, D_2, COS_T, SIN_T, T, O
    real(4) :: COS_MIN

    ! calculate trigonometry
    !╭────────────────────────────────────────────────╮!
    COS_MIN = 0.03

    COS1 = MU1
    COS2 = MU2

    if (COS1 .lt. COS_MIN) COS1 = COS_MIN
    if (COS2 .lt. COS_MIN) COS2 = COS_MIN

    SIN1 = sqrt(1.0 - COS1**2)
    SIN2 = sqrt(1.0 - COS2**2)
    TAN1 = SIN1/COS1
    TAN2 = SIN2/COS2
    COS_PHI = cos(PHI)
    SIN_PHI = sin(PHI)
    !╰────────────────────────────────────────────────╯!

    ! K_VOL
    !╭────────────────────────────────────────────────╮!
    COS_XI = (COS1*COS2) + (SIN1*SIN2*COS_PHI)
    if (COS_XI .lt. -1.0) COS_XI = -1.0
    if (COS_XI .gt.  1.0) COS_XI =  1.0

    SIN_XI = sqrt(1.0 - COS_XI**2)
    XI = acos(COS_XI)

    K_VOL = ((PI_/2.0 - XI)*COS_XI + SIN_XI)/(COS1+COS2) - (PI_/4.0)
    !╰────────────────────────────────────────────────╯!

    ! K_GEO
    !╭────────────────────────────────────────────────╮!
    B2R = 1.0
    H2B = 2.0

    THETA1_ = atan(B2R*TAN2)
    COS1_ = cos(THETA1_)
    if (COS1_ .lt. -1.0) COS1_ = -1.0
    if (COS1_ .gt.  1.0) COS1_ =  1.0
    SIN1_ = sqrt(1.0 - COS1_**2)
    TAN1_ = SIN1_/COS1_
    SEC1_ = 1.0/COS1_

    THETA2_ = atan(B2R*TAN1)
    COS2_ = cos(THETA2_)
    if (COS2_ .lt. -1.0) COS2_ = -1.0
    if (COS2_ .gt.  1.0) COS2_ =  1.0
    SIN2_ = sqrt(1.0 - COS2_**2)
    TAN2_ = SIN2_/COS2_
    SEC2_ = 1.0/COS2_

    COS_XI_ = (COS1_*COS2_) + (SIN1_*SIN2_*COS_PHI)

    D_2 = TAN1_**2 + TAN2_**2 - 2.0*TAN1_*TAN2_*COS_PHI
    if (D_2 .lt. 0.0) D_2 = 0.0

    COS_T = H2B * (sqrt(D_2 + (TAN1_*TAN2_*SIN_PHI)**2))/(SEC1_ + SEC2_)
    if (COS_T .lt. -1.0) COS_T = -1.0
    if (COS_T .gt.  1.0) COS_T =  1.0
    SIN_T = sqrt(1.0 - COS_T**2)
    T = acos(COS_T)

    O = (1.0/PI_) * (T - SIN_T*COS_T) * (SEC1_ + SEC2_)

    K_GEO = O - SEC1_ - SEC2_ + 0.5d0*(1.0+COS_XI_)*SEC1_*SEC2_
    !╰────────────────────────────────────────────────╯!

    ! get Jiao kernel value K_J with normalization
    !╭────────────────────────────────────────────────╮!
    call shd_sfc__brdfJiao(MU1, MU2, PHI, ALPHA, 1, K_J)
    !╰────────────────────────────────────────────────╯!


    REFLECT = F_ISO + F_VOL*K_VOL + F_GEO*K_GEO + F_J*K_J

    if (REFLECT .lt. 0.0) REFLECT = 0.0 ! relfectance cannot be negative

    return

end subroutine shd_sfc__brdfLSRTJ



subroutine shd_sfc__brdfJiao(MU1, MU2, PHI, ALPHA, NORMALIZE, K_J)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Computes the Jiao BRDF for snow.
    ! MU1: cosine of zenith angle of the incident light (solar zenith angle)
    ! MU2: cosine of zenith angle of the receiver (viewing zenith angle)
    ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
    ! ALPHA: parameter for the empirical function to control the forward scattering. Jiao et al. (2019) suggest ALPHA=0.2-0.3.
    ! NORMALIZE: whether to add a normalizer term to the reflectance (1 for yes, 0 for no). If combining this kernel with an isotropic kernel (such as the case with RTLSR/LSRT), set NORMALIZE=1.
    !
    ! When MU1>0 and the hot spot is MU2=MU1 and PHI=0
    !
    ! Code reference:
    !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
    !
    ! Formula reference:
    !    Jiao et al. (2019) https://doi.org/10.1016/j.rse.2018.11.001
    !╰────────────────────────────────────────────────────────────────────────────╯!
    real(4), intent(in)  :: MU1, MU2, PHI, ALPHA
    integer, intent(in)  :: NORMALIZE
    real(4), intent(out) :: K_J
    real(4) :: COS1, COS2, SIN1, SIN2, TAN1, TAN2, COS_PHI, SIN_PHI
    real(4) :: R0, F_G, K0
    real(4) :: COS_MIN
    real(4) :: COS_XI_J, XI_J, XI_J_DEG, P_PHA
    real(4) :: K1, K2, K3 ! called a, b, c in Kokhanovsky-Zege

    ! set constants for parameterized escape function
    K1 = 1.247
    K2 = 1.186
    K3 = 5.157

    ! calculate trigonometry
    COS_MIN = 0.03

    COS1 = MU1
    COS2 = MU2

    if (COS1 .lt. COS_MIN) COS1 = COS_MIN
    if (COS2 .lt. COS_MIN) COS2 = COS_MIN

    SIN1 = sqrt(1.0 - COS1**2)
    SIN2 = sqrt(1.0 - COS2**2)
    COS_PHI = cos(PHI)

    ! to calculate R0 (reflectance at 0 absorption by snow medium),
    ! we need phase function p of scattering angle xi.
    COS_XI_J = (SIN1*SIN2*COS_PHI) + (COS1*COS2)
    if (COS_XI_J .lt. -1.0) COS_XI_J = -1.0
    if (COS_XI_J .gt.  1.0) COS_XI_J =  1.0
    XI_J = acos(COS_XI_J)
    XI_J_DEG = XI_J * 180.0 / PI_ ! convert to degrees to calculate phase function

    ! calculate phase function p of scattering angle xi
    P_PHA = 11.1 * exp(-0.087 * (180.0 - XI_J_DEG)) + 1.1 * exp(-0.014 * (180.0 - XI_J_DEG))

    ! calculate R0 (reflectance at 0 absorption by snow medium)
    R0 = (K1 + (K2 * (COS1 + COS2)) + (K3 * COS1 * COS2) + P_PHA) / (4 * (COS1 + COS2))

    ! calculate f_g (approximated, empirical function)
    F_G = COS_XI_J * exp(-COS_XI_J)

    ! normalizer term
    if (NORMALIZE .eq. 1) then
        K0 = (0.4076 * ALPHA) - 1.1081
    else
        K0 = 0.0
    end if

    K_J = (R0 * (1.0 - (ALPHA * F_G))) + K0

    return

end subroutine shd_sfc__brdfJiao


subroutine shd_sfc__brdfOcean(&
    & pws, xsal, pcl, pwl, xmuo, xmu, xphi, xpaw, REFLECT)
    !
    ! INPUT:  pws=wind speed (in m/s) [value not allowed below 0.25 m/s]
    !         xsal=salinity (in ppt) [if xsal<0 then 34.3ppt by default]
    !         pcl=pigment concentration (in mg/m^3)
    !         pwl=wavelength of the computation (in micrometer)
    !         xmuo=cosine of the zenith angle of the incident radiation
    !         xmu=cosine of the zenith angle of the reflected radiation
    !         xphi=relative azimuth angle (between 0 and 2*pi)
    !         xpaw=azim. of sun - azim. of wind (in radian)
    !              For the direct solar radiation, the wind direction
    !              is fixed at the azimuth angle of 0 degree so that
    !              the azimuth angle of the sun = xpaw.
    !
    ! OUTPUT: REFLECT=the total reflectance of the sea water

    real, intent(in)  :: pws, xsal, pcl, pwl, xmuo, xmu, xphi, xpaw
    real, intent(out) :: REFLECT

    real Ref(39)
    real :: paw, wl, wspd, C, azw, phi
    real fac, nr, ni, n12
    real tetas, w, wlp, ref_i, rwc, rw, tds
    real tetav, tdv, fi, rog, a, rwb
    real angbnd(5), wsbnd(6), tdsbnd(5, 6), tdvbnd(5, 6)

    integer :: iwl, iang, iws, isz1, isz2, ivz1, ivz2, iws1, iws2

    ! effective reflectance of the whitecaps (Koepke, 1984)
    data Ref/&
        & 0.220, 0.220, 0.220, 0.220, 0.220, 0.220, 0.215, 0.210, 0.200, 0.190,&
        & 0.175, 0.155, 0.130, 0.080, 0.100, 0.105, 0.100, 0.080, 0.045, 0.055,&
        & 0.065, 0.060, 0.055, 0.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
        & 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/

    data angbnd/0.0, 45.0, 60.0, 75.0, 85.0/
    data wsbnd/1.0, 3.0, 5.0, 7.0, 9.0, 20.0/

    data((tdsbnd(iang, iws), iws=1, 6), iang=1, 5)/&
        & 0.9787803, 0.9787738, 0.9787626, 0.9787467, 0.9787264,&
        & 0.9785573, 0.9706900, 0.9698871, 0.9691746, 0.9685547,&
        & 0.9680276, 0.9666586, 0.9479931, 0.9404608, 0.9385692,&
        & 0.9381815, 0.9384519, 0.9430056, 0.9690591, 0.9275920,&
        & 0.9058769, 0.8951812, 0.8899654, 0.8892645, 0.9980542,&
        & 0.9602273, 0.9114283, 0.8713799, 0.8417820, 0.7800314/

    data((tdvbnd(iang, iws), iws=1, 6), iang=1, 5)/&
        & 0.9787764, 0.9787535, 0.9787106, 0.9786453, 0.9785548,&
        & 0.9775019, 0.9692680, 0.9637051, 0.9564344, 0.9495727,&
        & 0.9438773, 0.9288712, 0.9225163, 0.9069787, 0.9044844,&
        & 0.9052351, 0.9068328, 0.9153687, 0.8048478, 0.8479503,&
        & 0.8678726, 0.8797889, 0.8878716, 0.9091171, 0.7294627,&
        & 0.8137348, 0.8453338, 0.8629867, 0.8745421, 0.9036854/

    fac = PI_/180.
    ! paw= azim. of sun - azim. of wind (in deg.)
    paw = xpaw/fac
    ! Conversion of parameter
    if (xphi .lt. 0.0) then
        phi = -xphi
    else if (xphi .ge. 2.0*PI_) then
        phi = xphi - 2.0*PI_
    else
        phi = xphi
    end if

    if (xmuo .le. 0.028) then
        tetas = acos(0.028)/fac
    else
        tetas = acos(xmuo)/fac
    end if

    if (xmu .le. 0.028) then
        tetav = acos(0.028)/fac
    else
        tetav = acos(xmu)/fac
    end if

    if (pwl .lt. 0.4) then
        wl = 0.4
    else if (pwl .gt. 4.0) then
        wl = 4.0
    else
        wl = pwl
    end if

    fi = 180.0 - phi/fac

    C = pcl
    wspd = max(0.25, pws)
    azw = paw

    ! COMPUTE INDEX OF WATER
    call indwat(wl, xsal, nr, ni)
    n12 = sqrt(nr*nr + ni*ni)
    ! COMPUTE WHITECAPS REFLECTANCE (LAMBERTIAN)
    W = 2.95d-6*(wspd**3.52d0)
    iwl = 1 + int((wl - 0.2)/0.1)
    wlp = 0.5d0 + (iwl - 1)*0.1d0
    !Ref_i = ref(iwl + 1) + (wl - wlp)/0.1*(ref(iwl) - ref(iwl + 1))
    Ref_i = ref(iwl + 1) + (wl - wlp)*10.0*(ref(iwl) - ref(iwl + 1))
    Rwc = W*Ref_i

    ! COMPUTE BACKSCATTERED REFLECTANCE FROM THE SEA WATER (LAMBERTIAN)
    !  water reflectance below the sea surface
    call morcasiwat(wl, C, Rw)

    call getbound( wsbnd, 1, 6,  wspd, iws1, iws2)
    call getbound(angbnd, 1, 5, tetas, isz1, isz2)
    call getbound(angbnd, 1, 5, tetav, ivz1, ivz2)
    tds = tdsbnd(isz1, iws1)
    tdv = tdvbnd(ivz1, iws1)

    ! SUNGLINT REFLECTANCE
    call sunglint(wspd, nr, ni, azw, tetas, tetav, fi, rog)
    ! water reflectance above the sea surface
    ! for explanation on value of a see OCEAALBE.f
    a = 0.485d0
    ! add change in solid angle from under to above to surface
    ! that account for 1/(n12*n12) decrease in sea water directional
    ! reflectance
    Rwb = (1.0d0/(n12*n12))*tds*tdv*Rw/(1.0d0 - a*Rw)
    ! TOTAL REFLECTANCE OF SEA WATER
    REFLECT = Rwc + (1.0d0 - W)*Rog + (1.0d0 - Rwc)*Rwb
    return
end subroutine shd_sfc__brdfOcean


subroutine morcasiwat(wl, C, R2)
    ! Spectral diffuse attenuation coefficient of Case I Waters as Predicted
    ! by MOREL within the spectral range 400-700nm (1988, Journal of Geophys
    ! Research, Vol.93, No C9, pp 10749-10768)
    !
    ! input parameters: wl wavelength (IN MICROMETERS)
    !                   C  pigment concentration
    ! output parameter: R2  reflectance of water
    !
    ! According Morel,1988, we use:
    !
    ! Kd    spectral value of the attenuation coefficient for
    !        downwelling irradiance
    !        with: Kd=Kw+Xc*C**e
    ! Kw    spectral value of the diffuse attenuation coefficient
    !        for pure oceanic water
    ! Xc, e spectral coefficients to compute the diffuse attenuation
    !       coefficient for pigment
    ! bb  total backscattering coefficient
    !  with: bb=0.5*bw+bbt*b
    ! bw    spectral value of the molecular scattering coefficient of water
    ! bbt,b parameters to compute the scattering coefficients of pigments
    !
    ! R2    reflectance of water below the surface
    !    with: R2=(0.33/u)*(bb/Kd)  where u is depending of R2
    !
    real Kw, Kd
    real tKw(61), tXc(61), te(61), tbw(61)
    real wl, c, r2, xc, e, bw, bb, b, bbt, u1, r1, u2, err
    integer iwl

    data tKw/0.0209, 0.0200, 0.0196, 0.0189, 0.0183,&
        & 0.0182, 0.0171, 0.0170, 0.0168, 0.0166,&
        & 0.0168, 0.0170, 0.0173, 0.0174, 0.0175,&
        & 0.0184, 0.0194, 0.0203, 0.0217, 0.0240,&
        & 0.0271, 0.0320, 0.0384, 0.0445, 0.0490,&
        & 0.0505, 0.0518, 0.0543, 0.0568, 0.0615,&
        & 0.0640, 0.0640, 0.0717, 0.0762, 0.0807,&
        & 0.0940, 0.1070, 0.1280, 0.1570, 0.2000,&
        & 0.2530, 0.2790, 0.2960, 0.3030, 0.3100,&
        & 0.3150, 0.3200, 0.3250, 0.3300, 0.3400,&
        & 0.3500, 0.3700, 0.4050, 0.4180, 0.4300,&
        & 0.4400, 0.4500, 0.4700, 0.5000, 0.5500,&
        & 0.6500/
    data tXc/0.1100, 0.1110, 0.1125, 0.1135, 0.1126,&
        & 0.1104, 0.1078, 0.1065, 0.1041, 0.0996,&
        & 0.0971, 0.0939, 0.0896, 0.0859, 0.0823,&
        & 0.0788, 0.0746, 0.0726, 0.0690, 0.0660,&
        & 0.0636, 0.0600, 0.0578, 0.0540, 0.0498,&
        & 0.0475, 0.0467, 0.0450, 0.0440, 0.0426,&
        & 0.0410, 0.0400, 0.0390, 0.0375, 0.0360,&
        & 0.0340, 0.0330, 0.0328, 0.0325, 0.0330,&
        & 0.0340, 0.0350, 0.0360, 0.0375, 0.0385,&
        & 0.0400, 0.0420, 0.0430, 0.0440, 0.0445,&
        & 0.0450, 0.0460, 0.0475, 0.0490, 0.0515,&
        & 0.0520, 0.0505, 0.0440, 0.0390, 0.0340,&
        & 0.0300/
    data te/0.668, 0.672, 0.680, 0.687, 0.693,&
        & 0.701, 0.707, 0.708, 0.707, 0.704,&
        & 0.701, 0.699, 0.700, 0.703, 0.703,&
        & 0.703, 0.703, 0.704, 0.702, 0.700,&
        & 0.700, 0.695, 0.690, 0.685, 0.680,&
        & 0.675, 0.670, 0.665, 0.660, 0.655,&
        & 0.650, 0.645, 0.640, 0.630, 0.623,&
        & 0.615, 0.610, 0.614, 0.618, 0.622,&
        & 0.626, 0.630, 0.634, 0.638, 0.642,&
        & 0.647, 0.653, 0.658, 0.663, 0.667,&
        & 0.672, 0.677, 0.682, 0.687, 0.695,&
        & 0.697, 0.693, 0.665, 0.640, 0.620,&
        & 0.600/
    data tbw/0.0076, 0.0072, 0.0068, 0.0064, 0.0061,&
        & 0.0058, 0.0055, 0.0052, 0.0049, 0.0047,&
        & 0.0045, 0.0043, 0.0041, 0.0039, 0.0037,&
        & 0.0036, 0.0034, 0.0033, 0.0031, 0.0030,&
        & 0.0029, 0.0027, 0.0026, 0.0025, 0.0024,&
        & 0.0023, 0.0022, 0.0022, 0.0021, 0.0020,&
        & 0.0019, 0.0018, 0.0018, 0.0017, 0.0017,&
        & 0.0016, 0.0016, 0.0015, 0.0015, 0.0014,&
        & 0.0014, 0.0013, 0.0013, 0.0012, 0.0012,&
        & 0.0011, 0.0011, 0.0010, 0.0010, 0.0010,&
        & 0.0010, 0.0009, 0.0008, 0.0008, 0.0008,&
        & 0.0007, 0.0007, 0.0007, 0.0007, 0.0007,&
        & 0.0007/
    if (wl .lt. 0.400 .or. wl .gt. 0.700) then
        R2 = 0.000
        return
    end if

    iwl = 1 + nint((wl - 0.400)/0.005)
    Kw = tKw(iwl)
    Xc = tXc(iwl)
    e = te(iwl)
    bw = tbw(iwl)

    if (abs(C) .lt. 0.0001) then
        bb = 0.5d0*bw
        Kd = Kw
    else
        b = 0.30d0*C**0.62d0
        bbt = 0.002d0 + 0.02d0*(0.5d0 - 0.25d0*alog10(C))*0.550d0/wl
        bb = 0.5d0*bw + bbt*b
        Kd = Kw + Xc*C**e
    end if

    u1 = 0.75
    R1 = 0.33d0*bb/u1/Kd

50  u2 = 0.90d0*(1.0-R1)/(1.0+2.25*R1)
    R2 = 0.33d0*bb/u2/Kd
    err = abs((R2 - R1)/R2)
    if (err .lt. 0.0001) return
    R1 = R2
    goto 50
    return
end subroutine

subroutine indwat(wl, xsal, nr, ni)
    !
    ! input parameters:  wl=wavelength (in micrometers)
    !                    xsal=salinity (in ppt), if xsal<0 then 34.3ppt by d
    ! output parameters: nr=index of refraction of sea water
    !                    ni=extinction coefficient of sea water
    !
    real twl(62), tnr(62), tni(62)
    real nr, ni, wl, xwl, yr, yi, nrc, nic, xsal
    integer i
    ! Indices of refraction for pure water from Hale and Querry,
    ! Applied Optique, March 1973, Vol. 12,  No. 3, pp. 555-563
    data twl/&
        & 0.250, 0.275, 0.300, 0.325, 0.345, 0.375, 0.400, 0.425, 0.445, 0.475,&
        & 0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725,&
        & 0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975,&
        & 1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.650,&
        & 2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000, 3.050, 3.100, 3.150,&
        & 3.200, 3.250, 3.300, 3.350, 3.400, 3.450, 3.500, 3.600, 3.700, 3.800,&
        & 3.900, 4.000/
    data tnr/&
        & 1.362, 1.354, 1.349, 1.346, 1.343, 1.341, 1.339, 1.338, 1.337, 1.336,&
        & 1.335, 1.334, 1.333, 1.333, 1.332, 1.332, 1.331, 1.331, 1.331, 1.330,&
        & 1.330, 1.330, 1.329, 1.329, 1.329, 1.328, 1.328, 1.328, 1.327, 1.327,&
        & 1.327, 1.324, 1.321, 1.317, 1.312, 1.306, 1.296, 1.279, 1.242, 1.219,&
        & 1.188, 1.157, 1.142, 1.149, 1.201, 1.292, 1.371, 1.426, 1.467, 1.483,&
        & 1.478, 1.467, 1.450, 1.432, 1.420, 1.410, 1.400, 1.385, 1.374, 1.364,&
        & 1.357, 1.351/
    data tni/&
        & 3.35e-08, 2.35e-08, 1.60e-08, 1.08e-08, 6.50e-09,&
        & 3.50e-09, 1.86e-09, 1.30e-09, 1.02e-09, 9.35e-10,&
        & 1.00e-09, 1.32e-09, 1.96e-09, 3.60e-09, 1.09e-08,&
        & 1.39e-08, 1.64e-08, 2.23e-08, 3.35e-08, 9.15e-08,&
        & 1.56e-07, 1.48e-07, 1.25e-07, 1.82e-07, 2.93e-07,&
        & 3.91e-07, 4.86e-07, 1.06e-06, 2.93e-06, 3.48e-06,&
        & 2.89e-06, 9.89e-06, 1.38e-04, 8.55e-05, 1.15e-04,&
        & 1.10e-03, 2.89e-04, 9.56e-04, 3.17e-03, 6.70e-03,&
        & 1.90e-02, 5.90e-02, 1.15e-01, 1.85e-01, 2.68e-01,&
        & 2.98e-01, 2.72e-01, 2.40e-01, 1.92e-01, 1.35e-01,&
        & 9.24e-02, 6.10e-02, 3.68e-02, 2.61e-02, 1.95e-02,&
        & 1.32e-02, 9.40e-03, 5.15e-03, 3.60e-03, 3.40e-03,&
        & 3.80e-03, 4.60e-03/
    i = 2
10  if (wl .lt. twl(i)) goto 20
    if (i .lt. 62) then
        i = i + 1
        goto 10
    end if
20  xwl = twl(i) - twl(i - 1)
    yr = tnr(i) - tnr(i - 1)
    yi = tni(i) - tni(i - 1)
    nr = tnr(i - 1) + (wl - twl(i - 1))*yr/xwl
    ni = tni(i - 1) + (wl - twl(i - 1))*yi/xwl
    !
    ! Correction to be applied to the index of refraction and to the extinct
    ! coefficients of the pure water to obtain the ocean water one (see for
    ! example Friedman). By default, a typical sea water is assumed
    ! (Salinity=34.3ppt, Chlorinity=19ppt) as reported by Sverdrup.
    ! In that case there is no correction for the extinction coefficient bet
    ! 0.25 and 4 microns. For the index of refraction, a correction of +0.00
    ! has to be applied (McLellan). For a chlorinity of 19.0ppt the correcti
    ! is a linear function of the salt concentration. Then, in 6S users are
    ! to enter the salt concentration (in ppt).
    ! REFERENCES:
    ! Friedman D., Applied Optics, 1969, Vol.8, No.10, pp.2073-2078.
    ! McLellan H.J., Elements of physical Oceanography, Pergamon Press, Inc.
    !        New-York, 1965, p 129.
    ! Sverdrup H.V. et al., The Oceans (Prentice-Hall, Inc., Englewood Cliff
    !        N.J., 1942, p 173.

    nrc = 0.006
    nic = 0.000
    if (xsal .ge. 0.0) then
        nr = nr + nrc*(xsal/34.3)
        ni = ni + nic*(xsal/34.3)
    else
        nr = nr + nrc
        ni = ni + nic
    end if
    return
end subroutine

subroutine sunglint(wspd, nr, ni, azw, ts, tv, fi, rog)
    ! input parameters:   wspd=speed of the wind (in m/s)
    !                     nr=index of refraction of the sea water
    !                     ni=extinction coefficient of the sea water
    !                     azw=azim. of the sun - azim. of the wind (in deg.)
    !                     ts=solar zenith angle (in deg.)
    !                     tv=view zenith angle (in deg.)
    !                     fi=relative azimuth (sun-satellite)
    ! output parameters:  rog=reflectance of the sun glint
    !
    real fac
    real wspd, nr, ni, ts, tv, fi, rog, azw, phw
    real cs, cv, ss, sv, phi, zx, zy, tantilt, tilt, proba, xe, xn, xe2, xn2
    real coef, cos2chi, coschi, sinchi
    real r1, sigmaC, sigmaU, C21, C03, C40, C04, C22

    fac = PI_/180.0
    phw = azw*fac
    cs = cos(ts*fac)
    cv = cos(tv*fac)
    ss = sin(ts*fac)
    sv = sin(tv*fac)
    phi = fi*fac
    Zx = -sv*sin(phi)/(cs + cv)
    Zy = (ss + sv*cos(phi))/(cs + cv)
    tantilt = sqrt(zx*zx + zy*zy)
    tilt = atan(tantilt)
    !  Anisotropic Gaussian distribution
    !    phw=phi_sun-phi_wind
    sigmaC = 0.003d0 + 0.00192d0*wspd
    sigmaU = 0.00316d0*wspd
    C21 = 0.01d0 - 0.0086d0*wspd
    C03 = 0.04d0 - 0.033d0*wspd
    C40 = 0.40
    C22 = 0.12
    C04 = 0.23
    xe = ( cos(phw)*Zx + sin(phw)*Zy)/sqrt(SigmaC)
    xn = (-sin(phw)*Zx + cos(phw)*Zy)/sqrt(SigmaU)
    xe2 = xe*xe
    xn2 = xn*xn
    coef = 1.0 - C21/2.0*(xe2 - 1.0)*xn - C03/6.0*(xn2 - 3.0)*xn
    coef = coef + c40/24.0*(xe2*xe2 - 6.0*xe2 + 3.0)
    coef = coef + C04/24.0*(xn2*xn2 - 6.0*xn2 + 3.0)
    coef = coef + C22/4.0*(xe2 - 1)*(xn2 - 1)
    proba = coef/2.0/PI_/sqrt(sigmaU)/sqrt(sigmaC)*exp(-(xe2 + xn2)/2.0)
    ! Compute Fresnel's coefficient R1
    cos2chi = cv*cs + sv*ss*cos(phi)
    if (cos2chi .gt.  1.0) cos2chi = 0.99999999999
    if (cos2chi .lt. -1.0) cos2chi = -0.99999999999
    coschi = sqrt(0.5d0*(1 + cos2chi))
    sinchi = sqrt(0.5d0*(1 - cos2chi))
    call Fresnel(nr, ni, coschi, sinchi, R1)
    ! Compute Reflectance of the sun glint
    Rog = PI_*R1*proba/4.0/cs/cv/(cos(tilt)**4)
    return
end subroutine

subroutine Fresnel(nr, ni, coschi, sinchi, R1)
    !
    ! to compute the Fresnel's coefficient of reflection (see for
    ! example M. Born and E. Wolf, Principles of Optics, Pergamon Press, fif
    ! edition, 1975, pp 628
    ! input parameters: nr=index of refraction of the sea water
    !                   ni=extinction coefficient of the sea water
    !                   coschi & sinchi=cosine and sine of the incident radi
    !                                   with respect of the wave facet norma
    ! output parameter: R1=Fresnel's coefficient for reflection
    !
    real nr, ni, a1, a2, u, v, Rr2, Rl2, b1, b2, R1, coschi, sinchi
    ! absolute value for a1 to get v=0 when ni=0
    a1 = abs(nr*nr - ni*ni - sinchi*sinchi)
    a2 = sqrt((nr*nr - ni*ni - sinchi*sinchi)**2.+4*nr*nr*ni*ni)
    u = sqrt(0.5d0*(a1 + a2))
    v = sqrt(max(0.0, 0.5d0*(-a1 + a2)))
    Rr2 = ((coschi - u)**2 + v*v)/((coschi + u)**2 + v*v)
    b1 = (nr*nr - ni*ni)*coschi
    b2 = 2.0*nr*ni*coschi
    Rl2 = ((b1 - u)**2 + (b2 - v)**2)/((b1 + u)**2 + (b2 + v)**2)
    R1 = (Rr2 + Rl2)/2.0
    return
end subroutine

subroutine glitalbe(wspd, nr, ni, azw, rge)
    !
    ! To compute the spherical albedo of the sea water. See for example
    ! Masuda et al., Remote Sens. Environ., 24, 313-329, 1988.
    !
    ! input parameters: wsp=wind of speed
    !                   nr=index of refraction of the sea water
    !                   ni=extinction coefficient of the sea water
    !                   azw=azim. of sun - azim. of wind (in deg.)
    ! output parameter: rge=spherical albedo of the sun glint
    !
    real nr, ni, azw, phw, rge, q, wspd, prefl, proba, pr, pp, fac
    real sigma, sigmaC, sigmaU, C21, C03, C40, C04, C22
    real costt, hta, htb, hfa, cotb, cota, cofa, diff, coef
    real phin, cosphin, sinphin, costet, tet, sintet
    real costetn, sintetn, tantetn, coschi, sinchi
    real zx, zy, xe, xn, xe2, xn2, fonc0, pond, r1
    integer km, i, j
    integer :: nfa, nta, ntb

    fac = PI_/180.0
    sigma  = 0.003d0 + 0.00512d0*wspd
    sigmaC = 0.003d0 + 0.00192d0*wspd
    sigmaU = 0.00316d0*wspd
    C21 = 0.01d0 - 0.0086d0*wspd
    C03 = 0.04d0 -  0.033d0*wspd
    C40 = 0.40
    C22 = 0.12
    C04 = 0.23
    ! costt to minimize the time of the computation
    !     integration between 1 and costt instead of 1 and 0
    q = 50
    costt = 1.0/sqrt(1.0 + q*sigma/4.0)
    phw = azw*fac

    prefl = 0.0
    proba = 0.0

    ntb = 31
    htb = 1.0/float(ntb - 1)
    ! loops on the zenith angle of the emitted radiation
    do km = 1, ntb
        costet = (km - 1)*htb
        tet = acos(costet)
        sintet = sin(tet)
        tet = tet/fac
        ! Simpson's rules for the angle of the emitted radiation teta
        cotb = 2.0
        diff = abs(km/2 - km/2.0)
        if (diff .lt. 0.00001) cotb = 4.0
        if (km .eq. 1 .or. km .eq. ntb) cotb = 1.0
        !  loops step for phiN and tetaN (N is the facet unit normal vector)
        if (tet .lt. 91) nta = 801
        if (tet .lt. 81) nta = 301
        if (tet .lt. 75) nta = 101
        if (tet .lt. 65) nta = 31
        nfa = nta
        hta = (1.-costt)/float(nta - 1)
        hfa = PI_/float(nfa - 1)
        ! loops on phiN (azimuth angle of the facet normal vector)
        pr = 0.0
        pp = 0.0
        do i = 1, nfa
            phin = (i - 1)*hfa
            cosphin = cos(phin)
            sinphin = sin(phin)
            !  Simpson's rules for phin
            cofa = 2.0
            diff = abs(i/2 - i/2.0)
            if (diff .lt. 0.00001) cofa = 4.0
            if (i .eq. 1 .or. i .eq. nfa) cofa = 1.0
            ! loops on tetaN (zenith angle of the facet normal vector)
            do j = 1, nta
                costetn = costt + (j - 1)*hta
                sintetn = sqrt(1 - costetn*costetn)
                tantetn = sintetn/costetn
                !  Simpson's rules for tetaN
                cota = 2.0
                diff = abs(j/2 - j/2.0)
                if (diff .lt. 0.00001) cota = 4.0
                if (j .eq. 1 .or. j .eq. nta) cota = 1.0
                ! Fresnel's reflection coefficient R1
                coschi = costet*costetn + sintet*sintetn*cosphin
                if (coschi*coschi .gt. 1.0) coschi = 0.99999999999
                sinchi = sqrt(1 - coschi*coschi)
                if (coschi .lt. 0.0) then
                    r1 = 0.0
                    cota = 0.0
                else
                    call Fresnel(nr, ni, coschi, sinchi, r1)
                end if
                !  Anisotropic Gaussian distribution for wave facets slopes
                Zx = -tantetn*cosphin
                Zy = -tantetn*sinphin
                xe = (cos(phw)*Zx + sin(phw)*Zy)/sqrt(SigmaC)
                xn = (-sin(phw)*Zx + cos(phw)*Zy)/sqrt(SigmaU)
                xe2 = xe*xe
                xn2 = xn*xn
                coef = 1.0 - C21/2.0*(xe2 - 1.0)*xn - C03/6.0*(xn2 - 3.0)*xn
                coef = coef + c40/24.0*(xe2*xe2 - 6.0*xe2 + 3.0)
                coef = coef + C04/24.0*(xn2*xn2 - 6.0*xn2 + 3.0)
                coef = coef + C22/4.0*(xe2 - 1.0)*(xn2 - 1.0)
                fonc0 = 0.5d0*coschi*coef*exp(-(xe2 + xn2)/2.0)/(costetn**4)
                pr = pr + r1*fonc0*cofa*cota*cotb
                pp = pp + fonc0*cofa*cota*cotb
            end do
        end do

        !pond = 2.0*hta*hfa*htb/PI_/sqrt(sigmaC)/sqrt(sigmaU)/3./3./3.
        pond = 2.0*hta*hfa*htb/PI_/sqrt(sigmaC)/sqrt(sigmaU)/27.0
        prefl = prefl + pr*pond
        proba = proba + pp*pond
    end do
    rge = prefl/proba
    return
end subroutine

subroutine getbound(xvals, ifirst, ilast, x, ind1, ind2)
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    !
    !     PROGRAM: GETBOUND.FOR
    !
    !     DESCRIPTION: FINDS THE BOUNDS XVALS(IND1) AND XVALS(IND2)
    !                  FOR WHICH X IS BOUNDED.
    !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
    real xvals(*), x
    integer ind1, ind2, ifirst, ilast, i, imid

    imid = ilast/2 + 1

    if (xvals(ilast) .gt. xvals(ifirst)) then
        if (x .gt. xvals(imid)) then
            do i = imid, ilast - 1
                if (xvals(i) .le. x .and. xvals(i + 1) .ge. x) then
                    ind1 = i
                    ind2 = i + 1
                    return
                end if
            end do
        else
            do i = ifirst, imid
                if (xvals(i) .le. x .and. xvals(i + 1) .ge. x) then
                    ind1 = i
                    ind2 = i + 1
                    return
                end if
            end do
        end if
        if (x .lt. xvals(ifirst)) then
            ind1 = ifirst
            ind2 = ifirst + 1
        else
            ind1 = ilast - 1
            ind2 = ilast
        end if

    else

        if (x .lt. xvals(imid)) then
            do i = imid, ilast - 1
                if (xvals(i) .ge. x .and. xvals(i + 1) .le. x) then
                    ind1 = i
                    ind2 = i + 1
                    return
                end if
            end do
        else
            do i = ifirst, imid
                if (xvals(i) .ge. x .and. xvals(i + 1) .le. x) then
                    ind1 = i
                    ind2 = i + 1
                    return
                end if
            end do
        end if
        if (x .gt. xvals(ifirst)) then
            ind1 = ifirst
            ind2 = ifirst + 1
        else
            ind1 = ilast - 1
            ind2 = ilast
        end if
    end if
    return
end subroutine

end module shd_sfc
