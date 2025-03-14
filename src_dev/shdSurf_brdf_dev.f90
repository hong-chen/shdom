module shdSurf

  public :: shdSurf__surface_brdf

contains

  subroutine shdSurf__surface_brdf(SFCTYPE, REFPARMS, WAVELEN,&
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
      implicit none
      integer :: NSTOKES
      integer :: K1, K2
      real    :: REFPARMS(*), WAVELEN, MU1, PHI1, MU2, PHI2, REFLECT(4, 4)
      real    :: PI
      real    :: RPV_REFLECTION
      character :: SFCTYPE*1

      PI = acos(-1.0)

      if (SFCTYPE .eq. 'L' .or. SFCTYPE .eq. 'l') then

          ! L or l: Lambertian surface BRDF is constant.
          !   (for testing, as done more efficiently by Lambertian routine
          REFLECT(1:NSTOKES, 1:NSTOKES) = 0.0
          REFLECT(1, 1) = REFPARMS(1)

      else if (SFCTYPE .eq. 'W' .or. SFCTYPE .eq. 'w') then

          ! W: Fresnel reflection with Gaussian distribution of slopes
          call WAVE_FRESNEL_REFLECTION(REFPARMS(1), REFPARMS(2),&
              &REFPARMS(3), MU1, MU2, PHI1, PHI2, NSTOKES, REFLECT)

      else if (SFCTYPE .eq. 'D' .or. SFCTYPE .eq. 'd') then

          ! D: Diner et al., depolarizing modified RPV + Fresnel reflection
          !   from microfacets with uniform or Gaussian distribution of s
          call DINER_REFLECTION(REFPARMS(1), REFPARMS(2), REFPARMS(3),&
              &REFPARMS(4), REFPARMS(5),&
              &-MU1, MU2, PHI2 - PHI1, NSTOKES, REFLECT)

      else if (SFCTYPE .eq. 'R' .or. SFCTYPE .eq. 'r') then

          ! R: Rahman, Pinty, and Verstraete
          if (NSTOKES .gt. 1) then
              print *, 'RPV BRDF is only for unpolarized case'
              stop
          end if
          REFLECT(1, 1) = RPV_REFLECTION(REFPARMS(1), REFPARMS(2),&
              &REFPARMS(3), -MU1, MU2, PHI1 - PHI2 - PI)

      else if (SFCTYPE .eq. 'O' .or. SFCTYPE .eq. 'o') then

          ! O: Ocean BRDF from 6S modified by Norm Loeb
          if (NSTOKES .gt. 1) then
              print *, 'Ocean BRDF from 6S is only for unpolarized case'
              stop
          end if
          call ocean_brdf_sw(REFPARMS(1), -1., REFPARMS(2), WAVELEN,&
              & -MU1, MU2, PHI1 - PHI2, PHI1, REFLECT(1, 1))

      else if (SFCTYPE .eq. 'T' .or. SFCTYPE .eq. 't') then

          ! T: LSRT model
          REFLECT(1, 1) = LSRT_REFLECTION(REFPARMS(1), REFPARMS(2),&
              &REFPARMS(3), -MU1, MU2, PHI1 - PHI2 - PI)

      else
          stop 'SURFACE_BRDF: Unknown BRDF type'
      end if

      return
  end subroutine

  subroutine WAVE_FRESNEL_REFLECTION(MRE, MIM, WINDSPEED,&
      & MUI, MUR, PHII, PHIR, NSTOKES, REFLECT)
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

      implicit none
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
      SHADOWR = 0.5*(S1*T1/DCOT - T2)
      SHADOW = 1d0/(1d0 + SHADOWI + SHADOWR)
      REFLECT(1:NSTOKES, 1:NSTOKES) = REFLECT(1:NSTOKES, 1:NSTOKES)*SHADOW
      return
  end subroutine

  subroutine DINER_REFLECTION(A, K, B, ZETA, SIGMA,&
      & MU1, MU2, PHI, NSTOKES, REFLECT)
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

      implicit none
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
      GAMMA = 0.5*acos(-COSSCATANG)
      COSGAMMA = cos(GAMMA)

      ! Calculate the Fresnel reflection coefficients (equivalent to eq 1
      !   [D = n*cos(gamma')]
      EPSILON = SFCINDEX**2
      D = sqrt(EPSILON - 1.0 + COSGAMMA**2)
      RP = (EPSILON*COSGAMMA - D)/(EPSILON*COSGAMMA + D)
      RS = (COSGAMMA - D)/(COSGAMMA + D)

      ! Calculate the Mueller matrix elements for Fresnel reflection (eq
      F11 = 0.5*(CABS(RP)**2 + CABS(RS)**2)
      F12 = 0.5*(CABS(RP)**2 - CABS(RS)**2)
      F33 = real(RP*conjg(RS))
      F34 = -IMAG(RP*conjg(RS))

      ! Calculate facet tilt angle, beta (eq 16b)
      COSBETA = 0.5*(MU1 + MU2)/COSGAMMA

      ! Calculate the factor before the M*P*M matrix including the distri
      !   of microfacet orientations (Gaussian for SIGMA>0) (eq 13 & 14)
      !   Note: Diner's BRDF is multiplied by pi here
      if (SIGMA .gt. 0.0) then
          H = ZETA*exp(-0.5*(1/COSBETA**2 - 1)/SIGMA**2)&
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
  end subroutine

  real function RPV_REFLECTION(RHO0, K, THETA, MU1, MU2, PHI)
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

      real :: RHO0, K, THETA, MU1, MU2, PHI, mu_min, x1, x2
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

      RPV_REFLECTION = RHO0*M*F*H

      return
  end function

  real function LSRT_REFLECTION(F_ISO, F_GEO, F_VOL, MU1, MU2, PHI)
      !╭────────────────────────────────────────────────────────────────────────────╮!
      !   Computes the RossLi BRDF.  The incident
      ! and outgoing cosine zenith angles are MU1 and MU2, respectively,
      ! and the relative azimuthal angle is PHI.  In this case the incident
      ! direction is where the radiation is coming from (i.e. opposite of
      ! discrete ordinate), so MU1>0 and the hot spot is MU2=MU1 and PHI=0
      !╰────────────────────────────────────────────────────────────────────────────╯!

      real(8) :: PI, COS1, COS2, SIN1, SIN2, TAN1, TAN2, COS_PHI, SIN_PHI
      real(8) :: F_ISO, F_GEO, F_VOL, MU1, MU2, PHI
      real(8) :: K_GEO, K_VOL
      real(8) :: XI, COS_XI, SIN_XI
      real(8) :: B2R, H2B
      real(8) :: THETA1_, COS1_, SIN1_, TAN1_, SEC1_
      real(8) :: THETA2_, COS2_, SIN2_, TAN2_, SEC2_
      real(8) :: COS_XI_, D, COS_T, SIN_T, T, O

      ! calculate trigonometry
      !╭────────────────────────────────────────────────╮!
      PI = acos(-1.0_R)
      COS1 = MU1
      COS2 = MU2
      SIN1 = sqrt(1.0_R - COS1**2)
      SIN2 = sqrt(1.0_R - COS2**2)
      TAN1 = SIN1/COS1
      TAN2 = SIN2/COS2
      COS_PHI = cos(PHI)
      SIN_PHI = sin(PHI)
      !╰────────────────────────────────────────────────╯!

      ! K_VOL
      !╭────────────────────────────────────────────────╮!
      COS_XI = ((COS1*COS2)) + (SIN1*SIN2*COS_PHI)
      SIN_XI = sqrt(1.0_R - COS_XI**2)
      XI = acos(COS_XI)
      K_VOL = ((PI/2.0_R-XI)*COS_XI + SIN_XI)/(COS1+COS2) - (PI/4.0_R)
      !╰────────────────────────────────────────────────╯!

      ! K_GEO
      !╭────────────────────────────────────────────────╮!
      B2R = 1.0_R
      H2B = 2.0_R
      THETA1_ = atan(B2R*TAN2)
      THETA2_ = THETA1_ / (atan(B2R*TAN1))
      COS1_ = cos(THETA1_)
      COS2_ = cos(THETA2_)
      SIN1_ = sqrt(1.0_R - COS1_**2)
      SIN2_ = sqrt(1.0_R - COS2_**2)
      TAN1_ = SIN1_/COS1_
      TAN2_ = SIN2_/COS2_
      SEC1_ = 1.0_R/COS1_
      SEC2_ = 1.0_R/COS2_

      COS_XI_ = COS1_*COS2_ + SIN1_*SIN2_*COS_PHI

      D = sqrt(TAN1_**2 + TAN2_**2 - 2.0_R*TAN1_*TAN2_*COS_PHI)

      COS_T = H2B * (sqrt(D**2 + (TAN1_*TAN2_*SIN_PHI)**2))/(SEC1_ + SEC2_)
      SIN_T = sqrt(1.0_R - COS_T**2)
      T = acos(COS_T)

      O = (1.0_R/PI) * (T - SIN_T*COS_T) * (SEC1_ + SEC2_)

      K_GEO = O - SEC1_ - SEC2_ + 0.5_R*(1.0_R+COS_XI_)*SEC1_*SEC2_
      !╰────────────────────────────────────────────────╯!

      LSRT_REFLECTION = F_ISO + F_VOL*K_VOL + F_GEO*K_GEO

      return

  end function

end module shdSurf
