!╭────────────────────────────────────────────────────────────────────────────╮!
! @LICENSE_HEADER_START@
!
!   This file is part of SHDOMNG.
!
!   --
!   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
!   --
!
!   Copyright (C) 2025-
!   EaR³T Community: https://discord.gg/ntqsguwaWv
!
!   SHDOMNG is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   SHDOMNG is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SHDOMNG. If not, see <http://www.gnu.org/licenses/>.
!
! @LICENSE_HEADER_END@
!
! Acknowledgement:
!   SHDOM was originally developed by K. Franklin Evans and is available
!   at https://coloradolinux.com/shdom.
!
! History:
!   Created by:  K. Franklin Evans
!                evans@nit.colorado.edu
!
!   Modified by: Hong Chen
!                hong.chen@hchen.com
!
!                Vikas Nataraja
!                vikas.hanasogenataraja@lasp.colorado.edu
!╰────────────────────────────────────────────────────────────────────────────╯!

module shd_sfc

    use globals_const, only: pi_, Dpi_, deg2rad_, rad2deg_
    use globals_cfg,   only: isPol_, isMpi_
    use globals_obj,   only: prp_obj_, &
                             sfc_obj_

    use shd_util, only: shd_util__computePlanck, &
                        shd_util__getBound

    implicit none

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_sfc__interpSfc
    public :: shd_sfc__brdfSurface

    ! surface brdf (including ocean brdf)
    !╭────────────────────────────────────────────────╮!
    public :: shd_sfc__brdfFresnel
    public :: shd_sfc__brdfDiner
    public :: shd_sfc__brdfRPV
    public :: shd_sfc__brdfRTLSR
    public :: shd_sfc__brdfJiao
    public :: shd_sfc__brdfOcean
    !╰────────────────────────────────────────────────╯!

    ! ocean brdf sub-modules (originally from 6SV code,
    !   https://salsa.umd.edu/6spage.html)
    !╭────────────────────────────────────────────────╮!
    public :: computeRefWhitecaps
    public :: computeRefSunglint
    public :: computeRefBackscatter

    public :: getRefracIndexWater
    public :: computeRefFresnel
    public :: computeRefMorel
    !╰────────────────────────────────────────────────╯!

contains

subroutine shd_sfc__interpSfc(&
    NptsBot, BCPTR, GRIDPOS, &
    SRCTYPE, WAVENO, waveLen, UNITS, &
    brdf_Sfc, brdf_Rtm)

    !   Bilinearly interpolates the surface properties from the regular
    ! surface grid to the boundary points of the adaptive grid.
    ! The locations of the NptsBot grid points on the bottom boundary
    ! (pointed to by BCPTR) are in the GRIDPOS array.  The input regular
    ! surface info is specified by the number of grid lines (sfc_obj_%Nx,sfc_obj_%Ny),
    ! the grid spacing (sfc_obj_%dx,sfc_obj_%dy), the number of parameters at
    ! each grid point (sfc_obj_%Npar), and the surface parameters (brdf_Sfc).
    ! The output interpolated surface parameters at the NptsBot boundary
    ! points are in brdf_grid_Sfc.  The first parameter is input as
    ! temperature, but converted to Planck function for output.

    integer, intent(in) :: NptsBot, BCPTR(*)
    character(len=1), intent(in) :: SRCTYPE, UNITS

    real(kind=R__), intent(in)  :: brdf_Sfc(1:sfc_obj_%Npar, 1:sfc_obj_%Ny+1, 1:sfc_obj_%Nx+1)
    real(kind=R__), intent(out) :: brdf_Rtm(sfc_obj_%Npar, *)

    real :: WAVENO(2), waveLen
    real :: GRIDPOS(3, *)

    integer :: igrid, ibc, ix, iy, ipar, i1, i2, i3, i4
    real :: rx, ry, u, v, temp, planck
    logical :: doInterp

    doInterp = .true.
    if (sfc_obj_%Npar .gt. 6) then
        doInterp = .false.
    end if

    ! Loop over all bottom points
    do ibc = 1, NptsBot

        igrid = BCPTR(ibc)

        if  ((prp_obj_%Nx .gt. 1) .and. (prp_obj_%Ny .gt. 1) .and.&
            &(sfc_obj_%Nx .gt. 1) .and. (sfc_obj_%Ny .gt. 1) .and.&
            &((abs(prp_obj_%dx - sfc_obj_%dx) .gt. 1.0e-6) .or.&
             &(abs(prp_obj_%dy - sfc_obj_%dy) .gt. 1.0e-6))) then
            ! added half dx,dy (prp_obj_%dx/2.0-sfc_obj_%dx/2.0, prp_obj_%dy/2.0-sfc_obj_%dy/2.0) to
            ! find surface point that is closest to grid center
            !╭────────────────────────────────────────────────╮!
            rx = (GRIDPOS(1, igrid) + prp_obj_%dx/2.0 - sfc_obj_%dx/2.0 - sfc_obj_%xStart) / sfc_obj_%dx
            ry = (GRIDPOS(2, igrid) + prp_obj_%dy/2.0 - sfc_obj_%dy/2.0 - sfc_obj_%yStart) / sfc_obj_%dy
            !╰────────────────────────────────────────────────╯!
        else
            rx = (GRIDPOS(1, igrid) - sfc_obj_%xStart) / sfc_obj_%dx
            ry = (GRIDPOS(2, igrid) - sfc_obj_%yStart) / sfc_obj_%dy
        end if

        ix = max(1, min(sfc_obj_%Nx, int(rx) + 1))
        iy = max(1, min(sfc_obj_%Ny, int(ry) + 1))
        u  = max(0.0, min(1.0, rx - (ix - 1)))
        v  = max(0.0, min(1.0, ry - (iy - 1)))

        do ipar = 1, sfc_obj_%Npar
            ! only interpolate BRDF parameter, not surface type ID (currently stored at 8)
            if (doInterp) then
                ! bi-linear interpolation
                !╭────────────────────────────────────────────────╮!
                brdf_Rtm(ipar, ibc) =&
                    &  (1.0 - u)*(1.0 - v)*brdf_Sfc(ipar, iy,   ix)&
                    &+ (1.0 - u)*       v *brdf_Sfc(ipar, iy+1, ix)&
                    &+        u *(1.0 - v)*brdf_Sfc(ipar, iy,   ix+1)&
                    &+        u *       v *brdf_Sfc(ipar, iy+1, ix+1)
                !╰────────────────────────────────────────────────╯!
            else
                ! nearest interpolation
                !╭────────────────────────────────────────────────╮!
                brdf_Rtm(ipar, ibc) = brdf_Sfc(ipar, iy, ix)
                !╰────────────────────────────────────────────────╯!
            end if
        end do

        if (SRCTYPE .eq. 'S') then
            planck = 0.0
        else
            temp = brdf_Rtm(1, ibc)
            call shd_util__computePlanck(temp, UNITS, WAVENO, waveLen, planck)
        end if

        brdf_Rtm(1, ibc) = planck

    end do

    return

end subroutine shd_sfc__interpSfc

subroutine shd_sfc__brdfSurface(SFCTYPE, REFPARMS, waveLen,&
    &muOut, phiOut, muInc, phiInc, NSTOKES, refBRDF)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Returns the reflection matrix for the general bidirectional
    ! reflection distribution function of the specified type (SFCTYPE).
    ! The incident direction is (muInc,phiInc), and the outgoing direction
    ! is (muOut,phiOut) (MU is cosine of zenith angle, and PHI is the azimuth
    ! angle in radians).  The incident directions have mu<0 (downward),
    ! while the outgoing directions have mu>0 (upward). The reflection
    ! function is normalized so that for a Lambertian surface (uniform
    ! reflection) the returned value (refBRDF) is simply the albedo.
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
    ! T  RTLSR          k_iso, k_geo, k_vol
    ! J  RTLSR-Jiao     k_iso, k_geo, k_vol, k_j, alpha
    ! M  Mixed surface
    !╰────────────────────────────────────────────────────────────────────────────╯!

    character(len=1), intent(in) :: SFCTYPE
    integer :: NSTOKES
    integer :: K1, K2
    real    :: REFPARMS(*), waveLen, muInc, phiInc, muOut, phiOut, refBRDF(4, 4)
    real    :: dPhi
    integer, parameter :: iSFC_ID=7

    SFCTYPE_CASE: select case (SFCTYPE)

        case ('L', 'l')

            ! L or l: Lambertian surface BRDF is constant.
            !   (for testing, as done more efficiently by Lambertian routine
            refBRDF(1:NSTOKES, 1:NSTOKES) = 0.0
            refBRDF(1, 1) = REFPARMS(1)

        case ('O', 'o')

            ! O: Ocean BRDF from 6S modified by Norm Loeb
            if (isPol_) then
                stop 'Error [_brdfSurface]: Ocean BRDF from 6S is only for unpolarized case.'
            end if

            dPhi = phiInc-phiOut
            call shd_sfc__brdfOcean(waveLen, REFPARMS(1), REFPARMS(2), phiInc, -1.0, &
                & -muInc, muOut, dPhi, phiInc, refBRDF(1, 1))

        case ('R', 'r')

            ! R: Rahman, Pinty, and Verstraete
            if (isPol_) then
                stop 'Error [_brdfSurface]: RPV BRDF is only for unpolarized case.'
            end if

            dPhi = phiInc-phiOut-pi_
            call shd_sfc__brdfRPV(REFPARMS(1), REFPARMS(2),&
                &REFPARMS(3), -muInc, muOut, dPhi, refBRDF(1, 1))

        case ('T', 't')

            ! T: RTLSR model
            if (isPol_) then
                stop 'Error [_brdfSurface]: RTLSR BRDF is only for unpolarized case.'
            end if

            dPhi = phiInc-phiOut-pi_
            call shd_sfc__brdfRTLSR(REFPARMS(1), REFPARMS(2),&
                &REFPARMS(3), -muInc, muOut, dPhi, refBRDF(1, 1))

        case ('J', 'j')

            ! J: RTLSR-Jiao model
            if (isPol_) then
                stop 'Error [_brdfSurface]: RTLSR-Jiao BRDF is only for unpolarized case'
            end if

            !dPhi = phiInc-phiOut-pi_
            dPhi = phiInc-phiOut
            ! J: RTLSR-Jiao model
            call shd_sfc__brdfRTLSRJ(REFPARMS(1), REFPARMS(2),&
                &REFPARMS(3), REFPARMS(4), REFPARMS(5), -muInc, muOut, dPhi, refBRDF(1, 1))

        case ('W', 'w')

            ! W: Fresnel reflection with Gaussian distribution of slopes
            call shd_sfc__brdfFresnel(REFPARMS(1), REFPARMS(2),&
                &REFPARMS(3), muInc, muOut, phiInc, phiOut, NSTOKES, refBRDF)

        case ('D', 'd')

            ! D: Diner et al., depolarizing modified RPV + Fresnel reflection
            !   from microfacets with uniform or Gaussian distribution of s
            dPhi = phiOut - phiInc
            call shd_sfc__brdfDiner(REFPARMS(1), REFPARMS(2), REFPARMS(3),&
                &REFPARMS(4), REFPARMS(5),&
                &-muInc, muOut, dPhi, NSTOKES, refBRDF)

        case ('X', 'x')

            ! for Lambertian   (VL):IX, IY, TEMP,   ALB,    -999.0, -999.0, -999.0, -999.0, -999.0,  0.0
            ! for Ocean        (VO):IX, IY, TEMP,   WSPD,   PCL,    -999.0, -999.0, -999.0, -999.0,  1.0
            ! for RPV          (VR):IX, IY, TEMP,   RHO0,   KR,     THETA,  -999.0, -999.0, -999.0,  2.0
            ! for RTLSR        (VT):IX, IY, TEMP,   F_ISO,  F_GEO,  F_VOL,  -999.0, -999.0, -999.0,  3.0
            ! for RTLSR-J       (VJ):IX, IY, TEMP,   F_ISO,  F_GEO,  F_VOL,  F_JIAO, ALPHA,  -999.0,  4.0
            ! for Wave-Fresnel (VW):IX, IY, TEMP,   MRE,    MIM,    WSPD,   -999.0, -999.0, -999.0, -1.0
            ! for Diner        (VD):IX, IY, TEMP,   A,      KR,     B,      ZETA,   SIGMA,  -999.0, -2.0

            if ((REFPARMS(iSFC_ID) .gt. -0.01) .and. (NSTOKES .gt. 1)) then
                stop 'Error [_brdfSurface]: Mixed surface is only for unpolarized case.'
            end if

            MIXSFC_CASE: select case (int(REFPARMS(iSFC_ID)))

                ! 0.0
                case (0)

                    refBRDF(1:NSTOKES, 1:NSTOKES) = 0.0
                    refBRDF(1, 1) = REFPARMS(1)

                ! 1.0
                case (1)

                    dPhi = phiInc-phiOut
                    call shd_sfc__brdfOcean(waveLen, REFPARMS(1), REFPARMS(2), phiInc, -1.0, &
                        & -muInc, muOut, dPhi, phiInc, refBRDF(1, 1))

                ! 2.0
                case (2)

                    dPhi = phiInc-phiOut-pi_
                    call shd_sfc__brdfRPV(REFPARMS(1), REFPARMS(2),&
                        &REFPARMS(3), -muInc, muOut, dPhi, refBRDF(1, 1))

                ! 3.0
                case (3)

                    dPhi = phiInc-phiOut-pi_
                    call shd_sfc__brdfRTLSR(REFPARMS(1), REFPARMS(2),&
                        &REFPARMS(3), -muInc, muOut, dPhi, refBRDF(1, 1))

                ! 4.0
                case (4)

                    dPhi = phiInc-phiOut
                    ! J: RTLSR-Jiao model
                    call shd_sfc__brdfRTLSRJ(REFPARMS(1), REFPARMS(2),&
                        &REFPARMS(3), REFPARMS(4), REFPARMS(5), -muInc, muOut, dPhi, refBRDF(1, 1))

                ! -1.0
                case (-1)

                    call shd_sfc__brdfFresnel(REFPARMS(1), REFPARMS(2),&
                        &REFPARMS(3), muInc, muOut, phiInc, phiOut, NSTOKES, refBRDF)

                ! -2.0
                case (-2)

                    dPhi = phiOut - phiInc
                    call shd_sfc__brdfDiner(REFPARMS(1), REFPARMS(2), REFPARMS(3),&
                        &REFPARMS(4), REFPARMS(5),&
                        &-muInc, muOut, dPhi, NSTOKES, refBRDF)

                case default

                    write(6, *) 'SFC_ID=', REFPARMS(iSFC_ID)
                    stop 'Error [_brdfSurface]: Unknown BRDF type.'

                end select MIXSFC_CASE

        case default

            write(6, *) 'SFCTYPE=', SFCTYPE
            stop 'Error [_brdfSurface]: Unknown BRDF type.'

    end select SFCTYPE_CASE

    return
end subroutine shd_sfc__brdfSurface

subroutine shd_sfc__brdfFresnel(MRE, MIM, WINDSPEED,&
    &MUI, MUR, PHII, PHIR, NSTOKES, refBRDF)

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
    real MRE, MIM, WINDSPEED, MUI, MUR, PHII, PHIR, refBRDF(4, 4)
    integer I, J
    real*8 SIGMA2, DMUI, DMUR, DCOSI, DSINI, DCOSR, DSINR, DSI, DSR
    real*8 VI1, VI2, VI3, VR1, VR2, VR3, UNIT1, UNIT2, UNIT3
    real*8 FACT1, FACTOR, XI1
    real*8 TI1, TI2, TI3, TR1, TR2, TR3, PI1, PI2, PI3, PR1, PR2, PR3
    real*8 PIKR, PRKI, TIKR, TRKI, E1, E2, E3, E4
    real*8 VP1, VP2, VP3, DMOD, RDZ2, RDZ4, DCOEFF, DEX
    real*8 AF, AF11, AF12, AF21, AF22
    real*8 P, S1, S2, S3, XI, XXI, DCOT, T1, T2
    real*8 SHADOWI, SHADOWR, SHADOW
    complex*16 CN1, CN2, CXI2, C1, C2, CRPER, CRPAR
    complex*16 CF11, CF12, CF21, CF22
    complex*16 CI, C21, C22, CTTTP, CTTPT, CTTPP, CTPPT, CTPPP, CPTPP

    CN1 = (1.0d0, 0.0d0)
    CN2 = cmplx(MRE, MIM)
    SIGMA2 = max(0.0005d0, 0.0015d0 + 0.00256d0*WINDSPEED)
    DMUI = abs(dble(MUI))
    DMUR = dble(MUR)

    ! Cartesian components of the unit vectors of the incident and
    ! scattered beams
    if (abs(DMUI - 1d0) .lt. 1d-10) DMUI = 0.999999999999d0
    if (abs(DMUR - 1d0) .lt. 1d-10) DMUR = 0.999999999999d0
    DCOSI = cos(dble(PHII))
    DSINI = sin(dble(PHII))
    DCOSR = cos(dble(PHIR))
    DSINR = sin(dble(PHIR))
    DSI = sqrt(1d0 - DMUI*DMUI)
    DSR = sqrt(1d0 - DMUR*DMUR)
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
    CXI2 = sqrt(1d0 - (1d0 - XI1*XI1)*CN1*CN1/(CN2*CN2))
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
    TR1 =  DMUR*DCOSR
    TR2 =  DMUR*DSINR
    TR3 = -DSR
    PI1 = -DSINI
    PI2 =  DCOSI
    PI3 = 0d0
    PR1 = -DSINR
    PR2 =  DCOSR
    PR3 = 0d0
    PIKR = PI1*VR1 + PI2*VR2 + PI3*VR3
    PRKI = PR1*VI1 + PR2*VI2 + PR3*VI3
    TIKR = TI1*VR1 + TI2*VR2 + TI3*VR3
    TRKI = TR1*VI1 + TR2*VI2 + TR3*VI3
    E1 = PIKR*PRKI
    E2 = TIKR*TRKI
    E3 = TIKR*PRKI
    E4 = PIKR*TRKI
    CF11 =  E1*CRPER + E2*CRPAR
    CF12 = -E3*CRPER + E4*CRPAR
    CF21 = -E4*CRPER + E3*CRPAR
    CF22 =  E2*CRPER + E1*CRPAR

    ! Calculation of the Stokes reflection matrix
    VP1 = VI2*VR3 - VI3*VR2
    VP2 = VI3*VR1 - VI1*VR3
    VP3 = VI1*VR2 - VI2*VR1
    DMOD = (VP1*VP1 + VP2*VP2 + VP3*VP3)**2
    RDZ2 = UNIT3*UNIT3
    RDZ4 = RDZ2*RDZ2
    DEX = exp(-(UNIT1*UNIT1 + UNIT2*UNIT2)/(2*SIGMA2*RDZ2))
    DCOEFF = FACT1**2*DEX/(4*DMUI*DMUR*DMOD*RDZ4*2*SIGMA2)
    AF = 0.5d0*DCOEFF
    AF11 = abs(CF11)**2
    AF12 = abs(CF12)**2
    AF21 = abs(CF21)**2
    AF22 = abs(CF22)**2

    refBRDF(1, 1) = (AF11 + AF12 + AF21 + AF22)*AF
    if (NSTOKES .ge. 2) then
        refBRDF(1, 2) = (AF11 - AF12 + AF21 - AF22)*AF
        refBRDF(2, 1) = (AF11 - AF22 + AF12 - AF21)*AF
        refBRDF(2, 2) = (AF11 - AF12 - AF21 + AF22)*AF
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
        refBRDF(1, 3) = (-CTTTP - CPTPP)*DCOEFF
        refBRDF(2, 3) = (-CTTTP + CPTPP)*DCOEFF
        refBRDF(3, 1) = (-CTTPT - CTPPP)*DCOEFF
        refBRDF(3, 2) = (-CTTPT + CTPPP)*DCOEFF
        refBRDF(3, 3) = ( CTTPP + CTPPT)*DCOEFF
    end if

    if (NSTOKES .ge. 4) then
        refBRDF(1, 4) = -CI*(CTTTP + CPTPP)*DCOEFF
        refBRDF(2, 4) = -CI*(CTTTP - CPTPP)*DCOEFF
        refBRDF(3, 4) =  CI*(CTTPP - CTPPT)*DCOEFF
        refBRDF(4, 1) =  CI*(CTTPT + CTPPP)*DCOEFF
        refBRDF(4, 2) =  CI*(CTTPT - CTPPP)*DCOEFF
        refBRDF(4, 3) = -CI*(CTTPP + CTPPT)*DCOEFF
        refBRDF(4, 4) = (CTTPP - CTPPT)*DCOEFF
    end if

    ! Shadowing
    !╭────────────────────────────────────────────────────────────────────────────╮!
    P = Dpi_
    S1 = sqrt(2d0*SIGMA2/P)
    S3 = 1d0/(sqrt(2d0*SIGMA2))
    S2 = S3*S3
    DCOT = DMUI/sqrt(1d0 - DMUI**2)

    T1 = exp(-S2*DCOT**2)
    T2 = erfc(DCOT*S3)
    SHADOWI = 0.5d0*(S1*T1/DCOT - T2)

    DCOT = DMUR/sqrt(1d0 - DMUR**2)
    T1 = exp(-S2*DCOT**2)
    T2 = erfc(DCOT*S3)
    SHADOWR = 0.5d0*(S1*T1/DCOT - T2)

    SHADOW = 1d0/(1d0 + SHADOWI + SHADOWR)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    refBRDF(1:NSTOKES, 1:NSTOKES) = refBRDF(1:NSTOKES, 1:NSTOKES)*SHADOW

    return

end subroutine shd_sfc__brdfFresnel

subroutine shd_sfc__brdfDiner(A, K, B, ZETA, SIGMA,&
    &muInc, muOut, PHI, NSTOKES, refBRDF)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !   Computes polarized reflection from the Diner et al. (2012) model
    ! which two terms: 1) completely depolarizing volumetric scattering
    ! the modified RPV model, and 2) Fresnel reflection from randomly
    ! oriented microfacets with either a uniform or Gaussian distributio
    ! of slopes.  The three parameters of the modified-RPV model are A,
    ! and B.  The weight given to the Fresnel-microfacet term is ZETA (i
    ! ZETA=0 then only the refBRDF(1,1) element is non-zero with the
    ! modified-RPV model).  SIGMA is the standard deviation of microface
    ! slopes for the Gaussian orientation distribution, but if SIGMA<=0
    ! the uniform orientation distribution is used.  The complex index o
    ! refraction of the Fresnel microfacets is specified in this subrout
    ! The Diner et al. (2012) model leaves out the hotspot factor of the
    ! modified RPV model (and the original RPV), but the hotspot factor
    ! included here.
    !   The cosine of the incident zenith angle is muInc (muInc>0), the
    ! cosine of the outgoing zenith angle is muOut (muOut>0), and the relati
    ! azimuthal angle is PHI=phiOut-phiInc (radians).
    !   For all formulas and definitions, see the paper:
    ! Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
    ! V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
    ! Exploration of a polarized surface bidirectional reflectance model
    ! using the ground-based multiangle spectropolarimetric imager.
    ! Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer NSTOKES
    real A, K, B, ZETA, SIGMA, muInc, muOut, PHI
    real refBRDF(4, 4)
    integer K1, K2
    real SINTH1, SINTH2, COSPHI, COSSCATANG, TAN1, TAN2, CAPG, HOT
    real GAMMA, COSGAMMA, F11, F12, F33, F34
    real COSBETA, H, SINPHI, ALPHA1, ALPHA2
    real COS2ALPHA1, SIN2ALPHA1, COS2ALPHA2, SIN2ALPHA2
    complex SFCINDEX, EPSILON, D, RP, RS

    ! Another free parameter of the Fresnel microfacet model is the
    !   index of refraction specified here
    SFCINDEX = cmplx(1.5, 0.0)

    refBRDF(1:NSTOKES, 1:NSTOKES) = 0.0

    ! The modified RPV model (Martonchik et al., 1998, IEEE-TGARS, 36,
    !   is implemented here.  Equation 6 in Diner et al. leaves out the h
    SINTH1 = sqrt(1.0 - muInc**2)
    SINTH2 = sqrt(1.0 - muOut**2)
    COSPHI = cos(PHI)
    COSSCATANG = -muInc*muOut + SINTH1*SINTH2*COSPHI
    COSSCATANG = min(1.0, max(-1.0, COSSCATANG))
    TAN1 = SINTH1/muInc
    TAN2 = SINTH2/muOut
    CAPG = sqrt(abs(TAN1**2 + TAN2**2 + 2*TAN1*TAN2*COSPHI))
    HOT = 1 + (1 - A)/(1 + CAPG)
    refBRDF(1, 1) = A*((muInc + muOut)*muInc*muOut)**(K - 1)*exp(B*COSSCATANG)
    refBRDF(1, 1) = refBRDF(1, 1)*HOT
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
    F11 = 0.5d0*(abs(RP)**2 + abs(RS)**2)
    F12 = 0.5d0*(abs(RP)**2 - abs(RS)**2)
    F33 = real(RP*conjg(RS))
    F34 = -aimag(RP*conjg(RS))

    ! Calculate facet tilt angle, beta (eq 16b)
    COSBETA = 0.5d0*(muInc + muOut)/COSGAMMA

    ! Calculate the factor before the M*P*M matrix including the distri
    !   of microfacet orientations (Gaussian for SIGMA>0) (eq 13 & 14)
    !   Note: Diner's BRDF is multiplied by pi here
    if (SIGMA .gt. 0.0) then
        H = ZETA*exp(-0.5d0*(1/COSBETA**2 - 1)/SIGMA**2)&
            & /(8*SIGMA**2*muOut*muInc*COSBETA**4)
    else
        H = ZETA/(8*muOut*muInc*COSBETA)
    end if

    refBRDF(1, 1) = refBRDF(1, 1) + H*F11

    if (NSTOKES .ge. 2) then
        ! Calculate the polarization rotation angles (eq 20 & 21)
        SINPHI = sin(PHI)
        ALPHA1 = atan(SINTH2*SINPHI/(muOut*SINTH1 + SINTH2*muInc*COSPHI))
        ALPHA2 = atan(SINTH1*SINPHI/(SINTH2*muInc + muOut*SINTH1*COSPHI))
        COS2ALPHA1 = cos(2*ALPHA1)
        SIN2ALPHA1 = sin(2*ALPHA1)
        COS2ALPHA2 = cos(2*ALPHA2)
        SIN2ALPHA2 = sin(2*ALPHA2)

        ! Calculate the rotated Mueller matrix (M*P*M) for Fresnel reflec
        !   and multiply by the H factor to get the final BRDF
        refBRDF(1, 2) = H*F12*COS2ALPHA1
        refBRDF(2, 1) = H*F12*COS2ALPHA2
        refBRDF(2, 2) = H*(F11*COS2ALPHA1*COS2ALPHA2&
            & + F33*SIN2ALPHA1*SIN2ALPHA2)
    end if

    if (NSTOKES .ge. 3) then
        refBRDF(1, 3) = -H*F12*SIN2ALPHA1
        refBRDF(2, 3) = H*(-F11*SIN2ALPHA1*COS2ALPHA2&
            & + F33*COS2ALPHA1*SIN2ALPHA2)
        refBRDF(3, 1) = -H*F12*SIN2ALPHA2
        refBRDF(3, 2) = H*(-F11*COS2ALPHA1*SIN2ALPHA2&
            & + F33*SIN2ALPHA1*COS2ALPHA2)
        refBRDF(3, 3) = H*(F11*SIN2ALPHA1*SIN2ALPHA2&
            & + F33*COS2ALPHA1*COS2ALPHA2)
    end if

    if (NSTOKES .ge. 4) then
        refBRDF(1, 4) = 0.0
        refBRDF(2, 4) = H*F34*SIN2ALPHA2
        refBRDF(3, 4) = H*F34*COS2ALPHA2
        refBRDF(4, 1) = 0.0
        refBRDF(4, 2) = -H*F34*SIN2ALPHA1
        refBRDF(4, 3) = -H*F34*COS2ALPHA1
        refBRDF(4, 4) = H*F33
    end if

    return

end subroutine shd_sfc__brdfDiner

subroutine shd_sfc__brdfRPV(RHO0, K, THETA, muInc, muOut, PHI, refBRDF)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !   Computes the Rahman, Pinty, Verstraete BRDF.  The incident
    ! and outgoing cosine zenith angles are muInc and muOut, respectively,
    ! and the relative azimuthal angle is PHI.  In this case the incident
    ! direction is where the radiation is coming from (i.e. opposite of
    ! discrete ordinate), so muInc>0 and the hot spot is muOut=muInc and PHI=0
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

    real, intent(in)  :: RHO0, K, THETA, muInc, muOut, PHI
    real, intent(out) :: refBRDF
    real :: mu_min, x1, x2
    real :: M, F, H, COSPHI, SIN1, SIN2, COSG, TAN1, TAN2, CAPG

    ! calculate M
    !╭────────────────────────────────────────────────╮!
    mu_min = 0.03
    x1 = muInc
    if (x1 .lt. mu_min) x1 = mu_min
    x2 = muOut
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

    refBRDF = RHO0*M*F*H

    return

end subroutine shd_sfc__brdfRPV

subroutine shd_sfc__brdfRTLSR(f_iso, f_geo, f_vol, muInc, muOut, dPhi, refBRDF)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Computes the RossThick-LiSparseReciprocal BRDF used by MODIS and VIIRS BRDF/
    !   Albedo products.
    !
    ! Input:
    !   f_iso: isotropic scattering coefficient
    !   f_geo: geometric coefficient
    !   f_vol: volumetric scattering coefficient
    !   muInc: cosine of zenith angle of the incident light (solar zenith angle)
    !   muOut: cosine of zenith angle of the receiver (viewing zenith angle)
    !   dPhi: relative azimuth angle, 0 pointing east, positive counter-clockwise
    !
    ! Output:
    !   refBRDF: BRDF reflectance
    !
    ! Note:
    !   When muInc>0 and the hot spot is muOut=muInc and dPhi=0
    !
    ! Code reference:
    !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
    ! Formula reference:
    !   MODIS BRDF/Albedo Product ATBD v5.0, equations 37-44
    !╰────────────────────────────────────────────────────────────────────────────╯!
    real(R__), intent(in)  :: f_iso, f_geo, f_vol, muInc, muOut, dPhi
    real(R__), intent(out) :: refBRDF

    real(R__), parameter :: b2r=1.0, h2b=2.0

    real(R__) :: cosI, cosO, sinI, sinO, tanI, tanO, cos_phi, sin_phi
    real(R__) :: K_geo, K_vol
    real(R__) :: xi, cos_xi, sin_xi
    real(R__) :: cosI_p, sinI_p, tanI_p, secI_p
    real(R__) :: cosO_p, sinO_p, tanO_p, secO_p
    real(R__) :: cos_xi_p, D_2, cos_t, sin_t, t, Omicron
    real(R__) :: cosMin


    ! calculate trigonometry
    !╭────────────────────────────────────────────────╮!
    cosMin = 0.03

    cosI = muInc
    cosO = muOut

    if (cosI .lt. cosMin) cosI = cosMin
    if (cosO .lt. cosMin) cosO = cosMin

    sinI = sqrt(1.0 - cosI*cosI)
    sinO = sqrt(1.0 - cosO*cosO)
    tanI = sinI/cosI
    tanO = sinO/cosO
    cos_phi = cos(dPhi)
    sin_phi = sin(dPhi)
    !╰────────────────────────────────────────────────╯!

    ! K_vol
    !╭────────────────────────────────────────────────╮!
    cos_xi = (cosI*cosO) + (sinI*sinO*cos_phi)
    if (cos_xi .lt. -1.0) cos_xi = -1.0
    if (cos_xi .gt.  1.0) cos_xi =  1.0

    sin_xi = sqrt(1.0 - cos_xi*cos_xi)
    xi = acos(cos_xi)

    ! Eq. 38
    K_vol = ((pi_/2.0 - xi)*cos_xi + sin_xi)/(cosI+cosO) - (pi_/4.0)
    !╰────────────────────────────────────────────────╯!

    ! K_geo
    !╭────────────────────────────────────────────────╮!
    ! Eq. 44 (left) [used tan instead of atan, can save a few calculations]
    tanI_p = b2r*tanI
    cosI_p = 1.0/sqrt(1.0+tanI_p*tanI_p)
    sinI_p = tanI_p*cosI_p
    secI_p = 1.0/cosI_p

    ! Eq. 44 (right) [used tan instead of atan, can save a few calculations]
    tanO_p = b2r*tanO
    cosO_p = 1.0/sqrt(1.0+tanO_p*tanO_p)
    sinO_p = tanO_p*cosO_p
    secO_p = 1.0/cosO_p

    ! Eq. 43
    cos_xi_p = (cosI_p*cosO_p) + (sinI_p*sinO_p*cos_phi)

    ! Eq. 42 (without square root)
    D_2 = tanI_p*tanI_p + tanO_p*tanO_p - 2.0*tanI_p*tanO_p*cos_phi
    if (D_2 .lt. 0.0) D_2 = 0.0

    ! Eq. 41
    cos_t = h2b * (sqrt(D_2 + (tanI_p*tanO_p*sin_phi)**2))/(secI_p + secO_p)
    if (cos_t .lt. -1.0) cos_t = -1.0
    if (cos_t .gt.  1.0) cos_t =  1.0
    sin_t = sqrt(1.0 - cos_t*cos_t)
    t = acos(cos_t)

    ! Eq. 40
    Omicron = (1.0/pi_) * (t - sin_t*cos_t) * (secI_p + secO_p)

    ! Eq. 39
    K_geo = Omicron - secI_p - secO_p + 0.5d0*(1.0+cos_xi_p)*secI_p*secO_p
    !╰────────────────────────────────────────────────╯!

    ! Eq. 37
    refBRDF = f_iso + f_vol*K_vol + f_geo*K_geo

    if (refBRDF .lt. 0.0) refBRDF = 0.0

    return

end subroutine shd_sfc__brdfRTLSR

! subroutine shd_sfc__brdfRTLSRKZ(F_ISO, F_GEO, F_VOL, F_KZ, muInc, muOut, PHI, refBRDF)

!     !╭────────────────────────────────────────────────────────────────────────────╮!
!     ! Computes the RossLi-KokhanovskyZege BRDF (combines RTLSR model above with a fourth kernel for snow from Kokhanovsky-Zege).
!     ! muInc: cosine of zenith angle of the incident light (solar zenith angle)
!     ! muOut: cosine of zenith angle of the receiver (viewing zenith angle)
!     ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
!     !
!     ! When muInc>0 and the hot spot is muOut=muInc and PHI=0
!     !
!     ! Code reference:
!     !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
!     !
!     ! Formula reference:
!     !    KZ04: Kokhanovsky & Zege (2004), doi:10.1364/ao.43.001589
!     !    KB12: Kokhanovsky & Breon (2012), doi:10.1109/lgrs.2012.2185775
!     !╰────────────────────────────────────────────────────────────────────────────╯!

!     real(4), intent(in)  :: F_ISO, F_GEO, F_VOL, F_KZ, muInc, muOut, PHI
!     real(4), intent(out) :: refBRDF
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

!     COS1 = muInc
!     COS2 = muOut

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

!     K_VOL = ((pi_/2.0 - XI)*COS_XI + SIN_XI)/(COS1+COS2) - (pi_/4.0)
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

!     O = (1.0/pi_) * (T - SIN_T*COS_T) * (SEC1_ + SEC2_)

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
!     XI_KZ_DEG = XI_KZ * 180.0 / pi_ ! convert to degrees to calculate phase function

!     ! calculate phase function p of scattering angle xi
!     P_pha = 11.1 * exp(-0.087 * XI_KZ_DEG) + 1.1 * exp(-0.014 * XI_KZ_DEG)

!     ! calculate R0 (reflectance at 0 absorption by snow medium)
!     R0 = (a + (b * (muInc + muOut)) + (c * muInc * muOut) + P_pha) / (4 * (muInc + muOut))

!     ! calculate f (escape function)
!     ! note that this is a parameterization of the escape function
!     ! but one can run RT simulation to get the full escape function
!     K0_sza  = (3/7) * (1 + (2 * muInc))
!     K0_vza  = (3/7) * (1 + (2 * muOut))
!     f = K0_sza * K0_vza / R0


!     refBRDF = F_ISO + F_VOL*K_VOL + F_GEO*K_GEO + F_KZ*K_KZ

!     if (refBRDF .lt. 0.0) refBRDF = 0.0 ! relfectance cannot be negative

!     return

! end subroutine shd_sfc__brdfRTLSRKZ

subroutine shd_sfc__brdfRTLSRJ(F_ISO, F_GEO, F_VOL, F_J, ALPHA, muInc, muOut, PHI, refBRDF)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Computes the RossLi-Jiao BRDF (combines RTLSR model above with a fourth kernel for snow from Jiao et al.).
    ! muInc: cosine of zenith angle of the incident light (solar zenith angle)
    ! muOut: cosine of zenith angle of the receiver (viewing zenith angle)
    ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
    ! ALPHA: parameter for the empirical function to control the forward scattering. Jiao et al. (2019) suggest ALPHA=0.2-0.3.
    !
    ! When muInc>0 and the hot spot is muOut=muInc and PHI=0
    !
    ! Code reference:
    !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
    !
    ! Formula reference:
    !    Jiao et al. (2019) https://doi.org/10.1016/j.rse.2018.11.001
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(4), intent(in)  :: F_ISO, F_GEO, F_VOL, F_J, ALPHA, muInc, muOut, PHI
    real(4), intent(out) :: refBRDF
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

    COS1 = muInc
    COS2 = muOut

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

    K_VOL = ((pi_/2.0 - XI)*COS_XI + SIN_XI)/(COS1+COS2) - (pi_/4.0)
    !╰────────────────────────────────────────────────╯!

    ! K_GEO
    !╭────────────────────────────────────────────────╮!
    B2R = 1.0
    H2B = 2.0

    THETA1_ = atan(B2R*TAN1)
    COS1_ = cos(THETA1_)
    if (COS1_ .lt. -1.0) COS1_ = -1.0
    if (COS1_ .gt.  1.0) COS1_ =  1.0
    SIN1_ = sqrt(1.0 - COS1_**2)
    TAN1_ = SIN1_/COS1_
    SEC1_ = 1.0/COS1_

    THETA2_ = atan(B2R*TAN2)
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

    O = (1.0/pi_) * (T - SIN_T*COS_T) * (SEC1_ + SEC2_)

    K_GEO = O - SEC1_ - SEC2_ + 0.5d0*(1.0+COS_XI_)*SEC1_*SEC2_
    !╰────────────────────────────────────────────────╯!

    ! get Jiao kernel value K_J with normalization
    !╭────────────────────────────────────────────────╮!
    call shd_sfc__brdfJiao(muInc, muOut, PHI, ALPHA, 1, K_J)
    !╰────────────────────────────────────────────────╯!


    refBRDF = F_ISO + F_VOL*K_VOL + F_GEO*K_GEO + F_J*K_J

    if (refBRDF .lt. 0.0) refBRDF = 0.0 ! relfectance cannot be negative

    return

end subroutine shd_sfc__brdfRTLSRJ

subroutine shd_sfc__brdfJiao(muInc, muOut, PHI, ALPHA, NORMALIZE, K_J)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Computes the Jiao BRDF for snow.
    ! muInc: cosine of zenith angle of the incident light (solar zenith angle)
    ! muOut: cosine of zenith angle of the receiver (viewing zenith angle)
    ! PHI: relative azimuth angle, 0 pointing east, positive counter-clockwise
    ! ALPHA: parameter for the empirical function to control the forward scattering. Jiao et al. (2019) suggest ALPHA=0.2-0.3.
    ! NORMALIZE: whether to add a normalizer term to the reflectance (1 for yes, 0 for no). If combining this kernel with an isotropic kernel (such as the case with RTLSR/RTLSR), set NORMALIZE=1.
    !
    ! When muInc>0 and the hot spot is muOut=muInc and PHI=0
    !
    ! Code reference:
    !   https://github.com/vikasnataraja/sif-sat/blob/master/sifsat/kernels.py
    !
    ! Formula reference:
    !    Jiao et al. (2019) https://doi.org/10.1016/j.rse.2018.11.001
    !╰────────────────────────────────────────────────────────────────────────────╯!
    real(4), intent(in)  :: muInc, muOut, PHI, ALPHA
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

    COS1 = muInc
    COS2 = muOut

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
    XI_J_DEG = XI_J * rad2deg_ ! convert to degrees to calculate phase function

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
    wvl, wspd, pcl, phiW, xsal, &
    muInc, muOut, dPhi, phiS, &
    refBRDF)

    ! Description:
    !   Calculate ocean BRDF based on given windspeed, pigment concentration,
    !   salinity etc. The code was originally from 6S code (v2.1)
    !   https://salsa.umd.edu/6spage.html
    !
    ! Input:
    !   wvl   : wavelength of the computation (in micrometer)
    !   wspd  : wind speed (in m/s) [value not allowed below 0.25 m/s]
    !   pcl   : pigment concentration (in mg/m^3)
    !   phiW  : azimuth wind direction (in radians)
    !   xsal  : salinity (in ppt) [if xsal<0 then 34.3ppt by default]
    !   muInc : cosine of the zenith angle of the incident radiation
    !   muOut : cosine of the zenith angle of the reflected radiation
    !   dPhi  : relative azimuth angle (between 0 and 2*pi, in radians)
    !   phiS  : azimuth of sun
    !
    !   note:
    !   Originally, the following settings is enabled:
    !     For the direct solar radiation, the wind direction
    !     is fixed at the azimuth angle of 0 degree so that
    !     phiW = 0.0
    !   Now after benchmarking with libRadtran
    !     phiW is set to phiS, which means wind direction is
    !     the same as the direction of Sun
    !
    ! Output:
    !   refBRDF: the total reflectance of the sea water
    !
    ! Reference:
    ! Original Fortran 77 code was mostly adapted from 6S code (v2.1), which can be
    !   referenced to the following two papers:
    !
    ! Kotchenova, S. Y., Vermote, E. F., Matarrese, R., and Klemm, F. J., Jr.:
    !   Validation of a vector version of the 6S radiative transfer code
    !   for atmospheric correction of satellite data. Part I: Path radiance,
    !   Appl. Optics, 45, 6762–6774, doi:10.1364/ao.45.006762, 2006.
    ! Kotchenova, S. Y. and Vermote, E. F.: Validation of a vector version
    !   of the 6S radiative transfer code for atmospheric correction of
    !   satellite data. Part II. Homogeneous Lambertian and anisotropic
    !   surfaces, Appl. Optics, 46, 4455–4464, doi:10.1364/ao.46.004455, 2007.
    !
    ! The algorithms can be referenced to the following papers:
    !
    ! Cox, C. and Munk, W.: Some problems in optical oceanography,
    !   J. Mar. Res., 14, 63–78, 1954. 
    ! Vermote, E. F., Tanré, D., Deuze, J. L., Herman, M., and Morcette, J. J.:
    !   Second simulation of the satellite signal in the solar spectrum, 6S:
    !   An overview, IEEE T. Geosci. Remote, 35, 675–686, 1997. 
    ! Koepke, P.: Effective reflectance of oceanic whitecaps,
    !   J. Geophys. Res., 89(D3), 3669–3676,
    !   https://doi.org/10.1364/AO.23.001816, 1984.
    ! Austin, R. W.: The remote sensing of spectral radiance from below the
    !   ocean surface, in: Optical Aspects of Oceanography, Academic Press, 1974.

    real, intent(in)  :: wvl, wspd, pcl, phiW, xsal
    real, intent(in)  :: muInc, muOut, dPhi, phiS
    real, intent(out) :: refBRDF

    real :: nr, ni, windSpeed, waveLen, phiSV, phiSW
    real :: angZenView, angZenSun, angAziSunView, angAziSunWind
    real :: propWC, refWC, refGL, refSW

    ! preprocess input parameters
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! minimum and maximum wavelength cut (400 to 4000 nm)
    if (wvl .lt. 0.4) then
        waveLen = 0.4
    else if (wvl .gt. 4.0) then
        waveLen = 4.0
    else
        waveLen = wvl
    end if

    ! minimum wind speed cut (too low wind speed will cause extreme sunglint)
    windSpeed = max(wspd, 0.25)

    ! solar zenith angle
    if (muInc .le. 0.028) then
        angZenSun = acos(0.028)*rad2deg_
    else
        angZenSun = acos(muInc)*rad2deg_
    end if

    ! viewing zenith angle
    if (muOut .le. 0.028) then
        angZenView = acos(0.028)*rad2deg_
    else
        angZenView = acos(muOut)*rad2deg_
    end if

    ! relative azimuth between Sun and Viewer (0 to 2PI)
    if (dPhi .lt. 0.0) then
        phiSV = dPhi + 2.0*pi_
    else if (dPhi .ge. 2.0*pi_) then
        phiSV = dPhi - 2.0*pi_
    else
        phiSV = dPhi
    end if

    ! relative azimuth between Sun and Wind (0 to 2PI)
    ! note: disabled Sun-Wind azimuthal difference (phiW=phiS)
    !       after being disabled, calculations agree perfectly with libRadtran
    phiSW = phiS - phiW
    if (phiSW .lt. 0.0)     phiSW = phiSW + 2.0*pi_
    if (phiSW .ge. 2.0*pi_) phiSW = phiSW - 2.0*pi_
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! get water refractive indices based on wavelenth (micron) and salinity
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call getRefracIndexWater(waveLen, xsal, nr, ni)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Compute whitecaps reflectance (Lambertian)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call computeRefWhitecaps(waveLen, windSpeed, propWC, refWC)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! refGL: reflectance from sunglint (specular)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    angAziSunWind = phiSW*rad2deg_
    angAziSunView = 180.0 - phiSV*rad2deg_
    call computeRefSunglint(windSpeed, nr, ni, angZenSun, angZenView, angAziSunWind, angAziSunView, refGL)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! refSW: reflectance due to scattering of underwater (Lambertian)
    !        refSW1: from Morel et al., 1988
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call computeRefBackscatter(waveLen, windSpeed, pcl, nr, ni, angZenSun, angZenView, refSW)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Total reflectance of sea water
    !   corresponds to Eq. (C1, Appendix C) in Kotchenova and Vermote, 2007
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !   refWC: reflectance from whitecaps
    !   refGL: reflectance (specular) from sunglint
    !   refSW: reflectance due to scattering of underwater
    !   propWC: proportion of whitecaps
    refBRDF = refWC + (1.0 - propWC)*refGL + (1.0 - refWC)*refSW
    !╰────────────────────────────────────────────────────────────────────────────╯!

    return

end subroutine shd_sfc__brdfOcean

subroutine computeRefWhitecaps(waveLen, windSpeed, propWC, refWC)

    ! Description: Effective reflectance of the whitecaps
    ! input:
    !   waveLen: wavelength (in micrometer)
    !   windSpeed: wind speed (in meter/second)
    !
    ! output:
    !   refWC: effective reflectance of whitecaps
    !
    ! Reference:
    ! Koepke, P.: Effective reflectance of oceanic whitecaps,
    !   J. Geophys. Res., 89(D3), 3669–3676,
    !   https://doi.org/10.1364/AO.23.001816, 1984.
    ! Monahan, E. C. and O’Muircheartaigh, I.: Optimal power-law
    !   description of oceanic white cap coverage dependence on wind speed,
    !   J. Phys. Oceanogr., 10, 2094–2099, 1980.

    real, intent(in)  :: waveLen, windSpeed
    real, intent(out) :: propWC, refWC

    integer :: iwvl
    real :: wvlP, refEff, refTable(39)

    ! effective reflectance of the whitecaps (Koepke, 1984)
    data refTable/&
        0.220, 0.220, 0.220, 0.220, 0.220, 0.220, 0.215, 0.210, 0.200, 0.190,&
        0.175, 0.155, 0.130, 0.080, 0.100, 0.105, 0.100, 0.080, 0.045, 0.055,&
        0.065, 0.060, 0.055, 0.040, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000,&
        0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000, 0.000/

    iwvl = 1 + int((waveLen-0.2)/0.1)
    wvlP = 0.5 + (iwvl - 1)*0.1
    refEff = refTable(iwvl+1) + (waveLen-wvlP)/0.1 * (refTable(iwvl) - refTable(iwvl+1))

    ! corresponds to Eq. (3) in Koepke, 1984 (originally from Monahan and O’Muircheartaigh, 1980)
    propWC = 2.95e-6*(windSpeed**3.52)

    ! corresponds to Eq. (5) in Koepke, 1984
    refWC = refEff*propWC

    return

end subroutine computeRefWhitecaps

subroutine computeRefSunglint(&
    windSpeed, nr, ni, &
    angZenSun, angZenView, &
    angAziSunWind, angAziSunView, &
    refGL)

    ! Description:
    !   Calculate ocean reflectance from sun glint
    !
    ! Input:
    !   windSpeed: speed of the wind (in m/s)
    !   nr: index of refraction of the sea water
    !   ni: extinction coefficient of the sea water
    !   angZenSun: solar zenith angle (in deg.)
    !   angZenView: view zenith angle (in deg.)
    !   angAziSunWind: azimuth of the sun minus azimuth of the wind (in deg.)
    !   angAziSunView: relative azimuth (sun-satellite, in deg.)
    !
    ! Output:
    !   refGL: reflectance of the sun glint
    !
    ! Reference:
    ! C. Cox and W. Munk, “Measurement of the roughness of the sea surface
    !   from photographs of the sun’s glitter,” J. Opt. Soc. Am. 44, 838–850 (1954).
    ! C. Cox and W. Munk, “Statistics of the sea surface derived from sun
    !   glitter,” J. Mar. Res. 13, 198–227 (1954).
    ! C. Cox and W. Munk, “Some problems in optical oceanography,”
    !   J. Mar. Res. 14, 63–78 (1955).
    ! M. Born and E. Wolf, Principles of Optics, 5th ed. (Pergamon, 1975).

    real, intent(in)  :: windSpeed, nr, ni
    real, intent(in)  :: angZenSun, angZenView, angAziSunWind, angAziSunView
    real, intent(out) :: refGL

    real :: phiSV, sinThetaS, sinThetaV, cosThetaS, cosThetaV, cosBeta
    real :: Zx, Zy, tanTilt, betaTilt
    real :: cosChi, sinChi, cos2Chi, coefFR
    real :: sigma2C, sigma2U, sigmaC, sigmaU
    real :: C21, C03, C40, C04, C22, phiSW, xi, eta, xi2, eta2, coef, Pxy

    logical :: isClean = .true. ! True -> clean surface; False -> slick surface

    cosThetaS = cos(angZenSun *deg2rad_)
    cosThetaV = cos(angZenView*deg2rad_)
    sinThetaS = sin(angZenSun *deg2rad_)
    sinThetaV = sin(angZenView*deg2rad_)

    phiSV = angAziSunView*deg2rad_

    ! calculate tilt angle beta (in radians)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    Zx =            (-sinThetaV*sin(phiSV))/(cosThetaS + cosThetaV)
    Zy = (sinThetaS + sinThetaV*cos(phiSV))/(cosThetaS + cosThetaV)
    tanTilt  = sqrt(Zx*Zx + Zy*Zy)
    betaTilt = atan(tanTilt)
    cosBeta  = cos(betaTilt)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Anisotropic Gaussian distribution probability
    ! Cox and Munk, 1954 - Measurement of the Roughness of the Sea
    !   Surface from Photographs of the Sun's Glitter
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Sect. 6.3
    if (isClean) then
        sigma2C = 0.003 + 0.00192*windSpeed
        sigma2U =         0.00316*windSpeed
    else
        sigma2C = 0.003 + 0.00084*windSpeed
        sigma2U = 0.005 + 0.00078*windSpeed
    end if
    sigmaC = sqrt(sigma2C)
    sigmaU = sqrt(sigma2U)

    ! Sect. 6.4
    if (isClean) then
        C21 = 0.01 - 0.0086*windSpeed
        C03 = 0.04 - 0.033 *windSpeed
        C40 = 0.40
        C22 = 0.12
        C04 = 0.23
    else
        C21 = 0.00
        C03 = 0.02
        C40 = 0.36
        C22 = 0.10
        C04 = 0.26
    end if

    ! relative azimuth between Sun and Wind
    phiSW = angAziSunWind*deg2rad_

    xi   = ( cos(phiSW)*Zx + sin(phiSW)*Zy)/SigmaC
    xi2  = xi*xi

    eta  = (-sin(phiSW)*Zx + cos(phiSW)*Zy)/SigmaU
    eta2 = eta*eta

    ! square bracket in Eq. 18
    ! note: in Cox and Munk, 1954 Statistics paper,
    !   `eta` term is not included in [C21/2.0 *(xi2  - 1.0)*eta]
    !   in Eq. 3, I believe is a typo
    coef =  1.0 - &
            C21/2.0 *(xi2  - 1.0)*eta - &
            C03/6.0 *(eta2 - 3.0)*eta + &
            C40/24.0*(xi2*xi2 - 6.0*xi2 + 3.0) + &
            C22/4.0 *(xi2 - 1.0)*(eta2 - 1.0) + &
            C04/24.0*(eta2*eta2 - 6.0*eta2 + 3.0)

    ! Eq. 18
    Pxy = exp(-(xi2+eta2)/2.0)*coef/(2.0*pi_*sigmaU*sigmaC)

    ! note: this is added (Jan 20, 2026, by Hong Chen)
    !   to prevent negative probability
    !   occured when windSpeed=20 m/s and SunZenith=30 degree
    Pxy = max(0.0, Pxy)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Compute coefficient coefFR (Fresnel's reflection)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    cos2Chi = cosThetaV*cosThetaS + sinThetaV*sinThetaS*cos(phiSV)
    if (cos2Chi .gt.  1.0) cos2Chi =  0.99999999999
    if (cos2Chi .lt. -1.0) cos2Chi = -0.99999999999
    cosChi = sqrt(0.5*(1.0 + cos2Chi))
    sinChi = sqrt(0.5*(1.0 - cos2Chi))

    ! note: the following is added after reviewing 6S code (v2.1)
    if (cosChi .gt.  1.0) cosChi =  0.99999999999
    if (cosChi .lt. -1.0) cosChi = -0.99999999999
    if (sinChi .gt.  1.0) sinChi =  0.99999999999
    if (sinChi .lt. -1.0) sinChi = -0.99999999999

    call computeRefFresnel(nr, ni, cosChi, sinChi, coefFR)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Compute specular reflection from sunglint
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Eq. C3 (in Appendix C) in Kotchenova et al., 2007
    refGL = pi_*coefFR*Pxy/(4.0*cosThetaS*cosThetaV*(cosBeta**4))
    !╰────────────────────────────────────────────────────────────────────────────╯!

    return

end subroutine computeRefSunglint

subroutine computeRefBackscatter(&
    waveLen, windSpeed, pcl, nr, ni, &
    angZenSun, angZenView, &
    refSW)

    real, intent(in)  :: waveLen, windSpeed, pcl, nr, ni
    real, intent(in)  :: angZenSun, angZenView
    real, intent(out) :: refSW

    real :: refSW1, tds, tdv, a, n12
    integer :: iws1, iws2, isz1, isz2, ivz1, ivz2
    integer :: iang, iws

    real angbnd(5), wsbnd(6), tdsbnd(5, 6), tdvbnd(5, 6)

    data angbnd/&
        0.0, 45.0, 60.0, 75.0, 85.0/

    data wsbnd/&
        1.0, 3.0, 5.0, 7.0, 9.0, 20.0/

    data((tdsbnd(iang, iws), iws=1, 6), iang=1, 5)/&
        0.9787803, 0.9787738, 0.9787626, 0.9787467, 0.9787264,&
        0.9785573, 0.9706900, 0.9698871, 0.9691746, 0.9685547,&
        0.9680276, 0.9666586, 0.9479931, 0.9404608, 0.9385692,&
        0.9381815, 0.9384519, 0.9430056, 0.9690591, 0.9275920,&
        0.9058769, 0.8951812, 0.8899654, 0.8892645, 0.9980542,&
        0.9602273, 0.9114283, 0.8713799, 0.8417820, 0.7800314/

    data((tdvbnd(iang, iws), iws=1, 6), iang=1, 5)/&
        0.9787764, 0.9787535, 0.9787106, 0.9786453, 0.9785548,&
        0.9775019, 0.9692680, 0.9637051, 0.9564344, 0.9495727,&
        0.9438773, 0.9288712, 0.9225163, 0.9069787, 0.9044844,&
        0.9052351, 0.9068328, 0.9153687, 0.8048478, 0.8479503,&
        0.8678726, 0.8797889, 0.8878716, 0.9091171, 0.7294627,&
        0.8137348, 0.8453338, 0.8629867, 0.8745421, 0.9036854/

    ! Compute backscattered reflectance from the sea water (Lambertian)
    !   water reflectance below the sea surface
    call computeRefMorel(waveLen, pcl, refSW1)

    call shd_util__getBound( wsbnd, 1, 6,       windSpeed, iws1, iws2)
    call shd_util__getBound(angbnd, 1, 5,  angZenSun, isz1, isz2)
    call shd_util__getBound(angbnd, 1, 5, angZenView, ivz1, ivz2)
    tds = tdsbnd(isz1, iws1)
    tdv = tdvbnd(ivz1, iws1)

    ! albedo is a=re is taken from table 2 in Austin, 1974
    a=0.485

    ! Add change in solid angle from under to above to surface
    !   that account for 1/(n12*n12) decrease in sea water directional
    !   reflectance
    n12 = sqrt(nr*nr + ni*ni)

    refSW = (1.0/(n12*n12)) * refSW1*tds*tdv/(1.0 - a*refSW1)

end subroutine computeRefBackscatter

subroutine computeRefMorel(waveLen, pcl, refSW)

    ! Description:
    !   Spectral diffuse attenuation coefficient of Case I Waters as predicted
    !   by MOREL, 1998
    !
    ! Input:
    !   waveLen: wavelength (in micron)
    !   pcl: pigment concentration (in mg/m^3)
    !
    ! Output:
    !   refSW: reflectance due to scattering of underwater
    !
    ! Additional notes:
    ! according Morel,1988, we use:
    ! Kd    : spectral value of the attenuation coefficient for
    !         downwelling irradiance
    !         with: Kd=Kw+Xc*pcl**e
    ! Kw    : spectral value of the diffuse attenuation coefficient
    !         for pure oceanic water
    ! Xc, e : spectral coefficients to compute the diffuse attenuation
    !         coefficient for pigment
    ! bb    : total backscattering coefficient
    !         with: bb=0.5*bw+bbt*b
    ! bw    : spectral value of the molecular scattering coefficient of water
    ! bbt,b : parameters to compute the scattering coefficients of pigments
    ! refSW : reflectance of water below the surface
    !         with: refSW=(0.33/u)*(bb/Kd)  where u is depending of refSW
    !
    ! Reference:
    ! MOREL within the spectral range 400-700nm (1988, Journal of Geophys
    !   Research, Vol.93, No C9, pp 10749-10768)

    real, intent(in)  :: waveLen, pcl
    real, intent(out) :: refSW

    real :: Kw, Kd
    real :: Xc, e, bw, bb, b, bbt, u1, u2, err, refSW0
    real :: tKw(61), tXc(61), te(61), tbw(61)

    integer :: iwl

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

    ! refSW only applies to wavelength range 400-700 nm
    if ((waveLen .lt. 0.4) .or. (waveLen .gt. 0.7)) then
        refSW = 0.000
        return
    end if

    iwl = 1 + nint((waveLen - 0.4)/0.005)
    Kw = tKw(iwl)
    Xc = tXc(iwl)
    e  = te(iwl)
    bw = tbw(iwl)

    if (abs(pcl) .lt. 0.0001) then
        bb = 0.5*bw
        Kd = Kw
    else
        b   = 0.30*pcl**0.62
        bbt = 0.002 + 0.02*(0.5 - 0.25*log10(pcl))*0.550/waveLen
        bb  = 0.5*bw + bbt*b
        Kd  = Kw + Xc*pcl**e
    end if

    u1 = 0.75
    refSW0 = 0.33*bb/u1/Kd

    ! origianl code (uses goto) replaced by do
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! 50  u2 = 0.90d0*(1.0d0-refSW0)/(1.0d0+2.25d0*refSW0)
    !     refSW = 0.33d0*bb/u2/Kd
    !     err = abs((refSW - refSW0)/refSW)
    !     if (err .lt. 0.0001d0) return
    !     refSW0 = refSW
    !     goto 50
    !     return
    !╰────────────────────────────────────────────────────────────────────────────╯!

    do
        u2  = 0.90 * (1.0 - refSW0) / (1.0 + 2.25*refSW0)
        refSW  = 0.33*bb/u2/Kd
        err = abs((refSW-refSW0) / refSW)

        if (err .lt. 0.0001) return

        refSW0 = refSW
    end do

end subroutine computeRefMorel

subroutine computeRefFresnel(&
    nr, ni, cosChi, sinChi, &
    coefFR)

    ! Description:
    ! Compute the Fresnel's coefficient of reflection
    !
    ! Input:
    !   nr: index of refraction of the sea water
    !   ni: extinction coefficient of the sea water
    !   cosChi: cosine of the incident radiation with respect of the wave facet normal
    !   sinChi: sine of the incident radiation with respect of the wave facet normal
    !
    ! Output:
    !   coefFR: coefficient for Fresnel's reflection
    !
    ! Reference:
    ! M. Born and E. Wolf, Principles of Optics, Pergamon Press, fifth
    !   edition, 1975, pp 628

    real, intent(in)  :: nr, ni, cosChi, sinChi
    real, intent(out) :: coefFR

    real :: a1, a2, u, v, Rr2, Rl2, b1, b2
    real :: nr2, ni2

    ! original code before optimization
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !a1 = abs(nr*nr - ni*ni - sinChi*sinChi)
    !a2 = sqrt((nr*nr - ni*ni - sinChi*sinChi)**2 + 4.0d0*nr*nr*ni*ni)
    !u = sqrt(0.5d0*(a1 + a2))
    !v = sqrt(max(0.0, 0.5d0*(-a1 + a2)))
    !Rr2 = ((cosChi - u)**2 + v*v)/((cosChi + u)**2 + v*v)
    !b1 = (nr*nr - ni*ni)*cosChi
    !b2 = 2.0*nr*ni*cosChi
    !Rl2 = ((b1 - u)**2 + (b2 - v)**2)/((b1 + u)**2 + (b2 + v)**2)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    nr2=nr*nr
    ni2=ni*ni

    ! absolute value for a1 to get v=0 when ni=0
    a1 = abs(nr2 - ni2 - sinChi*sinChi)
    a2 = sqrt(a1*a1 + 4.0*nr2*ni2)

    u = sqrt(0.5*( a1 + a2))
    v = sqrt(0.5*(-a1 + a2))
    Rr2 = ((cosChi-u)*(cosChi-u) + v*v)/((cosChi+u)*(cosChi+u) + v*v)

    b1 = (nr2-ni2)*cosChi
    b2 = 2.0*nr*ni*cosChi
    Rl2 = ((b1-u)*(b1-u)+(b2-v)*(b2-v)) / ((b1+u)*(b1+u) + (b2+v)*(b2+v))

    coefFR = (Rr2 + Rl2) * 0.5d0

    return

end subroutine computeRefFresnel

subroutine getRefracIndexWater(waveLen, xsal, nr, ni)

    ! Description:
    !   Retrieve indices of refraction for pure water
    !
    ! Input:
    !   waveLen: wavelength (in micron)
    !   xsal: salinity (in ppt), if xsal<0 then 34.3ppt by default
    !
    ! Output:
    !   nr: index of refraction of sea water
    !   ni: extinction coefficient of sea water
    !
    ! Reference:
    ! Hale, G. M., and Querry, M. R.: Optical Constants of Water in the 200-nm to
    !   200-μm Wavelength Region, Appl. Opt. 12, 555-563,
    !   https://doi.org/10.1364/AO.12.000555, 1973.

    real, intent(in)  :: waveLen, xsal
    real, intent(out) :: nr, ni

    real :: xwl, yr, yi, nrc, nic
    integer :: i

    real :: twl(62), tnr(62), tni(62)

    data twl/&
        0.250, 0.275, 0.300, 0.325, 0.345, 0.375, 0.400, 0.425, 0.445, 0.475,&
        0.500, 0.525, 0.550, 0.575, 0.600, 0.625, 0.650, 0.675, 0.700, 0.725,&
        0.750, 0.775, 0.800, 0.825, 0.850, 0.875, 0.900, 0.925, 0.950, 0.975,&
        1.000, 1.200, 1.400, 1.600, 1.800, 2.000, 2.200, 2.400, 2.600, 2.650,&
        2.700, 2.750, 2.800, 2.850, 2.900, 2.950, 3.000, 3.050, 3.100, 3.150,&
        3.200, 3.250, 3.300, 3.350, 3.400, 3.450, 3.500, 3.600, 3.700, 3.800,&
        3.900, 4.000/
    data tnr/&
        1.362, 1.354, 1.349, 1.346, 1.343, 1.341, 1.339, 1.338, 1.337, 1.336,&
        1.335, 1.334, 1.333, 1.333, 1.332, 1.332, 1.331, 1.331, 1.331, 1.330,&
        1.330, 1.330, 1.329, 1.329, 1.329, 1.328, 1.328, 1.328, 1.327, 1.327,&
        1.327, 1.324, 1.321, 1.317, 1.312, 1.306, 1.296, 1.279, 1.242, 1.219,&
        1.188, 1.157, 1.142, 1.149, 1.201, 1.292, 1.371, 1.426, 1.467, 1.483,&
        1.478, 1.467, 1.450, 1.432, 1.420, 1.410, 1.400, 1.385, 1.374, 1.364,&
        1.357, 1.351/
    data tni/&
        3.35e-08, 2.35e-08, 1.60e-08, 1.08e-08, 6.50e-09,&
        3.50e-09, 1.86e-09, 1.30e-09, 1.02e-09, 9.35e-10,&
        1.00e-09, 1.32e-09, 1.96e-09, 3.60e-09, 1.09e-08,&
        1.39e-08, 1.64e-08, 2.23e-08, 3.35e-08, 9.15e-08,&
        1.56e-07, 1.48e-07, 1.25e-07, 1.82e-07, 2.93e-07,&
        3.91e-07, 4.86e-07, 1.06e-06, 2.93e-06, 3.48e-06,&
        2.89e-06, 9.89e-06, 1.38e-04, 8.55e-05, 1.15e-04,&
        1.10e-03, 2.89e-04, 9.56e-04, 3.17e-03, 6.70e-03,&
        1.90e-02, 5.90e-02, 1.15e-01, 1.85e-01, 2.68e-01,&
        2.98e-01, 2.72e-01, 2.40e-01, 1.92e-01, 1.35e-01,&
        9.24e-02, 6.10e-02, 3.68e-02, 2.61e-02, 1.95e-02,&
        1.32e-02, 9.40e-03, 5.15e-03, 3.60e-03, 3.40e-03,&
        3.80e-03, 4.60e-03/

    i = 2

    ! origianl code (uses goto) replaced by do-while
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! 10  if (waveLen .lt. twl(i)) goto 20
    !     if (i .lt. 62) then
    !         i = i + 1
    !         goto 10
    !     end if
    ! 20  xwl = twl(i) - twl(i - 1)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    do while ((waveLen .ge. twl(i)) .and. (i .lt. 62))
        i = i + 1
    end do

    xwl = twl(i) - twl(i-1)

    yr = tnr(i) - tnr(i-1)
    yi = tni(i) - tni(i-1)
    nr = tnr(i-1) + (waveLen - twl(i-1))*yr/xwl
    ni = tni(i-1) + (waveLen - twl(i-1))*yi/xwl

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

end subroutine getRefracIndexWater

end module shd_sfc
