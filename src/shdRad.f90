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
!                hong.chen-1@colorado.edu
!╰────────────────────────────────────────────────────────────────────────────╯!

module shd_rad

    use globals_const, only: RD_, PI_, DTOR_, &
                             FRAC12_, FRAC13_
    use globals_cfg,   only: verbose_, accSolv_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, &
                             is3Ddir1Ddif_, &
                             isPol_, isPol4_, &
                             isSolar_, isThermal_
    use globals_type,  only: shd_type_prp, shd_type_sfc, &
                             shd_type_ckd, shd_type_sen

    use shd_scat, only: shd_scat__computePhase

    use shd_shdo, only: shd_shdo__computeYlm

    use shd_src,  only: shd_src__computeSrcOneDir, &
                        shd_src__integrateSrc, &
                        shd_src__integrateSrcOneRay

    use shd_init, only: shd_init__makeTopRad, &
                        shd_init__makeBotRadUL, &
                        shd_init__makeBotRadVL, &
                        shd_init__makeBotRadV

    implicit none
    private

    public :: shd_rad__makeDirect
    public :: shd_rad__computeDirect

    public :: shd_rad__interpRad
    public :: shd_rad__computeRad
    public :: shd_rad__computeCamRad

    public :: shd_rad__computeNetFluxDiv

contains

subroutine shd_rad__makeDirect(&
    &NPTS, BCFLAG, IPFLAG, DELTAM,&
    &ML, SOLARFLUX, SOLARMU, SOLARAZ,&
    &GRIDPOS, DIRFLUX, prp_obj, ckd_obj)

    ! Makes the direct beam solar flux for the internal base grid.
    !   DIRFLUX is set to F*exp(-tau_sun).
    !   Actually calls shd_rad__computeDirect to do all the hard work.

    type(shd_type_prp), intent(in) :: prp_obj
    type(shd_type_ckd), intent(in) :: ckd_obj

    integer NPTS, BCFLAG, IPFLAG, ML
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), DIRFLUX(NPTS)
    integer SIDE, IP
    logical VALIDBEAM
    real UNIFZLEV, XO, YO, ZO, DIR, DIRPATH

    call shd_rad__computeDirect(&
        &1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
        &DELTAM, ML,&
        &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, prp_obj, ckd_obj)

    do IP = 1, NPTS
        DIRPATH = 0.0
        call shd_rad__computeDirect&
            &(0, GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &BCFLAG, IPFLAG, DELTAM, ML,&
            &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(IP),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, prp_obj, ckd_obj)
    end do

    return
end subroutine shd_rad__makeDirect

subroutine shd_rad__computeDirect(&
    &INIT, XI, YI, ZI, BCFLAG, IPFLAG,&
    &DELTAM, ML, SOLARFLUX, SOLARMU, SOLARAZ,&
    &DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM, prp_obj, ckd_obj)

    !╭──────────────────────────────────────────────────────────────╮!
    ! Computes the direct beam flux at point (XI,YI,ZI) by integrating
    !   the extinction through the property grid.  If called with
    !   INIT=1 then the property grid extinction array, solar direction
    !   terms, and lowest uniform level are computed and saved.  If called
    !   with INIT=2 then the input lowest uniform level (UNIFZLEV) is stored.
    !   Call with INIT=0 to do the path integration and return the direct
    !   beam flux on the horizontal (DIRFLUX) at XI,YI,ZI.
    !   The DELTAM flag is used to determine whether to delta-M scale the
    !   extinction for the direct beam computation.  The extinction includes
    !   the gaseous absorption.  If the IPFLAG has bit 2 set then the direct
    !   beam tracing is done in 3D, otherwise the lower two bits determine
    !   the type of tracing: 0 for 3D tracing, 1 for XZ only tracing,
    !   2 for YZ only tracing, and 3 for Z only tracing.  If BCFLAG bits 0
    !   or 1 are set then have open boundary conditions in X and/or Y.  In
    !   this case when the ray tracing back to the sun reaches the boundary
    !   then independent pixel mode is entered so that only Z grid
    !   intersections occur.
    ! For use with multiple processors (bit 2 or 3 set in BCFLAG),
    !   the VALIDBEAM flag is returned true if the ray to the sun made it
    !   to the top of the domain before hitting the subdomain side.  If the
    !   flag is false SIDE is returned with the boundary hit (1=-X, 2=+X,
    !   3=-Y, 4=+Y).  XE,YE,ZE returns the location of the exitting ray,
    !   and path is the optical path from XI,YI,ZI to the sun.
    !╰──────────────────────────────────────────────────────────────╯!

    type(shd_type_prp), intent(in) :: prp_obj
    type(shd_type_ckd), intent(in) :: ckd_obj

    integer INIT, BCFLAG, IPFLAG, ML, SIDE
    logical DELTAM, VALIDBEAM
    real XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
    real DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH

    integer IX, IY, IZ, JZ, IL, IM, IU
    integer I, J, K, L, IPH, IP, JP, I1, I2, I3, I4
    integer IPDIRECT, DI, DJ, DK
    logical CONSTX, CONSTY, HITBOUNDARY
    real(RD_) :: EXTINCT, ALBEDO, F
    real(RD_) :: SUNMU, SUNAZ, PATH
    real(RD_) :: EXTBEAMCUT, UNIFORMZLEV
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: EPSS, EPSZ, XDOMAIN, YDOMAIN
    real(RD_) :: XOFFS, YOFFS, DELXD, DELYD
    real(RD_) :: X, Y, Z, XE, YE, ZE, XP, YP, ZP
    real(RD_) :: X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
    real(RD_) :: U0, V0, W0, U1, V1, W1, AX, AY, AZ
    real(RD_) :: U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
    real(RD_) :: E1, E2, E3, E4, E5, E6, E7, E8, A, B, C, D
    real(RD_) :: B1, B2, B3, B4, B5, B6, B7, B8, C1, C2, C3, C4, C5, C6, C7, C8
    real(RD_) :: UV, UMV, UVM, UMVM, UW, UMW, UWM, UMWM, VW, VMW, VWM, VMWM
    real(RD_) :: VWU, VWUM, UWV, UWVM, UVW, UVWM
    real, save, allocatable :: GASEXT(:), EXTMIN(:), EXTMAX(:)

    save CX, CY, CZ, CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT
    save DELXD, DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV

    if (INIT .eq. 9) then
        deallocate (GASEXT, EXTMIN, EXTMAX)
        return
    end if

    if (INIT .eq. 1) then
        if (.not. allocated(GASEXT)) allocate (GASEXT(prp_obj%Nz), EXTMIN(prp_obj%Nz), EXTMAX(prp_obj%Nz))
        ! Get the gaseous extinction at the property grid levels
        do IZ = 1, prp_obj%Nz
            if (ckd_obj%Nz .gt. 0) then
                IL = 1
                IU = ckd_obj%Nz
                do while (IU - IL .gt. 1)
                    IM = (IU + IL)/2
                    if (prp_obj%z(IZ) .le. ckd_obj%z(IM)) then
                        IL = IM
                    else
                        IU = IM
                    end if
                end do
                I = min(max(IL, 1), ckd_obj%Nz - 1)
                W0 = (prp_obj%z(IZ) - ckd_obj%z(I))/(ckd_obj%z(I + 1) - ckd_obj%z(I))
                W0 = min(max(W0, 0.0d0), 1.0d0)
                GASEXT(IZ) = (1.0 - W0)*ckd_obj%kabs(I, ckd_obj%ig) + W0*ckd_obj%kabs(I + 1, ckd_obj%ig)
            else
                GASEXT(IZ) = 0.0
            end if
        end do
        ! First make the property grid extinction field, on which all
        !   the direct beam paths will be computed.
        IP = 0
        do IX = 1, prp_obj%Nx
            do IY = 1, prp_obj%Ny
                do IZ = 1, prp_obj%Nz
                    IP = IP + 1
                    EXTINCT = prp_obj%ext(IP)
                    ALBEDO = prp_obj%ssa(IP)
                    ! Add in the gaseous absorption to extinction
                    if (GASEXT(IZ) .gt. 0.0) then
                        ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                        EXTINCT = EXTINCT + GASEXT(IZ)
                    end if
                    ! Do the Delta-M scaling if needed
                    if (DELTAM) then
                        if (prp_obj%Npha .gt. 0) then
                            IPH = prp_obj%ipha(IP)
                        else
                            IPH = IP
                        end if
                        L = ML + 1
                        !F = prp_obj%leg(1 + prp_obj%NlegSet*(L + (prp_obj%Nleg+ 1)*(IPH - 1)))/(2*L + 1)
                        F = prp_obj%leg(1, L, IPH)/(2*L + 1)
                        EXTINCT = (1.0 - ALBEDO*F)*EXTINCT
                    end if
                    prp_obj%extDir(IP) = EXTINCT
                end do
            end do
        end do

        ! Bit 2 of IPFLAG means do the direct beam in 3D
        IPDIRECT = IPFLAG
        if (is3Ddir1Ddif_) IPDIRECT = 0

        ! Make the ray direction (opposite to the solar beam)
        !   SOLAR?? is direction of beam, SUN?? is direction to sun (opposite)
        SUNMU = -SOLARMU
        SUNAZ = SOLARAZ + acos(-1.0)
        CX = sqrt(1.0 - SUNMU**2)*cos(SUNAZ)
        CY = sqrt(1.0 - SUNMU**2)*sin(SUNAZ)
        CZ = abs(SUNMU)
        if (abs(CX) .gt. 1.0e-6) then
            CXINV = 1.0d0/CX
        else
            CX = 0.0
            CXINV = 1.0e20
        end if
        if (abs(CY) .gt. 1.0e-6) then
            CYINV = 1.0d0/CY
        else
            CY = 0.0
            CYINV = 1.0e20
        end if
        if (abs(CZ) .gt. 1.0e-6) then
            CZINV = 1.0d0/CZ
        else
            CZ = 0.0
            CZINV = 1.0e20
        end if
        DI = nint(sign(1.0d0, CX))
        DJ = nint(sign(1.0d0, CY))
        DK = nint(sign(1.0d0, CZ))
        EPSZ = 1.0e-6*(prp_obj%z(prp_obj%Nz) - prp_obj%z(1))
        EPSS = 1.0e-3*(prp_obj%z(prp_obj%Nz) - prp_obj%z(1))/prp_obj%Nz
        if (.not. btest(IPDIRECT, 0)) EPSS = max(EPSS, 1.0d-4*prp_obj%dx)
        if (.not. btest(IPDIRECT, 1)) EPSS = max(EPSS, 1.0d-4*prp_obj%dy)
        DELXD = dble(prp_obj%dx)
        DELYD = dble(prp_obj%dy)
        EPSS = max(0.001*DELXD, 0.001*DELYD, EPSS)

        if (isMpiX_ .and. (prp_obj%Nx .gt. 1)) then
            XDOMAIN = DELXD*(prp_obj%Nx - 1)
        else
            XDOMAIN = DELXD*prp_obj%Nx
        end if

        if (isMpiY_ .and. (prp_obj%Ny .gt. 1)) then
            YDOMAIN = DELYD*(prp_obj%Ny - 1)
        else
            YDOMAIN = DELYD*prp_obj%Ny
        end if

        ! Find the Z level above which the medium is plane-parallel
        EXTBEAMCUT = 1.0e-4
        do IZ = 1, prp_obj%Nz
            EXTMIN(IZ) = 1.0e20
            EXTMAX(IZ) = 0.0
        end do
        IP = 0
        do IX = 1, prp_obj%Nx
            do IY = 1, prp_obj%Ny
                do IZ = 1, prp_obj%Nz
                    IP = IP + 1
                    EXTMIN(IZ) = min(prp_obj%ext(IP), EXTMIN(IZ))
                    EXTMAX(IZ) = max(prp_obj%ext(IP), EXTMAX(IZ))
                end do
            end do
        end do
        JZ = 0
        do IZ = 1, prp_obj%Nz
            if (EXTMAX(IZ) - EXTMIN(IZ) .gt. EXTBEAMCUT) JZ = IZ
        end do
        JZ = min(prp_obj%Nz, JZ + 1)
        UNIFORMZLEV = prp_obj%z(JZ)
        UNIFZLEV = UNIFORMZLEV
        return
    ! Done with initialization
    end if

    if (INIT .eq. 2) then
        UNIFORMZLEV = UNIFZLEV
        return
    end if

    ! Here for computing the direct beam path for one starting point.
    Z = ZI
    X = XI - prp_obj%xStart
    Y = YI - prp_obj%yStart

    ! Find the grid location
    IL = 0
    IU = prp_obj%Nz
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. prp_obj%z(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    K = max(IL, 1)

    I = int(X/DELXD) + 1
    if (I .gt. prp_obj%Nx .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) I = prp_obj%Nx
    if (I .lt. 1 .or. I .gt. prp_obj%Nx) then
        write (6, *) '_computeDirect: Beyond X domain', I, XI, YI, ZI
        stop
    end if

    J = int(Y/DELYD) + 1
    if (J .gt. prp_obj%Ny .and. abs(Y - YDOMAIN) .lt. 0.001*DELYD) J = prp_obj%Ny
    if (J .lt. 1 .or. J .gt. prp_obj%Ny) then
        write (6, *) '_computeDirect: Beyond Y domain', J, XI, YI, ZI
        stop
    end if
    XE = X
    YE = Y
    ZE = Z
    XP = XE
    YP = YE
    ZP = ZE
    CONSTX = btest(IPDIRECT, 0)
    CONSTY = btest(IPDIRECT, 1)
    if (CX .eq. 0.0) CONSTX = .true.
    if (CY .eq. 0.0) CONSTY = .true.
    if (isOpenX_ .and. (abs(X) .lt. 0.01*DELXD&
        &.or. abs(X - (prp_obj%Nx - 1)*DELXD) .lt. 0.01*DELXD)) CONSTX = .true.
    if (isOpenYx_ .and. (abs(Y) .lt. 0.01*DELYD&
        &.or. abs(Y - (prp_obj%Ny - 1)*DELYD) .lt. 0.01*DELYD)) CONSTY = .true.

    ! If have multiple subdomains (processors) and ray is going outwards
    !   from a boundary then set the HITBOUNDARY flag
    HITBOUNDARY = .false.
    if (isMpiX_) then
        if (CX .gt. 0.0 .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) then
            SIDE = 2
            HITBOUNDARY = .true.
        else if (CX .lt. 0.0 .and. abs(X) .lt. 0.001*DELXD) then
            SIDE = 1
            HITBOUNDARY = .true.
        end if
    end if
    if (isMpiY_) then
        if (CY .gt. 0.0 .and. abs(Y - YDOMAIN) .lt. 0.001*DELYD) then
            SIDE = 4
            HITBOUNDARY = .true.
        end if
        if (CY .lt. 0.0 .and. abs(Y) .lt. 0.001*DELYD) then
            SIDE = 3
            HITBOUNDARY = .true.
        end if
    end if

    ! Grid cell loop begin
    PATH = DIRPATH
    do while (.not. HITBOUNDARY .and. abs(ZE - prp_obj%z(prp_obj%Nz)) .gt. EPSZ)
        IP = I + 1
        if (I .eq. prp_obj%Nx) then
            if (isOpenX_ .or. isMpiX_) then
                IP = prp_obj%Nx
            else
                IP = 1
            end if
        end if
        JP = J + 1
        if (J .eq. prp_obj%Ny) then
            if (isOpenYx_ .or. isMpiY_) then
                JP = prp_obj%Ny
            else
                JP = 1
            end if
        end if
        X0 = DELXD*(I - 1)
        X1 = X0 + DELXD
        Y0 = DELYD*(J - 1)
        Y1 = Y0 + DELYD
        Z0 = prp_obj%z(K)
        Z1 = prp_obj%z(K + 1)
        if  (I .lt. 1 .or. I .gt. prp_obj%Nx .or.&
            &J .lt. 1 .or. J .gt. prp_obj%Ny .or.&
            &K .lt. 1 .or. K .ge. prp_obj%Nz) then
            write (6, '(A,3I4)') '_computeDirect: beyond grid!', I, J, K
            write (6, '(1(2X,3F9.5))') X, Y, Z
            stop
        end if
        ! Get the eight corner extinction values
        I1 = K + prp_obj%Nz*(J - 1) + prp_obj%Nz*prp_obj%Ny*(I - 1)
        I2 = K + prp_obj%Nz*(J - 1) + prp_obj%Nz*prp_obj%Ny*(IP - 1)
        I3 = K + prp_obj%Nz*(JP - 1) + prp_obj%Nz*prp_obj%Ny*(I - 1)
        I4 = K + prp_obj%Nz*(JP - 1) + prp_obj%Nz*prp_obj%Ny*(IP - 1)

        ! zgrid correction (to address layer averaging issue in z direction)
        !╭────────────────────────────────────────────────╮!
        !E1 = prp_obj%extDir(I1)
        !E2 = prp_obj%extDir(I2)
        !E3 = prp_obj%extDir(I3)
        !E4 = prp_obj%extDir(I4)
        E1 = prp_obj%extDir(I1 + 1)
        E2 = prp_obj%extDir(I2 + 1)
        E3 = prp_obj%extDir(I3 + 1)
        E4 = prp_obj%extDir(I4 + 1)
        E5 = prp_obj%extDir(I1 + 1)
        E6 = prp_obj%extDir(I2 + 1)
        E7 = prp_obj%extDir(I3 + 1)
        E8 = prp_obj%extDir(I4 + 1)
        !╰────────────────────────────────────────────────╯!

        ! Compute the distance to the next grid plane in  X, Y, and Z
        !   If in horizontal uniform region or doing IP then fix X and/or Y.
        if (ZE .ge. UNIFORMZLEV) then
            CONSTX = .true.
            CONSTY = .true.
        end if
        if (CONSTX) then
            SOX = 1.0e30
        else if (CX .gt. 0.0) then
            SOX = (X1 - XE)*CXINV
            XP = X1
        else
            SOX = (X0 - XE)*CXINV
            XP = X0
        end if
        if (CONSTY) then
            SOY = 1.0e30
        else if (CY .gt. 0.0) then
            SOY = (Y1 - YE)*CYINV
            YP = Y1
        else
            SOY = (Y0 - YE)*CYINV
            YP = Y0
        end if
        if (CZ .gt. 0.0) then
            SOZ = (Z1 - ZE)*CZINV
            ZP = Z1
        else if (CZ .lt. 0.0) then
            SOZ = (Z0 - ZE)*CZINV
            ZP = Z0
        else
            SOZ = 1.0e30
        end if

        ! The shortest distance is the plane we stop at:
        !   get the exitting location and increment the cell
        XOFFS = 0.0
        YOFFS = 0.0
        if (SOZ .le. SOX .and. SOZ .le. SOY) then
            SO = SOZ
            if (.not. CONSTX) XP = XE + SO*CX
            if (.not. CONSTY) YP = YE + SO*CY
            K = K + DK
        else if (SOX .le. SOY) then
            SO = SOX
            if (.not. CONSTY) YP = YE + SO*CY
            ZP = ZE + SO*CZ
            I = I + DI
            ! If have reached a horizontal boundary then either wrap around
            !   (periodic) or go into IP mode (open boundaries).
            if (I .eq. 0) then
                if (isOpenX_) then
                    I = 1
                    CONSTX = .true.
                else if (isMpiX_) then
                    HITBOUNDARY = .true.
                    SIDE = 1
                else
                    I = prp_obj%Nx
                    XOFFS = XDOMAIN
                end if
            else if (I .ge. prp_obj%Nx .and. isMpiX_) then
                HITBOUNDARY = .true.
                SIDE = 2
            else if (I .eq. prp_obj%Nx + 1) then
                if (isOpenX_) then
                    I = prp_obj%Nx
                    CONSTX = .true.
                else
                    I = 1
                    XOFFS = -XDOMAIN
                end if
            end if
        else
            SO = SOY
            if (.not. CONSTX) XP = XE + SO*CX
            ZP = ZE + SO*CZ
            J = J + DJ
            if (J .eq. 0) then
                if (isOpenYx_) then
                    J = 1
                    CONSTY = .true.
                else if (isMpiY_) then
                    HITBOUNDARY = .true.
                    SIDE = 3
                else
                    J = prp_obj%Ny
                    YOFFS = YDOMAIN
                end if
            else if (J .ge. prp_obj%Ny .and. isMpiY_) then
                HITBOUNDARY = .true.
                SIDE = 4
            else if (J .eq. prp_obj%Ny + 1) then
                if (isOpenYx_) then
                    J = prp_obj%Ny
                    CONSTY = .true.
                else
                    J = 1
                    YOFFS = -YDOMAIN
                end if
            end if
        end if
        if (SO .lt. -EPSS) then
            write (6, *) '_computeDirect: SO<0', X, Y, Z,&
                &XE, YE, ZE, XP, YP, ZP, CX, CY, CZ, SOX, SOY, SOZ
            stop
        end if
        SO = max(SO, 0.0d0)
        ! Make the starting and ending interpolation factors
        AX = 1.0d0/(X1 - X0)
        AY = 1.0d0/(Y1 - Y0)
        AZ = 1.0d0/(Z1 - Z0)
        U0 = (XE - X0)*AX
        V0 = (YE - Y0)*AY
        W0 = (ZE - Z0)*AZ
        U1 = (XP - X0)*AX
        V1 = (YP - Y0)*AY
        W1 = (ZP - Z0)*AZ
        ! Compute the cubic polynomial extinction coefficients
        U0M = 1.0 - U0
        V0M = 1.0 - V0
        W0M = 1.0 - W0
        U1M = 1.0 - U1
        V1M = 1.0 - V1
        W1M = 1.0 - W1
        DU = U1 - U0
        DV = V1 - V0
        DW = W1 - W0
        UV = U0*V0
        UMV = U0M*V0
        UVM = U0*V0M
        UMVM = U0M*V0M
        UW = U0*W0
        UMW = U0M*W0
        UWM = U0*W0M
        UMWM = U0M*W0M
        VW = V0*W0
        VMW = V0M*W0
        VWM = V0*W0M
        VMWM = V0M*W0M
        A =    (E1*U0M + E2*U0)*VMWM&
            &+ (E3*U0M + E4*U0)*VWM&
            &+ (E5*U0M + E6*U0)*VMW&
            &+ (E7*U0M + E8*U0)*VW
        B1 = -DU*VMWM - DV*UMWM - DW*UMVM
        B2 = DU*VMWM - DV*UWM - DW*UVM
        B3 = -DU*VWM + DV*UMWM - DW*UMV
        B4 = DU*VWM + DV*UWM - DW*UV
        B5 = -DU*VMW - DV*UMW + DW*UMVM
        B6 = DU*VMW - DV*UW + DW*UVM
        B7 = -DU*VW + DV*UMW + DW*UMV
        B8 = DU*VW + DV*UW + DW*UV
        B = B1*E1 + B2*E2 + B3*E3 + B4*E4 + B5*E5 + B6*E6 + B7*E7 + B8*E8
        VW = DV*DW
        VWU = VW*U0
        VWUM = VW*U0M
        UW = DU*DW
        UWV = UW*V0
        UWVM = UW*V0M
        UV = DU*DV
        UVW = UV*W0
        UVWM = UV*W0M
        C1 = +VWUM + UWVM + UVWM
        C2 = +VWU - UWVM - UVWM
        C3 = -VWUM + UWV - UVWM
        C4 = -VWU - UWV + UVWM
        C5 = -VWUM - UWVM + UVW
        C6 = -VWU + UWVM - UVW
        C7 = +VWUM - UWV - UVW
        C8 = +VWU + UWV + UVW
        C = C1*E1 + C2*E2 + C3*E3 + C4*E4 + C5*E5 + C6*E6 + C7*E7 + C8*E8
        D = DU*DV*DW*(E2 + E3 + E5 + E8 - E1 - E4 - E6 - E7)
        ! Compute the path through the cell: integration of extinction
        PATH = PATH + SO*(A + 0.5d0*B + FRAC13_*C + 0.25d0*D)

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP
    end do
    DIRFLUX = SOLARFLUX*exp(-PATH)
    DIRPATH = PATH
    XO = XE + prp_obj%xStart; YO = YE + prp_obj%yStart; ZO = ZE
    VALIDBEAM = abs(ZE - prp_obj%z(prp_obj%Nz)) .lt. EPSZ
    if (VALIDBEAM) SIDE = 6

    return

end subroutine shd_rad__computeDirect

subroutine shd_rad__interpRad(&
    &NSTOKES, OLDNPTS, NPTS, RSHPTR,&
    &RADIANCE, NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)

    ! Interpolates the spherical harmonic radiance array to the new
    !   grid points (>OLDNPTS) from the old grid points.

    integer NSTOKES, OLDNPTS, NPTS, NBCELLS, NCELLS
    integer RSHPTR(NPTS + 1), GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    real RADIANCE(NSTOKES, *), GRIDPOS(3, *)
    integer IP, IC, ICELL, IPARENT, IDIR, I, J, K, DIR1, DIR2, IP1, IP2
    integer IR, IR1, IR2, NR, NR1, NR2
    integer GRIDCORNER(2, 4, 3)
    real RAD1(NSTOKES), RAD2(NSTOKES)
    data GRIDCORNER/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        &1, 5, 2, 6, 3, 7, 4, 8/

    if (OLDNPTS .le. 0) return

    ! Go through all the non-base grid cells
    do ICELL = NBCELLS + 1, NCELLS
        do IC = 1, 8
            IP = GRIDPTR(IC, ICELL)
            if (IP .gt. OLDNPTS) then
                IPARENT = TREEPTR(1, ICELL)
                ! Go through all 12 edges of parent cell to find which one h
                do IDIR = 1, 3
                    DIR1 = mod(IDIR - 1 + 1, 3) + 1
                    DIR2 = mod(IDIR - 1 + 2, 3) + 1
                    do I = 1, 4
                        IP1 = GRIDPTR(GRIDCORNER(1, I, IDIR), IPARENT)
                        IP2 = GRIDPTR(GRIDCORNER(2, I, IDIR), IPARENT)
                        if (GRIDPOS(DIR1, IP1) .eq. GRIDPOS(DIR1, IP) .and.&
                            &GRIDPOS(DIR2, IP1) .eq. GRIDPOS(DIR2, IP)) goto 120
                    end do
                end do
                stop 'Error [_interpRad]: Point not on an edge.'
120             continue
                ! Interpolate the data from the two parent points at the
                !   ends of the edge; assume the cell was split in half.
                IR1 = RSHPTR(IP1)
                IR2 = RSHPTR(IP2)
                NR1 = RSHPTR(IP1 + 1) - RSHPTR(IP1)
                NR2 = RSHPTR(IP2 + 1) - RSHPTR(IP2)
                NR = max(NR1, NR2)
                IR = RSHPTR(IP)
                RSHPTR(IP + 1) = IR + NR
                do J = 1, NR
                    if (J .le. NR1) then
                        RAD1(:) = RADIANCE(:, IR1 + J)
                    else
                        RAD1(:) = 0.0
                    end if
                    if (J .le. NR2) then
                        RAD2(:) = RADIANCE(:, IR2 + J)
                    else
                        RAD2(:) = 0.0
                    end if
                    RADIANCE(:, IR + J) = FRAC12_*(RAD1(:) + RAD2(:))
                end do
            end if
        end do
    end do

    return
end subroutine shd_rad__interpRad

subroutine shd_rad__computeRad(&
    &NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &SKYRAD, WAVENO, WAVELEN, UNITS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
    &SHPTR, SOURCE, SOURCE1, GRIDRAD,&
    &OUTPARMS, IRAD, RADOUT, prp_obj, sfc_obj, sen_obj)

    ! Computes the radiances for the specified locations and directions

    type(shd_type_prp), intent(in) :: prp_obj
    type(shd_type_sfc), intent(in) :: sfc_obj
    type(shd_type_sen), intent(in) :: sen_obj

    integer NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer ML, MM, NLM
    integer NMU, NPHI0MAX, NPHI0(NMU), IRAD
    integer MAXNBC, NTOPPTS, NBOTPTS
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ
    real SKYRAD, WAVENO(2), WAVELEN
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real BCRAD(NSTOKES, *)
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(prp_obj%NlegSet, 0:prp_obj%Nleg, NPTS)
    real DIRFLUX(NPTS), FLUXES(2, NPTS)
    real SOURCE1(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    real SOURCE(NSTOKES, *)
    real OUTPARMS(*), RADOUT(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IBC, IANG, JX, JY, K, SIDE
    integer NANGOUT, NXOUT, NYOUT
    logical VALIDRAD
    real(RD_) :: X0, Y0, Z0, XE, YE, ZE, TRANSMIT, RADIANCE(4)
    real MUOUT, PHIOUT, PHID
    real XDOMAIN, YDOMAIN, STARTX, STARTY

    ! added for fixing the "domain shift" issue
    !╭──────────────────────────────────────────────────────────────╮!
    real :: DS_Z0, DS_VZA_RAD, DS_XOFFSET, DS_YOFFSET
    !╰──────────────────────────────────────────────────────────────╯!

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))

    ! Make the bottom boundary radiances for the Lambertian surfaces.
    ! Compute the upwelling bottom radiances using the downwelling fluxes.
    if (sfc_obj%type .eq. 'FL') then

        call shd_init__makeBotRadUL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS), sfc_obj)

    else if (sfc_obj%type .eq. 'VL') then

        call shd_init__makeBotRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE,&
            &NSTOKES, BCRAD(1, 1 + NTOPPTS), sfc_obj)

    end if

    ! Setup the regularly spaced radiance locations
    XDOMAIN = XGRID(NX + 1 - ibits(BCFLAG, 0, 1)) - XGRID(1) + 2*OUTPARMS(4)
    STARTX = -OUTPARMS(4)
    if (OUTPARMS(2) .eq. 0.0) then
        NXOUT = 1
    else
        NXOUT = max(1, nint(XDOMAIN/OUTPARMS(2)))
    end if
    YDOMAIN = YGRID(NY + 1 - ibits(BCFLAG, 1, 1)) - YGRID(1) + 2*OUTPARMS(5)
    STARTY = -OUTPARMS(5)
    if (OUTPARMS(3) .eq. 0.0) then
        NYOUT = 1
    else
        NYOUT = max(1, nint(YDOMAIN/OUTPARMS(3)))
    end if
    Z0 = min(max(OUTPARMS(1), ZGRID(1)), ZGRID(NZ))
    NANGOUT = nint(OUTPARMS(6))

    ! Loop over the radiance directions
    do IANG = 1, NANGOUT
        if (sen_obj%fname(1:4) .ne. 'NONE') then
            MUOUT = sen_obj%view(2*IANG - 1)
            PHID  = sen_obj%view(2*IANG)
        else
            MUOUT = OUTPARMS(2*IANG + 5)
            PHID  = OUTPARMS(2*IANG + 6)
        end if
        PHIOUT = PHID*acos(-1.0)/180.0

        ! added for fixing the "domain shift" issue (only for 3D)
        !╭────────────────────────────────────────────────╮!
        if (IPFLAG .eq. 0) then
            DS_Z0 = Z0 - ZGRID(1)
            DS_VZA_RAD = acos(MUOUT)
            DS_XOFFSET = DS_Z0*tan(DS_VZA_RAD)*cos(PHIOUT)
            DS_YOFFSET = DS_Z0*tan(DS_VZA_RAD)*sin(PHIOUT)
            STARTX = -OUTPARMS(4) + DS_XOFFSET
            STARTY = -OUTPARMS(5) + DS_YOFFSET
        endif
        !╰────────────────────────────────────────────────╯!

        if (MUOUT .eq. 0.0) then
            if (verbose_ .ge. 4) write (6, *) 'Warning [_computeRad]: Bad MUOUT for radiance', IANG, MUOUT, PHID
            MUOUT = accSolv_
            write (6, *) '                       New MUOUT=', MUOUT
        end if

        if (abs(MUOUT) .gt. 1.0) then
            if (verbose_ .ge. 4) write (6, *) 'Warning [_computeRad]: Bad MUOUT for radiance', IANG, MUOUT, PHID
            if (MUOUT .gt.  1.0) MUOUT =  1.0
            if (MUOUT .lt. -1.0) MUOUT = -1.0
            write (6, *) '                       New MUOUT=', MUOUT
        end if

        ! Compute the source function throughout grid for this angle
        call shd_src__computeSrcOneDir(NSTOKES, ML, MM, NLM,&
            &NPTS, DELTAM, MUOUT, PHIOUT, SRCTYPE,&
            &SOLARMU, SOLARAZ, ALBEDO, LEGEN,&
            &IPHASE, DIRFLUX, SHPTR, SOURCE, SOURCE1, prp_obj)

        ! Set the radiance field to -1, so we can determine valid ra
        !   Also set the source array to the extinction times the so
        !   function, so the grid interpolation is correct.
        do I = 1, NPTS
            GRIDRAD(1, I) = -1.0
            SOURCE1(:, I) = SOURCE1(:, I)*EXTINCT(I)
        end do
        ! Get boundary radiances: either top or bottom
        ! Isotropic top boundary or Lambertian bottom boundary can u
        !   the previously computed boundary radiances in BCRAD,
        !   otherwise, compute the radiance for this angle
        !   by integrating over the stored downwelling radiances.
        if (MUOUT .lt. 0.0) then
            do IBC = 1, NTOPPTS
                I = BCPTR(IBC, 1)
                GRIDRAD(:, I) = BCRAD(:, IBC)
            end do
        else
            if (.not. sfc_obj%isLambertian) then
                call shd_init__makeBotRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT,&
                    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                    &NSTOKES,&
                    &BCRAD(1, 1 + NTOPPTS), sfc_obj)
            end if
            do IBC = 1, NBOTPTS
                I = BCPTR(IBC, 2)
                GRIDRAD(:, I) = BCRAD(:, NTOPPTS + IBC)
            end do
        end if

        ! Integrate backward from the location to get the radiance
        Y0 = STARTY
        do JY = 1, NYOUT
            X0 = STARTX
            do JX = 1, NXOUT
                IRAD = IRAD + 1
                TRANSMIT = 1.0d0; RADIANCE(:) = 0.0d0
                call shd_src__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
                    &NX, NY, NZ, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &XGRID, YGRID, ZGRID, GRIDPOS,&
                    &MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1,&
                    &X0, Y0, Z0, XE, YE, ZE, SIDE,&
                    &TRANSMIT, RADIANCE, VALIDRAD)
                RADOUT(:, IRAD) = RADIANCE(1:NSTOKES)
                X0 = X0 + OUTPARMS(2)
            end do
            Y0 = Y0 + OUTPARMS(3)
        end do

    end do

    return

end subroutine shd_rad__computeRad

subroutine shd_rad__computeCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &SKYRAD, WAVENO, WAVELEN, UNITS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
    &SHPTR, SOURCE, OUTPARMS, IVIS, VISOUT, prp_obj, sfc_obj, sen_obj)

    !   Computes Stokes radiances (output in VISOUT) for the visualization
    !  modes: 1) camera mode, and 2) cross track scanning.

    type(shd_type_prp), intent(in) :: prp_obj
    type(shd_type_sfc), intent(in) :: sfc_obj
    type(shd_type_sen), intent(in) :: sen_obj

    integer NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer ML, MM, NLM
    integer NMU, NPHI0MAX, NPHI0(NMU), IVIS
    integer MAXNBC, NTOPPTS, NBOTPTS
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ
    real SKYRAD, WAVENO(2), WAVELEN
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real BCRAD(NSTOKES, *)
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(prp_obj%NlegSet, 0:prp_obj%Nleg, *)
    real DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(NSTOKES, *)
    real OUTPARMS(*), VISOUT(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer NSCATANGLE, NSTPHASE, I, J, L, N, SIDE
    integer NL, NS, LINE, SAMP
    logical CAMERA_MODE, VALIDRAD
    real MURAY, PHIRAY
    real(RD_) :: XE, YE, ZE, TRANSMIT, VISRAD(NSTOKES)
    real(RD_) :: X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
    real(RD_) :: THETA0, THETA1, PHIR, PHI0
    real(RD_) :: COSTH, SINTH, SINTH0, SINTH1, MU2, PHI2
    real(RD_) :: U, V, UP, VP, COSDELPHI, ROTANG, R
    real(RD_) :: DIST, D, RX, RY, RZ, SCANANG
    integer MAXSCATANG
    parameter(MAXSCATANG=721)
    real, allocatable :: YLMSUN(:, :), PHASETAB(:, :, :)
    real              :: MEAN, STD1, STD2
    real, allocatable :: AOLP(:)

    allocate (YLMSUN(prp_obj%NlegSet, NLM))

    if (isSolar_) then
        call shd_shdo__computeYlm(.true., SOLARMU, SOLARAZ, ML, MM, prp_obj%NlegSet, YLMSUN)
        if (DELTAM .and. prp_obj%Npha .gt. 0) then
            NSCATANGLE = max(36, min(MAXSCATANG, 2*prp_obj%Nleg))
            NSTPHASE = min(prp_obj%NlegSet, 2)
            allocate (PHASETAB(NSTPHASE, prp_obj%Npha, NSCATANGLE))
            call shd_scat__computePhase(NSCATANGLE, prp_obj%Npha, NSTPHASE,&
                &NSTOKES, ML, prp_obj%NlegSet, prp_obj%Nleg, LEGEN, PHASETAB)
        end if
    end if

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces
    !   Compute the upwelling bottom radiances using the downwelling
    if (sfc_obj%type .eq. 'FL') then
        call shd_init__makeBotRadUL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS), sfc_obj)
    else if (sfc_obj%type .eq. 'VL') then
        call shd_init__makeBotRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE,&
            &NSTOKES, BCRAD(1, 1 + NTOPPTS), sfc_obj)
    end if

    CAMERA_MODE = nint(OUTPARMS(1)) .eq. 1
    if (CAMERA_MODE) then
        NL = nint(OUTPARMS(10))
        NS = nint(OUTPARMS(11))
    else
        NL = 1 + int(sqrt((OUTPARMS(4) - OUTPARMS(7))**2&
            &+ (OUTPARMS(5) - OUTPARMS(8))**2&
            &+ (OUTPARMS(6) - OUTPARMS(9))**2)/OUTPARMS(10))
        NS = 1 + int(abs(OUTPARMS(12) - OUTPARMS(11))/OUTPARMS(13))
    end if

    ! Loop over pixels in image
    do LINE = 1, NL
        do SAMP = 1, NS
            IVIS = IVIS + 1

            if (CAMERA_MODE) then
                ! Camera mode:
                !   1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline,
                !
                ! Use spherical trig to find the pixel direction (MURAY,PHIRAY
                !   from the camera center (THETA0,PHI0) and the relative pixe
                !   angles (U,V).
                UP = (SAMP - NS/2 - 1)*OUTPARMS(13)*DTOR_
                VP = (LINE - NL/2 - 1)*OUTPARMS(12)*DTOR_
                ROTANG = OUTPARMS(9)*DTOR_
                if (ROTANG .eq. 0.0) then
                    U = UP
                    V = VP
                else
                    U = cos(ROTANG)*UP - sin(ROTANG)*VP
                    V = sin(ROTANG)*UP + cos(ROTANG)*VP
                end if
                THETA0 = DTOR_*OUTPARMS(7)
                PHI0 = DTOR_*OUTPARMS(8)
                THETA1 = THETA0 + V
                if (V .eq. 0.0) then
                    COSTH = cos(U)*cos(THETA0)
                else
                    COSTH = cos(U)*(sin(THETA1)*cos(V) - sin(THETA0))/sin(V)
                    COSTH = min(+1.0d0, max(-1.0d0, COSTH))
                end if
                SINTH = sqrt(1 - COSTH**2)
                SINTH0 = sin(THETA0)
                SINTH1 = sin(THETA1)
                if (abs(SINTH) .lt. accSolv_) then
                    PHIR = 0.0
                else
                    if ((abs(SINTH0) .lt. accSolv_) .and. (abs(SINTH1) .le. accSolv_)) then
                        COSDELPHI = 0.0d0
                    else if (abs(SINTH1) .gt. accSolv_) then
                        COSDELPHI = (cos(U) - COSTH*cos(THETA1))/(SINTH1*SINTH)
                    else if (abs(SINTH0) .gt. accSolv_) then
                        COSDELPHI = (cos(U)*cos(V) - COSTH*cos(THETA0))&
                            &/(SINTH0*SINTH)
                    end if
                    COSDELPHI = min(+1.0d0, max(-1.0d0, COSDELPHI))
                    if (U .ge. 0.0) then
                        PHIR = PHI0 - acos(COSDELPHI)
                    else
                        PHIR = PHI0 + acos(COSDELPHI)
                    end if
                end if

                X0 = OUTPARMS(4)
                Y0 = OUTPARMS(5)
                Z0 = OUTPARMS(6)

            else

                ! Cross track scanning in the vertical plane:
                !   2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, d
                ! Start and end of scan angle (scan1, scan2) are +/- relative
                !   to nadir, and positive is on right side.
                X1 = OUTPARMS(4)
                Y1 = OUTPARMS(5)
                Z1 = OUTPARMS(6)
                X2 = OUTPARMS(7)
                Y2 = OUTPARMS(8)
                Z2 = OUTPARMS(9)
                DIST = sqrt((X2 - X1)**2 + (Y2 - Y1)**2 + (Z2 - Z1)**2)
                RX = (X2 - X1)/DIST
                RY = (Y2 - Y1)/DIST
                RZ = (Z2 - Z1)/DIST
                ! D = (LINE-1)*OUTPARMS(10)
                D = (NL - LINE)*OUTPARMS(10)
                X0 = X1 + D*RX
                Y0 = Y1 + D*RY
                Z0 = Z1 + D*RZ
                SCANANG = DTOR_*(OUTPARMS(11) + (SAMP - 1)*OUTPARMS(13))
                COSTH = cos(PI_ - abs(SCANANG))
                if (SCANANG .gt. 0.0) then
                    PHIR = atan2(-RX, RY)
                else
                    PHIR = atan2(RX, -RY)
                end if
            end if

            ! Extrapolate ray to domain top if above
            if (Z0 .gt. ZGRID(NZ)) then
                if (COSTH .ge. 0.0) then
                    VISRAD(:) = 0.0
                    goto 900
                end if
                R = (ZGRID(NZ) - Z0)/COSTH
                X0 = X0 + R*sqrt(1 - COSTH**2)*cos(PHIR)
                Y0 = Y0 + R*sqrt(1 - COSTH**2)*sin(PHIR)
                Z0 = ZGRID(NZ)
            else if (Z0 .lt. ZGRID(1)) then
                write (6, *) 'Error [_visualizeRad]: Level below domain.'
                stop
            end if

            ! MURAY,PHRAY is camera pixel viewing direction;
            ! MU2,PHI2 is radiance travel direction (the opposite)
            MURAY = sngl(COSTH)
            PHIRAY = sngl(PHIR)
            MU2 = -COSTH
            PHI2 = PHIR + PI_

            ! Integrate the extinction and source function along this ray
            !   to calculate the Stokes radiance vector for this pixel
            TRANSMIT = 1.0d0; VISRAD(:) = 0.0d0
            call shd_src__integrateSrcOneRay(BCFLAG, IPFLAG, NSTOKES,&
                &NSTPHASE, NSCATANGLE, PHASETAB,&
                &NX, NY, NZ, NPTS, NCELLS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                &XGRID, YGRID, ZGRID, GRIDPOS,&
                &ML, MM, NLM,&
                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                &DELTAM, SRCTYPE, WAVELEN, SOLARMU, SOLARAZ,&
                &EXTINCT, ALBEDO, LEGEN, IPHASE,&
                &DIRFLUX, SHPTR, SOURCE, YLMSUN,&
                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                &MU2, PHI2, X0, Y0, Z0,&
                &XE, YE, ZE, SIDE, TRANSMIT, VISRAD, VALIDRAD,&
                &prp_obj, sfc_obj)
900         continue
            ! write (6,'(1X,2F8.4,1X,2F11.7,4(1X,F11.6))')&
            !     &X0,Y0,MURAY,PHIRAY,VISRAD(:)
            VISOUT(1, IVIS) = VISRAD(1)
            if (VISRAD(1) .gt. 0.0) then
                if (isPol_) then
                    ! Output degree (0 to 1) and angle (-180 to 180) of linear polarization
                    VISOUT(2, IVIS) = sqrt(VISRAD(2)**2 + VISRAD(3)**2)/VISRAD(1)
                    VISOUT(3, IVIS) = (180.0d0/PI_)*FRAC12_*atan2(VISRAD(3), VISRAD(2))
                end if
                if (isPol4_) then
                    ! Output degree of circular polarization (-1 to 1)
                    VISOUT(4, IVIS) = VISRAD(4)/VISRAD(1)
                end if
            else
                VISOUT(2:, IVIS) = 0.0
            end if
        end do
    end do

    deallocate (YLMSUN)
    if (allocated(PHASETAB)) deallocate (PHASETAB)
    if (isPol_) then
        ! Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
        N = NL*NS
        allocate (AOLP(N))
        AOLP(:) = VISOUT(3, IVIS - N + 1:IVIS)
        MEAN = sum(AOLP(:))/N
        STD1 = sqrt(sum((AOLP(:) - MEAN)**2)/N)
        where (AOLP(:) < 0.0)
            AOLP(:) = AOLP(:) + 180.0
        end where
        MEAN = sum(AOLP(:))/N
        STD2 = sqrt(sum((AOLP(:) - MEAN)**2)/N)
        if (STD2 < STD1) then
            VISOUT(3, IVIS - N + 1:IVIS) = AOLP(:)
        end if
        deallocate (AOLP)
    end if
    return
end subroutine shd_rad__computeCamRad

subroutine shd_rad__computeNetFluxDiv(NSTOKES, NPTS, RSHPTR, SRCTYPE,&
    &SOLARMU, EXTINCT, ALBEDO, PLANCK, DIRFLUX,&
    &RADIANCE, NETFLUXDIV)

    ! Computes the net flux divergence at every grid point.

    integer NSTOKES, NPTS, RSHPTR(NPTS + 1)
    real PLANCK(NPTS), DIRFLUX(NPTS), SOLARMU
    real EXTINCT(NPTS), ALBEDO(NPTS)
    real RADIANCE(NSTOKES, *), NETFLUXDIV(NPTS)
    character(len=1) :: SRCTYPE
    integer I, IR
    real :: C, SECMU, HR

    C = sqrt(4.0*PI_)
    SECMU = 1.0/abs(SOLARMU)
    do I = 1, NPTS
        IR = RSHPTR(I)
        HR = C*EXTINCT(I)*(1 - ALBEDO(I))*RADIANCE(1, IR + 1)
        if (isSolar_) then
            HR = HR + EXTINCT(I)*(1.0 - ALBEDO(I))*DIRFLUX(I)*SECMU
        end if
        if (isThermal_) then
            HR = HR - 4*PI_*EXTINCT(I)*PLANCK(I)
        end if
        NETFLUXDIV(I) = -HR
    end do

    return

end subroutine shd_rad__computeNetFluxDiv

end module shd_rad
