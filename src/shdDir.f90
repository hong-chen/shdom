module shd_dir

    use globals

    implicit none

    public :: shd_dir__makeDirect
    public :: shd_dir__computeDirect

contains

subroutine shd_dir__makeDirect(&
    &NPTS, BCFLAG, IPFLAG, DELTAM,&
    &ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
    &GRIDPOS, DIRFLUX)

    ! Makes the direct beam solar flux for the internal base grid.
    !   DIRFLUX is set to F*exp(-tau_sun).
    !   Actually calls shd_dir__computeDirect to do all the hard work.

    integer NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), DIRFLUX(NPTS)
    integer SIDE, IP
    logical VALIDBEAM
    real UNIFZLEV, XO, YO, ZO, DIR, DIRPATH

    call shd_dir__computeDirect(&
        &1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
        &DELTAM, ML, NSTLEG, NLEG,&
        &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)

    do IP = 1, NPTS
        DIRPATH = 0.0
        call shd_dir__computeDirect&
            &(0, GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &BCFLAG, IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
            &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(IP),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end do

    return
end subroutine shd_dir__makeDirect

subroutine shd_dir__computeDirect(&
    &INIT, XI, YI, ZI, BCFLAG, IPFLAG,&
    &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
    &DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)

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

    use globals_shdom

    integer INIT, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, SIDE
    logical DELTAM, VALIDBEAM
    real XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
    real DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH

    integer IX, IY, IZ, JZ, IL, IM, IU
    integer I, J, K, L, IPH, IP, JP, I1, I2, I3, I4
    integer IPDIRECT, DI, DJ, DK
    logical CONSTX, CONSTY, HITBOUNDARY
    double precision EXTINCT, ALBEDO, F
    double precision SUNMU, SUNAZ, PATH
    double precision EXTBEAMCUT, UNIFORMZLEV
    double precision CX, CY, CZ, CXINV, CYINV, CZINV
    double precision EPSS, EPSZ, XDOMAIN, YDOMAIN
    double precision XOFFS, YOFFS, DELXD, DELYD
    double precision X, Y, Z, XE, YE, ZE, XP, YP, ZP
    double precision X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
    double precision U0, V0, W0, U1, V1, W1, AX, AY, AZ
    double precision U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
    double precision E1, E2, E3, E4, E5, E6, E7, E8, A, B, C, D
    double precision B1, B2, B3, B4, B5, B6, B7, B8, C1, C2, C3, C4, C5, C6, C7, C8
    double precision UV, UMV, UVM, UMVM, UW, UMW, UWM, UMWM, VW, VMW, VWM, VMWM
    double precision VWU, VWUM, UWV, UWVM, UVW, UVWM
    real, save, allocatable :: GASEXT(:), EXTMIN(:), EXTMAX(:)

    save CX, CY, CZ, CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT
    save DELXD, DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV

    if (INIT .eq. 9) then
        deallocate (GASEXT, EXTMIN, EXTMAX)
        return
    end if

    if (INIT .eq. 1) then
        if (.not. allocated(GASEXT)) allocate (GASEXT(NPZ_), EXTMIN(NPZ_), EXTMAX(NPZ_))
        ! Get the gaseous extinction at the property grid levels
        do IZ = 1, NPZ_
            if (NZCKD_ .gt. 0) then
                IL = 1
                IU = NZCKD_
                do while (IU - IL .gt. 1)
                    IM = (IU + IL)/2
                    if (ZLEVELS_(IZ) .le. ZCKD_(IM)) then
                        IL = IM
                    else
                        IU = IM
                    end if
                end do
                I = min(max(IL, 1), NZCKD_ - 1)
                W0 = (ZLEVELS_(IZ) - ZCKD_(I))/(ZCKD_(I + 1) - ZCKD_(I))
                W0 = min(max(W0, 0.0d0), 1.0d0)
                GASEXT(IZ) = (1.0 - W0)*GASABS_(I) + W0*GASABS_(I + 1)
            else
                GASEXT(IZ) = 0.0
            end if
        end do
        ! First make the property grid extinction field, on which all
        !   the direct beam paths will be computed.
        IP = 0
        do IX = 1, NPX_
            do IY = 1, NPY_
                do IZ = 1, NPZ_
                    IP = IP + 1
                    EXTINCT = EXTINCTP_(IP)
                    ALBEDO = ALBEDOP_(IP)
                    ! Add in the gaseous absorption to extinction
                    if (GASEXT(IZ) .gt. 0.0) then
                        ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                        EXTINCT = EXTINCT + GASEXT(IZ)
                    end if
                    ! Do the Delta-M scaling if needed
                    if (DELTAM) then
                        if (NUMPHASE_ .gt. 0) then
                            IPH = IPHASEP_(IP)
                        else
                            IPH = IP
                        end if
                        L = ML + 1
                        F = LEGENP_(1 + NSTLEG*(L + (NLEG + 1)*(IPH - 1)))/(2*L + 1)
                        EXTINCT = (1.0 - ALBEDO*F)*EXTINCT
                    end if
                    EXTDIRP_(IP) = EXTINCT
                end do
            end do
        end do

        ! Bit 2 of IPFLAG means do the direct beam in 3D
        IPDIRECT = IPFLAG
        if (btest(IPFLAG, 2)) IPDIRECT = 0

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
        EPSZ = 1.0e-6*(ZLEVELS_(NPZ_) - ZLEVELS_(1))
        EPSS = 1.0e-3*(ZLEVELS_(NPZ_) - ZLEVELS_(1))/NPZ_
        if (.not. btest(IPDIRECT, 0)) EPSS = max(EPSS, 1.0d-4*DELX_)
        if (.not. btest(IPDIRECT, 1)) EPSS = max(EPSS, 1.0d-4*DELY_)
        DELXD = dble(DELX_)
        DELYD = dble(DELY_)
        EPSS = max(0.001*DELXD, 0.001*DELYD, EPSS)

        if (btest(BCFLAG, 2) .and. (NPX_ .gt. 1)) then
            XDOMAIN = DELXD*(NPX_ - 1)
        else
            XDOMAIN = DELXD*NPX_
        end if

        if (btest(BCFLAG, 3) .and. (NPY_ .gt. 1)) then
            YDOMAIN = DELYD*(NPY_ - 1)
        else
            YDOMAIN = DELYD*NPY_
        end if

        ! Find the Z level above which the medium is plane-parallel
        EXTBEAMCUT = 1.0e-4
        do IZ = 1, NPZ_
            EXTMIN(IZ) = 1.0e20
            EXTMAX(IZ) = 0.0
        end do
        IP = 0
        do IX = 1, NPX_
            do IY = 1, NPY_
                do IZ = 1, NPZ_
                    IP = IP + 1
                    EXTMIN(IZ) = min(EXTINCTP_(IP), EXTMIN(IZ))
                    EXTMAX(IZ) = max(EXTINCTP_(IP), EXTMAX(IZ))
                end do
            end do
        end do
        JZ = 0
        do IZ = 1, NPZ_
            if (EXTMAX(IZ) - EXTMIN(IZ) .gt. EXTBEAMCUT) JZ = IZ
        end do
        JZ = min(NPZ_, JZ + 1)
        UNIFORMZLEV = ZLEVELS_(JZ)
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
    X = XI - XSTART_
    Y = YI - YSTART_

    ! Find the grid location
    IL = 0
    IU = NPZ_
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZLEVELS_(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    K = max(IL, 1)

    I = int(X/DELXD) + 1
    if (I .gt. NPX_ .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) I = NPX_
    if (I .lt. 1 .or. I .gt. NPX_) then
        write (6, *) '_computeDirectBeam: Beyond X domain', I, XI, YI, ZI
        stop
    end if

    J = int(Y/DELYD) + 1
    if (J .gt. NPY_ .and. abs(Y - YDOMAIN) .lt. 0.001*DELYD) J = NPY_
    if (J .lt. 1 .or. J .gt. NPY_) then
        write (6, *) '_computeDirectBeam: Beyond Y domain', J, XI, YI, ZI
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
    if (btest(BCFLAG, 0) .and. (abs(X) .lt. 0.01*DELXD&
        &.or. abs(X - (NPX_ - 1)*DELXD) .lt. 0.01*DELXD)) CONSTX = .true.
    if (btest(BCFLAG, 1) .and. (abs(Y) .lt. 0.01*DELYD&
        &.or. abs(Y - (NPY_ - 1)*DELYD) .lt. 0.01*DELYD)) CONSTY = .true.

    ! If have multiple subdomains (processors) and ray is going outwards
    !   from a boundary then set the HITBOUNDARY flag
    HITBOUNDARY = .false.
    if (btest(BCFLAG, 2)) then
        if (CX .gt. 0.0 .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) then
            SIDE = 2
            HITBOUNDARY = .true.
        else if (CX .lt. 0.0 .and. abs(X) .lt. 0.001*DELXD) then
            SIDE = 1
            HITBOUNDARY = .true.
        end if
    end if
    if (btest(BCFLAG, 3)) then
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
    do while (.not. HITBOUNDARY .and. abs(ZE - ZLEVELS_(NPZ_)) .gt. EPSZ)
        IP = I + 1
        if (I .eq. NPX_) then
            if (btest(BCFLAG, 0) .or. btest(BCFLAG, 2)) then
                IP = NPX_
            else
                IP = 1
            end if
        end if
        JP = J + 1
        if (J .eq. NPY_) then
            if (btest(BCFLAG, 1) .or. btest(BCFLAG, 3)) then
                JP = NPY_
            else
                JP = 1
            end if
        end if
        X0 = DELXD*(I - 1)
        X1 = X0 + DELXD
        Y0 = DELYD*(J - 1)
        Y1 = Y0 + DELYD
        Z0 = ZLEVELS_(K)
        Z1 = ZLEVELS_(K + 1)
        if  (I .lt. 1 .or. I .gt. NPX_ .or.&
            &J .lt. 1 .or. J .gt. NPY_ .or.&
            &K .lt. 1 .or. K .ge. NPZ_) then
            write (6, '(A,3I4)') '_computeDirectBeam: beyond grid!', I, J, K
            write (6, '(1(2X,3F9.5))') X, Y, Z
            stop
        end if
        ! Get the eight corner extinction values
        I1 = K + NPZ_*(J - 1) + NPZ_*NPY_*(I - 1)
        I2 = K + NPZ_*(J - 1) + NPZ_*NPY_*(IP - 1)
        I3 = K + NPZ_*(JP - 1) + NPZ_*NPY_*(I - 1)
        I4 = K + NPZ_*(JP - 1) + NPZ_*NPY_*(IP - 1)
        E1 = EXTDIRP_(I1)
        E2 = EXTDIRP_(I2)
        E3 = EXTDIRP_(I3)
        E4 = EXTDIRP_(I4)
        E5 = EXTDIRP_(I1 + 1)
        E6 = EXTDIRP_(I2 + 1)
        E7 = EXTDIRP_(I3 + 1)
        E8 = EXTDIRP_(I4 + 1)

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
                if (btest(BCFLAG, 0)) then
                    I = 1
                    CONSTX = .true.
                else if (btest(BCFLAG, 2)) then
                    HITBOUNDARY = .true.
                    SIDE = 1
                else
                    I = NPX_
                    XOFFS = XDOMAIN
                end if
            else if (I .ge. NPX_ .and. btest(BCFLAG, 2)) then
                HITBOUNDARY = .true.
                SIDE = 2
            else if (I .eq. NPX_ + 1) then
                if (btest(BCFLAG, 0)) then
                    I = NPX_
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
                if (btest(BCFLAG, 1)) then
                    J = 1
                    CONSTY = .true.
                else if (btest(BCFLAG, 3)) then
                    HITBOUNDARY = .true.
                    SIDE = 3
                else
                    J = NPY_
                    YOFFS = YDOMAIN
                end if
            else if (J .ge. NPY_ .and. btest(BCFLAG, 3)) then
                HITBOUNDARY = .true.
                SIDE = 4
            else if (J .eq. NPY_ + 1) then
                if (btest(BCFLAG, 1)) then
                    J = NPY_
                    CONSTY = .true.
                else
                    J = 1
                    YOFFS = -YDOMAIN
                end if
            end if
        end if
        if (SO .lt. -EPSS) then
            write (6, *) '_computeDirectBeam: SO<0', X, Y, Z,&
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
        PATH = PATH + SO*(A + 0.5d0*B + 0.3333333333333333d0*C + 0.25d0*D)

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP
    end do
    DIRFLUX = SOLARFLUX*exp(-PATH)
    DIRPATH = PATH
    XO = XE + XSTART_; YO = YE + YSTART_; ZO = ZE
    VALIDBEAM = abs(ZE - ZLEVELS_(NPZ_)) .lt. EPSZ
    if (VALIDBEAM) SIDE = 6

    return

end subroutine shd_dir__computeDirect

end module shd_dir
