module shd_comp

    use globals
    use shd_pha
    use shd_src
    use shd_init
    use shd_grid

    implicit none
    private

    public :: shd_comp__computeRad
    public :: shd_comp__computeNetFluxDiv
    public :: shd_comp__visualizeRad

    public :: shd_comp__sumOutput
    public :: shd_comp__interpOutput

    ! the following two subroutines are placed outside
    ! module for general access
    !╭────────────────────────────────────────────────╮!
    !public :: shd_comp__makeDirect
    !public :: shd_comp__computeDirect
    !╰────────────────────────────────────────────────╯!

contains

subroutine shd_comp__computeRad(&
    &NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
    &SHPTR, SOURCE, SOURCE1, GRIDRAD,&
    &OUTPARMS, IRAD, RADOUT)

    ! Computes the radiances for the specified locations and direction

    integer NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer NMU, NPHI0MAX, NPHI0(NMU), IRAD
    integer MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real DIRFLUX(NPTS), FLUXES(2, NPTS)
    real SOURCE1(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    real SOURCE(NSTOKES, *)
    real OUTPARMS(*), RADOUT(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    character(len=2) :: SFCTYPE
    integer I, IBC, IANG, JX, JY, K, SIDE
    integer NANGOUT, NXOUT, NYOUT
    logical LAMBERTIAN, VALIDRAD
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
    ! Compute the upwelling bottom radiances using the downwelling flu
    if (SFCTYPE .eq. 'FL') then

        call shd_init__makeSfcRadFL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))

    else if (SFCTYPE .eq. 'VL') then

        call shd_init__makeSfcRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
            &NSTOKES, BCRAD(1, 1 + NTOPPTS))

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

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'

    ! Loop over the radiance directions
    do IANG = 1, NANGOUT
        MUOUT = OUTPARMS(2*IANG + 5)
        PHID = OUTPARMS(2*IANG + 6)
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

        if (MUOUT .eq. 0.0 .or. abs(MUOUT) .gt. 1.0) then
            write (6, *) 'Warning [_computeRad]: Bad mu for radiance', MUOUT
        else

            ! Compute the source function throughout grid for this angle
            call shd_src__computeSrcOneDir(NSTOKES, ML, MM, NLM, NSTLEG, NLEG,&
                &NUMPHASE, NPTS, DELTAM, MUOUT, PHIOUT, SRCTYPE,&
                &SOLARMU, SOLARAZ, ALBEDO, LEGEN,&
                &IPHASE, DIRFLUX, SHPTR, SOURCE, SOURCE1)

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
                if (.not. LAMBERTIAN) then
                    call shd_init__makeSfcRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                        &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT,&
                        &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                        &SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES,&
                        &BCRAD(1, 1 + NTOPPTS))
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
        end if
    end do

    return

end subroutine shd_comp__computeRad

subroutine shd_comp__computeNetFluxDiv(NSTOKES, NPTS, RSHPTR, SRCTYPE,&
    &SOLARMU, EXTINCT, ALBEDO, PLANCK, DIRFLUX,&
    &RADIANCE, NETFLUXDIV)

    ! Computes the net flux divergence at every grid point.

    integer NSTOKES, NPTS, RSHPTR(NPTS + 1)
    real PLANCK(NPTS), DIRFLUX(NPTS), SOLARMU
    real EXTINCT(NPTS), ALBEDO(NPTS)
    real RADIANCE(NSTOKES, *), NETFLUXDIV(NPTS)
    character(len=1) :: SRCTYPE
    integer I, IR
    real PI, C, SECMU, HR

    PI = acos(-1.0)
    C = sqrt(4.0*PI)
    SECMU = 1.0/abs(SOLARMU)
    do I = 1, NPTS
        IR = RSHPTR(I)
        HR = C*EXTINCT(I)*(1 - ALBEDO(I))*RADIANCE(1, IR + 1)
        if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
            HR = HR + EXTINCT(I)*(1.0 - ALBEDO(I))*DIRFLUX(I)*SECMU
        end if
        if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
            HR = HR - 4*PI*EXTINCT(I)*PLANCK(I)
        end if
        NETFLUXDIV(I) = -HR
    end do

    return

end subroutine shd_comp__computeNetFluxDiv

subroutine shd_comp__visualizeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
    &SHPTR, SOURCE, OUTPARMS, IVIS, VISOUT)

    !   Computes Stokes radiances (output in VISOUT) for the visualization
    !  modes: 1) camera mode, and 2) cross track scanning.

    integer NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer NMU, NPHI0MAX, NPHI0(NMU), IVIS
    integer MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(NSTOKES, *)
    real OUTPARMS(*), VISOUT(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    character(len=2) :: SFCTYPE
    integer NSCATANGLE, NSTPHASE, I, J, L, N, SIDE
    integer NL, NS, LINE, SAMP
    logical CAMERA_MODE, VALIDRAD
    real MURAY, PHIRAY
    real(RD_) :: XE, YE, ZE, TRANSMIT, VISRAD(NSTOKES)
    real(RD_) :: X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
    real(RD_) :: THETA0, THETA1, PHIR, PHI0
    real(RD_) :: COSTH, SINTH, SINTH0, SINTH1, MU2, PHI2
    real(RD_) :: U, V, UP, VP, COSDELPHI, PI, ROTANG, DEGRAD, R
    real(RD_) :: DIST, D, RX, RY, RZ, SCANANG
    integer MAXSCATANG
    parameter(MAXSCATANG=721)
    real, allocatable :: YLMSUN(:, :), PHASETAB(:, :, :)
    real              :: MEAN, STD1, STD2
    real, allocatable :: AOLP(:)

    allocate (YLMSUN(NSTLEG, NLM))

    if (SRCTYPE .ne. 'T') then
        call shd_pha__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG, YLMSUN)
        if (DELTAM .and. NUMPHASE .gt. 0) then
            NSCATANGLE = max(36, min(MAXSCATANG, 2*NLEG))
            NSTPHASE = min(NSTLEG, 2)
            allocate (PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE))
            call shd_pha__computePhase(NSCATANGLE, NUMPHASE, NSTPHASE,&
                &NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
        end if
    end if

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces
    !   Compute the upwelling bottom radiances using the downwelling
    if (SFCTYPE .eq. 'FL') then
        call shd_init__makeSfcRadFL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
    else if (SFCTYPE .eq. 'VL') then
        call shd_init__makeSfcRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
            &NSTOKES, BCRAD(1, 1 + NTOPPTS))
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
    PI = acos(-1.0d0)
    DEGRAD = PI/180.

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
                UP = (SAMP - NS/2 - 1)*OUTPARMS(13)*DEGRAD
                VP = (LINE - NL/2 - 1)*OUTPARMS(12)*DEGRAD
                ROTANG = OUTPARMS(9)*DEGRAD
                if (ROTANG .eq. 0.0) then
                    U = UP
                    V = VP
                else
                    U = cos(ROTANG)*UP - sin(ROTANG)*VP
                    V = sin(ROTANG)*UP + cos(ROTANG)*VP
                end if
                THETA0 = DEGRAD*OUTPARMS(7)
                PHI0 = DEGRAD*OUTPARMS(8)
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
                if (abs(SINTH) .lt. 1.0e-6) then
                    PHIR = 0.0
                else
                    if (abs(SINTH0) .lt. 1e-6 .and. abs(SINTH1) .le. 1e-6) then
                        COSDELPHI = 0.0d0
                    else if (abs(SINTH1) .gt. 1.0e-6) then
                        COSDELPHI = (cos(U) - COSTH*cos(THETA1))/(SINTH1*SINTH)
                    else if (abs(SINTH0) .gt. 1.0e-6) then
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
                SCANANG = DEGRAD*(OUTPARMS(11) + (SAMP - 1)*OUTPARMS(13))
                COSTH = cos(PI - abs(SCANANG))
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
                write (6, *) 'Error [VISUALIZE_RADIANCE]: Level below domain.'
                stop
            end if

            ! MURAY,PHRAY is camera pixel viewing direction;
            ! MU2,PHI2 is radiance travel direction (the opposite)
            MURAY = SNGL(COSTH)
            PHIRAY = SNGL(PHIR)
            MU2 = -COSTH
            PHI2 = PHIR + PI

            ! Integrate the extinction and source function along this ray
            !   to calculate the Stokes radiance vector for this pixel
            TRANSMIT = 1.0d0; VISRAD(:) = 0.0d0
            call shd_src__integrateSrcOneRay(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
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
                &XE, YE, ZE, SIDE, TRANSMIT, VISRAD, VALIDRAD)
900         continue
            ! write (6,'(1X,2F8.4,1X,2F11.7,4(1X,F11.6))')&
            !     &X0,Y0,MURAY,PHIRAY,VISRAD(:)
            VISOUT(1, IVIS) = VISRAD(1)
            if (VISRAD(1) .gt. 0.0) then
                if (NSTOKES .gt. 1) then
                    ! Output degree (0 to 1) and angle (-180 to 180) of linear pol
                    VISOUT(2, IVIS) = sqrt(VISRAD(2)**2 + VISRAD(3)**2)/VISRAD(1)
                    VISOUT(3, IVIS) = (180/PI)*0.5*atan2(VISRAD(3), VISRAD(2))
                end if
                if (NSTOKES .eq. 4) then
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
    if (NSTOKES .gt. 1) then
        ! Choose the best range for the angle of linear polarization (-90
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
end subroutine shd_comp__visualizeRad

subroutine shd_comp__interpOutput(OLDNPTS, NPTS,&
    &FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
    &FLUXES, DIRFLUX, FLUXDIV, SHTERMS,&
    &NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)

    !   Interpolates the output data to the new grid points (those >OLDNPTS)
    ! from the old grid points.  The output data is in FLUXES, DIRFLUX,
    ! FLUXDIV, and SHTERMS.  OLDNPTS is set to zero to do no interpolation.
    ! Each new point is located in a cell, and the two parent grid points
    ! are found to interpolate the data from.

    integer OLDNPTS, NPTS, NBCELLS, NCELLS
    integer GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    logical FLUX_OUT, FLUXDIV_OUT, SH_OUT
    real FLUXES(2, NPTS), DIRFLUX(NPTS), FLUXDIV(NPTS)
    real SHTERMS(5, NPTS)
    real GRIDPOS(3, *)
    integer IP, IC, ICELL, IPARENT, IDIR, I, J, DIR1, DIR2, IP1, IP2
    integer GRIDCORNER(2, 4, 3)
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

                stop 'Error [_interpOutput]: Point not on an edge.'

120             continue

                ! Interpolate the data from the two parent points at the
                !   ends of the edge; assume the cell was split in half.
                if (FLUX_OUT) then
                    FLUXES(1, IP) = 0.5*(FLUXES(1, IP1) + FLUXES(1, IP2))
                    FLUXES(2, IP) = 0.5*(FLUXES(2, IP1) + FLUXES(2, IP2))
                    DIRFLUX(IP) = 0.5*(DIRFLUX(IP1) + DIRFLUX(IP2))
                end if
                if (FLUXDIV_OUT) then
                    FLUXDIV(IP) = 0.5*(FLUXDIV(IP1) + FLUXDIV(IP2))
                end if
                if (SH_OUT) then
                    do J = 1, 5
                        SHTERMS(J, IP) = 0.5*(SHTERMS(J, IP1) + SHTERMS(J, IP2))
                    end do
                end if
            end if
        end do
    end do

    return
end subroutine shd_comp__interpOutput

subroutine shd_comp__sumOutput(IFLAG, WT, N, data, SUMDATA)
    !     Sums an output array over the k-distribution.
    ! If IFLAG=0 then the output array is set to the input times the wei

    integer IFLAG, N, I
    real WT, data(N), SUMDATA(N)
    if (IFLAG .eq. 0) then
        do I = 1, N
            SUMDATA(I) = WT*data(I)
        end do
    else
        do I = 1, N
            SUMDATA(I) = SUMDATA(I) + WT*data(I)
        end do
    end if
    return

end subroutine shd_comp__sumOutput

end module shd_comp

subroutine shd_comp__makeDirect(&
    &NPTS, BCFLAG, IPFLAG, DELTAM,&
    &ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
    &GRIDPOS, DIRFLUX)

    ! Makes the direct beam solar flux for the internal base grid.
    !   DIRFLUX is set to F*exp(-tau_sun).
    !   Actually calls shd_comp__computeDirect to do all the hard work.

    implicit none

    integer NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), DIRFLUX(NPTS)
    integer SIDE, IP
    logical VALIDBEAM
    real UNIFZLEV, XO, YO, ZO, DIR, DIRPATH

    call shd_comp__computeDirect(&
        &1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
        &DELTAM, ML, NSTLEG, NLEG,&
        &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)

    do IP = 1, NPTS
        DIRPATH = 0.0
        call shd_comp__computeDirect&
            &(0, GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &BCFLAG, IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
            &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(IP),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end do

    return
end subroutine shd_comp__makeDirect

subroutine shd_comp__computeDirect(&
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

    use globals
    use globals_prp

    implicit none

    integer INIT, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, SIDE
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

end subroutine shd_comp__computeDirect
