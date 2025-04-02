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
    implicit none
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
    double precision X0, Y0, Z0, XE, YE, ZE, TRANSMIT, RADIANCE(4)
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
            call shd_comp__computeOneSrc(NSTOKES, ML, MM, NLM, NSTLEG, NLEG,&
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
                    call shd_comp__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
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
    implicit none
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

    implicit none
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
    double precision XE, YE, ZE, TRANSMIT, VISRAD(NSTOKES)
    double precision X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2
    double precision THETA0, THETA1, PHIR, PHI0
    double precision COSTH, SINTH, SINTH0, SINTH1, MU2, PHI2
    double precision U, V, UP, VP, COSDELPHI, PI, ROTANG, DEGRAD, R
    double precision DIST, D, RX, RY, RZ, SCANANG
    integer MAXSCATANG
    parameter(MAXSCATANG=721)
    real, allocatable :: YLMSUN(:, :), PHASETAB(:, :, :)
    real              :: MEAN, STD1, STD2
    real, allocatable :: AOLP(:)

    allocate (YLMSUN(NSTLEG, NLM))

    if (SRCTYPE .ne. 'T') then
        call shd_comp__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG, YLMSUN)
        if (DELTAM .and. NUMPHASE .gt. 0) then
            NSCATANGLE = max(36, min(MAXSCATANG, 2*NLEG))
            NSTPHASE = min(NSTLEG, 2)
            allocate (PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE))
            call shd_comp__computePhase(NSCATANGLE, NUMPHASE, NSTPHASE,&
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
            call shd_comp__integrateOneRay(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
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

subroutine shd_comp__integrateOneRay(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
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

    use shd_grid

    implicit none

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
    double precision MU2, PHI2, X0, Y0, Z0, XE, YE, ZE
    double precision TRANSMIT, RADIANCE(NSTOKES)
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
    double precision PI, CX, CY, CZ, CXINV, CYINV, CZINV
    double precision XN, YN, ZN, XI, YI, ZI
    double precision SO, SOX, SOY, SOZ, EPS
    double precision TAUTOL, TAUGRID, S, DELS, TRANSCUT
    double precision EXT, TAU, TRANSCELL, ABSCELL
    double precision SRC(NSTOKES)
    double precision U, V, W, DELX, DELY, DELZ, INVDELX, INVDELY, INVDELZ
    double precision COSSCAT
    real, allocatable :: YLMDIR(:, :), SINGSCAT(:, :)
    double precision, allocatable :: SUNDIRLEG(:)

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
    PI = acos(-1.0d0)

    ! Calculate the generalized spherical harmonics for this direction
    allocate (YLMDIR(NSTLEG, NLM))
    call shd_comp__computeYLM(.false., SNGL(MU2), SNGL(PHI2), ML, MM, NSTLEG, YLMDIR)

    if (SRCTYPE .ne. 'T' .and. DELTAM) then
        ! Get the solar single scattering Stokes vector for the outgoin
        !   direction my interpolating in scattering angle in the PHASETA
        !   table and then rotating the Q/U polarization to the outgoing
        COSSCAT = SOLARMU*MU2 + sqrt((1.0 - SOLARMU**2)*(1.0 - MU2**2))&
            &*cos(SOLARAZ - PHI2)
        if (NUMPHASE .gt. 0) then
            allocate (SINGSCAT(NSTOKES, NUMPHASE))
            F = (NSCATANGLE - 1)*(acos(COSSCAT)/PI) + 1
            J = min(NSCATANGLE - 1, int(F))
            F = F - J
            do I = 1, NUMPHASE
                SINGSCAT(1:NSTPHASE, I)&
                    &= (1 - F)*PHASETAB(:, I, J) + F*PHASETAB(:, I, J + 1)
                if (NSTOKES .gt. 1) then
                    call shd_comp__rotatePolarPlane(NSTOKES, COSSCAT, SOLARMU,&
                        &SNGL(MU2), SOLARAZ - SNGL(PHI2), SINGSCAT(:, I))
                end if
            end do
        else
            allocate (SUNDIRLEG(0:NLEG))
            call shd_comp__computeLegendre(COSSCAT, NLEG, SUNDIRLEG)
            do L = 0, NLEG
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)/(4*PI)
            end do
        end if
    end if

    ! Make the ray direction (opposite to the outgoing direction)
    CX = sqrt(1.0d0 - MU2**2)*cos(PHI2 - PI)
    CY = sqrt(1.0d0 - MU2**2)*sin(PHI2 - PI)
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
    XM = 0.5*(XGRID(1) + XGRID(NX))
    YM = 0.5*(YGRID(1) + YGRID(NY))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    ! ZE = MAX(MIN(Z0,DBLE(ZGRID(NZ))),DBLE(ZGRID(1)))
    call shd_grid__findGridCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
        &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
        &BCFLAG, IPFLAG, XE, YE, ZE, ICELL)
    ! print '(A,5(1X,F8.5),1X,I6)',&
    !     &'_integrateRay1:',X0,Y0,Z0,MU2,PHI2,ICELL
    IFACE = 0
    NGRID = 0

    ! Loop until reach a Z boundary or transmission is very small
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
        ! Make sure current cell is valid
        if (ICELL .le. 0) then
            write (6, *) 'Error [_integrateRay1]: ICELL=', ICELL,&
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
            call shd_comp__computeSrcOneCell_(ICELL, GRIDPTR,&
                &ML, MM, NLM, NLEG, NUMPHASE,&
                &NPTS, DELTAM, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,&
                &SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
                &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
                &EXTINCT8, SRCEXT8)
        else
            call shd_comp__computeSrcOneCell(ICELL, GRIDPTR,&
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
            write (6, *) 'Error [_integrateRay1]: SO<0 ',&
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
                ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - 0.33333333333*TAU))
                TRANSCELL = 1.0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:))&
                    &+ 0.08333333333*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*DELS&
                    &*(1.0 - 0.05*(EXT1 - EXT0)*DELS))/EXT
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
            call shd_grid__getNextCell(XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,&
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
            call shd_comp__findBoundaryRad(NSTOKES, XN, YN,&
                &SNGL(MU2), SNGL(PHI2),&
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
end subroutine shd_comp__integrateOneRay

subroutine shd_comp__findBoundaryRad(NSTOKES, XB, YB, MU2, PHI2,&
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
    implicit none
    integer NSTOKES, ICELL, KFACE, MAXNBC, NTOPPTS, NBOTPTS
    integer GRIDPTR(8, *), BCPTR(MAXNBC, 2)
    integer NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR
    real MU2, PHI2, RADBND(NSTOKES)
    double precision XB, YB
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
end subroutine shd_comp__findBoundaryRad

subroutine shd_comp__computeSrcOneCell(ICELL, GRIDPTR,&
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
    implicit none
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
    double precision SUNDIRLEG(0:NLEG)
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
end subroutine shd_comp__computeSrcOneCell

subroutine shd_comp__computeSrcOneCell_(ICELL, GRIDPTR,&
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

    implicit none
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
    double precision SUNDIRLEG(0:NLEG)
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
end subroutine shd_comp__computeSrcOneCell_

subroutine shd_comp__rotatePolarPlane(NSTOKES, COSSCAT,&
    &SOLARMU, MU, DELPHI, SCATVECT)

    ! Rotates the plane of polarization from the scattering plane
    ! to the outgoing meridional plane for the solar single
    ! scattered Stokes vector in SINGSCAT.

    implicit none
    integer NSTOKES
    double precision COSSCAT
    real SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    real B1
    double precision SIN_SCAT, SIN_THETA1, SIN_THETA2
    double precision SINPHI, COSPHI, SIN2, COS2, SIN22, COS22

    if (NSTOKES .gt. 1) then
        B1 = SCATVECT(2)
        SIN_SCAT = DSQRT(max(0.0d0, 1.d0 - COSSCAT**2))
        SIN_THETA1 = DSQRT(1.d0 - SOLARMU**2)
        SIN_THETA2 = DSQRT(1.d0 - MU**2)
        SINPHI = DSIN(dble(DELPHI))
        COSPHI = DCOS(dble(DELPHI))
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

end subroutine shd_comp__rotatePolarPlane

subroutine shd_comp__interpOutput(OLDNPTS, NPTS,&
    &FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
    &FLUXES, DIRFLUX, FLUXDIV, SHTERMS,&
    &NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)

    !   Interpolates the output data to the new grid points (those >OLDNPTS)
    ! from the old grid points.  The output data is in FLUXES, DIRFLUX,
    ! FLUXDIV, and SHTERMS.  OLDNPTS is set to zero to do no interpolation.
    ! Each new point is located in a cell, and the two parent grid points
    ! are found to interpolate the data from.

    implicit none
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
    implicit none
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
