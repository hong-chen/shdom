module shd_intg

    use globals
    use shd_sfc
    use shd_init
    use shd_grid
    use shd_util
    use shd_shdo
    use shd_mpi

    implicit none

    public :: shd_intg__integratePath

    private :: shd_intg__backIntegrate3D
    private :: shd_intg__backIntegrate3D_
    private :: shd_intg__backIntegrate2D
    private :: shd_intg__backIntegrate1D

contains

subroutine shd_intg__integratePath(&
    &NX, NY, NZ, NSTOKES, NSTLEG,&
    &ML, MM, NLM, NMU, NPHI0MAX, NPHI0, NANG,&
    &BCFLAG, IPFLAG,&
    &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN,&
    &UNITS, SFCTYPE, GNDTEMP, GNDALBEDO,&
    &NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
    &SFCGRIDPARMS, UNIFORM_SFC_BRDF, SFC_BRDF_DO,&
    &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &MU, PHI, WTDO,&
    &CMU1, CPHI1, CMU2, CPHI2, WSAVE, FFTFLAG,&
    &DIRFLUX, FLUXES, EXTINCT,&
    &SHPTR, SOURCE, RSHPTR, RADIANCE,&
    &WORK, SWEEPORD, GRIDRAD)

    !   Performs the path integrations through the medium specified by
    ! the extinction (EXTINCT) and source function (SOURCE) in
    ! spherical harmonic space.  The source function is transformed to
    ! discrete ordinates one zenith angle at a time.  The integrations
    ! are done from the top boundary (with no reflection) down to the
    ! bottom boundary for the downwelling angles, and then up from the
    ! bottom to the top after computing the surface reflection/emission.
    ! After the path integrations produce the discrete angle radiance
    ! it is transformed and added to output radiance spherical harmonic
    ! expansion.  The hemispheric flux FLUXES is also computed.
    ! The discrete ordinates in MU, PHI must be a complete set with the
    ! downwelling (mu<0) angles first.
    !   The surface reflection is handled differently for Lambertian
    ! surfaces and more general surfaces specified by bidirectional
    ! reflection distribution functions (BRDF).  For the general BRDF,
    ! the bottom boundary radiance must be computed for each upwelling
    ! ordinate, while for the Lambertian surface the boundary radiances
    ! can be computed just once, since they are isotropic.  For Lambertian
    ! cases the BCRAD array holds only the isotropic radiance for the
    ! boundaries at the top (downwelling) and bottom (upwelling) boundaries.
    ! For the general BRDF, the BCRAD array in addition holds the
    ! downwelling radiances at all discrete ordinates for the bottom
    ! boundary points.

    integer NX, NY, NZ, NSTOKES, NSTLEG, ML, MM, NLM
    integer NMU, NPHI0MAX, NPHI0(NMU), NANG
    integer BCFLAG, IPFLAG
    integer NPTS, NCELLS, MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS
    integer NXSFC, NYSFC, NSFCPAR
    integer SHPTR(NPTS + 1), RSHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    integer SWEEPORD(NPTS, *)
    logical FFTFLAG(NMU), UNIFORM_SFC_BRDF
    real SOLARMU, SOLARAZ, SKYRAD
    real GNDTEMP, GNDALBEDO
    real DELXSFC, DELYSFC
    real SFCPARMS(NSFCPAR, *), SFCGRIDPARMS(NSFCPAR, MAXNBC)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    real WAVENO(2), WAVELEN
    real XGRID(NX), YGRID(NY), ZGRID(NZ), GRIDPOS(3, NPTS)
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real CMU1(*), CPHI1(*), CMU2(*), CPHI2(*), WSAVE(*)
    real BCRAD(NSTOKES, MAXBCRAD)
    real DIRFLUX(NPTS), FLUXES(2, NPTS), EXTINCT(NPTS)
    real GRIDRAD(NSTOKES, NPTS), WORK(NSTOKES, NPHI0MAX, NPTS)
    real SOURCE(NSTOKES, *), RADIANCE(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    character(len=2) :: SFCTYPE

    integer I, I1, K, NR, IPHI, IMU, IBC, IANG, IUPDOWN, OLDNPTS
    integer joct, ipt, ip, iorder, ipcell, icorner
    logical :: LAMBERTIAN, theflag
    real A, TRANSMIN
    save OLDNPTS

    integer :: count_start, count_rate, count_max
    integer :: count_end
    real :: elapsed_time


    ! Set the minimum transmission for the cell tracing (1 for single cell)
    TRANSMIN = 1.00

    ! Make the new grid cell/point sweeping order if need to
    if (NPTS .ne. OLDNPTS) then
        call shd_grid__makeSweepOrder(NX, NY, NZ, NPTS, NCELLS,&
            &GRIDPTR, TREEPTR, CELLFLAGS,&
            &BCFLAG, IPFLAG, SWEEPORD, GRIDRAD)
    end if

    call shd_mpi__findBoundaryPnt(BCFLAG, IPFLAG, NPTS, SWEEPORD, GRIDPTR,&
        &GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'
    if (NPTS .ne. OLDNPTS) then
        ! Make the pointers to the grid points on the top and bottom boundaries
        call shd_grid__findBoundaryPnt(NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
            &ZGRID(1), ZGRID(NZ), GRIDPOS,&
            &NTOPPTS, NBOTPTS, BCPTR)
        if (SFCTYPE(1:1) .eq. 'V') then
            ! If this is a variable surface then bilinearly interpolate
            !   the surface temperature and reflection parameters to the
            !   bottom grid points.
            call shd_sfc__interpSfc(NBOTPTS, BCPTR(1, 2), GRIDPOS,&
                &SRCTYPE, WAVENO, WAVELEN, UNITS,&
                &NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
                &SFCGRIDPARMS)
        end if
    end if
    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))

    ! Zero the flux arrays and set temporary radiance to -1 for not valid
    do I = 1, NPTS
        FLUXES(1, I) = 0.0
        FLUXES(2, I) = 0.0
        GRIDRAD(1, I) = -1.0
    end do

    ! Zero the output spherical harmonic radiance
    NR = RSHPTR(NPTS + 1)
    do I = 1, NR
        RADIANCE(:, I) = 0.0
    end do

    ! print *, 'Path integration:'
    ! Loop over the zenith angles (downwelling angles must be first)
    IANG = 1

! -- loop mu
    do IMU = 1, NMU

        if (IMU .eq. NMU/2 + 1) then
            ! For the first upwelling angle, make the bottom boundary radiances
            !   for the Lambertian surfaces.  Compute the upwelling bottom
            !   radiances using the downwelling fluxes.
            if (SFCTYPE .eq. 'FL') then
                call shd_init__makeSfcRadFL(NBOTPTS, BCPTR(1, 2),&
                    &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
                    &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
            else if (SFCTYPE .eq. 'VL') then
                call shd_init__makeSfcRadVL(NBOTPTS, BCPTR(1, 2),&
                    &DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
                    &NSTOKES, BCRAD(1, 1 + NTOPPTS))
            end if
        end if

        if (NSTOKES .eq. 1) then
            call shd_shdo__transSHtoDO_(NPTS, ML, MM, NLM,&
                &NMU, NPHI0MAX, NPHI0(IMU),&
                &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
                &SHPTR, SOURCE, WORK)
        else
            call shd_shdo__transSHtoDO(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
                &NMU, NPHI0MAX, NPHI0(IMU),&
                &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
                &SHPTR, SOURCE, WORK)
        end if

        ! Enforce the discrete ordinate I source function to be non-negative
        do I = 1, NPTS
            do IPHI = 1, NPHI0(IMU)
                WORK(1, IPHI, I) = max(0.0, WORK(1, IPHI, I))
            end do
        end do

! -- loop phi
        ! Loop over the azimuthal angles
        do IPHI = 1, NPHI0(IMU)

            if (MU(IMU) .lt. 0.0) then

                ! For downward ordinates, initialize the top boundary radiances
                IUPDOWN = 1
                do IBC = 1, NTOPPTS
                    I = BCPTR(IBC, 1)
                    GRIDRAD(:, I) = BCRAD(:, IBC)
                end do

            else

                ! Initialize the bottom boundary radiances
                IUPDOWN = 2
                if (.not. LAMBERTIAN) then

                    ! If a Lambertian surface, use the radiances computed above
                    !   otherwise, compute the radiance for this angle
                    !   by integrating over the stored downwelling radiances.
                    ! If the surface reflection is uniform use the precomputed
                    !   BRDF matrix to do the integration.
                    if (UNIFORM_SFC_BRDF) then
                        call shd_init__makeSfcRadU(NBOTPTS, BCPTR(1, 2),&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, IMU, IPHI,&
                            &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS, SFC_BRDF_DO,&
                            &NSTOKES, BCRAD(1, 1 + NTOPPTS))
                    else
                        call shd_init__makeSfcRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &MU(IMU), PHI(IMU, IPHI),&
                            &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            &NSTOKES, BCRAD(1, 1 + NTOPPTS))
                    end if
                end if

                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    GRIDRAD(:, I) = BCRAD(:, IBC + NTOPPTS)
                end do

            end if

            ! Do backward path integration from each grid point not done yet
            if (btest(IPFLAG, 1) .and. btest(IPFLAG, 0)) then

                call shd_intg__backIntegrate1D(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &GRIDPOS, SWEEPORD, MU(IMU), PHI(IMU, IPHI), TRANSMIN,&
                    &EXTINCT, NSTOKES, WORK, IPHI, GRIDRAD)

            else if (btest(IPFLAG, 1) .and.&
                &.not. (btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) then

                call shd_intg__backIntegrate2D(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &GRIDPOS, SWEEPORD, MU(IMU), PHI(IMU, IPHI), TRANSMIN,&
                    &EXTINCT, NSTOKES, WORK, IPHI, GRIDRAD)

            else
                if (NSTOKES .eq. 1) then

                    call shd_intg__backIntegrate3D_(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                        &GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
                        &MU(IMU), PHI(IMU, IPHI), TRANSMIN, BCFLAG, IPFLAG,&
                        &EXTINCT, WORK, IPHI, GRIDRAD)

                else

                    call shd_intg__backIntegrate3D(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                        &GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
                        &MU(IMU), PHI(IMU, IPHI), TRANSMIN, BCFLAG, IPFLAG,&
                        &EXTINCT, NSTOKES, WORK, IPHI, GRIDRAD)

                end if
            end if

            ! Transfer the radiance for the (IMU,IPHI) discrete ordinate from
            !   GRIDRAD to the work array.
            do I = 1, NPTS
                if (GRIDRAD(1, I) .lt. -0.0001) then
                    write (6, *) 'Warning [_integratePath]: Grid point has no value',&
                        &I, MU(IMU), PHI(IMU, IPHI), GRIDRAD(1, I), XGRID(1), YGRID(1),&
                        &GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I)
                    ! stop
                    GRIDRAD(1, I) = 0.0
                else
                    WORK(:, IPHI, I) = GRIDRAD(:, I)
                    GRIDRAD(1, I) = -1.0
                end if
            end do

            ! Do the sum for the hemispheric fluxes for this angle
            A = abs(MU(IMU))*WTDO(IMU, IPHI)
            do I = 1, NPTS
                FLUXES(IUPDOWN, I) = FLUXES(IUPDOWN, I) + A*WORK(1, IPHI, I)
                !if (IUPDOWN .eq. 1) then
                !    print *, IMU, IPHI, I, FLUXES(IUPDOWN, I)
                !end if
            end do

            ! For a general BRDF surface save the downwelling radiance
            !   at the bottom boundary points.
            if (.not. LAMBERTIAN .and. IMU .le. NMU/2) then
                I1 = NTOPPTS + NBOTPTS*IANG
                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    BCRAD(:, IBC + I1) = WORK(:, IPHI, I)
                end do
            end if

            IANG = IANG + 1
        end do

        ! Transform the radiance to spherical harmonics for this zenith angle
        if (NSTOKES .eq. 1) then
            call shd_shdo__transDOtoSH_(NPTS, ML, MM, NLM,&
                &NMU, NPHI0MAX, NPHI0(IMU),&
                &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
                &RSHPTR, WORK, RADIANCE)
        else
            call shd_shdo__transDOtoSH(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
                &NMU, NPHI0MAX, NPHI0(IMU),&
                &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
                &RSHPTR, WORK, RADIANCE)
        end if

    end do

    OLDNPTS = NPTS

    return

end subroutine shd_intg__integratePath

subroutine shd_intg__backIntegrate3D(&
    &NX, NY, NZ, NA, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
    &MU, PHI, TRANSMIN, BCFLAG, IPFLAG,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)

    !   Sweeps through the spatial grid computing radiance vectors for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  The sweeping
    ! order through the grid points is in SWEEPORD. If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is traced
    ! back to a face containing known (valid) radiances.  As each cell is
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! at the end of the ray path is interpolated from the known grid point
    ! radiances on the face.  The resulting radiances are put in GRIDRAD.

    integer NX, NY, NZ, NA, NPTS, NCELLS, NSTOKES, KANG
    integer BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS), XGRID(*), YGRID(*), ZGRID(NZ)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITY, BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer I1, I2, I3, I4, IOPP, JFACE, KFACE, IC, IZ
    integer GRIDFACE(4, 6), OPPFACE(6)
    integer IORDER, ICORNER, JOCTORDER3(8), JOCTORDER2(8)
    logical VALIDRAD, VALIDFACE, IPINX, IPINY
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS, path
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(RD_) :: TAU, TRANSCELL, ABSCELL, TRANSMIT
    real(RD_) :: RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/
    data JOCTORDER3/1, 3, 5, 7, 2, 4, 6, 8/, JOCTORDER2/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
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

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negative.
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
        stop 'Error [_backIntegrate3D]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate goes to.
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    if (btest(IPFLAG, 1)) then
        JOCT = JOCTORDER2(IOCT)
    else
        JOCT = JOCTORDER3(IOCT)
    end if

    IZ = (NZ - 1)*(1 - BITZ) + 1
    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)

        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if ((GRIDPOS(3, IPT) .lt. ZGRID(IZ) .and. BITZ .eq. 0) .or.&
                &(GRIDPOS(3, IPT) .gt. ZGRID(IZ) .and. BITZ .eq. 1)) then
                IZ = IZ + 2*BITZ - 1
                call shd_mpi__calBoundaryRad(BCFLAG, IPFLAG, JOCT, IZ,&
                    &NX, NY, NZ, XGRID, YGRID, ZGRID,&
                    &NA, NPTS, NCELLS, GRIDPTR,&
                    &NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,&
                    &MU, PHI, EXTINCT, NSTOKES, SOURCE,&
                    &KANG, GRIDRAD)
            end if
        end if

        if (GRIDRAD(1, IPT) .lt. 0.0) then

            ! Start off at the grid point
            ICELL = IPCELL
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
            XE = GRIDPOS(1, IPT)
            YE = GRIDPOS(2, IPT)
            ZE = GRIDPOS(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)
                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate3D]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL, xe, ye, ze,&
                        &gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                    stop
                end if
                IPINX = btest(int(CELLFLAGS(ICELL)), 0)
                IPINY = btest(int(CELLFLAGS(ICELL)), 1)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = GRIDPTR(9 - IOCT, ICELL)
                ! Get the distances to the 3 planes and select the closest
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
                    write (6, *) 'Error [_backIntegrate3D]: SO<0 ',&
                        &MU, PHI, XE, YE, ZE, SOX, SOY, SOZ, SO, IPT, ICELL
                    stop
                end if
                ! Compute the coordinates of the cell exitting location
                XE = XE + SO*CX
                YE = YE + SO*CY
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                if (SOX .le. SOZ .and. SOX .le. SOY) then
                    IFACE = 2 - BITX
                    JFACE = 1
                else if (SOY .le. SOZ) then
                    IFACE = 4 - BITY
                    JFACE = 2
                else
                    IFACE = 6 - BITZ
                    JFACE = 3
                end if
                ! Get the next cell to go to
                INEXTCELL = NEIGHPTR(IFACE, ICELL)
                if (INEXTCELL .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
                end if

                ! Get the face grid pointers
                if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    KFACE = IFACE
                    IC = ICELL
                else
                    ! If going to smaller face then use next face (more accurate)
                    KFACE = OPPFACE(IFACE)
                    IC = INEXTCELL
                end if
                I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
                I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
                I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
                I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
                ! Compute the face interpolation factors
                if (JFACE .eq. 1) then
                    U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                    if (IPINY) then
                        V = 0.5
                    else
                        V = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I2) - GRIDPOS(2, I1))
                    end if
                else if (JFACE .eq. 2) then
                    U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                    end if
                else
                    if (IPINY) then
                        U = 0.5
                    else
                        U = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I3) - GRIDPOS(2, I1))
                    end if
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                    end if
                end if
                if (u .lt. -0.0001 .or. v .lt. -0.0001 .or.&
                    &u .gt. 1.0001 .or. v .gt. 1.0001) then
                    print *, '_backIntegrate3D: u,v<0 or u,v>1: ',&
                        &mu, phi, xe, ye, ze, u, v, iface, jface, kface,&
                        &gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    if (JFACE .eq. 1) then
                        XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
                    else if (JFACE .eq. 2) then
                        YE = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
                    else
                        ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                    end if
                end if

                ! Interpolate extinction and source function at face intersection
                F1 = (1 - U)*(1 - V)
                F2 = (1 - U)*V
                F3 = U*(1 - V)
                F4 = U*V
                EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)&
                    &+ F3*EXTINCT(I3) + F4*EXTINCT(I4)
                ! Correctly interpolate source using extinction*source
                SRCEXT0(:) = F1*SOURCE(:, KANG, I1)*EXTINCT(I1)&
                    &+ F2*SOURCE(:, KANG, I2)*EXTINCT(I2)&
                    &+ F3*SOURCE(:, KANG, I3)*EXTINCT(I3)&
                    &+ F4*SOURCE(:, KANG, I4)*EXTINCT(I4)
                ! Compute the cell radiance: integration of the source function
                EXT = 0.5*(EXT0 + EXT1)
                TAU = EXT*SO
                if (TAU .ge. 0.5) then
                    TRANSCELL = exp(-TAU)
                    ABSCELL = 1.0 - TRANSCELL
                else
                    ABSCELL = TAU*(1.0 - 0.5*TAU*&
                        &(1.0 - FRAC13_*TAU*(1 - 0.25*TAU)))
                    TRANSCELL = 1.0 - ABSCELL
                end if
                if (TAU .le. 2.0) then
                    if (EXT .eq. 0.0) then
                        SRC(:) = 0.0
                    else
                        ! Linear extinction, linear source*extinction, to first order
                        SRC(:) = (0.5d0*(SRCEXT0(:) + SRCEXT1(:)) + FRAC112_&
                            &*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                    end if
                else
                    ! Combined first order expansion and constant extinction formula
                    EXT0P = EXT0
                    SRCEXT0P(:) = SRCEXT0(:)
                    if (TAU .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0d0/TAU
                        if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                    end if
                    SRC(:) = 1.0d0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                        &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0d0/(EXT0P + EXT1)&
                        &*(1.0d0 - 2.0d0/TAU + 2.0d0*TRANSCELL/ABSCELL))
                end if
                SRC(1) = max(SRC(1), 0.0d0)

                ! Add in the cell radiance and update the transmission.
                ! If this is a boundary or the transmission is below the
                !   cutoff and we have a valid radiance then set the flag
                !   to stop the tracing and add in the interpolated face radiance.
                RAD(:) = RAD(:) + TRANSMIT*SRC(:)*ABSCELL
                TRANSMIT = TRANSMIT*TRANSCELL
                if (RAD(1) < -1.0e-5) then
                    print '(A,3(1X,F6.4),8(1X,E12.5))',&
                        &'_backIntegrate3D: RAD<0 ',&
                        &xe, ye, ze, ext0, ext1, so, tau, transmit, src, abscell, rad(:)
                    call shd_mpi__abort('The end')
                end if
                VALIDFACE = (GRIDRAD(1, I1) .ge. -0.1 .and. GRIDRAD(1, I2) .ge. -0.1&
                    &.and. GRIDRAD(1, I3) .ge. -0.1 .and. GRIDRAD(1, I4) .ge. -0.1)
                if (INEXTCELL .le. 0 .or.&
                    &(TRANSMIT .le. TRANSMIN .and. VALIDFACE)) then
                    if (VALIDFACE) then
                        VALIDRAD = .true.
                        RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)&
                            &+ F3*GRIDRAD(:, I3) + F4*GRIDRAD(:, I4)
                        RAD(:) = RAD(:) + TRANSMIT*RAD0(:)
                    else
                        print *, 'Error [_backIntegrate3D]: INEXTCELL=0: ', xe, ye, ze, icell
                        print *, iz, joct, mu, phi
                        print '(i7,4(1x,f7.4))', i1, GRIDRAD(1, I1),&
                            &gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                        print '(i7,4(1x,f7.4))', i2, GRIDRAD(1, I2),&
                            &gridpos(1, i2), gridpos(2, i2), gridpos(3, i2)
                        print '(i7,4(1x,f7.4))', i3, GRIDRAD(1, I3),&
                            &gridpos(1, i3), gridpos(2, i3), gridpos(3, i3)
                        print '(i7,4(1x,f7.4))', i4, GRIDRAD(1, I4),&
                            &gridpos(1, i4), gridpos(2, i4), gridpos(3, i4)
                        print '(i7,8x,3(1x,f7.4))', ipt,&
                            &gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                        stop
                    end if
                else
                    EXT1 = EXT0
                    SRCEXT1(:) = SRCEXT0(:)
                    ICELL = INEXTCELL
                end if
            end do
            GRIDRAD(:, IPT) = RAD(:)

        end if
    end do

    return
end subroutine shd_intg__backIntegrate3D

subroutine shd_intg__backIntegrate3D_(&
    &NX, NY, NZ, NA, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
    &MU, PHI, TRANSMIN, BCFLAG, IPFLAG,&
    &EXTINCT, SOURCE, KANG, GRIDRAD)

    !   Sweeps through the spatial grid computing radiances for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  The sweeping
    ! order through the grid points is in SWEEPORD. If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is traced
    ! back to a face containing known (valid) radiances.  As each cell is
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! at the end of the ray path is interpolated from the known grid point
    ! radiances on the face.  The resulting radiances are put in GRIDRAD.

    integer NX, NY, NZ, NA, NPTS, NCELLS, KANG, BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS), XGRID(*), YGRID(*), ZGRID(NZ)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NA, NPTS), GRIDRAD(NPTS)
    integer BITX, BITY, BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer I1, I2, I3, I4, IOPP, JFACE, KFACE, IC, IZ
    integer GRIDFACE(4, 6), OPPFACE(6)
    integer IORDER, ICORNER, JOCTORDER3(8), JOCTORDER2(8)
    logical :: VALIDRAD, VALIDFACE, IPINX, IPINY
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS, path
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
    real(RD_) :: EXT0P, SRCEXT0P
    real(RD_) :: TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANSMIT
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/
    data JOCTORDER3/1, 3, 5, 7, 2, 4, 6, 8/, JOCTORDER2/1, 3, 1, 3, 2, 4, 2, 4/

    integer :: iwhile, ibreak

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
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

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negative.
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
        stop 'Error [_backIntegrate3D_]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate goes to.
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    if (btest(IPFLAG, 1)) then
        JOCT = JOCTORDER2(IOCT)
    else
        JOCT = JOCTORDER3(IOCT)
    end if

    IZ = (NZ - 1)*(1 - BITZ) + 1

    ! Sweep through all the grid points
    do IORDER = 1, NPTS

        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)

        ! this is for mpi (e.g., shdMPI.f90, not OMP)
        !╭────────────────────────────────────────────────╮!
        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if ((GRIDPOS(3, IPT) .lt. ZGRID(IZ) .and. BITZ .eq. 0) .or.&
                &(GRIDPOS(3, IPT) .gt. ZGRID(IZ) .and. BITZ .eq. 1)) then
                IZ = IZ + 2*BITZ - 1
                call shd_mpi__calBoundaryRad(BCFLAG, IPFLAG, JOCT, IZ,&
                    &NX, NY, NZ, XGRID, YGRID, ZGRID,&
                    &NA, NPTS, NCELLS, GRIDPTR,&
                    &NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,&
                    &MU, PHI, EXTINCT, 1, SOURCE,&
                    &KANG, GRIDRAD)
            end if
        end if
        !╰────────────────────────────────────────────────╯!

        if (GRIDRAD(IPT) .lt. 0.0) then

            ! Start off at the grid point
            ICELL = IPCELL
            RAD = 0.0
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            SRCEXT1 = EXT1*SOURCE(KANG, IPT)
            XE = GRIDPOS(1, IPT)
            YE = GRIDPOS(2, IPT)
            ZE = GRIDPOS(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)

                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate3D_]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL, XE, YE, ZE,&
                        &GRIDPOS(1, IPT), GRIDPOS(2, IPT), GRIDPOS(3, IPT)
                    stop
                end if

                IPINX = btest(int(CELLFLAGS(ICELL)), 0)
                IPINY = btest(int(CELLFLAGS(ICELL)), 1)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = GRIDPTR(9 - IOCT, ICELL)

                ! Get the distances to the 3 planes and select the closest
                !╭────────────────────────────────────────────────╮!
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
                !╰────────────────────────────────────────────────╯!

                if (SO .lt. -EPS) then
                    write (6, *) 'Error [_backIntegrate3D_]: SO<0 ',&
                        &MU, PHI, XE, YE, ZE, SOX, SOY, SOZ, SO, IPT, ICELL
                    stop
                end if

                ! Compute the coordinates of the cell exitting location
                XE = XE + SO*CX
                YE = YE + SO*CY
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                if (SOX .le. SOZ .and. SOX .le. SOY) then
                    IFACE = 2 - BITX
                    JFACE = 1
                else if (SOY .le. SOZ) then
                    IFACE = 4 - BITY
                    JFACE = 2
                else
                    IFACE = 6 - BITZ
                    JFACE = 3
                end if

                ! Get the next cell to go to
                INEXTCELL = NEIGHPTR(IFACE, ICELL)
                if (INEXTCELL .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
                end if

                ! Get the face grid pointers
                if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    KFACE = IFACE
                    IC = ICELL
                else
                    ! If going to smaller face then use next face (more accurate)
                    KFACE = OPPFACE(IFACE)
                    IC = INEXTCELL
                end if

                I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
                I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
                I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
                I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)

                ! Compute the face interpolation factors
                if (JFACE .eq. 1) then
                    U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                    if (IPINY) then
                        V = 0.5
                    else
                        V = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I2) - GRIDPOS(2, I1))
                    end if
                else if (JFACE .eq. 2) then
                    U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                    end if
                else
                    if (IPINY) then
                        U = 0.5
                    else
                        U = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I3) - GRIDPOS(2, I1))
                    end if
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                    end if
                end if

                if (U .lt. -0.0001 .or. V .lt. -0.0001 .or.&
                    &U .gt. 1.0001 .or. V .gt. 1.0001) then
                    print *, '_backIntegrate3D_: U,V<0 or U,V>1: ',&
                        &MU, PHI, XE, YE, ZE, U, V, IFACE, JFACE, KFACE,&
                        &GRIDPOS(1, I1), GRIDPOS(2, I1), GRIDPOS(3, I1)
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    if (JFACE .eq. 1) then
                        XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
                    else if (JFACE .eq. 2) then
                        YE = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
                    else
                        ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                    end if
                end if

                ! Interpolate extinction and source function at face intersection
                F1 = (1 - U)*(1 - V)
                F2 = (1 - U)*V
                F3 = U*(1 - V)
                F4 = U*V
                EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)&
                    &+ F3*EXTINCT(I3) + F4*EXTINCT(I4)
                ! Correctly interpolate source using extinction*source
                SRCEXT0 = (F1*SOURCE(KANG, I1)*EXTINCT(I1)&
                    &+ F2*SOURCE(KANG, I2)*EXTINCT(I2)&
                    &+ F3*SOURCE(KANG, I3)*EXTINCT(I3)&
                    &+ F4*SOURCE(KANG, I4)*EXTINCT(I4))
                ! Compute the cell radiance: integration of the source function
                EXT = 0.5d0*(EXT0 + EXT1)

                TAU = EXT*SO
                if (TAU .ge. 0.5) then
                    TRANSCELL = exp(-TAU)
                    ABSCELL = 1.0 - TRANSCELL
                else
                    ABSCELL = TAU*(1.0d0 - 0.5d0*TAU*&
                        &(1.0d0 - 1.0d0/3.0d0*TAU*(1.0d0 - 0.25d0*TAU)))
                    TRANSCELL = 1.0 - ABSCELL
                end if

                if (TAU .le. 2.0) then
                    if (EXT .eq. 0.0) then
                        SRC = 0.0
                    else
                        ! Linear extinction, linear source*extinction, to first order
                        SRC = (0.5d0*(SRCEXT0 + SRCEXT1)&
                            &+ FRAC112_*(EXT0*SRCEXT1 - EXT1*SRCEXT0)*SO)/EXT
                    end if
                else
                    ! Combined first order expansion and constant extinction formula
                    ! SRC = 0.5/EXT *( SRCEXT0+SRCEXT1&
                    !     &+ (EXT0*SRCEXT1-EXT1*SRCEXT0)*SO&
                    !     &*(1-2/TAU+2*TRANSCELL/ABSCELL)/TAU )
                    EXT0P = EXT0
                    SRCEXT0P = SRCEXT0
                    if (TAU .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0d0/TAU
                        if (EXT0 .gt. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
                    end if
                    SRC = 1.0d0/(EXT0P + EXT1)*(SRCEXT0P + SRCEXT1&
                        &+ (EXT0P*SRCEXT1 - EXT1*SRCEXT0P)*2.0d0/(EXT0P + EXT1)&
                        &*(1.0d0 - 2.0d0/TAU + 2.0d0*TRANSCELL/ABSCELL))
                end if
                SRC = max(SRC, 0.0d0)

                ! Add in the cell radiance and update the transmission.
                ! If this is a boundary or the transmission is below the
                ! cutoff and we have a valid radiance then set the flag
                ! to stop the tracing and add in the interpolated face radiance.
                RAD = RAD + TRANSMIT*SRC*ABSCELL
                TRANSMIT = TRANSMIT*TRANSCELL
                if (RAD < -1.0e-5) then
                    print '(A,3(1X,F6.4),8(1X,E12.5))',&
                        &'_backIntegrate3D_: RAD<0 ',&
                        &XE, YE, ZE, EXT0, EXT1, SO, TAU, TRANSMIT, SRC, ABSCELL, RAD
                    call shd_mpi__abort('The end')
                end if

                VALIDFACE = (GRIDRAD(I1) .ge. -0.1 .and. GRIDRAD(I2) .ge. -0.1&
                    &.and. GRIDRAD(I3) .ge. -0.1 .and. GRIDRAD(I4) .ge. -0.1)
                if (INEXTCELL .le. 0 .or.&
                    &(TRANSMIT .le. TRANSMIN .and. VALIDFACE)) then
                    if (VALIDFACE) then
                        VALIDRAD = .true.
                        RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2)&
                            &+ F3*GRIDRAD(I3) + F4*GRIDRAD(I4)
                        RAD = RAD + TRANSMIT*RAD0
                    else
                        print *, 'Error [_backIntegrate3D_]: INEXTCELL=0: ', xe, ye, ze, icell
                        print *, iz, joct, mu, phi
                        print '(i7,4(1x,f7.4))', i1, GRIDRAD(I1),&
                            &gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                        print '(i7,4(1x,f7.4))', i2, GRIDRAD(I2),&
                            &gridpos(1, i2), gridpos(2, i2), gridpos(3, i2)
                        print '(i7,4(1x,f7.4))', i3, GRIDRAD(I3),&
                            &gridpos(1, i3), gridpos(2, i3), gridpos(3, i3)
                        print '(i7,4(1x,f7.4))', i4, GRIDRAD(I4),&
                            &gridpos(1, i4), gridpos(2, i4), gridpos(3, i4)
                        print '(i7,8x,3(1x,f7.4))', ipt,&
                            &gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                        stop
                    end if
                else
                    EXT1 = EXT0
                    SRCEXT1 = SRCEXT0
                    ICELL = INEXTCELL
                end if
            end do

            GRIDRAD(IPT) = RAD

        end if

    end do

    return

end subroutine shd_intg__backIntegrate3D_

subroutine shd_intg__backIntegrate2D(&
    &NX, NY, NZ, NA, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, SWEEPORD, MU, PHI, TRANSMIN,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)

    !   Sweeps through the spatial grid computing radiances for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is traced
    ! back to a face containing known (valid) radiances.  As each cell is
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! ad the end of the ray path is interpolated from the known grid point
    ! radiances on the face.  The resulting radiances are put in GRIDRAD.

    integer NX, NY, NZ, NA, NPTS, NCELLS, KANG, NSTOKES
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer I1, I2, IOPP, JFACE, KFACE, IC, IZ
    integer GRIDFACE(2, 6), OPPFACE(6)
    integer IORDER, ICORNER, JOCTORDER(8)
    logical :: VALIDRAD, VALIDFACE, IPINX
    real(RD_) :: CX, CZ, CXINV, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOZ, EPS
    real(RD_) :: U, F1, F2
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(RD_) :: TAU, TRANSCELL, ABSCELL, TRANSMIT
    real(RD_) :: RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/1, 5, 2, 6, 0, 0, 0, 0, 1, 2, 5, 6/
    data OPPFACE/2, 1, 4, 3, 6, 5/
    data JOCTORDER/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    CX = sqrt(1.0 - MU**2)*cos(PHI + DPI_)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negative.
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [_backIntegrate2D]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate goes to.
    IOCT = 1 + BITX + 4*BITZ
    JOCT = JOCTORDER(IOCT)

    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)
        if (GRIDRAD(1, IPT) .lt. 0.0) then

            ! Start off at the grid point
            ICELL = IPCELL
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
            XE = GRIDPOS(1, IPT)
            YE = GRIDPOS(2, IPT)
            ZE = GRIDPOS(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)
                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate2D]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL
                    stop
                end if
                IPINX = btest(int(CELLFLAGS(ICELL)), 0)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = GRIDPTR(9 - IOCT, ICELL)
                ! Get the distances to the 3 planes and select the closest
                if (IPINX) then
                    SOX = 1.0e20
                else
                    SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
                end if
                SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
                SO = min(SOX, SOZ)
                if (SO .lt. -EPS) then
                    write (6, *) 'Error [_backIntegrate2D]: SO<0 ',&
                        &MU, PHI, XE, ZE, SOX, SOZ, SO, IPT, ICELL
                    stop
                end if
                ! Compute the coordinates of the cell exitting location
                XE = XE + SO*CX
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                if (SOX .le. SOZ) then
                    IFACE = 2 - BITX
                    JFACE = 1
                else
                    IFACE = 6 - BITZ
                    JFACE = 3
                end if
                ! Get the next cell to go to
                INEXTCELL = NEIGHPTR(IFACE, ICELL)
                if (INEXTCELL .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
                end if

                ! Get the face grid pointers
                if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    KFACE = IFACE
                    IC = ICELL
                else
                    ! If going to smaller face then use next face (more accurate)
                    KFACE = OPPFACE(IFACE)
                    IC = INEXTCELL
                end if
                I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
                I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
                ! Compute the face interpolation factors
                if (JFACE .eq. 1) then
                    U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I2) - GRIDPOS(3, I1))
                else
                    if (IPINX) then
                        U = 0.5
                    else
                        U = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                    end if
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    if (JFACE .eq. 1) then
                        XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
                    else
                        ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                    end if
                end if

                ! Interpolate extinction and source function at face intersection
                F1 = 1 - U
                F2 = U
                EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)
                ! Correctly interpolate source using extinction*source
                SRCEXT0(:) = (F1*SOURCE(:, KANG, I1)*EXTINCT(I1)&
                    &+ F2*SOURCE(:, KANG, I2)*EXTINCT(I2))
                ! Compute the cell radiance: integration of the source function
                EXT = 0.5*(EXT0 + EXT1)
                TAU = EXT*SO
                if (TAU .ge. 0.5) then
                    TRANSCELL = exp(-TAU)
                    ABSCELL = 1.0 - TRANSCELL
                else
                    ABSCELL = TAU*(1.0 - 0.5*TAU*&
                        &(1.0 - FRAC13_*TAU*(1 - 0.25*TAU)))
                    TRANSCELL = 1.0 - ABSCELL
                end if
                if (TAU .le. 2.0) then
                    if (EXT .eq. 0.0) then
                        SRC(:) = 0.0
                    else
                        ! Linear extinction, linear source*extinction, to first order
                        SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + FRAC112_&
                            &*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                    end if
                else
                    ! Combined first order expansion and constant extinction formula
                    EXT0P = EXT0
                    SRCEXT0P(:) = SRCEXT0(:)
                    if (TAU .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                        if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                    end if
                    SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                        &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                        &*(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
                end if
                SRC(1) = max(SRC(1), 0.0d0)

                ! Add in the cell radiance and update the transmission.
                ! If this is a boundary or the transmission is below the
                ! cutoff and we have a valid radiance then set the flag
                ! to stop the tracing and add in the interpolated face radiance
                RAD(:) = RAD(:) + TRANSMIT*SRC(:)*ABSCELL
                TRANSMIT = TRANSMIT*TRANSCELL
                VALIDFACE = (GRIDRAD(1, I1) .ge. -0.1 .and. GRIDRAD(1, I2) .ge. -0.1)
                if (INEXTCELL .le. 0 .or.&
                    &(TRANSMIT .le. TRANSMIN .and. VALIDFACE)) then
                    if (VALIDFACE) then
                        VALIDRAD = .true.
                        RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)
                        RAD(:) = RAD(:) + TRANSMIT*RAD0(:)
                    else
                        print *, 'Error [_backIntegrate2D]: INEXTCELL=0: ', xe, ye, ze, icell
                        print *, iz, joct, mu, phi
                        print '(i7,4(1x,f7.4))', i1, GRIDRAD(1, I1),&
                            &gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                        print '(i7,4(1x,f7.4))', i2, GRIDRAD(1, I2),&
                            &gridpos(1, i2), gridpos(2, i2), gridpos(3, i2)
                        print '(i7,8x,3(1x,f7.4))', ipt,&
                            &gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                        stop
                    end if
                else
                    EXT1 = EXT0
                    SRCEXT1(:) = SRCEXT0(:)
                    ICELL = INEXTCELL
                end if
            end do

            GRIDRAD(:, IPT) = RAD(:)

        end if

    end do

    return

end subroutine shd_intg__backIntegrate2D

subroutine shd_intg__backIntegrate1D(&
    &NX, NY, NZ, NA, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, SWEEPORD, MU, PHI, TRANSMIN,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)

    !   Sweeps through the spatial grid computing radiances for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is traced
    ! back to a face containing known (valid) radiances.  As each cell is
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! ad the end of the ray path is interpolated from the known grid point
    ! radiances on the face.  The resulting radiances are put in GRIDRAD.

    integer NX, NY, NZ, NA, NPTS, NCELLS, KANG, NSTOKES
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)
    integer BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE, I1, IOPP
    integer GRIDFACE(6)
    integer IORDER, ICORNER, JOCTORDER(8)
    logical VALIDRAD
    real(RD_) :: CZ, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, EPS
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(RD_) :: TAU, TRANSCELL, ABSCELL, TRANSMIT
    real(RD_) :: RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/0, 0, 0, 0, 1, 5/
    data JOCTORDER/1, 1, 1, 1, 2, 2, 2, 2/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    CZ = -MU
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negative.
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [_backIntegrate1D]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate goes to.
    IOCT = 1 + 4*BITZ
    JOCT = JOCTORDER(IOCT)

    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)
        if (GRIDRAD(1, IPT) .lt. 0.0) then
            ! Start off at the grid point
            ICELL = IPCELL
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
            XE = GRIDPOS(1, IPT)
            YE = GRIDPOS(2, IPT)
            ZE = GRIDPOS(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)
                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate1D]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL
                    stop
                end if
                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = GRIDPTR(9 - IOCT, ICELL)
                ! Get the distances to the 3 planes and select the closest
                SO = (GRIDPOS(3, IOPP) - ZE)*CZINV
                if (SO .lt. -EPS) then
                    write (6, *) 'Error [_backIntegrate1D]: SO<0 ',&
                        &MU, PHI, ZE, SO, IPT, ICELL
                    stop
                end if
                ! Compute the coordinates of the cell exitting location
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                IFACE = 6 - BITZ
                ! Get the next cell to go to
                INEXTCELL = NEIGHPTR(IFACE, ICELL)

                ! Get the face grid pointers
                I1 = GRIDPTR(GRIDFACE(IFACE), ICELL)
                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                end if

                ! Get extinction and source function at face intersection
                EXT0 = EXTINCT(I1)
                SRCEXT0(:) = EXTINCT(I1)*SOURCE(:, KANG, I1)

                ! this is added to address layer averaging issue
                ! a more robust criteria should be developed in the future
                !   to distinguish the atmopheric medium intersection
                !╭────────────────────────────────────────────────╮!
                if (abs(EXT0-EXT1) .gt. min(EXT0, EXT1)*10.0) then
                    if (BITZ .eq. 0) then
                        EXT0 = EXT1
                        SRCEXT0(:) = SRCEXT1(:)
                    else
                        EXT1 = EXT0
                        SRCEXT1(:) = SRCEXT0(:)
                    endif
                end if
                !╰────────────────────────────────────────────────╯!

                ! Compute the cell radiance: integration of the source function
                EXT = 0.5*(EXT0 + EXT1)
                TAU = EXT*SO
                if (TAU .ge. 0.5) then
                    TRANSCELL = exp(-TAU)
                    ABSCELL = 1.0 - TRANSCELL
                else
                    ABSCELL = TAU*(1.0 - 0.5*TAU*&
                        &(1.0 - FRAC13_*TAU*(1 - 0.25*TAU)))
                    TRANSCELL = 1.0 - ABSCELL
                end if
                if (TAU .le. 2.0) then
                    if (EXT .eq. 0.0) then
                        SRC(:) = 0.0
                    else
                        SRC(:) = (0.5d0*(SRCEXT0(:) + SRCEXT1(:)) + FRAC112_&
                            &*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                    end if
                else
                    ! Combined first order expansion and constant extinction formula
                    EXT0P = EXT0
                    SRCEXT0P(:) = SRCEXT0(:)
                    if (TAU .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                        if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                    end if
                    SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                        &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                        &*(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
                end if
                SRC(1) = max(SRC(1), 0.0d0)

                ! See if there are valid radiances for the intersection face:
                ! If so, set flag, interpolate radiance, and add to
                !   radiance from source function in cell
                ! If not, only add in cell radiance and prepare for next cell
                if (GRIDRAD(1, I1) .ge. -0.1) then
                    VALIDRAD = .true.
                    RAD0(:) = GRIDRAD(:, I1)
                    RAD(:) = RAD(:) + TRANSMIT*(RAD0(:)*TRANSCELL + SRC(:)*ABSCELL)
                else
                    RAD(:) = RAD(:) + TRANSMIT*SRC(:)*ABSCELL
                    SRCEXT1(:) = SRCEXT0(:)
                    TRANSMIT = TRANSMIT*TRANSCELL
                    EXT1 = EXT0
                    ICELL = INEXTCELL
                end if
            end do

            GRIDRAD(:, IPT) = RAD(:)
        end if
    end do

    return

end subroutine shd_intg__backIntegrate1D

end module shd_intg
