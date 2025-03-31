subroutine NEW_GRIDS(&
    &BCFLAG, GRIDTYPE, NPX, NPY, NPZ, NX, NY, NZ,&
    &XSTART, YSTART, DELXP, DELYP, ZLEVELS,&
    &XGRID, YGRID, ZGRID)

    ! Makes the XGRID, YGRID, and ZGRID arrays from the input Z levels
    !   and X and Y spacing from the property file according to the GRIDTYPE.

    implicit none
    integer BCFLAG, NPX, NPY, NPZ, NX, NY, NZ
    real XSTART, YSTART, DELXP, DELYP, ZLEVELS(NPZ)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    character(len=1) :: GRIDTYPE
    integer IX, IY, IZ

    ! Make the X, Y, and Z grids, depending on type:
    !   E for evenly spaced, P from property file, F from zgrid.inp
    !   Only Even allowed horizontally
    if (DELXP .le. 0.0) DELXP = 1.0
    if (btest(BCFLAG, 2)) then
        do IX = 1, NX
            if (NPX .gt. 1) then
                XGRID(IX) = XSTART + (IX - 1)*(DELXP*(NPX - 1))/FLOAT(NX - 1)
            else
                XGRID(IX) = XSTART + (IX - 1)*(DELXP*NPX)/FLOAT(NX - 1)
            end if
        end do
    else
        do IX = 1, NX + 1
            XGRID(IX) = XSTART + (IX - 1)*(DELXP*NPX)/FLOAT(NX)
        end do
    end if

    if (DELYP .le. 0.0) DELYP = 1.0
    if (btest(BCFLAG, 3)) then
        do IY = 1, NY
            if (NPY .gt. 1) then
                YGRID(IY) = YSTART + (IY - 1)*(DELYP*(NPY - 1))/FLOAT(NY - 1)
            else
                YGRID(IY) = YSTART + (IY - 1)*(DELYP*NPY)/FLOAT(NY - 1)
            end if
        end do
    else
        do IY = 1, NY + 1
            YGRID(IY) = YSTART + (IY - 1)*(DELYP*NPY)/FLOAT(NY)
        end do
    end if

    if (GRIDTYPE .eq. 'P') then
        if (NZ .ne. NPZ) then
            write (6, *)&
                &'Error [NEW_GRIDS]: For internal Z grid same as property file grid,'
            write (6, *) '  must have same number of grid points (NZ=NPZ).'
            stop
        end if
        do IZ = 1, NZ
            ZGRID(IZ) = ZLEVELS(IZ)
        end do
    else if (GRIDTYPE .eq. 'E') then
        do IZ = 1, NZ
            ZGRID(IZ) = ZLEVELS(1)&
                &+ (IZ - 1)*(ZLEVELS(NPZ) - ZLEVELS(1))/FLOAT(NZ - 1)
        end do
    else if (GRIDTYPE .eq. 'F') then
        open (UNIT=3, FILE='zgrid.inp', STATUS='OLD')
        do IZ = 1, NZ
            read (3, *) ZGRID(IZ)
        end do
        close (3)
    else
        stop 'Error [NEW_GRIDS]: Illegal grid type in Z.'
    end if

    return
end subroutine

subroutine INIT_CELL_STRUCTURE(&
    &BCFLAG, IPFLAG,&
    &NX, NY, NZ, NX1, NY1, NPTS, NCELLS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)

    ! Make the initial grid cell data structure for the base grid.
    ! For periodic boundary conditions there are NX*NY*(NZ-1) grid cells
    ! and (NX+1)*(NY+1)*NZ grid points.  The IX=NX+1 column is the same
    ! as the IX=1 column, because of the periodicity, but the grid points
    ! are repeated.  The last grid cell on the "right" has the last column
    ! of grid points on the "left" and the repeated "virtual" set of grid
    ! points on the "right".  For open boundary conditions in X and Y
    ! there are (NX+1)*(NY+1)*(NZ-1) grid cells and NX*NY*NZ grid points.
    ! Open boundary conditions are done with the lowest "left most" and
    ! "right most" X and Y columns being independent pixels with
    ! their own grid cells (zero width), therefore one extra grid cell
    ! is required.   Multiple processor boundary conditions have
    ! (NX-1)*(NY-1)*(NZ-1) grid cells and NX*NY*NZ grid points.
    ! Independent pixel mode has NX*NY*(NZ-1) grid cells and NX*NY*NZ
    ! grid points
    !   Variable       Description
    ! IPFLAG         Independent pixel flag - 0 for none, 1 for X, 2 for Y,
    !                  3 for X&Y.
    ! BCFLAG         Boundary condition flag - 0 for periodic, 1 for open in X,
    !                  2 for open in Y, 3 for open in X&Y.  4 for multiple
    !                  processor in X, 8 for MP in Y, 12 for MP in both.
    ! NPTS           total number of grid points
    ! NCELLS         total number of cells
    ! GRIDPOS(3,IP)  x,y,z coordinate values of grid point locations
    ! GRIDPTR(8,IC)  pointer to 8 grid points belonging to a cell.
    !                  Order of pointers is based on binary digits zyx:
    !                  IOCT = 1 + BITX + 2*BITY + 4*BITZ, where
    !                  BITc is 0 for negative side and 1 for positive side.
    ! NEIGHPTR(6,IC) pointer to 6 neighboring cells of a cell.
    !                  Order of pointers are -X,+X,-Y,+Y,-Z,+Z (1-6).
    !                  Is 0 if there is no neighbor cell (boundary).
    !                  Is positive if there is just one cell across face
    !                  of current cell (unique neighbor).  Is negative
    !                  (-index) if there is more than one cell across
    !                  the face of current cell.  In this case the
    !                  cell pointed to is the parent cell that has a
    !                  has same or larger face adjoining current cell.
    ! TREEPTR(2,IC)  pointer to parent and children cells in tree structure.
    !                  First is pointer to parent cell, 0 for a base cell.
    !                  Second is pointer to first of the two children cells.
    !                  This pointer is 0 if there are no children cells.
    !                  The children cells are in successive order, with
    !                  the negative side first.  See CELLFLAGS for direction
    !                  of division of split.
    ! CELLFLAGS(IC)  a collection of binary flags for the cells:
    !                  bit  For bit set
    !                   0  Do independent pixel for X
    !                   1  Do independent pixel for Y
    !                  2,3 Direction cell is split (0=no split, 1=X, 2=Y, 3=Z)
    ! IP is a grid point index, IC is a cell index.

    implicit none
    integer IPFLAG, BCFLAG, NX, NY, NZ, NX1, NY1, NPTS, NCELLS
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real XGRID(NX1), YGRID(NY1), ZGRID(NZ), GRIDPOS(3, *)
    integer IX, IY, IZ, IX0, IX1, IY0, IY1, I, NXC, NYC
    logical BTEST

    ! First, make the grid point positions
    I = 0
    do IX = 1, NX1
        do IY = 1, NY1
            do IZ = 1, NZ
                I = I + 1
                GRIDPOS(1, I) = XGRID(IX)
                GRIDPOS(2, I) = YGRID(IY)
                GRIDPOS(3, I) = ZGRID(IZ)
            end do
        end do
    end do
    NPTS = I

    NXC = NX
    if (btest(BCFLAG, 0)) NXC = NX + 1
    if (btest(BCFLAG, 2) .and. .not. btest(IPFLAG, 0)) NXC = NX - 1
    NYC = NY
    if (btest(BCFLAG, 1)) NYC = NY + 1
    if (btest(BCFLAG, 3) .and. .not. btest(IPFLAG, 1)) NYC = NY - 1

    ! Then setup the cell structures
    I = 0
    do IX = 1, NXC
        if (btest(IPFLAG, 0)) then
            IX0 = IX
            IX1 = IX
        else if (btest(BCFLAG, 0)) then
            IX0 = max(1, IX - 1)
            IX1 = min(NX, IX)
        else
            IX0 = IX
            IX1 = IX + 1
        end if
        do IY = 1, NYC
            if (btest(IPFLAG, 1)) then
                IY0 = IY
                IY1 = IY
            else if (btest(BCFLAG, 1)) then
                IY0 = max(1, IY - 1)
                IY1 = min(NY, IY)
            else
                IY0 = IY
                IY1 = IY + 1
            end if
            do IZ = 1, NZ - 1
                I = I + 1
                GRIDPTR(1, I) = IZ + NZ*(IY0 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(2, I) = IZ + NZ*(IY0 - 1) + NZ*NY1*(IX1 - 1)
                GRIDPTR(3, I) = IZ + NZ*(IY1 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(4, I) = IZ + NZ*(IY1 - 1) + NZ*NY1*(IX1 - 1)
                GRIDPTR(5, I) = IZ + 1 + NZ*(IY0 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(6, I) = IZ + 1 + NZ*(IY0 - 1) + NZ*NY1*(IX1 - 1)
                GRIDPTR(7, I) = IZ + 1 + NZ*(IY1 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(8, I) = IZ + 1 + NZ*(IY1 - 1) + NZ*NY1*(IX1 - 1)
                CELLFLAGS(I) = 0
                ! Initialize the neighbor pointer; to self for IP mode or
                !   if open boundary conditions and at either end
                if (btest(IPFLAG, 0)) then
                    NEIGHPTR(1, I) = I
                    NEIGHPTR(2, I) = I
                    CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 0)
                else if (btest(BCFLAG, 0)) then
                    if (IX .eq. 1) then
                        NEIGHPTR(1, I) = I
                        NEIGHPTR(2, I) = I + (NZ - 1)*NYC
                        CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 0)
                    else if (IX .eq. NXC) then
                        NEIGHPTR(1, I) = I - (NZ - 1)*NYC
                        NEIGHPTR(2, I) = I
                        CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 0)
                    else
                        NEIGHPTR(1, I) = I - (NZ - 1)*NYC
                        NEIGHPTR(2, I) = I + (NZ - 1)*NYC
                    end if
                else if (btest(BCFLAG, 2)) then
                    if (IX .eq. 1) then
                        NEIGHPTR(1, I) = 0
                        NEIGHPTR(2, I) = I + (NZ - 1)*NYC
                    else if (IX .eq. NXC) then
                        NEIGHPTR(1, I) = I - (NZ - 1)*NYC
                        NEIGHPTR(2, I) = 0
                    else
                        NEIGHPTR(1, I) = I - (NZ - 1)*NYC
                        NEIGHPTR(2, I) = I + (NZ - 1)*NYC
                    end if
                else
                    if (IX .eq. 1) then
                        NEIGHPTR(1, I) = I + (NXC - 1)*(NZ - 1)*NYC
                    else
                        NEIGHPTR(1, I) = I - (NZ - 1)*NYC
                    end if
                    if (IX .eq. NX) then
                        NEIGHPTR(2, I) = I - (NXC - 1)*(NZ - 1)*NYC
                    else
                        NEIGHPTR(2, I) = I + (NZ - 1)*NYC
                    end if
                end if
                if (btest(IPFLAG, 1)) then
                    NEIGHPTR(3, I) = I
                    NEIGHPTR(4, I) = I
                    CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 1)
                else if (btest(BCFLAG, 1)) then
                    if (IY .eq. 1) then
                        NEIGHPTR(3, I) = I
                        NEIGHPTR(4, I) = I + (NZ - 1)
                        CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 1)
                    else if (IY .eq. NYC) then
                        NEIGHPTR(3, I) = I - (NZ - 1)
                        NEIGHPTR(4, I) = I
                        CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 1)
                    else
                        NEIGHPTR(3, I) = I - (NZ - 1)
                        NEIGHPTR(4, I) = I + (NZ - 1)
                    end if
                else if (btest(BCFLAG, 3)) then
                    if (IY .eq. 1) then
                        NEIGHPTR(3, I) = 0
                        NEIGHPTR(4, I) = I + (NZ - 1)
                    else if (IY .eq. NYC) then
                        NEIGHPTR(3, I) = I - (NZ - 1)
                        NEIGHPTR(4, I) = 0
                    else
                        NEIGHPTR(3, I) = I - (NZ - 1)
                        NEIGHPTR(4, I) = I + (NZ - 1)
                    end if
                else
                    if (IY .eq. 1) then
                        NEIGHPTR(3, I) = I + (NYC - 1)*(NZ - 1)
                    else
                        NEIGHPTR(3, I) = I - (NZ - 1)
                    end if
                    if (IY .eq. NY) then
                        NEIGHPTR(4, I) = I - (NYC - 1)*(NZ - 1)
                    else
                        NEIGHPTR(4, I) = I + (NZ - 1)
                    end if
                end if
                if (IZ .eq. 1) then
                    NEIGHPTR(5, I) = 0
                else
                    NEIGHPTR(5, I) = I - 1
                end if
                if (IZ .eq. NZ - 1) then
                    NEIGHPTR(6, I) = 0
                else
                    NEIGHPTR(6, I) = I + 1
                end if
                ! Base grid cells have no parents or children
                TREEPTR(1, I) = 0
                TREEPTR(2, I) = 0
            end do
        end do
    end do
    NCELLS = I
    return
end subroutine

subroutine INTERP_GRID(&
    &NPTS, NSTLEG, NLEG, GRIDPOS,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)

    ! Calls TRILIN_INTERP_PROP to interpolate the input arrays from
    !   the property grid to each internal grid point.

    implicit none
    integer NPTS, NSTLEG, NLEG
    integer*2 IPHASE(*)
    real GRIDPOS(3, NPTS)
    real TEMP(*), EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    integer IP

    ! Initialize: transfer the tabulated phase functions
    call TRILIN_INTERP_PROP(0.0, 0.0, 0.0, .true., NSTLEG, NLEG,&
        &TEMP(1), EXTINCT(1), ALBEDO(1), LEGEN(1, 0, 1), IPHASE(1))

    ! Trilinearly interpolate from the property grid to the adaptive
    do IP = 1, NPTS
        call TRILIN_INTERP_PROP&
            &(GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &.false., NSTLEG, NLEG, TEMP(IP), EXTINCT(IP),&
            &ALBEDO(IP), LEGEN(1, 0, IP), IPHASE(IP))
    end do

    return
end subroutine

subroutine MAKE_DIRECT(NPTS, BCFLAG, IPFLAG, DELTAM,&
    &ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
    &GRIDPOS, DIRFLUX)

    ! Makes the direct beam solar flux for the internal base grid.
    !   DIRFLUX is set to F*exp(-tau_sun).
    !   Actually calls DIRECT_BEAM_PROP to do all the hard work.

    implicit none
    integer NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), DIRFLUX(NPTS)
    integer SIDE, IP
    logical VALIDBEAM
    real UNIFZLEV, XO, YO, ZO, DIR, DIRPATH

    call DIRECT_BEAM_PROP(1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
        &DELTAM, ML, NSTLEG, NLEG,&
        &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)

    do IP = 1, NPTS
        DIRPATH = 0.0
        call DIRECT_BEAM_PROP&
            &(0, GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &BCFLAG, IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
            &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(IP),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end do

    return
end subroutine

subroutine PREPARE_PROP(ML, MM, NSTLEG, NLEG, NPTS, DELTAM,&
    &NUMPHASE, SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX,&
    &EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)

    ! Prepares the grid arrays for the iterative solution process.
    ! If doing Delta-M scaling then the extinction, albedo, and Legend
    ! terms are scaled first; only the 0 to ML LEGEN terms are scaled.
    ! Outputs PLANCK with (1-omega)*B(T) for thermal source, where B(T
    ! the Planck function (meaning depends on UNITS).
    ! TEMP array is unchanged.

    implicit none
    integer ML, MM, NSTLEG, NLEG, NPTS, NUMPHASE
    integer*2 IPHASE(NPTS)
    logical DELTAM
    real WAVENO(2), WAVELEN, ALBMAX
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real TEMP(NPTS), PLANCK(NPTS)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IPH, L
    real F, BB

    ALBMAX = 0.0

    if (DELTAM .and. NUMPHASE .gt. 0) then
        do IPH = 1, NUMPHASE
            F = LEGEN(1, ML + 1, IPH)
            do L = 0, ML
                ! Scale the diagonal and off-diagonal phase matrix elements d
                LEGEN(1, L, IPH) = (LEGEN(1, L, IPH) - F)/(1 - F)
                if (NSTLEG .gt. 1) then
                    LEGEN(2, L, IPH) = (LEGEN(2, L, IPH) - F)/(1 - F)
                    LEGEN(3, L, IPH) = (LEGEN(3, L, IPH) - F)/(1 - F)
                    LEGEN(4, L, IPH) = (LEGEN(4, L, IPH) - F)/(1 - F)
                    LEGEN(5, L, IPH) = LEGEN(5, L, IPH)/(1 - F)
                    LEGEN(6, L, IPH) = LEGEN(6, L, IPH)/(1 - F)
                end if
            end do
        end do
    end if

    do I = 1, NPTS
        if (DELTAM) then
            if (NUMPHASE .gt. 0) then
                F = LEGEN(1, ML + 1, IPHASE(I))
            else
                ! NUMPHASE=0 is standard property format, which is unpolariz
                F = LEGEN(1, ML + 1, I)
                do L = 0, ML
                    LEGEN(1, L, I) = (LEGEN(1, L, I) - F)/(1 - F)
                end do
            end if
            EXTINCT(I) = (1.0 - ALBEDO(I)*F)*EXTINCT(I)
            ALBEDO(I) = (1.0 - F)*ALBEDO(I)/(1.0 - ALBEDO(I)*F)
        end if
        ALBMAX = max(ALBMAX, ALBEDO(I))
        if (SRCTYPE .ne. 'S') then
            call PLANCK_FUNCTION(TEMP(I), UNITS, WAVENO, WAVELEN, BB)
            PLANCK(I) = (1.0 - ALBEDO(I))*BB
        end if
    end do
    return
end subroutine

subroutine INIT_RADIANCE(NSTOKES, NXY, NZ, NSTLEG, NLEG,&
    &RSHPTR, ZGRID,&
    &EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP,&
    &SKYRAD, UNITS, WAVENO, WAVELEN, RADIANCE)
    ! Initializes radiance field by solving plane-parallel two-stream.
    ! Solves the L=1 M=0 SH system by transforming the pentadiagonal
    ! system to tridiagonal and calling solver.
    ! Does the initialization for the NXY columns of the base grid.
    ! Only the I Stokes parameter is initialized, the rest are zeroed.
    implicit none
    integer NSTOKES, NXY, NZ, NSTLEG, NLEG, NUMPHASE, RSHPTR(*)
    integer*2 IPHASE(NZ, NXY)
    real GNDALBEDO, GNDTEMP, SKYRAD, SOLARFLUX, SOLARMU
    real WAVENO(2), WAVELEN
    real ZGRID(NZ)
    real EXTINCT(NZ, NXY), ALBEDO(NZ, NXY), LEGEN(NSTLEG, 0:NLEG, *)
    real TEMP(NZ, NXY)
    real RADIANCE(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer NLAYER, I, IZ, IR, J, K, L
    logical DELTAM
    real PI, C0, C1
    real EXT0, EXT1, SCAT0, SCAT1, G0, G1, GNDEMIS
    real, allocatable :: OPTDEPTHS(:), ALBEDOS(:), ASYMMETRIES(:)
    real, allocatable :: TEMPS(:), FLUXES(:, :)
    character(len=1) :: SRCT

    allocate (OPTDEPTHS(NZ), ALBEDOS(NZ), ASYMMETRIES(NZ))
    allocate (TEMPS(NZ), FLUXES(3, NZ))

    NLAYER = NZ - 1
    DELTAM = .false.
    SRCT = SRCTYPE
    if (SRCT .eq. 'B') SRCT = 'T'
    PI = acos(-1.0)
    C0 = sqrt(1.0/PI)
    C1 = sqrt(3.0/(4*PI))

    ! Loop over all the columns in the base grid
    IR = 0
    RSHPTR(1) = 0
    do I = 1, NXY
        ! Make layer properties for the Eddington routine
        do IZ = 1, NLAYER
            L = NZ - IZ
            EXT0 = EXTINCT(IZ, I)
            EXT1 = EXTINCT(IZ + 1, I)
            SCAT0 = ALBEDO(IZ, I)*EXT0
            SCAT1 = ALBEDO(IZ + 1, I)*EXT1
            OPTDEPTHS(L) = (ZGRID(IZ + 1) - ZGRID(IZ))*(EXT0 + EXT1)/2
            if (EXT0 + EXT1 .gt. 0.0) then
                ALBEDOS(L) = (SCAT0 + SCAT1)/(EXT0 + EXT1)
            else
                ALBEDOS(L) = 0.0
            end if
            if (NUMPHASE .gt. 0) then
                J = IPHASE(IZ, I)
                G0 = LEGEN(1, 1, J)
                J = IPHASE(IZ + 1, I)
                G1 = LEGEN(1, 1, J)
            else
                J = NZ*(I - 1) + IZ
                G0 = LEGEN(1, 1, J)
                G1 = LEGEN(1, 1, J + 1)
            end if
            if (SCAT0 + SCAT1 .gt. 0.0) then
                ASYMMETRIES(L) = (SCAT0*G0 + SCAT1*G1)/(SCAT0 + SCAT1)
            else
                ASYMMETRIES(L) = 0.0
            end if
            TEMPS(L + 1) = TEMP(IZ, I)
        end do
        TEMPS(1) = TEMP(NZ, I)
        GNDEMIS = 1.0 - GNDALBEDO
        ! Call the Eddington flux routine
        call EDDRTF(NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES,&
            &TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU,&
            &GNDTEMP, GNDEMIS, SKYRAD, UNITS, WAVENO, WAVELEN,&
            &FLUXES)
        ! Convert fluxes to first two moments of spherical harmonics
        do IZ = 1, NZ
            L = NZ + 1 - IZ
            RADIANCE(:, IR + 1) = 0.0
            RADIANCE(:, IR + 2) = 0.0
            RADIANCE(:, IR + 3) = 0.0
            RADIANCE(:, IR + 4) = 0.0
            RADIANCE(1, IR + 1) = C0*(FLUXES(1, L) + FLUXES(2, L))
            RADIANCE(1, IR + 3) = C1*(FLUXES(1, L) - FLUXES(2, L))
            IR = IR + 4
            RSHPTR(IZ + NZ*(I - 1) + 1) = IR
        end do
    end do
    deallocate (OPTDEPTHS, ALBEDOS, ASYMMETRIES, TEMPS, FLUXES)
    return
end subroutine

subroutine EDDRTF(NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES,&
    &TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU,&
    &GNDTEMP, GNDEMIS, SKYRAD, UNITS, WAVENO, WAVELEN,&
    &FLUXES)
    ! EDDRTF computes the layer interface fluxes for a plane-parallel
    ! atmosphere with either solar or thermal source of radiation using
    ! the Eddington approximation.  The medium is specified by a number
    ! of homogeneous layers.  For a thermal source the Planck function is
    ! linear with optical depth, while for a solar source it is exponential.
    ! The temperatures, optical depth, single scattering albedo, and
    ! asymmetry parameter are specified for each layer.  The boundary
    ! conditions such as the solar flux, and reflection and/or emission from
    ! ground surface are also specified. Delta Eddington scaling may be
    ! used.  The diffuse Eddington fluxes and the solar direct flux at
    ! each level are returned.
    ! The model works by calculating the reflection, transmission, and
    ! source terms for each layer from the input properties.  A
    ! tri-diagonal matrix solver is then used to compute the diffuse fluxes
    ! at each layer interface from the applied boundary conditions.
    !
    ! Parameters:
    !   Input:
    ! NLAYER         integer      Number of homogenous layers
    !                              (layers are specified from the top do
    ! OPTDEPTHS      real array   Optical thickness of layers
    ! ALBEDOS        real array   Single scattering albedos
    ! ASYMMETRIES    real array   Asymmetry parameters
    ! TEMPS          real array   Temperatures (K) at layer interfaces
    !                              (e.g. TEMPS(1) is at top of top layer
    !                               TEMPS(2) is at bottom of top layer).
    ! DELTAM         logical      True for delta-Eddington scaling
    ! SRCTYPE        character    'S' for solar source, 'T' for thermal source
    ! SOLARFLUX      real         Incident solar flux on horizonal plane
    ! SOLARMU        real         Cosine of the solar zenith angle
    ! GNDTEMP        real         Ground temperature (Kelvin)
    ! GNDEMIS        real         Ground emissivity (1-albedo)
    ! SKYRAD         real         Radiance (for solar) or brightness
    !                              temperature (for thermal) of isotropic
    !                              incident radiation from above
    ! UNITS          character    Units of flux for thermal source
    !                             'T' for brightness temperature,
    !                             'R' for W/(m^2 micron),
    !                             'B' for W/(m^2)  for band integration.
    ! WAVENO(2)      real         Wavenumber range (cm^-1) (UNITS='B')
    ! WAVELEN        real         Wavelength (micron) (UNITS='R')
    !
    !   Output:
    ! FLUXES         real         Eddington fluxes at layer interfaces.
    !                               FLUXES(1,L) is upwelling diffuse,
    !                               FLUXES(2,L) is downwelling diffuse,
    !                               FLUXES(3,L) is downwelling direct,
    !                               L=1 is top, L=NUML+1 is bottom
    implicit none
    integer NLAYER
    logical DELTAM
    real TEMPS(*)
    real OPTDEPTHS(*), ALBEDOS(*), ASYMMETRIES(*)
    real GNDTEMP, GNDEMIS, SKYRAD, SOLARFLUX, SOLARMU
    real WAVENO(2), WAVELEN
    character(len=1) :: SRCTYPE, UNITS
    real FLUXES(3, *)

    integer N, L, I
    double precision DELTAU, G, OMEGA, F
    double precision LAMBDA, R, T, D, CP, CM, A, B, X1, X2
    double precision REFLECT, TRANS, SOURCEP, SOURCEM
    double precision RADP1P, RADP1M, RADP2P, RADP2M
    double precision PI, MU0, SKYFLUX, GNDFLUX, PLANCK1, PLANCK2, C, TAU
    double precision EXLP, EXLM, V, DS, B1, B2, SOLPP, SOLPM
    double precision, allocatable :: LOWER(:), UPPER(:), DIAG(:), RHS(:)
    real BBRAD
    parameter(PI=3.1415926535)

    N = 2*NLAYER + 2
    allocate (LOWER(N), UPPER(N), DIAG(N), RHS(N))

    ! Compute the reflection, transmission, and source
    !   coefficients for each layer for the diffuse Eddington
    !   two stream problem.
    if (SRCTYPE .eq. 'T') then
        call PLANCK_FUNCTION(TEMPS(1), UNITS, WAVENO, WAVELEN, BBRAD)
        PLANCK1 = PI*BBRAD
    end if
    MU0 = abs(SOLARMU)
    TAU = 0.0
    I = 2
    do L = 1, NLAYER
        DELTAU = OPTDEPTHS(L)
        if (DELTAU .lt. 0.0) stop 'Error [EDDRTF]: TAU<0.'
        ! Special case for zero optical depth
        if (DELTAU .eq. 0.0) then
            TRANS = 1.0
            REFLECT = 0.0
            SOURCEP = 0.0
            SOURCEM = 0.0
        else
            OMEGA = ALBEDOS(L)
            G = ASYMMETRIES(L)
            if (DELTAM) then
                F = G**2
                DELTAU = (1 - OMEGA*F)*DELTAU
                OMEGA = (1 - F)*OMEGA/(1 - OMEGA*F)
                G = (G - F)/(1 - F)
            end if
            R = (1.0 - OMEGA*(4.0 - 3.0*G))/4.0
            T = (7.0 - OMEGA*(4.0 + 3.0*G))/4.0
            LAMBDA = sqrt(3.0*(1.0 - OMEGA)*(1.0 - OMEGA*G))
            ! Special case for conservative scattering (lambda=0)
            if (LAMBDA .eq. 0.0) then
                D = 1.0/(1.0 + T*DELTAU)
                TRANS = D
                REFLECT = -R*DELTAU*D
            else
                X1 = -R
                X2 = LAMBDA + T
                EXLP = DEXP(min(LAMBDA*DELTAU, 75.d0))
                EXLM = 1.0/EXLP
                TRANS = 2.*LAMBDA/(X2*EXLP + (LAMBDA - T)*EXLM)
                REFLECT = X1*(EXLP - EXLM)*TRANS/(2.*LAMBDA)
                D = 1.0/(X2**2*EXLP - X1**2*EXLM)
            end if

            if (SRCTYPE .eq. 'T') then
                ! Calculate thermal source terms
                call PLANCK_FUNCTION(TEMPS(L + 1), UNITS, WAVENO, WAVELEN, BBRAD)
                PLANCK2 = PI*BBRAD
                V = 2.0*(PLANCK2 - PLANCK1)/(3.0*(1.-OMEGA*G)*DELTAU)
                RADP1P = -V + PLANCK1
                RADP2M = V + PLANCK2
                RADP2P = -V + PLANCK2
                RADP1M = V + PLANCK1
                if (LAMBDA .eq. 0.0) then
                    A = (R*DELTAU*RADP1P - RADP2M)*D
                    B = -(R*RADP1P + T*RADP2M)*D
                    SOURCEP = (B - T*(A + B*DELTAU))/R + RADP2P
                    SOURCEM = A + RADP1M
                else
                    CP = (X1*EXLM*RADP1P - X2*RADP2M)*D
                    CM = (-X2*EXLP*RADP1P + X1*RADP2M)*D
                    SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
                    SOURCEM = X2*CP + X1*CM + RADP1M
                end if
                PLANCK1 = PLANCK2
                FLUXES(3, L) = 0.0
            else
                ! Calculate solar source terms
                FLUXES(3, L) = SOLARFLUX*exp(-TAU/MU0)
                DS = 1.0/(LAMBDA**2 - 1.0/MU0**2)
                B1 = 0.5*OMEGA*(SOLARFLUX/MU0)*exp(-TAU/MU0)*DS
                B2 = 0.5*OMEGA*(SOLARFLUX/MU0)*exp(-(TAU + DELTAU)/MU0)*DS
                SOLPP = 1.0 + 1.5*G*MU0
                SOLPM = -1.0 + 1.5*G*MU0
                RADP1P = ((T + 1.0/MU0)*SOLPP + R*SOLPM)*B1
                RADP2M = ((-T + 1.0/MU0)*SOLPM - R*SOLPP)*B2
                RADP2P = ((T + 1.0/MU0)*SOLPP + R*SOLPM)*B2
                RADP1M = ((-T + 1.0/MU0)*SOLPM - R*SOLPP)*B1
                if (LAMBDA .eq. 0.0) then
                    A = (R*DELTAU*RADP1P - RADP2M)*D
                    B = -(R*RADP1P + T*RADP2M)*D
                    SOURCEP = (B - T*(A + B*DELTAU))/R + RADP2P
                    SOURCEM = A + RADP1M
                else
                    CP = (X1*EXLM*RADP1P - X2*RADP2M)*D
                    CM = (-X2*EXLP*RADP1P + X1*RADP2M)*D
                    SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
                    SOURCEM = X2*CP + X1*CM + RADP1M
                end if
                TAU = TAU + DELTAU
            end if
        end if
        DIAG(I) = -REFLECT
        DIAG(I + 1) = -REFLECT
        LOWER(I) = 1.0
        LOWER(I + 1) = -TRANS
        UPPER(I) = -TRANS
        UPPER(I + 1) = 1.0
        RHS(I) = SOURCEM
        RHS(I + 1) = SOURCEP
        I = I + 2
    end do

    ! Set up boundary radiances
    if (SRCTYPE .eq. 'S') then
        FLUXES(3, NLAYER + 1) = SOLARFLUX*exp(-TAU/MU0)
        GNDFLUX = (1.0 - GNDEMIS)*SOLARFLUX*exp(-TAU/MU0)
        SKYFLUX = PI*SKYRAD
    else
        FLUXES(3, NLAYER + 1) = 0.0
        call PLANCK_FUNCTION(GNDTEMP, UNITS, WAVENO, WAVELEN, BBRAD)
        GNDFLUX = PI*BBRAD*GNDEMIS
        call PLANCK_FUNCTION(SKYRAD, UNITS, WAVENO, WAVELEN, BBRAD)
        SKYFLUX = PI*BBRAD
    end if
    ! Setup for and call the tri-diagonal matrix solver
    RHS(1) = SKYFLUX
    DIAG(1) = 0.0
    UPPER(1) = 1.0
    DIAG(N) = -(1.0 - GNDEMIS)
    LOWER(N) = 1.0
    RHS(N) = GNDFLUX
    call TRIDIAG(N, LOWER, DIAG, UPPER, RHS)

    ! Put the fluxes in the output array
    if (UNITS .eq. 'T') then
        C = 1.0/PI
    else
        C = 1.0
    end if
    I = 1
    do L = 1, NLAYER + 1
        FLUXES(1, L) = C*RHS(I)
        FLUXES(2, L) = C*RHS(I + 1)
        I = I + 2
    end do

    deallocate (LOWER, UPPER, DIAG, RHS)
    return
end subroutine

subroutine TRIDIAG(N, LOWER, DIAG, UPPER, RHS)

    ! Computes the solution to a tridiagonal system.
    ! N is order of the matrix.  LOWER(2..N) is the subdiagonal,
    ! DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the
    ! superdiagonal.  On input RHS is the right hand side, while
    ! on output it is the solution vector.  Everything is destroyed.
    ! Hacked from Linpack DGTSL.

    implicit none
    integer N
    double precision LOWER(*), DIAG(*), UPPER(*), RHS(*)
    integer K, KB
    double precision T

    if (N .eq. 1) then
        if (DIAG(1) .eq. 0.0) goto 990
        RHS(1) = RHS(1)/DIAG(1)
    end if

    LOWER(1) = DIAG(1)
    DIAG(1) = UPPER(1)
    UPPER(1) = 0.0
    UPPER(N) = 0.0
    do K = 1, N - 1
        ! Interchange this and next row to the get the largest pivo
        if (abs(LOWER(K + 1)) .ge. abs(LOWER(K))) then
            T = LOWER(K + 1)
            LOWER(K + 1) = LOWER(K)
            LOWER(K) = T
            T = DIAG(K + 1)
            DIAG(K + 1) = DIAG(K)
            DIAG(K) = T
            T = UPPER(K + 1)
            UPPER(K + 1) = UPPER(K)
            UPPER(K) = T
            T = RHS(K + 1)
            RHS(K + 1) = RHS(K)
            RHS(K) = T
        end if
        if (LOWER(K) .eq. 0.0) goto 990
        T = -LOWER(K + 1)/LOWER(K)
        LOWER(K + 1) = DIAG(K + 1) + T*DIAG(K)
        DIAG(K + 1) = UPPER(K + 1) + T*UPPER(K)
        UPPER(K + 1) = 0.0
        RHS(K + 1) = RHS(K + 1) + T*RHS(K)
    end do
    if (LOWER(N) .eq. 0.0) goto 990

    ! Back substitute
    RHS(N) = RHS(N)/LOWER(N)
    RHS(N - 1) = (RHS(N - 1) - DIAG(N - 1)*RHS(N))/LOWER(N - 1)
    do KB = 1, N - 2
        K = N - 2 - KB + 1
        RHS(K) = (RHS(K) - DIAG(K)*RHS(K + 1) - UPPER(K)*RHS(K + 2))/LOWER(K)
    end do
    return

990 continue
    stop 'Error [TRIDIAG]: Singular matrix in TRIDIAG.'
end subroutine

subroutine MAKE_ANGLE_SET(NMU, NPHI, ITYPE, NPHI0MAX,&
    &NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the set of angles for the discrete space representation.
    ! The number of mu's (cosine zenith angles) and maximum number of
    ! phi's (azimuth angles) is input.  The actual number of azimuthal
    ! angles for each mu is output (NPHI0).  There are three types of
    ! discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
    ! gaussian grid, and 3 is a reduced double gaussian set.
    ! The output is the NMU mu values, the NPHI0 phi values for each mu,
    ! and the integration weight for each ordinate. The first NMU/2 mu
    ! angles are the downwelling (mu<0) angles.  Also output are the
    ! maximum number of azimuthal angles (NPHI0MAX), the flags for doing
    ! an azimuthal FFT for each mu and the total number of angles.

    implicit none
    integer NMU, NPHI, ITYPE, NPHI0MAX, NPHI0(NMU), NANG
    logical FFTFLAG(NMU)
    real MU(NMU), PHI(NMU, *), WTMU(NMU), WTDO(NMU, *)
    integer MAXNPHI, J, K, MM
    real DELPHI
    parameter(MAXNPHI=256)
    integer GOODNFFT(MAXNPHI)
    data GOODNFFT/1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 16,&
        &18, 18, 20, 20, 24, 24, 24, 24, 25, 27, 27, 30, 30, 30, 32, 32,&
        &36, 36, 36, 36, 40, 40, 40, 40, 45, 45, 45, 45, 45, 48, 48, 48,&
        &50, 50, 54, 54, 54, 54, 60, 60, 60, 60, 60, 60, 64, 64, 64, 64,&
        &72, 72, 72, 72, 72, 72, 72, 72, 80, 80, 80, 80, 80, 80, 80, 80,&
        &81, 90, 90, 90, 90, 90, 90, 90, 90, 90, 96, 96, 96, 96, 96, 96,&
        &100, 100, 100, 100, 108, 108, 108, 108,&
        &108, 108, 108, 108, 120, 120, 120, 120,&
        &120, 120, 120, 120, 120, 120, 120, 120,&
        &128, 128, 128, 128, 128, 128, 128, 128,&
        &144, 144, 144, 144, 144, 144, 144, 144,&
        &144, 144, 144, 144, 144, 144, 144, 144,&
        &160, 160, 160, 160, 160, 160, 160, 160,&
        &160, 160, 160, 160, 160, 160, 160, 160,&
        &180, 180, 180, 180, 180, 180, 180, 180, 180, 180,&
        &180, 180, 180, 180, 180, 180, 180, 180, 180, 180,&
        &192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192, 192,&
        &200, 200, 200, 200, 200, 200, 200, 200,&
        &216, 216, 216, 216, 216, 216, 216, 216,&
        &216, 216, 216, 216, 216, 216, 216, 216,&
        &240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240,&
        &240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 240,&
        &256, 256, 256, 256, 256, 256, 256, 256,&
        &256, 256, 256, 256, 256, 256, 256, 256/

    MM = max(0, int(NPHI/2) - 1)
    NANG = 0
    if (ITYPE .ge. 1 .or. ITYPE .le. 3) then
        if (ITYPE .le. 2) then
            call GAUSQUADS(NMU, MU, WTMU)
        else if (ITYPE .eq. 3) then
            call DGAUSQUADS(NMU, MU, WTMU)
        end if
        do J = 1, NMU
            NPHI0MAX = NPHI
            if (ITYPE .eq. 1) then
                NPHI0(J) = NPHI0MAX
            else
                if (NPHI0MAX .gt. MAXNPHI)&
                    &stop 'Error [MAKE_ANGLE_SET]: Exceeded GOODNFFT.'
                NPHI0(J) = int(0.9 + NPHI0MAX*sqrt(1 - MU(J)**2))
                NPHI0(J) = min(NPHI0MAX, GOODNFFT(NPHI0(J)))
            end if
            DELPHI = 2.0*acos(-1.0)/NPHI0(J)
            do K = 1, NPHI0(J)
                PHI(J, K) = (K - 1)*DELPHI
                WTDO(J, K) = DELPHI*WTMU(J)
            end do
            NANG = NANG + NPHI0(J)
            FFTFLAG(J) = (NPHI0(J) .gt. 14) .or. (MM .gt. 15)
        end do

    else
        stop 'Error [MAKE_ANGLE_SET]: Invalid discrete ordinate type.'
    end if
    return
end subroutine

subroutine MAKE_SH_DO_COEF(NSTLEG, ML, MM, NLM, NMU, NPHI0,&
    &NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,&
    &CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

    ! Makes the transformation coefficients for the real generalized
    ! spherical harmonic transform.  The main coefficients are output in
    ! four arrays: 1 is for the SH_TO_DO forward transform, 2 is for the
    ! DO_TO_SH back transform, which contains the discrete angle integration
    ! weights.  The coefficients are divided up into the mu dependent set
    ! CMUn (function of l, m, mu_j), and the phi dependent set CPHIn
    ! (function of m, phi_k) for each mu_j.
    ! The FFTPACK phase coefficients for the FFT in azimuth are also
    ! output in WPHISAVE.

    implicit none
    integer NSTLEG, ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX
    logical FFTFLAG(NMU)
    real MU(NMU), PHI(NMU, *), WTMU(NMU), WTDO(NMU, *)
    real CMU1(NSTLEG, NLM, NMU), CMU2(NSTLEG, NMU, NLM)
    real CPHI1(-16:16, 32, NMU), CPHI2(32, -16:16, NMU)
    real WPHISAVE(3*NPHI0MAX + 15, NMU)
    integer I, J, K, M, Q
    real X, W
    real, allocatable :: PRC(:, :)

    ! Make the mu part of the real generalized spherical harmonics
    allocate (PRC(6, NLM))
    do I = 1, NMU
        X = MU(I)
        PRC(:, :) = 0.0
        call PLMALL(.false., X, ML, MM, PRC)
        do J = 1, NLM
            do Q = 1, NSTLEG
                CMU1(Q, J, I) = PRC(Q, J)
            end do
        end do
        PRC(:, :) = 0.0
        call PLMALL(.true., X, ML, MM, PRC)
        do J = 1, NLM
            do Q = 1, NSTLEG
                CMU2(Q, I, J) = PRC(Q, J)*WTMU(I)
            end do
        end do
    end do
    deallocate (PRC)

    ! Make the to and from Fourier coefficients for each mu
    do I = 1, NMU
        if (.not. FFTFLAG(I)) then
            ! If not doing an FFT for this mu then make the DFT coeffici
            W = 1.0/WTMU(I)
            do K = 1, NPHI0(I)
                if (NPHI0(I) .gt. 32 .or. MM .gt. 16) stop &
                    &'Error [MAKE_SH_DO_COEF]: Fourier coefficients array exceeded.'
                CPHI1(0, K, I) = 1.0
                CPHI2(K, 0, I) = WTDO(I, K)*W
                do M = 1, MM
                    CPHI1(M, K, I) = cos(M*PHI(I, K))
                    CPHI1(-M, K, I) = sin(M*PHI(I, K))
                    CPHI2(K, M, I) = CPHI1(M, K, I)*WTDO(I, K)*W
                    CPHI2(K, -M, I) = CPHI1(-M, K, I)*WTDO(I, K)*W
                end do
            end do
        else
            ! Precompute the phase factors for the FFTs
            call DFFTI(NPHI0(I), WPHISAVE(1, I))
        end if
    end do
    return
end subroutine

subroutine SUM_OUTPUT(IFLAG, WT, N, data, SUMDATA)
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
end subroutine

subroutine COMPUTE_NETFLUXDIV(NSTOKES, NPTS, RSHPTR, SRCTYPE,&
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
end subroutine

subroutine COMPUTE_SH(NSHOUT, NSTOKES, NPTS, SRCTYPE,&
    &SOLARMU, SOLARAZ, DIRFLUX,&
    &RSHPTR, ML, MM, RADIANCE, SHTERMS)

    !   Computes the quantities for the spherical harmonic output.
    ! At each grid point the mean radiance and net flux (Fx, Fy, Fz)
    ! are computed.  The output includes the diffuse and direct solar
    ! components.  If NSHOUT=5 then the HIGHORDERRAD flag is set and the
    ! rms of the higher order terms in the radiance series is computed.

    implicit none
    integer NSTOKES, NPTS, RSHPTR(NPTS + 1), ML, MM, NSHOUT
    real SOLARMU, SOLARAZ, DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *), SHTERMS(NSHOUT, NPTS)
    character(len=1) :: SRCTYPE
    integer I, IR, J, NLM
    real PI, A, C, C0, IMEAN, FX, FY, FZ, HORMS
    real SECSOL, SINSOL, SOLX, SOLY, SOLM, F0
    real, allocatable :: YLMSUN(:)

    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    allocate (YLMSUN(NLM))

    ! Spherical Harmonic output: individual SH terms
    PI = acos(-1.0)
    A = sqrt(2.0*PI/3.0)
    C = sqrt(4.0*PI/3.0)
    C0 = 1.0/sqrt(4.0*PI)
    SECSOL = 1.0/abs(SOLARMU)
    SINSOL = sqrt(1.0 - SOLARMU**2)
    SOLX = SINSOL*cos(SOLARAZ)*SECSOL
    SOLY = SINSOL*sin(SOLARAZ)*SECSOL
    SOLM = SECSOL/(4.0*PI)
    F0 = 0.0
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        call YLMALL(.true., SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
    end if
    do I = 1, NPTS
        IR = RSHPTR(I)
        IMEAN = C0*RADIANCE(1, IR + 1)
        FZ = C*RADIANCE(1, IR + 3)
        if (MM .gt. 0) then
            ! FX = -C*RADIANCE(1,IR+4)   ! old spherical harmonics basis
            ! real generaliz
            FX = -A*(RADIANCE(1, IR + 4) + RADIANCE(1, IR + 2))
        else
            FX = 0.0
        end if
        if (MM .gt. 0) then
            ! FY = -C*RADIANCE(1,IR+2)
            ! real generaliz
            FY = -A*(RADIANCE(1, IR + 2) - RADIANCE(1, IR + 4))
        else
            FY = 0.0
        end if
        if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
            IMEAN = IMEAN + DIRFLUX(I)*SOLM
            FX = FX + DIRFLUX(I)*SOLX
            FY = FY + DIRFLUX(I)*SOLY
            FZ = FZ - DIRFLUX(I)
        end if
        SHTERMS(1, I) = IMEAN
        SHTERMS(2, I) = FX
        SHTERMS(3, I) = FY
        SHTERMS(4, I) = FZ
        if (NSHOUT .eq. 5) then
            HORMS = 0.0
            if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
                F0 = SECSOL*DIRFLUX(I)
            end if
            do J = 5, RSHPTR(I + 1) - IR
                HORMS = HORMS + (F0*YLMSUN(J) + RADIANCE(1, IR + J))**2
            end do
            SHTERMS(5, I) = sqrt(HORMS)
        end if
    end do
    deallocate (YLMSUN)
    return
end subroutine

subroutine INTERP_OUTPUT(OLDNPTS, NPTS,&
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
                stop 'Error [INTERP_OUTPUT]: Point not on an edge.'
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
end subroutine

subroutine INTERP_RADIANCE(NSTOKES, OLDNPTS, NPTS, RSHPTR,&
    &RADIANCE, NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)

    ! Interpolates the spherical harmonic radiance array to the new
    !   grid points (>OLDNPTS) from the old grid points.

    implicit none
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
                stop 'Error [INTERP_OUTPUT]: Point not on an edge.'
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
                    RADIANCE(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
                end do
            end if
        end do
    end do

    return
end subroutine

subroutine VISUALIZE_RADIANCE(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
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
        call YLMALL(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG, YLMSUN)
        if (DELTAM .and. NUMPHASE .gt. 0) then
            NSCATANGLE = max(36, min(MAXSCATANG, 2*NLEG))
            NSTPHASE = min(NSTLEG, 2)
            allocate (PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE))
            call PRECOMPUTE_PHASE(NSCATANGLE, NUMPHASE, NSTPHASE,&
                &NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
        end if
    end if

    ! Make the isotropic radiances for the top boundary
    call COMPUTE_TOP_RADIANCES(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces
    !   Compute the upwelling bottom radiances using the downwelling
    if (SFCTYPE .eq. 'FL') then
        call FIXED_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
    else if (SFCTYPE .eq. 'VL') then
        call VARIABLE_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1, 2),&
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
            call INTEGRATE_1RAY(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
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
end subroutine

subroutine INTEGRATE_1RAY(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
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
    call YLMALL(.false., SNGL(MU2), SNGL(PHI2), ML, MM, NSTLEG, YLMDIR)

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
                    call ROTATE_POL_PLANE(NSTOKES, COSSCAT, SOLARMU,&
                        &SNGL(MU2), SOLARAZ - SNGL(PHI2), SINGSCAT(:, I))
                end if
            end do
        else
            allocate (SUNDIRLEG(0:NLEG))
            call LEGENDRE_ALL(COSSCAT, NLEG, SUNDIRLEG)
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
    call LOCATE_GRID_CELL(NX, NY, NZ, XGRID, YGRID, ZGRID,&
        &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
        &BCFLAG, IPFLAG, XE, YE, ZE, ICELL)
    ! print '(A,5(1X,F8.5),1X,I6)',&
    !     &'INTEGRATE_1RAY:',X0,Y0,Z0,MU2,PHI2,ICELL
    IFACE = 0
    NGRID = 0

    ! Loop until reach a Z boundary or transmission is very small
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
        ! Make sure current cell is valid
        if (ICELL .le. 0) then
            write (6, *) 'Error [INTEGRATE_1RAY]: ICELL=', ICELL,&
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
            call COMPUTE_SOURCE_1CELL_UNPOL(ICELL, GRIDPTR,&
                &ML, MM, NLM, NLEG, NUMPHASE,&
                &NPTS, DELTAM, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX,&
                &SHPTR, SOURCE, YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT,&
                &DONETHIS, OLDIPTS, OEXTINCT8, OSRCEXT8,&
                &EXTINCT8, SRCEXT8)
        else
            call COMPUTE_SOURCE_1CELL(ICELL, GRIDPTR,&
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
            write (6, *) 'Error [INTEGRATE_1RAY]: SO<0 ',&
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
            call NEXT_CELL(XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,&
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
            call FIND_BOUNDARY_RADIANCE(NSTOKES, XN, YN,&
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
end subroutine

subroutine FIND_BOUNDARY_RADIANCE(NSTOKES, XB, YB, MU2, PHI2,&
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
                &stop 'Error [FIND_BOUNDARY_RADIANCE]: Not at boundary.'
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
                &stop 'Error [FIND_BOUNDARY_RADIANCE]: Not at boundary.'

            if (.not. LAMBERTIAN) then
                call VARIABLE_BRDF_SURFACE(NBOTPTS, IBC, IBC, BCPTR(1, 2),&
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
end subroutine

subroutine COMPUTE_SOURCE_1CELL(ICELL, GRIDPTR,&
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
end subroutine

subroutine COMPUTE_SOURCE_1CELL_UNPOL(ICELL, GRIDPTR,&
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
end subroutine

subroutine PRECOMPUTE_PHASE(NSCATANGLE, NUMPHASE, NSTPHASE,&
    &NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
    !   Precomputes the phase matrix elements I-I and I-Q as a function
    ! of scattering angle for solar direct scattering for all the
    ! tabulated phase functions. Output is in PHASETAB.
    ! NSTPHASE=1 for NSTOKES=1 and NSTPHASE=2 for NSTOKES>1.
    implicit none
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
        call WIGNERFCT(X, NLEG, 0, 0, DMM1)
        if (NSTOKES .gt. 1) then
            call WIGNERFCT(X, NLEG, 2, 0, DMM2)
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
                write (6, *) 'Error [PRECOMPUTE_PHASE]: negative phase function',&
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
end subroutine

subroutine ROTATE_POL_PLANE(NSTOKES, COSSCAT,&
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
end subroutine

subroutine COMPUTE_RADIANCE(&
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
    call COMPUTE_TOP_RADIANCES(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))

    ! Make the bottom boundary radiances for the Lambertian surfaces.
    ! Compute the upwelling bottom radiances using the downwelling flu
    if (SFCTYPE .eq. 'FL') then

        call FIXED_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))

    else if (SFCTYPE .eq. 'VL') then

        call VARIABLE_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1, 2),&
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
            write (6, *) 'Warning [COMPUTE_RADIANCE]: Bad mu for radiance', MUOUT
        else

            ! Compute the source function throughout grid for this angle
            call COMPUTE_ONE_SOURCE(NSTOKES, ML, MM, NLM, NSTLEG, NLEG,&
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
                    call VARIABLE_BRDF_SURFACE(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
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
                    call INTEGRATE_SOURCE(BCFLAG, IPFLAG, NSTOKES,&
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

end subroutine

subroutine COMPUTE_ONE_SOURCE(NSTOKES,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NPTS, DELTAM, MU, PHI, SRCTYPE, SOLARMU, SOLARAZ,&
    &ALBEDO, LEGEN, IPHASE, DIRFLUX,&
    &SHPTR, SOURCE, SOURCE1)

    ! Computes the source function (SOURCE1) in the direction (MU,PHI)
    ! for the whole domain (NX,NZ).  The spherical harmonic source function
    ! series is input in SOURCE.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    implicit none
    integer NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer SHPTR(NPTS + 1)
    integer*2 IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ, MU, PHI
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS), DIRFLUX(NPTS)
    real SOURCE(NSTOKES, *), SOURCE1(NSTOKES, NPTS)
    character(len=1) :: SRCTYPE
    integer I, J, JT, K, L, M, N, IPH, IS, NS, ME, MS
    real OFOURPI, DA, F, A, A1, B1, SECMU0, COSSCAT
    real, allocatable :: YLMDIR(:, :), YLMSUN(:, :)
    double precision, allocatable :: SUNDIRLEG(:)
    real, allocatable :: UNSCLEGEN(:, :), SINGSCAT(:, :)
    logical TMS

    allocate (YLMDIR(NSTLEG, NLM), YLMSUN(NSTLEG, NLM))
    allocate (SUNDIRLEG(0:NLEG))
    allocate (UNSCLEGEN(2, 0:NLEG), SINGSCAT(NSTOKES, NUMPHASE))
    if (abs(MU) .gt. 1.0) stop 'Error [COMPUTE_ONE_SOURCE]: Bad mu.'

    TMS = DELTAM
    ! TMS = .FALSE.

    ! Precompute Ylm's for output direction and solar direction
    call YLMALL(.false., MU, PHI, ML, MM, NSTLEG, YLMDIR)
    if (SRCTYPE .ne. 'T') then
        OFOURPI = 1.0/(4.0*acos(-1.0))
        call YLMALL(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG,&
            &YLMSUN)
        SECMU0 = 1.0/abs(SOLARMU)
        COSSCAT = SOLARMU*MU + sqrt((1.0d0 - SOLARMU**2)*(1.0d0 - MU**2))&
            &*cos(SOLARAZ - PHI)
        if (NUMPHASE .eq. 0 .and. TMS) then
            ! If non-tabulate phase functions (unpolarized) then
            !   compute the Legendre polynomials for the scattering angle
            !   for the untruncated solar single scattering computation.
            call LEGENDRE_ALL(dble(COSSCAT), NLEG, SUNDIRLEG)
            do L = 0, NLEG
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)*OFOURPI
            end do
        end if
    end if

    ! If solar source and tabulated phase functions then precompute
    !   the single scattered solar Stokes vector to outgoing direction
    if (SRCTYPE .ne. 'T' .and. TMS .and. NUMPHASE .gt. 0) then
        do IPH = 1, NUMPHASE
            ! Unscaled the needed Legendre coefficients, but divide by 1-f
            !   because extinction is still scaled (for TMS method)
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
            ! Sum the Wigner series (for A1 and B1), rotate the plane of
            !   polarization to the outgoing direction, and return the Stoke
            call SOLAR_SINGSCAT_VECT(NLEG, UNSCLEGEN, COSSCAT,&
                &SOLARMU, MU, SOLARAZ - PHI, NSTOKES,&
                &SINGSCAT(1, IPH))
            if (SINGSCAT(1, IPH) .le. 0.0) then
                write (6, *) 'Error [COMPUTE_ONE_SOURCE]: negative source function',&
                    &' for tabulated phase function:', IPH, MU, PHI, SINGSCAT(1, IPH)
                stop
            end if
        end do
    end if

    ! Loop over all the grid points, computing the source function
    !   at the viewing angle from the spherical harmonic source func
    do I = 1, NPTS
        IS = SHPTR(I)
        NS = SHPTR(I + 1) - IS
        ! Sum over the real generalized spherical harmonic series
        !   of the source function
        SOURCE1(:, I) = 0.0
        do J = 1, NS
            SOURCE1(1, I) = SOURCE1(1, I) + SOURCE(1, IS + J)*YLMDIR(1, J)
        end do
        if (NSTOKES .gt. 1) then
            do J = 1, NS
                SOURCE1(2, I) = SOURCE1(2, I) + SOURCE(2, IS + J)*YLMDIR(2, J)&
                    &+ SOURCE(3, IS + J)*YLMDIR(5, J)
                SOURCE1(3, I) = SOURCE1(3, I) + SOURCE(2, IS + J)*YLMDIR(6, J)&
                    &+ SOURCE(3, IS + J)*YLMDIR(3, J)
            end do
        end if
        if (NSTOKES .eq. 4) then
            do J = 1, NS
                SOURCE1(4, I) = SOURCE1(4, I) + SOURCE(4, IS + J)*YLMDIR(4, J)
            end do
        end if

        ! Special case for solar source and Delta-M
        if (SRCTYPE .ne. 'T' .and. TMS) then
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            ! First subtract off the truncated single scattering
            DA = ALBEDO(I)*DIRFLUX(I)*SECMU0
            J = 1
            do L = 0, ML
                ME = min(L, MM)
                MS = -ME
                A1 = DA*LEGEN(1, L, IPH)
                B1 = DA*LEGEN(5, L, IPH)
                if (J .le. NS) then
                    JT = J
                    do M = MS, ME
                        SOURCE1(1, I) = SOURCE1(1, I) - A1*YLMDIR(1, J)*YLMSUN(1, J)
                        J = J + 1
                    end do
                    if (NSTOKES .gt. 1) then
                        J = JT
                        do M = MS, ME
                            SOURCE1(2, I) = SOURCE1(2, I) - B1*YLMDIR(2, J)*YLMSUN(1, J)
                            SOURCE1(3, I) = SOURCE1(3, I) - B1*YLMDIR(6, J)*YLMSUN(1, J)
                            J = J + 1
                        end do
                    end if
                end if
            end do

            ! Then add in the single scattering contribution for the
            !   original unscaled phase matrix.
            if (NUMPHASE .gt. 0) then
                SOURCE1(:, I) = SOURCE1(:, I) + DA*SINGSCAT(:, IPH)
            else if (NSTOKES .eq. 1) then
                ! If we have a phase function for each grid point (which
                !   should be the unpolarized case) then unscale the Legendre
                !   coefficients for L<=ML and sum over the Legendre series
                !   in cosine of the scattering angle.
                F = LEGEN(1, ML + 1, IPH)
                do L = 0, NLEG
                    if (L .le. ML) then
                        A = DA*(LEGEN(1, L, IPH) + F/(1 - F))
                    else
                        A = DA*LEGEN(1, L, IPH)/(1 - F)
                    end if
                    SOURCE1(1, I) = SOURCE1(1, I) + A*SUNDIRLEG(L)
                end do
            end if
        end if
    end do

    deallocate (YLMDIR, YLMSUN, SUNDIRLEG, UNSCLEGEN, SINGSCAT)
    return
end subroutine

subroutine SOLAR_SINGSCAT_VECT(NLEG, SINGLEGEN, COSSCAT,&
    &SOLARMU, MU, DELPHI, NSTOKES,&
    &SCATVECT)

    ! Calculates the Stokes vector for a single scattering of unpolarized
    ! sunlight from the solar direction (SOLARMU) to the outgoing
    ! direction (MU) and with cosine scattering angle (COSSCAT).
    ! The A1 and B1 elements of the phase matrix are calculated from
    ! the Wigner function expansion coefficients (SINGLEGEN).
    ! The plane of polarization is rotated from the scattering plane to
    ! the outgoing meridional plane.  Also divides the output vector by 4*pi.

    implicit none
    integer NLEG, NSTOKES
    real SINGLEGEN(2, 0:NLEG), COSSCAT, SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    integer I, J, K, L
    double precision OFOURPI, X, FCT, A1, B1
    double precision, allocatable :: DMM1(:)
    double precision SIN_SCAT, SIN_THETA1, SIN_THETA2
    double precision SINPHI, COSPHI, SIN2, COS2, SIN22, COS22
    parameter(OFOURPI=0.07957747151)

    allocate (DMM1(0:NLEG))

    ! Sum the first Wigner function series, which is actually a Legend
    X = dble(COSSCAT)
    call WIGNERFCT(X, NLEG, 0, 0, DMM1)
    A1 = 0.0d0
    do L = 0, NLEG
        FCT = 2.0d0*L + 1.0d0
        A1 = A1 + FCT*dble(SINGLEGEN(1, L))*DMM1(L)
    end do
    SCATVECT(1) = A1*OFOURPI

    if (NSTOKES .gt. 1) then
        ! If doing polarization, sum the second Wigner function series
        call WIGNERFCT(X, NLEG, 2, 0, DMM1)
        B1 = 0.0d0
        do L = 0, NLEG
            FCT = 2.0d0*L + 1.0d0
            B1 = B1 - FCT*dble(SINGLEGEN(2, L))*DMM1(L)
        end do

        ! Rotate the plane of polarization from the scattering plane
        !   to the outgoing meridional plane.
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
        SCATVECT(2) = B1*COS22*OFOURPI
        SCATVECT(3) = B1*SIN22*OFOURPI
    end if
    if (NSTOKES .eq. 4) then
        SCATVECT(4) = 0.0
    end if
    deallocate (DMM1)
    return
end subroutine

subroutine INTEGRATE_SOURCE(BCFLAG, IPFLAG, NSTOKES,&
    &NX, NY, NZ, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &MU, PHI, GRIDRAD, EXTINCT, SOURCE,&
    &X0, Y0, Z0, XE, YE, ZE, SIDE,&
    &TRANSMIT, RADIANCE, VALIDRAD)

    ! Integrates the source function (SOURCE) through the extinction
    ! field (EXTINCT) backward in the direction opposite to (MU,PHI)
    ! to find the outgoing radiance (RADIANCE) at the point X0,Y0,Z0.
    ! The radiances at the opposite boundary are input in GRIDRAD.
    ! The transmission and radiance of the ray so far (TRANSMIT, RADIANCE)
    ! are input and returned after the integration along with the exitting
    ! ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
    ! 5=-Z,6=+Z).

    implicit none
    integer BCFLAG, IPFLAG, NSTOKES, NX, NY, NZ, NPTS, NCELLS, SIDE
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    logical VALIDRAD
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real MU, PHI
    double precision X0, Y0, Z0, XE, YE, ZE, TRANSMIT, RADIANCE(NSTOKES)
    real EXTINCT(NPTS)
    real SOURCE(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
    integer I1, I2, I3, I4, IOPP, NTAU, IT, K
    logical IPINX, IPINY, OPENBCFACE, BTEST
    integer GRIDFACE(4, 6), OPPFACE(6), JFACE, KFACE, IC
    real EXT0(1), EXT1(1), EXTN(1)
    real XM, YM
    real SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RAD0(NSTOKES)
    double precision PI, CX, CY, CZ, CXINV, CYINV, CZINV
    double precision XN, YN, ZN, XI, YI, ZI
    double precision SO, SOX, SOY, SOZ, EPS
    double precision TAUTOL, TAUGRID, S, DELS
    double precision EXT, TAU, TRANSCELL, ABSCELL
    double precision SRC(NSTOKES), RAD(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/

    ! added for fixing the "domain shift" issue
    !╭──────────────────────────────────────────────────────────────╮!
    !real :: DS_Z0, DS_VZA_RAD, DS_XOFFSET, DS_YOFFSET

    !if (IPFLAG .eq. 0) then
    !  DS_Z0 = ZGRID(NZ) - ZGRID(1)
    !  DS_VZA_RAD = acos(MU)
    !  DS_XOFFSET = DS_Z0*tan(DS_VZA_RAD)*cos(PHI)
    !  DS_YOFFSET = DS_Z0*tan(DS_VZA_RAD)*sin(PHI)
    !else
    !  DS_XOFFSET = 0.0
    !  DS_YOFFSET = 0.0
    !endif
    !X0 = X0 + DS_XOFFSET
    !Y0 = Y0 + DS_YOFFSET
    !╰──────────────────────────────────────────────────────────────╯!

    ! TAUTOL is the maximum optical path for the subgrid intervals
    TAUTOL = 0.1
    EPS = 1.0e-4*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    PI = acos(-1.0)
    CX = sqrt(1.0 - MU**2)*cos(PHI + PI)
    CY = sqrt(1.0 - MU**2)*sin(PHI + PI)
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
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [INTEGRATE_SOURCE]: Bad MU.'
    end if
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    XM = 0.5*(XGRID(1) + XGRID(NX))
    YM = 0.5*(YGRID(1) + YGRID(NY))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    call LOCATE_GRID_CELL(NX, NY, NZ, XGRID, YGRID, ZGRID,&
        &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
        &BCFLAG, IPFLAG, XE, YE, ZE, ICELL)
    ! print '(A,5(1X,F8.5),1X,I6)', &
    !     &'INTEGRATE_SOURCE:', X0, Y0, Z0, MU, PHI, ICELL
    call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
        &XE, YE, ZE, 1, EXTINCT, 1, 1, EXT1)
    call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
        &XE, YE, ZE, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT1)
    SRCEXT1(1) = max(0.0, SRCEXT1(1))

    ! Loop until finding a face with known radiances
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
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
        ! from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9 - IOCT, ICELL)
        ! Get the distances to the 3 planes and select the closest
        !   (always need to deal with the cell that is wrapped)
        if (IPINX .or. CX .eq. 0.0) then
            SOX = 1.0e20
        else
            SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
        end if
        if (IPINY .or. CY .eq. 0.0) then
            SOY = 1.0e20
        else
            SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
        end if
        SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Error [INTEGRATE_SOURCE]: SO<0 ',&
                &MU, PHI, X0, Y0, Z0, XE, YE, ZE, SO, SOX, SOY, SOZ, XGRID(1), YGRID(1),&
                &GRIDPOS(1, IOPP), GRIDPOS(2, IOPP), GRIDPOS(3, IOPP)
            stop
        end if
        if (IPINX) then
            XN = XE
        else
            XN = XE + SO*CX
        end if
        if (IPINY) then
            YN = YE
        else
            YN = YE + SO*CY
        end if
        ZN = ZE + SO*CZ

        ! Find the optical path across the grid cell and figure how
        !   many subgrid intervals to use
        call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
            &XN, YN, ZN, 1, EXTINCT, 1, 1, EXTN)
        TAUGRID = SO*0.5*(EXT1(1) + EXTN(1))
        NTAU = max(1, 1 + int(TAUGRID/TAUTOL))
        DELS = SO/NTAU
        ! Loop over the subgrid cells
        do IT = 1, NTAU
            S = IT*DELS
            if (IPINX) then
                XI = XE
            else
                XI = XE + S*CX
            end if
            if (IPINY) then
                YI = YE
            else
                YI = YE + S*CY
            end if
            ZI = ZE + S*CZ
            ! Interpolate extinction and source function along path
            if (IT .eq. NTAU) then
                EXT0(1) = EXTN(1)
            else
                call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
                    &XI, YI, ZI, 1, EXTINCT, 1, 1, EXT0)
            end if
            call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
                &XI, YI, ZI, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT0)
            SRCEXT0(1) = max(0.0, SRCEXT0(1))

            ! Compute the subgrid radiance: integration of the source fun
            EXT = 0.5*(EXT0(1) + EXT1(1))
            if (EXT .ne. 0.0) then
                TAU = EXT*DELS
                ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - 0.33333333333*TAU))
                TRANSCELL = 1.0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + 0.08333333333*&
                    &(EXT0(1)*SRCEXT1(:) - EXT1(1)*SRCEXT0(:))*DELS&
                    &*(1.0 - 0.05*(EXT1(1) - EXT0(1))*DELS))/EXT
            else
                ABSCELL = 0.0
                TRANSCELL = 1.0
                SRC(:) = 0.0
            end if

            RADIANCE(:) = RADIANCE(:) + TRANSMIT*SRC(:)*ABSCELL
            SRCEXT1(:) = SRCEXT0(:)
            EXT1(1) = EXT0(1)
            TRANSMIT = TRANSMIT*TRANSCELL
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
            call NEXT_CELL(XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
        end if
        ! Get the four grid points on the intersection face for
        !   telling if at the boundary
        !   If going to same or larger face then use previous face
        if (NEIGHPTR(IFACE, ICELL) .ge. 0 .and. .not. OPENBCFACE) then
            KFACE = IFACE
            IC = ICELL
        else
            ! If going to smaller face then use next face (more accurate
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
        end if
        I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
        I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
        I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
        I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
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

        ! See if there are valid radiances for the intersection face:
        !   If so, set flag, interpolate radiance, and add in
        !   If not, only prepare for next cell
        if (TRANSMIT .lt. 1.0e-5) then
            VALIDRAD = .true.
        else if (GRIDRAD(1, I1) .ge. 0.0 .and. GRIDRAD(1, I2) .ge. 0.0 .and.&
            &GRIDRAD(1, I3) .ge. 0.0 .and. GRIDRAD(1, I4) .ge. 0.0) then
            VALIDRAD = .true.
            call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS, XI, YI, ZI,&
                &NSTOKES, GRIDRAD, 1, NSTOKES, RAD0)
            RADIANCE(:) = RADIANCE(:) + TRANSMIT*RAD0(:)
        else
            ICELL = INEXTCELL
        end if
        XE = XN
        YE = YN
        ZE = ZN

    end do
    SIDE = IFACE
    return
end subroutine

subroutine INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
    &X, Y, Z, ND, FIELD, ID1, NV, VALUES)

    !   Trilinearly interpolates the field (FIELD) given a grid cell
    ! pointer (ICELL), to compute the field values (VALUES) at the desired
    ! point (X,Y,Z). ND is the leading index of the FIELD (1 for a
    ! scalar field), ID1 is the first element of the vector to interpolate,
    ! and NV elements are output.

    implicit none
    integer ICELL, GRIDPTR(8, *), ND, ID1, NV
    double precision X, Y, Z
    real GRIDPOS(3, *), FIELD(ND, *), VALUES(NV)
    integer IPT1, IPT2, ID
    double precision U, V, W, UM, VM, WM, DELX, DELY, DELZ

    IPT1 = GRIDPTR(1, ICELL)
    IPT2 = GRIDPTR(8, ICELL)
    DELX = GRIDPOS(1, IPT2) - GRIDPOS(1, IPT1)
    if (DELX .le. 0.0) DELX = 1.0
    DELY = GRIDPOS(2, IPT2) - GRIDPOS(2, IPT1)
    if (DELY .le. 0.0) DELY = 1.0
    DELZ = GRIDPOS(3, IPT2) - GRIDPOS(3, IPT1)
    U = (X - GRIDPOS(1, IPT1))/DELX
    V = (Y - GRIDPOS(2, IPT1))/DELY
    W = (Z - GRIDPOS(3, IPT1))/DELZ
    UM = 1 - U
    VM = 1 - V
    WM = 1 - W

    do ID = ID1, ID1 + NV - 1
        VALUES(ID - ID1 + 1) = WM*(VM*(UM*FIELD(ID, GRIDPTR(1, ICELL))&
            &+ U*FIELD(ID, GRIDPTR(2, ICELL)))&
            &+ V*(UM*FIELD(ID, GRIDPTR(3, ICELL))&
            &+ U*FIELD(ID, GRIDPTR(4, ICELL))))&
            &+ W*(VM*(UM*FIELD(ID, GRIDPTR(5, ICELL))&
            &+ U*FIELD(ID, GRIDPTR(6, ICELL)))&
            &+ V*(UM*FIELD(ID, GRIDPTR(7, ICELL))&
            &+ U*FIELD(ID, GRIDPTR(8, ICELL))))
    end do
    return
end subroutine

subroutine LOCATE_GRID_CELL(NX, NY, NZ, XGRID, YGRID, ZGRID,&
    &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
    &BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)

    !   Locates the grid cell in the tree structure containing the
    ! specified point (X0,Y0,Z0), and returns the cell pointer ICELL.
    ! First the base grid cell is found (using NX,NY,NZ, XGRID,YGRID,ZGRID),
    ! and the tree is traced down to find the smallest cell containing
    ! the point.

    implicit none
    integer NX, NY, NZ, NCELLS, ICELL, BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    double precision X0, Y0, Z0
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real GRIDPOS(3, *)
    integer IL, IU, IM, IX, IY, IZ, NXC, NYC, IC, IPTR, DIR, IBITS
    logical BTEST
    double precision XDOMAIN, YDOMAIN

    ! Move point to within horizontal domain if periodic
    if (.not. (btest(BCFLAG, 0) .or. btest(BCFLAG, 2))) then
        XDOMAIN = XGRID(NX + 1) - XGRID(1)
        if (X0 .lt. XGRID(1)) then
            X0 = X0 - XDOMAIN*(int((X0 - XGRID(1))/XDOMAIN) - 1)
        else if (X0 .gt. XGRID(NX + 1)) then
            X0 = X0 - XDOMAIN*int((X0 - XGRID(1))/XDOMAIN)
        end if
    end if
    if (.not. (btest(BCFLAG, 1) .or. btest(BCFLAG, 3))) then
        YDOMAIN = YGRID(NY + 1) - YGRID(1)
        if (Y0 .lt. YGRID(1)) then
            Y0 = Y0 - YDOMAIN*(int((Y0 - YGRID(1))/YDOMAIN) - 1)
        else if (Y0 .gt. YGRID(NY + 1)) then
            Y0 = Y0 - YDOMAIN*int((Y0 - YGRID(1))/YDOMAIN)
        end if
    end if

    ! Find the base grid cell (base grid does not have to be even sp
    !   so do binary search in X, Y, and Z)
    !   Special case for independent pixel: get closest grid column.
    IL = 0
    if (btest(IPFLAG, 0)) then
        IU = NX
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (X0 .ge. 0.5*(XGRID(IM) + XGRID(IM + 1))) then
                IL = IM
            else
                IU = IM
            end if
        end do
        IL = IL + 1
    else
        IU = NX + 1
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (X0 .ge. XGRID(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
    end if
    IX = max(IL, 1)

    IL = 0
    if (btest(IPFLAG, 1)) then
        IU = NY
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (Y0 .ge. 0.5*(YGRID(IM) + YGRID(IM + 1))) then
                IL = IM
            else
                IU = IM
            end if
        end do
        IL = IL + 1
    else
        IU = NY + 1
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (Y0 .ge. YGRID(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
    end if
    IY = max(IL, 1)

    IL = 0
    IU = NZ
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z0 .ge. ZGRID(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)

    ! Special case for open boundary conditions: if outside the 1 to
    !   domain then use the IP grid cells at either end
    NXC = NX
    if (btest(BCFLAG, 0)) then
        NXC = NX + 1
        if (X0 .lt. XGRID(1)) then
            ! X0 = XGRID(1)
            IX = 1
        else if (X0 .gt. XGRID(NX)) then
            ! X0 = XGRID(NX)
            IX = NX + 1
        else
            IX = IX + 1
        end if
    end if
    if (btest(BCFLAG, 2) .and. .not. btest(IPFLAG, 0)) then
        NXC = NX - 1
        IX = min(IX, NXC)
    end if
    NYC = NY
    if (btest(BCFLAG, 1)) then
        NYC = NY + 1
        if (Y0 .lt. YGRID(1)) then
            ! Y0 = YGRID(1)
            IY = 1
        else if (Y0 .gt. YGRID(NY)) then
            ! Y0 = YGRID(NY)
            IY = NY + 1
        else
            IY = IY + 1
        end if
    end if
    if (btest(BCFLAG, 3) .and. .not. btest(IPFLAG, 1)) then
        NYC = NY - 1
        IY = min(IY, NYC)
    end if

    ! Get the base grid cell pointer
    ICELL = IZ + (NZ - 1)*(IY - 1) + (NZ - 1)*NYC*(IX - 1)

    !   Trace down the tree: point to the positive child cell, get the
    ! grid point that is on the negative (X,Y,Z) side, use the
    ! flags to find which dimension (X,Y,Z) the parent cell was
    ! split in, then compare the grid point with the test point
    ! (X0,Y0,Z0) to find which child cell (pos or neg) the test
    ! point is in.
    do while (TREEPTR(2, ICELL) .gt. 0)
        DIR = ibits(int(CELLFLAGS(ICELL)), 2, 2)
        IC = TREEPTR(2, ICELL) + 1
        IPTR = GRIDPTR(1, IC)
        if (DIR .eq. 1) then
            if (X0 .lt. GRIDPOS(1, IPTR)) IC = IC - 1
        else if (DIR .eq. 2) then
            if (Y0 .lt. GRIDPOS(2, IPTR)) IC = IC - 1
        else if (DIR .eq. 3) then
            if (Z0 .lt. GRIDPOS(3, IPTR)) IC = IC - 1
        end if
        ICELL = IC
    end do

    return
end subroutine

subroutine LEGENDRE_ALL(COSSCAT, NLEG, P)
    !   This subroutine computes a set of Legendre polynomials for
    ! a particular scattering angle COSSCAT.  NLEG is the maximum term.
    ! The Legendre functions evaluated at COSSCAT are returned in
    ! P, starting at l=0 and ending with l=NLEG  (NLEG+1 terms).
    implicit none
    integer NLEG
    double precision COSSCAT, P(0:NLEG)
    integer L
    double precision X, PL, PL1, PL2

    X = dble(COSSCAT)
    if (X*X .gt. 1.) stop 'Error [LEGENDRE_ALL]: |COSSCAT| larger than 1.'
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
end subroutine

subroutine YLMALL(TRANSPOSE, MU, PHI, ML, MM, NSTLEG, YR)
    !   This subroutine computes a set of normalized real generalized
    ! spherical harmonic functions (YR) for a particular direction mu,phi.
    ! (Doicu et al. 2013; http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mode.
    ! NSTLEG is the number of matrix elements returned, either 1 for the
    ! 1,1 term or 6 for the 6 unique elements of the 4x4 matrix.
    ! J = 2*(L*(L+1))/2 + M+1  for L<=MM
    ! J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM
    implicit none
    integer ML, MM, NSTLEG
    real MU, PHI
    logical TRANSPOSE
    real YR(NSTLEG, *)
    integer J, M, MABS, L
    double precision X, PI, FCT, P1, P2, P3, COSM, SINM, SIGN
    double precision, allocatable :: DM0(:), DM2P(:), DM2M(:)

    if (NSTLEG .eq. 1) then
        call YLMALL_UNPOL(MU, PHI, ML, MM, YR)
        return
    end if

    allocate (DM0(0:ML), DM2P(0:ML), DM2M(0:ML))
    X = dble(MU)

    PI = DACOS(-1.d0)
    FCT = 1.d0/DSQRT(2.d0*PI)

    ! sign accounting for transpose matrices
    if (.not. TRANSPOSE) then
        SIGN = 1.d0
    else
        SIGN = -1.d0
    end if

    ! ..................................................................
    !                               M = 0
    M = 0
    call WIGNERFCT02P2M_NORMALIZED(X, ML, M, DM0, DM2P, DM2M)
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
        call WIGNERFCT02P2M_NORMALIZED(X, ML, MABS, DM0, DM2P, DM2M)
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
end subroutine

subroutine WIGNERFCT02P2M_NORMALIZED(MIU, NRANK, M,&
    &DM0, DM2P, DM2M)
    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for
    ! |               M >= 0; M1 = 0, 2, -2; N = 0,...,NRANK, and
    ! |              -1 <= MIU <= 1
    ! |  For M1 = -2, we use the symmetry relation
    ! |               DM-M1(MIU) = (-1)**(N+M) * DMM1(-MIU)
    ! ----------------------------------------------------------------
    implicit none
    integer NRANK, M
    double precision MIU
    double precision DM0(0:NRANK), DM2P(0:NRANK), DM2M(0:NRANK)
    integer N0, N
    double precision XP, XM, FACT1, FACT2, FACTP, FACTM, DMM1_N0

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
            DM0(M) = DMM1_N0(XP, M, 0)
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
        DM2P(N0) = DMM1_N0(XP, M, 2)
        DM2M(N0) = DMM1_N0(XM, M, 2)
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
end subroutine

function DMM1_N0(X, M, M1) result(DMM1N0)
    ! ---------------------------------------------------------------
    ! | THE ROUTINE COMPUTES THE WIGNER FUNCTIONS FOR               |
    ! |     M >= 0, M1 >= 0, AND N0 = MAX(M,M1)                     |
    ! ---------------------------------------------------------------
    implicit none
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
end function

subroutine YLMALL_UNPOL(MU, PHI, ML, MM, YR)
    !   This subroutine computes a set of the unpolarized (I-I) element
    ! normalized real generalized spherical harmonic functions (YR) for
    ! a particular direction mu,phi.  (Doicu et al. 2013;
    ! http://dx.doi.org/10.1016/j.jqsrt.2012.12.009)
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mod
    !   J = 2*(L*(L+1))/2 + M+1  for L<=MM
    !   J = (2*MM+1)*L-MM*(2+2*(MM-1))/2 + M+1  for L>MM
    implicit none
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
        call WIGNERFCT_DM0(X, ML, M, DM0)
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
end subroutine

subroutine WIGNERFCT_DM0(X, NRANK, M, DM0)
    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for      |
    ! |               M >= 0; M1 = 0; N = 0,...,NRANK, and           |
    ! |              -1 <= X <= 1                                    |
    ! ----------------------------------------------------------------
    implicit none
    integer NRANK, M
    double precision X
    double precision DM0(0:NRANK)
    integer N
    double precision DM_M10_N0

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
            DM0(M) = DM_M10_N0(X, M)
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
end subroutine

function DM_M10_N0(X, M) result(DMM0N0)
    ! The routine computes the Wigner functions for M >= 0, M1 = 0, a
    ! N0 = M
    implicit none
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
end function

subroutine WIGNERFCT(X, NRANK, M, M1, DMM1)
    ! ----------------------------------------------------------------
    ! |  The routine computes the DMM1N vector coefficients for
    ! |               M >= 0, M1 >= 0,  N = 0,...,NRANK, and
    ! |               -1 < X= COS(BETA) < 1
    ! ----------------------------------------------------------------
    implicit none
    integer NRANK, M, M1
    double precision X, DMM1(0:NRANK)

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
                &*DSQRT(dble(N**2 - M1**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M1**2))
            FACT2 = FACT2*dble(N + 1)/dble(N)

            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    end if
    return
end subroutine

subroutine PLMALL(TRANSPOSE, MU, ML, MM, PRC)
    implicit none
    integer ML, MM
    real MU
    logical TRANSPOSE
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
    if (.not. TRANSPOSE) then
        SIGN = 1.d0
    else
        SIGN = -1.d0
    end if
    ! ..................................................................
    !                            M = 0
    ! ..................................................................
    M = 0
    call WIGNERFCT02P2M_NORMALIZED&
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
        call WIGNERFCT02P2M_NORMALIZED&
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
end subroutine

subroutine PLANCK_FUNCTION(TEMP, UNITS, WAVENO, WAVELEN, PLANCK)
    !    Calculates the Planck blackbody radiance. If UNITS='T' then
    ! using brightness temperature units and the temperature is simply
    ! returned. If UNITS='B' then doing a band integration and the
    ! Planck blackbody radiance in [Watts /(meter^2 ster)] over a
    ! wavenumber range [cm^-1] is returned. Otherwise, the Planck
    ! blackbody radiance in [Watts /(meter^2 ster micron)] for a
    ! temperature in [Kelvins] at a wavelength in [microns] is returned.
    implicit none
    real TEMP, WAVENO(2), WAVELEN, PLANCK
    character(len=1) :: UNITS
    double precision X1, X2, F

    if (UNITS .eq. 'T') then
        PLANCK = TEMP
    else if (UNITS .eq. 'B') then
        if (TEMP .gt. 0.0) then
            X1 = 1.4388d0*WAVENO(1)/TEMP
            X2 = 1.4388d0*WAVENO(2)/TEMP
            call INTEGRATE_PLANCK(X1, X2, F)
            PLANCK = 1.1911d-8*(TEMP/1.4388d0)**4*F
        else
            PLANCK = 0.0
        end if
    else
        if (TEMP .gt. 0.0) then
            PLANCK = 1.1911e8/WAVELEN**5&
                &/(exp(1.4388e4/(WAVELEN*TEMP)) - 1)
        else
            PLANCK = 0.0
        end if
    end if
    return
end subroutine

subroutine INTEGRATE_PLANCK(X1, X2, F)
    ! Returns integral of x**3/(exp(x)-1) from x1 to x2.
    ! Accurate to better than 1 part in 10**9 for integral from
    ! 0 to infinity (pi**4/15).
    implicit none
    double precision X1, X2, F
    double precision INTLOW, INTHIGH
    double precision C

    C = 1.0d0
    if (X1 .lt. C .and. X2 .lt. C) then
        F = INTLOW(X2) - INTLOW(X1)
    else if (X1 .lt. C .and. X2 .ge. C) then
        F = INTLOW(C) - INTLOW(X1) + INTHIGH(C) - INTHIGH(X2)
    else if (X1 .ge. C .and. X2 .ge. C) then
        F = INTHIGH(X1) - INTHIGH(X2)
    else
        stop 'Error [INTEGRATE_PLANCK]: X1 and X2 out of order.'
    end if
    return
end subroutine

double precision function INTLOW(X)
    ! Integral of x**3/(exp(x)-1) from 0 to x.
    ! Accurate for x less than about 1.
    ! Uses Taylor series expansion around x=0.
    implicit none
    double precision X
    integer N
    double precision SUM, F, A(29)
    data A/0.0d0,&
        &0.0d0, 0.333333333333333333333d0,&
        &-0.125d0, 0.016666666666666666667d0,&
        &0.0d0, -0.0001984126984126984127d0,&
        &0.0d0, 0.36743092298647854203d-5,&
        &0.0d0, -0.75156325156325156325d-7,&
        &0.0d0, 0.16059043836821614599d-8,&
        &0.0d0, -0.35227934257916621232d-10,&
        &0.0d0, 0.78720803121674581370d-12,&
        &0.0d0, -0.17840422612224120352d-13,&
        &0.0d0, 0.40886009791799259829d-15,&
        &0.0d0, -0.94559508632959211871d-17,&
        &0.0d0, 0.22036011313440918061d-18,&
        &0.0d0, -0.51683202540046382743d-20,&
        &0.0d0, 0.12188644964239543006d-20/

    SUM = A(4)*X**4
    F = X
    do 100 N = 3, 29, 2
        F = X*X*F
        SUM = SUM + A(N)*F
100 end do
    INTLOW = SUM
    return
end function

double precision function INTHIGH(X)
    ! Integral of x**3/(exp(x)-1) from x to infinity.
    ! Accurate for x greater than about 1.
    implicit none
    double precision X
    integer N
    double precision SUM

    SUM = 0.0d0
    do 100 N = 1, 15
        SUM = SUM + exp(-N*X)&
            &*(X**3/N + 3*X**2/N**2 + 6*X/N**3 + 6.0d0/N**4)
100 end do
    INTHIGH = SUM
    return
end function

subroutine GAUSQUADS(N, XA, WT)
    !   Generates the abscissas (X) and weights (W) for an N point
    ! Gauss-Legendre quadrature.  The XA are returned in this order:
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    implicit none
    integer N
    real XA(*), WT(*)
    integer K, I, J, L
    double precision X, XP, PL, PL1, PL2, DPL, TINY
    parameter(TINY=3.0d-7)

    K = N/2
    if (2*K .ne. N) stop 'Error [GAUSQUADS]: N must be even.'
    do J = 1, K
        X = cos(3.141592654*(J - .25)/(N + .5))
        I = 0
100     continue
        PL1 = 1
        PL = X
        do L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        end do
        DPL = N*(X*PL - PL1)/(X*X - 1)
        XP = X
        X = XP - PL/DPL
        I = I + 1
        if (abs(X - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J) = -X
        XA(J + K) = X
        WT(J) = 2.0d0/((1.0d0 - X*X)*DPL*DPL)
        WT(J + K) = WT(J)
    end do
    return
end subroutine

subroutine DGAUSQUADS(N, XA, WT)
    !   Generates the abscissas (X) and weights (W) for an N point
    ! Double-Gauss-Legendre quadrature.  The XA are returned in this o
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    implicit none
    integer N
    real XA(*), WT(*)
    integer K, I, J, L, N2
    double precision X, XP, PL, PL1, PL2, DPL, TINY
    parameter(TINY=3.0d-7)

    N2 = N/2
    if (2*N2 .ne. N) stop 'Error [DGAUSQUADS]: N must be even.'
    K = (N2 + 1)/2
    do J = 1, K
        X = cos(3.141592654*(J - .25)/(N2 + .5))
        I = 0
100     continue
        PL1 = 1
        PL = X
        do L = 2, N2
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        end do
        DPL = N2*(X*PL - PL1)/(X*X - 1)
        XP = X
        X = XP - PL/DPL
        I = I + 1
        if (abs(X - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J + N2) = (1 + X)/2
        XA(N2 + 1 - J + N2) = (1 - X)/2
        WT(J + N2) = 1/((1 - X*X)*DPL*DPL)
        WT(N2 + 1 - J + N2) = WT(J + N2)
        XA(J) = -XA(J + N2)
        XA(N2 + 1 - J) = -XA(N2 + 1 - J + N2)
        WT(J) = WT(J + N2)
        WT(N2 + 1 - J) = WT(J + N2)
    end do
    return
end subroutine

subroutine TRILIN_INTERP_PROP(X, Y, Z, INIT, NSTLEG, NLEG,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
    !  Trilinearly interpolates the quantities on the input property
    ! grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
    ! ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the
    ! volume coefficients.  Also adds in the separate gaseous absorption.
    ! Divides the phase function Legendre coefficients LEGEN by 2*l+1.
    ! The phase function pointer IPHASE (for tabulated phase functions)
    ! is that of the maximum weighted scattering property grid point.
    ! If INIT=.TRUE. then transfers the tabulated phase functions.

    use globals_shdom

    implicit none
    integer NSTLEG, NLEG
    integer*2 IPHASE
    logical INIT
    real X, Y, Z, TEMP, EXTINCT, ALBEDO, LEGEN(NSTLEG, 0:NLEG, *)
    integer IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
    integer I1, I2, I3, I4, I5, I6, I7, I8, I
    double precision U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
    double precision SCAT1, SCAT2, SCAT3, SCAT4, SCAT5, SCAT6, SCAT7, SCAT8
    double precision SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN
    save EXTMIN, SCATMIN

    if (INIT) then
       ! If there are tabulated phase functions, then transfer them
        do I = 1, NUMPHASE
            do L = 0, NLEG
                do J = 1, NSTLEG
                    LEGEN(J, L, I) = LEGENP(J + NSTLEG*(L + (NLEG + 1)*(I - 1)))/(2*L + 1)
                end do
            end do
        end do
        EXTMIN = 1.0e-5/((ZLEVELS(NPZ) - ZLEVELS(1))/NPZ)
        SCATMIN = 0.1*EXTMIN
        return
    end if

    ! Find the grid location and compute the interpolation factors
    IL = 0
    IU = NPZ
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZLEVELS(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)
    W = dble(Z - ZLEVELS(IZ))/(ZLEVELS(IZ + 1) - ZLEVELS(IZ))
    W = max(min(W, 1.0d0), 0.0d0)
    IX = int((X - XSTART)/DELX) + 1
    if (abs(X - XSTART - NPX*DELX) .lt. 0.01*DELX) IX = NPX
    if (IX .lt. 1 .or. IX .gt. NPX) then
        write (6, *) 'TRILIN: Beyond X domain'
        stop
    end if
    IXP = mod(IX, NPX) + 1
    U = dble(X - XSTART - DELX*(IX - 1))/DELX
    U = max(min(U, 1.0d0), 0.0d0)
    if (U .lt. 1.0d-5) U = 0.0d0
    if (U .gt. 1.0d0 - 1.0d-5) U = 1.0d0
    IY = int((Y - YSTART)/DELY) + 1
    if (abs(Y - YSTART - NPY*DELY) .lt. 0.01*DELY) IY = NPY
    if (IY .lt. 1 .or. IY .gt. NPY) then
        write (6, *) 'TRILIN: Beyond Y domain'
        stop
    end if
    IYP = mod(IY, NPY) + 1
    V = dble(Y - YSTART - DELY*(IY - 1))/DELY
    V = max(min(V, 1.0d0), 0.0d0)
    if (V .lt. 1.0d-5) V = 0.0d0
    if (V .gt. 1.0d0 - 1.0d-5) V = 1.0d0

    F1 = (1 - U)*(1 - V)*(1 - W)
    F2 = U*(1 - V)*(1 - W)
    F3 = (1 - U)*V*(1 - W)
    F4 = U*V*(1 - W)
    F5 = (1 - U)*(1 - V)*W
    F6 = U*(1 - V)*W
    F7 = (1 - U)*V*W
    F8 = U*V*W
    I1 = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
    I2 = IZ + NPZ*(IY - 1) + NPZ*NPY*(IXP - 1)
    I3 = IZ + NPZ*(IYP - 1) + NPZ*NPY*(IX - 1)
    I4 = IZ + NPZ*(IYP - 1) + NPZ*NPY*(IXP - 1)
    I5 = I1 + 1
    I6 = I2 + 1
    I7 = I3 + 1
    I8 = I4 + 1

    ! Trilinearly interpolate the temperature, extinction, scattering
    TEMP = F1*TEMPP(I1) + F2*TEMPP(I2) + F3*TEMPP(I3) + F4*TEMPP(I4)&
        &+ F5*TEMPP(I5) + F6*TEMPP(I6) + F7*TEMPP(I7) + F8*TEMPP(I8)
    EXTINCT = F1*EXTINCTP(I1) + F2*EXTINCTP(I2) + F3*EXTINCTP(I3)&
        &+ F4*EXTINCTP(I4) + F5*EXTINCTP(I5) + F6*EXTINCTP(I6)&
        &+ F7*EXTINCTP(I7) + F8*EXTINCTP(I8)
    SCAT1 = F1*EXTINCTP(I1)*ALBEDOP(I1)
    SCAT2 = F2*EXTINCTP(I2)*ALBEDOP(I2)
    SCAT3 = F3*EXTINCTP(I3)*ALBEDOP(I3)
    SCAT4 = F4*EXTINCTP(I4)*ALBEDOP(I4)
    SCAT5 = F5*EXTINCTP(I5)*ALBEDOP(I5)
    SCAT6 = F6*EXTINCTP(I6)*ALBEDOP(I6)
    SCAT7 = F7*EXTINCTP(I7)*ALBEDOP(I7)
    SCAT8 = F8*EXTINCTP(I8)*ALBEDOP(I8)
    SCATTER = SCAT1 + SCAT2 + SCAT3 + SCAT4 + SCAT5 + SCAT6 + SCAT7 + SCAT8
    if (EXTINCT .gt. EXTMIN) then
        ALBEDO = SCATTER/EXTINCT
    else
        ALBEDO = SCATTER/EXTMIN
    end if
    ! For tabulated phase functions pick the one we are on top of
    !   or the one with the most scattering weight.
    if (NUMPHASE .gt. 0) then
        MAXSCAT = -1.0
        if (SCAT1 .gt. MAXSCAT .or. abs(F1 - 1) .lt. 0.001) then
            MAXSCAT = SCAT1
            IPHASE = IPHASEP(I1)
        end if
        if (SCAT2 .gt. MAXSCAT .or. abs(F2 - 1) .lt. 0.001) then
            MAXSCAT = SCAT2
            IPHASE = IPHASEP(I2)
        end if
        if (SCAT3 .gt. MAXSCAT .or. abs(F3 - 1) .lt. 0.001) then
            MAXSCAT = SCAT3
            IPHASE = IPHASEP(I3)
        end if
        if (SCAT4 .gt. MAXSCAT .or. abs(F4 - 1) .lt. 0.001) then
            MAXSCAT = SCAT4
            IPHASE = IPHASEP(I4)
        end if
        if (SCAT5 .gt. MAXSCAT .or. abs(F5 - 1) .lt. 0.001) then
            MAXSCAT = SCAT5
            IPHASE = IPHASEP(I5)
        end if
        if (SCAT6 .gt. MAXSCAT .or. abs(F6 - 1) .lt. 0.001) then
            MAXSCAT = SCAT6
            IPHASE = IPHASEP(I6)
        end if
        if (SCAT7 .gt. MAXSCAT .or. abs(F7 - 1) .lt. 0.001) then
            MAXSCAT = SCAT7
            IPHASE = IPHASEP(I7)
        end if
        if (SCAT8 .gt. MAXSCAT .or. abs(F8 - 1) .lt. 0.001) then
            MAXSCAT = SCAT8
            IPHASE = IPHASEP(I8)
        end if
    else
        ! Standard property file format: unpolarized
        LEGEN(1, 0, 1) = 1.0
        do L = 1, NLEG
            LEGEN(1, L, 1) = SCAT1*LEGENP(L + (NLEG + 1)*(I1 - 1))&
                &+ SCAT2*LEGENP(L + (NLEG + 1)*(I2 - 1))&
                &+ SCAT3*LEGENP(L + (NLEG + 1)*(I3 - 1))&
                &+ SCAT4*LEGENP(L + (NLEG + 1)*(I4 - 1))&
                &+ SCAT5*LEGENP(L + (NLEG + 1)*(I5 - 1))&
                &+ SCAT6*LEGENP(L + (NLEG + 1)*(I6 - 1))&
                &+ SCAT7*LEGENP(L + (NLEG + 1)*(I7 - 1))&
                &+ SCAT8*LEGENP(L + (NLEG + 1)*(I8 - 1))
            if (SCATTER .gt. SCATMIN) then
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATTER
            else
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATMIN
            end if
            LEGEN(1, L, 1) = LEGEN(1, L, 1)/(2*L + 1)
        end do
    end if

    ! Add in the gaseous absorption to extinction and albedo
    if (NZCKD .gt. 0) then
        IL = 1
        IU = NZCKD
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (Z .le. ZCKD(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
        I = min(max(IL, 1), NZCKD - 1)
        F = (Z - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
        F = min(max(F, 0.0d0), 1.0d0)
        KG = (1.0 - F)*GASABS(I) + F*GASABS(I + 1)
        if (EXTINCT + KG .gt. 0.0) then
            ALBEDO = ALBEDO*EXTINCT/(EXTINCT + KG)
        else
            ALBEDO = 0.0
        end if
        EXTINCT = EXTINCT + KG
    end if
    return
end subroutine

subroutine DIRECT_BEAM_PROP(INIT, XI, YI, ZI, BCFLAG, IPFLAG,&
    &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
    &DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)

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

    use globals_shdom

    implicit none
    integer INIT, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, SIDE
    logical DELTAM, VALIDBEAM
    real XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
    real DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH

    integer IX, IY, IZ, JZ, IL, IM, IU
    integer I, J, K, L, IPH, IP, JP, I1, I2, I3, I4
    integer IPDIRECT, DI, DJ, DK
    logical CONSTX, CONSTY, HITBOUNDARY, BTEST
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
        if (.not. allocated(GASEXT)) allocate (GASEXT(NPZ), EXTMIN(NPZ), EXTMAX(NPZ))
        ! Get the gaseous extinction at the property grid levels
        do IZ = 1, NPZ
            if (NZCKD .gt. 0) then
                IL = 1
                IU = NZCKD
                do while (IU - IL .gt. 1)
                    IM = (IU + IL)/2
                    if (ZLEVELS(IZ) .le. ZCKD(IM)) then
                        IL = IM
                    else
                        IU = IM
                    end if
                end do
                I = min(max(IL, 1), NZCKD - 1)
                W0 = (ZLEVELS(IZ) - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
                W0 = min(max(W0, 0.0d0), 1.0d0)
                GASEXT(IZ) = (1.0 - W0)*GASABS(I) + W0*GASABS(I + 1)
            else
                GASEXT(IZ) = 0.0
            end if
        end do
        ! First make the property grid extinction field, on which all
        !   the direct beam paths will be computed.
        IP = 0
        do IX = 1, NPX
            do IY = 1, NPY
                do IZ = 1, NPZ
                    IP = IP + 1
                    EXTINCT = EXTINCTP(IP)
                    ALBEDO = ALBEDOP(IP)
                    ! Add in the gaseous absorption to extinction
                    if (GASEXT(IZ) .gt. 0.0) then
                        ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                        EXTINCT = EXTINCT + GASEXT(IZ)
                    end if
                    ! Do the Delta-M scaling if needed
                    if (DELTAM) then
                        if (NUMPHASE .gt. 0) then
                            IPH = IPHASEP(IP)
                        else
                            IPH = IP
                        end if
                        L = ML + 1
                        F = LEGENP(1 + NSTLEG*(L + (NLEG + 1)*(IPH - 1)))/(2*L + 1)
                        EXTINCT = (1.0 - ALBEDO*F)*EXTINCT
                    end if
                    EXTDIRP(IP) = EXTINCT
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
        EPSZ = 1.0e-6*(ZLEVELS(NPZ) - ZLEVELS(1))
        EPSS = 1.0e-3*(ZLEVELS(NPZ) - ZLEVELS(1))/NPZ
        if (.not. btest(IPDIRECT, 0)) EPSS = max(EPSS, 1.0d-4*DELX)
        if (.not. btest(IPDIRECT, 1)) EPSS = max(EPSS, 1.0d-4*DELY)
        DELXD = dble(DELX)
        DELYD = dble(DELY)
        EPSS = max(0.001*DELXD, 0.001*DELYD, EPSS)
        XDOMAIN = DELXD*NPX
        if (btest(BCFLAG, 2)) XDOMAIN = DELXD*(NPX - 1)
        YDOMAIN = DELYD*NPY
        if (btest(BCFLAG, 3)) YDOMAIN = DELYD*(NPY - 1)

        ! Find the Z level above which the medium is plane-parallel
        EXTBEAMCUT = 1.0e-4
        do IZ = 1, NPZ
            EXTMIN(IZ) = 1.0e20
            EXTMAX(IZ) = 0.0
        end do
        IP = 0
        do IX = 1, NPX
            do IY = 1, NPY
                do IZ = 1, NPZ
                    IP = IP + 1
                    EXTMIN(IZ) = min(EXTINCTP(IP), EXTMIN(IZ))
                    EXTMAX(IZ) = max(EXTINCTP(IP), EXTMAX(IZ))
                end do
            end do
        end do
        JZ = 0
        do IZ = 1, NPZ
            if (EXTMAX(IZ) - EXTMIN(IZ) .gt. EXTBEAMCUT) JZ = IZ
        end do
        JZ = min(NPZ, JZ + 1)
        UNIFORMZLEV = ZLEVELS(JZ)
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
    X = XI - XSTART
    Y = YI - YSTART

    ! Find the grid location
    IL = 0
    IU = NPZ
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZLEVELS(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    K = max(IL, 1)
    I = int(X/DELXD) + 1
    if (I .gt. NPX .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) I = NPX
    if (I .lt. 1 .or. I .gt. NPX) then
        write (6, *) 'DIRECT_BEAM_PROP: Beyond X domain', I, XI, YI, ZI
        stop
    end if
    J = int(Y/DELYD) + 1
    if (J .gt. NPY .and. abs(Y - YDOMAIN) .lt. 0.001*DELYD) J = NPY
    if (J .lt. 1 .or. J .gt. NPY) then
        write (6, *) 'DIRECT_BEAM_PROP: Beyond Y domain', J, XI, YI, ZI
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
        &.or. abs(X - (NPX - 1)*DELXD) .lt. 0.01*DELXD)) CONSTX = .true.
    if (btest(BCFLAG, 1) .and. (abs(Y) .lt. 0.01*DELYD&
        &.or. abs(Y - (NPY - 1)*DELYD) .lt. 0.01*DELYD)) CONSTY = .true.

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
    do while (.not. HITBOUNDARY .and. abs(ZE - ZLEVELS(NPZ)) .gt. EPSZ)
        IP = I + 1
        if (I .eq. NPX) then
            if (btest(BCFLAG, 0) .or. btest(BCFLAG, 2)) then
                IP = NPX
            else
                IP = 1
            end if
        end if
        JP = J + 1
        if (J .eq. NPY) then
            if (btest(BCFLAG, 1) .or. btest(BCFLAG, 3)) then
                JP = NPY
            else
                JP = 1
            end if
        end if
        X0 = DELXD*(I - 1)
        X1 = X0 + DELXD
        Y0 = DELYD*(J - 1)
        Y1 = Y0 + DELYD
        Z0 = ZLEVELS(K)
        Z1 = ZLEVELS(K + 1)
        if  (I .lt. 1 .or. I .gt. NPX .or.&
            &J .lt. 1 .or. J .gt. NPY .or.&
            &K .lt. 1 .or. K .ge. NPZ) then
            write (6, '(A,3I4)') 'DIRECT_BEAM_PROP: beyond grid!', I, J, K
            write (6, '(1(2X,3F9.5))') X, Y, Z
            stop
        end if
        ! Get the eight corner extinction values
        I1 = K + NPZ*(J - 1) + NPZ*NPY*(I - 1)
        I2 = K + NPZ*(J - 1) + NPZ*NPY*(IP - 1)
        I3 = K + NPZ*(JP - 1) + NPZ*NPY*(I - 1)
        I4 = K + NPZ*(JP - 1) + NPZ*NPY*(IP - 1)
        E1 = EXTDIRP(I1)
        E2 = EXTDIRP(I2)
        E3 = EXTDIRP(I3)
        E4 = EXTDIRP(I4)
        E5 = EXTDIRP(I1 + 1)
        E6 = EXTDIRP(I2 + 1)
        E7 = EXTDIRP(I3 + 1)
        E8 = EXTDIRP(I4 + 1)

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
                    I = NPX
                    XOFFS = XDOMAIN
                end if
            else if (I .ge. NPX .and. btest(BCFLAG, 2)) then
                HITBOUNDARY = .true.
                SIDE = 2
            else if (I .eq. NPX + 1) then
                if (btest(BCFLAG, 0)) then
                    I = NPX
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
                    J = NPY
                    YOFFS = YDOMAIN
                end if
            else if (J .ge. NPY .and. btest(BCFLAG, 3)) then
                HITBOUNDARY = .true.
                SIDE = 4
            else if (J .eq. NPY + 1) then
                if (btest(BCFLAG, 1)) then
                    J = NPY
                    CONSTY = .true.
                else
                    J = 1
                    YOFFS = -YDOMAIN
                end if
            end if
        end if
        if (SO .lt. -EPSS) then
            write (6, *) 'DIRECT_BEAM_PROP: SO<0', X, Y, Z,&
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
    XO = XE + XSTART; YO = YE + YSTART; ZO = ZE
    VALIDBEAM = abs(ZE - ZLEVELS(NPZ)) .lt. EPSZ
    if (VALIDBEAM) SIDE = 6
    return
end subroutine

subroutine SSORT(X, Y, N, KFLAG)
    !***BEGIN PROLOGUE  SSORT
    !***PURPOSE  Sort an array and optionally make the same interchanges in
    !            an auxiliary array.  The array may be sorted in increasing
    !            or decreasing order.  A slightly modified QUICKSORT
    !            algorithm is used.
    !***LIBRARY   SLATEC
    !***CATEGORY  N6A2B
    !***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
    !***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
    !***AUTHOR  Jones, R. E., (SNLA)
    !           Wisniewski, J. A., (SNLA)
    !***DESCRIPTION
    !
    !   SSORT sorts array X and optionally makes the same interchanges in
    !   array Y.  The array X may be sorted in increasing order or
    !   decreasing order.  A slightly modified quicksort algorithm is used.
    !
    !   Description of Parameters
    !      X - array of values to be sorted   (usually abscissas)
    !      Y - array to be (optionally) carried along
    !      N - number of values in array X to be sorted
    !      KFLAG - control parameter
    !            =  2  means sort X in increasing order and carry Y along.
    !            =  1  means sort X in increasing order (ignoring Y)
    !            = -1  means sort X in decreasing order (ignoring Y)
    !            = -2  means sort X in decreasing order and carry Y along.
    !
    !***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
    !                 for sorting with minimal storage, Communications of
    !                 the ACM, 12, 3 (1969), pp. 185-187.
    !***END PROLOGUE  SSORT
    !     .. Scalar Arguments ..
    integer KFLAG, N
    ! .. Array Arguments ..
    ! real X(*), Y(*)
    real X(*)
    integer Y(*)
    ! .. Local Scalars ..
    ! real R, T, TT, TTY, TY
    real R, T, TT
    integer TY, TTY
    integer I, IJ, J, K, KK, L, M, NN
    ! .. Local Arrays ..
    integer IL(51), IU(51)
    ! .. External Subroutines ..
    ! .. Intrinsic Functions ..
    intrinsic ABS, INT
    !***First executable statement  SSORT
    NN = N
    if (NN .lt. 1) then
        stop 'Error [SSORT]: The number of values to be sorted is not positive.'
    end if

    KK = abs(KFLAG)
    if (KK .ne. 1 .and. KK .ne. 2) then
        stop 'Error [SSORT]: The sort control parameter, K, is not 2, 1, -1, or -2.'
    end if

    ! Alter array X to get decreasing order if needed
    if (KFLAG .le. -1) then
        do 10 I = 1, NN
            X(I) = -X(I)
10      continue
     end if

    if (KK .eq. 2) GO TO 100

    ! Sort X only
    M = 1
    I = 1
    J = NN
    R = 0.375e0

20  if (I .eq. J) GO TO 60
    if (R .le. 0.5898437e0) then
        R = R + 3.90625e-2
    else
        R = R - 0.21875e0
    end if

30  K = I

    ! Select a central element of the array and save it in location T
    IJ = I + int((J - I)*R)
    T = X(IJ)

    ! If first element of array is greater than T, interchange with T
    if (X(I) .gt. T) then
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
    end if
    L = J

    ! If last element of array is less than than T, interchange with T
    if (X(J) .lt. T) then
        X(IJ) = X(J)
        X(J) = T
        T = X(IJ)

        ! If first element of array is greater than T, interchange with T
        if (X(I) .gt. T) then
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
        end if
    end if

    ! Find an element in the second half of the array which is smaller
    ! than T
40  L = L - 1
    if (X(L) .gt. T) GO TO 40

    ! Find an element in the first half of the array which is greater
    ! than T
50  K = K + 1
    if (X(K) .lt. T) GO TO 50

    ! Interchange these elements
    if (K .le. L) then
        TT = X(L)
        X(L) = X(K)
        X(K) = TT
        GO TO 40
    end if

    ! Save upper and lower subscripts of the array yet to be sorted
    if (L - I .gt. J - K) then
        IL(M) = I
        IU(M) = L
        I = K
        M = M + 1
    else
        IL(M) = K
        IU(M) = J
        J = L
        M = M + 1
    end if
    GO TO 70

    ! Begin again on another portion of the unsorted array
60  M = M - 1
    if (M .eq. 0) GO TO 190
    I = IL(M)
    J = IU(M)

70  if (J - I .ge. 1) GO TO 30
    if (I .eq. 1) GO TO 20
    I = I - 1

80  I = I + 1
    if (I .eq. J) GO TO 60
    T = X(I + 1)
    if (X(I) .le. T) GO TO 80
    K = I

90  X(K + 1) = X(K)
    K = K - 1
    if (T .lt. X(K)) GO TO 90
    X(K + 1) = T
    GO TO 80

    ! Sort X and carry Y along
100 M = 1
    I = 1
    J = NN
    R = 0.375e0

110 if (I .eq. J) GO TO 150
    if (R .le. 0.5898437e0) then
        R = R + 3.90625e-2
    else
        R = R - 0.21875e0
    end if

120 K = I

    ! Select a central element of the array and save it in location T
    IJ = I + int((J - I)*R)
    T = X(IJ)
    TY = Y(IJ)

    ! If first element of array is greater than T, interchange with T
    if (X(I) .gt. T) then
        X(IJ) = X(I)
        X(I) = T
        T = X(IJ)
        Y(IJ) = Y(I)
        Y(I) = TY
        TY = Y(IJ)
    end if
    L = J

    ! If last element of array is less than T, interchange with T
    if (X(J) .lt. T) then
        X(IJ) = X(J)
        X(J) = T
        T = X(IJ)
        Y(IJ) = Y(J)
        Y(J) = TY
        TY = Y(IJ)

    ! If first element of array is greater than T, interchange with T
        if (X(I) .gt. T) then
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
        end if
    end if

    ! Find an element in the second half of the array which is smaller
    ! than T
130 L = L - 1
    if (X(L) .gt. T) GO TO 130

    ! Find an element in the first half of the array which is greater
    ! than T
140 K = K + 1
    if (X(K) .lt. T) GO TO 140

    ! Interchange these elements
    if (K .le. L) then
        TT = X(L)
        X(L) = X(K)
        X(K) = TT
        TTY = Y(L)
        Y(L) = Y(K)
        Y(K) = TTY
        GO TO 130
    end if

    ! Save upper and lower subscripts of the array yet to be sorted
    if (L - I .gt. J - K) then
        IL(M) = I
        IU(M) = L
        I = K
        M = M + 1
    else
        IL(M) = K
        IU(M) = J
        J = L
        M = M + 1
    end if
    GO TO 160

    ! Begin again on another portion of the unsorted array
150 M = M - 1
    if (M .eq. 0) GO TO 190
    I = IL(M)
    J = IU(M)

160 if (J - I .ge. 1) GO TO 120
    if (I .eq. 1) GO TO 110
    I = I - 1

170 I = I + 1
    if (I .eq. J) GO TO 150
    T = X(I + 1)
    TY = Y(I + 1)
    if (X(I) .le. T) GO TO 170
    K = I

180 X(K + 1) = X(K)
    Y(K + 1) = Y(K)
    K = K - 1
    if (T .lt. X(K)) GO TO 180
    X(K + 1) = T
    Y(K + 1) = TY
    GO TO 170

    ! Clean up
190 if (KFLAG .le. -1) then
        do 200 I = 1, NN
            X(I) = -X(I)
200     continue
    end if
    return
end subroutine
