module shd_grid

  use globals

  implicit none
  private

  public :: shd_grid__makeNewGridPnt
  public :: shd_grid__getNextCell
  public :: shd_grid__testGridSmooth
  public :: shd_grid__testCellSplit
  public :: shd_grid__splitGrid
  public :: shd_grid__makeSweepOrder

contains

subroutine shd_grid__makeNewGridPnt(IDIR, ICELL, NEWCELL, NPTS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, NEWPOINTS)

    !   Sets up the pointers to the grid points for the two new cells.
    ! Creates the four new points if they don't already exist, changing
    ! NPTS and GRIDPOS accordingly.  Returns the NEWPOINTS(3,4) arrays
    ! that gives information about the new points: NEWPOINTS(3,) is the
    ! point address of the new points (or zero), and NEWPOINTS(1,) and
    ! NEWPOINTS(2,) are the "parent" grid points these should be
    ! interpolated from.

    integer IDIR, ICELL, NEWCELL, NPTS, NEWPOINTS(3, 4)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer I, K, I1, I2, IP1, IP2, IPT, IPMATCH
    integer IFACE1, IFACE2, ICELL2
    integer GRIDCORNER(2, 4, 3), FACEGRID(2, 4, 3)
    real XP, YP, ZP
    data GRIDCORNER/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        &1, 5, 2, 6, 3, 7, 4, 8/
    data FACEGRID/3, 5, 4, 5, 3, 6, 4, 6, 1, 5, 2, 5, 1, 6, 2, 6,&
        &1, 3, 2, 3, 1, 4, 2, 4/

    !   In general, when subdividing a cell four new points may be made.
    ! Each of these points belongs to two faces of the original cell.
    ! Use the shd_comp__matchGridPnt routine to try to find matches in cells
    ! across the two faces that contain the point; also cross the face of
    ! one of these to get the third bounding cell.  If any of the appropriate
    ! gridpoints of these three cells has the same coordinates then the point
    ! already exists, otherwise create the point.

    ! Temporarily set new cells' grid pointer to the parent cell's ones
    do K = 1, 8
        GRIDPTR(K, NEWCELL) = GRIDPTR(K, ICELL)
        GRIDPTR(K, NEWCELL + 1) = GRIDPTR(K, ICELL)
    end do

    ! Loop over the four edges with potentially new points
    do I = 1, 4
        ! Get the grid corners and coordinates for this edge and split dir
        I1 = GRIDCORNER(1, I, IDIR)
        I2 = GRIDCORNER(2, I, IDIR)
        IP1 = GRIDPTR(I1, ICELL)
        IP2 = GRIDPTR(I2, ICELL)
        XP = (GRIDPOS(1, IP1) + GRIDPOS(1, IP2))/2
        YP = (GRIDPOS(2, IP1) + GRIDPOS(2, IP2))/2
        ZP = (GRIDPOS(3, IP1) + GRIDPOS(3, IP2))/2

        ! Look for gridpoints matching this new point in the three cell
        !   that border this edge
        IFACE1 = FACEGRID(1, I, IDIR)
        call shd_comp__matchGridPnt(XP, YP, ZP, ICELL, IFACE1, GRIDPOS,&
            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
        if (IPMATCH .eq. 0) then
            IFACE2 = FACEGRID(2, I, IDIR)
            call shd_comp__matchGridPnt(XP, YP, ZP, ICELL, IFACE2, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
        end if
        if (IPMATCH .eq. 0) then
            ICELL2 = abs(NEIGHPTR(IFACE1, ICELL))
            if (ICELL2 .gt. 0) then
                call shd_comp__matchGridPnt(XP, YP, ZP, ICELL2, IFACE2, GRIDPOS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
            end if
        end if

        ! Set some new cells' grid pointers to the new or matching one
        if (IPMATCH .eq. 0) then
            NPTS = NPTS + 1
            GRIDPTR(I2, NEWCELL) = NPTS
            GRIDPTR(I1, NEWCELL + 1) = NPTS
            GRIDPOS(1, NPTS) = XP
            GRIDPOS(2, NPTS) = YP
            GRIDPOS(3, NPTS) = ZP
            NEWPOINTS(1, I) = IP1
            NEWPOINTS(2, I) = IP2
            NEWPOINTS(3, I) = NPTS
        else
            GRIDPTR(I2, NEWCELL) = IPMATCH
            GRIDPTR(I1, NEWCELL + 1) = IPMATCH
            NEWPOINTS(3, I) = 0
        end if
    end do

    return
end subroutine shd_grid__makeNewGridPnt

subroutine shd_grid__getNextCell(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXT)

    !   Gets the next cell (INEXT) to go to that is adjacent to this face
    ! IFACE (1-6) has which face (-X,+X,-Y,+Y,-Z,+Z) of the current cell
    ! (ICELL) we are at, and JFACE (1-3) is the direction (X,Y,Z).
    ! The current location (XE,YE,ZE) on the face is used to find the cell
    ! if there are several adjacent to the current cell.

    integer IFACE, JFACE, ICELL, INEXT
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    double precision XE, YE, ZE
    integer IC, IC1, DIR

    ! Get the pointer to the neighbor cell
    INEXT = NEIGHPTR(IFACE, ICELL)
    if (INEXT .lt. 0) then
        ! If neighbor pointer is negative then we are going into a cell
        !   that has subcells, and the absolute value of the pointer is the
        !   cell to start the tree search.
        IC = abs(INEXT)
        ! Go down the tree until there is no children
        do while (TREEPTR(2, IC) .gt. 0)
            ! Get the direction cell is split
            DIR = ibits(int(CELLFLAGS(IC)), 2, 2)
            if (DIR .eq. 0) stop 'Error [_moveNextCell]: No split direction.'
            ! If split is in direction of face then choose the near child cell
            !   which is the negative side cell if this is a positive face.
            IC1 = TREEPTR(2, IC)
            if (DIR .eq. JFACE) then
                IC = IC1 + 1 - mod(IFACE - 1, 2)
            else
                ! If split is not in direction of face then find which of two cells
                !   has a face that contains the current location.
                ! Compare the appropriate coordinate (X,Y,Z) of the grid point
                !   on the split line with the current location.
                IC = IC1
                if (DIR .eq. 1) then
                    if (XE .gt. GRIDPOS(1, GRIDPTR(8, IC1))) IC = IC + 1
                else if (DIR .eq. 2) then
                    if (YE .gt. GRIDPOS(2, GRIDPTR(8, IC1))) IC = IC + 1
                else
                    if (ZE .gt. GRIDPOS(3, GRIDPTR(8, IC1))) IC = IC + 1
                end if
            end if
        end do
        INEXT = IC
    end if

    return
end subroutine shd_grid__getNextCell

subroutine shd_grid__testGridSmooth(GRIDPTR, GRIDPOS, NEIGHPTR,&
    &TREEPTR, CELLFLAGS, SHPTR, ICELL, IDIR)
    !   Tests the grid cell ICELL for splitting in order to improve the
    ! adaptive grid smoothness.  The direction to split the cell is returned
    ! in IDIR (1-X, 2-Y, 3-Z, 0-don't).
    ! The criterion for splitting the cell is
    !   1) cell has finer grid cells on either side
    !   2) a neighboring cell has grid spacing more than 3 times finer

    integer ICELL, IDIR
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *), SHPTR(*)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer ID, IFACE, IP, IP1, IP2, IN, I, J, ISUM, IBITS
    integer INCELL(2), DIRSPLIT(2)
    integer FACECORNER(4, 6), EDGECORNER(2, 3)
    double precision XE, YE, ZE
    real RATIO, SIZERATIO, CURGRIDSIZE(3), GRIDSIZE(2), GRIDSIZEINV
    data FACECORNER/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/
    data EDGECORNER/1, 2, 1, 3, 1, 5/

    IDIR = 0
    SIZERATIO = 1.0
    ! Loop over the three possible directions to do the split (X,Y,Z)
    do ID = 1, 3
        ! Get the grid points
        IP1 = GRIDPTR(EDGECORNER(1, ID), ICELL)
        IP2 = GRIDPTR(EDGECORNER(2, ID), ICELL)
        CURGRIDSIZE(ID) = 1.0e20
        if (IP1 .ne. IP2) then
            CURGRIDSIZE(ID) = abs(GRIDPOS(ID, IP1) - GRIDPOS(ID, IP2))
            GRIDSIZEINV = 1.0/CURGRIDSIZE(ID)
            ! Loop over two opposite faces
            do J = 1, 2
                DIRSPLIT(J) = 0
                IFACE = 2*(ID - 1) + J
                ! Get location of center of cell's face
                XE = 0.0
                YE = 0.0
                ZE = 0.0
                do I = 1, 4
                    IP = GRIDPTR(FACECORNER(I, IFACE), ICELL)
                    XE = XE + GRIDPOS(1, IP)*0.25
                    YE = YE + GRIDPOS(2, IP)*0.25
                    ZE = ZE + GRIDPOS(3, IP)*0.25
                end do
                ! Get the neighboring cell and its size in this direction
                call shd_grid__getNextCell(XE, YE, ZE, IFACE, ID, ICELL, GRIDPOS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INCELL(J))
                ! If the neighbor does not exist or is a base cell or
                !   is an independent pixel cell in X or Y then set
                !   its grid size the same as current cell (so don't split
                if (INCELL(J) .eq. 0) then
                    GRIDSIZE(J) = CURGRIDSIZE(ID)
                else
                    GRIDSIZE(J) = abs(GRIDPOS(ID, GRIDPTR(1, INCELL(J)))&
                        &- GRIDPOS(ID, GRIDPTR(8, INCELL(J))))

                    if (TREEPTR(1, INCELL(J)) .eq. 0) then
                        GRIDSIZE(J) = CURGRIDSIZE(ID)
                    end if
                    if (ID .le. 2 .and. btest(CELLFLAGS(INCELL(J)), ID - 1)) then
                        GRIDSIZE(J) = CURGRIDSIZE(ID)
                    end if
                end if
                ! If there are multiple neighboring cells then get
                !   direction the neighbor cell was split
                IN = NEIGHPTR(IFACE, ICELL)
                if (IN .lt. 0) then
                    DIRSPLIT(J) = ibits(int(CELLFLAGS(abs(IN))), 2, 2)
                end if
            end do

            ! Split cell this way if cells on opposite side are smaller in
            !   this direction OR opposite neighbor cells are split in
            !   same other direction.
            if (GRIDSIZE(1)*GRIDSIZEINV .lt. 0.75 .and.&
                &GRIDSIZE(2)*GRIDSIZEINV .lt. 0.75) then
                IDIR = ID
            end if
            if (DIRSPLIT(1) .gt. 0 .and. DIRSPLIT(1) .ne. ID .and.&
                &DIRSPLIT(1) .eq. DIRSPLIT(2)) then
                IDIR = DIRSPLIT(1)
            end if

            ! Split cell this way if ratio of neighboring grid cells' size
            !   to this cells size (in ID direction) is less than 0.4.
            RATIO = min(GRIDSIZE(1)*GRIDSIZEINV, GRIDSIZE(2)*GRIDSIZEINV)
            if (RATIO .lt. 0.4 .and. RATIO .lt. SIZERATIO) then
                ! Check if there is any source function for this cell, before
                !   deciding to split it (no point in dividing clear sky).
                ISUM = 0
                do I = 1, 8
                    IP = GRIDPTR(I, ICELL)
                    ISUM = ISUM + SHPTR(IP + 1) - SHPTR(IP)
                end do
                if (ISUM .gt. 0) then
                    IDIR = ID
                    SIZERATIO = RATIO
                end if
            end if

        end if
    end do

    return
end subroutine shd_grid__testGridSmooth

subroutine shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
    &SHPTR, NSTOKES, SOURCE, ICELL,&
    &ADAPTCRIT, MAXADAPT, IDIR)

    !   Computes the adaptive cell splitting criterion for the X, Y, and Z
    ! directions for cell (ICELL).  The three criterion are returned in
    ! ADAPTCRIT, the max criterion is returned in MAXADAPT, and the
    ! direction of max criterion in IDIR.
    ! The adaptive cell criterion is computed by averaging the division
    ! criterion over the four pairs of grid points that cross the potential
    ! split direction. The cell division criterion is based on the
    ! integration of the source function across the cell (including the
    ! optical depth/transmission factor) that gives the delta radiance.
    ! The criterion is actually the part of the integral that depends on
    ! the difference in the extinction times source function across a cell.
    ! The rms over all angles of the extinction/source function difference
    ! is computed in spherical harmonic space by summing the square of
    ! the difference over all terms.

    integer NSTOKES, ICELL, IDIR
    integer SHPTR(*), GRIDPTR(8, *)
    real GRIDPOS(3, *), EXTINCT(*), SOURCE(NSTOKES, *)
    real ADAPTCRIT(3), MAXADAPT
    integer ID, IE, IP1, IP2, IS1, IS2, NS1, NS2, NS
    integer J, NUM
    integer GRIDCORNER(2, 4, 3)
    real SPLIT, SUM1, JAY, TAU, EXT
    real C0
    data GRIDCORNER/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        &1, 5, 2, 6, 3, 7, 4, 8/
    data C0/0.282095/

    MAXADAPT = -1.0
    ! Loop over the three possible directions to do the split (X,Y,Z)
    do ID = 1, 3
        SUM1 = 0.0
        NUM = 0
        ! Loop over four edges that cross the ID direction
        do IE = 1, 4
            ! Get the grid points
            IP1 = GRIDPTR(GRIDCORNER(1, IE, ID), ICELL)
            IP2 = GRIDPTR(GRIDCORNER(2, IE, ID), ICELL)
            if (IP1 .ne. IP2) then
                NUM = NUM + 1
                IS1 = SHPTR(IP1)
                IS2 = SHPTR(IP2)
                NS1 = SHPTR(IP1 + 1) - IS1
                NS2 = SHPTR(IP2 + 1) - IS2
                JAY = 0.0
                NS = min(NS1, NS2)
                !   Compute the source function rms difference integrated
                ! over all angles, which is sum of squares of spherical
                ! harmonic term differences.
                do J = 1, NS
                    JAY = JAY + (EXTINCT(IP2)*SOURCE(1, IS2 + J)&
                        &- EXTINCT(IP1)*SOURCE(1, IS1 + J))**2
                end do
                do J = NS + 1, NS1
                    JAY = JAY + (EXTINCT(IP1)*SOURCE(1, IS1 + J))**2
                end do
                do J = NS + 1, NS2
                    JAY = JAY + (EXTINCT(IP2)*SOURCE(1, IS2 + J))**2
                end do
                EXT = 0.5*(EXTINCT(IP1) + EXTINCT(IP2))
                if (EXT .gt. 0.0) then
                    JAY = C0*sqrt(JAY)/EXT
                else
                    JAY = 0.0
                end if
                TAU = abs(EXT*(GRIDPOS(ID, IP2) - GRIDPOS(ID, IP1)))
                SPLIT = abs(JAY)*(1 - exp(-TAU))
                SUM1 = SUM1 + SPLIT
            end if
        end do

        if (NUM .gt. 0) then
            ADAPTCRIT(ID) = SUM1/NUM
        else
            ADAPTCRIT(ID) = 0.0
        end if
        if (ADAPTCRIT(ID) .gt. MAXADAPT) then
            MAXADAPT = ADAPTCRIT(ID)
            IDIR = ID
        end if
    end do

    return
end subroutine shd_grid__testCellSplit

subroutine shd_grid__splitGrid(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
    &DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
    &NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, NX, XGRID, NY, YGRID,&
    &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
    &TEMP, PLANCK, DIRFLUX, SRCTYPE,&
    &SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    &UNITS, WAVENO, WAVELEN, ADAPTCRIT, ADAPTIND)

    !   Splits the cells that have a cell dividing criterion greater than
    ! CURSPLITACC if DOSPLIT is true.  The current splitting criterion
    ! achieved is returned in SPLITCRIT.  If we are at the end of the
    ! grid point, grid cell, or spherical harmonic arrays (sizes given
    ! by MAXIG, MAXIC, MAXIV, MAXIDO) then the OUTOFMEM flag is returned true.

    use shd_ext_mpi

    integer MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX
    integer NPTS, NCELLS, BCFLAG, IPFLAG
    integer NX, NY, NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    real ADAPTIND(*)
    integer*2 CELLFLAGS(*), IPHASE(*)
    logical DELTAM, ACCELFLAG, DOSPLIT, OUTOFMEM
    real CURSPLITACC, SPLITCRIT
    real SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(NSTLEG, NLM)
    real WAVENO(2), WAVELEN
    real GRIDPOS(3, *), XGRID(*), YGRID(*), ADAPTCRIT(*)
    real EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    real RADIANCE(NSTOKES, *), SOURCE(NSTOKES, *)
    real TEMP(*), PLANCK(*), DIRFLUX(*)
    character(len=1) :: SRCTYPE, UNITS
    integer ICELL, ICELL1, IDIR, I, N, NEWPOINTS(3, 4), OLDNPTS
    integer MAXCELLS, MAXPTS, MAXWORK, MAXSH, MAXRAD
    logical OUTOFMEM0
    real ADAPT(3), FRAC

    if (DOSPLIT) then
        OLDNPTS = NPTS
        ! Do the cells in batchs, so the newly split cells can be tested
        !   for further splitting
        ICELL1 = 1
        do while (ICELL1 .le. NCELLS)
            ! For the current batch of cells, compute the adaptive grid
            !   splitting criterion for 3 directions for each cell.
            N = 0
            do ICELL = ICELL1, NCELLS
                if (TREEPTR(2, ICELL) .eq. 0) then
                    ! If this is a end node cell then compute adaptive criterion
                    call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
                        &SHPTR, NSTOKES, SOURCE, ICELL,&
                        &ADAPT, ADAPTCRIT(N + 1), IDIR)
                    ADAPTIND(N+1) = 4*ICELL + IDIR
                    N = N + 1
                end if
            end do
            ICELL1 = NCELLS + 1

            ! Sort the adaptive criterion in decreasing order
            call shd_comp__ssort(ADAPTCRIT, ADAPTIND, N, -2)

            ! Figure how many cells we can split, leaving some room
            !   for cells split during grid smoothing
            FRAC = 0.03
            MAXCELLS = MAXIC - FRAC*(MAXIC - NCELLS) - 2
            MAXPTS = MAXIG - FRAC*(MAXIG - NPTS) - 4
            MAXWORK = MAXIDO - FRAC*(MAXIDO - NPHI0MAX*NPTS) - 4*NPHI0MAX
            MAXSH = MAXIV - FRAC*(MAXIV - SHPTR(NPTS + 1)) - 4*NLM
            MAXRAD = MAXIV + MAXIG - FRAC*(MAXIV + MAXIG - RSHPTR(NPTS + 1)) - 4*NLM

            ! Go down the list from highest adaptive cell criterion on down
            OUTOFMEM0 = OUTOFMEM
            I = 1
            do while (ADAPTCRIT(I) .gt. CURSPLITACC .and. .not. OUTOFMEM0&
                &.and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    &NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    &SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    &RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM0 = .true.
                else
                    ICELL = int(ADAPTIND(I))/4
                    IDIR = iand(int(ADAPTIND(I)), 3)
                    ! Create the new grid points and make pointers
                    call shd_comp__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                    ! Interpolate the medium properties, radiance, and source
                    call shd_comp__interpPnt(NEWPOINTS, NPTS,&
                        &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                        &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                        &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                        &UNITS, WAVENO, WAVELEN,&
                        &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                        &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                end if
                I = I + 1
            end do

            MAXCELLS = MAXIC - 2
            MAXPTS = MAXIG - 4
            MAXWORK = MAXIDO - 4*NPHI0MAX
            MAXSH = MAXIV - 4*NLM
            MAXRAD = MAXIV + MAXIG - 4*NLM

            ! Split some more cells to make the adaptive grid "smoother"
            do while (.not. OUTOFMEM .and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    &NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    &SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    &RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM = .true.
                else
                    ! See if cell should be split to make better grid:
                    !   IDIR=0 if no, IDIR=direction (1,2,3) otherwise
                    ICELL = int(ADAPTIND(I))/4
                    call shd_grid__testGridSmooth(GRIDPTR, GRIDPOS, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, SHPTR, ICELL, IDIR)
                    ! Set IDIR=0 to stop adaptive grid smoothing
                    ! IDIR = 0
                    if (IDIR .gt. 0) then
                        ! Create the new grid points and make pointers
                        call shd_comp__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                            &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                        ! Interpolate the medium properties, radiance, and source
                        call shd_comp__interpPnt(NEWPOINTS, NPTS,&
                            &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                            &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                            &UNITS, WAVENO, WAVELEN,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                            &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                    end if
                end if
                I = I + 1
            end do

        end do
        if (OUTOFMEM0) then
            OUTOFMEM = .true.
            write (6, '(2A,I8,A,I8)')&
                &'Message [_splitGrid]: Further cell splitting requires larger array allocations: ',&
                &'Ncells = ', NCELLS, ', Npnts = ', NPTS
        end if

        ! If in parallel mode, calculate all the direct beam fluxes at once
        if (SRCTYPE .ne. 'T' .and.&
            &(btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) then
            call shd_ext__mpiMakeDirect(OLDNPTS + 1, NPTS, BCFLAG,&
                &IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
                &SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
                &NX, XGRID, NY, YGRID, DIRFLUX)
        end if
    end if

    ! Find the max adaptive cell criterion after the cell divisions
    SPLITCRIT = 0.0
    do ICELL = 1, NCELLS
        if (TREEPTR(2, ICELL) .eq. 0) then
            call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
                &SHPTR, NSTOKES, SOURCE, ICELL,&
                &ADAPT, ADAPTCRIT(1), IDIR)
            SPLITCRIT = max(SPLITCRIT, ADAPTCRIT(1))
        end if
    end do

    return
end subroutine shd_grid__splitGrid

logical function sweepCellBase(BCFLAG, NXC, NYC, NZ, IOCT, ICELL)

    !   Returns the next base grid cell.  The octant the ray direction
    ! is coming from (IOCT) determines the order of sweeping.  ICELL is 0
    ! for the first call to set things up.  The function returns false if
    ! we are done sweeping the grid.  Do loops over IZ, IY, and IX are
    ! simulated.  For open boundary condition do the endpoint cells
    ! (in X,Y) first, so the boundary grid points will be done with
    ! independent pixel radiative transfer.

    integer BCFLAG, NXC, NYC, NZ, IOCT, ICELL
    integer IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
    save IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ

    sweepCellBase = .true.

    ! If first time here then set up the loops
    if (ICELL .eq. 0) then
        if (btest(IOCT - 1, 0)) then
            DIX = +1
            if (btest(BCFLAG, 0)) then
                SIX = NXC
                EIX = max(NXC - 1, 1)
            else
                SIX = 1
                EIX = NXC
            end if
        else
            DIX = -1
            if (btest(BCFLAG, 0)) then
                SIX = 1
                EIX = min(2, NXC)
            else
                SIX = NXC
                EIX = 1
            end if
        end if
        if (btest(IOCT - 1, 1)) then
            DIY = +1
            if (btest(BCFLAG, 1)) then
                SIY = NYC
                EIY = max(NYC - 1, 1)
            else
                SIY = 1
                EIY = NYC
            end if
        else
            DIY = -1
            if (btest(BCFLAG, 1)) then
                SIY = 1
                EIY = min(2, NYC)
            else
                SIY = NYC
                EIY = 1
            end if
        end if
        if (btest(IOCT - 1, 2)) then
            DIZ = +1
            SIZ = 1
            EIZ = NZ - 1
        else
            DIZ = -1
            SIZ = NZ - 1
            EIZ = 1
        end if
        IX = SIX
        IY = SIY
        IZ = SIZ
    else
        ! Sweep through grid, first in X, in Y, then in Z
        if (IX .eq. EIX) then
            if (IY .eq. EIY) then
                if (IZ .eq. EIZ) then
                    sweepCellBase = .false.
                    return
                else
                    IZ = IZ + DIZ
                end if
                IY = SIY
            else
                IY = mod(IY + DIY + NYC - 1, NYC) + 1
            end if
            IX = SIX
        else
            IX = mod(IX + DIX + NXC - 1, NXC) + 1
        end if
    end if

    ! Compute the base grid cell
    ICELL = IZ + (NZ - 1)*(IY - 1) + (NZ - 1)*NYC*(IX - 1)
    return
end function sweepCellBase

logical function sweepCellNext(BCFLAG, NXC, NYC, NZ, IOCT,&
    &ICELL, TREEPTR, CELLFLAGS)

    !   Returns the next grid cell given the last cell (ICELL).  The
    ! octant the ray direction is coming from (IOCT) determines the
    ! order of sweeping.  ICELL is 0 for the first call to set things up
    ! The function returns false if we are done sweeping the grid.
    ! Within a base cell the tree structure is traced using a stack until
    ! we get to a new end cell.

    integer BCFLAG, NXC, NYC, NZ, IOCT, ICELL, TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    integer IC, IDIR, SP, STACK(50)
    logical DONE, sweepCellBase
    save SP, STACK

    IC = ICELL

    ! Special case for first time here
    if (ICELL .eq. 0) then
        IC = 0
        sweepCellNext = sweepCellBase(BCFLAG, NXC, NYC, NZ, IOCT, IC)
        SP = 0
    end if

    ! The scanning direction is the direction (pos or neg) for each
    !   of X,Y,Z that the ray is going in (indicated by IOCT).

    ! Go through the tree checking the cells:
    ! The stack holds the second child cells, which are the
    !   next cell to go to if at an end cell.
    DONE = .false.
    do while (.not. DONE)
        ! If this is an end cell then we are either at the desired cell
        !   or if this is the old cell then we need to pop the stack
        if (TREEPTR(2, IC) .eq. 0) then
            if (IC .ne. ICELL) then
                DONE = .true.
            else
                ! If no more to pop then we are at the base cell, so go to
                !   next base cell if there is one.
                if (SP .eq. 0) then
                    sweepCellNext=&
                        &sweepCellBase(BCFLAG, NXC, NYC, NZ, IOCT, IC)
                    if (.not. sweepCellNext) return
                else
                    ! Otherwise, pop the stack to get the next cell to go to.
                    IC = STACK(SP)
                    SP = SP - 1
                end if
            end if
        else
            ! Not an end cell, so this cell has children.
            ! Depending on the ordinate octant, decide which child to
            !   do first and put the second child on the stack.
            SP = SP + 1
            if (SP .gt. 50) stop 'Error [sweepCellNext]: Stack exceeded.'
            IDIR = ibits(int(CELLFLAGS(IC)), 2, 2) - 1
            if (btest(IOCT - 1, IDIR)) then
                STACK(SP) = TREEPTR(2, IC) + 1
                IC = TREEPTR(2, IC)
            else
                STACK(SP) = TREEPTR(2, IC)
                IC = TREEPTR(2, IC) + 1
            end if
        end if
    end do

    ICELL = IC
    sweepCellNext = .true.
    return
end function sweepCellNext

subroutine shd_grid__makeSweepOrder(NX, NY, NZ, NPTS, NCELLS,&
    &GRIDPTR, TREEPTR, CELLFLAGS,&
    &BCFLAG, IPFLAG, SWEEPORD, GRIDRAD)

    !   Make the array of grid point sweeping order for each discrete
    ! ordinate octant (SWEEPORD).  Rather than pointing to the grid points
    ! directly, the elements of SWEEPORD point to the cell (bits 3-30) and
    ! the corner of the cell (bits 0-2).  The first index is the grid
    ! point visiting order (1-NPTS), and the second index refers to the
    ! octant, but is not the octant number.  The octant index is translated
    ! to the octant number with IOCTORDER - this complication arises because
    ! of the desire to use just 2 or 4 octants in 1D or 2D.

    integer NX, NY, NZ, NPTS, NCELLS, BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDRAD(NPTS)
    integer NOCT, IOCT, JOCT, NXC, NYC
    integer IPCELL, IPT, IORDER
    integer INDEXCORN, ICORNER, CORNDOG(8, 8), IOCTORDER(8)
    logical sweepCellNext
    data CORNDOG/8, 7, 6, 5, 4, 3, 2, 1, 7, 8, 5, 6, 3, 4, 1, 2, 6, 5, 8, 7, 2, 1, 4, 3,&
        &5, 6, 7, 8, 1, 2, 3, 4, 4, 3, 2, 1, 8, 7, 6, 5, 3, 4, 1, 2, 7, 8, 5, 6,&
        &2, 1, 4, 3, 6, 5, 8, 7, 1, 2, 3, 4, 5, 6, 7, 8/
    data IOCTORDER/1, 5, 2, 6, 3, 7, 4, 8/

    ! IOCT is the octant the ray come from or the discrete ordinate
    !   IOCT = 1 + BITX + 2*BITY + 4*BITZ

    if (btest(IPFLAG, 1) .and. btest(IPFLAG, 0)) then
        NOCT = 2
    else if (btest(IPFLAG, 1)) then
        NOCT = 4
    else
        NOCT = 8
    end if
    NXC = NX
    if (btest(BCFLAG, 0)) NXC = NX + 1
    if (btest(BCFLAG, 2) .and. .not. btest(IPFLAG, 0)) NXC = NX - 1
    NYC = NY
    if (btest(BCFLAG, 1)) NYC = NY + 1
    if (btest(BCFLAG, 3) .and. .not. btest(IPFLAG, 1)) NYC = NY - 1

    ! Loop over the number of octants needed; get the octant number
    do JOCT = 1, NOCT
        IOCT = IOCTORDER(JOCT)

        ! Set each grid point flag to indicate unvisited.
        do IPT = 1, NPTS
            GRIDRAD(IPT) = -1.0
        end do

        ! Sweep through all the grid points
        IORDER = 1
        IPCELL = 0
        INDEXCORN = 8
100     continue
        ! INDEXCORN is a logical ordering index to the 8 corners in a cell;
        !   it is not the actual corner. ICORNER is the actual corner and
        !   depends on the ray direction (IOCT) and INDEXCORN.
        !   Scan the 8 points in a cell according to the ray direction:
        !   Do the ones furthest "upstream" first; loop until we find
        !   a point that is not already done.
        if (INDEXCORN .eq. 8) then
            ! If done with the current cell, then go to next
            INDEXCORN = 1
            if (.not. sweepCellNext(BCFLAG, NXC, NYC, NZ, IOCT,&
                &IPCELL, TREEPTR, CELLFLAGS)) goto 900
        else
            INDEXCORN = INDEXCORN + 1
        end if
        ! Compute the corresponding point
        !   ICORNER = IAND(NOT(IEOR(INDEXCORN-1,IOCT-1)),7) + 1
        ICORNER = CORNDOG(INDEXCORN, IOCT)
        IPT = GRIDPTR(ICORNER, IPCELL)
        if (GRIDRAD(IPT) .ge. 0.0) goto 100

        GRIDRAD(IPT) = 1.0
        ! Fill in the next entry with the grid cell and corner of this point
        SWEEPORD(IORDER, JOCT) = ior(ishft(IPCELL, 3), ICORNER - 1)
        IORDER = IORDER + 1
        goto 100

900     continue

    end do

    return
end subroutine shd_grid__makeSweepOrder

end module shd_grid
