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

module shd_init

    use globals_const, only: RD_, PI_, OPI_
    use globals_cfg,   only: VERBOSE_, PLANET_RADIUS_, &
                             isOpenX_, isOpenYorOpen_, isMpiX_, isMpiY_, &
                             is1DX_, is1DYor1D_

    use shd_util, only: shd_util__computePlanck, &
                        shd_util__solveTriDiag

    use shd_sfc,  only: shd_sfc__brdfSurface

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_init__makeGrid
    public :: shd_init__makeCell

    public :: shd_init__makePrp
    public :: shd_init__makeRad
    public :: shd_init__makeFlux

    public :: shd_init__makeTopRad
    public :: shd_init__makeSfcRadFL
    public :: shd_init__makeSfcRadVL
    public :: shd_init__makeSfcRadU
    public :: shd_init__makeSfcRadV

    public :: shd_init__checkSfcU

contains

subroutine shd_init__makeGrid(&
    &BCFLAG, GRIDTYPE, NPX, NPY, NPZ, NX, NY, NZ,&
    &XSTART, YSTART, DELXP, DELYP, ZLEVELS,&
    &XGRID, YGRID, ZGRID)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Makes the XGRID, YGRID, and ZGRID arrays from the input Z levels
    !   and X and Y spacing from the property file according to the GRIDTYPE.
    !
    ! Input:
    !   BCFLAG  : integer, bit flags to specify the horizontal boundary conditions
    !             0: periodic in X&Y; 1: open in x; 2: open in y; 3: open in X&Y
    !             4: for MPI in X (perodic)
    !             8: for MPI in Y (perodic)
    !            12: for MPI in both (perodic)
    !   GRIDTYPE: character(len=1), grid type
    !             "E": even Z base grid between bottom and top;
    !             "P": Z base grid levels from property file;
    !             "F": Z base grid levels from file: shdom-zgrid.txt
    !   NPX, NPY, NPZ : integer, size of the input property arrays
    !   NX, NY, NZ    : integer, base grid size in X, Y and Z
    !   XSTART, YSTART: real, starting X and Y
    !   DELXP, DELYP  : real, horizontal grid spacing of the property arrays
    !   ZLEVELS       : real array, input Z levels from the property file
    !
    ! Output:
    !   XGRID, YGRID, ZGRID: real array, grid locations for each dimension.
    !                        To facilitate the periodic horizontal boundaries,
    !                        there are NX+1 values in XGRID, NY+1 in YGRID,
    !                        NZ values in ZGRID. Currently the grid locations
    !                        are evenly spaced in X and Y.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=I__), intent(in) :: BCFLAG, NPX, NPY, NPZ, NX, NY, NZ
    character(len=1) , intent(in) :: GRIDTYPE
    real(kind=R__), intent(in)    :: XSTART, YSTART, ZLEVELS(NPZ)
    real(kind=R__), intent(inout) :: DELXP, DELYP
    real(kind=R__), intent(out)   :: XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)

    integer(kind=I__) :: IX, IY, IZ

    ! Make the X, Y, and Z grids, depending on type:
    !   E for evenly spaced, P from property file, F from zgrid.inp
    !   Only Even allowed horizontally
    if (DELXP .le. 0.0) then
        DELXP = 1.0
        if (VERBOSE_ .ge. 4) write(*, '(A)') 'Warning [_makeGrid]: Found DELXP<=0.0, reset DELXP to 1.0.'
    end if

    if (isMpiX_) then ! open in X
        do IX = 1, NX
            if (NPX .gt. 1) then
                XGRID(IX) = XSTART + (IX - 1)*(DELXP*(NPX - 1))/float(NX - 1)
            else
                XGRID(IX) = XSTART + (IX - 1)*(DELXP*NPX)/float(NX - 1)
            end if
        end do
    else ! periodic in X
        do IX = 1, NX + 1
            XGRID(IX) = XSTART + (IX - 1)*(DELXP*NPX)/float(NX)
        end do
    end if

    if (DELYP .le. 0.0) then
        DELYP = 1.0
        if (VERBOSE_ .ge. 4) write(*, '(A)') 'Warning [_makeGrid]: Found DELYP<=0.0, reset DELYP to 1.0.'
    end if

    if (isMpiY_) then ! open in Y
        do IY = 1, NY
            if (NPY .gt. 1) then
                YGRID(IY) = YSTART + (IY - 1)*(DELYP*(NPY - 1))/float(NY - 1)
            else
                YGRID(IY) = YSTART + (IY - 1)*(DELYP*NPY)/float(NY - 1)
            end if
        end do
    else ! periodic in Y
        do IY = 1, NY + 1
            YGRID(IY) = YSTART + (IY - 1)*(DELYP*NPY)/float(NY)
        end do
    end if

    if (GRIDTYPE .eq. 'P') then
        if (NZ .ne. NPZ) then
            write (6, *)&
                &'Error [_makeGrid]: For internal Z grid same as property file grid,'
            write (6, *) '  must have same number of grid points (NZ=NPZ).'
            stop
        end if
        do IZ = 1, NZ
            ZGRID(IZ) = ZLEVELS(IZ)
        end do
    else if (GRIDTYPE .eq. 'E') then
        do IZ = 1, NZ
            ZGRID(IZ) = ZLEVELS(1)&
                + (IZ - 1)*(ZLEVELS(NPZ) - ZLEVELS(1))/float(NZ - 1)
        end do
    else if (GRIDTYPE .eq. 'F') then
        if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_makeGrid]: Reading ZGRID from <shdom-zgrid.txt> ...'
        open (UNIT=3, FILE='shdom-zgrid.txt', STATUS='OLD')
        do IZ = 1, NZ
            read (3, *) ZGRID(IZ)
        end do
        close (3)
    else
        stop 'Error [_makeGrid]: Illegal grid type in Z.'
    end if

    return

end subroutine shd_init__makeGrid

subroutine shd_init__makeCell(&
    &BCFLAG, IPFLAG, NX, NY, NZ, NX1, NY1,&
    &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Make the initial grid cell data structure for the base grid.
    !   For periodic boundary conditions there are NX*NY*(NZ-1) grid cells
    !   and (NX+1)*(NY+1)*NZ grid points. The IX=NX+1 column is the same
    !   as the IX=1 column, because of the periodicity, but the grid points
    !   are repeated.  The last grid cell on the "right" has the last column
    !   of grid points on the "left" and the repeated "virtual" set of grid
    !   points on the "right". For open boundary conditions in X and Y
    !   there are (NX+1)*(NY+1)*(NZ-1) grid cells and NX*NY*NZ grid points.
    !   Open boundary conditions are done with the lowest "left most" and
    !   "right most" X and Y columns being independent pixels with
    !   their own grid cells (zero width), therefore one extra grid cell
    !   is required.   Multiple processor boundary conditions have
    !   (NX-1)*(NY-1)*(NZ-1) grid cells and NX*NY*NZ grid points.
    !   Independent pixel mode has NX*NY*(NZ-1) grid cells and NX*NY*NZ
    !   grid points
    !
    ! Input:
    !   IPFLAG    : Independent pixel flag
    !                0: 3D, 1: IPA in X, 2: IPA in Y, 3: IPA in X&Y.
    !   BCFLAG    : integer, bit flags to specify the horizontal boundary conditions
    !               0: periodic in X&Y; 1: open in x; 2: open in y; 3: open in X&Y
    !               4: for MPI in X (perodic)
    !               8: for MPI in Y (perodic)
    !              12: for MPI in both (perodic)
    !   NX, NY, NZ: integer, base grid size in X, Y and Z
    !   NX1, NY1  : integer, base grid actual size in X, Y
    !   XGRID, YGRID, ZGRID: real array, grid locations for each dimension.
    !                        To facilitate the periodic horizontal boundaries,
    !                        there are NX+1 values in XGRID, NY+1 in YGRID,
    !                        NZ values in ZGRID. Currently the grid locations
    !                        are evenly spaced in X and Y.
    !
    ! Output:
    !   NPTS          : total number of grid points
    !   NCELLS        : total number of cells
    !   GRIDPOS(3,IP) : X,Y,Z coordinate values of grid point locations
    !   GRIDPTR(8,IC) : pointer to 8 grid points belonging to a cell.
    !                   Order of pointers is based on binary digits Z-Y-X:
    !                   IOCT = 1 + BITX + 2*BITY + 4*BITZ, where
    !                   BITc is 0 for negative side and 1 for positive side.
    !                   1: (-X, -Y, -Z), 2: (+X, -Y, -Z), 3: (-X, +Y, -Z), 4: (+X, +Y, -Z)
    !                   5: (-X, -Y, +Z), 6: (+X, -Y, +Z), 7: (-X, +Y, +Z), 8: (+X, +Y, +Z)
    !   NEIGHPTR(6,IC): pointer to 6 neighboring cells of a cell.
    !                   Order of pointers are 1:-X, 2:+X, 3:-Y, 4:+Y, 5:-Z, 6:+Z.
    !                   Is 0 if there is no neighbor cell (boundary).
    !                   Is positive if there is just one cell across face
    !                   of current cell (unique neighbor).  Is negative
    !                   (-index) if there is more than one cell across
    !                   the face of current cell.  In this case the
    !                   cell pointed to is the parent cell that has a
    !                   has same or larger face adjoining current cell.
    !   TREEPTR(2,IC) : pointer to parent and children cells in tree structure.
    !                   First is pointer to parent cell, 0 for a base cell.
    !                   Second is pointer to first of the two children cells.
    !                   This pointer is 0 if there are no children cells.
    !                   The children cells are in successive order, with
    !                   the negative side first. See CELLFLAGS for direction
    !                   of division of split.
    !   CELLFLAGS(IC) : a collection of binary flags for the cells:
    !                   bit  For bit set
    !                   0:  Do independent pixel for X
    !                   1:  Do independent pixel for Y
    !                   2,3: Direction cell is split (0=no split, 1=X, 2=Y, 3=Z)
    ! IP is a grid point index, IC is a cell index.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=I__), intent(in)    :: IPFLAG, BCFLAG, NX, NY, NZ, NX1, NY1
    real(kind=R__)   , intent(in)    :: XGRID(NX1), YGRID(NY1), ZGRID(NZ)

    real(kind=R__)   , intent(out)   :: GRIDPOS(3, *)
    integer(kind=I__), intent(inout) :: NPTS, NCELLS
    integer(kind=I__), intent(out)   :: GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer(kind=2)  , intent(out)   :: CELLFLAGS(*)

    integer(kind=I__) :: IX, IY, IZ, IX0, IX1, IY0, IY1, I, NXC, NYC

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
    if (isOpenX_) NXC = NX + 1
    if (isMpiX_ .and. (.not. is1DX_)) NXC = NX - 1
    NYC = NY
    if (isOpenYorOpen_) NYC = NY + 1
    if (isMpiY_ .and. (.not. is1DYor1D_)) NYC = NY - 1

    ! Then setup the cell structures
    I = 0
    do IX = 1, NXC
        if (is1DX_) then
            IX0 = IX
            IX1 = IX
        else if (isOpenX_) then
            IX0 = max(1, IX - 1)
            IX1 = min(NX, IX)
        else
            IX0 = IX
            IX1 = IX + 1
        end if
        do IY = 1, NYC
            if (is1DYor1D_) then
                IY0 = IY
                IY1 = IY
            else if (isOpenYorOpen_) then
                IY0 = max(1, IY - 1)
                IY1 = min(NY, IY)
            else
                IY0 = IY
                IY1 = IY + 1
            end if
            do IZ = 1, NZ - 1
                I = I + 1
                GRIDPTR(1, I) = IZ +     NZ*(IY0 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(2, I) = IZ +     NZ*(IY0 - 1) + NZ*NY1*(IX1 - 1)
                GRIDPTR(3, I) = IZ +     NZ*(IY1 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(4, I) = IZ +     NZ*(IY1 - 1) + NZ*NY1*(IX1 - 1)
                GRIDPTR(5, I) = IZ + 1 + NZ*(IY0 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(6, I) = IZ + 1 + NZ*(IY0 - 1) + NZ*NY1*(IX1 - 1)
                GRIDPTR(7, I) = IZ + 1 + NZ*(IY1 - 1) + NZ*NY1*(IX0 - 1)
                GRIDPTR(8, I) = IZ + 1 + NZ*(IY1 - 1) + NZ*NY1*(IX1 - 1)
                CELLFLAGS(I) = 0
                ! Initialize the neighbor pointer; to self for IP mode or
                !   if open boundary conditions and at either end
                if (is1DX_) then
                    NEIGHPTR(1, I) = I
                    NEIGHPTR(2, I) = I
                    CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 0)
                else if (isOpenX_) then
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
                else if (isMpiX_) then
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
                if (is1DYor1D_) then
                    NEIGHPTR(3, I) = I
                    NEIGHPTR(4, I) = I
                    CELLFLAGS(I) = ibset(int(CELLFLAGS(I)), 1)
                else if (isOpenYorOpen_) then
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
                else if (isMpiY_) then
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

end subroutine shd_init__makeCell

subroutine shd_init__makePrp(&
    &ML, MM, NSTLEG, NLEG, NPTS, DELTAM,&
    &NUMPHASE, SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX,&
    &EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)

    ! Prepares the grid arrays for the iterative solution process.
    ! If doing Delta-M scaling then the extinction, albedo, and Legend
    ! terms are scaled first; only the 0 to ML LEGEN terms are scaled.
    ! Outputs PLANCK with (1-omega)*B(T) for thermal source, where B(T
    ! the Planck function (meaning depends on UNITS).
    ! TEMP array is unchanged.

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
                ! NUMPHASE=0 is standard property format, which is unpolarized
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
            call shd_util__computePlanck(TEMP(I), UNITS, WAVENO, WAVELEN, BB)
            PLANCK(I) = (1.0 - ALBEDO(I))*BB
        end if
    end do

    return

end subroutine shd_init__makePrp

subroutine shd_init__makeRad(&
    &NSTOKES, NXY, NZ, NSTLEG, NLEG,&
    &RSHPTR, ZGRID,&
    &EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP,&
    &SKYRAD, UNITS, WAVENO, WAVELEN, RADIANCE)

    ! Initializes radiance field by solving plane-parallel two-stream.
    ! Solves the L=1 M=0 SH system by transforming the pentadiagonal
    ! system to tridiagonal and calling solver.
    ! Does the initialization for the NXY columns of the base grid.
    ! Only the I Stokes parameter is initialized, the rest are zeroed.

    integer NSTOKES, NXY, NZ, NSTLEG, NLEG, NUMPHASE, RSHPTR(*)
    integer*2 IPHASE(NZ, NXY)
    real GNDALBEDO, GNDTEMP, SKYRAD, SOLARFLUX, SOLARMU
    real WAVENO(2), WAVELEN
    real ZGRID(NZ)
    real EXTINCT(NZ, NXY), ALBEDO(NZ, NXY), LEGEN(NSTLEG, 0:NLEG, *)
    real TEMP(NZ, NXY)
    real RADIANCE(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer NLAYER, I, IZ, IR, J, K, L, ZINDEX
    logical DELTAM
    real C0, C1
    real EXT0, ALB0, SCAT0, G0, GNDEMIS
    real, allocatable :: OPTDEPTHS(:), ALBEDOS(:), ASYMMETRIES(:)
    real, allocatable :: TEMPS(:), FLUXES(:, :)
    character(len=1) :: SRCT

    NLAYER = NZ - 1

    allocate (OPTDEPTHS(NZ), ALBEDOS(NZ), ASYMMETRIES(NZ))
    allocate (TEMPS(NZ), FLUXES(3, NZ))

    !real, allocatable :: CH(:, :)
    !allocate (CH(NLAYER, NLAYER))
    !call shd_init__chapman(NLAYER, SOLARMU, PLANET_RADIUS_, ZGRID(NZ:1:-1), CH)

    DELTAM = .false.

    SRCT = SRCTYPE
    if (SRCT .eq. 'B') SRCT = 'T'

    ! SH coefficients
    C0 = sqrt(1.0d0/PI_)
    C1 = sqrt(3.0d0/(4.0d0*PI_))

    ! Loop over all the columns in the base grid
    IR = 0
    RSHPTR(1) = 0
    do I = 1, NXY

        ! this is modified to address layer averaging issue
        !╭────────────────────────────────────────────────╮!
        do IZ = 1, NLAYER

            L = NZ - IZ

            ZINDEX = IZ + 1

            EXT0 = EXTINCT(ZINDEX, I)
            ALB0 = ALBEDO(ZINDEX, I)
            SCAT0 = ALB0*EXT0

            OPTDEPTHS(L) = (ZGRID(IZ + 1) - ZGRID(IZ))*EXT0

            if (EXT0 .gt. 0.0) then
                ALBEDOS(L) = ALB0
            else
                ALBEDOS(L) = 0.0
            end if

            if (NUMPHASE .gt. 0) then
                J = IPHASE(ZINDEX, I)
            else
                J = NZ*(I - 1) + IZ
            end if
            G0 = LEGEN(1, 1, J)

            if (SCAT0 .gt. 0.0) then
                ASYMMETRIES(L) = G0
            else
                ASYMMETRIES(L) = 0.0
            end if
            TEMPS(L + 1) = TEMP(IZ, I)
        end do
        !╰────────────────────────────────────────────────╯!

        TEMPS(1) = TEMP(NZ, I)

        GNDEMIS = 1.0 - GNDALBEDO

        ! Call the Eddington flux routine
        call shd_init__makeFlux(NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES,&
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

end subroutine shd_init__makeRad

subroutine shd_init__makeFlux(&
    &NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES,&
    &TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU,&
    &GNDTEMP, GNDEMIS, SKYRAD, UNITS, WAVENO, WAVELEN,&
    &FLUXES)

    ! shd_init__makeFlux (EDDRTF)
    ! Computes the layer interface fluxes for a plane-parallel
    ! atmosphere with either solar or thermal source of radiation using
    ! the Eddington approximation.  The medium is specified by a number
    ! of homogeneous layers. For a thermal source the Planck function is
    ! linear with optical depth, while for a solar source it is exponential.
    ! The temperatures, optical depth, single scattering albedo, and
    ! asymmetry parameter are specified for each layer. The boundary
    ! conditions such as the solar flux, and reflection and/or emission from
    ! ground surface are also specified. Delta Eddington scaling may be
    ! used. The diffuse Eddington fluxes and the solar direct flux at
    ! each level are returned.
    ! The model works by calculating the reflection, transmission, and
    ! source terms for each layer from the input properties. A
    ! tri-diagonal matrix solver is then used to compute the diffuse fluxes
    ! at each layer interface from the applied boundary conditions.
    !
    ! Parameters:
    !   Input:
    ! NLAYER         integer      Number of homogenous layers
    !                             (layers are specified from the top down)
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

    integer NLAYER
    logical DELTAM
    real TEMPS(*)
    real OPTDEPTHS(*), ALBEDOS(*), ASYMMETRIES(*)
    real GNDTEMP, GNDEMIS, SKYRAD, SOLARFLUX, SOLARMU
    real WAVENO(2), WAVELEN
    character(len=1) :: SRCTYPE, UNITS
    real FLUXES(3, *)

    integer N, L, I
    real(RD_) :: DELTAU, G, OMEGA, F
    real(RD_) :: LAMBDA, R, T, D, CP, CM, A, B, X1, X2
    real(RD_) :: REFLECT, TRANS, SOURCEP, SOURCEM
    real(RD_) :: RADP1P, RADP1M, RADP2P, RADP2M
    real(RD_) :: MU0, SKYFLUX, GNDFLUX, PLANCK1, PLANCK2, C, TAU
    real(RD_) :: EXLP, EXLM, V, DS, B1, B2, SOLPP, SOLPM
    real(RD_),  allocatable :: LOWER(:), UPPER(:), DIAG(:), RHS(:)
    real BBRAD

    N = 2*NLAYER + 2
    allocate (LOWER(N), UPPER(N), DIAG(N), RHS(N))

    ! Compute the reflection, transmission, and source
    !   coefficients for each layer for the diffuse Eddington
    !   two stream problem.
    if (SRCTYPE .eq. 'T') then
        call shd_util__computePlanck(TEMPS(1), UNITS, WAVENO, WAVELEN, BBRAD)
        PLANCK1 = PI_*BBRAD
    end if
    MU0 = abs(SOLARMU)
    TAU = 0.0
    I = 2
    do L = 1, NLAYER
        DELTAU = OPTDEPTHS(L)
        if (DELTAU .lt. 0.0) stop 'Error [_makeFlux]: TAU<0.'
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
                call shd_util__computePlanck(TEMPS(L + 1), UNITS, WAVENO, WAVELEN, BBRAD)
                PLANCK2 = PI_*BBRAD
                V = 2.0*(PLANCK2 - PLANCK1)/(3.0*(1.0-OMEGA*G)*DELTAU)
                RADP1P = -V + PLANCK1
                RADP2M =  V + PLANCK2
                RADP2P = -V + PLANCK2
                RADP1M =  V + PLANCK1
                if (LAMBDA .eq. 0.0) then
                    A = (R*DELTAU*RADP1P - RADP2M)*D
                    B = -(R*RADP1P + T*RADP2M)*D
                    SOURCEP = (B - T*(A + B*DELTAU))/R + RADP2P
                    SOURCEM = A + RADP1M
                else
                    CP = ( X1*EXLM*RADP1P - X2*RADP2M)*D
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
                SOLPP =  1.0 + 1.5*G*MU0
                SOLPM = -1.0 + 1.5*G*MU0
                RADP1P = (( T + 1.0/MU0)*SOLPP + R*SOLPM)*B1
                RADP2M = ((-T + 1.0/MU0)*SOLPM - R*SOLPP)*B2
                RADP2P = (( T + 1.0/MU0)*SOLPP + R*SOLPM)*B2
                RADP1M = ((-T + 1.0/MU0)*SOLPM - R*SOLPP)*B1
                if (LAMBDA .eq. 0.0) then
                    A =  (R*DELTAU*RADP1P -   RADP2M)*D
                    B = -(R*RADP1P        + T*RADP2M)*D
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
        SKYFLUX = PI_*SKYRAD
    else
        FLUXES(3, NLAYER + 1) = 0.0
        call shd_util__computePlanck(GNDTEMP, UNITS, WAVENO, WAVELEN, BBRAD)
        GNDFLUX = PI_*BBRAD*GNDEMIS
        call shd_util__computePlanck(SKYRAD, UNITS, WAVENO, WAVELEN, BBRAD)
        SKYFLUX = PI_*BBRAD
    end if
    ! Setup for and call the tri-diagonal matrix solver
    RHS(1) = SKYFLUX
    DIAG(1) = 0.0
    UPPER(1) = 1.0
    DIAG(N) = -(1.0 - GNDEMIS)
    LOWER(N) = 1.0
    RHS(N) = GNDFLUX
    call shd_util__solveTriDiag(N, LOWER, DIAG, UPPER, RHS)

    ! Put the fluxes in the output array
    if (UNITS .eq. 'T') then
        C = 1.0/PI_
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

end subroutine shd_init__makeFlux

subroutine shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
    &UNITS, NTOPPTS, NSTOKES, BCRAD)

    !   Returns the beginning of the BCRAD array with the downwelling
    ! radiance for the NTOPPTS top boundary points.  Currently all
    ! points have the same unpolarized isotropic radiance.

    integer NTOPPTS, NSTOKES
    real SKYRAD, WAVENO(2), WAVELEN
    real BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer IBC, K
    real SKYRAD2

    ! At top, boundary radiance is any isotropic sky radiation
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        SKYRAD2 = SKYRAD
    else
        call shd_util__computePlanck(SKYRAD, UNITS, WAVENO, WAVELEN, SKYRAD2)
    end if
    ! Loop over all points assigning the uniform radiance
    do IBC = 1, NTOPPTS
        BCRAD(1, IBC) = SKYRAD2
        BCRAD(2:, IBC) = 0.0
    end do

    return

end subroutine shd_init__makeTopRad

subroutine shd_init__makeSfcRadFL(NBOTPTS, BCPTR,&
    &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
    &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD)

    !   Returns the BCRAD array with the upwelling unpolarized radiance for
    ! a Lambertian surface for the bottom points (pointed to by BCPTR).
    ! If there is reflection the downwelling hemispheric flux at the
    ! bottom points is gotten (from FLUXES). There may also be thermal
    ! emission (depending on temperature and emissivity) or reflected
    ! collimated solar radiation (DIRFLUX).

    integer NBOTPTS, BCPTR(*), NSTOKES
    real DIRFLUX(*), FLUXES(2, *)
    real GNDTEMP, GNDALBEDO, WAVENO(2), WAVELEN
    real BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IBC, K
    real ALB, GNDRAD

    ALB = GNDALBEDO/acos(-1.0)
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        call shd_util__computePlanck(GNDTEMP, UNITS, WAVENO, WAVELEN, GNDRAD)
        GNDRAD = GNDRAD*(1.0 - GNDALBEDO)
    end if

    ! Loop over all bottom points finding the upwelling radiance
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        if (SRCTYPE .eq. 'S') then
            BCRAD(1, IBC) = ALB*(DIRFLUX(I) + FLUXES(1, I))
        else if (SRCTYPE .eq. 'T') then
            BCRAD(1, IBC) = GNDRAD + ALB*FLUXES(1, I)
        else if (SRCTYPE .eq. 'B') then
            BCRAD(1, IBC) = ALB*(DIRFLUX(I) + FLUXES(1, I)) + GNDRAD
        end if
        BCRAD(2:, IBC) = 0.0
    end do

    return

end subroutine shd_init__makeSfcRadFL

subroutine shd_init__makeSfcRadVL(NBOTPTS, BCPTR,&
    &DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &NSTOKES, BCRAD)

    !   Returns the BCRAD array with the upwelling unpolarized radiance for a
    ! variable Lambertian surface for the bottom points (pointed to by BCPTR).
    ! SFCGRIDPARMS contains the interpolated surface properties for the
    ! bottom boundary points; the first is Planck function and the second is
    ! albedo.  The reflection is computed from the downwelling hemispheric
    ! flux (in FLUXES).  There may also be thermal emission (depending on
    ! temperature and emissivity) or reflected collimated solar radiation
    ! (from DIRFLUX).

    integer NBOTPTS, BCPTR(*), NSFCPAR, NSTOKES
    real DIRFLUX(*), FLUXES(2, *)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS)
    real BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE
    integer I, IBC, K
    real ALB, GNDRAD

    ! Loop over all bottom points
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        ALB = SFCGRIDPARMS(2, IBC)
        if (SRCTYPE .eq. 'S') then
            BCRAD(1, IBC) = OPI_*ALB*(DIRFLUX(I) + FLUXES(1, I))
        else
            GNDRAD = SFCGRIDPARMS(1, IBC)*(1 - ALB)
            if (SRCTYPE .eq. 'T') then
                BCRAD(1, IBC) = GNDRAD + OPI_*ALB*FLUXES(1, I)
            else if (SRCTYPE .eq. 'B') then
                BCRAD(1, IBC) = OPI_*ALB*(DIRFLUX(I) + FLUXES(1, I)) + GNDRAD
            end if
        end if
        BCRAD(2:, IBC) = 0.0
    end do
    return

end subroutine shd_init__makeSfcRadVL

subroutine shd_init__makeSfcRadU(NBOTPTS, BCPTR,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, IMU, IPHI,&
    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS, SFC_BRDF_DO,&
    &NSTOKES, BCRAD)

    !   Computes the upwelling Stokes vector of the bottom boundary points
    ! for one outgoing discrete ordinate using the precomputed bidirectional
    ! reflectance distribution matrix for the discrete ordinates.
    ! The upwelling Stokes vector includes the reflection of the incident
    ! radiance, the thermal emission (emissivity*Planck function), and the
    ! reflected direct solar flux (if applicable). The upwelling radiance
    ! vector is the integral over all incident directions of the BRDF times
    ! the downwelling radiance, so a discrete sum is done using the SFC_BRDF_DO
    ! matrix, which includes the integral weights (WTDO). The emissivity
    ! is computed implicitly from the integral of the BRDF. The outgoing
    ! direction is specified with (IMU,IPHI), and the BRDF is computed for
    ! all incident directions (loop over JMU,JPHI). The incident downwelling
    ! radiance vectors are input in BCRAD(*,*,2...) and the outgoing
    ! upwelling radiance vectors are output in BCRAD(*,*,1).

    integer NBOTPTS, BCPTR(NBOTPTS)
    integer NMU, NPHI0MAX, NPHI0(NMU), IMU, IPHI, NSFCPAR, NSTOKES
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    real BCRAD(NSTOKES, NBOTPTS, *)
    character(len=1) :: SRCTYPE
    character(len=2) :: SFCTYPE
    integer IMU1, JMU, JPHI, IBC, I, JANG, K, K1, K2
    real REFLECT(4, 4), W

    IMU1 = IMU - NMU/2
    do IBC = 1, NBOTPTS

        ! Initialize the upwelling boundary radiances to zero or to
        !   the reflection of direct unpolarized solar flux.
        BCRAD(:, IBC, 1) = 0.0
        if (SRCTYPE .ne. 'T') then
            I = BCPTR(IBC)
            call shd_sfc__brdfSurface(SFCTYPE(2:2), SFCGRIDPARMS(2, IBC), WAVELEN,&
                &MU(IMU), PHI(IMU, IPHI), SOLARMU, SOLARAZ,&
                &NSTOKES, REFLECT)
            BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                &+ OPI_*REFLECT(1:NSTOKES, 1)*DIRFLUX(I)
        end if

        ! Integrate over the incident discrete ordinate directions (JMU, JPHI)
        JANG = 1
        do JMU = 1, NMU/2
            do JPHI = 1, NPHI0(JMU)
                ! Add in the polarized reflection
                do K1 = 1, NSTOKES
                    BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1) + BCRAD(K1, IBC, JANG + 1)&
                        &*SFC_BRDF_DO(:, IMU1, IPHI, K1, JMU, JPHI)
                end do
                ! Add in the polarized thermal emission
                W = OPI_*abs(MU(JMU))*WTDO(JMU, JPHI)
                BCRAD(1, IBC, 1) = BCRAD(1, IBC, 1) + SFCGRIDPARMS(1, IBC)&
                    &*(W - SFC_BRDF_DO(1, IMU1, IPHI, 1, JMU, JPHI))
                BCRAD(2:, IBC, 1) = BCRAD(2:, IBC, 1) - SFCGRIDPARMS(1, IBC)&
                    &*SFC_BRDF_DO(2:, IMU1, IPHI, 1, JMU, JPHI)
                JANG = JANG + 1
            end do
        end do
    end do

    return

end subroutine shd_init__makeSfcRadU

subroutine shd_init__makeSfcRadV(&
    &NBOTPTS, IBEG, IEND, BCPTR,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,&
    &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES, BCRAD)

    !   Computes the upwelling Stokes vector of the bottom boundary points
    ! for one outgoing direction using the specified bidirectional
    ! reflectance distribution function.  The upwelling Stokes vector
    ! includes the reflection of the incident radiance, the thermal
    ! emission (emissivity*Planck function), and the reflected direct
    ! solar flux (if applicable).  The upwelling radiance vector is the
    ! integral over all incident directions of the BRDF times the
    ! downwelling radiance, so a discrete sum is done and the integral
    ! weights (WTDO) are included. The general BRDF function is called
    ! to compute the Stokes reflectance matrix for the particular BRDF
    ! type (SFCTYPE) with parameters (SFCGRIDPARMS) for the incident and
    ! outgoing directions.  The emissivity is computed implicitly from the
    ! integral of the BRDF.  The outgoing direction is specified with
    ! (MU2,PHI2), and the BRDF is computed for all incident directions
    ! (loop over JMU,JPHI).  The incident downwelling radiance vectors are
    ! input in BCRAD(*,*,2...) and the outgoing upwelling radiance vectors
    ! are output in BCRAD(*,*,1).

    integer NBOTPTS, IBEG, IEND, BCPTR(NBOTPTS)
    integer NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR, NSTOKES
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real MU2, PHI2, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS), BCRAD(NSTOKES, NBOTPTS, *)
    character(len=1) :: SRCTYPE
    character(len=2) :: SFCTYPE
    integer JMU, JPHI, IBC, I, JANG, K, K1, K2
    real REFLECT(4, 4), W, SUM0, SUM1, RADSPEC(4)

    do IBC = IBEG, IEND

        ! Initialize the upwelling boundary radiances to zero or to
        !   the reflection of direct unpolarized solar flux.
        BCRAD(:, IBC, 1) = 0.0
        if (SRCTYPE .ne. 'T') then
            I = BCPTR(IBC)
            call shd_sfc__brdfSurface(SFCTYPE(2:2), SFCGRIDPARMS(2, IBC), WAVELEN,&
                &MU2, PHI2, SOLARMU, SOLARAZ, NSTOKES, REFLECT)
            BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                &+ OPI_*REFLECT(1:NSTOKES, 1)*DIRFLUX(I)
        end if

        ! Integrate over the incident discrete ordinate directions (JMU,JPHI)
        JANG = 1
        do JMU = 1, NMU/2
            do JPHI = 1, NPHI0(JMU)
                call shd_sfc__brdfSurface(SFCTYPE(2:2), SFCGRIDPARMS(2, IBC),&
                    &WAVELEN, MU2, PHI2, MU(JMU), PHI(JMU, JPHI),&
                    &NSTOKES, REFLECT)
                W = OPI_*abs(MU(JMU))*WTDO(JMU, JPHI)
                ! Add in the polarized reflection
                do K1 = 1, NSTOKES
                    BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                        &+ W*REFLECT(1:NSTOKES, K1)*BCRAD(K1, IBC, JANG + 1)
                end do
                ! Add in the polarized thermal emission
                BCRAD(1, IBC, 1) = BCRAD(1, IBC, 1)&
                    &+ W*(1 - REFLECT(1, 1))*SFCGRIDPARMS(1, IBC)
                BCRAD(2:, IBC, 1) = BCRAD(2:, IBC, 1)&
                    &- W*REFLECT(2:NSTOKES, 1)*SFCGRIDPARMS(1, IBC)
                JANG = JANG + 1
            end do
        end do

    end do
    return
end subroutine shd_init__makeSfcRadV

subroutine shd_init__checkSfcU(&
    &NXSFC, NYSFC, NSFCPAR, SFCPARMS,&
    &UNIFORM_SFC_BRDF)

    ! Description:
    !   Checks to see if the surface reflection parameters are uniform,
    !   and if so sets the UNIFORM_SFC_BRDF flag.  The temperature need
    !   not be uniform.
    !
    ! Input:
    !   NXSFC, NYSFC: integer, size of surface arrays (internally NXSFC+1 by NYSFC+1)
    !   NSFCPAR: integer, number of surface parameters for modeling different BRDF
    !   SFCPARMS: real array, surface parameters for modeling different BRDF
    !
    ! Output:
    !   UNIFORM_SFC_BRDF: logical, whether surface reflection parameters are uniform

    integer(kind=I__), intent(in) :: NXSFC, NYSFC, NSFCPAR
    real(kind=R__)   , intent(in) :: SFCPARMS(NSFCPAR, NXSFC + 1, NYSFC + 1)

    logical, intent(out) :: UNIFORM_SFC_BRDF

    integer(kind=I__) :: IX, IY, J

    UNIFORM_SFC_BRDF = .true.

    do IY = 1, NYSFC
        do IX = 1, NXSFC
            do J = 2, NSFCPAR
                if (SFCPARMS(J, IX, IY) .ne. SFCPARMS(J, 1, 1)) then
                    UNIFORM_SFC_BRDF = .false.
                    return
                end if
            end do
        end do
    end do

    return

end subroutine shd_init__checkSfcU

subroutine shd_init__chapman(&
    &NLAYER, SOLARMU, PLANET_RADIUS, ZLEVELS,&
    &CH)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Apply Chapman function for psuedo spherical correction
    !
    ! Input:
    !   NLAYER : integer, number of layers
    !   SOLARMU: real, cosine of incident angle
    !   PLANET_RADIUS: real, planet radius in km
    !   ZLEVELS: real array, layer height, go from top to surface
    !
    ! Output:
    !   CH: Chapman function (see Eq. 9 in Dahlback and Stamnes, 1990)
    !
    ! Notes:
    !   This subroutine is adopted from DISORT 4.0.99 (http://www.rtatmocn.com/disort/)
    !   subroutine <chapman> in <disort.f>
    !   Reference: Dahlback and Stamnes, 1990 - A NEW SPHERICAL MODEL FOR COMPUTING
    !              THE RADIATION FIELD AVAILABLE FOR PHOTOLYSIS AND HEATING AT TWILIGHT
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=I__), intent(in) :: NLAYER
    real(kind=R__)   , intent(in) :: PLANET_RADIUS, ZLEVELS(0:NLAYER), SOLARMU

    real(kind=R__), intent(out) :: CH(NLAYER, NLAYER)

    integer(kind=I__) :: lc, n
    real(kind=R__)    :: cos_theta0, sin_theta0, cos_theta, sin_theta
    real(kind=R__)    :: s1, s2
    real(kind=R__)    :: h(0:NLAYER)

    if (h(0) .lt. h(NLAYER)) stop 'Error [_chapman]: Input ZLEVELS should go from top to bottom.'

    cos_theta0 = SOLARMU
    sin_theta0 = sqrt(1.0 - SOLARMU*SOLARMU)
    h = ZLEVELS + PLANET_RADIUS

    CH(:, :) = -1.0

    do lc = 1, NLAYER

        sin_theta = sin_theta0 * h(lc) / h(0)
        cos_theta = sqrt(1.0 - sin_theta*sin_theta)
        s1 = 0.0

        do n = 1, lc

            s2 = h(0)*cos_theta - sqrt(h(n)*h(n) - h(0)*h(0)*sin_theta*sin_theta)

            !consider special case for delta h = 0
            if ((s1 .ne. s2) .and. (h(n-1) .ne. h(n))) then
                CH(lc,n) = (h(n-1) - h(n)) / (s2 - s1)
            else
                if (n .gt. 1) then
                    CH(lc,n) = CH(lc, n-1)
                else
                    CH(lc,n) = SOLARMU
                end if
            endif
            s1 = s2

        end do

    end do

    do n = 1, NLAYER
        print *, n
        print *, CH(n, :)
        print *
    end do

    return

end subroutine shd_init__chapman


end module shd_init
