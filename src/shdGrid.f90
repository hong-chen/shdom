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

module shd_grid

    use globals_const, only: RD_
    use globals_cfg,   only: accSolv_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, isMpi_, &
                             is1DX_, is1DYx_, is1D_, &
                             isPol_, isSolar_, isThermal_
    use globals_type,  only: shd_type_rtm, &
                             shd_type_grid

    use globals_obj, only: shdo_obj_, iradI_Shdo, radI_Shdo, isrcJ_Shdo, srcJ_Shdo, &
                            src_obj_, &
                            prp_obj_, temp_Prp, ext_Prp, ssa_Prp, leg_Prp, ipha_Prp, &
                            sfc_obj_, &
                            ckd_obj_, kabs_Ckd

    use shd_util, only: shd_util__findNextCell, &
                        shd_util__computePlanck, &
                        shd_util__ssort

    use shd_rad,  only: shd_rad__computeDirect

    use shd_src,  only: shd_src__computeSrcOnePnt_, &
                        shd_src__computeSrcOnePnt

    use shd_mpi,  only: shd_mpi__makeDirect

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    integer(kind=I__) :: gridCorner(2, 4, 3)
    data gridCorner /1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        &1, 5, 2, 6, 3, 7, 4, 8/

    integer(kind=I__) :: oppFace(6)
    data oppFace /2, 1, 4, 3, 6, 5/

    public :: shd_grid__findBndryPnt
    public :: shd_grid__matchGridPnt
    public :: shd_grid__makeNewGridPnt

    public :: shd_grid__makeSweepOrder
    public :: shd_grid__divideCell
    public :: shd_grid__matchNeighCellFace
    public :: shd_grid__inheritNeighCell

    public :: shd_grid__testGridSmooth
    public :: shd_grid__testCellSplit

    public :: shd_grid__interpGrid
    public :: shd_grid__interpNewGridPnt

    public :: shd_grid__splitGrid

    private :: shd_grid__InterpTriLinGridPrp
    private :: sweepCellNext
    private :: sweepCellBase

contains

subroutine shd_grid__findBndryPnt(&
    &NANG, MAXNBC, MAXBCRAD, &
    &ZBOT, ZTOP,&
    &NTOPPTS, NBOTPTS, BCPTR,&
    &grid_obj)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Returns the BCPTR pointer array for the NTOPPTS top boundary
    !   points (BCPTR(*,1)) and the NBOTPTS bottom boundary points (BCPTR(*,2)).
    !   All the grid points are looked at to find the boundary ones.
    !
    ! Input:
    !   grid_obj%Npts : integer, total number of grid points
    !   NANG      : integer, number of discrete ordinate angles (NMU and NPHI)
    !   MAXNBC    : integer, maximum number of boundary grid points
    !   MAXBCRAD  : integer, maximum number of radiances at boundary grid points
    !   sfc_obj_%isLambertian: logical, lambertian surface or not
    !   ZBOT      : real, altitude of bottom boundary
    !   ZTOP      : real, altitude of bottom boundary
    !   grid_obj%xyz(3, grid_obj%Npts): real array, (x,y,z) coordinate values of grid point locations
    !
    ! Output:
    !   NTOPPTS   : integer, number of top boundary points
    !   NBOTPTS   : integer, number of bottom boundary points
    !   BCPTR(MAXNBC, 2): integer array, pointer array for NTOPPTS (top boundary points)
    !               and NBOTPTS (bottom boundary points)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    type(shd_type_grid), intent(in) :: grid_obj

    integer(kind=I__), intent(in) :: NANG, MAXNBC, MAXBCRAD
    real(kind=R__)   , intent(in) :: ZBOT, ZTOP

    integer(kind=I__), intent(out) :: NTOPPTS, NBOTPTS
    integer(kind=I__), intent(out) :: BCPTR(MAXNBC, 2)

    integer(kind=I__) :: I, IT, IB, NA

    if (sfc_obj_%isLambertian) then
        NA = 1
    else
        NA = NANG/2 + 1
    end if

    IT = 0
    IB = 0
    ! Loop over all points finding the top and bottom ones
    do I = 1, grid_obj%Npts

        if (grid_obj%xyz(3, I) .ge. ZTOP) then
            IT = IT + 1
            if (IT .gt. MAXNBC .or. IT .gt. MAXBCRAD) then
                stop 'Error [_findBndryPnt]: MAXNBC exceeded - top.'
            end if
            BCPTR(IT, 1) = I
        end if

        if (grid_obj%xyz(3, I) .le. ZBOT) then
            IB = IB + 1
            if ((IB .gt. MAXNBC) .or. (IT + IB*NA .gt. MAXBCRAD)) then
                stop 'Error [_findBndryPnt]: MAXNBC exceeded - bottom.'
            end if
            BCPTR(IB, 2) = I
        end if

    end do

    NTOPPTS = IT
    NBOTPTS = IB

    return

end subroutine shd_grid__findBndryPnt

subroutine shd_grid__matchGridPnt(&
    &XP, YP, ZP, icell, iface,&
    &iptsMatch,&
    &grid_obj)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Determines whether the neighboring cells to icell adjacent to
    !   face iface (1-6: -X,+X,-Y,+Y,-Z,+Z) have grid points that match
    !   the point at XP,YP,ZP on face iface of cell icell.
    !   If the face of icell has several neighboring cells, then the tree
    !   is traced down to find those cells that might contain XP,YP,ZP.
    !   iptsMatch is returned 0 if there is no match, otherwise it has the
    !   number of the matching existing gridpoint.
    !
    ! Input:
    !   XP        : real, total number of grid points
    !   grid_obj%xyz : real array, (x,y,z) coordinate values of grid point locations
    !
    ! Output:
    !   iptsMatch : integer, pointer array for NTOPPTS (top boundary points)
    !               and NBOTPTS (bottom boundary points)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    type(shd_type_grid), intent(in) :: grid_obj

    real(kind=R__)   , intent(in)  :: XP, YP, ZP
    integer(kind=I__), intent(in)  :: icell, iface
    integer(kind=I__), intent(out) :: iptsMatch

    integer(kind=2), parameter :: MAXSTACK=64
    integer(kind=I__) :: cellStack(MAXSTACK), SP

    integer(kind=I__) :: idir, DIR, IC, IC1, KFACE, ipts, I

    integer(kind=I__) :: gridFace(4, 6)
    data gridFace /1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8/

    iptsMatch = 0
    idir = (iface + 1)/2
    KFACE = oppFace(iface)

    ! Get the pointer to the neighbor cell
    IC = abs(grid_obj%i6Cell(iface, icell))
    ! Don't process if the neighbor cell is a boundary
    if (IC .eq. 0) return

    ! Go down all relevant branches of the tree until there are no children
    SP = 0
    do while (.true.)

        ! If at an end cell then check for a gridpoint match at the four
        !   corners of the cell's face
        do while (grid_obj%i2Tree(2, IC) .eq. 0)
            do I = 1, 4
                ipts = grid_obj%i8Grid(gridFace(I, KFACE), IC)
                if (XP .eq. grid_obj%xyz(1, ipts) .and. YP .eq. grid_obj%xyz(2, ipts)&
                    &.and. ZP .eq. grid_obj%xyz(3, ipts)) then
                    iptsMatch = ipts
                    return
                end if
            end do
            ! Return if there are no more branches on the stack, otherwise pop
            if (SP .eq. 0) return
            IC = cellStack(SP)
            SP = SP - 1
        end do

        ! Get the direction current cell is split
        DIR = ibits(int(grid_obj%flag(IC)), 2, 2)

        if (DIR .eq. 0) stop 'Error [_matchGridPnt]: No split direction.'

        ! If split is in direction of face then choose the near child cell
        !   which is the negative side cell if this is a positive face.
        IC1 = grid_obj%i2Tree(2, IC)
        if (DIR .eq. idir) then
            IC = IC1 + 1 - mod(iface - 1, 2)
        else
            !   If split is not in direction of face then find which of two cells
            ! has a face that contains the current location by comparing the
            ! appropriate coordinate (X,Y,Z) of the grid point on the split line
            ! with the current location.  If the split line matches a coordinate
            ! then we have to trace down both cells so put one on the stack.
            IC = IC1
            if (DIR .eq. 1) then
                if (XP .eq. grid_obj%xyz(1, grid_obj%i8Grid(8, IC1))) then
                    SP = SP + 1
                    cellStack(SP) = IC + 1
                else if (XP .gt. grid_obj%xyz(1, grid_obj%i8Grid(8, IC1))) then
                    IC = IC + 1
                end if
            else if (DIR .eq. 2) then
                if (YP .eq. grid_obj%xyz(2, grid_obj%i8Grid(8, IC1))) then
                    SP = SP + 1
                    cellStack(SP) = IC + 1
                else if (YP .gt. grid_obj%xyz(2, grid_obj%i8Grid(8, IC1))) then
                    IC = IC + 1
                end if
            else
                if (ZP .eq. grid_obj%xyz(3, grid_obj%i8Grid(8, IC1))) then
                    SP = SP + 1
                    cellStack(SP) = IC + 1
                else if (ZP .gt. grid_obj%xyz(3, grid_obj%i8Grid(8, IC1))) then
                    IC = IC + 1
                end if
            end if
            if (SP .ge. MAXSTACK)&
                &stop 'Error [_matchGridPnt]: MAXSTACK exceeded.'
        end if
    end do

    return

end subroutine shd_grid__matchGridPnt

subroutine shd_grid__makeNewGridPnt(idir, icell, newCell,&
    &newPts,&
    &grid_obj)

    !   Sets up the pointers to the grid points for the two new cells.
    ! Creates the four new points if they don't already exist, changing
    ! grid_obj%Npts and grid_obj%xyz accordingly.  Returns the newPts(3,4) arrays
    ! that gives information about the new points: newPts(3,) is the
    ! point address of the new points (or zero), and newPts(1,) and
    ! newPts(2,) are the "parent" grid points these should be
    ! interpolated from.

    type(shd_type_grid), intent(inout) :: grid_obj

    integer idir, icell, newCell, newPts(3, 4)
    integer I, K, I1, I2, IP1, IP2, ipts, iptsMatch
    integer iface1, iface2, ICELL2
    real XP, YP, ZP

    integer(kind=I__) :: FACEGRID(2, 4, 3)
    data FACEGRID /3, 5, 4, 5, 3, 6, 4, 6, 1, 5, 2, 5, 1, 6, 2, 6,&
        &1, 3, 2, 3, 1, 4, 2, 4/

    !   In general, when subdividing a cell four new points may be made.
    ! Each of these points belongs to two faces of the original cell.
    ! Use the shd_grid__matchGridPnt routine to try to find matches in cells
    ! across the two faces that contain the point; also cross the face of
    ! one of these to get the third bounding cell.  If any of the appropriate
    ! gridpoints of these three cells has the same coordinates then the point
    ! already exists, otherwise create the point.

    ! Temporarily set new cells' grid pointer to the parent cell's ones
    do K = 1, 8
        grid_obj%i8Grid(K, newCell) = grid_obj%i8Grid(K, icell)
        grid_obj%i8Grid(K, newCell + 1) = grid_obj%i8Grid(K, icell)
    end do

    ! Loop over the four edges with potentially new points
    do I = 1, 4
        ! Get the grid corners and coordinates for this edge and split direction
        I1 = gridCorner(1, I, idir)
        I2 = gridCorner(2, I, idir)
        IP1 = grid_obj%i8Grid(I1, icell)
        IP2 = grid_obj%i8Grid(I2, icell)
        XP = (grid_obj%xyz(1, IP1) + grid_obj%xyz(1, IP2))/2
        YP = (grid_obj%xyz(2, IP1) + grid_obj%xyz(2, IP2))/2
        ZP = (grid_obj%xyz(3, IP1) + grid_obj%xyz(3, IP2))/2

        ! Look for gridpoints matching this new point in the three cell
        !   that border this edge
        iface1 = FACEGRID(1, I, idir)
        call shd_grid__matchGridPnt(XP, YP, ZP, icell, iface1,&
            &iptsMatch, grid_obj)
        if (iptsMatch .eq. 0) then
            iface2 = FACEGRID(2, I, idir)
            call shd_grid__matchGridPnt(XP, YP, ZP, icell, iface2,&
                &iptsMatch, grid_obj)
        end if
        if (iptsMatch .eq. 0) then
            ICELL2 = abs(grid_obj%i6Cell(iface1, icell))
            if (ICELL2 .gt. 0) then
                call shd_grid__matchGridPnt(XP, YP, ZP, ICELL2, iface2,&
                    &iptsMatch, grid_obj)
            end if
        end if

        ! Set some new cells' grid pointers to the new or matching one
        if (iptsMatch .eq. 0) then
            grid_obj%Npts = grid_obj%Npts + 1
            grid_obj%i8Grid(I2, newCell) = grid_obj%Npts
            grid_obj%i8Grid(I1, newCell + 1) = grid_obj%Npts
            grid_obj%xyz(1, grid_obj%Npts) = XP
            grid_obj%xyz(2, grid_obj%Npts) = YP
            grid_obj%xyz(3, grid_obj%Npts) = ZP
            newPts(1, I) = IP1
            newPts(2, I) = IP2
            newPts(3, I) = grid_obj%Npts
        else
            grid_obj%i8Grid(I2, newCell) = iptsMatch
            grid_obj%i8Grid(I1, newCell + 1) = iptsMatch
            newPts(3, I) = 0
        end if
    end do

    return

end subroutine shd_grid__makeNewGridPnt

subroutine shd_grid__makeSweepOrder(&
    &SWEEPORD, GRIDRAD,&
    &rtm_obj, grid_obj)

    ! Description:
    !   Make the array of grid point sweeping order for each discrete
    !   ordinate octant (SWEEPORD). Rather than pointing to the grid points
    !   directly, the elements of SWEEPORD point to the cell (bits 3-30) and
    !   the corner of the cell (bits 0-2). The first index is the grid
    !   point visiting order (1-grid_obj%Npts), and the second index refers to the
    !   octant, but is not the octant number. The octant index is translated
    !   to the octant number with IOCTORDER - this complication arises because
    !   of the desire to use just 2 or 4 octants in 1D or 2D.
    !
    ! Input:
    !

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer(kind=I__), intent(out) :: SWEEPORD(grid_obj%Npts, *)
    real(kind=R__)   , intent(out) :: GRIDRAD(grid_obj%Npts)

    integer(kind=I__) :: NOCT, IOCT, JOCT, NXC, NYC
    integer(kind=I__) :: IPCELL, ipts, IORDER
    integer(kind=I__) :: INDEXCORN, ICORNER

    integer(kind=I__) :: CORNDOG(8, 8)
    data CORNDOG /&
        &8, 7, 6, 5, 4, 3, 2, 1,&
        &7, 8, 5, 6, 3, 4, 1, 2,&
        &6, 5, 8, 7, 2, 1, 4, 3,&
        &5, 6, 7, 8, 1, 2, 3, 4,&
        &4, 3, 2, 1, 8, 7, 6, 5,&
        &3, 4, 1, 2, 7, 8, 5, 6,&
        &2, 1, 4, 3, 6, 5, 8, 7,&
        &1, 2, 3, 4, 5, 6, 7, 8/

    integer(kind=I__) :: IOCTORDER(8)
    data IOCTORDER /1, 5, 2, 6, 3, 7, 4, 8/


    ! IOCT is the octant the ray come from or the discrete ordinate
    !   IOCT = 1 + BITX + 2*BITY + 4*BITZ

    if (is1D_) then
        NOCT = 2
    else if (is1DYx_) then
        NOCT = 4
    else
        NOCT = 8
    end if

    NXC = rtm_obj%Nx
    if (isOpenX_) NXC = rtm_obj%Nx + 1
    if (isMpiX_ .and. (.not. is1DX_)) NXC = rtm_obj%Nx - 1

    NYC = rtm_obj%Ny
    if (isOpenYx_) NYC = rtm_obj%Ny + 1
    if (isMpiY_ .and. (.not. is1DYx_)) NYC = rtm_obj%Ny - 1

    ! Loop over the number of octants needed; get the octant number
    do JOCT = 1, NOCT
        IOCT = IOCTORDER(JOCT)

        ! Set each grid point flag to indicate unvisited.
        do ipts = 1, grid_obj%Npts
            GRIDRAD(ipts) = -1.0
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
            if (.not. sweepCellNext(rtm_obj%flag_bc, NXC, NYC, rtm_obj%Nz, IOCT,&
                &IPCELL, grid_obj)) goto 900
        else
            INDEXCORN = INDEXCORN + 1
        end if

        ! Compute the corresponding point
        ! ICORNER = IAND(NOT(IEOR(INDEXCORN-1,IOCT-1)),7) + 1
        ICORNER = CORNDOG(INDEXCORN, IOCT)
        ipts = grid_obj%i8Grid(ICORNER, IPCELL)
        if (GRIDRAD(ipts) .ge. 0.0) goto 100

        GRIDRAD(ipts) = 1.0

        ! Fill in the next entry with the grid cell and corner of this point
        SWEEPORD(IORDER, JOCT) = ior(ishft(IPCELL, 3), ICORNER - 1)
        IORDER = IORDER + 1
        goto 100

900     continue

    end do

    return

end subroutine shd_grid__makeSweepOrder

subroutine shd_grid__divideCell(&
    &icell, idir, newPts,&
    &grid_obj)

    !   Divides a single cell (icell) in the idir direction (1=Z,2=Y,3=X)
    ! Makes the pointers to the two new cells (grid_obj%i2Tree), set the flags
    ! (grid_obj%flag), creates the new grid points if they don't exist already,
    ! and updates the neighbor pointers (grid_obj%i6Cell) for the new cells and
    ! their neighboring cells.

    type(shd_type_grid), intent(inout) :: grid_obj

    integer icell, idir, newPts(3, 4)
    integer newCell, iface, I

    ! Only allow to divide an end cell (i.e. one not yet divided)
    if (grid_obj%i2Tree(2, icell) .ne. 0) then
        stop 'Error [_divideCell]: Cannot divide already split cell.'
    end if

    ! Do the tree pointers: parent points to negative child,
    !   children point back to parent.
    newCell = grid_obj%Ncell + 1
    grid_obj%Ncell = grid_obj%Ncell + 2
    grid_obj%i2Tree(2, icell) = newCell
    grid_obj%i2Tree(1, newCell) = icell
    grid_obj%i2Tree(2, newCell) = 0
    grid_obj%i2Tree(1, newCell + 1) = icell
    grid_obj%i2Tree(2, newCell + 1) = 0
    ! Do the flags: put the splitting direction in the parent and
    !   inherit the boundary flags for the children.
    grid_obj%flag(icell) = ior(int(grid_obj%flag(icell)), ishft(idir, 2))
    grid_obj%flag(newCell) = 0
    grid_obj%flag(newCell + 1) = 0
    do I = 0, 1
        if (btest(int(grid_obj%flag(icell)), I)) then
            grid_obj%flag(newCell) = ibset(int(grid_obj%flag(newCell)), I)
            grid_obj%flag(newCell + 1) = ibset(int(grid_obj%flag(newCell + 1)), I)
        end if
    end do

    ! Do the grid points:
    call shd_grid__makeNewGridPnt(idir, icell, newCell,&
        &newPts, grid_obj)

    ! Do the neighbor pointers:
    do iface = 1, 6
        ! Special case if the neighbor is the cell itself
        if (grid_obj%i6Cell(iface, icell) .eq. icell) then
            grid_obj%i6Cell(iface, newCell) = newCell
            grid_obj%i6Cell(iface, newCell + 1) = newCell + 1
        ! For the faces in the split direction, the children either
        !   point to the other child or inherit the parent's neighbor.
        !   Also update the neighbor's neighbor pointer.
        else if (iface .eq. 2*idir) then
            grid_obj%i6Cell(iface, newCell) = newCell + 1
            grid_obj%i6Cell(iface, newCell + 1) = grid_obj%i6Cell(iface, icell)
        else if (iface .eq. 2*idir - 1) then
            grid_obj%i6Cell(iface, newCell + 1) = newCell
            grid_obj%i6Cell(iface, newCell) = grid_obj%i6Cell(iface, icell)
        else
            ! For the other four faces:
            !   The children inherit the parent's neighbor pointers
            grid_obj%i6Cell(iface, newCell) = grid_obj%i6Cell(iface, icell)
            grid_obj%i6Cell(iface, newCell + 1) = grid_obj%i6Cell(iface, icell)
        end if
        ! Then see if we can match faces further down the neighbor's tree;
        !   also update the neighboring cells pointers
        call shd_grid__matchNeighCellFace(iface, newCell,&
            &grid_obj)
        call shd_grid__matchNeighCellFace(iface, newCell + 1,&
            &grid_obj)
    end do

    return

end subroutine shd_grid__divideCell

subroutine shd_grid__matchNeighCellFace(iface, IC,&
    &grid_obj)

    !   Resolves the neighbor pointers for a new cell IC (from a subdivision)
    ! for the iface face.  Based on the size of this new face and
    ! the corresponding faces of the neighboring cells (down its tree),
    ! the new cell's grid_obj%i6Cell and those of the neighboring cells are updated.

    type(shd_type_grid), intent(in) :: grid_obj

    integer iface, IC
    integer icellNeigh, INN, JFACE, IC1, IC8, IN1, IN8, DIR, DIR1, DIR2
    logical DONE
    real POS(3)

    icellNeigh = abs(grid_obj%i6Cell(iface, IC))
    ! If no neighbor then return
    if (icellNeigh .eq. 0) return
    JFACE = (iface + 1)/2
    ! Get the center of the new face
    IC1 = grid_obj%i8Grid(1, IC)
    IC8 = grid_obj%i8Grid(8, IC)
    DIR1 = mod(JFACE - 1 + 1, 3) + 1
    DIR2 = mod(JFACE - 1 + 2, 3) + 1
    POS(DIR1) = (grid_obj%xyz(DIR1, IC1) + grid_obj%xyz(DIR1, IC8))/2
    POS(DIR2) = (grid_obj%xyz(DIR2, IC1) + grid_obj%xyz(DIR2, IC8))/2
    ! Go down the neighbors cell tree
    DONE = .false.
    do while (.not. DONE .and. IC .ne. icellNeigh)
        ! Compare the coordinates of this neighbor cell to see
        !   if the new face is contained within the neighbor face.
        IN1 = grid_obj%i8Grid(1, icellNeigh)
        IN8 = grid_obj%i8Grid(8, icellNeigh)
        if  (grid_obj%xyz(DIR1, IC1) .ge. grid_obj%xyz(DIR1, IN1) .and.&
            &grid_obj%xyz(DIR1, IC8) .le. grid_obj%xyz(DIR1, IN8) .and.&
            &grid_obj%xyz(DIR2, IC1) .ge. grid_obj%xyz(DIR2, IN1) .and.&
            &grid_obj%xyz(DIR2, IC8) .le. grid_obj%xyz(DIR2, IN8)) then
            ! If this is the bottom of the tree then we have a face match
            if (grid_obj%i2Tree(2, icellNeigh) .eq. 0) then
                grid_obj%i6Cell(iface, IC) =  icellNeigh
            else
                grid_obj%i6Cell(iface, IC) = -icellNeigh
            end if
        end if
        ! Compare the size of the neighbor face with the new face.
        !   Update the neighbor's neighbor pointer to the new face
        !   if the neighbor face is smaller, otherwise make the
        !   pointer negative, because it no longer points to an end node.
        if  (grid_obj%xyz(DIR1, IN1) .ge. grid_obj%xyz(DIR1, IC1) .and.&
            &grid_obj%xyz(DIR1, IN8) .le. grid_obj%xyz(DIR1, IC8) .and.&
            &grid_obj%xyz(DIR2, IN1) .ge. grid_obj%xyz(DIR2, IC1) .and.&
            &grid_obj%xyz(DIR2, IN8) .le. grid_obj%xyz(DIR2, IC8)) then
            call shd_grid__inheritNeighCell(icellNeigh, oppFace(iface), IC,&
                &grid_obj)
        else
            grid_obj%i6Cell(oppFace(iface), icellNeigh) = -abs(grid_obj%i6Cell(oppFace(iface), icellNeigh))
        end if

        ! Go down the neighbors cell tree until there is no children,
        !   finding the cell that contains the center of the new face.
        if (grid_obj%i2Tree(2, icellNeigh) .eq. 0) then
            DONE = .true.
        else
            DIR = ibits(int(grid_obj%flag(icellNeigh)), 2, 2)
            INN = grid_obj%i2Tree(2, icellNeigh)
            if (DIR .eq. JFACE) then
                icellNeigh = INN + 1 - mod(iface - 1, 2)
            else
                if (POS(DIR) .gt. grid_obj%xyz(DIR, grid_obj%i8Grid(8, INN))) then
                    icellNeigh = INN + 1
                else
                    icellNeigh = INN
                end if
            end if
        end if
    end do

    return

end subroutine shd_grid__matchNeighCellFace

subroutine shd_grid__inheritNeighCell(icell, iface, icellNeigh,&
    &grid_obj)

    !   Updates the icell cell's new neighbor cell (icellNeigh) for all
    ! children cells that border the face (iface).

    type(shd_type_grid), intent(in) :: grid_obj

    integer icell, iface, icellNeigh
    integer JFACE, IC, DIR, SP, STACK(50)
    logical DONE

    ! Go through the tree checking the cells:
    ! The stack holds the second child cells, which are the
    !   next cell to go to if at an end cell.
    JFACE = (iface + 1)/2
    IC = icell
    SP = 0
    DONE = .false.
    do while (.not. DONE)
        grid_obj%i6Cell(iface, IC) = icellNeigh
        ! If this is an end cell then pop the stack
        if (grid_obj%i2Tree(2, IC) .eq. 0) then
            ! If no more to pop then we are done
            if (SP .eq. 0) then
                DONE = .true.
            else
                ! Otherwise, pop the stack to get the next cell to go to
                IC = STACK(SP)
                SP = SP - 1
            end if
        else
            ! Not and end cell, so this cell has children.
            ! If split is in same direction as face then do only close child
            DIR = ibits(int(grid_obj%flag(IC)), 2, 2)
            if (DIR .eq. JFACE) then
                IC = grid_obj%i2Tree(2, IC) + mod(iface - 1, 2)
            else
                ! Otherwise need to go to both children so put the second child
                !   on the stack and set the pointer to go to the first child next.
                SP = SP + 1
                if (SP .gt. 50) stop 'Error [_inheritNeighCell]: Stack exceeded.'
                STACK(SP) = grid_obj%i2Tree(2, IC) + 1
                IC = grid_obj%i2Tree(2, IC)
            end if
        end if
    end do

    return

end subroutine shd_grid__inheritNeighCell

subroutine shd_grid__testGridSmooth(&
    &icell, idir,&
    &grid_obj)
    !   Tests the grid cell icell for splitting in order to improve the
    ! adaptive grid smoothness.  The direction to split the cell is returned
    ! in idir (1-X, 2-Y, 3-Z, 0-don't).
    ! The criterion for splitting the cell is
    !   1) cell has finer grid cells on either side
    !   2) a neighboring cell has grid spacing more than 3 times finer

    type(shd_type_grid), intent(in) :: grid_obj

    integer icell, idir
    integer ID, iface, IP, IP1, IP2, icellNeigh, I, J, ISUM
    integer INCELL(2), DIRSPLIT(2)
    integer FACECORNER(4, 6), edgeCorner(2, 3)
    real(RD_) :: XE, YE, ZE
    real RATIO, SIZERATIO, CURGRIDSIZE(3), GRIDSIZE(2), GRIDSIZEINV
    data FACECORNER/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/
    data edgeCorner/1, 2, 1, 3, 1, 5/

    idir = 0
    SIZERATIO = 1.0
    ! Loop over the three possible directions to do the split (X,Y,Z)
    do ID = 1, 3
        ! Get the grid points
        IP1 = grid_obj%i8Grid(edgeCorner(1, ID), icell)
        IP2 = grid_obj%i8Grid(edgeCorner(2, ID), icell)
        CURGRIDSIZE(ID) = 1.0e20
        if (IP1 .ne. IP2) then
            CURGRIDSIZE(ID) = abs(grid_obj%xyz(ID, IP1) - grid_obj%xyz(ID, IP2))
            GRIDSIZEINV = 1.0/CURGRIDSIZE(ID)
            ! Loop over two opposite faces
            do J = 1, 2
                DIRSPLIT(J) = 0
                iface = 2*(ID - 1) + J
                ! Get location of center of cell's face
                XE = 0.0
                YE = 0.0
                ZE = 0.0
                do I = 1, 4
                    IP = grid_obj%i8Grid(FACECORNER(I, iface), icell)
                    XE = XE + grid_obj%xyz(1, IP)*0.25
                    YE = YE + grid_obj%xyz(2, IP)*0.25
                    ZE = ZE + grid_obj%xyz(3, IP)*0.25
                end do
                ! Get the neighboring cell and its size in this direction
                call shd_util__findNextCell(XE, YE, ZE, iface, ID, icell,&
                    &INCELL(J), grid_obj)
                ! If the neighbor does not exist or is a base cell or
                !   is an independent pixel cell in X or Y then set
                !   its grid size the same as current cell (so don't split
                if (INCELL(J) .eq. 0) then
                    GRIDSIZE(J) = CURGRIDSIZE(ID)
                else
                    GRIDSIZE(J) = abs(grid_obj%xyz(ID, grid_obj%i8Grid(1, INCELL(J)))&
                        &- grid_obj%xyz(ID, grid_obj%i8Grid(8, INCELL(J))))

                    if (grid_obj%i2Tree(1, INCELL(J)) .eq. 0) then
                        GRIDSIZE(J) = CURGRIDSIZE(ID)
                    end if
                    if (ID .le. 2 .and. btest(grid_obj%flag(INCELL(J)), ID - 1)) then
                        GRIDSIZE(J) = CURGRIDSIZE(ID)
                    end if
                end if
                ! If there are multiple neighboring cells then get
                !   direction the neighbor cell was split
                icellNeigh = grid_obj%i6Cell(IFACE, icell)
                if (icellNeigh .lt. 0) then
                    DIRSPLIT(J) = ibits(int(grid_obj%flag(abs(icellNeigh))), 2, 2)
                end if
            end do

            ! Split cell this way if cells on opposite side are smaller in
            !   this direction OR opposite neighbor cells are split in
            !   same other direction.
            if  (GRIDSIZE(1)*GRIDSIZEINV .lt. 0.75 .and.&
                &GRIDSIZE(2)*GRIDSIZEINV .lt. 0.75) then
                idir = ID
            end if
            if (DIRSPLIT(1) .gt. 0 .and. DIRSPLIT(1) .ne. ID .and.&
                &DIRSPLIT(1) .eq. DIRSPLIT(2)) then
                idir = DIRSPLIT(1)
            end if

            ! Split cell this way if ratio of neighboring grid cells' size
            !   to this cells size (in ID direction) is less than 0.4.
            RATIO = min(GRIDSIZE(1)*GRIDSIZEINV, GRIDSIZE(2)*GRIDSIZEINV)
            if (RATIO .lt. 0.4 .and. RATIO .lt. SIZERATIO) then
                ! Check if there is any source function for this cell, before
                !   deciding to split it (no point in dividing clear sky).
                ISUM = 0
                do I = 1, 8
                    IP = grid_obj%i8Grid(I, icell)
                    ISUM = ISUM + isrcJ_Shdo(IP + 1) - isrcJ_Shdo(IP)
                end do
                if (ISUM .gt. 0) then
                    idir = ID
                    SIZERATIO = RATIO
                end if
            end if

        end if
    end do

    return

end subroutine shd_grid__testGridSmooth

subroutine shd_grid__testCellSplit(icell,&
    &ADAPTCRIT, MAXADAPT, idir,&
    &rtm_obj, grid_obj)

    !   Computes the adaptive cell splitting criterion for the X, Y, and Z
    ! directions for cell (icell).  The three criterion are returned in
    ! ADAPTCRIT, the max criterion is returned in MAXADAPT, and the
    ! direction of max criterion in idir.
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

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer icell, idir
    real ADAPTCRIT(3), MAXADAPT
    integer ID, IE, IP1, IP2, IS1, IS2, NS1, NS2, NS
    integer J, NUM
    real TOSPLIT, SUM1, JAY, TAU, EXT
    real, parameter :: C0 = 0.28209479177387814 ! sqrt(1.0/(4.0*PI_)), Y00 (l=0, m=0)

    MAXADAPT = -1.0
    ! Loop over the three possible directions to do the split (X,Y,Z)
    do ID = 1, 3
        SUM1 = 0.0
        NUM = 0
        ! Loop over four edges that cross the ID direction
        do IE = 1, 4
            ! Get the grid points
            IP1 = grid_obj%i8Grid(gridCorner(1, IE, ID), icell)
            IP2 = grid_obj%i8Grid(gridCorner(2, IE, ID), icell)
            if (IP1 .ne. IP2) then
                NUM = NUM + 1
                IS1 = isrcJ_Shdo(IP1)
                IS2 = isrcJ_Shdo(IP2)
                NS1 = isrcJ_Shdo(IP1 + 1) - IS1
                NS2 = isrcJ_Shdo(IP2 + 1) - IS2
                JAY = 0.0
                NS = min(NS1, NS2)
                !   Compute the source function rms difference integrated
                ! over all angles, which is sum of squares of spherical
                ! harmonic term differences.
                do J = 1, NS
                    JAY = JAY + (rtm_obj%ext(IP2)*srcJ_Shdo(1, IS2 + J)&
                        &- rtm_obj%ext(IP1)*srcJ_Shdo(1, IS1 + J))**2
                end do
                do J = NS + 1, NS1
                    JAY = JAY + (rtm_obj%ext(IP1)*srcJ_Shdo(1, IS1 + J))**2
                end do
                do J = NS + 1, NS2
                    JAY = JAY + (rtm_obj%ext(IP2)*srcJ_Shdo(1, IS2 + J))**2
                end do
                EXT = 0.5*(rtm_obj%ext(IP1) + rtm_obj%ext(IP2))
                if (EXT .gt. 0.0) then
                    JAY = C0*sqrt(JAY)/EXT
                else
                    JAY = 0.0
                end if
                TAU = abs(EXT*(grid_obj%xyz(ID, IP2) - grid_obj%xyz(ID, IP1)))
                TOSPLIT = abs(JAY)*(1 - exp(-TAU))
                SUM1 = SUM1 + TOSPLIT
            end if
        end do

        if (NUM .gt. 0) then
            ADAPTCRIT(ID) = SUM1/NUM
        else
            ADAPTCRIT(ID) = 0.0
        end if
        if (ADAPTCRIT(ID) .gt. MAXADAPT) then
            MAXADAPT = ADAPTCRIT(ID)
            idir = ID
        end if
    end do

    return

end subroutine shd_grid__testCellSplit

logical function sweepCellBase(BCFLAG, NXC, NYC, NZ, IOCT, icell)

    !   Returns the next base grid cell.  The octant the ray direction
    ! is coming from (IOCT) determines the order of sweeping.  icell is 0
    ! for the first call to set things up.  The function returns false if
    ! we are done sweeping the grid.  Do loops over IZ, IY, and IX are
    ! simulated.  For open boundary condition do the endpoint cells
    ! (in X,Y) first, so the boundary grid points will be done with
    ! independent pixel radiative transfer.

    integer BCFLAG, NXC, NYC, NZ, IOCT, icell
    integer IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
    save IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ

    sweepCellBase = .true.

    ! If first time here then set up the loops
    if (icell .eq. 0) then
        if (btest(IOCT - 1, 0)) then
            DIX = +1
            if (isOpenX_) then
                SIX = NXC
                EIX = max(NXC - 1, 1)
            else
                SIX = 1
                EIX = NXC
            end if
        else
            DIX = -1
            if (isOpenX_) then
                SIX = 1
                EIX = min(2, NXC)
            else
                SIX = NXC
                EIX = 1
            end if
        end if
        if (btest(IOCT - 1, 1)) then
            DIY = +1
            if (isOpenYx_) then
                SIY = NYC
                EIY = max(NYC - 1, 1)
            else
                SIY = 1
                EIY = NYC
            end if
        else
            DIY = -1
            if (isOpenYx_) then
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
    icell = IZ + (NZ - 1)*(IY - 1) + (NZ - 1)*NYC*(IX - 1)

    return

end function sweepCellBase

logical function sweepCellNext(BCFLAG, NXC, NYC, NZ, IOCT,&
    &icell,&
    &grid_obj)

    !   Returns the next grid cell given the last cell (icell).  The
    ! octant the ray direction is coming from (IOCT) determines the
    ! order of sweeping.  icell is 0 for the first call to set things up
    ! The function returns false if we are done sweeping the grid.
    ! Within a base cell the tree structure is traced using a stack until
    ! we get to a new end cell.

    type(shd_type_grid), intent(in) :: grid_obj

    integer BCFLAG, NXC, NYC, NZ, IOCT, icell
    integer IC, idir, SP, STACK(50)
    logical DONE
    save SP, STACK

    IC = icell

    ! Special case for first time here
    if (icell .eq. 0) then
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
        if (grid_obj%i2Tree(2, IC) .eq. 0) then
            if (IC .ne. icell) then
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
            idir = ibits(int(grid_obj%flag(IC)), 2, 2) - 1
            if (btest(IOCT - 1, idir)) then
                STACK(SP) = grid_obj%i2Tree(2, IC) + 1
                IC = grid_obj%i2Tree(2, IC)
            else
                STACK(SP) = grid_obj%i2Tree(2, IC)
                IC = grid_obj%i2Tree(2, IC) + 1
            end if
        end if
    end do

    icell = IC
    sweepCellNext = .true.

    return

end function sweepCellNext

subroutine shd_grid__interpGrid(rtm_obj, grid_obj)

    ! Calls shd_grid__InterpTriLinGridPrp to interpolate the input arrays from
    !   the property grid to each internal grid point.

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer(kind=I__) :: ipts

    ! Initialize: transfer the tabulated phase functions
    call shd_grid__InterpTriLinGridPrp(&
        &0.0, 0.0, 0.0,&
        &.true., rtm_obj%temp(1), rtm_obj%ext(1), rtm_obj%ssa(1),&
        &rtm_obj%leg(1:prp_obj_%NlegSet, 0:prp_obj_%Nleg, 1), rtm_obj%ipha(1))

    ! Trilinearly interpolate from the property grid to the adaptive
    do ipts = 1, grid_obj%Npts
        call shd_grid__InterpTriLinGridPrp(&
            &grid_obj%xyz(1, ipts), grid_obj%xyz(2, ipts), grid_obj%xyz(3, ipts),&
            &.false., rtm_obj%temp(ipts), rtm_obj%ext(ipts), rtm_obj%ssa(ipts),&
            &rtm_obj%leg(1:prp_obj_%NlegSet, 0:prp_obj_%Nleg, ipts), rtm_obj%ipha(ipts))
    end do

    return
end subroutine shd_grid__interpGrid

subroutine shd_grid__InterpTriLinGridPrp(&
    &x0, y0, z0, doInit,&
    &temp0, ext0, ssa0, leg0, ipha0)

    !  Trilinearly interpolates the quantities on the input property
    ! grid at the single point (x0,y0,z0) to get the output temp0,ext0,
    ! ssa0, and leg0 or ipha0.  Interpolation is done on the
    ! volume coefficients.  Also adds in the separate gaseous absorption.
    ! Divides the phase function Legendre coefficients leg0 by 2*l+1.
    ! The phase function pointer ipha0 (for tabulated phase functions)
    ! is that of the maximum weighted scattering property grid point.
    ! If doInit=.TRUE. then transfers the tabulated phase functions.


    real(kind=R__), intent(in)  :: x0, y0, z0
    logical, intent(in) :: doInit

    integer(kind=2), intent(out) :: ipha0
    real(kind=R__),  intent(out) :: temp0, ext0, ssa0, leg0(prp_obj_%NlegSet, 0:prp_obj_%Nleg, *)

    integer(kind=I__) :: IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
    integer(kind=I__) :: I1, I2, I3, I4, I5, I6, I7, I8, I

    real(kind=RD_) :: U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
    real(kind=RD_) :: scat1, scat2, scat3, scat4, scat5, scat6, scat7, scat8
    real(kind=RD_) :: scat0, minScat, maxScat
    real(kind=RD_), save :: kabs0, minExt

    integer :: iSelect

    if (doInit) then
       ! If there are tabulated phase functions, then transfer them
        do I = 1, prp_obj_%Npha
            do L = 0, prp_obj_%Nleg
                do J = 1, prp_obj_%NlegSet
                    !leg0(J, L, I) = leg_Prp(J + prp_obj_%NlegSet*(L + (prp_obj_%Nleg+ 1)*(I - 1)))/(2*L + 1)
                    leg0(J, L, I) = leg_Prp(J, L, I)/(2*L + 1)
                end do
            end do
        end do
        !minExt = 1.0d-5/((prp_obj_%z(prp_obj_%Nz) - prp_obj_%z(1))/prp_obj_%Nz)
        minExt = accSolv_/((prp_obj_%z(prp_obj_%Nz) - prp_obj_%z(1))/prp_obj_%Nz) ! attempt to improve agreement between libRadtran in water vapor abosprtion bands
        minScat = 0.1d0*minExt
        return
    end if

    ! Find the grid location and compute the interpolation factors
    IL = 0
    IU = prp_obj_%Nz
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (z0 .ge. prp_obj_%z(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)
    W = dble(z0 - prp_obj_%z(IZ))/(prp_obj_%z(IZ + 1) - prp_obj_%z(IZ))
    W = max(min(W, 1.0d0), 0.0d0)
    IX = int((x0 - prp_obj_%xStart)/prp_obj_%dx) + 1
    if (abs(x0 - prp_obj_%xStart - prp_obj_%Nx*prp_obj_%dx) .lt. 0.01*prp_obj_%dx) IX = prp_obj_%Nx
    if (IX .lt. 1 .or. IX .gt. prp_obj_%Nx) then
        write (6, *) 'Error [_InterpTriLinGridPrp]: Beyond X domain.'
        stop
    end if
    IXP = mod(IX, prp_obj_%Nx) + 1
    U = dble(x0 - prp_obj_%xStart - prp_obj_%dx*(IX - 1))/prp_obj_%dx
    U = max(min(U, 1.0d0), 0.0d0)
    if (U .lt. 1.0d-5) U = 0.0d0
    if (U .gt. 1.0d0 - 1.0d-5) U = 1.0d0
    IY = int((y0 - prp_obj_%yStart)/prp_obj_%dy) + 1
    if (abs(y0 - prp_obj_%yStart - prp_obj_%Ny*prp_obj_%dy) .lt. 0.01*prp_obj_%dy) IY = prp_obj_%Ny
    if (IY .lt. 1 .or. IY .gt. prp_obj_%Ny) then
        write (6, *) 'Error [_InterpTriLinGridPrp]: Beyond Y domain.'
        stop
    end if
    IYP = mod(IY, prp_obj_%Ny) + 1
    V = dble(y0 - prp_obj_%yStart - prp_obj_%dy*(IY - 1))/prp_obj_%dy
    V = max(min(V, 1.0d0), 0.0d0)
    if (V .lt. 1.0d-5) V = 0.0d0
    if (V .gt. 1.0d0 - 1.0d-5) V = 1.0d0

    F1 = (1.0d0 - U)*(1.0d0 - V)*(1.0d0 - W)
    F2 =          U *(1.0d0 - V)*(1.0d0 - W)
    F3 = (1.0d0 - U)*         V *(1.0d0 - W)
    F4 =          U *         V *(1.0d0 - W)
    F5 = (1.0d0 - U)*(1.0d0 - V)*         W
    F6 =          U *(1.0d0 - V)*         W
    F7 = (1.0d0 - U)*         V *         W
    F8 =          U *         V *         W
    I1 = IZ + prp_obj_%Nz*( IY - 1) + prp_obj_%Nz*prp_obj_%Ny*( IX - 1)
    I2 = IZ + prp_obj_%Nz*( IY - 1) + prp_obj_%Nz*prp_obj_%Ny*(IXP - 1)
    I3 = IZ + prp_obj_%Nz*(IYP - 1) + prp_obj_%Nz*prp_obj_%Ny*( IX - 1)
    I4 = IZ + prp_obj_%Nz*(IYP - 1) + prp_obj_%Nz*prp_obj_%Ny*(IXP - 1)

    I5 = I1 + 1
    I6 = I2 + 1
    I7 = I3 + 1
    I8 = I4 + 1

    ! zgrid correction (to address layer averaging issue in z direction)
    !╭────────────────────────────────────────────────╮!
    I1 = I5
    I2 = I6
    I3 = I7
    I4 = I8
    !╰────────────────────────────────────────────────╯!

    scat1 = F1*ext_Prp(I1)*ssa_Prp(I1)
    scat2 = F2*ext_Prp(I2)*ssa_Prp(I2)
    scat3 = F3*ext_Prp(I3)*ssa_Prp(I3)
    scat4 = F4*ext_Prp(I4)*ssa_Prp(I4)
    scat5 = F5*ext_Prp(I5)*ssa_Prp(I5)
    scat6 = F6*ext_Prp(I6)*ssa_Prp(I6)
    scat7 = F7*ext_Prp(I7)*ssa_Prp(I7)
    scat8 = F8*ext_Prp(I8)*ssa_Prp(I8)

    maxScat = -1.0
    if (scat1 .gt. maxScat .or. abs(F1 - 1) .lt. 0.001) then
        maxScat = scat1
        iSelect = I1
    end if
    if (scat2 .gt. maxScat .or. abs(F2 - 1) .lt. 0.001) then
        maxScat = scat2
        iSelect = I2
    end if
    if (scat3 .gt. maxScat .or. abs(F3 - 1) .lt. 0.001) then
        maxScat = scat3
        iSelect = I3
    end if
    if (scat4 .gt. maxScat .or. abs(F4 - 1) .lt. 0.001) then
        maxScat = scat4
        iSelect = I4
    end if
    if (scat5 .gt. maxScat .or. abs(F5 - 1) .lt. 0.001) then
        maxScat = scat5
        iSelect = I5
    end if
    if (scat6 .gt. maxScat .or. abs(F6 - 1) .lt. 0.001) then
        maxScat = scat6
        iSelect = I6
    end if
    if (scat7 .gt. maxScat .or. abs(F7 - 1) .lt. 0.001) then
        maxScat = scat7
        iSelect = I7
    end if
    if (scat8 .gt. maxScat .or. abs(F8 - 1) .lt. 0.001) then
        maxScat = scat8
        iSelect = I8
    end if

    !temp0 = temp_Prp(iSelect)
    !ext0 = ext_Prp(iSelect)
    !scat0 = ext_Prp(iSelect)*ssa_Prp(iSelect)

    temp0 = F1*temp_Prp(I1) + F2*temp_Prp(I2) + F3*temp_Prp(I3) + F4*temp_Prp(I4)&
        &+ F5*temp_Prp(I5) + F6*temp_Prp(I6) + F7*temp_Prp(I7) + F8*temp_Prp(I8)
    ext0 = F1*ext_Prp(I1) + F2*ext_Prp(I2) + F3*ext_Prp(I3)&
        &+ F4*ext_Prp(I4) + F5*ext_Prp(I5) + F6*ext_Prp(I6)&
        &+ F7*ext_Prp(I7) + F8*ext_Prp(I8)
    scat0 = scat1 + scat2 + scat3 + scat4 + scat5 + scat6 + scat7 + scat8

    if (ext0 .gt. minExt) then
        ssa0 = scat0/ext0
    else
        ssa0 = scat0/minExt
    end if

    if (prp_obj_%Npha .gt. 0) then
        ipha0 = ipha_Prp(iSelect)
    else
        ! Standard property file format: unpolarized
        leg0(1, 0, 1) = 1.0
        do L = 1, prp_obj_%Nleg
            !leg0(1, L, 1) = scat1*leg_Prp(L + (prp_obj_%Nleg + 1)*(I1 - 1))&
            !    &+ scat2*leg_Prp(L + (prp_obj_%Nleg + 1)*(I2 - 1))&
            !    &+ scat3*leg_Prp(L + (prp_obj_%Nleg + 1)*(I3 - 1))&
            !    &+ scat4*leg_Prp(L + (prp_obj_%Nleg + 1)*(I4 - 1))&
            !    &+ scat5*leg_Prp(L + (prp_obj_%Nleg + 1)*(I5 - 1))&
            !    &+ scat6*leg_Prp(L + (prp_obj_%Nleg + 1)*(I6 - 1))&
            !    &+ scat7*leg_Prp(L + (prp_obj_%Nleg + 1)*(I7 - 1))&
            !    &+ scat8*leg_Prp(L + (prp_obj_%Nleg + 1)*(I8 - 1))
            leg0(1, L, 1) = scat1*leg_Prp(1, L, I1)&
                         &+ scat2*leg_Prp(1, L, I2)&
                         &+ scat3*leg_Prp(1, L, I3)&
                         &+ scat4*leg_Prp(1, L, I4)&
                         &+ scat5*leg_Prp(1, L, I5)&
                         &+ scat6*leg_Prp(1, L, I6)&
                         &+ scat7*leg_Prp(1, L, I7)&
                         &+ scat8*leg_Prp(1, L, I8)
            if (scat0 .gt. minScat) then
                leg0(1, L, 1) = leg0(1, L, 1)/scat0
            else
                leg0(1, L, 1) = leg0(1, L, 1)/minScat
            end if
            leg0(1, L, 1) = leg0(1, L, 1)/(2*L + 1)
        end do
    end if

    ! Trilinearly interpolate the temperature, extinction, scattering
    ! original
    !╭────────────────────────────────────────────────╮!
    !temp0 = F1*temp_Prp(I1) + F2*temp_Prp(I2) + F3*temp_Prp(I3) + F4*temp_Prp(I4)&
    !    &+ F5*temp_Prp(I5) + F6*temp_Prp(I6) + F7*temp_Prp(I7) + F8*temp_Prp(I8)
    !ext0 = F1*ext_Prp(I1) + F2*ext_Prp(I2) + F3*ext_Prp(I3)&
    !    &+ F4*ext_Prp(I4) + F5*ext_Prp(I5) + F6*ext_Prp(I6)&
    !    &+ F7*ext_Prp(I7) + F8*ext_Prp(I8)
    !scat1 = F1*ext_Prp(I1)*ssa_Prp(I1)
    !scat2 = F2*ext_Prp(I2)*ssa_Prp(I2)
    !scat3 = F3*ext_Prp(I3)*ssa_Prp(I3)
    !scat4 = F4*ext_Prp(I4)*ssa_Prp(I4)
    !scat5 = F5*ext_Prp(I5)*ssa_Prp(I5)
    !scat6 = F6*ext_Prp(I6)*ssa_Prp(I6)
    !scat7 = F7*ext_Prp(I7)*ssa_Prp(I7)
    !scat8 = F8*ext_Prp(I8)*ssa_Prp(I8)
    !scat0 = scat1 + scat2 + scat3 + scat4 + scat5 + scat6 + scat7 + scat8
    !if (ext0 .gt. minExt) then
    !    ssa0 = scat0/ext0
    !else
    !    ssa0 = scat0/minExt
    !end if
    !╰────────────────────────────────────────────────╯!

    ! For tabulated phase functions pick the one we are on top of
    !   or the one with the most scattering weight.
    ! original
    !╭────────────────────────────────────────────────╮!
    !if (prp_obj_%Npha .gt. 0) then
    !    maxScat = -1.0
    !    if (scat1 .gt. maxScat .or. abs(F1 - 1) .lt. 0.001) then
    !        maxScat = scat1
    !        ipha0 = ipha_Prp(I1)
    !    end if
    !    if (scat2 .gt. maxScat .or. abs(F2 - 1) .lt. 0.001) then
    !        maxScat = scat2
    !        ipha0 = ipha_Prp(I2)
    !    end if
    !    if (scat3 .gt. maxScat .or. abs(F3 - 1) .lt. 0.001) then
    !        maxScat = scat3
    !        ipha0 = ipha_Prp(I3)
    !    end if
    !    if (scat4 .gt. maxScat .or. abs(F4 - 1) .lt. 0.001) then
    !        maxScat = scat4
    !        ipha0 = ipha_Prp(I4)
    !    end if
    !    if (scat5 .gt. maxScat .or. abs(F5 - 1) .lt. 0.001) then
    !        maxScat = scat5
    !        ipha0 = ipha_Prp(I5)
    !    end if
    !    if (scat6 .gt. maxScat .or. abs(F6 - 1) .lt. 0.001) then
    !        maxScat = scat6
    !        ipha0 = ipha_Prp(I6)
    !    end if
    !    if (scat7 .gt. maxScat .or. abs(F7 - 1) .lt. 0.001) then
    !        maxScat = scat7
    !        ipha0 = ipha_Prp(I7)
    !    end if
    !    if (scat8 .gt. maxScat .or. abs(F8 - 1) .lt. 0.001) then
    !        maxScat = scat8
    !        ipha0 = ipha_Prp(I8)
    !    end if
    !else
    !    ! Standard property file format: unpolarized
    !    leg0(1, 0, 1) = 1.0
    !    do L = 1, prp_obj_%Nleg
    !        leg0(1, L, 1) = scat1*leg_Prp(L + (prp_obj_%Nleg + 1)*(I1 - 1))&
    !            &+ scat2*leg_Prp(L + (prp_obj_%Nleg + 1)*(I2 - 1))&
    !            &+ scat3*leg_Prp(L + (prp_obj_%Nleg + 1)*(I3 - 1))&
    !            &+ scat4*leg_Prp(L + (prp_obj_%Nleg + 1)*(I4 - 1))&
    !            &+ scat5*leg_Prp(L + (prp_obj_%Nleg + 1)*(I5 - 1))&
    !            &+ scat6*leg_Prp(L + (prp_obj_%Nleg + 1)*(I6 - 1))&
    !            &+ scat7*leg_Prp(L + (prp_obj_%Nleg + 1)*(I7 - 1))&
    !            &+ scat8*leg_Prp(L + (prp_obj_%Nleg + 1)*(I8 - 1))
    !        if (scat0 .gt. minScat) then
    !            leg0(1, L, 1) = leg0(1, L, 1)/scat0
    !        else
    !            leg0(1, L, 1) = leg0(1, L, 1)/minScat
    !        end if
    !        leg0(1, L, 1) = leg0(1, L, 1)/(2*L + 1)
    !    end do
    !end if
    !╰────────────────────────────────────────────────╯!

    ! Add in the gaseous absorption to extinction and albedo
    if (ckd_obj_%Nz .gt. 0) then
        IL = 1
        IU = ckd_obj_%Nz
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (z0 .le. ckd_obj_%z(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
        I = min(max(IL, 1), ckd_obj_%Nz - 1)
        F = (z0 - ckd_obj_%z(I))/(ckd_obj_%z(I + 1) - ckd_obj_%z(I))
        F = min(max(F, 0.0d0), 1.0d0)
        kabs0 = (1.0 - F)*kabs_Ckd(I, ckd_obj_%ig) + F*kabs_Ckd(I + 1, ckd_obj_%ig)
        kabs0 = (1.0 - F)*kabs_Ckd(I, ckd_obj_%ig) + F*kabs_Ckd(I + 1, ckd_obj_%ig)
        if ((ext0 + kabs0) .gt. 0.0) then
            ssa0 = ssa0*ext0/(ext0 + kabs0)
        else
            ssa0 = 0.0
        end if
        ext0 = ext0 + kabs0
    end if

    return

end subroutine shd_grid__InterpTriLinGridPrp

subroutine shd_grid__interpNewGridPnt(newPts,&
    &doDeltaM, ACCELFLAG,&
    &SOLARFLUX, YLMSUN,&
    &DIRFLUX,&
    &OSHPTR,&
    &rtm_obj, grid_obj)

    !   Interpolates the medium properties, radiance, and source function for
    ! new grid points (specified in newPts).  The medium properties are
    ! interpolated from the property grid, since interpolating from the
    ! internal grid may not faithfully represent the true medium.
    ! The delta-M scaling is done to the new grid point properties if needed.
    ! The direct beam flux and the Planck function for the new points is
    ! recomputed (rather than interpolated).  The radiance is averaged
    ! from the two parent grid points, and used to compute the source
    ! function for the new points.

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer newPts(3, 4)
    integer OSHPTR(*)
    logical doDeltaM, ACCELFLAG
    real SOLARFLUX, YLMSUN(prp_obj_%NlegSet, shdo_obj_%Nlm)
    real DIRFLUX(grid_obj%Npts)
    integer I, IP1, IP2, IP, IPH, J, K, L, M, ME, MS
    integer IS, NS, NS1, NS2
    integer IR, IR1, IR2, NR, NR1, NR2, SIDE
    logical VALIDBEAM
    real X, Y, Z, F, BB, C, D, SECMU0
    real RAD1(shdo_obj_%Nstokes), RAD2(shdo_obj_%Nstokes)
    real UNIFZLEV, XO, YO, ZO, DIRPATH
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)

    allocate(LOFJ(shdo_obj_%Nlm), SOURCET(shdo_obj_%Nstokes, shdo_obj_%Nlm))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, shdo_obj_%maxNl
        ME = min(L, shdo_obj_%maxNm)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    C = sqrt(4.0*acos(-1.0))
    if (isSolar_) then
        SECMU0 = 1.0/abs(src_obj_%avgMu)
    end if

    ! Loop over the (potentially) four new grid points
    do I = 1, 4
        if (newPts(3, I) .gt. 0) then
            ! Get the parent and new grid pointers
            IP1 = newPts(1, I)
            IP2 = newPts(2, I)
            IP  = newPts(3, I)

            ! Interpolate the medium properties from the property grid
            X = grid_obj%xyz(1, IP)
            Y = grid_obj%xyz(2, IP)
            Z = grid_obj%xyz(3, IP)
            call shd_grid__InterpTriLinGridPrp(X, Y, Z, .false.,&
                &rtm_obj%temp(IP), rtm_obj%ext(IP), rtm_obj%ssa(IP),&
                &rtm_obj%leg(1:prp_obj_%NlegSet, 0:prp_obj_%Nleg, 1), rtm_obj%ipha(IP))
            ! Do the Delta-M scaling of extinction and albedo for this point
            if (doDeltaM) then
                if (prp_obj_%Npha .gt. 0) then
                    F = rtm_obj%leg(1, shdo_obj_%maxNl + 1, rtm_obj%ipha(IP))
                else
                    F = rtm_obj%leg(1, shdo_obj_%maxNl + 1, IP)
                    do L = 0, shdo_obj_%maxNl
                        rtm_obj%leg(1, L, IP) = (rtm_obj%leg(1, L, IP) - F)/(1 - F)
                    end do
                end if
                rtm_obj%ext(IP) = (1.0 - rtm_obj%ssa(IP)*F)*rtm_obj%ext(IP)
                rtm_obj%ssa(IP) = (1.0 - F)*rtm_obj%ssa(IP)/(1.0 - rtm_obj%ssa(IP)*F)
            end if

            ! Compute the new Planck source function (if needed)
            if (isThermal_) then
                call shd_util__computePlanck(rtm_obj%temp(IP), src_obj_%unit, src_obj_%waveNum, src_obj_%avgWaveLen, BB)
                rtm_obj%planck(IP) = (1.0 - rtm_obj%ssa(IP))*BB
            end if

            ! Compute the new direct beam flux (if needed)
            if (isSolar_) then
                if (isMpi_) then
                    ! If MPI mode then calculate a temporary approximation
                    !   to the direct flux (for source function calculation below)
                    DIRFLUX(IP) = exp(0.5*(log(DIRFLUX(IP1)) + log(DIRFLUX(IP2))))
                else
                    ! Otherwise, calculate the exact direct beam from property grid
                    DIRPATH = 0.0
                    call shd_rad__computeDirect(0, X, Y, Z, rtm_obj%flag_bc, rtm_obj%flag_ip, doDeltaM,&
                        &shdo_obj_%maxNl, SOLARFLUX, src_obj_%avgMu, src_obj_%avgAz,&
                        &DIRFLUX(IP),&
                        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                end if
            end if

            ! Interpolate the radiance
            IR1 = iradI_Shdo(IP1)
            IR2 = iradI_Shdo(IP2)
            NR1 = iradI_Shdo(IP1 + 1) - iradI_Shdo(IP1)
            NR2 = iradI_Shdo(IP2 + 1) - iradI_Shdo(IP2)
            NR = max(NR1, NR2)
            IR = iradI_Shdo(IP)
            iradI_Shdo(IP + 1) = IR + NR
            do J = 1, NR
                if (J .le. NR1) then
                    RAD1(:) = radI_Shdo(:, IR1 + J)
                else
                    RAD1(:) = 0.0
                end if
                if (J .le. NR2) then
                    RAD2(:) = radI_Shdo(:, IR2 + J)
                else
                    RAD2(:) = 0.0
                end if
                ! this is modified to address the layer averaging issue
                !╭────────────────────────────────────────────────╮!
                !if ((abs(grid_obj%xyz(1, IP1)-grid_obj%xyz(1, IP2)) .lt. 1.0e-6) .and.&
                !    &(abs(grid_obj%xyz(2, IP1)-grid_obj%xyz(2, IP2)) .lt. 1.0e-6)) then
                !    if (grid_obj%xyz(3, IP2) .gt. grid_obj%xyz(3, IP1)) then
                !        radI_Shdo(:, IR + J) = RAD2(:)
                !    else
                !        radI_Shdo(:, IR + J) = RAD1(:)
                !    end if
                !else
                    !radI_Shdo(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
                !end if
                !╰────────────────────────────────────────────────╯!
                radI_Shdo(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
            end do

            ! Compute the source function for the new points
            NS1 = isrcJ_Shdo(IP1 + 1) - isrcJ_Shdo(IP1)
            NS2 = isrcJ_Shdo(IP2 + 1) - isrcJ_Shdo(IP2)
            NS = max(NS1, NS2)
            IS = isrcJ_Shdo(IP)
            isrcJ_Shdo(IP + 1) = IS + NS
            if (ACCELFLAG) OSHPTR(IP + 1) = OSHPTR(IP)
            if (prp_obj_%Npha .gt. 0) then
                IPH = rtm_obj%ipha(IP)
            else
                IPH = IP
            end if
            if (isPol_) then
                call shd_src__computeSrcOnePnt(shdo_obj_%Nstokes, shdo_obj_%Nlm, LOFJ,&
                    &src_obj_%type, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    &rtm_obj%planck(IP), rtm_obj%ssa(IP), IPH, rtm_obj%leg,&
                    &NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), SOURCET)
            else
                call shd_src__computeSrcOnePnt_(shdo_obj_%Nlm, LOFJ,&
                    &src_obj_%type, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    &rtm_obj%planck(IP), rtm_obj%ssa(IP), IPH, rtm_obj%leg,&
                    &NR, radI_Shdo(1:shdo_obj_%Nstokes, IR+1:NR), SOURCET)
            end if
            do J = 1, NS
                srcJ_Shdo(:, IS + J) = SOURCET(:, J)
            end do
        end if
    end do

    deallocate(LOFJ, SOURCET)

    return

end subroutine shd_grid__interpNewGridPnt

subroutine shd_grid__splitGrid(&
    &doSplit, isOutMem, CURSPLITACC, SPLITCRIT,&
    &doDeltaM, ACCELFLAG,&
    &OSHPTR,&
    &DIRFLUX,&
    &SOLARFLUX, YLMSUN,&
    &ADAPTCRIT, ADAPTIND,&
    &rtm_obj, grid_obj)

    !   Splits the cells that have a cell dividing criterion greater than
    ! CURSPLITACC if doSplit is true.  The current splitting criterion
    ! achieved is returned in SPLITCRIT.  If we are at the end of the
    ! grid point, grid cell, or spherical harmonic arrays (sizes given
    ! by grid_obj%maxNpts, grid_obj%maxNcell, shdo_obj_%maxNptsSH, shdo_obj_%maxNptsDO) then the isOutMem flag is returned true.

    type(shd_type_rtm),  intent(in)    :: rtm_obj
    type(shd_type_grid), intent(inout) :: grid_obj

    integer OSHPTR(*)
    real ADAPTIND(*)
    logical doDeltaM, ACCELFLAG, doSplit, isOutMem
    real CURSPLITACC, SPLITCRIT
    real SOLARFLUX, YLMSUN(prp_obj_%NlegSet, shdo_obj_%Nlm)
    real ADAPTCRIT(*)
    real DIRFLUX(*)
    integer icell, ICELL1, idir, I, N, newPts(3, 4), OLDNPTS
    integer MAXCELLS, MAXPTS, MAXWORK, MAXSH, MAXRAD
    logical isOutMem0
    real ADAPT(3), FRAC

    if (doSplit) then
        OLDNPTS = grid_obj%Npts
        ! Do the cells in batchs, so the newly split cells can be tested
        !   for further splitting
        ICELL1 = 1
        do while (ICELL1 .le. grid_obj%Ncell)
            ! For the current batch of cells, compute the adaptive grid
            !   splitting criterion for 3 directions for each cell.
            N = 0
            do icell = ICELL1, grid_obj%Ncell
                if (grid_obj%i2Tree(2, icell) .eq. 0) then
                    ! If this is a end node cell then compute adaptive criterion
                    call shd_grid__testCellSplit(icell,&
                        &ADAPT, ADAPTCRIT(N + 1), idir,&
                        &rtm_obj, grid_obj)
                    ADAPTIND(N+1) = 4*icell + idir
                    N = N + 1
                end if
            end do
            ICELL1 = grid_obj%Ncell + 1

            ! Sort the adaptive criterion in decreasing order
            call shd_util__ssort(ADAPTCRIT, ADAPTIND, N, -2)

            ! Figure how many cells we can split, leaving some room
            !   for cells split during grid smoothing
            FRAC = 0.03
            MAXCELLS = grid_obj%maxNcell - FRAC*(grid_obj%maxNcell - grid_obj%Ncell) - 2
            MAXPTS   = grid_obj%maxNpts  - FRAC*(grid_obj%maxNpts  - grid_obj%Npts)  - 4
            MAXWORK  = shdo_obj_%maxNptsDO - FRAC*(shdo_obj_%maxNptsDO - shdo_obj_%maxNphi0*grid_obj%Npts) - 4*shdo_obj_%maxNphi0
            MAXSH    = shdo_obj_%maxNptsSH - FRAC*(shdo_obj_%maxNptsSH - isrcJ_Shdo(grid_obj%Npts+1)) - 4*shdo_obj_%Nlm
            MAXRAD   = shdo_obj_%maxNptsSH + grid_obj%maxNpts - FRAC*(shdo_obj_%maxNptsSH + grid_obj%maxNpts - iradI_Shdo(grid_obj%Npts + 1)) - 4*shdo_obj_%Nlm

            ! Go down the list from highest adaptive cell criterion on down
            isOutMem0 = isOutMem
            I = 1
            do while ((.not. isOutMem0) .and. (ADAPTCRIT(I) .gt. CURSPLITACC)&
                &.and. (I .le. N))
                if  ((grid_obj%Ncell .gt. MAXCELLS) .or.&
                    &(grid_obj%Npts  .gt. MAXPTS) .or.&
                    &(grid_obj%Npts*shdo_obj_%maxNphi0   .gt. MAXWORK) .or.&
                    &(isrcJ_Shdo(grid_obj%Npts + 1) .gt. MAXSH)   .or.&
                    &(iradI_Shdo(grid_obj%Npts + 1) .gt. MAXRAD)) then
                    isOutMem0 = .true.
                else
                    icell = int(ADAPTIND(I))/4
                    idir = iand(int(ADAPTIND(I)), 3)
                    ! Create the new grid points and make pointers
                    call shd_grid__divideCell(icell, idir, newPts, grid_obj)

                    ! Interpolate the medium properties, radiance, and source
                    call shd_grid__interpNewGridPnt(newPts,&
                        &doDeltaM, ACCELFLAG,&
                        &SOLARFLUX, YLMSUN,&
                        &DIRFLUX,&
                        &OSHPTR,&
                        &rtm_obj, grid_obj)
                end if
                I = I + 1
            end do

            MAXCELLS = grid_obj%maxNcell - 2
            MAXPTS   = grid_obj%maxNpts  - 4
            MAXWORK  = shdo_obj_%maxNptsDO - 4*shdo_obj_%maxNphi0
            MAXSH    = shdo_obj_%maxNptsSH - 4*shdo_obj_%Nlm
            MAXRAD   = shdo_obj_%maxNptsSH + grid_obj%maxNpts - 4*shdo_obj_%Nlm

            ! Split some more cells to make the adaptive grid "smoother"
            do while ((.not. isOutMem) .and. (I .le. N))
                if  ((grid_obj%Ncell .gt. MAXCELLS) .or.&
                    &(grid_obj%Npts  .gt. MAXPTS) .or.&
                    &(grid_obj%Npts * shdo_obj_%maxNphi0 .gt. MAXWORK) .or.&
                    &(isrcJ_Shdo(grid_obj%Npts + 1) .gt. MAXSH) .or.&
                    &(iradI_Shdo(grid_obj%Npts + 1) .gt. MAXRAD)) then
                    isOutMem = .true.
                else
                    ! See if cell should be split to make better grid:
                    !   idir=0 if no, idir=direction (1,2,3) otherwise
                    icell = int(ADAPTIND(I))/4
                    call shd_grid__testGridSmooth(icell, idir, grid_obj)

                    ! Set idir=0 to stop adaptive grid smoothing
                    ! idir = 0
                    if (idir .gt. 0) then
                        ! Create the new grid points and make pointers
                        call shd_grid__divideCell(icell, idir, newPts, grid_obj)

                        ! Interpolate the medium properties, radiance, and source
                        call shd_grid__interpNewGridPnt(newPts,&
                            &doDeltaM, ACCELFLAG,&
                            &SOLARFLUX, YLMSUN,&
                            &DIRFLUX,&
                            &OSHPTR,&
                            &rtm_obj, grid_obj)
                    end if
                end if
                I = I + 1
            end do

        end do
        if (isOutMem0) then
            isOutMem = .true.
            !write (6, '(2A,I8,A,I8)')&
            !    &'Message [_splitGrid]: Further cell splitting requires larger array allocations: ',&
            !    &'grid_obj%Ncell= ', grid_obj%Ncell, ', Npnts = ', grid_obj%Npts
        end if

        ! If in parallel mode, calculate all the direct beam fluxes at once
        if (isSolar_) then
            if (isMpi_) then
                call shd_mpi__makeDirect(OLDNPTS + 1,&
                    &doDeltaM, shdo_obj_%maxNl,&
                    &SOLARFLUX, src_obj_%avgMu, src_obj_%avgAz,&
                    &DIRFLUX,&
                    &rtm_obj, grid_obj)
            end if
        end if
    end if

    ! Find the max adaptive cell criterion after the cell divisions
    SPLITCRIT = 0.0
    do icell = 1, grid_obj%Ncell
        if (grid_obj%i2Tree(2, icell) .eq. 0) then
            call shd_grid__testCellSplit(icell,&
                &ADAPT, ADAPTCRIT(1), idir,&
                &rtm_obj, grid_obj)
            SPLITCRIT = max(SPLITCRIT, ADAPTCRIT(1))
        end if
    end do

    return

end subroutine shd_grid__splitGrid

end module shd_grid
