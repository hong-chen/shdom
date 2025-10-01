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
!   Hong Chen, hong.chen-1@colorado.edu
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
!╰────────────────────────────────────────────────────────────────────────────╯!

module shd_grid

    use globals_const
    use globals_cfg
    use shd_util
    use shd_rad
    use shd_src
    use shd_mpi

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_grid__findBoundaryPnt
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

subroutine shd_grid__findBoundaryPnt(&
    &NPTS, NANG, MAXNBC, MAXBCRAD, LAMBERTIAN, &
    &ZBOT, ZTOP, GRIDPOS,&
    &NTOPPTS, NBOTPTS, BCPTR)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Returns the BCPTR pointer array for the NTOPPTS top boundary
    !   points (BCPTR(*,1)) and the NBOTPTS bottom boundary points (BCPTR(*,2)).
    !   All the grid points are looked at to find the boundary ones.
    !
    ! Input:
    !   NPTS      : integer, total number of grid points
    !   NANG      : integer, number of discrete ordinate angles (NMU and NPHI)
    !   MAXNBC    : integer, maximum number of boundary grid points
    !   MAXBCRAD  : integer, maximum number of radiances at boundary grid points
    !   LAMBERTIAN: logical, lambertian surface or not
    !   ZBOT      : real, altitude of bottom boundary
    !   ZTOP      : real, altitude of bottom boundary
    !   GRIDPOS(3, NPTS): real array, (x,y,z) coordinate values of grid point locations
    !
    ! Output:
    !   NTOPPTS   : integer, number of top boundary points
    !   NBOTPTS   : integer, number of bottom boundary points
    !   BCPTR(MAXNBC, 2): integer array, pointer array for NTOPPTS (top boundary points)
    !               and NBOTPTS (bottom boundary points)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=I__), intent(in) :: NPTS, NANG, MAXNBC, MAXBCRAD
    real(kind=R__)   , intent(in) :: ZBOT, ZTOP, GRIDPOS(3, NPTS)
    logical, intent(in) :: LAMBERTIAN

    integer(kind=I__), intent(out) :: NTOPPTS, NBOTPTS
    integer(kind=I__), intent(out) :: BCPTR(MAXNBC, 2)

    integer(kind=I__) :: I, IT, IB, NA

    if (LAMBERTIAN) then
        NA = 1
    else
        NA = NANG/2 + 1
    end if

    IT = 0
    IB = 0
    ! Loop over all points finding the top and bottom ones
    do I = 1, NPTS

        if (GRIDPOS(3, I) .ge. ZTOP) then
            IT = IT + 1
            if (IT .gt. MAXNBC .or. IT .gt. MAXBCRAD) then
                stop 'Error [_findBoundaryPnt]: MAXNBC exceeded - top.'
            end if
            BCPTR(IT, 1) = I
        end if

        if (GRIDPOS(3, I) .le. ZBOT) then
            IB = IB + 1
            if ((IB .gt. MAXNBC) .or. (IT + IB*NA .gt. MAXBCRAD)) then
                stop 'Error [_findBoundaryPnt]: MAXNBC exceeded - bottom.'
            end if
            BCPTR(IB, 2) = I
        end if

    end do

    NTOPPTS = IT
    NBOTPTS = IB

    return

end subroutine shd_grid__findBoundaryPnt

subroutine shd_grid__matchGridPnt(&
    &XP, YP, ZP, ICELL, IFACE,&
    &GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &IPMATCH)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Determines whether the neighboring cells to ICELL adjacent to
    !   face IFACE (1-6: -X,+X,-Y,+Y,-Z,+Z) have grid points that match
    !   the point at XP,YP,ZP on face IFACE of cell ICELL.
    !   If the face of ICELL has several neighboring cells, then the tree
    !   is traced down to find those cells that might contain XP,YP,ZP.
    !   IPMATCH is returned 0 if there is no match, otherwise it has the
    !   number of the matching existing gridpoint.
    !
    ! Input:
    !   XP        : real, total number of grid points
    !   GRIDPOS   : real array, (x,y,z) coordinate values of grid point locations
    !
    ! Output:
    !   IPMATCH   : integer, pointer array for NTOPPTS (top boundary points)
    !               and NBOTPTS (bottom boundary points)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(kind=R__)   , intent(in)  :: XP, YP, ZP
    integer(kind=I__), intent(in)  :: ICELL, IFACE
    real(kind=R__)   , intent(in)  :: GRIDPOS(3, *)
    integer(kind=I__), intent(in)  :: GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer(kind=2)  , intent(in)  :: CELLFLAGS(*)
    integer(kind=I__), intent(out) :: IPMATCH

    ! if MAXSTACK is greater than 128, kind=2 is needed
    integer(kind=1), parameter :: MAXSTACK=64
    integer(kind=I__) :: CELLSTACK(MAXSTACK), SP
    integer(kind=1)   :: OPPFACE(6), GRIDFACE(4, 6)
    data OPPFACE /2, 1, 4, 3, 6, 5/
    data GRIDFACE /1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8/

    integer(kind=I4_) :: INEXT, IDIR, DIR, IC, IC1, KFACE, IPT, I

    IPMATCH = 0
    IDIR = (IFACE + 1)/2
    KFACE = OPPFACE(IFACE)

    ! Get the pointer to the neighbor cell
    IC = abs(NEIGHPTR(IFACE, ICELL))
    ! Don't process if the neighbor cell is a boundary
    if (IC .eq. 0) return

    ! Go down all relevant branches of the tree until there are no children
    SP = 0
    do while (.true.)

        ! If at an end cell then check for a gridpoint match at the four
        !   corners of the cell's face
        do while (TREEPTR(2, IC) .eq. 0)
            do I = 1, 4
                IPT = GRIDPTR(GRIDFACE(I, KFACE), IC)
                if (XP .eq. GRIDPOS(1, IPT) .and. YP .eq. GRIDPOS(2, IPT)&
                    &.and. ZP .eq. GRIDPOS(3, IPT)) then
                    IPMATCH = IPT
                    return
                end if
            end do
            ! Return if there are no more branches on the stack, otherwise pop
            if (SP .eq. 0) return
            IC = CELLSTACK(SP)
            SP = SP - 1
        end do

        ! Get the direction current cell is split
        DIR = ibits(int(CELLFLAGS(IC)), 2, 2)

        if (DIR .eq. 0) stop 'Error [_matchGridPnt]: No split direction.'

        ! If split is in direction of face then choose the near child cell
        !   which is the negative side cell if this is a positive face.
        IC1 = TREEPTR(2, IC)
        if (DIR .eq. IDIR) then
            IC = IC1 + 1 - mod(IFACE - 1, 2)
        else
            !   If split is not in direction of face then find which of two cells
            ! has a face that contains the current location by comparing the
            ! appropriate coordinate (X,Y,Z) of the grid point on the split line
            ! with the current location.  If the split line matches a coordinate
            ! then we have to trace down both cells so put one on the stack.
            IC = IC1
            if (DIR .eq. 1) then
                if (XP .eq. GRIDPOS(1, GRIDPTR(8, IC1))) then
                    SP = SP + 1
                    CELLSTACK(SP) = IC + 1
                else if (XP .gt. GRIDPOS(1, GRIDPTR(8, IC1))) then
                    IC = IC + 1
                end if
            else if (DIR .eq. 2) then
                if (YP .eq. GRIDPOS(2, GRIDPTR(8, IC1))) then
                    SP = SP + 1
                    CELLSTACK(SP) = IC + 1
                else if (YP .gt. GRIDPOS(2, GRIDPTR(8, IC1))) then
                    IC = IC + 1
                end if
            else
                if (ZP .eq. GRIDPOS(3, GRIDPTR(8, IC1))) then
                    SP = SP + 1
                    CELLSTACK(SP) = IC + 1
                else if (ZP .gt. GRIDPOS(3, GRIDPTR(8, IC1))) then
                    IC = IC + 1
                end if
            end if
            if (SP .ge. MAXSTACK)&
                &stop 'Error [_matchGridPnt]: MAXSTACK exceeded.'
        end if
    end do

    return

end subroutine shd_grid__matchGridPnt

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
    ! Use the shd_grid__matchGridPnt routine to try to find matches in cells
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
        call shd_grid__matchGridPnt(XP, YP, ZP, ICELL, IFACE1, GRIDPOS,&
            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
        if (IPMATCH .eq. 0) then
            IFACE2 = FACEGRID(2, I, IDIR)
            call shd_grid__matchGridPnt(XP, YP, ZP, ICELL, IFACE2, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
        end if
        if (IPMATCH .eq. 0) then
            ICELL2 = abs(NEIGHPTR(IFACE1, ICELL))
            if (ICELL2 .gt. 0) then
                call shd_grid__matchGridPnt(XP, YP, ZP, ICELL2, IFACE2, GRIDPOS,&
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

subroutine shd_grid__makeSweepOrder(&
    &NX, NY, NZ, NPTS, NCELLS,&
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

subroutine shd_grid__divideCell(&
    &NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
    &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)

    !   Divides a single cell (ICELL) in the IDIR direction (1=Z,2=Y,3=X)
    ! Makes the pointers to the two new cells (TREEPTR), set the flags
    ! (CELLFLAGS), creates the new grid points if they don't exist already,
    ! and updates the neighbor pointers (NEIGHPTR) for the new cells and
    ! their neighboring cells.

    integer NPTS, NCELLS, ICELL, IDIR, NEWPOINTS(3, 4)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer NEWCELL, IFACE, I

    ! Only allow to divide an end cell (i.e. one not yet divided)
    if (TREEPTR(2, ICELL) .ne. 0) then
        stop 'Error [_divideCell]: Cannot divide already split cell.'
    end if

    ! Do the tree pointers: parent points to negative child,
    !   children point back to parent.
    NEWCELL = NCELLS + 1
    NCELLS = NCELLS + 2
    TREEPTR(2, ICELL) = NEWCELL
    TREEPTR(1, NEWCELL) = ICELL
    TREEPTR(2, NEWCELL) = 0
    TREEPTR(1, NEWCELL + 1) = ICELL
    TREEPTR(2, NEWCELL + 1) = 0
    ! Do the flags: put the splitting direction in the parent and
    !   inherit the boundary flags for the children.
    CELLFLAGS(ICELL) = ior(int(CELLFLAGS(ICELL)), ishft(IDIR, 2))
    CELLFLAGS(NEWCELL) = 0
    CELLFLAGS(NEWCELL + 1) = 0
    do I = 0, 1
        if (btest(int(CELLFLAGS(ICELL)), I)) then
            CELLFLAGS(NEWCELL) = ibset(int(CELLFLAGS(NEWCELL)), I)
            CELLFLAGS(NEWCELL + 1) = ibset(int(CELLFLAGS(NEWCELL + 1)), I)
        end if
    end do

    ! Do the grid points:
    call shd_grid__makeNewGridPnt(IDIR, ICELL, NEWCELL, NPTS,&
        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
        &GRIDPOS, NEWPOINTS)

    ! Do the neighbor pointers:
    do IFACE = 1, 6
        ! Special case if the neighbor is the cell itself
        if (NEIGHPTR(IFACE, ICELL) .eq. ICELL) then
            NEIGHPTR(IFACE, NEWCELL) = NEWCELL
            NEIGHPTR(IFACE, NEWCELL + 1) = NEWCELL + 1
        ! For the faces in the split direction, the children either
        !   point to the other child or inherit the parent's neighbor.
        !   Also update the neighbor's neighbor pointer.
        else if (IFACE .eq. 2*IDIR) then
            NEIGHPTR(IFACE, NEWCELL) = NEWCELL + 1
            NEIGHPTR(IFACE, NEWCELL + 1) = NEIGHPTR(IFACE, ICELL)
        else if (IFACE .eq. 2*IDIR - 1) then
            NEIGHPTR(IFACE, NEWCELL + 1) = NEWCELL
            NEIGHPTR(IFACE, NEWCELL) = NEIGHPTR(IFACE, ICELL)
        else
            ! For the other four faces:
            !   The children inherit the parent's neighbor pointers
            NEIGHPTR(IFACE, NEWCELL) = NEIGHPTR(IFACE, ICELL)
            NEIGHPTR(IFACE, NEWCELL + 1) = NEIGHPTR(IFACE, ICELL)
        end if
        ! Then see if we can match faces further down the neighbor's tree;
        !   also update the neighboring cells pointers
        call shd_grid__matchNeighCellFace(IFACE, NEWCELL, NEIGHPTR,&
            &CELLFLAGS, TREEPTR, GRIDPTR, GRIDPOS)
        call shd_grid__matchNeighCellFace(IFACE, NEWCELL + 1, NEIGHPTR,&
            &CELLFLAGS, TREEPTR, GRIDPTR, GRIDPOS)
    end do

    return

end subroutine shd_grid__divideCell

subroutine shd_grid__matchNeighCellFace(IFACE, IC, NEIGHPTR, CELLFLAGS,&
    &TREEPTR, GRIDPTR, GRIDPOS)

    !   Resolves the neighbor pointers for a new cell IC (from a subdivision)
    ! for the IFACE face.  Based on the size of this new face and
    ! the corresponding faces of the neighboring cells (down its tree),
    ! the new cell's NEIGHPTR and those of the neighboring cells are updated.

    integer IFACE, IC, NEIGHPTR(6, *), TREEPTR(2, *), GRIDPTR(8, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer INeigh, INN, JFACE, IC1, IC8, IN1, IN8, DIR, DIR1, DIR2
    logical DONE
    real POS(3)
    integer OPPFACE(6)
    data OPPFACE/2, 1, 4, 3, 6, 5/

    INeigh = abs(NEIGHPTR(IFACE, IC))
    ! If no neighbor then return
    if (INeigh .eq. 0) return
    JFACE = (IFACE + 1)/2
    ! Get the center of the new face
    IC1 = GRIDPTR(1, IC)
    IC8 = GRIDPTR(8, IC)
    DIR1 = mod(JFACE - 1 + 1, 3) + 1
    DIR2 = mod(JFACE - 1 + 2, 3) + 1
    POS(DIR1) = (GRIDPOS(DIR1, IC1) + GRIDPOS(DIR1, IC8))/2
    POS(DIR2) = (GRIDPOS(DIR2, IC1) + GRIDPOS(DIR2, IC8))/2
    ! Go down the neighbors cell tree
    DONE = .false.
    do while (.not. DONE .and. IC .ne. INeigh)
        ! Compare the coordinates of this neighbor cell to see
        !   if the new face is contained within the neighbor face.
        IN1 = GRIDPTR(1, INeigh)
        IN8 = GRIDPTR(8, INeigh)
        if (GRIDPOS(DIR1, IC1) .ge. GRIDPOS(DIR1, IN1) .and.&
            &GRIDPOS(DIR1, IC8) .le. GRIDPOS(DIR1, IN8) .and.&
            &GRIDPOS(DIR2, IC1) .ge. GRIDPOS(DIR2, IN1) .and.&
            &GRIDPOS(DIR2, IC8) .le. GRIDPOS(DIR2, IN8)) then
            ! If this is the bottom of the tree then we have a face match
            if (TREEPTR(2, INeigh) .eq. 0) then
                NEIGHPTR(IFACE, IC) = INeigh
            else
                NEIGHPTR(IFACE, IC) = -INeigh
            end if
        end if
        ! Compare the size of the neighbor face with the new face.
        !   Update the neighbor's neighbor pointer to the new face
        !   if the neighbor face is smaller, otherwise make the
        !   pointer negative, because it no longer points to an end node.
        if (GRIDPOS(DIR1, IN1) .ge. GRIDPOS(DIR1, IC1) .and.&
            &GRIDPOS(DIR1, IN8) .le. GRIDPOS(DIR1, IC8) .and.&
            &GRIDPOS(DIR2, IN1) .ge. GRIDPOS(DIR2, IC1) .and.&
            &GRIDPOS(DIR2, IN8) .le. GRIDPOS(DIR2, IC8)) then
            call shd_grid__inheritNeighCell(INeigh, OPPFACE(IFACE), IC,&
                &NEIGHPTR, TREEPTR, CELLFLAGS)
        else
            NEIGHPTR(OPPFACE(IFACE), INeigh) = -abs(NEIGHPTR(OPPFACE(IFACE), INeigh))
        end if

        ! Go down the neighbors cell tree until there is no children,
        !   finding the cell that contains the center of the new face.
        if (TREEPTR(2, INeigh) .eq. 0) then
            DONE = .true.
        else
            DIR = ibits(int(CELLFLAGS(INeigh)), 2, 2)
            INN = TREEPTR(2, INeigh)
            if (DIR .eq. JFACE) then
                INeigh = INN + 1 - mod(IFACE - 1, 2)
            else
                if (POS(DIR) .gt. GRIDPOS(DIR, GRIDPTR(8, INN))) then
                    INeigh = INN + 1
                else
                    INeigh = INN
                end if
            end if
        end if
    end do

    return

end subroutine shd_grid__matchNeighCellFace

subroutine shd_grid__inheritNeighCell(ICELL, IFACE, INeigh,&
    &NEIGHPTR, TREEPTR, CELLFLAGS)

    !   Updates the ICELL cell's new neighbor cell (INeigh) for all
    ! children cells that border the face (IFACE).

    integer ICELL, IFACE, INeigh, NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    integer JFACE, IC, DIR, SP, STACK(50)
    logical DONE

    ! Go through the tree checking the cells:
    ! The stack holds the second child cells, which are the
    !   next cell to go to if at an end cell.
    JFACE = (IFACE + 1)/2
    IC = ICELL
    SP = 0
    DONE = .false.
    do while (.not. DONE)
        NEIGHPTR(IFACE, IC) = INeigh
        ! If this is an end cell then pop the stack
        if (TREEPTR(2, IC) .eq. 0) then
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
            DIR = ibits(int(CELLFLAGS(IC)), 2, 2)
            if (DIR .eq. JFACE) then
                IC = TREEPTR(2, IC) + mod(IFACE - 1, 2)
            else
                ! Otherwise need to go to both children so put the second child
                !   on the stack and set the pointer to go to the first child next.
                SP = SP + 1
                if (SP .gt. 50) stop 'Error [_inheritNeighCell]: Stack exceeded.'
                STACK(SP) = TREEPTR(2, IC) + 1
                IC = TREEPTR(2, IC)
            end if
        end if
    end do

    return

end subroutine shd_grid__inheritNeighCell

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
    integer ID, IFACE, IP, IP1, IP2, INeigh, I, J, ISUM
    integer INCELL(2), DIRSPLIT(2)
    integer FACECORNER(4, 6), EDGECORNER(2, 3)
    real(RD_) :: XE, YE, ZE
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
                call shd_util__findNextCell(XE, YE, ZE, IFACE, ID, ICELL, GRIDPOS,&
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
                INeigh = NEIGHPTR(IFACE, ICELL)
                if (INeigh .lt. 0) then
                    DIRSPLIT(J) = ibits(int(CELLFLAGS(abs(INeigh))), 2, 2)
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
    real TOSPLIT, SUM1, JAY, TAU, EXT
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
            IDIR = ID
        end if
    end do

    return

end subroutine shd_grid__testCellSplit

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
    logical DONE
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

subroutine shd_grid__interpGrid(&
    &NPTS, NSTLEG, NLEG, GRIDPOS,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)

    ! Calls shd_grid__InterpTriLinGridPrp to interpolate the input arrays from
    !   the property grid to each internal grid point.

    integer NPTS, NSTLEG, NLEG
    integer*2 IPHASE(*)
    real GRIDPOS(3, NPTS)
    real TEMP(*), EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    integer IP

    ! Initialize: transfer the tabulated phase functions
    call shd_grid__InterpTriLinGridPrp(0.0, 0.0, 0.0, .true., NSTLEG, NLEG,&
        &TEMP(1), EXTINCT(1), ALBEDO(1), LEGEN(1, 0, 1), IPHASE(1))

    ! Trilinearly interpolate from the property grid to the adaptive
    do IP = 1, NPTS
        call shd_grid__InterpTriLinGridPrp(&
            &GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &.false., NSTLEG, NLEG, TEMP(IP), EXTINCT(IP),&
            &ALBEDO(IP), LEGEN(1, 0, IP), IPHASE(IP))
    end do

    return
end subroutine shd_grid__interpGrid

subroutine shd_grid__InterpTriLinGridPrp(X, Y, Z, INIT, NSTLEG, NLEG,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
    !  Trilinearly interpolates the quantities on the input property
    ! grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
    ! ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the
    ! volume coefficients.  Also adds in the separate gaseous absorption.
    ! Divides the phase function Legendre coefficients LEGEN by 2*l+1.
    ! The phase function pointer IPHASE (for tabulated phase functions)
    ! is that of the maximum weighted scattering property grid point.
    ! If INIT=.TRUE. then transfers the tabulated phase functions.

    use globals_prp

    integer NSTLEG, NLEG
    integer*2 IPHASE
    logical INIT
    real X, Y, Z, TEMP, EXTINCT, ALBEDO, LEGEN(NSTLEG, 0:NLEG, *)
    integer IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
    integer I1, I2, I3, I4, I5, I6, I7, I8, I
    real(RD_) :: U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
    real(RD_) :: SCAT1, SCAT2, SCAT3, SCAT4, SCAT5, SCAT6, SCAT7, SCAT8
    real(RD_) :: SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN
    save EXTMIN, SCATMIN

    integer :: ISELECT

    if (INIT) then
       ! If there are tabulated phase functions, then transfer them
        do I = 1, NUMPHASE_
            do L = 0, NLEG
                do J = 1, NSTLEG
                    LEGEN(J, L, I) = LEGENP_(J + NSTLEG*(L + (NLEG + 1)*(I - 1)))/(2*L + 1)
                end do
            end do
        end do
        !EXTMIN = 1.0d-5/((ZLEVELS_(NPZ_) - ZLEVELS_(1))/NPZ_)
        EXTMIN  = RACC_/((ZLEVELS_(NPZ_) - ZLEVELS_(1))/NPZ_) ! attempt to improve agreement between libRadtran in water vapor abosprtion bands
        SCATMIN = 0.1d0*EXTMIN
        return
    end if

    ! Find the grid location and compute the interpolation factors
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
    IZ = max(IL, 1)
    W = dble(Z - ZLEVELS_(IZ))/(ZLEVELS_(IZ + 1) - ZLEVELS_(IZ))
    W = max(min(W, 1.0d0), 0.0d0)
    IX = int((X - XSTART_)/DELX_) + 1
    if (abs(X - XSTART_ - NPX_*DELX_) .lt. 0.01*DELX_) IX = NPX_
    if (IX .lt. 1 .or. IX .gt. NPX_) then
        write (6, *) 'TRILIN: Beyond X domain'
        stop
    end if
    IXP = mod(IX, NPX_) + 1
    U = dble(X - XSTART_ - DELX_*(IX - 1))/DELX_
    U = max(min(U, 1.0d0), 0.0d0)
    if (U .lt. 1.0d-5) U = 0.0d0
    if (U .gt. 1.0d0 - 1.0d-5) U = 1.0d0
    IY = int((Y - YSTART_)/DELY_) + 1
    if (abs(Y - YSTART_ - NPY_*DELY_) .lt. 0.01*DELY_) IY = NPY_
    if (IY .lt. 1 .or. IY .gt. NPY_) then
        write (6, *) 'TRILIN: Beyond Y domain'
        stop
    end if
    IYP = mod(IY, NPY_) + 1
    V = dble(Y - YSTART_ - DELY_*(IY - 1))/DELY_
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
    I1 = IZ + NPZ_*( IY - 1) + NPZ_*NPY_*( IX - 1)
    I2 = IZ + NPZ_*( IY - 1) + NPZ_*NPY_*(IXP - 1)
    I3 = IZ + NPZ_*(IYP - 1) + NPZ_*NPY_*( IX - 1)
    I4 = IZ + NPZ_*(IYP - 1) + NPZ_*NPY_*(IXP - 1)

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

    SCAT1 = F1*EXTINCTP_(I1)*ALBEDOP_(I1)
    SCAT2 = F2*EXTINCTP_(I2)*ALBEDOP_(I2)
    SCAT3 = F3*EXTINCTP_(I3)*ALBEDOP_(I3)
    SCAT4 = F4*EXTINCTP_(I4)*ALBEDOP_(I4)
    SCAT5 = F5*EXTINCTP_(I5)*ALBEDOP_(I5)
    SCAT6 = F6*EXTINCTP_(I6)*ALBEDOP_(I6)
    SCAT7 = F7*EXTINCTP_(I7)*ALBEDOP_(I7)
    SCAT8 = F8*EXTINCTP_(I8)*ALBEDOP_(I8)

    MAXSCAT = -1.0
    if (SCAT1 .gt. MAXSCAT .or. abs(F1 - 1) .lt. 0.001) then
        MAXSCAT = SCAT1
        ISELECT = I1
    end if
    if (SCAT2 .gt. MAXSCAT .or. abs(F2 - 1) .lt. 0.001) then
        MAXSCAT = SCAT2
        ISELECT = I2
    end if
    if (SCAT3 .gt. MAXSCAT .or. abs(F3 - 1) .lt. 0.001) then
        MAXSCAT = SCAT3
        ISELECT = I3
    end if
    if (SCAT4 .gt. MAXSCAT .or. abs(F4 - 1) .lt. 0.001) then
        MAXSCAT = SCAT4
        ISELECT = I4
    end if
    if (SCAT5 .gt. MAXSCAT .or. abs(F5 - 1) .lt. 0.001) then
        MAXSCAT = SCAT5
        ISELECT = I5
    end if
    if (SCAT6 .gt. MAXSCAT .or. abs(F6 - 1) .lt. 0.001) then
        MAXSCAT = SCAT6
        ISELECT = I6
    end if
    if (SCAT7 .gt. MAXSCAT .or. abs(F7 - 1) .lt. 0.001) then
        MAXSCAT = SCAT7
        ISELECT = I7
    end if
    if (SCAT8 .gt. MAXSCAT .or. abs(F8 - 1) .lt. 0.001) then
        MAXSCAT = SCAT8
        ISELECT = I8
    end if

    !TEMP = TEMPP_(ISELECT)
    !EXTINCT = EXTINCTP_(ISELECT)
    !SCATTER = EXTINCTP_(ISELECT)*ALBEDOP_(ISELECT)

    TEMP = F1*TEMPP_(I1) + F2*TEMPP_(I2) + F3*TEMPP_(I3) + F4*TEMPP_(I4)&
        &+ F5*TEMPP_(I5) + F6*TEMPP_(I6) + F7*TEMPP_(I7) + F8*TEMPP_(I8)
    EXTINCT = F1*EXTINCTP_(I1) + F2*EXTINCTP_(I2) + F3*EXTINCTP_(I3)&
        &+ F4*EXTINCTP_(I4) + F5*EXTINCTP_(I5) + F6*EXTINCTP_(I6)&
        &+ F7*EXTINCTP_(I7) + F8*EXTINCTP_(I8)
    SCATTER = SCAT1 + SCAT2 + SCAT3 + SCAT4 + SCAT5 + SCAT6 + SCAT7 + SCAT8

    if (EXTINCT .gt. EXTMIN) then
        ALBEDO = SCATTER/EXTINCT
    else
        ALBEDO = SCATTER/EXTMIN
    end if

    if (NUMPHASE_ .gt. 0) then
        IPHASE = IPHASEP_(ISELECT)
    else
        ! Standard property file format: unpolarized
        LEGEN(1, 0, 1) = 1.0
        do L = 1, NLEG
            LEGEN(1, L, 1) = SCAT1*LEGENP_(L + (NLEG + 1)*(I1 - 1))&
                &+ SCAT2*LEGENP_(L + (NLEG + 1)*(I2 - 1))&
                &+ SCAT3*LEGENP_(L + (NLEG + 1)*(I3 - 1))&
                &+ SCAT4*LEGENP_(L + (NLEG + 1)*(I4 - 1))&
                &+ SCAT5*LEGENP_(L + (NLEG + 1)*(I5 - 1))&
                &+ SCAT6*LEGENP_(L + (NLEG + 1)*(I6 - 1))&
                &+ SCAT7*LEGENP_(L + (NLEG + 1)*(I7 - 1))&
                &+ SCAT8*LEGENP_(L + (NLEG + 1)*(I8 - 1))
            if (SCATTER .gt. SCATMIN) then
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATTER
            else
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATMIN
            end if
            LEGEN(1, L, 1) = LEGEN(1, L, 1)/(2*L + 1)
        end do
    end if

    ! Trilinearly interpolate the temperature, extinction, scattering
    ! original
    !╭────────────────────────────────────────────────╮!
    !TEMP = F1*TEMPP_(I1) + F2*TEMPP_(I2) + F3*TEMPP_(I3) + F4*TEMPP_(I4)&
    !    &+ F5*TEMPP_(I5) + F6*TEMPP_(I6) + F7*TEMPP_(I7) + F8*TEMPP_(I8)
    !EXTINCT = F1*EXTINCTP_(I1) + F2*EXTINCTP_(I2) + F3*EXTINCTP_(I3)&
    !    &+ F4*EXTINCTP_(I4) + F5*EXTINCTP_(I5) + F6*EXTINCTP_(I6)&
    !    &+ F7*EXTINCTP_(I7) + F8*EXTINCTP_(I8)
    !SCAT1 = F1*EXTINCTP_(I1)*ALBEDOP_(I1)
    !SCAT2 = F2*EXTINCTP_(I2)*ALBEDOP_(I2)
    !SCAT3 = F3*EXTINCTP_(I3)*ALBEDOP_(I3)
    !SCAT4 = F4*EXTINCTP_(I4)*ALBEDOP_(I4)
    !SCAT5 = F5*EXTINCTP_(I5)*ALBEDOP_(I5)
    !SCAT6 = F6*EXTINCTP_(I6)*ALBEDOP_(I6)
    !SCAT7 = F7*EXTINCTP_(I7)*ALBEDOP_(I7)
    !SCAT8 = F8*EXTINCTP_(I8)*ALBEDOP_(I8)
    !SCATTER = SCAT1 + SCAT2 + SCAT3 + SCAT4 + SCAT5 + SCAT6 + SCAT7 + SCAT8
    !if (EXTINCT .gt. EXTMIN) then
    !    ALBEDO = SCATTER/EXTINCT
    !else
    !    ALBEDO = SCATTER/EXTMIN
    !end if
    !╰────────────────────────────────────────────────╯!

    ! For tabulated phase functions pick the one we are on top of
    !   or the one with the most scattering weight.
    ! original
    !╭────────────────────────────────────────────────╮!
    !if (NUMPHASE_ .gt. 0) then
    !    MAXSCAT = -1.0
    !    if (SCAT1 .gt. MAXSCAT .or. abs(F1 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT1
    !        IPHASE = IPHASEP_(I1)
    !    end if
    !    if (SCAT2 .gt. MAXSCAT .or. abs(F2 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT2
    !        IPHASE = IPHASEP_(I2)
    !    end if
    !    if (SCAT3 .gt. MAXSCAT .or. abs(F3 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT3
    !        IPHASE = IPHASEP_(I3)
    !    end if
    !    if (SCAT4 .gt. MAXSCAT .or. abs(F4 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT4
    !        IPHASE = IPHASEP_(I4)
    !    end if
    !    if (SCAT5 .gt. MAXSCAT .or. abs(F5 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT5
    !        IPHASE = IPHASEP_(I5)
    !    end if
    !    if (SCAT6 .gt. MAXSCAT .or. abs(F6 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT6
    !        IPHASE = IPHASEP_(I6)
    !    end if
    !    if (SCAT7 .gt. MAXSCAT .or. abs(F7 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT7
    !        IPHASE = IPHASEP_(I7)
    !    end if
    !    if (SCAT8 .gt. MAXSCAT .or. abs(F8 - 1) .lt. 0.001) then
    !        MAXSCAT = SCAT8
    !        IPHASE = IPHASEP_(I8)
    !    end if
    !else
    !    ! Standard property file format: unpolarized
    !    LEGEN(1, 0, 1) = 1.0
    !    do L = 1, NLEG
    !        LEGEN(1, L, 1) = SCAT1*LEGENP_(L + (NLEG + 1)*(I1 - 1))&
    !            &+ SCAT2*LEGENP_(L + (NLEG + 1)*(I2 - 1))&
    !            &+ SCAT3*LEGENP_(L + (NLEG + 1)*(I3 - 1))&
    !            &+ SCAT4*LEGENP_(L + (NLEG + 1)*(I4 - 1))&
    !            &+ SCAT5*LEGENP_(L + (NLEG + 1)*(I5 - 1))&
    !            &+ SCAT6*LEGENP_(L + (NLEG + 1)*(I6 - 1))&
    !            &+ SCAT7*LEGENP_(L + (NLEG + 1)*(I7 - 1))&
    !            &+ SCAT8*LEGENP_(L + (NLEG + 1)*(I8 - 1))
    !        if (SCATTER .gt. SCATMIN) then
    !            LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATTER
    !        else
    !            LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATMIN
    !        end if
    !        LEGEN(1, L, 1) = LEGEN(1, L, 1)/(2*L + 1)
    !    end do
    !end if
    !╰────────────────────────────────────────────────╯!

    ! Add in the gaseous absorption to extinction and albedo
    if (NZCKD_ .gt. 0) then
        IL = 1
        IU = NZCKD_
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (Z .le. ZCKD_(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
        I = min(max(IL, 1), NZCKD_ - 1)
        F = (Z - ZCKD_(I))/(ZCKD_(I + 1) - ZCKD_(I))
        F = min(max(F, 0.0d0), 1.0d0)
        KG = (1.0 - F)*GASABS_(I) + F*GASABS_(I + 1)
        if (EXTINCT + KG .gt. 0.0) then
            ALBEDO = ALBEDO*EXTINCT/(EXTINCT + KG)
        else
            ALBEDO = 0.0
        end if
        EXTINCT = EXTINCT + KG
    end if

    return
end subroutine shd_grid__InterpTriLinGridPrp

subroutine shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
    &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    &UNITS, WAVENO, WAVELEN,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)

    !   Interpolates the medium properties, radiance, and source function for
    ! new grid points (specified in NEWPOINTS).  The medium properties are
    ! interpolated from the property grid, since interpolating from the
    ! internal grid may not faithfully represent the true medium.
    ! The delta-M scaling is done to the new grid point properties if needed.
    ! The direct beam flux and the Planck function for the new points is
    ! recomputed (rather than interpolated).  The radiance is averaged
    ! from the two parent grid points, and used to compute the source
    ! function for the new points.

    integer NEWPOINTS(3, 4), NPTS
    integer NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer BCFLAG, IPFLAG
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    integer*2 IPHASE(NPTS)
    logical DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(NSTLEG, NLM)
    real WAVENO(2), WAVELEN
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real TEMP(NPTS), PLANCK(NPTS), DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *), SOURCE(NSTOKES, *)
    real GRIDPOS(3, *)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IP1, IP2, IP, IPH, J, K, L, M, ME, MS
    integer IS, NS, NS1, NS2
    integer IR, IR1, IR2, NR, NR1, NR2, SIDE
    logical VALIDBEAM
    real X, Y, Z, F, BB, C, D, SECMU0
    real RAD1(NSTOKES), RAD2(NSTOKES)
    real UNIFZLEV, XO, YO, ZO, DIRPATH
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    C = sqrt(4.0*acos(-1.0))
    if (SRCTYPE .ne. 'T') then
        SECMU0 = 1.0/abs(SOLARMU)
    end if

    ! Loop over the (potentially) four new grid points
    do I = 1, 4
        if (NEWPOINTS(3, I) .gt. 0) then
            ! Get the parent and new grid pointers
            IP1 = NEWPOINTS(1, I)
            IP2 = NEWPOINTS(2, I)
            IP = NEWPOINTS(3, I)

            ! Interpolate the medium properties from the property grid
            X = GRIDPOS(1, IP)
            Y = GRIDPOS(2, IP)
            Z = GRIDPOS(3, IP)
            call shd_grid__InterpTriLinGridPrp(X, Y, Z, .false., NSTLEG, NLEG,&
                &TEMP(IP), EXTINCT(IP), ALBEDO(IP),&
                &LEGEN(1, 0, IP), IPHASE(IP))
            ! Do the Delta-M scaling of extinction and albedo for this point
            if (DELTAM) then
                if (NUMPHASE .gt. 0) then
                    F = LEGEN(1, ML + 1, IPHASE(IP))
                else
                    F = LEGEN(1, ML + 1, IP)
                    do L = 0, ML
                        LEGEN(1, L, IP) = (LEGEN(1, L, IP) - F)/(1 - F)
                    end do
                end if
                EXTINCT(IP) = (1.0 - ALBEDO(IP)*F)*EXTINCT(IP)
                ALBEDO(IP) = (1.0 - F)*ALBEDO(IP)/(1.0 - ALBEDO(IP)*F)
            end if

            ! Compute the new Planck source function (if needed)
            if (SRCTYPE .ne. 'S') then
                call shd_util__computePlanck(TEMP(IP), UNITS, WAVENO, WAVELEN, BB)
                PLANCK(IP) = (1.0 - ALBEDO(IP))*BB
            end if

            ! Compute the new direct beam flux (if needed)
            if (SRCTYPE .ne. 'T') then
                if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
                    ! If parallel mode then calculate a temporary approximation
                    !   to the direct flux (for source function calculation below)
                    DIRFLUX(IP) = exp(0.5*(log(DIRFLUX(IP1)) + log(DIRFLUX(IP2))))
                else
                    ! Otherwise, calculate the exact direct beam from property grid
                    DIRPATH = 0.0
                    call shd_rad__computeDirect(0, X, Y, Z, BCFLAG, IPFLAG, DELTAM,&
                        &ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                        &DIRFLUX(IP),&
                        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                end if
            end if

            ! Interpolate the radiance
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
                ! this is modified to address the layer averaging issue
                !╭────────────────────────────────────────────────╮!
                !if ((abs(GRIDPOS(1, IP1)-GRIDPOS(1, IP2)) .lt. 1.0e-6) .and.&
                !    &(abs(GRIDPOS(2, IP1)-GRIDPOS(2, IP2)) .lt. 1.0e-6)) then
                !    if (GRIDPOS(3, IP2) .gt. GRIDPOS(3, IP1)) then
                !        RADIANCE(:, IR + J) = RAD2(:)
                !    else
                !        RADIANCE(:, IR + J) = RAD1(:)
                !    end if
                !else
                    !RADIANCE(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
                !end if
                !╰────────────────────────────────────────────────╯!
                RADIANCE(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
            end do

            ! Compute the source function for the new points
            NS1 = SHPTR(IP1 + 1) - SHPTR(IP1)
            NS2 = SHPTR(IP2 + 1) - SHPTR(IP2)
            NS = max(NS1, NS2)
            IS = SHPTR(IP)
            SHPTR(IP + 1) = IS + NS
            if (ACCELFLAG) OSHPTR(IP + 1) = OSHPTR(IP)
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(IP)
            else
                IPH = IP
            end if
            if (NSTOKES .eq. 1) then
                call shd_src__computeSrcOnePnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    &PLANCK(IP), ALBEDO(IP), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_src__computeSrcOnePnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    &PLANCK(IP), ALBEDO(IP), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if
            do J = 1, NS
                SOURCE(:, IS + J) = SOURCET(:, J)
            end do
        end if
    end do
    deallocate (LOFJ, SOURCET)

    return

end subroutine shd_grid__interpNewGridPnt

subroutine shd_grid__splitGrid(&
    &MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
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
            call shd_util__ssort(ADAPTCRIT, ADAPTIND, N, -2)

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
                    call shd_grid__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                    ! Interpolate the medium properties, radiance, and source
                    call shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
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
                        call shd_grid__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                            &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                        ! Interpolate the medium properties, radiance, and source
                        call shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
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
            !write (6, '(2A,I8,A,I8)')&
            !    &'Message [_splitGrid]: Further cell splitting requires larger array allocations: ',&
            !    &'Ncells = ', NCELLS, ', Npnts = ', NPTS
        end if

        ! If in parallel mode, calculate all the direct beam fluxes at once
        if (SRCTYPE .ne. 'T') then
            if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
                call shd_mpi__makeDirect(OLDNPTS + 1, NPTS, BCFLAG,&
                    &IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
                    &NX, XGRID, NY, YGRID, DIRFLUX)
            end if
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

end module shd_grid
