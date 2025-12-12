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

    use globals_const, only: RD_, PI_, OPI_, DTOR_
    use globals_cfg,   only: verbose_, doGridPtsPrp_, &
                             planetRadius_, A2_, B2_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, &
                             is1DX_, is1DYx_

    use globals_type,  only: shd_type_rtm, &
                             shd_type_cell, &
                             shd_type_prp, &
                             shd_type_sfc

    use globals_obj, only:  rtm_obj_, xyz_Rtm_, &
                           shdo_obj_, &
                            src_obj_, &
                            prp_obj_, &
                            sfc_obj_, brdf_Rtm_

    use shd_util, only: shd_util__computePlanck, &
                        shd_util__computeSolarAngle, &
                        shd_util__solveTriDiag

    use shd_sfc,  only: shd_sfc__brdfSurface


    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_init__makeSolarVect

    public :: shd_init__makeGrid
    public :: shd_init__makeCell

    public :: shd_init__makePrp
    public :: shd_init__makeRad
    public :: shd_init__makeFlux

    public :: shd_init__makeTopRad
    public :: shd_init__makeBotRadUL
    public :: shd_init__makeBotRadVL
    public :: shd_init__makeBotRadU
    public :: shd_init__makeBotRadV

    public :: shd_init__checkSfcU

    real(kind=R__) :: lon0, lat0, mu0, phi0

contains

subroutine shd_init__makeSolarVect(&
    lon_deg, lat_deg, mu0, phi0, &
    solarVect)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Initialize Sun vector at center, converting local ENU to ECEF
    !   mu0 is the SHDOM muSun
    !   phi0 is the SHDOM phiSun
    !╰────────────────────────────────────────────────────────────────────────────╯!

    implicit none
    real(kind=R__), intent(in)  :: lon_deg, lat_deg, mu0, phi0
    real(kind=R__), intent(out) :: solarVect(3)
    real(kind=R__) :: mu, phi, lon, lat, Se, Sn, sinmu
    real(kind=R__) :: east_x, east_y
    real(kind=R__) :: up_x, up_y, up_z
    real(kind=R__) :: north_x, north_y, north_z
    real(kind=R__) :: nx, ny, nz, norm

    ! convert MU (cos(SZA)) and PHI from SHDOM convension to legacy convension
    mu  = abs(mu0)
    phi = 1.5*PI_ - phi0

    ! convert longitude and latitude from degree to radian
    lon = lon_deg * DTOR_
    lat = lat_deg * DTOR_

    ! Local ENU Sun vector at center
    sinmu = sqrt(1.0 - mu*mu)

    Se = sinmu * sin(phi)
    Sn = sinmu * cos(phi)

    ! Ellipsoidal surface normal (UP vector)
    nx = cos(lat)*cos(lon)/A2_
    ny = cos(lat)*sin(lon)/A2_
    nz = sin(lat)/B2_
    norm = sqrt(nx*nx + ny*ny + nz*nz)

    up_x = nx / norm
    up_y = ny / norm
    up_z = nz / norm

    ! East vector (same as spherical)
    east_x = -sin(lon)
    east_y =  cos(lon)

    ! North = Up × East
    north_x = -up_z*east_y
    north_y =  up_z*east_x
    north_z =  up_x*east_y - up_y*east_x

    ! Convert local ENU Sun vector -> global ECEF
    solarVect(1) = Se*east_x  + Sn*north_x + mu*up_x
    solarVect(2) = Se*east_y  + Sn*north_y + mu*up_y
    solarVect(3) = Sn*north_z + mu*up_z

    return

end subroutine shd_init__makeSolarVect

subroutine shd_init__makeGrid(rtm_OBJ, &
                              prp_OBJ)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Makes the rtm_obj%x, rtm_obj%y, and rtm_obj%z arrays from the input Z levels
    !   and X and Y spacing from the property file according to the rtm_obj%type.
    !
    ! Input:
    !   rtm_obj%type: character(len=1), grid type
    !             "E": even Z base grid between bottom and top;
    !             "P": Z base grid levels from property file;
    !             "F": Z base grid levels from file: shdom-zgrid.txt
    !   prp_obj%Nx, prp_obj%Ny, prp_obj%Nz: integer, size of the input property arrays
    !   rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz: integer, base grid size in X, Y and Z
    !   prp_obj%xStart, prp_obj%yStart: real, starting X and Y
    !   prp_obj%dx, prp_obj%dy: real, horizontal grid spacing of the property arrays
    !   prp_obj%z: real array, input Z levels from the property file
    !
    ! Output:
    !   rtm_obj%x, rtm_obj%y, rtm_obj%z: real array, grid locations for each dimension.
    !                        To facilitate the periodic horizontal boundaries,
    !                        there are rtm_obj%Nx+1 values in rtm_obj%x, rtm_obj%Ny+1 in rtm_obj%y,
    !                        rtm_obj%Nz values in rtm_obj%z. Currently the grid locations
    !                        are evenly spaced in X and Y.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    type(shd_type_rtm),  intent(inout) :: rtm_obj
    type(shd_type_prp) , intent(inout) :: prp_obj

    integer(kind=I__) :: ix, iy, iz

    ! Make the X, Y, and Z grids, depending on type:
    !   E for evenly spaced, P from property file, F from shdom-zgrid.txt
    !   Only Even allowed horizontally

    ! x
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (prp_obj%dx .le. 0.0) then
        prp_obj%dx = 1.0
        if (verbose_ .ge. 4) write(*, '(A)') 'Warning [_makeGrid]: Found prp_obj%dx<=0.0, reset prp_obj%dx to 1.0.'
    end if

    if (isMpiX_) then ! open in X
        do ix = 1, rtm_obj%Nx
            if (prp_obj%Nx .gt. 1) then
                rtm_obj%x(ix) = prp_obj%xStart + (ix - 1)*(prp_obj%dx*(prp_obj%Nx - 1))/float(rtm_obj%Nx - 1)
            else
                rtm_obj%x(ix) = prp_obj%xStart + (ix - 1)*(prp_obj%dx*prp_obj%Nx)/float(rtm_obj%Nx - 1)
            end if
        end do
    else ! periodic in X
        do ix = 1, rtm_obj%Nx + 1
            rtm_obj%x(ix) = prp_obj%xStart + (ix - 1)*(prp_obj%dx*prp_obj%Nx)/float(rtm_obj%Nx)
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! y
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (prp_obj%dy .le. 0.0) then
        prp_obj%dy = 1.0
        if (verbose_ .ge. 4) write(*, '(A)') 'Warning [_makeGrid]: Found prp_obj%dy<=0.0, reset prp_obj%dy to 1.0.'
    end if

    if (isMpiY_) then ! open in Y
        do iy = 1, rtm_obj%Ny
            if (prp_obj%Ny .gt. 1) then
                rtm_obj%y(iy) = prp_obj%yStart + (iy - 1)*(prp_obj%dy*(prp_obj%Ny - 1))/float(rtm_obj%Ny - 1)
            else
                rtm_obj%y(iy) = prp_obj%yStart + (iy - 1)*(prp_obj%dy*prp_obj%Ny)/float(rtm_obj%Ny - 1)
            end if
        end do
    else ! periodic in Y
        do iy = 1, rtm_obj%Ny + 1
            rtm_obj%y(iy) = prp_obj%yStart + (iy - 1)*(prp_obj%dy*prp_obj%Ny)/float(rtm_obj%Ny)
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! z
    !╭────────────────────────────────────────────────────────────────────────────╮!
    GRIDTYPE_CASE: select case (rtm_obj%type)

        case ('P', 'p', 'S', 's')

            if (rtm_obj%Nz .ne. prp_obj%Nz) then
                write(6, *)&
                    &'Error [_makeGrid]: For internal Z grid same as property file grid,'
                write(6, *) '  must have same number of grid points (rtm_obj%Nz=prp_obj%Nz).'
                stop
            end if
            do iz = 1, rtm_obj%Nz
                rtm_obj%z(iz) = prp_obj%z(iz)
            end do

        case ('E', 'e')

            do iz = 1, rtm_obj%Nz
                rtm_obj%z(iz) = prp_obj%z(1)&
                    + (iz - 1)*(prp_obj%z(prp_obj%Nz) - prp_obj%z(1))/float(rtm_obj%Nz - 1)
            end do

        case ('F', 'f')

            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_makeGrid]: Reading rtm_obj%z from <shdom-zgrid.txt> ...'
            open(UNIT=3, FILE='shdom-zgrid.txt', STATUS='OLD')
            do iz = 1, rtm_obj%Nz
                read(3, *) rtm_obj%z(iz)
            end do
            close(3)

        case default

            stop 'Error [_makeGrid]: Illegal grid type in Z.'

    end select GRIDTYPE_CASE
    !╰────────────────────────────────────────────────────────────────────────────╯!

    return

end subroutine shd_init__makeGrid

subroutine shd_init__makeCell(rtm_obj, &
                              cell_obj, &
                              xyz_Rtm, i8Grid_Cell, i6Cell_Cell, i2Tree_Cell, flag_Cell)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Make the initial grid cell data structure for the base grid.
    !   For periodic boundary conditions there are rtm_obj%Nx*rtm_obj%Ny*(rtm_obj%Nz-1) grid cells
    !   and (rtm_obj%Nx+1)*(rtm_obj%Ny+1)*rtm_obj%Nz grid points. The ix=rtm_obj%Nx+1 column is the same
    !   as the ix=1 column, because of the periodicity, but the grid points
    !   are repeated.  The last grid cell on the "right" has the last column
    !   of grid points on the "left" and the repeated "virtual" set of grid
    !   points on the "right". For open boundary conditions in X and Y
    !   there are (rtm_obj%Nx+1)*(rtm_obj%Ny+1)*(rtm_obj%Nz-1) grid cells and rtm_obj%Nx*rtm_obj%Ny*rtm_obj%Nz grid points.
    !   Open boundary conditions are done with the lowest "left most" and
    !   "right most" X and Y columns being independent pixels with
    !   their own grid cells (zero width), therefore one extra grid cell
    !   is required.   Multiple processor boundary conditions have
    !   (rtm_obj%Nx-1)*(rtm_obj%Ny-1)*(rtm_obj%Nz-1) grid cells and rtm_obj%Nx*rtm_obj%Ny*rtm_obj%Nz grid points.
    !   Independent pixel mode has rtm_obj%Nx*rtm_obj%Ny*(rtm_obj%Nz-1) grid cells and rtm_obj%Nx*rtm_obj%Ny*rtm_obj%Nz
    !   grid points
    !
    ! Input:
    !   rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz: integer, base grid size in X, Y and Z
    !   rtm_obj%NxA, rtm_obj%NyA: integer, base grid actual size in X, Y
    !   rtm_obj%x, rtm_obj%y, rtm_obj%z: real array, grid locations for each dimension.
    !                        To facilitate the periodic horizontal boundaries,
    !                        there are rtm_obj%Nx+1 values in rtm_obj%x, rtm_obj%Ny+1 in rtm_obj%y,
    !                        rtm_obj%Nz values in rtm_obj%z. Currently the grid locations
    !                        are evenly spaced in X and Y.
    !
    ! Output:
    !   rtm_obj%Npts : total number of grid points
    !   cell_obj%Npts: total number of cells
    !   xyz_Rtm(3,IP) : X,Y,Z coordinate values of grid point locations
    !   i8Grid_Cell(8,IC) : pointer to 8 grid points belonging to a cell.
    !                   Order of pointers is based on binary digits Z-Y-X:
    !                   iOctant = 1 + bitX + 2*bitY + 4*bitZ, where
    !                   bit is 0 for negative side and 1 for positive side.
    !                   1: (-X, -Y, -Z), 2: (+X, -Y, -Z), 3: (-X, +Y, -Z), 4: (+X, +Y, -Z)
    !                   5: (-X, -Y, +Z), 6: (+X, -Y, +Z), 7: (-X, +Y, +Z), 8: (+X, +Y, +Z)
    !   i6Cell_Cell(6,IC): pointer to 6 neighboring cells of a cell.
    !                   Order of pointers are 1:-X, 2:+X, 3:-Y, 4:+Y, 5:-Z, 6:+Z.
    !                   Is 0 if there is no neighbor cell (boundary).
    !                   Is positive if there is just one cell across face
    !                   of current cell (unique neighbor).  Is negative
    !                   (-index) if there is more than one cell across
    !                   the face of current cell.  In this case the
    !                   cell pointed to is the parent cell that has a
    !                   has same or larger face adjoining current cell.
    !   i2Tree_Cell(2,IC) : pointer to parent and children cells in tree structure.
    !                   First is pointer to parent cell, 0 for a base cell.
    !                   Second is pointer to first of the two children cells.
    !                   This pointer is 0 if there are no children cells.
    !                   The children cells are in successive order, with
    !                   the negative side first. See flag_Cell for direction
    !                   of division of split.
    !   flag_Cell(IC) : a collection of binary flags for the cells:
    !                   bit  For bit set
    !                   0:  Do independent pixel for X
    !                   1:  Do independent pixel for Y
    !                   2,3: Direction cell is split (0=no split, 1=X, 2=Y, 3=Z)
    ! IP is a grid point index, IC is a cell index.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    type(shd_type_rtm),  intent(inout) :: rtm_obj
    type(shd_type_cell), intent(inout) :: cell_obj

    real(kind=R__),    intent(out) :: xyz_Rtm(3, *)
    integer(kind=I__), intent(out) :: i8Grid_Cell(8, *)
    integer(kind=I__), intent(out) :: i6Cell_Cell(6, *)
    integer(kind=I__), intent(out) :: i2Tree_Cell(2, *)
    integer(kind=2),   intent(out) :: flag_Cell(*)

    integer(kind=I__) :: ix, iy, iz, ix0, ix1, iy0, iy1, i, Nxc, Nyc

    ! First, make the grid point positions
    i = 0
    do ix = 1, rtm_obj%NxA
        do iy = 1, rtm_obj%NyA
            do iz = 1, rtm_obj%Nz
                i = i + 1
                xyz_Rtm(1, i) = rtm_obj%x(ix)
                xyz_Rtm(2, i) = rtm_obj%y(iy)
                xyz_Rtm(3, i) = rtm_obj%z(iz)
            end do
        end do
    end do
    rtm_obj%Npts = i

    Nxc = rtm_obj%Nx
    if (isOpenX_) Nxc = rtm_obj%Nx + 1
    if (isMpiX_ .and. (.not. is1DX_)) Nxc = rtm_obj%Nx - 1
    Nyc = rtm_obj%Ny
    if (isOpenYx_) Nyc = rtm_obj%Ny + 1
    if (isMpiY_ .and. (.not. is1DYx_)) Nyc = rtm_obj%Ny - 1

    ! Then setup the cell structures
    i = 0
    do ix = 1, Nxc
        if (is1DX_) then
            ix0 = ix
            ix1 = ix
        else if (isOpenX_) then
            ix0 = max(1, ix - 1)
            ix1 = min(rtm_obj%Nx, ix)
        else
            ix0 = ix
            ix1 = ix + 1
        end if
        do iy = 1, Nyc
            if (is1DYx_) then
                iy0 = iy
                iy1 = iy
            else if (isOpenYx_) then
                iy0 = max(1, iy - 1)
                iy1 = min(rtm_obj%Ny, iy)
            else
                iy0 = iy
                iy1 = iy + 1
            end if
            do iz = 1, rtm_obj%Nz - 1
                i = i + 1
                i8Grid_Cell(1, i) = iz +     rtm_obj%Nz*(iy0 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix0 - 1)
                i8Grid_Cell(2, i) = iz +     rtm_obj%Nz*(iy0 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix1 - 1)
                i8Grid_Cell(3, i) = iz +     rtm_obj%Nz*(iy1 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix0 - 1)
                i8Grid_Cell(4, i) = iz +     rtm_obj%Nz*(iy1 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix1 - 1)
                i8Grid_Cell(5, i) = iz + 1 + rtm_obj%Nz*(iy0 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix0 - 1)
                i8Grid_Cell(6, i) = iz + 1 + rtm_obj%Nz*(iy0 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix1 - 1)
                i8Grid_Cell(7, i) = iz + 1 + rtm_obj%Nz*(iy1 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix0 - 1)
                i8Grid_Cell(8, i) = iz + 1 + rtm_obj%Nz*(iy1 - 1) + rtm_obj%Nz*rtm_obj%NyA*(ix1 - 1)
                flag_Cell(I) = 0
                ! Initialize the neighbor pointer; to self for IP mode or
                !   if open boundary conditions and at either end
                if (is1DX_) then
                    i6Cell_Cell(1, i) = i
                    i6Cell_Cell(2, i) = i
                    flag_Cell(I) = ibset(int(flag_Cell(I)), 0)
                else if (isOpenX_) then
                    if (ix .eq. 1) then
                        i6Cell_Cell(1, i) = i
                        i6Cell_Cell(2, i) = i + (rtm_obj%Nz - 1)*Nyc
                        flag_Cell(I) = ibset(int(flag_Cell(I)), 0)
                    else if (ix .eq. Nxc) then
                        i6Cell_Cell(1, i) = i - (rtm_obj%Nz - 1)*Nyc
                        i6Cell_Cell(2, i) = i
                        flag_Cell(I) = ibset(int(flag_Cell(I)), 0)
                    else
                        i6Cell_Cell(1, i) = i - (rtm_obj%Nz - 1)*Nyc
                        i6Cell_Cell(2, i) = i + (rtm_obj%Nz - 1)*Nyc
                    end if
                else if (isMpiX_) then
                    if (ix .eq. 1) then
                        i6Cell_Cell(1, i) = 0
                        i6Cell_Cell(2, i) = i + (rtm_obj%Nz - 1)*Nyc
                    else if (ix .eq. Nxc) then
                        i6Cell_Cell(1, i) = i - (rtm_obj%Nz - 1)*Nyc
                        i6Cell_Cell(2, i) = 0
                    else
                        i6Cell_Cell(1, i) = i - (rtm_obj%Nz - 1)*Nyc
                        i6Cell_Cell(2, i) = i + (rtm_obj%Nz - 1)*Nyc
                    end if
                else
                    if (ix .eq. 1) then
                        i6Cell_Cell(1, i) = i + (Nxc - 1)*(rtm_obj%Nz - 1)*Nyc
                    else
                        i6Cell_Cell(1, i) = i - (rtm_obj%Nz - 1)*Nyc
                    end if
                    if (ix .eq. rtm_obj%Nx) then
                        i6Cell_Cell(2, i) = i - (Nxc - 1)*(rtm_obj%Nz - 1)*Nyc
                    else
                        i6Cell_Cell(2, i) = i + (rtm_obj%Nz - 1)*Nyc
                    end if
                end if
                if (is1DYx_) then
                    i6Cell_Cell(3, i) = i
                    i6Cell_Cell(4, i) = i
                    flag_Cell(I) = ibset(int(flag_Cell(I)), 1)
                else if (isOpenYx_) then
                    if (iy .eq. 1) then
                        i6Cell_Cell(3, i) = i
                        i6Cell_Cell(4, i) = i + (rtm_obj%Nz - 1)
                        flag_Cell(I) = ibset(int(flag_Cell(I)), 1)
                    else if (iy .eq. Nyc) then
                        i6Cell_Cell(3, i) = i - (rtm_obj%Nz - 1)
                        i6Cell_Cell(4, i) = i
                        flag_Cell(I) = ibset(int(flag_Cell(I)), 1)
                    else
                        i6Cell_Cell(3, i) = i - (rtm_obj%Nz - 1)
                        i6Cell_Cell(4, i) = i + (rtm_obj%Nz - 1)
                    end if
                else if (isMpiY_) then
                    if (iy .eq. 1) then
                        i6Cell_Cell(3, i) = 0
                        i6Cell_Cell(4, i) = i + (rtm_obj%Nz - 1)
                    else if (iy .eq. Nyc) then
                        i6Cell_Cell(3, i) = i - (rtm_obj%Nz - 1)
                        i6Cell_Cell(4, i) = 0
                    else
                        i6Cell_Cell(3, i) = i - (rtm_obj%Nz - 1)
                        i6Cell_Cell(4, i) = i + (rtm_obj%Nz - 1)
                    end if
                else
                    if (iy .eq. 1) then
                        i6Cell_Cell(3, i) = i + (Nyc - 1)*(rtm_obj%Nz - 1)
                    else
                        i6Cell_Cell(3, i) = i - (rtm_obj%Nz - 1)
                    end if
                    if (iy .eq. rtm_obj%Ny) then
                        i6Cell_Cell(4, i) = i - (Nyc - 1)*(rtm_obj%Nz - 1)
                    else
                        i6Cell_Cell(4, i) = i + (rtm_obj%Nz - 1)
                    end if
                end if
                if (iz .eq. 1) then
                    i6Cell_Cell(5, i) = 0
                else
                    i6Cell_Cell(5, i) = i - 1
                end if
                if (iz .eq. rtm_obj%Nz - 1) then
                    i6Cell_Cell(6, i) = 0
                else
                    i6Cell_Cell(6, i) = i + 1
                end if
                ! Base grid cells have no parents or children
                i2Tree_Cell(1, i) = 0
                i2Tree_Cell(2, i) = 0
            end do
        end do
    end do
    cell_obj%Npts = i

    return

end subroutine shd_init__makeCell

subroutine shd_init__makePrp(&
    doDeltaM, maxSsa,&
    rtm_obj, &
    temp_Rtm, planck_Rtm, ext_Rtm, ssa_Rtm, leg_Rtm, ipha_Rtm)

    ! Prepares the grid arrays for the iterative solution process.
    ! If doing Delta-M scaling then the extinction, albedo, and Legend
    ! terms are scaled first; only the 0 to shdo_obj_%maxNl leg_Rtm terms are scaled.
    ! Outputs planck_Rtm with (1-omega)*B(T) for thermal source, where B(T
    ! the function (meaning depends on src_obj_%unit).
    ! temp_Rtm array is unchanged.

    type(shd_type_rtm),  intent(inout) :: rtm_obj

    real(kind=R__), intent(in)    ::   temp_Rtm(*)
    real(kind=R__), intent(out)   :: planck_Rtm(*)
    real(kind=R__), intent(inout) ::    ext_Rtm(*)
    real(kind=R__), intent(inout) ::    ssa_Rtm(*)
    real(kind=R__), intent(inout) ::    leg_Rtm(prp_obj_%NlegSet, &
                                              0:prp_obj_%maxNleg, *)
    integer(kind=2), intent(in)   ::   ipha_Rtm(*)

    logical, intent(in) :: doDeltaM
    real(kind=R__), intent(inout) :: maxSsa
    integer I, IPH, L
    real F, BB

    maxSsa = 0.0

    if (doDeltaM .and. (prp_obj_%Npha .gt. 0)) then
        do IPH = 1, prp_obj_%Npha
            F = leg_Rtm(1, shdo_obj_%maxNl + 1, IPH)
            do L = 0, shdo_obj_%maxNl
                ! Scale the diagonal and off-diagonal phase matrix elements d
                leg_Rtm(1, L, IPH) = (leg_Rtm(1, L, IPH) - F)/(1 - F)
                if (prp_obj_%NlegSet .gt. 1) then
                    leg_Rtm(2, L, IPH) = (leg_Rtm(2, L, IPH) - F)/(1 - F)
                    leg_Rtm(3, L, IPH) = (leg_Rtm(3, L, IPH) - F)/(1 - F)
                    leg_Rtm(4, L, IPH) = (leg_Rtm(4, L, IPH) - F)/(1 - F)
                    leg_Rtm(5, L, IPH) =  leg_Rtm(5, L, IPH)     /(1 - F)
                    leg_Rtm(6, L, IPH) =  leg_Rtm(6, L, IPH)     /(1 - F)
                end if
            end do
        end do
    end if

    do I = 1, rtm_obj_%Npts
        if (doDeltaM) then
            if (prp_obj_%Npha .gt. 0) then
                F = leg_Rtm(1, shdo_obj_%maxNl + 1, ipha_Rtm(I))
            else
                ! prp_obj_%Npha =0 is standard property format, which is unpolarized
                F = leg_Rtm(1, shdo_obj_%maxNl + 1, I)
                do L = 0, shdo_obj_%maxNl
                    leg_Rtm(1, L, I) = (leg_Rtm(1, L, I) - F)/(1 - F)
                end do
            end if

            ext_Rtm(I) = (1.0 - ssa_Rtm(I)*F)*ext_Rtm(I)
            ssa_Rtm(I) = (1.0 - F)*ssa_Rtm(I)/(1.0 - ssa_Rtm(I)*F)

        end if

        maxSsa = max(maxSsa, ssa_Rtm(I))
        if (src_obj_%doTherm) then
            call shd_util__computePlanck(temp_Rtm(I), src_obj_%unit, src_obj_%waveNum, src_obj_%avgWaveLen, BB)
            planck_Rtm(I) = (1.0 - ssa_Rtm(I))*BB
        end if
    end do

    return

end subroutine shd_init__makePrp

subroutine shd_init__makeRad(&
    &NSTOKES, NX, NY, NZ,&
    &iradI_Shdo, ZGRID,&
    &ext, ssa, leg, temp, ipha,&
    &SRCTYPE, SOLARFLUX, muSun,&
    &SKYRAD, UNITS, WAVENO, WAVELEN, radI_Shdo)

    ! Initializes radiance field by solving plane-parallel two-stream.
    ! Solves the L=1 M=0 SH system by transforming the pentadiagonal
    ! system to tridiagonal and calling solver.
    ! Does the initialization for the NX*NY columns of the base grid.
    ! Only the I Stokes parameter is initialized, the rest are zeroed.

    integer NSTOKES, NX, NY, NZ, iradI_Shdo(*)

    integer(kind=2), intent(in) :: ipha(NZ, NY, NX)

    real SKYRAD, SOLARFLUX, muSun
    real WAVENO(2), WAVELEN
    real ZGRID(NZ)
    real ext(NZ, NY, NX), ssa(NZ, NY, NX), leg(prp_obj_%NlegSet, 0:prp_obj_%Nleg, *)
    real temp(NZ, NY, NX)
    real radI_Shdo(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer NLAYER, IX, IY, IZ, IR, I, J, K, L, ZINDEX
    logical DELTAM
    real C0, C1
    real EXT0, ALB0, SCAT0, G0, GNDEMIS
    real EXT1, SCAT1, G1
    real, allocatable :: tauZ(:), ssaZ(:), asyZ(:)
    real, allocatable :: tempZ(:), FLUXES(:, :)
    character(len=1) :: SRCT

    NLAYER = NZ - 1

    allocate(tauZ(NZ), ssaZ(NZ), asyZ(NZ))
    allocate(tempZ(NZ), FLUXES(3, NZ))

    !real, allocatable :: CH(:, :)
    !allocate(CH(NLAYER, NLAYER))
    !call shd_init__chapman(NLAYER, muSun, planetRadius_, ZGRID(NZ:1:-1), CH)

    DELTAM = .false.

    SRCT = SRCTYPE
    if (SRCT .eq. 'B') SRCT = 'T'

    ! SH coefficients
    C0 = sqrt(1.0d0/PI_)
    C1 = sqrt(3.0d0/(4.0d0*PI_))

    ! Loop over all the columns in the base grid
    IR = 0
    iradI_Shdo(1) = 0
    do IX = 1, NX
        do IY = 1, NY

            I = IY + NY*(IX - 1)

            if (doGridPtsPrp_) then

                do IZ = 1, NLAYER

                    L = NZ - IZ
                    EXT0  = ext(IZ,     IY, IX)
                    EXT1  = ext(IZ + 1, IY, IX)
                    SCAT0 = ssa(IZ,     IY, IX)*EXT0
                    SCAT1 = ssa(IZ + 1, IY, IX)*EXT1

                    tauZ(L) = (ZGRID(IZ + 1) - ZGRID(IZ))*(EXT0 + EXT1)/2.0

                    if ((EXT0 + EXT1) .gt. 0.0) then
                        ssaZ(L) = (SCAT0 + SCAT1)/(EXT0 + EXT1)
                    else
                        ssaZ(L) = 0.0
                    end if

                    if (prp_obj_%Npha .gt. 0) then
                        J  = ipha(IZ, IY, IX)
                        G0 = leg(1, 1, J)
                        J  = ipha(IZ + 1, IY, IX)
                        G1 = leg(1, 1, J)
                    else
                        J = NZ*(I - 1) + IZ
                        G0 = leg(1, 1, J)
                        G1 = leg(1, 1, J + 1)
                    end if

                    if (SCAT0 + SCAT1 .gt. 0.0) then
                        asyZ(L) = (SCAT0*G0 + SCAT1*G1)/(SCAT0 + SCAT1)
                    else
                        asyZ(L) = 0.0
                    end if

                    tempZ(L + 1) = temp(IZ, IY, IX)

                end do

            else

                ! zgrid correction (to address layer averaging issue in z direction)
                !╭────────────────────────────────────────────────╮!
                do IZ = 1, NLAYER

                    L = NZ - IZ

                    ZINDEX = IZ + 1

                    EXT0 = ext(ZINDEX, IY, IX)
                    ALB0 = ssa(ZINDEX, IY, IX)
                    SCAT0 = ALB0*EXT0

                    tauZ(L) = (ZGRID(IZ + 1) - ZGRID(IZ))*EXT0

                    if (EXT0 .gt. 0.0) then
                        ssaZ(L) = ALB0
                    else
                        ssaZ(L) = 0.0
                    end if

                    if (prp_obj_%Npha .gt. 0) then
                        J = ipha(ZINDEX, IY, IX)
                    else
                        J = NZ*(I - 1) + IZ
                    end if
                    G0 = leg(1, 1, J)

                    if (SCAT0 .gt. 0.0) then
                        asyZ(L) = G0
                    else
                        asyZ(L) = 0.0
                    end if
                    tempZ(L + 1) = temp(IZ, IY, IX)
                end do
                !╰────────────────────────────────────────────────╯!

            end if

            tempZ(1) = temp(NZ, IY, IX)

            GNDEMIS = 1.0 - sfc_obj_%avgAlb

            if (src_obj_%doVariableSolar) then
                lon0 = prp_obj_%xStart + (IX-1)*rtm_obj_%dx + rtm_obj_%xOffset
                lat0 = prp_obj_%yStart + (IY-1)*rtm_obj_%dy + rtm_obj_%yOffset
                call shd_util__computeSolarAngle(lon0, lat0, mu0, phi0)
            else
                mu0 = muSun
            end if

            ! Call the Eddington flux routine
            call shd_init__makeFlux(NLAYER, tauZ, ssaZ, asyZ,&
                &tempZ, DELTAM, SRCTYPE, SOLARFLUX, mu0,&
                &sfc_obj_%avgTemp, GNDEMIS, SKYRAD, UNITS, WAVENO, WAVELEN,&
                &FLUXES)

            ! Convert fluxes to first two moments of spherical harmonics
            do IZ = 1, NZ
                L = NZ + 1 - IZ
                radI_Shdo(:, IR + 1) = 0.0
                radI_Shdo(:, IR + 2) = 0.0
                radI_Shdo(:, IR + 3) = 0.0
                radI_Shdo(:, IR + 4) = 0.0
                radI_Shdo(1, IR + 1) = C0*(FLUXES(1, L) + FLUXES(2, L))
                radI_Shdo(1, IR + 3) = C1*(FLUXES(1, L) - FLUXES(2, L))
                IR = IR + 4
                iradI_Shdo(IZ + NZ*(I - 1) + 1) = IR
            end do

        end do
    end do

    deallocate(tauZ, ssaZ, asyZ, tempZ, FLUXES)

    return

end subroutine shd_init__makeRad

subroutine shd_init__makeFlux(&
    &NLAYER, tauZ, ssaZ, asyZ,&
    &tempZ, DELTAM, SRCTYPE, SOLARFLUX, muSun,&
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
    ! tauZ           real array   Optical thickness of layers
    ! ssaZ           real array   Single scattering albedos
    ! asyZ           real array   Asymmetry parameters
    ! tempZ          real array   Temperatures (K) at layer interfaces
    !                              (e.g. tempZ(1) is at top of top layer
    !                               tempZ(2) is at bottom of top layer).
    ! DELTAM         logical      True for delta-Eddington scaling
    ! SRCTYPE        character    'S' for solar source, 'T' for thermal source
    ! SOLARFLUX      real         Incident solar flux on horizonal plane
    ! muSun        real         Cosine of the solar zenith angle
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
    real tempZ(*)
    real tauZ(*), ssaZ(*), asyZ(*)
    real GNDTEMP, GNDEMIS, SKYRAD, SOLARFLUX, muSun
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
    allocate(LOWER(N), UPPER(N), DIAG(N), RHS(N))

    ! Compute the reflection, transmission, and source
    !   coefficients for each layer for the diffuse Eddington
    !   two stream problem.
    if (SRCTYPE .eq. 'T') then
        call shd_util__computePlanck(tempZ(1), UNITS, WAVENO, WAVELEN, BBRAD)
        PLANCK1 = PI_*BBRAD
    end if
    MU0 = abs(muSun)
    TAU = 0.0
    I = 2
    do L = 1, NLAYER
        DELTAU = tauZ(L)
        if (DELTAU .lt. 0.0) stop 'Error [_makeFlux]: TAU<0.'
        ! Special case for zero optical depth
        if (DELTAU .eq. 0.0) then
            TRANS = 1.0
            REFLECT = 0.0
            SOURCEP = 0.0
            SOURCEM = 0.0
        else
            OMEGA = ssaZ(L)
            G = asyZ(L)
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
                call shd_util__computePlanck(tempZ(L + 1), UNITS, WAVENO, WAVELEN, BBRAD)
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

    deallocate(LOWER, UPPER, DIAG, RHS)

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
    if (src_obj_%doSolar) then
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

subroutine shd_init__makeBotRadUL(NBOTPTS, BCPTR,&
    &fluxDir, FLUXES, SRCTYPE,&
    &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD)

    !   Returns the BCRAD array with the upwelling unpolarized radiance for
    ! a Lambertian surface for the bottom points (pointed to by BCPTR).
    ! If there is reflection the downwelling hemispheric flux at the
    ! bottom points is gotten (from FLUXES). There may also be thermal
    ! emission (depending on temperature and emissivity) or reflected
    ! collimated solar radiation (fluxDir).

    integer NBOTPTS, BCPTR(*), NSTOKES
    real fluxDir(*), FLUXES(2, *)
    real WAVENO(2), WAVELEN
    real BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IBC, K
    real ALB, GNDRAD

    ALB = sfc_obj_%avgAlb * OPI_
    if (src_obj_%doTherm) then
        call shd_util__computePlanck(sfc_obj_%avgTemp, UNITS, WAVENO, WAVELEN, GNDRAD)
        GNDRAD = GNDRAD*(1.0 - sfc_obj_%avgAlb)
    end if

    ! Loop over all bottom points finding the upwelling radiance
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        if (SRCTYPE .eq. 'S') then
            BCRAD(1, IBC) = ALB*(fluxDir(I) + FLUXES(1, I))
        else if (SRCTYPE .eq. 'T') then
            BCRAD(1, IBC) = GNDRAD + ALB*FLUXES(1, I)
        else if (SRCTYPE .eq. 'B') then
            BCRAD(1, IBC) = ALB*(fluxDir(I) + FLUXES(1, I)) + GNDRAD
        end if
        BCRAD(2:, IBC) = 0.0
    end do

    return

end subroutine shd_init__makeBotRadUL

subroutine shd_init__makeBotRadVL(NBOTPTS, BCPTR,&
    &fluxDir, FLUXES, SRCTYPE,&
    &NSTOKES, BCRAD)

    !   Returns the BCRAD array with the upwelling unpolarized radiance for a
    ! variable Lambertian surface for the bottom points (pointed to by BCPTR).
    ! brdf_Rtm_ contains the interpolated surface properties for the
    ! bottom boundary points; the first is Planck function and the second is
    ! albedo.  The reflection is computed from the downwelling hemispheric
    ! flux (in FLUXES).  There may also be thermal emission (depending on
    ! temperature and emissivity) or reflected collimated solar radiation
    ! (from fluxDir).

    integer NBOTPTS, BCPTR(*), NSTOKES
    real fluxDir(*), FLUXES(2, *)
    real BCRAD(NSTOKES, *)
    character(len=1) :: SRCTYPE
    integer I, IBC, K
    real ALB, GNDRAD

    ! Loop over all bottom points
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        ALB = brdf_Rtm_(2, IBC)
        if (SRCTYPE .eq. 'S') then
            BCRAD(1, IBC) = OPI_*ALB*(fluxDir(I) + FLUXES(1, I))
        else
            GNDRAD = brdf_Rtm_(1, IBC)*(1 - ALB)
            if (SRCTYPE .eq. 'T') then
                BCRAD(1, IBC) = GNDRAD + OPI_*ALB*FLUXES(1, I)
            else if (SRCTYPE .eq. 'B') then
                BCRAD(1, IBC) = OPI_*ALB*(fluxDir(I) + FLUXES(1, I)) + GNDRAD
            end if
        end if
        BCRAD(2:, IBC) = 0.0
    end do

    return

end subroutine shd_init__makeBotRadVL

subroutine shd_init__makeBotRadU(NBOTPTS, BCPTR,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, IMU, IPHI,&
    &SRCTYPE, WAVELEN, muSun, phiSun, fluxDir,&
    &brdf_Shdo,&
    &NSTOKES, BCRAD)

    !   Computes the upwelling Stokes vector of the bottom boundary points
    ! for one outgoing discrete ordinate using the precomputed bidirectional
    ! reflectance distribution matrix for the discrete ordinates.
    ! The upwelling Stokes vector includes the reflection of the incident
    ! radiance, the thermal emission (emissivity*Planck function), and the
    ! reflected direct solar flux (if applicable). The upwelling radiance
    ! vector is the integral over all incident directions of the BRDF times
    ! the downwelling radiance, so a discrete sum is done using the brdf_Shdo
    ! matrix, which includes the integral weights (WTDO). The emissivity
    ! is computed implicitly from the integral of the BRDF. The outgoing
    ! direction is specified with (IMU,IPHI), and the BRDF is computed for
    ! all incident directions (loop over JMU,JPHI). The incident downwelling
    ! radiance vectors are input in BCRAD(*,*,2...) and the outgoing
    ! upwelling radiance vectors are output in BCRAD(*,*,1).

    integer NBOTPTS, BCPTR(NBOTPTS)
    integer NMU, NPHI0MAX, NPHI0(NMU), IMU, IPHI, NSTOKES
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, muSun, phiSun, fluxDir(*)
    real brdf_Shdo(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    real BCRAD(NSTOKES, NBOTPTS, *)
    character(len=1) :: SRCTYPE
    integer IMU1, JMU, JPHI, IBC, I, JANG, K, K1, K2
    real REFLECT(4, 4), W

    IMU1 = IMU - NMU/2
    do IBC = 1, NBOTPTS

        ! Initialize the upwelling boundary radiances to zero or to
        !   the reflection of direct unpolarized solar flux.
        BCRAD(:, IBC, 1) = 0.0
        if (src_obj_%doSolar) then
            I = BCPTR(IBC)
            call shd_sfc__brdfSurface(sfc_obj_%type(2:2), brdf_Rtm_(2:sfc_obj_%Npar, IBC), WAVELEN,&
                &MU(IMU), PHI(IMU, IPHI), muSun, phiSun,&
                &NSTOKES, REFLECT)
            BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                &+ OPI_*REFLECT(1:NSTOKES, 1)*fluxDir(I)
        end if

        ! Integrate over the incident discrete ordinate directions (JMU, JPHI)
        JANG = 1
        do JMU = 1, NMU/2
            do JPHI = 1, NPHI0(JMU)
                ! Add in the polarized reflection
                do K1 = 1, NSTOKES
                    BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1) + BCRAD(K1, IBC, JANG + 1)&
                        &*brdf_Shdo(:, IMU1, IPHI, K1, JMU, JPHI)
                end do
                ! Add in the polarized thermal emission
                W = OPI_*abs(MU(JMU))*WTDO(JMU, JPHI)
                BCRAD(1, IBC, 1) = BCRAD(1, IBC, 1) + brdf_Rtm_(1, IBC)&
                    &*(W - brdf_Shdo(1, IMU1, IPHI, 1, JMU, JPHI))
                BCRAD(2:, IBC, 1) = BCRAD(2:, IBC, 1) - brdf_Rtm_(1, IBC)&
                    &*brdf_Shdo(2:, IMU1, IPHI, 1, JMU, JPHI)
                JANG = JANG + 1
            end do
        end do
    end do

    return

end subroutine shd_init__makeBotRadU

subroutine shd_init__makeBotRadV(&
    &NBOTPTS, IBEG, IEND, BCPTR,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,&
    &SRCTYPE, WAVELEN, muSun, phiSun, fluxDir,&
    &NSTOKES, BCRAD)

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
    ! type (sfc_obj_%type) with parameters (brdf_Rtm_) for the incident and
    ! outgoing directions.  The emissivity is computed implicitly from the
    ! integral of the BRDF.  The outgoing direction is specified with
    ! (MU2,PHI2), and the BRDF is computed for all incident directions
    ! (loop over JMU,JPHI).  The incident downwelling radiance vectors are
    ! input in BCRAD(*,*,2...) and the outgoing upwelling radiance vectors
    ! are output in BCRAD(*,*,1).

    integer NBOTPTS, IBEG, IEND, BCPTR(NBOTPTS)
    integer NMU, NPHI0MAX, NPHI0(NMU), NSTOKES
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real MU2, PHI2, WAVELEN, muSun, phiSun, fluxDir(*)
    real BCRAD(NSTOKES, NBOTPTS, *)
    character(len=1) :: SRCTYPE
    integer JMU, JPHI, IBC, I, JANG, K, K1, K2
    real REFLECT(4, 4), W, SUM0, SUM1, RADSPEC(4)

    do IBC = IBEG, IEND

        ! Initialize the upwelling boundary radiances to zero or to
        !   the reflection of direct unpolarized solar flux.
        BCRAD(:, IBC, 1) = 0.0
        if (src_obj_%doSolar) then
            I = BCPTR(IBC)

            if (src_obj_%doVariableSolar) then
                lon0 = xyz_Rtm_(1, I) + rtm_obj_%xOffset
                lat0 = xyz_Rtm_(2, I) + rtm_obj_%yOffset
                call shd_util__computeSolarAngle(lon0, lat0, mu0, phi0)
            else
                mu0  = muSun
                phi0 = phiSun
            end if

            call shd_sfc__brdfSurface(sfc_obj_%type(2:2), brdf_Rtm_(2:sfc_obj_%Npar, IBC), WAVELEN,&
                &MU2, PHI2, mu0, phi0, NSTOKES, REFLECT)

            BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1) + &
                OPI_*REFLECT(1:NSTOKES, 1)*fluxDir(I)
        end if

        ! Integrate over the incident discrete ordinate directions (JMU,JPHI)
        JANG = 1
        do JMU = 1, NMU/2
            do JPHI = 1, NPHI0(JMU)
                call shd_sfc__brdfSurface(sfc_obj_%type(2:2), brdf_Rtm_(2:sfc_obj_%Npar, IBC),&
                    &WAVELEN, MU2, PHI2, MU(JMU), PHI(JMU, JPHI),&
                    &NSTOKES, REFLECT)
                W = OPI_*abs(MU(JMU))*WTDO(JMU, JPHI)
                ! Add in the polarized reflection
                do K1 = 1, NSTOKES
                    BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1) + &
                        W*REFLECT(1:NSTOKES, K1)*BCRAD(K1, IBC, JANG + 1)
                end do
                ! Add in the polarized thermal emission
                BCRAD(1, IBC, 1) = BCRAD(1, IBC, 1) + &
                    W*(1 - REFLECT(1, 1))*brdf_Rtm_(1, IBC)
                BCRAD(2:, IBC, 1) = BCRAD(2:, IBC, 1) - &
                    W*REFLECT(2:NSTOKES, 1)*brdf_Rtm_(1, IBC)
                JANG = JANG + 1
            end do
        end do

    end do

    return

end subroutine shd_init__makeBotRadV

subroutine shd_init__checkSfcU(sfc_obj, brdf_Sfc)

    ! Description:
    !   Checks to see if the surface reflection parameters are uniform,
    !   and if so sets the sfc_obj%isUniform flag.  The temperature need
    !   not be uniform.
    !
    ! Input:
    !   sfc_obj%Nx, sfc_obj%Ny: integer, size of surface arrays (internally sfc_obj%Nx+1 by sfc_obj%Ny+1)
    !   sfc_obj%Npar: integer, number of surface parameters for modeling different BRDF
    !   brdf_Sfc: real array, surface parameters for modeling different BRDF
    !
    ! Output:
    !   sfc_obj%isUniform: logical, whether surface reflection parameters are uniform

    type(shd_type_sfc), intent(inout) :: sfc_obj
    real(R__),          intent(in)    :: brdf_Sfc(sfc_obj%Npar, *)

    integer(kind=I__) :: ix, iy, J, idata

    sfc_obj%isUniform = .true.

    do iy = 1, sfc_obj%Ny
        do ix = 1, sfc_obj%Nx
            idata = iy + (sfc_obj%Ny + 1)*(ix - 1)
            do J = 2, sfc_obj%Npar
                if (brdf_Sfc(J, idata) .ne. brdf_Sfc(J, 1)) then

                    sfc_obj%isUniform = .false.

                    return

                end if
            end do
        end do
    end do

    return

end subroutine shd_init__checkSfcU

subroutine shd_init__chapman(&
    &NLAYER, muSun, planetRadius, ZLEVELS,&
    &CH)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Apply Chapman function for psuedo spherical correction
    !
    ! Input:
    !   NLAYER : integer, number of layers
    !   muSun: real, cosine of incident angle
    !   planetRadius: real, planet radius in km
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
    real(kind=R__)   , intent(in) :: planetRadius, ZLEVELS(0:NLAYER), muSun

    real(kind=R__), intent(out) :: CH(NLAYER, NLAYER)

    integer(kind=I__) :: lc, n
    real(kind=R__)    :: cos_theta0, sin_theta0, cos_theta, sin_theta
    real(kind=R__)    :: s1, s2
    real(kind=R__)    :: h(0:NLAYER)

    if (h(0) .lt. h(NLAYER)) stop 'Error [_chapman]: Input ZLEVELS should go from top to bottom.'

    cos_theta0 = muSun
    sin_theta0 = sqrt(1.0 - muSun*muSun)
    h = ZLEVELS + planetRadius

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
                    CH(lc,n) = muSun
                end if
            endif
            s1 = s2

        end do

    end do

    return

end subroutine shd_init__chapman


end module shd_init
