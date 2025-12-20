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
!
!                Vikas Nataraja
!                vikas.hanasogenataraja@lasp.colorado.edu
!╰────────────────────────────────────────────────────────────────────────────╯!

module shd_util

    use globals_const, only: RD_, PI_, DPI_, deg2rad_
    use globals_cfg,   only: charMaxLen_, planetRadiusA_, planetRadiusB_, &
                             A2_, B2_, E2_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, &
                             is1DX_, is1DYx_, verbose_

    use globals_obj,   only: rtm_obj_, xyz_Rtm_, &
                            cell_obj_, i8Grid_Cell_, &
                                       i6Cell_Cell_, &
                                       i2Tree_Cell_, &
                                         flag_Cell_

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_util__sumData
    public :: shd_util__interpData

    public :: shd_util__findCell
    public :: shd_util__findNextCell

    public :: shd_util__interpTriLin
    public :: shd_util__solveTriDiag
    public :: shd_util__computeGausQuad_new
    public :: shd_util__computeGausQuad
    public :: shd_util__computeGausQuadDble

    public :: shd_util__computePlanck
    public :: shd_util__integratePlanck

    public :: shd_util__ssort
    public :: shd_util__ssort_i
    public :: shd_util__ssort_r

    public :: indexAt
    public :: distKm

    public :: shd_util__computeRayVect
    public :: shd_util__computeDegKmRatio
    public :: shd_util__computeVecECEF
    public :: shd_util__computeDirENU

    public :: shd_util__tickTock

    private :: integralLow
    private :: integralHigh

    interface shd_util__ssort
        module procedure shd_util__ssort_i
        module procedure shd_util__ssort_r
    end interface

contains

subroutine shd_util__sumData(&
    iFlag, wgt, Ndata, data_in, &
    data_sum)

    ! Sums an output array over the k-distribution.
    ! If iFlag=0 then the output array is set to the input times the weight.

    integer, intent(in) :: iFlag, Ndata
    real,    intent(in) :: wgt, data_in(Ndata)

    real, intent(out) :: data_sum(Ndata)

    integer :: iData

    if (iFlag .eq. 0) then
        do iData = 1, Ndata
            data_sum(iData) = wgt*data_in(iData)
        end do
    else
        do iData = 1, Ndata
            data_sum(iData) = data_sum(iData) + wgt*data_in(iData)
        end do
    end if

    return

end subroutine shd_util__sumData

subroutine shd_util__interpData(&
    isFluxOut, isHROut, isSHOut, &
    NptsOld, &
    fluxTot, fluxDir, fluxDiv, SHterms)

    !   Interpolates the output data to the new grid points (those >NptsOld)
    ! from the old grid points.  The output data is in fluxTot, fluxDir,
    ! fluxDiv, and SHterms.  NptsOld is set to zero to do no interpolation.
    ! Each new point is located in a cell, and the two parent grid points
    ! are found to interpolate the data from.

    logical, intent(in) :: isFluxOut, isHROut, isSHOut
    integer(kind=I__), intent(in) :: NptsOld

    real(kind=R__), intent(out) :: fluxTot(2, rtm_obj_%Npts), &
                                   fluxDir(rtm_obj_%Npts), &
                                   fluxDiv(rtm_obj_%Npts), &
                                   SHterms(5, rtm_obj_%Npts)

    integer IP, IC, iCell, iParent, iDir, I, J, DIR1, DIR2, IP1, IP2

    integer cornerGrid(2, 4, 3)
    data cornerGrid/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        &1, 5, 2, 6, 3, 7, 4, 8/

    if (NptsOld .le. 0) return

    ! Go through all the non-base grid cells
    do iCell = cell_obj_%Npts0+1, cell_obj_%Npts
        do IC = 1, 8
            IP = i8Grid_Cell_(IC, iCell)
            if (IP .gt. NptsOld) then
                iParent = i2Tree_Cell_(1, iCell)
                ! Go through all 12 edges of parent cell to find which one has point
                do iDir = 1, 3
                    DIR1 = mod(iDir - 1 + 1, 3) + 1
                    DIR2 = mod(iDir - 1 + 2, 3) + 1
                    do I = 1, 4
                        IP1 = i8Grid_Cell_(cornerGrid(1, I, iDir), iParent)
                        IP2 = i8Grid_Cell_(cornerGrid(2, I, iDir), iParent)
                        if  (xyz_Rtm_(DIR1, IP1) .eq. xyz_Rtm_(DIR1, IP) .and.&
                            &xyz_Rtm_(DIR2, IP1) .eq. xyz_Rtm_(DIR2, IP)) goto 120
                    end do
                end do

                stop 'Error [_interpOutput]: Point not on an edge.'

120             continue

                ! Interpolate the data from the two parent points at the
                !   ends of the edge; assume the cell was split in half.
                if (isFluxOut) then
                    fluxTot(1, IP) = 0.5*(fluxTot(1, IP1) + fluxTot(1, IP2))
                    fluxTot(2, IP) = 0.5*(fluxTot(2, IP1) + fluxTot(2, IP2))
                    fluxDir(IP) = 0.5*(fluxDir(IP1) + fluxDir(IP2))
                end if
                if (isHROut) then
                    fluxDiv(IP) = 0.5*(fluxDiv(IP1) + fluxDiv(IP2))
                end if
                if (isSHOut) then
                    do J = 1, 5
                        SHterms(J, IP) = 0.5*(SHterms(J, IP1) + SHterms(J, IP2))
                    end do
                end if
            end if
        end do
    end do

    return

end subroutine shd_util__interpData

subroutine shd_util__findCell(&
    x, y, z, &
    iCell)

    !   Locates the grid cell in the tree structure containing the
    ! specified point (x,y,z), and returns the cell pointer iCell.
    ! First the base grid cell is found (using
    ! rtm_obj_%Nx,rtm_obj_%Ny,rtm_obj_%Nz,
    ! rtm_obj_%x,rtm_obj_%y,rtm_obj_%z),
    ! and the tree is traced down to find the smallest cell containing
    ! the point.

    real(kind=RD_), intent(inout) :: x, y, z

    integer(kind=I__), intent(out) :: iCell

    integer :: IL, IU, IM, IX, IY, IZ, NXC, NYC, IC, IPTR, DIR
    real(kind=RD_) :: xDomain, yDomain

    ! Move point to within horizontal domain if periodic
    if (.not. (isOpenX_ .or. isMpiX_)) then
        xDomain = rtm_obj_%x(rtm_obj_%Nx + 1) - rtm_obj_%x(1)
        if (x .lt. rtm_obj_%x(1)) then
            x = x - xDomain*(int((x - rtm_obj_%x(1))/xDomain) - 1)
        else if (x .gt. rtm_obj_%x(rtm_obj_%Nx + 1)) then
            x = x - xDomain*int((x - rtm_obj_%x(1))/xDomain)
        end if
    end if

    if (.not. (isOpenYx_ .or. isMpiY_)) then
        yDomain = rtm_obj_%y(rtm_obj_%Ny + 1) - rtm_obj_%y(1)
        if (y .lt. rtm_obj_%y(1)) then
            y = y - yDomain*(int((y - rtm_obj_%y(1))/yDomain) - 1)
        else if (y .gt. rtm_obj_%y(rtm_obj_%Ny + 1)) then
            y = y - yDomain*int((y - rtm_obj_%y(1))/yDomain)
        end if
    end if

    ! Find the base grid cell (base grid does not have to be even sp
    !   so do binary search in x, y, and z)
    !   Special case for independent pixel: get closest grid column.
    IL = 0
    if (is1DX_) then
        IU = rtm_obj_%Nx
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (x .ge. 0.5*(rtm_obj_%x(IM) + rtm_obj_%x(IM + 1))) then
                IL = IM
            else
                IU = IM
            end if
        end do
        IL = IL + 1
    else
        IU = rtm_obj_%Nx + 1
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (x .ge. rtm_obj_%x(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
    end if
    IX = max(IL, 1)

    IL = 0
    if (is1DYx_) then
        IU = rtm_obj_%Ny
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (y .ge. 0.5*(rtm_obj_%y(IM) + rtm_obj_%y(IM + 1))) then
                IL = IM
            else
                IU = IM
            end if
        end do
        IL = IL + 1
    else
        IU = rtm_obj_%Ny + 1
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (y .ge. rtm_obj_%y(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
    end if
    IY = max(IL, 1)

    IL = 0
    IU = rtm_obj_%Nz
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (z .ge. rtm_obj_%z(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)

    ! Special case for open boundary conditions: if outside the 1 to
    !   domain then use the IP grid cells at either end
    NXC = rtm_obj_%Nx
    if (isOpenX_) then
        NXC = rtm_obj_%Nx + 1
        if (x .lt. rtm_obj_%x(1)) then
            ! x = rtm_obj_%x(1)
            IX = 1
        else if (x .gt. rtm_obj_%x(rtm_obj_%Nx)) then
            ! x = rtm_obj_%x(rtm_obj_%Nx)
            IX = rtm_obj_%Nx + 1
        else
            IX = IX + 1
        end if
    end if

    if (isMpiX_ .and. .not. is1DX_) then
        NXC = rtm_obj_%Nx - 1
        IX = min(IX, NXC)
    end if

    NYC = rtm_obj_%Ny
    if (isOpenYx_) then
        NYC = rtm_obj_%Ny + 1
        if (y .lt. rtm_obj_%y(1)) then
            ! y = rtm_obj_%y(1)
            IY = 1
        else if (y .gt. rtm_obj_%y(rtm_obj_%Ny)) then
            ! y = rtm_obj_%y(rtm_obj_%Ny)
            IY = rtm_obj_%Ny + 1
        else
            IY = IY + 1
        end if
    end if

    if (isMpiY_ .and. .not. is1DYx_) then
        NYC = rtm_obj_%Ny - 1
        IY = min(IY, NYC)
    end if

    ! Get the base grid cell pointer
    iCell = IZ + (rtm_obj_%Nz - 1)*(IY - 1) + (rtm_obj_%Nz - 1)*NYC*(IX - 1)

    !   Trace down the tree: point to the positive child cell, get the
    ! grid point that is on the negative (x,y,z) side, use the
    ! flags to find which dimension (x,y,z) the parent cell was
    ! split in, then compare the grid point with the test point
    ! (x,y,z) to find which child cell (pos or neg) the test
    ! point is in.
    do while (i2Tree_Cell_(2, iCell) .gt. 0)
        DIR = ibits(int(flag_Cell_(iCell)), 2, 2)
        IC = i2Tree_Cell_(2, iCell) + 1
        IPTR = i8Grid_Cell_(1, IC)
        if (DIR .eq. 1) then
            if (x .lt. xyz_Rtm_(1, IPTR)) IC = IC - 1
        else if (DIR .eq. 2) then
            if (y .lt. xyz_Rtm_(2, IPTR)) IC = IC - 1
        else if (DIR .eq. 3) then
            if (z .lt. xyz_Rtm_(3, IPTR)) IC = IC - 1
        end if
        iCell = IC
    end do

    return

end subroutine shd_util__findCell

subroutine shd_util__findNextCell(&
    x, y, z, &
    iFace, jFace, iCell, &
    iCellNext)

    !   Gets the next cell (iCellNext) to go to that is adjacent to this face
    ! iFace (1-6) has which face (-x,+x,-y,+y,-z,+z) of the current cell
    ! (iCell) we are at, and jFace (1-3) is the direction (x,y,z).
    ! The current location (x,y,z) on the face is used to find the cell
    ! if there are several adjacent to the current cell.

    real(kind=RD_),    intent(in) :: x, y, z
    integer(kind=I__), intent(in) :: iFace, jFace, iCell

    integer(kind=I__), intent(out) :: iCellNext

    integer :: iCellTmp0, iCellTmp1, iDir

    ! Get the pointer to the neighbor cell
    iCellNext = i6Cell_Cell_(iFace, iCell)

    if (iCellNext .lt. 0) then
        ! If neighbor pointer is negative then we are going into a cell
        !   that has subcells, and the absolute value of the pointer is the
        !   cell to start the tree search.
        iCellTmp0 = abs(iCellNext)

        ! Go down the tree until there is no children
        do while (i2Tree_Cell_(2, iCellTmp0) .gt. 0)
            ! Get the direction cell is split
            iDir = ibits(int(flag_Cell_(iCellTmp0)), 2, 2)
            if (iDir .eq. 0) &
                stop 'Error [_findNextCell]: No split direction.'

            ! If split is in direction of face then choose the near child cell
            !   which is the negative side cell if this is a positive face.
            iCellTmp1 = i2Tree_Cell_(2, iCellTmp0)

            if (iDir .eq. jFace) then
                iCellTmp0 = iCellTmp1 + 1 - mod(iFace - 1, 2)
            else
                ! If split is not in direction of face then find which of two cells
                !   has a face that contains the current location.
                ! Compare the appropriate coordinate (x,y,z) of the grid point
                !   on the split line with the current location.
                iCellTmp0 = iCellTmp1
                if (iDir .eq. 1) then
                    if (x .gt. xyz_Rtm_(1, i8Grid_Cell_(8, iCellTmp1))) &
                        iCellTmp0 = iCellTmp0 + 1
                else if (iDir .eq. 2) then
                    if (y .gt. xyz_Rtm_(2, i8Grid_Cell_(8, iCellTmp1))) &
                        iCellTmp0 = iCellTmp0 + 1
                else
                    if (z .gt. xyz_Rtm_(3, i8Grid_Cell_(8, iCellTmp1))) &
                        iCellTmp0 = iCellTmp0 + 1
                end if
            end if
        end do
        iCellNext = iCellTmp0
    end if

    return

end subroutine shd_util__findNextCell

subroutine shd_util__interpTriLin(&
    iCell, x, y, z, &
    Npar, data_in, &
    iParStart, NparOut, &
    data_out)

    !   Trilinearly interpolates the field (data_in) given a grid cell
    ! pointer (iCell), to compute the field values (data_out) at the desired
    ! point (x,y,z). Npar is the leading index of the data_in (1 for a
    ! scalar field), iParStart is the first element of the vector to interpolate,
    ! and NparOut elements are output.

    integer, intent(in) :: iCell, Npar, iParStart, NparOut
    real(kind=RD_), intent(in) :: x, y, z
    real(kind=R__), intent(in) :: data_in(Npar, *)

    real(kind=R__), intent(out) :: data_out(NparOut)

    integer :: iPts1, iPts2, iPar
    real(RD_) :: U, V, W, UM, VM, WM, dx, dy, dz

    iPts1 = i8Grid_Cell_(1, iCell)
    iPts2 = i8Grid_Cell_(8, iCell)

    dx = xyz_Rtm_(1, iPts2) - xyz_Rtm_(1, iPts1)
    if (dx .le. 0.0) dx = 1.0

    dy = xyz_Rtm_(2, iPts2) - xyz_Rtm_(2, iPts1)
    if (dy .le. 0.0) dy = 1.0

    dz = xyz_Rtm_(3, iPts2) - xyz_Rtm_(3, iPts1)

    U = (x - xyz_Rtm_(1, iPts1))/dx
    V = (y - xyz_Rtm_(2, iPts1))/dy
    W = (z - xyz_Rtm_(3, iPts1))/dz
    UM = 1 - U
    VM = 1 - V
    WM = 1 - W

    do iPar = iParStart, iParStart + NparOut - 1

        data_out(iPar - iParStart + 1) = &
            WM*(VM*(UM*data_in(iPar, i8Grid_Cell_(1, iCell)) + &
                     U*data_in(iPar, i8Grid_Cell_(2, iCell))) + &
                 V*(UM*data_in(iPar, i8Grid_Cell_(3, iCell)) + &
                     U*data_in(iPar, i8Grid_Cell_(4, iCell)))) + &
             W*(VM*(UM*data_in(iPar, i8Grid_Cell_(5, iCell)) + &
                     U*data_in(iPar, i8Grid_Cell_(6, iCell))) + &
                 V*(UM*data_in(iPar, i8Grid_Cell_(7, iCell)) + &
                     U*data_in(iPar, i8Grid_Cell_(8, iCell))))

    end do

    return

end subroutine shd_util__interpTriLin

subroutine shd_util__solveTriDiag(&
    &N, LOWER, DIAG, UPPER, RHS)

    ! Computes the solution to a tridiagonal system.
    ! N is order of the matrix.  LOWER(2..N) is the subdiagonal,
    ! DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the
    ! superdiagonal.  On input RHS is the right hand side, while
    ! on output it is the solution vector. Everything is destroyed.
    ! Hacked from Linpack DGTSL.

    integer, intent(in) :: N
    real(RD_) :: LOWER(*), DIAG(*), UPPER(*), RHS(*)
    integer K, KB
    real(RD_) :: T

    if (N .eq. 1) then
        if (DIAG(1) .eq. 0.0) stop 'Error [_solveTriDiag]: Singular matrix.'
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
        if (LOWER(K) .eq. 0.0) stop 'Error [_solveTriDiag]: Singular matrix.'
        T = -LOWER(K + 1)/LOWER(K)
        LOWER(K + 1) = DIAG(K + 1) + T*DIAG(K)
        DIAG(K + 1) = UPPER(K + 1) + T*UPPER(K)
        UPPER(K + 1) = 0.0
        RHS(K + 1) = RHS(K + 1) + T*RHS(K)
    end do

    if (LOWER(N) .eq. 0.0) stop 'Error [_solveTriDiag]: Singular matrix.'

    ! Back substitute
    RHS(N) = RHS(N)/LOWER(N)
    RHS(N - 1) = (RHS(N - 1) - DIAG(N - 1)*RHS(N))/LOWER(N - 1)
    do KB = 1, N - 2
        K = N - 2 - KB + 1
        RHS(K) = (RHS(K) - DIAG(K)*RHS(K + 1) - UPPER(K)*RHS(K + 2))/LOWER(K)
    end do

    return

end subroutine shd_util__solveTriDiag

subroutine shd_util__computeGausQuad_new(m, gmu, gwt)

    ! this is adapted from libRadtran/libsrc_c/cdisort.c <c_gaussian_quadrature> for testing

    integer, intent(in) :: m
    real, intent(out) :: gmu(*), gwt(*)
    integer :: iter, k, lim, nn, np1, MAXIT
    double precision :: cona, t, en, nnp1, p, p2pri, pm1, pm2, ppr, prod, tmp, x, xi
    double precision :: tol
    logical :: initialized

    MAXIT = 1000
    initialized = .false.

    if (.not. initialized) then
        tol = 10.0d0 * epsilon(1.0d0)
        initialized = .true.
    end if

    lim = m / 2
    if (2*lim .ne. m) stop 'Error [_computeGausQuad_new]: N must be even.'

    en = real(m)
    np1 = m + 1
    nnp1 = m * np1
    cona = real(m - 1) / (8.0d0 * m**3)

    do k = 1, lim
        t = real(4 * k - 1) * PI_ / (4.0d0 * m + 2.0d0)
        x = cos(t + cona / tan(t))

        do iter = 1, MAXIT

            if (iter > MAXIT) then
                stop "Error [_computeGausQuad_new]: Max iteration count."
            end if

            pm2 = 1.0d0
            pm1 = x
            do nn = 2, m
                p = (real(2 * nn - 1) * x * pm1 - real(nn - 1) * pm2) / nn
                pm2 = pm1
                pm1 = p
            end do

            tmp = 1.0d0 / (1.0d0 - x**2)
            ppr = en * (pm2 - x * p) * tmp
            p2pri = (2.0d0 * x * ppr - nnp1 * p) * tmp
            xi = x - p / ppr * (1.0d0 + p / ppr * p2pri / (2.0d0 * ppr))

            if (abs(xi - x) .le. tol) then
                exit
            else
                x = xi
            end if
        end do

        gmu(k) = -x
        gwt(k) = 2.0d0 / (tmp * en**2 * pm2**2)
        gmu(lim + k) = -gmu(k)
        gwt(lim + k) = gwt(k)
    end do

    return

end subroutine shd_util__computeGausQuad_new

subroutine shd_util__computeGausQuad(N, XA, WT)

    !╭──────────────────────────────────────────────────────────────╮!
    !   Generates the abscissas (x) and weights (W) for an N point
    ! Gauss-Legendre quadrature.  The XA are returned in this order:
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    !╰──────────────────────────────────────────────────────────────╯!

    integer, intent(in) :: N
    real, intent(inout) :: XA(*), WT(*)
    integer :: K, I, J, L
    real(RD_) :: x, XP, PL, PL1, PL2, DPL
    real(RD_), parameter :: TINY=3.0d-7

    K = N/2
    if (2*K .ne. N) stop 'Error [_computeGausQuad]: N must be even.'
    do J = 1, K
        x = dcos(PI_*(J - 0.25d0)/(N + 0.5d0))
        I = 0
100     continue
        PL1 = 1
        PL = x
        do L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*x*PL1 - (L - 1)*PL2)/L
        end do
        DPL = N*(x*PL - PL1)/(x*x - 1)
        XP = x
        x = XP - PL/DPL
        I = I + 1
        if (abs(x - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J) = -x
        XA(J + K) = x
        WT(J) = 2.0d0/((1.0d0 - x*x)*DPL*DPL)
        WT(J + K) = WT(J)
    end do

    return

end subroutine shd_util__computeGausQuad

subroutine shd_util__computeGausQuadDble(N, XA, WT)

    !╭──────────────────────────────────────────────────────────────╮!
    !   Generates the abscissas (x) and weights (W) for an N point
    ! Double-Gauss-Legendre quadrature. The XA are returned in this order
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    !╰──────────────────────────────────────────────────────────────╯!

    integer, intent(in) :: N
    real, intent(inout) :: XA(*), WT(*)

    integer :: K, I, J, L, N2
    real(RD_) :: x, XP, PL, PL1, PL2, DPL
    real(RD_), parameter :: TINY=3.0d-7

    N2 = N/2
    if (2*N2 .ne. N) stop 'Error [_computeGausQuadDble]: N must be even.'

    K = (N2 + 1)/2
    do J = 1, K
        x = cos(DPI_*(J - 0.25d0)/(N2 + 0.5d0))
        I = 0
100     continue
        PL1 = 1
        PL = x
        do L = 2, N2
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*x*PL1 - (L - 1)*PL2)/L
        end do
        DPL = N2*(x*PL - PL1)/(x*x - 1)
        XP = x
        x = XP - PL/DPL
        I = I + 1
        if (abs(x - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J + N2) = (1 + x)/2
        XA(N2 + 1 - J + N2) = (1 - x)/2
        WT(J + N2) = 1/((1 - x*x)*DPL*DPL)
        WT(N2 + 1 - J + N2) = WT(J + N2)
        XA(J) = -XA(J + N2)
        XA(N2 + 1 - J) = -XA(N2 + 1 - J + N2)
        WT(J) = WT(J + N2)
        WT(N2 + 1 - J) = WT(J + N2)
    end do

    return

end subroutine shd_util__computeGausQuadDble

subroutine shd_util__computePlanck(TEMP, UNITS, WAVENO, WAVELEN, PLANCK)

    !    Calculates the Planck blackbody radiance. If UNITS='T' then
    ! using brightness temperature units and the temperature is simply
    ! returned. If UNITS='B' then doing a band integration and the
    ! Planck blackbody radiance in [Watts /(meter^2 ster)] over a
    ! wavenumber range [cm^-1] is returned. Otherwise, the Planck
    ! blackbody radiance in [Watts /(meter^2 ster micron)] for a
    ! temperature in [Kelvins] at a wavelength in [microns] is returned.

    character(len=1), intent(in) :: UNITS
    real, intent(in) :: TEMP, WAVENO(2), WAVELEN
    real, intent(out) :: PLANCK
    real(RD_) :: X1, X2, F

    if (UNITS .eq. 'T') then
        PLANCK = TEMP
    else if (UNITS .eq. 'B') then
        if (TEMP .gt. 0.0) then
            X1 = 1.4388d0*WAVENO(1)/TEMP
            X2 = 1.4388d0*WAVENO(2)/TEMP
            call shd_util__integratePlanck(X1, X2, F)
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

end subroutine shd_util__computePlanck

subroutine shd_util__integratePlanck(X1, X2, F)

    ! Returns integral of x**3/(exp(x)-1) from x1 to x2.
    ! Accurate to better than 1 part in 10**9 for integral from
    ! 0 to infinity (pi**4/15).

    real(RD_), intent(in)  :: X1, X2
    real(RD_), intent(out) :: F
    real(RD_) :: C

    C = 1.0d0

    if (X1 .lt. C .and. X2 .lt. C) then
        F = integralLow(X2) - integralLow(X1)
    else if (X1 .lt. C .and. X2 .ge. C) then
        F = integralLow(C) - integralLow(X1) + integralHigh(C) - integralHigh(X2)
    else if (X1 .ge. C .and. X2 .ge. C) then
        F = integralHigh(X1) - integralHigh(X2)
    else
        stop 'Error [_integratePlanck]: X1 and X2 out of order.'
    end if

    return

end subroutine shd_util__integratePlanck

real(RD_) function integralLow(x)

    ! Integral of x**3/(exp(x)-1) from 0 to x.
    ! Accurate for x less than about 1.
    ! Uses Taylor series expansion around x=0.

    real(RD_) :: x
    integer :: N
    real(RD_) :: SUM, F, A(29)
    data A/ 0.0d0,&
        &   0.0d0,  0.333333333333333333333d0,&
        &-0.125d0,  0.016666666666666666667d0,&
        &   0.0d0, -0.0001984126984126984127d0,&
        &   0.0d0,  0.36743092298647854203d-5,&
        &   0.0d0, -0.75156325156325156325d-7,&
        &   0.0d0,  0.16059043836821614599d-8,&
        &   0.0d0, -0.35227934257916621232d-10,&
        &   0.0d0,  0.78720803121674581370d-12,&
        &   0.0d0, -0.17840422612224120352d-13,&
        &   0.0d0,  0.40886009791799259829d-15,&
        &   0.0d0, -0.94559508632959211871d-17,&
        &   0.0d0,  0.22036011313440918061d-18,&
        &   0.0d0, -0.51683202540046382743d-20,&
        &   0.0d0,  0.12188644964239543006d-20/

    SUM = A(4)*x**4
    F = x
    do N = 3, 29, 2
        F = x*x*F
        SUM = SUM + A(N)*F
    end do

    integralLow = SUM

    return

end function integralLow

real(RD_) function integralHigh(x)

    ! Integral of x**3/(exp(x)-1) from x to infinity.
    ! Accurate for x greater than about 1.

    real(RD_) :: x
    integer :: N
    real(RD_) :: SUM

    SUM = 0.0d0
    do N = 1, 15
        SUM = SUM + exp(-N*x)&
            &*(x**3/N + 3*x**2/N**2 + 6*x/N**3 + 6.0d0/N**4)
    end do

    integralHigh = SUM

    return

end function integralHigh

subroutine shd_util__ssort_i(x, y, N, kFlag)

    !╭──────────────────────────────────────────────────────────────╮!
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
    !   SSORT sorts array x and optionally makes the same interchanges in
    !   array y.  The array x may be sorted in increasing order or
    !   decreasing order.  A slightly modified quicksort algorithm is used.
    !
    !   Description of Parameters
    !      x - array of values to be sorted   (usually abscissas)
    !      y - array to be (optionally) carried along
    !      N - number of values in array x to be sorted
    !      kFlag - control parameter
    !            =  2  means sort x in increasing order and carry y along.
    !            =  1  means sort x in increasing order (ignoring y)
    !            = -1  means sort x in decreasing order (ignoring y)
    !            = -2  means sort x in decreasing order and carry y along.
    !
    !***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
    !                 for sorting with minimal storage, Communications of
    !                 the ACM, 12, 3 (1969), pp. 185-187.
    !╰──────────────────────────────────────────────────────────────╯!

    ! Scalar Arguments
    integer kFlag, N

    ! Array Arguments
    !real x(*), y(*)
    real, intent(inout) :: x(*)
    integer, intent(inout) :: y(*)

    ! Local Scalars
    ! real R, T, TT, TTY, TY
    real R, T, TT
    integer TY, TTY
    integer I, IJ, J, K, KK, L, M, NN

    ! Local Arrays
    integer IL(51), IU(51)

    ! External Subroutines
    ! Intrinsic Functions
    intrinsic abs, int

    !***First executable statement  SSORT
    NN = N
    if (NN .lt. 1) then
        stop 'Error [_ssort_i]: The number of values to be sorted is not positive.'
    end if

    KK = abs(kFlag)
    if (KK .ne. 1 .and. KK .ne. 2) then
        stop 'Error [_ssort_i]: The sort control parameter, K, is not 2, 1, -1, or -2.'
    end if

    ! Alter array x to get decreasing order if needed
    if (kFlag .le. -1) then
        do 10 I = 1, NN
            x(I) = -x(I)
10      continue
     end if

    if (KK .eq. 2) GO TO 100

    ! Sort x only
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
    T = x(IJ)

    ! If first element of array is greater than T, interchange with T
    if (x(I) .gt. T) then
        x(IJ) = x(I)
        x(I) = T
        T = x(IJ)
    end if
    L = J

    ! If last element of array is less than than T, interchange with T
    if (x(J) .lt. T) then
        x(IJ) = x(J)
        x(J) = T
        T = x(IJ)

        ! If first element of array is greater than T, interchange with T
        if (x(I) .gt. T) then
            x(IJ) = x(I)
            x(I) = T
            T = x(IJ)
        end if
    end if

    ! Find an element in the second half of the array which is smaller
    ! than T
40  L = L - 1
    if (x(L) .gt. T) GO TO 40

    ! Find an element in the first half of the array which is greater
    ! than T
50  K = K + 1
    if (x(K) .lt. T) GO TO 50

    ! Interchange these elements
    if (K .le. L) then
        TT = x(L)
        x(L) = x(K)
        x(K) = TT
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
    T = x(I + 1)
    if (x(I) .le. T) GO TO 80
    K = I

90  x(K + 1) = x(K)
    K = K - 1
    if (T .lt. x(K)) GO TO 90
    x(K + 1) = T
    GO TO 80

    ! Sort x and carry y along
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
    T = x(IJ)
    TY = y(IJ)

    ! If first element of array is greater than T, interchange with T
    if (x(I) .gt. T) then
        x(IJ) = x(I)
        x(I) = T
        T = x(IJ)
        y(IJ) = y(I)
        y(I) = TY
        TY = y(IJ)
    end if
    L = J

    ! If last element of array is less than T, interchange with T
    if (x(J) .lt. T) then
        x(IJ) = x(J)
        x(J) = T
        T = x(IJ)
        y(IJ) = y(J)
        y(J) = TY
        TY = y(IJ)

    ! If first element of array is greater than T, interchange with T
        if (x(I) .gt. T) then
            x(IJ) = x(I)
            x(I) = T
            T = x(IJ)
            y(IJ) = y(I)
            y(I) = TY
            TY = y(IJ)
        end if
    end if

    ! Find an element in the second half of the array which is smaller
    ! than T
130 L = L - 1
    if (x(L) .gt. T) GO TO 130

    ! Find an element in the first half of the array which is greater
    ! than T
140 K = K + 1
    if (x(K) .lt. T) GO TO 140

    ! Interchange these elements
    if (K .le. L) then
        TT = x(L)
        x(L) = x(K)
        x(K) = TT
        TTY = y(L)
        y(L) = y(K)
        y(K) = TTY
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
    T = x(I + 1)
    TY = y(I + 1)
    if (x(I) .le. T) GO TO 170
    K = I

180 x(K + 1) = x(K)
    y(K + 1) = y(K)
    K = K - 1
    if (T .lt. x(K)) GO TO 180
    x(K + 1) = T
    y(K + 1) = TY
    GO TO 170

    ! Clean up
190 if (kFlag .le. -1) then
        do 200 I = 1, NN
            x(I) = -x(I)
200     continue
    end if

    return

end subroutine shd_util__ssort_i

subroutine shd_util__ssort_r(x, y, N, kFlag)

    !╭──────────────────────────────────────────────────────────────╮!
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
    !   SSORT sorts array x and optionally makes the same interchanges in
    !   array y.  The array x may be sorted in increasing order or
    !   decreasing order.  A slightly modified quicksort algorithm is used.
    !
    !   Description of Parameters
    !      x - array of values to be sorted   (usually abscissas)
    !      y - array to be (optionally) carried along
    !      N - number of values in array x to be sorted
    !      kFlag - control parameter
    !            =  2  means sort x in increasing order and carry y along.
    !            =  1  means sort x in increasing order (ignoring y)
    !            = -1  means sort x in decreasing order (ignoring y)
    !            = -2  means sort x in decreasing order and carry y along.
    !
    !***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
    !                 for sorting with minimal storage, Communications of
    !                 the ACM, 12, 3 (1969), pp. 185-187.
    !╰──────────────────────────────────────────────────────────────╯!

    ! Scalar Arguments
    integer kFlag, N

    ! Array Arguments
    !real x(*), y(*)
    real, intent(inout) :: x(*)
    real, intent(inout) :: y(*)

    ! Local Scalars
    ! real R, T, TT, TTY, TY
    real R, T, TT
    integer TY, TTY
    integer I, IJ, J, K, KK, L, M, NN

    ! Local Arrays
    integer IL(51), IU(51)

    ! External Subroutines
    ! Intrinsic Functions
    intrinsic abs, int

    !***First executable statement  SSORT
    NN = N
    if (NN .lt. 1) then
        stop 'Error [_ssort_r]: The number of values to be sorted is not positive.'
    end if

    KK = abs(kFlag)
    if (KK .ne. 1 .and. KK .ne. 2) then
        stop 'Error [_ssort_r]: The sort control parameter, K, is not 2, 1, -1, or -2.'
    end if

    ! Alter array x to get decreasing order if needed
    if (kFlag .le. -1) then
        do 10 I = 1, NN
            x(I) = -x(I)
10      continue
     end if

    if (KK .eq. 2) GO TO 100

    ! Sort x only
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
    T = x(IJ)

    ! If first element of array is greater than T, interchange with T
    if (x(I) .gt. T) then
        x(IJ) = x(I)
        x(I) = T
        T = x(IJ)
    end if
    L = J

    ! If last element of array is less than than T, interchange with T
    if (x(J) .lt. T) then
        x(IJ) = x(J)
        x(J) = T
        T = x(IJ)

        ! If first element of array is greater than T, interchange with T
        if (x(I) .gt. T) then
            x(IJ) = x(I)
            x(I) = T
            T = x(IJ)
        end if
    end if

    ! Find an element in the second half of the array which is smaller
    ! than T
40  L = L - 1
    if (x(L) .gt. T) GO TO 40

    ! Find an element in the first half of the array which is greater
    ! than T
50  K = K + 1
    if (x(K) .lt. T) GO TO 50

    ! Interchange these elements
    if (K .le. L) then
        TT = x(L)
        x(L) = x(K)
        x(K) = TT
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
    T = x(I + 1)
    if (x(I) .le. T) GO TO 80
    K = I

90  x(K + 1) = x(K)
    K = K - 1
    if (T .lt. x(K)) GO TO 90
    x(K + 1) = T
    GO TO 80

    ! Sort x and carry y along
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
    T = x(IJ)
    TY = y(IJ)

    ! If first element of array is greater than T, interchange with T
    if (x(I) .gt. T) then
        x(IJ) = x(I)
        x(I) = T
        T = x(IJ)
        y(IJ) = y(I)
        y(I) = TY
        TY = y(IJ)
    end if
    L = J

    ! If last element of array is less than T, interchange with T
    if (x(J) .lt. T) then
        x(IJ) = x(J)
        x(J) = T
        T = x(IJ)
        y(IJ) = y(J)
        y(J) = TY
        TY = y(IJ)

    ! If first element of array is greater than T, interchange with T
        if (x(I) .gt. T) then
            x(IJ) = x(I)
            x(I) = T
            T = x(IJ)
            y(IJ) = y(I)
            y(I) = TY
            TY = y(IJ)
        end if
    end if

    ! Find an element in the second half of the array which is smaller
    ! than T
130 L = L - 1
    if (x(L) .gt. T) GO TO 130

    ! Find an element in the first half of the array which is greater
    ! than T
140 K = K + 1
    if (x(K) .lt. T) GO TO 140

    ! Interchange these elements
    if (K .le. L) then
        TT = x(L)
        x(L) = x(K)
        x(K) = TT
        TTY = y(L)
        y(L) = y(K)
        y(K) = TTY
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
    T = x(I + 1)
    TY = y(I + 1)
    if (x(I) .le. T) GO TO 170
    K = I

180 x(K + 1) = x(K)
    y(K + 1) = y(K)
    K = K - 1
    if (T .lt. x(K)) GO TO 180
    x(K + 1) = T
    y(K + 1) = TY
    GO TO 170

    ! Clean up
190 if (kFlag .le. -1) then
        do 200 I = 1, NN
            x(I) = -x(I)
200     continue
    end if

    return

end subroutine shd_util__ssort_r

function indexAt(x0, x) result(i)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Locates the position of 'x0' in the sorted array 'x' using binary search.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(kind=R__), intent(in) :: x0     ! value to find
    real(kind=R__), intent(in) :: x(:)   ! sorted array
    integer(kind=I__) :: i               ! return index
    integer :: low, high, mid
    integer :: n

    n = size(x)
    low = 1
    high = n

    do while (low .le. high)
        mid = low + (high - low) / 2  ! Avoids potential integer overflow

        if (x(mid) .eq. x0) then
            i = mid
            return
        else if (x(mid) .lt. x0) then
            low = mid + 1
        else
            high = mid - 1
        end if
    end do

    ! If the loop terminates, the value was not found.
    ! The 'low' index is the insertion point (the index of the first
    ! element greater than x0, or n+1 if all elements are smaller).
    i = low

end function indexAt

elemental pure function distKm(&
    lon1_deg, lat1_deg, h1_km, &
    lon2_deg, lat2_deg, h2_km) result(dist_km)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Calculate distance between two points, each point defined by
    !   (longitude [degree], latitude [degree], altitude [km])
    !   on oblate spheroid Earth
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(kind=R__), intent(in) :: lon1_deg, lat1_deg, h1_km
    real(kind=R__), intent(in) :: lon2_deg, lat2_deg, h2_km
    real(kind=R__) :: dist_km

    real(kind=R__) :: lon1, lat1, lon2, lat2
    real(kind=R__) :: sinlat1, coslat1, sinlon1, coslon1
    real(kind=R__) :: sinlat2, coslat2, sinlon2, coslon2
    real(kind=R__) :: N1, N2
    real(kind=R__) :: x1, y1, z1, x2, y2, z2

    ! convert degrees to radians
    lon1 = lon1_deg * deg2rad_
    lat1 = lat1_deg * deg2rad_
    lon2 = lon2_deg * deg2rad_
    lat2 = lat2_deg * deg2rad_

    ! compute sin/cos efficiently
    call sin_cos(lat1, sinlat1, coslat1)
    call sin_cos(lon1, sinlon1, coslon1)
    call sin_cos(lat2, sinlat2, coslat2)
    call sin_cos(lon2, sinlon2, coslon2)

    ! radii of curvature
    N1 = planetRadiusA_ / sqrt(1.0 - E2_ * sinlat1 * sinlat1)
    N2 = planetRadiusA_ / sqrt(1.0 - E2_ * sinlat2 * sinlat2)

    ! ECEF coordinates (km, Earth-Centered, Earth-Fixed coordinates)
    x1 = (N1 + h1_km) * coslat1 * coslon1
    y1 = (N1 + h1_km) * coslat1 * sinlon1
    z1 = ((1.0 - E2_) * N1 + h1_km) * sinlat1

    x2 = (N2 + h2_km) * coslat2 * coslon2
    y2 = (N2 + h2_km) * coslat2 * sinlon2
    z2 = ((1.0 - E2_) * N2 + h2_km) * sinlat2

    ! Euclidean distance
    dist_km = sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)

end function distKm

pure subroutine sin_cos(x, s, c)
    real(kind=R__), intent(in)  :: x
    real(kind=R__), intent(out) :: s, c
    s = sin(x)
    c = cos(x)
end subroutine sin_cos

subroutine shd_util__computeRayVect(&
    muRay, phiRay, &
    cx, cy, cz, &
    cx_inv, cy_inv, cz_inv, &
    ratio_x, ratio_y, &
    xyz_info)

    real(kind=RD_), intent(in)  :: muRay, phiRay
    real(kind=RD_), intent(out) :: cx, cy, cz
    real(kind=RD_), intent(out) :: cx_inv, cy_inv, cz_inv
    real(kind=RD_), intent(out) :: ratio_x, ratio_y

    real(kind=R__), intent(in), optional :: xyz_info(5)

    cx = dsqrt(dble(1.0d0 - muRay**2))*dcos(dble(phiRay + DPI_))
    cy = dsqrt(dble(1.0d0 - muRay**2))*dsin(dble(phiRay + DPI_))
    cz = dble(-muRay)

    if (abs(cx) .gt. 1.0d-7) then
        cx_inv = 1.0d0/cx
    else
        if (verbose_ .ge. 7) &
            write(6, *) 'Message [_computeRayVect]: cx too small, reset cx to 0, cx_inv to 1e8.'
        cx = 0.0d0
        cx_inv = 1.0d8
    end if

    if (abs(cy) .gt. 1.0d-7) then
        cy_inv = 1.0d0/cy
    else
        if (verbose_ .ge. 7) &
            write(6, *) 'Message [_computeRayVect]: cy too small, reset cy to 0, cy_inv to 1e8.'
        cy = 0.0d0
        cy_inv = 1.0d8
    end if

    if (abs(cz) .gt. 1.0d-7) then
        cz_inv = 1.0d0/cz
    else
        if (verbose_ .ge. 7) &
            write(6, *) 'Message [_computeRayVect]: cz too small, reset cz to 0, cz_inv to 1e8.'
        cz = 0.0d0
        cz_inv = 1.0d8
    end if

    if (present(xyz_info)) then
        call shd_util__computeDegKmRatio(&
            xyz_info(1), xyz_info(2), xyz_info(3), &
            xyz_info(4), xyz_info(5), &
            ratio_x, ratio_y)
    else
        ratio_x = 1.0d0
        ratio_y = 1.0d0
    end if

    return

end subroutine shd_util__computeRayVect

subroutine shd_util__computeDegKmRatio(&
    lon_deg, lat_deg, h_km, dlon, dlat, &
    ratio_x, ratio_y)

    real(kind=R__), intent(in)  :: lon_deg, lat_deg, h_km, dlon, dlat
    real(kind=RD_), intent(out) :: ratio_x, ratio_y

    ratio_x = distKm(lon_deg-0.5*dlon, lat_deg, h_km, &
                   lon_deg+0.5*dlon, lat_deg, h_km) / dlon

    ratio_y = distKm(lon_deg, lat_deg-0.5*dlat, h_km, &
                   lon_deg, lat_deg+0.5*dlat, h_km) / dlat

    return

end subroutine shd_util__computeDegKmRatio

subroutine shd_util__computeVecECEF(&
    lon_deg, lat_deg, muIn, phiIn, &
    vectECEF)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Initialize angle vector at center, converting local ENU to ECEF
    !   muIn is the SHDOM mu
    !   phiIn is the SHDOM phi
    !╰────────────────────────────────────────────────────────────────────────────╯!

    implicit none
    real(kind=R__), intent(in)  :: lon_deg, lat_deg, muIn, phiIn
    real(kind=R__), intent(out) :: vectECEF(3)
    real(kind=R__) :: mu, phi, lon, lat, Se, Sn
    real(kind=R__) :: sinLon, cosLon, sinLat, cosLat, sinZen
    real(kind=R__) :: east_x, east_y
    real(kind=R__) :: up_x, up_y, up_z
    real(kind=R__) :: north_x, north_y, north_z
    real(kind=R__) :: nx, ny, nz, norm

    ! convert MU (cos(zenith)) and PHI from SHDOM convension to legacy convension
    mu  = abs(muIn)
    phi = 1.5*PI_ - phiIn

    ! convert longitude and latitude from degree to radian
    lon = lon_deg * deg2rad_
    lat = lat_deg * deg2rad_

    ! Local ENU vector at center
    sinZen = sqrt(1.0 - mu*mu)

    Se = sinZen * sin(phi)
    Sn = sinZen * cos(phi)

    cosLon = cos(lon)
    sinLon = sin(lon)
    cosLat = cos(lat)
    sinLat = sin(lat)

    ! Ellipsoidal surface normal (UP vector)
    nx = cosLat*cosLon
    ny = cosLat*sinLon
    nz = sinLat/(1.0-E2_)
    norm = sqrt(nx*nx + ny*ny + nz*nz)

    up_x = nx / norm
    up_y = ny / norm
    up_z = nz / norm

    ! East vector (same as spherical)
    east_x = -sinLon
    east_y =  cosLon

    ! North = Up × East
    north_x = -up_z*east_y
    north_y =  up_z*east_x
    north_z =  up_x*east_y - up_y*east_x

    ! Convert local ENU vector -> global ECEF
    vectECEF(1) = Se*east_x  + Sn*north_x + mu*up_x
    vectECEF(2) = Se*east_y  + Sn*north_y + mu*up_y
    vectECEF(3) =              Sn*north_z + mu*up_z

    return

end subroutine shd_util__computeVecECEF

subroutine shd_util__computeDirENU(&
    vectECEF, &
    lon_deg, lat_deg, &
    muOut, phiOut)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Compute solar angles at any lat/lon on oblate ellipsoid
    !╰────────────────────────────────────────────────────────────────────────────╯!

    implicit none
    real(kind=R__), intent(in)  :: vectECEF(3)
    real(kind=R__), intent(in)  :: lon_deg, lat_deg
    real(kind=R__), intent(out) :: muOut, phiOut

    real(kind=R__) :: mu, phi, lon, lat, east_x, east_y
    real(kind=R__) :: up_x, up_y, up_z
    real(kind=R__) :: north_x, north_y, north_z
    real(kind=R__) :: nx, ny, nz, norm
    real(kind=R__) :: Se, Sn

    lon = lon_deg * deg2rad_
    lat = lat_deg * deg2rad_

    ! East (unchanged)
    east_x = -sin(lon)
    east_y =  cos(lon)

    ! Ellipsoidal UP normal
    nx = cos(lat)*cos(lon)/A2_
    ny = cos(lat)*sin(lon)/A2_
    nz = sin(lat)/B2_
    norm = sqrt(nx*nx + ny*ny + nz*nz)

    up_x = nx / norm
    up_y = ny / norm
    up_z = nz / norm

    ! North = Up × East
    north_x = -up_z*east_y
    north_y =  up_z*east_x
    north_z =  up_x*east_y - up_y*east_x

    ! Project Sun vector into local ENU (East, North, Up) frame
    Se =  east_x*vectECEF(1) +  east_y*vectECEF(2)
    Sn = north_x*vectECEF(1) + north_y*vectECEF(2) + north_z*vectECEF(3)
    mu =    up_x*vectECEF(1) +    up_y*vectECEF(2) +    up_z*vectECEF(3)

    if ((mu .ge. 0.0) .and. (mu .le. 1.0)) then
        phi = atan2(Se, Sn)

        ! convert to SHDOM convenction (opposite to sun direction)
        ! phi: from  [0 at north, positive clockwise] to
        !            [0 at east,  positive counterclockwise]
        !╭────────────────────────────────────────────────────────────────────────────╮!
        muOut  = -abs(mu)
        phiOut = 1.5*PI_ - phi

        if (phiOut < 0.0)     phiOut = phiOut + 2.0*PI_
        if (phiOut > 2.0*PI_) phiOut = phiOut - 2.0*PI_
        !╰────────────────────────────────────────────────────────────────────────────╯!
    else
        muOut  = 0.0
        phiOut = 0.0
    end if

    return

end subroutine shd_util__computeDirENU

subroutine shd_util__tickTock(message, count_start, more_digits)

    character(len=*), intent(inout) :: message
    integer(I__), intent(inout) :: count_start
    logical, intent(in), optional :: more_digits

    integer(I__) :: count_rate, count_max, count_end
    real(R__) :: elapsed_time

    call system_clock(count_end, count_rate, count_max)
    elapsed_time = (real(count_end) - real(count_start)) / real(count_rate)
    count_start = count_end

    if (elapsed_time .lt. 60.0) then
        if (present(more_digits)) then
            if (more_digits) then
                write(message, '(A,1X,F0.4,A)') trim(adjustl(message)), elapsed_time, ' sec.'
            else
                write(message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time, ' sec.'
            end if
        else
            write(message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time, ' sec.'
        end if
    else
        if (present(more_digits)) then
            if (more_digits) then
                write(message, '(A,1X,F0.4,A)') trim(adjustl(message)), elapsed_time/60.0, ' min.'
            else
                write(message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time/60.0, ' min.'
            end if
        else
            write(message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time/60.0, ' min.'
        end if
    end if

    return

end subroutine shd_util__tickTock

function change_range(angle, min_value, max_value) result(output_angle)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !***PURPOSE  Change the range of an angle to a specified range.
    !            This function is similar to Python's modulo operator
    !            but works with real numbers and allows for a specified
    !            minimum and maximum value.
    !***LIBRARY   SHD_UTIL
    !***CATEGORY  N6A2B
    !***TYPE      REAL
    !╰────────────────────────────────────────────────────────────────────────────╯!

    IMPLICIT NONE

    REAL, INTENT(IN) :: angle
    REAL, INTENT(IN) :: min_value
    REAL, INTENT(IN) :: max_value
    REAL             :: output_angle

    REAL :: current_range_span

    current_range_span = max_value - min_value

    ! This implementation assumes current_range_span is not zero.
    ! If current_range_span is zero, Python's % operator would raise a ZeroDivisionError.
    ! Fortran's MODULO intrinsic function would similarly result in a runtime error.
    ! The MODULO intrinsic behaves like Python's % operator for real numbers:
    ! result = x - FLOOR(x / y) * y
    output_angle = modulo(angle - min_value, current_range_span) + min_value

end function change_range

subroutine read_ice_refractive_index(filepath, wavelengths, n_real, n_imag, num_entries)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    !***PURPOSE  Read ice refractive index data from a file.
    !            The file is expected to contain three columns (delimted by commas only, no spaces):
    !            wavelength, real part of refractive index, and imaginary part.
    !***LIBRARY   SHD_UTIL
    !***CATEGORY  N6A2B
    !***TYPE      REAL
    !***KEYWORDS  Refractive Index, Ice, File I/O
    !***DESCRIPTION
    !   This subroutine reads ice refractive index data from a specified file.
    !   The file reference used is from Warren and Brandt (2008), doi:10.1029/2007JD009744
    !╰────────────────────────────────────────────────────────────────────────────╯!

    implicit none

    character(len=*),  intent(in)  :: filepath
    real, allocatable, intent(out) :: wavelengths(:)
    real, allocatable, intent(out) :: n_real(:)
    real, allocatable, intent(out) :: n_imag(:)
    integer, intent(out) :: num_entries

    integer :: unit_num, io_status, count
    character(len=charMaxLen_) :: line_buffer

    ! Find a free unit number
    call shd_util__get_free_unit(unit_num)

    ! First pass: Count the number of lines/entries
    open(unit=unit_num, file=filepath, status='OLD', action='READ', iostat=io_status)
    if (io_status .ne. 0) then
        write(*,*) 'Error opening file for counting: ', filepath
        num_entries = 0
        return
    end if

    count = 0
    do
        read(unit_num, fmt='(A)', iostat=io_status) line_buffer
        if (io_status .ne. 0) exit ! Exit on EOF or error
        count = count + 1
    end do
    close(unit_num)
    num_entries = count

    if (num_entries .eq. 0) then
        write(*,*) 'No data found in file: ', filepath
        return
    end if

    ! Allocate arrays
    allocate(wavelengths(num_entries), &
             n_real(num_entries), &
             n_imag(num_entries), &
             stat=io_status)

    if (io_status .ne. 0) then
        write(*,*) 'Error allocating memory for refractive index data.'
        num_entries = 0
        return
    end if

    ! Second pass: Read the data
    open(unit=unit_num, file=filepath, status='OLD', action='READ', iostat=io_status)
    if (io_status .ne. 0) then
        write(*,*) 'Error opening file for reading: ', filepath
        deallocate(wavelengths, n_real, n_imag)
        num_entries = 0
        return
    end if

    do count = 1, num_entries
        read(unit_num, *, iostat=io_status) wavelengths(count), n_real(count), n_imag(count)
        if (io_status .ne. 0) THEN
            write(*,*) 'Error reading data from file at line: ', count
            close(unit_num)
            deallocate(wavelengths, n_real, n_imag)
            num_entries = 0
            return
        end if
    end do

    close(unit_num)

end subroutine read_ice_refractive_index

subroutine shd_util__get_free_unit(unit)

    ! Subroutine to find a free file unit number

    integer, intent(out) :: unit

    logical :: unit_exists
    integer :: test_unit

    unit = 10 ! Starting unit number to check

    do
        inquire(unit=unit, opened=unit_exists)
        if (.not. unit_exists) then
            inquire(unit=unit, iostat=test_unit) ! Check if unit number is valid at all
            if (test_unit == 0) exit
        end if
        unit = unit + 1
        if (unit > 99) then ! Avoid standard units and very high numbers
            write(*,*) "Error [_get_free_unit]: Could not find a free file unit."
            stop
        end if
    end do

end subroutine shd_util__get_free_unit

end module shd_util
