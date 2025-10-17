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

    use globals_const, only: RD_, PI_, DPI_
    use globals_cfg,   only: isOpenBndryX_, isOpenBndryY_, isMpiX_, isMpiY_, &
                             is1DX_, is1DYor1D_

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

    public :: shd_util__tickTock

    private :: integralLow
    private :: integralHigh

    interface shd_util__ssort
        module procedure shd_util__ssort_i
        module procedure shd_util__ssort_r
    end interface

contains

subroutine shd_util__sumData(IFLAG, WT, N, data, SUMDATA)

    ! Sums an output array over the k-distribution.
    ! If IFLAG=0 then the output array is set to the input times the weight.

    integer, intent(in) :: IFLAG, N
    real, intent(in)  :: WT, data(N)
    real, intent(out) :: SUMDATA(N)
    integer :: I
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

end subroutine shd_util__sumData

subroutine shd_util__interpData(OLDNPTS, NPTS,&
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
                ! Go through all 12 edges of parent cell to find which one has point
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

end subroutine shd_util__interpData

subroutine shd_util__findCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
    &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
    &BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)

    !   Locates the grid cell in the tree structure containing the
    ! specified point (X0,Y0,Z0), and returns the cell pointer ICELL.
    ! First the base grid cell is found (using NX,NY,NZ, XGRID,YGRID,ZGRID),
    ! and the tree is traced down to find the smallest cell containing
    ! the point.

    integer NX, NY, NZ, NCELLS, ICELL, BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    real(RD_) :: X0, Y0, Z0
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real GRIDPOS(3, *)
    integer IL, IU, IM, IX, IY, IZ, NXC, NYC, IC, IPTR, DIR
    real(RD_) :: XDOMAIN, YDOMAIN

    ! Move point to within horizontal domain if periodic
    if (.not. (isOpenBndryX_ .or. isMpiX_)) then
        XDOMAIN = XGRID(NX + 1) - XGRID(1)
        if (X0 .lt. XGRID(1)) then
            X0 = X0 - XDOMAIN*(int((X0 - XGRID(1))/XDOMAIN) - 1)
        else if (X0 .gt. XGRID(NX + 1)) then
            X0 = X0 - XDOMAIN*int((X0 - XGRID(1))/XDOMAIN)
        end if
    end if
    if (.not. (isOpenBndryY_ .or. isMpiY_)) then
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
    if (is1DX_) then
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
    if (is1DYor1D_) then
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
    if (isOpenBndryX_) then
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
    if (isMpiX_ .and. .not. is1DX_) then
        NXC = NX - 1
        IX = min(IX, NXC)
    end if
    NYC = NY
    if (isOpenBndryY_) then
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
    if (isMpiY_ .and. .not. is1DYor1D_) then
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

end subroutine shd_util__findCell

subroutine shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
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
    real(RD_) :: XE, YE, ZE
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
            if (DIR .eq. 0) stop 'Error [_findNextCell]: No split direction.'
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

end subroutine shd_util__findNextCell

subroutine shd_util__interpTriLin(&
    &ICELL, GRIDPTR, GRIDPOS,&
    &X, Y, Z, ND, FIELD, ID1, NV, VALUES)

    !   Trilinearly interpolates the field (FIELD) given a grid cell
    ! pointer (ICELL), to compute the field values (VALUES) at the desired
    ! point (X,Y,Z). ND is the leading index of the FIELD (1 for a
    ! scalar field), ID1 is the first element of the vector to interpolate,
    ! and NV elements are output.

    integer, intent(in) :: ICELL, GRIDPTR(8, *), ND, ID1, NV
    real(RD_), intent(in) :: X, Y, Z
    real :: GRIDPOS(3, *), FIELD(ND, *), VALUES(NV)
    integer :: IPT1, IPT2, ID
    real(RD_) :: U, V, W, UM, VM, WM, DELX, DELY, DELZ

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
    !   Generates the abscissas (X) and weights (W) for an N point
    ! Gauss-Legendre quadrature.  The XA are returned in this order:
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    !╰──────────────────────────────────────────────────────────────╯!

    integer, intent(in) :: N
    real, intent(inout) :: XA(*), WT(*)
    integer :: K, I, J, L
    real(RD_) :: X, XP, PL, PL1, PL2, DPL
    real(RD_), parameter :: TINY=3.0d-7

    K = N/2
    if (2*K .ne. N) stop 'Error [_computeGausQuad]: N must be even.'
    do J = 1, K
        X = dcos(PI_*(J - 0.25d0)/(N + 0.5d0))
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

end subroutine shd_util__computeGausQuad

subroutine shd_util__computeGausQuadDble(N, XA, WT)

    !╭──────────────────────────────────────────────────────────────╮!
    !   Generates the abscissas (X) and weights (W) for an N point
    ! Double-Gauss-Legendre quadrature. The XA are returned in this order
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    !╰──────────────────────────────────────────────────────────────╯!

    integer, intent(in) :: N
    real, intent(inout) :: XA(*), WT(*)

    integer :: K, I, J, L, N2
    real(RD_) :: X, XP, PL, PL1, PL2, DPL
    real(RD_), parameter :: TINY=3.0d-7

    N2 = N/2
    if (2*N2 .ne. N) stop 'Error [_computeGausQuadDble]: N must be even.'

    K = (N2 + 1)/2
    do J = 1, K
        X = cos(DPI_*(J - 0.25d0)/(N2 + 0.5d0))
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

real(RD_) function integralLow(X)

    ! Integral of x**3/(exp(x)-1) from 0 to x.
    ! Accurate for x less than about 1.
    ! Uses Taylor series expansion around x=0.

    real(RD_) :: X
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

    SUM = A(4)*X**4
    F = X
    do N = 3, 29, 2
        F = X*X*F
        SUM = SUM + A(N)*F
    end do

    integralLow = SUM

    return

end function integralLow

real(RD_) function integralHigh(X)

    ! Integral of x**3/(exp(x)-1) from x to infinity.
    ! Accurate for x greater than about 1.

    real(RD_) :: X
    integer :: N
    real(RD_) :: SUM

    SUM = 0.0d0
    do N = 1, 15
        SUM = SUM + exp(-N*X)&
            &*(X**3/N + 3*X**2/N**2 + 6*X/N**3 + 6.0d0/N**4)
    end do

    integralHigh = SUM

    return

end function integralHigh

subroutine shd_util__ssort_i(X, Y, N, KFLAG)

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
    !╰──────────────────────────────────────────────────────────────╯!

    ! Scalar Arguments
    integer KFLAG, N

    ! Array Arguments
    !real X(*), Y(*)
    real, intent(inout) :: X(*)
    integer, intent(inout) :: Y(*)

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

    KK = abs(KFLAG)
    if (KK .ne. 1 .and. KK .ne. 2) then
        stop 'Error [_ssort_i]: The sort control parameter, K, is not 2, 1, -1, or -2.'
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

end subroutine shd_util__ssort_i

subroutine shd_util__ssort_r(X, Y, N, KFLAG)

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
    !╰──────────────────────────────────────────────────────────────╯!

    ! Scalar Arguments
    integer KFLAG, N

    ! Array Arguments
    !real X(*), Y(*)
    real, intent(inout) :: X(*)
    real, intent(inout) :: Y(*)

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

    KK = abs(KFLAG)
    if (KK .ne. 1 .and. KK .ne. 2) then
        stop 'Error [_ssort_r]: The sort control parameter, K, is not 2, 1, -1, or -2.'
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

end subroutine shd_util__ssort_r

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
                write (message, '(A,1X,F0.4,A)') trim(adjustl(message)), elapsed_time, ' sec.'
            else
                write (message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time, ' sec.'
            end if
        else
            write (message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time, ' sec.'
        end if
    else
        if (present(more_digits)) then
            if (more_digits) then
                write (message, '(A,1X,F0.4,A)') trim(adjustl(message)), elapsed_time/60.0, ' min.'
            else
                write (message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time/60.0, ' min.'
            end if
        else
            write (message, '(A,1X,F0.2,A)') trim(adjustl(message)), elapsed_time/60.0, ' min.'
        end if
    end if

    return

end subroutine shd_util__tickTock

FUNCTION change_range(angle, min_value, max_value) RESULT(output_angle)
    !╭──────────────────────────────────────────────────────────────╮!
    !***PURPOSE  Change the range of an angle to a specified range.
    !            This function is similar to Python's modulo operator
    !            but works with real numbers and allows for a specified
    !            minimum and maximum value.
    !***LIBRARY   SHD_UTIL
    !***CATEGORY  N6A2B
    !***TYPE      REAL

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
    output_angle = MODULO(angle - min_value, current_range_span) + min_value

END FUNCTION change_range

SUBROUTINE read_ice_refractive_index(filepath, wavelengths, n_real, n_imag, num_entries)
    !╭──────────────────────────────────────────────────────────────╮!
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

    IMPLICIT NONE

    CHARACTER(LEN=*), INTENT(IN) :: filepath
    REAL, ALLOCATABLE, INTENT(OUT) :: wavelengths(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: n_real(:)
    REAL, ALLOCATABLE, INTENT(OUT) :: n_imag(:)
    INTEGER, INTENT(OUT) :: num_entries

    INTEGER :: unit_num, io_status, count
    CHARACTER(LEN=200) :: line_buffer

    ! Find a free unit number
    CALL shd_util__get_free_unit(unit_num)

    ! First pass: Count the number of lines/entries
    OPEN(UNIT=unit_num, FILE=filepath, STATUS='OLD', ACTION='READ', IOSTAT=io_status)
    IF (io_status /= 0) THEN
        WRITE(*,*) 'Error opening file for counting: ', filepath
        num_entries = 0
        RETURN
    END IF

    count = 0
    DO
        READ(unit_num, FMT='(A)', IOSTAT=io_status) line_buffer
        IF (io_status /= 0) EXIT ! Exit on EOF or error
        count = count + 1
    END DO
    CLOSE(unit_num)
    num_entries = count

    IF (num_entries == 0) THEN
        WRITE(*,*) 'No data found in file: ', filepath
        RETURN
    END IF

    ! Allocate arrays
    ALLOCATE(wavelengths(num_entries), n_real(num_entries), n_imag(num_entries), STAT=io_status)
    IF (io_status /= 0) THEN
        WRITE(*,*) 'Error allocating memory for refractive index data.'
        num_entries = 0
        RETURN
    END IF

    ! Second pass: Read the data
    OPEN(UNIT=unit_num, FILE=filepath, STATUS='OLD', ACTION='READ', IOSTAT=io_status)
    IF (io_status /= 0) THEN
        WRITE(*,*) 'Error opening file for reading: ', filepath
        DEALLOCATE(wavelengths, n_real, n_imag)
        num_entries = 0
        RETURN
    END IF

    DO count = 1, num_entries
        READ(unit_num, *, IOSTAT=io_status) wavelengths(count), n_real(count), n_imag(count)
        IF (io_status /= 0) THEN
            WRITE(*,*) 'Error reading data from file at line: ', count
            CLOSE(unit_num)
            DEALLOCATE(wavelengths, n_real, n_imag)
            num_entries = 0
            RETURN
        END IF
    END DO

    CLOSE(unit_num)

CONTAINS
    ! Subroutine to find a free file unit number
    SUBROUTINE shd_util__get_free_unit(unit)
        INTEGER, INTENT(OUT) :: unit
        LOGICAL :: unit_exists
        INTEGER :: test_unit
        unit = 10 ! Starting unit number to check
        DO
            INQUIRE(UNIT=unit, OPENED=unit_exists)
            IF (.NOT. unit_exists) THEN
                INQUIRE(UNIT=unit, IOSTAT=test_unit) ! Check if unit number is valid at all
                IF (test_unit == 0) EXIT
            END IF
            unit = unit + 1
            IF (unit > 99) THEN ! Avoid standard units and very high numbers
                WRITE(*,*) "Could not find a free file unit."
                STOP
            END IF
        END DO
    END SUBROUTINE shd_util__get_free_unit

END SUBROUTINE read_ice_refractive_index

end module shd_util
