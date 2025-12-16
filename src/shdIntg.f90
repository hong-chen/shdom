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

module shd_intg

    use globals_const, only: RD_, DPI_, &
                             FRAC13_, FRAC112_
    use globals_cfg,   only: verbose_, accSolv_, doGridPtsPrp_, &
                             isMpi_, is1DX_, is1DYx_, is1D_, &
                             isPol_

    use globals_obj,   only: rtm_obj_, xyz_Rtm_, &
                                       ext_Rtm_, &
                                      brdf_Rtm_, &
                             cell_obj_, i8Grid_Cell_, &
                                        i6Cell_Cell_, &
                                          flag_Cell_, &
                            shdo_obj_, &
                             src_obj_, &
                             prp_obj_, &
                             sfc_obj_, brdf_Sfc_

    use shd_sfc,  only: shd_sfc__interpSfc

    use shd_init, only: shd_init__makeTopRad, &
                        shd_init__makeBotRadUL, &
                        shd_init__makeBotRadVL, &
                        shd_init__makeBotRadU, &
                        shd_init__makeBotRadV

    use shd_grid, only: shd_grid__makeSweepOrder, &
                        shd_grid__findBndryPts

    use shd_util, only: shd_util__findNextCell, &
                        shd_util__computeRayVect

    use shd_shdo, only: shd_shdo__transSHtoDO_, &
                        shd_shdo__transSHtoDO, &
                        shd_shdo__transDOtoSH_, &
                        shd_shdo__transDOtoSH

    use shd_mpi,  only: shd_mpi__findBndryPts, &
                        shd_mpi__calBndryRad, &
                        shd_mpi__abort


    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_intg__integratePath

    private :: shd_intg__backIntegrate3D
    private :: shd_intg__backIntegrate3D_
    private :: shd_intg__backIntegrate2D
    private :: shd_intg__backIntegrate1D

contains

subroutine shd_intg__integratePath(&
    &brdf_Shdo,&
    &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &CMU1, CPHI1, CMU2, CPHI2, WSAVE, FFTFLAG,&
    &fluxDir, fluxTot,&
    &WORK, sweepOrder, GRIDRAD,&
    &iradI_Shdo, radI_Shdo, isrcJ_Shdo, srcJ_Shdo)

    !   Performs the path integrations through the medium specified by
    ! the extinction (ext_Rtm_) and source function (srcJ_Shdo) in
    ! spherical harmonic space.  The source function is transformed to
    ! discrete ordinates one zenith angle at a time.  The integrations
    ! are done from the top boundary (with no reflection) down to the
    ! bottom boundary for the downwelling angles, and then up from the
    ! bottom to the top after computing the surface reflection/emission.
    ! After the path integrations produce the discrete angle radiance
    ! it is transformed and added to output radiance spherical harmonic
    ! expansion.  The hemispheric flux fluxTot is also computed.
    ! The discrete ordinates in shdo_obj_%mu, shdo_obj_%phi must be a complete set with the
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

    integer(I__), intent(inout) :: iradI_Shdo(*)
    real(R__),    intent(inout) ::  radI_Shdo(shdo_obj_%Nstokes, *)

    integer(I__), intent(inout) :: isrcJ_Shdo(*)
    real(R__),    intent(inout) ::  srcJ_Shdo(shdo_obj_%Nstokes, *)

    integer MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS
    integer BCPTR(MAXNBC, 2)
    integer sweepOrder(rtm_obj_%Npts, *)
    logical FFTFLAG(shdo_obj_%Nmu)
    real brdf_Shdo(shdo_obj_%Nstokes, shdo_obj_%Nmu/2, shdo_obj_%maxNphiA, shdo_obj_%Nstokes, shdo_obj_%Nmu/2, shdo_obj_%maxNphiA)
    real CMU1(*), CPHI1(*), CMU2(*), CPHI2(*), WSAVE(*)
    real BCRAD(shdo_obj_%Nstokes, MAXBCRAD)
    real fluxDir(rtm_obj_%Npts), fluxTot(2, rtm_obj_%Npts)
    real GRIDRAD(shdo_obj_%Nstokes, rtm_obj_%Npts), WORK(shdo_obj_%Nstokes, shdo_obj_%maxNphiA, rtm_obj_%Npts)

    integer I, I1, K, NR, iPhiDO, iMuDO, IBC, iAngDO, iUpDn, NptsOld
    integer jOctant, iPts, ip, iorder, iCellPts, iCorner
    real A, minTrans
    save NptsOld

    ! Set the minimum transmission for the cell tracing (1 for single cell)
    minTrans = 1.00

    ! Make the new grid cell/point sweeping order if need to
    if (rtm_obj_%Npts .ne. NptsOld) then
        call shd_grid__makeSweepOrder(sweepOrder, GRIDRAD)
    end if

    call shd_mpi__findBndryPts(rtm_obj_%flag_bc, rtm_obj_%flag_ip, sweepOrder,&
        &rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz)

    if (rtm_obj_%Npts .ne. NptsOld) then
        ! Make the pointers to the grid points on the top and bottom boundaries
        call shd_grid__findBndryPts(shdo_obj_%Nang, MAXNBC, MAXBCRAD, &
            &rtm_obj_%z(1), rtm_obj_%z(rtm_obj_%Nz),&
            &NTOPPTS, NBOTPTS, BCPTR)
        if (sfc_obj_%type(1:1) .eq. 'V') then
            ! If this is a variable surface then bilinearly interpolate
            !   the surface temperature and reflection parameters to the
            !   bottom grid points.
            call shd_sfc__interpSfc(NBOTPTS, BCPTR(1, 2), xyz_Rtm_, &
                src_obj_%type, src_obj_%waveNum, src_obj_%avgWaveLen, src_obj_%unit, &
                brdf_Sfc_, brdf_Rtm_)
        end if
    end if
    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(src_obj_%type, src_obj_%skyRad, src_obj_%waveNum, src_obj_%avgWaveLen,&
        &src_obj_%unit, NTOPPTS, shdo_obj_%Nstokes, BCRAD(1, 1))

    ! Zero the flux arrays and set temporary radiance to -1 for not valid
    do I = 1, rtm_obj_%Npts
        fluxTot(1, I) = 0.0
        fluxTot(2, I) = 0.0
        GRIDRAD(1, I) = -1.0
    end do

    ! Zero the output spherical harmonic radiance
    NR = iradI_Shdo(rtm_obj_%Npts + 1)
    do I = 1, NR
        radI_Shdo(:, I) = 0.0
    end do

    ! write(*, *) 'Path integration:'
    ! Loop over the zenith angles (downwelling angles must be first)
    iAngDO = 1
    MU_DO_LOOP: do iMuDO = 1, shdo_obj_%Nmu

        if (iMuDO .eq. shdo_obj_%Nmu/2 + 1) then
            ! For the first upwelling angle, make the bottom boundary radiances
            !   for the Lambertian surfaces.  Compute the upwelling bottom
            !   radiances using the downwelling fluxTot.
            if (sfc_obj_%type .eq. 'FL') then
                call shd_init__makeBotRadUL(NBOTPTS, BCPTR(1, 2),&
                    &fluxDir, fluxTot, src_obj_%type,&
                    &src_obj_%waveNum, src_obj_%avgWaveLen, src_obj_%unit, shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
            else if (sfc_obj_%type .eq. 'VL') then
                call shd_init__makeBotRadVL(NBOTPTS, BCPTR(1, 2),&
                    &fluxDir, fluxTot, src_obj_%type,&
                    &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
            end if
        end if

        if (isPol_) then
            call shd_shdo__transSHtoDO(rtm_obj_%Npts, shdo_obj_%Nstokes, prp_obj_%NlegSet, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA(iMuDO),&
                &iMuDO, FFTFLAG, CMU1, CPHI1, WSAVE,&
                &isrcJ_Shdo, srcJ_Shdo, WORK)
        else
            call shd_shdo__transSHtoDO_(rtm_obj_%Npts, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA(iMuDO),&
                &iMuDO, FFTFLAG, CMU1, CPHI1, WSAVE,&
                &isrcJ_Shdo, srcJ_Shdo, WORK)
        end if

        ! Enforce the discrete ordinate I source function to be non-negative
        do I = 1, rtm_obj_%Npts
            do iPhiDO = 1, shdo_obj_%NphiA(iMuDO)
                WORK(1, iPhiDO, I) = max(0.0, WORK(1, iPhiDO, I))
            end do
        end do

        ! Loop over the azimuthal angles
        PHI_DO_LOOP: do iPhiDO = 1, shdo_obj_%NphiA(iMuDO)

            if (shdo_obj_%mu(iMuDO) .lt. 0.0) then

                ! For downward ordinates, initialize the top boundary radiances
                iUpDn = 1
                do IBC = 1, NTOPPTS
                    I = BCPTR(IBC, 1)
                    GRIDRAD(:, I) = BCRAD(:, IBC)
                end do

            else

                ! Initialize the bottom boundary radiances
                iUpDn = 2
                if (.not. sfc_obj_%isLambertian) then

                    ! If a Lambertian surface, use the radiances computed above
                    !   otherwise, compute the radiance for this angle
                    !   by integrating over the stored downwelling radiances.
                    ! If the surface reflection is uniform use the precomputed
                    !   BRDF matrix to do the integration.
                    if (sfc_obj_%isUniform) then
                        call shd_init__makeBotRadU(NBOTPTS, BCPTR(1, 2),&
                            &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt, iMuDO, iPhiDO,&
                            &src_obj_%type, src_obj_%avgWaveLen, src_obj_%mu0, src_obj_%phi0, fluxDir,&
                            &brdf_Shdo,&
                            &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
                    else
                        call shd_init__makeBotRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                            &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt,&
                            &shdo_obj_%mu(iMuDO), shdo_obj_%phi(iMuDO, iPhiDO),&
                            &src_obj_%type, src_obj_%avgWaveLen, src_obj_%mu0, src_obj_%phi0, fluxDir,&
                            &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
                    end if
                end if

                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    GRIDRAD(:, I) = BCRAD(:, IBC + NTOPPTS)
                end do

            end if

            ! Do backward path integration from each grid point not done yet
            if (is1D_) then

                call shd_intg__backIntegrate1D(rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz, shdo_obj_%maxNphiA,&
                    &sweepOrder, shdo_obj_%mu(iMuDO), shdo_obj_%phi(iMuDO, iPhiDO), minTrans,&
                    &ext_Rtm_, shdo_obj_%Nstokes, WORK, iPhiDO, GRIDRAD)

            else if (is1DYx_ .and. (.not. isMpi_)) then

                call shd_intg__backIntegrate2D(rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz, shdo_obj_%maxNphiA,&
                    &sweepOrder, shdo_obj_%mu(iMuDO), shdo_obj_%phi(iMuDO, iPhiDO), minTrans,&
                    &ext_Rtm_, shdo_obj_%Nstokes, WORK, iPhiDO, GRIDRAD)

            else
                if (isPol_) then

                    call shd_intg__backIntegrate3D(rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz, shdo_obj_%maxNphiA,&
                        &sweepOrder,&
                        &shdo_obj_%mu(iMuDO), shdo_obj_%phi(iMuDO, iPhiDO), minTrans, rtm_obj_%flag_bc, rtm_obj_%flag_ip,&
                        &ext_Rtm_, shdo_obj_%Nstokes, WORK, iPhiDO, GRIDRAD)

                else

                    call shd_intg__backIntegrate3D_(rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz, shdo_obj_%maxNphiA,&
                        &sweepOrder,&
                        &shdo_obj_%mu(iMuDO), shdo_obj_%phi(iMuDO, iPhiDO), minTrans, rtm_obj_%flag_bc, rtm_obj_%flag_ip,&
                        &ext_Rtm_, WORK, iPhiDO, GRIDRAD)

                end if
            end if

            ! Transfer the radiance for the (iMuDO,iPhiDO) discrete ordinate from
            !   GRIDRAD to the work array.
            do I = 1, rtm_obj_%Npts
                if (GRIDRAD(1, I) .lt. -accSolv_) then
                    if (verbose_ .ge. 4) write(6, *) 'Warning [_integratePath]: Grid point has no value',&
                        &I, shdo_obj_%mu(iMuDO), shdo_obj_%phi(iMuDO, iPhiDO), GRIDRAD(1, I), rtm_obj_%x(1), rtm_obj_%y(1),&
                        &xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I)
                    ! stop
                    GRIDRAD(1, I) = 0.0
                else
                    WORK(:, iPhiDO, I) = GRIDRAD(:, I)
                    GRIDRAD(1, I) = -1.0
                end if
            end do

            ! Do the sum for the hemispheric fluxTot for this angle
            A = abs(shdo_obj_%mu(iMuDO))*shdo_obj_%wgt(iMuDO, iPhiDO)
            do I = 1, rtm_obj_%Npts
                fluxTot(iUpDn, I) = fluxTot(iUpDn, I) + A*WORK(1, iPhiDO, I)
            end do

            ! For a general BRDF surface save the downwelling radiance
            !   at the bottom boundary points.
            if ((.not. sfc_obj_%isLambertian) .and. (iMuDO .le. shdo_obj_%Nmu/2)) then
                I1 = NTOPPTS + NBOTPTS*iAngDO
                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    BCRAD(:, IBC + I1) = WORK(:, iPhiDO, I)
                end do
            end if

            iAngDO = iAngDO + 1
        end do PHI_DO_LOOP

        ! Transform the radiance to spherical harmonics for this zenith angle
        if (isPol_) then
            call shd_shdo__transDOtoSH(rtm_obj_%Npts, shdo_obj_%Nstokes, prp_obj_%NlegSet, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA(iMuDO),&
                &iMuDO, FFTFLAG, CMU2, CPHI2, WSAVE,&
                &iradI_Shdo, WORK, radI_Shdo)
        else
            call shd_shdo__transDOtoSH_(rtm_obj_%Npts, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA(iMuDO),&
                &iMuDO, FFTFLAG, CMU2, CPHI2, WSAVE,&
                &iradI_Shdo, WORK, radI_Shdo)
        end if

    end do MU_DO_LOOP

    NptsOld = rtm_obj_%Npts

    return

end subroutine shd_intg__integratePath

subroutine shd_intg__backIntegrate3D(&
    &NX, NY, NZ, NA,&
    &sweepOrder,&
    &MU, PHI, minTrans, BCFLAG, IPFLAG,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)

    !   Sweeps through the spatial grid computing radiance vectors for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  The sweeping
    ! order through the grid points is in sweepOrder. If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is traced
    ! back to a face containing known (valid) radiances.  As each cell is
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! at the end of the ray path is interpolated from the known grid point
    ! radiances on the face.  The resulting radiances are put in GRIDRAD.

    integer NX, NY, NZ, NA, NSTOKES, KANG
    integer BCFLAG, IPFLAG
    integer sweepOrder(rtm_obj_%Npts, *)
    real MU, PHI, minTrans
    real EXTINCT(rtm_obj_%Npts), SOURCE(NSTOKES, NA, rtm_obj_%Npts)
    real GRIDRAD(NSTOKES, rtm_obj_%Npts)
    integer bitX, bitY, bitZ, iOctant, jOctant
    integer iCellPts, iCell, iCellNext, iPts, iFace
    integer I1, I2, I3, I4, iOpp, jFace, kFace, IC, IZ
    integer faceGrid(4, 6), faceOpp(6)
    integer iOrder, iCorner, jOctantOrder3(8), jOctantOrder2(8)
    logical isValidRad, isValidFace, is1DX, is1DY
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV, ratioX, ratioY
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS, path
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(RD_) :: tau, transCell, absCell, trans
    real(RD_) :: RAD(NSTOKES), RAD0(NSTOKES)
    data faceGrid/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, faceOpp/2, 1, 4, 3, 6, 5/
    data jOctantOrder3/1, 3, 5, 7, 2, 4, 6, 8/, jOctantOrder2/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    call shd_util__computeRayVect(dble(MU), dble(PHI), &
        CX, CY, CZ, CXINV, CYINV, CZINV, &
        ratioX, ratioY)

    ! Setup the sweeping direction and the gridpoint location
    ! bitc is 0 for positive direction trace back ray, 1 for negative.
    if (CX .lt. 0.0) then
        bitX = 1
    else
        bitX = 0
    end if
    if (CY .lt. 0.0) then
        bitY = 1
    else
        bitY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        bitZ = 1
    else if (CZ .gt. 1.0e-3) then
        bitZ = 0
    else
        stop 'Error [_backIntegrate3D]: Bad MU.'
    end if

    ! iOctant is the octant the ray come from or the discrete ordinate goes to.
    iOctant = 1 + bitX + 2*bitY + 4*bitZ
    if (is1DYx_) then
        jOctant = jOctantOrder2(iOctant)
    else
        jOctant = jOctantOrder3(iOctant)
    end if

    IZ = (NZ - 1)*(1 - bitZ) + 1
    ! Sweep through all the grid points
    do iOrder = 1, rtm_obj_%Npts
        iCellPts = ishft(sweepOrder(iOrder, jOctant), -3)
        iCorner = ibits(sweepOrder(iOrder, jOctant), 0, 3) + 1
        iPts = i8Grid_Cell_(iCorner, iCellPts)

        if (isMpi_) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if  ((xyz_Rtm_(3, iPts) .lt. rtm_obj_%z(IZ) .and. bitZ .eq. 0) .or.&
                &(xyz_Rtm_(3, iPts) .gt. rtm_obj_%z(IZ) .and. bitZ .eq. 1)) then
                IZ = IZ + 2*bitZ - 1
                call shd_mpi__calBndryRad(BCFLAG, IPFLAG, jOctant, IZ,&
                    &NX, NY, NZ,&
                    &NA,&
                    &MU, PHI, EXTINCT, NSTOKES, SOURCE,&
                    &KANG, GRIDRAD)
            end if
        end if

        if (GRIDRAD(1, iPts) .lt. 0.0) then

            ! Start off at the grid point
            iCell = iCellPts
            trans = 1.0
            EXT1 = EXTINCT(iPts)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, iPts)
            XE = xyz_Rtm_(1, iPts)
            YE = xyz_Rtm_(2, iPts)
            ZE = xyz_Rtm_(3, iPts)

            ! Loop until finding a face with known radiances
            isValidRad = .false.
            do while (.not. isValidRad)
                ! Make sure current cell is valid
                if (iCell .le. 0) then
                    write(6, *) 'Error [_backIntegrate3D]: iCell=0 ',&
                        &MU, PHI, iPts, iCellNext, xe, ye, ze,&
                        &xyz_Rtm_(1, iPts), xyz_Rtm_(2, iPts), xyz_Rtm_(3, iPts)
                    stop
                end if
                is1DX = btest(int(flag_Cell_(iCell)), 0)
                is1DY = btest(int(flag_Cell_(iCell)), 1)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                iOpp = i8Grid_Cell_(9 - iOctant, iCell)
                ! Get the distances to the 3 planes and select the closest
                if (is1DX) then
                    SOX = 1.0e20
                else
                    SOX = (xyz_Rtm_(1, iOpp) - XE)*CXINV
                end if
                if (is1DY) then
                    SOY = 1.0e20
                else
                    SOY = (xyz_Rtm_(2, iOpp) - YE)*CYINV
                end if
                SOZ = (xyz_Rtm_(3, iOpp) - ZE)*CZINV

                SO = min(SOX, SOY, SOZ)
                if (SO .lt. -EPS) then
                    write(6, *) 'Error [_backIntegrate3D]: SO<0 ',&
                        &MU, PHI, XE, YE, ZE, SOX, SOY, SOZ, SO, iPts, iCell
                    stop
                end if
                ! Compute the coordinates of the cell exitting location
                XE = XE + SO*CX
                YE = YE + SO*CY
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                if (SOX .le. SOZ .and. SOX .le. SOY) then
                    iFace = 2 - bitX
                    jFace = 1
                else if (SOY .le. SOZ) then
                    iFace = 4 - bitY
                    jFace = 2
                else
                    iFace = 6 - bitZ
                    jFace = 3
                end if
                ! Get the next cell to go to
                iCellNext = i6Cell_Cell_(iFace, iCell)
                if (iCellNext .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, iFace, jFace, iCell,&
                        &iCellNext)
                end if

                ! Get the face grid pointers
                if (i6Cell_Cell_(iFace, iCell) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    kFace = iFace
                    IC = iCell
                else
                    ! If going to smaller face then use next face (more accurate)
                    kFace = faceOpp(iFace)
                    IC = iCellNext
                end if
                I1 = i8Grid_Cell_(faceGrid(1, kFace), IC)
                I2 = i8Grid_Cell_(faceGrid(2, kFace), IC)
                I3 = i8Grid_Cell_(faceGrid(3, kFace), IC)
                I4 = i8Grid_Cell_(faceGrid(4, kFace), IC)
                ! Compute the face interpolation factors
                if (jFace .eq. 1) then
                    U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
                    if (is1DY) then
                        V = 0.5
                    else
                        V = (YE - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I2) - xyz_Rtm_(2, I1))
                    end if
                else if (jFace .eq. 2) then
                    U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
                    if (is1DX) then
                        V = 0.5
                    else
                        V = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
                    end if
                else
                    if (is1DY) then
                        U = 0.5
                    else
                        U = (YE - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I3) - xyz_Rtm_(2, I1))
                    end if
                    if (is1DX) then
                        V = 0.5
                    else
                        V = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
                    end if
                end if
                if (u .lt. -0.0001 .or. v .lt. -0.0001 .or.&
                    &u .gt. 1.0001 .or. v .gt. 1.0001) then
                    write(*, *) '_backIntegrate3D: u,v<0 or u,v>1: ',&
                        &mu, phi, xe, ye, ze, u, v, iFace, jFace, kFace,&
                        &xyz_Rtm_(1, i1), xyz_Rtm_(2, i1), xyz_Rtm_(3, i1)
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (iCellNext .gt. 0) then
                    if (jFace .eq. 1) then
                        XE = xyz_Rtm_(1, i8Grid_Cell_(iOctant, iCellNext))
                    else if (jFace .eq. 2) then
                        YE = xyz_Rtm_(2, i8Grid_Cell_(iOctant, iCellNext))
                    else
                        ZE = xyz_Rtm_(3, i8Grid_Cell_(iOctant, iCellNext))
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
                tau = EXT*SO
                if (tau .ge. 0.5) then
                    transCell = exp(-tau)
                    absCell = 1.0 - transCell
                else
                    absCell = tau*(1.0 - 0.5*tau*&
                        &(1.0 - FRAC13_*tau*(1 - 0.25*tau)))
                    transCell = 1.0 - absCell
                end if
                if (tau .le. 2.0) then
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
                    if (tau .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0/tau
                        if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                    end if
                    SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                        &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                        &*(1.0 - 2.0/tau + 2.0*transCell/absCell))
                end if
                SRC(1) = max(SRC(1), 0.0d0)

                ! Add in the cell radiance and update the transmission.
                ! If this is a boundary or the transmission is below the
                !   cutoff and we have a valid radiance then set the flag
                !   to stop the tracing and add in the interpolated face radiance.
                RAD(:) = RAD(:) + trans*SRC(:)*absCell
                trans = trans*transCell
                !if (RAD(1) < -1.0e-5) then
                if (RAD(1) < -accSolv_) then
                    write(*, '(A,3(1X,F6.4),8(1X,E12.5))') &
                        &'_backIntegrate3D: RAD<0 ',&
                        &xe, ye, ze, ext0, ext1, so, tau, trans, src, absCell, rad(:)
                    call shd_mpi__abort('The end')
                end if
                isValidFace = (GRIDRAD(1, I1) .ge. -0.1 .and. GRIDRAD(1, I2) .ge. -0.1&
                    &.and. GRIDRAD(1, I3) .ge. -0.1 .and. GRIDRAD(1, I4) .ge. -0.1)
                if (iCellNext .le. 0 .or.&
                    &(trans .le. minTrans .and. isValidFace)) then
                    if (isValidFace) then
                        isValidRad = .true.
                        RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)&
                            &+ F3*GRIDRAD(:, I3) + F4*GRIDRAD(:, I4)
                        RAD(:) = RAD(:) + trans*RAD0(:)
                    else
                        write(*, *) 'Error [_backIntegrate3D]: iCellNext=0: ', xe, ye, ze, iCell
                        write(*, *) iz, jOctant, mu, phi
                        write(*, '(i7,4(1x,f7.4))') i1, GRIDRAD(1, I1),&
                            &xyz_Rtm_(1, i1), xyz_Rtm_(2, i1), xyz_Rtm_(3, i1)
                        write(*, '(i7,4(1x,f7.4))') i2, GRIDRAD(1, I2),&
                            &xyz_Rtm_(1, i2), xyz_Rtm_(2, i2), xyz_Rtm_(3, i2)
                        write(*, '(i7,4(1x,f7.4))') i3, GRIDRAD(1, I3),&
                            &xyz_Rtm_(1, i3), xyz_Rtm_(2, i3), xyz_Rtm_(3, i3)
                        write(*, '(i7,4(1x,f7.4))') i4, GRIDRAD(1, I4),&
                            &xyz_Rtm_(1, i4), xyz_Rtm_(2, i4), xyz_Rtm_(3, i4)
                        write(*, '(i7,8x,3(1x,f7.4))') iPts,&
                            &xyz_Rtm_(1, iPts), xyz_Rtm_(2, iPts), xyz_Rtm_(3, iPts)
                        stop
                    end if
                else
                    EXT1 = EXT0
                    SRCEXT1(:) = SRCEXT0(:)
                    iCell = iCellNext
                end if
            end do
            GRIDRAD(:, iPts) = RAD(:)

        end if
    end do

    return
end subroutine shd_intg__backIntegrate3D

subroutine shd_intg__backIntegrate3D_( &
    NX, NY, NZ, NA, &
    sweepOrder, &
    MU, PHI, minTrans, BCFLAG, IPFLAG, &
    EXTINCT, SOURCE, KANG, GRIDRAD)

    ! Description:
    !   Sweeps through the spatial grid computing radiances for the
    !   discrete ordinate direction (MU,PHI) by integrating the source
    !   function and extinction backward from each point. The sweeping
    !   order through the grid points is in sweepOrder. If a point is
    !   already done (nonnegative in GRIDRAD) then it is skipped.
    !   For each grid point a ray opposite to the ordinate direction is traced
    !   back to a face containing known (valid) radiances. As each cell is
    !   crossed (usually only one) the integration of the source function
    !   across the cell is done and added to the radiance. The radiance
    !   at the end of the ray path is interpolated from the known grid point
    !   radiances on the face. The resulting radiances are put in GRIDRAD.
    !
    ! Input:
    !   NX, NY, NZ: integer, base grid size in X, Y and Z
    !   NA        : integer, total number of angles (MU and PHI)
    !   rtm_obj_%Npts : integer, total number of grid points
    !   cell_obj_%Npts : integer, total number of cells
    !   i8Grid_Cell_(8,cell_obj_%Npts) : integer array, pointer to 8 grid points belonging to a cell.
    !                       Order of pointers is based on binary digits Z-Y-X:
    !                       iOctant = 1 + bitX + 2*bitY + 4*bitZ, where
    !                       bitc is 0 for negative side and 1 for positive side.
    !                       1: (-X, -Y, -Z), 2: (+X, -Y, -Z), 3: (-X, +Y, -Z), 4: (+X, +Y, -Z)
    !                       5: (-X, -Y, +Z), 6: (+X, -Y, +Z), 7: (-X, +Y, +Z), 8: (+X, +Y, +Z)
    !   i6Cell_Cell_(6,cell_obj_%Npts): integer array, pointer to 6 neighboring cells of a cell.
    !                       Order of pointers are 1:-X, 2:+X, 3:-Y, 4:+Y, 5:-Z, 6:+Z.
    !                       Is 0 if there is no neighbor cell (boundary).
    !                       Is positive if there is just one cell across face
    !                       of current cell (unique neighbor).  Is negative
    !                       (-index) if there is more than one cell across
    !                       the face of current cell.  In this case the
    !                       cell pointed to is the parent cell that has a
    !                       has same or larger face adjoining current cell.
    !   i2Tree_Cell_(2,cell_obj_%Npts) : integer array, pointer to parent and children cells in tree structure.
    !                       First is pointer to parent cell, 0 for a base cell.
    !                       Second is pointer to first of the two children cells.
    !                       This pointer is 0 if there are no children cells.
    !                       The children cells are in successive order, with
    !                       the negative side first. See flag_Cell_ for direction
    !                       of division of split.
    !   flag_Cell_(cell_obj_%Npts) : integer array, a collection of binary flags for the cells:
    !                       bit  For bit set
    !                       0:  Do independent pixel for X
    !                       1:  Do independent pixel for Y
    !                       2,3: Direction cell is split (0=no split, 1=X, 2=Y, 3=Z)
    !   xyz_Rtm_(3,rtm_obj_%Npts)    : real array, X,Y,Z coordinate values of grid point locations
    !   rtm_obj_%x, rtm_obj_%y, rtm_obj_%z: real array, grid locations for each dimension.
    !                        To facilitate the periodic horizontal boundaries,
    !                        there are NX+1 values in rtm_obj_%x, NY+1 in rtm_obj_%y,
    !                        NZ values in rtm_obj_%z. Currently the grid locations
    !                        are evenly spaced in X and Y.
    !   sweepOrder(rtm_obj_%Npts,*)   : integer array, grid point sweeping order

    integer(kind=I__), intent(in) :: NX, NY, NZ, NA, KANG, BCFLAG, IPFLAG
    integer(kind=I__), intent(in) :: sweepOrder(rtm_obj_%Npts, *)
    real(kind=R__)   , intent(in) :: MU, PHI, minTrans
    real(kind=R__)   , intent(in) :: EXTINCT(rtm_obj_%Npts), SOURCE(NA, rtm_obj_%Npts)

    real   (kind=R__), intent(inout) :: GRIDRAD(rtm_obj_%Npts)

    integer(kind=I__) :: bitX, bitY, bitZ, iOctant, jOctant
    integer(kind=I__) :: iCellPts, iCell, iCellNext, iPts, iFace
    integer(kind=I__) :: I1, I2, I3, I4, iOpp, jFace, kFace, IC, IZ
    integer(kind=I__) :: iOrder, iCorner
    logical :: isValidRad, isValidFace, is1DX, is1DY

    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV, ratioX, ratioY
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS, path
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
    real(RD_) :: EXT0P, SRCEXT0P
    real(RD_) :: tau, transCell, absCell, RAD, RAD0, trans

    integer(kind=I__) :: faceGrid(4, 6), faceOpp(6), jOctantOrder_1D(8), jOctantOrder_3D(8)
    data faceGrid /  &
        1, 3, 5, 7, &
        2, 4, 6, 8, &
        1, 2, 5, 6, &
        3, 4, 7, 8, &
        1, 2, 3, 4, &
        5, 6, 7, 8/
    data faceOpp /2, 1, 4, 3, 6, 5/
    data jOctantOrder_1D /1, 3, 1, 3, 2, 4, 2, 4/
    data jOctantOrder_3D /1, 3, 5, 7, 2, 4, 6, 8/

    integer :: iwhile, ibreak

    EPS = 1.0e-3*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    call shd_util__computeRayVect(dble(MU), dble(PHI), &
        CX, CY, CZ, CXINV, CYINV, CZINV, &
        ratioX, ratioY)

    ! Setup the sweeping direction and the gridpoint location
    ! bitc is 0 for positive direction trace back ray, 1 for negative.
    if (CX .lt. 0.0) then
        bitX = 1
    else
        bitX = 0
    end if

    if (CY .lt. 0.0) then
        bitY = 1
    else
        bitY = 0
    end if

    if (CZ .lt. -1.0e-3) then
        bitZ = 1
    else if (CZ .gt. 1.0e-3) then
        bitZ = 0
    else
        stop 'Error [_backIntegrate3D_]: Bad MU.'
    end if

    ! iOctant is the octant the ray come from or the discrete ordinate goes to.
    iOctant = 1 + bitX + 2*bitY + 4*bitZ
    if (is1DYx_) then
        jOctant = jOctantOrder_1D(iOctant)
    else
        jOctant = jOctantOrder_3D(iOctant)
    end if

    IZ = (NZ - 1)*(1 - bitZ) + 1

    ! Sweep through all the grid points
    do iOrder = 1, rtm_obj_%Npts

        iCellPts  = ishft(sweepOrder(iOrder, jOctant), -3)
        iCorner = ibits(sweepOrder(iOrder, jOctant), 0, 3) + 1
        iPts = i8Grid_Cell_(iCorner, iCellPts)

        ! this is for mpi (e.g., shdMpi.f90, not OMP)
        !╭────────────────────────────────────────────────╮!
        if (isMpi_) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if  ((xyz_Rtm_(3, iPts) .lt. rtm_obj_%z(IZ) .and. bitZ .eq. 0) .or. &
                 (xyz_Rtm_(3, iPts) .gt. rtm_obj_%z(IZ) .and. bitZ .eq. 1)) then
                IZ = IZ + 2*bitZ - 1
                call shd_mpi__calBndryRad(BCFLAG, IPFLAG, jOctant, IZ, &
                    NX, NY, NZ, &
                    NA, &
                    MU, PHI, EXTINCT, 1, SOURCE, &
                    KANG, GRIDRAD)
            end if
        end if
        !╰────────────────────────────────────────────────╯!

        if (GRIDRAD(iPts) .lt. 0.0) then

            ! Start off at the grid point
            iCell = iCellPts
            RAD = 0.0
            trans = 1.0
            EXT1 = EXTINCT(iPts)
            SRCEXT1 = EXT1*SOURCE(KANG, iPts)
            XE = xyz_Rtm_(1, iPts)
            YE = xyz_Rtm_(2, iPts)
            ZE = xyz_Rtm_(3, iPts)

            ! Loop until finding a face with known radiances
            isValidRad = .false.
            do while (.not. isValidRad)

                ! Make sure current cell is valid
                if (iCell .le. 0) then
                    write(6, *) 'Error [_backIntegrate3D_]: iCell=0 ', &
                        MU, PHI, iPts, iCellNext, XE, YE, ZE, &
                        xyz_Rtm_(1, iPts), xyz_Rtm_(2, iPts), xyz_Rtm_(3, iPts)
                    stop
                end if

                is1DX = btest(int(flag_Cell_(iCell)), 0)
                is1DY = btest(int(flag_Cell_(iCell)), 1)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                iOpp = i8Grid_Cell_(9 - iOctant, iCell)

                ! Get the distances to the 3 planes and select the closest
                !╭────────────────────────────────────────────────╮!
                if (is1DX) then
                    SOX = 1.0e20
                else
                    SOX = (xyz_Rtm_(1, iOpp) - XE)*CXINV
                end if

                if (is1DY) then
                    SOY = 1.0e20
                else
                    SOY = (xyz_Rtm_(2, iOpp) - YE)*CYINV
                end if

                SOZ = (xyz_Rtm_(3, iOpp) - ZE)*CZINV

                SO = min(SOX, SOY, SOZ)
                !╰────────────────────────────────────────────────╯!

                if (SO .lt. -EPS) then
                    write(6, *) 'Error [_backIntegrate3D_]: SO<0 ', &
                        MU, PHI, XE, YE, ZE, SOX, SOY, SOZ, SO, iPts, iCell
                    stop
                end if

                ! Compute the coordinates of the cell exitting location
                XE = XE + SO*CX
                YE = YE + SO*CY
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                if (SOX .le. SOZ .and. SOX .le. SOY) then
                    iFace = 2 - bitX
                    jFace = 1
                else if (SOY .le. SOZ) then
                    iFace = 4 - bitY
                    jFace = 2
                else
                    iFace = 6 - bitZ
                    jFace = 3
                end if

                ! Get the next cell to go to
                iCellNext = i6Cell_Cell_(iFace, iCell)
                if (iCellNext .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, iFace, jFace, iCell, &
                        iCellNext)
                end if

                ! Get the face grid pointers
                if (i6Cell_Cell_(iFace, iCell) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    kFace = iFace
                    IC = iCell
                else
                    ! If going to smaller face then use next face (more accurate)
                    kFace = faceOpp(iFace)
                    IC = iCellNext
                end if

                I1 = i8Grid_Cell_(faceGrid(1, kFace), IC)
                I2 = i8Grid_Cell_(faceGrid(2, kFace), IC)
                I3 = i8Grid_Cell_(faceGrid(3, kFace), IC)
                I4 = i8Grid_Cell_(faceGrid(4, kFace), IC)

                ! Compute the face interpolation factors
                if (jFace .eq. 1) then
                    U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
                    if (is1DY) then
                        V = 0.5
                    else
                        V = (YE - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I2) - xyz_Rtm_(2, I1))
                    end if
                else if (jFace .eq. 2) then
                    U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
                    if (is1DX) then
                        V = 0.5
                    else
                        V = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
                    end if
                else
                    if (is1DY) then
                        U = 0.5
                    else
                        U = (YE - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I3) - xyz_Rtm_(2, I1))
                    end if
                    if (is1DX) then
                        V = 0.5
                    else
                        V = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
                    end if
                end if

                if (U .lt. -0.0001 .or. V .lt. -0.0001 .or. &
                    U .gt.  1.0001 .or. V .gt.  1.0001) then
                    write(*, *) '_backIntegrate3D_: U,V<0 or U,V>1: ', &
                        MU, PHI, XE, YE, ZE, U, V, iFace, jFace, kFace, &
                        xyz_Rtm_(1, I1), xyz_Rtm_(2, I1), xyz_Rtm_(3, I1)
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (iCellNext .gt. 0) then
                    if (jFace .eq. 1) then
                        XE = xyz_Rtm_(1, i8Grid_Cell_(iOctant, iCellNext))
                    else if (jFace .eq. 2) then
                        YE = xyz_Rtm_(2, i8Grid_Cell_(iOctant, iCellNext))
                    else
                        ZE = xyz_Rtm_(3, i8Grid_Cell_(iOctant, iCellNext))
                    end if
                end if

                ! Interpolate extinction and source function at face intersection
                F1 = (1.0 - U)*(1.0 - V)
                F2 = (1.0 - U)*       V
                F3 =        U *(1.0 - V)
                F4 =        U *       V
                EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + &
                       F3*EXTINCT(I3) + F4*EXTINCT(I4)
                ! Correctly interpolate source using extinction*source
                SRCEXT0 = (F1*SOURCE(KANG, I1)*EXTINCT(I1) + &
                           F2*SOURCE(KANG, I2)*EXTINCT(I2) + &
                           F3*SOURCE(KANG, I3)*EXTINCT(I3) + &
                           F4*SOURCE(KANG, I4)*EXTINCT(I4))
                ! Compute the cell radiance: integration of the source function
                EXT = 0.5*(EXT0 + EXT1)

                tau = EXT*SO
                if (tau .ge. 0.5) then
                    transCell = exp(-tau)
                    absCell = 1.0 - transCell
                else
                    absCell = tau*(1.0 - 0.5*tau*&
                        &(1.0 - FRAC13_*tau*(1.0 - 0.25*tau)))
                    transCell = 1.0 - absCell
                end if

                if (tau .le. 2.0) then
                    if (EXT .eq. 0.0) then
                        SRC = 0.0
                    else
                        ! Linear extinction, linear source*extinction, to first order
                        SRC = (0.5*(SRCEXT0 + SRCEXT1)&
                            &+ FRAC112_*(EXT0*SRCEXT1 - EXT1*SRCEXT0)*SO)/EXT
                    end if
                else
                    ! Combined first order expansion and constant extinction formula
                    ! SRC = 0.5/EXT *( SRCEXT0+SRCEXT1&
                    !     &+ (EXT0*SRCEXT1-EXT1*SRCEXT0)*SO&
                    !     &*(1-2/tau+2*transCell/absCell)/tau )
                    EXT0P = EXT0
                    SRCEXT0P = SRCEXT0
                    if (tau .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0/tau
                        if (EXT0 .gt. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
                    end if
                    SRC = 1.0/(EXT0P + EXT1)*(SRCEXT0P + SRCEXT1&
                        &+ (EXT0P*SRCEXT1 - EXT1*SRCEXT0P)*2.0/(EXT0P + EXT1)&
                        &*(1.0 - 2.0/tau + 2.0*transCell/absCell))
                end if
                SRC = max(SRC, 0.0d0)

                ! Add in the cell radiance and update the transmission.
                ! If this is a boundary or the transmission is below the
                ! cutoff and we have a valid radiance then set the flag
                ! to stop the tracing and add in the interpolated face radiance.
                RAD = RAD + trans*SRC*absCell
                trans = trans*transCell
                !if (RAD < -1.0e-5) then
                if (RAD < -accSolv_) then
                    write(*, '(A,3(1X,F6.4),8(1X,E12.5))') &
                        &'_backIntegrate3D_: RAD<0 ',&
                        &XE, YE, ZE, EXT0, EXT1, SO, tau, trans, SRC, absCell, RAD
                    call shd_mpi__abort('The end')
                end if

                isValidFace = (GRIDRAD(I1) .ge. -0.1 .and. GRIDRAD(I2) .ge. -0.1&
                      &.and. GRIDRAD(I3) .ge. -0.1 .and. GRIDRAD(I4) .ge. -0.1)
                if ((iCellNext .le. 0) .or. &
                    ((trans .le. minTrans) .and. isValidFace)) then
                    if (isValidFace) then
                        isValidRad = .true.
                        RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2)&
                            &+ F3*GRIDRAD(I3) + F4*GRIDRAD(I4)
                        RAD = RAD + trans*RAD0
                    else
                        write(*, *) 'Error [_backIntegrate3D_]: iCellNext=0: ', xe, ye, ze, iCell
                        write(*, *) iz, jOctant, mu, phi
                        write(*, '(i7,4(1x,f7.4))') i1, GRIDRAD(I1),&
                            &xyz_Rtm_(1, i1), xyz_Rtm_(2, i1), xyz_Rtm_(3, i1)
                        write(*, '(i7,4(1x,f7.4))') i2, GRIDRAD(I2),&
                            &xyz_Rtm_(1, i2), xyz_Rtm_(2, i2), xyz_Rtm_(3, i2)
                        write(*, '(i7,4(1x,f7.4))') i3, GRIDRAD(I3),&
                            &xyz_Rtm_(1, i3), xyz_Rtm_(2, i3), xyz_Rtm_(3, i3)
                        write(*, '(i7,4(1x,f7.4))') i4, GRIDRAD(I4),&
                            &xyz_Rtm_(1, i4), xyz_Rtm_(2, i4), xyz_Rtm_(3, i4)
                        write(*, '(i7,8x,3(1x,f7.4))') iPts,&
                            &xyz_Rtm_(1, iPts), xyz_Rtm_(2, iPts), xyz_Rtm_(3, iPts)
                        stop
                    end if
                else
                    EXT1 = EXT0
                    SRCEXT1 = SRCEXT0
                    iCell = iCellNext
                end if
            end do

            GRIDRAD(iPts) = RAD

        end if

    end do

    return

end subroutine shd_intg__backIntegrate3D_

subroutine shd_intg__backIntegrate2D(&
    &NX, NY, NZ, NA,&
    &sweepOrder, MU, PHI, minTrans,&
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

    integer NX, NY, NZ, NA, KANG, NSTOKES
    integer sweepOrder(rtm_obj_%Npts, *)
    real MU, PHI, minTrans
    real EXTINCT(rtm_obj_%Npts), SOURCE(NSTOKES, NA, rtm_obj_%Npts)
    real GRIDRAD(NSTOKES, rtm_obj_%Npts)
    integer bitX, bitZ, iOctant, jOctant
    integer iCellPts, iCell, iCellNext, iPts, iFace
    integer I1, I2, iOpp, jFace, kFace, IC, IZ
    integer faceGrid(2, 6), faceOpp(6)
    integer iOrder, iCorner, jOctantOrder(8)
    logical :: isValidRad, isValidFace, is1DX
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV, ratioX, ratioY
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOZ, EPS
    real(RD_) :: U, F1, F2
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(RD_) :: tau, transCell, absCell, trans
    real(RD_) :: RAD(NSTOKES), RAD0(NSTOKES)
    data faceGrid/1, 5, 2, 6, 0, 0, 0, 0, 1, 2, 5, 6/
    data faceOpp/2, 1, 4, 3, 6, 5/
    data jOctantOrder/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    call shd_util__computeRayVect(dble(MU), dble(PHI), &
        CX, CY, CZ, CXINV, CYINV, CZINV, &
        ratioX, ratioY)

    ! Setup the sweeping direction and the gridpoint location
    ! bitc is 0 for positive direction trace back ray, 1 for negative.
    if (CX .lt. 0.0) then
        bitX = 1
    else
        bitX = 0
    end if
    if (CZ .lt. -1.0e-3) then
        bitZ = 1
    else if (CZ .gt. 1.0e-3) then
        bitZ = 0
    else
        stop 'Error [_backIntegrate2D]: Bad MU.'
    end if
    ! iOctant is the octant the ray come from or the discrete ordinate goes to.
    iOctant = 1 + bitX + 4*bitZ
    jOctant = jOctantOrder(iOctant)

    ! Sweep through all the grid points
    do iOrder = 1, rtm_obj_%Npts
        iCellPts = ishft(sweepOrder(iOrder, jOctant), -3)
        iCorner = ibits(sweepOrder(iOrder, jOctant), 0, 3) + 1
        iPts = i8Grid_Cell_(iCorner, iCellPts)
        if (GRIDRAD(1, iPts) .lt. 0.0) then

            ! Start off at the grid point
            iCell = iCellPts
            trans = 1.0
            EXT1 = EXTINCT(iPts)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, iPts)
            XE = xyz_Rtm_(1, iPts)
            YE = xyz_Rtm_(2, iPts)
            ZE = xyz_Rtm_(3, iPts)

            ! Loop until finding a face with known radiances
            isValidRad = .false.
            do while (.not. isValidRad)
                ! Make sure current cell is valid
                if (iCell .le. 0) then
                    write(6, *) 'Error [_backIntegrate2D]: iCell=0 ',&
                        &MU, PHI, iPts, iCellNext
                    stop
                end if
                is1DX = btest(int(flag_Cell_(iCell)), 0)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                iOpp = i8Grid_Cell_(9 - iOctant, iCell)
                ! Get the distances to the 3 planes and select the closest
                if (is1DX) then
                    SOX = 1.0e20
                else
                    SOX = (xyz_Rtm_(1, iOpp) - XE)*CXINV
                end if
                SOZ = (xyz_Rtm_(3, iOpp) - ZE)*CZINV
                SO = min(SOX, SOZ)
                if (SO .lt. -EPS) then
                    write(6, *) 'Error [_backIntegrate2D]: SO<0 ',&
                        &MU, PHI, XE, ZE, SOX, SOZ, SO, iPts, iCell
                    stop
                end if
                ! Compute the coordinates of the cell exitting location
                XE = XE + SO*CX
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                if (SOX .le. SOZ) then
                    iFace = 2 - bitX
                    jFace = 1
                else
                    iFace = 6 - bitZ
                    jFace = 3
                end if
                ! Get the next cell to go to
                iCellNext = i6Cell_Cell_(iFace, iCell)
                if (iCellNext .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, iFace, jFace, iCell,&
                        &iCellNext)
                end if

                ! Get the face grid pointers
                if (i6Cell_Cell_(iFace, iCell) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    kFace = iFace
                    IC = iCell
                else
                    ! If going to smaller face then use next face (more accurate)
                    kFace = faceOpp(iFace)
                    IC = iCellNext
                end if
                I1 = i8Grid_Cell_(faceGrid(1, kFace), IC)
                I2 = i8Grid_Cell_(faceGrid(2, kFace), IC)
                ! Compute the face interpolation factors
                if (jFace .eq. 1) then
                    U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I2) - xyz_Rtm_(3, I1))
                else
                    if (is1DX) then
                        U = 0.5
                    else
                        U = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
                    end if
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (iCellNext .gt. 0) then
                    if (jFace .eq. 1) then
                        XE = xyz_Rtm_(1, i8Grid_Cell_(iOctant, iCellNext))
                    else
                        ZE = xyz_Rtm_(3, i8Grid_Cell_(iOctant, iCellNext))
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
                tau = EXT*SO
                if (tau .ge. 0.5) then
                    transCell = exp(-tau)
                    absCell = 1.0 - transCell
                else
                    absCell = tau*(1.0 - 0.5*tau*&
                        &(1.0 - FRAC13_*tau*(1 - 0.25*tau)))
                    transCell = 1.0 - absCell
                end if
                if (tau .le. 2.0) then
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
                    if (tau .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0/tau
                        if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                    end if
                    SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                        &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                        &*(1 - 2/tau + 2*transCell/absCell))
                end if
                SRC(1) = max(SRC(1), 0.0d0)

                ! Add in the cell radiance and update the transmission.
                ! If this is a boundary or the transmission is below the
                ! cutoff and we have a valid radiance then set the flag
                ! to stop the tracing and add in the interpolated face radiance
                RAD(:) = RAD(:) + trans*SRC(:)*absCell
                trans = trans*transCell
                isValidFace = (GRIDRAD(1, I1) .ge. -0.1 .and. GRIDRAD(1, I2) .ge. -0.1)
                if (iCellNext .le. 0 .or.&
                    &(trans .le. minTrans .and. isValidFace)) then
                    if (isValidFace) then
                        isValidRad = .true.
                        RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)
                        RAD(:) = RAD(:) + trans*RAD0(:)
                    else
                        write(*, *) 'Error [_backIntegrate2D]: iCellNext=0: ', xe, ye, ze, iCell
                        write(*, *) iz, jOctant, mu, phi
                        write(*, '(i7,4(1x,f7.4))') i1, GRIDRAD(1, I1),&
                            &xyz_Rtm_(1, i1), xyz_Rtm_(2, i1), xyz_Rtm_(3, i1)
                        write(*, '(i7,4(1x,f7.4))') i2, GRIDRAD(1, I2),&
                            &xyz_Rtm_(1, i2), xyz_Rtm_(2, i2), xyz_Rtm_(3, i2)
                        write(*, '(i7,8x,3(1x,f7.4))') iPts,&
                            &xyz_Rtm_(1, iPts), xyz_Rtm_(2, iPts), xyz_Rtm_(3, iPts)
                        stop
                    end if
                else
                    EXT1 = EXT0
                    SRCEXT1(:) = SRCEXT0(:)
                    iCell = iCellNext
                end if
            end do

            GRIDRAD(:, iPts) = RAD(:)

        end if

    end do

    return

end subroutine shd_intg__backIntegrate2D

subroutine shd_intg__backIntegrate1D(&
    &NX, NY, NZ, NA,&
    &sweepOrder, MU, PHI, minTrans,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Sweeps through the spatial grid computing radiances for the
    !   discrete ordinate direction (MU,PHI) by integrating the source
    !   function and extinction backward from each point. If a point is
    !   already done (nonnegative in GRIDRAD) then it is skipped.
    !   For each grid point a ray opposite to the ordinate direction is traced
    !   back to a face containing known (valid) radiances. As each cell is
    !   crossed (usually only one) the integration of the source function
    !   across the cell is done and added to the radiance. The radiance
    !   at the end of the ray path is interpolated from the known grid point
    !   radiances on the face. The resulting radiances are put in GRIDRAD.
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=I__), intent(in) :: NX, NY, NZ, NA, KANG, NSTOKES
    integer(kind=I__), intent(in) :: sweepOrder(rtm_obj_%Npts, *)
    real(kind=R__)   , intent(in) :: MU, PHI, minTrans
    real(kind=R__)   , intent(in) :: EXTINCT(rtm_obj_%Npts), SOURCE(NSTOKES, NA, rtm_obj_%Npts)

    real(kind=R__), intent(inout) :: GRIDRAD(NSTOKES, rtm_obj_%Npts)

    real(kind=RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV, ratioX, ratioY
    real(kind=RD_) :: XE, YE, ZE
    real(kind=RD_) :: SO, EPS
    real(kind=RD_) :: EXT, EXT0, EXT1
    real(kind=RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(kind=RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(kind=RD_) :: tau, transCell, absCell, trans
    real(kind=RD_) :: RAD(NSTOKES), RAD0(NSTOKES)

    integer(kind=I__) :: bitZ, iOctant, jOctant
    integer(kind=I__) :: iCellPts, iCell, iCellNext, iPts, iFace, I1, iOpp
    integer(kind=I__) :: iOrder, iCorner
    logical :: isValidRad

    integer(kind=I__) :: faceGrid(6), jOctantOrder(8)
    data faceGrid /0, 0, 0, 0, 1, 5/
    data jOctantOrder /1, 1, 1, 1, 2, 2, 2, 2/

    EPS = 1.0e-3*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    call shd_util__computeRayVect(dble(MU), dble(PHI), &
        CX, CY, CZ, CXINV, CYINV, CZINV, &
        ratioX, ratioY)

    ! Setup the sweeping direction and the gridpoint location
    ! bitc is 0 for positive direction trace back ray, 1 for negative.
    if (CZ .lt. -1.0e-3) then
        bitZ = 1
    else if (CZ .gt. 1.0e-3) then
        bitZ = 0
    else
        stop 'Error [_backIntegrate1D]: Bad MU.'
    end if

    ! iOctant is the octant the ray come from or the discrete ordinate goes to.
    iOctant = 1 + 4*bitZ
    jOctant = jOctantOrder(iOctant)

    ! Sweep through all the grid points
    do iOrder = 1, rtm_obj_%Npts
        iCellPts = ishft(sweepOrder(iOrder, jOctant), -3)
        iCorner = ibits(sweepOrder(iOrder, jOctant), 0, 3) + 1
        iPts = i8Grid_Cell_(iCorner, iCellPts)
        if (GRIDRAD(1, iPts) .lt. 0.0) then
            ! Start off at the grid point
            iCell = iCellPts
            trans = 1.0
            EXT1 = EXTINCT(iPts)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, iPts)
            XE = xyz_Rtm_(1, iPts)
            YE = xyz_Rtm_(2, iPts)
            ZE = xyz_Rtm_(3, iPts)

            ! Loop until finding a face with known radiances
            isValidRad = .false.
            do while (.not. isValidRad)
                ! Make sure current cell is valid
                if (iCell .le. 0) then
                    write(6, *) 'Error [_backIntegrate1D]: iCell=0 ',&
                        &MU, PHI, iPts, iCellNext
                    stop
                end if
                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                iOpp = i8Grid_Cell_(9 - iOctant, iCell)
                ! Get the distances to the 3 planes and select the closest
                SO = (xyz_Rtm_(3, iOpp) - ZE)*CZINV
                if (SO .lt. -EPS) then
                    write(6, *) 'Error [_backIntegrate1D]: SO<0 ',&
                        &MU, PHI, ZE, SO, iPts, iCell
                    stop
                end if
                ! Compute the coordinates of the cell exitting location
                ZE = ZE + SO*CZ

                ! Get the intersection face number (i.e. neighptr index)
                iFace = 6 - bitZ
                ! Get the next cell to go to
                iCellNext = i6Cell_Cell_(iFace, iCell)

                ! Get the face grid pointers
                I1 = i8Grid_Cell_(faceGrid(iFace), iCell)
                ! Get the location coordinate (does the boundary wrapping)
                if (iCellNext .gt. 0) then
                    ZE = xyz_Rtm_(3, i8Grid_Cell_(iOctant, iCellNext))
                end if

                ! Get extinction and source function at face intersection
                EXT0 = EXTINCT(I1)
                SRCEXT0(:) = EXTINCT(I1)*SOURCE(:, KANG, I1)

                if (.not. doGridPtsPrp_) then
                    ! zgrid correction (to address layer averaging issue in z direction)
                    ! a more robust criteria should be developed in the future
                    !   to distinguish the atmopheric medium intersection
                    !╭────────────────────────────────────────────────╮!
                    if (abs(EXT0-EXT1) .gt. min(EXT0, EXT1)*10.0) then
                        if (bitZ .eq. 0) then
                            EXT0 = EXT1
                            SRCEXT0(:) = SRCEXT1(:)
                        else
                            EXT1 = EXT0
                            SRCEXT1(:) = SRCEXT0(:)
                        endif
                    end if
                    !╰────────────────────────────────────────────────╯!
                end if

                ! Compute the cell radiance: integration of the source function
                EXT = 0.5*(EXT0 + EXT1)
                tau = EXT*SO
                if (tau .ge. 0.5) then
                    transCell = exp(-tau)
                    absCell = 1.0 - transCell
                else
                    absCell = tau*(1.0 - 0.5*tau*&
                        &(1.0 - FRAC13_*tau*(1 - 0.25*tau)))
                    transCell = 1.0 - absCell
                end if
                if (tau .le. 2.0) then
                    if (EXT .eq. 0.0) then
                        SRC(:) = 0.0
                    else
                        SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + FRAC112_&
                            &*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                    end if
                else
                    ! Combined first order expansion and constant extinction formula
                    EXT0P = EXT0
                    SRCEXT0P(:) = SRCEXT0(:)
                    if (tau .gt. 4.0) then
                        EXT0P = EXT1 + (EXT0 - EXT1)*4.0/tau
                        if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                    end if
                    SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                        &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                        &*(1 - 2/tau + 2*transCell/absCell))
                end if

                SRC(1) = max(SRC(1), 0.0)

                ! See if there are valid radiances for the intersection face:
                ! If so, set flag, interpolate radiance, and add to
                !   radiance from source function in cell
                ! If not, only add in cell radiance and prepare for next cell

                !if (GRIDRAD(1, I1) .ge. -0.1) then
                if (GRIDRAD(1, I1) .ge. -accSolv_) then
                    isValidRad = .true.
                    RAD0(:) = GRIDRAD(:, I1)
                    RAD(:) = RAD(:) + trans*(RAD0(:)*transCell + SRC(:)*absCell)
                else
                    RAD(:) = RAD(:) + trans*SRC(:)*absCell
                    SRCEXT1(:) = SRCEXT0(:)
                    trans = trans*transCell
                    EXT1 = EXT0
                    iCell = iCellNext
                end if
            end do

            GRIDRAD(:, iPts) = RAD(:)
        end if
    end do

    return

end subroutine shd_intg__backIntegrate1D

end module shd_intg
