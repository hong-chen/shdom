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
    use globals_cfg,   only: verbose_, accSolv_, &
                             isMpi_, is1DX_, is1DYx_, is1D_, &
                             isPol_
    use globals_type,  only: shd_type_rtm, &
                             shd_type_grid

    use globals_obj, only: shdo_obj_, iradI_Shdo, radI_Shdo, isrcJ_Shdo, srcJ_Shdo, &
                            src_obj_, &
                            prp_obj_, &
                            sfc_obj_, brdf_Rtm

    use shd_sfc,  only: shd_sfc__interpSfc

    use shd_init, only: shd_init__makeTopRad, &
                        shd_init__makeBotRadUL, &
                        shd_init__makeBotRadVL, &
                        shd_init__makeBotRadU, &
                        shd_init__makeBotRadV

    use shd_grid, only: shd_grid__makeSweepOrder, &
                        shd_grid__findBndryPnt

    use shd_util, only: shd_util__findNextCell

    use shd_shdo, only: shd_shdo__transSHtoDO_, &
                        shd_shdo__transSHtoDO, &
                        shd_shdo__transDOtoSH_, &
                        shd_shdo__transDOtoSH

    use shd_mpi,  only: shd_mpi__findBndryPnt, &
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
    &SFC_BRDF_DO,&
    &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &CMU1, CPHI1, CMU2, CPHI2, WSAVE, FFTFLAG,&
    &DIRFLUX, FLUXES,&
    &WORK, SWEEPORD, GRIDRAD,&
    &rtm_obj, grid_obj)

    !   Performs the path integrations through the medium specified by
    ! the extinction (rtm_obj%ext) and source function (srcJ_Shdo) in
    ! spherical harmonic space.  The source function is transformed to
    ! discrete ordinates one zenith angle at a time.  The integrations
    ! are done from the top boundary (with no reflection) down to the
    ! bottom boundary for the downwelling angles, and then up from the
    ! bottom to the top after computing the surface reflection/emission.
    ! After the path integrations produce the discrete angle radiance
    ! it is transformed and added to output radiance spherical harmonic
    ! expansion.  The hemispheric flux FLUXES is also computed.
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

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS
    integer BCPTR(MAXNBC, 2)
    integer SWEEPORD(grid_obj%Npts, *)
    logical FFTFLAG(shdo_obj_%Nmu)
    real SFC_BRDF_DO(shdo_obj_%Nstokes, shdo_obj_%Nmu/2, shdo_obj_%maxNphi0, shdo_obj_%Nstokes, shdo_obj_%Nmu/2, shdo_obj_%maxNphi0)
    real CMU1(*), CPHI1(*), CMU2(*), CPHI2(*), WSAVE(*)
    real BCRAD(shdo_obj_%Nstokes, MAXBCRAD)
    real DIRFLUX(grid_obj%Npts), FLUXES(2, grid_obj%Npts)
    real GRIDRAD(shdo_obj_%Nstokes, grid_obj%Npts), WORK(shdo_obj_%Nstokes, shdo_obj_%maxNphi0, grid_obj%Npts)

    integer I, I1, K, NR, IPHI, IMU, IBC, IANG, IUPDOWN, OLDNPTS
    integer joct, ipt, ip, iorder, ipcell, icorner
    logical :: theflag
    real A, TRANSMIN
    save OLDNPTS

    integer :: count_start, count_rate, count_max
    integer :: count_end
    real :: elapsed_time

    ! Set the minimum transmission for the cell tracing (1 for single cell)
    TRANSMIN = 1.00

    ! Make the new grid cell/point sweeping order if need to
    if (grid_obj%Npts .ne. OLDNPTS) then
        call shd_grid__makeSweepOrder(SWEEPORD, GRIDRAD, rtm_obj, grid_obj)
    end if

    call shd_mpi__findBndryPnt(rtm_obj%flag_bc, rtm_obj%flag_ip, SWEEPORD,&
        &rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, rtm_obj, grid_obj)

    if (grid_obj%Npts .ne. OLDNPTS) then
        ! Make the pointers to the grid points on the top and bottom boundaries
        call shd_grid__findBndryPnt(shdo_obj_%Nang, MAXNBC, MAXBCRAD, &
            &rtm_obj%z(1), rtm_obj%z(rtm_obj%Nz),&
            &NTOPPTS, NBOTPTS, BCPTR, grid_obj)
        if (sfc_obj_%type(1:1) .eq. 'V') then
            ! If this is a variable surface then bilinearly interpolate
            !   the surface temperature and reflection parameters to the
            !   bottom grid points.
            call shd_sfc__interpSfc(NBOTPTS, BCPTR(1, 2), grid_obj%xyz,&
                &src_obj_%type, src_obj_%waveNum, src_obj_%avgWaveLen, src_obj_%unit, brdf_Rtm)
        end if
    end if
    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(src_obj_%type, src_obj_%skyRad, src_obj_%waveNum, src_obj_%avgWaveLen,&
        &src_obj_%unit, NTOPPTS, shdo_obj_%Nstokes, BCRAD(1, 1))

    ! Zero the flux arrays and set temporary radiance to -1 for not valid
    do I = 1, grid_obj%Npts
        FLUXES(1, I) = 0.0
        FLUXES(2, I) = 0.0
        GRIDRAD(1, I) = -1.0
    end do

    ! Zero the output spherical harmonic radiance
    NR = iradI_Shdo(grid_obj%Npts + 1)
    do I = 1, NR
        radI_Shdo(:, I) = 0.0
    end do

    ! print *, 'Path integration:'
    ! Loop over the zenith angles (downwelling angles must be first)
    IANG = 1
    DO_MU_LOOP: do IMU = 1, shdo_obj_%Nmu

        if (IMU .eq. shdo_obj_%Nmu/2 + 1) then
            ! For the first upwelling angle, make the bottom boundary radiances
            !   for the Lambertian surfaces.  Compute the upwelling bottom
            !   radiances using the downwelling fluxes.
            if (sfc_obj_%type .eq. 'FL') then
                call shd_init__makeBotRadUL(NBOTPTS, BCPTR(1, 2),&
                    &DIRFLUX, FLUXES, src_obj_%type,&
                    &src_obj_%waveNum, src_obj_%avgWaveLen, src_obj_%unit, shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
            else if (sfc_obj_%type .eq. 'VL') then
                call shd_init__makeBotRadVL(NBOTPTS, BCPTR(1, 2),&
                    &DIRFLUX, FLUXES, src_obj_%type,&
                    &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
            end if
        end if

        if (isPol_) then
            call shd_shdo__transSHtoDO(grid_obj%Npts, shdo_obj_%Nstokes, prp_obj_%NlegSet, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphi0, shdo_obj_%Nphi0(IMU),&
                &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
                &isrcJ_Shdo, srcJ_Shdo, WORK)
        else
            call shd_shdo__transSHtoDO_(grid_obj%Npts, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphi0, shdo_obj_%Nphi0(IMU),&
                &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
                &isrcJ_Shdo, srcJ_Shdo, WORK)
        end if

        ! Enforce the discrete ordinate I source function to be non-negative
        do I = 1, grid_obj%Npts
            do IPHI = 1, shdo_obj_%Nphi0(IMU)
                WORK(1, IPHI, I) = max(0.0, WORK(1, IPHI, I))
            end do
        end do

        ! Loop over the azimuthal angles
        DO_PHI_LOOP: do IPHI = 1, shdo_obj_%Nphi0(IMU)

            if (shdo_obj_%mu(IMU) .lt. 0.0) then

                ! For downward ordinates, initialize the top boundary radiances
                IUPDOWN = 1
                do IBC = 1, NTOPPTS
                    I = BCPTR(IBC, 1)
                    GRIDRAD(:, I) = BCRAD(:, IBC)
                end do

            else

                ! Initialize the bottom boundary radiances
                IUPDOWN = 2
                if (.not. sfc_obj_%isLambertian) then

                    ! If a Lambertian surface, use the radiances computed above
                    !   otherwise, compute the radiance for this angle
                    !   by integrating over the stored downwelling radiances.
                    ! If the surface reflection is uniform use the precomputed
                    !   BRDF matrix to do the integration.
                    if (sfc_obj_%isUniform) then
                        call shd_init__makeBotRadU(NBOTPTS, BCPTR(1, 2),&
                            &shdo_obj_%Nmu, shdo_obj_%maxNphi0, shdo_obj_%Nphi0, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt, IMU, IPHI,&
                            &src_obj_%type, src_obj_%avgWaveLen, src_obj_%avgMu, src_obj_%avgAz, DIRFLUX,&
                            &SFC_BRDF_DO,&
                            &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
                    else
                        call shd_init__makeBotRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                            &shdo_obj_%Nmu, shdo_obj_%maxNphi0, shdo_obj_%Nphi0, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt,&
                            &shdo_obj_%mu(IMU), shdo_obj_%phi(IMU, IPHI),&
                            &src_obj_%type, src_obj_%avgWaveLen, src_obj_%avgMu, src_obj_%avgAz, DIRFLUX,&
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

                call shd_intg__backIntegrate1D(rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, shdo_obj_%maxNphi0,&
                    &SWEEPORD, shdo_obj_%mu(IMU), shdo_obj_%phi(IMU, IPHI), TRANSMIN,&
                    &rtm_obj%ext, shdo_obj_%Nstokes, WORK, IPHI, GRIDRAD, grid_obj)

            else if (is1DYx_ .and. (.not. isMpi_)) then

                call shd_intg__backIntegrate2D(rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, shdo_obj_%maxNphi0,&
                    &SWEEPORD, shdo_obj_%mu(IMU), shdo_obj_%phi(IMU, IPHI), TRANSMIN,&
                    &rtm_obj%ext, shdo_obj_%Nstokes, WORK, IPHI, GRIDRAD, grid_obj)

            else
                if (isPol_) then

                    call shd_intg__backIntegrate3D(rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, shdo_obj_%maxNphi0,&
                        &SWEEPORD,&
                        &shdo_obj_%mu(IMU), shdo_obj_%phi(IMU, IPHI), TRANSMIN, rtm_obj%flag_bc, rtm_obj%flag_ip,&
                        &rtm_obj%ext, shdo_obj_%Nstokes, WORK, IPHI, GRIDRAD, rtm_obj, grid_obj)

                else

                    call shd_intg__backIntegrate3D_(rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, shdo_obj_%maxNphi0,&
                        &SWEEPORD,&
                        &shdo_obj_%mu(IMU), shdo_obj_%phi(IMU, IPHI), TRANSMIN, rtm_obj%flag_bc, rtm_obj%flag_ip,&
                        &rtm_obj%ext, WORK, IPHI, GRIDRAD, rtm_obj, grid_obj)

                end if
            end if

            ! Transfer the radiance for the (IMU,IPHI) discrete ordinate from
            !   GRIDRAD to the work array.
            do I = 1, grid_obj%Npts
                if (GRIDRAD(1, I) .lt. -accSolv_) then
                    if (verbose_ .ge. 4) write (6, *) 'Warning [_integratePath]: Grid point has no value',&
                        &I, shdo_obj_%mu(IMU), shdo_obj_%phi(IMU, IPHI), GRIDRAD(1, I), rtm_obj%x(1), rtm_obj%y(1),&
                        &grid_obj%xyz(1, I), grid_obj%xyz(2, I), grid_obj%xyz(3, I)
                    ! stop
                    GRIDRAD(1, I) = 0.0
                else
                    WORK(:, IPHI, I) = GRIDRAD(:, I)
                    GRIDRAD(1, I) = -1.0
                end if
            end do

            ! Do the sum for the hemispheric fluxes for this angle
            A = abs(shdo_obj_%mu(IMU))*shdo_obj_%wgt(IMU, IPHI)
            do I = 1, grid_obj%Npts
                FLUXES(IUPDOWN, I) = FLUXES(IUPDOWN, I) + A*WORK(1, IPHI, I)
            end do

            ! For a general BRDF surface save the downwelling radiance
            !   at the bottom boundary points.
            if ((.not. sfc_obj_%isLambertian) .and. (IMU .le. shdo_obj_%Nmu/2)) then
                I1 = NTOPPTS + NBOTPTS*IANG
                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    BCRAD(:, IBC + I1) = WORK(:, IPHI, I)
                end do
            end if

            IANG = IANG + 1
        end do DO_PHI_LOOP

        ! Transform the radiance to spherical harmonics for this zenith angle
        if (isPol_) then
            call shd_shdo__transDOtoSH(grid_obj%Npts, shdo_obj_%Nstokes, prp_obj_%NlegSet, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphi0, shdo_obj_%Nphi0(IMU),&
                &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
                &iradI_Shdo, WORK, radI_Shdo)
        else
            call shd_shdo__transDOtoSH_(grid_obj%Npts, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                &shdo_obj_%Nmu, shdo_obj_%maxNphi0, shdo_obj_%Nphi0(IMU),&
                &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
                &iradI_Shdo, WORK, radI_Shdo)
        end if

    end do DO_MU_LOOP

    OLDNPTS = grid_obj%Npts

    return

end subroutine shd_intg__integratePath

subroutine shd_intg__backIntegrate3D(&
    &NX, NY, NZ, NA,&
    &SWEEPORD,&
    &MU, PHI, TRANSMIN, BCFLAG, IPFLAG,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD,&
    &rtm_obj, grid_obj)

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

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer NX, NY, NZ, NA, NSTOKES, KANG
    integer BCFLAG, IPFLAG
    integer SWEEPORD(grid_obj%Npts, *)
    real MU, PHI, TRANSMIN
    real EXTINCT(grid_obj%Npts), SOURCE(NSTOKES, NA, grid_obj%Npts)
    real GRIDRAD(NSTOKES, grid_obj%Npts)
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

    EPS = 1.0e-3*(grid_obj%xyz(3, grid_obj%i8Grid(8, 1)) - grid_obj%xyz(3, grid_obj%i8Grid(1, 1)))

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
    if (is1DYx_) then
        JOCT = JOCTORDER2(IOCT)
    else
        JOCT = JOCTORDER3(IOCT)
    end if

    IZ = (NZ - 1)*(1 - BITZ) + 1
    ! Sweep through all the grid points
    do IORDER = 1, grid_obj%Npts
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = grid_obj%i8Grid(ICORNER, IPCELL)

        if (isMpi_) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if  ((grid_obj%xyz(3, IPT) .lt. rtm_obj%z(IZ) .and. BITZ .eq. 0) .or.&
                &(grid_obj%xyz(3, IPT) .gt. rtm_obj%z(IZ) .and. BITZ .eq. 1)) then
                IZ = IZ + 2*BITZ - 1
                call shd_mpi__calBndryRad(BCFLAG, IPFLAG, JOCT, IZ,&
                    &NX, NY, NZ,&
                    &NA,&
                    &MU, PHI, EXTINCT, NSTOKES, SOURCE,&
                    &KANG, GRIDRAD, rtm_obj, grid_obj)
            end if
        end if

        if (GRIDRAD(1, IPT) .lt. 0.0) then

            ! Start off at the grid point
            ICELL = IPCELL
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
            XE = grid_obj%xyz(1, IPT)
            YE = grid_obj%xyz(2, IPT)
            ZE = grid_obj%xyz(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)
                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate3D]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL, xe, ye, ze,&
                        &grid_obj%xyz(1, ipt), grid_obj%xyz(2, ipt), grid_obj%xyz(3, ipt)
                    stop
                end if
                IPINX = btest(int(grid_obj%flag(ICELL)), 0)
                IPINY = btest(int(grid_obj%flag(ICELL)), 1)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = grid_obj%i8Grid(9 - IOCT, ICELL)
                ! Get the distances to the 3 planes and select the closest
                if (IPINX) then
                    SOX = 1.0e20
                else
                    SOX = (grid_obj%xyz(1, IOPP) - XE)*CXINV
                end if
                if (IPINY) then
                    SOY = 1.0e20
                else
                    SOY = (grid_obj%xyz(2, IOPP) - YE)*CYINV
                end if
                SOZ = (grid_obj%xyz(3, IOPP) - ZE)*CZINV
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
                INEXTCELL = grid_obj%i6Cell(IFACE, ICELL)
                if (INEXTCELL .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL,&
                        &INEXTCELL, grid_obj)
                end if

                ! Get the face grid pointers
                if (grid_obj%i6Cell(IFACE, ICELL) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    KFACE = IFACE
                    IC = ICELL
                else
                    ! If going to smaller face then use next face (more accurate)
                    KFACE = OPPFACE(IFACE)
                    IC = INEXTCELL
                end if
                I1 = grid_obj%i8Grid(GRIDFACE(1, KFACE), IC)
                I2 = grid_obj%i8Grid(GRIDFACE(2, KFACE), IC)
                I3 = grid_obj%i8Grid(GRIDFACE(3, KFACE), IC)
                I4 = grid_obj%i8Grid(GRIDFACE(4, KFACE), IC)
                ! Compute the face interpolation factors
                if (JFACE .eq. 1) then
                    U = (ZE - grid_obj%xyz(3, I1))/(grid_obj%xyz(3, I3) - grid_obj%xyz(3, I1))
                    if (IPINY) then
                        V = 0.5
                    else
                        V = (YE - grid_obj%xyz(2, I1))/(grid_obj%xyz(2, I2) - grid_obj%xyz(2, I1))
                    end if
                else if (JFACE .eq. 2) then
                    U = (ZE - grid_obj%xyz(3, I1))/(grid_obj%xyz(3, I3) - grid_obj%xyz(3, I1))
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - grid_obj%xyz(1, I1))/(grid_obj%xyz(1, I2) - grid_obj%xyz(1, I1))
                    end if
                else
                    if (IPINY) then
                        U = 0.5
                    else
                        U = (YE - grid_obj%xyz(2, I1))/(grid_obj%xyz(2, I3) - grid_obj%xyz(2, I1))
                    end if
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - grid_obj%xyz(1, I1))/(grid_obj%xyz(1, I2) - grid_obj%xyz(1, I1))
                    end if
                end if
                if (u .lt. -0.0001 .or. v .lt. -0.0001 .or.&
                    &u .gt. 1.0001 .or. v .gt. 1.0001) then
                    print *, '_backIntegrate3D: u,v<0 or u,v>1: ',&
                        &mu, phi, xe, ye, ze, u, v, iface, jface, kface,&
                        &grid_obj%xyz(1, i1), grid_obj%xyz(2, i1), grid_obj%xyz(3, i1)
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    if (JFACE .eq. 1) then
                        XE = grid_obj%xyz(1, grid_obj%i8Grid(IOCT, INEXTCELL))
                    else if (JFACE .eq. 2) then
                        YE = grid_obj%xyz(2, grid_obj%i8Grid(IOCT, INEXTCELL))
                    else
                        ZE = grid_obj%xyz(3, grid_obj%i8Grid(IOCT, INEXTCELL))
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
                !if (RAD(1) < -1.0e-5) then
                if (RAD(1) < -accSolv_) then
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
                            &grid_obj%xyz(1, i1), grid_obj%xyz(2, i1), grid_obj%xyz(3, i1)
                        print '(i7,4(1x,f7.4))', i2, GRIDRAD(1, I2),&
                            &grid_obj%xyz(1, i2), grid_obj%xyz(2, i2), grid_obj%xyz(3, i2)
                        print '(i7,4(1x,f7.4))', i3, GRIDRAD(1, I3),&
                            &grid_obj%xyz(1, i3), grid_obj%xyz(2, i3), grid_obj%xyz(3, i3)
                        print '(i7,4(1x,f7.4))', i4, GRIDRAD(1, I4),&
                            &grid_obj%xyz(1, i4), grid_obj%xyz(2, i4), grid_obj%xyz(3, i4)
                        print '(i7,8x,3(1x,f7.4))', ipt,&
                            &grid_obj%xyz(1, ipt), grid_obj%xyz(2, ipt), grid_obj%xyz(3, ipt)
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
    &NX, NY, NZ, NA,&
    &SWEEPORD,&
    &MU, PHI, TRANSMIN, BCFLAG, IPFLAG,&
    &EXTINCT, SOURCE, KANG, GRIDRAD,&
    &rtm_obj, grid_obj)

    ! Description:
    !   Sweeps through the spatial grid computing radiances for the
    !   discrete ordinate direction (MU,PHI) by integrating the source
    !   function and extinction backward from each point. The sweeping
    !   order through the grid points is in SWEEPORD. If a point is
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
    !   grid_obj%Npts : integer, total number of grid points
    !   grid_obj%Ncell : integer, total number of cells
    !   grid_obj%i8Grid(8,grid_obj%Ncell) : integer array, pointer to 8 grid points belonging to a cell.
    !                       Order of pointers is based on binary digits Z-Y-X:
    !                       IOCT = 1 + BITX + 2*BITY + 4*BITZ, where
    !                       BITc is 0 for negative side and 1 for positive side.
    !                       1: (-X, -Y, -Z), 2: (+X, -Y, -Z), 3: (-X, +Y, -Z), 4: (+X, +Y, -Z)
    !                       5: (-X, -Y, +Z), 6: (+X, -Y, +Z), 7: (-X, +Y, +Z), 8: (+X, +Y, +Z)
    !   grid_obj%i6Cell(6,grid_obj%Ncell): integer array, pointer to 6 neighboring cells of a cell.
    !                       Order of pointers are 1:-X, 2:+X, 3:-Y, 4:+Y, 5:-Z, 6:+Z.
    !                       Is 0 if there is no neighbor cell (boundary).
    !                       Is positive if there is just one cell across face
    !                       of current cell (unique neighbor).  Is negative
    !                       (-index) if there is more than one cell across
    !                       the face of current cell.  In this case the
    !                       cell pointed to is the parent cell that has a
    !                       has same or larger face adjoining current cell.
    !   grid_obj%i2Tree(2,grid_obj%Ncell) : integer array, pointer to parent and children cells in tree structure.
    !                       First is pointer to parent cell, 0 for a base cell.
    !                       Second is pointer to first of the two children cells.
    !                       This pointer is 0 if there are no children cells.
    !                       The children cells are in successive order, with
    !                       the negative side first. See grid_obj%flag for direction
    !                       of division of split.
    !   grid_obj%flag(grid_obj%Ncell) : integer array, a collection of binary flags for the cells:
    !                       bit  For bit set
    !                       0:  Do independent pixel for X
    !                       1:  Do independent pixel for Y
    !                       2,3: Direction cell is split (0=no split, 1=X, 2=Y, 3=Z)
    !   grid_obj%xyz(3,grid_obj%Npts)    : real array, X,Y,Z coordinate values of grid point locations
    !   rtm_obj%x, rtm_obj%y, rtm_obj%z: real array, grid locations for each dimension.
    !                        To facilitate the periodic horizontal boundaries,
    !                        there are NX+1 values in rtm_obj%x, NY+1 in rtm_obj%y,
    !                        NZ values in rtm_obj%z. Currently the grid locations
    !                        are evenly spaced in X and Y.
    !   SWEEPORD(grid_obj%Npts,*)   : integer array, grid point sweeping order

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_grid), intent(in) :: grid_obj

    integer(kind=I__), intent(in) :: NX, NY, NZ, NA, KANG, BCFLAG, IPFLAG
    integer(kind=I__), intent(in) :: SWEEPORD(grid_obj%Npts, *)
    real(kind=R__)   , intent(in) :: MU, PHI, TRANSMIN
    real(kind=R__)   , intent(in) :: EXTINCT(grid_obj%Npts), SOURCE(NA, grid_obj%Npts)

    real   (kind=R__), intent(inout) :: GRIDRAD(grid_obj%Npts)

    integer(kind=I__) :: BITX, BITY, BITZ, IOCT, JOCT
    integer(kind=I__) :: IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer(kind=I__) :: I1, I2, I3, I4, IOPP, JFACE, KFACE, IC, IZ
    integer(kind=I__) :: IORDER, ICORNER
    logical :: VALIDRAD, VALIDFACE, IPINX, IPINY

    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS, path
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
    real(RD_) :: EXT0P, SRCEXT0P
    real(RD_) :: TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANSMIT

    integer(kind=I__) :: GRIDFACE(4, 6), OPPFACE(6), JOCTORDER_1D(8), JOCTORDER_3D(8)
    data GRIDFACE /&
        &1, 3, 5, 7,&
        &2, 4, 6, 8,&
        &1, 2, 5, 6,&
        &3, 4, 7, 8,&
        &1, 2, 3, 4,&
        &5, 6, 7, 8/
    data OPPFACE /2, 1, 4, 3, 6, 5/
    data JOCTORDER_1D /1, 3, 1, 3, 2, 4, 2, 4/
    data JOCTORDER_3D /1, 3, 5, 7, 2, 4, 6, 8/

    integer :: iwhile, ibreak

    EPS = 1.0e-3*(grid_obj%xyz(3, grid_obj%i8Grid(8, 1)) - grid_obj%xyz(3, grid_obj%i8Grid(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    CX = dsqrt(dble(1.0d0 - MU**2))*dcos(dble(PHI + DPI_))
    CY = dsqrt(dble(1.0d0 - MU**2))*dsin(dble(PHI + DPI_))
    CZ = dble(-MU)

    if (abs(CX) .gt. 1.0e-7) then
        CXINV = 1.0d0/CX
    else
        if (verbose_ .ge. 7) write (6, *) &
            'Message [_backIntegrate3D_]: CX too small, reset CX to 0, CXINV to 1e8.'
        CX = 0.0
        CXINV = 1.0e8
    end if

    if (abs(CY) .gt. 1.0e-7) then
        CYINV = 1.0d0/CY
    else
        if (verbose_ .ge. 7) write (6, *) &
            'Message [_backIntegrate3D_]: CY too small, reset CY to 0, CYINV to 1e8.'
        CY = 0.0
        CYINV = 1.0e8
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
    if (is1DYx_) then
        JOCT = JOCTORDER_1D(IOCT)
    else
        JOCT = JOCTORDER_3D(IOCT)
    end if

    IZ = (NZ - 1)*(1 - BITZ) + 1

    ! Sweep through all the grid points
    do IORDER = 1, grid_obj%Npts

        IPCELL  = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = grid_obj%i8Grid(ICORNER, IPCELL)

        ! this is for mpi (e.g., shdMpi.f90, not OMP)
        !╭────────────────────────────────────────────────╮!
        if (isMpi_) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if  ((grid_obj%xyz(3, IPT) .lt. rtm_obj%z(IZ) .and. BITZ .eq. 0) .or.&
                &(grid_obj%xyz(3, IPT) .gt. rtm_obj%z(IZ) .and. BITZ .eq. 1)) then
                IZ = IZ + 2*BITZ - 1
                call shd_mpi__calBndryRad(BCFLAG, IPFLAG, JOCT, IZ,&
                    &NX, NY, NZ,&
                    &NA,&
                    &MU, PHI, EXTINCT, 1, SOURCE,&
                    &KANG, GRIDRAD, rtm_obj, grid_obj)
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
            XE = grid_obj%xyz(1, IPT)
            YE = grid_obj%xyz(2, IPT)
            ZE = grid_obj%xyz(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)

                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate3D_]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL, XE, YE, ZE,&
                        &grid_obj%xyz(1, IPT), grid_obj%xyz(2, IPT), grid_obj%xyz(3, IPT)
                    stop
                end if

                IPINX = btest(int(grid_obj%flag(ICELL)), 0)
                IPINY = btest(int(grid_obj%flag(ICELL)), 1)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = grid_obj%i8Grid(9 - IOCT, ICELL)

                ! Get the distances to the 3 planes and select the closest
                !╭────────────────────────────────────────────────╮!
                if (IPINX) then
                    SOX = 1.0e20
                else
                    SOX = (grid_obj%xyz(1, IOPP) - XE)*CXINV
                end if

                if (IPINY) then
                    SOY = 1.0e20
                else
                    SOY = (grid_obj%xyz(2, IOPP) - YE)*CYINV
                end if

                SOZ = (grid_obj%xyz(3, IOPP) - ZE)*CZINV

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
                INEXTCELL = grid_obj%i6Cell(IFACE, ICELL)
                if (INEXTCELL .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL,&
                        &INEXTCELL, grid_obj)
                end if

                ! Get the face grid pointers
                if (grid_obj%i6Cell(IFACE, ICELL) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    KFACE = IFACE
                    IC = ICELL
                else
                    ! If going to smaller face then use next face (more accurate)
                    KFACE = OPPFACE(IFACE)
                    IC = INEXTCELL
                end if

                I1 = grid_obj%i8Grid(GRIDFACE(1, KFACE), IC)
                I2 = grid_obj%i8Grid(GRIDFACE(2, KFACE), IC)
                I3 = grid_obj%i8Grid(GRIDFACE(3, KFACE), IC)
                I4 = grid_obj%i8Grid(GRIDFACE(4, KFACE), IC)

                ! Compute the face interpolation factors
                if (JFACE .eq. 1) then
                    U = (ZE - grid_obj%xyz(3, I1))/(grid_obj%xyz(3, I3) - grid_obj%xyz(3, I1))
                    if (IPINY) then
                        V = 0.5
                    else
                        V = (YE - grid_obj%xyz(2, I1))/(grid_obj%xyz(2, I2) - grid_obj%xyz(2, I1))
                    end if
                else if (JFACE .eq. 2) then
                    U = (ZE - grid_obj%xyz(3, I1))/(grid_obj%xyz(3, I3) - grid_obj%xyz(3, I1))
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - grid_obj%xyz(1, I1))/(grid_obj%xyz(1, I2) - grid_obj%xyz(1, I1))
                    end if
                else
                    if (IPINY) then
                        U = 0.5
                    else
                        U = (YE - grid_obj%xyz(2, I1))/(grid_obj%xyz(2, I3) - grid_obj%xyz(2, I1))
                    end if
                    if (IPINX) then
                        V = 0.5
                    else
                        V = (XE - grid_obj%xyz(1, I1))/(grid_obj%xyz(1, I2) - grid_obj%xyz(1, I1))
                    end if
                end if

                if (U .lt. -0.0001 .or. V .lt. -0.0001 .or.&
                    &U .gt. 1.0001 .or. V .gt. 1.0001) then
                    print *, '_backIntegrate3D_: U,V<0 or U,V>1: ',&
                        &MU, PHI, XE, YE, ZE, U, V, IFACE, JFACE, KFACE,&
                        &grid_obj%xyz(1, I1), grid_obj%xyz(2, I1), grid_obj%xyz(3, I1)
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    if (JFACE .eq. 1) then
                        XE = grid_obj%xyz(1, grid_obj%i8Grid(IOCT, INEXTCELL))
                    else if (JFACE .eq. 2) then
                        YE = grid_obj%xyz(2, grid_obj%i8Grid(IOCT, INEXTCELL))
                    else
                        ZE = grid_obj%xyz(3, grid_obj%i8Grid(IOCT, INEXTCELL))
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
                !if (RAD < -1.0e-5) then
                if (RAD < -accSolv_) then
                    print '(A,3(1X,F6.4),8(1X,E12.5))',&
                        &'_backIntegrate3D_: RAD<0 ',&
                        &XE, YE, ZE, EXT0, EXT1, SO, TAU, TRANSMIT, SRC, ABSCELL, RAD
                    call shd_mpi__abort('The end')
                end if

                VALIDFACE = (GRIDRAD(I1) .ge. -0.1 .and. GRIDRAD(I2) .ge. -0.1&
                      &.and. GRIDRAD(I3) .ge. -0.1 .and. GRIDRAD(I4) .ge. -0.1)
                if ( (INEXTCELL .le. 0) .or.&
                    &((TRANSMIT .le. TRANSMIN) .and. VALIDFACE)) then
                    if (VALIDFACE) then
                        VALIDRAD = .true.
                        RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2)&
                            &+ F3*GRIDRAD(I3) + F4*GRIDRAD(I4)
                        RAD = RAD + TRANSMIT*RAD0
                    else
                        print *, 'Error [_backIntegrate3D_]: INEXTCELL=0: ', xe, ye, ze, icell
                        print *, iz, joct, mu, phi
                        print '(i7,4(1x,f7.4))', i1, GRIDRAD(I1),&
                            &grid_obj%xyz(1, i1), grid_obj%xyz(2, i1), grid_obj%xyz(3, i1)
                        print '(i7,4(1x,f7.4))', i2, GRIDRAD(I2),&
                            &grid_obj%xyz(1, i2), grid_obj%xyz(2, i2), grid_obj%xyz(3, i2)
                        print '(i7,4(1x,f7.4))', i3, GRIDRAD(I3),&
                            &grid_obj%xyz(1, i3), grid_obj%xyz(2, i3), grid_obj%xyz(3, i3)
                        print '(i7,4(1x,f7.4))', i4, GRIDRAD(I4),&
                            &grid_obj%xyz(1, i4), grid_obj%xyz(2, i4), grid_obj%xyz(3, i4)
                        print '(i7,8x,3(1x,f7.4))', ipt,&
                            &grid_obj%xyz(1, ipt), grid_obj%xyz(2, ipt), grid_obj%xyz(3, ipt)
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
    &NX, NY, NZ, NA,&
    &SWEEPORD, MU, PHI, TRANSMIN,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD,&
    &grid_obj)

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

    type(shd_type_grid), intent(in) :: grid_obj
    integer NX, NY, NZ, NA, KANG, NSTOKES
    integer SWEEPORD(grid_obj%Npts, *)
    real MU, PHI, TRANSMIN
    real EXTINCT(grid_obj%Npts), SOURCE(NSTOKES, NA, grid_obj%Npts)
    real GRIDRAD(NSTOKES, grid_obj%Npts)
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

    EPS = 1.0e-3*(grid_obj%xyz(3, grid_obj%i8Grid(8, 1)) - grid_obj%xyz(3, grid_obj%i8Grid(1, 1)))

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
    do IORDER = 1, grid_obj%Npts
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = grid_obj%i8Grid(ICORNER, IPCELL)
        if (GRIDRAD(1, IPT) .lt. 0.0) then

            ! Start off at the grid point
            ICELL = IPCELL
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
            XE = grid_obj%xyz(1, IPT)
            YE = grid_obj%xyz(2, IPT)
            ZE = grid_obj%xyz(3, IPT)

            ! Loop until finding a face with known radiances
            VALIDRAD = .false.
            do while (.not. VALIDRAD)
                ! Make sure current cell is valid
                if (ICELL .le. 0) then
                    write (6, *) 'Error [_backIntegrate2D]: ICELL=0 ',&
                        &MU, PHI, IPT, INEXTCELL
                    stop
                end if
                IPINX = btest(int(grid_obj%flag(ICELL)), 0)

                ! Find boundaries of the current cell
                ! Find the three possible intersection planes (X,Y,Z)
                !   from the coordinates of the opposite corner grid point
                IOPP = grid_obj%i8Grid(9 - IOCT, ICELL)
                ! Get the distances to the 3 planes and select the closest
                if (IPINX) then
                    SOX = 1.0e20
                else
                    SOX = (grid_obj%xyz(1, IOPP) - XE)*CXINV
                end if
                SOZ = (grid_obj%xyz(3, IOPP) - ZE)*CZINV
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
                INEXTCELL = grid_obj%i6Cell(IFACE, ICELL)
                if (INEXTCELL .lt. 0) then
                    call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL,&
                        &INEXTCELL, grid_obj)
                end if

                ! Get the face grid pointers
                if (grid_obj%i6Cell(IFACE, ICELL) .ge. 0) then
                    ! If going to same or larger face then use previous face
                    KFACE = IFACE
                    IC = ICELL
                else
                    ! If going to smaller face then use next face (more accurate)
                    KFACE = OPPFACE(IFACE)
                    IC = INEXTCELL
                end if
                I1 = grid_obj%i8Grid(GRIDFACE(1, KFACE), IC)
                I2 = grid_obj%i8Grid(GRIDFACE(2, KFACE), IC)
                ! Compute the face interpolation factors
                if (JFACE .eq. 1) then
                    U = (ZE - grid_obj%xyz(3, I1))/(grid_obj%xyz(3, I2) - grid_obj%xyz(3, I1))
                else
                    if (IPINX) then
                        U = 0.5
                    else
                        U = (XE - grid_obj%xyz(1, I1))/(grid_obj%xyz(1, I2) - grid_obj%xyz(1, I1))
                    end if
                end if

                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    if (JFACE .eq. 1) then
                        XE = grid_obj%xyz(1, grid_obj%i8Grid(IOCT, INEXTCELL))
                    else
                        ZE = grid_obj%xyz(3, grid_obj%i8Grid(IOCT, INEXTCELL))
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
                            &grid_obj%xyz(1, i1), grid_obj%xyz(2, i1), grid_obj%xyz(3, i1)
                        print '(i7,4(1x,f7.4))', i2, GRIDRAD(1, I2),&
                            &grid_obj%xyz(1, i2), grid_obj%xyz(2, i2), grid_obj%xyz(3, i2)
                        print '(i7,8x,3(1x,f7.4))', ipt,&
                            &grid_obj%xyz(1, ipt), grid_obj%xyz(2, ipt), grid_obj%xyz(3, ipt)
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
    &NX, NY, NZ, NA,&
    &SWEEPORD, MU, PHI, TRANSMIN,&
    &EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD,&
    &grid_obj)

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

    type(shd_type_grid), intent(in) :: grid_obj

    integer(kind=I__), intent(in) :: NX, NY, NZ, NA, KANG, NSTOKES
    integer(kind=I__), intent(in) :: SWEEPORD(grid_obj%Npts, *)
    real(kind=R__)   , intent(in) :: MU, PHI, TRANSMIN
    real(kind=R__)   , intent(in) :: EXTINCT(grid_obj%Npts), SOURCE(NSTOKES, NA, grid_obj%Npts)

    real(kind=R__), intent(inout) :: GRIDRAD(NSTOKES, grid_obj%Npts)

    real(kind=RD_) :: CZ, CZINV
    real(kind=RD_) :: XE, YE, ZE
    real(kind=RD_) :: SO, EPS
    real(kind=RD_) :: EXT, EXT0, EXT1
    real(kind=RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(kind=RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(kind=RD_) :: TAU, TRANSCELL, ABSCELL, TRANSMIT
    real(kind=RD_) :: RAD(NSTOKES), RAD0(NSTOKES)

    integer(kind=I__) :: BITZ, IOCT, JOCT
    integer(kind=I__) :: IPCELL, ICELL, INEXTCELL, IPT, IFACE, I1, IOPP
    integer(kind=I__) :: IORDER, ICORNER
    logical :: VALIDRAD

    integer(kind=I__) :: GRIDFACE(6), JOCTORDER(8)
    data GRIDFACE /0, 0, 0, 0, 1, 5/
    data JOCTORDER /1, 1, 1, 1, 2, 2, 2, 2/

    EPS = 1.0e-3*(grid_obj%xyz(3, grid_obj%i8Grid(8, 1)) - grid_obj%xyz(3, grid_obj%i8Grid(1, 1)))

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
    do IORDER = 1, grid_obj%Npts
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = grid_obj%i8Grid(ICORNER, IPCELL)
        if (GRIDRAD(1, IPT) .lt. 0.0) then
            ! Start off at the grid point
            ICELL = IPCELL
            TRANSMIT = 1.0
            EXT1 = EXTINCT(IPT)
            RAD(:) = 0.0
            SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
            XE = grid_obj%xyz(1, IPT)
            YE = grid_obj%xyz(2, IPT)
            ZE = grid_obj%xyz(3, IPT)

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
                IOPP = grid_obj%i8Grid(9 - IOCT, ICELL)
                ! Get the distances to the 3 planes and select the closest
                SO = (grid_obj%xyz(3, IOPP) - ZE)*CZINV
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
                INEXTCELL = grid_obj%i6Cell(IFACE, ICELL)

                ! Get the face grid pointers
                I1 = grid_obj%i8Grid(GRIDFACE(IFACE), ICELL)
                ! Get the location coordinate (does the boundary wrapping)
                if (INEXTCELL .gt. 0) then
                    ZE = grid_obj%xyz(3, grid_obj%i8Grid(IOCT, INEXTCELL))
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

                SRC(1) = max(SRC(1), 0.0)

                ! See if there are valid radiances for the intersection face:
                ! If so, set flag, interpolate radiance, and add to
                !   radiance from source function in cell
                ! If not, only add in cell radiance and prepare for next cell

                !if (GRIDRAD(1, I1) .ge. -0.1) then
                if (GRIDRAD(1, I1) .ge. -accSolv_) then
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
