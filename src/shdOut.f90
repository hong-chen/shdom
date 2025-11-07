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

module shd_out

    use globals_const, only: RD_, PI_
    use globals_cfg,   only: verbose_, charMaxLen_, is1DX_, &
                             isPol_, isPol4_

    use globals_type,  only: shd_type_rtm

    use globals_obj,   only: rtm_obj_, xyz_Rtm_, &
                                      temp_Rtm_, &
                                       ext_Rtm_, &
                                       ssa_Rtm_, &
                                       leg_Rtm_, &
                                      ipha_Rtm_, &
                            cell_obj_, i8Grid_Cell_, &
                                       i6Cell_Cell_, &
                                       i2Tree_Cell_, &
                                         flag_Cell_, &
                            shdo_obj_, iradI_Shdo_, &
                                        radI_Shdo_, &
                                       isrcJ_Shdo_, &
                                        srcJ_Shdo_, &
                             src_obj_, &
                             prp_obj_, &
                             sfc_obj_, &
                             ckd_obj_, &
                             sen_obj_, view_Sen_

    use shd_grid, only: shd_grid__testCellSplit

    use shd_util, only: shd_util__findCell, &
                        shd_util__interpTriLin

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_io__outputResult
    public :: shd_io__outputResultMpi
    public :: shd_io__outputCamRad
    public :: shd_io__outputImg

    public :: shd_io__outputCellSplit
    public :: shd_io__visualizeCell

    public :: shd_io__saveState

    private :: shd_io__outputDataBin

contains

subroutine shd_io__outputResult(&
    &INSAVEFILE, OUTSAVEFILE,&
    &DELTAM,&
    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    &FLUXES, DIRFLUX, FLUXDIV,&
    &IRAD, RADOUT, NSHOUT, SHTERMS, SOURCE1OUT,&
    &OUTTYPE, OUTPARMS, OUTFILE)

    !   Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are five types (OUTTYPE) of output: 'R' - for radiance,
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms, and 'M' for
    ! medium properties.
    ! See the overall program documentation for output parameters.

    integer, intent(in) :: MAXITER, TOTITER, NSHOUT
    logical  , intent(in) :: DELTAM
    real, intent(in) :: SOLACC, SPLITACC, SHACC
    real, intent(inout) :: OUTPARMS(*)
    character(len=1), intent(in) :: OUTTYPE
    character(len=*), intent(in) :: INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real, intent(in) :: FLUXES(2, rtm_obj_%Npts), DIRFLUX(rtm_obj_%Npts), FLUXDIV(rtm_obj_%Npts)
    real, intent(in) :: SOURCE1OUT(shdo_obj_%Nstokes, *)

    real, intent(in) :: RADOUT(shdo_obj_%Nstokes, *), SHTERMS(NSHOUT, rtm_obj_%Npts)
    integer, intent(inout) :: IRAD

    integer IZ, I, K, N, NX2, NY2, J, JX, JY, icell
    integer NANGOUT, NXOUT, NYOUT
    logical ALLGRID
    real(RD_) :: MUOUT, PHID
    real :: C, SUM
    real XDOMAIN, YDOMAIN, STARTX, STARTY
    real, allocatable :: ASYM(:)
    real, allocatable :: FUP(:), FDN(:), FDIR(:)
    real, allocatable :: SUM1(:), SUM2(:), SUM3(:)
    real, allocatable :: SHTERM5(:)
    real(RD_) :: X0, Y0, Z0
    character(len=128) :: FORM
    character(len=128) :: GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

    character(len=charMaxLen_) :: OUTFILE_DATA
    character(len=20)  :: DATA_POSTFIX
    integer :: NVAR, NDATA, IDATA

    DATA_POSTFIX = '.sHdOmNG-out'
    OUTFILE_DATA = trim(adjustl(OUTFILE))//trim(adjustl(DATA_POSTFIX))

    GRIDNAME = 'EVEN (X,Y)  '
    if (rtm_obj_%type .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (rtm_obj_%type .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (src_obj_%type .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (src_obj_%type .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (src_obj_%type .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (sfc_obj_%type .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (sfc_obj_%type .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (sfc_obj_%type .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (sfc_obj_%type .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (sfc_obj_%type .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (sfc_obj_%type .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (sfc_obj_%type .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (sfc_obj_%type(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (src_obj_%unit .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (src_obj_%unit .eq. 'B') then
        if (OUTTYPE .eq. 'F') then
            UNITSNAME = 'WATTS/(M^2)'
        else if (OUTTYPE .eq. 'H') then
            UNITSNAME = 'WATTS/(M^2 KM)'
        else
            UNITSNAME = 'WATTS/(M^2 STER)'
        end if
    else
        if (OUTTYPE .eq. 'F' .or. OUTTYPE .eq. 'S') then
            UNITSNAME = 'WATTS/(M^2 MICRON)'
        else if (OUTTYPE .eq. 'H') then
            UNITSNAME = 'WATTS/(M^2 MICRON KM)'
        else
            UNITSNAME = 'WATTS/(M^2 MICRON STER)'
        end if
    end if
    if (OUTTYPE .eq. 'R') OUTNAME = 'RADIANCE'
    if (OUTTYPE .eq. 'F') OUTNAME = 'FLUX'
    if (OUTTYPE .eq. 'H') OUTNAME = 'NET_FLUX_DIV'
    if (OUTTYPE .eq. 'S') OUTNAME = 'SPHERICAL-HARMONIC'
    if (OUTTYPE .eq. 'J') OUTNAME = 'SOURCE_FUNCTION'
    if (OUTTYPE .eq. 'M') OUTNAME = 'MEDIUM-PROPERTIES'

    open (UNIT=2, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN')
    write (2, '(A,A)')&
        &'! Polarized Spherical Harmonic Discrete Ordinate',&
        &' Radiative Transfer Output'
    write (2, '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)')&
        &'!  L=', shdo_obj_%maxNl, '  M=', shdo_obj_%maxNm, '  NLM=', shdo_obj_%Nlm,&
        &'   NMU=', shdo_obj_%Nmu, '  NPHI=', shdo_obj_%Nphi, '  NANG=', shdo_obj_%Nang,&
        &'   NSH=', shdo_obj_%Nsh
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', shdo_obj_%Nstokes,&
        &'   NX=', rtm_obj_%Nx, '   NY=', rtm_obj_%Ny, '   NZ=', rtm_obj_%Nz,&
        &'    NPTS=', rtm_obj_%Npts, '   NCELLS=', cell_obj_%Npts
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(prp_obj_%fname)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',&
        &trim(ckd_obj_%fname), '   NUM_G=', ckd_obj_%Ng
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)

    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,&
            &'      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if

    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,&
        &'   INDEPENDENT_PIXEL=',  rtm_obj_%flag_ip
    write (2, '(A,A22,A,I1)') '!  SURFACE_TYPE=', SFCNAME,&
        &'   HORIZ_BOUNDARY_COND=', rtm_obj_%flag_bc

    if (src_obj_%doSolar) then

        if (sfc_obj_%type(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(sfc_obj_%fname),&
                &'  SKY_RAD=', src_obj_%skyRad
        else
            if (src_obj_%type .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', sfc_obj_%avgAlb, '  GROUND_TEMP=', sfc_obj_%avgTemp,&
                    &'  SKY_RAD=', src_obj_%skyRad
            else
                write (2, '(A,F9.7,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', sfc_obj_%avgAlb, '  SKY_RAD=', src_obj_%skyRad
            end if
        end if

        write (2, '(A,E13.6,A,F10.7,A,F8.3)')&
            &'!  SOLAR_FLUX=', src_obj_%avgSolar, '   SOLAR_MU=', src_obj_%avgMu,&
            &'   SOLAR_AZ=', src_obj_%avgAz*180.0/acos(-1.0)

    else

        if (sfc_obj_%type(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(sfc_obj_%fname)
        else
            write (2, '(A,F8.3,A,F9.7)')&
                &'!  GROUND_TEMP=', sfc_obj_%avgTemp, '  GROUND_EMIS=', 1.0 - sfc_obj_%avgAlb
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', src_obj_%skyRad

    end if

    if (src_obj_%unit .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (src_obj_%unit .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,&
            &'   WAVENUMBER_RANGE=', src_obj_%waveNum(1), src_obj_%waveNum(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,&
            &'   WAVELENGTH=', src_obj_%avgWaveLen
    end if

    write (2, '(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,&
        &'   SPHERICAL_HARMONIC_ACCURACY=', SHACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER,&
        &'   NUMBER_ITERATIONS=', TOTITER
    write (2, '(A,A)') '!  OUTPUT_TYPE=', trim(adjustl(OUTNAME))

    if (OUTTYPE .eq. 'R') then

        ! Radiance output
        ! NXOUT
        !╭────────────────────────────────────────────────╮!
        XDOMAIN = rtm_obj_%x(rtm_obj_%Nx + 1 - ibits(rtm_obj_%flag_bc, 0, 1)) - rtm_obj_%x(1) + 2*OUTPARMS(4)
        STARTX = -OUTPARMS(4)
        if (OUTPARMS(2) .eq. 0.0) then
            NXOUT = 1
        else
            NXOUT = max(1, nint(XDOMAIN/OUTPARMS(2)))
        end if
        !╰────────────────────────────────────────────────╯!

        ! NYOUT
        !╭────────────────────────────────────────────────╮!
        YDOMAIN = rtm_obj_%y(rtm_obj_%Ny + 1 - ibits(rtm_obj_%flag_bc, 1, 1)) - rtm_obj_%y(1) + 2*OUTPARMS(5)
        STARTY = -OUTPARMS(5)
        if (OUTPARMS(3) .eq. 0.0) then
            NYOUT = 1
        else
            NYOUT = max(1, nint(YDOMAIN/OUTPARMS(3)))
        end if
        !╰────────────────────────────────────────────────╯!

        ! Z0
        !╭────────────────────────────────────────────────╮!
        Z0 = min(max(OUTPARMS(1), rtm_obj_%z(1)), rtm_obj_%z(rtm_obj_%Nz))
        !╰────────────────────────────────────────────────╯!

        ! NANGOUT
        !╭────────────────────────────────────────────────╮!
        NANGOUT = nint(OUTPARMS(6))
        !╰────────────────────────────────────────────────╯!

        write (2, '(A,F7.3,3(A,I10))') '!    RADIANCE AT Z=', Z0,&
            &'    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (src_obj_%unit .eq. 'T') then
            FORM = '(2(1X,F7.3),1X,4(1X,F6.2))'
            write (2, '(A)')&
                &'!   X       Y       RADIANCE (STOKES I, Q, U, V)'
        else
            FORM = '(2(1X,F7.3),1X,4(1X,E12.5))'
            write (2, '(A)')&
                &'!   X       Y       RADIANCE (STOKES I, Q, U, V)'
        end if

        ! data processing
        !╭────────────────────────────────────────────────╮!
        do K = 1, NANGOUT
            if (sen_obj_%fname(1:4) .ne. 'NONE') then
                MUOUT = view_Sen_(1, K)
                PHID  = view_Sen_(2, K)
            else
                MUOUT = OUTPARMS(2*K + 5)
                PHID  = OUTPARMS(2*K + 6)
            end if

            if ((K .lt. 5) .or. (K .gt. (NANGOUT-5))) then
                write (2, '(A,1X,I8,F8.5,1X,F7.2,2X,A)')&
                    &'! ', K, MUOUT, PHID, '<- (mu,phi)'
            else if (K .eq. 11) then
                write (2, '(A)') '! ...'
            end if

        end do
        !╰────────────────────────────────────────────────╯!

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following provides information for interpreting binary data:'
        write (2, '(A)') '! R'
        write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
        write (2, '(A, 5(I10,A))') '! ', NXOUT, ',', NYOUT, ',', 1, ',', NANGOUT, ',', shdo_obj_%Nstokes, ''
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <R> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        N = NANGOUT*NXOUT*NYOUT
        call shd_io__outputDataBin(RADOUT, shdo_obj_%Nstokes, N, IRAD, OUTFILE_DATA)
        !╰────────────────────────────────────────────────╯!

    else if (OUTTYPE .eq. 'F') then

        if (src_obj_%unit .eq. 'T') then
            C = 1.0 / PI_
        else
            C = 1.0
        end if

        if (src_obj_%doSolar) then
            NVAR = 3
        else
            NVAR = 2
        end if

        ! Hemispheric flux output, F1
        ! Format 1 is flux at top and bottom of medium at grid point
        ! NVAR=2 or 3, NSET=1, NX=rtm_obj_%Nx, NY=rtm_obj_%Ny, NZ=1
        if (nint(OUTPARMS(1)) .le. 1) then

            write (2, '(A,F7.3)') '!    UPWELLING FLUX:     Z=', rtm_obj_%z(rtm_obj_%Nz)
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', rtm_obj_%z(1)

            if (src_obj_%unit .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
                    write (2, '(A)')&
                        & '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)') '!   X       Y           UP           DOWN'
                end if
            end if

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! F1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtm_obj_%NxA, ',', rtm_obj_%NyA, ',', 1, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = rtm_obj_%NxA*rtm_obj_%NyA
            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, rtm_obj_%Nz:rtm_obj_%Npts0:rtm_obj_%Nz), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1, 1:rtm_obj_%Npts0:rtm_obj_%Nz), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(DIRFLUX(1:rtm_obj_%Npts0:rtm_obj_%Nz), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do I = 1, rtm_obj_%Npts0, rtm_obj_%Nz
            !    if  (xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
            !        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) then
            !        if (src_obj_%type .ne. 'T') then
            !            write (2, FORM) xyz_Rtm_(1, I), xyz_Rtm_(2, I),&
            !                & C*FLUXES(2, I + rtm_obj_%Nz - 1), C*FLUXES(1, I), DIRFLUX(I)
            !        else
            !            write (2, FORM) xyz_Rtm_(1, I), xyz_Rtm_(2, I),&
            !                & C*FLUXES(2, I + rtm_obj_%Nz - 1), C*FLUXES(1, I)
            !        end if
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 2 is flux at a given level at regular locations
        ! NVAR=2 or 3, NSET=1, NX=NXOUT, NY=NYOUT, NZ=1
        else if (nint(OUTPARMS(1)) .eq. 2) then

            Z0 = min(max(OUTPARMS(2), rtm_obj_%z(1)), rtm_obj_%z(rtm_obj_%Nz))
            if (OUTPARMS(3) .eq. 0.0) OUTPARMS(3) = 1.0
            NX2 = rtm_obj_%Nx + 1 - ibits(rtm_obj_%flag_bc, 0, 1)
            NXOUT = max(1, nint((rtm_obj_%x(NX2) - rtm_obj_%x(1))/OUTPARMS(3)))
            if (OUTPARMS(4) .eq. 0.0) OUTPARMS(4) = 1.0
            NY2 = rtm_obj_%Ny + 1 - ibits(rtm_obj_%flag_bc, 1, 1)
            NYOUT = max(1, nint((rtm_obj_%y(NY2) - rtm_obj_%y(1))/OUTPARMS(4)))
            write (2, '(A,F7.3,2(A,I4))') '!    UPWELLING FLUX:     Z=', Z0,&
                & '    NXO=', NXOUT, '    NYO=', NYOUT
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', Z0
            if (src_obj_%unit .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
                    write (2, '(A)')&
                        & '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)') '!   X       Y           UP           DOWN'
                end if
            end if

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(FUP(NXOUT*NYOUT))
            allocate(FDN(NXOUT*NYOUT))
            allocate(FDIR(NXOUT*NYOUT))

            IDATA = 0
            Y0 = rtm_obj_%y(1)
            do JY = 1, NYOUT
                X0 = rtm_obj_%x(1)
                do JX = 1, NXOUT
                    IDATA = IDATA + 1
                    call shd_util__findCell(X0, Y0, Z0, icell)
                    call shd_util__interpTriLin(icell,&
                        & X0, Y0, Z0, 2, FLUXES, 1, 1, FDN(IDATA))
                    call shd_util__interpTriLin(icell,&
                        & X0, Y0, Z0, 2, FLUXES, 2, 1, FUP(IDATA))
                    call shd_util__interpTriLin(icell,&
                        & X0, Y0, Z0, 1, DIRFLUX, 1, 1, FDIR(IDATA))
                    X0 = X0 + OUTPARMS(3)
                end do
                Y0 = Y0 + OUTPARMS(4)
            end do
            !╰────────────────────────────────────────────────╯!

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! F2'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', NXOUT, ',', NYOUT, ',', 1, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NXOUT*NYOUT
            IDATA = 0
            call shd_io__outputDataBin(C*FUP, 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FDN, 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(FDIR, 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(FUP)
            deallocate(FDN)
            deallocate(FDIR)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !Y0 = rtm_obj_%y(1)
            !do JY = 1, NYOUT
            !    X0 = rtm_obj_%x(1)
            !    do JX = 1, NXOUT
            !        call shd_util__findCell(rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz,&
            !            & rtm_obj_%flag_bc, rtm_obj_%flag_ip, X0, Y0, Z0, icell, cell_obj_)
            !        call shd_util__interpTriLin(icell,&
            !            & X0, Y0, Z0, 2, FLUXES, 1, 1, FDN, cell_obj_)
            !        call shd_util__interpTriLin(icell,&
            !            & X0, Y0, Z0, 2, FLUXES, 2, 1, FUP, cell_obj_)
            !        call shd_util__interpTriLin(icell,&
            !            & X0, Y0, Z0, 1, DIRFLUX, 1, 1, FDIR, cell_obj_)
            !        if (src_obj_%type .ne. 'T') then
            !            write (2, FORM) X0, Y0, C*FUP, C*FDN, FDIR
            !        else
            !            write (2, FORM) X0, Y0, C*FUP, C*FDN
            !        end if
            !        X0 = X0 + OUTPARMS(3)
            !    end do
            !    Y0 = Y0 + OUTPARMS(4)
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 3 is domain averaged vertical profile
        ! NVAR=2 or 3, NSET=1, NX=rtm_obj_%NxA, NY=rtm_obj_%NyA, NZ=1
        else if (nint(OUTPARMS(1)) .eq. 3) then

            if (src_obj_%unit .eq. 'T') then
                FORM = '(1X,F7.3,1X,2(2X,F6.2))'
                write (2, '(A)') '!   Z        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(1X,F7.3,1X,3(2X,E12.5))'
                    write (2, '(A)')&
                        & '!    Z          UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(1X,F7.3,1X,2(2X,E12.5))'
                    write (2, '(A)') '!    Z          UP           DOWN'
                end if
            end if

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(SUM1(rtm_obj_%Nz))
            allocate(SUM2(rtm_obj_%Nz))
            allocate(SUM3(rtm_obj_%Nz))

            do IZ = 1, rtm_obj_%Nz
                SUM1(IZ) = 0.0
                SUM2(IZ) = 0.0
                SUM3(IZ) = 0.0
                N = 0
                do I = IZ, rtm_obj_%Npts0, rtm_obj_%Nz
                    if  (xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
                        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) then
                        SUM1(IZ) = SUM1(IZ) + FLUXES(1, I)
                        SUM2(IZ) = SUM2(IZ) + FLUXES(2, I)
                        SUM3(IZ) = SUM3(IZ) + DIRFLUX(I)
                        N = N + 1
                    end if
                end do
            end do
            !╰────────────────────────────────────────────────╯!

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! F3'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F3> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            IDATA = 0
            call shd_io__outputDataBin(C*SUM2/N, 1, rtm_obj_%Nz, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*SUM1/N, 1, rtm_obj_%Nz, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(SUM3/N, 1, rtm_obj_%Nz, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(SUM1)
            deallocate(SUM2)
            deallocate(SUM3)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do IZ = 1, rtm_obj_%Nz
            !    SUM1 = 0.0
            !    SUM2 = 0.0
            !    SUM3 = 0.0
            !    N = 0
            !    do I = IZ, rtm_obj_%Npts0, rtm_obj_%Nz
            !        if  (xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
            !            &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) then
            !            SUM1 = SUM1 + FLUXES(1, I)
            !            SUM2 = SUM2 + FLUXES(2, I)
            !            SUM3 = SUM3 + DIRFLUX(I)
            !            N = N + 1
            !        end if
            !    end do
            !    if (src_obj_%type .ne. 'T') then
            !        write (2, FORM) rtm_obj_%z(IZ), C*SUM2/N, C*SUM1/N, SUM3/N
            !    else
            !        write (2, FORM) rtm_obj_%z(IZ), C*SUM2/N, C*SUM1/N
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 4 is fluxes at every base grid point
        ! Format 5 is fluxes at every grid point
        ! NVAR=2 or 3, NSET=1, NX=rtm_obj_%NxA, NY=rtm_obj_%NyA, NZ=rtm_obj_%Nz
        else if (nint(OUTPARMS(1)) .ge. 4) then

            if (nint(OUTPARMS(1)) .eq. 5) then
                N = rtm_obj_%Npts
                ALLGRID = .true.
            else
                N = rtm_obj_%Npts0
                ALLGRID = .false.
            end if

            if (src_obj_%unit .eq. 'T') then
                FORM = '(3(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y       Z        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(3(1X,F7.3),1X,3(2X,E12.5))'
                    write (2, '(A,A)') '!   X       Y       Z',&
                        & '           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(3(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)')&
                        & '!   X       Y       Z           UP           DOWN'
                end if
            end if

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            if (ALLGRID) then
                write (2, '(A)') '! F5'
                write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
                ! added x, y, z (NVAR+3)
                write (2, '(A, 2(I10,A))') '! ', N, ',', NVAR+3, ''
            else
                write (2, '(A)') '! F4'
                write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
                write (2, '(A, 5(I10,A))') '! ', rtm_obj_%NxA, ',', rtm_obj_%NyA, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
            end if
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (ALLGRID) then
                if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F5> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            else
                if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F4> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            end if

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, :), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(DIRFLUX(:), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if

            ! add grid positions (or binary data cannot be reconstructed)
            if (ALLGRID) then
                IDATA = 0
                call shd_io__outputDataBin(xyz_Rtm_(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do I = 1, N
            !    if ((xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
            !        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) .or. ALLGRID) then
            !        if (src_obj_%type .ne. 'T') then
            !            write (2, FORM) xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I),&
            !                & C*FLUXES(2, I), C*FLUXES(1, I), DIRFLUX(I)
            !        else
            !            write (2, FORM) xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I),&
            !                & C*FLUXES(2, I), C*FLUXES(1, I)
            !        end if
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        end if

    ! Heating output: net flux convergence
    else if (OUTTYPE .eq. 'H') then

        NVAR = 1

        ! Format 1 is domain averaged vertical profile
        if (nint(OUTPARMS(1)) .eq. 1) then

            write (2, '(A,A)') '!    Z       -DIV(Fnet)'

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(SUM1(rtm_obj_%Nz))

            do IZ = 1, rtm_obj_%Nz
                SUM1(IZ) = 0.0
                N = 0
                do I = IZ, rtm_obj_%Npts0, rtm_obj_%Nz
                    if  (xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
                        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) then
                        SUM1(IZ) = SUM1(IZ) + FLUXDIV(I)
                        N = N + 1
                    end if
                end do
            end do
            !╰────────────────────────────────────────────────╯!

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! H1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <H1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            IDATA = 0
            call shd_io__outputDataBin(-SUM1/N, 1, rtm_obj_%Nz, IDATA, OUTFILE_DATA)

            deallocate(SUM1)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do IZ = 1, rtm_obj_%Nz
            !    SUM = 0.0
            !    N = 0
            !    do I = IZ, rtm_obj_%Npts0, rtm_obj_%Nz
            !        if  (xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
            !            &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) then
            !            SUM = SUM + FLUXDIV(I)
            !            N = N + 1
            !        end if
            !    end do
            !    write (2, '(1X,F7.3,2X,E13.5)') rtm_obj_%z(IZ), -SUM/N
            !end do
            !╰────────────────────────────────────────────────╯!

        else

            ! Format 2 is flux convergence for every base grid point
            ! Format 3 is flux convergence for every grid point
            if (nint(OUTPARMS(1)) .eq. 3) then
                N = rtm_obj_%Npts
                ALLGRID = .true.
            else
                N = rtm_obj_%Npts0
                ALLGRID = .false.
            end if

            write (2, '(A)') '!    X       Y        Z       -DIV(Fnet)'

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            if (ALLGRID) then
                write (2, '(A)') '! H3'
                write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
                ! added x, y, z (NVAR+3)
                write (2, '(A, 2(I10,A))') '! ', N, ',', NVAR+3, ''
            else
                write (2, '(A)') '! H2'
                write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
                write (2, '(A, 5(I10,A))') '! ', rtm_obj_%NxA, ',', rtm_obj_%NyA, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
            end if
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (ALLGRID) then
                if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <H3> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            else
                if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <H2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            end if

            IDATA = 0
            call shd_io__outputDataBin(-FLUXDIV, 1, N, IDATA, OUTFILE_DATA)

            ! add grid positions (or binary data cannot be reconstructed)
            if (ALLGRID) then
                IDATA = 0
                call shd_io__outputDataBin(xyz_Rtm_(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do I = 1, N
            !    if ((xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
            !        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) .or. ALLGRID) then
            !        write (2, '(3(1X,F7.3),2X,E13.5)')&
            !            & xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I), -FLUXDIV(I)
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        end if

    ! Spherical Harmonic output:
    !   Output mean intensity and net flux (x,y,z components).
    !   and maybe normalized rms of higher order terms
    else if (OUTTYPE .eq. 'S') then

        if (nint(OUTPARMS(1)) .eq. 2) then
            N = rtm_obj_%Npts
            ALLGRID = .true.
        else
            N = rtm_obj_%Npts0
            ALLGRID = .false.
        end if

        ! data processing
        !╭────────────────────────────────────────────────╮!
        allocate(SHTERM5(N))

        if (NSHOUT .eq. 5) then
            NVAR = 5
            write (2, '(A,A)') '!    X       Y        Z         Imean',&
                &'          Fx           Fy           Fz       HOrms'
            do I = 1, N
                SHTERM5(I) = SHTERMS(5, I)/max(1e-20, SHTERMS(1, I))
            end do
        else
            NVAR = 4
            write (2, '(A,A)') '!    X       Y        Z         Imean',&
                &'          Fx           Fy           Fz '
        end if
        !╰────────────────────────────────────────────────╯!

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following contains information for binary data intepretation:'
        if (ALLGRID) then
            write (2, '(A)') '! S2'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            ! added x, y, z (NVAR+3)
            write (2, '(A, 2(I10,A))') '! ', N, ',', NVAR+3, ''
        else
            write (2, '(A)') '! S1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtm_obj_%NxA, ',', rtm_obj_%NyA, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
        end if
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (ALLGRID) then
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <S2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        else
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <S1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        end if

        IDATA = 0
        call shd_io__outputDataBin(SHTERMS(1:4, :), 4, N, IDATA, OUTFILE_DATA)

        if (NSHOUT .eq. 5) then
            IDATA = 0
            call shd_io__outputDataBin(SHTERM5, 1, N, IDATA, OUTFILE_DATA, 'ADD')
        end if

        ! add grid positions (or binary data cannot be reconstructed)
        if (ALLGRID) then
            IDATA = 0
            call shd_io__outputDataBin(xyz_Rtm_(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
        end if

        deallocate(SHTERM5)
        !╰────────────────────────────────────────────────╯!

        ! original
        !╭────────────────────────────────────────────────╮!
        !if (NSHOUT .eq. 5) then
        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz       HOrms'
        !    do I = 1, N
        !        if ((xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
        !            &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) .or. ALLGRID) then
        !            write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
        !                &xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I),&
        !                &SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I),&
        !                &SHTERMS(5, I)/max(1e-20, SHTERMS(1, I))
        !        end if
        !    end do
        !else
        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz '
        !    do I = 1, N
        !        if ((xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
        !            &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) .or. ALLGRID) then
        !            write (2, '(3(1X,F7.3),1X,4(1X,E12.5))')&
        !                &xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I),&
        !                &SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I)
        !        end if
        !    end do
        !end if
        !╰────────────────────────────────────────────────╯!

    else if (OUTTYPE .eq. 'J') then

        NVAR = 5

        ! J source function output: source function for a specified ang
        !   at every grid point
        if (nint(OUTPARMS(1)) .eq. 2) then
            N = rtm_obj_%Npts
            ALLGRID = .true.
        else
            N = rtm_obj_%Npts0
            ALLGRID = .false.
        end if

        MUOUT = OUTPARMS(2)
        PHID = OUTPARMS(3)
        write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')&
            &'! ', MUOUT, PHID, '<- (mu,phi)'
        write (2, '(A,A)') '!    X       Y       Z    ',&
            &' Extinction   Source Function (Stokes I, Q, U, V)'

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following contains information for binary data intepretation:'
        if (ALLGRID) then
            write (2, '(A)') '! J2'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            ! added x, y, z (NVAR+3)
            write (2, '(A, 2(I10,A))') '! ', N, ',', NVAR+3, ''
        else
            write (2, '(A)') '! J1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtm_obj_%NxA, ',', rtm_obj_%NyA, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
        end if
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (ALLGRID) then
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <J2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        else
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <J1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        end if

        IDATA = 0
        call shd_io__outputDataBin(ext_Rtm_, 1, N, IDATA, OUTFILE_DATA)
        IDATA = 0
        call shd_io__outputDataBin(SOURCE1OUT, 4, N, IDATA, OUTFILE_DATA, 'ADD')

        ! add grid positions (or binary data cannot be reconstructed)
        if (ALLGRID) then
            IDATA = 0
            call shd_io__outputDataBin(xyz_Rtm_(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
        end if
        !╰────────────────────────────────────────────────╯!

        ! original
        !╭────────────────────────────────────────────────╮!
        !do I = 1, N
        !    if ((xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
        !        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) .or. ALLGRID) then
        !        write (2, '(3(1X,F7.3),5(1X,E12.5))')&
        !            & xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I),&
        !            & ext_Rtm_(I), SOURCE1OUT(:, I)
        !    end if
        !end do
        !╰────────────────────────────────────────────────╯!

    ! Medium output: grid point medium properties -
    !   Output extinction, single scattering albedo,
    !   asymmetry parameter, and temperature.
    else if (OUTTYPE .eq. 'M') then

        NVAR = 4

        if (nint(OUTPARMS(1)) .eq. 2) then
            N = rtm_obj_%Npts
            ALLGRID = .true.
        else
            N = rtm_obj_%Npts0
            ALLGRID = .false.
        end if

        write (2, '(A,A)') '!    X       Y       Z',&
            &'      Extinct  Albedo  Asymmetry Temp'

        ! data processing
        !╭────────────────────────────────────────────────╮!
        allocate(ASYM(N))
        do I = 1, N
            if (prp_obj_%Npha .gt. 0) then
                ASYM(I) = leg_Rtm_(1, 1, ipha_Rtm_(I))
            else
                ASYM(I) = leg_Rtm_(1, 1, I)
            end if
        end do
        !╰────────────────────────────────────────────────╯!

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following contains information for binary data intepretation:'
        if (ALLGRID) then
            write (2, '(A)') '! M2'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            ! added x, y, z (NVAR+3)
            write (2, '(A, 2(I10,A))') '! ', N, ',', NVAR+3, ''
        else
            write (2, '(A)') '! M1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtm_obj_%NxA, ',', rtm_obj_%NyA, ',', rtm_obj_%Nz, ',', 1, ',', NVAR, ''
        end if
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (ALLGRID) then
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <M2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        else
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <M1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        end if

        IDATA = 0
        call shd_io__outputDataBin(ext_Rtm_, 1, N, IDATA, OUTFILE_DATA)
        IDATA = 0
        call shd_io__outputDataBin(ssa_Rtm_, 1, N, IDATA, OUTFILE_DATA, 'ADD')
        IDATA = 0
        call shd_io__outputDataBin(ASYM, 1, N, IDATA, OUTFILE_DATA, 'ADD')
        IDATA = 0
        call shd_io__outputDataBin(temp_Rtm_, 1, N, IDATA, OUTFILE_DATA, 'ADD')

        ! add grid positions (or binary data cannot be reconstructed)
        if (ALLGRID) then
            IDATA = 0
            call shd_io__outputDataBin(xyz_Rtm_(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
        end if

        deallocate(ASYM)
        !╰────────────────────────────────────────────────╯!


        ! original
        !╭────────────────────────────────────────────────╮!
        !do I = 1, N
        !    if (prp_obj_%Npha .gt. 0) then
        !        ASYM = leg_Rtm_(1, 1, ipha_Rtm_(I))
        !    else
        !        ASYM = leg_Rtm_(1, 1, I)
        !    end if
        !    if ((xyz_Rtm_(1, I) .lt. rtm_obj_%x(rtm_obj_%Nx + 1) .and.&
        !        &xyz_Rtm_(2, I) .lt. rtm_obj_%y(rtm_obj_%Ny + 1)) .or. ALLGRID) then
        !        write (2, '(3(1X,F7.3),1X,F10.6,1X,F8.6,1X,F8.5,1X,F7.2)')&
        !            & xyz_Rtm_(1, I), xyz_Rtm_(2, I), xyz_Rtm_(3, I),&
        !            & ext_Rtm_(I), ssa_Rtm_(I), ASYM, temp_Rtm_(I)
        !    end if
        !end do
        !╰────────────────────────────────────────────────╯!

    end if

    close (2)
    return
end subroutine shd_io__outputResult

subroutine shd_io__radShiftCorr(&
    &DATA, NVAR, NSET, NX, NY, DX, DY,&
    &XDOMAIN, YDOMAIN,&
    &DS_Z, MU_SET, PHI_SET)

    integer, intent(in) :: NVAR, NSET, NX, NY
    real, intent(inout) :: DATA(NVAR, *)
    real, intent(in) :: DX, DY, XDOMAIN, YDOMAIN
    real, intent(in) :: MU_SET(*), PHI_SET(*)

    real(RD_) :: DS_Z
    real(RD_) :: DS_VZA_RAD, DS_XOFFSET, DS_YOFFSET
    real ::MUOUT, PHIOUT
    real(RD_) :: STARTX, STARTY
    integer :: XSHIFT, YSHIFT, ISET, IVAR
    real, allocatable :: DATA_COPY(:, :, :, :)

    allocate(DATA_COPY(NX, NY, NSET, NVAR))

    DATA_COPY = reshape(DATA(:, 1:NX*NY*NSET), (/NX, NY, NSET, NVAR/))

    do ISET = 1, NSET

        MUOUT = MU_SET(ISET)
        PHIOUT = PHI_SET(ISET)

        DS_VZA_RAD = acos(MUOUT)
        DS_XOFFSET = DS_Z*dtan(DS_VZA_RAD)*cos(PHIOUT)
        DS_YOFFSET = DS_Z*dtan(DS_VZA_RAD)*sin(PHIOUT)

        STARTX = DS_XOFFSET - DX*nint(DS_XOFFSET/DX)
        XSHIFT = nint(DS_XOFFSET/DX)
        if (STARTX .lt. 0.0d0) then
            XSHIFT = XSHIFT - 1
        end if

        STARTY = DS_YOFFSET - DY*nint(DS_YOFFSET/DY)
        YSHIFT = nint(DS_YOFFSET/DY)
        if (STARTY .lt. 0.0d0) then
            YSHIFT = YSHIFT - 1
        end if

        if (XSHIFT .ne. 0) then
            do IVAR = 1, NVAR
                DATA_COPY(:, :, ISET, IVAR) = cshift(DATA_COPY(:, :, ISET, IVAR), SHIFT=XSHIFT, DIM=1)
            end do
        end if

        if (YSHIFT .ne. 0) then
            do IVAR = 1, NVAR
                DATA_COPY(:, :, ISET, IVAR) = cshift(DATA_COPY(:, :, ISET, IVAR), SHIFT=YSHIFT, DIM=2)
            end do
        end if

    end do

    DATA(:, 1:NX*NY*NSET) = reshape(DATA_COPY, (/NVAR, NX*NY*NSET/))

    deallocate(DATA_COPY)

    return

end subroutine shd_io__radShiftCorr

subroutine shd_io__outputResultMpi(&
    &Npts_tot, Ncell_tot, Nsh_tot,&
    &INSAVEFILE, OUTSAVEFILE,&
    &DELTAM,&
    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    &FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT,&
    &OUTTYPE, OUTPARMS, OUTFILE,&
    &rtmAll_obj)

    ! Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are four supported types (OUTTYPE) of output: 'R' - for radi
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms
    ! See the overall program documentation for output parameters.

    type(shd_type_rtm),  intent(in) :: rtmAll_obj

    integer Nsh_tot, Npts_tot, Ncell_tot
    integer MAXITER, TOTITER, IRAD, NSHOUT
    logical DELTAM
    real SOLACC, SPLITACC, SHACC
    real OUTPARMS(*)
    character(len=1), intent(in) :: OUTTYPE
    character(len=*), intent(in) :: INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real FLUXES(3, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx), FLUXDIV(rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx)
    real RADOUT(shdo_obj_%Nstokes, *), SHTERMS(NSHOUT, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx)

    integer I, J, K, N, JX, JY, IZ
    integer NANGOUT, NXOUT, NYOUT
    real :: MUOUT, PHID
    real :: C, SUM
    real(RD_) :: STARTX, STARTY
    real(RD_) :: X0, Y0, Z0
    character(len=128) :: FORM
    character(len=64) :: GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

    real, allocatable :: MU_SET(:), PHI_SET(:)
    character(len=charMaxLen_) :: OUTFILE_DATA
    character(len=20)  :: DATA_POSTFIX
    integer :: NVAR, NDATA, IDATA
    real, allocatable :: SUM1(:), SUM2(:), SUM3(:)
    real, allocatable :: SHTERM5(:, :, :)

    DATA_POSTFIX = '.mpi.sHdOmNG-out'
    OUTFILE_DATA = trim(adjustl(OUTFILE))//trim(adjustl(DATA_POSTFIX))

    GRIDNAME = 'EVEN (X,Y)  '
    if (rtmAll_obj%type .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (rtmAll_obj%type .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (src_obj_%type .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (src_obj_%type .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (src_obj_%type .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (sfc_obj_%type .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (sfc_obj_%type .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (sfc_obj_%type .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (sfc_obj_%type .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (sfc_obj_%type .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (sfc_obj_%type .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (sfc_obj_%type .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (sfc_obj_%type(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (src_obj_%unit .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (src_obj_%unit .eq. 'B') then
        if (OUTTYPE .eq. 'F') then
            UNITSNAME = 'WATTS/(M^2)'
        else if (OUTTYPE .eq. 'H') then
            UNITSNAME = 'WATTS/(M^2 KM)'
        else
            UNITSNAME = 'WATTS/(M^2 STER)'
        end if
    else
        if (OUTTYPE .eq. 'F' .or. OUTTYPE .eq. 'S') then
            UNITSNAME = 'WATTS/(M^2 MICRON)'
        else if (OUTTYPE .eq. 'H') then
            UNITSNAME = 'WATTS/(M^2 MICRON KM)'
        else
            UNITSNAME = 'WATTS/(M^2 MICRON STER)'
        end if
    end if
    if (OUTTYPE .eq. 'R') OUTNAME = 'RADIANCE'
    if (OUTTYPE .eq. 'F') OUTNAME = 'FLUX'
    if (OUTTYPE .eq. 'H') OUTNAME = 'NET_FLUX_DIV'
    if (OUTTYPE .eq. 'S') OUTNAME = 'SPHERICAL-HARMONIC'

    open (UNIT=2, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN')
    write (2, '(A,A)')&
        & '! Polarized Spherical Harmonic Discrete Ordinate',&
        & ' Radiative Transfer Output'
    write (2, '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)')&
        & '!  L=', shdo_obj_%maxNl, '  M=', shdo_obj_%maxNm, '  NLM=', shdo_obj_%Nlm,&
        & '   NMU=', shdo_obj_%Nmu, '  NPHI=', shdo_obj_%Nphi, '  NANG=', shdo_obj_%Nang,&
        & '   NSH=', Nsh_tot
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', shdo_obj_%Nstokes,&
        & '   NX=', rtmAll_obj%Nx, '   NY=', rtmAll_obj%Ny, '   NZ=', rtmAll_obj%Nz,&
        & '   NPTS=', Npts_tot, '   NCELLS=', Ncell_tot
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(prp_obj_%fname)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',&
        & trim(ckd_obj_%fname), '   NUM_G=', ckd_obj_%Ng
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)
    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,&
            & '      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if
    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,&
        & '   INDEPENDENT_PIXEL=', rtmAll_obj%flag_ip
    write (2, '(A,A22,A,I2)') '!  SURFACE_TYPE=', SFCNAME,&
        & '   HORIZ_BOUNDARY_COND=', rtmAll_obj%flag_bc
    if (src_obj_%doSolar) then
        if (sfc_obj_%type(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(sfc_obj_%fname),&
                & '  SKY_RAD=', src_obj_%skyRad
        else
            if (src_obj_%type .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')&
                    & '!  GROUND_ALBEDO=', sfc_obj_%avgAlb, '  GROUND_TEMP=', sfc_obj_%avgTemp,&
                    & '  SKY_RAD=', src_obj_%skyRad
            else
                write (2, '(A,F9.7,A,E12.5)')&
                    & '!  GROUND_ALBEDO=', sfc_obj_%avgAlb, '  SKY_RAD=', src_obj_%skyRad
            end if
        end if
        write (2, '(A,E13.6,A,F10.7,A,F8.3)')&
            & '!  SOLAR_FLUX=', src_obj_%avgSolar, '   SOLAR_MU=', src_obj_%avgMu,&
            & '   SOLAR_AZ=', src_obj_%avgAz*180.0/acos(-1.0)
    else
        if (sfc_obj_%type(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(sfc_obj_%fname)
        else
            write (2, '(A,F8.3,A,F9.7)')&
                & '!  GROUND_TEMP=', sfc_obj_%avgTemp, '  GROUND_EMIS=', 1.0 - sfc_obj_%avgAlb
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', src_obj_%skyRad
    end if

    if (src_obj_%unit .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (src_obj_%unit .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,&
            & '   WAVENUMBER_RANGE=', src_obj_%waveNum(1), src_obj_%waveNum(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,&
            & '   WAVELENGTH=', src_obj_%avgWaveLen
    end if
    write (2, '(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,&
        & '   SPHERICAL_HARMONIC_ACCURACY=', SHACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER,&
        & '   NUMBER_ITERATIONS=', TOTITER
    write (2, '(A,A)') '!  OUTPUT_TYPE=', trim(adjustl(OUTNAME))

    if (OUTTYPE .eq. 'R') then
        ! Radiance output
        STARTX = 0.0
        if (OUTPARMS(2) .le. 0.0) then
            NXOUT = 1
        else
            NXOUT = max(1, nint(prp_obj_%xDomain/OUTPARMS(2)))
        end if
        STARTY = 0.0
        if (OUTPARMS(3) .le. 0.0) then
            NYOUT = 1
        else
            NYOUT = max(1, nint(prp_obj_%yDomain/OUTPARMS(3)))
        end if
        Z0 = min(max(OUTPARMS(1), rtmAll_obj%z(1)), rtmAll_obj%z(rtmAll_obj%Nz))
        NANGOUT = nint(OUTPARMS(6))

        write (2, '(A,F7.3,3(A,I10))') '!    RADIANCE AT Z=', Z0,&
            & '    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (src_obj_%unit .eq. 'T') then
            FORM = '(2(1X,F7.3),1X,4(1X,F6.2))'
            write (2, '(A)')&
                & '!   X       Y      RADIANCE (STOKES I, Q, U, V)'
        else
            FORM = '(2(1X,F7.3),1X,4(1X,E12.5))'
            write (2, '(A)')&
                & '!   X       Y      RADIANCE (STOKES I, Q, U, V)'
        end if

        ! data processing
        !╭────────────────────────────────────────────────╮!
        allocate(MU_SET(NANGOUT))
        allocate(PHI_SET(NANGOUT))

        do K = 1, NANGOUT
            if (sen_obj_%fname(1:4) .ne. 'NONE') then
                MUOUT = view_Sen_(1, K)
                PHID  = view_Sen_(2, K)
            else
                MUOUT = OUTPARMS(2*K + 5)
                PHID  = OUTPARMS(2*K + 6)
            endif

            if ((K .lt. 5) .or. (K .gt. (NANGOUT-5))) then
                write (2, '(A,1X,I8,F8.5,1X,F7.2,2X,A)')&
                    &'! ', K, MUOUT, PHID, '<- (mu,phi)'
            else if (K .eq. 11) then
                write (2, '(A)') '! ...'
            end if

            MU_SET(K) = MUOUT
            PHI_SET(K) = PHID*acos(-1.0)/180.0
        end do

        if (rtmAll_obj%flag_ip .eq. 0) then
            call shd_io__radShiftCorr(RADOUT(:, 1:NXOUT*NYOUT*NANGOUT), shdo_obj_%Nstokes, NANGOUT, NXOUT, NYOUT,&
                OUTPARMS(2), OUTPARMS(3), prp_obj_%xDomain, prp_obj_%yDomain, Z0-rtmAll_obj%z(1), MU_SET, PHI_SET)
        end if

        deallocate(MU_SET, PHI_SET)
        !╰────────────────────────────────────────────────╯!

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following contains information for binary data intepretation:'
        write (2, '(A)') '! R'
        write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
        write (2, '(A, 5(I10,A))') '! ', NXOUT, ',', NYOUT, ',', 1, ',', NANGOUT, ',', shdo_obj_%Nstokes, ''
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <R> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        N = NANGOUT*NXOUT*NYOUT
        call shd_io__outputDataBin(RADOUT, shdo_obj_%Nstokes, N, IRAD, OUTFILE_DATA)
        !╰────────────────────────────────────────────────╯!

        ! original
        !╭────────────────────────────────────────────────╮!
        !do K = 1, NANGOUT
        !    MUOUT = OUTPARMS(2*K + 5)
        !    PHID = OUTPARMS(2*K + 6)
        !    write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')&
        !        & '! ', MUOUT, PHID, '<- (mu,phi)'
        !    Y0 = STARTY
        !    do JY = 1, NYOUT
        !        X0 = STARTX
        !        do JX = 1, NXOUT
        !            IRAD = IRAD + 1
        !            write (2, FORM) X0, Y0, RADOUT(:, IRAD)
        !            X0 = X0 + OUTPARMS(2)
        !        end do
        !        Y0 = Y0 + OUTPARMS(3)
        !    end do
        !end do
        !╰────────────────────────────────────────────────╯!

    ! Hemispheric flux output
    ! Format 1 is flux at top and bottom of medium at grid point
    else if (OUTTYPE .eq. 'F') then

        if (src_obj_%unit .eq. 'T') then
            C = 1.0 / PI_
        else
            C = 1.0
        end if

        if (src_obj_%doSolar) then
            NVAR = 3
        else
            NVAR = 2
        end if

        if (nint(OUTPARMS(1)) .le. 1) then

            write (2, '(A,F7.3)') '!    UPWELLING FLUX:     Z=', rtmAll_obj%z(rtmAll_obj%Nz)
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', rtmAll_obj%z(1)
            if (src_obj_%unit .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
                    write (2, '(A)')&
                        & '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)') '!   X       Y           UP           DOWN'
                end if
            end if

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! F1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtmAll_obj%Nx, ',', rtmAll_obj%Ny, ',', 1, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <F1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = rtmAll_obj%Nx*rtmAll_obj%Ny
            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, rtmAll_obj%Nz, :, :), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1,  1, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(FLUXES(3,  1, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !C = 1.0
            !if (src_obj_%unit .eq. 'T') C = 1.0 / PI_
            !do I = 1, rtmAll_obj%Nx
            !    do J = 1, rtmAll_obj%Ny
            !        if (src_obj_%type .ne. 'T') then
            !            write (2, FORM) rtmAll_obj%x(I), rtmAll_obj%y(J),&
            !                & C*FLUXES(2, rtmAll_obj%Nz, J, I), C*FLUXES(1, 1, J, I), FLUXES(3, 1, J, I)
            !        else
            !            write (2, FORM) rtmAll_obj%x(I), rtmAll_obj%y(J),&
            !                & C*FLUXES(2, rtmAll_obj%Nz, J, I), C*FLUXES(1, 1, J, I)
            !        end if
            !    end do
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 2 is flux at a given level at regular locations
        ! Not supported for multiple processors
        else if (nint(OUTPARMS(1)) .eq. 2) then

            write(*, '(A)') 'Error [_outputResultMpi]: MPI does NOT support <F2>.'

        ! Format 3 is domain averaged vertical profile
        else if (nint(OUTPARMS(1)) .eq. 3) then

            if (src_obj_%unit .eq. 'T') then
                FORM = '(1X,F7.3,1X,2(2X,F6.2))'
                write (2, '(A)') '!   Z        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(1X,F7.3,1X,3(2X,E12.5))'
                    write (2, '(A)')&
                        & '!    Z          UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(1X,F7.3,1X,2(2X,E12.5))'
                    write (2, '(A)') '!    Z          UP           DOWN'
                end if
            end if

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(SUM1(rtmAll_obj%Nz))
            allocate(SUM2(rtmAll_obj%Nz))
            allocate(SUM3(rtmAll_obj%Nz))

            do IZ = 1, rtmAll_obj%Nz
                SUM1(IZ) = 0.0
                SUM2(IZ) = 0.0
                SUM3(IZ) = 0.0

                do I = 1, rtmAll_obj%Nx
                    do J = 1, rtmAll_obj%Ny
                        SUM1(IZ) = SUM1(IZ) + FLUXES(1, IZ, J, I)
                        SUM2(IZ) = SUM2(IZ) + FLUXES(2, IZ, J, I)
                        SUM3(IZ) = SUM3(IZ) + FLUXES(3, IZ, J, I)
                    end do
                end do

            end do
            !╰────────────────────────────────────────────────╯!

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! F3'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', rtmAll_obj%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <F3> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = rtmAll_obj%Nx*rtmAll_obj%Ny
            IDATA = 0
            call shd_io__outputDataBin(C*SUM2/N, 1, rtmAll_obj%Nz, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*SUM1/N, 1, rtmAll_obj%Nz, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(SUM3/N, 1, rtmAll_obj%Nz, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(SUM1)
            deallocate(SUM2)
            deallocate(SUM3)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !C = 1.0
            !if (src_obj_%unit .eq. 'T') C = 1.0 / PI_
            !N = rtmAll_obj%Nx*rtmAll_obj%Ny
            !do K = 1, rtmAll_obj%Nz
            !    SUM1 = 0.0
            !    SUM2 = 0.0
            !    SUM3 = 0.0
            !    do I = 1, rtmAll_obj%Nx
            !        do J = 1, rtmAll_obj%Ny
            !            SUM1 = SUM1 + FLUXES(1, K, J, I)
            !            SUM2 = SUM2 + FLUXES(2, K, J, I)
            !            SUM3 = SUM3 + FLUXES(3, K, J, I)
            !        end do
            !    end do
            !    if (src_obj_%type .ne. 'T') then
            !        write (2, FORM) rtmAll_obj%z(K), C*SUM2/N, C*SUM1/N, SUM3/N
            !    else
            !        write (2, FORM) rtmAll_obj%z(K), C*SUM2/N, C*SUM1/N
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 4 is fluxes at every base grid point
        else if (nint(OUTPARMS(1)) .eq. 4) then

            if (src_obj_%unit .eq. 'T') then
                FORM = '(3(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y       Z        UP     DOWN'
            else
                if (src_obj_%doSolar) then
                    FORM = '(3(1X,F7.3),1X,3(2X,E12.5))'
                    write (2, '(A,A)') '!   X       Y       Z',&
                        & '           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(3(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)')&
                        & '!   X       Y       Z           UP           DOWN'
                end if
            end if

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! F4'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtmAll_obj%Nx, ',', rtmAll_obj%Ny, ',', rtmAll_obj%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <F4> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = rtmAll_obj%Nx*rtmAll_obj%Ny*rtmAll_obj%Nz
            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, :, :, :), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1, :, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (src_obj_%doSolar) then
                IDATA = 0
                call shd_io__outputDataBin(FLUXES(3, :, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !C = 1.0
            !if (src_obj_%unit .eq. 'T') C = 1.0 / PI_
            !do I = 1, rtmAll_obj%Nx
            !    do J = 1, rtmAll_obj%Ny
            !        do K = 1, rtmAll_obj%Nz
            !            if (src_obj_%type .ne. 'T') then
            !                write (2, FORM) rtmAll_obj%x(I), rtmAll_obj%y(J), rtmAll_obj%z(K),&
            !                    & C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I), FLUXES(3, K, J, I)
            !            else
            !                write (2, FORM) rtmAll_obj%x(I), rtmAll_obj%y(J), rtmAll_obj%z(K),&
            !                    & C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I)
            !            end if
            !        end do
            !    end do
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 5 is fluxes at every grid point (not supported)
        else

            write(*, '(A)') 'Error [_outputResultMpi]: MPI does NOT support <F5>.'

        end if

    else if (OUTTYPE .eq. 'H') then

        NVAR = 1

        ! Format 1 is domain averaged vertical profile
        if (nint(OUTPARMS(1)) .eq. 1) then

            write (2, '(A,A)') '!    Z       -DIV(Fnet)'

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(SUM1(rtmAll_obj%Nz))

            do IZ = 1, rtmAll_obj%Nz
                SUM1(IZ) = 0.0
                do I = 1, rtmAll_obj%Nx
                    do J = 1, rtmAll_obj%Ny
                        SUM1(IZ) = SUM1(IZ) + FLUXDIV(IZ, J, I)
                    end do
                end do
            end do
            !╰────────────────────────────────────────────────╯!

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! H1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', rtmAll_obj%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <H1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = rtmAll_obj%Nx*rtmAll_obj%Ny
            IDATA = 0
            call shd_io__outputDataBin(-SUM1/N, 1, rtmAll_obj%Nz, IDATA, OUTFILE_DATA)

            deallocate(SUM1)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do K = 1, rtmAll_obj%Nz
            !    SUM = 0.0
            !    do I = 1, rtmAll_obj%Nx
            !        do J = 1, rtmAll_obj%Ny
            !            SUM = SUM + FLUXDIV(K, J, I)
            !        end do
            !    end do
            !    write (2, '(1X,F7.3,2X,E13.5)') rtmAll_obj%z(K), -SUM/(rtmAll_obj%Nx*rtmAll_obj%Ny)
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 2 is flux convergence for every base grid point
        else if (nint(OUTPARMS(1)) .eq. 2) then

            write (2, '(A)') '!    X       Y        Z       -DIV(Fnet)'

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! H2'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtmAll_obj%Nx, ',', rtmAll_obj%Ny, ',', rtmAll_obj%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <H2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = rtmAll_obj%Nx*rtmAll_obj%Ny*rtmAll_obj%Nz
            IDATA = 0
            call shd_io__outputDataBin(-FLUXDIV, 1, N, IDATA, OUTFILE_DATA)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !    do I = 1, rtmAll_obj%Nx
            !        do J = 1, rtmAll_obj%Ny
            !            do K = 1, rtmAll_obj%Nz
            !                write (2, '(3(1X,F7.3),2X,E13.5)')&
            !                    & rtmAll_obj%x(I), rtmAll_obj%y(J), rtmAll_obj%z(K), -FLUXDIV(K, J, I)
            !            end do
            !        end do
            !    end do
            !end if
            !╰────────────────────────────────────────────────╯!

        ! Format 3 is flux convergence for every grid point (not supported)
        else

            write(*, '(A)') 'Error [_outputResultMpi]: MPI does NOT support <H3>.'

        end if


    ! Spherical Harmonic output:
    !   Output mean intensity and net flux (x,y,z components).
    !   and maybe normalized rms of higher order terms
    else if (OUTTYPE .eq. 'S') then

        if (nint(OUTPARMS(1)) .eq. 1) then

            N = rtmAll_obj%Nx*rtmAll_obj%Ny*rtmAll_obj%Nz

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(SHTERM5(rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))

            if (NSHOUT .eq. 5) then

                NVAR = 5
                write (2, '(A,A)') '!    X       Y        Z         Imean',&
                    &'          Fx           Fy           Fz       HOrms'
                do IZ = 1, rtmAll_obj%Nz
                    do I = 1, rtmAll_obj%Nx
                        do J = 1, rtmAll_obj%Ny
                            SHTERM5(IZ, J, I) = SHTERMS(5, IZ, J, I)/max(1e-20, SHTERMS(1, IZ, J, I))
                        end do
                    end do
                end do

            else

                NVAR = 4
                write (2, '(A,A)') '!    X       Y        Z         Imean',&
                    &'          Fx           Fy           Fz '

            end if
            !╰────────────────────────────────────────────────╯!

            ! description for binary data
            !╭────────────────────────────────────────────────╮!
            write (2, '(A)') '! The following contains information for binary data intepretation:'
            write (2, '(A)') '! S1'
            write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
            write (2, '(A, 5(I10,A))') '! ', rtmAll_obj%Nx, ',', rtmAll_obj%Ny, ',', rtmAll_obj%Nz, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <S1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            IDATA = 0
            call shd_io__outputDataBin(SHTERMS(1:4, :, :, :), 4, N, IDATA, OUTFILE_DATA)

            if (NSHOUT .eq. 5) then
                IDATA = 0
                call shd_io__outputDataBin(SHTERM5, 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(SHTERM5)
            !╰────────────────────────────────────────────────╯!

        else

            write(*, '(A)') 'Error [_outputResultMpi]: MPI does NOT support <S2>.'

        end if

        ! original
        !╭────────────────────────────────────────────────╮!
        !if (NSHOUT .eq. 5) then

        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz       HOrms'
        !    do I = 1, rtmAll_obj%Nx
        !        do J = 1, rtmAll_obj%Ny
        !            do K = 1, rtmAll_obj%Nz
        !                write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
        !                    &rtmAll_obj%x(I), rtmAll_obj%y(J), rtmAll_obj%z(K), SHTERMS(1, K, J, I),&
        !                    &SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I),&
        !                    &SHTERMS(5, K, J, I)/max(1e-20, SHTERMS(1, K, J, I))
        !            end do
        !        end do
        !    end do

        !else

        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz '
        !    do I = 1, rtmAll_obj%Nx
        !        do J = 1, rtmAll_obj%Ny
        !            do K = 1, rtmAll_obj%Nz
        !                write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
        !                    &rtmAll_obj%x(I), rtmAll_obj%y(J), rtmAll_obj%z(K), SHTERMS(1, K, J, I),&
        !                    &SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I)
        !            end do
        !        end do
        !    end do

        !end if
        !╰────────────────────────────────────────────────╯!

    ! J source function output: source function for a specified ang
    !   at every grid point (not supported)
    else if (OUTTYPE .eq. 'J') then

        write(*, '(A)') 'Error [_outputResultMpi]: MPI does NOT support <J1> or <J2>.'

    ! Medium output: grid point medium properties -
    !   Output extinction, single scattering albedo,
    !   asymmetry parameter, and temperature.  (not supported)
    else if (OUTTYPE .eq. 'M') then

        write(*, '(A)') 'Error [_outputResultMpi]: MPI does NOT support <M1> or <M2>.'

    end if

    close (2)
    return
end subroutine shd_io__outputResultMpi

subroutine shd_io__outputCamRad(&
    &INSAVEFILE, OUTSAVEFILE,&
    &DELTAM,&
    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    &IVIS, VISOUT, OUTPARMS, OUTFILE)

    ! Writes the visualization image in PGM binary format with either
    ! one byte or two byte pixels.  The PGM binary format (P5) has an
    ! ascii header which has the image parameters amd can include comments
    ! (here the SHDOM input parameters), followed by the byte or two byte
    ! image pixels.

    integer MAXITER, TOTITER, IVIS
    logical, intent(in) :: DELTAM
    real SOLACC, SPLITACC, SHACC
    real OUTPARMS(*)
    character(len=*), intent(in) :: INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real VISOUT(shdo_obj_%Nstokes, *)

    character(len=128) :: GRIDNAME, SOURCENAME, UNITSNAME, SFCNAME, OUTNAME
    character(len=charMaxLen_) :: OUTFILE_DATA
    character(len=20)  :: DATA_POSTFIX
    integer :: N, IDATA, NL, NS, CAMID

    DATA_POSTFIX = '.sHdOmNG-out'
    OUTFILE_DATA = trim(adjustl(OUTFILE))//trim(adjustl(DATA_POSTFIX))

    GRIDNAME = 'EVEN (X,Y)  '
    if (rtm_obj_%type .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (rtm_obj_%type .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (src_obj_%type .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (src_obj_%type .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (src_obj_%type .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (sfc_obj_%type .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (sfc_obj_%type .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (sfc_obj_%type .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (sfc_obj_%type .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (sfc_obj_%type .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (sfc_obj_%type .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (sfc_obj_%type .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (sfc_obj_%type(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (src_obj_%unit .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (src_obj_%unit .eq. 'B') then
        UNITSNAME = 'WATTS/(M^2 STER)'
    else
        UNITSNAME = 'WATTS/(M^2 MICRON STER)'
    end if


    CAMID = nint(OUTPARMS(1))
    if (CAMID .eq. 1) then ! camera mode
        NL = nint(OUTPARMS(10))
        NS = nint(OUTPARMS(11))
        OUTNAME = 'CAMERA RADIANCE (V1 Camera Mode)'
    else if (CAMID .eq. 2) then ! cross track mode
        NL = 1 + sqrt((OUTPARMS(4) - OUTPARMS(7))**2&
                  & + (OUTPARMS(5) - OUTPARMS(8))**2&
                  & + (OUTPARMS(6) - OUTPARMS(9))**2)/OUTPARMS(10)
        NS = 1 + abs(OUTPARMS(12) - OUTPARMS(11))/OUTPARMS(13)
        OUTNAME = 'CAMERA RADIANCE (V2 Crosstrack Mode)'
    else if (CAMID .eq. 3) then ! user-specified sensor mode
        NL = 1
        NS = nint(OUTPARMS(2))
        OUTNAME = 'CAMERA RADIANCE (V3 User-Specified Sensor Mode)'
    else
        write(*, '(A)') 'Error [_outputCamRad]: Unkown camera mode.'
        stop
    end if

    open (UNIT=2, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN')
    write (2, '(A,A)')&
        &'! Polarized Spherical Harmonic Discrete Ordinate',&
        &' Radiative Transfer Output'
    write (2, '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)')&
        &'!  L=', shdo_obj_%maxNl, '  M=', shdo_obj_%maxNm, '  NLM=', shdo_obj_%Nlm,&
        &'   NMU=', shdo_obj_%Nmu, '  NPHI=', shdo_obj_%Nphi, '  NANG=', shdo_obj_%Nang,&
        &'   NSH=', shdo_obj_%Nsh
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', shdo_obj_%Nstokes,&
        &'   NX=', rtm_obj_%Nx, '   NY=', rtm_obj_%Ny, '   NZ=', rtm_obj_%Nz,&
        &'    NPTS=', rtm_obj_%Npts, '   NCELLS=', cell_obj_%Npts
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(prp_obj_%fname)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',&
        &trim(ckd_obj_%fname), '   NUM_G=', ckd_obj_%Ng
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)

    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,&
            &'      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if

    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,&
        &'   INDEPENDENT_PIXEL=', rtm_obj_%flag_ip
    write (2, '(A,A22,A,I1)') '!  SURFACE_TYPE=', SFCNAME,&
        &'   HORIZ_BOUNDARY_COND=', rtm_obj_%flag_bc

    if (src_obj_%doSolar) then

        if (sfc_obj_%type(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(sfc_obj_%fname),&
                &'  SKY_RAD=', src_obj_%skyRad
        else
            if (src_obj_%type .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', sfc_obj_%avgAlb, '  GROUND_TEMP=', sfc_obj_%avgTemp,&
                    &'  SKY_RAD=', src_obj_%skyRad
            else
                write (2, '(A,F9.7,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', sfc_obj_%avgAlb, '  SKY_RAD=', src_obj_%skyRad
            end if
        end if

        write (2, '(A,E13.6,A,F10.7,A,F8.3)')&
            &'!  SOLAR_FLUX=', src_obj_%avgSolar, '   SOLAR_MU=', src_obj_%avgMu,&
            &'   SOLAR_AZ=', src_obj_%avgAz*180.0/acos(-1.0)

    else

        if (sfc_obj_%type(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(sfc_obj_%fname)
        else
            write (2, '(A,F8.3,A,F9.7)')&
                &'!  GROUND_TEMP=', sfc_obj_%avgTemp, '  GROUND_EMIS=', 1.0 - sfc_obj_%avgAlb
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', src_obj_%skyRad

    end if

    if (src_obj_%unit .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (src_obj_%unit .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,&
            &'   WAVENUMBER_RANGE=', src_obj_%waveNum(1), src_obj_%waveNum(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,&
            &'   WAVELENGTH=', src_obj_%avgWaveLen
    end if

    write (2, '(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,&
        &'   SPHERICAL_HARMONIC_ACCURACY=', SHACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER,&
        &'   NUMBER_ITERATIONS=', TOTITER
    write (2, '(A,A)') '!  OUTPUT_TYPE=', trim(adjustl(OUTNAME))

    if (OUTPARMS(1) .gt. 0.0) then ! placeholder for treating camera mode differently if needed

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following provides information for interpreting binary data:'
        write (2, '(A)') '! R'
        write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
        write (2, '(A, 5(I10,A))') '! ', NS, ',', NL, ',', 1, ',', 1, ',', shdo_obj_%Nstokes, ''
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (verbose_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <V> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        N = NL*NS
        IDATA = 0
        call shd_io__outputDataBin(VISOUT, shdo_obj_%Nstokes, N, IDATA, OUTFILE_DATA)
        !╰────────────────────────────────────────────────╯!

    end if

    close (2)
    return
end subroutine shd_io__outputCamRad

subroutine shd_io__outputImg(&
    &INSAVEFILE, OUTSAVEFILE,&
    &DELTAM,&
    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    &IVIS, VISOUT, OUTPARMS, OUTFILE)

    ! Writes the visualization image in PGM binary format with either
    ! one byte or two byte pixels.  The PGM binary format (P5) has an
    ! ascii header which has the image parameters amd can include comments
    ! (here the SHDOM input parameters), followed by the byte or two byte
    ! image pixels.

    integer MAXITER, TOTITER, IVIS
    logical, intent(in) :: DELTAM
    real SOLACC, SPLITACC, SHACC
    real OUTPARMS(*)
    character(len=*), intent(in) :: INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real VISOUT(shdo_obj_%Nstokes, *)

    integer MAXHEAD, NL, NS, NLOUT, LABRECS, RECLEN, NBYTES
    integer DN, I, J, K, L, M
    real MAXLINPOL, MINANGPOL, MAXANGPOL, MAXCIRCPOL
    real VISMIN(4), VISMAX(4), VISSCALE(4)
    logical CAMERA_MODE, LSBFIRST
    integer*1, allocatable :: BUFFER1(:)
    character(len=128) :: GRIDNAME, SOURCENAME, UNITSNAME, SFCNAME
    character(len=128) :: POLNAMES(4), VISFORM(4)
    character(len=1) :: LF
    character(len=2000) :: BUF

    MAXHEAD = 2000
    do I = 1, MAXHEAD
        BUF(I:I) = ' '
    end do

    ! Find the image size
    CAMERA_MODE = nint(OUTPARMS(1)) .eq. 1
    if (CAMERA_MODE) then
        NL = nint(OUTPARMS(10))
        NS = nint(OUTPARMS(11))
    else if (OUTPARMS(1) .eq. 2) then
        NL = 1 + sqrt((OUTPARMS(4) - OUTPARMS(7))**2&
                  & + (OUTPARMS(5) - OUTPARMS(8))**2&
                  & + (OUTPARMS(6) - OUTPARMS(9))**2)/OUTPARMS(10)
        NS = 1 + abs(OUTPARMS(12) - OUTPARMS(11))/OUTPARMS(13)
    end if
    NBYTES = min(2, max(1, nint(abs(OUTPARMS(2)))))
    LSBFIRST = OUTPARMS(2) .lt. 0.0
    allocate(BUFFER1(NBYTES*NS))
    NLOUT = shdo_obj_%Nstokes*NL

    ! Calculate the output range for each Stokes parameter and setup fo
    POLNAMES(1) = 'RADIANCE'
    POLNAMES(2) = 'DEGREE_OF_LINEAR_POLARIZATION'
    POLNAMES(3) = 'ANGLE_OF_LINEAR_POLARIZATION'
    POLNAMES(4) = 'DEGREE_OF_CIRCULAR_POLARIZATION'
    VISFORM(:) = '(3A,2(1X,F8.5),A1)'
    VISFORM(3) = '(3A,2(1X,F8.3),A1)'
    VISMIN(1) = 0.0
    VISSCALE(1) = OUTPARMS(3)
    if (NBYTES .eq. 1) then
        VISMAX(1) = 255./OUTPARMS(3)
        if (isPol_) then
            MAXLINPOL = max(0.01, maxval(VISOUT(2, IVIS + 1:IVIS + NL*NS)))
            MINANGPOL = minval(VISOUT(3, IVIS + 1:IVIS + NL*NS))
            MAXANGPOL = maxval(VISOUT(3, IVIS + 1:IVIS + NL*NS))
            VISMIN(2) = 0.0
            VISMAX(2) = MAXLINPOL
            VISSCALE(2) = 255./MAXLINPOL
            VISMIN(3) = MINANGPOL
            VISMAX(3) = max(MAXANGPOL, MINANGPOL + 1.0)
            VISSCALE(3) = 255./(VISMAX(3) - VISMIN(3))
        end if
        if (isPol4_) then
            MAXCIRCPOL = max(0.01, maxval(abs(VISOUT(4, IVIS + 1:IVIS + NL*NS))))
            VISMIN(4) = -MAXCIRCPOL
            VISMAX(4) = MAXCIRCPOL
            VISSCALE(4) = 255./(2*MAXCIRCPOL)
        end if
    else
        VISMAX(1) = 65535./OUTPARMS(3)
        VISMIN(2) = 0.0
        VISMAX(2) = 1.0
        VISSCALE(2) = 65535.
        MINANGPOL = minval(VISOUT(3, IVIS + 1:IVIS + NL*NS))
        if (MINANGPOL .lt. -1.0) then
            VISMIN(3) = -90.
            VISMAX(3) = 90.
        else
            VISMIN(3) = 0.
            VISMAX(3) = 180.
        end if
        VISSCALE(3) = 65535./180.
        VISMIN(4) = -1.0
        VISMAX(4) = 1.0
        VISSCALE(4) = 65535./2.0
    end if

    ! Open the PGM image, which has to be direct access to output the
    ! clean binary pixel data
    RECLEN = NS*NBYTES
    open (UNIT=2, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN',&
        & ACCESS='DIRECT', FORM='UNFORMATTED', RECL=RECLEN)

    ! Start the PGM header in the buffer
    LF = char(10)
    I = 1
    write (BUF(I:I + 2), '(A2,A1)') 'P5', LF
    I = I + 3

    ! Make the SHDOM part of the header
    GRIDNAME = 'EVEN (X,Y)  '
    if (rtm_obj_%type .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (rtm_obj_%type .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (src_obj_%type .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (src_obj_%type .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (src_obj_%type .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (sfc_obj_%type .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (sfc_obj_%type .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (sfc_obj_%type .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (sfc_obj_%type .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (sfc_obj_%type .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (sfc_obj_%type .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (sfc_obj_%type .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (sfc_obj_%type .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (sfc_obj_%type(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (src_obj_%unit .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (src_obj_%unit .eq. 'B') then
        UNITSNAME = 'WATTS/(M^2 STER)'
    else
        UNITSNAME = 'WATTS/(M^2 MICRON STER)'
    end if

    write (BUF(I:I + 79), '(A,A1)')&
        & '# Polarized SHDOM Visualization Output', LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10,A1)')&
        & '# L=', shdo_obj_%maxNl, '  M=', shdo_obj_%maxNm, '  NLM=', shdo_obj_%Nlm,&
        & '   NMU=', shdo_obj_%Nmu, '  NPHI=', shdo_obj_%Nphi, '  NANG=', shdo_obj_%Nang,&
        & '   NSH=', shdo_obj_%Nsh, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 90), '(A,I1,3(A,I4),2(A,I8),A1)')&
        & '# NSTOKES=', shdo_obj_%Nstokes, '  NX=', rtm_obj_%Nx, '   NY=', rtm_obj_%Ny, '  NZ=', rtm_obj_%Nz,&
        & '    NPTS=', rtm_obj_%Npts, '   NCELLS=', cell_obj_%Npts, LF
    I = I + index(BUF(I:I + 79), LF)
    !write (BUF(I:I + 79), '(A,A,A1)')&
    !    & '# PROPERTY_FILE=', trim(prp_obj_%fname), LF
    I = I + index(BUF(I:I + 79), LF)
    !write (BUF(I:I + 79), '(A,A,A,I2,A1)')&
    !    & '# CORRELATED_K-DIST_FILE=', trim(ckd_obj_%fname), '   NUM_G=', ckd_obj_%Ng, LF
    I = I + index(BUF(I:I + 79), LF)
    !write (BUF(I:I + 79), '(A,A,A1)')&
    !    & '# INPUT_SAVE_FILE=', trim(INSAVEFILE), LF
    I = I + index(BUF(I:I + 79), LF)
    !write (BUF(I:I + 79), '(A,A,A1)')&
    !    & '# OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE), LF
    I = I + index(BUF(I:I + 79), LF)
    if (DELTAM) then
        write (BUF(I:I + 79), '(A,A14,A,A1)')&
            & '# SOURCE_TYPE=', SOURCENAME, '    DELTA-M METHOD', LF
        I = I + index(BUF(I:I + 79), LF)
    else
        write (BUF(I:I + 79), '(A,A14,A1)')&
            & '# SOURCE_TYPE=', SOURCENAME, LF
        I = I + index(BUF(I:I + 79), LF)
    end if
    write (BUF(I:I + 79), '(A,A32,A,I1,A1)') '# GRID_TYPE=', GRIDNAME,&
        & '   INDEPENDENT_PIXEL=', rtm_obj_%flag_ip, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A22,A,I2,A1)')&
        & '# SURFACE_TYPE=', SFCNAME, '   HORIZ_BOUNDARY_COND=', rtm_obj_%flag_bc, LF
    I = I + index(BUF(I:I + 79), LF)
    if (src_obj_%doSolar) then
        if (sfc_obj_%type(1:1) .eq. 'V') then
            !write (BUF(I:I + 79), '(A,A,A,E12.5,A1)')&
            !    & '# SURFACE_FILE=', trim(sfc_obj_%fname),&
            !    & '  SKY_RAD=', src_obj_%skyRad, LF
            I = I + index(BUF(I:I + 79), LF)
        else
            if (src_obj_%type .eq. 'B') then
                write (BUF(I:I + 79), '(A,F9.7,A,F8.3,A,E12.5,A1)')&
                    &'# GROUND_ALBEDO=', sfc_obj_%avgAlb, '  GROUND_TEMP=', sfc_obj_%avgTemp,&
                    &'  SKY_RAD=', src_obj_%skyRad, LF
                I = I + index(BUF(I:I + 79), LF)
            else
                write (BUF(I:I + 79), '(A,F9.7,A,E12.5,A1)')&
                    &'# GROUND_ALBEDO=', sfc_obj_%avgAlb, '  SKY_RAD=', src_obj_%skyRad, LF
                I = I + index(BUF(I:I + 79), LF)
            end if
        end if
        write (BUF(I:I + 79), '(A,E13.6,A,F10.7,A,F8.3,A1)')&
            &'# SOLAR_FLUX=', src_obj_%avgSolar, '   SOLAR_MU=', src_obj_%avgMu,&
            &'   SOLAR_AZ=', src_obj_%avgAz*180.0/acos(-1.0), LF
        I = I + index(BUF(I:I + 79), LF)
    else
        if (sfc_obj_%type(1:1) .eq. 'V') then
            !write (BUF(I:I + 79), '(A,A,A1)')&
            !    &'# SURFACE_FILE=', trim(sfc_obj_%fname), LF
            I = I + index(BUF(I:I + 79), LF)
        else
            write (BUF(I:I + 79), '(A,F8.3,A,F9.7,A1)')&
                &'# GROUND_TEMP=', sfc_obj_%avgTemp,&
                &'  GROUND_EMIS=', 1.0 - sfc_obj_%avgAlb, LF
            I = I + index(BUF(I:I + 79), LF)
        end if
        write (BUF(I:I + 79), '(A,F8.3,A1)') '# SKY_TEMP=', src_obj_%skyRad, LF
        I = I + index(BUF(I:I + 79), LF)
    end if

    if (src_obj_%unit .eq. 'T') then
        write (BUF(I:I + 79), '(A,A24,A1)') '# UNITS=', UNITSNAME, LF
        I = I + index(BUF(I:I + 79), LF)
    else if (src_obj_%unit .eq. 'B') then
        write (BUF(I:I + 79), '(A,A24,A,2(1X,F10.2),A1)')&
            &'# UNITS=', UNITSNAME,&
            &'   WAVENUMBER_RANGE=', src_obj_%waveNum(1), src_obj_%waveNum(2), LF
        I = I + index(BUF(I:I + 79), LF)
    else
        write (BUF(I:I + 79), '(A,A24,A,F10.2,A1)')&
            &'# UNITS=', UNITSNAME, '   WAVELENGTH=', src_obj_%avgWaveLen, LF
        I = I + index(BUF(I:I + 79), LF)
    end if
    write (BUF(I:I + 79), '(2(A,E10.3),A1)')&
        &'# SPLITTING_ACCURACY=', SPLITACC,&
        &'   SPHERICAL_HARMONIC_ACCURACY=', SHACC, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,E10.3,A1)')&
        &'# SOLUTION_ACCURACY=', SOLACC, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(2(A,I4),A1)')&
        &'# MAXIMUM_ITERATIONS=', MAXITER,&
        &'   NUMBER_ITERATIONS=', TOTITER, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 7), '(A,A1)') '# ', LF
    I = I + index(BUF(I:I + 7), LF)

    if (CAMERA_MODE) then
        write (BUF(I:I + 79), '(A,A1)') '# CAMERA MODE:', LF
        I = I + index(BUF(I:I + 79), LF)
        ! 1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline,
        write (BUF(I:I + 79), '(A,F9.2,A1)')&
            &'# RADIANCE_SCALE=', OUTPARMS(3), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,3(1X,F9.3),A1)') '# CAMERA_LOCATION=',&
            &OUTPARMS(4), OUTPARMS(5), OUTPARMS(6), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F10.5,A1)')&
            &'# CAMERA_ZENITH_ANGLE=', OUTPARMS(7), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F10.5,A1)')&
            &'# CAMERA_AZIMUTH_ANGLE=', OUTPARMS(8), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F10.5,A1)')&
            &'# CAMERA_ROTATION_ANGLE=', OUTPARMS(9), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,I5,A1)')&
            &'# CAMERA_NUMBER_LINES=', nint(OUTPARMS(10)), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,I5,A1)')&
            &'# CAMERA_NUMBER_SAMPLES=', nint(OUTPARMS(11)), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.5,A1)')&
            &'# CAMERA_LINE_SPACING_ANGLE=', OUTPARMS(12), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.5,A1)')&
            &'# CAMERA_SAMPLE_SPACING_ANGLE=', OUTPARMS(13), LF
        I = I + index(BUF(I:I + 79), LF)

    else
        write (BUF(I:I + 79), '(A,A1)')&
            &'# CROSS TRACK SCANNING MODE:', LF
        I = I + index(BUF(I:I + 79), LF)
        ! 2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, d
        write (BUF(I:I + 79), '(A,F9.2,A1)')&
            &'# RADIANCE_SCALE=', OUTPARMS(3), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,3(1X,F9.3),A1)') '# START_LOCATION=',&
            &OUTPARMS(4), OUTPARMS(5), OUTPARMS(6), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,3(1X,F9.3),A1)') '# END_LOCATION=',&
            &OUTPARMS(7), OUTPARMS(8), OUTPARMS(9), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.4,A1)')&
            &'# SPACING_ON_TRACK=', OUTPARMS(10), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.3,A1)')&
            &'# START_SCAN_ANGLE=', OUTPARMS(11), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.3,A1)')&
            &'# END_SCAN_ANGLE=', OUTPARMS(12), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F7.4,A1)')&
            &'# SCAN_ANGLE_SPACING=', OUTPARMS(13), LF
        I = I + index(BUF(I:I + 79), LF)
    end if
    do M = 1, shdo_obj_%Nstokes
        write (BUF(I:I + 79), VISFORM(M)) '# ', trim(POLNAMES(M)),&
            &'_RANGE_IN_PIXEL=', VISMIN(M), VISMAX(M), LF
        I = I + index(BUF(I:I + 79), LF)
    end do

    ! Make the final part of PGM header at the end of the last label r
    write (BUF(I:I + 1), '(A2)') '# '
    I = I + 2
    LABRECS = int((I + 19 - 1)/RECLEN) + 1
    I = RECLEN*LABRECS - 19 + 1
    write (BUF(I:I + 12), '(A1,I5,1X,I5,A1)') LF, NS, NLOUT, LF
    I = I + 13
    if (NBYTES .eq. 1) then
        write (BUF(I:I + 5), '(I5,A1)') 255, LF
    else
        write (BUF(I:I + 5), '(I5,A1)') 65535, LF
    end if

    ! Write the header records
    do I = 1, LABRECS
        write (2, REC=I) BUF(1 + RECLEN*(I - 1):RECLEN*I)
    end do

    ! Output the image, converting from real to byte or integer
    L = 0
    do M = 1, shdo_obj_%Nstokes
        K = IVIS + 1
        do J = 1, NL
            L = L + 1
            if (NBYTES .eq. 1) then
                do I = 1, NS
                    DN = nint(VISSCALE(M)*(VISOUT(M, K) - VISMIN(M)))
                    DN = max(0, min(255, DN))
                    BUFFER1(I) = DN
                    ! if (DN .GT. 127) then
                    !     BUFFER1(I) = DN - 256
                    ! else
                    !     BUFFER1(I) = DN
                    ! endif
                    K = K + 1
                end do
                write (2, REC=LABRECS + L) BUFFER1(1:NS)
            else
                do I = 1, NS
                    DN = nint(VISSCALE(M)*(VISOUT(M, K) - VISMIN(M)))
                    DN = max(0, min(65535, DN))
                    if (LSBFIRST) then
                        BUFFER1(2*I - 1) = iand(DN, 255)
                        BUFFER1(2*I) = ishft(DN, -8)
                    else
                        BUFFER1(2*I - 1) = ishft(DN, -8)
                        BUFFER1(2*I) = iand(DN, 255)
                    end if
                    K = K + 1
                end do
                write (2, REC=LABRECS + L) BUFFER1(1:2*NS)
            end if
        end do
    end do
    IVIS = IVIS + NL*NS

    deallocate(BUFFER1)
    close (2)
    return
end subroutine shd_io__outputImg

subroutine shd_io__outputCellSplit(&
    &OUTFILE,&
    &SHPTR, SOURCE)

    ! Outputs the final cell splitting criterion for each grid cell
    ! and direction (X,Y,Z).  Calls CELL_SPLIT_TEST.
    ! Mainly for debugging purposes.

    integer SHPTR(*)
    real SOURCE(shdo_obj_%Nstokes, *)
    character(len=*), intent(in) :: OUTFILE
    integer IP1, IP2, icell, IDIR
    real ADAPTCRIT(3), MAXADAPT, XC, YC, ZC, EXT

    open (UNIT=9, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN')
    write (9, *) '!  icell   Xc     Yc     Zc  ',&
        &'Extinct SplitX  SplitY  SplitZ'

    ! Loop over all the cells
    do icell = 1, cell_obj_%Npts
        ! Get the adaptive grid criterion for the 3 directions
        call shd_grid__testCellSplit(icell,&
            &ADAPTCRIT, MAXADAPT, IDIR)
        ! Compute the cell center location and mean extinction
        IP1 = i8Grid_Cell_(1, icell)
        IP2 = i8Grid_Cell_(8, icell)
        XC = 0.5*(xyz_Rtm_(1, IP1) + xyz_Rtm_(1, IP2))
        YC = 0.5*(xyz_Rtm_(2, IP1) + xyz_Rtm_(2, IP2))
        ZC = 0.5*(xyz_Rtm_(3, IP1) + xyz_Rtm_(3, IP2))
        EXT =  (ext_Rtm_(i8Grid_Cell_(1, icell)) + ext_Rtm_(i8Grid_Cell_(2, icell))&
             &+ ext_Rtm_(i8Grid_Cell_(3, icell)) + ext_Rtm_(i8Grid_Cell_(4, icell))&
             &+ ext_Rtm_(i8Grid_Cell_(5, icell)) + ext_Rtm_(i8Grid_Cell_(6, icell))&
             &+ ext_Rtm_(i8Grid_Cell_(7, icell)) + ext_Rtm_(i8Grid_Cell_(8, icell)))/8
        write (9, '(1X,I5,3(1X,F6.3),1X,E11.4,3(1X,F7.5))')&
            &icell, XC, YC, ZC, EXT,&
            &ADAPTCRIT(1), ADAPTCRIT(2), ADAPTCRIT(3)
    end do

    close (9)
    return
end subroutine shd_io__outputCellSplit

subroutine shd_io__visualizeCell(&
    &OUTFILE, YVAL, OUTTYPE,&
    &XDOMAIN, ZDOMAIN)

    ! Visualizes the cell structure by outputting GLE code to draw
    ! a diagram of the cells in the XZ plane.  Only the end cells are
    ! output.  The outlines of the cells are drawn, and arrows between
    ! the cell centers show the neighbor relationships.
    ! OUTTYPE=0 for cell outlines only, OUTTYPE=1 for cells and arrows.
    ! There is a special output for IP mode.

    integer OUTTYPE
    real YVAL, XDOMAIN, ZDOMAIN
    character(len=*), intent(in) :: OUTFILE
    integer IC
    integer J, INeigh
    real XSIZE, ZSIZE, S, DELX
    real X1, Z1, X2, Z2, X3, Z3, X4, Z4, XC, ZC
    real DX, DZ, XN, ZN, XA, ZA

    S = min(24.0/XDOMAIN, 18.0/ZDOMAIN)
    XSIZE = S*XDOMAIN
    ZSIZE = S*ZDOMAIN
    open (UNIT=9, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN')
    write (9, '(A,2(1X,F6.2))') 'size ', XSIZE, ZSIZE
    write (9, '(A)') 'set just cc'
    write (9, '(A)') 'amove 0 0'
    write (9, '(A)') 'begin origin'
    write (9, '(A)') 'begin scale 1 1'
    write (9, '(A,1X,F7.3)') '! Y value: ', YVAL
    write (9, '(A)')
    if (is1DX_) then
        DELX = 0.5*XDOMAIN/rtm_obj_%Nx
        write (9, '(A)') 'sub drawbox x1 y1 x3 y3 xc yc n'
        write (9, '(A)') ' set color black'
        write (9, '(A)') ' amove x1 y1'
        write (9, '(A,2F7.3)') ' rline ', S*DELX, 0
        write (9, '(A)') ' amove x3 y3'
        write (9, '(A,2F7.3)') ' rline ', S*DELX, 0
        write (9, '(A)') ' amove xc yc'
        write (9, '(A)') ' set hei 0.3'
        write (9, '(A)') ' write num$(n)'
        write (9, '(A)') 'end sub'
    else
        write (9, '(A)') 'sub drawbox x1 y1 x2 y2 x3 y3 x4 y4  xc yc n'
        write (9, '(A)') ' set color black'
        write (9, '(A)') ' amove x1 y1'
        write (9, '(A)') ' aline x2 y2'
        write (9, '(A)') ' aline x3 y3'
        write (9, '(A)') ' aline x4 y4'
        write (9, '(A)') ' aline x1 y1'
        write (9, '(A)') ' amove xc yc'
        write (9, '(A)') ' set hei 0.3'
        write (9, '(A)') ' write num$(n)'
        write (9, '(A)') 'end sub'
    end if
    write (9, '(A)')
    if (OUTTYPE .eq. 1) then
        write (9, '(A)') 'sub drawarrow x1 y1 x2 y2 i'
        write (9, '(A)') ' if (i<0) then'
        write (9, '(A)') '   set color red'
        write (9, '(A)') ' else'
        write (9, '(A)') '   set color green'
        write (9, '(A)') ' end if'
        write (9, '(A)') ' amove x1 y1'
        write (9, '(A)') ' set hei 0.2'
        write (9, '(A)') ' aline x2 y2 arrow end'
        write (9, '(A)') ' set color black'
        write (9, '(A)') 'end sub'
        write (9, '(A)')
    end if

    do IC = 1, cell_obj_%Npts
        ! Output if it is an end cell and crosses the Y value
        if (i2Tree_Cell_(2, IC) .eq. 0 .and.&
            & xyz_Rtm_(2, i8Grid_Cell_(1, IC)) .le. YVAL .and.&
            & xyz_Rtm_(2, i8Grid_Cell_(3, IC)) .ge. YVAL) then
            ! Get the center of the cell
            XC = (xyz_Rtm_(1, i8Grid_Cell_(1, IC)) + xyz_Rtm_(1, i8Grid_Cell_(8, IC)))/2
            ZC = (xyz_Rtm_(3, i8Grid_Cell_(1, IC)) + xyz_Rtm_(3, i8Grid_Cell_(8, IC)))/2
            ! Draw the outline in XZ (four lines and four grid points)
            !   and the cell number at the center using the gle subrouti
            X1 = xyz_Rtm_(1, i8Grid_Cell_(1, IC))
            Z1 = xyz_Rtm_(3, i8Grid_Cell_(1, IC))
            X2 = xyz_Rtm_(1, i8Grid_Cell_(2, IC))
            Z2 = xyz_Rtm_(3, i8Grid_Cell_(2, IC))
            X3 = xyz_Rtm_(1, i8Grid_Cell_(6, IC))
            Z3 = xyz_Rtm_(3, i8Grid_Cell_(6, IC))
            X4 = xyz_Rtm_(1, i8Grid_Cell_(5, IC))
            Z4 = xyz_Rtm_(3, i8Grid_Cell_(5, IC))
            if (is1DX_) then
                write (9, '(A,6(1X,F6.3),I8)') '@drawbox ',&
                    & S*X1, S*Z1, S*X3, S*Z3, S*XC, S*ZC, IC
            else
                write (9, '(A,10(1X,F6.3),I8)') '@drawbox ',&
                    & S*X1, S*Z1, S*X2, S*Z2, S*X3, S*Z3, S*X4, S*Z4, S*XC, S*ZC, IC
            end if
            if (OUTTYPE .eq. 1) then
                ! Draw the neighbor arrows
                DX = 0.1*(X2 - X1)
                DZ = 0.1*(Z4 - Z1)
                do J = 1, 6
                    if (i6Cell_Cell_(J, IC) .ne. 0 .and. J .ne. 3 .and. J .ne. 4) then
                        INeigh = abs(i6Cell_Cell_(J, IC))
                        XN = (xyz_Rtm_(1, i8Grid_Cell_(1, INeigh)) + xyz_Rtm_(1, i8Grid_Cell_(8, INeigh)))/2
                        ZN = (xyz_Rtm_(3, i8Grid_Cell_(1, INeigh)) + xyz_Rtm_(3, i8Grid_Cell_(8, INeigh)))/2
                        XA = XC + 1.0*(XN - XC)
                        ZA = ZC + 1.0*(ZN - ZC)
                        if (.not. (XN - XC .gt. 0.0 .and. J .eq. 1 .or.&
                            & XN - XC .lt. 0.0 .and. J .eq. 2)) then
                            write (9, '(A,4(1X,F6.3),1X,I8)') '@drawarrow ',&
                                & S*XC, S*ZC, S*XA, S*ZA, i6Cell_Cell_(J, IC)
                        end if
                    end if
                end do
            end if

        end if
    end do
    write (9, '(A)') 'end scale'
    write (9, '(A)') 'end origin'

    close (9)
    return
end subroutine shd_io__visualizeCell

subroutine shd_io__saveState(&
    &OUTSAVEFILE, WORK, FLUXES)

    ! Write a binary save file containing the current state of
    ! the SHDOM solution. This includes the cell tree structure,
    ! the spherical harmonic expansion of the source function,
    ! the low order terms of the radiance expansion, and the flux array.
    ! If the file name is 'NONE' then no file is written.

    character(len=*) , intent(in) :: OUTSAVEFILE
    real(kind=R__)   , intent(inout) :: WORK(*)
    real(kind=R__)   , intent(in) :: FLUXES(2, rtm_obj_%Npts)

    integer(kind=4) :: I, J, K, IR, NS, NR, istoke
    integer(kind=4) :: DATA_BYTE
    real(kind=4) :: DATA_SIZE

    if (OUTSAVEFILE(1:4) .eq. 'NONE') return

    open (unit=91, file=trim(adjustl(OUTSAVEFILE)), status='REPLACE',&
        &form='UNFORMATTED')

    write (91) shdo_obj_%Nstokes, rtm_obj_%Nx, rtm_obj_%Ny, rtm_obj_%Nz, shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm, rtm_obj_%Npts, cell_obj_%Npts

    write (91) (rtm_obj_%x(I), I=1, rtm_obj_%Nx+1)
    write (91) (rtm_obj_%y(I), I=1, rtm_obj_%Ny+1)
    write (91) (rtm_obj_%z(I), I=1, rtm_obj_%Nz)

    write (91) ((i8Grid_Cell_(J, I), J=1, 8), I=1, cell_obj_%Npts)
    write (91) ((i6Cell_Cell_(J, I), J=1, 6), I=1, cell_obj_%Npts)
    write (91) ((i2Tree_Cell_(J, I), J=1, 2), I=1, cell_obj_%Npts)
    write (91)  (flag_Cell_(I), I=1, cell_obj_%Npts)

    write (91) ((xyz_Rtm_(J, I), J=1, 3), I=1, rtm_obj_%Npts)
    write (91) ((FLUXES(J, I) , J=1, 2), I=1, rtm_obj_%Npts)

    write (91) (isrcJ_Shdo_(I), I=1, rtm_obj_%Npts + 1)

    NS = isrcJ_Shdo_(rtm_obj_%Npts + 1)
    write (91) ((srcJ_Shdo_(J, I), J=1, shdo_obj_%Nstokes), I=1, NS)

    !K = 1
    !do I = 1, rtm_obj_%Npts
    !    IR = iradI_Shdo_(I)
    !    do J = 1, shdo_obj_%Nstokes*4

    !        ! the following is an attempt to fix the segmentation issue
    !        !╭────────────────────────────────────────────────╮!
    !        ! iradI_Shdo_(1) = 0 (from shd_init__makeRad)
    !        ! originally it was:
    !        ! WORK(K) = radI_Shdo_(shdo_obj_%Nstokes*(IR-1)+J)
    !        ! when I=1, IR=0, shdo_obj_%Nstokes*(IR-1)+J is 0 when J=1 (shdo_obj_%Nstokes=1)
    !        WORK(K) = radI_Shdo_(shdo_obj_%Nstokes*(IR)+J)
    !        !╰────────────────────────────────────────────────╯!

    !        K = K + 1
    !    end do
    !end do

    K = 1
    do I = 1, rtm_obj_%Npts
        IR = iradI_Shdo_(I)   !radI(NSTOKES, NPTS)
        do J = 1, 4
            do istoke = 1, shdo_obj_%Nstokes
                WORK(K) = radI_Shdo_(istoke, IR+J)
                K = K + 1
            end do
        end do
    end do

    NR = 4*rtm_obj_%Npts
    write (91) (WORK(I), I=1, shdo_obj_%Nstokes*NR)

    close (unit=91)

    ! display file size
    !╭────────────────────────────────────────────────╮!
    inquire(FILE=trim(adjustl(OUTSAVEFILE)), size=DATA_BYTE)
    DATA_SIZE = DATA_BYTE/1000000.0
    if (verbose_ .ge. 5) write(*, '(A,F6.2,A)') 'Message [_saveState]: Data saved at <'//trim(adjustl(OUTSAVEFILE))//'> (',&
        &DATA_SIZE, 'MB).'
    !╰────────────────────────────────────────────────╯!

    return

end subroutine shd_io__saveState

subroutine shd_io__outputDataBin(&
    &DATA_OUT, NVAR, NDATA,&
    &IDATA, OUTFILE_DATA, MODE)

    ! Writes data into binary file

    real, intent(in) :: DATA_OUT(NVAR, *)
    integer :: NVAR, NDATA
    integer, intent(inout) :: IDATA
    character(len=*), intent(in) :: OUTFILE_DATA
    character(len=*), intent(in), optional :: MODE
    integer :: ivar, icount, idata_start, istat
    integer :: DATA_BYTE
    real :: DATA_SIZE

    if ((present(MODE)) .and. (MODE .eq. 'ADD')) then
        open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='OLD', position='APPEND',&
            &form='UNFORMATTED', access='STREAM', iostat=istat)
    else
        open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='REPLACE',&
            &form='UNFORMATTED', access='STREAM', iostat=istat)
    end if

    idata_start = IDATA

    if (istat .eq. 0) then

        do ivar = 1, NVAR
            ! rewind for each column variable
            IDATA = idata_start
            do icount = 1, NDATA
                 IDATA = IDATA + 1
                 write (88) DATA_OUT(ivar, IDATA)
            end do
        end do

        close(88, iostat=istat)

        inquire(FILE=trim(adjustl(OUTFILE_DATA)), SIZE=DATA_BYTE)
        DATA_SIZE = DATA_BYTE/1000000.0

        if ((present(MODE)) .and. (MODE .eq. 'ADD')) then
            if (verbose_ .ge. 5) write(*, '(A,F0.2,A)') 'Message [_outputDataBin]: Data added to <'//trim(adjustl(OUTFILE_DATA))//'> (',&
                &DATA_SIZE, 'MB).'
        else
            if (verbose_ .ge. 5) write(*, '(A,F0.2,A)') 'Message [_outputDataBin]: Data saved at <'//trim(adjustl(OUTFILE_DATA))//'> (',&
                &DATA_SIZE, 'MB).'
        end if

    else

        close(88, iostat=istat)

        write(*, '(A)') 'Error [_outputDataBin]: Failed to open file <'//trim(adjustl(OUTFILE_DATA))//'>.'
        stop

    end if

    return

end subroutine shd_io__outputDataBin

end module shd_out
