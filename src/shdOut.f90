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
    use globals_cfg,   only: VERBOSE_, FNAMEMAXLEN_, is1DX_, &
                             isPol_, isPol4_

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
    &NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NCELLS,&
    &NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NMU, NPHI, NANG, NG,&
    &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
    &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
    &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &TREEPTR, GRIDPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,&
    &FLUXES, DIRFLUX, FLUXDIV,&
    &IRAD, RADOUT, NSHOUT, SHTERMS, SOURCE1OUT,&
    &OUTTYPE, OUTPARMS, OUTFILE, SENPARMS)

    !   Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are five types (OUTTYPE) of output: 'R' - for radiance,
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms, and 'M' for
    ! medium properties.
    ! See the overall program documentation for output parameters.

    integer, intent(in) :: NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NCELLS, NSH
    integer, intent(in) :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, NMU, NPHI, NANG, NG
    integer, intent(in) :: BCFLAG, IPFLAG
    integer, intent(in) :: MAXITER, TOTITER, NSHOUT
    integer, intent(in) :: TREEPTR(2, NCELLS), GRIDPTR(8, NCELLS)
    integer*2, intent(in) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical  , intent(in) :: DELTAM
    real, intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ
    real, intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD
    real, intent(in) :: SOLACC, SPLITACC, SHACC
    real, intent(in) :: WAVENO(2), WAVELEN
    real, intent(in) :: XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real, intent(inout) :: OUTPARMS(*)
    real, optional :: SENPARMS(*)
    character(len=1), intent(in) :: SRCTYPE, UNITS, GRIDTYPE, OUTTYPE
    character(len=2), intent(in) :: SFCTYPE
    character(len=*), intent(in) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real, intent(in) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real, intent(in) :: TEMP(NPTS)
    real, intent(in) :: FLUXES(2, NPTS), DIRFLUX(NPTS), FLUXDIV(NPTS)
    real, intent(in) :: SOURCE1OUT(NSTOKES, *)

    real, intent(in) :: RADOUT(NSTOKES, *), SHTERMS(NSHOUT, NPTS)
    integer, intent(inout) :: IRAD

    integer IZ, I, K, N, NX2, NY2, J, JX, JY, ICELL
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

    character(len=FNAMEMAXLEN_) :: OUTFILE_DATA
    character(len=20)  :: DATA_POSTFIX
    integer :: NVAR, NDATA, IDATA

    DATA_POSTFIX = '.sHdOmNG-out'
    OUTFILE_DATA = trim(adjustl(OUTFILE))//trim(adjustl(DATA_POSTFIX))

    GRIDNAME = 'EVEN (X,Y)  '
    if (GRIDTYPE .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (GRIDTYPE .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (SRCTYPE .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (SFCTYPE .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (SFCTYPE .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (SFCTYPE .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (SFCTYPE .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (SFCTYPE .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (SFCTYPE .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (SFCTYPE .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (SFCTYPE .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (SFCTYPE .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (SFCTYPE(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (UNITS .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (UNITS .eq. 'B') then
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
        &'!  L=', ML, '  M=', MM, '  NLM=', NLM,&
        &'   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,&
        &'   NSH=', NSH
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', NSTOKES,&
        &'   NX=', NX, '   NY=', NY, '   NZ=', NZ,&
        &'    NPTS=', NPTS, '   NCELLS=', NCELLS
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(PROPFILE)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',&
        &trim(CKDFILE), '   NUM_G=', NG
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)

    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,&
            &'      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if

    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,&
        &'   INDEPENDENT_PIXEL=', IPFLAG
    write (2, '(A,A22,A,I1)') '!  SURFACE_TYPE=', SFCNAME,&
        &'   HORIZ_BOUNDARY_COND=', BCFLAG

    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(SFCFILE),&
                &'  SKY_RAD=', SKYRAD
        else
            if (SRCTYPE .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,&
                    &'  SKY_RAD=', SKYRAD
            else
                write (2, '(A,F9.7,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD
            end if
        end if

        write (2, '(A,E13.6,A,F10.7,A,F8.3)')&
            &'!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,&
            &'   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0)

    else

        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(SFCFILE)
        else
            write (2, '(A,F8.3,A,F9.7)')&
                &'!  GROUND_TEMP=', GNDTEMP, '  GROUND_EMIS=', 1.0 - GNDALBEDO
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', SKYRAD

    end if

    if (UNITS .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (UNITS .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,&
            &'   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,&
            &'   WAVELENGTH=', WAVELEN
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
        XDOMAIN = XGRID(NX + 1 - ibits(BCFLAG, 0, 1)) - XGRID(1) + 2*OUTPARMS(4)
        STARTX = -OUTPARMS(4)
        if (OUTPARMS(2) .eq. 0.0) then
            NXOUT = 1
        else
            NXOUT = max(1, nint(XDOMAIN/OUTPARMS(2)))
        end if
        !╰────────────────────────────────────────────────╯!

        ! NYOUT
        !╭────────────────────────────────────────────────╮!
        YDOMAIN = YGRID(NY + 1 - ibits(BCFLAG, 1, 1)) - YGRID(1) + 2*OUTPARMS(5)
        STARTY = -OUTPARMS(5)
        if (OUTPARMS(3) .eq. 0.0) then
            NYOUT = 1
        else
            NYOUT = max(1, nint(YDOMAIN/OUTPARMS(3)))
        end if
        !╰────────────────────────────────────────────────╯!

        ! Z0
        !╭────────────────────────────────────────────────╮!
        Z0 = min(max(OUTPARMS(1), ZGRID(1)), ZGRID(NZ))
        !╰────────────────────────────────────────────────╯!

        ! NANGOUT
        !╭────────────────────────────────────────────────╮!
        NANGOUT = nint(OUTPARMS(6))
        !╰────────────────────────────────────────────────╯!

        write (2, '(A,F7.3,3(A,I10))') '!    RADIANCE AT Z=', Z0,&
            &'    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (UNITS .eq. 'T') then
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
            if (present(SENPARMS)) then
                MUOUT = SENPARMS(2*K - 1)
                PHID  = SENPARMS(2*K)
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
        write (2, '(A, 5(I10,A))') '! ', NXOUT, ',', NYOUT, ',', 1, ',', NANGOUT, ',', NSTOKES, ''
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <R> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        N = NANGOUT*NXOUT*NYOUT
        call shd_io__outputDataBin(RADOUT, NSTOKES, N, IRAD, OUTFILE_DATA)
        !╰────────────────────────────────────────────────╯!

    else if (OUTTYPE .eq. 'F') then

        if (UNITS .eq. 'T') then
            C = 1.0 / PI_
        else
            C = 1.0
        end if

        if (SRCTYPE .ne. 'T') then
            NVAR = 3
        else
            NVAR = 2
        end if

        ! Hemispheric flux output, F1
        ! Format 1 is flux at top and bottom of medium at grid point
        ! NVAR=2 or 3, NSET=1, NX=NX, NY=NY, NZ=1
        if (nint(OUTPARMS(1)) .le. 1) then

            write (2, '(A,F7.3)') '!    UPWELLING FLUX:     Z=', ZGRID(NZ)
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', ZGRID(1)

            if (UNITS .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
            write (2, '(A, 5(I10,A))') '! ', NX1, ',', NY1, ',', 1, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NX1*NY1
            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, NZ:NBPTS:NZ), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1, 1:NBPTS:NZ), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(DIRFLUX(1:NBPTS:NZ), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do I = 1, NBPTS, NZ
            !    if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
            !        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
            !        if (SRCTYPE .ne. 'T') then
            !            write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I),&
            !                & C*FLUXES(2, I + NZ - 1), C*FLUXES(1, I), DIRFLUX(I)
            !        else
            !            write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I),&
            !                & C*FLUXES(2, I + NZ - 1), C*FLUXES(1, I)
            !        end if
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 2 is flux at a given level at regular locations
        ! NVAR=2 or 3, NSET=1, NX=NXOUT, NY=NYOUT, NZ=1
        else if (nint(OUTPARMS(1)) .eq. 2) then

            Z0 = min(max(OUTPARMS(2), ZGRID(1)), ZGRID(NZ))
            if (OUTPARMS(3) .eq. 0.0) OUTPARMS(3) = 1.0
            NX2 = NX + 1 - ibits(BCFLAG, 0, 1)
            NXOUT = max(1, nint((XGRID(NX2) - XGRID(1))/OUTPARMS(3)))
            if (OUTPARMS(4) .eq. 0.0) OUTPARMS(4) = 1.0
            NY2 = NY + 1 - ibits(BCFLAG, 1, 1)
            NYOUT = max(1, nint((YGRID(NY2) - YGRID(1))/OUTPARMS(4)))
            write (2, '(A,F7.3,2(A,I4))') '!    UPWELLING FLUX:     Z=', Z0,&
                & '    NXO=', NXOUT, '    NYO=', NYOUT
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', Z0
            if (UNITS .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
            Y0 = YGRID(1)
            do JY = 1, NYOUT
                X0 = XGRID(1)
                do JX = 1, NXOUT
                    IDATA = IDATA + 1
                    call shd_util__findCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
                        & NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
                        & BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)
                    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
                        & X0, Y0, Z0, 2, FLUXES, 1, 1, FDN(IDATA))
                    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
                        & X0, Y0, Z0, 2, FLUXES, 2, 1, FUP(IDATA))
                    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
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
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NXOUT*NYOUT
            IDATA = 0
            call shd_io__outputDataBin(C*FUP, 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FDN, 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(FDIR, 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(FUP)
            deallocate(FDN)
            deallocate(FDIR)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !Y0 = YGRID(1)
            !do JY = 1, NYOUT
            !    X0 = XGRID(1)
            !    do JX = 1, NXOUT
            !        call shd_util__findCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
            !            & NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
            !            & BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)
            !        call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
            !            & X0, Y0, Z0, 2, FLUXES, 1, 1, FDN)
            !        call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
            !            & X0, Y0, Z0, 2, FLUXES, 2, 1, FUP)
            !        call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
            !            & X0, Y0, Z0, 1, DIRFLUX, 1, 1, FDIR)
            !        if (SRCTYPE .ne. 'T') then
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
        ! NVAR=2 or 3, NSET=1, NX=NX1, NY=NY1, NZ=1
        else if (nint(OUTPARMS(1)) .eq. 3) then

            if (UNITS .eq. 'T') then
                FORM = '(1X,F7.3,1X,2(2X,F6.2))'
                write (2, '(A)') '!   Z        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
            allocate(SUM1(NZ))
            allocate(SUM2(NZ))
            allocate(SUM3(NZ))

            do IZ = 1, NZ
                SUM1(IZ) = 0.0
                SUM2(IZ) = 0.0
                SUM3(IZ) = 0.0
                N = 0
                do I = IZ, NBPTS, NZ
                    if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
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
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F3> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            IDATA = 0
            call shd_io__outputDataBin(C*SUM2/N, 1, NZ, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*SUM1/N, 1, NZ, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(SUM3/N, 1, NZ, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(SUM1)
            deallocate(SUM2)
            deallocate(SUM3)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do IZ = 1, NZ
            !    SUM1 = 0.0
            !    SUM2 = 0.0
            !    SUM3 = 0.0
            !    N = 0
            !    do I = IZ, NBPTS, NZ
            !        if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
            !            & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
            !            SUM1 = SUM1 + FLUXES(1, I)
            !            SUM2 = SUM2 + FLUXES(2, I)
            !            SUM3 = SUM3 + DIRFLUX(I)
            !            N = N + 1
            !        end if
            !    end do
            !    if (SRCTYPE .ne. 'T') then
            !        write (2, FORM) ZGRID(IZ), C*SUM2/N, C*SUM1/N, SUM3/N
            !    else
            !        write (2, FORM) ZGRID(IZ), C*SUM2/N, C*SUM1/N
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 4 is fluxes at every base grid point
        ! Format 5 is fluxes at every grid point
        ! NVAR=2 or 3, NSET=1, NX=NX1, NY=NY1, NZ=NZ
        else if (nint(OUTPARMS(1)) .ge. 4) then

            if (nint(OUTPARMS(1)) .eq. 5) then
                N = NPTS
                ALLGRID = .true.
            else
                N = NBPTS
                ALLGRID = .false.
            end if

            if (UNITS .eq. 'T') then
                FORM = '(3(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y       Z        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
                write (2, '(A, 5(I10,A))') '! ', NX1, ',', NY1, ',', NZ, ',', 1, ',', NVAR, ''
            end if
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (ALLGRID) then
                if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F5> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            else
                if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <F4> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            end if

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, :), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(DIRFLUX(:), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if

            ! add grid positions (or binary data cannot be reconstructed)
            if (ALLGRID) then
                IDATA = 0
                call shd_io__outputDataBin(GRIDPOS(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do I = 1, N
            !    if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
            !        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
            !        if (SRCTYPE .ne. 'T') then
            !            write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
            !                & C*FLUXES(2, I), C*FLUXES(1, I), DIRFLUX(I)
            !        else
            !            write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
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
            allocate(SUM1(NZ))

            do IZ = 1, NZ
                SUM1(IZ) = 0.0
                N = 0
                do I = IZ, NBPTS, NZ
                    if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
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
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <H1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            IDATA = 0
            call shd_io__outputDataBin(-SUM1/N, 1, NZ, IDATA, OUTFILE_DATA)

            deallocate(SUM1)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do IZ = 1, NZ
            !    SUM = 0.0
            !    N = 0
            !    do I = IZ, NBPTS, NZ
            !        if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
            !            & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
            !            SUM = SUM + FLUXDIV(I)
            !            N = N + 1
            !        end if
            !    end do
            !    write (2, '(1X,F7.3,2X,E13.5)') ZGRID(IZ), -SUM/N
            !end do
            !╰────────────────────────────────────────────────╯!

        else

            ! Format 2 is flux convergence for every base grid point
            ! Format 3 is flux convergence for every grid point
            if (nint(OUTPARMS(1)) .eq. 3) then
                N = NPTS
                ALLGRID = .true.
            else
                N = NBPTS
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
                write (2, '(A, 5(I10,A))') '! ', NX1, ',', NY1, ',', NZ, ',', 1, ',', NVAR, ''
            end if
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (ALLGRID) then
                if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <H3> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            else
                if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <H2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
            end if

            IDATA = 0
            call shd_io__outputDataBin(-FLUXDIV, 1, N, IDATA, OUTFILE_DATA)

            ! add grid positions (or binary data cannot be reconstructed)
            if (ALLGRID) then
                IDATA = 0
                call shd_io__outputDataBin(GRIDPOS(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do I = 1, N
            !    if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
            !        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
            !        write (2, '(3(1X,F7.3),2X,E13.5)')&
            !            & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I), -FLUXDIV(I)
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        end if

    ! Spherical Harmonic output:
    !   Output mean intensity and net flux (x,y,z components).
    !   and maybe normalized rms of higher order terms
    else if (OUTTYPE .eq. 'S') then

        if (nint(OUTPARMS(1)) .eq. 2) then
            N = NPTS
            ALLGRID = .true.
        else
            N = NBPTS
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
            write (2, '(A, 5(I10,A))') '! ', NX1, ',', NY1, ',', NZ, ',', 1, ',', NVAR, ''
        end if
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (ALLGRID) then
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <S2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        else
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <S1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
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
            call shd_io__outputDataBin(GRIDPOS(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
        end if

        deallocate(SHTERM5)
        !╰────────────────────────────────────────────────╯!

        ! original
        !╭────────────────────────────────────────────────╮!
        !if (NSHOUT .eq. 5) then
        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz       HOrms'
        !    do I = 1, N
        !        if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
        !            &GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
        !            write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
        !                &GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
        !                &SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I),&
        !                &SHTERMS(5, I)/max(1e-20, SHTERMS(1, I))
        !        end if
        !    end do
        !else
        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz '
        !    do I = 1, N
        !        if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
        !            &GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
        !            write (2, '(3(1X,F7.3),1X,4(1X,E12.5))')&
        !                &GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
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
            N = NPTS
            ALLGRID = .true.
        else
            N = NBPTS
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
            write (2, '(A, 5(I10,A))') '! ', NX1, ',', NY1, ',', NZ, ',', 1, ',', NVAR, ''
        end if
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (ALLGRID) then
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <J2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        else
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <J1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        end if

        IDATA = 0
        call shd_io__outputDataBin(EXTINCT, 1, N, IDATA, OUTFILE_DATA)
        IDATA = 0
        call shd_io__outputDataBin(SOURCE1OUT, 4, N, IDATA, OUTFILE_DATA, 'ADD')

        ! add grid positions (or binary data cannot be reconstructed)
        if (ALLGRID) then
            IDATA = 0
            call shd_io__outputDataBin(GRIDPOS(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
        end if
        !╰────────────────────────────────────────────────╯!

        ! original
        !╭────────────────────────────────────────────────╮!
        !do I = 1, N
        !    if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
        !        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
        !        write (2, '(3(1X,F7.3),5(1X,E12.5))')&
        !            & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
        !            & EXTINCT(I), SOURCE1OUT(:, I)
        !    end if
        !end do
        !╰────────────────────────────────────────────────╯!

    ! Medium output: grid point medium properties -
    !   Output extinction, single scattering albedo,
    !   asymmetry parameter, and temperature.
    else if (OUTTYPE .eq. 'M') then

        NVAR = 4

        if (nint(OUTPARMS(1)) .eq. 2) then
            N = NPTS
            ALLGRID = .true.
        else
            N = NBPTS
            ALLGRID = .false.
        end if

        write (2, '(A,A)') '!    X       Y       Z',&
            &'      Extinct  Albedo  Asymmetry Temp'

        ! data processing
        !╭────────────────────────────────────────────────╮!
        allocate(ASYM(N))
        do I = 1, N
            if (NUMPHASE .gt. 0) then
                ASYM(I) = LEGEN(1, 1, IPHASE(I))
            else
                ASYM(I) = LEGEN(1, 1, I)
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
            write (2, '(A, 5(I10,A))') '! ', NX1, ',', NY1, ',', NZ, ',', 1, ',', NVAR, ''
        end if
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (ALLGRID) then
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <M2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        else
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <M1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        end if

        IDATA = 0
        call shd_io__outputDataBin(EXTINCT, 1, N, IDATA, OUTFILE_DATA)
        IDATA = 0
        call shd_io__outputDataBin(ALBEDO, 1, N, IDATA, OUTFILE_DATA, 'ADD')
        IDATA = 0
        call shd_io__outputDataBin(ASYM, 1, N, IDATA, OUTFILE_DATA, 'ADD')
        IDATA = 0
        call shd_io__outputDataBin(TEMP, 1, N, IDATA, OUTFILE_DATA, 'ADD')

        ! add grid positions (or binary data cannot be reconstructed)
        if (ALLGRID) then
            IDATA = 0
            call shd_io__outputDataBin(GRIDPOS(1:3, :), 3, N, IDATA, OUTFILE_DATA, 'ADD')
        end if

        deallocate(ASYM)
        !╰────────────────────────────────────────────────╯!


        ! original
        !╭────────────────────────────────────────────────╮!
        !do I = 1, N
        !    if (NUMPHASE .gt. 0) then
        !        ASYM = LEGEN(1, 1, IPHASE(I))
        !    else
        !        ASYM = LEGEN(1, 1, I)
        !    end if
        !    if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
        !        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
        !        write (2, '(3(1X,F7.3),1X,F10.6,1X,F8.6,1X,F8.5,1X,F7.2)')&
        !            & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
        !            & EXTINCT(I), ALBEDO(I), ASYM, TEMP(I)
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
    &NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
    &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
    &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
    &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    &XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID,&
    &FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT,&
    &OUTTYPE, OUTPARMS, OUTFILE, SENPARMS)
    ! Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are four supported types (OUTTYPE) of output: 'R' - for radi
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms
    ! See the overall program documentation for output parameters.

    integer NSTOKES, NX, NY, NZ, NPTS, NCELLS, NSH
    integer ML, MM, NLM, NMU, NPHI, NANG, NG
    integer BCFLAG, IPFLAG
    integer MAXITER, TOTITER, IRAD, NSHOUT
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real SOLACC, SPLITACC, SHACC
    real WAVENO(2), WAVELEN
    real XDOMAIN, YDOMAIN, XGRID(NX), YGRID(NY), ZGRID(NZ)
    real OUTPARMS(*)
    real, optional :: SENPARMS(*)
    character(len=1), intent(in) :: SRCTYPE, UNITS, GRIDTYPE, OUTTYPE
    character(len=2), intent(in) :: SFCTYPE
    character(len=*), intent(in) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real FLUXES(3, NZ, NY, NX), FLUXDIV(NZ, NY, NX)
    real RADOUT(NSTOKES, *), SHTERMS(NSHOUT, NZ, NY, NX)

    integer I, J, K, N, JX, JY, IZ
    integer NANGOUT, NXOUT, NYOUT
    real :: MUOUT, PHID
    real :: C, SUM
    real(RD_) :: STARTX, STARTY
    real(RD_) :: X0, Y0, Z0
    character(len=128) :: FORM
    character(len=64) :: GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

    real, allocatable :: MU_SET(:), PHI_SET(:)
    character(len=FNAMEMAXLEN_) :: OUTFILE_DATA
    character(len=20)  :: DATA_POSTFIX
    integer :: NVAR, NDATA, IDATA
    real, allocatable :: SUM1(:), SUM2(:), SUM3(:)
    real, allocatable :: SHTERM5(:, :, :)

    DATA_POSTFIX = '.mpi.sHdOmNG-out'
    OUTFILE_DATA = trim(adjustl(OUTFILE))//trim(adjustl(DATA_POSTFIX))

    GRIDNAME = 'EVEN (X,Y)  '
    if (GRIDTYPE .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (GRIDTYPE .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (SRCTYPE .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (SFCTYPE .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (SFCTYPE .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (SFCTYPE .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (SFCTYPE .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (SFCTYPE .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (SFCTYPE .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (SFCTYPE .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (SFCTYPE .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (SFCTYPE .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (SFCTYPE(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (UNITS .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (UNITS .eq. 'B') then
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
        & '!  L=', ML, '  M=', MM, '  NLM=', NLM,&
        & '   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,&
        & '   NSH=', NSH
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', NSTOKES,&
        & '   NX=', NX, '   NY=', NY, '   NZ=', NZ,&
        & '   NPTS=', NPTS, '   NCELLS=', NCELLS
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(PROPFILE)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',&
        & trim(CKDFILE), '   NUM_G=', NG
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)
    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,&
            & '      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if
    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,&
        & '   INDEPENDENT_PIXEL=', IPFLAG
    write (2, '(A,A22,A,I2)') '!  SURFACE_TYPE=', SFCNAME,&
        & '   HORIZ_BOUNDARY_COND=', BCFLAG
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(SFCFILE),&
                & '  SKY_RAD=', SKYRAD
        else
            if (SRCTYPE .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')&
                    & '!  GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,&
                    & '  SKY_RAD=', SKYRAD
            else
                write (2, '(A,F9.7,A,E12.5)')&
                    & '!  GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD
            end if
        end if
        write (2, '(A,E13.6,A,F10.7,A,F8.3)')&
            & '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,&
            & '   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0)
    else
        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(SFCFILE)
        else
            write (2, '(A,F8.3,A,F9.7)')&
                & '!  GROUND_TEMP=', GNDTEMP, '  GROUND_EMIS=', 1.0 - GNDALBEDO
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', SKYRAD
    end if

    if (UNITS .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (UNITS .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,&
            & '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,&
            & '   WAVELENGTH=', WAVELEN
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
            NXOUT = max(1, nint(XDOMAIN/OUTPARMS(2)))
        end if
        STARTY = 0.0
        if (OUTPARMS(3) .le. 0.0) then
            NYOUT = 1
        else
            NYOUT = max(1, nint(YDOMAIN/OUTPARMS(3)))
        end if
        Z0 = min(max(OUTPARMS(1), ZGRID(1)), ZGRID(NZ))
        NANGOUT = nint(OUTPARMS(6))

        write (2, '(A,F7.3,3(A,I10))') '!    RADIANCE AT Z=', Z0,&
            & '    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (UNITS .eq. 'T') then
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
            if (present(SENPARMS)) then
                MUOUT = SENPARMS(2*K - 1)
                PHID  = SENPARMS(2*K)
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

        if (IPFLAG .eq. 0) then
            call shd_io__radShiftCorr(RADOUT(:, 1:NXOUT*NYOUT*NANGOUT), NSTOKES, NANGOUT, NXOUT, NYOUT,&
                OUTPARMS(2), OUTPARMS(3), XDOMAIN, YDOMAIN, Z0-ZGRID(1), MU_SET, PHI_SET)
        end if

        deallocate(MU_SET, PHI_SET)
        !╰────────────────────────────────────────────────╯!

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (2, '(A)') '! The following contains information for binary data intepretation:'
        write (2, '(A)') '! R'
        write (2, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
        write (2, '(A, 5(I10,A))') '! ', NXOUT, ',', NYOUT, ',', 1, ',', NANGOUT, ',', NSTOKES, ''
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <R> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        N = NANGOUT*NXOUT*NYOUT
        call shd_io__outputDataBin(RADOUT, NSTOKES, N, IRAD, OUTFILE_DATA)
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

        if (UNITS .eq. 'T') then
            C = 1.0 / PI_
        else
            C = 1.0
        end if

        if (SRCTYPE .ne. 'T') then
            NVAR = 3
        else
            NVAR = 2
        end if

        if (nint(OUTPARMS(1)) .le. 1) then

            write (2, '(A,F7.3)') '!    UPWELLING FLUX:     Z=', ZGRID(NZ)
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', ZGRID(1)
            if (UNITS .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
            write (2, '(A, 5(I10,A))') '! ', NX, ',', NY, ',', 1, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <F1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NX*NY
            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, NZ, :, :), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1,  1, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(FLUXES(3,  1, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !C = 1.0
            !if (UNITS .eq. 'T') C = 1.0 / PI_
            !do I = 1, NX
            !    do J = 1, NY
            !        if (SRCTYPE .ne. 'T') then
            !            write (2, FORM) XGRID(I), YGRID(J),&
            !                & C*FLUXES(2, NZ, J, I), C*FLUXES(1, 1, J, I), FLUXES(3, 1, J, I)
            !        else
            !            write (2, FORM) XGRID(I), YGRID(J),&
            !                & C*FLUXES(2, NZ, J, I), C*FLUXES(1, 1, J, I)
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

            if (UNITS .eq. 'T') then
                FORM = '(1X,F7.3,1X,2(2X,F6.2))'
                write (2, '(A)') '!   Z        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
            allocate(SUM1(NZ))
            allocate(SUM2(NZ))
            allocate(SUM3(NZ))

            do IZ = 1, NZ
                SUM1(IZ) = 0.0
                SUM2(IZ) = 0.0
                SUM3(IZ) = 0.0

                do I = 1, NX
                    do J = 1, NY
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
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <F3> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NX*NY
            IDATA = 0
            call shd_io__outputDataBin(C*SUM2/N, 1, NZ, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*SUM1/N, 1, NZ, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(SUM3/N, 1, NZ, IDATA, OUTFILE_DATA, 'ADD')
            end if

            deallocate(SUM1)
            deallocate(SUM2)
            deallocate(SUM3)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !C = 1.0
            !if (UNITS .eq. 'T') C = 1.0 / PI_
            !N = NX*NY
            !do K = 1, NZ
            !    SUM1 = 0.0
            !    SUM2 = 0.0
            !    SUM3 = 0.0
            !    do I = 1, NX
            !        do J = 1, NY
            !            SUM1 = SUM1 + FLUXES(1, K, J, I)
            !            SUM2 = SUM2 + FLUXES(2, K, J, I)
            !            SUM3 = SUM3 + FLUXES(3, K, J, I)
            !        end do
            !    end do
            !    if (SRCTYPE .ne. 'T') then
            !        write (2, FORM) ZGRID(K), C*SUM2/N, C*SUM1/N, SUM3/N
            !    else
            !        write (2, FORM) ZGRID(K), C*SUM2/N, C*SUM1/N
            !    end if
            !end do
            !╰────────────────────────────────────────────────╯!

        ! Format 4 is fluxes at every base grid point
        else if (nint(OUTPARMS(1)) .eq. 4) then

            if (UNITS .eq. 'T') then
                FORM = '(3(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y       Z        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
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
            write (2, '(A, 5(I10,A))') '! ', NX, ',', NY, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <F4> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NX*NY*NZ
            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(2, :, :, :), 1, N, IDATA, OUTFILE_DATA)

            IDATA = 0
            call shd_io__outputDataBin(C*FLUXES(1, :, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')

            if (SRCTYPE .ne. 'T') then
                IDATA = 0
                call shd_io__outputDataBin(FLUXES(3, :, :, :), 1, N, IDATA, OUTFILE_DATA, 'ADD')
            end if
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !C = 1.0
            !if (UNITS .eq. 'T') C = 1.0 / PI_
            !do I = 1, NX
            !    do J = 1, NY
            !        do K = 1, NZ
            !            if (SRCTYPE .ne. 'T') then
            !                write (2, FORM) XGRID(I), YGRID(J), ZGRID(K),&
            !                    & C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I), FLUXES(3, K, J, I)
            !            else
            !                write (2, FORM) XGRID(I), YGRID(J), ZGRID(K),&
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
            allocate(SUM1(NZ))

            do IZ = 1, NZ
                SUM1(IZ) = 0.0
                do I = 1, NX
                    do J = 1, NY
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
            write (2, '(A, 5(I10,A))') '! ', 1, ',', 1, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <H1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NX*NY
            IDATA = 0
            call shd_io__outputDataBin(-SUM1/N, 1, NZ, IDATA, OUTFILE_DATA)

            deallocate(SUM1)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !do K = 1, NZ
            !    SUM = 0.0
            !    do I = 1, NX
            !        do J = 1, NY
            !            SUM = SUM + FLUXDIV(K, J, I)
            !        end do
            !    end do
            !    write (2, '(1X,F7.3,2X,E13.5)') ZGRID(K), -SUM/(NX*NY)
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
            write (2, '(A, 5(I10,A))') '! ', NX, ',', NY, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <H2> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

            N = NX*NY*NZ
            IDATA = 0
            call shd_io__outputDataBin(-FLUXDIV, 1, N, IDATA, OUTFILE_DATA)
            !╰────────────────────────────────────────────────╯!

            ! original
            !╭────────────────────────────────────────────────╮!
            !    do I = 1, NX
            !        do J = 1, NY
            !            do K = 1, NZ
            !                write (2, '(3(1X,F7.3),2X,E13.5)')&
            !                    & XGRID(I), YGRID(J), ZGRID(K), -FLUXDIV(K, J, I)
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

            N = NX*NY*NZ

            ! data processing
            !╭────────────────────────────────────────────────╮!
            allocate(SHTERM5(NZ, NY, NX))

            if (NSHOUT .eq. 5) then

                NVAR = 5
                write (2, '(A,A)') '!    X       Y        Z         Imean',&
                    &'          Fx           Fy           Fz       HOrms'
                do IZ = 1, NZ
                    do I = 1, NX
                        do J = 1, NY
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
            write (2, '(A, 5(I10,A))') '! ', NX, ',', NY, ',', NZ, ',', 1, ',', NVAR, ''
            !╰────────────────────────────────────────────────╯!

            ! binary data
            !╭────────────────────────────────────────────────╮!
            if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResultMpi]: <S1> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'

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
        !    do I = 1, NX
        !        do J = 1, NY
        !            do K = 1, NZ
        !                write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
        !                    &XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1, K, J, I),&
        !                    &SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I),&
        !                    &SHTERMS(5, K, J, I)/max(1e-20, SHTERMS(1, K, J, I))
        !            end do
        !        end do
        !    end do

        !else

        !    write (2, '(A,A)') '!    X       Y        Z         Imean',&
        !        &'          Fx           Fy           Fz '
        !    do I = 1, NX
        !        do J = 1, NY
        !            do K = 1, NZ
        !                write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
        !                    &XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1, K, J, I),&
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

subroutine shd_io__outputCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    & NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
    & PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
    & BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
    & SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    & SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
    & SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    & IVIS, VISOUT, OUTPARMS, OUTFILE, SENPARMS)

    ! Writes the visualization image in PGM binary format with either
    ! one byte or two byte pixels.  The PGM binary format (P5) has an
    ! ascii header which has the image parameters amd can include comments
    ! (here the SHDOM input parameters), followed by the byte or two byte
    ! image pixels.

    integer NSTOKES, NX, NY, NZ, NPTS, NCELLS, NSH
    integer ML, MM, NLM, NMU, NPHI, NANG, NG
    integer BCFLAG, IPFLAG
    integer MAXITER, TOTITER, IVIS
    logical, intent(in) :: DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real SOLACC, SPLITACC, SHACC
    real WAVENO(2), WAVELEN
    real OUTPARMS(*)
    real, optional :: SENPARMS(*)
    character(len=1), intent(in) :: SRCTYPE, UNITS, GRIDTYPE
    character(len=2), intent(in) :: SFCTYPE
    character(len=*), intent(in) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real VISOUT(NSTOKES, *)

    character(len=128) :: GRIDNAME, SOURCENAME, UNITSNAME, SFCNAME, OUTNAME
    character(len=FNAMEMAXLEN_) :: OUTFILE_DATA
    character(len=20)  :: DATA_POSTFIX
    integer :: N, IDATA, NL, NS, CAMID

    DATA_POSTFIX = '.sHdOmNG-out'
    OUTFILE_DATA = trim(adjustl(OUTFILE))//trim(adjustl(DATA_POSTFIX))

    GRIDNAME = 'EVEN (X,Y)  '
    if (GRIDTYPE .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (GRIDTYPE .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (SRCTYPE .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (SFCTYPE .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (SFCTYPE .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (SFCTYPE .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (SFCTYPE .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (SFCTYPE .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (SFCTYPE .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (SFCTYPE .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (SFCTYPE .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (SFCTYPE .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (SFCTYPE(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (UNITS .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (UNITS .eq. 'B') then
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
        &'!  L=', ML, '  M=', MM, '  NLM=', NLM,&
        &'   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,&
        &'   NSH=', NSH
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', NSTOKES,&
        &'   NX=', NX, '   NY=', NY, '   NZ=', NZ,&
        &'    NPTS=', NPTS, '   NCELLS=', NCELLS
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(PROPFILE)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',&
        &trim(CKDFILE), '   NUM_G=', NG
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)

    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,&
            &'      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if

    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,&
        &'   INDEPENDENT_PIXEL=', IPFLAG
    write (2, '(A,A22,A,I1)') '!  SURFACE_TYPE=', SFCNAME,&
        &'   HORIZ_BOUNDARY_COND=', BCFLAG

    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(SFCFILE),&
                &'  SKY_RAD=', SKYRAD
        else
            if (SRCTYPE .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,&
                    &'  SKY_RAD=', SKYRAD
            else
                write (2, '(A,F9.7,A,E12.5)')&
                    &'!  GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD
            end if
        end if

        write (2, '(A,E13.6,A,F10.7,A,F8.3)')&
            &'!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,&
            &'   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0)

    else

        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(SFCFILE)
        else
            write (2, '(A,F8.3,A,F9.7)')&
                &'!  GROUND_TEMP=', GNDTEMP, '  GROUND_EMIS=', 1.0 - GNDALBEDO
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', SKYRAD

    end if

    if (UNITS .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (UNITS .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,&
            &'   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,&
            &'   WAVELENGTH=', WAVELEN
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
        write (2, '(A, 5(I10,A))') '! ', NS, ',', NL, ',', 1, ',', 1, ',', NSTOKES, ''
        !╰────────────────────────────────────────────────╯!

        ! binary data
        !╭────────────────────────────────────────────────╮!
        if (VERBOSE_ .ge. 5) write(*, '(A)') 'Message [_outputResult]: <V> Data info saved at <'//trim(adjustl(OUTFILE))//'>.'
        N = NL*NS
        IDATA = 0
        call shd_io__outputDataBin(VISOUT, NSTOKES, N, IDATA, OUTFILE_DATA)
        !╰────────────────────────────────────────────────╯!

    end if

    close (2)
    return
end subroutine shd_io__outputCamRad

subroutine shd_io__outputImg(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    & NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
    & PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
    & BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
    & SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    & SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
    & SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    & IVIS, VISOUT, OUTPARMS, OUTFILE)

    ! Writes the visualization image in PGM binary format with either
    ! one byte or two byte pixels.  The PGM binary format (P5) has an
    ! ascii header which has the image parameters amd can include comments
    ! (here the SHDOM input parameters), followed by the byte or two byte
    ! image pixels.

    integer NSTOKES, NX, NY, NZ, NPTS, NCELLS, NSH
    integer ML, MM, NLM, NMU, NPHI, NANG, NG
    integer BCFLAG, IPFLAG
    integer MAXITER, TOTITER, IVIS
    logical, intent(in) :: DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real SOLACC, SPLITACC, SHACC
    real WAVENO(2), WAVELEN
    real OUTPARMS(*)
    character(len=1), intent(in) :: SRCTYPE, UNITS, GRIDTYPE
    character(len=2), intent(in) :: SFCTYPE
    character(len=*), intent(in) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=*), intent(in) :: OUTFILE
    real VISOUT(NSTOKES, *)

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
    allocate (BUFFER1(NBYTES*NS))
    NLOUT = NSTOKES*NL

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
    if (GRIDTYPE .eq. 'P') then
        GRIDNAME = GRIDNAME(1:14)//'PROPERTY-FILE (Z) '
    else if (GRIDTYPE .eq. 'F') then
        GRIDNAME = GRIDNAME(1:14)//'INPUT-FILE (Z)    '
    else
        GRIDNAME = GRIDNAME(1:14)//'EVEN (Z)          '
    end if

    if (SRCTYPE .eq. 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE .eq. 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE .eq. 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (SFCTYPE .eq. 'FL') then
        SFCNAME = 'FIXED LAMBERTIAN'
    else if (SFCTYPE .eq. 'VL') then
        SFCNAME = 'VARIABLE LAMBERTIAN'
    else if (SFCTYPE .eq. 'VW') then
        SFCNAME = 'VARIABLE WAVE-FRESNEL'
    else if (SFCTYPE .eq. 'VD') then
        SFCNAME = 'VARIABLE DINER'
    else if (SFCTYPE .eq. 'VO') then
        SFCNAME = 'VARIABLE OCEAN'
    else if (SFCTYPE .eq. 'VR') then
        SFCNAME = 'VARIABLE RPV'
    else if (SFCTYPE .eq. 'VT') then
        SFCNAME = 'VARIABLE LSRT'
    else if (SFCTYPE .eq. 'VJ') then
        SFCNAME = 'VARIABLE LSRT-J'
    else if (SFCTYPE .eq. 'VX') then
        SFCNAME = 'VARIABLE MIXED'
    else if (SFCTYPE(1:1) .eq. 'V') then
        SFCNAME = 'VARIABLE OTHER'
    end if

    if (UNITS .eq. 'T') then
        UNITSNAME = 'KELVIN'
    else if (UNITS .eq. 'B') then
        UNITSNAME = 'WATTS/(M^2 STER)'
    else
        UNITSNAME = 'WATTS/(M^2 MICRON STER)'
    end if

    write (BUF(I:I + 79), '(A,A1)')&
        & '# Polarized SHDOM Visualization Output', LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10,A1)')&
        & '# L=', ML, '  M=', MM, '  NLM=', NLM,&
        & '   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,&
        & '   NSH=', NSH, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 90), '(A,I1,3(A,I4),2(A,I8),A1)')&
        & '# NSTOKES=', NSTOKES, '  NX=', NX, '   NY=', NY, '  NZ=', NZ,&
        & '    NPTS=', NPTS, '   NCELLS=', NCELLS, LF
    I = I + index(BUF(I:I + 79), LF)
    !write (BUF(I:I + 79), '(A,A,A1)')&
    !    & '# PROPERTY_FILE=', trim(PROPFILE), LF
    I = I + index(BUF(I:I + 79), LF)
    !write (BUF(I:I + 79), '(A,A,A,I2,A1)')&
    !    & '# CORRELATED_K-DIST_FILE=', trim(CKDFILE), '   NUM_G=', NG, LF
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
        & '   INDEPENDENT_PIXEL=', IPFLAG, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A22,A,I2,A1)')&
        & '# SURFACE_TYPE=', SFCNAME, '   HORIZ_BOUNDARY_COND=', BCFLAG, LF
    I = I + index(BUF(I:I + 79), LF)
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        if (SFCTYPE(1:1) .eq. 'V') then
            !write (BUF(I:I + 79), '(A,A,A,E12.5,A1)')&
            !    & '# SURFACE_FILE=', trim(SFCFILE),&
            !    & '  SKY_RAD=', SKYRAD, LF
            I = I + index(BUF(I:I + 79), LF)
        else
            if (SRCTYPE .eq. 'B') then
                write (BUF(I:I + 79), '(A,F9.7,A,F8.3,A,E12.5,A1)')&
                    &'# GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,&
                    &'  SKY_RAD=', SKYRAD, LF
                I = I + index(BUF(I:I + 79), LF)
            else
                write (BUF(I:I + 79), '(A,F9.7,A,E12.5,A1)')&
                    &'# GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD, LF
                I = I + index(BUF(I:I + 79), LF)
            end if
        end if
        write (BUF(I:I + 79), '(A,E13.6,A,F10.7,A,F8.3,A1)')&
            &'# SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,&
            &'   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0), LF
        I = I + index(BUF(I:I + 79), LF)
    else
        if (SFCTYPE(1:1) .eq. 'V') then
            !write (BUF(I:I + 79), '(A,A,A1)')&
            !    &'# SURFACE_FILE=', trim(SFCFILE), LF
            I = I + index(BUF(I:I + 79), LF)
        else
            write (BUF(I:I + 79), '(A,F8.3,A,F9.7,A1)')&
                &'# GROUND_TEMP=', GNDTEMP,&
                &'  GROUND_EMIS=', 1.0 - GNDALBEDO, LF
            I = I + index(BUF(I:I + 79), LF)
        end if
        write (BUF(I:I + 79), '(A,F8.3,A1)') '# SKY_TEMP=', SKYRAD, LF
        I = I + index(BUF(I:I + 79), LF)
    end if

    if (UNITS .eq. 'T') then
        write (BUF(I:I + 79), '(A,A24,A1)') '# UNITS=', UNITSNAME, LF
        I = I + index(BUF(I:I + 79), LF)
    else if (UNITS .eq. 'B') then
        write (BUF(I:I + 79), '(A,A24,A,2(1X,F10.2),A1)')&
            &'# UNITS=', UNITSNAME,&
            &'   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2), LF
        I = I + index(BUF(I:I + 79), LF)
    else
        write (BUF(I:I + 79), '(A,A24,A,F10.2,A1)')&
            &'# UNITS=', UNITSNAME, '   WAVELENGTH=', WAVELEN, LF
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
    do M = 1, NSTOKES
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
    do M = 1, NSTOKES
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

    deallocate (BUFFER1)
    close (2)
    return
end subroutine shd_io__outputImg

subroutine shd_io__outputCellSplit(&
    &OUTFILE, GRIDPTR, GRIDPOS, EXTINCT,&
    &SHPTR, SOURCE, NCELLS, NSTOKES)

    ! Outputs the final cell splitting criterion for each grid cell
    ! and direction (X,Y,Z).  Calls CELL_SPLIT_TEST.
    ! Mainly for debugging purposes.

    integer NCELLS, NSTOKES
    integer SHPTR(*), GRIDPTR(8, *)
    real GRIDPOS(3, *), EXTINCT(*), SOURCE(NSTOKES, *)
    character(len=*), intent(in) :: OUTFILE
    integer IP1, IP2, ICELL, IDIR
    real ADAPTCRIT(3), MAXADAPT, XC, YC, ZC, EXT

    open (UNIT=9, FILE=trim(adjustl(OUTFILE)), STATUS='UNKNOWN')
    write (9, *) '!  ICELL   Xc     Yc     Zc  ',&
        &'Extinct SplitX  SplitY  SplitZ'

    ! Loop over all the cells
    do ICELL = 1, NCELLS
        ! Get the adaptive grid criterion for the 3 directions
        call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
            &SHPTR, NSTOKES, SOURCE, ICELL,&
            &ADAPTCRIT, MAXADAPT, IDIR)
        ! Compute the cell center location and mean extinction
        IP1 = GRIDPTR(1, ICELL)
        IP2 = GRIDPTR(8, ICELL)
        XC = 0.5*(GRIDPOS(1, IP1) + GRIDPOS(1, IP2))
        YC = 0.5*(GRIDPOS(2, IP1) + GRIDPOS(2, IP2))
        ZC = 0.5*(GRIDPOS(3, IP1) + GRIDPOS(3, IP2))
        EXT =  (EXTINCT(GRIDPTR(1, ICELL)) + EXTINCT(GRIDPTR(2, ICELL))&
            &+ EXTINCT(GRIDPTR(3, ICELL)) + EXTINCT(GRIDPTR(4, ICELL))&
            &+ EXTINCT(GRIDPTR(5, ICELL)) + EXTINCT(GRIDPTR(6, ICELL))&
            &+ EXTINCT(GRIDPTR(7, ICELL)) + EXTINCT(GRIDPTR(8, ICELL)))/8
        write (9, '(1X,I5,3(1X,F6.3),1X,E11.4,3(1X,F7.5))')&
            &ICELL, XC, YC, ZC, EXT,&
            &ADAPTCRIT(1), ADAPTCRIT(2), ADAPTCRIT(3)
    end do

    close (9)
    return
end subroutine shd_io__outputCellSplit

subroutine shd_io__visualizeCell(&
    & OUTFILE, YVAL, OUTTYPE,&
    & IPFLAG, NX, XDOMAIN, ZDOMAIN,&
    & NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS)

    ! Visualizes the cell structure by outputting GLE code to draw
    ! a diagram of the cells in the XZ plane.  Only the end cells are
    ! output.  The outlines of the cells are drawn, and arrows between
    ! the cell centers show the neighbor relationships.
    ! OUTTYPE=0 for cell outlines only, OUTTYPE=1 for cells and arrows.
    ! There is a special output for IP mode.

    integer OUTTYPE, IPFLAG, NX, NCELLS
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real YVAL, XDOMAIN, ZDOMAIN, GRIDPOS(3, *)
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
        DELX = 0.5*XDOMAIN/NX
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

    do IC = 1, NCELLS
        ! Output if it is an end cell and crosses the Y value
        if (TREEPTR(2, IC) .eq. 0 .and.&
            & GRIDPOS(2, GRIDPTR(1, IC)) .le. YVAL .and.&
            & GRIDPOS(2, GRIDPTR(3, IC)) .ge. YVAL) then
            ! Get the center of the cell
            XC = (GRIDPOS(1, GRIDPTR(1, IC)) + GRIDPOS(1, GRIDPTR(8, IC)))/2
            ZC = (GRIDPOS(3, GRIDPTR(1, IC)) + GRIDPOS(3, GRIDPTR(8, IC)))/2
            ! Draw the outline in XZ (four lines and four grid points)
            !   and the cell number at the center using the gle subrouti
            X1 = GRIDPOS(1, GRIDPTR(1, IC))
            Z1 = GRIDPOS(3, GRIDPTR(1, IC))
            X2 = GRIDPOS(1, GRIDPTR(2, IC))
            Z2 = GRIDPOS(3, GRIDPTR(2, IC))
            X3 = GRIDPOS(1, GRIDPTR(6, IC))
            Z3 = GRIDPOS(3, GRIDPTR(6, IC))
            X4 = GRIDPOS(1, GRIDPTR(5, IC))
            Z4 = GRIDPOS(3, GRIDPTR(5, IC))
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
                    if (NEIGHPTR(J, IC) .ne. 0 .and. J .ne. 3 .and. J .ne. 4) then
                        INeigh = abs(NEIGHPTR(J, IC))
                        XN = (GRIDPOS(1, GRIDPTR(1, INeigh)) + GRIDPOS(1, GRIDPTR(8, INeigh)))/2
                        ZN = (GRIDPOS(3, GRIDPTR(1, INeigh)) + GRIDPOS(3, GRIDPTR(8, INeigh)))/2
                        XA = XC + 1.0*(XN - XC)
                        ZA = ZC + 1.0*(ZN - ZC)
                        if (.not. (XN - XC .gt. 0.0 .and. J .eq. 1 .or.&
                            & XN - XC .lt. 0.0 .and. J .eq. 2)) then
                            write (9, '(A,4(1X,F6.3),1X,I8)') '@drawarrow ',&
                                & S*XC, S*ZC, S*XA, S*ZA, NEIGHPTR(J, IC)
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
    &OUTSAVEFILE, NSTOKES, NX, NY, NZ,&
    &ML, MM, NLM, WORK, XGRID, YGRID, ZGRID,&
    &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
    &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

    ! Write a binary save file containing the current state of
    ! the SHDOM solution. This includes the cell tree structure,
    ! the spherical harmonic expansion of the source function,
    ! the low order terms of the radiance expansion, and the flux array.
    ! If the file name is 'NONE' then no file is written.

    character(len=*) , intent(in) :: OUTSAVEFILE
    integer(kind=I__), intent(in) :: NSTOKES, NX, NY, NZ, ML, MM, NLM, NPTS, NCELLS
    integer(kind=I__), intent(in) :: SHPTR(NPTS+1), RSHPTR(NPTS)
    integer(kind=I__), intent(in) :: GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer(kind=2)  , intent(in) :: CELLFLAGS(NCELLS)
    real(kind=R__)   , intent(in) :: XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
    real(kind=R__)   , intent(inout) :: WORK(*)
    real(kind=R__)   , intent(in) :: GRIDPOS(3, NPTS), FLUXES(2, NPTS), SOURCE(*), RADIANCE(*)

    integer(kind=4) :: I, J, K, IR, NS, NR
    integer(kind=4) :: DATA_BYTE
    real(kind=4) :: DATA_SIZE

    if (OUTSAVEFILE(1:4) .eq. 'NONE') return

    open (unit=91, file=trim(adjustl(OUTSAVEFILE)), status='REPLACE',&
        &form='UNFORMATTED')

    write (91) NSTOKES, NX, NY, NZ, ML, MM, NLM, NPTS, NCELLS

    write (91) (XGRID(I), I=1, NX+1)
    write (91) (YGRID(I), I=1, NY+1)
    write (91) (ZGRID(I), I=1, NZ)

    write (91) ((GRIDPTR(J, I) , J=1, 8), I=1, NCELLS)
    write (91) ((NEIGHPTR(J, I), J=1, 6), I=1, NCELLS)
    write (91) ((TREEPTR(J, I) , J=1, 2), I=1, NCELLS)
    write (91) (CELLFLAGS(I)   , I=1, NCELLS)

    write (91) ((GRIDPOS(J, I), J=1, 3), I=1, NPTS)
    write (91) ((FLUXES(J, I) , J=1, 2), I=1, NPTS)

    write (91) (SHPTR(I), I=1, NPTS + 1)

    NS = SHPTR(NPTS + 1)
    write (91) (SOURCE(I), I=1, NSTOKES*NS)

    K = 1
    do I = 1, NPTS
        IR = RSHPTR(I)
        do J = 1, NSTOKES*4

            ! the following is an attempt to fix the segmentation issue
            !╭────────────────────────────────────────────────╮!
            ! RSHPTR(1) = 0 (from shd_init__makeRad)
            ! originally it was:
            ! WORK(K) = RADIANCE(NSTOKES*(IR-1)+J)
            ! when I=1, IR=0, NSTOKES*(IR-1)+J is 0 when J=1 (NSTOKES=1)
            WORK(K) = RADIANCE(NSTOKES*(IR)+J)
            !╰────────────────────────────────────────────────╯!

            K = K + 1
        end do
    end do

    NR = 4*NPTS
    write (91) (WORK(I), I=1, NSTOKES*NR)

    close (unit=91)

    ! display file size
    !╭────────────────────────────────────────────────╮!
    inquire(FILE=trim(adjustl(OUTSAVEFILE)), size=DATA_BYTE)
    DATA_SIZE = DATA_BYTE/1000000.0
    if (VERBOSE_ .ge. 5) write(*, '(A,F6.2,A)') 'Message [_saveState]: Data saved at <'//trim(adjustl(OUTSAVEFILE))//'> (',&
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
            if (VERBOSE_ .ge. 5) write(*, '(A,F0.2,A)') 'Message [_outputDataBin]: Data added to <'//trim(adjustl(OUTFILE_DATA))//'> (',&
                &DATA_SIZE, 'MB).'
        else
            if (VERBOSE_ .ge. 5) write(*, '(A,F0.2,A)') 'Message [_outputDataBin]: Data saved at <'//trim(adjustl(OUTFILE_DATA))//'> (',&
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
