! Input/Output for SHDOM
!╭──────────────────────────────────────────────────────────────╮!
!╰──────────────────────────────────────────────────────────────╯!

module shdUtil_io

  implicit none
  private

  ! Output types
  !╭────────────────────────────────────────────────╮!
  ! R : radiance, Nang, Nx, Ny, Nz=1
  ! V1: camera mode
  ! V2: cross track
  ! F1: irradiance,
  !╰────────────────────────────────────────────────╯!

  ! public modules
  !╭────────────────────────────────────────────────╮!
  public :: shdUtil__outputData
  public :: shdUtil__outputDesc
  !╰────────────────────────────────────────────────╯!

contains

subroutine shdUtil__outputDataRad()
    integer :: i, j, k

    ! R Nx, Ny, Nz=1, NANGOUT

    ! Radiance output
    XDOMAIN = XGRID(NX + 1 - ibits(BCFLAG, 0, 1)) - XGRID(1) + 2*OUTPARMS(4)
    STARTX = -OUTPARMS(4)
    if (OUTPARMS(2) .eq. 0.0) then
        NXOUT = 1
    else
        NXOUT = max(1, nint(XDOMAIN/OUTPARMS(2)))
    end if
    YDOMAIN = YGRID(NY + 1 - ibits(BCFLAG, 1, 1)) - YGRID(1) + 2*OUTPARMS(5)
    STARTY = -OUTPARMS(5)
    if (OUTPARMS(3) .eq. 0.0) then
        NYOUT = 1
    else
        NYOUT = max(1, nint(YDOMAIN/OUTPARMS(3)))
    end if
    Z0 = min(max(OUTPARMS(1), ZGRID(1)), ZGRID(NZ))
    NANGOUT = nint(OUTPARMS(6))

    write (2, '(A,F7.3,3(A,I4))') '!    RADIANCE AT Z=', Z0,&
        & '    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
    if (UNITS .eq. 'T') then
        FORM = '(2(1X,F7.3),1X,4(1X,F6.2))'
        write (2, '(A)')&
            & '!   X       Y       RADIANCE (STOKES I, Q, U, V)'
    else
        FORM = '(2(1X,F7.3),1X,4(1X,E12.5))'
        write (2, '(A)')&
            & '!   X       Y        RADIANCE (STOKES I, Q, U, V)'
    end if

    do K = 1, NANGOUT
        MUOUT = OUTPARMS(2*K + 5)
        PHID = OUTPARMS(2*K + 6)
        write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')&
            & '! ', MUOUT, PHID, '<- (mu,phi)'
        Y0 = STARTY
        do JY = 1, NYOUT
            X0 = STARTX
            do JX = 1, NXOUT
                IRAD = IRAD + 1
                write (2, FORM) X0, Y0, RADOUT(:, IRAD)
                X0 = X0 + OUTPARMS(2)
            end do
            Y0 = Y0 + OUTPARMS(3)
        end do
    end do
    return

end subroutine shdUtil__outputDataRad()

subroutine shdUtil__outputDataFlx()
    integer :: i, j, k

    ! Hemispheric flux output
    !   Format 1 is flux at top and bottom of medium at grid point
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
        C = 1.0
        if (UNITS .eq. 'T') C = 1.0/PI
        do I = 1, NBPTS, NZ
            if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
                if (SRCTYPE .ne. 'T') then
                    write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I),&
                        & C*FLUXES(2, I + NZ - 1), C*FLUXES(1, I), DIRFLUX(I)
                else
                    write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I),&
                        & C*FLUXES(2, I + NZ - 1), C*FLUXES(1, I)
                end if
            end if
        end do

    else if (nint(OUTPARMS(1)) .eq. 2) then
        ! Format 2 is flux at a given level at regular locations
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
        C = 1.0
        if (UNITS .eq. 'T') C = 1.0/PI
        Y0 = YGRID(1)
        do JY = 1, NYOUT
            X0 = XGRID(1)
            do JX = 1, NXOUT
                call LOCATE_GRID_CELL(NX, NY, NZ, XGRID, YGRID, ZGRID,&
                    & NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
                    & BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)
                call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
                    & X0, Y0, Z0, 2, FLUXES, 1, 1, FDN)
                call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
                    & X0, Y0, Z0, 2, FLUXES, 2, 1, FUP)
                call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,&
                    & X0, Y0, Z0, 1, DIRFLUX, 1, 1, FDIR)
                if (SRCTYPE .ne. 'T') then
                    write (2, FORM) X0, Y0, C*FUP, C*FDN, FDIR
                else
                    write (2, FORM) X0, Y0, C*FUP, C*FDN
                end if
                X0 = X0 + OUTPARMS(3)
            end do
            Y0 = Y0 + OUTPARMS(4)
        end do

    else if (nint(OUTPARMS(1)) .eq. 3) then
        ! Format 3 is domain averaged vertical profile
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
        C = 1.0
        if (UNITS .eq. 'T') C = 1.0/PI
        do IZ = 1, NZ
            SUM1 = 0.0
            SUM2 = 0.0
            SUM3 = 0.0
            N = 0
            do I = IZ, NBPTS, NZ
                if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                    & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
                    SUM1 = SUM1 + FLUXES(1, I)
                    SUM2 = SUM2 + FLUXES(2, I)
                    SUM3 = SUM3 + DIRFLUX(I)
                    N = N + 1
                end if
            end do
            if (SRCTYPE .ne. 'T') then
                write (2, FORM) ZGRID(IZ), C*SUM2/N, C*SUM1/N, SUM3/N
            else
                write (2, FORM) ZGRID(IZ), C*SUM2/N, C*SUM1/N
            end if
        end do

    else if (nint(OUTPARMS(1)) .ge. 4) then
        ! Format 4 is fluxes at every base grid point
        ! Format 5 is fluxes at every grid point
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
        C = 1.0
        if (UNITS .eq. 'T') C = 1.0/PI
        do I = 1, N
            if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                if (SRCTYPE .ne. 'T') then
                    write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
                        & C*FLUXES(2, I), C*FLUXES(1, I), DIRFLUX(I)
                else
                    write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
                        & C*FLUXES(2, I), C*FLUXES(1, I)
                end if
            end if
        end do
    end if

    return
end subroutine shdUtil__outputDataFlx()

subroutine shdUtil__outputData(&
    & NSTOKES, NX, NY, NZ, NBPTS, NPTS, NCELLS,&
    & NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    & NMU, NPHI, NANG, NG,&
    & PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
    & BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
    & SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    & SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
    & SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    & XGRID, YGRID, ZGRID, GRIDPOS,&
    & TREEPTR, GRIDPTR, CELLFLAGS,&
    & EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,&
    & FLUXES, DIRFLUX, FLUXDIV,&
    & IRAD, RADOUT, NSHOUT, SHTERMS, SOURCE1OUT,&
    & OUTTYPE, OUTPARMS, OUTFILE)

    ! Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are five types (OUTTYPE) of output: 'R' - for radiance,
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms, and 'M' for
    ! medium properties.
    ! See the overall program documentation for output parameters.

    implicit none
    integer :: NSTOKES, NX, NY, NZ, NBPTS, NPTS, NCELLS, NSH
    integer :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, NMU, NPHI, NANG, NG
    integer :: BCFLAG, IPFLAG
    integer :: MAXITER, TOTITER, IRAD, NSHOUT
    integer :: TREEPTR(2, NCELLS), GRIDPTR(8, NCELLS)
    integer(2) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical :: DELTAM
    real :: SOLARFLUX, SOLARMU, SOLARAZ
    real :: GNDTEMP, GNDALBEDO, SKYRAD
    real :: SOLACC, SPLITACC, SHACC
    real :: WAVENO(2), WAVELEN
    real :: XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real :: OUTPARMS(*)
    character :: SRCTYPE(1), SFCTYPE(2), UNITS(1), GRIDTYPE(1), OUTTYPE(1)
    character(256) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 OUTFILE
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real TEMP(NPTS)
    real FLUXES(2, NPTS), DIRFLUX(NPTS), FLUXDIV(NPTS)
    real RADOUT(NSTOKES, *), SHTERMS(NSHOUT, NPTS)
    real SOURCE1OUT(NSTOKES, *)

    integer IZ, I, K, N, NX2, NY2, J, JX, JY, ICELL
    integer NANGOUT, NXOUT, NYOUT
    logical ALLGRID
    real :: MUOUT, PHID, C, PI, SUM, SUM1, SUM2, SUM3
    real :: XDOMAIN, YDOMAIN, STARTX, STARTY
    real :: FUP, FDN, FDIR, ASYM
    real(8) :: X0, Y0, Z0
    character(88) :: FORM
    character(88) :: GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

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
        & '    NPTS=', NPTS, '   NCELLS=', NCELLS
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
    write (2, '(A,A22,A,I1)') '!  SURFACE_TYPE=', SFCNAME,&
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

    PI = acos(-1.0)

    if (OUTTYPE .eq. 'R') then

      call shdUtil__outputDataRad()

    else if (OUTTYPE .eq. 'F') then

      call shdUtil__outputDataFlx()

    else if (OUTTYPE .eq. 'H') then

        ! Heating output: net flux convergence
        if (nint(OUTPARMS(1)) .eq. 1) then
            ! Format 1 is domain averaged vertical profile
            write (2, '(A,A)') '!    Z       -DIV(Fnet)'
            do IZ = 1, NZ
                SUM = 0.0
                N = 0
                do I = IZ, NBPTS, NZ
                    if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                        & GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
                        SUM = SUM + FLUXDIV(I)
                        N = N + 1
                    end if
                end do
                write (2, '(1X,F7.3,2X,E13.5)') ZGRID(IZ), -SUM/N
            end do

        else
            if (nint(OUTPARMS(1)) .eq. 3) then
                N = NPTS
                ALLGRID = .true.
            else
                N = NBPTS
                ALLGRID = .false.
            end if
            ! Format 2 is flux convergence for every base grid point
            ! Format 3 is flux convergence for every grid point
            write (2, '(A)') '!    X       Y        Z       -DIV(Fnet)'
            do I = 1, N
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                    & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    write (2, '(3(1X,F7.3),2X,E13.5)')&
                        & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I), -FLUXDIV(I)
                end if
            end do
        end if

    else if (OUTTYPE .eq. 'S') then

        ! Spherical Harmonic output:
        !   Output mean intensity and net flux (x,y,z components).
        !   and maybe normalized rms of higher order terms
        if (nint(OUTPARMS(1)) .eq. 2) then
            N = NPTS
            ALLGRID = .true.
        else
            N = NBPTS
            ALLGRID = .false.
        end if
        if (NSHOUT .eq. 5) then
            write (2, '(A,A)') '!    X       Y        Z         Imean',&
                & '          Fx           Fy           Fz       HOrms'
            do I = 1, N
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                    & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
                        & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
                        & SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I),&
                        & SHTERMS(5, I)/max(1e-20, SHTERMS(1, I))
                end if
            end do
        else
            write (2, '(A,A)') '!    X       Y        Z         Imean',&
                & '          Fx           Fy           Fz '
            do I = 1, N
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                    & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    write (2, '(3(1X,F7.3),1X,4(1X,E12.5))')&
                        & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
                        & SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I)
                end if
            end do
        end if

    else if (OUTTYPE .eq. 'J') then

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
            & '! ', MUOUT, PHID, '<- (mu,phi)'
        write (2, '(A,A)') '!    X       Y       Z    ',&
            & ' Extinction   Source Function (Stokes I, Q, U, V)'
        do I = 1, N
            if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                write (2, '(3(1X,F7.3),5(1X,E12.5))')&
                    & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
                    & EXTINCT(I), SOURCE1OUT(:, I)
            end if
        end do

    else if (OUTTYPE .eq. 'M') then

        ! Medium output: grid point medium properties -
        !   Output extinction, single scattering albedo,
        !   asymmetry parameter, and temperature.
        if (nint(OUTPARMS(1)) .eq. 2) then
            N = NPTS
            ALLGRID = .true.
        else
            N = NBPTS
            ALLGRID = .false.
        end if
        write (2, '(A,A)') '!    X       Y       Z',&
            & '      Extinct  Albedo  Asymmetry Temp'
        do I = 1, N
            if (NUMPHASE .gt. 0) then
                ASYM = LEGEN(1, 1, IPHASE(I))
            else
                ASYM = LEGEN(1, 1, I)
            end if
            if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.&
                & GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                write (2, '(3(1X,F7.3),1X,F10.6,1X,F8.6,1X,F8.5,1X,F7.2)')&
                    & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),&
                    & EXTINCT(I), ALBEDO(I), ASYM, TEMP(I)
            end if
        end do

    end if

    close (2)
    return
end subroutine

subroutine shdUtil__outputDesc(&
    & NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    & NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
    & PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
    & BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
    & SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    & SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
    & SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
    & XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID,&
    & FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT,&
    & OUTTYPE, OUTPARMS, OUTFILE)
    ! Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are four supported types (OUTTYPE) of output: 'R' - for radi
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms
    ! See the overall program documentation for output parameters.
    implicit none
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
    character SRCTYPE*1, SFCTYPE*2, UNITS*1, GRIDTYPE*1, OUTTYPE*1
    character*256 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 OUTFILE
    real FLUXES(3, NZ, NY, NX), FLUXDIV(NZ, NY, NX)
    real RADOUT(NSTOKES, *), SHTERMS(NSHOUT, NZ, NY, NX)

    integer I, J, K, N, JX, JY
    integer NANGOUT, NXOUT, NYOUT
    real MUOUT, PHID, C, PI, SUM, SUM1, SUM2, SUM3
    real STARTX, STARTY
    double precision X0, Y0, Z0
    character*128 FORM
    character*64 GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

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

    PI = acos(-1.0)

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

        write (2, '(A,F7.3,3(A,I4))') '!    RADIANCE AT Z=', Z0,&
            & '    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (UNITS .eq. 'T') then
            FORM = '(2(1X,F7.3),1X,4(1X,F6.2))'
            write (2, '(A)')&
                & '!   X       Y      RADIANCE (STOKES I, Q, U, V)'
        else
            FORM = '(2(1X,F7.3),1X,4(1X,E12.5))'
            write (2, '(A)')&
                & '!   X       Y        RADIANCE (STOKES I, Q, U, V)'
        end if

        do K = 1, NANGOUT
            MUOUT = OUTPARMS(2*K + 5)
            PHID = OUTPARMS(2*K + 6)
            write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')&
                & '! ', MUOUT, PHID, '<- (mu,phi)'
            Y0 = STARTY
            do JY = 1, NYOUT
                X0 = STARTX
                do JX = 1, NXOUT
                    IRAD = IRAD + 1
                    write (2, FORM) X0, Y0, RADOUT(:, IRAD)
                    X0 = X0 + OUTPARMS(2)
                end do
                Y0 = Y0 + OUTPARMS(3)
            end do
        end do

    else if (OUTTYPE .eq. 'F') then
        ! Hemispheric flux output
        ! Format 1 is flux at top and bottom of medium at grid point
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
            C = 1.0
            if (UNITS .eq. 'T') C = 1.0/PI
            do I = 1, NX
                do J = 1, NY
                    if (SRCTYPE .ne. 'T') then
                        write (2, FORM) XGRID(I), YGRID(J),&
                            & C*FLUXES(2, NZ, J, I), C*FLUXES(1, 1, J, I), FLUXES(3, 1, J, I)
                    else
                        write (2, FORM) XGRID(I), YGRID(J),&
                            & C*FLUXES(2, NZ, J, I), C*FLUXES(1, 1, J, I)
                    end if
                end do
            end do

        else if (nint(OUTPARMS(1)) .eq. 2) then
        ! Format 2 is flux at a given level at regular locations
        ! Not supported for multiple processors

        else if (nint(OUTPARMS(1)) .eq. 3) then
            ! Format 3 is domain averaged vertical profile
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
            C = 1.0
            if (UNITS .eq. 'T') C = 1.0/PI
            N = NX*NY
            do K = 1, NZ
                SUM1 = 0.0
                SUM2 = 0.0
                SUM3 = 0.0
                do I = 1, NX
                    do J = 1, NY
                        SUM1 = SUM1 + FLUXES(1, K, J, I)
                        SUM2 = SUM2 + FLUXES(2, K, J, I)
                        SUM3 = SUM3 + FLUXES(3, K, J, I)
                    end do
                end do
                if (SRCTYPE .ne. 'T') then
                    write (2, FORM) ZGRID(K), C*SUM2/N, C*SUM1/N, SUM3/N
                else
                    write (2, FORM) ZGRID(K), C*SUM2/N, C*SUM1/N
                end if
            end do

        else
            ! Format 4 is fluxes at every base grid point
            ! Format 5 is fluxes at every grid point (not supported)
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
            C = 1.0
            if (UNITS .eq. 'T') C = 1.0/PI
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        if (SRCTYPE .ne. 'T') then
                            write (2, FORM) XGRID(I), YGRID(J), ZGRID(K),&
                                & C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I), FLUXES(3, K, J, I)
                        else
                            write (2, FORM) XGRID(I), YGRID(J), ZGRID(K),&
                                & C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I)
                        end if
                    end do
                end do
            end do
        end if

    else if (OUTTYPE .eq. 'H') then

        ! Heating output: net flux convergence
        if (nint(OUTPARMS(1)) .eq. 1) then
            ! Format 1 is domain averaged vertical profile
            write (2, '(A,A)') '!    Z       -DIV(Fnet)'
            do K = 1, NZ
                SUM = 0.0
                do I = 1, NX
                    do J = 1, NY
                        SUM = SUM + FLUXDIV(K, J, I)
                    end do
                end do
                write (2, '(1X,F7.3,2X,E13.5)') ZGRID(K), -SUM/(NX*NY)
            end do

        else

          ! Format 2 is flux convergence for every base grid point
          ! Format 3 is flux convergence for every grid point (not sup
            write (2, '(A)') '!    X       Y        Z       -DIV(Fnet)'
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        write (2, '(3(1X,F7.3),2X,E13.5)')&
                            & XGRID(I), YGRID(J), ZGRID(K), -FLUXDIV(K, J, I)
                    end do
                end do
            end do
        end if

    else if (OUTTYPE .eq. 'S') then

        ! Spherical Harmonic output:
        !   Output mean intensity and net flux (x,y,z components).
        !   and maybe normalized rms of higher order terms
        if (NSHOUT .eq. 5) then
            write (2, '(A,A)') '!    X       Y        Z         Imean',&
                & '          Fx           Fy           Fz       HOrms'
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
                            & XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1, K, J, I),&
                            & SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I),&
                            & SHTERMS(5, K, J, I)/max(1e-20, SHTERMS(1, K, J, I))
                    end do
                end do
            end do
        else
            write (2, '(A,A)') '!    X       Y        Z         Imean',&
                & '          Fx           Fy           Fz '
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')&
                            & XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1, K, J, I),&
                            & SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I)
                    end do
                end do
            end do
        end if

    else if (OUTTYPE .eq. 'J') then

    ! J source function output: source function for a specified ang
    !   at every grid point (not supported)

    else if (OUTTYPE .eq. 'M') then

    ! Medium output: grid point medium properties -
    !   Output extinction, single scattering albedo,
    !   asymmetry parameter, and temperature.  (not supported)

    end if

    close (2)
    return
end subroutine

end module shdUtil_io
