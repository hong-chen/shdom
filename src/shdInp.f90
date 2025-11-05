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

module shd_inp

    use globals_const, only: PI_
    use globals_cfg,   only: verbose_, charMaxLen_
    use globals_type,  only: shd_type_rtm, &
                             shd_type_grid, &
                             shd_type_shdo, &
                             shd_type_src, &
                             shd_type_prp, &
                             shd_type_sfc, &
                             shd_type_ckd, &
                             shd_type_sen

    use shd_util, only: shd_util__computePlanck

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public :: shd_io__userInp
    public :: shd_io__nmlInp
    public :: shd_io__checkInp

    public :: shd_io__readPrpSize
    public :: shd_io__readPrp
    public :: shd_io__checkPrp

    public :: shd_io__readSfcSize
    public :: shd_io__readSfc

    public :: shd_io__readCkdSize
    public :: shd_io__readCkd

    public :: shd_io__readSenSize
    public :: shd_io__readSen

    public :: shd_io__restoreState

contains

subroutine shd_io__userInp(&
    &RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
    &INSAVEFILE, OUTSAVEFILE, SENFILE,&
    &NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
    &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
    &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
    &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
    &OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, VERBOSE)

    ! Obtains the input parameters for the program by writing prompts
    ! to and reading responses from the terminal.  Echos the inputs to
    ! make a useful log file when running non-interactively.  See the
    ! overall program documentation for the list of input parameters.

    integer NSTOKES, NX, NY, NZ, NMU, NPHI
    integer BCFLAG, IPFLAG
    integer MAXITER, MAXOUT, MAXPAR, NUMOUT
    logical KDIST, DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real SOLACC, SPLITACC, SHACC
    real OUTPARMS(MAXPAR, MAXOUT)
    real MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(len=1) :: SRCTYPE, GRIDTYPE, UNITS
    character(len=1) :: OUTTYPES(*)
    character(len=charMaxLen_) :: OUTFILES(*)
    character(len=charMaxLen_) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, SENFILE
    character(len=charMaxLen_) :: OutFileNC, RUNNAME
    integer(kind=1) :: VERBOSE
    integer i, j, N
    real GNDEMIS

    SOLARMU = -1.0
    GNDALBEDO = 0.0
    WAVELEN = 1.0
    WAVENO(1) = 10000.
    WAVENO(2) = 10001.
    SENFILE = 'NONE'
    VERBOSE = verbose_

    write (*, *)
    write (*, '(1X,A)') 'Name of this SHDOM run:'
    read (*, '(A)') RUNNAME
    write (*, '(A)') RUNNAME

    write (*, '(1X,A)') 'Atmospheric optical properties file name'
    read (*, '(A)') PROPFILE
    write (*, *) PROPFILE
    write (*, '(1X,A)') 'Surface properties file name (or NONE)'
    read (*, '(A)') SFCFILE
    write (*, *) SFCFILE
    write (*, '(1X,A)')&
        &'Correlated k-distribution file name (or NONE)'
    read (*, '(A)') CKDFILE
    write (*, *) CKDFILE
    KDIST = CKDFILE(1:4) .ne. 'NONE'
    write (*, '(1X,A)') 'Input binary save file name'
    read (*, '(A)') INSAVEFILE
    write (*, *) INSAVEFILE
    write (*, '(1X,A)') 'Output binary save file name'
    read (*, '(A)') OUTSAVEFILE
    write (*, *) OUTSAVEFILE

    write (*, '(1X,A)')&
        &'Number of Stokes parameters (1, 3, or 4; 1=unpolarized)'
    read (*, *) NSTOKES
    write (*, *) NSTOKES
    write (*, '(1X,A)') 'Internal base grid size (Nx, Ny, Nz)'
    read (*, *) NX, NY, NZ
    write (*, *) NX, NY, NZ
    write (*, '(1X,A)') 'Number of discrete ordinates in mu and phi'
    read (*, *) NMU, NPHI
    write (*, *) NMU, NPHI
    write (*, '(1X,A)')&
        &'Boundary conditions (0=periodic, 3=open)'
    read (*, *) BCFLAG
    write (*, *) BCFLAG

    write (*, '(1X,A)') 'Independent pixel approximation (0,1,2,3)'
    read (*, *) IPFLAG
    write (*, *) IPFLAG
    write (*, '(1X,A)') 'Delta-M method (T,F)'
    read (*, *) DELTAM
    write (*, *) DELTAM
    write (*, '(1X,A)') 'Z grid type (E-even, P-property file)'
    read (*, '(A)') GRIDTYPE
    write (*, *) GRIDTYPE

    write (*, '(1X,A)') 'Thermal, solar, or both source (T, S, B)'
    read (*, '(A)') SRCTYPE
    write (*, *) SRCTYPE
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        write (*, '(1X,A)') 'Solar flux and direction (F, mu, phi[deg])'
        read (*, *) SOLARFLUX, SOLARMU, SOLARAZ
        write (*, *) SOLARFLUX, SOLARMU, SOLARAZ
        write (*, '(1X,A)') 'Isotropic sky radiance'
        read (*, *) SKYRAD
        write (*, *) SKYRAD
        if (SFCFILE(1:4) .eq. 'NONE') then
            write (*, '(1X,A)') 'Ground albedo'
            read (*, *) GNDALBEDO
            write (*, *) GNDALBEDO
        end if
        SOLARMU = -abs(SOLARMU)
        SOLARAZ = SOLARAZ*acos(-1.0)/180.0
    end if
    GNDTEMP = 300.
    if (SFCFILE(1:4) .eq. 'NONE') then
        if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
            write (*, '(1X,A)') 'Ground temperature'
            read (*, *) GNDTEMP
            write (*, *) GNDTEMP
        end if
        if (SRCTYPE .eq. 'T') then
            write (*, '(1X,A)') 'Ground emissivity'
            read (*, *) GNDEMIS
            write (*, *) GNDEMIS
            GNDALBEDO = 1.0 - GNDEMIS
        end if
    end if
    if (SRCTYPE .eq. 'T') then
        write (*, '(1X,A)') 'Sky temperature'
        read (*, *) SKYRAD
        write (*, *) SKYRAD
    end if

    if (KDIST) then
        write (*, '(1X,A)') 'Wavenumber range (cm^-1)'
        read (*, *) WAVENO(1), WAVENO(2)
        write (*, *) WAVENO(1), WAVENO(2)
        if (WAVENO(2) .gt. WAVENO(1)) then
            UNITS = 'B'
        else
            UNITS = 'R'
        end if
        WAVELEN = 1.0e4/(0.5*(WAVENO(1) + WAVENO(2)))
    else
        if (SRCTYPE .eq. 'T') then
            write (*, '(1X,A)') 'Radiance units (R-radiance, T-temp)'
            read (*, '(A)') UNITS
            write (*, *) UNITS
            if (UNITS .ne. 'T') then
                write (*, '(1X,A)') 'Wavelength (micron)'
                read (*, *) WAVELEN
                write (*, *) WAVELEN
            end if
        else
            UNITS = 'R'
            write (*, '(1X,A)') 'Wavelength (micron)'
            read (*, *) WAVELEN
            write (*, *) WAVELEN
        end if
    end if

    write (*, '(1X,A)')&
        &'Cell splitting and spherical harmonic series accuracy'
    read (*, *) SPLITACC, SHACC
    write (*, *) SPLITACC, SHACC

    write (*, '(1X,A,A)') 'Acceleration flag (T or F), ',&
        &'solution accuracy, and max iterations'
    read (*, *) ACCELFLAG, SOLACC, MAXITER
    write (*, *) ACCELFLAG, SOLACC, MAXITER

    write (*, '(1X,A)') 'Number of output files'
    read (*, *) NUMOUT
    write (*, *) NUMOUT
    if (NUMOUT .gt. MAXOUT) stop 'Error [_userInp]: MAXOUT exceeded.'
    do i = 1, NUMOUT
        write (*, '(1X,A,A)') 'Output type (R-radiance, ',&
            &'V-visualization, F-flux, H-heating,'
        write (*, '(14X,A)') 'S-spherical harmonic, J-source, M-medium)'
        read (*, '(A)') OUTTYPES(i)
        write (*, *) OUTTYPES(i)
        if (OUTTYPES(i) (1:1) .eq. 'R') then
            write (*, '(1X,A,A)')&
                &'Z level, X and Y spacing, X and Y offset, ',&
                &'number of angles, and mu, phi pairs'
            read (*, *) OUTPARMS(1, i), OUTPARMS(2, i), OUTPARMS(3, i),&
                &OUTPARMS(4, i), OUTPARMS(5, i),&
                &N, (OUTPARMS(j + 6, i), j=1, min(2*N, MAXPAR - 6))
            if (N .gt. (MAXPAR - 6)/2) then
                write (*, *) 'Error [_userInp]: Number of radiance angles is too large for MAXPAR'
                stop
            end if
            N = min(N, (MAXPAR - 6)/2)
            OUTPARMS(6, i) = float(N)
            write (*, *) (OUTPARMS(j, i), j=1, 2*N + 6)
        else if (OUTTYPES(i) (1:1) .eq. 'V') then
            write (*, '(1X,A,A)')&
                &'1, bytes,scale, X,Y,Z, theta,phi,rotang, NL,NS, delline,delsamp'
            write (*, '(1X,A,A)')&
                &'2, bytes,scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1,scan2,delscan'
            read (*, *) (OUTPARMS(j, i), j=1, 13)
            write (*, *) (OUTPARMS(j, i), j=1, 13)
        else if (OUTTYPES(i) (1:1) .eq. 'F') then
            write (*, '(1X,A)')&
                &'Format (1,2,3,4,5), [Z, del X, del Y for format 2]'
            read (*, *) N, (OUTPARMS(1 + j, i), j=1, 3*max(0, 1 - abs(N - 2)))
            OUTPARMS(1, i) = float(N)
            if (N .eq. 2) then
                write (*, *) N, (OUTPARMS(1 + j, i), j=1, 3)
            else
                write (*, *) N
            end if
        else if (OUTTYPES(i) (1:1) .eq. 'H') then
            write (*, '(1X,A)') 'Format (1,2,3)'
            read (*, *) OUTPARMS(1, i)
            write (*, *) OUTPARMS(1, i)
        else if (OUTTYPES(i) (1:1) .eq. 'S') then
            write (*, '(1X,A)') 'Format (1,2)'
            read (*, *) OUTPARMS(1, i)
            write (*, *) OUTPARMS(1, i)
        else if (OUTTYPES(i) (1:1) .eq. 'J') then
            write (*, '(1X,A)') 'Format (1,2),  mu, phi '
            read (*, *) OUTPARMS(1, i), OUTPARMS(2, i), OUTPARMS(3, i)
            write (*, *) OUTPARMS(1, i), OUTPARMS(2, i), OUTPARMS(3, i)
        else if (OUTTYPES(i) (1:1) .eq. 'M') then
            write (*, '(1X,A)') 'Format (1,2)'
            read (*, *) OUTPARMS(1, i)
            write (*, *) OUTPARMS(1, i)
        end if
        write (*, '(1X,A,A)') 'Output file name'
        read (*, '(A)') OUTFILES(i)
        write (*, '(A)') OUTFILES(i)
    end do

    write (*, '(1X,A)') 'Output netcdf file name (or NONE or "")'
    read (*, '(A)') OutFileNC
    write (*, '(A)') OutFileNC
    if (OutFileNC(1:4) .ne. 'NONE') then
        if (verbose_ .ge. 4) write(*, '(A)') 'Warning [_userInp]: Removed support for netCDF, set <OUTFILENC="NONE">.'
        OutFileNC = 'NONE'
    end if

    write (*, '(1X,A)')&
        &'Approximate maximum memory to use (MB for 4 byte reals)'
    read (*, *) MAX_TOTAL_MB
    write (*, *) MAX_TOTAL_MB
    write (*, '(1X,A)')&
        &'Ratio of total grid points to base grid points (1.0 to ~10)'
    read (*, *) ADAPT_GRID_FACTOR
    write (*, *) ADAPT_GRID_FACTOR
    write (*, '(1X,A,A)') 'Ratio of average number of spherical',&
        &' harmonic terms to total possible (<=1.0)'
    read (*, *) NUM_SH_TERM_FACTOR
    write (*, *) NUM_SH_TERM_FACTOR
    write (*, '(1X,A)')&
        &'Ratio of number of grid cells to grid points (1.0 to 2.0)'
    read (*, *) CELL_TO_POINT_RATIO
    write (*, *) CELL_TO_POINT_RATIO

    write (*, *)
    return
end subroutine shd_io__userInp

subroutine shd_io__nmlInp(&
    &dRUNNAME, dPROPFILE, dSFCFILE, dCKDFILE,&
    &dINSAVEFILE, dOUTSAVEFILE, dSENFILE,&
    &dNSTOKES, dNX, dNY, dNZ, dNMU, dNPHI,&
    &dBCFLAG, dIPFLAG,&
    &dKDIST, dDELTAM, dGRIDTYPE,&
    &dSRCTYPE, dSOLARFLUX, dSOLARMU, dSOLARAZ, dSKYRAD,&
    &dGNDTEMP, dGNDALBEDO, dUNITS, dWAVENO, dWAVELEN,&
    &dACCELFLAG, dSOLACC, dMAXITER, dSPLITACC, dSHACC,&
    &dMAXOUT, dMAXPAR, dNUMOUT, dOUTTYPES, dOUTPARMS, dOUTFILES,&
    &dOutFileNC, dMAX_TOTAL_MB, dADAPT_GRID_FACTOR,&
    &dNUM_SH_TERM_FACTOR, dCELL_TO_POINT_RATIO,&
    &dVERBOSE)

    ! Obtains the input parameters for the program from a namelist
    ! read from stdin.  See the overall program documentation for the
    ! list of input parameters.  The subroutine dummy variables can't
    ! be in the namelist, hence the two sets of variables.

    integer dMAXOUT, dMAXPAR
    integer dNSTOKES, dNX, dNY, dNZ, dNMU, dNPHI
    integer dBCFLAG, dIPFLAG
    integer dMAXITER, dNUMOUT
    logical dKDIST, dDELTAM, dACCELFLAG
    real dSOLARFLUX, dSOLARMU, dSOLARAZ
    real dGNDTEMP, dGNDALBEDO, dSKYRAD
    real dWAVENO(2), dWAVELEN, dSOLACC, dSPLITACC, dSHACC
    real dOUTPARMS(dMAXPAR, dMAXOUT)
    real dMAX_TOTAL_MB, dADAPT_GRID_FACTOR
    real dNUM_SH_TERM_FACTOR, dCELL_TO_POINT_RATIO
    character(len=1) :: dSRCTYPE, dGRIDTYPE, dUNITS
    character(len=1) ::  dOUTTYPES(*)
    character(len=charMaxLen_) :: dOUTFILES(*)
    character(len=charMaxLen_) :: dOutFileNC, dRUNNAME
    character(len=charMaxLen_) :: dPROPFILE, dSFCFILE, dCKDFILE
    character(len=charMaxLen_) :: dINSAVEFILE, dOUTSAVEFILE, dSENFILE
    integer(kind=1) :: dVERBOSE
    integer i, j, N

    integer MAXOUT, MAXPAR
    parameter(MAXOUT=20, MAXPAR=100)
    integer NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    integer MAXITER, NUMOUT
    logical DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real WAVENO(2), WAVELEN, SOLACC, SPLITACC, SHACC
    real MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real OUTPARMS(MAXPAR, MAXOUT)
    character(len=1) :: SRCTYPE, GRIDTYPE, UNITS
    character(len=1) :: OUTTYPES(MAXOUT)
    character(len=charMaxLen_) :: OUTFILES(MAXOUT)
    character(len=charMaxLen_) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=charMaxLen_) :: SENFILE, OutFileNC, RUNNAME
    integer(kind=1) :: VERBOSE
    namelist /SHDOMINPUT/&
        &RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
        &INSAVEFILE, OUTSAVEFILE, SENFILE,&
        &NSTOKES, NX, NY, NZ, NMU, NPHI,&
        &BCFLAG, IPFLAG, DELTAM, GRIDTYPE,&
        &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
        &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
        &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
        &NUMOUT, OUTTYPES, OUTPARMS, OUTFILES, OutFileNC,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO,&
        &VERBOSE

    character(len=charMaxLen_) :: NMLFILE

    if (dMAXOUT .ne. MAXOUT) stop 'Error [_nmlInp]: MAXOUT not the same.'
    if (dMAXPAR .ne. MAXPAR) stop 'Error [_nmlInp]: MAXPAR not the same.'

    SOLARMU = -1.0
    GNDALBEDO = 0.0
    WAVELEN = 1.0
    WAVENO(1) = 10000.
    WAVENO(2) = 10001.
    SENFILE = 'NONE'
    VERBOSE = verbose_

    ! Read in the namelist
    call get_command_argument(1, NMLFILE)
    if (verbose_ .ge. 6) write (*, '(A,"<",A,">")') 'Message [_nmlInp]: ', trim(adjustl(NMLFILE))
    open (UNIT=5, FILE=trim(adjustl(NMLFILE)), STATUS='OLD')
    read (5, NML=SHDOMINPUT)
    close (5)

    ! Make the assignments to the dummy variables
    dRUNNAME = RUNNAME
    dPROPFILE = PROPFILE
    dSFCFILE = SFCFILE
    dCKDFILE = CKDFILE
    dKDIST = CKDFILE(1:4) .ne. 'NONE'
    dINSAVEFILE = INSAVEFILE
    dOUTSAVEFILE = OUTSAVEFILE
    dSENFILE = SENFILE ! this is added for user-specified sensor option, output type [R] and [V3]
    dNSTOKES = NSTOKES
    dNX = NX
    dNY = NY
    dNZ = NZ
    dNMU = NMU
    dNPHI = NPHI
    dIPFLAG = IPFLAG
    dBCFLAG = BCFLAG
    dDELTAM = DELTAM
    dGRIDTYPE = GRIDTYPE

    dSRCTYPE = SRCTYPE
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        dSOLARFLUX = SOLARFLUX
        dSOLARMU = -abs(SOLARMU)
        dSOLARAZ = SOLARAZ*acos(-1.0)/180.0
    else
        dSOLARFLUX = 0.0
        dSOLARMU = 1.0
        dSOLARAZ = 0.0
    end if
    dSKYRAD = SKYRAD
    dGNDTEMP = GNDTEMP
    dGNDALBEDO = GNDALBEDO

    if (dKDIST) then
        dUNITS = 'B'
        dWAVENO(1) = WAVENO(1)
        dWAVENO(2) = WAVENO(2)
        dWAVELEN = 1.0d4/(0.5d0*(WAVENO(1) + WAVENO(2)))
    else
        dUNITS = UNITS
        dWAVELEN = WAVELEN
        if (SRCTYPE .ne. 'T') then
            dUNITS = 'R'
        end if
    end if

    dSPLITACC = SPLITACC
    dSHACC = SHACC
    dACCELFLAG = ACCELFLAG
    dSOLACC = SOLACC
    dMAXITER = MAXITER

    dNUMOUT = NUMOUT
    if (NUMOUT .gt. MAXOUT) stop 'Error [_nmlInp]: MAXOUT exceeded.'
    do i = 1, NUMOUT
        dOUTTYPES(i) = OUTTYPES(i)
        if (OUTTYPES(i) (1:1) .eq. 'R') then
            ! new binary angle input do not need the loop above
            if (SENFILE(1:4) .ne. 'NONE') then
                do j = 1, 6
                    dOUTPARMS(j, i) = OUTPARMS(j, i)
                end do
            else
                N = nint(OUTPARMS(6, i))
                do j = 1, min(2*N + 6, MAXPAR)
                    dOUTPARMS(j, i) = OUTPARMS(j, i)
                end do
            end if
        else if (OUTTYPES(i) (1:1) .eq. 'V') then
            dOUTPARMS(1, i) = OUTPARMS(1, i)
            dOUTPARMS(2, i) = OUTPARMS(2, i)
            dOUTPARMS(3, i) = OUTPARMS(3, i)
            ! new user-specified sensor mode V3 only need 3 parameters
            ! read more parameters if for V1 and V2
            if (nint(OUTPARMS(1, i)) .le. 2) then
                do j = 4, 13
                    dOUTPARMS(j, i) = OUTPARMS(j, i)
                end do
            end if
        else if (OUTTYPES(i) (1:1) .eq. 'F') then
            dOUTPARMS(1, i) = OUTPARMS(1, i)
            if (nint(OUTPARMS(1, i)) .eq. 2) then
                do j = 2, 4
                    dOUTPARMS(j, i) = OUTPARMS(j, i)
                end do
            end if
        else if (OUTTYPES(i) (1:1) .eq. 'H') then
            dOUTPARMS(1, i) = OUTPARMS(1, i)
        else if (OUTTYPES(i) (1:1) .eq. 'S') then
            dOUTPARMS(1, i) = OUTPARMS(1, i)
        else if (OUTTYPES(i) (1:1) .eq. 'J') then
            dOUTPARMS(1, i) = OUTPARMS(1, i)
            dOUTPARMS(2, i) = OUTPARMS(2, i)
            dOUTPARMS(3, i) = OUTPARMS(3, i)
        else if (OUTTYPES(i) (1:1) .eq. 'M') then
            dOUTPARMS(1, i) = OUTPARMS(1, i)
        end if
        dOUTFILES(i) = OUTFILES(i)
    end do

    dOutFileNC = OutFileNC

    if (dOutFileNC(1:4) .ne. 'NONE') then
        if (verbose_ .ge. 4) write(*, '(A)') 'Warning [_nmlInp]: Removed support for netCDF, set <OUTFILENC="NONE">.'
        dOutFileNC = 'NONE'
    end if

    dMAX_TOTAL_MB = MAX_TOTAL_MB
    dADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR
    dNUM_SH_TERM_FACTOR = NUM_SH_TERM_FACTOR
    dCELL_TO_POINT_RATIO = CELL_TO_POINT_RATIO
    dVERBOSE = VERBOSE

    return
end subroutine shd_io__nmlInp

subroutine shd_io__checkInp(&
    &DELTAM,&
    &SPLITACC, SHACC,&
    &rtm_obj, shdo_obj, src_obj, prp_obj, sfc_obj)

    ! Checks for validity of the input parameters. Gives warning for
    ! parameters that might not be correct.
    ! Check the letter flags.
    ! Delta-M should be used for solar transfer with forward peaked
    ! phase functions.
    ! Check for reasonable ranges on SPLITACC and SHACC.
    ! Check for reasonable range of sfc_obj%avgTemp and valid range of sfc_obj%avgAlb.
    ! Check for src_obj%avgMu=0.
    ! Check for shdo_obj%Nphi<shdo_obj%Nmu.

    type(shd_type_rtm),  intent(in) :: rtm_obj
    type(shd_type_shdo), intent(in) :: shdo_obj
    type(shd_type_src),  intent(in) :: src_obj
    type(shd_type_prp),  intent(in) :: prp_obj
    type(shd_type_sfc),  intent(in) :: sfc_obj

    logical DELTAM
    real SPLITACC, SHACC
    real FLUX

    if ((src_obj%type .ne. 'S') .and. (src_obj%type .ne. 'T') .and. (src_obj%type .ne. 'B'))&
        &stop 'Error [_checkInp]: Illegal source type.'
    if ((src_obj%unit .ne. 'R') .and. (src_obj%unit .ne. 'T') .and. (src_obj%unit .ne. 'B'))&
        &stop 'Error [_checkInp]: Illegal units.'
    if ((rtm_obj%type .ne. 'E') .and. (rtm_obj%type .ne. 'P') .and. (rtm_obj%type .ne. 'F'))&
        &stop 'Error [_checkInp]: Illegal grid type.'

    if ((src_obj%type .ne. 'T') .and. (src_obj%avgMu .eq. 0.0)) then
        stop 'Error [_checkInp]: Cosine solar zenith angle&
            & (src_obj%avgMu) cannot be 0.'
    end if

    if (.not. DELTAM .and. (src_obj%type .ne. 'T') .and. prp_obj%maxAsym .gt. 0.5) then
        if (verbose_ .ge. 4) write (*, '(A)') 'Warning [_checkInp]: Delta-M should be used&
            & (DELTAM=.TRUE.) for solar radiative transfer&
            & with highly peaked phase functions.'
    end if

    ! Check for SPLITACC and SHACC magnitudes
    if ((src_obj%type .eq. 'T') .or. (src_obj%type .eq. 'B')) then
        call shd_util__computePlanck(sfc_obj%avgTemp, src_obj%unit, src_obj%waveNum, src_obj%avgWaveLen, FLUX)
        FLUX = PI_*FLUX
    else
        FLUX = 0.0
    end if
    if ((src_obj%type .eq. 'S') .or. (src_obj%type .eq. 'B')) then
        FLUX = FLUX + src_obj%inpSolar
    end if
    if (SPLITACC .gt. 0.0 .and.&
        &(SPLITACC .lt. 0.001*FLUX .or. SPLITACC .gt. 0.1*FLUX)) then
        if (verbose_ .ge. 4) write (*, '(A)') 'Warning [_checkInp]: Splitting accuracy parameter&
            & is not a fraction, but scales with fluxes.'
    end if
    if (SHACC .gt. 0.0 .and. SHACC .gt. 0.03*FLUX) then
        if (verbose_ .ge. 4) write (*, '(A)') 'Warning [_checkInp]: Spherical harmonic accuracy&
            & (SHACC) is not a fraction, but scales with fluxes.'
    end if

    if ((src_obj%type .ne. 'S') .and.&
        &(sfc_obj%avgTemp .lt. 150. .or. sfc_obj%avgTemp .gt. 350.)) then
        if (verbose_ .ge. 4) write (*, '(A)') 'Warning [_checkInp]:&
            & ground temperature out of Earth range.'
    end if
    if (sfc_obj%avgAlb .lt. 0.0 .or. sfc_obj%avgAlb .gt. 1.0) then
        stop 'Error [_checkInp]: Ground albedo must be 0 to 1.'
    end if

    if (shdo_obj%Nphi .lt. shdo_obj%Nmu) then
      if (verbose_ .ge. 4) write (*, '(A)') 'Warning [_checkInp]: Usually want shdo_obj%Nphi>= shdo_obj%Nmu.'
    end if

    return
end subroutine shd_io__checkInp

subroutine shd_io__readPrpSize(prp_obj)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Reads parts of the property file to get the maximum array sizes
    !   needed for allocatable arrays.  For extinction only and tabulated
    !   phase function formats, only the header is read, while for the
    !   standard format, the whole file must be read to determine prp_obj%maxNleg.
    !
    ! Input:
    !   prp_obj%fname: string, file path for the property file (contains extinction, phase function etc.
    !          generated by <propgen>)
    !   prp_obj%NlegSet: integer, the number of phase function elements used, =1 when unpolarized, =6 when polarized
    !   prp_obj%Nleg: integer, the number of Legendre series phase function terms used in the solution procedure
    !
    ! Output:
    !   prp_obj%Nx: integer, number of grid points in X
    !   prp_obj%Ny: integer, number of grid points in Y
    !   prp_obj%Nz: integer, number of grid points in Z
    !   prp_obj%dx: real, X grid resolution
    !   prp_obj%dy: real, Y grid resolution
    !   prp_obj%Npha: integer, number of tabulated phase functions
    !   prp_obj%maxNleg: integer, maximum Legendre series order for phase functions
    !   prp_obj%maxNptsLeg: integer, maximum data points to store all phase functions
    !╰────────────────────────────────────────────────────────────────────────────╯!

    type(shd_type_prp), intent(inout) :: prp_obj

    integer(kind=I__) :: Nleg, ix, iy, iz, i, j, k, l
    real(kind=R__)    :: z0, temp0, ext0, ssa0, chi0

    ! Open the file, figure out the type, and get the grid size
    open (unit=11, file=trim(adjustl(prp_obj%fname)), status='old')

    read (11, '(A1)') prp_obj%type

    if ((prp_obj%type .ne. 'E') .and. (prp_obj%type .ne. 'T')&
        &.and. (prp_obj%type .ne. 'P')) rewind (1)

    if ((prp_obj%type .ne. 'P') .and. (prp_obj%NlegSet .gt. 1)) then
        stop 'Error [_readPrpSize]: Nstokes>1 requires polarized phase function property file.'
    end if

    read (11, *) prp_obj%Nx, prp_obj%Ny, prp_obj%Nz
    prp_obj%maxNz = prp_obj%Nz
    prp_obj%maxNpts = prp_obj%Nx * prp_obj%Ny * prp_obj%Nz

    read (11, *) prp_obj%dx, prp_obj%dy, (z0, k=1, prp_obj%Nz)
    prp_obj%xDomain = prp_obj%Nx * prp_obj%dx
    prp_obj%yDomain = prp_obj%Ny * prp_obj%dy


    ! Property file type E is for extinction only format
    if (prp_obj%type .eq. 'E') then

        prp_obj%Npha = 1
        read (11, *) (temp0, iz=1, prp_obj%Nz)
        read (11, *) ssa0, Nleg

        prp_obj%Nleg = max(prp_obj%Nleg, Nleg)
        prp_obj%maxNleg = prp_obj%Nleg
        prp_obj%maxNptsLeg = prp_obj%maxNleg + 1

    ! Property file type T is for tabulated phase function format
    else if (prp_obj%type .eq. 'T') then

        prp_obj%maxNleg = prp_obj%Nleg

        read (11, *) prp_obj%Npha
        do i = 1, prp_obj%Npha
            read (11, *) Nleg, (chi0, l=1, Nleg)
            prp_obj%maxNleg = max(Nleg, prp_obj%maxNleg)
        end do

        prp_obj%maxNptsLeg = prp_obj%Npha * (prp_obj%maxNleg + 1)

    ! Property file type P is for polarized tabulated phase functi
    else if (prp_obj%type .eq. 'P') then

        prp_obj%maxNleg = prp_obj%Nleg
        read (11, *) prp_obj%Npha
        do i = 1, prp_obj%Npha
            do j = 1, 6
                read (11, *) k, Nleg, (chi0, l=0, Nleg)
                prp_obj%maxNleg = max(Nleg, prp_obj%maxNleg)
            end do
        end do
        prp_obj%maxNptsLeg = prp_obj%NlegSet * prp_obj%Npha * (prp_obj%maxNleg + 1)

    else

        prp_obj%Npha = 0

        ! Standard property file has everything variable
        prp_obj%maxNleg = prp_obj%Nleg
        do while (.true.)
            read (11, *, end=190) ix, iy, iz, temp0, ext0, ssa0, Nleg,&
                &(chi0, l=1, Nleg)
            prp_obj%maxNleg = max(Nleg, prp_obj%maxNleg)
        end do

190     continue

        prp_obj%maxNptsLeg = prp_obj%Nx * prp_obj%Ny * prp_obj%Nz * (prp_obj%maxNleg + 1)

    end if

    close (unit=11)

    return

end subroutine shd_io__readPrpSize

subroutine shd_io__readPrp(isDeltaM, &
                           prp_obj, temp_Prp, ext_Prp, ssa_Prp, leg_Prp, ipha_Prp)

    ! Reads the medium properties from the file into the property arrays.
    ! There are four formats: the standard format specifies the extinction,
    ! scattering albedo, temperature, and Legendre phase function expansion
    ! at every grid point; the extinction only format only specifies
    ! extinction at every grid point, while the other quantities are given
    ! in the header part; and the tabulated and polarized phase function
    ! format specifies all properties except the phase function for each
    ! grid point, while the phase functions are tabulated in the header.
    ! The polarized format specifies the six unique elements of the phase
    ! (Mueller) matrix in the header and includes the l=0 term.
    ! If the format is extinction only or tabulated phase function then
    ! prp_obj%Npha is the number of phase functions, otherwise prp_obj%Npha=0.
    ! If doing delta-M then prp_obj%Nleg is the minimum number of Legendre terms on
    ! input and the actual maximum number of terms on output, except that
    ! it may not exceed prp_obj%maxNleg; otherwise prp_obj%Nleg is the number of Legendre
    ! terms to be used (regardless of what is in property file).
    ! See the overall program documentation for the four file formats.

    logical, intent(in) :: isDeltaM

    type(shd_type_prp), intent(inout) :: prp_obj
    real(kind=R__), intent(out) :: temp_Prp(prp_obj%maxNpts), &
                                   ext_Prp(prp_obj%maxNpts), &
                                   ssa_Prp(prp_obj%maxNpts), &
                                   leg_Prp(prp_obj%NlegSet, &
                                           0:prp_obj%maxNleg, *)
    integer(kind=2), intent(out) :: ipha_Prp(prp_obj%maxNpts)

    !integer i, ix, iy, iz, j, js, k, kl, l, il, Nleg, ipha0, N
    integer i, ix, iy, iz, j, js, k, l, Nleg, ipha0, N
    real :: temp0, ext0, ssa0, chi0

    character(len=charMaxLen_) :: fnameBin
    character(len=charMaxLen_) :: lineTmp
    character(len=charMaxLen_) :: lineDim
    character(len=20) :: dataPostfix
    integer :: Nxx, Nyy, Nzz, Npp

    ! Open the file, figure out the type, and get the grid size
    open (unit=1, file=trim(adjustl(prp_obj%fname)), status='old')

    read (1, '(A1)') prp_obj%type

    if (prp_obj%type .ne. 'E' .and. prp_obj%type .ne. 'T'&
        &.and. prp_obj%type .ne. 'P') rewind (1)

    if (prp_obj%type .ne. 'P' .and. prp_obj%NlegSet .gt. 1) then
        write (*, '(A,A)') 'Error [_readPrp]: Nstokes>1 requires ',&
            &'Polarized phase function property file.'
        stop
    end if

    read (1, *) prp_obj%Nx, prp_obj%Ny, prp_obj%Nz

    if (prp_obj%Nz .gt. prp_obj%maxNz) stop 'Error [_readPrp]: prp_obj%maxNz exceeded.'

    read (1, *) prp_obj%dx, prp_obj%dy, (prp_obj%z(k), k=1, prp_obj%Nz)

    ! Property file type E is for extinction only format
    if (prp_obj%type .eq. 'E') then

        prp_obj%Npha = 1

        read (1, *) (temp0, iz=1, prp_obj%Nz)
        read (1, *) ssa0, Nleg

        if (isDeltaM) then
            prp_obj%Nleg = max(prp_obj%Nleg, Nleg)
            prp_obj%Nleg = min(prp_obj%Nleg, prp_obj%maxNleg)
        end if

        if (prp_obj%Nx * prp_obj%Ny * prp_obj%Nz .gt. prp_obj%maxNpts) then
            stop 'Error [_readPrp]: prp_obj%maxNpts exceeded.'
        end if

        if (prp_obj%Nleg + 1 .gt. prp_obj%maxNptsLeg)&
            &stop 'Error [_readPrp]: prp_obj%maxNptsLeg exceeded.'

        rewind (1)
        read (1, *)
        read (1, *) prp_obj%Nx, prp_obj%Ny, prp_obj%Nz
        read (1, *) prp_obj%dx, prp_obj%dy, (prp_obj%z(k), k=1, prp_obj%Nz)
        read (1, *) (temp_Prp(iz), iz=1, prp_obj%Nz)

        !read (1, *) ssa0, Nleg, (leg_Prp(l + 1), l=1, min(Nleg, prp_obj%Nleg)),&
        !    &(chi0, l=prp_obj%Nleg + 1, Nleg)
        read (1, *) ssa0, Nleg, (leg_Prp(1, l, 1), l=1, min(Nleg, prp_obj%Nleg)),&
            &(chi0, l=prp_obj%Nleg + 1, Nleg)

        !leg_Prp(1) = 1.0
        leg_Prp(1, 0, 1) = 1.0

        !prp_obj%maxAsym = leg_Prp(1 + 1)/3.0
        prp_obj%maxAsym = leg_Prp(1, 1, 1)/3.0
        do l = Nleg + 1, prp_obj%Nleg
            !leg_Prp(l) = 0.0
            leg_Prp(1, l, 1) = 0.0
        end do

        iy = 1
        do while (.true.)
            if (prp_obj%Ny .eq. 1) then
                read (1, *, end=290) ix, iz, ext0
            else
                read (1, *, end=290) ix, iy, iz, ext0
            end if
            if (ix .ge. 1 .and. ix .le. prp_obj%Nx .and.&
               &iy .ge. 1 .and. iy .le. prp_obj%Ny .and.&
               &iz .ge. 1 .and. iz .le. prp_obj%Nz) then
                k = iz + prp_obj%Nz*(iy - 1) + prp_obj%Nz*prp_obj%Ny*(ix - 1)
                temp_Prp(k) = temp_Prp(iz)
                ext_Prp(k)  = ext0
                ssa_Prp(k)  = ssa0
                ipha_Prp(k) = 1
            end if
        end do

    ! Property file type T or P is for tabulated phase function format
    else if (prp_obj%type .eq. 'T' .or. prp_obj%type .eq. 'P') then

        read (1, *) prp_obj%Npha

        ! If delta-M then find the largest number of Legendre terms
        if (isDeltaM) then

            do i = 1, prp_obj%Npha
                if (prp_obj%type .eq. 'T') then
                    read (1, *) Nleg, (chi0, l=1, Nleg)
                else
                    do j = 1, 6
                        read (1, *) js, Nleg, (chi0, l=0, Nleg)
                    end do
                end if
                prp_obj%Nleg = max(Nleg, prp_obj%Nleg)
            end do

            if (prp_obj%Nleg .gt. prp_obj%maxNleg) then
                if (verbose_ .ge. 4) write (*, '(A,I4,A,I4,A)')&
                    &'Warning [_readPrp]: Truncating phase function from ',&
                    &prp_obj%Nleg, ' to ', prp_obj%maxNleg, ' Legendre terms.'
            end if
            prp_obj%Nleg = min(prp_obj%Nleg, prp_obj%maxNleg)

        end if

        N = prp_obj%Nx * prp_obj%Ny * prp_obj%Nz

        if (N .gt. prp_obj%maxNpts) stop 'Error [_readPrp]: prp_obj%maxNpts exceeded.'

        if (prp_obj%Npha*(prp_obj%Nleg+1)*prp_obj%NlegSet .gt. prp_obj%maxNptsLeg)&
            &stop 'Error [_readPrp]: prp_obj%maxNptsLeg exceeded.'

        rewind (1)
        read (1, *)
        read (1, *) prp_obj%Nx, prp_obj%Ny, prp_obj%Nz
        read (1, *) prp_obj%dx, prp_obj%dy, (prp_obj%z(k), k=1, prp_obj%Nz)
        read (1, *) prp_obj%Npha

        prp_obj%maxAsym = -1.0
        do i = 1, prp_obj%Npha

            !il = (prp_obj%Nleg + 1)*(i - 1)

            if (prp_obj%type .eq. 'T') then

                !leg_Prp(il + 1) = 1.0
                leg_Prp(1, 0, i) = 1.0

                !read (1, *) Nleg, (leg_Prp(il + l + 1), l=1, min(Nleg, prp_obj%Nleg)),&
                !    &(chi0, l=prp_obj%Nleg + 1, Nleg)
                read (1, *) Nleg, (leg_Prp(1, l, i), l=1, min(Nleg, prp_obj%Nleg)),&
                    &(chi0, l=prp_obj%Nleg + 1, Nleg)

                do l = Nleg + 1, prp_obj%Nleg
                    !leg_Prp(l + il + 1) = 0.0
                    leg_Prp(1, l, i) = 0.0
                end do
            else
                do j = 1, 6
                    if (j .le. prp_obj%NlegSet) then
                        !read (1, *) js, Nleg,&
                        !    &(leg_Prp((il + l)*prp_obj%NlegSet + j), l=0, min(Nleg, prp_obj%Nleg)),&
                        !    &(chi0, l=prp_obj%Nleg + 1, Nleg)
                        read (1, *) js, Nleg,&
                            &(leg_Prp(j, l, i), l=0, min(Nleg, prp_obj%Nleg)),&
                            &(chi0, l=prp_obj%Nleg + 1, Nleg)
                        do l = Nleg + 1, prp_obj%Nleg
                            !leg_Prp((il + l)*prp_obj%NlegSet + j) = 0.0
                            leg_Prp(j, l, i) = 0.0
                        end do
                    else
                        read (1, *) js, Nleg, (chi0, l=0, Nleg)
                    end if
                end do
            end if
            !prp_obj%maxAsym = max(prp_obj%maxAsym, leg_Prp((il + 1)*prp_obj%NlegSet + 1)/3.0)
            prp_obj%maxAsym = max(prp_obj%maxAsym, leg_Prp(1, 1, i)/3.0)
        end do

        read (1, '(A)') lineTmp ! skip one line

        if (lineTmp(1:1) .eq. '!') then

            read (1, '(A)') dataPostfix
            read (1, '(A)') lineDim

            read (lineDim( 3:12), '(I10)') Nxx
            read (lineDim(14:23), '(I10)') Nyy
            read (lineDim(25:34), '(I10)') Nzz
            read (lineDim(36:45), '(I10)') Npp

            fnameBin = trim(adjustl(prp_obj%fname))//trim(adjustl(dataPostfix(2:)))

            open (UNIT=88, FILE=trim(adjustl(fnameBin)), STATUS='OLD',&
                &ACCESS='STREAM', FORM='UNFORMATTED')

            ! temp_Prp: real(4)
            read (88) (((temp_Prp(iz+prp_obj%Nz*(iy-1)+prp_obj%Nz*prp_obj%Ny*(ix-1)), iz=1,Nzz), iy=1,Nyy), ix=1,Nxx)
            ! ext_Prp: real(4)
            read (88) ((( ext_Prp(iz+prp_obj%Nz*(iy-1)+prp_obj%Nz*prp_obj%Ny*(ix-1)), iz=1,Nzz), iy=1,Nyy), ix=1,Nxx)
            ! ssa_Prp: real(4)
            read (88) ((( ssa_Prp(iz+prp_obj%Nz*(iy-1)+prp_obj%Nz*prp_obj%Ny*(ix-1)), iz=1,Nzz), iy=1,Nyy), ix=1,Nxx)
            ! ipha_Prp: integer(2)
            read (88) (((ipha_Prp(iz+prp_obj%Nz*(iy-1)+prp_obj%Nz*prp_obj%Ny*(ix-1)), iz=1,Nzz), iy=1,Nyy), ix=1,Nxx)

            close (88)

        else

            backspace(1)

            ! original
            !╭────────────────────────────────────────────────╮!
            do while (.true.)
                read (1, *, end=290) ix, iy, iz, temp0, ext0, ssa0, ipha0
                if (ix .ge. 1 .and. ix .le. prp_obj%Nx .and.&
                   &iy .ge. 1 .and. iy .le. prp_obj%Ny .and.&
                   &iz .ge. 1 .and. iz .le. prp_obj%Nz) then
                    k = iz + prp_obj%Nz*(iy - 1) + prp_obj%Nz*prp_obj%Ny*(ix - 1)
                    temp_Prp(k) = temp0
                    ext_Prp(k)  = ext0
                    ssa_Prp(k)  = ssa0
                    ipha_Prp(k) = ipha0
                end if
            end do
            !╰────────────────────────────────────────────────╯!

        end if

    else

        ! Standard property file has everything variable
        !   If delta-M then find the largest number of Legendre terms
        prp_obj%Npha = 0

        if (isDeltaM) then
            do while (.true.)
                read (1, *, end=190) ix, iy, iz, temp0, ext0, ssa0, Nleg,&
                    &(chi0, l=1, Nleg)
                prp_obj%Nleg = max(Nleg, prp_obj%Nleg)
            end do

190         continue

            prp_obj%Nleg = min(prp_obj%Nleg, prp_obj%maxNleg)

            rewind (1)

            read (1, *) prp_obj%Nx, prp_obj%Ny, prp_obj%Nz
            read (1, *) prp_obj%dx, prp_obj%dy, (prp_obj%z(k), k=1, prp_obj%Nz)
        end if

        N = prp_obj%Nx * prp_obj%Ny * prp_obj%Nz

        if (N .gt. prp_obj%maxNpts) stop 'Error [_readPrp]: prp_obj%maxNpts exceeded.'

        if (N*(prp_obj%Nleg + 1) .gt. prp_obj%maxNptsLeg)&
            &stop 'Error [_readPrp]: prp_obj%maxNptsLeg exceeded.'

        prp_obj%maxAsym = -1.0
        do while (.true.)

            read (1, *, end=290) ix, iy, iz, temp0, ext0, ssa0, Nleg

            if  (ix .ge. 1 .and. ix .le. prp_obj%Nx .and.&
                &iy .ge. 1 .and. iy .le. prp_obj%Ny.and.&
                &iz .ge. 1 .and. iz .le. prp_obj%Nz) then

                k = iz + prp_obj%Nz*(iy - 1) + prp_obj%Nz*prp_obj%Ny*(ix - 1)

                !kl = (prp_obj%Nleg + 1)*(k - 1)

                backspace (1)

                !read (1, *, end=290) ix, iy, iz, temp_Prp(k), ext_Prp(k),&
                !    &ssa_Prp(k), Nleg, (leg_Prp(l + kl), l=1, min(Nleg, prp_obj%Nleg)),&
                !    &(chi0, l=prp_obj%Nleg + 1, Nleg)
                read (1, *, end=290) ix, iy, iz, temp_Prp(k), ext_Prp(k),&
                    &ssa_Prp(k), Nleg, (leg_Prp(1, l, k), l=1, min(Nleg, prp_obj%Nleg)),&
                    &(chi0, l=prp_obj%Nleg + 1, Nleg)

                !prp_obj%maxAsym = max(prp_obj%maxAsym, leg_Prp(kl + 1 + 1)/3.0)
                !leg_Prp(1 + kl) = 1.0

                prp_obj%maxAsym = max(prp_obj%maxAsym, leg_Prp(1, 1, k)/3.0)
                leg_Prp(1, 0, k) = 1.0

                do l = Nleg + 1, prp_obj%Nleg
                    !leg_Prp(l + kl) = 0.0
                    leg_Prp(1, l, k) = 0.0
                end do

            end if
        end do

    end if

290 continue
    close (1)

    return

end subroutine shd_io__readPrp

subroutine shd_io__checkPrp(prp_obj, temp_Prp, ext_Prp, ssa_Prp, leg_Prp, ipha_Prp)

    ! Checks the input medium property values for correctness.
    ! Warning if temperature below 150 K or above 350 K.
    ! Stop if negative extinction, single scattering albedo below 0 or
    ! above 1, asymmetry parameter below -1 or above 1, or IPHASE
    ! out of range (of tabulated phase functions).
    ! Warning if Chi1=1, for people specifying Chi0.
    ! Warning if optical depth across a grid cell exceeds 2.

    type(shd_type_prp), intent(inout) :: prp_obj

    real(kind=R__), intent(in)  :: temp_Prp(prp_obj%maxNpts), &
                                   ext_Prp(prp_obj%maxNpts), &
                                   ssa_Prp(prp_obj%maxNpts), &
                                   leg_Prp(prp_obj%NlegSet, &
                                           0:prp_obj%maxNleg, *)
    integer(kind=2), intent(in) :: ipha_Prp(prp_obj%maxNpts)

    integer(kind=I__) :: i, iz, N, Nthick
    real(kind=R__)    :: tmin, tmax, tau, maxTau, dz
    logical :: legWarn

    ! prp_obj%Nx and prp_obj%Ny must be 1 or greater, prp_obj%Nz must be 2 or greater
    if (prp_obj%Nx .lt. 1 .or. prp_obj%Ny .lt. 1 .or. prp_obj%Nz .lt. 2) then
        write (*, *) 'Error [_checkPrp]: Bad grid size - ',&
            &prp_obj%Nx, prp_obj%Ny, prp_obj%Nz
        stop
    end if

    if (prp_obj%dx .le. 0.0 .or. prp_obj%dy .le. 0.0) then
        write (*, *) 'Error [_checkPrp]: Bad grid spacing - ',&
            &'prp_obj%dx,prp_obj%dy: ', prp_obj%dx, prp_obj%dy
        stop
    end if

    ! The vertical levels must be increasing
    do iz = 1, prp_obj%Nz - 1
        if (prp_obj%z(iz + 1) .le. prp_obj%z(iz)) then
            write (*, *) 'Error [_checkPrp]: prp_obj%z must increase : ',&
                &iz, prp_obj%z(iz), prp_obj%z(iz + 1)
            stop
        end if
    end do

    ! Loop over the property grid points, checking parameters
    !   Warn about temperatures outside of Earth atmosphere temperat
    N = prp_obj%Nx*prp_obj%Ny*prp_obj%Nz
    tmin = 1000.0
    tmax = -1000.0
    do i = 1, N
        tmin = min(tmin, temp_Prp(i))
        tmax = max(tmax, temp_Prp(i))
    end do

    if (tmin .lt. 150.0 .or. tmax .gt. 350.0) then
        if (verbose_ .ge. 4) write (*, *) 'Warning [_checkPrp]: ',&
            &'temperature min and max (Kelvin): ', tmin, tmax
        write (*, *)
    end if

    do i = 1, N
        if (ext_Prp(i) .lt. 0.0) then
            write (*, *) 'Error [_checkPrp]: Negative extinction - ',&
                &ext_Prp(i), ' at grid point ', i
            stop
        end if
        if (ssa_Prp(i) .lt. 0.0 .or. ssa_Prp(i) .gt. 1.0) then
            write (*, *) 'Error [_checkPrp]: Single scattering ',&
                &'albedo out of range - ', ssa_Prp(i), ' at grid point ', i
            stop
        end if
    end do

    if (prp_obj%Npha .gt. 0) N = prp_obj%Npha
    legWarn = .false.
    do i = 1, N
        if (leg_Prp(1, 1, i) .lt. -3.0 .or. leg_Prp(1, 1, i) .gt. 3.0) then
            write (*, *) 'Error [_checkPrp]: Asymmetry parameter ',&
                &'out of range - ', leg_Prp(1, 1, i)/3.0, ' at grid point ', i
            stop
        end if
        if (leg_Prp(1, 1, i) .eq. 1.0) legWarn = .true.
    end do

    if (legWarn) then
        if (verbose_ .ge. 4) write (*, *) 'Warning [_checkPrp]: ',&
            &'an asymmetry parameter is 1/3 (Chi1=1); ',&
            &'  Chi0 should not be specified for Chi1.'
        write (*, *)
    end if

    if (prp_obj%Npha .gt. 0) then
        N = prp_obj%Nx*prp_obj%Ny*prp_obj%Nz
        do i = 1, N
            if (ipha_Prp(i) .lt. 1 .or. ipha_Prp(i) .gt. prp_obj%Npha) then
                write (*, *) 'Error [_checkPrp]: Phase function pointer',&
                    &' out of range - ', ipha_Prp(i), ' at grid point ', i
                stop
            end if
        end do
    end if

    N = prp_obj%Nx*prp_obj%Ny*prp_obj%Nz
    Nthick = 0
    maxTau = 0.0
    do iz = 1, prp_obj%Nz - 1
        dz = prp_obj%z(iz + 1) - prp_obj%z(iz)
        do i = 0, N - 1, prp_obj%Nz
            tau = 0.5*(ext_Prp(i + iz) + ext_Prp(i + iz + 1))*dz
            if (tau .gt. 2.0) then
                Nthick = Nthick + 1
                maxTau = max(maxTau, tau)
            end if
        end do
    end do

    if (Nthick .gt. 0) then
        if (verbose_ .ge. 5) write (*, '(A)') 'Message [_checkPrp]:'
        if (verbose_ .ge. 5) write (*, '(A,I0)')   '  Npoints (optical depth > 2): ', Nthick
        if (verbose_ .ge. 5) write (*, '(A,F0.2)') '       Max cell optical depth: ', maxTau
    end if

    return

end subroutine shd_io__checkPrp

subroutine shd_io__readSfcSize(sfc_obj)

    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Description:
    !   Gets the size of the surface property inputs (temperature is
    !   always included).
    !
    ! Input:
    !   sfc_obj:
    !     %fname: character, file path for the surface file (contains BRDF coefficients)
    !
    ! Output:
    !   sfc_obj:
    !     %Nx: integer, number of grid points in X
    !     %Ny: integer, number of grid points in Y
    !     %dx: real, X grid resolution
    !     %dy: real, Y grid resolution
    !     %maxNpts: integer, maximum surface points (including extra one line of grids for both X and Y)
    !     %maxNpar: integer, maximum surface parameters
    !╰────────────────────────────────────────────────────────────────────────────╯!

    type(shd_type_sfc), intent(inout) :: sfc_obj

    ! Open the file, figure out the type, and get the grid size
    sfc_obj%type(1:1) = 'V'

    open (unit=21, file=trim(adjustl(sfc_obj%fname)), status='old')

    read (21, '(A1)') sfc_obj%type(2:2)
    read (21, *) sfc_obj%Nx, sfc_obj%Ny, sfc_obj%dx, sfc_obj%dy

    sfc_obj%maxNpts = (sfc_obj%Nx + 1)*(sfc_obj%Ny + 1)

    ! 0
    if ((sfc_obj%type .eq. 'VL') .or. (sfc_obj%type .eq. 'Vl')) then
        ! Surface file type L is for variable Lambertian surface
        sfc_obj%maxNpar = 2
        sfc_obj%isLambertian = .true.
    ! 1
    else if ((sfc_obj%type .eq. 'VO') .or. (sfc_obj%type .eq. 'Vo')) then
        ! Surface file type O is for Ocean BRDF surface
        sfc_obj%maxNpar = 3
    ! 2
    else if ((sfc_obj%type .eq. 'VR') .or. (sfc_obj%type .eq. 'Vr')) then
        ! Surface file type R is for RPV surface
        sfc_obj%maxNpar = 4
    ! 3
    else if ((sfc_obj%type .eq. 'VT') .or. (sfc_obj%type .eq. 'Vt')) then
        ! Surface file type T is for LSRT surface
        sfc_obj%maxNpar = 4
    ! 4
    else if ((sfc_obj%type .eq. 'VJ') .or. (sfc_obj%type .eq. 'Vj')) then
        ! Surface file type J is for LSRT-Jiao surface
        sfc_obj%maxNpar = 6
    ! -1
    else if ((sfc_obj%type .eq. 'VW') .or. (sfc_obj%type .eq. 'Vw')) then
        ! Surface file type W is for Wave-Fresnel surface
        sfc_obj%maxNpar = 4
    ! -2
    else if ((sfc_obj%type .eq. 'VD') .or. (sfc_obj%type .eq. 'Vd')) then
        ! Surface file type D is for Diner polarized RPV surface:
        sfc_obj%maxNpar = 6
    else if ((sfc_obj%type .eq. 'VX') .or. (sfc_obj%type .eq. 'Vx')) then
        ! Surface file type X is for Mixed-surface type:
        sfc_obj%maxNpar = 8
    else
        stop 'Error [_readSfcSize]: Unknown BRDF type.'
    end if

    close (unit=21)

    return

end subroutine shd_io__readSfcSize

subroutine shd_io__readSfc(sfc_obj, brdf_Sfc)

    ! Reads the surface properties from the file into the surface arrays.
    ! There are currently six surface types (Lambertian, Wave-Fresnel,
    ! Ocean, RPV, and Diner polarized RPV), though it is simple to add more
    ! In all cases the surface properties are specified on a regular, evenly.
    ! spaced grid of size sfc_obj%Nx by sfc_obj%Ny and spacing sfc_obj%dx by sfc_obj%dy.
    ! The properties are returned on a sfc_obj%Nx+1 by sfc_obj%Ny+1 grid that is
    ! periodic (NX+1 same as 1, etc).  The surface properties at each
    ! grid point are returned in brdf_Sfc, and the first parameter must
    ! be temperature.  The number of parameters for each grid point
    ! (including temperature) are returned in sfc_obj%Npar.
    ! Also returned are the domain average ground albedo and temperature
    ! (sfc_obj%avgAlb, sfc_obj%avgTemp) for use with initialization (hence they can be
    ! highly approximate).
    !        Type        Parameters
    !   L  Lambertian    albedo
    !   O  Ocean         [albedo], Wind speed (m/s), Pigmentation (mg/m^3) {albedo is added for easily retrieve average albedo}
    !   R  RPV           rho0, k, Theta
    !   T  LSRT          f_iso, f_geo, f_vol
    !   J  LSRT-Jiao     f_iso, f_geo, f_vol, f_jiao, alpha
    !   W  Wave-Fresnel  [albedo], Real, Imaginary index of refraction, wind speed (m/s) {albedo is added for easily retrieve
    !                    average albedo}
    !   D  Diner et al   a, k, b, zeta, sigma
    !   M  Mixed

    type(shd_type_sfc), intent(inout) :: sfc_obj
    real(kind=R__), intent(out) :: brdf_Sfc(sfc_obj%maxNpar, sfc_obj%maxNpts)

    real :: temp_sfc, alb_sfc

    integer :: N, i, i0, ix, iy, j
    real :: alb0, temp0, MRE, MIM, RHO0, KR, THETA, WSPD, PCL
    real :: param1, param2, param3, param4, param5, param6, param7, param8
    real F_ISO, F_VOL, F_GEO, F_JIAO
    real A, B, ZETA, SIGMA, ALPHA

    character(len=charMaxLen_) :: fnameBin
    character(len=charMaxLen_) :: lineTmp
    character(len=charMaxLen_) :: lineDim
    character(len=20) :: dataPostfix
    integer :: Nxx, Nyy, Npp, ipar, idata, Ndata

    ! Surface type is variable something
    sfc_obj%type(1:1) = 'V'

    ! Open the file, figure out the type, and get the grid size
    open (unit=1, file=trim(adjustl(sfc_obj%fname)), status='old')

    read (1, '(A1)') sfc_obj%type(2:2)
    read (1, *) sfc_obj%Nx, sfc_obj%Ny, sfc_obj%dx, sfc_obj%dy

    if ((sfc_obj%Nx + 1)*(sfc_obj%Ny + 1) .gt. sfc_obj%maxNpts) then
        write (*, *) 'Error [_readSfc]: sfc_obj%maxNpts exceeded.', (sfc_obj%Nx + 1)*(sfc_obj%Ny + 1), sfc_obj%maxNpts
        stop
    end if

    temp_sfc = 0.0
    alb_sfc  = 0.0
    N = 0

    ! Surface file type L is for variable Lambertian surface
    if ((sfc_obj%type .eq. 'VL') .or. (sfc_obj%type .eq. 'Vl')) then

        sfc_obj%Npar = 2
        sfc_obj%isLambertian = .true.

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = alb0

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! albedo
                    if (alb0 .lt. 0.0 .or. alb0 .gt. 1.0)&
                        &stop 'Error [_readSfc]: Illegal surface albedo.'
                    if (temp0 .lt. 150. .or. temp0 .gt. 350.) then
                        if (verbose_ .ge. 4) write (*, *) 'Warning [_readSfc]: ',&
                            &' surface temperature: ', ix, iy, temp0
                    end if

                    temp_sfc = temp_sfc + param1
                    alb_sfc  = alb_sfc  + param2

                    N = N + 1

                end if
            end do
        end if

    ! Surface file type O is for Ocean surface BRDF
    !   Norm Loeb's modification of 6S ocean reflectance module.
    else if ((sfc_obj%type .eq. 'VO') .or. (sfc_obj%type .eq. 'Vo')) then

        sfc_obj%Npar = 3

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2, param3
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = WSPD
                    !brdf_Sfc(i + 3) = PCL

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! 10m wind speed, units: m/s
                    brdf_Sfc(3, i) = param3 ! pigment concentration, units: mg/m^3

                    temp_sfc = temp_sfc + param1

                    N = N + 1
                end if
            end do
        end if

    ! Surface file type R is for RPV surface:
    !   (Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphe
    !   Reflectance (CSAR) Model. 2. Semiempirical Surface Model
    !   Usable With NOAA Advanced Very High Resolution Radiometer
    !   Data, J. Geophys. Res., 98, 20791-20801.)
    else if ((sfc_obj%type .eq. 'VR') .or. (sfc_obj%type .eq. 'Vr')) then

        sfc_obj%Npar = 4

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2, param3, param4
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = RHO0
                    !brdf_Sfc(i + 3) = KR
                    !brdf_Sfc(i + 4) = THETA

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! rho0
                    brdf_Sfc(3, i) = param3 ! kr
                    brdf_Sfc(4, i) = param4 ! theta

                    temp_sfc = temp_sfc + param1
                    alb_sfc  = alb_sfc  + param2

                    N = N + 1

                end if
            end do
        end if

    ! Surface file type T is for LSRT surface:
    else if ((sfc_obj%type .eq. 'VT') .or. (sfc_obj%type .eq. 'Vt')) then

        sfc_obj%Npar = 4

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2, param3, param4
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = F_ISO
                    !brdf_Sfc(i + 3) = F_GEO
                    !brdf_Sfc(i + 4) = F_VOL

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! f_iso
                    brdf_Sfc(3, i) = param3 ! f_geo
                    brdf_Sfc(4, i) = param4 ! f_vol

                    temp_sfc = temp_sfc + param1
                    alb_sfc  = alb_sfc  + param2

                    N = N + 1

                end if
            end do
        end if

    ! Surface file type J is for LSRT-Jiao surface:
    else if ((sfc_obj%type .eq. 'VJ') .or. (sfc_obj%type .eq. 'Vj')) then

        sfc_obj%Npar = 6

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2, param3, param4, param5, param6
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = F_ISO
                    !brdf_Sfc(i + 3) = F_GEO
                    !brdf_Sfc(i + 4) = F_VOL
                    !brdf_Sfc(i + 5) = F_JIAO
                    !brdf_Sfc(i + 6) = ALPHA

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! f_iso
                    brdf_Sfc(3, i) = param3 ! f_geo
                    brdf_Sfc(4, i) = param4 ! f_vol
                    brdf_Sfc(5, i) = param5 ! f_jiao
                    brdf_Sfc(6, i) = param6 ! alpha

                    temp_sfc = temp_sfc + param1
                    alb_sfc  = alb_sfc  + param2

                    N = N + 1

                end if
            end do
        end if

    ! Surface file type W is for Wave-Fresnel surface
    else if ((sfc_obj%type .eq. 'VW') .or. (sfc_obj%type .eq. 'Vw')) then

        sfc_obj%Npar = 4

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                 ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2, param3, param4
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = MRE
                    !brdf_Sfc(i + 3) = MIM
                    !brdf_Sfc(i + 4) = WSPD

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! M_real
                    brdf_Sfc(3, i) = param3 ! M_img
                    brdf_Sfc(4, i) = param4 ! 10m wind speed

                    temp_sfc = temp_sfc + param1

                    N = N + 1

                end if
            end do
        end if

    ! Surface file type D is for Diner et al surface reflection:
    ! (Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
    ! V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
    ! Exploration of a polarized surface bidirectional reflectance model
    ! using the ground-based multiangle spectropolarimetric imager.
    ! Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.)
    else if ((sfc_obj%type .eq. 'VD') .or. (sfc_obj%type .eq. 'Vd')) then

        sfc_obj%Npar = 6

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                 ! if not, go back one line
            do while (.true.)
                read (1, *, end=290) ix, iy, param1, param2, param3, param4, param5, param6
                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = temp0
                    !brdf_Sfc(i + 2) = A
                    !brdf_Sfc(i + 3) = KR
                    !brdf_Sfc(i + 4) = B
                    !brdf_Sfc(i + 5) = ZETA
                    !brdf_Sfc(i + 6) = SIGMA

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)
                    brdf_Sfc(1, i) = param1 ! temperature
                    brdf_Sfc(2, i) = param2 ! A
                    brdf_Sfc(3, i) = param3 ! Kr
                    brdf_Sfc(4, i) = param4 ! B
                    brdf_Sfc(5, i) = param5 ! Zeta
                    brdf_Sfc(6, i) = param6 ! Sigma

                    temp_sfc = temp_sfc + param1
                    alb_sfc  = alb_sfc  + param2

                    N = N + 1

                end if
            end do
        end if

    ! Surface file type X is for Mixed-type surface:
    else if ((sfc_obj%type .eq. 'VX') .or. (sfc_obj%type .eq. 'Vx')) then

        sfc_obj%Npar = 8

        if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
            write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.'
            stop
        end if

        read (1, '(A)') lineTmp
        if (lineTmp(1:1) .ne. '!') then ! check if binary input,
            backspace(1)                 ! if not, go back one line
            do while (.true.)
                read (1, *, end=290)    ix, iy, param1, param2, param3, param4, param5, param6, param7, param8
                ! for Lambertian   (VL):ix, iy, temp0,   alb0,    -999.0, -999.0, -999.0, -999.0, -999.0,  0.0
                ! for Ocean        (VO):ix, iy, temp0,   WSPD,   PCL,    -999.0, -999.0, -999.0, -999.0,  1.0
                ! for RPV          (VR):ix, iy, temp0,   RHO0,   KR,     THETA,  -999.0, -999.0, -999.0,  2.0
                ! for LSRT         (VT):ix, iy, temp0,   F_ISO,  F_GEO,  F_VOL,  -999.0, -999.0, -999.0,  3.0
                ! for LSRT-J       (VJ):ix, iy, temp0,   F_ISO,  F_GEO,  F_VOL,  F_JIAO, ALPHA,  -999.0,  4.0
                ! for Wave-Fresnel (VW):ix, iy, temp0,   MRE,    MIM,    WSPD,   -999.0, -999.0, -999.0, -1.0
                ! for Diner        (VD):ix, iy, temp0,   A,      KR,     B,      ZETA,   SIGMA,  -999.0, -2.0

                if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                    &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                    !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))
                    !brdf_Sfc(i + 1) = param1
                    !brdf_Sfc(i + 2) = param2
                    !brdf_Sfc(i + 3) = param3
                    !brdf_Sfc(i + 4) = param4
                    !brdf_Sfc(i + 5) = param5
                    !brdf_Sfc(i + 6) = param6
                    !brdf_Sfc(i + 7) = param7
                    !brdf_Sfc(i + 8) = param8 ! surface identification

                    i = iy + (sfc_obj%Ny + 1)*(ix - 1)
                    brdf_Sfc(1, i) = param1
                    brdf_Sfc(2, i) = param2
                    brdf_Sfc(3, i) = param3
                    brdf_Sfc(4, i) = param4
                    brdf_Sfc(5, i) = param5
                    brdf_Sfc(6, i) = param6
                    brdf_Sfc(7, i) = param7
                    brdf_Sfc(8, i) = param8 ! surface identification

                    temp_sfc = temp_sfc + param1

                    ! no albedo for polarized (W [-1]) and unpolarized (O [1]) Ocean parameterization
                    if ((abs(param8)-1.0) .gt. 0.01) then
                        alb_sfc = alb_sfc + param2
                    end if

                    N = N + 1
                end if
            end do
        end if

    else

        close (1)
        stop 'Error [_readSfc]: Unknown BRDF type <'//trim(adjustl(sfc_obj%type))//'>.'

    end if

    ! read binary
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (sfc_obj%Npar .gt. sfc_obj%maxNpar) then
        write (*, *) 'Error [_readSfc]: sfc_obj%maxNpar exceeded.', sfc_obj%Npar, sfc_obj%maxNpar
        stop
    end if

    read (1, '(A)') dataPostfix
    read (1, '(A)') lineDim

    read (lineDim( 3:12), '(I10)') Nxx
    read (lineDim(14:23), '(I10)') Nyy
    read (lineDim(25:34), '(I10)') Npp

    fnameBin = trim(adjustl(sfc_obj%fname))//trim(adjustl(dataPostfix(2:)))
    open (UNIT=88, FILE=trim(adjustl(fnameBin)), STATUS='OLD',&
        &ACCESS='STREAM', FORM='UNFORMATTED')
    Ndata = Nxx*Nyy
    read (88) ((brdf_Sfc(ipar, idata), ipar=1, Npp), idata=1, Ndata)
    close (88)

    do iy = 1, sfc_obj%Ny + 1
        do ix = 1, sfc_obj%Nx + 1

            if  (ix .ge. 1 .and. ix .le. sfc_obj%Nx .and.&
                &iy .ge. 1 .and. iy .le. sfc_obj%Ny) then

                !i = sfc_obj%Npar * (iy - 1 + (sfc_obj%Ny + 1)*(ix - 1))

                i = iy + (sfc_obj%Ny + 1)*(ix - 1)

                temp_sfc = temp_sfc + brdf_Sfc(1, i)

                if  (((sfc_obj%type .ne. 'VO') .and. (sfc_obj%type .ne. 'Vo')) .and.&
                    &((sfc_obj%type .ne. 'VW') .and. (sfc_obj%type .ne. 'Vw'))) then
                    if ((sfc_obj%type .eq. 'VX') .or. (sfc_obj%type .eq. 'Vx')) then
                        if ((abs(brdf_Sfc(8, i))-1.0) .gt. 0.01) then ! exclude [O] and [W]
                            alb_sfc = alb_sfc + brdf_Sfc(2, i)
                        end if
                    else
                        alb_sfc = alb_sfc + brdf_Sfc(2, i)
                    end if
                end if

                N = N + 1

            end if

        end do
    end do
    !╰────────────────────────────────────────────────────────────────────────────╯!

290 continue

    close (1)

    ! Copy the edge points to the opposite side for periodic boundar
    ! for x, copy bottom line to top line
    do ix = 1, sfc_obj%Nx
        !i0 = sfc_obj%Npar * (ix - 1 + (sfc_obj%Nx + 1)*(1 - 1))
        !i  = sfc_obj%Npar * (ix - 1 + (sfc_obj%Nx + 1)*(sfc_obj%Ny + 1 - 1))
        i0 = ix
        i  = ix + (sfc_obj%Nx + 1)*(sfc_obj%Ny)
        do j = 1, sfc_obj%Npar
            brdf_Sfc(j, i) = brdf_Sfc(j, i0)
        end do
    end do

    ! for y, copy left line to right line
    do iy = 1, sfc_obj%Ny + 1
        !i0 = sfc_obj%Npar * ((sfc_obj%Nx + 1)*(iy - 1))
        !i  = sfc_obj%Npar * (sfc_obj%Nx + (sfc_obj%Nx + 1)*(iy - 1))
        i0 = 1 + (sfc_obj%Nx + 1)*(iy - 1)
        i  = 1 + sfc_obj%Nx + (sfc_obj%Nx + 1)*(iy - 1)
        do j = 1, sfc_obj%Npar
            brdf_Sfc(j, i) = brdf_Sfc(j, i0)
        end do
    end do

    sfc_obj%avgTemp = temp_sfc/N
    sfc_obj%avgAlb  = alb_sfc/N
    return

end subroutine shd_io__readSfc

subroutine shd_io__readCkdSize(waveNum, ckd_obj)

    ! Find the size of the k-distribution arrays for this wavenumber band

    real(kind=R__), intent(in) :: waveNum(2)
    type(shd_type_ckd), intent(inout) :: ckd_obj

    integer(kind=I__) :: Nband, iband, ibandCheck, ibandSelect

    open (unit=1, file=trim(adjustl(ckd_obj%fname)), status='old')
    read (1, *)
    read (1, *) Nband
    read (1, *)

    ! Read band information until the right one is found
    ibandSelect = 0
    do iband = 1, Nband
        read (1, *) ibandCheck, ckd_obj%waveNum(1), ckd_obj%waveNum(2), ckd_obj%Ng
        if  (abs(ckd_obj%waveNum(1) - waveNum(1)) .lt. 1.0 .and.&
            &abs(ckd_obj%waveNum(2) - waveNum(2)) .lt. 1.0) then
            backspace (1)
            read (1, *) ibandSelect, ckd_obj%waveNum(1), ckd_obj%waveNum(2), ckd_obj%Ng
        end if
    end do

    if (ibandSelect .eq. 0) then
        write (*, *) 'Error [_readCkdSize]: Wavenumber range not found in CKD file : ',&
            &waveNum(1), waveNum(2), ckd_obj%fname
        stop
    end if

    read (1, *) ckd_obj%Nz
    close (1)

    return

end subroutine shd_io__readCkdSize

subroutine shd_io__readCkd(waveNum, &
    ckd_obj, solar_Ckd, wgt_Ckd, kabs_Ckd)

    ! Reads the information appropriate for one band from a
    ! correlated k-distribution file.   The wavenumber range (cm^-1)
    ! in the file must match the desired in waveNum.  The solar flux
    ! is returned in solar_Ckd, the number of "g"'s in ckd_obj%Ng, the weights
    ! or delta g's in wgt_Ckd, the number of levels in ckd_obj%Nz, the
    ! levels from the top down in ckd_obj%z, and the absorption coefficients
    ! for each g and level in kabs_Ckd.
    ! The array sizes ckd_obj%maxNg and ckd_obj%maxNz are checked.

    real(kind=R__), intent(inout) :: waveNum(2)
    type(shd_type_ckd), intent(inout) :: ckd_obj
    real(kind=R__), intent(inout) :: solar_Ckd(ckd_obj%Ng), &
                                     wgt_Ckd(ckd_obj%Ng), &
                                     kabs_Ckd(ckd_obj%Nz, ckd_obj%Ng)

    integer(kind=I__) :: Nband, iband, ibandCheck, ibandSelect
    integer(kind=I__) :: ig, iz, itmp

    open (unit=1, file=trim(adjustl(ckd_obj%fname)), status='old')

    read (1, *)
    read (1, *) Nband
    read (1, *)

    ! Read band information until the right one is found
    ibandSelect = 0
    do iband = 1, Nband
        read (1, *) ibandCheck, ckd_obj%waveNum(1), ckd_obj%waveNum(2), ckd_obj%Ng
        if  (abs(ckd_obj%waveNum(1) - waveNum(1)) .lt. 1.0 .and.&
            &abs(ckd_obj%waveNum(2) - waveNum(2)) .lt. 1.0) then
            backspace (1)
            if (ckd_obj%Ng .gt. ckd_obj%maxNg) stop 'Error [_readCkd]: ckd_obj%maxNg exceeded.'
            read (1, *) ibandSelect, waveNum(1), waveNum(2),&
                &ckd_obj%Ng, (solar_Ckd(ig), ig=1, ckd_obj%Ng), (wgt_Ckd(ig), ig=1, ckd_obj%Ng)
        end if
    end do

    if (ibandSelect .eq. 0) then
        write (*, *) 'Error [_readCkd]: Wavenumber range not found in CKD file : ',&
            &waveNum(1), waveNum(2), ckd_obj%fname
        stop
    end if

    ! Read the Z level heights
    read (1, *) ckd_obj%Nz
    if (ckd_obj%Nz .gt. ckd_obj%maxNz) stop 'Error [_readCkd]: ckd_obj%maxNz exceeded.'

    read (1, *)
    read (1, *)
    do iz = 1, ckd_obj%maxNz
        read (1, *) ckd_obj%z(iz)
    end do

    ! Skip over the irrelevant bands and read the absorption coeffic
    read (1, *)
    do itmp = 1, (ibandSelect - 1)*ckd_obj%maxNz
        read (1, *)
    end do

    do iz = 1, ckd_obj%maxNz
        read (1, *) iband, itmp, (kabs_Ckd(iz, ig), ig=1, ckd_obj%Ng)
        do ig = 1, ckd_obj%Ng
            if (kabs_Ckd(iz, ig) .lt. 0.0) then
                write (*, '(1X,A,A,I3,A,I3,A,I2)')&
                    &'Negative absorption in CKD file.  ',&
                    &'Band: ', iband, '   Level: ', iz, '   k: ', ig
            end if
        end do
    end do
    close (1)

    return

end subroutine shd_io__readCkd

subroutine shd_io__readSenSize(sen_obj)

    ! Reads the user-specified sensor data size information

    type(shd_type_sen), intent(inout) :: sen_obj

    character(len=charMaxLen_) :: lineTmp
    character(len=charMaxLen_) :: lineDim

    ! Open the file, figure out the type, and get the grid size
    open (unit=1, file=trim(adjustl(sen_obj%fname)), status='old')
    read (1, '(A)') lineTmp ! skip one line
    read (1, '(A)') lineTmp ! skip one line
    read (1, '(A)') lineTmp ! skip one line
    read (1, '(A)') lineDim

    read (lineDim( 3:12), '(I10)') sen_obj%Npar
    read (lineDim(14:23), '(I10)') sen_obj%Npts

    close (1)

    return

end subroutine shd_io__readSenSize

subroutine shd_io__readSen(sen_obj, view_Sen)

    ! Reads the user-specified sensor information

    type(shd_type_sen), intent(inout) :: sen_obj
    real(kind=R__),     intent(out)   :: view_Sen(sen_obj%Npar, sen_obj%Npts)

    integer :: ipar, ipts

    character(len=charMaxLen_) :: fnameBin
    character(len=charMaxLen_) :: lineTmp
    character(len=charMaxLen_) :: lineDim
    character(len=20) :: dataPostfix

    ! Open the file, figure out the type, and get the grid size
    open (unit=1, file=trim(adjustl(sen_obj%fname)), status='old')
    read (1, '(A)') lineTmp ! skip one line
    read (1, '(A)') lineTmp ! skip one line
    read (1, '(A)') dataPostfix
    read (1, '(A)') lineDim

    read (lineDim( 3:12), '(I10)') sen_obj%Npar
    read (lineDim(14:23), '(I10)') sen_obj%Npts

    fnameBin = trim(adjustl(sen_obj%fname))//trim(adjustl(dataPostfix(2:)))

    open (unit=88, file=trim(adjustl(fnameBin)), status='old',&
        &access='stream', form='unformatted')

    read (88) ((view_Sen(ipar, ipts), ipar=1, sen_obj%Npar), ipts=1, sen_obj%Npts)

    close (88)

    close (1)

    return

end subroutine shd_io__readSen

subroutine shd_io__restoreState(&
    INSAVEFILE, &
    INRADFLAG, NEWGRIDFLAG, &
    FLUXES, &
    rtm_obj, &
    grid_obj, xyz_Grid, &
    shdo_obj, iradI_Shdo, radI_Shdo, isrcJ_Shdo, srcJ_Shdo)

    ! Reads the binary save file containing a previous state of
    ! the SHDOM solution. This includes the cell tree structure,
    ! the spherical harmonic expansion of the source function,
    ! the low order terms of the radiance expansion, and the flux array.
    ! If the file name is 'NONE' then the INRADFLAG is set to false,
    ! otherwise it is true.  The input file may have any spherical
    ! harmonic truncation smaller than the current (shdo_obj%maxNl,shdo_obj%maxNm,shdo_obj%Nlm).
    ! The flag NEWGRIDFLAG is set to false if reading succeeds.

    type(shd_type_rtm),  intent(inout) :: rtm_obj
    type(shd_type_grid), intent(inout) :: grid_obj
    type(shd_type_shdo), intent(inout) :: shdo_obj

    real(kind=R__),    intent(out) :: xyz_Grid(3, *)
    integer(kind=I__), intent(out) :: iradI_Shdo(*)
    real(kind=R__),    intent(out) ::  radI_Shdo(shdo_obj%Nstokes, *)
    integer(kind=I__), intent(out) :: isrcJ_Shdo(*)
    real(kind=R__),    intent(out) ::  srcJ_Shdo(shdo_obj%Nstokes, *)


    character(*), intent(in) :: INSAVEFILE

    logical, intent(out) :: INRADFLAG, NEWGRIDFLAG
    real(kind=R__), intent(out) :: FLUXES(2, *)

    integer :: NSTI, NXI, NYI, NZI, MLI, MMI, NLMI
    integer :: i, j, NS, NR

    INRADFLAG = .false.
    NEWGRIDFLAG = .true.
    if (INSAVEFILE(1:4) .eq. 'NONE') return

    open (UNIT=3, FILE=trim(adjustl(INSAVEFILE)), STATUS='OLD',&
        &FORM='UNFORMATTED', ERR=900)

    read (3) NSTI, NXI, NYI, NZI, MLI, MMI, NLMI, grid_obj%Npts, grid_obj%Ncell

    if ((NXI .ne. rtm_obj%Nx) .or. (NYI .ne. rtm_obj%Ny) .or. (NZI .ne. rtm_obj%Nz)) then
        if (verbose_ .ge. 4) write (6, '(A, 3I6)') 'Warning [_restoreState]: wrong grid size:', NXI, NYI, NZI
        close (3)
        return
    end if

    if ((NSTI .ne. shdo_obj%Nstokes) .or. &
        (MLI .gt. shdo_obj%maxNl) .or. &
        (MMI .gt. shdo_obj%maxNm)) then
        if (verbose_ .ge. 4) write (6, '(A, 3I6)') 'Warning [_restoreState]: incompatible SH:', NSTI, MLI, MMI
        close (3)
        return
    end if

    read (3) (rtm_obj%x(i), i=1, rtm_obj%Nx+1)
    read (3) (rtm_obj%y(i), i=1, rtm_obj%Ny+1)
    read (3) (rtm_obj%z(i), i=1, rtm_obj%Nz)
    read (3) ((grid_obj%i8Grid(j, i), j=1, 8), i=1, grid_obj%Ncell)
    read (3) ((grid_obj%i6Cell(j, i), j=1, 6), i=1, grid_obj%Ncell)
    read (3) ((grid_obj%i2Tree(j, i), j=1, 2), i=1, grid_obj%Ncell)
    read (3) (grid_obj%flag(i), i=1, grid_obj%Ncell)

    read (3) ((xyz_Grid(j, i), j=1, 3), i=1, grid_obj%Npts)
    read (3) ((FLUXES(j, i), j=1, 2), i=1, grid_obj%Npts)

    read (3) (isrcJ_Shdo(i), i=1, grid_obj%Npts + 1)

    NS = isrcJ_Shdo(grid_obj%Npts + 1)
    read (3) ((srcJ_Shdo(j, i), j=1, shdo_obj%Nstokes), i=1, NS)

    NR = 4*grid_obj%Npts
    read (3) ((radI_Shdo(j, i), j=1, shdo_obj%Nstokes), i=1, NR)

    iradI_Shdo(1) = 0
    do i = 1, grid_obj%Npts
        iradI_Shdo(i+1) = iradI_Shdo(i) + shdo_obj%Nstokes*4
    end do

    close (3)

    if (verbose_ .ge. 5) write (6, '(A)') 'Message [_restoreState]: State restored from <'//trim(adjustl(INSAVEFILE))//'> ...'

    INRADFLAG = .true.
    NEWGRIDFLAG = .false.
    return

900 continue

    if (verbose_ .ge. 4) write (6, '(A)') 'Warning [_restoreState]: Failed to restore state file from <'//trim(adjustl(INSAVEFILE))//'>.'
    INRADFLAG = .false.
    NEWGRIDFLAG = .true.

    return

end subroutine shd_io__restoreState

end module shd_inp
