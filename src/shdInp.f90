module shd_inp

  use globals
  use shd_util

  implicit none
  private

  public :: shd_io__userInput
  public :: shd_io__nmlInput

  public :: shd_io__checkInputParam

  public :: shd_io__readPrpSize
  public :: shd_io__readPrp
  public :: shd_io__checkPrpInput

  public :: shd_io__readSfcSize
  public :: shd_io__readSfc

  public :: shd_io__readCkdSize
  public :: shd_io__readCkd

  public :: shd_io__readSenSize
  public :: shd_io__readSen

  public :: shd_io__restoreState

contains

subroutine shd_io__userInput(&
    &RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
    &INSAVEFILE, OUTSAVEFILE,&
    &NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
    &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
    &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
    &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
    &OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

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
    character(len=FNAMEMAXLEN_) :: OUTFILES(*)
    character(len=FNAMEMAXLEN_) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=FNAMEMAXLEN_) :: OutFileNC, RUNNAME
    integer I, J, N
    real GNDEMIS

    SOLARMU = -1.0
    GNDALBEDO = 0.0
    WAVELEN = 1.0
    WAVENO(1) = 10000.
    WAVENO(2) = 10001.

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
    KDIST = CKDFILE(1:2) .ne. 'NO'
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
        if (SFCFILE(1:2) .eq. 'NO') then
            write (*, '(1X,A)') 'Ground albedo'
            read (*, *) GNDALBEDO
            write (*, *) GNDALBEDO
        end if
        SOLARMU = -abs(SOLARMU)
        SOLARAZ = SOLARAZ*acos(-1.0)/180.0
    end if
    GNDTEMP = 300.
    if (SFCFILE(1:2) .eq. 'NO') then
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
    if (NUMOUT .gt. MAXOUT) stop 'Error [_userInput]: MAXOUT exceeded.'
    do I = 1, NUMOUT
        write (*, '(1X,A,A)') 'Output type (R-radiance, ',&
            &'V-visualization, F-flux, H-heating,'
        write (*, '(14X,A)') 'S-spherical harmonic, J-source, M-medium)'
        read (*, '(A)') OUTTYPES(I)
        write (*, *) OUTTYPES(I)
        if (OUTTYPES(I) (1:1) .eq. 'R') then
            write (*, '(1X,A,A)')&
                &'Z level, X and Y spacing, X and Y offset, ',&
                &'number of angles, and mu, phi pairs'
            read (*, *) OUTPARMS(1, I), OUTPARMS(2, I), OUTPARMS(3, I),&
                &OUTPARMS(4, I), OUTPARMS(5, I),&
                &N, (OUTPARMS(J + 6, I), J=1, min(2*N, MAXPAR - 6))
            if (N .gt. (MAXPAR - 6)/2) then
                print *, 'Error [_userInput]: Number of radiance angles is too large for MAXPAR'
                stop
            end if
            N = min(N, (MAXPAR - 6)/2)
            OUTPARMS(6, I) = float(N)
            write (*, *) (OUTPARMS(J, I), J=1, 2*N + 6)
        else if (OUTTYPES(I) (1:1) .eq. 'V') then
            write (*, '(1X,A,A)')&
                &'1, bytes,scale, X,Y,Z, theta,phi,rotang, NL,NS, delline,delsamp'
            write (*, '(1X,A,A)')&
                &'2, bytes,scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1,scan2,delscan'
            read (*, *) (OUTPARMS(J, I), J=1, 13)
            write (*, *) (OUTPARMS(J, I), J=1, 13)
        else if (OUTTYPES(I) (1:1) .eq. 'F') then
            write (*, '(1X,A)')&
                &'Format (1,2,3,4,5), [Z, del X, del Y for format 2]'
            read (*, *) N, (OUTPARMS(1 + J, I), J=1, 3*max(0, 1 - abs(N - 2)))
            OUTPARMS(1, I) = float(N)
            if (N .eq. 2) then
                write (*, *) N, (OUTPARMS(1 + J, I), J=1, 3)
            else
                write (*, *) N
            end if
        else if (OUTTYPES(I) (1:1) .eq. 'H') then
            write (*, '(1X,A)') 'Format (1,2,3)'
            read (*, *) OUTPARMS(1, I)
            write (*, *) OUTPARMS(1, I)
        else if (OUTTYPES(I) (1:1) .eq. 'S') then
            write (*, '(1X,A)') 'Format (1,2)'
            read (*, *) OUTPARMS(1, I)
            write (*, *) OUTPARMS(1, I)
        else if (OUTTYPES(I) (1:1) .eq. 'J') then
            write (*, '(1X,A)') 'Format (1,2),  mu, phi '
            read (*, *) OUTPARMS(1, I), OUTPARMS(2, I), OUTPARMS(3, I)
            write (*, *) OUTPARMS(1, I), OUTPARMS(2, I), OUTPARMS(3, I)
        else if (OUTTYPES(I) (1:1) .eq. 'M') then
            write (*, '(1X,A)') 'Format (1,2)'
            read (*, *) OUTPARMS(1, I)
            write (*, *) OUTPARMS(1, I)
        end if
        write (*, '(1X,A,A)') 'Output file name'
        read (*, '(A)') OUTFILES(I)
        write (*, '(A)') OUTFILES(I)
    end do

    write (*, '(1X,A)') 'Output netcdf file name (or NONE or "")'
    read (*, '(A)') OutFileNC
    write (*, '(A)') OutFileNC
    if (OutFileNC .ne. 'NONE') then
        write(*, '(A)') 'Warning [_userInput]: Removed support for netCDF, set <OUTFILENC="NONE">.'
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
end subroutine shd_io__userInput

subroutine shd_io__nmlInput(&
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
    &dNUM_SH_TERM_FACTOR, dCELL_TO_POINT_RATIO)

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
    character(len=FNAMEMAXLEN_) :: dOUTFILES(*)
    character(len=FNAMEMAXLEN_) :: dOutFileNC, dRUNNAME
    character(len=FNAMEMAXLEN_) :: dPROPFILE, dSFCFILE, dCKDFILE
    character(len=FNAMEMAXLEN_) :: dINSAVEFILE, dOUTSAVEFILE, dSENFILE
    integer I, J, N

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
    character(len=FNAMEMAXLEN_) :: OUTFILES(MAXOUT)
    character(len=FNAMEMAXLEN_) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=FNAMEMAXLEN_) :: SENFILE, OutFileNC, RUNNAME
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
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

    character(len=FNAMEMAXLEN_) :: NMLFILE

    if (dMAXOUT .ne. MAXOUT) stop 'Error [_nmlInput]: MAXOUT not the same.'
    if (dMAXPAR .ne. MAXPAR) stop 'Error [_nmlInput]: MAXPAR not the same.'

    SOLARMU = -1.0
    GNDALBEDO = 0.0
    WAVELEN = 1.0
    WAVENO(1) = 10000.
    WAVENO(2) = 10001.
    SENFILE = 'NONE'

    ! Read in the namelist
    call get_command_argument(1, NMLFILE)
    write (*, '(A,"<",A,">")') 'Message [_nmlInput]: ', trim(adjustl(NMLFILE))
    open (UNIT=5, FILE=trim(adjustl(NMLFILE)), STATUS='OLD')
    read (5, NML=SHDOMINPUT)
    close (5)

    ! Make the assignments to the dummy variables
    dRUNNAME = RUNNAME
    dPROPFILE = PROPFILE
    dSFCFILE = SFCFILE
    dCKDFILE = CKDFILE
    dKDIST = CKDFILE(1:2) .ne. 'NO'
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
        dWAVELEN = 1.0e4/(0.5*(WAVENO(1) + WAVENO(2)))
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
    if (NUMOUT .gt. MAXOUT) stop 'Error [_nmlInput]: MAXOUT exceeded.'
    do I = 1, NUMOUT
        dOUTTYPES(I) = OUTTYPES(I)
        if (OUTTYPES(I) (1:1) .eq. 'R') then
            ! new binary angle input do not need the loop above
            if (SENFILE(1:2) .ne. 'NO') then
                do J = 1, 6
                    dOUTPARMS(J, I) = OUTPARMS(J, I)
                end do
            else
                N = nint(OUTPARMS(6, I))
                do J = 1, min(2*N + 6, MAXPAR)
                    dOUTPARMS(J, I) = OUTPARMS(J, I)
                end do
            end if
        else if (OUTTYPES(I) (1:1) .eq. 'V') then
            dOUTPARMS(1, I) = OUTPARMS(1, I)
            dOUTPARMS(2, I) = OUTPARMS(2, I)
            dOUTPARMS(3, I) = OUTPARMS(3, I)
            ! new user-specified sensor mode V3 only need 3 parameters
            ! read more parameters if for V1 and V2
            if (nint(OUTPARMS(1, I)) .le. 2) then
                do J = 4, 13
                    dOUTPARMS(J, I) = OUTPARMS(J, I)
                end do
            end if
        else if (OUTTYPES(I) (1:1) .eq. 'F') then
            dOUTPARMS(1, I) = OUTPARMS(1, I)
            if (nint(OUTPARMS(1, I)) .eq. 2) then
                do J = 2, 4
                    dOUTPARMS(J, I) = OUTPARMS(J, I)
                end do
            end if
        else if (OUTTYPES(I) (1:1) .eq. 'H') then
            dOUTPARMS(1, I) = OUTPARMS(1, I)
        else if (OUTTYPES(I) (1:1) .eq. 'S') then
            dOUTPARMS(1, I) = OUTPARMS(1, I)
        else if (OUTTYPES(I) (1:1) .eq. 'J') then
            dOUTPARMS(1, I) = OUTPARMS(1, I)
            dOUTPARMS(2, I) = OUTPARMS(2, I)
            dOUTPARMS(3, I) = OUTPARMS(3, I)
        else if (OUTTYPES(I) (1:1) .eq. 'M') then
            dOUTPARMS(1, I) = OUTPARMS(1, I)
        end if
        dOUTFILES(I) = OUTFILES(I)
    end do

    dOutFileNC = OutFileNC

    if (dOutFileNC .ne. 'NONE') then
        write(*, '(A)') 'Warning [_nmlInput]: Removed support for netCDF, set <OUTFILENC="NONE">.'
        dOutFileNC = 'NONE'
    end if

    dMAX_TOTAL_MB = MAX_TOTAL_MB
    dADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR
    dNUM_SH_TERM_FACTOR = NUM_SH_TERM_FACTOR
    dCELL_TO_POINT_RATIO = CELL_TO_POINT_RATIO

    return
end subroutine shd_io__nmlInput

subroutine shd_io__checkInputParam(&
    &NMU, NPHI, DELTAM, MAXASYM,&
    &GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU,&
    &WAVENO, WAVELEN, GNDTEMP, GNDALBEDO,&
    &SPLITACC, SHACC)

    ! Checks for validity of the input parameters. Gives warning for
    ! parameters that might not be correct.
    ! Check the letter flags.
    ! Delta-M should be used for solar transfer with forward peaked
    ! phase functions.
    ! Check for reasonable ranges on SPLITACC and SHACC.
    ! Check for reasonable range of GNDTEMP and valid range of GNDALBEDO.
    ! Check for SOLARMU=0.
    ! Check for NPHI<NMU.

    integer NMU, NPHI
    logical DELTAM
    real MAXASYM, SOLARFLUX, SOLARMU, WAVENO(2), WAVELEN
    real GNDTEMP, GNDALBEDO, SPLITACC, SHACC
    character(len=1) :: GRIDTYPE, SRCTYPE, UNITS
    real FLUX

    if (SRCTYPE .ne. 'S' .and. SRCTYPE .ne. 'T' .and. SRCTYPE .ne. 'B')&
        &stop 'Error [_checkInputParam]: Illegal source type.'
    if (UNITS .ne. 'R' .and. UNITS .ne. 'T' .and. UNITS .ne. 'B')&
        &stop 'Error [_checkInputParam]: Illegal units.'
    if (GRIDTYPE .ne. 'E' .and. GRIDTYPE .ne. 'P' .and. GRIDTYPE .ne. 'F')&
        &stop 'Error [_checkInputParam]: Illegal grid type.'

    if (SRCTYPE .ne. 'T' .and. SOLARMU .eq. 0.0) then
        write (*, *) 'Error [_checkInputParam]: Cosine solar zenith angle'
        write (*, *) ' (SOLARMU) cannot be zero.'
        stop
    end if

    if (.not. DELTAM .and. SRCTYPE .ne. 'T' .and. MAXASYM .gt. 0.5) then
        write (*, *) 'Warning [_checkInputParam]: Delta-M should be used'
        write (*, *) ' (DELTAM=.TRUE.) for solar radiative transfer'
        write (*, *) ' with highly peaked phase functions.'
        write (*, *)
    end if

    ! Check for SPLITACC and SHACC magnitudes
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        call shd_util__computePlanck(GNDTEMP, UNITS, WAVENO, WAVELEN, FLUX)
        FLUX = PI_*FLUX
    else
        FLUX = 0.0
    end if
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        FLUX = FLUX + SOLARFLUX
    end if
    if (SPLITACC .gt. 0.0 .and.&
        &(SPLITACC .lt. 0.001*FLUX .or. SPLITACC .gt. 0.1*FLUX)) then
        write (*, *) 'Warning [_checkInputParam]: Splitting accuracy parameter'
        write (*, *) ' is not a fraction, but scales with fluxes'
        write (*, *)
    end if
    if (SHACC .gt. 0.0 .and. SHACC .gt. 0.03*FLUX) then
        write (*, *) 'Warning [_checkInputParam]: Spherical harmonic accuracy'
        write (*, *) ' (SHACC) is not a fraction, but scales with fluxes'
        write (*, *)
    end if

    if (SRCTYPE .ne. 'S' .and.&
        &(GNDTEMP .lt. 150. .or. GNDTEMP .gt. 350.)) then
        write (*, *) 'Warning [_checkInputParam]: '
        write (*, *) '  ground temperature out of Earth range.'
        write (*, *)
    end if
    if (GNDALBEDO .lt. 0.0 .or. GNDALBEDO .gt. 1.0) then
        write (*, *) 'Error [_checkInputParam]: Ground albedo must be 0 to 1.'
        stop
    end if

    if (NPHI .lt. NMU) then
        write (*, *) 'Warning [_checkInputParam]: Usually want NPHI >= NMU.'
        write (*, *)
    end if

    return
end subroutine shd_io__checkInputParam

subroutine shd_io__readPrpSize(PROPFILE, NSTLEG, NLEG, NPX, NPY, NPZ,&
    &NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)

    ! Reads parts of the property file to get the maximum array sizes
    ! needed for allocatable arrays.  For extinction only and tabulated
    ! phase function formats, only the header is read, while for the
    ! standard format, the whole file must be read to determine MAXLEG.

    integer NSTLEG, NPX, NPY, NPZ
    integer NUMPHASE, NLEG, MAXLEG, MAXPGL
    real DELX, DELY
    character(len=*), intent(in) :: PROPFILE
    integer NUML, IX, IY, IZ, I, J, K, L
    real ZLEVELS, TMP, EXT, ALB, CHI
    character(len=1) :: PROPTYPE

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=trim(adjustl(PROPFILE)), STATUS='OLD')
    read (1, '(A1)') PROPTYPE
    if (PROPTYPE .ne. 'E' .and. PROPTYPE .ne. 'T'&
        &.and. PROPTYPE .ne. 'P') rewind (1)
    if (PROPTYPE .ne. 'P' .and. NSTLEG .gt. 1) then
        write (*, '(A,A)') 'Error [_readPrpSize]: NSTOKES>1 requires ',&
            &'Polarized phase function property file'
        stop
    end if
    read (1, *) NPX, NPY, NPZ
    read (1, *) DELX, DELY, (ZLEVELS, K=1, NPZ)

    ! Property file type E is for extinction only format
    if (PROPTYPE .eq. 'E') then
        NUMPHASE = 1
        read (1, *) (TMP, IZ=1, NPZ)
        read (1, *) ALB, NUML
        MAXLEG = max(NLEG, NUML)
        MAXPGL = MAXLEG + 1

    ! Property file type T is for tabulated phase function format
    else if (PROPTYPE .eq. 'T') then
        MAXLEG = NLEG
        read (1, *) NUMPHASE
        do I = 1, NUMPHASE
            read (1, *) NUML, (CHI, L=1, NUML)
            MAXLEG = max(NUML, MAXLEG)
        end do
        MAXPGL = NUMPHASE*(MAXLEG + 1)

    ! Property file type P is for polarized tabulated phase functi
    else if (PROPTYPE .eq. 'P') then
        MAXLEG = NLEG
        read (1, *) NUMPHASE
        do I = 1, NUMPHASE
            do J = 1, 6
                read (1, *) K, NUML, (CHI, L=0, NUML)
                MAXLEG = max(NUML, MAXLEG)
            end do
        end do
        MAXPGL = NSTLEG*NUMPHASE*(MAXLEG + 1)

    else
        ! Standard property file has everything variable
        MAXLEG = NLEG
        do while (.true.)
            read (1, *, end=190) IX, IY, IZ, TMP, EXT, ALB, NUML,&
                &(CHI, L=1, NUML)
            MAXLEG = max(NUML, MAXLEG)
        end do
190     continue
        MAXPGL = NPX*NPY*NPZ*(MAXLEG + 1)
    end if
    close (1)
    return
end subroutine shd_io__readPrpSize

subroutine shd_io__readPrp(&
    &PROPFILE, NPX, NPY, NPZ,&
    &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,&
    &PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM,&
    &TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)

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
    ! NUMPHASE is the number of phase functions, otherwise NUMPHASE=0.
    ! If doing delta-M then NLEG is the minimum number of Legendre terms on
    ! input and the actual maximum number of terms on output, except that
    ! it may not exceed MAXLEG; otherwise NLEG is the number of Legendre
    ! terms to be used (regardless of what is in property file).
    ! See the overall program documentation for the four file formats.

    integer NPX, NPY, NPZ, NSTLEG, MAXLEG, NLEG, NUMPHASE
    integer MAXNZ, MAXPG, MAXPGL
    integer*2 IPHASEP(*)
    logical DELTAM
    real DELX, DELY, ZLEVELS(*)
    real TEMPP(*), EXTINCTP(*)
    real ALBEDOP(*), LEGENP(*), MAXASYM
    character(len=*), intent(in) :: PROPFILE
    character(len=1) :: PROPTYPE
    integer I, IX, IY, IZ, J, JS, K, KL, L, IL, NUML, IPH, N
    real :: TMP, EXT, ALB, CHI

    character(len=FNAMEMAXLEN_) :: OUTFILE_DATA
    character(len=FNAMEMAXLEN_) :: TMP_CHAR
    character(len=FNAMEMAXLEN_) :: DIMS_CHAR
    character(len=20) :: DATA_POSTFIX
    integer :: NXX, NYY, NZZ, NPP

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=trim(adjustl(PROPFILE)), STATUS='OLD')
    read (1, '(A1)') PROPTYPE
    if (PROPTYPE .ne. 'E' .and. PROPTYPE .ne. 'T'&
        &.and. PROPTYPE .ne. 'P') rewind (1)
    if (PROPTYPE .ne. 'P' .and. NSTLEG .gt. 1) then
        write (*, '(A,A)') 'Error [_readPrp]: NSTOKES>1 requires ',&
            &'Polarized phase function property file.'
        stop
    end if
    read (1, *) NPX, NPY, NPZ
    if (NPZ .gt. MAXNZ) stop 'Error [_readPrp]: MAXNZ exceeded.'
    read (1, *) DELX, DELY, (ZLEVELS(K), K=1, NPZ)

    ! Property file type E is for extinction only format
    if (PROPTYPE .eq. 'E') then
        NUMPHASE = 1
        read (1, *) (TMP, IZ=1, NPZ)
        read (1, *) ALB, NUML
        if (DELTAM) then
            NLEG = max(NLEG, NUML)
            NLEG = min(NLEG, MAXLEG)
        end if
        if (NPX*NPY*NPZ .gt. MAXPG)&
            &stop 'Error [_readPrp]: MAXPG exceeded.'
        if (NLEG + 1 .gt. MAXPGL)&
            &stop 'Error [_readPrp]: MAXPGL exceeded.'
        rewind (1)
        read (1, *)
        read (1, *) NPX, NPY, NPZ
        read (1, *) DELX, DELY, (ZLEVELS(K), K=1, NPZ)
        read (1, *) (TEMPP(IZ), IZ=1, NPZ)
        read (1, *) ALB, NUML, (LEGENP(L + 1), L=1, min(NUML, NLEG)),&
            &(CHI, L=NLEG + 1, NUML)

        LEGENP(1) = 1.0
        MAXASYM = LEGENP(1 + 1)/3.0
        do L = NUML + 1, NLEG
            LEGENP(L) = 0.0
        end do

        IY = 1
        do while (.true.)
            if (NPY .eq. 1) then
                read (1, *, end=290) IX, IZ, EXT
            else
                read (1, *, end=290) IX, IY, IZ, EXT
            end if
            if (IX .ge. 1 .and. IX .le. NPX .and.&
               &IY .ge. 1 .and. IY .le. NPY .and.&
               &IZ .ge. 1 .and. IZ .le. NPZ) then
                K = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
                TEMPP(K) = TEMPP(IZ)
                EXTINCTP(K) = EXT
                ALBEDOP(K) = ALB
                IPHASEP(K) = 1
            end if
        end do

    ! Property file type T or P is for tabulated phase function format
    else if (PROPTYPE .eq. 'T' .or. PROPTYPE .eq. 'P') then

        read (1, *) NUMPHASE
        ! If delta-M then find the largest number of Legendre terms
        if (DELTAM) then

            do I = 1, NUMPHASE
                if (PROPTYPE .eq. 'T') then
                    read (1, *) NUML, (CHI, L=1, NUML)
                else
                    do J = 1, 6
                        read (1, *) JS, NUML, (CHI, L=0, NUML)
                    end do
                end if
                NLEG = max(NUML, NLEG)
            end do

            if (NLEG .gt. MAXLEG) then
                write (*, '(A,I4,A,I4,A)')&
                    &'Warning [_readPrp]: Truncating phase function from ',&
                    &NLEG, ' to ', MAXLEG, ' Legendre terms.'
            end if
            NLEG = min(NLEG, MAXLEG)

        end if

        N = NPX*NPY*NPZ
        if (N .gt. MAXPG) stop 'Error [_readPrp]: MAXPG exceeded.'
        if (NUMPHASE*(NLEG + 1)*NSTLEG .gt. MAXPGL)&
            &stop 'Error [_readPrp]: MAXPGL exceeded.'
        rewind (1)
        read (1, *)
        read (1, *) NPX, NPY, NPZ
        read (1, *) DELX, DELY, (ZLEVELS(K), K=1, NPZ)
        read (1, *) NUMPHASE
        MAXASYM = -1.0
        do I = 1, NUMPHASE
            IL = (NLEG + 1)*(I - 1)
            if (PROPTYPE .eq. 'T') then
                LEGENP(IL + 1) = 1.0
                read (1, *) NUML, (LEGENP(IL + L + 1), L=1, min(NUML, NLEG)),&
                    &(CHI, L=NLEG + 1, NUML)
                do L = NUML + 1, NLEG
                    LEGENP(L + IL + 1) = 0.0
                end do
            else
                do J = 1, 6
                    if (J .le. NSTLEG) then
                        read (1, *) JS, NUML,&
                            &(LEGENP((IL + L)*NSTLEG + J), L=0, min(NUML, NLEG)),&
                            &(CHI, L=NLEG + 1, NUML)
                        do L = NUML + 1, NLEG
                            LEGENP((IL + L)*NSTLEG + J) = 0.0
                        end do
                    else
                        read (1, *) JS, NUML, (CHI, L=0, NUML)
                    end if
                end do
            end if
            MAXASYM = max(MAXASYM, LEGENP((IL + 1)*NSTLEG + 1)/3.0)
        end do

        read (1, '(A)') TMP_CHAR ! skip one line

        if (TMP_CHAR(1:1) .eq. '!') then

            read (1, '(A)') DATA_POSTFIX
            read (1, '(A)') DIMS_CHAR

            read (DIMS_CHAR( 3:12), '(I10)') NXX
            read (DIMS_CHAR(14:23), '(I10)') NYY
            read (DIMS_CHAR(25:34), '(I10)') NZZ
            read (DIMS_CHAR(36:45), '(I10)') NPP

            OUTFILE_DATA = trim(adjustl(PROPFILE))//trim(adjustl(DATA_POSTFIX(2:)))

            open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='OLD',&
                &ACCESS='STREAM', FORM='UNFORMATTED')

            ! TEMPP: real(4)
            read (88) (((   TEMPP(IZ+NPZ*(IY-1)+NPZ*NPY*(IX-1)), IZ=1,NZZ), IY=1,NYY), IX=1,NXX)
            ! EXTINCTP: real(4)
            read (88) (((EXTINCTP(IZ+NPZ*(IY-1)+NPZ*NPY*(IX-1)), IZ=1,NZZ), IY=1,NYY), IX=1,NXX)
            ! ALBEDOP: real(4)
            read (88) ((( ALBEDOP(IZ+NPZ*(IY-1)+NPZ*NPY*(IX-1)), IZ=1,NZZ), IY=1,NYY), IX=1,NXX)
            ! IPHASEP: integer(2)
            read (88) ((( IPHASEP(IZ+NPZ*(IY-1)+NPZ*NPY*(IX-1)), IZ=1,NZZ), IY=1,NYY), IX=1,NXX)

            close (88)

        else

            backspace(1)

            ! original
            !╭────────────────────────────────────────────────╮!
            do while (.true.)
                read (1, *, end=290) IX, IY, IZ, TMP, EXT, ALB, IPH
                if (IX .ge. 1 .and. IX .le. NPX .and.&
                   &IY .ge. 1 .and. IY .le. NPY .and.&
                   &IZ .ge. 1 .and. IZ .le. NPZ) then
                    K = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
                    TEMPP(K) = TMP
                    EXTINCTP(K) = EXT
                    ALBEDOP(K) = ALB
                    IPHASEP(K) = IPH
                end if
            end do
            !╰────────────────────────────────────────────────╯!

        end if

    else

        ! Standard property file has everything variable
        !   If delta-M then find the largest number of Legendre terms
        NUMPHASE = 0
        if (DELTAM) then
            do while (.true.)
                read (1, *, end=190) IX, IY, IZ, TMP, EXT, ALB, NUML,&
                    &(CHI, L=1, NUML)
                NLEG = max(NUML, NLEG)
            end do
190         continue
            NLEG = min(NLEG, MAXLEG)
            rewind (1)
            read (1, *) NPX, NPY, NPZ
            read (1, *) DELX, DELY, (ZLEVELS(K), K=1, NPZ)
        end if
        N = NPX*NPY*NPZ
        if (N .gt. MAXPG) stop 'Error [_readPrp]: MAXPG exceeded.'
        if (N*(NLEG + 1) .gt. MAXPGL)&
            &stop 'Error [_readPrp]: MAXPGL exceeded.'
        MAXASYM = -1.0
        do while (.true.)
            read (1, *, end=290) IX, IY, IZ, TMP, EXT, ALB, NUML
            if (IX .ge. 1 .and. IX .le. NPX .and.&
                &IY .ge. 1 .and. IY .le. NPY .and.&
                &IZ .ge. 1 .and. IZ .le. NPZ) then
                K = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
                KL = (NLEG + 1)*(K - 1)
                backspace (1)
                read (1, *, end=290) IX, IY, IZ, TEMPP(K), EXTINCTP(K),&
                    &ALBEDOP(K), NUML, (LEGENP(L + KL), L=1, min(NUML, NLEG)),&
                    &(CHI, L=NLEG + 1, NUML)
                MAXASYM = max(MAXASYM, LEGENP(KL + 1 + 1)/3.0)
                LEGENP(1 + KL) = 1.0
                do L = NUML + 1, NLEG
                    LEGENP(L + KL) = 0.0
                end do
            end if
        end do

    end if

290 continue
    close (1)
    return

end subroutine shd_io__readPrp

subroutine shd_io__checkPrpInput(&
    &NPX, NPY, NPZ, NSTLEG, NLEG,&
    &DELX, DELY, ZLEVELS,&
    &TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)

    ! Checks the input medium property values for correctness.
    ! Warning if temperature below 150 K or above 350 K.
    ! Stop if negative extinction, single scattering albedo below 0 or
    ! above 1, asymmetry parameter below -1 or above 1, or IPHASE
    ! out of range (of tabulated phase functions).
    ! Warning if Chi1=1, for people specifying Chi0.
    ! Warning if optical depth across a grid cell exceeds 2.

    integer NPX, NPY, NPZ, NSTLEG, NLEG, NUMPHASE
    integer*2 IPHASEP(*)
    real DELX, DELY, ZLEVELS(*)
    real TEMPP(*), EXTINCTP(*)
    real ALBEDOP(*), LEGENP(NSTLEG, 0:NLEG, *)
    integer I, IZ, N, NTHICK
    logical LEGWARN
    real TMIN, TMAX, TAU, MAXTAU, DZ

    ! NPX and NPY must be 1 or greater, NPZ must be 2 or greater
    if (NPX .lt. 1 .or. NPY .lt. 1 .or. NPZ .lt. 2) then
        write (*, *) 'Error [_checkPrpInput]: Bad grid size - ',&
            &NPX, NPY, NPZ
        stop
    end if

    if (DELX .le. 0.0 .or. DELY .le. 0.0) then
        write (*, *) 'Error [_checkPrpInput]: Bad grid spacing - ',&
            &'DELX,DELY: ', DELX, DELY
        stop
    end if

    ! The vertical levels must be increasing
    do IZ = 1, NPZ - 1
        if (ZLEVELS(IZ + 1) .le. ZLEVELS(IZ)) then
            write (*, *) 'Error [_checkPrpInput]: ZLEVELS must increase : ',&
                &IZ, ZLEVELS(IZ), ZLEVELS(IZ + 1)
            stop
        end if
    end do

    ! Loop over the property grid points, checking parameters
    !   Warn about temperatures outside of Earth atmosphere temperat
    N = NPX*NPY*NPZ
    TMIN = 1000.0
    TMAX = -1000.0
    do I = 1, N
        TMIN = min(TMIN, TEMPP(I))
        TMAX = max(TMAX, TEMPP(I))
    end do

    if (TMIN .lt. 150.0 .or. TMAX .gt. 350.0) then
        write (*, *) 'Warning [_checkPrpInput]: ',&
            &'temperature min and max (Kelvin): ', TMIN, TMAX
        write (*, *)
    end if

    do I = 1, N
        if (EXTINCTP(I) .lt. 0.0) then
            write (*, *) 'Error [_checkPrpInput]: Negative extinction - ',&
                &EXTINCTP(I), ' at grid point ', I
            stop
        end if
        if (ALBEDOP(I) .lt. 0.0 .or. ALBEDOP(I) .gt. 1.0) then
            write (*, *) 'Error [_checkPrpInput]: Single scattering ',&
                &'albedo out of range - ', ALBEDOP(I), ' at grid point ', I
            stop
        end if
    end do

    if (NUMPHASE .gt. 0) N = NUMPHASE
    LEGWARN = .false.
    do I = 1, N
        if (LEGENP(1, 1, I) .lt. -3.0 .or. LEGENP(1, 1, I) .gt. 3.0) then
            write (*, *) 'Error [_checkPrpInput]: Asymmetry parameter ',&
                &'out of range - ', LEGENP(1, 1, I)/3.0, ' at grid point ', I
            stop
        end if
        if (LEGENP(1, 1, I) .eq. 1.0) LEGWARN = .true.
    end do

    if (LEGWARN) then
        write (*, *) 'Warning [_checkPrpInput]: ',&
            &'an asymmetry parameter is 1/3 (Chi1=1); ',&
            &'  Chi0 should not be specified for Chi1.'
        write (*, *)
    end if

    if (NUMPHASE .gt. 0) then
        N = NPX*NPY*NPZ
        do I = 1, N
            if (IPHASEP(I) .lt. 1 .or. IPHASEP(I) .gt. NUMPHASE) then
                write (*, *) 'Error [_checkPrpInput]: Phase function pointer',&
                    &' out of range - ', IPHASEP(I), ' at grid point ', I
                stop
            end if
        end do
    end if

    N = NPX*NPY*NPZ
    NTHICK = 0
    MAXTAU = 0.0
    do IZ = 1, NPZ - 1
        DZ = ZLEVELS(IZ + 1) - ZLEVELS(IZ)
        do I = 0, N - 1, NPZ
            TAU = 0.5*(EXTINCTP(I + IZ) + EXTINCTP(I + IZ + 1))*DZ
            if (TAU .gt. 2.0) then
                NTHICK = NTHICK + 1
                MAXTAU = max(MAXTAU, TAU)
            end if
        end do
    end do

    if (NTHICK .gt. 0) then
        write (*, '(A)') 'Message [_checkPrpInput]:'
        write (*, '(A,I0)')   '  Npoints (optical depth > 2): ', NTHICK
        write (*, '(A,F0.2)') '       Max cell optical depth: ', MAXTAU
    end if

    return

end subroutine shd_io__checkPrpInput

subroutine shd_io__readSfcSize(SFCFILE, MAXSFCPTS, MAXSFCPARS)

    ! Gets the size of the surface property inputs (temperature is
    ! always included).

    integer MAXSFCPTS, MAXSFCPARS
    character(len=*), intent(in) :: SFCFILE
    integer NXSFC, NYSFC
    character(len=2) :: SFCTYPE

    ! Open the file, figure out the type, and get the grid size
    SFCTYPE(1:1) = 'V'
    open (UNIT=1, FILE=trim(adjustl(SFCFILE)), STATUS='OLD')
    read (1, '(A1)') SFCTYPE(2:2)
    read (1, *) NXSFC, NYSFC
    MAXSFCPTS = (NXSFC + 1)*(NYSFC + 1)

    ! 0
    if ((SFCTYPE .eq. 'VL') .or. (SFCTYPE .eq. 'Vl')) then
        ! Surface file type L is for variable Lambertian surface
        MAXSFCPARS = 2
    ! 1
    else if ((SFCTYPE .eq. 'VO') .or. (SFCTYPE .eq. 'Vo')) then
        ! Surface file type O is for Ocean BRDF surface
        MAXSFCPARS = 3
    ! 2
    else if ((SFCTYPE .eq. 'VR') .or. (SFCTYPE .eq. 'Vr')) then
        ! Surface file type R is for RPV surface
        MAXSFCPARS = 4
    ! 3
    else if ((SFCTYPE .eq. 'VT') .or. (SFCTYPE .eq. 'Vt')) then
        ! Surface file type T is for LSRT surface
        MAXSFCPARS = 4
    ! 4
    else if ((SFCTYPE .eq. 'VJ') .or. (SFCTYPE .eq. 'Vj')) then
        ! Surface file type J is for LSRT-Jiao surface
        MAXSFCPARS = 6
    ! -1
    else if ((SFCTYPE .eq. 'VW') .or. (SFCTYPE .eq. 'Vw')) then
        ! Surface file type W is for Wave-Fresnel surface
        MAXSFCPARS = 4
    ! -2
    else if ((SFCTYPE .eq. 'VD') .or. (SFCTYPE .eq. 'Vd')) then
        ! Surface file type D is for Diner polarized RPV surface:
        MAXSFCPARS = 6
    else if ((SFCTYPE .eq. 'VX') .or. (SFCTYPE .eq. 'Vx')) then
        ! Surface file type X is for Mixed-surface type:
        MAXSFCPARS = 8
    else
        stop 'Error [_readSfcSize]: Unknown BRDF type.'
    end if
    close (1)
    return

end subroutine shd_io__readSfcSize

subroutine shd_io__readSfc(&
    &SFCFILE, MAXSFCPTS, MAXSFCPARS,&
    &SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,&
    &NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)

    ! Reads the surface properties from the file into the surface arrays.
    ! There are currently six surface types (Lambertian, Wave-Fresnel,
    ! Ocean, RPV, and Diner polarized RPV), though it is simple to add more
    ! In all cases the surface properties are specified on a regular, evenly.
    ! spaced grid of size NXSFC by NYSFC and spacing DELXSFC by DELYSFC.
    ! The properties are returned on a NXSFC+1 by NYSFC+1 grid that is
    ! periodic (NX+1 same as 1, etc).  The surface properties at each
    ! grid point are returned in SFCPARMS, and the first parameter must
    ! be temperature.  The number of parameters for each grid point
    ! (including temperature) are returned in NSFCPAR.
    ! Also returned are the domain average ground albedo and temperature
    ! (GNDALBEDO, GNDTEMP) for use with initialization (hence they can be
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

    integer :: MAXSFCPTS, MAXSFCPARS, NXSFC, NYSFC, NSFCPAR
    real :: DELXSFC, DELYSFC
    real :: SFCPARMS(*), GNDTEMP, GNDALBEDO
    character(len=*), intent(in) :: SFCFILE
    character(len=2) :: SFCTYPE
    integer :: N, I, I0, IX, IY, J
    real :: ALB, TEMP, MRE, MIM, RHO0, KR, THETA, WSPD, PCL
    real :: PARAM1, PARAM2, PARAM3, PARAM4, PARAM5, PARAM6, PARAM7, PARAM8
    real F_ISO, F_VOL, F_GEO, F_JIAO
    real A, B, ZETA, SIGMA, ALPHA

    character(len=FNAMEMAXLEN_) :: OUTFILE_DATA
    character(len=FNAMEMAXLEN_) :: TMP_CHAR
    character(len=FNAMEMAXLEN_) :: DIMS_CHAR
    character(len=20) :: DATA_POSTFIX
    integer :: NXX, NYY, NPP, IDATA, NDATA

    ! Surface type is variable something
    SFCTYPE(1:1) = 'V'
    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=trim(adjustl(SFCFILE)), STATUS='OLD')
    read (1, '(A1)') SFCTYPE(2:2)
    read (1, *) NXSFC, NYSFC, DELXSFC, DELYSFC
    if ((NXSFC + 1)*(NYSFC + 1) .gt. MAXSFCPTS) then
        write (*, *) 'Error [_readSfc]: MAXSFCPTS exceeded.', (NXSFC + 1)*(NYSFC + 1), MAXSFCPTS
        stop
    end if

    GNDTEMP = 0.0
    GNDALBEDO = 0.0
    N = 0

    ! Surface file type L is for variable Lambertian surface
    if ((SFCTYPE .eq. 'VL') .or. (SFCTYPE .eq. 'Vl')) then

        NSFCPAR = 2

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, ALB
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = ALB
        !        if (ALB .lt. 0.0 .or. ALB .gt. 1.0)&
        !            &stop 'Error [_readSfc]: Illegal surface albedo.'
        !        if (TEMP .lt. 150. .or. TEMP .gt. 350.)&
        !            &write (*, *) 'Warning [_readSfc]: ',&
        !            &' surface temperature: ', IX, IY, TEMP
        !        GNDTEMP = GNDTEMP + TEMP
        !        GNDALBEDO = GNDALBEDO + ALB
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type O is for Ocean surface BRDF
    !   Norm Loeb's modification of 6S ocean reflectance module.
    else if ((SFCTYPE .eq. 'VO') .or. (SFCTYPE .eq. 'Vo')) then

        NSFCPAR = 3

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, ALB, WSPD, PCL
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = ALB
        !        SFCPARMS(I + 3) = WSPD
        !        SFCPARMS(I + 4) = PCL
        !        GNDTEMP = GNDTEMP + TEMP
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type R is for RPV surface:
    !   (Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphe
    !   Reflectance (CSAR) Model. 2. Semiempirical Surface Model
    !   Usable With NOAA Advanced Very High Resolution Radiometer
    !   Data, J. Geophys. Res., 98, 20791-20801.)
    else if ((SFCTYPE .eq. 'VR') .or. (SFCTYPE .eq. 'Vr')) then

        NSFCPAR = 4

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, RHO0, KR, THETA
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = RHO0
        !        SFCPARMS(I + 3) = KR
        !        SFCPARMS(I + 4) = THETA
        !        GNDTEMP = GNDTEMP + TEMP
        !        GNDALBEDO = GNDALBEDO + RHO0
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type T is for LSRT surface:
    else if ((SFCTYPE .eq. 'VT') .or. (SFCTYPE .eq. 'Vt')) then

        NSFCPAR = 4

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, F_ISO, F_GEO, F_VOL
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = F_ISO
        !        SFCPARMS(I + 3) = F_GEO
        !        SFCPARMS(I + 4) = F_VOL
        !        GNDTEMP = GNDTEMP + TEMP
        !        GNDALBEDO = GNDALBEDO + F_ISO
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type J is for LSRT-Jiao surface:
    else if ((SFCTYPE .eq. 'VJ') .or. (SFCTYPE .eq. 'Vj')) then

        NSFCPAR = 6

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, F_ISO, F_GEO, F_VOL, F_JIAO, ALPHA
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = F_ISO
        !        SFCPARMS(I + 3) = F_GEO
        !        SFCPARMS(I + 4) = F_VOL
        !        SFCPARMS(I + 5) = F_JIAO
        !        SFCPARMS(I + 6) = ALPHA
        !        GNDTEMP = GNDTEMP + TEMP
        !        GNDALBEDO = GNDALBEDO + F_ISO
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type W is for Wave-Fresnel surface
    else if ((SFCTYPE .eq. 'VW') .or. (SFCTYPE .eq. 'Vw')) then

        NSFCPAR = 4

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, ALB, MRE, MIM, WSPD
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = ALB
        !        SFCPARMS(I + 3) = MRE
        !        SFCPARMS(I + 4) = MIM
        !        SFCPARMS(I + 5) = WSPD
        !        GNDTEMP = GNDTEMP + TEMP
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type D is for Diner et al surface reflection:
    ! (Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
    ! V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
    ! Exploration of a polarized surface bidirectional reflectance model
    ! using the ground-based multiangle spectropolarimetric imager.
    ! Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.)
    else if ((SFCTYPE .eq. 'VD') .or. (SFCTYPE .eq. 'Vd')) then

        NSFCPAR = 6

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290) IX, IY, TEMP, A, KR, B, ZETA, SIGMA
        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = TEMP
        !        SFCPARMS(I + 2) = A
        !        SFCPARMS(I + 3) = KR
        !        SFCPARMS(I + 4) = B
        !        SFCPARMS(I + 5) = ZETA
        !        SFCPARMS(I + 6) = SIGMA
        !        GNDTEMP = GNDTEMP + TEMP
        !        GNDALBEDO = GNDALBEDO + A
        !        N = N + 1
        !    end if
        !end do

    ! Surface file type X is for Mixed-type surface:
    else if ((SFCTYPE .eq. 'VX') .or. (SFCTYPE .eq. 'Vx')) then

        NSFCPAR = 8

        !if (NSFCPAR .gt. MAXSFCPARS) then
        !    write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.'
        !    stop
        !end if
        !do while (.true.)
        !    read (1, *, end=290)    IX, IY, PARAM1, PARAM2, PARAM3, PARAM4, PARAM5, PARAM6, PARAM7, PARAM8
        !    ! for Lambertian   (VL):IX, IY, TEMP,   ALB,    -999.0, -999.0, -999.0, -999.0, -999.0,  0.0
        !    ! for Ocean        (VO):IX, IY, TEMP,   WSPD,   PCL,    -999.0, -999.0, -999.0, -999.0,  1.0
        !    ! for RPV          (VR):IX, IY, TEMP,   RHO0,   KR,     THETA,  -999.0, -999.0, -999.0,  2.0
        !    ! for LSRT         (VT):IX, IY, TEMP,   F_ISO,  F_GEO,  F_VOL,  -999.0, -999.0, -999.0,  3.0
        !    ! for LSRT-J       (VJ):IX, IY, TEMP,   F_ISO,  F_GEO,  F_VOL,  F_JIAO, ALPHA,  -999.0,  4.0
        !    ! for Wave-Fresnel (VW):IX, IY, TEMP,   MRE,    MIM,    WSPD,   -999.0, -999.0, -999.0, -1.0
        !    ! for Diner        (VD):IX, IY, TEMP,   A,      KR,     B,      ZETA,   SIGMA,  -999.0, -2.0

        !    if (IX .ge. 1 .and. IX .le. NXSFC .and.&
        !        &IY .ge. 1 .and. IY .le. NYSFC) then
        !        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
        !        SFCPARMS(I + 1) = PARAM1
        !        SFCPARMS(I + 2) = PARAM2
        !        SFCPARMS(I + 3) = PARAM3
        !        SFCPARMS(I + 4) = PARAM4
        !        SFCPARMS(I + 5) = PARAM5
        !        SFCPARMS(I + 6) = PARAM6
        !        SFCPARMS(I + 7) = PARAM7
        !        SFCPARMS(I + 8) = PARAM8 ! surface identification

        !        GNDTEMP = GNDTEMP + PARAM1

        !        ! no albedo for polarized (W) and unpolarized (O) Ocean parameterization
        !        if ((abs(PARAM8)-1.0) .gt. 0.01) then
        !            GNDALBEDO = GNDALBEDO + PARAM2
        !        end if

        !        N = N + 1
        !    end if
        !end do

    else

        stop 'Error [_readSfc]: Unknown BRDF type <'//trim(adjustl(SFCTYPE))//'>.'

    end if

    if (NSFCPAR .gt. MAXSFCPARS) then
        write (*, *) 'Error [_readSfc]: MAXSFCPARS exceeded.', NSFCPAR, MAXSFCPARS
        stop
    end if

    ! read binary
    !╭────────────────────────────────────────────────────────────────────────────╮!
    read (1, '(A)') TMP_CHAR ! skip one line
    read (1, '(A)') DATA_POSTFIX
    read (1, '(A)') DIMS_CHAR

    read (DIMS_CHAR( 3:12), '(I10)') NXX
    read (DIMS_CHAR(14:23), '(I10)') NYY
    read (DIMS_CHAR(25:34), '(I10)') NPP

    OUTFILE_DATA = trim(adjustl(SFCFILE))//trim(adjustl(DATA_POSTFIX(2:)))
    open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='OLD',&
        &ACCESS='STREAM', FORM='UNFORMATTED')
    NDATA = NXX*NYY*NPP
    read (88) (SFCPARMS(IDATA), IDATA=1,NDATA)
    close (88)
    !╰────────────────────────────────────────────────────────────────────────────╯!

290 continue

    close (1)

    ! Copy the edge points to the opposite side for periodic boundar
    do IX = 1, NXSFC
        I0 = NSFCPAR*(IX - 1 + (NXSFC + 1)*(1 - 1))
        I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(NYSFC + 1 - 1))
        do J = 1, NSFCPAR
            SFCPARMS(I + J) = SFCPARMS(I0 + J)
        end do
    end do

    do IY = 1, NYSFC + 1
        I0 = NSFCPAR*((NXSFC + 1)*(IY - 1))
        I = NSFCPAR*(NXSFC + (NXSFC + 1)*(IY - 1))
        do J = 1, NSFCPAR
            SFCPARMS(I + J) = SFCPARMS(I0 + J)
        end do
    end do

    do IY = 1, NYSFC + 1
        do IX = 1, NXSFC + 1

            if (IX .ge. 1 .and. IX .le. NXSFC .and.&
                &IY .ge. 1 .and. IY .le. NYSFC) then

                I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))

                GNDTEMP = GNDTEMP + SFCPARMS(I + 1)

                if (((SFCTYPE .ne. 'VO') .and. (SFCTYPE .ne. 'Vo')) .and.&
                    &((SFCTYPE .ne. 'VW') .and. (SFCTYPE .ne. 'Vw'))) then
                    if ((SFCTYPE .eq. 'VX') .or. (SFCTYPE .eq. 'Vx')) then
                        if ((abs(SFCPARMS(I+8))-1.0) .gt. 0.01) then ! exclude [O] and [W]
                            GNDALBEDO = GNDALBEDO + SFCPARMS(I + 2)
                        end if
                    else
                        GNDALBEDO = GNDALBEDO + SFCPARMS(I + 2)
                    end if
                end if

                N = N + 1

            end if

        end do
    end do

    GNDTEMP   = GNDTEMP/N
    GNDALBEDO = GNDALBEDO/N
    return

end subroutine shd_io__readSfc

subroutine shd_io__readCkdSize(CKDFILE, WAVENO, NG, NZCKD)

    ! Find the size of the k-distribution arrays for this wavenumber band

    integer NG, NZCKD
    real WAVENO(2)
    character(len=*), intent(in) :: CKDFILE
    integer NB, IB, JB, KB, N
    real WAVENUM1, WAVENUM2

    open (UNIT=1, FILE=trim(adjustl(CKDFILE)), STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)

    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, N
        if (abs(WAVENUM1 - WAVENO(1)) .lt. 1.0 .and.&
            &abs(WAVENUM2 - WAVENO(2)) .lt. 1.0) then
            backspace (1)
            NG = N
            read (1, *) JB, WAVENUM1, WAVENUM2, N
        end if
    end do

    if (JB .eq. 0) then
        write (*, *) 'Error [_readCkdSize]: Wavenumber range not found in CKD file : ',&
            &WAVENO(1), WAVENO(2), CKDFILE
        stop
    end if

    read (1, *) NZCKD
    close (1)
    return

end subroutine shd_io__readCkdSize

subroutine shd_io__readCkd(&
    &MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,&
    &NG, DELG, NZCKD, ZCKD, KABS)

    ! Reads the information appropriate for one band from a
    ! correlated k-distribution file.   The wavenumber range (cm^-1)
    ! in the file must match the desired in WAVENO.  The solar flux
    ! is returned in SOLFLUX, the number of "g"'s in NG, the weights
    ! or delta g's in DELG, the number of levels in NZCKD, the
    ! levels from the top down in ZCKD, and the absorption coefficients
    ! for each g and level in KABS.  The array sizes MAXNG and MAXNZ
    ! are checked.

    integer MAXNG, MAXNZ, NG, NZCKD
    real WAVENO(2), SOLFLUX(*), DELG(*), ZCKD(*), KABS(*)
    character(len=*), intent(in) :: CKDFILE
    integer NB, IB, JB, KB, IG, I, J, N
    real WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=trim(adjustl(CKDFILE)), STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, N
        if (abs(WAVENUM1 - WAVENO(1)) .lt. 1.0 .and.&
            &abs(WAVENUM2 - WAVENO(2)) .lt. 1.0) then
            backspace (1)
            NG = N
            if (NG .gt. MAXNG) stop 'Error [_readCkd]: MAXNG exceeded.'
            read (1, *) JB, WAVENO(1), WAVENO(2),&
                &N, (SOLFLUX(IG), IG=1, NG), (DELG(IG), IG=1, NG)
        end if
    end do
    if (JB .eq. 0) then
        write (*, *) 'Error [_readCkd]: Wavenumber range not found in CKD file : ',&
            &WAVENO(1), WAVENO(2), CKDFILE
        stop
    end if
    ! Read the Z level heights
    read (1, *) NZCKD
    if (NZCKD .gt. MAXNZ) stop 'Error [_readCkd]: MAXNZ exceeded.'
    read (1, *)
    read (1, *)
    do I = 1, NZCKD
        read (1, *) ZCKD(I)
    end do
    ! Skip over the irrelevant bands and read the absorption coeffic
    read (1, *)
    do I = 1, (JB - 1)*NZCKD
        read (1, *)
    end do
    do I = 1, NZCKD
        read (1, *) IB, J, (KABS(I + (IG - 1)*NZCKD), IG=1, NG)
        do IG = 1, NG
            if (KABS(I + (IG - 1)*NZCKD) .lt. 0.0) then
                write (*, '(1X,A,A,I3,A,I3,A,I2)')&
                    &'Negative absorption in CKD file.  ',&
                    &'Band: ', IB, '   Level: ', I, '   k: ', IG
            end if
        end do
    end do
    close (1)
    return
end subroutine shd_io__readCkd

subroutine shd_io__readSenSize(&
    &SENFILE, NSENDATA, NSENPAR)

    ! Reads the user-specified sensor data size information

    character(len=*), intent(in) :: SENFILE
    integer, intent(out) :: NSENDATA, NSENPAR

    character(len=FNAMEMAXLEN_) :: TMP_CHAR
    character(len=FNAMEMAXLEN_) :: DIMS_CHAR

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=trim(adjustl(SENFILE)), STATUS='OLD')
    read (1, '(A)') TMP_CHAR ! skip one line
    read (1, '(A)') TMP_CHAR ! skip one line
    read (1, '(A)') TMP_CHAR ! skip one line
    read (1, '(A)') DIMS_CHAR

    read (DIMS_CHAR( 3:12), '(I10)') NSENPAR
    read (DIMS_CHAR(14:23), '(I10)') NSENDATA

    close (1)

    return

end subroutine shd_io__readSenSize

subroutine shd_io__readSen(&
    &SENFILE, NSENDATA, NSENPAR,&
    &SENPARMS)

    ! Reads the user-specified sensor information

    character(len=*), intent(in) :: SENFILE
    real, intent(inout) :: SENPARMS(*)
    integer :: NSENDATA, NSENPAR, ISENDATA

    character(len=FNAMEMAXLEN_) :: OUTFILE_DATA
    character(len=FNAMEMAXLEN_) :: TMP_CHAR
    character(len=FNAMEMAXLEN_) :: DIMS_CHAR
    character(len=20) :: DATA_POSTFIX

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=trim(adjustl(SENFILE)), STATUS='OLD')
    read (1, '(A)') TMP_CHAR ! skip one line
    read (1, '(A)') TMP_CHAR ! skip one line
    read (1, '(A)') DATA_POSTFIX
    read (1, '(A)') DIMS_CHAR

    read (DIMS_CHAR( 3:12), '(I10)') NSENPAR
    read (DIMS_CHAR(14:23), '(I10)') NSENDATA

    OUTFILE_DATA = trim(adjustl(SENFILE))//trim(adjustl(DATA_POSTFIX(2:)))

    open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='OLD',&
        &ACCESS='STREAM', FORM='UNFORMATTED')

    read (88) (SENPARMS(ISENDATA), ISENDATA=1,NSENDATA*NSENPAR)

    close (88)

    close (1)

    return

end subroutine shd_io__readSen

subroutine shd_io__restoreState(&
    &INSAVEFILE, NSTOKES, NX, NY, NZ,&
    &ML, MM, NLM, INRADFLAG, NEWGRIDFLAG, XGRID, YGRID, ZGRID,&
    &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
    &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

    ! Reads the binary save file containing a previous state of
    ! the SHDOM solution. This includes the cell tree structure,
    ! the spherical harmonic expansion of the source function,
    ! the low order terms of the radiance expansion, and the flux array.
    ! If the file name is 'NONE' then the INRADFLAG is set to false,
    ! otherwise it is true.  The input file may have any spherical
    ! harmonic truncation smaller than the current (ML,MM,NLM).
    ! The flag NEWGRIDFLAG is set to false if reading succeeds.

    character(*), intent(in) :: INSAVEFILE
    integer, intent(in) :: NSTOKES, NX, NY, NZ, ML, MM, NLM

    integer, intent(out) :: NPTS, NCELLS
    integer, intent(out) :: SHPTR(*), RSHPTR(*)
    integer, intent(out) :: GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer(2), intent(out) :: CELLFLAGS(*)
    logical, intent(out) :: INRADFLAG, NEWGRIDFLAG
    real, intent(out) :: XGRID(NX+1), YGRID(NY+1), ZGRID(NZ)
    real, intent(out) :: GRIDPOS(3, *), FLUXES(2, *), SOURCE(*), RADIANCE(*)

    integer :: NSTI, NXI, NYI, NZI, MLI, MMI, NLMI
    integer :: I, J, NS, NR

    INRADFLAG = .false.
    NEWGRIDFLAG = .true.
    if (INSAVEFILE(1:4) .eq. 'NONE') return

    open (UNIT=3, FILE=trim(adjustl(INSAVEFILE)), STATUS='OLD',&
        &FORM='UNFORMATTED', ERR=900)

    read (3) NSTI, NXI, NYI, NZI, MLI, MMI, NLMI, NPTS, NCELLS

    if (NXI .ne. NX .or. NYI .ne. NY .or. NZI .ne. NZ) then
        write (6, '(A, 3I6)') 'Warning [_restoreState]: wrong grid size:', NXI, NYI, NZI
        close (3)
        return
    end if

    if (NSTI .ne. NSTOKES .or.&
        &MLI .gt. ML .or. MMI .gt. MM) then
        write (6, '(A, 3I6)') 'Warning [_restoreState]: incompatible SH:', NSTI, MLI, MMI
        close (3)
        return
    end if

    read (3) (XGRID(I), I=1, NX+1)
    read (3) (YGRID(I), I=1, NY+1)
    read (3) (ZGRID(I), I=1, NZ)
    read (3) ((GRIDPTR(J, I), J=1, 8), I=1, NCELLS)
    read (3) ((NEIGHPTR(J, I), J=1, 6), I=1, NCELLS)
    read (3) ((TREEPTR(J, I), J=1, 2), I=1, NCELLS)
    read (3) (CELLFLAGS(I), I=1, NCELLS)

    read (3) ((GRIDPOS(J, I), J=1, 3), I=1, NPTS)
    read (3) ((FLUXES(J, I), J=1, 2), I=1, NPTS)

    read (3) (SHPTR(I), I=1, NPTS + 1)
    NS = SHPTR(NPTS + 1)
    read (3) (SOURCE(I), I=1, NSTOKES*NS)

    NR = 4*NPTS
    read (3) (RADIANCE(I), I=1, NSTOKES*NR)
    RSHPTR(1) = 0
    do I = 1, NPTS
        RSHPTR(I+1) = RSHPTR(I) + NSTOKES*4
    end do

    close (3)

    write (6, '(A)') 'Message [_restoreState]: State restored from <'//trim(adjustl(INSAVEFILE))//'> ...'

    INRADFLAG = .true.
    NEWGRIDFLAG = .false.
    return

900 continue

    write (6, '(A)') 'Warning [_restoreState]: Failed to restore state file from <'//trim(adjustl(INSAVEFILE))//'>.'
    INRADFLAG = .false.
    NEWGRIDFLAG = .true.

    return

end subroutine shd_io__restoreState

end module shd_inp
