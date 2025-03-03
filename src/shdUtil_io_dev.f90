subroutine USER_INPUT( &
    & RUNNAME, PROPFILE, SFCFILE, CKDFILE, &
    & INSAVEFILE, OUTSAVEFILE, &
    & NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, &
    & IPFLAG, KDIST, DELTAM, GRIDTYPE, &
    & SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
    & GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN, &
    & ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC, &
    & MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES, &
    & OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR, &
    & NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

    ! Obtains the input parameters for the program by writing prompts
    ! to and reading responses from the terminal.  Echos the inputs to
    ! make a useful log file when running non-interactively.  See the ov
    ! program documentation for the list of input parameters.
    implicit none
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
    character SRCTYPE*1, GRIDTYPE*1, UNITS*1
    character OUTTYPES(*)*1, OUTFILES(*)*256
    character*256 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 OutFileNC, RUNNAME
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
    write (*, '(1X,A)') &
        & 'Correlated k-distribution file name (or NONE)'
    read (*, '(A)') CKDFILE
    write (*, *) CKDFILE
    KDIST = CKDFILE(1:2) .ne. 'NO'
    write (*, '(1X,A)') 'Input binary save file name'
    read (*, '(A)') INSAVEFILE
    write (*, *) INSAVEFILE
    write (*, '(1X,A)') 'Output binary save file name'
    read (*, '(A)') OUTSAVEFILE
    write (*, *) OUTSAVEFILE

    write (*, '(1X,A)') &
        & 'Number of Stokes parameters (1, 3, or 4; 1=unpolarized)'
    read (*, *) NSTOKES
    write (*, *) NSTOKES
    write (*, '(1X,A)') 'Internal base grid size (Nx, Ny, Nz)'
    read (*, *) NX, NY, NZ
    write (*, *) NX, NY, NZ
    write (*, '(1X,A)') 'Number of discrete ordinates in mu and phi'
    read (*, *) NMU, NPHI
    write (*, *) NMU, NPHI
    write (*, '(1X,A)') &
        & 'Boundary conditions (0=periodic, 3=open)'
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

    write (*, '(1X,A)') &
        & 'Cell splitting and spherical harmonic series accuracy'
    read (*, *) SPLITACC, SHACC
    write (*, *) SPLITACC, SHACC

    write (*, '(1X,A,A)') 'Acceleration flag (T or F), ', &
        & 'solution accuracy, and max iterations'
    read (*, *) ACCELFLAG, SOLACC, MAXITER
    write (*, *) ACCELFLAG, SOLACC, MAXITER

    write (*, '(1X,A)') 'Number of output files'
    read (*, *) NUMOUT
    write (*, *) NUMOUT
    if (NUMOUT .gt. MAXOUT) stop 'USER_INPUT: MAXOUT exceeded'
    do I = 1, NUMOUT
        write (*, '(1X,A,A)') 'Output type (R-radiance, ', &
            & 'V-visualization, F-flux, H-heating,'
        write (*, '(14X,A)') 'S-spherical harmonic, J-source, M-medium)'
        read (*, '(A)') OUTTYPES(I)
        write (*, *) OUTTYPES(I)
        if (OUTTYPES(I) (1:1) .eq. 'R') then
            write (*, '(1X,A,A)') &
                & 'Z level, X and Y spacing, X and Y offset, ', &
                & 'number of angles, and mu, phi pairs'
            read (*, *) OUTPARMS(1, I), OUTPARMS(2, I), OUTPARMS(3, I), &
                & OUTPARMS(4, I), OUTPARMS(5, I), &
                & N, (OUTPARMS(J + 6, I), J=1, min(2*N, MAXPAR - 6))
            if (N .gt. (MAXPAR - 6)/2) then
                print *, 'Number of radiance angles is too large for MAXPAR'
                stop
            end if
            N = min(N, (MAXPAR - 6)/2)
            OUTPARMS(6, I) = FLOAT(N)
            write (*, *) (OUTPARMS(J, I), J=1, 2*N + 6)
        else if (OUTTYPES(I) (1:1) .eq. 'V') then
            write (*, '(1X,A,A)') &
                & '1, bytes,scale, X,Y,Z, theta,phi,rotang, NL,NS, delline,delsamp'
            write (*, '(1X,A,A)') &
                & '2, bytes,scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1,scan2,delscan'
            read (*, *) (OUTPARMS(J, I), J=1, 13)
            write (*, *) (OUTPARMS(J, I), J=1, 13)
        else if (OUTTYPES(I) (1:1) .eq. 'F') then
            write (*, '(1X,A)') &
                & 'Format (1,2,3,4,5), [Z, del X, del Y for format 2]'
            read (*, *) N, (OUTPARMS(1 + J, I), J=1, 3*max(0, 1 - abs(N - 2)))
            OUTPARMS(1, I) = FLOAT(N)
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

    write (*, '(1X,A)') &
        & 'Approximate maximum memory to use (MB for 4 byte reals)'
    read (*, *) MAX_TOTAL_MB
    write (*, *) MAX_TOTAL_MB
    write (*, '(1X,A)') &
        & 'Ratio of total grid points to base grid points (1.0 to ~10)'
    read (*, *) ADAPT_GRID_FACTOR
    write (*, *) ADAPT_GRID_FACTOR
    write (*, '(1X,A,A)') 'Ratio of average number of spherical', &
        & ' harmonic terms to total possible (<=1.0)'
    read (*, *) NUM_SH_TERM_FACTOR
    write (*, *) NUM_SH_TERM_FACTOR
    write (*, '(1X,A)') &
        & 'Ratio of number of grid cells to grid points (1.0 to 2.0)'
    read (*, *) CELL_TO_POINT_RATIO
    write (*, *) CELL_TO_POINT_RATIO

    write (*, *)
    return
end subroutine

subroutine NAMELIST_INPUT ( &
    & dRUNNAME, dPROPFILE, dSFCFILE, dCKDFILE, &
    & dINSAVEFILE, dOUTSAVEFILE, &
    & dNSTOKES, dNX, dNY, dNZ, dNMU, dNPHI, &
    & dBCFLAG, dIPFLAG, &
    & dKDIST, dDELTAM, dGRIDTYPE, &
    & dSRCTYPE, dSOLARFLUX, dSOLARMU, dSOLARAZ, dSKYRAD, &
    & dGNDTEMP, dGNDALBEDO, dUNITS, dWAVENO, dWAVELEN, &
    & dACCELFLAG, dSOLACC, dMAXITER, dSPLITACC, dSHACC, &
    & dMAXOUT, dMAXPAR, dNUMOUT, dOUTTYPES, dOUTPARMS, dOUTFILES, &
    & dOutFileNC, dMAX_TOTAL_MB, dADAPT_GRID_FACTOR, &
    & dNUM_SH_TERM_FACTOR, dCELL_TO_POINT_RATIO)
    ! Obtains the input parameters for the program from a namelist
    ! read from stdin.  See the overall program documentation for the
    ! list of input parameters.  The subroutine dummy variables can't
    ! be in the namelist, hence the two sets of variables.
    implicit none
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
    character dSRCTYPE*1, dGRIDTYPE*1, dUNITS*1
    character dOUTTYPES(*)*1, dOUTFILES(*)*256
    character*256 dOutFileNC, dRUNNAME
    character*256 dPROPFILE, dSFCFILE, dCKDFILE
    character*256 dINSAVEFILE, dOUTSAVEFILE
    integer I, J, N

    integer MAXOUT, MAXPAR
    parameter(MAXOUT=20, MAXPAR=500)
    integer NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    integer MAXITER, NUMOUT
    logical DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real WAVENO(2), WAVELEN, SOLACC, SPLITACC, SHACC
    real MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real OUTPARMS(MAXPAR, MAXOUT)
    character SRCTYPE*1, GRIDTYPE*1, UNITS*1
    character OUTTYPES(MAXOUT)*1, OUTFILES(MAXOUT)*256
    character*256 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 OutFileNC, RUNNAME
    namelist /SHDOMINPUT/ &
        & RUNNAME, PROPFILE, SFCFILE, CKDFILE, &
        & INSAVEFILE, OUTSAVEFILE, &
        & NSTOKES, NX, NY, NZ, NMU, NPHI, &
        & BCFLAG, IPFLAG, DELTAM, GRIDTYPE, &
        & SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
        & GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN, &
        & ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC, &
        & NUMOUT, OUTTYPES, OUTPARMS, OUTFILES, OutFileNC, &
        & MAX_TOTAL_MB, ADAPT_GRID_FACTOR, &
        & NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

    character*256 NMLFILE

    if (dMAXOUT .ne. MAXOUT) stop 'NAMELIST_INPUT: MAXOUT not the same.'
    if (dMAXPAR .ne. MAXPAR) stop 'NAMELIST_INPUT: MAXPAR not the same.'

    SOLARMU = -1.0
    GNDALBEDO = 0.0
    WAVELEN = 1.0
    WAVENO(1) = 10000.
    WAVENO(2) = 10001.

    ! Read in the namelist
    write (*, '(1X,A)') 'Namelist file name'
    read (*, '(A)') NMLFILE
    open (UNIT=5, FILE=NMLFILE, STATUS='OLD')
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
    if (NUMOUT .gt. MAXOUT) stop 'NAMELIST_INPUT: MAXOUT exceeded'
    do I = 1, NUMOUT
        dOUTTYPES(I) = OUTTYPES(I)
        if (OUTTYPES(I) (1:1) .eq. 'R') then
            N = nint(OUTPARMS(6, I))
            do J = 1, min(2*N + 6, MAXPAR)
                dOUTPARMS(J, I) = OUTPARMS(J, I)
            end do
        else if (OUTTYPES(I) (1:1) .eq. 'V') then
            do J = 1, 13
                dOUTPARMS(J, I) = OUTPARMS(J, I)
            end do
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

    dMAX_TOTAL_MB = MAX_TOTAL_MB
    dADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR
    dNUM_SH_TERM_FACTOR = NUM_SH_TERM_FACTOR
    dCELL_TO_POINT_RATIO = CELL_TO_POINT_RATIO

    return
end subroutine

subroutine CHECK_INPUT_PARAMETERS(NMU, NPHI, DELTAM, MAXASYM,    &
    &               GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU,      &
    &               WAVENO, WAVELEN, GNDTEMP, GNDALBEDO,               &
    &               SPLITACC, SHACC)
    ! Checks for validity of the input parameters. Gives warning for
    ! parameters that might not be correct.
    ! Check the letter flags.
    ! Delta-M should be used for solar transfer with forward peaked
    ! phase functions.
    ! Check for reasonable ranges on SPLITACC and SHACC.
    ! Check for reasonable range of GNDTEMP and valid range of GNDALBEDO
    ! Check for SOLARMU=0.
    ! Check for NPHI<NMU.
    implicit none
    integer NMU, NPHI
    logical DELTAM
    real MAXASYM, SOLARFLUX, SOLARMU, WAVENO(2), WAVELEN
    real GNDTEMP, GNDALBEDO, SPLITACC, SHACC
    character*1 GRIDTYPE, SRCTYPE, UNITS
    real FLUX

    if (SRCTYPE .ne. 'S' .and. SRCTYPE .ne. 'T' .and. SRCTYPE .ne. 'B')  &
   &    stop 'CHECK_INPUT_PARMETERS: Illegal source type'
    if (UNITS .ne. 'R' .and. UNITS .ne. 'T' .and. UNITS .ne. 'B')     &
   &    stop 'CHECK_INPUT_PARMETERS: Illegal units'
    if (GRIDTYPE .ne. 'E' .and. GRIDTYPE .ne. 'P' .and. GRIDTYPE .ne. 'F')  &
   &    stop 'CHECK_INPUT_PARMETERS: Illegal grid type'

    if (SRCTYPE .ne. 'T' .and. SOLARMU .eq. 0.0) then
        write (*, *) 'CHECK_INPUT_PARMETERS: cosine solar zenith angle'
        write (*, *) ' (SOLARMU) cannot be zero.'
        stop
    end if

    if (.not. DELTAM .and. SRCTYPE .ne. 'T' .and. MAXASYM .gt. 0.5) then
        write (*, *) 'CHECK_INPUT_PARMETERS: Delta-M should be used'
        write (*, *) ' (DELTAM=.TRUE.) for solar radiative transfer'
        write (*, *) ' with highly peaked phase functions.'
        write (*, *)
    end if

    ! Check for SPLITACC and SHACC magnitudes
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        call PLANCK_FUNCTION(GNDTEMP, UNITS, WAVENO, WAVELEN, FLUX)
        FLUX = 3.14*FLUX
    else
        FLUX = 0.0
    end if
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        FLUX = FLUX + SOLARFLUX
    end if
    if (SPLITACC .gt. 0.0 .and.                                       &
   &    (SPLITACC .lt. 0.001*FLUX .or. SPLITACC .gt. 0.1*FLUX)) then
        write (*, *) 'CHECK_INPUT_PARMETERS: Splitting accuracy parameter'
        write (*, *) ' is not a fraction, but scales with fluxes'
        write (*, *)
    end if
    if (SHACC .gt. 0.0 .and. SHACC .gt. 0.03*FLUX) then
        write (*, *) 'CHECK_INPUT_PARMETERS: Spherical harmonic accuracy'
        write (*, *) ' (SHACC) is not a fraction, but scales with fluxes'
        write (*, *)
    end if

    if (SRCTYPE .ne. 'S' .and.                                        &
   &    (GNDTEMP .lt. 150. .or. GNDTEMP .gt. 350.)) then
        write (*, *) 'CHECK_INPUT_PARMETERS: Warning - '
        write (*, *) '  ground temperature out of Earth range.'
        write (*, *)
    end if
    if (GNDALBEDO .lt. 0.0 .or. GNDALBEDO .gt. 1.0) then
        write (*, *) 'CHECK_INPUT_PARMETERS: Ground albedo must be 0 to 1'
        stop
    end if

    if (NPHI .lt. NMU) then
        write (*, *) 'CHECK_INPUT_PARMETERS: Usually want NPHI >= NMU'
        write (*, *)
    end if

    return
end

subroutine READ_PROPERTY_SIZE(PROPFILE, NSTLEG, NLEG, NPX, NPY, NPZ, &
&                            NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
    ! Reads parts of the property file to get the maximum array sizes
    ! needed for allocatable arrays.  For extinction only and tabulated
    ! phase function formats, only the header is read, while for the
    ! standard format, the whole file must be read to determine MAXLEG.
    implicit none
    integer NSTLEG, NPX, NPY, NPZ
    integer NUMPHASE, NLEG, MAXLEG, MAXPGL
    real DELX, DELY
    character PROPFILE*256
    integer NUML, IX, IY, IZ, I, J, K, L
    real ZLEVELS, TMP, EXT, ALB, CHI
    character PROPTYPE*1

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=PROPFILE, STATUS='OLD')
    read (1, '(A1)') PROPTYPE
    if (PROPTYPE .ne. 'E' .and. PROPTYPE .ne. 'T'                     &
   &    .and. PROPTYPE .ne. 'P') rewind (1)
    if (PROPTYPE .ne. 'P' .and. NSTLEG .gt. 1) then
        write (*, '(A,A)') 'Error: NSTOKES>1 requires ',                 &
     &                    'Polarized phase function property file'
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
            read (1, *, end=190) IX, IY, IZ, TMP, EXT, ALB, NUML,           &
       &                         (CHI, L=1, NUML)
            MAXLEG = max(NUML, MAXLEG)
        end do
190     continue
        MAXPGL = NPX*NPY*NPZ*(MAXLEG + 1)
    end if
    close (1)
    return
end

subroutine READ_PROPERTIES(PROPFILE, NPX, NPY, NPZ,              &
    &             NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,  &
    &             PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM,              &
    &             TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
    ! Reads the medium properties from the file into the property arra
    ! There are four formats: the standard format specifies the extincti
    ! scattering albedo, temperature, and Legendre phase function expans
    ! at every grid point; the extinction only format only specifies
    ! extinction at every grid point, while the other quantities are giv
    ! in the header part; and the tabulated and polarized phase function
    ! format specifies all properties except the phase function for each
    ! grid point, while the phase functions are tabulated in the header.
    ! The polarized format specifies the six unique elements of the phas
    ! (Mueller) matrix in the header and includes the l=0 term.
    ! If the format is extinction only or tabulated phase function then
    ! NUMPHASE is the number of phase functions, otherwise NUMPHASE=0.
    ! If doing delta-M then NLEG is the minimum number of Legendre terms
    ! input and the actual maximum number of terms on output, except tha
    ! it may not exceed MAXLEG; otherwise NLEG is the number of Legendre
    ! terms to be used (regardless of what is in property file).
    ! See the overall program documentation for the four file formats.
    implicit none
    integer NPX, NPY, NPZ, NSTLEG, MAXLEG, NLEG, NUMPHASE
    integer MAXNZ, MAXPG, MAXPGL
    integer*2 IPHASEP(*)
    logical DELTAM
    real DELX, DELY, ZLEVELS(*)
    real TEMPP(*), EXTINCTP(*)
    real ALBEDOP(*), LEGENP(*), MAXASYM
    character PROPFILE*256
    character PROPTYPE*1
    integer I, IX, IY, IZ, J, JS, K, KL, L, IL, NUML, IPH, N
    real TMP, EXT, ALB, CHI

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=PROPFILE, STATUS='OLD')
    read (1, '(A1)') PROPTYPE
    if (PROPTYPE .ne. 'E' .and. PROPTYPE .ne. 'T'                     &
   &    .and. PROPTYPE .ne. 'P') rewind (1)
    if (PROPTYPE .ne. 'P' .and. NSTLEG .gt. 1) then
        write (*, '(A,A)') 'Error: NSTOKES>1 requires ',                 &
     &                    'Polarized phase function property file'
        stop
    end if
    read (1, *) NPX, NPY, NPZ
    if (NPZ .gt. MAXNZ) stop 'READ_PROPERTIES: MAXNZ exceeded'
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
        if (NPX*NPY*NPZ .gt. MAXPG)                                     &
     &    stop 'READ_PROPERTIES: MAXPG exceeded'
        if (NLEG + 1 .gt. MAXPGL)                                         &
     &    stop 'READ_PROPERTIES: MAXPGL exceeded'
        rewind (1)
        read (1, *)
        read (1, *) NPX, NPY, NPZ
        read (1, *) DELX, DELY, (ZLEVELS(K), K=1, NPZ)
        read (1, *) (TEMPP(IZ), IZ=1, NPZ)
        read (1, *) ALB, NUML, (LEGENP(L + 1), L=1, min(NUML, NLEG)),        &
     &                        (CHI, L=NLEG + 1, NUML)
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
            if (IX .ge. 1 .and. IX .le. NPX .and.                         &
       &        IY .ge. 1 .and. IY .le. NPY .and.                         &
       &        IZ .ge. 1 .and. IZ .le. NPZ) then
                K = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
                TEMPP(K) = TEMPP(IZ)
                EXTINCTP(K) = EXT
                ALBEDOP(K) = ALB
                IPHASEP(K) = 1
            end if
        end do

    ! Property file type T or P is for tabulated phase function fo
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
                write (*, '(A,I4,A,I4,A)')                                   &
         &        'Warning: truncating phase function from ',               &
         &         NLEG, ' to ', MAXLEG, ' Legendre terms.'
            end if
            NLEG = min(NLEG, MAXLEG)
        end if
        N = NPX*NPY*NPZ
        if (N .gt. MAXPG) stop 'READ_PROPERTIES: MAXPG exceeded'
        if (NUMPHASE*(NLEG + 1)*NSTLEG .gt. MAXPGL)                       &
     &    stop 'READ_PROPERTIES: MAXPGL exceeded'
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
                read (1, *) NUML, (LEGENP(IL + L + 1), L=1, min(NUML, NLEG)),      &
         &                     (CHI, L=NLEG + 1, NUML)
                do L = NUML + 1, NLEG
                    LEGENP(L + IL + 1) = 0.0
                end do
            else
                do J = 1, 6
                    if (J .le. NSTLEG) then
                        read (1, *) JS, NUML,                                    &
             &            (LEGENP((IL + L)*NSTLEG + J), L=0, min(NUML, NLEG)),        &
             &                     (CHI, L=NLEG + 1, NUML)
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

        do while (.true.)
            read (1, *, end=290) IX, IY, IZ, TMP, EXT, ALB, IPH
            if (IX .ge. 1 .and. IX .le. NPX .and.                         &
       &        IY .ge. 1 .and. IY .le. NPY .and.                         &
       &        IZ .ge. 1 .and. IZ .le. NPZ) then
                K = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
                TEMPP(K) = TMP
                EXTINCTP(K) = EXT
                ALBEDOP(K) = ALB
                IPHASEP(K) = IPH
            end if
        end do

    else

        ! Standard property file has everything variable
        !   If delta-M then find the largest number of Legendre terms
        NUMPHASE = 0
        if (DELTAM) then
            do while (.true.)
                read (1, *, end=190) IX, IY, IZ, TMP, EXT, ALB, NUML,         &
         &                         (CHI, L=1, NUML)
                NLEG = max(NUML, NLEG)
            end do
190         continue
            NLEG = min(NLEG, MAXLEG)
            rewind (1)
            read (1, *) NPX, NPY, NPZ
            read (1, *) DELX, DELY, (ZLEVELS(K), K=1, NPZ)
        end if
        N = NPX*NPY*NPZ
        if (N .gt. MAXPG) stop 'READ_PROPERTIES: MAXPG exceeded'
        if (N*(NLEG + 1) .gt. MAXPGL)                                     &
     &    stop 'READ_PROPERTIES: MAXPGL exceeded'
        MAXASYM = -1.0
        do while (.true.)
            read (1, *, end=290) IX, IY, IZ, TMP, EXT, ALB, NUML
            if (IX .ge. 1 .and. IX .le. NPX .and.                         &
       &        IY .ge. 1 .and. IY .le. NPY .and.                         &
       &        IZ .ge. 1 .and. IZ .le. NPZ) then
                K = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
                KL = (NLEG + 1)*(K - 1)
                backspace (1)
                read (1, *, end=290) IX, IY, IZ, TEMPP(K), EXTINCTP(K),     &
       &            ALBEDOP(K), NUML, (LEGENP(L + KL), L=1, min(NUML, NLEG)), &
       &                              (CHI, L=NLEG + 1, NUML)
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
end

subroutine CHECK_PROPERTY_INPUT(NPX, NPY, NPZ, NSTLEG, NLEG,     &
    &             DELX, DELY, ZLEVELS,                                 &
    &             TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
    ! Checks the input medium property values for correctness.
    ! Warning if temperature below 150 K or above 350 K.
    ! Stop if negative extinction, single scattering albedo below 0 or
    ! above 1, asymmetry parameter below -1 or above 1, or IPHASE
    ! out of range (of tabulated phase functions).
    ! Warning if Chi1=1, for people specifying Chi0.
    ! Warning if optical depth across a grid cell exceeds 2.
    implicit none
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
        write (*, *) 'CHECK_PROPERTY_INPUT: Bad grid size - ',           &
     &         NPX, NPY, NPZ
        stop
    end if
    if (DELX .le. 0.0 .or. DELY .le. 0.0) then
        write (*, *) 'CHECK_PROPERTY_INPUT: Bad grid spacing - ',        &
     &         'DELX,DELY: ', DELX, DELY
        stop
    end if
    ! The vertical levels must be increasing
    do IZ = 1, NPZ - 1
        if (ZLEVELS(IZ + 1) .le. ZLEVELS(IZ)) then
            write (*, *) 'CHECK_PROPERTY_INPUT: ZLEVELS must increase : ', &
       &         IZ, ZLEVELS(IZ), ZLEVELS(IZ + 1)
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
        write (*, *) 'CHECK_PROPERTY_INPUT: Warning - ',                 &
     &        'temperature min and max (Kelvin): ', TMIN, TMAX
        write (*, *)
    end if
    do I = 1, N
        if (EXTINCTP(I) .lt. 0.0) then
            write (*, *) 'CHECK_PROPERTY_INPUT: Negative extinction - ',   &
       &        EXTINCTP(I), ' at grid point ', I
            stop
        end if
        if (ALBEDOP(I) .lt. 0.0 .or. ALBEDOP(I) .gt. 1.0) then
            write (*, *) 'CHECK_PROPERTY_INPUT: Single scattering ',       &
       &      'albedo out of range - ', ALBEDOP(I), ' at grid point ', I
            stop
        end if
    end do

    if (NUMPHASE .gt. 0) N = NUMPHASE
    LEGWARN = .false.
    do I = 1, N
        if (LEGENP(1, 1, I) .lt. -3.0 .or. LEGENP(1, 1, I) .gt. 3.0) then
            write (*, *) 'CHECK_PROPERTY_INPUT: Asymmetry parameter ',     &
       &      'out of range - ', LEGENP(1, 1, I)/3.0, ' at grid point ', I
            stop
        end if
        if (LEGENP(1, 1, I) .eq. 1.0) LEGWARN = .true.
    end do
    if (LEGWARN) then
        write (*, *) 'CHECK_PROPERTY_INPUT: Warning - ',                 &
     &        'an asymmetry parameter is 1/3 (Chi1=1); ',               &
     &        '  Chi0 should not be specified for Chi1.'
        write (*, *)
    end if

    if (NUMPHASE .gt. 0) then
        N = NPX*NPY*NPZ
        do I = 1, N
            if (IPHASEP(I) .lt. 1 .or. IPHASEP(I) .gt. NUMPHASE) then
                write (*, *) 'CHECK_PROPERTY_INPUT: Phase function pointer', &
         &      ' out of range - ', IPHASEP(I), ' at grid point ', I
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
        write (*, *) 'CHECK_PROPERTY_INPUT: Warning -'
        write (*, '(A,A,I7)') '  Number of property grid cells',         &
     &     ' with optical depth greater than 2 : ', NTHICK
        write (*, *) ' Max cell optical depth : ', MAXTAU
        write (*, *)
    end if

    return
end

subroutine READ_SURFACE_SIZE(SFCFILE, MAXSFCPTS, MAXSFCPARS)
    ! Gets the size of the surface property inputs (temperature is
    ! always included).
    implicit none
    integer MAXSFCPTS, MAXSFCPARS
    character SFCFILE*256
    integer NXSFC, NYSFC
    character SFCTYPE*2

    ! Open the file, figure out the type, and get the grid size
    SFCTYPE(1:1) = 'V'
    open (UNIT=1, FILE=SFCFILE, STATUS='OLD')
    read (1, '(A1)') SFCTYPE(2:2)
    read (1, *) NXSFC, NYSFC
    MAXSFCPTS = (NXSFC + 1)*(NYSFC + 1)

    if (SFCTYPE .eq. 'VL' .or. SFCTYPE .eq. 'Vl') then
        ! Surface file type L is for variable Lambertian surface
        MAXSFCPARS = 2
    else if (SFCTYPE .eq. 'VW') then
        ! Surface file type W is for Wave-Fresnel surface
        MAXSFCPARS = 4
    else if (SFCTYPE .eq. 'VD') then
        ! Surface file type D is for Diner polarized reflection model
        MAXSFCPARS = 6
    else if (SFCTYPE .eq. 'VO') then
        ! Surface file type O is for Ocean BRDF surface
        MAXSFCPARS = 3
    else if (SFCTYPE .eq. 'VR') then
        ! Surface file type R is for RPV surface
        MAXSFCPARS = 4
    else if (SFCTYPE .eq. 'VD') then
        ! Surface file type D is for Diner polarized RPV surface:
        MAXSFCPARS = 6
    else
        stop 'READ_SURFACE_SIZE: Unknown BRDF type'
    end if
    close (1)
    return
end

subroutine READ_SURFACE(SFCFILE, MAXSFCPTS, MAXSFCPARS,          &
    &                     SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,     &
    &                     NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
    ! Reads the surface properties from the file into the surface arra
    ! There are currently six surface types (Lambertian, Wave-Fresnel,
    ! Ocean, RPV, and Diner polarized RPV), though it is simple to add m
    ! In all cases the surface properties are specified on a regular, ev
    ! spaced grid of size NXSFC by NYSFC and spacing DELXSFC by DELYSFC.
    ! The properties are returned on a NXSFC+1 by NYSFC+1 grid that is
    ! periodic (NX+1 same as 1, etc).  The surface properties at each
    ! grid point are returned in SFCPARMS, and the first parameter must
    ! be temperature.  The number of parameters for each grid point
    ! (including temperature) are returned in NSFCPAR.
    ! Also returned are the domain average ground albedo and temperature
    ! (GNDALBEDO, GNDTEMP) for use with initialization (hence they can b
    ! highly approximate).
    !        Type        Parameters
    !   L  Lambertian    albedo
    !   W  Wave-Fresnel  Real, Imaginary index of refraction, wind speed
    !   O  Ocean         Wind speed (m/s), Pigmentation (mg/m^3)
    !   D  Diner et al   a, k, b, zeta, sigma
    !   R  RPV           rho0, k, Theta
    !
    implicit none
    integer MAXSFCPTS, MAXSFCPARS, NXSFC, NYSFC, NSFCPAR
    real DELXSFC, DELYSFC
    real SFCPARMS(*), GNDTEMP, GNDALBEDO
    character SFCFILE*256, SFCTYPE*2
    integer N, I, I0, IX, IY, J
    real ALB, TEMP, MRE, MIM, RHO0, KR, THETA, WSPD, PCL
    real A, B, ZETA, SIGMA

    ! Surface type is variable something
    SFCTYPE(1:1) = 'V'
    ! Open the file, figure out the type, and get the grid size
    open (UNIT=1, FILE=SFCFILE, STATUS='OLD')
    read (1, '(A1)') SFCTYPE(2:2)
    read (1, *) NXSFC, NYSFC, DELXSFC, DELYSFC
    if ((NXSFC + 1)*(NYSFC + 1) .gt. MAXSFCPTS) then
        write (*, *) 'READ_SURFACE: MAXSFCPTS exceeded'
        stop
    end if

    GNDTEMP = 0.0
    GNDALBEDO = 0.0
    N = 0

    ! Surface file type L is for variable Lambertian surface
    if (SFCTYPE .eq. 'VL' .or. SFCTYPE .eq. 'Vl') then
        NSFCPAR = 2
        if (NSFCPAR .gt. MAXSFCPARS) then
            write (*, *) 'READ_SURFACE: MAXSFCPARS exceeded'
            stop
        end if
        do while (.true.)
            read (1, *, end=290) IX, IY, TEMP, ALB
            if (IX .ge. 1 .and. IX .le. NXSFC .and.                       &
       &        IY .ge. 1 .and. IY .le. NYSFC) then
                I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
                SFCPARMS(I + 1) = TEMP
                SFCPARMS(I + 2) = ALB
                if (ALB .lt. 0.0 .or. ALB .gt. 1.0)                         &
         &        stop 'READ_SURFACE: Illegal surface albedo'
                if (TEMP .lt. 150. .or. TEMP .gt. 350.)                     &
         &        write (*, *) 'READ_SURFACE: Warning -',                    &
         &          ' surface temperature: ', IX, IY, TEMP
                GNDTEMP = GNDTEMP + TEMP
                GNDALBEDO = GNDALBEDO + ALB
                N = N + 1
            end if
        end do

    ! Surface file type W is for Wave-Fresnel surface
    else if (SFCTYPE .eq. 'VW') then
        NSFCPAR = 4
        if (NSFCPAR .gt. MAXSFCPARS) then
            write (*, *) 'READ_SURFACE: MAXSFCPARS exceeded'
            stop
        end if
        do while (.true.)
            read (1, *, end=290) IX, IY, TEMP, MRE, MIM, WSPD
            if (IX .ge. 1 .and. IX .le. NXSFC .and.                       &
       &        IY .ge. 1 .and. IY .le. NYSFC) then
                I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
                SFCPARMS(I + 1) = TEMP
                SFCPARMS(I + 2) = MRE
                SFCPARMS(I + 3) = MIM
                SFCPARMS(I + 4) = WSPD
                GNDTEMP = GNDTEMP + TEMP
                N = N + 1
            end if
        end do

    ! Surface file type O is for Ocean surface BRDF
    !   Norm Loeb's modification of 6S ocean reflectance module.
    else if (SFCTYPE .eq. 'VO') then
        NSFCPAR = 3
        if (NSFCPAR .gt. MAXSFCPARS) then
            write (*, *) 'READ_SURFACE: MAXSFCPARS exceeded'
            stop
        end if
        do while (.true.)
            read (1, *, end=290) IX, IY, TEMP, WSPD, PCL
            if (IX .ge. 1 .and. IX .le. NXSFC .and.                       &
       &        IY .ge. 1 .and. IY .le. NYSFC) then
                I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
                SFCPARMS(I + 1) = TEMP
                SFCPARMS(I + 2) = WSPD
                SFCPARMS(I + 3) = PCL
                GNDTEMP = GNDTEMP + TEMP
                N = N + 1
            end if
        end do

    ! Surface file type R is for RPV surface:
    !   (Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphe
    !   Reflectance (CSAR) Model. 2. Semiempirical Surface Model
    !   Usable With NOAA Advanced Very High Resolution Radiometer
    !   Data, J. Geophys. Res., 98, 20791-20801.)
    else if (SFCTYPE .eq. 'VR') then
        NSFCPAR = 4
        if (NSFCPAR .gt. MAXSFCPARS) then
            write (*, *) 'READ_SURFACE: MAXSFCPARS exceeded'
            stop
        end if
        do while (.true.)
            read (1, *, end=290) IX, IY, TEMP, RHO0, KR, THETA
            if (IX .ge. 1 .and. IX .le. NXSFC .and.                       &
       &        IY .ge. 1 .and. IY .le. NYSFC) then
                I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
                SFCPARMS(I + 1) = TEMP
                SFCPARMS(I + 2) = RHO0
                SFCPARMS(I + 3) = KR
                SFCPARMS(I + 4) = THETA
                GNDTEMP = GNDTEMP + TEMP
                GNDALBEDO = GNDALBEDO + RHO0
                N = N + 1
            end if
        end do

    ! Surface file type D is for Diner et al surface reflection:
    ! (Diner, D. J., F. Xu, J. V. Martonchik, B. E. Rheingans, S. Geier,
    ! V. M. Jovanovic, A. Davis, R. A. Chipman, S. C. McClain, 2012:
    ! Exploration of a polarized surface bidirectional reflectance model
    ! using the ground-based multiangle spectropolarimetric imager.
    ! Atmosphere 2012, 3, 591-619; doi:10.3390/atmos3040591.)
    else if (SFCTYPE .eq. 'VD') then
        NSFCPAR = 6
        if (NSFCPAR .gt. MAXSFCPARS) then
            write (*, *) 'READ_SURFACE: MAXSFCPARS exceeded'
            stop
        end if
        do while (.true.)
            read (1, *, end=290) IX, IY, TEMP, A, KR, B, ZETA, SIGMA
            if (IX .ge. 1 .and. IX .le. NXSFC .and.                       &
       &        IY .ge. 1 .and. IY .le. NYSFC) then
                I = NSFCPAR*(IX - 1 + (NXSFC + 1)*(IY - 1))
                SFCPARMS(I + 1) = TEMP
                SFCPARMS(I + 2) = A
                SFCPARMS(I + 3) = KR
                SFCPARMS(I + 4) = B
                SFCPARMS(I + 5) = ZETA
                SFCPARMS(I + 6) = SIGMA
                GNDTEMP = GNDTEMP + TEMP
                GNDALBEDO = GNDALBEDO + A
                N = N + 1
            end if
        end do

    else
        stop 'READ_SURFACE: Unknown BRDF type'
    end if

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
    GNDALBEDO = GNDALBEDO/N
    GNDTEMP = GNDTEMP/N
    return
end

subroutine READ_CKD_SIZE(CKDFILE, WAVENO, NG, NZCKD)
    ! Find the size of the k-distribution arrays for this wavenumber b
    implicit none
    integer NG, NZCKD
    real WAVENO(2)
    character CKDFILE*256
    integer NB, IB, JB, KB, N
    real WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=CKDFILE, STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, SF, N
        if (abs(WAVENUM1 - WAVENO(1)) .lt. 1.0 .and.                     &
     &       abs(WAVENUM2 - WAVENO(2)) .lt. 1.0) then
            backspace (1)
            NG = N
            read (1, *) JB, WAVENUM1, WAVENUM2, SF, N
        end if
    end do
    if (JB .eq. 0) then
        write (*, *) ' Wavenumber range not found in CKD file : ',       &
     &              WAVENO(1), WAVENO(2), CKDFILE
        stop
    end if
    read (1, *) NZCKD
    close (1)
    return
end

subroutine READ_CKD(MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,      &
    &                     NG, DELG, NZCKD, ZCKD, KABS)
    ! Reads the information appropriate for one band from a
    ! correlated k-distribution file.   The wavenumber range (cm^-1)
    ! in the file must match the desired in WAVENO.  The solar flux
    ! is returned in SOLFLUX, the number of "g"'s in NG, the weights
    ! or delta g's in DELG, the number of levels in NZCKD, the
    ! levels from the top down in ZCKD, and the absorption coefficients
    ! for each g and level in KABS.  The array sizes MAXNG and MAXNZ
    ! are checked.
    implicit none
    integer MAXNG, MAXNZ, NG, NZCKD
    real WAVENO(2), SOLFLUX, DELG(*), ZCKD(*), KABS(*)
    character CKDFILE*256
    integer NB, IB, JB, KB, IG, I, J, N
    real WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=CKDFILE, STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, SF, N
        if (abs(WAVENUM1 - WAVENO(1)) .lt. 1.0 .and.                     &
     &       abs(WAVENUM2 - WAVENO(2)) .lt. 1.0) then
            backspace (1)
            NG = N
            if (NG .gt. MAXNG) stop 'READ_CKD: MAXNG exceeded'
            read (1, *) JB, WAVENO(1), WAVENO(2), SOLFLUX,                 &
       &               N, (DELG(IG), IG=1, NG)
        end if
    end do
    if (JB .eq. 0) then
        write (*, *) ' Wavenumber range not found in CKD file : ',       &
     &              WAVENO(1), WAVENO(2), CKDFILE
        stop
    end if
    ! Read the Z level heights
    read (1, *) NZCKD
    if (NZCKD .gt. MAXNZ) stop 'READ_CKD: MAXNZ exceeded'
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
                write (*, '(1X,A,A,I3,A,I3,A,I2)')                           &
         &        'Negative absorption in CKD file.  ',                     &
         &        'Band: ', IB, '   Level: ', I, '   k: ', IG
            end if
        end do
    end do
    close (1)
    return
end

subroutine RESTORE_STATE(INSAVEFILE, NSTOKES, NX, NY, NZ,        &
    &      ML, MM, NLM, INRADFLAG, NEWGRIDFLAG, XGRID, YGRID, ZGRID,   &
    &      NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,          &
    &      CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)
    ! Reads the binary save file containing a previous state of
    ! the SHDOM solution. This includes the cell tree structure,
    ! the spherical harmonic expansion of the source function,
    ! the low order terms of the radiance expansion, and the flux array.
    ! If the file name is 'NONE' then the INRADFLAG is set to false,
    ! otherwise it is true.  The input file may have any spherical
    ! harmonic truncation smaller than the current (ML,MM,NLM).
    ! The flag NEWGRIDFLAG is set to false if reading succeeds.
    implicit none
    integer NSTOKES, NX, NY, NZ, ML, MM, NLM, NPTS, NCELLS
    integer SHPTR(*), RSHPTR(*)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    logical INRADFLAG, NEWGRIDFLAG
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real GRIDPOS(3, *), FLUXES(2, *), SOURCE(*), RADIANCE(*)
    character INSAVEFILE*256
    integer NSTI, NXI, NYI, NZI, MLI, MMI, NLMI
    integer I, J, NS, NR

    INRADFLAG = .false.
    NEWGRIDFLAG = .true.
    if (INSAVEFILE(1:4) .eq. 'NONE') return

    open (UNIT=3, FILE=INSAVEFILE, STATUS='OLD',                      &
   &              FORM='UNFORMATTED', ERR=900)
    read (3) NSTI, NXI, NYI, NZI, MLI, MMI, NLMI, NPTS, NCELLS
    if (NXI .ne. NX .or. NYI .ne. NY .or. NZI .ne. NZ) then
        write (6, *) 'RESTORE_STATE: wrong grid size:', NXI, NYI, NZI
        close (3)
        return
    end if
    if (NSTI .ne. NSTOKES .or.                                        &
   &    MLI .gt. ML .or. MMI .gt. MM) then
        write (6, *) 'RESTORE_STATE: incompatible SH:', NSTI, MLI, MMI
        close (3)
        return
    end if
    read (3) (XGRID(I), I=1, NX + 1)
    read (3) (YGRID(I), I=1, NY + 1)
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
        RSHPTR(I + 1) = RSHPTR(I) + NSTOKES*4
    end do

    close (3)
    INRADFLAG = .true.
    NEWGRIDFLAG = .false.
    return

900 continue
    write (6, *) 'Input radiance vector file not found.'
    INRADFLAG = .false.
    NEWGRIDFLAG = .true.
    return
end

subroutine SAVE_STATE(OUTSAVEFILE, NSTOKES, NX, NY, NZ,          &
    &       ML, MM, NLM, WORK, XGRID, YGRID, ZGRID,                    &
    &       NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,         &
    &       CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)
    ! Write a binary save file containing the current state of
    ! the SHDOM solution. This includes the cell tree structure,
    ! the spherical harmonic expansion of the source function,
    ! the low order terms of the radiance expansion, and the flux array.
    ! If the file name is 'NONE' then no file is written.
    implicit none
    integer NSTOKES, NX, NY, NZ, ML, MM, NLM, NPTS, NCELLS
    integer SHPTR(NPTS + 1), RSHPTR(NPTS)
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), WORK(*)
    real GRIDPOS(3, NPTS), FLUXES(2, NPTS), SOURCE(*), RADIANCE(*)
    character OUTSAVEFILE*256
    integer I, J, K, IR, NS, NR

    if (OUTSAVEFILE(1:4) .eq. 'NONE') return

    open (UNIT=3, FILE=OUTSAVEFILE, STATUS='UNKNOWN',                 &
   &              FORM='UNFORMATTED')
    write (3) NSTOKES, NX, NY, NZ, ML, MM, NLM, NPTS, NCELLS

    write (3) (XGRID(I), I=1, NX + 1)
    write (3) (YGRID(I), I=1, NY + 1)
    write (3) (ZGRID(I), I=1, NZ)
    write (3) ((GRIDPTR(J, I), J=1, 8), I=1, NCELLS)
    write (3) ((NEIGHPTR(J, I), J=1, 6), I=1, NCELLS)
    write (3) ((TREEPTR(J, I), J=1, 2), I=1, NCELLS)
    write (3) (CELLFLAGS(I), I=1, NCELLS)

    write (3) ((GRIDPOS(J, I), J=1, 3), I=1, NPTS)
    write (3) ((FLUXES(J, I), J=1, 2), I=1, NPTS)

    write (3) (SHPTR(I), I=1, NPTS + 1)
    NS = SHPTR(NPTS + 1)
    write (3) (SOURCE(I), I=1, NSTOKES*NS)

    K = 1
    do I = 1, NPTS
        IR = RSHPTR(I)
        do J = 1, NSTOKES*4
            WORK(K) = RADIANCE(NSTOKES*(IR - 1) + J)
            K = K + 1
        end do
    end do
    NR = 4*NPTS
    write (3) (WORK(I), I=1, NSTOKES*NR)

    close (3)
    return
end

subroutine OUTPUT_CELL_SPLIT(OUTFILE, GRIDPTR, GRIDPOS, EXTINCT,  &
    &                              SHPTR, SOURCE, NCELLS)
    ! Outputs the final cell splitting criterion for each grid cell
    ! and direction (X,Y,Z).  Calls CELL_SPLIT_TEST.
    ! Mainly for debugging purposes.
    implicit none
    integer NCELLS
    integer SHPTR(*), GRIDPTR(8, *)
    real GRIDPOS(3, *), EXTINCT(*), SOURCE(*)
    character OUTFILE*256
    integer IP1, IP2, ICELL, IDIR
    real ADAPTCRIT(3), MAXADAPT, XC, YC, ZC, EXT

    open (UNIT=9, FILE=OUTFILE, STATUS='UNKNOWN')
    write (9, *) '!  ICELL   Xc     Yc     Zc  ',                      &
   &            'Extinct SplitX  SplitY  SplitZ'

    ! Loop over all the cells
    do ICELL = 1, NCELLS
        ! Get the adaptive grid criterion for the 3 directions
        call CELL_SPLIT_TEST(GRIDPTR, GRIDPOS, EXTINCT,                &
     &                        SHPTR, SOURCE, ICELL,                     &
     &                        ADAPTCRIT, MAXADAPT, IDIR)
        ! Compute the cell center location and mean extinction
        IP1 = GRIDPTR(1, ICELL)
        IP2 = GRIDPTR(8, ICELL)
        XC = 0.5*(GRIDPOS(1, IP1) + GRIDPOS(1, IP2))
        YC = 0.5*(GRIDPOS(2, IP1) + GRIDPOS(2, IP2))
        ZC = 0.5*(GRIDPOS(3, IP1) + GRIDPOS(3, IP2))
        EXT = (EXTINCT(GRIDPTR(1, ICELL)) + EXTINCT(GRIDPTR(2, ICELL))    &
     &       + EXTINCT(GRIDPTR(3, ICELL)) + EXTINCT(GRIDPTR(4, ICELL))    &
     &       + EXTINCT(GRIDPTR(5, ICELL)) + EXTINCT(GRIDPTR(6, ICELL))    &
     &       + EXTINCT(GRIDPTR(7, ICELL)) + EXTINCT(GRIDPTR(8, ICELL)))/8
        write (9, '(1X,I5,3(1X,F6.3),1X,E11.4,3(1X,F7.5))')              &
     &        ICELL, XC, YC, ZC, EXT,                                   &
     &        ADAPTCRIT(1), ADAPTCRIT(2), ADAPTCRIT(3)
    end do

    close (9)
    return
end

subroutine VISUALIZE_CELLS(OUTFILE, YVAL, OUTTYPE,               &
    &      IPFLAG, NX, XDOMAIN, ZDOMAIN,                               &
    &      NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS)
    ! Visualizes the cell structure by outputting GLE code to draw
    ! a diagram of the cells in the XZ plane.  Only the end cells are
    ! output.  The outlines of the cells are drawn, and arrows between
    ! the cell centers show the neighbor relationships.
    ! OUTTYPE=0 for cell outlines only, OUTTYPE=1 for cells and arrows.
    ! There is a special output for IP mode.
    implicit none
    integer OUTTYPE, IPFLAG, NX, NCELLS
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real YVAL, XDOMAIN, ZDOMAIN, GRIDPOS(3, *)
    character OUTFILE*256
    integer IC
    integer J, IN
    logical BTEST
    real XSIZE, ZSIZE, S, DELX
    real X1, Z1, X2, Z2, X3, Z3, X4, Z4, XC, ZC
    real DX, DZ, XN, ZN, XA, ZA

    S = min(24.0/XDOMAIN, 18.0/ZDOMAIN)
    XSIZE = S*XDOMAIN
    ZSIZE = S*ZDOMAIN
    open (UNIT=9, FILE=OUTFILE, STATUS='UNKNOWN')
    write (9, '(A,2(1X,F6.2))') 'size ', XSIZE, ZSIZE
    write (9, '(A)') 'set just cc'
    write (9, '(A)') 'amove 0 0'
    write (9, '(A)') 'begin origin'
    write (9, '(A)') 'begin scale 1 1'
    write (9, '(A,1X,F7.3)') '! Y value: ', YVAL
    write (9, '(A)')
    if (btest(IPFLAG, 0)) then
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
        if (TREEPTR(2, IC) .eq. 0 .and.                                  &
     &      GRIDPOS(2, GRIDPTR(1, IC)) .le. YVAL .and.                    &
     &      GRIDPOS(2, GRIDPTR(3, IC)) .ge. YVAL) then
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
            if (btest(IPFLAG, 0)) then
                write (9, '(A,6(1X,F6.3),I8)') '@drawbox ',                  &
         &       S*X1, S*Z1, S*X3, S*Z3, S*XC, S*ZC, IC
            else
                write (9, '(A,10(1X,F6.3),I8)') '@drawbox ',                 &
         &       S*X1, S*Z1, S*X2, S*Z2, S*X3, S*Z3, S*X4, S*Z4, S*XC, S*ZC, IC
            end if
            if (OUTTYPE .eq. 1) then
                ! Draw the neighbor arrows
                DX = 0.1*(X2 - X1)
                DZ = 0.1*(Z4 - Z1)
                do J = 1, 6
                    if (NEIGHPTR(J, IC) .ne. 0 .and. J .ne. 3 .and. J .ne. 4) then
                        IN = abs(NEIGHPTR(J, IC))
                        XN = (GRIDPOS(1, GRIDPTR(1, IN)) + GRIDPOS(1, GRIDPTR(8, IN)))/2
                        ZN = (GRIDPOS(3, GRIDPTR(1, IN)) + GRIDPOS(3, GRIDPTR(8, IN)))/2
                        XA = XC + 1.0*(XN - XC)
                        ZA = ZC + 1.0*(ZN - ZC)
                        if (.not. (XN - XC .gt. 0.0 .and. J .eq. 1 .or.             &
               &                   XN - XC .lt. 0.0 .and. J .eq. 2)) then
                            write (9, '(A,4(1X,F6.3),1X,I8)') '@drawarrow ',         &
                 &            S*XC, S*ZC, S*XA, S*ZA, NEIGHPTR(J, IC)
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
end

subroutine OUTPUT_RESULTS(NSTOKES, NX, NY, NZ, NBPTS, NPTS, NCELLS, &
    &             NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,              &
    &             NMU, NPHI, NANG, NG,                                 &
    &             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,  &
    &             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,           &
    &             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,                 &
    &             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
    &             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,           &
    &             XGRID, YGRID, ZGRID, GRIDPOS,                        &
    &             TREEPTR, GRIDPTR, CELLFLAGS,                         &
    &             EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,                &
    &             FLUXES, DIRFLUX, FLUXDIV,                            &
    &             IRAD, RADOUT, NSHOUT, SHTERMS, SOURCE1OUT,           &
    &             OUTTYPE, OUTPARMS, OUTFILE)
    ! Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are five types (OUTTYPE) of output: 'R' - for radiance,
    ! 'F' - for hemispheric flux, 'H' - for heating rate (net flux
    ! convergence), 'S' - for spherical harmonic terms, and 'M' for
    ! medium properties.
    ! See the overall program documentation for output parameters.
    implicit none
    integer NSTOKES, NX, NY, NZ, NBPTS, NPTS, NCELLS, NSH
    integer ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, NMU, NPHI, NANG, NG
    integer BCFLAG, IPFLAG
    integer MAXITER, TOTITER, IRAD, NSHOUT
    integer TREEPTR(2, NCELLS), GRIDPTR(8, NCELLS)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real SOLACC, SPLITACC, SHACC
    real WAVENO(2), WAVELEN
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real OUTPARMS(*)
    character SRCTYPE*1, SFCTYPE*2, UNITS*1, GRIDTYPE*1, OUTTYPE*1
    character*256 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 OUTFILE
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real TEMP(NPTS)
    real FLUXES(2, NPTS), DIRFLUX(NPTS), FLUXDIV(NPTS)
    real RADOUT(NSTOKES, *), SHTERMS(NSHOUT, NPTS)
    real SOURCE1OUT(NSTOKES, *)

    integer IZ, I, K, N, NX2, NY2, J, JX, JY, ICELL
    integer NANGOUT, NXOUT, NYOUT
    logical ALLGRID
    real MUOUT, PHID, C, PI, SUM, SUM1, SUM2, SUM3
    real XDOMAIN, YDOMAIN, STARTX, STARTY
    real FUP, FDN, FDIR, ASYM
    double precision X0, Y0, Z0
    character*128 FORM
    character*128 GRIDNAME, SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

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

    open (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
    write (2, '(A,A)')                                                 &
   &  '! Polarized Spherical Harmonic Discrete Ordinate',             &
   &  ' Radiative Transfer Output'
    write (2, '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)')                     &
   &    '!  L=', ML, '  M=', MM, '  NLM=', NLM,                           &
   &    '   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,               &
   &    '   NSH=', NSH
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', NSTOKES,         &
   &  '   NX=', NX, '   NY=', NY, '   NZ=', NZ,                            &
   &  '    NPTS=', NPTS, '   NCELLS=', NCELLS
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(PROPFILE)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',              &
   &                    trim(CKDFILE), '   NUM_G=', NG
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)
    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,           &
     &                        '      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if
    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,               &
   &                         '   INDEPENDENT_PIXEL=', IPFLAG
    write (2, '(A,A22,A,I1)') '!  SURFACE_TYPE=', SFCNAME,             &
   &                         '   HORIZ_BOUNDARY_COND=', BCFLAG
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(SFCFILE),   &
       &                                 '  SKY_RAD=', SKYRAD
        else
            if (SRCTYPE .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')                         &
         &        '!  GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,  &
         &        '  SKY_RAD=', SKYRAD
            else
                write (2, '(A,F9.7,A,E12.5)')                                &
         &        '!  GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD
            end if
        end if
        write (2, '(A,E13.6,A,F10.7,A,F8.3)')                            &
     &         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,    &
     &         '   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0)
    else
        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(SFCFILE)
        else
            write (2, '(A,F8.3,A,F9.7)')                                   &
       &      '!  GROUND_TEMP=', GNDTEMP, '  GROUND_EMIS=', 1.0 - GNDALBEDO
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', SKYRAD
    end if

    if (UNITS .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (UNITS .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,       &
     &           '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,             &
     &           '   WAVELENGTH=', WAVELEN
    end if
    write (2, '(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,      &
   &                         '   SPHERICAL_HARMONIC_ACCURACY=', SHACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER,       &
   &                         '   NUMBER_ITERATIONS=', TOTITER
    write (2, '(A,A20)') '!  OUTPUT_TYPE=', OUTNAME

    PI = acos(-1.0)

    if (OUTTYPE .eq. 'R') then
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

        write (2, '(A,F7.3,3(A,I4))') '!    RADIANCE AT Z=', Z0,         &
     &      '    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (UNITS .eq. 'T') then
            FORM = '(2(1X,F7.3),1X,4(1X,F6.2))'
            write (2, '(A)')                                               &
       &      '!   X       Y       RADIANCE (STOKES I, Q, U, V)'
        else
            FORM = '(2(1X,F7.3),1X,4(1X,E12.5))'
            write (2, '(A)')                                               &
       &      '!   X       Y        RADIANCE (STOKES I, Q, U, V)'
        end if

        do K = 1, NANGOUT
            MUOUT = OUTPARMS(2*K + 5)
            PHID = OUTPARMS(2*K + 6)
            write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')                          &
       &          '! ', MUOUT, PHID, '<- (mu,phi)'
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
                    write (2, '(A)')                                           &
           &     '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(2(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)') '!   X       Y           UP           DOWN'
                end if
            end if
            C = 1.0
            if (UNITS .eq. 'T') C = 1.0/PI
            do I = 1, NBPTS, NZ
                if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                     &
         &          GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
                    if (SRCTYPE .ne. 'T') then
                        write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I),              &
             &            C*FLUXES(2, I + NZ - 1), C*FLUXES(1, I), DIRFLUX(I)
                    else
                        write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I),              &
             &              C*FLUXES(2, I + NZ - 1), C*FLUXES(1, I)
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
            write (2, '(A,F7.3,2(A,I4))') '!    UPWELLING FLUX:     Z=', Z0, &
       &        '    NXO=', NXOUT, '    NYO=', NYOUT
            write (2, '(A,F7.3)') '!    DOWNWELLING FLUX:   Z=', Z0
            if (UNITS .eq. 'T') then
                FORM = '(2(1X,F7.3),1X,2(2X,F6.2))'
                write (2, '(A)') '!   X       Y        UP     DOWN'
            else
                if (SRCTYPE .ne. 'T') then
                    FORM = '(2(1X,F7.3),1X,3(2X,E12.5))'
                    write (2, '(A)')                                           &
           &     '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
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
                    call LOCATE_GRID_CELL(NX, NY, NZ, XGRID, YGRID, ZGRID,   &
           &                  NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,   &
           &                  BCFLAG, IPFLAG, X0, Y0, Z0, ICELL)
                    call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,          &
           &                             X0, Y0, Z0, 2, FLUXES, 1, 1, FDN)
                    call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,          &
           &                             X0, Y0, Z0, 2, FLUXES, 2, 1, FUP)
                    call INTERPOLATE_FIELD(ICELL, GRIDPTR, GRIDPOS,          &
           &                             X0, Y0, Z0, 1, DIRFLUX, 1, 1, FDIR)
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
                    write (2, '(A)')                                           &
           &            '!    Z          UP       DOWN_DIFFUSE   DOWN_DIRECT'
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
                    if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                   &
           &            GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
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
                    write (2, '(A,A)') '!   X       Y       Z',                &
           &          '           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(3(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)')                                           &
           &       '!   X       Y       Z           UP           DOWN'
                end if
            end if
            C = 1.0
            if (UNITS .eq. 'T') C = 1.0/PI
            do I = 1, N
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                    &
         &           GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    if (SRCTYPE .ne. 'T') then
                        write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),  &
             &            C*FLUXES(2, I), C*FLUXES(1, I), DIRFLUX(I)
                    else
                        write (2, FORM) GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),  &
             &            C*FLUXES(2, I), C*FLUXES(1, I)
                    end if
                end if
            end do
        end if

    else if (OUTTYPE .eq. 'H') then

        ! Heating output: net flux convergence
        if (nint(OUTPARMS(1)) .eq. 1) then
            ! Format 1 is domain averaged vertical profile
            write (2, '(A,A)') '!    Z       -DIV(Fnet)'
            do IZ = 1, NZ
                SUM = 0.0
                N = 0
                do I = IZ, NBPTS, NZ
                    if (GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                   &
           &            GRIDPOS(2, I) .lt. YGRID(NY + 1)) then
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
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                    &
         &           GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    write (2, '(3(1X,F7.3),2X,E13.5)')                         &
           &          GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I), -FLUXDIV(I)
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
            write (2, '(A,A)') '!    X       Y        Z         Imean',    &
       &     '          Fx           Fy           Fz       HOrms'
            do I = 1, N
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                    &
         &         GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')           &
           &          GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),               &
           &          SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I),   &
           &          SHTERMS(5, I)/max(1e-20, SHTERMS(1, I))
                end if
            end do
        else
            write (2, '(A,A)') '!    X       Y        Z         Imean',    &
       &     '          Fx           Fy           Fz '
            do I = 1, N
                if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                    &
         &         GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                    write (2, '(3(1X,F7.3),1X,4(1X,E12.5))')                   &
           &          GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),               &
           &          SHTERMS(1, I), SHTERMS(2, I), SHTERMS(3, I), SHTERMS(4, I)
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
        write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')                            &
     &          '! ', MUOUT, PHID, '<- (mu,phi)'
        write (2, '(A,A)') '!    X       Y       Z    ',                 &
     &      ' Extinction   Source Function (Stokes I, Q, U, V)'
        do I = 1, N
            if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                      &
       &         GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                write (2, '(3(1X,F7.3),5(1X,E12.5))')                        &
         &          GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),               &
         &          EXTINCT(I), SOURCE1OUT(:, I)
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
        write (2, '(A,A)') '!    X       Y       Z',                     &
       &         '      Extinct  Albedo  Asymmetry Temp'
        do I = 1, N
            if (NUMPHASE .gt. 0) then
                ASYM = LEGEN(1, 1, IPHASE(I))
            else
                ASYM = LEGEN(1, 1, I)
            end if
            if ((GRIDPOS(1, I) .lt. XGRID(NX + 1) .and.                      &
           &     GRIDPOS(2, I) .lt. YGRID(NY + 1)) .or. ALLGRID) then
                write (2, '(3(1X,F7.3),1X,F10.6,1X,F8.6,1X,F8.5,1X,F7.2)')   &
               &       GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I),              &
               &       EXTINCT(I), ALBEDO(I), ASYM, TEMP(I)
            end if
        end do

    end if

    close (2)
    return
end

subroutine OUTPUT_RESULTS_PAR(NSTOKES, NX, NY, NZ, NPTS, NCELLS, &
&             NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,                   &
&             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,  &
&             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,           &
&             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,                 &
&             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
&             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,           &
&             XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID,               &
&             FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT,      &
&             OUTTYPE, OUTPARMS, OUTFILE)
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

    open (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
    write (2, '(A,A)')                                                 &
   &  '! Polarized Spherical Harmonic Discrete Ordinate',             &
   &  ' Radiative Transfer Output'
    write (2, '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10)')                     &
   &    '!  L=', ML, '  M=', MM, '  NLM=', NLM,                           &
   &    '   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,               &
   &    '   NSH=', NSH
    write (2, '(A,I1,3(A,I4),2(A,I8))') '!  NSTOKES=', NSTOKES,         &
   &  '   NX=', NX, '   NY=', NY, '   NZ=', NZ,                            &
   &  '   NPTS=', NPTS, '   NCELLS=', NCELLS
    write (2, '(A,A)') '!  PROPERTY_FILE=', trim(PROPFILE)
    write (2, '(A,A,A,I2)') '!  CORRELATED_K-DIST_FILE=',              &
   &                     trim(CKDFILE), '   NUM_G=', NG
    write (2, '(A,A)') '!  INPUT_SAVE_FILE=', trim(INSAVEFILE)
    write (2, '(A,A)') '!  OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE)
    if (DELTAM) then
        write (2, '(A,A14,A)') '!  SOURCE_TYPE=', SOURCENAME,           &
     &                        '      DELTA-M METHOD'
    else
        write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    end if
    write (2, '(A,A32,A,I1)') '!  GRID_TYPE=', GRIDNAME,               &
   &                         '   INDEPENDENT_PIXEL=', IPFLAG
    write (2, '(A,A22,A,I2)') '!  SURFACE_TYPE=', SFCNAME,             &
   &                         '   HORIZ_BOUNDARY_COND=', BCFLAG
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A,A,E12.5)') '!  SURFACE_FILE=', trim(SFCFILE),   &
       &                                 '  SKY_RAD=', SKYRAD
        else
            if (SRCTYPE .eq. 'B') then
                write (2, '(A,F9.7,A,F8.3,A,E12.5)')                         &
         &        '!  GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,  &
         &        '  SKY_RAD=', SKYRAD
            else
                write (2, '(A,F9.7,A,E12.5)')                                &
         &        '!  GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD
            end if
        end if
        write (2, '(A,E13.6,A,F10.7,A,F8.3)')                            &
     &         '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,    &
     &         '   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0)
    else
        if (SFCTYPE(1:1) .eq. 'V') then
            write (2, '(A,A)') '!  SURFACE_FILE=', trim(SFCFILE)
        else
            write (2, '(A,F8.3,A,F9.7)')                                   &
       &      '!  GROUND_TEMP=', GNDTEMP, '  GROUND_EMIS=', 1.0 - GNDALBEDO
        end if
        write (2, '(A,F8.3)') '!  SKY_TEMP=', SKYRAD
    end if

    if (UNITS .eq. 'T') then
        write (2, '(A,A24)') '!  UNITS=', UNITSNAME
    else if (UNITS .eq. 'B') then
        write (2, '(A,A24,A,2(1X,F10.2))') '!  UNITS=', UNITSNAME,       &
     &           '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2)
    else
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME,             &
     &           '   WAVELENGTH=', WAVELEN
    end if
    write (2, '(2(A,E10.3))') '!  SPLITTING_ACCURACY=', SPLITACC,      &
   &                         '   SPHERICAL_HARMONIC_ACCURACY=', SHACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER,       &
   &                         '   NUMBER_ITERATIONS=', TOTITER
    write (2, '(A,A20)') '!  OUTPUT_TYPE=', OUTNAME

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

        write (2, '(A,F7.3,3(A,I4))') '!    RADIANCE AT Z=', Z0,         &
     &      '    NXO=', NXOUT, '    NYO=', NYOUT, '    NDIR=', NANGOUT
        if (UNITS .eq. 'T') then
            FORM = '(2(1X,F7.3),1X,4(1X,F6.2))'
            write (2, '(A)')                                               &
       &       '!   X       Y      RADIANCE (STOKES I, Q, U, V)'
        else
            FORM = '(2(1X,F7.3),1X,4(1X,E12.5))'
            write (2, '(A)')                                               &
       &       '!   X       Y        RADIANCE (STOKES I, Q, U, V)'
        end if

        do K = 1, NANGOUT
            MUOUT = OUTPARMS(2*K + 5)
            PHID = OUTPARMS(2*K + 6)
            write (2, '(A,1X,F8.5,1X,F7.2,2X,A)')                          &
       &          '! ', MUOUT, PHID, '<- (mu,phi)'
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
                    write (2, '(A)')                                           &
           &     '!   X       Y           UP       DOWN_DIFFUSE   DOWN_DIRECT'
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
                        write (2, FORM) XGRID(I), YGRID(J),                      &
             &            C*FLUXES(2, NZ, J, I), C*FLUXES(1, 1, J, I), FLUXES(3, 1, J, I)
                    else
                        write (2, FORM) XGRID(I), YGRID(J),                      &
             &              C*FLUXES(2, NZ, J, I), C*FLUXES(1, 1, J, I)
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
                    write (2, '(A)')                                           &
           &            '!    Z          UP       DOWN_DIFFUSE   DOWN_DIRECT'
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
                    write (2, '(A,A)') '!   X       Y       Z',                &
           &          '           UP       DOWN_DIFFUSE   DOWN_DIRECT'
                else
                    FORM = '(3(1X,F7.3),1X,2(2X,E12.5))'
                    write (2, '(A)')                                           &
           &       '!   X       Y       Z           UP           DOWN'
                end if
            end if
            C = 1.0
            if (UNITS .eq. 'T') C = 1.0/PI
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        if (SRCTYPE .ne. 'T') then
                            write (2, FORM) XGRID(I), YGRID(J), ZGRID(K),              &
                 &            C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I), FLUXES(3, K, J, I)
                        else
                            write (2, FORM) XGRID(I), YGRID(J), ZGRID(K),              &
                 &            C*FLUXES(2, K, J, I), C*FLUXES(1, K, J, I)
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
                        write (2, '(3(1X,F7.3),2X,E13.5)')                         &
               &             XGRID(I), YGRID(J), ZGRID(K), -FLUXDIV(K, J, I)
                    end do
                end do
            end do
        end if

    else if (OUTTYPE .eq. 'S') then

        ! Spherical Harmonic output:
        !   Output mean intensity and net flux (x,y,z components).
        !   and maybe normalized rms of higher order terms
        if (NSHOUT .eq. 5) then
            write (2, '(A,A)') '!    X       Y        Z         Imean',    &
       &     '          Fx           Fy           Fz       HOrms'
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')           &
               &          XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1, K, J, I),         &
               &          SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I),     &
               &          SHTERMS(5, K, J, I)/max(1e-20, SHTERMS(1, K, J, I))
                    end do
                end do
            end do
        else
            write (2, '(A,A)') '!    X       Y        Z         Imean',    &
       &     '          Fx           Fy           Fz '
            do I = 1, NX
                do J = 1, NY
                    do K = 1, NZ
                        write (2, '(3(1X,F7.3),1X,4(1X,E12.5),1X,F6.3)')           &
               &          XGRID(I), YGRID(J), ZGRID(K), SHTERMS(1, K, J, I),         &
               &          SHTERMS(2, K, J, I), SHTERMS(3, K, J, I), SHTERMS(4, K, J, I)
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
end

subroutine OUTPUT_IMAGE(NSTOKES, NX, NY, NZ, NPTS, NCELLS,       &
    &             NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,                 &
    &             PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,  &
    &             BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,           &
    &             SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,                 &
    &             SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
    &             SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,           &
    &             IVIS, VISOUT, OUTPARMS, OUTFILE)
    ! Writes the visualization image in PGM binary format with either
    ! one byte or two byte pixels.  The PGM binary format (P5) has an
    ! ascii header which has the image parameters amd can include commen
    ! (here the SHDOM input parameters), followed by the byte or two byt
    ! image pixels.
    implicit none
    integer NSTOKES, NX, NY, NZ, NPTS, NCELLS, NSH
    integer ML, MM, NLM, NMU, NPHI, NANG, NG
    integer BCFLAG, IPFLAG
    integer MAXITER, TOTITER, IVIS
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD
    real SOLACC, SPLITACC, SHACC
    real WAVENO(2), WAVELEN
    real OUTPARMS(*)
    character SRCTYPE*1, SFCTYPE*2, UNITS*1, GRIDTYPE*1
    character*256 PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 OUTFILE
    real VISOUT(NSTOKES, *)

    integer MAXHEAD, NL, NS, NLOUT, LABRECS, RECLEN, NBYTES
    integer DN, I, J, K, L, M
    real MAXLINPOL, MINANGPOL, MAXANGPOL, MAXCIRCPOL
    real VISMIN(4), VISMAX(4), VISSCALE(4)
    logical CAMERA_MODE, LSBFIRST
    integer*1, allocatable :: BUFFER1(:)
    character*128 GRIDNAME, SOURCENAME, UNITSNAME, SFCNAME
    character*128 POLNAMES(4), VISFORM(4)
    character*1 LF
    character BUF*2000

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
        NL = 1 + sqrt((OUTPARMS(4) - OUTPARMS(7))**2                      &
     &               + (OUTPARMS(5) - OUTPARMS(8))**2                      &
     &               + (OUTPARMS(6) - OUTPARMS(9))**2)/OUTPARMS(10)
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
        if (NSTOKES .gt. 1) then
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
        if (NSTOKES .eq. 4) then
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
    open (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN',                     &
   &      ACCESS='DIRECT', FORM='UNFORMATTED', RECL=RECLEN)

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

    write (BUF(I:I + 79), '(A,A1)')                                      &
   &  '# Polarized SHDOM Visualization Output', LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(2(A,I3),A,I5,2(A,I3),A,I5,A,I10,A1)')        &
   &    '# L=', ML, '  M=', MM, '  NLM=', NLM,                            &
   &    '   NMU=', NMU, '  NPHI=', NPHI, '  NANG=', NANG,               &
   &    '   NSH=', NSH, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 90), '(A,I1,3(A,I4),2(A,I8),A1)')                   &
   &  '# NSTOKES=', NSTOKES, '  NX=', NX, '   NY=', NY, '  NZ=', NZ,        &
   &  '    NPTS=', NPTS, '   NCELLS=', NCELLS, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A,A1)')                                    &
   &  '# PROPERTY_FILE=', trim(PROPFILE), LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A,A,I2,A1)')                               &
   &  '# CORRELATED_K-DIST_FILE=', trim(CKDFILE), '   NUM_G=', NG, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A,A1)')                                    &
   &  '# INPUT_SAVE_FILE=', trim(INSAVEFILE), LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A,A1)')                                    &
   &  '# OUTPUT_SAVE_FILE=', trim(OUTSAVEFILE), LF
    I = I + index(BUF(I:I + 79), LF)
    if (DELTAM) then
        write (BUF(I:I + 79), '(A,A14,A,A1)')                              &
     &    '# SOURCE_TYPE=', SOURCENAME, '    DELTA-M METHOD', LF
        I = I + index(BUF(I:I + 79), LF)
    else
        write (BUF(I:I + 79), '(A,A14,A1)')                                &
     &    '# SOURCE_TYPE=', SOURCENAME, LF
        I = I + index(BUF(I:I + 79), LF)
    end if
    write (BUF(I:I + 79), '(A,A32,A,I1,A1)') '# GRID_TYPE=', GRIDNAME,   &
   &                         '   INDEPENDENT_PIXEL=', IPFLAG, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,A22,A,I2,A1)')                             &
   & '# SURFACE_TYPE=', SFCNAME, '   HORIZ_BOUNDARY_COND=', BCFLAG, LF
    I = I + index(BUF(I:I + 79), LF)
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        if (SFCTYPE(1:1) .eq. 'V') then
            write (BUF(I:I + 79), '(A,A,A,E12.5,A1)')                        &
       &        '# SURFACE_FILE=', trim(SFCFILE),                         &
       &        '  SKY_RAD=', SKYRAD, LF
            I = I + index(BUF(I:I + 79), LF)
        else
            if (SRCTYPE .eq. 'B') then
                write (BUF(I:I + 79), '(A,F9.7,A,F8.3,A,E12.5,A1)')            &
         &        '# GROUND_ALBEDO=', GNDALBEDO, '  GROUND_TEMP=', GNDTEMP,   &
         &        '  SKY_RAD=', SKYRAD, LF
                I = I + index(BUF(I:I + 79), LF)
            else
                write (BUF(I:I + 79), '(A,F9.7,A,E12.5,A1)')                   &
         &        '# GROUND_ALBEDO=', GNDALBEDO, '  SKY_RAD=', SKYRAD, LF
                I = I + index(BUF(I:I + 79), LF)
            end if
        end if
        write (BUF(I:I + 79), '(A,E13.6,A,F10.7,A,F8.3,A1)')               &
     &         '# SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU,     &
     &         '   SOLAR_AZ=', SOLARAZ*180.0/acos(-1.0), LF
        I = I + index(BUF(I:I + 79), LF)
    else
        if (SFCTYPE(1:1) .eq. 'V') then
            write (BUF(I:I + 79), '(A,A,A1)')                                &
       &      '# SURFACE_FILE=', trim(SFCFILE), LF
            I = I + index(BUF(I:I + 79), LF)
        else
            write (BUF(I:I + 79), '(A,F8.3,A,F9.7,A1)')                      &
       &      '# GROUND_TEMP=', GNDTEMP,                                   &
       &      '  GROUND_EMIS=', 1.0 - GNDALBEDO, LF
            I = I + index(BUF(I:I + 79), LF)
        end if
        write (BUF(I:I + 79), '(A,F8.3,A1)') '# SKY_TEMP=', SKYRAD, LF
        I = I + index(BUF(I:I + 79), LF)
    end if

    if (UNITS .eq. 'T') then
        write (BUF(I:I + 79), '(A,A24,A1)') '# UNITS=', UNITSNAME, LF
        I = I + index(BUF(I:I + 79), LF)
    else if (UNITS .eq. 'B') then
        write (BUF(I:I + 79), '(A,A24,A,2(1X,F10.2),A1)')                  &
     &    '# UNITS=', UNITSNAME,                                        &
     &    '   WAVENUMBER_RANGE=', WAVENO(1), WAVENO(2), LF
        I = I + index(BUF(I:I + 79), LF)
    else
        write (BUF(I:I + 79), '(A,A24,A,F10.2,A1)')                        &
     &    '# UNITS=', UNITSNAME, '   WAVELENGTH=', WAVELEN, LF
        I = I + index(BUF(I:I + 79), LF)
    end if
    write (BUF(I:I + 79), '(2(A,E10.3),A1)')                             &
   &  '# SPLITTING_ACCURACY=', SPLITACC,                              &
   &  '   SPHERICAL_HARMONIC_ACCURACY=', SHACC, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(A,E10.3,A1)')                                &
   &  '# SOLUTION_ACCURACY=', SOLACC, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 79), '(2(A,I4),A1)')                                &
   &  '# MAXIMUM_ITERATIONS=', MAXITER,                               &
   &  '   NUMBER_ITERATIONS=', TOTITER, LF
    I = I + index(BUF(I:I + 79), LF)
    write (BUF(I:I + 7), '(A,A1)') '# ', LF
    I = I + index(BUF(I:I + 7), LF)

    if (CAMERA_MODE) then
        write (BUF(I:I + 79), '(A,A1)') '# CAMERA MODE:', LF
        I = I + index(BUF(I:I + 79), LF)
!          1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline,
        write (BUF(I:I + 79), '(A,F9.2,A1)')                               &
     &      '# RADIANCE_SCALE=', OUTPARMS(3), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,3(1X,F9.3),A1)') '# CAMERA_LOCATION=',   &
     &      OUTPARMS(4), OUTPARMS(5), OUTPARMS(6), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F10.5,A1)')                              &
     &      '# CAMERA_ZENITH_ANGLE=', OUTPARMS(7), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F10.5,A1)')                              &
     &      '# CAMERA_AZIMUTH_ANGLE=', OUTPARMS(8), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F10.5,A1)')                              &
     &      '# CAMERA_ROTATION_ANGLE=', OUTPARMS(9), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,I5,A1)')                                 &
     &      '# CAMERA_NUMBER_LINES=', nint(OUTPARMS(10)), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,I5,A1)')                                 &
     &      '# CAMERA_NUMBER_SAMPLES=', nint(OUTPARMS(11)), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.5,A1)')                               &
     &      '# CAMERA_LINE_SPACING_ANGLE=', OUTPARMS(12), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.5,A1)')                               &
     &      '# CAMERA_SAMPLE_SPACING_ANGLE=', OUTPARMS(13), LF
        I = I + index(BUF(I:I + 79), LF)

    else
        write (BUF(I:I + 79), '(A,A1)')                                    &
     &      '# CROSS TRACK SCANNING MODE:', LF
        I = I + index(BUF(I:I + 79), LF)
!          2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, d
        write (BUF(I:I + 79), '(A,F9.2,A1)')                               &
     &      '# RADIANCE_SCALE=', OUTPARMS(3), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,3(1X,F9.3),A1)') '# START_LOCATION=',    &
     &      OUTPARMS(4), OUTPARMS(5), OUTPARMS(6), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,3(1X,F9.3),A1)') '# END_LOCATION=',      &
     &      OUTPARMS(7), OUTPARMS(8), OUTPARMS(9), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.4,A1)')                               &
     &      '# SPACING_ON_TRACK=', OUTPARMS(10), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.3,A1)')                               &
     &      '# START_SCAN_ANGLE=', OUTPARMS(11), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F8.3,A1)')                               &
     &      '# END_SCAN_ANGLE=', OUTPARMS(12), LF
        I = I + index(BUF(I:I + 79), LF)
        write (BUF(I:I + 79), '(A,F7.4,A1)')                               &
     &      '# SCAN_ANGLE_SPACING=', OUTPARMS(13), LF
        I = I + index(BUF(I:I + 79), LF)
    end if
    do M = 1, NSTOKES
        write (BUF(I:I + 79), VISFORM(M)) '# ', trim(POLNAMES(M)),          &
     &        '_RANGE_IN_PIXEL=', VISMIN(M), VISMAX(M), LF
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
!              IF (DN .GT. 127) THEN
!                BUFFER1(I) = DN - 256
!              ELSE
!                BUFFER1(I) = DN
!              ENDIF
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
end
