!   SHDOM: Spherical harmonic discrete ordinate radiative transfer method.
!      See shdom.txt for documentation.

module shd_common_arr
    integer, save :: NPX, NPY, NPZ, NUMPHASE
    real, save    :: DELX, DELY, XSTART, YSTART
    real, save, allocatable :: ZLEVELS(:)
    real, save, allocatable :: TEMPP(:), EXTINCTP(:), ALBEDOP(:)
    real, save, allocatable :: LEGENP(:), EXTDIRP(:)
    integer*2, save, allocatable :: IPHASEP(:)
    integer, save :: NZCKD
    real, save, allocatable :: ZCKD(:), GASABS(:)
end module shd_common_arr

program shdom_main
    use shd_common_arr
    use shd_io_inp
    use shd_io_out
    implicit none
    ! Number of output files and inputs parameter per file
    integer, parameter :: MAXOUT = 20, MAXPAR = 500
    real      :: OUTPARMS(MAXPAR, MAXOUT)
    character :: OUTTYPES(MAXOUT)*1, OUTFILES(MAXOUT)*256

    real, allocatable :: DELG(:), KABS(:, :)
    real, allocatable :: XGRID(:), YGRID(:), ZGRID(:)
    real, allocatable :: TEMP(:), PLANCK(:), EXTINCT(:), ALBEDO(:)
    real, allocatable :: LEGEN(:)
    integer*2, allocatable:: IPHASE(:)
    integer, allocatable :: GRIDPTR(:, :), NEIGHPTR(:, :), TREEPTR(:, :)
    integer*2, allocatable:: CELLFLAGS(:)
    real, allocatable :: GRIDPOS(:, :)
    integer, allocatable :: RSHPTR(:), SHPTR(:), OSHPTR(:)
    real, allocatable :: RADIANCE(:, :), SOURCE(:, :), DELSOURCE(:, :)
    integer, allocatable :: WORK1(:)
    real, allocatable :: WORK(:), WORK2(:)
    real, allocatable :: FLUXES(:), DIRFLUX(:)
    integer, allocatable :: NPHI0(:)
    real, allocatable :: MU(:), WTDO(:), PHI(:)
    integer, allocatable :: BCPTR(:)
    real, allocatable :: BCRAD(:, :)
    real, allocatable :: SFCPARMS(:)
    real, allocatable :: SFCGRIDPARMS(:)
    real, allocatable :: SUMFLUXES(:), SUMDIRFLUX(:), SUMFLUXDIV(:)
    real, allocatable :: SUMSHTERMS(:), SOURCE1OUT(:)
    real, allocatable :: SUMRADOUT(:), SUMVISOUT(:)
    real, allocatable :: TEMPPT(:), EXTINCTPT(:), ALBEDOPT(:)
    real, allocatable :: LEGENPT(:)
    integer*2, allocatable :: IPHASEPT(:)
    real, allocatable :: ALLFLUXES(:, :, :, :)
    real, allocatable :: ALLFLUXDIV(:, :, :), ALLSHTERMS(:, :, :, :)
    real, allocatable :: XALLGRID(:), YALLGRID(:)

    integer :: MAXPG, MAXIG, MAXIC, MAXIV, MAXIDO, MAXNZ, MAXNG
    integer :: MAXLEG, MAXPGL, MAXIGL
    integer :: MAXNBC, MAXBCRAD
    integer :: MAXSFCPTS, MAXSFCPARS, BIG_ARRAYS
    real :: MEMWORD, WANTMEM
    real :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer :: NSTOKES, NX, NY, NZ, NX1, NY1
    integer :: NXT, NYT, NPXT, NPYT, MAXPGT, MAXPGLT
    integer :: NXSFC, NYSFC, NSFCPAR, NTOPPTS, NBOTPTS
    integer :: NG, SIDE
    integer :: ML, MM, NLM, NSTLEG, NLEG, NANG, NMU, NPHI, NPHI0MAX
    integer :: NBPTS, NPTS, NCELLS, NSH, NBCELLS, OLDNPTS
    integer :: MAXITER, ITER, TOTITER, IG, I, NUMOUT, BCFLAG, IPFLAG
    integer :: NRAD, IRAD, NXOUT, NYOUT, NVIS, IVIS, NSCAN, NPIX, NSHOUT
    integer :: NCELLSTOT, NPTSTOT, NSHTOT
    logical :: KDIST, INRADFLAG, NEWGRIDFLAG, DELTAM, ACCELFLAG
    logical :: BASEOUT, HIGHORDERRAD, MASTERPROC, VALIDBEAM, BTEST
    logical :: RADIANCE_OUT, VISUAL_OUT
    logical :: FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
    real :: SOLARFLUX, SOLARMU, SOLARAZ, SOLFLUX
    real :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real :: SOLACC, SOLCRIT, SPLITACC, SHACC
    real :: DELXSFC, DELYSFC, MAXASYM, REDUCE
    real :: MUOUT, PHIOUT
    real :: UNIFZLEV, XO, YO, ZO, DIRPATH
    real :: CPUTIME1, CPUTIME2, CPUTIMETOTAL, SUM_CPU_TIME
    character :: SRCTYPE*1, UNITS*1
    character :: GRIDTYPE*1, PROPTYPE*1, SFCTYPE*2
    character*256 :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character*256 :: OutFileNC, RUNNAME
    external :: SUM_CPU_TIME

    logical :: nmlMode = .true.

    ! Choose the option to compute higher order radiance SH terms
    HIGHORDERRAD = .false.
    NSHOUT = 4
    if (HIGHORDERRAD) NSHOUT = 5
    ! Choose type of initialization/output for k-distributions:
    !   True means initiatize with and output only the base grid.
    BASEOUT = .false.

    call cpu_time(CPUTIME1)

    ! Start up MPI for multiple processors if available, and
    !   return flag for master process (procnum=0)
    call START_MPI(MASTERPROC)

    if (MASTERPROC) then
        if (nmlMode) then
            call shd_io__nmlInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                INSAVEFILE, OUTSAVEFILE,&
                NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
                IPFLAG, KDIST, DELTAM, GRIDTYPE,&
                SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
                ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
                MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
                OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
        else
            call shd_io__userInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                INSAVEFILE, OUTSAVEFILE,&
                NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
                IPFLAG, KDIST, DELTAM, GRIDTYPE,&
                SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
                ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
                MAXOUT,MAXPAR, NUMOUT,OUTTYPES,OUTPARMS,OUTFILES,&
                OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
        end if
    end if

    call BROADCAST_USER_INPUT(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
        &INSAVEFILE, OUTSAVEFILE,&
        &NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
        &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
        &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
        &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
        &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
        &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    NSTOKES = max(1, min(4, NSTOKES))
    if (NSTOKES .eq. 2) NSTOKES = 3

    NX = max(1, NX); NY = max(1, NY); NZ = max(2, NZ)
    ! If single plane or column then force IP mode
    if (NX .eq. 1) IPFLAG = ior(IPFLAG, 1)
    if (NY .eq. 1) IPFLAG = ior(IPFLAG, 2)

    ! Make ML and MM from NMU and NPHI
    NMU = max(2, 2*int((NMU + 1)/2))
    NPHI = max(1, NPHI)
    ML = NMU - 1
    MM = max(0, int(NPHI/2) - 1)
    NLEG = ML
    if (DELTAM) NLEG = ML + 1

    ! Make ML and MM from NMU and NPHI
    ! Compute NLM, NPHI0MAX.
    if (SRCTYPE .eq. 'T') SOLARAZ = 0.0
    NPHI0MAX = NPHI
    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    MEMWORD = NMU*(2 + 2*NPHI + 2*NLM + 2*33*32)

    ! NSTLEG is the number of phase function elements used (1 or 6)
    ! Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (NSTOKES .eq. 1) then
        NSTLEG = 1
    else
        NSTLEG = 6
    end if

    allocate (MU(NMU), WTDO(NMU*NPHI), PHI(NMU*NPHI), NPHI0(NMU))

    ! Find out the size of the input property file grid
    if (MASTERPROC) then
        call shd_io__readPropSize(PROPFILE, NSTLEG, NLEG, NPX, NPY, NPZ,&
                                NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
    end if
    call BROADCAST_PROPERTY_SIZE(NPX, NPY, NPZ, DELX, DELY,&
                                 NUMPHASE, MAXLEG, MAXPGL)

    ! Set the full domain grid sizes (variables end in T)
    NXT = NX; NYT = NY; NPXT = NPX; NPYT = NPY

    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and NX,NY to the subdomain size, and output
    !   the X and Y starting position of this subdomain.  Set bits
    !   2 and 3 in BCFLAG to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    call MAP_SHDOM_MPI(BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY, PROPFILE,&
        &XSTART, YSTART, RUNNAME)

    ! Set up base grid point actual size (NX1xNY1xNZ)
    NX1 = NX + 1
    if (btest(IPFLAG, 0) .or. btest(BCFLAG, 0) .or. btest(BCFLAG, 2)) NX1 = NX
    NY1 = NY + 1
    if (btest(IPFLAG, 1) .or. btest(BCFLAG, 1) .or. btest(BCFLAG, 3)) NY1 = NY
    if (NZ .le. 1) then
        write (6, '(A)') 'Error [shdom_main]: NZ < 2.'
        stop
    end if
    NBPTS = NX1*NY1*NZ
    ! Calculate the number of base grid cells depending on the BCFLAG
    NBCELLS = (NZ - 1)*(NX + ibits(BCFLAG, 0, 1) - ibits(BCFLAG, 2, 1))&
        *(NY + ibits(BCFLAG, 1, 1) - ibits(BCFLAG, 3, 1))

    ! Determine the size of the property grid arrays so can allocate
    MAXPG = NPX*NPY*NPZ
    if (NUMPHASE .eq. 0) MAXPGL = MAXPG*(MAXLEG + 1)*NSTLEG

    allocate (ZLEVELS(NPZ))
    allocate (TEMPP(MAXPG), EXTINCTP(MAXPG), ALBEDOP(MAXPG))
    allocate (EXTDIRP(MAXPG), IPHASEP(MAXPG), LEGENP(MAXPGL))

    MEMWORD = MEMWORD + 4.5*MAXPG + MAXPGL + NUMPHASE*NSTLEG*(MAXLEG + 1)

    MAXNZ = NPZ
    if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then

        ! Read the properties of the medium to the master processor and
        !  distribute to all the processors
        if (MASTERPROC) then
            MAXPGT = NPXT*NPYT*NPZ
            MAXPGLT = MAXPGL
            if (NUMPHASE .eq. 0) MAXPGLT = MAXPGT*(MAXLEG + 1)*NSTLEG

            allocate (TEMPPT(MAXPGT), EXTINCTPT(MAXPGT), ALBEDOPT(MAXPGT))
            allocate (IPHASEPT(MAXPGT), LEGENPT(MAXPGLT))

            call shd_io__readProp(PROPFILE, NPXT, NPYT, NPZ,&
                &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPGT, MAXPGLT, DELTAM,&
                &PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM,&
                &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE, IPHASEPT)

            if (NUMPHASE == 0) then
                write (6, '(A)') 'Error [shdom_main]: Standard property file (with phase function Legendre series at every point)'
                write (6, '(A)') '  is not allowed with multiple processors.'
                stop
            end if
            ! Check values of medium properties
            call shd_io__checkPropInput(NPXT, NPYT, NPZ, NSTLEG, NLEG,&
                &DELX, DELY, ZLEVELS,&
                &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE, IPHASEPT)
        end if

        call SCATTER_PROPERTIES(NPXT, NPYT, NPX, NPY, NPZ,&
            &NSTLEG, NLEG, NUMPHASE,&
            &ZLEVELS, MAXASYM,&
            &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,&
            &TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)

        if (MASTERPROC) deallocate (TEMPPT, EXTINCTPT, ALBEDOPT, IPHASEPT)

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readProp(PROPFILE, NPX, NPY, NPZ,&
            &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,&
            &PROPTYPE, DELX, DELY, ZLEVELS, MAXASYM,&
            &TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)

        ! Check values of medium properties
        call shd_io__checkPropInput(NPX, NPY, NPZ, NSTLEG, NLEG,&
            &DELX, DELY, ZLEVELS,&
            &TEMPP, EXTINCTP, ALBEDOP, LEGENP, NUMPHASE, IPHASEP)
    end if

    if (MASTERPROC) then
        ! Check input parameters
        call shd_io__checkInputParams(NMU, NPHI, DELTAM, MAXASYM,&
            &GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU,&
            &WAVENO, WAVELEN, GNDTEMP, GNDALBEDO,&
            &SPLITACC, SHACC)
    end if

    if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
        ! Read the memory parameters from a file, if it exists, and send
        !  the parameters to all processors. Otherwise use the input ones.
        call READ_BROADCAST_MEM_PARAMS(MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
            &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, RUNNAME)
    end if

    ! Guess maximum number of grid points, cells, SH vector size needed
    ! but don't let MAX_TOTAL_MB be exceeded
    if (MAX_TOTAL_MB*1024.0**2 > 1.80*huge(MAXIV)) then
        MAX_TOTAL_MB = 1.80*huge(MAXIV)/1024.0**2
        write (6, '(A,F10.2)') 'Message [shdom_main]: MAX_TOTAL_MB reduced to fit memory model:', MAX_TOTAL_MB
    end if

    if (SPLITACC .le. 0.0 .and. INSAVEFILE(1:4) .eq. 'NONE') then
        ADAPT_GRID_FACTOR = 1.0
    end if

    BIG_ARRAYS = 3
    if (.not. ACCELFLAG) BIG_ARRAYS = 2

    WANTMEM = ADAPT_GRID_FACTOR*NBPTS*(28 + 16.5*CELL_TO_POINT_RATIO&
        &+ NPHI0MAX*NSTOKES + NUM_SH_TERM_FACTOR*NSTOKES*NLM*BIG_ARRAYS)
    REDUCE = min(1.0, ((MAX_TOTAL_MB*1024**2)/4 - MEMWORD)/WANTMEM)

    ADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR*REDUCE
    if (ADAPT_GRID_FACTOR < 1.0) then
        write (6, '(A)') 'Error [shdom_main]: MAX_TOTAL_MB memory limit exceeded with just base grid.'
        stop
    end if

    if (REDUCE < 1.0) then
        write (6, '(A,F8.2)') 'Message [shdom_main]: ADAPT_GRID_FACTOR reduced to:', ADAPT_GRID_FACTOR
    end if

    WANTMEM = REDUCE*WANTMEM
    if (WANTMEM > huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom_main]: Number of words of memory exceeds max integer size:', WANTMEM
        stop
    end if

    MAXIG = nint(ADAPT_GRID_FACTOR*NBPTS)
    MAXIC = nint(CELL_TO_POINT_RATIO*MAXIG)
    MAXIV = nint(NUM_SH_TERM_FACTOR*NLM*MAXIG)
    MAXIDO = MAXIG*NPHI0MAX
    if (NUMPHASE .gt. 0) then
        MAXIGL = NUMPHASE*NSTLEG*(MAXLEG + 1)
    else
        MAXIGL = MAXIG*NSTLEG*(MAXLEG + 1)
    end if
    if (4.0*(MAXIV + MAXIG)*NSTOKES > huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom_main]: Size of big SH arrays (MAXIV) probably exceeds max integer number of bytes:', MAXIV
        stop
    end if
    if (4.0*8.0*MAXIC > huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom_main]: Size of GRIDPTR array (8*MAXIC) probably exceeds max integer number of bytes:', 8*MAXIC
        stop
    end if

    ! Allocate the internal adaptable grid arrays, including the BIG 3
    allocate (XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ))
    allocate (TEMP(MAXIG), PLANCK(MAXIG))
    allocate (EXTINCT(MAXIG), ALBEDO(MAXIG), IPHASE(MAXIG))
    allocate (LEGEN(MAXIGL))
    allocate (GRIDPTR(8, MAXIC), NEIGHPTR(6, MAXIC), TREEPTR(2, MAXIC))
    allocate (CELLFLAGS(MAXIC))
    allocate (GRIDPOS(3, MAXIG))
    allocate (RSHPTR(MAXIG + 2), SHPTR(MAXIG + 1), OSHPTR(MAXIG + 1))
    allocate (RADIANCE(NSTOKES, MAXIV + MAXIG), SOURCE(NSTOKES, MAXIV))
    if (ACCELFLAG) allocate (DELSOURCE(NSTOKES, MAXIV))
    allocate (WORK(NSTOKES*MAXIDO), WORK1(8*MAXIG))
    allocate (FLUXES(2*MAXIG), DIRFLUX(MAXIG))

    MAXNBC = MAXIG*3/NZ
    ! Read in the surface properties
    if (SFCFILE(1:2) .eq. 'NO') then
        SFCTYPE = 'FL'
    else
        if (MASTERPROC) then
            call shd_io__readSurfSize(SFCFILE, MAXSFCPTS, MAXSFCPARS)
        end if
        ! Broadcast the surface array sizes
        call BROADCAST_SURFACE_SIZE(MAXSFCPTS, MAXSFCPARS)
        allocate (SFCPARMS(MAXSFCPARS*MAXSFCPTS))
        if (MASTERPROC) then
            call shd_io__readSurf(SFCFILE, MAXSFCPTS, MAXSFCPARS,&
                &SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,&
                &NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
        end if
        ! Broadcast the surface array data
        call BROADCAST_SURFACE_PARMS(SFCTYPE, NXSFC, NYSFC, NSFCPAR,&
            &DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
        allocate (SFCGRIDPARMS(MAXSFCPARS*MAXNBC))
    end if

    if (SFCTYPE(2:2) .eq. 'L') then
        MAXBCRAD = 2*MAXNBC
    else
        MAXBCRAD = (2 + NMU*NPHI0MAX/2)*MAXNBC
    end if
    allocate (BCPTR(2*MAXNBC), BCRAD(NSTOKES, MAXBCRAD))

    ! If doing a k-distribution then get the band info from the CKD file
    if (KDIST) then
        if (MASTERPROC) then
            call shd_io__readCkdSize(CKDFILE, WAVENO, NG, NZCKD)
        end if
        ! Broadcast the k-distribution array sizes
        call BROADCAST_KDIST_SIZE(NG, NZCKD)
        allocate (ZCKD(NZCKD), GASABS(NZCKD), DELG(NG), KABS(NZCKD, NG))
        if (MASTERPROC) then
            MAXNG = NG; MAXNZ = NZCKD
            call shd_io__readCkd(MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,&
                &NG, DELG, NZCKD, ZCKD, KABS)
        end if
        call BROADCAST_KDIST_PARMS(SOLFLUX, NG, DELG, NZCKD, ZCKD, KABS)
        SOLARFLUX = SOLFLUX*SOLARFLUX*abs(SOLARMU)
    else
        NG = 1
        allocate (DELG(NG))
        DELG(1) = 1.0
        NZCKD = 0
        BASEOUT = .false.
    end if

    ! Figure out what type of output we will be generating
    NRAD = 0
    NVIS = 0
    RADIANCE_OUT = .false.
    VISUAL_OUT = .false.
    SOURCE_OUT = .false.
    FLUX_OUT = .false.
    FLUXDIV_OUT = .false.
    SH_OUT = .false.
    do I = 1, NUMOUT
        if (OUTTYPES(I) .eq. 'R') then
            RADIANCE_OUT = .true.
            if (OUTPARMS(2, I) .le. 0.0) then
                NXOUT = 1
            else
                if (btest(BCFLAG, 2)) then
                    NXOUT = max(1, nint(NPXT*DELX/OUTPARMS(2, I)))
                else
                    NXOUT = max(1, nint((NPX*DELX + 2*OUTPARMS(4, I))/OUTPARMS(2, I)))
                end if
            end if
            if (OUTPARMS(3, I) .le. 0.0) then
                NYOUT = 1
            else
                if (btest(BCFLAG, 3)) then
                    NYOUT = max(1, nint(NPYT*DELY/OUTPARMS(3, I)))
                else
                    NYOUT = max(1, nint((NPY*DELY + 2*OUTPARMS(5, I))/OUTPARMS(3, I)))
                end if
            end if
            NRAD = NRAD + NSTOKES*NXOUT*NYOUT*nint(OUTPARMS(6, I))
        else if (OUTTYPES(I) .eq. 'V') then
            VISUAL_OUT = .true.
            if (nint(OUTPARMS(1, I)) .eq. 1) then
                NVIS = NVIS + NSTOKES*nint(OUTPARMS(10, I))*nint(OUTPARMS(11, I))
            else
                NSCAN = 1 + sqrt((OUTPARMS(4, I) - OUTPARMS(7, I))**2&
                    &+ (OUTPARMS(5, I) - OUTPARMS(8, I))**2&
                    &+ (OUTPARMS(6, I) - OUTPARMS(9, I))**2)/OUTPARMS(10, I)
                NPIX = 1 + abs(OUTPARMS(12, I) - OUTPARMS(11, I))/OUTPARMS(13, I)
                NVIS = NVIS + NSTOKES*NSCAN*NPIX
            end if
        else if (OUTTYPES(I) .eq. 'F') then
            if (.not. FLUX_OUT) then
                allocate (SUMFLUXES(2*MAXIG), SUMDIRFLUX(MAXIG))
            end if
            FLUX_OUT = .true.
        else if (OUTTYPES(I) .eq. 'H') then
            if (.not. FLUXDIV_OUT) then
                allocate (SUMFLUXDIV(MAXIG))
            end if
            FLUXDIV_OUT = .true.
        else if (OUTTYPES(I) .eq. 'S') then
            if (.not. SH_OUT) then
                allocate (SUMSHTERMS(NSHOUT*MAXIG))
            end if
            SH_OUT = .true.
        else if (OUTTYPES(I) .eq. 'J') then
            if (.not. SOURCE_OUT) then
                allocate (SOURCE1OUT(MAXIG*NSTOKES))
            end if
            SOURCE_OUT = .true.
        end if
    end do

    allocate (SUMRADOUT(NRAD), SUMVISOUT(NVIS))
    allocate (WORK2(max(NVIS, NRAD, NSTOKES*MAXIG)))

    ! Notify user that certain outputs are not supported with multiple processors
    if ((btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) .and. MASTERPROC) then
        do I = 1, NUMOUT
            if (OUTTYPES(I) .eq. 'F') then
                if (nint(OUTPARMS(1, I)) .eq. 2 .or. nint(OUTPARMS(1, I)) .eq. 5) then
                    write (6, '(A,I1,A,A)') 'F output format ', nint(OUTPARMS(1, I)),&
                        &' not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
                end if
            else if (OUTTYPES(I) .eq. 'H') then
                if (nint(OUTPARMS(1, I)) .eq. 3) then
                    write (6, '(A,I1,A,A)') 'H output format ', nint(OUTPARMS(1, I)),&
                        &' not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
                end if
            else if (OUTTYPES(I) .eq. 'S') then
                if (nint(OUTPARMS(1, I)) .eq. 2) then
                    write (6, '(A,I1,A,A)') 'S output format ', nint(OUTPARMS(1, I)),&
                        &' not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
                end if
            else if (OUTTYPES(I) .eq. 'J') then
                write (6, '(A,A)') 'J output format not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
            else if (OUTTYPES(I) .eq. 'M') then
                write (6, '(A,A)') 'M output format not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
            end if
        end do
    end if

    ! Read in previous run from binary file if desired
    call shd_io__restoreState(INSAVEFILE, NSTOKES, NX, NY, NZ, ML, MM, NLM,&
        &INRADFLAG, NEWGRIDFLAG, XGRID, YGRID, ZGRID,&
        &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
        &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

    ! check if we can use pre-calculated "state", if not, generate new grids
    ! NEWGRIDFLAG is from shd_io__restoreState call
    if (NEWGRIDFLAG) then
        ! Make the internal base grid lines
        call NEW_GRIDS(BCFLAG, GRIDTYPE, NPX, NPY, NPZ, NX, NY, NZ,&
            &XSTART, YSTART, DELX, DELY, ZLEVELS,&
            &XGRID, YGRID, ZGRID)
    end if

    write (6, '(A,"<",A,">")') 'Message [shdom_main]: Starting solution iterations for SHDOM run: ', trim(adjustl(RUNNAME))
    TOTITER = 0
    OLDNPTS = 0

    ! Loop over the k-distribution g's from high to low absorption
    do IG = NG, 1, -1
        ! If doing k-distribution then get the gas absorption profile.
        if (KDIST) then
            do I = 1, NZCKD
                GASABS(I) = KABS(I, IG)
            end do
        end if

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (BASEOUT .or. (IG .eq. NG .and. NEWGRIDFLAG)) then
            call INIT_CELL_STRUCTURE(BCFLAG, IPFLAG,&
                &NX, NY, NZ, NX1, NY1, NPTS, NCELLS,&
                &XGRID, YGRID, ZGRID, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
            NBCELLS = NCELLS
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in SOURCE (SHPTR), RADIANCE (RSHPTR), FLUXES, DIRFLUX.
        call shdSolv__solveRte(NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
            &ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE,&
            &NPHI0, MU, PHI, WTDO,&
            &MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG,&
            &BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD,&
            &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
            &SFCTYPE, GNDTEMP, GNDALBEDO,&
            &NXSFC, NYSFC, DELXSFC, DELYSFC,&
            &NSFCPAR, SFCPARMS, SFCGRIDPARMS,&
            &UNITS, WAVENO, WAVELEN,&
            &ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,&
            &SPLITACC, SHACC, XGRID, YGRID, ZGRID,&
            &TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,&
            &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
            &NPTS, GRIDPOS,&
            &NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
            &RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,&
            &SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
        INRADFLAG = .false.
        TOTITER = TOTITER + ITER

        ! Compute the output and put in arrays.  First the radiance.
        if (RADIANCE_OUT) then
            IRAD = 0
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'R') then
                    if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
                        call COMPUTE_RADIANCE_PAR(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &BCFLAG, NPXT*DELX, NPYT*DELY, IPFLAG,&
                            &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                            &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                            &SHPTR, SOURCE, WORK, WORK(1 + NPTS*NSTOKES),&
                            &OUTPARMS(1, I), IRAD, WORK2)
                    else
                        call COMPUTE_RADIANCE(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                            &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                            &SHPTR, SOURCE, WORK(1), WORK(1 + NPTS*NSTOKES),&
                            &OUTPARMS(1, I), IRAD, WORK2)
                    end if
                end if
            end do
            call SUM_OUTPUT(NG - IG, DELG(IG), NRAD, WORK2, SUMRADOUT)
        end if

        if (VISUAL_OUT) then
            IVIS = 0
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'V') then
                    if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
                        call VISUALIZE_RADIANCE_PAR(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &BCFLAG, NPXT*DELX, NPYT*DELY, IPFLAG,&
                            &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                            &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                            &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2)
                    else
                        call VISUALIZE_RADIANCE(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                            &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                            &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2)
                    end if
                end if
            end do
            call SUM_OUTPUT(NG-IG, DELG(IG), NVIS, WORK2, SUMVISOUT)
        end if

        if (SOURCE_OUT) then
            ! Compute source function throughout grid for this one angle
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'J') then
                    MUOUT = OUTPARMS(2, I)
                    PHIOUT = OUTPARMS(3, I)*acos(-1.0)/180.0
                    call COMPUTE_ONE_SOURCE(NSTOKES, ML, MM, NLM,&
                        &NSTLEG, NLEG, NUMPHASE,&
                        &NPTS, DELTAM, MUOUT, PHIOUT,&
                        &SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE,&
                        &DIRFLUX, SHPTR, SOURCE, SOURCE1OUT)
                end if
            end do
        end if

        if (BASEOUT) then
            ! If outputting only the base grid then throw out the other points
            call INIT_CELL_STRUCTURE(BCFLAG, IPFLAG,&
                &NX, NY, NZ, NX1, NY1,&
                &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
        else
            ! For the output from the grid point arrays (Flux, Net flux, SH)
            !   we may need to interpolate the running k-distribution sum
            !   to the newly created grid points.
            call INTERP_OUTPUT(OLDNPTS, NPTS, FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
                SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
                NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        end if

        if (FLUX_OUT) then
            call SUM_OUTPUT(NG - IG, DELG(IG), 2*NPTS, FLUXES, SUMFLUXES)
            call SUM_OUTPUT(NG - IG, DELG(IG), NPTS, DIRFLUX, SUMDIRFLUX)
        end if
        if (FLUXDIV_OUT) then
            ! Compute the net flux divergence at every point
            call COMPUTE_NETFLUXDIV(NSTOKES, NPTS, RSHPTR, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, PLANCK, DIRFLUX, RADIANCE, WORK)
            call SUM_OUTPUT(NG - IG, DELG(IG), NPTS, WORK, SUMFLUXDIV)
        end if

        if (SH_OUT) then
            call COMPUTE_SH(NSHOUT, NSTOKES, NPTS, SRCTYPE, SOLARMU, SOLARAZ,&
                &DIRFLUX, RSHPTR, ML, MM, RADIANCE, WORK)
            call SUM_OUTPUT(NG - IG, DELG(IG), NSHOUT*NPTS, WORK, SUMSHTERMS)
        end if

        OLDNPTS = NPTS
    end do

    call shd_io__saveState(OUTSAVEFILE, NSTOKES, NX, NY, NZ, ML, MM, NLM,&
        &WORK, XGRID, YGRID, ZGRID,&
        &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
        &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

    call cpu_time(CPUTIME2)
    CPUTIMETOTAL = SUM_CPU_TIME(CPUTIME2 - CPUTIME1)
    if (MASTERPROC) then
        if (CPUTIMETOTAL .lt. 60.0) then
            write (6, '(A,F6.2,A)') 'Message [shdom_main]: Total CPU time used:', CPUTIMETOTAL, ' sec'
        else
            write (6, '(A,F6.2,A)') 'Message [shdom_main]: Total CPU time used:', CPUTIMETOTAL/60.0, ' min'
        endif
    end if

    NSH = SHPTR(NPTS + 1)
    IRAD = 0
    IVIS = 0
    if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then

        if (MASTERPROC) then
            if (FLUX_OUT) allocate (ALLFLUXES(3, NZ, NYT, NXT))
            if (FLUXDIV_OUT) allocate (ALLFLUXDIV(NZ, NYT, NXT))
            if (SH_OUT) allocate (ALLSHTERMS(NSHOUT, NZ, NYT, NXT))
            allocate (XALLGRID(NXT), YALLGRID(NYT))
        end if

        call GATHER_OUTPUT(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
            &IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,&
            &SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
            &NPXT, NPYT, DELX, DELY, XALLGRID, YALLGRID,&
            &ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,&
            &NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)

        if (MASTERPROC) then
            do I = 1, NUMOUT
                if (OUTTYPES(I) /= 'V' .and. OUTFILES(I) /= 'NONE') then
                    call shd_io__outputResultMpi(NSTOKES, NXT, NYT, NZ,&
                        &NPTSTOT, NCELLSTOT, NSHTOT,&
                        &ML, MM, NLM, NMU, NPHI, NANG, NG,&
                        &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                        &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                        &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                        &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                        &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                        &NPXT*DELX, NPYT*DELY, XALLGRID, YALLGRID, ZGRID,&
                        &ALLFLUXES, ALLFLUXDIV,&
                        &NSHOUT, ALLSHTERMS, IRAD, SUMRADOUT,&
                        &OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I))
                end if
            end do

            if (FLUX_OUT) deallocate (ALLFLUXES)
            if (FLUXDIV_OUT) deallocate (ALLFLUXDIV)
            if (SH_OUT) deallocate (ALLSHTERMS)
            deallocate (XALLGRID, YALLGRID)
        end if

    else

        do I = 1, NUMOUT
            if (OUTTYPES(I) /= 'V' .and. OUTFILES(I) /= 'NONE') then
                call shd_io__outputResult(NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NCELLS,&
                    NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, NMU, NPHI, NANG, NG,&
                    PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                    BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                    SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                    SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                    SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                    XGRID, YGRID, ZGRID, GRIDPOS,&
                    TREEPTR, GRIDPTR, CELLFLAGS,&
                    EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,&
                    SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,&
                    IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT,&
                    OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I))
            end if
        end do
    end if

    ! Output visualization images if there are any
    if (MASTERPROC) then
        do I = 1, NUMOUT
            if (OUTTYPES(I) == 'V') then
                call shd_io__outputImg(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                    &NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
                    &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                    &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                    &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                    &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                    &IVIS, SUMVISOUT, OUTPARMS(1, I), OUTFILES(I))
            end if
        end do
    end if

    ! Calculate the actual memory parameters and output them to the log file
    MAXMB_OUT = 4*(NMU*(2 + 2*NPHI + 2*NLM + 2*33*32)&
        &+ 4.5*MAXPG + MAXPGL + NUMPHASE*(MAXLEG + 1)&
        &+ 16.5*NCELLS + NPTS*(28 + NPHI0MAX*NSTOKES)&
        &+ NSTOKES*NSH*BIG_ARRAYS)/1024**2
    ADAPTGRIDFAC_OUT = float(NPTS)/NBPTS
    SHTERM_FAC_OUT = float(NSH)/(NLM*NPTS)
    CELL_POINT_OUT = float(NCELLS)/NPTS
    write (6, '(A)')'Message [shdom_main]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    write (6, '(28X,F7.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    call END_SHDOM_MPI(NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART, DELX, DELY,&
        &NPXT, NPYT, PROPFILE, RUNNAME)

    if (SRCTYPE .ne. 'T') then
        call DIRECT_BEAM_PROP(9, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
            &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end if

    deallocate (SUMRADOUT, SUMVISOUT, WORK2)
    if (FLUX_OUT) deallocate (SUMFLUXES, SUMDIRFLUX)
    if (FLUXDIV_OUT) deallocate (SUMFLUXDIV)
    if (SH_OUT) deallocate (SUMSHTERMS)
    if (SOURCE_OUT) deallocate (SOURCE1OUT)
    if (KDIST) deallocate (ZCKD, GASABS, KABS)
    deallocate (BCPTR, BCRAD, DELG)
    if (allocated(SFCPARMS)) deallocate (SFCPARMS, SFCGRIDPARMS)
    deallocate (RADIANCE, SOURCE)
    if (ACCELFLAG) deallocate (DELSOURCE)
    deallocate (RSHPTR, SHPTR, OSHPTR, WORK, WORK1, FLUXES, DIRFLUX)
    deallocate (GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS)
    deallocate (TEMP, PLANCK, EXTINCT, ALBEDO, IPHASE, LEGEN)
    deallocate (XGRID, YGRID, ZGRID)
    deallocate (ZLEVELS, TEMPP, EXTINCTP, ALBEDOP, EXTDIRP, IPHASEP, LEGENP)
    deallocate (MU, WTDO, PHI, NPHI0)

end program shdom_main

subroutine TRILIN_INTERP_PROP(X, Y, Z, INIT, NSTLEG, NLEG,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
    !  Trilinearly interpolates the quantities on the input property
    ! grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
    ! ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the
    ! volume coefficients.  Also adds in the separate gaseous absorption.
    ! Divides the phase function Legendre coefficients LEGEN by 2*l+1.
    ! The phase function pointer IPHASE (for tabulated phase functions)
    ! is that of the maximum weighted scattering property grid point.
    ! If INIT=.TRUE. then transfers the tabulated phase functions.

    use shd_common_arr
    implicit none
    integer NSTLEG, NLEG
    integer*2 IPHASE
    logical INIT
    real X, Y, Z, TEMP, EXTINCT, ALBEDO, LEGEN(NSTLEG, 0:NLEG, *)
    integer IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
    integer I1, I2, I3, I4, I5, I6, I7, I8, I
    double precision U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
    double precision SCAT1, SCAT2, SCAT3, SCAT4, SCAT5, SCAT6, SCAT7, SCAT8
    double precision SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN
    save EXTMIN, SCATMIN

    if (INIT) then
       ! If there are tabulated phase functions, then transfer them
        do I = 1, NUMPHASE
            do L = 0, NLEG
                do J = 1, NSTLEG
                    LEGEN(J, L, I) = LEGENP(J + NSTLEG*(L + (NLEG + 1)*(I - 1)))/(2*L + 1)
                end do
            end do
        end do
        EXTMIN = 1.0e-5/((ZLEVELS(NPZ) - ZLEVELS(1))/NPZ)
        SCATMIN = 0.1*EXTMIN
        return
    end if

    ! Find the grid location and compute the interpolation factors
    IL = 0
    IU = NPZ
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZLEVELS(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)
    W = dble(Z - ZLEVELS(IZ))/(ZLEVELS(IZ + 1) - ZLEVELS(IZ))
    W = max(min(W, 1.0d0), 0.0d0)
    IX = int((X - XSTART)/DELX) + 1
    if (abs(X - XSTART - NPX*DELX) .lt. 0.01*DELX) IX = NPX
    if (IX .lt. 1 .or. IX .gt. NPX) then
        write (6, *) 'TRILIN: Beyond X domain'
        stop
    end if
    IXP = mod(IX, NPX) + 1
    U = dble(X - XSTART - DELX*(IX - 1))/DELX
    U = max(min(U, 1.0d0), 0.0d0)
    if (U .lt. 1.0d-5) U = 0.0d0
    if (U .gt. 1.0d0 - 1.0d-5) U = 1.0d0
    IY = int((Y - YSTART)/DELY) + 1
    if (abs(Y - YSTART - NPY*DELY) .lt. 0.01*DELY) IY = NPY
    if (IY .lt. 1 .or. IY .gt. NPY) then
        write (6, *) 'TRILIN: Beyond Y domain'
        stop
    end if
    IYP = mod(IY, NPY) + 1
    V = dble(Y - YSTART - DELY*(IY - 1))/DELY
    V = max(min(V, 1.0d0), 0.0d0)
    if (V .lt. 1.0d-5) V = 0.0d0
    if (V .gt. 1.0d0 - 1.0d-5) V = 1.0d0

    F1 = (1 - U)*(1 - V)*(1 - W)
    F2 = U*(1 - V)*(1 - W)
    F3 = (1 - U)*V*(1 - W)
    F4 = U*V*(1 - W)
    F5 = (1 - U)*(1 - V)*W
    F6 = U*(1 - V)*W
    F7 = (1 - U)*V*W
    F8 = U*V*W
    I1 = IZ + NPZ*(IY - 1) + NPZ*NPY*(IX - 1)
    I2 = IZ + NPZ*(IY - 1) + NPZ*NPY*(IXP - 1)
    I3 = IZ + NPZ*(IYP - 1) + NPZ*NPY*(IX - 1)
    I4 = IZ + NPZ*(IYP - 1) + NPZ*NPY*(IXP - 1)
    I5 = I1 + 1
    I6 = I2 + 1
    I7 = I3 + 1
    I8 = I4 + 1

    ! Trilinearly interpolate the temperature, extinction, scattering
    TEMP = F1*TEMPP(I1) + F2*TEMPP(I2) + F3*TEMPP(I3) + F4*TEMPP(I4)&
        &+ F5*TEMPP(I5) + F6*TEMPP(I6) + F7*TEMPP(I7) + F8*TEMPP(I8)
    EXTINCT = F1*EXTINCTP(I1) + F2*EXTINCTP(I2) + F3*EXTINCTP(I3)&
        &+ F4*EXTINCTP(I4) + F5*EXTINCTP(I5) + F6*EXTINCTP(I6)&
        &+ F7*EXTINCTP(I7) + F8*EXTINCTP(I8)
    SCAT1 = F1*EXTINCTP(I1)*ALBEDOP(I1)
    SCAT2 = F2*EXTINCTP(I2)*ALBEDOP(I2)
    SCAT3 = F3*EXTINCTP(I3)*ALBEDOP(I3)
    SCAT4 = F4*EXTINCTP(I4)*ALBEDOP(I4)
    SCAT5 = F5*EXTINCTP(I5)*ALBEDOP(I5)
    SCAT6 = F6*EXTINCTP(I6)*ALBEDOP(I6)
    SCAT7 = F7*EXTINCTP(I7)*ALBEDOP(I7)
    SCAT8 = F8*EXTINCTP(I8)*ALBEDOP(I8)
    SCATTER = SCAT1 + SCAT2 + SCAT3 + SCAT4 + SCAT5 + SCAT6 + SCAT7 + SCAT8
    if (EXTINCT .gt. EXTMIN) then
        ALBEDO = SCATTER/EXTINCT
    else
        ALBEDO = SCATTER/EXTMIN
    end if
    ! For tabulated phase functions pick the one we are on top of
    !   or the one with the most scattering weight.
    if (NUMPHASE .gt. 0) then
        MAXSCAT = -1.0
        if (SCAT1 .gt. MAXSCAT .or. abs(F1 - 1) .lt. 0.001) then
            MAXSCAT = SCAT1
            IPHASE = IPHASEP(I1)
        end if
        if (SCAT2 .gt. MAXSCAT .or. abs(F2 - 1) .lt. 0.001) then
            MAXSCAT = SCAT2
            IPHASE = IPHASEP(I2)
        end if
        if (SCAT3 .gt. MAXSCAT .or. abs(F3 - 1) .lt. 0.001) then
            MAXSCAT = SCAT3
            IPHASE = IPHASEP(I3)
        end if
        if (SCAT4 .gt. MAXSCAT .or. abs(F4 - 1) .lt. 0.001) then
            MAXSCAT = SCAT4
            IPHASE = IPHASEP(I4)
        end if
        if (SCAT5 .gt. MAXSCAT .or. abs(F5 - 1) .lt. 0.001) then
            MAXSCAT = SCAT5
            IPHASE = IPHASEP(I5)
        end if
        if (SCAT6 .gt. MAXSCAT .or. abs(F6 - 1) .lt. 0.001) then
            MAXSCAT = SCAT6
            IPHASE = IPHASEP(I6)
        end if
        if (SCAT7 .gt. MAXSCAT .or. abs(F7 - 1) .lt. 0.001) then
            MAXSCAT = SCAT7
            IPHASE = IPHASEP(I7)
        end if
        if (SCAT8 .gt. MAXSCAT .or. abs(F8 - 1) .lt. 0.001) then
            MAXSCAT = SCAT8
            IPHASE = IPHASEP(I8)
        end if
    else
        ! Standard property file format: unpolarized
        LEGEN(1, 0, 1) = 1.0
        do L = 1, NLEG
            LEGEN(1, L, 1) = SCAT1*LEGENP(L + (NLEG + 1)*(I1 - 1))&
                &+ SCAT2*LEGENP(L + (NLEG + 1)*(I2 - 1))&
                &+ SCAT3*LEGENP(L + (NLEG + 1)*(I3 - 1))&
                &+ SCAT4*LEGENP(L + (NLEG + 1)*(I4 - 1))&
                &+ SCAT5*LEGENP(L + (NLEG + 1)*(I5 - 1))&
                &+ SCAT6*LEGENP(L + (NLEG + 1)*(I6 - 1))&
                &+ SCAT7*LEGENP(L + (NLEG + 1)*(I7 - 1))&
                &+ SCAT8*LEGENP(L + (NLEG + 1)*(I8 - 1))
            if (SCATTER .gt. SCATMIN) then
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATTER
            else
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATMIN
            end if
            LEGEN(1, L, 1) = LEGEN(1, L, 1)/(2*L + 1)
        end do
    end if

    ! Add in the gaseous absorption to extinction and albedo
    if (NZCKD .gt. 0) then
        IL = 1
        IU = NZCKD
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (Z .le. ZCKD(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
        I = min(max(IL, 1), NZCKD - 1)
        F = (Z - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
        F = min(max(F, 0.0d0), 1.0d0)
        KG = (1.0 - F)*GASABS(I) + F*GASABS(I + 1)
        if (EXTINCT + KG .gt. 0.0) then
            ALBEDO = ALBEDO*EXTINCT/(EXTINCT + KG)
        else
            ALBEDO = 0.0
        end if
        EXTINCT = EXTINCT + KG
    end if
    return
end subroutine

subroutine DIRECT_BEAM_PROP(INIT, XI, YI, ZI, BCFLAG, IPFLAG,&
    &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
    &DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    ! Computes the direct beam flux at point (XI,YI,ZI) by integrating
    !   the extinction through the property grid.  If called with
    !   INIT=1 then the property grid extinction array, solar direction
    !   terms, and lowest uniform level are computed and saved.  If called
    !   with INIT=2 then the input lowest uniform level (UNIFZLEV) is stored.
    !   Call with INIT=0 to do the path integration and return the direct
    !   beam flux on the horizontal (DIRFLUX) at XI,YI,ZI.
    !   The DELTAM flag is used to determine whether to delta-M scale the
    !   extinction for the direct beam computation.  The extinction includes
    !   the gaseous absorption.  If the IPFLAG has bit 2 set then the direct
    !   beam tracing is done in 3D, otherwise the lower two bits determine
    !   the type of tracing: 0 for 3D tracing, 1 for XZ only tracing,
    !   2 for YZ only tracing, and 3 for Z only tracing.  If BCFLAG bits 0
    !   or 1 are set then have open boundary conditions in X and/or Y.  In
    !   this case when the ray tracing back to the sun reaches the boundary
    !   then independent pixel mode is entered so that only Z grid
    !   intersections occur.
    ! For use with multiple processors (bit 2 or 3 set in BCFLAG),
    !   the VALIDBEAM flag is returned true if the ray to the sun made it
    !   to the top of the domain before hitting the subdomain side.  If the
    !   flag is false SIDE is returned with the boundary hit (1=-X, 2=+X,
    !   3=-Y, 4=+Y).  XE,YE,ZE returns the location of the exitting ray,
    !   and path is the optical path from XI,YI,ZI to the sun.
    use shd_common_arr
    implicit none
    integer INIT, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, SIDE
    logical DELTAM, VALIDBEAM
    real XI, YI, ZI, SOLARFLUX, SOLARMU, SOLARAZ
    real DIRFLUX, UNIFZLEV, XO, YO, ZO, DIRPATH

    integer IX, IY, IZ, JZ, IL, IM, IU
    integer I, J, K, L, IPH, IP, JP, I1, I2, I3, I4
    integer IPDIRECT, DI, DJ, DK
    logical CONSTX, CONSTY, HITBOUNDARY, BTEST
    double precision EXTINCT, ALBEDO, F
    double precision SUNMU, SUNAZ, PATH
    double precision EXTBEAMCUT, UNIFORMZLEV
    double precision CX, CY, CZ, CXINV, CYINV, CZINV
    double precision EPSS, EPSZ, XDOMAIN, YDOMAIN
    double precision XOFFS, YOFFS, DELXD, DELYD
    double precision X, Y, Z, XE, YE, ZE, XP, YP, ZP
    double precision X0, X1, Y0, Y1, Z0, Z1, SO, SOX, SOY, SOZ
    double precision U0, V0, W0, U1, V1, W1, AX, AY, AZ
    double precision U0M, V0M, W0M, U1M, V1M, W1M, DU, DV, DW
    double precision E1, E2, E3, E4, E5, E6, E7, E8, A, B, C, D
    double precision B1, B2, B3, B4, B5, B6, B7, B8, C1, C2, C3, C4, C5, C6, C7, C8
    double precision UV, UMV, UVM, UMVM, UW, UMW, UWM, UMWM, VW, VMW, VWM, VMWM
    double precision VWU, VWUM, UWV, UWVM, UVW, UVWM
    real, save, allocatable :: GASEXT(:), EXTMIN(:), EXTMAX(:)

    save CX, CY, CZ, CXINV, CYINV, CZINV, DI, DJ, DK, IPDIRECT
    save DELXD, DELYD, XDOMAIN, YDOMAIN, EPSS, EPSZ, UNIFORMZLEV

    if (INIT .eq. 9) then
        deallocate (GASEXT, EXTMIN, EXTMAX)
        return
    end if

    if (INIT .eq. 1) then
        if (.not. allocated(GASEXT)) allocate (GASEXT(NPZ), EXTMIN(NPZ), EXTMAX(NPZ))
        ! Get the gaseous extinction at the property grid levels
        do IZ = 1, NPZ
            if (NZCKD .gt. 0) then
                IL = 1
                IU = NZCKD
                do while (IU - IL .gt. 1)
                    IM = (IU + IL)/2
                    if (ZLEVELS(IZ) .le. ZCKD(IM)) then
                        IL = IM
                    else
                        IU = IM
                    end if
                end do
                I = min(max(IL, 1), NZCKD - 1)
                W0 = (ZLEVELS(IZ) - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
                W0 = min(max(W0, 0.0d0), 1.0d0)
                GASEXT(IZ) = (1.0 - W0)*GASABS(I) + W0*GASABS(I + 1)
            else
                GASEXT(IZ) = 0.0
            end if
        end do
        ! First make the property grid extinction field, on which all
        !   the direct beam paths will be computed.
        IP = 0
        do IX = 1, NPX
            do IY = 1, NPY
                do IZ = 1, NPZ
                    IP = IP + 1
                    EXTINCT = EXTINCTP(IP)
                    ALBEDO = ALBEDOP(IP)
                    ! Add in the gaseous absorption to extinction
                    if (GASEXT(IZ) .gt. 0.0) then
                        ALBEDO = ALBEDO*EXTINCT/(EXTINCT + GASEXT(IZ))
                        EXTINCT = EXTINCT + GASEXT(IZ)
                    end if
                    ! Do the Delta-M scaling if needed
                    if (DELTAM) then
                        if (NUMPHASE .gt. 0) then
                            IPH = IPHASEP(IP)
                        else
                            IPH = IP
                        end if
                        L = ML + 1
                        F = LEGENP(1 + NSTLEG*(L + (NLEG + 1)*(IPH - 1)))/(2*L + 1)
                        EXTINCT = (1.0 - ALBEDO*F)*EXTINCT
                    end if
                    EXTDIRP(IP) = EXTINCT
                end do
            end do
        end do

        ! Bit 2 of IPFLAG means do the direct beam in 3D
        IPDIRECT = IPFLAG
        if (btest(IPFLAG, 2)) IPDIRECT = 0

        ! Make the ray direction (opposite to the solar beam)
        !   SOLAR?? is direction of beam, SUN?? is direction to sun (opposite)
        SUNMU = -SOLARMU
        SUNAZ = SOLARAZ + acos(-1.0)
        CX = sqrt(1.0 - SUNMU**2)*cos(SUNAZ)
        CY = sqrt(1.0 - SUNMU**2)*sin(SUNAZ)
        CZ = abs(SUNMU)
        if (abs(CX) .gt. 1.0e-6) then
            CXINV = 1.0d0/CX
        else
            CX = 0.0
            CXINV = 1.0e20
        end if
        if (abs(CY) .gt. 1.0e-6) then
            CYINV = 1.0d0/CY
        else
            CY = 0.0
            CYINV = 1.0e20
        end if
        if (abs(CZ) .gt. 1.0e-6) then
            CZINV = 1.0d0/CZ
        else
            CZ = 0.0
            CZINV = 1.0e20
        end if
        DI = nint(sign(1.0d0, CX))
        DJ = nint(sign(1.0d0, CY))
        DK = nint(sign(1.0d0, CZ))
        EPSZ = 1.0e-6*(ZLEVELS(NPZ) - ZLEVELS(1))
        EPSS = 1.0e-3*(ZLEVELS(NPZ) - ZLEVELS(1))/NPZ
        if (.not. btest(IPDIRECT, 0)) EPSS = max(EPSS, 1.0d-4*DELX)
        if (.not. btest(IPDIRECT, 1)) EPSS = max(EPSS, 1.0d-4*DELY)
        DELXD = dble(DELX)
        DELYD = dble(DELY)
        EPSS = max(0.001*DELXD, 0.001*DELYD, EPSS)
        XDOMAIN = DELXD*NPX
        if (btest(BCFLAG, 2)) XDOMAIN = DELXD*(NPX - 1)
        YDOMAIN = DELYD*NPY
        if (btest(BCFLAG, 3)) YDOMAIN = DELYD*(NPY - 1)

        ! Find the Z level above which the medium is plane-parallel
        EXTBEAMCUT = 1.0e-4
        do IZ = 1, NPZ
            EXTMIN(IZ) = 1.0e20
            EXTMAX(IZ) = 0.0
        end do
        IP = 0
        do IX = 1, NPX
            do IY = 1, NPY
                do IZ = 1, NPZ
                    IP = IP + 1
                    EXTMIN(IZ) = min(EXTINCTP(IP), EXTMIN(IZ))
                    EXTMAX(IZ) = max(EXTINCTP(IP), EXTMAX(IZ))
                end do
            end do
        end do
        JZ = 0
        do IZ = 1, NPZ
            if (EXTMAX(IZ) - EXTMIN(IZ) .gt. EXTBEAMCUT) JZ = IZ
        end do
        JZ = min(NPZ, JZ + 1)
        UNIFORMZLEV = ZLEVELS(JZ)
        UNIFZLEV = UNIFORMZLEV
        return
    ! Done with initialization
    end if

    if (INIT .eq. 2) then
        UNIFORMZLEV = UNIFZLEV
        return
    end if

    ! Here for computing the direct beam path for one starting point.
    Z = ZI
    X = XI - XSTART
    Y = YI - YSTART

    ! Find the grid location
    IL = 0
    IU = NPZ
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZLEVELS(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    K = max(IL, 1)
    I = int(X/DELXD) + 1
    if (I .gt. NPX .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) I = NPX
    if (I .lt. 1 .or. I .gt. NPX) then
        write (6, *) 'DIRECT_BEAM_PROP: Beyond X domain', I, XI, YI, ZI
        stop
    end if
    J = int(Y/DELYD) + 1
    if (J .gt. NPY .and. abs(Y - YDOMAIN) .lt. 0.001*DELYD) J = NPY
    if (J .lt. 1 .or. J .gt. NPY) then
        write (6, *) 'DIRECT_BEAM_PROP: Beyond Y domain', J, XI, YI, ZI
        stop
    end if
    XE = X
    YE = Y
    ZE = Z
    XP = XE
    YP = YE
    ZP = ZE
    CONSTX = btest(IPDIRECT, 0)
    CONSTY = btest(IPDIRECT, 1)
    if (CX .eq. 0.0) CONSTX = .true.
    if (CY .eq. 0.0) CONSTY = .true.
    if (btest(BCFLAG, 0) .and. (abs(X) .lt. 0.01*DELXD&
        &.or. abs(X - (NPX - 1)*DELXD) .lt. 0.01*DELXD)) CONSTX = .true.
    if (btest(BCFLAG, 1) .and. (abs(Y) .lt. 0.01*DELYD&
        &.or. abs(Y - (NPY - 1)*DELYD) .lt. 0.01*DELYD)) CONSTY = .true.

    ! If have multiple subdomains (processors) and ray is going outwards
    !   from a boundary then set the HITBOUNDARY flag
    HITBOUNDARY = .false.
    if (btest(BCFLAG, 2)) then
        if (CX .gt. 0.0 .and. abs(X - XDOMAIN) .lt. 0.001*DELXD) then
            SIDE = 2
            HITBOUNDARY = .true.
        else if (CX .lt. 0.0 .and. abs(X) .lt. 0.001*DELXD) then
            SIDE = 1
            HITBOUNDARY = .true.
        end if
    end if
    if (btest(BCFLAG, 3)) then
        if (CY .gt. 0.0 .and. abs(Y - YDOMAIN) .lt. 0.001*DELYD) then
            SIDE = 4
            HITBOUNDARY = .true.
        end if
        if (CY .lt. 0.0 .and. abs(Y) .lt. 0.001*DELYD) then
            SIDE = 3
            HITBOUNDARY = .true.
        end if
    end if

    ! Grid cell loop begin
    PATH = DIRPATH
    do while (.not. HITBOUNDARY .and. abs(ZE - ZLEVELS(NPZ)) .gt. EPSZ)
        IP = I + 1
        if (I .eq. NPX) then
            if (btest(BCFLAG, 0) .or. btest(BCFLAG, 2)) then
                IP = NPX
            else
                IP = 1
            end if
        end if
        JP = J + 1
        if (J .eq. NPY) then
            if (btest(BCFLAG, 1) .or. btest(BCFLAG, 3)) then
                JP = NPY
            else
                JP = 1
            end if
        end if
        X0 = DELXD*(I - 1)
        X1 = X0 + DELXD
        Y0 = DELYD*(J - 1)
        Y1 = Y0 + DELYD
        Z0 = ZLEVELS(K)
        Z1 = ZLEVELS(K + 1)
        if  (I .lt. 1 .or. I .gt. NPX .or.&
            &J .lt. 1 .or. J .gt. NPY .or.&
            &K .lt. 1 .or. K .ge. NPZ) then
            write (6, '(A,3I4)') 'DIRECT_BEAM_PROP: beyond grid!', I, J, K
            write (6, '(1(2X,3F9.5))') X, Y, Z
            stop
        end if
        ! Get the eight corner extinction values
        I1 = K + NPZ*(J - 1) + NPZ*NPY*(I - 1)
        I2 = K + NPZ*(J - 1) + NPZ*NPY*(IP - 1)
        I3 = K + NPZ*(JP - 1) + NPZ*NPY*(I - 1)
        I4 = K + NPZ*(JP - 1) + NPZ*NPY*(IP - 1)
        E1 = EXTDIRP(I1)
        E2 = EXTDIRP(I2)
        E3 = EXTDIRP(I3)
        E4 = EXTDIRP(I4)
        E5 = EXTDIRP(I1 + 1)
        E6 = EXTDIRP(I2 + 1)
        E7 = EXTDIRP(I3 + 1)
        E8 = EXTDIRP(I4 + 1)

        ! Compute the distance to the next grid plane in  X, Y, and Z
        !   If in horizontal uniform region or doing IP then fix X and/or Y.
        if (ZE .ge. UNIFORMZLEV) then
            CONSTX = .true.
            CONSTY = .true.
        end if
        if (CONSTX) then
            SOX = 1.0e30
        else if (CX .gt. 0.0) then
            SOX = (X1 - XE)*CXINV
            XP = X1
        else
            SOX = (X0 - XE)*CXINV
            XP = X0
        end if
        if (CONSTY) then
            SOY = 1.0e30
        else if (CY .gt. 0.0) then
            SOY = (Y1 - YE)*CYINV
            YP = Y1
        else
            SOY = (Y0 - YE)*CYINV
            YP = Y0
        end if
        if (CZ .gt. 0.0) then
            SOZ = (Z1 - ZE)*CZINV
            ZP = Z1
        else if (CZ .lt. 0.0) then
            SOZ = (Z0 - ZE)*CZINV
            ZP = Z0
        else
            SOZ = 1.0e30
        end if

        ! The shortest distance is the plane we stop at:
        !   get the exitting location and increment the cell
        XOFFS = 0.0
        YOFFS = 0.0
        if (SOZ .le. SOX .and. SOZ .le. SOY) then
            SO = SOZ
            if (.not. CONSTX) XP = XE + SO*CX
            if (.not. CONSTY) YP = YE + SO*CY
            K = K + DK
        else if (SOX .le. SOY) then
            SO = SOX
            if (.not. CONSTY) YP = YE + SO*CY
            ZP = ZE + SO*CZ
            I = I + DI
            ! If have reached a horizontal boundary then either wrap around
            !   (periodic) or go into IP mode (open boundaries).
            if (I .eq. 0) then
                if (btest(BCFLAG, 0)) then
                    I = 1
                    CONSTX = .true.
                else if (btest(BCFLAG, 2)) then
                    HITBOUNDARY = .true.
                    SIDE = 1
                else
                    I = NPX
                    XOFFS = XDOMAIN
                end if
            else if (I .ge. NPX .and. btest(BCFLAG, 2)) then
                HITBOUNDARY = .true.
                SIDE = 2
            else if (I .eq. NPX + 1) then
                if (btest(BCFLAG, 0)) then
                    I = NPX
                    CONSTX = .true.
                else
                    I = 1
                    XOFFS = -XDOMAIN
                end if
            end if
        else
            SO = SOY
            if (.not. CONSTX) XP = XE + SO*CX
            ZP = ZE + SO*CZ
            J = J + DJ
            if (J .eq. 0) then
                if (btest(BCFLAG, 1)) then
                    J = 1
                    CONSTY = .true.
                else if (btest(BCFLAG, 3)) then
                    HITBOUNDARY = .true.
                    SIDE = 3
                else
                    J = NPY
                    YOFFS = YDOMAIN
                end if
            else if (J .ge. NPY .and. btest(BCFLAG, 3)) then
                HITBOUNDARY = .true.
                SIDE = 4
            else if (J .eq. NPY + 1) then
                if (btest(BCFLAG, 1)) then
                    J = NPY
                    CONSTY = .true.
                else
                    J = 1
                    YOFFS = -YDOMAIN
                end if
            end if
        end if
        if (SO .lt. -EPSS) then
            write (6, *) 'DIRECT_BEAM_PROP: SO<0', X, Y, Z,&
                &XE, YE, ZE, XP, YP, ZP, CX, CY, CZ, SOX, SOY, SOZ
            stop
        end if
        SO = max(SO, 0.0d0)
        ! Make the starting and ending interpolation factors
        AX = 1.0d0/(X1 - X0)
        AY = 1.0d0/(Y1 - Y0)
        AZ = 1.0d0/(Z1 - Z0)
        U0 = (XE - X0)*AX
        V0 = (YE - Y0)*AY
        W0 = (ZE - Z0)*AZ
        U1 = (XP - X0)*AX
        V1 = (YP - Y0)*AY
        W1 = (ZP - Z0)*AZ
        ! Compute the cubic polynomial extinction coefficients
        U0M = 1.0 - U0
        V0M = 1.0 - V0
        W0M = 1.0 - W0
        U1M = 1.0 - U1
        V1M = 1.0 - V1
        W1M = 1.0 - W1
        DU = U1 - U0
        DV = V1 - V0
        DW = W1 - W0
        UV = U0*V0
        UMV = U0M*V0
        UVM = U0*V0M
        UMVM = U0M*V0M
        UW = U0*W0
        UMW = U0M*W0
        UWM = U0*W0M
        UMWM = U0M*W0M
        VW = V0*W0
        VMW = V0M*W0
        VWM = V0*W0M
        VMWM = V0M*W0M
        A =    (E1*U0M + E2*U0)*VMWM&
            &+ (E3*U0M + E4*U0)*VWM&
            &+ (E5*U0M + E6*U0)*VMW&
            &+ (E7*U0M + E8*U0)*VW
        B1 = -DU*VMWM - DV*UMWM - DW*UMVM
        B2 = DU*VMWM - DV*UWM - DW*UVM
        B3 = -DU*VWM + DV*UMWM - DW*UMV
        B4 = DU*VWM + DV*UWM - DW*UV
        B5 = -DU*VMW - DV*UMW + DW*UMVM
        B6 = DU*VMW - DV*UW + DW*UVM
        B7 = -DU*VW + DV*UMW + DW*UMV
        B8 = DU*VW + DV*UW + DW*UV
        B = B1*E1 + B2*E2 + B3*E3 + B4*E4 + B5*E5 + B6*E6 + B7*E7 + B8*E8
        VW = DV*DW
        VWU = VW*U0
        VWUM = VW*U0M
        UW = DU*DW
        UWV = UW*V0
        UWVM = UW*V0M
        UV = DU*DV
        UVW = UV*W0
        UVWM = UV*W0M
        C1 = +VWUM + UWVM + UVWM
        C2 = +VWU - UWVM - UVWM
        C3 = -VWUM + UWV - UVWM
        C4 = -VWU - UWV + UVWM
        C5 = -VWUM - UWVM + UVW
        C6 = -VWU + UWVM - UVW
        C7 = +VWUM - UWV - UVW
        C8 = +VWU + UWV + UVW
        C = C1*E1 + C2*E2 + C3*E3 + C4*E4 + C5*E5 + C6*E6 + C7*E7 + C8*E8
        D = DU*DV*DW*(E2 + E3 + E5 + E8 - E1 - E4 - E6 - E7)
        ! Compute the path through the cell: integration of extinction
        PATH = PATH + SO*(A + 0.5d0*B + 0.3333333333333333d0*C + 0.25d0*D)

        XE = XP + XOFFS
        YE = YP + YOFFS
        ZE = ZP
    end do
    DIRFLUX = SOLARFLUX*exp(-PATH)
    DIRPATH = PATH
    XO = XE + XSTART; YO = YE + YSTART; ZO = ZE
    VALIDBEAM = abs(ZE - ZLEVELS(NPZ)) .lt. EPSZ
    if (VALIDBEAM) SIDE = 6
    return
end subroutine
