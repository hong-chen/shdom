!   SHDOM: Spherical harmonic discrete ordinate radiative transfer method.
!      See shdom.txt for documentation.

program shdom_main

    use globals
    use globals_shdom
    use shd_init
    use shd_io_inp
    use shd_io_out
    use shd_ext_mpi
    use shd_dir
    use shd_comp
    use shd_shdo
    use shd_src

    implicit none
    ! Number of output files and inputs parameter per file
    integer, parameter :: MAXOUT = 20, MAXPAR = 500
    real      :: OUTPARMS(MAXPAR, MAXOUT)
    character(len=1) :: OUTTYPES(MAXOUT)
    character(len=FNAMEMAXLEN_) :: OUTFILES(MAXOUT)

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
    logical :: BASEOUT, HIGHORDERRAD, MASTERPROC, VALIDBEAM
    logical :: RADIANCE_OUT, VISUAL_OUT
    logical :: FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
    real :: SOLARFLUX, SOLARMU, SOLARAZ, SOLFLUX
    real :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real :: SOLACC, SOLCRIT, SPLITACC, SHACC
    real :: DELXSFC, DELYSFC, MAXASYM, REDUCE_RATIO
    real :: MUOUT, PHIOUT
    real :: UNIFZLEV, XO, YO, ZO, DIRPATH
    character(len=1) :: SRCTYPE, UNITS
    character(len=1) :: GRIDTYPE, PROPTYPE
    character(len=2) :: SFCTYPE
    character(len=FNAMEMAXLEN_) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=FNAMEMAXLEN_) :: OutFileNC, RUNNAME

    integer :: count_start, count_rate, count_max
    integer :: count_end
    real :: elapsed_time
    real :: DELX0, DELY0
    logical :: mpiX, mpiY, mpiMode
    logical :: nmlMode = .true.

    ! Choose the option to compute higher order radiance SH terms
    HIGHORDERRAD = .false.
    NSHOUT = 4
    if (HIGHORDERRAD) NSHOUT = 5
    ! Choose type of initialization/output for k-distributions:
    !   True means initiatize with and output only the base grid.
    BASEOUT = .false.

    !call cpu_time(CPUTIME1)
    call system_clock(count_start, count_rate, count_max)

    ! Start up MPI for multiple processors if available, and
    !   return flag for master process (procnum=0)
    call shd_ext__mpiStart(MASTERPROC)

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

    call shd_ext__mpiBcastInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
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
        call shd_io__readPropSize(PROPFILE, NSTLEG, NLEG, NPX_, NPY_, NPZ_,&
                                NUMPHASE_, MAXLEG, MAXPGL, DELX_, DELY_)
    end if
    call shd_ext__mpiBcastPropSize(NPX_, NPY_, NPZ_, DELX_, DELY_,&
                                 NUMPHASE_, MAXLEG, MAXPGL)

    ! Set the full domain grid sizes (variables end in T)
    NXT = NX; NYT = NY; NPXT = NPX_; NPYT = NPY_

    XDOMAIN_ = NPXT*DELX_
    YDOMAIN_ = NPYT*DELY_

    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and NX,NY to the subdomain size, and output
    !   the X and Y starting position of this subdomain.  Set bits
    !   2 and 3 in BCFLAG to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    call shd_ext__mpiMap(BCFLAG, NPX_, NPY_, NX, NY, NZ, DELX_, DELY_, PROPFILE,&
        &XSTART_, YSTART_, RUNNAME)

    if (btest(BCFLAG, 2)) then
        mpiX = .true.
    else
        mpiX = .false.
    end if

    if (btest(BCFLAG, 3)) then
        mpiY = .true.
    else
        mpiY = .false.
    end if

    if (mpiX .or. mpiY) then
        mpiMode = .true.
    else
        mpiMode = .false.
    end if

    ! Set up base grid point actual size (NX1xNY1xNZ)
    NX1 = NX + 1
    if (btest(IPFLAG, 0) .or. btest(BCFLAG, 0) .or. mpiX) NX1 = NX
    NY1 = NY + 1
    if (btest(IPFLAG, 1) .or. btest(BCFLAG, 1) .or. mpiY) NY1 = NY
    if (NZ .le. 1) then
        write (6, '(A)') 'Error [shdom_main]: NZ < 2.'
        stop
    end if
    NBPTS = NX1*NY1*NZ
    ! Calculate the number of base grid cells depending on the BCFLAG
    NBCELLS = (NZ - 1)*(NX + ibits(BCFLAG, 0, 1) - ibits(BCFLAG, 2, 1))&
        *(NY + ibits(BCFLAG, 1, 1) - ibits(BCFLAG, 3, 1))

    ! Determine the size of the property grid arrays so can allocate
    MAXPG = NPX_*NPY_*NPZ_
    if (NUMPHASE_ .eq. 0) MAXPGL = MAXPG*(MAXLEG + 1)*NSTLEG

    allocate (ZLEVELS_(NPZ_))
    allocate (TEMPP_(MAXPG), EXTINCTP_(MAXPG), ALBEDOP_(MAXPG))
    allocate (EXTDIRP_(MAXPG), IPHASEP_(MAXPG), LEGENP_(MAXPGL))

    MEMWORD = MEMWORD + 4.5*MAXPG + MAXPGL + NUMPHASE_*NSTLEG*(MAXLEG + 1)

    MAXNZ = NPZ_

    if (mpiMode) then

        ! Read the properties of the medium to the master processor and
        !  distribute to all the processors
        if (MASTERPROC) then
            MAXPGT = NPXT*NPYT*NPZ_
            MAXPGLT = MAXPGL
            if (NUMPHASE_ .eq. 0) MAXPGLT = MAXPGT*(MAXLEG + 1)*NSTLEG

            allocate (TEMPPT(MAXPGT), EXTINCTPT(MAXPGT), ALBEDOPT(MAXPGT))
            allocate (IPHASEPT(MAXPGT), LEGENPT(MAXPGLT))

            ! only use DELX and DELY from initial read, DELX0 and DELY0 are discarded
            call shd_io__readProp(PROPFILE, NPXT, NPYT, NPZ_,&
                &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPGT, MAXPGLT, DELTAM,&
                &PROPTYPE, DELX0, DELY0, ZLEVELS_, MAXASYM,&
                &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE_, IPHASEPT)

            if (NUMPHASE_ .eq. 0) then
                write (6, '(A)') 'Error [shdom_main]: Standard property file (with phase function Legendre series at every point)'
                write (6, '(A)') '  is not allowed with multiple processors.'
                stop
            end if

            ! Check values of medium properties
            call shd_io__checkPropInput(NPXT, NPYT, NPZ_, NSTLEG, NLEG,&
                &DELX_, DELY_, ZLEVELS_,&
                &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE_, IPHASEPT)

        end if

        call shd_ext__mpiScatProp(NPXT, NPYT, NPX_, NPY_, NPZ_,&
            &NSTLEG, NLEG, NUMPHASE_,&
            &ZLEVELS_, MAXASYM,&
            &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, IPHASEP_)

        if (MASTERPROC) deallocate (TEMPPT, EXTINCTPT, ALBEDOPT, IPHASEPT)

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readProp(PROPFILE, NPX_, NPY_, NPZ_,&
            &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,&
            &PROPTYPE, DELX_, DELY_, ZLEVELS_, MAXASYM,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, NUMPHASE_, IPHASEP_)

        ! Check values of medium properties
        call shd_io__checkPropInput(NPX_, NPY_, NPZ_, NSTLEG, NLEG,&
            &DELX_, DELY_, ZLEVELS_,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, NUMPHASE_, IPHASEP_)

    end if

    if (MASTERPROC) then
        ! Check input parameters
        call shd_io__checkInputParam(NMU, NPHI, DELTAM, MAXASYM,&
            &GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX, SOLARMU,&
            &WAVENO, WAVELEN, GNDTEMP, GNDALBEDO,&
            &SPLITACC, SHACC)
    end if

    if (mpiMode) then
        ! Read the memory parameters from a file, if it exists, and send
        !  the parameters to all processors. Otherwise use the input ones.
        call shd_ext__mpiReadBcastMemParam(MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
            &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, RUNNAME)
    end if

    ! Guess maximum number of grid points, cells, SH vector size needed
    ! but don't let MAX_TOTAL_MB be exceeded
    if (MAX_TOTAL_MB*1024.0**2 .gt. 1.80*huge(MAXIV)) then
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
    REDUCE_RATIO = min(1.0, ((MAX_TOTAL_MB*1024**2)/4 - MEMWORD)/WANTMEM)

    ADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR*REDUCE_RATIO
    if (ADAPT_GRID_FACTOR .lt. 1.0) then
        write (6, '(A)') 'Error [shdom_main]: MAX_TOTAL_MB memory limit exceeded with just base grid.'
        stop
    end if

    if (REDUCE_RATIO .lt. 1.0) then
        write (6, '(A,F8.2)') 'Message [shdom_main]: ADAPT_GRID_FACTOR reduced to:', ADAPT_GRID_FACTOR
    end if

    WANTMEM = REDUCE_RATIO*WANTMEM
    if (WANTMEM .gt. huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom_main]: Number of words of memory exceeds max integer size:', WANTMEM
        stop
    end if

    MAXIG = nint(ADAPT_GRID_FACTOR*NBPTS)
    MAXIC = nint(CELL_TO_POINT_RATIO*MAXIG)
    MAXIV = nint(NUM_SH_TERM_FACTOR*NLM*MAXIG)
    MAXIDO = MAXIG*NPHI0MAX
    if (NUMPHASE_ .gt. 0) then
        MAXIGL = NUMPHASE_*NSTLEG*(MAXLEG + 1)
    else
        MAXIGL = MAXIG*NSTLEG*(MAXLEG + 1)
    end if
    if (4.0*(MAXIV + MAXIG)*NSTOKES .gt. huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom_main]: Size of big SH arrays (MAXIV) probably exceeds max integer number of bytes:', MAXIV
        stop
    end if
    if (4.0*8.0*MAXIC .gt. huge(MAXIV)) then
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

    ! print dimension parameters
    !╭────────────────────────────────────────────────╮!
    !print *, 'NSTOKES                :', NSTOKES
    !print *, 'MAXIG (grid)           :', MAXIG
    !print *, 'MAXIC (cell)           :', MAXIC
    !print *, 'MAXIGL (Legendre)      :', MAXIGL
    !print *, 'MAXIDO (grid x nphi)   :', MAXIDO
    !print *, 'MAXIV (harmonic terms) :', MAXIV
    !╰────────────────────────────────────────────────╯!

    MAXNBC = MAXIG*3/NZ
    ! Read in the surface properties
    if (SFCFILE(1:2) .eq. 'NO') then
        SFCTYPE = 'FL'
    else
        if (MASTERPROC) then
            call shd_io__readSurfSize(SFCFILE, MAXSFCPTS, MAXSFCPARS)
        end if
        ! Broadcast the surface array sizes
        call shd_ext__mpiBcastSurfSize(MAXSFCPTS, MAXSFCPARS)
        allocate (SFCPARMS(MAXSFCPARS*MAXSFCPTS))
        if (MASTERPROC) then
            call shd_io__readSurf(SFCFILE, MAXSFCPTS, MAXSFCPARS,&
                &SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,&
                &NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
        end if
        ! Broadcast the surface array data
        call shd_ext__mpiBcastSurfParam(SFCTYPE, NXSFC, NYSFC, NSFCPAR,&
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
            call shd_io__readCkdSize(CKDFILE, WAVENO, NG, NZCKD_)
        end if
        ! Broadcast the k-distribution array sizes
        call shd_ext__mpiBcastKdistSize(NG, NZCKD_)
        allocate (ZCKD_(NZCKD_), GASABS_(NZCKD_), DELG(NG), KABS(NZCKD_, NG))
        if (MASTERPROC) then
            MAXNG = NG; MAXNZ = NZCKD_
            call shd_io__readCkd(MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,&
                &NG, DELG, NZCKD_, ZCKD_, KABS)
        end if
        call shd_ext__mpiBcastKdistParam(SOLFLUX, NG, DELG, NZCKD_, ZCKD_, KABS)
        SOLARFLUX = SOLFLUX*SOLARFLUX*abs(SOLARMU)
    else
        NG = 1
        allocate (DELG(NG))
        DELG(1) = 1.0
        NZCKD_ = 0
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
                if (mpiX) then
                    NXOUT = max(1, nint(XDOMAIN_/OUTPARMS(2, I)))
                else
                    NXOUT = max(1, nint((XDOMAIN_ + 2*OUTPARMS(4, I))/OUTPARMS(2, I)))
                end if
            end if
            if (OUTPARMS(3, I) .le. 0.0) then
                NYOUT = 1
            else
                if (mpiY) then
                    NYOUT = max(1, nint(YDOMAIN_/OUTPARMS(3, I)))
                else
                    NYOUT = max(1, nint((YDOMAIN_ + 2*OUTPARMS(5, I))/OUTPARMS(3, I)))
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
    if (mpiMode .and. MASTERPROC) then
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
        call shd_init__makeGrid(BCFLAG, GRIDTYPE, NPX_, NPY_, NPZ_, NX, NY, NZ,&
            &XSTART_, YSTART_, DELX_, DELY_, ZLEVELS_,&
            &XGRID, YGRID, ZGRID)
    end if

    write (6, '(A,"<",A,">")') 'Message [shdom_main]: Starting solution iterations for SHDOM run: ', trim(adjustl(RUNNAME))
    TOTITER = 0
    OLDNPTS = 0

    ! Loop over the k-distribution g's from high to low absorption
    do IG = NG, 1, -1
        ! If doing k-distribution then get the gas absorption profile.
        if (KDIST) then
            do I = 1, NZCKD_
                GASABS_(I) = KABS(I, IG)
            end do
        end if

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (BASEOUT .or. (IG .eq. NG .and. NEWGRIDFLAG)) then
            call shd_init__makeCell(BCFLAG, IPFLAG,&
                &NX, NY, NZ, NX1, NY1, NPTS, NCELLS,&
                &XGRID, YGRID, ZGRID, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
            NBCELLS = NCELLS
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in SOURCE (SHPTR), RADIANCE (RSHPTR), FLUXES, DIRFLUX.
        call shd_solv__solveRte(NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
            &ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE_,&
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
                    if (mpiMode) then
                        call shd_ext__mpiComputeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &BCFLAG, XDOMAIN_, YDOMAIN_, IPFLAG,&
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
                        call shd_comp__computeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
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
            call shd_comp__sumOutput(NG - IG, DELG(IG), NRAD, WORK2, SUMRADOUT)
        end if

        if (VISUAL_OUT) then
            IVIS = 0
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'V') then
                    if (mpiMode) then
                        call shd_ext__mpiVisualizeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                            &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            &BCFLAG, XDOMAIN_, YDOMAIN_, IPFLAG,&
                            &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                            &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                            &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                            &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2)
                    else
                        call shd_comp__visualizeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                            &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
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
            call shd_comp__sumOutput(NG-IG, DELG(IG), NVIS, WORK2, SUMVISOUT)
        end if

        if (SOURCE_OUT) then
            ! Compute source function throughout grid for this one angle
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'J') then
                    MUOUT = OUTPARMS(2, I)
                    PHIOUT = OUTPARMS(3, I)*acos(-1.0)/180.0
                    call shd_src__computeOneSrc(NSTOKES, ML, MM, NLM,&
                        &NSTLEG, NLEG, NUMPHASE_,&
                        &NPTS, DELTAM, MUOUT, PHIOUT,&
                        &SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE,&
                        &DIRFLUX, SHPTR, SOURCE, SOURCE1OUT)
                end if
            end do
        end if

        if (BASEOUT) then
            ! If outputting only the base grid then throw out the other points
            call shd_init__makeCell(BCFLAG, IPFLAG,&
                &NX, NY, NZ, NX1, NY1,&
                &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
        else
            ! For the output from the grid point arrays (Flux, Net flux, SH)
            !   we may need to interpolate the running k-distribution sum
            !   to the newly created grid points.
            call shd_comp__interpOutput(OLDNPTS, NPTS, FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
                SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
                NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        end if

        if (FLUX_OUT) then
            call shd_comp__sumOutput(NG - IG, DELG(IG), 2*NPTS, FLUXES, SUMFLUXES)
            call shd_comp__sumOutput(NG - IG, DELG(IG), NPTS, DIRFLUX, SUMDIRFLUX)
        end if
        if (FLUXDIV_OUT) then
            ! Compute the net flux divergence at every point
            call shd_comp__computeNetFluxDiv(NSTOKES, NPTS, RSHPTR, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, PLANCK, DIRFLUX, RADIANCE, WORK)
            call shd_comp__sumOutput(NG - IG, DELG(IG), NPTS, WORK, SUMFLUXDIV)
        end if

        if (SH_OUT) then
            call shd_shdo__computeSH(NSHOUT, NSTOKES, NPTS, SRCTYPE, SOLARMU, SOLARAZ,&
                &DIRFLUX, RSHPTR, ML, MM, RADIANCE, WORK)
            call shd_comp__sumOutput(NG - IG, DELG(IG), NSHOUT*NPTS, WORK, SUMSHTERMS)
        end if

        OLDNPTS = NPTS
    end do

    call shd_io__saveState(OUTSAVEFILE, NSTOKES, NX, NY, NZ, ML, MM, NLM,&
        &WORK, XGRID, YGRID, ZGRID,&
        &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
        &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)

    call system_clock(count_end, count_rate, count_max)
    elapsed_time = (real(count_end) - real(count_start)) / real(count_rate)
    if (MASTERPROC) then
        if (elapsed_time .lt. 60.0) then
            write (6, '(A,F6.2,A)') 'Message [shdom_main]: Total time used:', elapsed_time, ' sec'
        else
            write (6, '(A,F6.2,A)') 'Message [shdom_main]: Total time used:', elapsed_time/60.0, ' min'
        endif
    end if

    NSH = SHPTR(NPTS + 1)
    IRAD = 0
    IVIS = 0
    if (mpiMode) then

        if (MASTERPROC) then
            if (FLUX_OUT) allocate (ALLFLUXES(3, NZ, NYT, NXT))
            if (FLUXDIV_OUT) allocate (ALLFLUXDIV(NZ, NYT, NXT))
            if (SH_OUT) allocate (ALLSHTERMS(NSHOUT, NZ, NYT, NXT))
            allocate (XALLGRID(NXT), YALLGRID(NYT))
        end if

        call shd_ext__mpiGatherOutput(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
            &IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,&
            &SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
            &XDOMAIN_, YDOMAIN_, XALLGRID, YALLGRID,&
            &ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,&
            &NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)

        if (MASTERPROC) then
            do I = 1, NUMOUT
                if (OUTTYPES(I) .ne. 'V' .and. OUTFILES(I) .ne. 'NONE') then
                    call shd_io__outputResultMpi(NSTOKES, NXT, NYT, NZ,&
                        &NPTSTOT, NCELLSTOT, NSHTOT,&
                        &ML, MM, NLM, NMU, NPHI, NANG, NG,&
                        &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                        &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                        &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                        &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                        &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                        &XDOMAIN_, YDOMAIN_, XALLGRID, YALLGRID, ZGRID,&
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
            if (OUTTYPES(I) .ne. 'V' .and. OUTFILES(I) .ne. 'NONE') then
                call shd_io__outputResult(NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NCELLS,&
                    NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_, NMU, NPHI, NANG, NG,&
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
            if (OUTTYPES(I) .eq. 'V') then
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
        &+ 4.5*MAXPG + MAXPGL + NUMPHASE_*(MAXLEG + 1)&
        &+ 16.5*NCELLS + NPTS*(28 + NPHI0MAX*NSTOKES)&
        &+ NSTOKES*NSH*BIG_ARRAYS)/1024**2
    ADAPTGRIDFAC_OUT = float(NPTS)/NBPTS
    SHTERM_FAC_OUT = float(NSH)/(NLM*NPTS)
    CELL_POINT_OUT = float(NCELLS)/NPTS
    write (6, '(A)')'Message [shdom_main]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    write (6, '(28X,F7.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    call shd_ext__mpiEnd(NPTS, GRIDPOS, NPX_, NPY_, XSTART_, YSTART_, XDOMAIN_, YDOMAIN_,&
        &NPXT, NPYT, PROPFILE, RUNNAME)

    if (SRCTYPE .ne. 'T') then
        call shd_dir__computeDirect(9, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
            &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end if

    deallocate (SUMRADOUT, SUMVISOUT, WORK2)
    if (FLUX_OUT) deallocate (SUMFLUXES, SUMDIRFLUX)
    if (FLUXDIV_OUT) deallocate (SUMFLUXDIV)
    if (SH_OUT) deallocate (SUMSHTERMS)
    if (SOURCE_OUT) deallocate (SOURCE1OUT)
    if (KDIST) deallocate (ZCKD_, GASABS_, KABS)
    deallocate (BCPTR, BCRAD, DELG)
    if (allocated(SFCPARMS)) deallocate (SFCPARMS, SFCGRIDPARMS)
    deallocate (RADIANCE, SOURCE)
    if (ACCELFLAG) deallocate (DELSOURCE)
    deallocate (RSHPTR, SHPTR, OSHPTR, WORK, WORK1, FLUXES, DIRFLUX)
    deallocate (GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS)
    deallocate (TEMP, PLANCK, EXTINCT, ALBEDO, IPHASE, LEGEN)
    deallocate (XGRID, YGRID, ZGRID)
    deallocate (ZLEVELS_, TEMPP_, EXTINCTP_, ALBEDOP_, EXTDIRP_, IPHASEP_, LEGENP_)
    deallocate (MU, WTDO, PHI, NPHI0)

end program shdom_main
