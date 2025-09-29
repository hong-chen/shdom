! SHDOMPP: Spherical harmonic discrete ordinate radiative transfer method
!          for plane-parallel atmospheres.
!     See shdompp.doc for documentation.

program SHDOMPP

    use globals
    use globals_prp
    use shd_inp ! shd_io__
    use shd_out ! shd_io__

    implicit none

    !  MAXDELTAU  maximum scattering optical depth allowed in initial sublayers
    !  SPLITTING_FACTOR  max number of internal grid points is this factor
    !                    times the number of input layers + optical depth

    real, parameter :: MAXDELTAU = 1.0
    real, parameter :: SPLITTING_FACTOR = 10.0
    logical, parameter :: PRNT = .true. ! True for iteration printing
    !LOGICAL, PARAMETER :: PRNT=.FALSE.

    integer, parameter :: MAXPAR = 100 ! Max array size for user input arrays
    integer, parameter :: MAXOUT = 20
    real :: OUTPARMS(MAXPAR, MAXOUT)
    character(len=1) :: OUTTYPES(MAXOUT)
    character(len=FNAMEMAXLEN_) :: OUTFILES(MAXOUT)

    real :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
    real, allocatable :: SFCPARMS(:), SFCPARMS_TMP(:)

    integer :: NLAY
    integer :: NUMITER
    integer :: NZOUT, NMUOUT, NPHIOUT, FLUXPARM
    integer :: NCASES, ICASE
    real    :: SFCTEMP

    real, allocatable :: SOLFLUX(:)
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
    real, allocatable :: FLUXES(:, :), DIRFLUX(:), DIRFLUX_TMP(:)
    integer, allocatable :: NPHI0(:)
    real, allocatable :: MU(:), WTDO(:), PHI(:)
    integer, allocatable :: BCPTR(:)
    real, allocatable :: BCRAD(:, :)
    real, allocatable :: SFCGRIDPARMS(:)
    real, allocatable :: SENPARMS(:)
    real, allocatable :: SUMFLUXES(:, :), SUMDIRFLUX(:), SUMFLUXDIV(:)
    real, allocatable :: SUMSHTERMS(:), SOURCE1OUT(:)
    real, allocatable :: SUMRADOUT(:), SUMVISOUT(:)
    real, allocatable :: EXTINCTP(:)
    real, allocatable :: ALLFLUXES(:, :, :, :)
    real, allocatable :: ALLFLUXDIV(:, :, :), ALLSHTERMS(:, :, :, :)
    real, allocatable :: XALLGRID(:), YALLGRID(:)
    real, allocatable :: TEMPPT(:), EXTINCTPT(:), ALBEDOPT(:), LEGENPT(:)
    integer*2, allocatable :: IPHASEPT(:)

    integer :: MAXPG, MAXIG, MAXIC, MAXIV, MAXIDO, MAXNZ, MAXNG
    integer :: MAXLEG, MAXPGL, MAXIGL
    integer :: MAXNBC, MAXBCRAD
    integer :: MAXSFCPTS, MAXSFCPARS, BIG_ARRAYS
    real :: MEMWORD, WANTMEM
    real :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer :: NXSFC, NYSFC, NSFCPAR, NTOPPTS, NBOTPTS
    integer :: NG, SIDE, IG, ITER, I
    integer :: ML, MM, NLM, NSTLEG, NLEG, NANG, NPHI0MAX, NX1, NY1
    integer :: NBPTS, NPTS, NCELLS, NSH, NBCELLS, OLDNPTS
    integer :: NRAD, IRAD, NXOUT, NYOUT, NVIS, IVIS, NSCAN, NPIX, NSHOUT
    integer :: NCELLSTOT, NPTSTOT, NSHTOT
    logical :: INRADFLAG, NEWGRIDFLAG
    logical :: BASEOUT, HIGHORDERRAD, MASTERPROC, VALIDBEAM
    logical :: RADIANCE_OUT, VISUAL_OUT
    logical :: FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
    real :: DELXSFC, DELYSFC, MAXASYM, REDUCE_RATIO
    real :: UNIFZLEV, XO, YO, ZO, DIRPATH
    character(len=2) :: SFCTYPE
    integer :: NSENDATA, NSENPAR

    character(len=FNAMEMAXLEN_) :: RUNNAME, PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, SENFILE, OutFileNC
    integer :: NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    logical :: KDIST, DELTAM, ACCELFLAG
    character(len=1) :: GRIDTYPE, SRCTYPE, UNITS
    character(len=1) :: PROPTYPE
    real :: SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD
    real :: GNDTEMP, GNDALBEDO, WAVENO(2), WAVELEN
    real :: SOLACC, SOLCRIT, SPLITACC, SHACC
    integer :: MAXITER, NUMOUT
    real :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

    real :: SOLARFLUX

    integer :: count_start, count_rate, count_max
    integer :: count_end
    real :: elapsed_time
    real :: DELX0, DELY0
    logical :: mpiX, mpiY, mpiMode
    logical :: nmlMode = .true.

    integer, allocatable :: NLEGP(:), IXP(:, :)
    real, allocatable :: TAUP(:), ALBEDOP(:), LEGENP(:, :)
    real, allocatable :: TAUPM(:), ALBEDOPM(:)
    real, allocatable :: HEIGHTP(:), TEMPP(:)
    real, allocatable :: DELG(:), ZCKD(:), KABS(:, :)
    real, allocatable :: HEIGHTG(:), TAUG(:)
    real, allocatable :: FLUXUP(:), FLUXDN(:), MEANRAD(:), RADOUT(:, :, :)
    real, allocatable :: FLUXUPM(:), FLUXDNM(:), MEANRADM(:), RADOUTM(:, :, :)

    ! namelist input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call shd_io__nmlInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
        INSAVEFILE, OUTSAVEFILE, SENFILE,&
        NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
        IPFLAG, KDIST, DELTAM, GRIDTYPE,&
        SRCTYPE, SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD,&
        GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
        ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
        MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
        OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! process input parameters
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! set calculation accuracy (1.0e-4 to 1.0e-8)
    RACC_ = max(REPS_, min(SOLACC*0.1, 1.0e-4))

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    NSTOKES = max(1, min(4, NSTOKES))
    if (NSTOKES .eq. 2) NSTOKES = 3

    ! X, Y, Z grids (force NX&NY to 1 and IPFLAG to 3)
    NX = 1; NY = 1; NZ = max(2, NZ)
    IPFLAG = 3

    ! Make ML and MM from NMU and NPHI
    !   NMU : number of discrete ordinate zenith angles in both hemispheres
    !   NPHI: number of discrete ordinate azimuth angles in 2*pi
    !   ML  : The maximum L of the spherical harmonic truncation
    !   MM  : The maximum M of the spherical harmonic truncation
    NMU  = max(2, 2*int((NMU + 1)/2)) ! force NMU an even number
    NPHI = max(1, NPHI)

    ML = NMU - 1
    MM = max(0, int(NPHI/2) - 1)

    ! NLEG: The number of Legendre series phase function terms used in the solution procedure
    NLEG = ML
    if (DELTAM) NLEG = ML + 1

    ! NSTLEG is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (NSTOKES .eq. 1) then
        NSTLEG = 1
    else
        NSTLEG = 6
    end if

    NLEG = ML
    if (DELTAM) NLEG = ML + 1
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! read PROPFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call shd_io__readPrpSize(PROPFILE, NSTLEG, NLEG, NPX_, NPY_, NPZ_,&
        &NUMPHASE_, MAXLEG, MAXPGL, DELX_, DELY_)

    MAXPG = NPZ_
    if (NUMPHASE_ .eq. 0) MAXPGL = MAXPG*(MAXLEG + 1)*NSTLEG
    MAXNZ = NPZ_

    allocate (ZLEVELS_(NPZ_))
    allocate (TEMPPT(MAXPG), EXTINCTPT(MAXPG), ALBEDOPT(MAXPG))
    allocate (IPHASEPT(MAXPG), LEGENPT(MAXPGL))

    call shd_io__readPrp(PROPFILE, NPX_, NPY_, NPZ_,&
        &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,&
        &PROPTYPE, DELX_, DELY_, ZLEVELS_, MAXASYM,&
        &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE_, IPHASEPT)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! process layer property
    !╭────────────────────────────────────────────────────────────────────────────╮!
    NLAY = NPZ_ - 1
    allocate (HEIGHTP(NLAY+1), TEMPP(NLAY+1))
    HEIGHTP(1:NPZ_) = ZLEVELS_(1:NPZ_)
    TEMPP(1:NPZ_)   = TEMPPT(1:NPZ_)

    allocate (TAUP(NLAY), ALBEDOP(NLAY))
    TAUP(1:NLAY) = EXTINCTPT(1:NLAY)*(HEIGHTP(2:NPZ_)-HEIGHTP(1:NLAY))
    ALBEDOP(1:NLAY) = ALBEDOPT(1:NLAY)

    allocate (NLEGP(NLAY), LEGENP(0:MAXLEG, NLAY))
    !╰────────────────────────────────────────────────────────────────────────────╯!
    allocate (TAUPM(NLAY), ALBEDOPM(NLAY), IXP(2, NLAY))

    !call USER_INPUT(PROPFILE, CKDFILE, KDIST, NMU, NPHI, &
    !                SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
    !                SFCTEMP, SFCTYPE, SFCPARMS, WAVELEN, WAVENO, &
    !                ACCELFLAG, SOLACC, MAXITER, SPLITACC, &
    !                FLUXFILE, FLUXPARM, RADFILE, MAXPAR, NZOUT, ZOUT, &
    !                NMUOUT, MUOUT, NPHIOUT, PHIOUT)

    ! Find out the size of the property grid arrays for allocating arrays
    !call READ_PROPERTY_SIZE(PROPFILE, NLAY, MAXLEG)

    ! Read in the properties of the medium
    !call READ_PROPERTIES(PROPFILE, NLAY, MAXLEG, &
    !                     HEIGHTP, TEMPP, TAUP, ALBEDOP, NLEGP, LEGENP)

    !allocate (RADOUT(NZOUT, NMUOUT, NPHIOUT), RADOUTM(NZOUT, NMUOUT, NPHIOUT))

    ! Guess the maximum size of the adaptive grid and allocate arrays
    MAXIG = nint(SPLITTING_FACTOR*(NLAY + sum(TAUP(:)*ALBEDOP(:))/MAXDELTAU))
    allocate (HEIGHTG(MAXIG), TAUG(MAXIG))
    allocate (FLUXUP(MAXIG), FLUXDN(MAXIG), MEANRAD(MAXIG))
    allocate (FLUXUPM(MAXIG), FLUXDNM(MAXIG), MEANRADM(MAXIG))

    ! read SFCFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (SFCFILE(1:4) .eq. 'NONE') then
        SFCTYPE = 'FL'
    else
        call shd_io__readSfcSize(SFCFILE, MAXSFCPTS, MAXSFCPARS)
        allocate (SFCPARMS_TMP(MAXSFCPARS*MAXSFCPTS), SFCPARMS(MAXSFCPARS))
        call shd_io__readSfc(SFCFILE, MAXSFCPTS, MAXSFCPARS,&
            &SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,&
            &NSFCPAR, SFCPARMS_TMP, GNDTEMP, GNDALBEDO)
        SFCPARMS(1:MAXSFCPARS) = SFCPARMS_TMP(1:MAXSFCPARS)
        deallocate (SFCPARMS_TMP)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! read CKDFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (KDIST) then
        call shd_io__readCkdSize(CKDFILE, WAVENO, NG, NZCKD_)

        allocate (ZCKD_(NZCKD_), GASABS_(NZCKD_), SOLFLUX(NG), DELG(NG), KABS(NZCKD_, NG))

        MAXNG = NG
        MAXNZ = NZCKD_

        call shd_io__readCkd(MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,&
            &NG, DELG, NZCKD_, ZCKD_, KABS)
    else
        NG = 1
        allocate (SOLFLUX(NG), DELG(NG))
        SOLFLUX(1) = 1.0
        DELG(1) = 1.0
        NZCKD_ = 0
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! read SENFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (SENFILE(1:4) .ne. 'NONE') then
        call shd_io__readSenSize(SENFILE, NSENDATA, NSENPAR)
        allocate (SENPARMS(NSENDATA*NSENPAR))
        call shd_io__readSen(SENFILE, NSENDATA, NSENPAR, SENPARMS)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!

        !call READ_CKD_SIZE(CKDFILE, WAVENO, NG, NZCKD)
        !call READ_CKD(CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)

    !WRITE (*,*) 'Number of times to run SHDOMPP routine for timing'
    !READ (*,*) NCASES
    !WRITE (*,*) NCASES
    !DO ICASE = 1, NCASES

    ! Loop over the k-distribution g's from low to high absorption
    NUMITER = 0
    do IG = 1, NG
        SOLARFLUX = SOLFLUX(IG)*SOLARFLUX0*abs(SOLARMU)
        write (6, '(A,I2,A,I2,A)') 'Message [shdompp]: Starting ', IG, '/', NG, ' (ig/Ng) ...'
    !    ! Add in the layer gaseous absorption if doing a k-distribution
    !    call ADD_IN_KDIST(NLAY, HEIGHTP, TAUP, ALBEDOP, &
    !                      NZCKD, ZCKD, KABS(:, IG), TAUPM, ALBEDOPM)

    !    call SOLVE_SHDOMPP(NLAY, MAXLEG, HEIGHTP, TEMPP, TAUPM, ALBEDOPM, &
    !                       NLEGP, LEGENP, SRCTYPE, SOLARFLUX, SOLARMU, &
    !                       WAVELEN, WAVENO, SKYRAD, SFCTYPE, SFCPARMS, SFCTEMP, &
    !                       NMU, NPHI, MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
    !                       MAXITER, ACCELFLAG, PRNT, SOLCRIT, ITER, &
    !                       NPTS, HEIGHTG, TAUG, IXP, FLUXUPM, FLUXDNM, MEANRADM, &
    !                       NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
    !                       RADOUTM)

        NUMITER = NUMITER + ITER
    !    if (KDIST) then
    !        call SUM_KDIST_OUTPUT(IG, DELG(IG), NLAY, TAUP, IXP, NPTS, TAUG, &
    !                              HEIGHTP, HEIGHTG, FLUXUPM, FLUXUP, &
    !                              FLUXDNM, FLUXDN, MEANRADM, MEANRAD, &
    !                              NZOUT, NMUOUT, NPHIOUT, RADOUTM, RADOUT)
    !    else
    !        FLUXUP(1:NPTS) = FLUXUPM(1:NPTS)
    !        FLUXDN(1:NPTS) = FLUXDNM(1:NPTS)
    !        MEANRAD(1:NPTS) = MEANRADM(1:NPTS)
    !        RADOUT(:, :, :) = RADOUTM(:, :, :)
    !    end if
    end do
    print *, 'hahaha'
    stop


    !end do

    !call OUTPUT_RESULTS(NLAY, NPTS, NMU, NPHI, NG, &
    !                    PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
    !                    SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, &
    !                    SPLITACC, SOLACC, MAXITER, NUMITER, &
    !                    IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
    !                    NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
    !                    'F', FLUXFILE)

    !if (NZOUT*NMUOUT*NPHIOUT > 0) then
    !    call OUTPUT_RESULTS(NLAY, NPTS, NMU, NPHI, NG, &
    !                        PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
    !                        SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, &
    !                        SPLITACC, SOLACC, MAXITER, NUMITER, &
    !                        IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
    !                        NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
    !                        'R', RADFILE)
    !end if

    deallocate (HEIGHTP, TEMPP, TAUP, ALBEDOP, TAUPM, ALBEDOPM, NLEGP, LEGENP)
    deallocate (IXP, HEIGHTG, TAUG, FLUXUP, FLUXDN, MEANRAD, RADOUT)
    deallocate (FLUXUPM, FLUXDNM, MEANRADM, RADOUTM)
    deallocate (DELG, ZCKD, KABS)

    if (allocated(SFCPARMS)) deallocate (SFCPARMS)
    if (allocated(SENPARMS)) deallocate (SENPARMS)
    deallocate (TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT)
    deallocate (ZLEVELS_, TEMPP_, EXTINCTP_, ALBEDOP_, EXTDIRP_, IPHASEP_, LEGENP_)

end program SHDOMPP

subroutine SOLVE_SHDOMPP(NLAY, MAXLEG, HEIGHTP, TEMPP, TAUP, ALBEDOP, &
                         NLEGP, LEGENP, SRCTYPE, SOLARFLUX, SOLARMU, &
                         WAVELEN, WAVENO, SKYRAD, SFCTYPE, SFCPARMS, SFCTEMP, &
                         NMU, NPHI, MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
                         MAXITER, ACCELFLAG, PRNT, SOLCRIT, ITER, &
                         NPTS, HEIGHTG, TAUG, IXP, FLUXUP, FLUXDN, MEANRAD, &
                         NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                         RADOUT)
    ! Performs the SHDOMPP solution procedure to calculate plane-parallel
    ! monochromatic, unpolarized radiative transfer.  The optical properties,
    ! uniform in each layer, are input.  The output are the upwelling and
    ! downwelling fluxes and mean radiance at the adaptive grid optical depth
    ! levels (TAUG), and the radiances in the NMUOUT*NPHIOUT specified
    ! directions at each of NZOUT radiance output levels.
    !
    ! Inputs:
    !  NLAY      number of input property layers (TAUP, ALBEDOP, etc.)
    !  MAXLEG    maximum Legendre series order for phase functions
    !  HEIGHTP   heights (e.g. km) of input layer interfaces (NLAY+1 values)
    !  TEMPP     temperatures (K) of input layer interfaces (NLAY+1 values)
    !              temperature is linearly interpolated in tau within each layer
    !  TAUP      optical depths of the input layers (each layer from top down)
    !  ALBEDOP   single scattering albedo of the inputs layers
    !  NLEGP     order of each layer Legendre series phase function
    !  LEGENP    phase function Legendre series coefficients for each layer
    !  SRCTYPE   'S' for solar source, 'T' for thermal source, 'B' for both
    !  SOLARFLUX solar flux on a *horizontal* surface
    !  SOLARMU   negative of cosine of solar zenith angle (< 0)
    !              direct beam is assumed to be going in phi=0 azimuth
    !  WAVELEN   wavelength (microns) for thermal emission
    !  WAVENO(2) starting and ending wavenumber (cm^-1) for thermal (if WAVELEN=0)
    !  SKYRAD    top of atmosphere downward sky radiance (or temp for SRCTYPE='T')
    !  SFCTYPE   surface type (currently 'L' for Lambertian, 'R' for RPV
    !  SFCPARMS  array of surface parameters (first is albedo for Lambertian)
    !  SFCTEMP   surface temperature (K)
    !  NMU       number of discrete ordinate zenith angles in both hemispheres
    !  NPHI      number of discrete ordinate azimuth angles in 2\pi
    !  MAXIG     maximum possible size of adaptive grid arrays
    !  MAXDELTAU max sublayer scattering optical depth in adaptive grid
    !  SPLITACC  adaptive layer splitting accuracy
    !  SOLACC    desired solution accuracy (try 1E-4 or less)
    !  MAXITER   maximum number of iterations allowed
    !  ACCELFLAG logical flag: true to perform iteration acceleration
    !  PRNT      logical flat: true to printer iteration information
    !  NZOUT     number of radiance output height levels
    !  ZOUT      height levels of radiance output (set to closest adaptive level)
    !  NMUOUT    number of cosine zenith angles for radiance output
    !  MUOUT     cosine zenith angles for radiance output
    !  NPHIOUT   number of azimuth angles for radiance output
    !  PHIOUT    azimuth angles (degrees) for radiance output
    !
    ! Outputs:
    !  SOLCRIT   actual solution accuracy achieved
    !  ITER      actual number of iterations
    !  NPTS      number of adaptive grid points
    !  HEIGHTG   heights of the adaptive grid levels (same units as HEIGHTP)
    !  TAUG      sorted grid point optical depth array for the adaptive grid
    !  IXP       index of starting and ending point in TAUG for each input layer
    !  FLUXUP    upward hemispheric flux at each adaptive grid level
    !  FLUXDN    downward hemispheric fluxes (includes direct solar) at each level
    !  MEANRAD   mean radiance (includes direct solar) at each level
    !  RADOUT    array of output radiances
    !
    !  TAUG(IXP(1,layer)) is the optical depth at the top of the input layer
    !  TAUG(IXP(2,layer)) is the optical depth at the bottom of the input layer
    !  FLUXUP(IXP(1,layer)) is the upwelling flux at the top of the layer
    !  FLUXDN(IXP(2,layer)) is the downwelling flux at the bottom of the layer
    !
    implicit none
    integer, intent(IN)  :: NLAY, MAXLEG, NLEGP(NLAY), MAXIG
    integer, intent(IN)  :: NMU, NPHI, MAXITER
    integer, intent(IN)  :: NZOUT, NMUOUT, NPHIOUT
    logical, intent(IN)  :: ACCELFLAG, PRNT
    real, intent(IN)  :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
    real, intent(IN)  :: HEIGHTP(NLAY + 1), TEMPP(NLAY + 1)
    real, intent(IN)  :: SOLARFLUX, SOLARMU, WAVELEN, WAVENO(2), SKYRAD
    real, intent(IN)  :: SFCPARMS(*), SFCTEMP
    real, intent(IN)  :: MAXDELTAU, SPLITACC, SOLACC
    real, intent(INOUT)  :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
    character(LEN=1), intent(IN) :: SRCTYPE, SFCTYPE
    integer, intent(OUT) :: ITER, NPTS, IXP(2, NLAY)
    real, intent(OUT) :: SOLCRIT
    real, intent(OUT) :: HEIGHTG(MAXIG), TAUG(MAXIG)
    real, intent(OUT) :: FLUXUP(MAXIG), FLUXDN(MAXIG), MEANRAD(MAXIG)
    real, intent(OUT) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)

    integer :: ML, MM, NLM, NPHI0MAX
    integer :: ORDINATESET, OLDNPTS, PREVNPTS, LASTNPTS, NANG, I
    logical, parameter :: DELTAM = .true.
    logical :: FFTFLAG(NMU)
    real, parameter :: MINDELTAUSPLIT = 1.0e-3
    real    :: DELJDOT, DELJOLD, DELJNEW
    real    :: SOLARAZ, PI, GNDALB, DIRFLUX
    real    :: TRUESOLCRIT, CURSPLITACC
    integer, allocatable :: IXG(:), NPHI0(:)
    real, allocatable :: TAUPSC(:), ALBEDOSC(:), LEGENSC(:, :)
    real, allocatable :: TAUGSC(:), RADIANCE(:, :, :), SOURCE(:, :)
    real, allocatable :: MU(:), WTDO(:, :), PHI(:, :), WTMU(:), YLMSUN(:)
    real, allocatable :: CMU1(:, :), CMU2(:, :), CPHI1(:, :, :), CPHI2(:, :, :)
    real, allocatable :: WPHISAVE(:, :), SRCINT(:, :, :)
    real, allocatable :: CONSTSOURCE(:, :), SHRADIANCE(:, :), DOSOURCE(:, :, :)
    real, allocatable :: DELSOURCE(:, :)

    ! Make NPHI0MAX, ML, MM, and NLM from NMU and NPHI
    if (NMU /= max(2, 2*int((NMU + 1)/2))) stop 'SOLVE_SHDOMPP: Bad NMU'
    if (NPHI < 1) stop 'SOLVE_SHDOMPP: NPHI must be greater than 0'
    NPHI0MAX = int((NPHI + 2)/2)
    ML = NMU - 1
    MM = max(0, int(NPHI/2) - 1)
    NLM = (MM + 1)*(ML + 1) - (MM*(MM + 1))/2

    allocate (TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY))
    allocate (TAUGSC(MAXIG), IXG(MAXIG))
    allocate (MU(NMU), PHI(NPHI0MAX, NMU), NPHI0(NMU), WTDO(NMU, NPHI0MAX))
    allocate (YLMSUN(NLM), WTMU(NMU), CMU1(NLM, NMU), CMU2(NMU, NLM))
    allocate (CPHI1(0:16, 32, NMU), CPHI2(32, 0:16, NMU), WPHISAVE(3*NPHI0MAX + 15, NMU))
    allocate (SOURCE(NLM, MAXIG), RADIANCE(NPHI0MAX, NMU, MAXIG))
    allocate (CONSTSOURCE(NLM, MAXIG), DELSOURCE(NLM, MAXIG))
    allocate (SRCINT(0:4, NMU, MAXIG))

    ! Set up some things before solution loop

    ! If Delta-M then scale the optical depth, albedo, and Legendre terms.
    if (DELTAM) then
        call DELTA_SCALE(NLAY, ML, TAUP, ALBEDOP, MAXLEG, NLEGP, LEGENP, &
                         TAUPSC, ALBEDOSC, LEGENSC)
    else
        TAUPSC(:) = TAUP(:)
        ALBEDOSC(:) = ALBEDOP(:)
        LEGENSC(0:ML, :) = LEGENP(0:ML, :)
    end if

    call MAKE_TAU_GRID(NLAY, TAUPSC, ALBEDOSC, MAXDELTAU, MINDELTAUSPLIT, &
                       NPTS, TAUGSC, IXP, IXG)

    ! Precompute Ylm's for solar direction
    if (SRCTYPE /= 'T') then
        SOLARAZ = 0.0
        call YLMALL(SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
    end if
    ! Prepare the constant part of the level source function, which are
    !   the Planck blackbody source and the solar attenuated single scattered
    !   radiance, as a function of level.
    CONSTSOURCE(1, 1:NPTS) = -1.0
    call PREPARE_SOURCE(NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                        SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                        SOLARFLUX, SOLARMU, YLMSUN, WAVELEN, WAVENO, &
                        CONSTSOURCE)

    ! Make the discrete ordinates (angles)
    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
    ORDINATESET = 3
    call MAKE_ANGLE_SET(NMU, NPHI, 1, NPHI0MAX, ORDINATESET, &
                        NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the Ylm transform coefficients
    call MAKE_SH_DO_COEF(ML, MM, NLM, NMU, NPHI0, &
                         NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
                         CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

    ! Initialize the radiance on the base grid using Eddington
    !   two-stream plane-parallel
    if (SFCTYPE == 'L') then
        GNDALB = SFCPARMS(1)
    else
        GNDALB = 0.2
    end if
    allocate (SHRADIANCE(NLM, NPTS))
    call INIT_RADIANCE(NLAY, ML, TAUPSC, ALBEDOSC, LEGENSC, TEMPP, &
                       SRCTYPE, SOLARFLUX, SOLARMU, GNDALB, SFCTEMP, &
                       SKYRAD, WAVELEN, WAVENO, NPTS, TAUGSC, IXP, NLM, &
                       SHRADIANCE)

    ! Initialize the SH source function from the SH radiance field
    OLDNPTS = 0
    SOURCE(:, :) = 0.0
    call COMPUTE_SOURCE(ML, MM, NLM, NLAY, NPTS, IXP, IXG, &
                        ALBEDOSC, LEGENSC, &
                        SHRADIANCE, CONSTSOURCE, SOURCE, &
                        OLDNPTS, DELSOURCE, &
                        DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
    deallocate (SHRADIANCE)

    SOLCRIT = 1.0
    TRUESOLCRIT = 1.0
    CURSPLITACC = 1.0
    LASTNPTS = -1
    OLDNPTS = NPTS
    RADIANCE(:, :, :) = 0.0
    DELSOURCE(:, :) = 0.0
    ITER = 0
    if (MAXITER <= ITER) return

    if (PRNT) write (*, *) '! Iter Log(Sol)  Log(True)  Npoints'

    ! Main solution loop
    ! Iterate until the iterations exceed the limit, or the desired
    ! solution criterion is reached
    do while (ITER < MAXITER .and. (SOLCRIT > SOLACC .or. NPTS /= OLDNPTS))
        ITER = ITER + 1

        ! Transform the source from spherical harmonics to discrete ordinates
        allocate (DOSOURCE(NPHI0MAX, NMU, MAXIG))
        call SH_TO_DO(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                      FFTFLAG, CMU1, CPHI1, WPHISAVE, SOURCE, DOSOURCE)
        ! Enforce the discrete ordinate source function to be non-negative
        DOSOURCE(:, :, :) = max(0.0, DOSOURCE(:, :, :))

        ! Keep on doing integrations as long as new levels are added
        OLDNPTS = NPTS
        PREVNPTS = 0
        do while (NPTS > PREVNPTS)
            ! Integrate the source function along discrete ordinates to
            ! compute the radiance field.
            call PATH_INTEGRATION(NLAY, NPTS, IXP, IXG, TAUGSC, &
                                  NMU, NPHI0MAX, NPHI0, &
                                  SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                                  SFCTYPE, SFCPARMS, SFCTEMP, WAVELEN, WAVENO, &
                                  MU, PHI, WTDO, SRCINT, LASTNPTS, &
                                  FLUXUP, FLUXDN, MEANRAD, DOSOURCE, RADIANCE)
            PREVNPTS = NPTS
            ! Do the adaptive grid cell splitting stuff if desired
            call SPLIT_CELLS(NLAY, NPTS, MAXIG, IXP, IXG, TAUGSC, &
                             NMU, NPHI0MAX, NPHI0, MU, DOSOURCE, RADIANCE, &
                             SPLITACC, TRUESOLCRIT, CURSPLITACC)
            if (NPTS > PREVNPTS) then
                CONSTSOURCE(1, PREVNPTS + 1:NPTS) = -1.0
                call PREPARE_SOURCE(NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                                    SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                                    SOLARFLUX, SOLARMU, YLMSUN, WAVELEN, WAVENO, &
                                    CONSTSOURCE)
            end if
        end do
        deallocate (DOSOURCE)

        ! Transform the radiance from discrete ordinates to spherical harmonics
        allocate (SHRADIANCE(NLM, NPTS))
        call DO_TO_SH(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                      FFTFLAG, CMU2, CPHI2, WPHISAVE, RADIANCE, SHRADIANCE)

        ! Calculate the source function from the radiance in spherical harmonics.
        !   Also computes the solution criterion and the series acceleration
        !   stuff (DELJDOT, DELJOLD, DELJNEW, DELSOURCE).
        call COMPUTE_SOURCE(ML, MM, NLM, NLAY, NPTS, IXP, IXG, &
                            ALBEDOSC, LEGENSC, &
                            SHRADIANCE, CONSTSOURCE, SOURCE, &
                            OLDNPTS, DELSOURCE, &
                            DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
        deallocate (SHRADIANCE)

        ! Accelerate the convergence of the source function vector
        if (ACCELFLAG) then
            call ACCELERATE_SOLUTION(NLM, NPTS, OLDNPTS, SOURCE, DELSOURCE, &
                                     DELJDOT, DELJOLD, DELJNEW)
        end if
        if (DELJNEW < DELJOLD) &
            TRUESOLCRIT = SOLCRIT/(1.0 - sqrt(DELJNEW/DELJOLD))

        ! Print out the iteration, log solution criterion, and number of points
        if (PRNT) write (*, '(2X,I4,2F8.3,1X,I6)') ITER, &
            log10(max(SOLCRIT, 1.0e-20)), log10(max(TRUESOLCRIT, 1.0e-20)), NPTS
    end do

    if (PRNT) write (*, '(1X,A,I6,A,F9.6)') &
        '! Iterations: ', ITER, '     Final Criterion: ', SOLCRIT

    ! Add the direct solar component to the downwelling flux and mean radiance
    if (SRCTYPE /= 'T') then
        PI = acos(-1.0)
        do I = 1, NPTS
            DIRFLUX = SOLARFLUX*exp(-TAUGSC(I)/abs(SOLARMU))
            FLUXDN(I) = FLUXDN(I) + DIRFLUX
            MEANRAD(I) = MEANRAD(I) + DIRFLUX/(4*PI*abs(SOLARMU))
        end do
    end if

    ! Make an unscaled optical depth grid for the output file
    if (DELTAM) then
        call UNDELTA_SCALE(NLAY, TAUPSC, TAUP, IXP, NPTS, TAUGSC, TAUG)
    else
        TAUG(:) = TAUGSC(:)
    end if
    call OUTPUT_HEIGHT_GRID(NLAY, TAUP, IXP, NPTS, TAUG, HEIGHTP, HEIGHTG)

    if (NZOUT*NMUOUT*NPHIOUT > 0) then
        ! Compute the radiance output and put in array
        call COMPUTE_RADIANCE(NLAY, ALBEDOSC, LEGENSC, &
                              MAXLEG, NLEGP, LEGENP, DELTAM, &
                              NPTS, HEIGHTG, TAUG, TAUGSC, IXP, IXG, &
                              ML, MM, NLM, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                              SRCTYPE, SOLARFLUX, SOLARMU, &
                              SKYRAD, WAVELEN, WAVENO, SFCTEMP, &
                              SFCTYPE, SFCPARMS, SOURCE, RADIANCE, &
                              NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                              RADOUT)
    end if

    deallocate (TAUPSC, ALBEDOSC, LEGENSC, TAUGSC)
    deallocate (MU, PHI, NPHI0, WTDO, YLMSUN, WTMU)
    deallocate (CPHI1, CPHI2, WPHISAVE, CMU1, CMU2)
    deallocate (SOURCE, RADIANCE, CONSTSOURCE, DELSOURCE, SRCINT)
end subroutine SOLVE_SHDOMPP

subroutine DELTA_SCALE(NLAY, ML, TAUP, ALBEDOP, MAXLEG, NLEGP, LEGENP, &
                       TAUPSC, ALBEDOSC, LEGENSC)
    ! Delta scales the layer optical depth, albedo, and Legendre terms;
    ! only the 0 to ML LEGEN terms are scaled.
    implicit none
    integer, intent(IN) :: NLAY, ML, MAXLEG, NLEGP(NLAY)
    real, intent(IN) ::  TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
    real, intent(OUT) ::  TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    integer :: I, L
    real    :: F

    do I = 1, NLAY
        if (MAXLEG < ML + 1) then
            F = 0.0
            LEGENSC(:, I) = 0.0
        else
            F = LEGENP(ML + 1, I)/(2*(ML + 1) + 1)
        end if
        do L = 0, min(ML, MAXLEG)
            LEGENSC(L, I) = (2*L + 1)*(LEGENP(L, I)/(2*L + 1) - F)/(1 - F)
        end do
        TAUPSC(I) = (1.0 - ALBEDOP(I)*F)*TAUP(I)
        ALBEDOSC(I) = (1.0 - F)*ALBEDOP(I)/(1.0 - ALBEDOP(I)*F)
    end do
end subroutine DELTA_SCALE

subroutine UNDELTA_SCALE(NLAY, TAUPSC, TAUP, IXP, NPTS, TAUGSC, TAUG)
    ! Undoes the delta scaling to get the optical depth grid in the
    ! original scale.  Uses the scaled and unscaled input layer optical depths.
    implicit none
    integer, intent(IN) :: NLAY, IXP(2, NLAY), NPTS
    real, intent(IN) :: TAUPSC(NLAY), TAUP(NLAY), TAUGSC(NPTS)
    real, intent(OUT) :: TAUG(NPTS)
    integer :: LAY, I
    real    :: R

    TAUG(1) = 0.0
    do LAY = 1, NLAY
        if (TAUPSC(LAY) > 0.0) then
            R = TAUP(LAY)/TAUPSC(LAY)
        else
            R = 1.0
        end if
        do I = IXP(1, LAY), IXP(2, LAY) - 1
            TAUG(I + 1) = TAUG(I) + R*(TAUGSC(I + 1) - TAUGSC(I))
        end do
        if (LAY < NLAY) TAUG(IXP(1, LAY + 1)) = TAUG(IXP(2, LAY))
    end do
end subroutine UNDELTA_SCALE

subroutine OUTPUT_HEIGHT_GRID(NLAY, TAUP, IXP, NPTS, TAUG, HEIGHTP, HEIGHTG)
    ! Makes the output height grid (HEIGHTG) by interpolating between the
    ! levels in the input height grid (HEIGHTP) according to the optical depth
    ! grid (TAUG) assuming constant extinction in each input layer.
    ! The levels in HEIGHTG correspond to the optical depths in TAUG.
    implicit none
    integer, intent(IN) :: NLAY, IXP(2, NLAY), NPTS
    real, intent(IN) :: TAUP(NLAY), TAUG(NPTS), HEIGHTP(NLAY + 1)
    real, intent(OUT) :: HEIGHTG(NPTS)
    integer :: LAY, I
    real    :: E

    do LAY = 1, NLAY
        HEIGHTG(IXP(1, LAY)) = HEIGHTP(LAY)
        HEIGHTG(IXP(2, LAY)) = HEIGHTP(LAY + 1)
        if (TAUP(LAY) > 0.0) then
            E = (HEIGHTP(LAY) - HEIGHTP(LAY + 1))/TAUP(LAY)
            do I = IXP(1, LAY) + 1, IXP(2, LAY) - 1
                HEIGHTG(I) = HEIGHTP(LAY) - E*(TAUG(I) - TAUG(IXP(1, LAY)))
            end do
        else
            do I = IXP(1, LAY) + 1, IXP(2, LAY) - 1
                HEIGHTG(I) = HEIGHTP(LAY) + (HEIGHTP(LAY + 1) - HEIGHTP(LAY))* &
                             (I - IXP(1, LAY))/FLOAT(IXP(2, LAY) - IXP(1, LAY))
            end do
        end if
    end do
end subroutine OUTPUT_HEIGHT_GRID

subroutine MAKE_TAU_GRID(NLAY, TAUPSC, ALBEDOSC, MAXDELTAU, MINDELTAUSPLIT, &
                         NPTS, TAUGSC, IXP, IXG)
    ! Makes the initial optical depth grid, TAUGSC, of increasing values
    ! (representing decreasing altitudes in the atmosphere). TAUGSC(1) is zero.
    ! The "base" tau grid is the input property grid boundaries and even
    ! subdivision of each of these layers so that the sublayers are less
    ! than MAXDELTAU scattering optical depth.   All layers with scaled
    ! optical depth greater than MINDELTAUSPLIT are split into two sublayers.
    ! The IXP index array points to the top and bottom TAUGSC levels for each
    ! property layer. The IXG index array points to the source and
    ! radiance arrays for each TAUGSC level.
    implicit none
    integer, intent(IN)  :: NLAY
    real, intent(IN)  :: TAUPSC(NLAY), ALBEDOSC(NLAY)
    real, intent(IN)  :: MAXDELTAU, MINDELTAUSPLIT
    real, intent(OUT) :: TAUGSC(*)
    integer, intent(OUT) :: NPTS, IXP(2, NLAY), IXG(*)
    integer :: I, L, N
    real    :: TAU

    NPTS = 0
    TAU = 0.0
    do L = 1, NLAY
        NPTS = NPTS + 1
        TAUGSC(NPTS) = TAU
        IXP(1, L) = NPTS
        if (TAUPSC(L) > MINDELTAUSPLIT) then
            N = max(1, int(ALBEDOSC(L)*TAUPSC(L)/MAXDELTAU))
        else
            N = 0
        end if
        do I = 1, N
            NPTS = NPTS + 1
            TAUGSC(NPTS) = TAUGSC(IXP(1, L)) + TAUPSC(L)*FLOAT(I)/(N + 1)
        end do
        NPTS = NPTS + 1
        IXP(2, L) = NPTS
        TAUGSC(NPTS) = TAUGSC(IXP(1, L)) + TAUPSC(L)
        TAU = TAUGSC(NPTS)
    end do
    IXG(1:NPTS) = (/(I, I=1, NPTS)/)
end subroutine MAKE_TAU_GRID

subroutine PREPARE_SOURCE(NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                          SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                          SOLARFLUX, SOLARMU, YLMSUN, WAVELEN, WAVENO, &
                          CONSTSOURCE)
    ! Prepares the constant part of the source function for each level.
    ! This is either or both thermal emission or the single scattering
    ! solar pseudo-source.  Only those levels (points) with CONSTSOURCE(1,:)
    ! less than zero are changed, which allows this to be run for each time
    ! adaptive grid points are added.
    ! The thermal emission source is \sqrt{4\pi}(1-\omega)B(T), where the
    ! temperature T is assumed to vary linearly in optical depth.
    ! The solar source is
    !   S_0/\mu_0 e^{-\tau/\mu_0} Y_{lm}(\mu_0,\phi_0) {\omega \chi_l \over 2l+1}
    ! where S_0 is the solar flux on a horizontal plane, \chi_l are the Legendre
    ! coefficients, and \tau is the optical depth at the level.
    implicit none
    integer, intent(IN) :: NLAY, ML, MM, NLM, NPTS
    integer, intent(IN) :: IXP(2, NLAY), IXG(NPTS)
    real, intent(IN) :: SOLARFLUX, SOLARMU, YLMSUN(NLM), WAVELEN, WAVENO(2)
    real, intent(IN) :: TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(IN) :: TEMPP(NLAY + 1), TAUGSC(NPTS)
    real, intent(INOUT) :: CONSTSOURCE(NLM, NPTS)
    character(LEN=1), intent(IN) :: SRCTYPE
    integer :: I, J, K, L, LAY, ME, M, LOFJ(NLM)
    real    :: S0, C, TAU1, U, TEMP, BB, D, PLANCKFUNC

    S0 = SOLARFLUX/abs(SOLARMU)
    C = sqrt(4.0*acos(-1.0))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        do M = 0, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    do LAY = 1, NLAY
        TAU1 = TAUGSC(IXP(1, LAY))
        do I = IXP(1, LAY), IXP(2, LAY)
            K = IXG(I)
            if (CONSTSOURCE(1, K) < 0) then
                CONSTSOURCE(:, K) = 0.0
                if (SRCTYPE /= 'S') then
                    if (TAUPSC(LAY) > 0.0) then
                        U = (TAUGSC(I) - TAU1)/TAUPSC(LAY)
                    else
                        U = (I - IXP(1, LAY))/FLOAT(IXP(2, LAY) - IXP(1, LAY))
                    end if
                    TEMP = TEMPP(LAY) + U*(TEMPP(LAY + 1) - TEMPP(LAY))
                    CONSTSOURCE(1, K) = C*(1.0 - ALBEDOSC(LAY)) &
                                        *PLANCKFUNC(WAVELEN, WAVENO, TEMP)
                end if
                if (SRCTYPE /= 'T') then
                    D = S0*exp(-TAUGSC(I)/abs(SOLARMU))
                    do J = 1, NLM
                        L = LOFJ(J)
                        CONSTSOURCE(J, K) = CONSTSOURCE(J, K) &
                                            + D*YLMSUN(J)*ALBEDOSC(LAY)*LEGENSC(L, LAY)/(2*L + 1)
                    end do
                end if
            end if
        end do
    end do
end subroutine PREPARE_SOURCE

function PLANCKFUNC(WAVELEN, WAVENO, TEMP)
    ! Calculates the Planck blackbody radiance. If WAVELEN<=0 then doing a
    ! spectral band integration and the Planck blackbody radiance in
    ! [Watts /(meter^2 ster)] over a wavenumber range [cm^-1] (WAVENO) is
    ! returned.  Otherwise, the Planck blackbody radiance in
    ! [Watts /(meter^2 ster micron)] for a temperature in [Kelvins] at
    ! a wavelength (WAVELEN) in [microns] is returned.
    implicit none
    real, intent(IN) :: WAVELEN, WAVENO(2), TEMP
    real :: PLANCKFUNC
    double precision :: X1, X2, F

    if (TEMP > 0.0) then
        if (WAVELEN > 0) then
            PLANCKFUNC = 1.1911e8/WAVELEN**5/(exp(1.4388e4/(WAVELEN*TEMP)) - 1)
        else
            X1 = 1.4388d0*WAVENO(1)/TEMP
            X2 = 1.4388d0*WAVENO(2)/TEMP
            call INTEGRATE_PLANCK(X1, X2, F)
            PLANCKFUNC = 1.1911d-8*(TEMP/1.4388d0)**4*F
        end if
    else
        PLANCKFUNC = 0.0
    end if
end function PLANCKFUNC

subroutine MAKE_ANGLE_SET(NMU, NPHI, NCS, NPHI0MAX, ITYPE, &
                          NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)
    ! Make the set of angles for the discrete space representation.
    ! The number of mu's (cosine zenith angles) and maximum number of
    ! phi's (azimuth angles) is input.  The actual number of azimuthal
    ! angles for each mu is output (NPHI0).  There are three types of
    ! discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
    ! gaussian grid, and 3 is a reduced double gaussian set.
    ! If NCS=1 then only about half the azimuthal angles (from 0 to pi)
    ! are used because the radiance is even in phi (cosine terms).
    ! The output is the NMU mu values, the NPHI0 phi values for each mu,
    ! and the integration weight for each ordinate. The first NMU/2 mu
    ! angles are the downwelling (mu<0) angles.  Also output are the
    ! the flags for doing an azimuthal FFT for each mu and the total
    ! number of angles (NANG).
    implicit none
    integer, intent(IN)  :: NMU, NPHI, NCS, NPHI0MAX, ITYPE
    integer, intent(OUT) :: NPHI0(NMU), NANG
    logical, intent(OUT) :: FFTFLAG(NMU)
    real, intent(OUT) :: MU(NMU), WTMU(NMU)
    real, intent(OUT) :: PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    integer :: J, K, MM
    real    :: DELPHI
    integer, parameter :: MAXNPHI = 1024
    integer :: GOODNFFT(MAXNPHI)

    call GOOD_FFT_LENGTHS(MAXNPHI, GOODNFFT)
    GOODNFFT(1:14) = (/(K, K=1, 14)/)

    MM = max(0, int(NPHI/2) - 1)
    NANG = 0
    if (ITYPE .ge. 1 .or. ITYPE .le. 3) then
        if (ITYPE .le. 2) then
            call GAUSQUADS(NMU, MU, WTMU)
        else if (ITYPE .eq. 3) then
            call DGAUSQUADS(NMU, MU, WTMU)
        end if
        do J = 1, NMU
            if (NCS .eq. 1) then
                if (NPHI0MAX .ne. int((NPHI + 2)/2)) stop 'MAKE_ANGLE_SET: bad NPHI0MAX'
                if (ITYPE .eq. 1) then
                    NPHI0(J) = NPHI0MAX
                else
                    ! For the reduced gaussian set, make the smaller NPHI0 values
                    ! still be good for the FFT (but don't let NPHI0MAX be exceeded).
                    if (NPHI0MAX .gt. MAXNPHI) stop 'MAKE_ANGLE_SET: exceeded GOODNFFT'
                    NPHI0(J) = int(0.9 + 1 + (NPHI0MAX - 1)*sqrt(1 - MU(J)**2))
                    if (NPHI0(J) .gt. 1) NPHI0(J) = min(NPHI0MAX, GOODNFFT(NPHI0(J) - 1) + 1)
                end if
                ! Compute the azimuth angles and weights
                DELPHI = acos(-1.0)/max(1, NPHI0(J) - 1)
                do K = 1, NPHI0(J)
                    PHI(K, J) = (K - 1)*DELPHI
                    if ((K .eq. 1 .or. K .eq. NPHI0(J)) .and. NPHI0(J) .ne. 1) then
                        WTDO(J, K) = DELPHI*WTMU(J)
                    else
                        WTDO(J, K) = 2.0*DELPHI*WTMU(J)
                    end if
                end do
            else
                if (NPHI0MAX .ne. NPHI) stop 'MAKE_ANGLE_SET: bad NPHI0MAX'
                if (ITYPE .eq. 1) then
                    NPHI0(J) = NPHI0MAX
                else
                    if (NPHI0MAX .gt. MAXNPHI) stop 'MAKE_ANGLE_SET: exceeded GOODNFFT'
                    NPHI0(J) = int(0.9 + NPHI0MAX*sqrt(1 - MU(J)**2))
                    NPHI0(J) = min(NPHI0MAX, GOODNFFT(NPHI0(J)))
                end if
                DELPHI = 2.0*acos(-1.0)/NPHI0(J)
                do K = 1, NPHI0(J)
                    PHI(K, J) = (K - 1)*DELPHI
                    WTDO(J, K) = DELPHI*WTMU(J)
                end do
            end if
            NANG = NANG + NPHI0(J)
            FFTFLAG(J) = (NPHI0(J) .gt. 14) .or. (MM .gt. 16)
        end do

    else
        stop 'MAKE_ANGLE_SET: invalid discrete ordinate type'
    end if
end subroutine MAKE_ANGLE_SET

subroutine GOOD_FFT_LENGTHS(NMAX, GOODNFFT)
    ! Returns a list in GOODNFFT(I) that rounds up to the nearest I that
    ! has only divisors of 2 or 3.
    implicit none
    integer, intent(IN)  :: NMAX
    integer, intent(OUT) :: GOODNFFT(NMAX)
    integer :: I, J, ILAST
    integer, parameter :: PRIMES(310) = &
                          (/2, 3, 5, 7, 11, 13, 17, 19, 23, 29, &
                            31, 37, 41, 43, 47, 53, 59, 61, 67, 71, &
                            73, 79, 83, 89, 97, 101, 103, 107, 109, 113, &
                            127, 131, 137, 139, 149, 151, 157, 163, 167, 173, &
                            179, 181, 191, 193, 197, 199, 211, 223, 227, 229, &
                            233, 239, 241, 251, 257, 263, 269, 271, 277, 281, &
                            283, 293, 307, 311, 313, 317, 331, 337, 347, 349, &
                            353, 359, 367, 373, 379, 383, 389, 397, 401, 409, &
                            419, 421, 431, 433, 439, 443, 449, 457, 461, 463, &
                            467, 479, 487, 491, 499, 503, 509, 521, 523, 541, &
                            547, 557, 563, 569, 571, 577, 587, 593, 599, 601, &
                            607, 613, 617, 619, 631, 641, 643, 647, 653, 659, &
                            661, 673, 677, 683, 691, 701, 709, 719, 727, 733, &
                            739, 743, 751, 757, 761, 769, 773, 787, 797, 809, &
                            811, 821, 823, 827, 829, 839, 853, 857, 859, 863, &
                            877, 881, 883, 887, 907, 911, 919, 929, 937, 941, &
                            947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, &
                            1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
                            1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
                            1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, &
                            1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
                            1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
                            1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
                            1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
                            1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
                            1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
                            1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
                            1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
                            1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
                            1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, &
                            1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053/)

    GOODNFFT(:) = (/(I, I=1, NMAX)/)
    do I = 3, size(PRIMES)
        do J = PRIMES(I), NMAX, PRIMES(I)
            GOODNFFT(J) = 0
        end do
    end do
    ILAST = NMAX
    if (GOODNFFT(ILAST) == 0) stop 'GOOD_FFT_LENGTHS: NMAX not good'
    do I = NMAX, 1, -1
        if (GOODNFFT(I) == 0) then
            GOODNFFT(I) = ILAST
        else
            ILAST = I
        end if
    end do
end subroutine GOOD_FFT_LENGTHS

subroutine MAKE_SH_DO_COEF(ML, MM, NLM, NMU, NPHI0, &
                           NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
                           CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)
    ! Makes the transformation coefficients for the spherical harmonic
    ! transform.  The main coefficients are output in four arrays: 1 is for
    ! the SH_TO_DO forward transform, 2 is for the DO_TO_SH back transform,
    ! which contains the discrete angle integration weights.
    ! The coefficients are divided up into the mu dependent set CMUn
    ! (function of l, m, mu_j), and the phi dependent set CPHIn
    ! (function of m, phi_k) for each mu_j.
    ! The FFTPACK phase coefficients for the FFT in azimuth are also
    ! output in WPHISAVE.
    implicit none
    integer, intent(IN) :: ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX
    logical, intent(IN) :: FFTFLAG(NMU)
    real, intent(IN) :: MU(NMU), WTMU(NMU)
    real, intent(IN) :: PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(OUT) :: CMU1(NLM, NMU), CMU2(NMU, NLM)
    real, intent(OUT) :: CPHI1(0:16, 32, NMU), CPHI2(32, 0:16, NMU)
    real, intent(OUT) :: WPHISAVE(3*NPHI0MAX + 15, NMU)
    integer :: I, K, M, ier
    real    :: W

    ! Make the to and from associate Legendre coefficients
    do I = 1, NMU
        call YLMALL(MU(I), 0.0, ML, MM, -1, CMU1(1, I))
    end do
    do I = 1, NMU
        CMU2(I, 1:NLM) = CMU1(1:NLM, I)*WTMU(I)
    end do

    ! Make the to and from Fourier coefficients for each mu
    do I = 1, NMU
        if (.not. FFTFLAG(I)) then
            ! If not doing an FFT for this mu then make the DFT coefficient
            W = 1.0/WTMU(I)
            do K = 1, NPHI0(I)
                if (NPHI0(I) .gt. 32 .or. MM .gt. 16) &
                    stop 'MAKE_SH_DO_COEF: Fourier coefficients array exceeded'
                do M = 0, MM
                    CPHI1(M, K, I) = cos(M*PHI(K, I))
                    CPHI2(K, M, I) = CPHI1(M, K, I)*WTDO(I, K)*W
                end do
            end do
        else
            ! Precompute the phase factors for the FFTs
            call cost1i(NPHI0(I), WPHISAVE(1, I), 3*NPHI0MAX + 15, ier)
            if (ier /= 0) stop
        end if
    end do
end subroutine MAKE_SH_DO_COEF

subroutine SH_TO_DO(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                    FFTFLAG, CMU1, CPHI1, WSAVE, INDATA, OUTDATA)
    ! Transforms the input data from spherical harmonic space to discrete
    ! ordinate space for all points in the array.  Successive transforms are
    ! done in zenith and then azimuth angle.  The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small primes
    ! (powers of 2 are best).  The number of floating point operations for
    ! the DFT is 2*Nmu*Nphi0*Nm. The number of floating point operations for
    ! the zenith angle transform is 2*Nmu*Nlm.
    implicit none
    integer, intent(IN) :: NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
    logical, intent(IN) :: FFTFLAG(NMU)
    real, intent(IN) :: INDATA(NLM, NPTS)
    real, intent(IN) :: CMU1(NLM, NMU), CPHI1(0:16, 32, NMU)
    real, intent(IN) :: WSAVE(3*NPHI0MAX + 15, NMU)
    real, intent(OUT) :: OUTDATA(NPHI0MAX, NMU, NPTS)
    integer :: I, J, L, M, ME, IMU, IPHI, ier
    integer :: MOFJ(NLM)
    real    :: S(0:MM), WORK(NPHI0MAX)

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = 0, min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        do IMU = 1, NMU
            ! Figure the max Fourier azimuthal mode we can do for this Nphi0
            ME = max(0, min(NPHI0(IMU) - 2, MOFJ(NLM)))
            ! First do Legendre transform by summing over l for each m.
            S(:) = 0.0
            do J = 1, NLM
                M = MOFJ(J)
                S(M) = S(M) + CMU1(J, IMU)*INDATA(J, I)
            end do
            ! Then do the Fourier transform from m to phi for each mu
            if (FFTFLAG(IMU)) then
                ! For cosine modes only use a cosine FFT
                OUTDATA(1:ME + 1, IMU, I) = S(0:ME)
                OUTDATA(ME + 2:NPHI0(IMU), IMU, I) = 0.0
                call cost1b(NPHI0(IMU), 1, OUTDATA(:, IMU, I), NPHI0MAX, &
                            WSAVE(:, IMU), 3*NPHI0MAX + 15, WORK, NPHI0MAX, ier)
                if (ier /= 0) stop
            else
                ! Else do a slow DFT
                do IPHI = 1, NPHI0(IMU)
                    OUTDATA(IPHI, IMU, I) = sum(CPHI1(0:ME, IPHI, IMU)*S(0:ME))
                end do
            end if
        end do
    end do
end subroutine SH_TO_DO

subroutine DO_TO_SH(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                    FFTFLAG, CMU2, CPHI2, WSAVE, INDATA, OUTDATA)
    ! Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for all the points.  (See SH_TO_DO for more).
    ! The SELECTPNTS logical array indicates whether a point should be
    ! transformed or to use the previous points (for the input layer interfaces).
    implicit none
    integer, intent(IN) :: NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
    logical, intent(IN) :: FFTFLAG(NMU)
    real, intent(IN) :: INDATA(NPHI0MAX, NMU, NPTS)
    real, intent(IN) :: CMU2(NMU, NLM), CPHI2(32, 0:16, NMU)
    real, intent(IN) :: WSAVE(3*NPHI0MAX + 15, NMU)
    real, intent(OUT) :: OUTDATA(NLM, NPTS)
    integer :: I, IMU, IPHI, J, L, M, ME, N, ier
    real    :: PI, DELPHI
    integer :: MOFJ(NLM)
    real    :: S(0:NPHI0MAX - 1), WORK(NPHI0MAX)

    PI = acos(-1.0)
    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = 0, min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        OUTDATA(:, I) = 0.0
        do IMU = 1, NMU
            N = NPHI0(IMU)
            ! Figure the max Fourier azimuthal mode we can do for this Nphi0
            ME = max(0, min(N - 2, MOFJ(NLM)))
            ! First do Fourier transform from phi to m for each mu
            if (FFTFLAG(IMU)) then
                ! For cosine modes only use a cosine FFT
                S(0:N - 1) = INDATA(1:N, IMU, I)
                call cost1f(N, 1, S(0:N - 1), NPHI0MAX, &
                            WSAVE(:, IMU), 3*NPHI0MAX + 15, WORK, NPHI0MAX, ier)
                if (ier /= 0) stop
                S(0:ME) = PI*S(0:ME)
                S(0) = 2.0*S(0)
                S(ME) = 2.0*S(ME)
            else
                ! Else do a slow DFT
                S(:) = 0.0
                do M = 0, ME
                    S(M) = sum(CPHI2(1:N, M, IMU)*INDATA(1:N, IMU, I))
                end do
            end if
            ! Then do Legendre transform by adding to output for each l, m
            do J = 1, NLM
                M = MOFJ(J)
                OUTDATA(J, I) = OUTDATA(J, I) + CMU2(IMU, J)*S(M)
            end do
        end do
    end do
end subroutine DO_TO_SH

subroutine COMPUTE_SOURCE(ML, MM, NLM, NLAY, NPTS, IXP, IXG, &
                          ALBEDOSC, LEGENSC, &
                          SHRADIANCE, CONSTSOURCE, SOURCE, &
                          OLDNPTS, DELSOURCE, &
                          DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
    ! Computes the source function (SOURCE) in spherical harmonic space
    ! for all the grid levels in all the input layers.
    ! The thermal source and/or solar pseudo-sources in CONSTSOURCE is
    ! added to the scattering source (computed from LEGENP and the spherical
    ! harmonic expansion in SHRADIANCE).
    ! To save memory this routine also computes some things for the sequence
    ! acceleration: the dot products of the difference in source function
    ! between successive iterations and the new difference source function
    ! vector (DELSOURCE).
    ! Computes the new solution criterion (SOLCRIT, the RMS difference
    ! in succesive source function fields normalized by the RMS of the field
    implicit none
    integer, intent(IN) :: ML, MM, NLM, NLAY, NPTS, OLDNPTS
    integer, intent(IN) :: IXP(2, NLAY), IXG(NPTS)
    real, intent(IN) :: ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(IN) :: SHRADIANCE(NLM, NPTS), CONSTSOURCE(NLM, NPTS)
    real, intent(INOUT) :: SOURCE(NLM, NPTS), DELSOURCE(NLM, OLDNPTS)
    real, intent(OUT) :: DELJDOT, DELJOLD, DELJNEW, SOLCRIT
    integer :: LOFJ(NLM), I, J, K, L, LAY, M, ME
    real    :: SOURCET, NORM

    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        do M = 0, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    DELJDOT = 0.0
    DELJOLD = 0.0
    DELJNEW = 0.0
    NORM = 0.0
    do LAY = 1, NLAY
        do I = IXP(1, LAY), IXP(2, LAY)
            K = IXG(I)
            do J = 1, NLM
                L = LOFJ(J)
                SOURCET = CONSTSOURCE(J, K) &
                          + SHRADIANCE(J, K)*ALBEDOSC(LAY)*LEGENSC(L, LAY)/(2*L + 1)
                if (K <= OLDNPTS) then
                    DELJDOT = DELJDOT + (SOURCET - SOURCE(J, K))*DELSOURCE(J, K)
                    DELJOLD = DELJOLD + DELSOURCE(J, K)**2
                    DELSOURCE(J, K) = SOURCET - SOURCE(J, K)
                    DELJNEW = DELJNEW + (SOURCET - SOURCE(J, K))**2
                    NORM = NORM + SOURCET**2
                end if
                SOURCE(J, K) = SOURCET
            end do
        end do
    end do
    if (NORM > 0.0) then
        SOLCRIT = sqrt(DELJNEW/NORM)
    else
        SOLCRIT = 1.0
    end if
end subroutine COMPUTE_SOURCE

subroutine ACCELERATE_SOLUTION(NLM, NPTS, OLDNPTS, SOURCE, DELSOURCE, &
                               DELJDOT, DELJOLD, DELJNEW)
    ! Accelerates the successive order of scattering series.  Uses
    ! information about the ratio in lengths and angle between successive
    ! source function differences to decide how far along the current
    ! source function difference to extrapolate.
    implicit none
    integer, intent(IN) :: NLM, NPTS, OLDNPTS
    real, intent(IN) :: DELSOURCE(NLM, OLDNPTS)
    real, intent(INOUT) :: SOURCE(NLM, NPTS)
    real, intent(IN) :: DELJDOT, DELJOLD, DELJNEW
    integer :: I, J, IS, ISD, NS, NSD, NSC
    real, save :: A, R, THETA

    ! Accelerate if the grids are the same, didn't last time, and
    !   things are converging.
    if (NPTS == OLDNPTS .and. A == 0.0 .and. DELJNEW < DELJOLD) then
        ! Compute the acceleration extrapolation factor and apply it.
        R = sqrt(DELJNEW/DELJOLD)
        THETA = acos(DELJDOT/sqrt(DELJOLD*DELJNEW))
        A = (1 - R*cos(THETA) + R**(1 + 0.5*3.14159/THETA)) &
            /(1 + R**2 - 2*R*cos(THETA)) - 1.0
        A = max(0.0, A)
        ! WRITE (*,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
        SOURCE(:, :) = SOURCE(:, :) + A*DELSOURCE(:, :)
    else
        A = 0.0
    end if
end subroutine ACCELERATE_SOLUTION

subroutine PATH_INTEGRATION(NLAY, NPTS, IXP, IXG, TAUGSC, &
                            NMU, NPHI0MAX, NPHI0, &
                            SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                            SFCTYPE, SFCPARMS, SFCTEMP, WAVELEN, WAVENO, &
                            MU, PHI, WTDO, SRCINT, LASTNPTS, &
                            FLUXUP, FLUXDN, MEANRAD, DOSOURCE, RADIANCE)
    ! Performs the path integrations through the medium specified by
    ! the optical depth (TAUGSC) and source function (DOSOURCE) in
    ! discrete ordinates.  The integrations are done from the top boundary
    ! (with no reflection) down to the bottom boundary for the downwelling
    ! angles, and then up from the bottom to the top after computing the
    ! surface reflection/emission.   The hemispheric flux (FLUXES) and mean
    ! radiance (MEANRAD) are also computed.
    ! The discrete ordinates in MU, PHI must be a complete set with the
    ! downwelling (mu<0) angles first and downwelling and upwelling matching
    ! [MU(J) = -MU(NMU/2+J)].
    ! The surface reflection is handled differently for Lambertian surfaces
    ! and more general surfaces specified by bidirectional reflection
    ! distribution functions (BRDF).  For the general BRDF, the bottom
    ! boundary radiance must be computed for each upwelling ordinate, while
    ! for the Lambertian surface the boundary radiances can be computed just
    ! once, since they are isotropic.
    implicit none
    integer, intent(IN) :: NLAY, NPTS, IXP(2, NLAY), IXG(NPTS)
    integer, intent(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(INOUT) :: LASTNPTS
    real, intent(IN) :: SOLARFLUX, SOLARMU, SKYRAD
    real, intent(IN) :: SFCPARMS(*), SFCTEMP, WAVELEN, WAVENO(2)
    real, intent(IN) :: MU(NMU), PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(INOUT) :: SRCINT(0:4, NMU, NPTS)
    real, intent(IN) :: TAUGSC(NPTS), DOSOURCE(NPHI0MAX, NMU, NPTS)
    real, intent(OUT) :: FLUXUP(NPTS), FLUXDN(NPTS), MEANRAD(NPTS)
    real, intent(OUT) :: RADIANCE(NPHI0MAX, NMU, NPTS)
    character(LEN=1), intent(IN) :: SRCTYPE, SFCTYPE

    integer :: LAY, I, I1, I2, J, K, IPHI, IMU, M, N
    real    :: PI, TEMP, SFCPLANCK, PLANCKFUNC, ALB
    real    :: DIRFLUXSFC, TOPRAD, BOTRAD
    double precision, allocatable :: A(:, :), B(:, :), C(:, :), D(:, :), E(:)
    double precision, allocatable :: T(:), T1(:), T3(:), T4(:), TR(:)

    PI = acos(-1.0)

    M = NMU/2

    if (LASTNPTS /= NPTS) then
        ! Calculate and store the source function integration coefficients
        allocate (A(4, M), B(4, M), C(4, M), D(4, M), E(M), T(M), T1(M), T3(M), T4(M), TR(M))
        do LAY = 1, NLAY
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then ! special case for layer with just two points
                T(:) = (TAUGSC(I1 + 1) - TAUGSC(I1))/abs(MU(1:M))
                TR(:) = exp(-T(:)) ! sublayer transmission along ordinate
                if (T(1) < 1.0e-4) then
                    D(1, :) = 0.5*T(:) - (1.0d0/3)*T(:)**2
                    D(2, :) = 0.5*T(:) - (1.0d0/6)*T(:)**2
                else
                    D(1, :) = (1 - TR(:)*(1 + T(:)))/T(:)
                    D(2, :) = (T(:) - (1 - TR(:)))/T(:)
                end if
                SRCINT(0, 1:M, I1) = TR(:)
                SRCINT(0, M + 1:NMU, I1) = TR(:)
                SRCINT(1, 1:M, I1) = D(2, :)
                SRCINT(2, 1:M, I1) = D(1, :)
                SRCINT(1, M + 1:NMU, I1) = D(1, :)
                SRCINT(2, M + 1:NMU, I1) = D(2, :)
            else
                do I = I1, I2
                    A(:, :) = 0.0; B(:, :) = 0.0; C(:, :) = 0.0; D(:, :) = 0.0
                    T3(:) = (TAUGSC(I + 1) - TAUGSC(I))/abs(MU(1:M))
                    if (I == I1) then
                        T4(:) = (TAUGSC(I + 2) - TAUGSC(I))/abs(MU(1:M))
                        E = T3*T4; B(2, :) = 1/E; C(2, :) = -(T3 + T4)/E; D(2, :) = 1.0
                        E = T3*(T3 - T4); B(3, :) = 1/E; C(3, :) = -T4/E; D(3, :) = 0.0
                        E = T4*(T4 - T3); B(4, :) = 1/E; C(4, :) = -T3/E; D(4, :) = 0.0
                    else if (I < I2) then
                        T1(:) = (TAUGSC(I - 1) - TAUGSC(I))/abs(MU(1:M))
                        T4(:) = (TAUGSC(I + 2) - TAUGSC(I))/abs(MU(1:M))
                        E = T1*(T1 - T3)*(T1 - T4); A(1, :) = 1/E
                        B(1, :) = -(T3 + T4)/E; C(1, :) = T3*T4/E; D(1, :) = 0.0
                        E = -T1*T3*T4; A(2, :) = 1/E
                        B(2, :) = -(T1 + T3 + T4)/E; C(2, :) = (T1*T3 + T1*T4 + T3*T4)/E; D(2, :) = 1.0
                        E = T3*(T3 - T1)*(T3 - T4); A(3, :) = 1/E
                        B(3, :) = -(T1 + T4)/E; C(3, :) = T1*T4/E; D(3, :) = 0.0
                        E = T4*(T4 - T1)*(T4 - T3); A(4, :) = 1/E
                        B(4, :) = -(T1 + T3)/E; C(4, :) = T1*T3/E; D(4, :) = 0.0
                    else
                        T1(:) = (TAUGSC(I - 1) - TAUGSC(I))/abs(MU(1:M))
                        E = T1*(T1 - T3); B(1, :) = 1/E; C(1, :) = -T3/E; D(1, :) = 0.0
                        E = T1*T3; B(2, :) = 1/E; C(2, :) = -(T1 + T3)/E; D(2, :) = 1.0
                        E = T3*(T3 - T1); B(3, :) = 1/E; C(3, :) = -T1/E; D(3, :) = 0.0
                    end if
                    T(:) = T3(:)
                    TR(:) = exp(-T(:))
                    SRCINT(0, 1:M, I) = TR(:)
                    SRCINT(0, M + 1:NMU, I) = TR(:)
                    do J = 1, 4
                        SRCINT(J, 1:M, I) = A(J, :)*(6*TR + T**3 - 3*T**2 + 6*T - 6) &
                                            + B(J, :)*(-2*TR + T**2 - 2*T + 2) &
                                            + C(J, :)*(TR + T - 1) + D(J, :)*(1 - TR)
                        SRCINT(J, M + 1:NMU, I) = A(J, :)*(6 - TR*(T**3 + 3*T**2 + 6*T + 6)) &
                                                  + B(J, :)*(2 - TR*(T**2 + 2*T + 2)) &
                                                  + C(J, :)*(1 - TR*(T + 1)) + D(J, :)*(1 - TR)
                    end do
                end do
            end if
        end do
        deallocate (A, B, C, D, E, T, T1, T3, T4, TR)
    end if

    ! Zero the fluxes and mean radiance arrays
    FLUXUP(:) = 0.0; FLUXDN(:) = 0.0; MEANRAD(:) = 0.0

    ! Make the isotropic radiances for the top boundary
    if (SRCTYPE .eq. 'T') then
        TEMP = SKYRAD
        TOPRAD = PLANCKFUNC(WAVELEN, WAVENO, TEMP)
    else
        TOPRAD = SKYRAD
    end if
    RADIANCE(:, 1:NMU/2, IXG(IXP(1, 1))) = TOPRAD

    ! First integrate the downwelling zenith angles (these angles must be first)
    ! For downwelling ordinates integrate source function from the top down
    !   Integrations are performed over the vector of azimuths
    do IMU = 1, NMU/2
        N = NPHI0(IMU)
        do LAY = 1, NLAY
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADIANCE(1:N, IMU, IXG(I1 + 1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(1, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1))
            else
                RADIANCE(1:N, IMU, IXG(I1 + 1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(3, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1)) &
                                                  + SRCINT(4, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 2))
                do I = I1 + 1, I2 - 1
                    RADIANCE(1:N, IMU, IXG(I + 1)) = SRCINT(0, IMU, I)*RADIANCE(1:N, IMU, IXG(I)) &
                                                     + SRCINT(1, IMU, I)*DOSOURCE(1:N, IMU, IXG(I - 1)) &
                                                     + SRCINT(2, IMU, I)*DOSOURCE(1:N, IMU, IXG(I)) &
                                                     + SRCINT(3, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 1)) &
                                                     + SRCINT(4, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 2))
                end do
                RADIANCE(1:N, IMU, IXG(I2 + 1)) = SRCINT(0, IMU, I2)*RADIANCE(1:N, IMU, IXG(I2)) &
                                                  + SRCINT(1, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 - 1)) &
                                                  + SRCINT(2, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2)) &
                                                  + SRCINT(3, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 + 1))
            end if
            if (LAY < NLAY) RADIANCE(:, IMU, IXG(IXP(1, LAY + 1))) = RADIANCE(:, IMU, IXG(IXP(2, LAY)))
        end do

        ! Integrate over the DO radiances to get the downwelling fluxes
        do I = 1, NPTS
            do K = 1, NPHI0(IMU)
                FLUXDN(I) = FLUXDN(I) + abs(MU(IMU))*WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))
                MEANRAD(I) = MEANRAD(I) + WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))/(4*PI)
            end do
        end do
    end do

    ! Now do all the surface reflection stuff
    if (SRCTYPE /= 'T') then
        DIRFLUXSFC = SOLARFLUX*exp(-TAUGSC(IXP(2, NLAY))/abs(SOLARMU))
    end if
    SFCPLANCK = 0.0
    if (SRCTYPE /= 'S') then
        SFCPLANCK = PLANCKFUNC(WAVELEN, WAVENO, SFCTEMP)
    end if

    ! For a Lambertian surface calculate the isotropic surface radiance
    !   from the thermal emission or direct reflection and the reflected
    !   downwelling flux at the surface.
    if (SFCTYPE == 'L') then
        ALB = SFCPARMS(1)
        if (SRCTYPE == 'T') then
            BOTRAD = (1.0 - ALB)*SFCPLANCK
        else if (SRCTYPE == 'S') then
            BOTRAD = ALB*DIRFLUXSFC/PI
        else if (SRCTYPE .eq. 'B') then
            BOTRAD = (1.0 - ALB)*SFCPLANCK + ALB*DIRFLUXSFC/PI
        end if
        BOTRAD = BOTRAD + FLUXDN(IXP(2, NLAY))*ALB/PI
        RADIANCE(:, NMU/2 + 1:NMU, IXG(IXP(2, NLAY))) = BOTRAD
    else
        ! If not a Lambertian surface, compute the radiance for each of the
        !   upwelling ordinates by integrating the BRDF over the stored
        !   downwelling radiances.
        do IMU = NMU/2 + 1, NMU
            do IPHI = 1, NPHI0(IMU)
                call VARIABLE_BRDF_SURFACE(NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                           MU(IMU), PHI(IPHI, IMU), &
                                           SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                           SFCTYPE, SFCPARMS, SFCPLANCK, &
                                           RADIANCE(:, :, IXG(IXP(2, NLAY))), BOTRAD)
                RADIANCE(IPHI, IMU, IXG(IXP(2, NLAY))) = BOTRAD
            end do
        end do
    end if

    ! Then integrate the upwelling zenith angles
    ! For upwelling ordinates integrate source function from the bottom up
    do IMU = NMU/2 + 1, NMU
        N = NPHI0(IMU)
        do LAY = NLAY, 1, -1
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADIANCE(1:N, IMU, IXG(I1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1 + 1)) &
                                              + SRCINT(1, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                              + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1))
            else
                RADIANCE(1:N, IMU, IXG(I2)) = SRCINT(0, IMU, I2)*RADIANCE(1:N, IMU, IXG(I2 + 1)) &
                                              + SRCINT(1, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 - 1)) &
                                              + SRCINT(2, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2)) &
                                              + SRCINT(3, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 + 1))
                do I = I2 - 1, I1 + 1, -1
                    RADIANCE(1:N, IMU, IXG(I)) = SRCINT(0, IMU, I)*RADIANCE(1:N, IMU, IXG(I + 1)) &
                                                 + SRCINT(1, IMU, I)*DOSOURCE(1:N, IMU, IXG(I - 1)) &
                                                 + SRCINT(2, IMU, I)*DOSOURCE(1:N, IMU, IXG(I)) &
                                                 + SRCINT(3, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 1)) &
                                                 + SRCINT(4, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 2))
                end do
                RADIANCE(1:N, IMU, IXG(I1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1 + 1)) &
                                              + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                              + SRCINT(3, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1)) &
                                              + SRCINT(4, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 2))
            end if
            if (LAY > 1) RADIANCE(:, IMU, IXG(IXP(2, LAY - 1))) = RADIANCE(:, IMU, IXG(IXP(1, LAY)))
        end do

        ! Integrate over the DO radiances to get the upwelling fluxes
        do I = 1, NPTS
            do K = 1, NPHI0(IMU)
                FLUXUP(I) = FLUXUP(I) + abs(MU(IMU))*WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))
                MEANRAD(I) = MEANRAD(I) + WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))/(4*PI)
            end do
        end do
    end do

    LASTNPTS = NPTS
end subroutine PATH_INTEGRATION

subroutine SPLIT_CELLS(NLAY, NPTS, MAXIG, IXP, IXG, TAUGSC, &
                       NMU, NPHI0MAX, NPHI0, MU, DOSOURCE, RADIANCE, &
                       SPLITACC, TRUESOLCRIT, CURSPLITACC)
    ! Determines the current cell splitting accuracy (CURSPLITACC) from the
    ! solution accuracy (SOLCRIT) and the desired final cell splitting accuracy
    ! (SPLITACC).  Calculates the cell splitting criterion for each sublayer,
    ! The cell splitting criterion for a sublayer is the rms difference over
    ! outgoing ordinates in radiance between the original cubic source function
    ! interpolation integration (in RADIANCE) and a linear interpolation
    ! integration, normalized by the rms of RADIANCE over all ordinates and
    ! levels.  Splits cells with criterion greater than the current splitting
    ! accuracy.  Interpolates the discrete ordinate source function to the new
    ! levels.  The number of levels (NPTS), the adaptive index arrays (IXP, IXG),
    ! and the optical depth grid array (TAUGSC) are also updated.
    ! The boundary condition information is not needed for this routine because
    ! it has already been computed by PATH_INTEGRATION and stored in RADIANCE.
    ! SPLITCRITRATIO is the factor needed to divide the cell splitting criterion
    ! to convert it to a typical radiance accuracy of the split sublayer with
    ! the cubic interpolation integration.
    implicit none
    integer, intent(IN) :: NLAY, MAXIG
    integer, intent(INOUT) :: NPTS, IXP(2, NLAY), IXG(*)
    integer, intent(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
    real, intent(IN) :: MU(NMU)
    real, intent(IN) :: SPLITACC, TRUESOLCRIT
    real, intent(IN) :: RADIANCE(NPHI0MAX, NMU, NPTS)
    real, intent(INOUT) :: CURSPLITACC
    real, intent(INOUT) :: TAUGSC(*), DOSOURCE(NPHI0MAX, NMU, *)
    integer :: LAY, I, IMU, IMU2, N
    integer :: IXPN(2, NLAY), NEWNPTS, I1, I2, I3, I4
    real    :: SPLITCRIT, NORM, T1, T2, T3, T4
    double precision  :: TAU(NMU/2), TR(NMU/2)
    double precision, allocatable :: SRCINT1(:, :, :), DIFFRAD(:, :)

    NORM = sqrt(sum(RADIANCE(:, :, 1:NPTS)**2)/(NPTS*sum(NPHI0(:))))
    ! Leave immediately if splitting is not desired
    if (SPLITACC <= 0.0 .or. NORM == 0.0) return

    allocate (SRCINT1(0:2, NMU/2, NPTS), DIFFRAD(NPHI0MAX, NMU))
    DIFFRAD(:, :) = 0.0

    CURSPLITACC = max(SPLITACC, min(CURSPLITACC, TRUESOLCRIT))

    ! Calculate and store the linear source function integration coefficients
    do LAY = 1, NLAY
        do I = IXP(1, LAY), IXP(2, LAY) - 1
            TAU(:) = (TAUGSC(I + 1) - TAUGSC(I))/abs(MU(1:NMU/2))
            TR(:) = exp(-TAU(:)) ! sublayer transmission along ordinate
            SRCINT1(0, :, I) = TR(:)
            if (TAU(1) < 1.0e-4) then
                SRCINT1(1, :, I) = 0.5*TAU(:) - (1.0d0/3)*TAU(:)**2
                SRCINT1(2, :, I) = 0.5*TAU(:) - (1.0d0/6)*TAU(:)**2
            else
                SRCINT1(1, :, I) = (1 - TR(:)*(1 + TAU(:)))/TAU(:)
                SRCINT1(2, :, I) = (TAU(:) - (1 - TR(:)))/TAU(:)
            end if
        end do
    end do

    ! For each sublayer integrate the downwelling zenith angles from top to
    ! bottom and the upwelling angles from bottom to top using linear
    ! source function interpolation and difference from the previous
    ! cubic/quadratic interpolation results.
    NEWNPTS = NPTS
    IXPN(:, :) = IXP(:, :)
    do LAY = 1, NLAY
        if (IXP(1, LAY) < IXP(2, LAY) - 1) then
            do I = IXP(1, LAY), IXP(2, LAY) - 1
                do IMU = 1, NMU/2
                    IMU2 = IMU + NMU/2
                    N = NPHI0(IMU)
                    DIFFRAD(1:N, IMU) = -RADIANCE(1:N, IMU, IXG(I + 1)) &
                                        + SRCINT1(0, IMU, I)*RADIANCE(1:N, IMU, IXG(I)) &
                                        + SRCINT1(1, IMU, I)*DOSOURCE(1:N, IMU, IXG(I)) &
                                        + SRCINT1(2, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 1))
                    DIFFRAD(1:N, IMU2) = -RADIANCE(1:N, IMU2, IXG(I)) &
                                         + SRCINT1(0, IMU, I)*RADIANCE(1:N, IMU2, IXG(I + 1)) &
                                         + SRCINT1(1, IMU, I)*DOSOURCE(1:N, IMU2, IXG(I + 1)) &
                                         + SRCINT1(2, IMU, I)*DOSOURCE(1:N, IMU2, IXG(I))
                end do
                SPLITCRIT = sqrt(sum(DIFFRAD(:, :)**2)/sum(NPHI0(:)))/NORM
                ! If the splitting criterion is large enough then add a level in
                !   the center of the current sublayer.  Put the new level at the
                !   end of the radiance array and update the pointers.
                if (SPLITCRIT > CURSPLITACC) then
                    if (NEWNPTS >= MAXIG) then
                        write (*, *)
                        write (*, *) 'Adaptive grid exceeded MAXIG=', MAXIG
                        write (*, *) 'Increase SPLITTING_FACTOR'
                        stop
                    end if
                    NEWNPTS = NEWNPTS + 1
                    TAUGSC(NEWNPTS) = 0.5*(TAUGSC(I) + TAUGSC(I + 1))
                    ! print *, newnpts,taug(newnpts), splitcrit
                    IXG(NEWNPTS) = NEWNPTS
                    IXPN(2, LAY:NLAY) = IXPN(2, LAY:NLAY) + 1
                    IXPN(1, LAY + 1:NLAY) = IXPN(1, LAY + 1:NLAY) + 1
                    I1 = max(I - 1, IXP(1, LAY)); I2 = I
                    I3 = I + 1; I4 = min(I + 2, IXP(2, LAY))
                    call INTERP_SOURCE(NPHI0MAX, NMU, TAUGSC(I1), DOSOURCE(:, :, IXG(I1)), &
                                       TAUGSC(I2), DOSOURCE(:, :, IXG(I2)), &
                                       TAUGSC(I3), DOSOURCE(:, :, IXG(I3)), &
                                       TAUGSC(I4), DOSOURCE(:, :, IXG(I4)), &
                                       TAUGSC(NEWNPTS), DOSOURCE(:, :, IXG(NEWNPTS)))
                end if
            end do
        end if
    end do

    if (NEWNPTS > NPTS) then
        call ssort(TAUGSC, IXG, NEWNPTS, 2)
        NPTS = NEWNPTS
        IXP(:, :) = IXPN(:, :)
        ! Make sure that the IXG pointers are in the correct order
        do I = 1, NPTS - 1
            if (TAUGSC(I) == TAUGSC(I + 1) .and. IXG(I) > IXG(I + 1)) then
                I1 = IXG(I)
                IXG(I) = IXG(I + 1)
                IXG(I + 1) = I1
            end if
        end do
    end if
    deallocate (SRCINT1, DIFFRAD)
end subroutine SPLIT_CELLS

subroutine ssort (X, Y, N, KFLAG)
    ! ***BEGIN PROLOGUE  SSORT
    ! ***PURPOSE  Sort an array and optionally make the same interchanges in
    !             an auxiliary array.  The array may be sorted in increasing
    !             or decreasing order.  A slightly modified QUICKSORT
    !             algorithm is used.
    ! ***LIBRARY   SLATEC
    ! ***CATEGORY  N6A2B
    ! ***TYPE      SINGLE PRECISION (SSORT-S, DSORT-D, ISORT-I)
    ! ***KEYWORDS  SINGLETON QUICKSORT, SORT, SORTING
    ! ***AUTHOR  Jones, R. E., (SNLA)
    !            Wisniewski, J. A., (SNLA)
    ! ***DESCRIPTION
    !
    !    SSORT sorts array X and optionally makes the same interchanges in
    !    array Y.  The array X may be sorted in increasing order or
    !    decreasing order.  A slightly modified quicksort algorithm is used.
    !
    !    Description of Parameters
    !       X - array of values to be sorted   (usually abscissas)
    !       Y - array to be (optionally) carried along
    !       N - number of values in array X to be sorted
    !       KFLAG - control parameter
    !             =  2  means sort X in increasing order and carry Y along.
    !             =  1  means sort X in increasing order (ignoring Y)
    !             = -1  means sort X in decreasing order (ignoring Y)
    !             = -2  means sort X in decreasing order and carry Y along.
    !
    ! ***REFERENCES  R. C. Singleton, Algorithm 347, An efficient algorithm
    !                  for sorting with minimal storage, Communications of
    !                  the ACM, 12, 3 (1969), pp. 185-187.
    ! ***END PROLOGUE  SSORT
!      .. Scalar Arguments ..
      INTEGER KFLAG, N
!      .. Array Arguments ..
!       REAL X(*), Y(*)
      REAL X(*)
      INTEGER Y(*)
!      .. Local Scalars ..
!       REAL R, T, TT, TTY, TY
      REAL R, T, TT
      INTEGER TY, TTY
      INTEGER I, IJ, J, K, KK, L, M, NN
!      .. Local Arrays ..
      INTEGER IL(51), IU(51)
!      .. External Subroutines ..
!      .. Intrinsic Functions ..
      INTRINSIC ABS, INT
! ***First executable statement  SSORT
      NN = N
      IF (NN .LT. 1) THEN
         STOP 'The number of values to be sorted is not positive.'
      ENDIF

      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        STOP 'The sort control parameter, K, is not 2, 1, -1, or -2.'
      ENDIF

!      Alter array X to get decreasing order if needed
      IF (KFLAG .LE. -1) THEN
         DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF

      IF (KK .EQ. 2) GO TO 100

!      Sort X only
      M = 1
      I = 1
      J = NN
      R = 0.375E0

   20 IF (I .EQ. J) GO TO 60
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF

   30 K = I

!      Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)

!      If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J

!      If last element of array is less than than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)

!         If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF

!      Find an element in the second half of the array which is smaller
!      than T
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40

!      Find an element in the first half of the array which is greater
!      than T
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50

!      Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF

!      Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 70

!      Begin again on another portion of the unsorted array
   60 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

   70 IF (J-I .GE. 1) GO TO 30
      IF (I .EQ. 1) GO TO 20
      I = I-1

   80 I = I+1
      IF (I .EQ. J) GO TO 60
      T = X(I+1)
      IF (X(I) .LE. T) GO TO 80
      K = I

   90 X(K+1) = X(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 90
      X(K+1) = T
      GO TO 80

!      Sort X and carry Y along
  100 M = 1
      I = 1
      J = NN
      R = 0.375E0

  110 IF (I .EQ. J) GO TO 150
      IF (R .LE. 0.5898437E0) THEN
         R = R+3.90625E-2
      ELSE
         R = R-0.21875E0
      ENDIF
!
  120 K = I

!      Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)

!      If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J

!      If last element of array is less than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)

!         If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF

!      Find an element in the second half of the array which is smaller
!      than T
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130

!      Find an element in the first half of the array which is greater
!      than T
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140

!      Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF

!      Save upper and lower subscripts of the array yet to be sorted
      IF (L-I .GT. J-K) THEN
         IL(M) = I
         IU(M) = L
         I = K
         M = M+1
      ELSE
         IL(M) = K
         IU(M) = J
         J = L
         M = M+1
      ENDIF
      GO TO 160

!      Begin again on another portion of the unsorted array
  150 M = M-1
      IF (M .EQ. 0) GO TO 190
      I = IL(M)
      J = IU(M)

  160 IF (J-I .GE. 1) GO TO 120
      IF (I .EQ. 1) GO TO 110
      I = I-1

  170 I = I+1
      IF (I .EQ. J) GO TO 150
      T = X(I+1)
      TY = Y(I+1)
      IF (X(I) .LE. T) GO TO 170
      K = I

  180 X(K+1) = X(K)
      Y(K+1) = Y(K)
      K = K-1
      IF (T .LT. X(K)) GO TO 180
      X(K+1) = T
      Y(K+1) = TY
      GO TO 170

!      Clean up
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
      RETURN
end subroutine ssort



subroutine INTERP_SOURCE(NPH, NMU, T1, S1, T2, S2, T3, S3, T4, S4, T, S)
    ! Does either quadratic or cubic polynomial interpolation of the discrete
    ! ordinate source function as a function of optical depth.  The source
    ! function (S1,S2,S3,S4) is specified at three or four optical depths
    ! (T1,T2,T3,T4), interpolated to T, and output in S.
    implicit none
    integer, intent(IN) :: NPH, NMU
    real, intent(IN) :: T1, T2, T3, T4, T
    real, intent(IN) :: S1(NPH, NMU), S2(NPH, NMU), S3(NPH, NMU), S4(NPH, NMU)
    real, intent(OUT) :: S(NPH, NMU)

    if (abs(T1 - T2) < 4*spacing(T1)) then
        S = (T - T3)*(T - T4)*S2/((T2 - T3)*(T2 - T4)) &
            + (T - T2)*(T - T4)*S3/((T3 - T2)*(T3 - T4)) &
            + (T - T2)*(T - T3)*S4/((T4 - T2)*(T4 - T3))
    else if (abs(T3 - T4) < 4*spacing(T3)) then
        S = (T - T2)*(T - T3)*S1/((T1 - T2)*(T1 - T3)) &
            + (T - T1)*(T - T3)*S2/((T2 - T1)*(T2 - T3)) &
            + (T - T1)*(T - T2)*S3/((T3 - T1)*(T3 - T2))
    else
        S = (T - T2)*(T - T3)*(T - T4)*S1/((T1 - T2)*(T1 - T3)*(T1 - T4)) &
            + (T - T1)*(T - T3)*(T - T4)*S2/((T2 - T1)*(T2 - T3)*(T2 - T4)) &
            + (T - T1)*(T - T2)*(T - T4)*S3/((T3 - T1)*(T3 - T2)*(T3 - T4)) &
            + (T - T1)*(T - T2)*(T - T3)*S4/((T4 - T1)*(T4 - T2)*(T4 - T3))
    end if
end subroutine INTERP_SOURCE

subroutine VARIABLE_BRDF_SURFACE(NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                 MU2, PHI2, SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                 SFCTYPE, SFCPARMS, SFCPLANCK, &
                                 SFCRADIANCE, BOTRAD)
    ! Computes the upwelling radiance at the bottom boundary for one
    ! outgoing direction using the specified bidirectional reflectance
    ! distribution function.  The upwelling radiance includes the reflection of
    ! the incident radiance, the thermal emission (emissivity*Planck function),
    ! and the reflected direct solar flux (if applicable).  The upwelling
    ! radiance is the integral over all incident directions of the BRDF times
    ! the downwelling radiance, so a discrete sum is done and the integral
    ! weights (WTDO) are included.  The general BRDF function is called
    ! to compute the reflectance for the particular BRDF type (SFCTYPE)
    ! with parameters (SFCPARMS) for the incident and outgoing directions.
    ! The emissivity is computed implicitly from the integral of the BRDF.
    ! The outgoing direction is specified with (MU2,PHI2), and the BRDF is
    ! computed for all incident directions (loop over JMU,JPHI).  The incident
    ! downwelling radiances are input in SFCRADIANCE(:,:) and the outgoing
    ! upwelling radiance is output in BOTRAD.
    implicit none
    integer, intent(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
    real, intent(IN) :: MU(NMU), PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(IN) :: MU2, PHI2, SOLARMU, DIRFLUXSFC, SFCPLANCK
    real, intent(IN) :: SFCPARMS(*), SFCRADIANCE(NPHI0MAX, NMU)
    real, intent(OUT):: BOTRAD
    character(LEN=1), intent(IN) :: SRCTYPE, SFCTYPE
    integer :: JMU, JPHI
    real    :: OPI, SOLARAZ, REFLECT, W

    OPI = 1.0/acos(-1.0)
    SOLARAZ = 0.0

    ! Initialize the upwelling boundary radiances to zero (for thermal)
    !   or to the reflected direct solar flux.
    if (SRCTYPE == 'T') then
        BOTRAD = 0.0
    else
        call SURFACE_BRDF(SFCTYPE, SFCPARMS, MU2, PHI2, SOLARMU, SOLARAZ, &
                          REFLECT)
        BOTRAD = OPI*REFLECT*DIRFLUXSFC
    end if

    ! Integrate over the incident discrete ordinate directions (JMU,JPHI)
    do JMU = 1, NMU/2
        do JPHI = 1, NPHI0(JMU)
            call SURFACE_BRDF(SFCTYPE, SFCPARMS, MU2, PHI2, &
                              MU(JMU), PHI(JPHI, JMU), REFLECT)
            W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
            BOTRAD = BOTRAD + W*REFLECT*SFCRADIANCE(JPHI, JMU) &
                     + W*(1 - REFLECT)*SFCPLANCK
        end do
    end do

end subroutine VARIABLE_BRDF_SURFACE

subroutine SURFACE_BRDF(SFCTYPE, REFPARMS, &
                        MU2, PHI2, MU1, PHI1, REFLECT)
    !  Returns the reflection coefficient for the general bidirectional
    ! reflection distribution function of the specified type (SFCTYPE).
    ! The incident direction is (MU1,PHI1), and the outgoing direction
    ! is (MU2,PHI2) (MU is cosine of zenith angle, and PHI is the azimuthal
    ! angle in radians).  The incident directions have mu<0 (downward),
    ! while the outgoing directions have mu>0 (upward). The reflection
    ! function is normalized so that for a Lambertian surface (uniform
    ! reflection) the returned value (REFLECT) is simply the albedo.
    ! This routine calls the desired BRDF function, passing the
    ! appropriate parameters.  More BRDF surface types may be added easily
    ! by putting in the appropriate function calls.
    !        Type        Parameters
    !    R  Rahman et al  rho0, k, Theta
    implicit none
    real, intent(IN) :: REFPARMS(*), MU1, PHI1, MU2, PHI2
    real, intent(OUT):: REFLECT
    character(LEN=1), intent(IN) :: SFCTYPE
    real  :: PI
    real  :: RPV_REFLECTION

    PI = acos(-1.0)
    if (SFCTYPE == 'R') then
        ! R: Rahman, Pinty, and Verstraete
        REFLECT = RPV_REFLECTION(REFPARMS(1), REFPARMS(2), REFPARMS(3), &
                                 -MU1, MU2, PHI1 - PHI2 - PI)
    else
        stop 'SURFACE_BRDF: Unknown BRDF type'
    end if
end subroutine SURFACE_BRDF

function RPV_REFLECTION(RHO0, K, THETA, MU1, MU2, PHI)
    !  Computes the Rahman, Pinty, Verstraete BRDF.  The incident
    ! and outgoing cosine zenith angles are MU1 and MU2, respectively,
    ! and the relative azimuthal angle is PHI.  In this case the incident
    ! direction is where the radiation is coming from (i.e. opposite of the
    ! discrete ordinate), so MU1>0 and the hot spot is MU2=MU1 and PHI=0.
    ! The reference is:
    !    Rahman, Pinty, Verstraete, 1993: Coupled Surface-Atmosphere
    !    Reflectance (CSAR) Model. 2. Semiempirical Surface Model Usable
    !    With NOAA Advanced Very High Resolution Radiometer Data,
    !    J. Geophys. Res., 98, 20791-20801.
    implicit none
    real, intent(IN) :: RHO0, K, THETA, MU1, MU2, PHI
    real :: RPV_REFLECTION
    real M, F, H, COSPHI, SIN1, SIN2, COSG, TAN1, TAN2, CAPG

    M = MU1**(K - 1)*MU2**(K - 1)/(MU1 + MU2)**(1 - K)
    COSPHI = cos(PHI)
    SIN1 = sqrt(1.0 - MU1**2)
    SIN2 = sqrt(1.0 - MU2**2)
    COSG = MU1*MU2 + SIN1*SIN2*COSPHI
    F = (1 - THETA**2)/(1 + 2*THETA*COSG + THETA**2)**1.5
    TAN1 = SIN1/MU1
    TAN2 = SIN2/MU2
    CAPG = sqrt(TAN1**2 + TAN2**2 - 2*TAN1*TAN2*COSPHI)
    H = 1 + (1 - RHO0)/(1 + CAPG)
    RPV_REFLECTION = RHO0*M*F*H
end function RPV_REFLECTION

subroutine COMPUTE_RADIANCE(NLAY, ALBEDOSC, LEGENSC, &
                            MAXLEG, NLEGP, LEGENP, DELTAM, &
                            NPTS, HEIGHTG, TAUG, TAUGSC, IXP, IXG, &
                            ML, MM, NLM, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                            SRCTYPE, SOLARFLUX, SOLARMU, &
                            SKYRAD, WAVELEN, WAVENO, SFCTEMP, &
                            SFCTYPE, SFCPARMS, SOURCE, RADIANCE, &
                            NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                            RADOUT)
    ! Computes the radiances for the specified tau's, mu's, and phi's.
    implicit none
    integer, intent(IN) :: NLAY, MAXLEG, NLEGP(NLAY)
    integer, intent(IN) :: NPTS, IXP(2, NLAY), IXG(NPTS)
    integer, intent(IN) :: ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(IN) :: NZOUT, NMUOUT, NPHIOUT
    logical, intent(IN) :: DELTAM
    character(LEN=1), intent(IN) :: SRCTYPE, SFCTYPE
    real, intent(IN) :: ALBEDOSC(NLAY), LEGENSC(0:Ml, NLAY)
    real, intent(IN) :: LEGENP(0:MAXLEG, NLAY)
    real, intent(IN) :: HEIGHTG(NPTS), TAUG(NPTS), TAUGSC(NPTS)
    real, intent(IN) :: MU(NMU), PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(IN) :: SOLARFLUX, SOLARMU, SKYRAD
    real, intent(IN) :: SFCPARMS(*), SFCTEMP, WAVELEN, WAVENO(2)
    real, intent(IN) :: SOURCE(NLM, NPTS), RADIANCE(NPHI0MAX, NMU, NPTS)
    real, intent(INOUT) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
    real, intent(OUT) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)

    integer :: I, J, K, L
    real    :: ZDIF, RADDEG, PHIRAD
    real    :: DIRFLUXSFC, SFCPLANCK, BNDRAD, PLANCKFUNC
    integer, allocatable :: LZOUT(:)
    real, allocatable :: DIRFLUX(:), SOURCE1(:), RADVEC(:)

    RADDEG = acos(-1.0)/180.0
    allocate (DIRFLUX(NPTS), SOURCE1(NPTS), RADVEC(NPTS), LZOUT(NZOUT))

    ! Find the closest level in height to the output heights
    do I = 1, NZOUT
        ZDIF = 1.0e10
        do L = 1, NPTS
            if (abs(HEIGHTG(L) - ZOUT(I)) < ZDIF) then
                LZOUT(I) = L
                ZDIF = abs(HEIGHTG(L) - ZOUT(I))
            end if
        end do
        ZOUT(I) = HEIGHTG(LZOUT(I))
    end do

    do L = 1, NPTS
        ! Make the perpendicular direct beam flux for all the levels
        DIRFLUX(IXG(L)) = SOLARFLUX/abs(SOLARMU)*exp(-TAUGSC(L)/abs(SOLARMU))
    end do

    if (SFCTYPE /= 'L') then
        if (SRCTYPE /= 'T') then
            DIRFLUXSFC = SOLARFLUX*exp(-TAUGSC(IXP(2, NLAY))/abs(SOLARMU))
            SFCPLANCK = 0.0
        else
            SFCPLANCK = PLANCKFUNC(WAVELEN, WAVENO, SFCTEMP)
        end if
    end if

    ! Loop over the output radiance directions
    do J = 1, NMUOUT
        do K = 1, NPHIOUT
            PHIRAD = PHIOUT(K)*RADDEG
            if (MUOUT(J) == 0.0 .or. abs(MUOUT(J)) > 1.0) then
                write (*, *) 'COMPUTE_RADIANCE: Bad mu for radiance', MUOUT(J)
                RADOUT(:, J, K) = -999.0
                cycle
            end if

            ! Compute the source function throughout grid for this angle
            call COMPUTE_ONE_SOURCE(NLAY, ALBEDOSC, LEGENSC, &
                                    MAXLEG, NLEGP, LEGENP, &
                                    NPTS, IXP, IXG, ML, MM, NLM, &
                                    MUOUT(J), PHIRAD, SRCTYPE, DELTAM, &
                                    SOLARMU, DIRFLUX, SOURCE, SOURCE1)

            ! Get boundary radiances: either top or bottom
            !  Isotropic top boundary or Lambertian bottom boundary can use
            !  the previously computed boundary radiances in RADIANCE,
            !  otherwise, compute the radiance for this angle by integrating
            !  over the stored downwelling radiances.
            if (MUOUT(J) < 0.0) then
                BNDRAD = RADIANCE(1, 1, IXG(IXP(1, 1)))
            else
                if (SFCTYPE == 'L') then
                    BNDRAD = RADIANCE(1, NMU/2 + 1, IXG(IXP(2, NLAY)))
                else
                    call VARIABLE_BRDF_SURFACE(NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                               MUOUT(J), PHIRAD, &
                                               SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                               SFCTYPE, SFCPARMS, SFCPLANCK, &
                                               RADIANCE(:, :, IXG(IXP(2, NLAY))), BNDRAD)
                end if
            end if

            ! Integrate the source function for this direction to get the radiances
            !   at all the levels (RADVEC), then pull out the desired one.
            call INTEGRATE_SOURCE(NLAY, NPTS, IXP, IXG, TAUGSC, &
                                  MUOUT(J), BNDRAD, SOURCE1, RADVEC)
            do I = 1, NZOUT
                RADOUT(I, J, K) = RADVEC(LZOUT(I))
            end do
        end do
    end do

    deallocate (RADVEC, SOURCE1, DIRFLUX, LZOUT)
end subroutine COMPUTE_RADIANCE

subroutine COMPUTE_ONE_SOURCE(NLAY, ALBEDOSC, LEGENSC, &
                              MAXLEG, NLEGP, LEGENP, &
                              NPTS, IXP, IXG, ML, MM, NLM, &
                              MUO, PHIO, SRCTYPE, DELTAM, &
                              SOLARMU, DIRFLUX, SOURCE, SOURCE1)
    ! Computes the source function (SOURCE1) in the direction (MUO,PHIO)
    ! for all levels (NPTS).  The spherical harmonic source function
    ! series is input in SOURCE.  DIRFLUX is the attenuated direct beam
    ! perpendicular flux.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.
    implicit none
    integer, intent(IN) :: NLAY, MAXLEG, NLEGP(NLAY), NPTS
    integer, intent(IN) :: ML, MM, NLM, IXP(2, NLAY), IXG(NPTS)
    logical, intent(IN) :: DELTAM
    character(LEN=1), intent(IN) :: SRCTYPE
    real, intent(IN) :: ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(IN) :: LEGENP(0:MAXLEG, NLAY)
    real, intent(IN) :: MUO, PHIO, SOLARMU, DIRFLUX(NPTS)
    real, intent(IN) :: SOURCE(NLM, NPTS)
    real, intent(OUT) :: SOURCE1(NPTS)
    integer :: I, J, K, L, LAY, M, ME
    double precision  :: F, A, COSSCAT
    real, allocatable :: YLMDIR(:), YLMSUN(:)
    double precision, allocatable :: SUNDIRLEG(:), SINGSCAT(:)

    allocate (YLMDIR(NLM), YLMSUN(NLM), SUNDIRLEG(0:MAXLEG), SINGSCAT(NLAY))

    ! Precompute Ylm's for output direction
    call YLMALL(MUO, PHIO, ML, MM, 1, YLMDIR)

    ! Special case for solar source and Delta-M
    if (SRCTYPE /= 'T' .and. DELTAM) then
        ! Compute the Legendre polynomials for the scattering angle
        !   for the untruncated solar single scattering computation.
        COSSCAT = SOLARMU*MUO + sqrt((1.0d0 - SOLARMU**2)*(1.0d0 - MUO**2))*cos(PHIO)
        call LEGENDRE_ALL(COSSCAT, MAXLEG, SUNDIRLEG)
        SUNDIRLEG(:) = SUNDIRLEG(:)/(4.0d0*acos(-1.0d0))
        call YLMALL(SOLARMU, 0.0, ML, MM, 1, YLMSUN)

        ! Precompute the single scattering contribution for the original
        ! unscaled phase function (except for the direct beam transmission).
        ! Need to put the inverse of the tau scaling in the source function
        ! because tau is still scaled.
        do LAY = 1, NLAY
            SINGSCAT(LAY) = 0.0
            if (ML + 1 > NLEGP(LAY)) then
                F = 0
            else
                F = LEGENP(ML + 1, LAY)/(2*(ML + 1) + 1)
            end if
            do L = 0, NLEGP(LAY)
                SINGSCAT(LAY) = SINGSCAT(LAY) &
                                + SUNDIRLEG(L)*ALBEDOSC(LAY)*LEGENP(L, LAY)/(1 - F)
            end do
            if (SINGSCAT(LAY) < 0.0) then
                write (*, *) 'COMPUTE_ONE_SOURCE: negative source function for layer: ', LAY
                stop
            end if
        end do

        ! Loop over all the grid points, computing the source function
        !   at the viewing angle from the spherical harmonic source function.
        do LAY = 1, NLAY
            do I = IXP(1, LAY), IXP(2, LAY)
                K = IXG(I)
                ! Sum over the spherical harmonic series of the source function
                SOURCE1(K) = sum(SOURCE(:, K)*YLMDIR(:))
                ! First subtract off the truncated single scattering
                J = 1
                do L = 0, ML
                    ME = min(L, MM)
                    A = DIRFLUX(K)*ALBEDOSC(LAY)*LEGENSC(L, LAY)/(2*L + 1)
                    do M = 0, ME
                        SOURCE1(K) = SOURCE1(K) - A*YLMDIR(J)*YLMSUN(J)
                        J = J + 1
                    end do
                end do
                ! Then add in the single scattering contribution
                SOURCE1(K) = SOURCE1(K) + DIRFLUX(K)*SINGSCAT(LAY)
            end do
        end do

    else ! for thermal sources or no delta-M it is very simple:
        ! Sum over the spherical harmonic series of the source function
        do K = 1, NPTS
            SOURCE1(K) = sum(SOURCE(:, K)*YLMDIR(:))
        end do
    end if

    deallocate (YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT)
end subroutine COMPUTE_ONE_SOURCE

subroutine INTEGRATE_SOURCE(NLAY, NPTS, IXP, IXG, TAUGSC, &
                            MUO, BNDRAD, SOURCE1, RADVEC)
    ! Integrates the radiative transfer equation through the medium specified
    ! by the optical depth (TAUGSC) and source function (SOURCE1) in direction
    ! cosine MUO.  The radiance at the approprate boundary (top or bottom
    ! depending on the sign of MUO) is input (BNDRAD).  The radiance for
    ! the direction the source function was made at all the adaptive grid
    ! levels is output (RADVEC with indexing corresponding to TAUGSC).
    ! Uses the same quadratic/cubic polynomial source function interpolation
    ! as in PATH_INTEGRATION.
    implicit none
    integer, intent(IN) :: NLAY, NPTS, IXP(2, NLAY), IXG(NPTS)
    real, intent(IN) :: MUO, BNDRAD
    real, intent(IN) :: TAUGSC(NPTS), SOURCE1(NPTS)
    real, intent(OUT) :: RADVEC(NPTS)
    integer :: LAY, I, I1, I2
    real, allocatable :: SRCINT(:, :)
    double precision :: A(4), B(4), C(4), D(4), E, T, T1, T3, T4, TR

    ! Calculate and store the source function integration coefficients
    allocate (SRCINT(0:4, NPTS))
    do LAY = 1, NLAY
        I1 = IXP(1, LAY)
        I2 = IXP(2, LAY) - 1
        if (I1 == I2) then ! special case for layer with just two points
            T = (TAUGSC(I1 + 1) - TAUGSC(I1))/abs(MUO)
            TR = exp(-T) ! sublayer transmission along ordinate
            if (T < 1.0e-4) then
                D(1) = 0.5*T - (1.0d0/3)*T**2
                D(2) = 0.5*T - (1.0d0/6)*T**2
            else
                D(1) = (1 - TR*(1 + T))/T
                D(2) = (T - (1 - TR))/T
            end if
            SRCINT(0, I1) = TR
            if (MUO < 0.0) then
                SRCINT(1, I1) = D(2)
                SRCINT(2, I1) = D(1)
            else
                SRCINT(1, I1) = D(1)
                SRCINT(2, I1) = D(2)
            end if
        else
            do I = I1, I2
                A(:) = 0.0; B(:) = 0.0; C(:) = 0.0; D(:) = 0.0
                T3 = (TAUGSC(I + 1) - TAUGSC(I))/abs(MUO)
                if (I == I1) then
                    T4 = (TAUGSC(I + 2) - TAUGSC(I))/abs(MUO)
                    E = T3*T4; B(2) = 1/E; C(2) = -(T3 + T4)/E; D(2) = 1.0
                    E = T3*(T3 - T4); B(3) = 1/E; C(3) = -T4/E; D(3) = 0.0
                    E = T4*(T4 - T3); B(4) = 1/E; C(4) = -T3/E; D(4) = 0.0
                else if (I < I2) then
                    T1 = (TAUGSC(I - 1) - TAUGSC(I))/abs(MUO)
                    T4 = (TAUGSC(I + 2) - TAUGSC(I))/abs(MUO)
                    E = T1*(T1 - T3)*(T1 - T4); A(1) = 1/E
                    B(1) = -(T3 + T4)/E; C(1) = T3*T4/E; D(1) = 0.0
                    E = -T1*T3*T4; A(2) = 1/E
                    B(2) = -(T1 + T3 + T4)/E; C(2) = (T1*T3 + T1*T4 + T3*T4)/E; D(2) = 1.0
                    E = T3*(T3 - T1)*(T3 - T4); A(3) = 1/E
                    B(3) = -(T1 + T4)/E; C(3) = T1*T4/E; D(3) = 0.0
                    E = T4*(T4 - T1)*(T4 - T3); A(4) = 1/E
                    B(4) = -(T1 + T3)/E; C(4) = T1*T3/E; D(4) = 0.0
                else
                    T1 = (TAUGSC(I - 1) - TAUGSC(I))/abs(MUO)
                    E = T1*(T1 - T3); B(1) = 1/E; C(1) = -T3/E; D(1) = 0.0
                    E = T1*T3; B(2) = 1/E; C(2) = -(T1 + T3)/E; D(2) = 1.0
                    E = T3*(T3 - T1); B(3) = 1/E; C(3) = -T1/E; D(3) = 0.0
                end if
                T = T3
                TR = exp(-T)
                SRCINT(0, I) = TR
                if (MUO < 0.0) then
                    SRCINT(1:4, I) = A(:)*(6*TR + T**3 - 3*T**2 + 6*T - 6) + D(:)*(1 - TR) &
                                     + B(:)*(-2*TR + T**2 - 2*T + 2) + C(:)*(TR + T - 1)
                else
                    SRCINT(1:4, I) = A(:)*(6 - TR*(T**3 + 3*T**2 + 6*T + 6)) + D(:)*(1 - TR) &
                                     + B(:)*(2 - TR*(T**2 + 2*T + 2)) + C(:)*(1 - TR*(T + 1))
                end if
            end do
        end if
    end do

    if (MUO < 0.0) then
        ! For downwelling direction integrate source function from the top down
        RADVEC(IXP(1, 1)) = BNDRAD
        do LAY = 1, NLAY
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADVEC(I1 + 1) = SRCINT(0, I1)*RADVEC(I1) &
                                 + SRCINT(1, I1)*SOURCE1(IXG(I1)) + SRCINT(2, I1)*SOURCE1(IXG(I1 + 1))
            else
                RADVEC(I1 + 1) = SRCINT(0, I1)*RADVEC(I1) + SRCINT(2, I1)*SOURCE1(IXG(I1)) &
                                 + SRCINT(3, I1)*SOURCE1(IXG(I1 + 1)) + SRCINT(4, I1)*SOURCE1(IXG(I1 + 2))
                do I = I1 + 1, I2 - 1
                    RADVEC(I + 1) = SRCINT(0, I)*RADVEC(I) + SRCINT(1, I)*SOURCE1(IXG(I - 1)) &
                                    + SRCINT(2, I)*SOURCE1(IXG(I)) + SRCINT(3, I)*SOURCE1(IXG(I + 1)) &
                                    + SRCINT(4, I)*SOURCE1(IXG(I + 2))
                end do
                RADVEC(I2 + 1) = SRCINT(0, I2)*RADVEC(I2) + SRCINT(1, I2)*SOURCE1(IXG(I2 - 1)) &
                                 + SRCINT(2, I2)*SOURCE1(IXG(I2)) + SRCINT(3, I2)*SOURCE1(IXG(I2 + 1))
            end if
            if (LAY < NLAY) RADVEC(IXP(1, LAY + 1)) = RADVEC(IXP(2, LAY))
        end do
    else
        ! For upwelling direction integrate source function from the bottom up
        RADVEC(IXP(2, NLAY)) = BNDRAD
        do LAY = NLAY, 1, -1
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADVEC(I1) = SRCINT(0, I1)*RADVEC(I1 + 1) &
                             + SRCINT(1, I1)*SOURCE1(IXG(I1)) + SRCINT(2, I1)*SOURCE1(IXG(I1 + 1))
            else
                RADVEC(I2) = SRCINT(0, I2)*RADVEC(I2 + 1) + SRCINT(1, I2)*SOURCE1(IXG(I2 - 1)) &
                             + SRCINT(2, I2)*SOURCE1(IXG(I2)) + SRCINT(3, I2)*SOURCE1(IXG(I2 + 1))
                do I = I2 - 1, I1 + 1, -1
                    RADVEC(I) = SRCINT(0, I)*RADVEC(I + 1) + SRCINT(1, I)*SOURCE1(IXG(I - 1)) &
                                + SRCINT(2, I)*SOURCE1(IXG(I)) + SRCINT(3, I)*SOURCE1(IXG(I + 1)) &
                                + SRCINT(4, I)*SOURCE1(IXG(I + 2))
                end do
                RADVEC(I1) = SRCINT(0, I1)*RADVEC(I1 + 1) + SRCINT(2, I1)*SOURCE1(IXG(I1)) &
                             + SRCINT(3, I1)*SOURCE1(IXG(I1 + 1)) + SRCINT(4, I1)*SOURCE1(IXG(I1 + 2))
            end if
            if (LAY > 1) RADVEC(IXP(2, LAY - 1)) = RADVEC(IXP(1, LAY))
        end do
    end if
    deallocate (SRCINT)
end subroutine INTEGRATE_SOURCE

subroutine INIT_RADIANCE(NLAY, ML, TAUPSC, ALBEDOSC, LEGENSC, TEMPP, &
                         SRCTYPE, SOLARFLUX, SOLARMU, GNDALB, SFCTEMP, &
                         SKYRAD, WAVELEN, WAVENO, NPTS, TAUGSC, IXP, NLM, &
                         SHRADIANCE)
    ! Initializes the spherical harmonic radiance field by solving the
    ! plane-parallel Eddington system to give the first two terms of the
    ! spherical harmonic series at the NPTS sublayer interfaces.
    implicit none
    integer, intent(IN) :: NLAY, ML, NPTS, IXP(2, NLAY), NLM
    real, intent(IN) :: TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(IN) :: TEMPP(NLAY + 1)
    real, intent(IN) :: SOLARFLUX, SOLARMU, GNDALB, SFCTEMP
    real, intent(IN) :: SKYRAD, WAVELEN, WAVENO(2), TAUGSC(NPTS)
    character(LEN=1), intent(IN) :: SRCTYPE
    real, intent(OUT):: SHRADIANCE(NLM, NPTS)

    integer :: NLAYER, LAY, I, L
    logical, parameter :: DELTAM = .false.
    real    :: PI, C0, C1, GNDEMIS, TAU1, U
    real, allocatable :: OPTDEPTHS(:), ALBEDOS(:), ASYMMETRIES(:)
    real, allocatable :: TEMPS(:), FLUXES(:, :)
    character(LEN=1) :: SRCTYPE2

    NLAYER = sum(IXP(2, :) - IXP(1, :))
    allocate (OPTDEPTHS(NLAYER), ALBEDOS(NLAYER), ASYMMETRIES(NLAYER))
    allocate (TEMPS(NLAYER + 1), FLUXES(3, NLAYER + 1))

    SRCTYPE2 = SRCTYPE
    if (SRCTYPE == 'B') SRCTYPE2 = 'T'
    PI = acos(-1.0)
    C0 = sqrt(1.0/PI)
    C1 = sqrt(3.0/(4*PI))

    ! Make layer properties for the Eddington routine
    L = 0
    do LAY = 1, NLAY
        TAU1 = TAUGSC(IXP(1, LAY))
        do I = IXP(1, LAY), IXP(2, LAY) - 1
            L = L + 1
            OPTDEPTHS(L) = TAUGSC(I + 1) - TAUGSC(I)
            ALBEDOS(L) = ALBEDOSC(LAY)
            ASYMMETRIES(L) = LEGENSC(1, LAY)/3.0
            if (TAUPSC(LAY) > 0.0) then
                U = (TAUGSC(I) - TAU1)/TAUPSC(LAY)
            else
                U = (I - IXP(1, LAY))/FLOAT(IXP(2, LAY) - IXP(1, LAY))
            end if
            TEMPS(L) = TEMPP(LAY) + U*(TEMPP(LAY + 1) - TEMPP(LAY))
        end do
    end do
    if (L /= NLAYER) stop 'INIT_RADIANCE: initial grid error'
    TEMPS(NLAYER + 1) = TEMPP(NLAY + 1)

    ! Call the Eddington flux routine
    GNDEMIS = 1.0 - GNDALB
    call EDDRTF(NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, &
                TEMPS, DELTAM, SRCTYPE2, SOLARFLUX, SOLARMU, &
                SFCTEMP, GNDEMIS, SKYRAD, WAVELEN, WAVENO, FLUXES)

    ! Convert fluxes to first two moments of spherical harmonics
    SHRADIANCE(:, :) = 0.0
    L = 0
    do LAY = 1, NLAY
        do I = IXP(1, LAY), IXP(2, LAY)
            L = L + 1
            SHRADIANCE(1, I) = C0*(FLUXES(1, L) + FLUXES(2, L))
            SHRADIANCE(2, I) = C1*(FLUXES(1, L) - FLUXES(2, L))
        end do
        L = L - 1
    end do

    deallocate (OPTDEPTHS, ALBEDOS, ASYMMETRIES, TEMPS, FLUXES)
end subroutine INIT_RADIANCE

subroutine EDDRTF(NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, &
                  TEMPS, DELTAM, SRCTYPE, SOLARFLUX, SOLARMU, &
                  SFCTEMP, GNDEMIS, SKYRAD, WAVELEN, WAVENO, FLUXES)
    ! EDDRTF computes the layer interface fluxes for a plane-parallel
    ! atmosphere with either solar or thermal source of radiation using
    ! the Eddington approximation.  The medium is specified by a number
    ! of homogeneous layers.  For a thermal source the Planck function is
    ! linear with optical depth, while for a solar source it is exponential.
    ! The temperatures, optical depth, single scattering albedo, and
    ! asymmetry parameter are specified for each layer.  The boundary
    ! conditions such as the solar flux, and reflection and/or emission from
    ! ground surface are also specified. Delta Eddington scaling may be
    ! used.  The diffuse Eddington fluxes and the solar direct flux at
    ! each level are returned.
    ! The model works by calculating the reflection, transmission, and
    ! source terms for each layer from the input properties.  A
    ! tri-diagonal matrix solver is then used to compute the diffuse fluxes
    ! at each layer interface from the applied boundary conditions.
    !
    ! Parameters:
    !   Input:
    ! NLAYER         integer      Number of homogenous layers
    !                              (layers are specified from the top down)
    ! OPTDEPTHS      real array   Optical thickness of layers
    ! ALBEDOS        real array   Single scattering albedos
    ! ASYMMETRIES    real array   Asymmetry parameters
    ! TEMPS          real array   Temperatures (K) at layer interfaces
    !                              (e.g. TEMPS(1) is at top of top layer,
    !                               TEMPS(2) is at bottom of top layer).
    ! DELTAM         logical      True for delta-Eddington scaling
    ! SRCTYPE        character    'S' for solar source, 'T' for thermal source
    ! SOLARFLUX      real         Incident solar flux on horizonal plane
    ! SOLARMU        real         Cosine of the solar zenith angle
    ! SFCTEMP        real         Surface temperature (Kelvin)
    ! GNDEMIS        real         Ground emissivity (1-albedo)
    ! SKYRAD         real         Radiance (for solar) or brightness
    !                              temperature (for thermal) of isotropic
    !                              incident radiation from above
    ! WAVELEN        real         Wavelength (micron)
    ! WAVENO(2)      real         Wavenumber range (cm^-1) (if WAVELEN=0)
    !
    !   Output:
    ! FLUXES         real         Eddington fluxes at layer interfaces.
    !                               FLUXES(1,L) is upwelling diffuse,
    !                               FLUXES(2,L) is downwelling diffuse,
    !                               FLUXES(3,L) is downwelling direct,
    !                               L=1 is top, L=NLAYER+1 is bottom
    implicit none
    integer, intent(IN) :: NLAYER
    logical, intent(IN) :: DELTAM
    real, intent(IN) :: TEMPS(NLAYER + 1)
    real, intent(IN) :: OPTDEPTHS(NLAYER), ALBEDOS(NLAYER), ASYMMETRIES(NLAYER)
    real, intent(IN) :: SFCTEMP, GNDEMIS, SKYRAD, SOLARFLUX, SOLARMU
    real, intent(IN) :: WAVELEN, WAVENO(2)
    character(LEN=1), intent(IN) ::  SRCTYPE
    real, intent(OUT) :: FLUXES(3, NLAYER + 1)

    integer :: N, L, I
    double precision :: DELTAU, G, OMEGA, F
    double precision :: LAMBDA, R, T, D, CP, CM, A, B, X1, X2
    double precision :: REFLECT, TRANS, SOURCEP, SOURCEM
    double precision :: RADP1P, RADP1M, RADP2P, RADP2M
    double precision :: PI, MU0, SKYFLUX, GNDFLUX, PLANCK1, PLANCK2, C, TAU
    double precision :: EXLP, EXLM, V, DS, B1, B2, SOLPP, SOLPM
    double precision, allocatable :: LOWER(:), UPPER(:), DIAG(:), RHS(:)
    parameter(PI=3.1415926535)
    real :: PLANCKFUNC

    ! Compute the reflection, transmission, and source coefficients
    !   for each layer for the diffuse Eddington two stream problem.
    N = 2*NLAYER + 2
    allocate (LOWER(N), UPPER(N), DIAG(N), RHS(N))
    if (SRCTYPE == 'T') then
        PLANCK1 = PI*PLANCKFUNC(WAVELEN, WAVENO, TEMPS(1))
    end if
    MU0 = abs(SOLARMU)
    TAU = 0.0
    I = 2
    do L = 1, NLAYER
        DELTAU = OPTDEPTHS(L)
        if (DELTAU < 0.0) stop 'EDDRTF: TAU<0'
        ! Special case for zero optical depth
        if (DELTAU == 0.0) then
            TRANS = 1.0
            REFLECT = 0.0
            SOURCEP = 0.0
            SOURCEM = 0.0
        else
            OMEGA = ALBEDOS(L)
            G = ASYMMETRIES(L)
            if (DELTAM) then
                F = G**2
                DELTAU = (1 - OMEGA*F)*DELTAU
                OMEGA = (1 - F)*OMEGA/(1 - OMEGA*F)
                G = (G - F)/(1 - F)
            end if
            R = (1.0 - OMEGA*(4.0 - 3.0*G))/4.0
            T = (7.0 - OMEGA*(4.0 + 3.0*G))/4.0
            LAMBDA = sqrt(3.0*(1.0 - OMEGA)*(1.0 - OMEGA*G))
            ! Special case for conservative scattering (lambda=0)
            if (LAMBDA == 0.0) then
                D = 1.0/(1.0 + T*DELTAU)
                TRANS = D
                REFLECT = -R*DELTAU*D
            else
                X1 = -R
                X2 = LAMBDA + T
                EXLP = DEXP(min(LAMBDA*DELTAU, 75.d0))
                EXLM = 1.0/EXLP
                TRANS = 2.*LAMBDA/(X2*EXLP + (LAMBDA - T)*EXLM)
                REFLECT = X1*(EXLP - EXLM)*TRANS/(2.*LAMBDA)
                D = 1.0/(X2**2*EXLP - X1**2*EXLM)
            end if

            if (SRCTYPE == 'T') then
                ! Calculate thermal source terms
                PLANCK2 = PI*PLANCKFUNC(WAVELEN, WAVENO, TEMPS(L + 1))
                V = 2.0*(PLANCK2 - PLANCK1)/(3.0*(1.-OMEGA*G)*DELTAU)
                RADP1P = -V + PLANCK1
                RADP2M = V + PLANCK2
                RADP2P = -V + PLANCK2
                RADP1M = V + PLANCK1
                if (LAMBDA .eq. 0.0) then
                    A = (R*DELTAU*RADP1P - RADP2M)*D
                    B = -(R*RADP1P + T*RADP2M)*D
                    SOURCEP = (B - T*(A + B*DELTAU))/R + RADP2P
                    SOURCEM = A + RADP1M
                else
                    CP = (X1*EXLM*RADP1P - X2*RADP2M)*D
                    CM = (-X2*EXLP*RADP1P + X1*RADP2M)*D
                    SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
                    SOURCEM = X2*CP + X1*CM + RADP1M
                end if
                PLANCK1 = PLANCK2
                FLUXES(3, L) = 0.0
            else
                ! Calculate solar source terms
                FLUXES(3, L) = SOLARFLUX*exp(-TAU/MU0)
                DS = 1.0/(LAMBDA**2 - 1.0/MU0**2)
                B1 = 0.5*OMEGA*(SOLARFLUX/MU0)*exp(-TAU/MU0)*DS
                B2 = 0.5*OMEGA*(SOLARFLUX/MU0)*exp(-(TAU + DELTAU)/MU0)*DS
                SOLPP = 1.0 + 1.5*G*MU0
                SOLPM = -1.0 + 1.5*G*MU0
                RADP1P = ((T + 1.0/MU0)*SOLPP + R*SOLPM)*B1
                RADP2M = ((-T + 1.0/MU0)*SOLPM - R*SOLPP)*B2
                RADP2P = ((T + 1.0/MU0)*SOLPP + R*SOLPM)*B2
                RADP1M = ((-T + 1.0/MU0)*SOLPM - R*SOLPP)*B1
                if (LAMBDA .eq. 0.0) then
                    A = (R*DELTAU*RADP1P - RADP2M)*D
                    B = -(R*RADP1P + T*RADP2M)*D
                    SOURCEP = (B - T*(A + B*DELTAU))/R + RADP2P
                    SOURCEM = A + RADP1M
                else
                    CP = (X1*EXLM*RADP1P - X2*RADP2M)*D
                    CM = (-X2*EXLP*RADP1P + X1*RADP2M)*D
                    SOURCEP = X1*CP*EXLP + X2*CM*EXLM + RADP2P
                    SOURCEM = X2*CP + X1*CM + RADP1M
                end if
                TAU = TAU + DELTAU
            end if
        end if
        DIAG(I) = -REFLECT
        DIAG(I + 1) = -REFLECT
        LOWER(I) = 1.0
        LOWER(I + 1) = -TRANS
        UPPER(I) = -TRANS
        UPPER(I + 1) = 1.0
        RHS(I) = SOURCEM
        RHS(I + 1) = SOURCEP
        I = I + 2
    end do

    ! Set up boundary radiances
    if (SRCTYPE == 'S') then
        FLUXES(3, NLAYER + 1) = SOLARFLUX*exp(-TAU/MU0)
        GNDFLUX = (1.0 - GNDEMIS)*SOLARFLUX*exp(-TAU/MU0)
        SKYFLUX = PI*SKYRAD
    else
        FLUXES(3, NLAYER + 1) = 0.0
        GNDFLUX = PI*GNDEMIS*PLANCKFUNC(WAVELEN, WAVENO, SFCTEMP)
        SKYFLUX = PI*PLANCKFUNC(WAVELEN, WAVENO, SKYRAD)
    end if
    ! Setup for and call the tri-diagonal matrix solver
    RHS(1) = SKYFLUX
    DIAG(1) = 0.0
    UPPER(1) = 1.0
    DIAG(N) = -(1.0 - GNDEMIS)
    LOWER(N) = 1.0
    RHS(N) = GNDFLUX
    call TRIDIAG(N, LOWER, DIAG, UPPER, RHS)
    ! Put the fluxes in the output array
    I = 1
    do L = 1, NLAYER + 1
        FLUXES(1, L) = RHS(I)
        FLUXES(2, L) = RHS(I + 1)
        I = I + 2
    end do

    deallocate (LOWER, UPPER, DIAG, RHS)
end subroutine EDDRTF

subroutine TRIDIAG(N, LOWER, DIAG, UPPER, RHS)
    ! Computes the solution to a tridiagonal system.
    ! N is order of the matrix.  LOWER(2..N) is the subdiagonal,
    ! DIAG(1..N) is the diagonal, and UPPER(1..N-1) is the
    ! superdiagonal.  On input RHS is the right hand side, while
    ! on output it is the solution vector.  Everything is destroyed.
    ! Hacked from Linpack DGTSL.
    implicit none
    integer :: N
    double precision :: LOWER(*), DIAG(*), UPPER(*), RHS(*)
    integer :: K, KB
    double precision :: T

    if (N .eq. 1) then
        if (DIAG(1) .eq. 0.0) goto 990
        RHS(1) = RHS(1)/DIAG(1)
    end if
    LOWER(1) = DIAG(1)
    DIAG(1) = UPPER(1)
    UPPER(1) = 0.0
    UPPER(N) = 0.0
    do K = 1, N - 1
        ! Interchange this and next row to the get the largest pivot.
        if (abs(LOWER(K + 1)) .ge. abs(LOWER(K))) then
            T = LOWER(K + 1)
            LOWER(K + 1) = LOWER(K)
            LOWER(K) = T
            T = DIAG(K + 1)
            DIAG(K + 1) = DIAG(K)
            DIAG(K) = T
            T = UPPER(K + 1)
            UPPER(K + 1) = UPPER(K)
            UPPER(K) = T
            T = RHS(K + 1)
            RHS(K + 1) = RHS(K)
            RHS(K) = T
        end if
        if (LOWER(K) .eq. 0.0) goto 990
        T = -LOWER(K + 1)/LOWER(K)
        LOWER(K + 1) = DIAG(K + 1) + T*DIAG(K)
        DIAG(K + 1) = UPPER(K + 1) + T*UPPER(K)
        UPPER(K + 1) = 0.0
        RHS(K + 1) = RHS(K + 1) + T*RHS(K)
    end do
    if (LOWER(N) .eq. 0.0) goto 990
    ! Back substitute
    RHS(N) = RHS(N)/LOWER(N)
    RHS(N - 1) = (RHS(N - 1) - DIAG(N - 1)*RHS(N))/LOWER(N - 1)
    do KB = 1, N - 2
        K = N - 2 - KB + 1
        RHS(K) = (RHS(K) - DIAG(K)*RHS(K + 1) - UPPER(K)*RHS(K + 2))/LOWER(K)
    end do
    return
990 continue
    stop 'Singular matrix in TRIDIAG'
end subroutine TRIDIAG

subroutine LEGENDRE_ALL(COSSCAT, NLEG, P)
    ! This subroutine computes a set of Legendre polynomials for
    ! a particular scattering angle COSSCAT.  NLEG is the maximum term.
    ! The Legendre functions evaluated at COSSCAT are returned in
    ! P, starting at l=0 and ending with l=NLEG  (NLEG+1 terms).
    implicit none
    integer NLEG
    double precision COSSCAT, P(0:NLEG)
    integer L
    double precision X, PL, PL1, PL2

    X = dble(COSSCAT)
    if (X*X .gt. 1.) stop 'LEGENDRE_ALL: |COSSCAT| larger than 1'
    ! Use the stable upward recursion on l, starting from P_0
    PL2 = 1.0d0
    P(0) = PL2
    if (NLEG .gt. 1) then
        PL1 = X
        P(1) = X
    end if
    do L = 2, NLEG
        PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        P(L) = PL
        PL2 = PL1
        PL1 = PL
    end do
end subroutine LEGENDRE_ALL

subroutine YLMALL(MU, PHI, ML, MM, NCS, P)
    ! This subroutine computes a set of normalized spherical harmonic
    ! functions, P(J), for a particular direction mu,phi.
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mode,
    ! and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for
    ! sines and cosines).  Returns normalized associated Legendre functions
    ! only if NCS<0. The set is returned for triangular truncation:
    ! J = NCS*(L*(L+1))/2 + M+1  for L<=MM
    ! J = (NCS*MM+1)*L-MM*(2+NCS*(MM-1))/2 + M+1  for L>MM
    implicit none
    integer ML, MM, NCS
    real MU, PHI, P(*)
    integer J, L, M, C
    double precision X, Y, A, PMM, PL, PL1, PL2, PHI8

    C = abs(NCS)
    if (C .ne. 1 .and. C .ne. 2) stop 'YLMALL: bad NCS'
    if (MM .gt. ML) stop 'YLMALL: MM greater than LM'
    if (MU*MU .gt. 1.) stop 'YLMALL: |MU| larger than 1'
    X = dble(MU)
    Y = sqrt(1.0d0 - X*X)
    ! Use the stable upward recursion on l, starting from P^m_m
    ! Put in the spherical harmonic normalization as it goes
    PMM = 1.0d0/sqrt(4.0d0*acos(-1.0d0))
    do M = 0, MM
        if (M .gt. 0) PMM = -PMM*Y*sqrt((2*M + 1.0d0)/(2.0d0*M))
        J = C*(M*(M + 1))/2 + M + 1
        P(J) = PMM
        PL2 = PMM
        if (M .lt. ML) then
            if (M + 1 .le. MM) J = C*((M + 1)*(M + 2))/2 + M + 1
            if (M + 1 .gt. MM) J = (C*MM + 1)*(M + 1) - (MM*(2 + C*(MM - 1)))/2 + M + 1
            PL1 = sqrt(2*M + 3.0d0)*X*PMM
            P(J) = PL1
        end if
        do L = M + 1, ML - 1
            if (L + 1 .le. MM) J = C*((L + 1)*(L + 2))/2 + M + 1
            if (L + 1 .gt. MM) J = (C*MM + 1)*(L + 1) - (MM*(2 + C*(MM - 1)))/2 + M + 1
            A = 1.0d0/((L + M + 1.d0)*(L - M + 1.d0))
            PL = sqrt((2*L + 1)*A*(2*L + 3))*X*PL1 &
                 - sqrt((2*L + 3)*A*(L + M)*(L - M)/(2*L - 1.))*PL2
            P(J) = PL
            PL2 = PL1
            PL1 = PL
        end do
        if (M .eq. 0) PMM = PMM*sqrt(2.0d0)
    end do
    ! If there are M<0 terms then fill them in
    if (C .eq. 2) then
        do L = 0, ML
            do M = 1, min(L, MM)
                if (L .le. MM) J = L*(L + 1) + M + 1
                if (L .gt. MM) J = MM*(2*L - MM) + L + M + 1
                P(J - 2*M) = P(J)
            end do
        end do
    end if
    ! Put in the azimuthal dependence
    PHI8 = dble(PHI)
    if (NCS .gt. 0) then
        do M = (1 - NCS)*MM, MM
            if (M .lt. 0) then
                A = sin(-M*PHI8)
            else if (M .gt. 0) then
                A = cos(M*PHI8)
            else
                A = 1.0d0
            end if
            do L = abs(M), ML
                if (L .le. MM) then
                    J = C*(L*(L + 1))/2 + M + 1
                else
                    J = (C*MM + 1)*L - (MM*(2 + C*(MM - 1)))/2 + M + 1
                end if
                P(J) = A*P(J)
            end do
        end do
    end if
end subroutine YLMALL

subroutine GAUSQUADS(N, XA, WT)
    ! Generates the abscissas (X) and weights (W) for an N point
    ! Gauss-Legendre quadrature.  The XA are returned in this order:
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    implicit none
    integer N
    real XA(*), WT(*)
    integer K, I, J, L
    double precision X, XP, PL, PL1, PL2, DPL, TINY
    parameter(TINY=3.0d-7)

    K = N/2
    if (2*K .ne. N) stop 'GAUSQUADS: N must be even'
    do J = 1, K
        X = cos(3.141592654*(J - .25)/(N + .5))
        I = 0
100     continue
        PL1 = 1
        PL = X
        do L = 2, N
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        end do
        DPL = N*(X*PL - PL1)/(X*X - 1)
        XP = X
        X = XP - PL/DPL
        I = I + 1
        if (abs(X - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J) = -X
        XA(J + K) = X
        WT(J) = 2.0d0/((1.0d0 - X*X)*DPL*DPL)
        WT(J + K) = WT(J)
    end do
end subroutine GAUSQUADS

subroutine DGAUSQUADS(N, XA, WT)
    ! Generates the abscissas (X) and weights (W) for an N point
    ! Double-Gauss-Legendre quadrature.  The XA are returned in this order:
    ! -mu1, -mu2, ..., -muK, mu1, mu2, ..., muK  (mu1 > mu2 > muK,
    ! K=N/2, N must be even).
    implicit none
    integer N
    real XA(*), WT(*)
    integer K, I, J, L, N2
    double precision X, XP, PL, PL1, PL2, DPL, TINY
    parameter(TINY=3.0d-7)

    N2 = N/2
    if (2*N2 .ne. N) stop 'DGAUSQUADS: N must be even'
    K = (N2 + 1)/2
    do J = 1, K
        X = cos(3.141592654*(J - .25)/(N2 + .5))
        I = 0
100     continue
        PL1 = 1
        PL = X
        do L = 2, N2
            PL2 = PL1
            PL1 = PL
            PL = ((2*L - 1)*X*PL1 - (L - 1)*PL2)/L
        end do
        DPL = N2*(X*PL - PL1)/(X*X - 1)
        XP = X
        X = XP - PL/DPL
        I = I + 1
        if (abs(X - XP) .gt. TINY .and. I .lt. 10) GO TO 100
        XA(J + N2) = (1 + X)/2
        XA(N2 + 1 - J + N2) = (1 - X)/2
        WT(J + N2) = 1/((1 - X*X)*DPL*DPL)
        WT(N2 + 1 - J + N2) = WT(J + N2)
        XA(J) = -XA(J + N2)
        XA(N2 + 1 - J) = -XA(N2 + 1 - J + N2)
        WT(J) = WT(J + N2)
        WT(N2 + 1 - J) = WT(J + N2)
    end do
end subroutine DGAUSQUADS

subroutine INTEGRATE_PLANCK(X1, X2, F)
    ! Returns integral of x**3/(exp(x)-1) from x1 to x2.  Accurate to better
    ! than 1 part in 10**9 for integral from 0 to infinity (pi**4/15).
    implicit none
    double precision :: X1, X2, F
    double precision :: INTLOW, INTHIGH, C

    C = 1.0d0
    if (X1 < C .and. X2 < C) then
        F = INTLOW(X2) - INTLOW(X1)
    else if (X1 < C .and. X2 >= C) then
        F = INTLOW(C) - INTLOW(X1) + INTHIGH(C) - INTHIGH(X2)
    else if (X1 >= C .and. X2 >= C) then
        F = INTHIGH(X1) - INTHIGH(X2)
    else
        stop 'X1 and X2 out of order'
    end if
end subroutine INTEGRATE_PLANCK

function INTLOW(X)
    ! Integral of x**3/(exp(x)-1) from 0 to x.  Accurate for x less than about 1.
    ! Uses Taylor series expansion around x=0.
    implicit none
    double precision :: X, INTLOW
    integer :: N
    double precision :: SUM, F, A(29)
    data A/0.0d0, 0.0d0, 0.333333333333333333333d0, &
        -0.125d0, 0.016666666666666666667d0, 0.0d0, -0.0001984126984126984127d0, &
        0.0d0, 0.36743092298647854203d-5, 0.0d0, -0.75156325156325156325d-7, &
        0.0d0, 0.16059043836821614599d-8, 0.0d0, -0.35227934257916621232d-10, &
        0.0d0, 0.78720803121674581370d-12, 0.0d0, -0.17840422612224120352d-13, &
        0.0d0, 0.40886009791799259829d-15, 0.0d0, -0.94559508632959211871d-17, &
        0.0d0, 0.22036011313440918061d-18, 0.0d0, -0.51683202540046382743d-20, &
        0.0d0, 0.12188644964239543006d-20/

    SUM = A(4)*X**4
    F = X
    do N = 3, 29, 2
        F = X*X*F
        SUM = SUM + A(N)*F
    end do
    INTLOW = SUM
end function INTLOW

function INTHIGH(X)
    ! Integral of x**3/(exp(x)-1) from x to infinity.
    ! Accurate for x greater than about 1.
    implicit none
    double precision :: X, INTHIGH
    integer :: N
    double precision :: SUM

    SUM = 0.0d0
    do N = 1, 15
        SUM = SUM + exp(-N*X)*(X**3/N + 3*X**2/N**2 + 6*X/N**3 + 6.0d0/N**4)
    end do
    INTHIGH = SUM
end function INTHIGH

subroutine ADD_IN_KDIST(NLAY, HEIGHTP, TAUP, ALBEDOP, &
                        NZCKD, ZCKD, GASABS, TAUPM, ALBEDOPM)
    ! Adds the gaseous absorption in with the particle optical properties.
    ! The gaseous absorption optical depth is found for each of the NLAY
    ! optical property layers by interpolated the absorption linearly
    ! between the levels in ZCKD.
    integer, intent(IN) :: NLAY, NZCKD
    real, intent(IN) :: HEIGHTP(NLAY + 1), TAUP(NLAY), ALBEDOP(NLAY)
    real, intent(IN) :: ZCKD(NZCKD), GASABS(NZCKD)
    real, intent(OUT) :: TAUPM(NLAY), ALBEDOPM(NLAY)
    integer :: LAY, I
    real    :: GASTAU, Z, Z1, Z2, F, ABS

    if (NZCKD <= 0) then
        TAUPM(1:NLAY) = TAUP(1:NLAY)
        ALBEDOPM(1:NLAY) = ALBEDOP(1:NLAY)
    else
        if (HEIGHTP(1) > ZCKD(1) .or. HEIGHTP(NLAY + 1) < ZCKD(NZCKD)) then
            write (*, *) 'Layer heights outside of CKD file height range.'
            stop
        end if
        I = 1
        do LAY = 1, NLAY
            ! First, use trapezoidal integration to get optical depth of
            ! gaseous absorption assuming linear interpolation.
            GASTAU = 0.0
            Z1 = HEIGHTP(LAY)
            Z2 = HEIGHTP(LAY + 1)
            do while (ZCKD(I) >= Z1)
                I = I + 1
            end do
            I = I - 1
            Z = Z1
            do while (ZCKD(I) > Z2)
                F = (Z - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
                ABS = (1 - F)*GASABS(I) + F*GASABS(I + 1)
                GASTAU = GASTAU + 0.5*(ABS + GASABS(I + 1))*(Z - ZCKD(I + 1))
                Z = ZCKD(I + 1)
                I = I + 1
            end do
            I = I - 1
            F = (Z2 - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
            ABS = (1 - F)*GASABS(I) + F*GASABS(I + 1)
            GASTAU = GASTAU - 0.5*(ABS + GASABS(I + 1))*(Z2 - ZCKD(I + 1))
            ! Then, add the gaseous optical depth to the particle properties
            TAUPM(LAY) = TAUP(LAY) + GASTAU
            if (TAUPM(LAY) > 0.0) then
                ALBEDOPM(LAY) = ALBEDOP(LAY)*TAUP(LAY)/TAUPM(LAY)
            else
                ALBEDOPM(LAY) = 0.0
            end if
            ! print *, 'add_kdist: ',LAY,GASTAU,TAUPM(LAY),ALBEDOPM(LAY)
        end do
    end if
end subroutine ADD_IN_KDIST

subroutine SUM_KDIST_OUTPUT(IG, WT, NLAY, TAUP, IXP, NPTS, TAUG, &
                            HEIGHTP, HEIGHTG, FLUXUPM, FLUXUP, &
                            FLUXDNM, FLUXDN, MEANRADM, MEANRAD, &
                            NZOUT, NMUOUT, NPHIOUT, RADOUTM, RADOUT)
    ! Sums the output arrays over the k-distribution.  Only sums the output
    ! at the input layer interfaces because the adaptive grid changes for
    ! each k.  Also fixes up the tau grid and height grid accordingly.
    integer, intent(IN) :: IG, NLAY, NZOUT, NMUOUT, NPHIOUT
    integer, intent(INOUT) :: NPTS, IXP(2, NLAY)
    real, intent(IN) :: WT, TAUP(NLAY), HEIGHTP(NLAY + 1)
    real, intent(INOUT) :: TAUG(*), HEIGHTG(*)
    real, intent(IN) :: FLUXUPM(*), FLUXDNM(*), MEANRADM(*)
    real, intent(IN) :: RADOUTM(NZOUT, NMUOUT, NPHIOUT)
    real, intent(OUT) :: FLUXUP(*), FLUXDN(*), MEANRAD(*)
    real, intent(OUT) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)

    if (IG == 1) then
        FLUXUP(1:NLAY + 1) = 0.0
        FLUXDN(1:NLAY + 1) = 0.0
        MEANRAD(1:NLAY + 1) = 0.0
        RADOUT(:, :, :) = 0.0
    end if

    do L = 1, NLAY
        FLUXUP(L) = FLUXUP(L) + WT*FLUXUPM(IXP(1, L))
        FLUXDN(L) = FLUXDN(L) + WT*FLUXDNM(IXP(1, L))
        MEANRAD(L) = MEANRAD(L) + WT*MEANRADM(IXP(1, L))
    end do
    FLUXUP(NLAY + 1) = FLUXUP(NLAY + 1) + WT*FLUXUPM(IXP(2, NLAY))
    FLUXDN(NLAY + 1) = FLUXDN(NLAY + 1) + WT*FLUXDNM(IXP(2, NLAY))
    MEANRAD(NLAY + 1) = MEANRAD(NLAY + 1) + WT*MEANRADM(IXP(2, NLAY))
    RADOUT(:, :, :) = RADOUT(:, :, :) + WT*RADOUTM(:, :, :)

    HEIGHTG(1:NLAY + 1) = HEIGHTP(1:NLAY + 1)
    TAUG(1) = 0.0
    do L = 1, NLAY
        TAUG(L + 1) = TAUG(L) + TAUP(L)
        IXP(1, L) = L
        IXP(2, L) = L + 1
    end do
    NPTS = NLAY + 1
end subroutine SUM_KDIST_OUTPUT

subroutine USER_INPUT(PROPFILE, CKDFILE, KDIST, NMU, NPHI, &
                      SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                      SFCTEMP, SFCTYPE, SFCPARMS, WAVELEN, WAVENO, &
                      ACCELFLAG, SOLACC, MAXITER, SPLITACC, &
                      FLUXFILE, FLUXPARM, RADFILE, MAXPAR, NZOUT, ZOUT, &
                      NMUOUT, MUOUT, NPHIOUT, PHIOUT)
    !   Obtains the input parameters for the program by writing prompts
    ! to and reading responses from the terminal.  Echos the inputs to
    ! make a useful log file when running non-interactively.  See the overall
    ! program documentation for the list of input parameters.
    implicit none
    integer, intent(IN)  :: MAXPAR
    integer, intent(OUT) :: NMU, NPHI
    integer, intent(OUT) :: MAXITER, FLUXPARM, NZOUT, NMUOUT, NPHIOUT
    logical, intent(OUT) :: KDIST, ACCELFLAG
    real, intent(OUT) :: SOLARFLUX, SOLARMU
    real, intent(OUT) :: SFCTEMP, SKYRAD, WAVELEN, WAVENO(2)
    real, intent(OUT) :: SOLACC, SPLITACC
    real, intent(OUT) :: SFCPARMS(MAXPAR)
    real, intent(OUT) :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
    character(LEN=1), intent(OUT) :: SRCTYPE, SFCTYPE
    character(LEN=80), intent(OUT) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE

    write (*, *)
    write (*, '(1X,A)') 'Atmospheric optical properties file name'
    read (*, '(A)') PROPFILE
    write (*, *) PROPFILE

    write (*, '(1X,A)') 'Correlated k-distribution file name (or NONE)'
    read (*, '(A)') CKDFILE
    write (*, *) CKDFILE
    KDIST = CKDFILE(1:2) .ne. 'NO'

    write (*, '(1X,A)') 'Number of discrete ordinates in mu and phi'
    read (*, *) NMU, NPHI
    write (*, *) NMU, NPHI

    write (*, '(1X,A)') 'Thermal, solar, or both source (T, S, B)'
    read (*, '(A)') SRCTYPE
    write (*, *) SRCTYPE
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        write (*, '(1X,A)') 'Solar flux and direction (F, mu0)'
        read (*, *) SOLARFLUX, SOLARMU
        write (*, *) SOLARFLUX, SOLARMU
        write (*, '(1X,A)') 'Isotropic sky radiance'
        read (*, *) SKYRAD
        write (*, *) SKYRAD
        SOLARMU = -abs(SOLARMU)
    end if
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        write (*, '(1X,A)') 'Surface temperature'
        read (*, *) SFCTEMP
        write (*, *) SFCTEMP
    end if
    if (SRCTYPE .eq. 'T') then
        write (*, '(1X,A)') 'Sky temperature'
        read (*, *) SKYRAD
        write (*, *) SKYRAD
    end if

    write (*, '(1X,A)') 'Surface type (L-Lambertian, R-RPV vegetation BRDF)'
    read (*, '(A1)') SFCTYPE
    write (*, *) SFCTYPE
    if (SFCTYPE == 'L') then
        write (*, '(1X,A)') 'Surface albedo'
        read (*, *) SFCPARMS(1)
        write (*, *) SFCPARMS(1)
    else if (SFCTYPE == 'R') then
        write (*, '(1X,A)') 'RPV parameters (rho0, k, Theta)'
        read (*, *) SFCPARMS(1:3)
        write (*, *) SFCPARMS(1:3)
    else
        write (*, *) 'Unknown surface type'
        stop
    end if

    write (*, '(1X,A)') 'Wavelength (micron) (<=0 for wavenumber range)'
    read (*, *) WAVELEN
    write (*, *) WAVELEN
    if (KDIST .and. WAVELEN > 0.0) then
        write (*, *) 'Must use wavenumber range for k-distribution'
        stop
    end if
    if (WAVELEN <= 0.0) then
        write (*, '(1X,A)') 'Wavenumber range (cm^-1)'
        read (*, *) WAVENO(1:2)
        write (*, *) WAVENO(1:2)
    end if

    write (*, '(1X,A)') 'Layer splitting accuracy'
    read (*, *) SPLITACC
    write (*, *) SPLITACC

    write (*, '(1X,A,A)') 'Solution accuracy'
    read (*, *) SOLACC
    write (*, *) SOLACC

    write (*, '(1X,A,A)') 'Acceleration flag (T or F) and max iterations'
    read (*, *) ACCELFLAG, MAXITER
    write (*, *) ACCELFLAG, MAXITER

    write (*, '(1X,A,A)') 'Flux output file name'
    read (*, '(A)') FLUXFILE
    write (*, '(A)') FLUXFILE
    write (*, '(1X,A)') 'Fluxes output at input layer boundaries or adaptive levels (1,2)'
    read (*, *) FLUXPARM
    write (*, *) FLUXPARM

    write (*, '(1X,A,A)') 'Radiance output file name'
    read (*, '(A)') RADFILE
    write (*, '(A)') RADFILE

    write (*, *) 'Number of heights for output radiance'
    read (*, *) NZOUT
    write (*, *) NZOUT
    write (*, *) 'Height levels for output radiance'
    read (*, *) ZOUT(1:NZOUT)
    write (*, *) ZOUT(1:NZOUT)

    write (*, *) 'Number of cosine zenith angles for output radiance'
    read (*, *) NMUOUT
    write (*, *) NMUOUT
    write (*, *) 'Cosine zenith angles for output radiance'
    read (*, *) MUOUT(1:NMUOUT)
    write (*, *) MUOUT(1:NMUOUT)

    write (*, *) 'Number of azimuth angles for output radiance'
    read (*, *) NPHIOUT
    write (*, *) NPHIOUT
    write (*, *) 'Azimuth angles for output radiance (degrees)'
    read (*, *) PHIOUT(1:NPHIOUT)
    write (*, *) PHIOUT(1:NPHIOUT)

    write (*, *)
end subroutine USER_INPUT

subroutine READ_PROPERTY_SIZE(PROPFILE, NLAY, MAXLEG)
    !   Reads parts of the property file to get the maximum array sizes
    ! needed for allocatable arrays.
    implicit none
    integer, intent(OUT) :: NLAY, MAXLEG
    character(LEN=*), intent(IN) :: PROPFILE
    integer :: I, NLEG
    real    :: TAU, Z, TTOP, TBOT, SSALB

    ! Open the file, get the number of layers and the max number Legendre terms
    open (UNIT=1, FILE=PROPFILE, STATUS='OLD')
    read (1, *) NLAY
    read (1, *) Z, TTOP
    MAXLEG = 1
    do I = 1, NLAY
        read (1, *) Z, TBOT, TAU, SSALB, NLEG
        MAXLEG = max(MAXLEG, NLEG)
    end do
    close (1)
end subroutine READ_PROPERTY_SIZE

subroutine READ_PROPERTIES(PROPFILE, NLAY, MAXLEG, &
                           HEIGHTP, TEMPP, TAUP, ALBEDOP, NLEGP, LEGENP)
    ! Reads the plane-parallel medium properties from the file into the
    ! property arrays.  The optical properties are uniform within each layer.
    ! The heights must decrease (top of atmosphere downward).
    implicit none
    integer, intent(IN)  :: NLAY, MAXLEG
    integer, intent(OUT) :: NLEGP(NLAY)
    real, intent(OUT) :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
    real, intent(OUT) :: HEIGHTP(NLAY + 1), TEMPP(NLAY + 1)
    character(LEN=*), intent(IN) :: PROPFILE
    integer :: NL, I, L

    ! Open the file, get the number of layers and the max number Legendre terms
    open (UNIT=1, FILE=PROPFILE, STATUS='OLD')
    read (1, *) NL
    read (1, *) HEIGHTP(1), TEMPP(1)
    do I = 1, NLAY
        LEGENP(:, I) = 0.0
        read (1, *) HEIGHTP(I + 1), TEMPP(I + 1), TAUP(I), ALBEDOP(I), &
            NLEGP(I), (LEGENP(L, I), L=1, NLEGP(I))
        if (HEIGHTP(I + 1) > HEIGHTP(I)) then
            write (*, *) 'Heights must decrease in property file.'
            stop
        end if
        LEGENP(0, I) = 1.0
    end do
    close (1)
end subroutine READ_PROPERTIES

subroutine READ_CKD_SIZE(CKDFILE, WAVENO, NG, NZCKD)
    ! Find the size of the k-distribution arrays for this wavenumber band.
    integer, intent(OUT) :: NG, NZCKD
    real, intent(IN)  :: WAVENO(2)
    character(LEN=*), intent(IN) :: CKDFILE
    integer :: NB, IB, JB, KB, N
    real    :: WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=CKDFILE, STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, SF, N
        if (abs(WAVENUM1 - WAVENO(1)) < 1.0 .and. &
            abs(WAVENUM2 - WAVENO(2)) < 1.0) then
            backspace (1)
            NG = N
            read (1, *) JB, WAVENUM1, WAVENUM2, SF, N
        end if
    end do
    if (JB == 0) then
        write (*, *) ' Wavenumber range not found in CKD file : ', &
            WAVENO(1), WAVENO(2), CKDFILE
        stop
    end if
    read (1, *) NZCKD
    close (1)
end subroutine READ_CKD_SIZE

subroutine READ_CKD(CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)
    ! Reads the information appropriate for one band from a
    ! correlated k-distribution file.   The number of "g"'s are input in NG,
    ! and the number of levels are input in NZCKD. The wavenumber range (cm^-1)
    ! in the file must match the desired in WAVENO.  The solar flux
    ! is returned in SOLFLUX, the weights or delta g's in DELG, the
    ! levels from the top down in ZCKD, and the absorption coefficients
    ! for each g and level in KABS.
    integer, intent(IN) :: NG, NZCKD
    real, intent(INOUT) :: WAVENO(2)
    real, intent(OUT) :: SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG)
    character(LEN=*), intent(IN) :: CKDFILE
    integer :: NB, IB, JB, KB, IG, I, J, N
    real    :: WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=CKDFILE, STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, SF, N
        if (abs(WAVENUM1 - WAVENO(1)) < 1.0 .and. &
            abs(WAVENUM2 - WAVENO(2)) < 1.0) then
            backspace (1)
            read (1, *) JB, WAVENO(1), WAVENO(2), SOLFLUX, &
                N, (DELG(IG), IG=1, NG)
        end if
    end do
    ! Read the Z level heights
    read (1, *) N
    read (1, *)
    read (1, *)
    do I = 1, NZCKD
        read (1, *) ZCKD(I)
    end do
    ! Skip over the irrelevant bands and read the absorption coefficients
    read (1, *)
    do I = 1, (JB - 1)*NZCKD
        read (1, *)
    end do
    do I = 1, NZCKD
        read (1, *) IB, J, (KABS(I, IG), IG=1, NG)
        do IG = 1, NG
            if (KABS(I, IG) < 0.0) then
                write (*, '(1X,A,A,I3,A,I3,A,I2)') 'Negative absorption in CKD file. ', &
                    'Band: ', IB, '   Level: ', I, '   k: ', IG
            end if
        end do
    end do
    close (1)
end subroutine READ_CKD

subroutine OUTPUT_RESULTS(NLAY, NPTS, NMU, NPHI, NG, &
                          PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
                          SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, &
                          SPLITACC, SOLACC, MAXITER, ITER, &
                          IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
                          NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                          RADOUT, OUTTYPE, OUTFILE)
    !   Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are two types (OUTTYPE) of output: 'R' - for radiance,
    ! 'F' - for hemispheric/actinic flux.
    ! See the overall program documentation for output parameters.
    implicit none
    integer, intent(IN) :: NLAY, NPTS, NMU, NPHI, NG
    integer, intent(IN) :: NZOUT, NMUOUT, NPHIOUT, FLUXPARM
    integer, intent(IN) :: MAXITER, ITER, IXP(2, NLAY)
    real, intent(IN) :: SOLARFLUX, SOLARMU, SFCTEMP, SFCPARMS(*), SKYRAD
    real, intent(IN) :: SOLACC, SPLITACC, WAVELEN, WAVENO(2)
    real, intent(IN) :: HEIGHTG(NPTS)
    real, intent(IN) :: FLUXUP(NPTS), FLUXDN(NPTS), MEANRAD(NPTS)
    real, intent(IN) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
    real, intent(IN) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)
    character(LEN=1), intent(IN) :: SRCTYPE, SFCTYPE, OUTTYPE
    character(LEN=*), intent(IN) :: PROPFILE, CKDFILE, OUTFILE
    integer :: I, J, K, L, N
    real    :: FOURPI
    character(LEN=32) :: SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

    FOURPI = 4*acos(-1.0)
    if (SRCTYPE == 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE == 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE == 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (SFCTYPE == 'L') then
        SFCNAME = 'LAMBERTIAN'
    else if (SFCTYPE == 'R') then
        SFCNAME = 'RPV-BRDF'
    end if

    if (WAVELEN > 0.0) then
        if (OUTTYPE == 'R') then
            UNITSNAME = 'WATTS/(M^2 MICRON STER)'
        else
            UNITSNAME = 'WATTS/(M^2 MICRON)'
        end if
    else
        if (OUTTYPE == 'R') then
            UNITSNAME = 'WATTS/(M^2 STER)'
        else
            UNITSNAME = 'WATTS/(M^2)'
        end if
    end if
    if (OUTTYPE == 'R') OUTNAME = 'RADIANCE'
    if (OUTTYPE == 'F') OUTNAME = 'FLUX'

    open (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
    write (2, '(A)') '! SHDOM Plane-parallel Radiative Transfer Output'
    write (2, '(2(A,I3))') '!  NMU=', NMU, '  NPHI=', NPHI
    write (2, '(2(A,I4))') '!  NLAYER(input)=', NLAY, '    NPTS=', NPTS
    write (2, '(A,A32)') '!  PROPERTY_FILE=', PROPFILE
    write (2, '(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE, &
        '   NUM_G=', NG
    write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    write (2, '(A,A20)') '!  SURFACE_TYPE=', SFCNAME
    if (SFCTYPE == 'L') then
        write (2, '(A,F9.7)') '!  SURFACE_ALBEDO=', SFCPARMS(1)
    else if (SFCTYPE == 'R') then
        write (2, '(3(A,F8.5))') '!  SURFACE_RHO0=', SFCPARMS(1), &
            '   SURFACE_K=', SFCPARMS(2), '   SURFACE_THETA=', SFCPARMS(3)
    end if
    write (2, '(A,F8.3,A,E12.5)') &
        '!  SURFACE_TEMP=', SFCTEMP, '  SKY_RAD=', SKYRAD

    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        write (2, '(A,E13.6,A,F10.7)') &
            '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU
    end if
    if (WAVELEN > 0.0) then
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME, &
            '   WAVELENGTH(MICRON)=', WAVELEN
    else
        write (2, '(A,A18,A,2F10.3)') '!  UNITS=', UNITSNAME, &
            '   WAVENUMBER_RANGE(CM^-1)=', WAVENO(1:2)
    end if
    write (2, '(A,E10.3)') '!  SPLITTING_ACCURACY=', SPLITACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER, &
        '   NUMBER_ITERATIONS=', ITER
    write (2, '(A,A20)') '!  OUTPUT_TYPE=', OUTNAME

    if (OUTTYPE .eq. 'R') then
        ! Radiance output
        write (2, '(3(A,I3))') '!   NZOUT=', NZOUT, &
            '  NMUOUT=', NMUOUT, '  NPHIOUT=', NPHIOUT
        write (2, '(A)') '!    Z       MU     PHI     RADIANCE'
        do I = 1, NZOUT
            do J = 1, NMUOUT
                do K = 1, NPHIOUT
                    write (2, '(1X,F7.3,1X,F8.5,1X,F6.2,1X,E13.6)') &
                        ZOUT(I), MUOUT(J), PHIOUT(K), RADOUT(I, J, K)
                end do
            end do
        end do

    else if (OUTTYPE .eq. 'F') then
        ! Hemispheric flux output
        !   Format 1 is fluxes at every property grid point
        !   Format 2 is fluxes at every grid point
        write (2, '(A)') '!    Z         F_UP        F_DOWN       F_ACTINIC'
        if (FLUXPARM == 1) then
            do L = 1, NLAY
                I = IXP(1, L)
                write (2, '(1X,F7.3,3(2X,E12.5))') &
                    HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
            end do
            I = IXP(2, NLAY)
            write (2, '(1X,F7.3,3(2X,E12.5))') &
                HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
        else
            do L = 1, NLAY
                do I = IXP(1, L), IXP(2, L) - 1
                    write (2, '(1X,F7.3,3(2X,E12.5))') &
                        HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
                end do
            end do
            I = IXP(2, NLAY)
            write (2, '(1X,F7.3,3(2X,E12.5))') &
                HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
        end if

    end if
    close (2)
end subroutine OUTPUT_RESULTS
