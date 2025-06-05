module shd_solv

    use shd_init
    use shd_pha
    use shd_src
    use shd_rad
    use shd_sfc
    use shd_shdo
    use shd_grid
    use shd_intg
    use shd_util
    use shd_mpi

    implicit none

    public :: shd_solv__solveRte

    private :: shd_solv__accelSol
    private :: shd_solv__truncateRad

contains

subroutine shd_solv__solveRte(&
    &NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
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

    ! Performs the SHDOM solution procedure.
    ! Output is returned in SOURCE, RADIANCE, FLUXES, DIRFLUX.
    ! The adaptive grid stucture (NCELLS,GRIDPTR,NEIGHPTR,TREEPTR,
    ! CELLFLAGS,NPTS,GRIDPOS) is input and modified for new grid
    ! cells and points on output.

    integer, intent(in) :: NSTOKES, NX, NY, NZ, NX1, NY1, NXSFC, NYSFC, NSFCPAR
    integer, intent(in) :: ML, MM, NLM, NMU, NPHI, NANG, NLEG, NSTLEG, NUMPHASE
    integer NPHI0(*), MAXITER, ITER, BCFLAG, IPFLAG
    integer MAXIV, MAXIC, MAXIG, MAXIDO
    integer MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR(MAXNBC, 2)
    integer NPTS, NCELLS
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*), IPHASE(*)
    integer WORK1(*)
    logical INRADFLAG, DELTAM, ACCELFLAG, HIGHORDERRAD
    real SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD
    real GNDTEMP, GNDALBEDO
    real DELXSFC, DELYSFC, SFCPARMS(*), SFCGRIDPARMS(*)
    real WAVENO(2), WAVELEN
    real SOLACC, SOLCRIT, SPLITACC, SHACC
    real MU(*), WTDO(*), PHI(*)
    real TEMP(*), PLANCK(*)
    real EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    real XGRID(NX1), YGRID(NY1), ZGRID(NZ)
    real GRIDPOS(3, *)
    real BCRAD(*)
    real FLUXES(2, *), DIRFLUX(*)
    real SOURCE(NSTOKES, *), DELSOURCE(NSTOKES, *)
    real RADIANCE(NSTOKES, *)
    real WORK(*), WORK2(*)
    character(len=1) :: SRCTYPE, UNITS
    character(len=2) :: SFCTYPE

    integer NPHI0MAX, I, ORDINATESET, NBPTS, NBCELLS
    integer MAXNANGBND, NANGBND(4, 2), IERR
    logical FIXSH, SPLITTESTING, DOSPLIT, OUTOFMEM, LAMBERTIAN
    logical UNIFORM_SFC_BRDF
    logical, allocatable :: FFTFLAG(:)
    real ALBMAX
    real STARTADAPTSOL, ENDADAPTSOL, ADAPTRANGE, SPLITCRIT
    real STARTSPLITACC, CURSPLITACC, AVGSOLCRIT, BETA
    real DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR
    real, allocatable :: CMU1(:, :, :), CMU2(:, :, :)
    real, allocatable :: CPHI1(:, :, :), CPHI2(:, :, :)
    real, allocatable :: YLMSUN(:, :), WTMU(:), WPHISAVE(:)
    real, allocatable :: SFC_BRDF_DO(:, :, :, :, :, :)

    integer :: count_start, count_rate, count_max
    integer :: count_end
    real :: elapsed_time

    ! Allocate local arrays
    ! NSTLEG=1 for NSTOKES=1, otherwise NSTLEG=6
    allocate (CMU1(NSTLEG, NLM, NMU), CMU2(NSTLEG, NLM, NMU))
    allocate (CPHI1(-16:16, 32, NMU), CPHI2(32, -16:16, NMU))
    allocate (YLMSUN(NSTLEG, NLM))
    allocate (FFTFLAG(NMU), WTMU(NMU), WPHISAVE(NMU*(3*NPHI + 15)))

    ! Set up some things before solution loop

    ! Transfer the medium properties to the internal grid and add gas abs
    call shd_grid__interpGrid(NPTS, NSTLEG, NLEG, GRIDPOS,&
        &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)

    ! If Delta-M then scale the extinction, albedo, and Legendre terms
    ! Put the Planck blackbody source in PLANCK.
    call shd_init__makePrp(ML, MM, NSTLEG, NLEG, NPTS, DELTAM, NUMPHASE,&
        &SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX,&
        &EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)

    ! Get the maximum single scattering albedo over all processors
    call shd_mpi__findAlbMax(ALBMAX)

    ! Compute the solar transmission in DIRFLUX.
    ! Precompute Ylm's for solar direction
    if (SRCTYPE .ne. 'T') then
        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            call shd_mpi__makeDirect(1, NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                &NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                &GRIDPOS, NX, XGRID, NY, YGRID, DIRFLUX)
        else
            call shd_rad__makeDirect(NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                &NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                &GRIDPOS, DIRFLUX)
        end if
        call shd_pha__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG,&
            &YLMSUN)
    end if

    ! Make the discrete ordinates (angles)
    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
    ORDINATESET = 2
    call shd_init__makeAngle(NMU, NPHI, ORDINATESET, NPHI0MAX,&
        &NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the Ylm transform coefficients
    call shd_shdo__makeSHDO(NSTLEG, ML, MM, NLM, NMU, NPHI0,&
        &NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,&
        &CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

    FIXSH = .false.
    if (.not. INRADFLAG) then
        ! Initialize the radiance on the base grid using Eddington
        ! two-stream plane-parallel
        call shd_init__makeRad(NSTOKES, NX1*NY1, NZ, NSTLEG, NLEG,&
            &RSHPTR, ZGRID,&
            &EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,&
            &SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP, SKYRAD,&
            &UNITS, WAVENO, WAVELEN, RADIANCE)
        ! Interpolate the radiance on the non-base grid points from the
        !   base grid points
        NBPTS = NX1*NY1*NZ
        NBCELLS = (NZ - 1)*(NX + ibits(BCFLAG, 0, 1) - ibits(BCFLAG, 2, 1))&
            &*(NY + ibits(BCFLAG, 1, 1) - ibits(BCFLAG, 3, 1))
        call shd_rad__interpRad(NSTOKES, NBPTS, NPTS, RSHPTR, RADIANCE,&
            &NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        ! Initialize the source function from the radiance field
        call shd_src__computeSrc(NSTOKES, ML, MM, NLM,&
            &NSTLEG, NLEG, NUMPHASE, NPTS,&
            &FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            &IPHASE, PLANCK, DIRFLUX, SHACC,&
            &MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            &.true., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)
    end if

    if (ACCELFLAG) then
        OSHPTR(1:NPTS + 1) = SHPTR(1:NPTS + 1)
        DELSOURCE(:, 1:OSHPTR(NPTS + 1)) = 0.0
    end if
    SOLCRIT = 1.0

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'
    ! Make the pointers to the grid points on the top and bottom boundaries
    call shd_grid__findBoundaryPnt(NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
        &ZGRID(1), ZGRID(NZ), GRIDPOS,&
        &NTOPPTS, NBOTPTS, BCPTR)
    if (SFCTYPE(1:1) .eq. 'V') then
        ! If this is a variable surface then bilinearly interpolate
        !   the surface temperature and reflection parameters to the
        !   bottom grid points.
        call shd_sfc__interpSfc(NBOTPTS, BCPTR(1, 2), GRIDPOS,&
            &SRCTYPE, WAVENO, WAVELEN, UNITS,&
            &NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
            &SFCGRIDPARMS)
        call shd_init__checkSfcU(NXSFC, NYSFC, NSFCPAR, SFCPARMS,&
            &UNIFORM_SFC_BRDF)
        if (UNIFORM_SFC_BRDF) then
            if (NSTOKES*(NMU/2)*NPHI0MAX .gt.&
                &sqrt(1.0*huge(NMU)/kind(NMU))) then
                IERR = 1
            else
                allocate (SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX,&
                    &NSTOKES, NMU/2, NPHI0MAX), STAT=IERR)
            end if
            if (IERR .eq. 0) then
                call shd_shdo__makeSfcDO(SFCTYPE, WAVELEN,&
                    &NXSFC, NYSFC, NSFCPAR, SFCPARMS, NSTOKES,&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                    &SFC_BRDF_DO)
            else
                print '(A,A)', 'Not enough memory to allocate precomputed',&
                    &' surface BRDF discrete ordinate matrix: '
                print *, '  Proceeding more slowly.'
                UNIFORM_SFC_BRDF = .false.
            end if
        end if
    end if

    ITER = 0

    if (MAXITER .le. ITER) return
    ! Starting values for the adaptive cell splitting controlling method.
    ! The splitting accuracy used at a particular iteration is CURSPLITACC,
    ! and it is started out high an lowered to the desired SPLITACC.
    ! The cell splitting starts when the solution criterion (SOLCRIT)
    ! goes below STARTADAPTSOL and ends when the solution criterion
    ! goes below ENDADAPTSOL, though the current splitting accuracy
    ! is designed to get to the final SPLITACC over a smaller range
    ! (ADAPTRANGE).  The currect splitting accuracy is computed from
    ! a geometric average of the solution criterion (AVGSOLCRIT).
    ! The shd_grid__splitGrid routine is called and the splitting criterion returned
    ! as long as SPLITTESTING is true.
    SPLITTESTING = .true.
    OUTOFMEM = .false.
    ENDADAPTSOL = 0.001
    STARTADAPTSOL = 0.1
    ADAPTRANGE = STARTADAPTSOL/(3.0*ENDADAPTSOL)
    CURSPLITACC = SPLITACC*ADAPTRANGE
    STARTSPLITACC = CURSPLITACC
    AVGSOLCRIT = SOLCRIT
    SPLITCRIT = 0.0

    write (6, '(A)') 'Message [_solveRte]:'
    write (6, '(A)') '+ Iter  Log(Sol)  SplitCrit  Npoints   Ncells  Nsh(avg)'

    ! Main solution loop
    !   Iterate until the iterations exceed the limit, or the desired
    !   solution criterion is reached and the current splitting accuracy
    !   is down to the desired level.
    do while (ITER .lt. MAXITER .and. (SOLCRIT .gt. SOLACC&
        &.or. (SPLITCRIT .gt. SPLITACC .and.&
        &CURSPLITACC .gt. SPLITACC .and. .not. OUTOFMEM)))

        call system_clock(count_start, count_rate, count_max)

        ITER = ITER + 1

        ! Don't do the adaptive grid cell splitting stuff unless wanted
        if (SPLITACC .gt. 0.0) then
            ! If in the right solution range then try to split some cells, and
            !   interpolate the medium and source functions to the new points.
            ! The current splitting accuracy is started high and then lowered.
            AVGSOLCRIT = sqrt(AVGSOLCRIT*SOLCRIT)
            DOSPLIT = SOLCRIT .le. STARTADAPTSOL .and.&
                &(SOLCRIT .gt. ENDADAPTSOL .or. CURSPLITACC .gt. SPLITACC)&
                &.and. .not. OUTOFMEM
            ! Make sure all processors are going to split cells if any want to
            call shd_mpi__unifySplit(DOSPLIT, STARTSPLITACC)
            BETA = log(STARTSPLITACC/SPLITACC)/log(ADAPTRANGE)
            CURSPLITACC = min(CURSPLITACC, max(SPLITACC,&
                &SPLITACC*(AVGSOLCRIT/(3.0*ENDADAPTSOL))**BETA))
            if (SOLCRIT .le. ENDADAPTSOL) CURSPLITACC = SPLITACC
            if (SPLITTESTING) then
                call shd_grid__splitGrid(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
                    &DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
                    &NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &GRIDPOS, NX, XGRID, NY, YGRID,&
                    &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
                    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
                    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
                    &TEMP, PLANCK, DIRFLUX, SRCTYPE,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                    &UNITS, WAVENO, WAVELEN, WORK, WORK2)
                if (SOLCRIT .gt. STARTADAPTSOL) STARTSPLITACC = SPLITCRIT
            end if
            ! Find the maximum splitting criterion over all processors
            call shd_mpi__findTotSplitCritMax(SPLITCRIT)
            if (SOLCRIT .le. ENDADAPTSOL) SPLITTESTING = .false.
        else
            DOSPLIT = .false.
        end if

        ! Find the SH truncation to use for the radiance (RSHPTR)
        call shd_solv__truncateRad(HIGHORDERRAD,&
            &NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE,&
            &NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
            &MAXIV + MAXIG, FIXSH, SHACC, RSHPTR)

        ! Integrate the source function along discrete ordinates
        !   to compute the radiance field.  The input source function
        !   and output radiance are in spherical harmonic expansions,
        !   and the SH/DO transforms are done during the loop over
        !   discrete ordinates.
        call shd_intg__integratePath(NX, NY, NZ, NSTOKES, NSTLEG, ML, MM, NLM,&
            &NMU, NPHI0MAX, NPHI0, NANG, BCFLAG, IPFLAG,&
            &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
            &SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN,&
            &UNITS, SFCTYPE, GNDTEMP, GNDALBEDO,&
            &NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
            &SFCGRIDPARMS, UNIFORM_SFC_BRDF, SFC_BRDF_DO,&
            &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
            &MU, PHI, WTDO,&
            &CMU1, CPHI1, CMU2, CPHI2, WPHISAVE, FFTFLAG,&
            &DIRFLUX, FLUXES, EXTINCT,&
            &SHPTR, SOURCE, RSHPTR, RADIANCE,&
            &WORK, WORK1, WORK2)

        ! Compute the source function from the radiance field,
        !   do the adaptive spherical harmonics truncation, compute
        !   the solution criterion, and dot products for acceleration
        if (((SOLCRIT .lt. ENDADAPTSOL) .or. (ITER .gt. 30)) .and.&
            &((NX .gt. 10) .or. (NY .gt. 10))) FIXSH = .true.
        call shd_src__computeSrc(NSTOKES, ML, MM, NLM,&
            &NSTLEG, NLEG, NUMPHASE, NPTS,&
            &FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            &IPHASE, PLANCK, DIRFLUX, SHACC,&
            &MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            &.false., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

        ! Calculate the acceleration parameter and solution criterion
        !   from all the processors
        call shd_mpi__calAccelSolCrit(ACCELFLAG, DELJDOT, DELJOLD, DELJNEW,&
            &JNORM, ACCELPAR, SOLCRIT)

        ! Accelerate the convergence of the source function vector
        call shd_solv__accelSol(ACCELPAR, NPTS, SHPTR, OSHPTR,&
            &NSTOKES, SOURCE, DELSOURCE)

        ! If it is a not scattering medium then do only one iteration
        if (ALBMAX .lt. SOLACC) SOLCRIT = SOLACC

        call system_clock(count_end, count_rate, count_max)
        elapsed_time = (real(count_end) - real(count_start)) / real(count_rate)

        ! Print out the log solution criterion, number of points,
        !   and average SH truncation.
        if (DOSPLIT) then
            write (6, '("S",1X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                &SPLITCRIT,&
                &NPTS,&
                &NCELLS,&
                &float(SHPTR(NPTS + 1))/NPTS,&
                &elapsed_time
                !&float(SHPTR(NPTS + 1))/(NPTS*NLM)
        else
            if (OUTOFMEM) then
                write (6, '("M",1X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                    &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                    &SPLITCRIT,&
                    &NPTS,&
                    &NCELLS,&
                    &float(SHPTR(NPTS + 1))/NPTS,&
                    &elapsed_time
                    !&float(SHPTR(NPTS + 1))/(NPTS*NLM)
            else
                write (6, '(2X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                    &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                    &SPLITCRIT,&
                    &NPTS,&
                    &NCELLS,&
                    &float(SHPTR(NPTS + 1))/NPTS,&
                    &elapsed_time
                    !&float(SHPTR(NPTS + 1))/(NPTS*NLM)
            end if
        end if
    end do

    write (6, '(A,I4,A,F9.6)') '- Iterations: ', ITER,&
        &', Final Criterion: ', SOLCRIT

    if (SFCTYPE(1:1) .eq. 'V' .and. UNIFORM_SFC_BRDF) then
        deallocate (SFC_BRDF_DO)
    end if
    deallocate (YLMSUN, FFTFLAG, WTMU, WPHISAVE)
    deallocate (CMU1, CMU2, CPHI1, CPHI2)

    ! Comment in for GLE graphics output of cell structure (see routine)
    ! call shdUtil_io__plotCell('cellfinal.gle', YGRID(1), 1, IPFLAG,&
    !     &NX, XGRID(NX+1)-XGRID(1), ZGRID(NZ)-ZGRID(1),&
    !     &NCELLS, GRIDPTR, NEIGHPTR, TREEPTR,&
    !     &CELLFLAGS, GRIDPOS)

    ! Comment in for output of cell splitting criterion for every cell
    ! call shdUtil_io_outputCellSplit ('cellsplit.dat',GRIDPTR,&
    !     &GRIDPOS, EXTINCT, SHPTR, SOURCE, NCELLS)

    return

end subroutine shd_solv__solveRte

subroutine shd_solv__accelSol(&
    &ACCELPAR, NPTS, SHPTR, OSHPTR,&
    &NSTOKES, SOURCE, DELSOURCE)

    !   Accelerates the successive order of scattering series using the
    ! input acceleration parameter (ACCELPAR).

    integer, intent(in) :: NPTS, SHPTR(*), OSHPTR(*), NSTOKES
    real, intent(in) :: ACCELPAR, DELSOURCE(NSTOKES, *)
    real, intent(inout) :: SOURCE(NSTOKES, *)
    integer I, J, IS, ISD, K, NS, NSD, NSC

    if (ACCELPAR .gt. 0.0) then
        do I = 1, NPTS
            IS = SHPTR(I)
            NS = SHPTR(I + 1) - IS
            ISD = OSHPTR(I)
            NSD = OSHPTR(I + 1) - ISD
            NSC = min(NS, NSD)
            do J = 1, NSC
                SOURCE(:, IS + J) = SOURCE(:, IS + J)&
                    &+ ACCELPAR*DELSOURCE(:, ISD + J)
            end do
        end do
    end if

    return
end subroutine shd_solv__accelSol

subroutine shd_solv__truncateRad(HIGHORDERRAD,&
    &NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE,&
    &NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
    &MAXIR, FIXSH, SHACC, RSHPTR)

    ! Computes the next radiance spherical harmonic series truncation
    ! to use based on the current source function truncation and the
    ! scattering properties.  The radiance truncation is the smaller
    ! of that based on a minimum SH term radiance and that needed for
    ! the computation of the source function.  The truncation is always
    ! made for a complete set in the L index.  The radiance SH truncation
    ! is set to a few higher L index than the source function, so that
    ! the source function truncation may grow (since initialized with
    ! a low order truncation).  The minimum radiance truncation is
    ! L=1 so that the mean radiance and net flux terms may be computed,
    ! unless HIGHORDERRAD is on, in which case all terms are kept.

    integer NPTS, NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE, MAXIR
    integer SHPTR(NPTS + 1), RSHPTR(NPTS + 2)
    integer*2 IPHASE(NPTS)
    logical HIGHORDERRAD, FIXSH
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *), RADIANCE(NSTOKES, *)
    real SHACC
    integer IRO, IR, NRO, NR, I, J, K, L, LR, LS, M, MS, ME, NS, NLM
    integer, allocatable :: LOFJ(:)
    logical NOTEND
    real RADMIN

    if (FIXSH) goto 190

    RADMIN = SHACC
    NOTEND = .true.
    ! Make the l index as a function of SH term (J)
    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    allocate (LOFJ(NLM))
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    IRO = 0
    IR = 0
    RSHPTR(1) = 0
    do I = 1, NPTS
        ! Find how many terms wanted for the radiance based on the
        !   radiance threshold (RADMIN) (use the mean radiance as an
        !   upper limit to the SH term radiance).
        !   We don't have radiances for any new points, so don't use them
        NRO = RSHPTR(I + 1) - IRO
        if (NRO .eq. 0) NOTEND = .false.
        if (NOTEND) then
            if (NUMPHASE .gt. 0) then
                K = IPHASE(I)
            else
                K = I
            end if
            LR = 1
            do L = 1, ML
                if (ALBEDO(I)*LEGEN(1, L, K)*RADIANCE(1, IRO + 1) .gt. RADMIN) LR = L
            end do
        else
            LR = ML
        end if
        IRO = IRO + NRO

        ! The radiance L truncation is the smaller of the radiance cutoff
        !   from above and a few more than the source truncation, but
        !   it is at least 1 (so low order terms are always available).
        !   For higher order radiance option, use all terms.
        !   But if out of memory, use the source function truncation.
        NS = max(1, SHPTR(I + 1) - SHPTR(I))
        LS = LOFJ(NS)
        LR = min(LR, LS + ML/8 + 2)
        if (HIGHORDERRAD) LR = ML
        if (LR .le. MM) then
            NR = (LR*(LR + 1)) + LR + 1
        else
            NR = (2*MM + 1)*LR - (MM*(1 + (MM - 1))) + MM + 1
        end if
        IR = IR + NR
        RSHPTR(I + 1) = IR
        if (IR .gt. MAXIR) then
            write (6, *)&
                &'Message [_truncateRad]: Out of memory for more radiance terms.'
            goto 190
        end if
    end do
    RSHPTR(NPTS + 2) = IR
    deallocate (LOFJ)
    return

    ! If out of memory or the SH truncation is fixed then
    ! just set the radiance truncation to that of the source.
190 continue

    RSHPTR(1) = 0
    IR = 0
    do I = 1, NPTS
        NR = max(4, SHPTR(I + 1) - SHPTR(I))
        if (HIGHORDERRAD) then
            if (ML .le. MM) then
                NR = (ML*(ML + 1)) + ML + 1
            else
                NR = (2*MM + 1)*ML - (MM*(1 + (MM - 1))) + MM + 1
            end if
        end if
        IR = IR + NR
        if (IR .gt. MAXIR) then
            write (6, *) 'Error [_truncateRad]: Really out of memory ',&
                &'for more radiance terms.'
            write (6, *) '  Increase MAXIV.'
            stop
        end if
        RSHPTR(I + 1) = IR
    end do
    RSHPTR(NPTS + 2) = IR
    return
end subroutine shd_solv__truncateRad

end module shd_solv
