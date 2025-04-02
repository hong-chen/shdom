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

    use shd_init
    use shd_ext_mpi

    implicit none
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
    call shd_comp__interpGrid(NPTS, NSTLEG, NLEG, GRIDPOS,&
        &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)

    ! If Delta-M then scale the extinction, albedo, and Legendre terms
    ! Put the Planck blackbody source in PLANCK.
    call shd_init__makeProp(ML, MM, NSTLEG, NLEG, NPTS, DELTAM, NUMPHASE,&
        &SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX,&
        &EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)
    ! Get the maximum single scattering albedo over all processors
    call shd_ext__mpiFindAlbMax(ALBMAX)

    ! Compute the solar transmission in DIRFLUX.
    ! Precompute Ylm's for solar direction
    if (SRCTYPE .ne. 'T') then
        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            call shd_ext__mpiMakeDirect(1, NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                &NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                &GRIDPOS, NX, XGRID, NY, YGRID, DIRFLUX)
        else
            call shd_comp__makeDirect(NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                &NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                &GRIDPOS, DIRFLUX)
        end if
        call shd_comp__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG,&
            &YLMSUN)
    end if

    ! Make the discrete ordinates (angles)
    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
    ORDINATESET = 2
    call shd_init__makeAngle(NMU, NPHI, ORDINATESET, NPHI0MAX,&
        &NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the Ylm transform coefficients
    call shd_solv__makeSHDO(NSTLEG, ML, MM, NLM, NMU, NPHI0,&
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
        call shd_comp__interpRad(NSTOKES, NBPTS, NPTS, RSHPTR, RADIANCE,&
            &NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        ! Initialize the source function from the radiance field
        call shd_comp__computeSrc(NSTOKES, ML, MM, NLM,&
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
    call shd_solv__findBoundaryPnt(NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
        &ZGRID(1), ZGRID(NZ), GRIDPOS,&
        &NTOPPTS, NBOTPTS, BCPTR)
    if (SFCTYPE(1:1) .eq. 'V') then
        ! If this is a variable surface then bilinearly interpolate
        !   the surface temperature and reflection parameters to the
        !   bottom grid points.
        call shd_comp__interpSurfParams(NBOTPTS, BCPTR(1, 2), GRIDPOS,&
            &SRCTYPE, WAVENO, WAVELEN, UNITS,&
            &NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
            &SFCGRIDPARMS)
        call shd_init__checkSurfUnif(NXSFC, NYSFC, NSFCPAR, SFCPARMS,&
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
                call shd_solv__makeSurfDO(SFCTYPE, WAVELEN,&
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
    ! The shd_solv__splitGrid routine is called and the splitting criterion returned
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
    write (6, '(A)') '+ Iter  Log(Sol)  SplitCrit  Npoints  Nsh(avg)'

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
            call shd_ext__mpiUnifySplit(DOSPLIT, STARTSPLITACC)
            BETA = log(STARTSPLITACC/SPLITACC)/log(ADAPTRANGE)
            CURSPLITACC = min(CURSPLITACC, max(SPLITACC,&
                &SPLITACC*(AVGSOLCRIT/(3.0*ENDADAPTSOL))**BETA))
            if (SOLCRIT .le. ENDADAPTSOL) CURSPLITACC = SPLITACC
            if (SPLITTESTING) then
                call shd_solv__splitGrid(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
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
            call shd_ext__mpiFindTotSplitCritMax(SPLITCRIT)
            if (SOLCRIT .le. ENDADAPTSOL) SPLITTESTING = .false.
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
        call shd_solv__integratePath(NX, NY, NZ, NSTOKES, NSTLEG, ML, MM, NLM,&
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
        if (SOLCRIT .lt. ENDADAPTSOL .or. ITER .gt. 30) FIXSH = .true.
        call shd_comp__computeSrc(NSTOKES, ML, MM, NLM,&
            &NSTLEG, NLEG, NUMPHASE, NPTS,&
            &FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            &IPHASE, PLANCK, DIRFLUX, SHACC,&
            &MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            &.false., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

        ! Calculate the acceleration parameter and solution criterion
        !   from all the processors
        call shd_ext__mpiCalAccelSolCrit(ACCELFLAG, DELJDOT, DELJOLD, DELJNEW,&
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
        write (6, '(2X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,F8.2,2X,F6.2," sec")')&
            &ITER, log10(max(SOLCRIT, 1.0e-20)), SPLITCRIT,&
            &NPTS, FLOAT(SHPTR(NPTS + 1))/NPTS,&
            &elapsed_time
            !&FLOAT(SHPTR(NPTS + 1))/(NPTS*NLM)
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

subroutine shd_solv__calSrcPnt(&
    &NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    &SOURCET)

    ! Calculates the polarized source function in real generalized
    ! spherical harmonics at one grid point from the radiance in
    ! real generalized spherical harmonics and thermal and/or solar
    ! source information.

    implicit none

    integer, intent(in) :: NSTOKES, NSTLEG, NLEG, NLM, LOFJ(NLM), IPH, NR
    real, intent(in) :: FLUX0, YLMSUN(NSTLEG, NLM), PLANCK, ALBEDO
    real, intent(in) :: LEGEN(NSTLEG, 0:NLEG, *)
    real, intent(in) :: RADIANCE(NSTOKES, NR)
    character(len=1), intent(in) :: SRCTYPE

    real, intent(out) :: SOURCET(NSTOKES, NLM)

    real, parameter :: C = 3.544907703
    integer J

    SOURCET(:, :) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.  Do different parts of the
    ! matrix multiply depending on NSTOKES.
    do J = 1, NR
        SOURCET(1, J) = SOURCET(1, J) + LEGEN(1, LOFJ(J), IPH)*RADIANCE(1, J)
    end do

    if (NSTOKES .gt. 1) then

        do J = 1, NR
            SOURCET(1, J) = SOURCET(1, J) + LEGEN(5, LOFJ(J), IPH)*RADIANCE(2, J)
        end do

        do J = 5, NR
            SOURCET(2, J) = SOURCET(2, J) + LEGEN(5, LOFJ(J), IPH)*RADIANCE(1, J)&
                &+ LEGEN(2, LOFJ(J), IPH)*RADIANCE(2, J)
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(3, LOFJ(J), IPH)*RADIANCE(3, J)
        end do

    end if

    if (NSTOKES .eq. 4) then

        do J = 1, NR
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(6, LOFJ(J), IPH)*RADIANCE(4, J)
            SOURCET(4, J) = SOURCET(4, J) - LEGEN(6, LOFJ(J), IPH)*RADIANCE(3, J)&
                &+ LEGEN(4, LOFJ(J), IPH)*RADIANCE(4, J)
        end do

    end if

    ! Do the final multiplication by the single scattering albedo
    SOURCET(:, :) = ALBEDO*SOURCET(:, :)

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        do J = 1, NLM
            SOURCET(1, J) = SOURCET(1, J)&
                &+ FLUX0*ALBEDO*LEGEN(1, LOFJ(J), IPH)*YLMSUN(1, J)
        end do

        if (NSTOKES .gt. 1) then
            do J = 5, NLM
                SOURCET(2, J) = SOURCET(2, J)&
                    &+ FLUX0*ALBEDO*LEGEN(5, LOFJ(J), IPH)*YLMSUN(1, J)
            end do
        end if

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1, 1) = SOURCET(1, 1) + C*PLANCK
    end if

    return

end subroutine shd_solv__calSrcPnt

subroutine shd_solv__calSrcPnt_(&
    &NLM, NLEG, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    &SOURCET)

    ! Calculates the scalar source function in spherical harmonics at
    ! one grid point from the radiance in spherical harmonics and
    ! thermal and/or solar source information.

    implicit none

    integer, intent(in) :: NLEG, NLM, LOFJ(NLM), IPH, NR
    character(len=1), intent(in) :: SRCTYPE
    real, intent(in) :: FLUX0, YLMSUN(NLM), PLANCK, ALBEDO
    real, intent(in) :: LEGEN(0:NLEG, *)
    real, intent(in) :: RADIANCE(NR)

    real, intent(out) :: SOURCET(NLM)

    real, parameter :: C=3.544907703
    integer J


    SOURCET(:) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        do J = 1, NLM
            SOURCET(J) = FLUX0*ALBEDO*LEGEN(LOFJ(J), IPH)*YLMSUN(J)
        end do

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1) = SOURCET(1) + C*PLANCK
    end if

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.
    do J = 1, NR
        SOURCET(J) = SOURCET(J) + ALBEDO*LEGEN(LOFJ(J), IPH)*RADIANCE(J)
    end do

    return

end subroutine shd_solv__calSrcPnt_

subroutine shd_comp__computeSrc(&
    &NSTOKES, ML, MM, NLM,&
    &NSTLEG, NLEG, NUMPHASE, NPTS, FIXSH,&
    &SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
    &IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
    &FIRST, ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

    ! Computes the source function (SOURCE) in spherical harmonic space
    !   for all the grid points.  The thermal source and/or solar
    !   pseudo-source (in PLANCK or DIRFLUX) is added to the scattering source
    !   (computed from LEGEN and the spherical harmonic expansion in RADIAIANCE).
    !   Also computes the adaptive spherical harmonic truncation for
    !   the source function, unless FIXSH is true.  The SH truncation is set
    !   so that the higher terms of the source function are smaller in
    !   absolute value than SHACC. The truncation is to the nearest l,
    !   so that a triangular (for M=L) truncation is preserved.  If there
    !   is little or no scattering then the SH truncation can be set to
    !   zero terms.  The source function SH truncation is output in terms
    !   of the pointer SHPTR (note that RSHPTR refers to the RADIANCE array
    !   and OSHPTR to the DELSOURCE array).
    ! To save memory this routine also computes some things for the
    !   series acceleration: the dot products of the difference in source
    !   function between successive iterations and the new difference source
    !   function vector (DELSOURCE).  To avoid using another array and
    !   because of the adaptive truncation pointer stuff, this requires
    !   two extra calculations of the source function.
    ! Computes the sum of squares of the new source function for calculation
    !   of the solution criterion (SOLCRIT, the RMS difference in succesive
    !   source function fields normalized by the RMS of the field).

    implicit none
    integer NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG
    integer NUMPHASE, MAXIV
    integer RSHPTR(NPTS + 1), SHPTR(NPTS + 1), OSHPTR(NPTS + 1)
    integer*2 IPHASE(NPTS)
    logical FIXSH, FIRST, ACCELFLAG
    real SHACC, YLMSUN(NSTLEG, NLM), SOLARMU
    real DELJDOT, DELJOLD, DELJNEW, JNORM
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real PLANCK(NPTS), DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *)
    real SOURCE(NSTOKES, *), DELSOURCE(NSTOKES, *)
    character(len=1) :: SRCTYPE
    integer IS, ISO, IR, I, IPH, J, JS, K, L, LS, M, MS, ME, NS, NR
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)
    real SRCMIN, C, SECMU0, D

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))

    SRCMIN = SHACC
    ! Make source function from scattering integral and real sources:
    !   Solar source: DIRFLUX has F0*exp(-tau)
    !   Thermal source: PLANCK has (1-omega)*B(T)
    SECMU0 = 1.0/abs(SOLARMU)
    C = sqrt(4.0*acos(-1.0))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    ! Compute the dot products of the successive source differences
    if (.not. FIRST) then
        DELJDOT = 0.0
        DELJOLD = 0.0
        DELJNEW = 0.0
        JNORM = 0.0
        do I = 1, NPTS
            ! Compute the temporary source function for this point
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            IR = RSHPTR(I)
            NR = RSHPTR(I + 1) - IR
            if (NR .gt. NLM) then
                write (6, *) 'Error [_computeSrc]: NR>NLM 1', I, IR, RSHPTR(I + 1)
                stop
            end if
            if (NSTOKES .eq. 1) then
                call shd_solv__calSrcPnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_solv__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Compute the dot product of the new and old source difference,
            !   and the solution criterion (don't use new points for Jnorm).
            if (ACCELFLAG) then
                IS = SHPTR(I)
                ISO = OSHPTR(I)
                NS = min(SHPTR(I + 1) - IS, OSHPTR(I + 1) - ISO)
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJDOT = DELJDOT +&
                            &(SOURCET(K, J) - SOURCE(K, IS + J))*DELSOURCE(K, ISO + J)
                        DELJOLD = DELJOLD + DELSOURCE(K, ISO + J)**2
                        DELJNEW = DELJNEW + (SOURCET(K, J) - SOURCE(K, IS + J))**2
                        JNORM = JNORM + SOURCE(K, IS + J)**2
                    end do
                end do
            else
                IS = SHPTR(I)
                NS = SHPTR(I + 1) - IS
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJNEW = DELJNEW + (SOURCET(K, J) - SOURCE(K, IS + J))**2
                        JNORM = JNORM + SOURCE(K, IS + J)**2
                    end do
                end do
            end if
        end do
    end if

    ! Compute and store the new source function difference vector
    if (.not. FIRST .and. ACCELFLAG) then
        do I = 1, NPTS
            ! Compute the temporary source function for this point
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            IR = RSHPTR(I)
            NR = RSHPTR(I + 1) - IR
            if (NSTOKES .eq. 1) then
                call shd_solv__calSrcPnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_solv__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Make the source difference vector
            IS = SHPTR(I)
            NS = SHPTR(I + 1) - IS
            OSHPTR(I) = IS
            do J = 1, NS
                DELSOURCE(:, IS + J) = SOURCET(:, J) - SOURCE(:, IS + J)
            end do
        end do
        OSHPTR(NPTS + 1) = SHPTR(NPTS + 1)
    end if

    ! Compute the source function in the temporary array.
    IS = 0
    do I = 1, NPTS
        if (NUMPHASE .gt. 0) then
            IPH = IPHASE(I)
        else
            IPH = I
        end if
        IR = RSHPTR(I)
        NR = RSHPTR(I + 1) - IR
        if (NR .gt. NLM) then
            write (6, *) 'Error [_computeSrc]: NR>NLM 3', I, IR, RSHPTR(I + 1)
            stop
        end if
        if (NSTOKES .eq. 1) then
            call shd_solv__calSrcPnt_(NLM, NLEG, LOFJ,&
                &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                &NR, RADIANCE(1, IR + 1), SOURCET)
        else
            call shd_solv__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                &NR, RADIANCE(1, IR + 1), SOURCET)
        end if
        if (FIXSH) then
            NS = SHPTR(I + 1) - IS
        else
            ! Then do the adaptive spherical harmonics:
            !   Find how many source function terms to use, and how many
            !   terms wanted for the radiance.
            if (SRCTYPE .eq. 'S') then
                JS = 0
            else
                JS = 1
            end if
            do J = 1, NLM
                do K = 1, NSTOKES
                    if (abs(SOURCET(K, J)) .gt. SRCMIN) JS = J
                end do
            end do
            ! Use all the m's for the last l for the source function.
            if (JS .eq. 0) then
                NS = 0
            else
                LS = LOFJ(JS)
                if (LS .le. MM) then
                    NS = (LS*(LS + 1)) + LS + 1
                else
                    NS = (2*MM + 1)*LS - (MM*(1 + (MM - 1))) + MM + 1
                end if
            end if
            SHPTR(I) = IS
        end if
        if (IS + NS .gt. MAXIV) then
            write (6, *) 'Error [_computeSrc]: MAXIV exceeded.', MAXIV
            write (6, *) '  Out of memory for more spherical harmonic terms.'
            stop
        end if
        ! Transfer the new source function
        do J = 1, NS
            SOURCE(:, IS + J) = SOURCET(:, J)
        end do
        IS = IS + NS
    end do
    SHPTR(NPTS + 1) = IS
    deallocate (LOFJ, SOURCET)
    return
end subroutine shd_comp__computeSrc

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

    implicit none
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

subroutine shd_solv__accelSol(&
    &ACCELPAR, NPTS, SHPTR, OSHPTR,&
    &NSTOKES, SOURCE, DELSOURCE)

    !   Accelerates the successive order of scattering series using the
    ! input acceleration parameter (ACCELPAR).

    implicit none
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

subroutine shd_solv__findBoundaryPnt(&
    &NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
    &ZBOT, ZTOP, GRIDPOS,&
    &NTOPPTS, NBOTPTS, BCPTR)

    !   Returns the BCPTR pointer array for the NTOPPTS top boundary
    ! points (BCPTR(*,1)) and the NBOTPTS bottom boundary points (BCPTR(*,2)).
    ! All the grid points are looked at to find the boundary ones.

    implicit none
    integer, intent(in) :: NPTS, NANG, MAXNBC, MAXBCRAD
    integer, intent(inout) :: NTOPPTS, NBOTPTS
    logical, intent(in) :: LAMBERTIAN
    real, intent(in) :: ZBOT, ZTOP, GRIDPOS(3, NPTS)
    integer, intent(out) :: BCPTR(MAXNBC, 2)
    integer :: I, IT, IB, NA

    ! Loop over all points finding the top and bottom ones
    if (LAMBERTIAN) then
        NA = 1
    else
        NA = NANG/2 + 1
    end if

    IT = 0
    IB = 0
    do I = 1, NPTS

        if (GRIDPOS(3, I) .ge. ZTOP) then
            IT = IT + 1
            if (IT .gt. MAXNBC .or. IT .gt. MAXBCRAD) then
                write (6, *) 'Error [_findBoundaryPnt]: MAXNBC exceeded - top.'
                stop
            end if
            BCPTR(IT, 1) = I
        end if

        if (GRIDPOS(3, I) .le. ZBOT) then
            IB = IB + 1
            if (IB .gt. MAXNBC .or. IT + IB*NA .gt. MAXBCRAD) then
                write (6, *) 'Error [_findBoundaryPnt]: MAXNBC exceeded - bottom.'
                stop
            end if
            BCPTR(IB, 2) = I
        end if

    end do

    NTOPPTS = IT
    NBOTPTS = IB

    return

end subroutine shd_solv__findBoundaryPnt

subroutine shd_solv__makeSurfDO(SFCTYPE, WAVELEN,&
    &NXSFC, NYSFC, NSFCPAR, SFCPARMS, NSTOKES,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &SFC_BRDF_DO)

    !   Makes the precomputed BRDF surface reflection matrix for
    ! discrete ordinate incident and outgoing directions.  For each
    ! pair of discrete ordinates shd_surf__brdfSurface is called and the
    ! integration weights are included.

    use shd_surf

    implicit none

    integer NXSFC, NYSFC, NSFCPAR
    integer NSTOKES, NMU, NPHI0MAX, NPHI0(NMU)
    real SFCPARMS(NSFCPAR, NXSFC + 1, NYSFC + 1), WAVELEN
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    character(len=2) :: SFCTYPE
    integer IMU, IMU1, IPHI, JMU, JPHI, K1, K2
    real OPI, REFLECT(4, 4), W

    OPI = 1.0/acos(-1.0)
    ! Loop over outgoing discrete ordinates (IMU,IPHI)
    do IMU = NMU/2 + 1, NMU
        IMU1 = IMU - NMU/2
        do IPHI = 1, NPHI0(IMU)
            ! Loop over incident discrete ordinates (JMU,JPHI)
            do JMU = 1, NMU/2
                do JPHI = 1, NPHI0(JMU)
                    ! Get the BRDF for this
                    call shd_surf__brdfSurface(SFCTYPE(2:2), SFCPARMS(2, 1, 1),&
                        &WAVELEN, MU(IMU), PHI(IMU, IPHI),&
                        &MU(JMU), PHI(JMU, JPHI), NSTOKES, REFLECT)
                    W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
                    SFC_BRDF_DO(:, IMU1, IPHI, :, JMU, JPHI)&
                        &= W*REFLECT(1:NSTOKES, 1:NSTOKES)
                end do
            end do
        end do
    end do
    return
end subroutine shd_solv__makeSurfDO

subroutine shd_solv__SHtoDO(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
    &SHPTR, INDATA, OUTDATA)

    !   Transforms the input data from spherical harmonic space to discrete
    ! ordinate space for one zenith angle (IMU).  Successive transforms are
    ! done in zenith and then azimuth angle.  The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(2*MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small primes
    ! (powers of 2 are best).  The spherical harmonic truncation of the
    ! input array is adaptive (with pointer SHPTR); the last term for
    ! each grid point must be the maximum m for the particular l term.

    implicit none
    integer NPTS, NSTOKES, NSTLEG
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, SHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NSTOKES, *), OUTDATA(NSTOKES, NPHI0MAX, NPTS)
    real CMU1(NSTLEG, NLM, NMU), CPHI1(-16:16, 32, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, J, K, L, M, M2, N, IPHI, MS, ME, MR, IS, NSH
    real SUM1
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:, :), SUMCS(:, :), TMP(:)

    allocate (MOFJ(NLM), TMP(2*MM + 2))
    allocate (SUMUV(NSTOKES, -MM:MM), SUMCS(NSTOKES, -MM:MM))

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        IS = SHPTR(I)
        NSH = SHPTR(I + 1) - IS
        if (NSH .eq. 0) then
            do K = 1, NPHI0
                OUTDATA(:, K, I) = 0.0
            end do
        else
            ! Figure the max Fourier azimuthal mode we can do for this Nphi0
            ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
            MS = -ME
            ! First do mu part of real generalized spherical harmonics
            !   transform by summing over l for each m.
            SUMUV(:, :) = 0.0
            do J = 1, NSH
                M = MOFJ(J)
                SUMUV(1, M) = SUMUV(1, M) + CMU1(1, J, IMU)*INDATA(1, IS + J)
                if (NSTOKES .gt. 1) then
                    SUMUV(2, M) = SUMUV(2, M)&
                        &+ CMU1(2, J, IMU)*INDATA(2, IS + J)&
                        &+ CMU1(5, J, IMU)*INDATA(3, IS + J)
                    SUMUV(3, M) = SUMUV(3, M)&
                        &+ CMU1(6, J, IMU)*INDATA(2, IS + J)&
                        &+ CMU1(3, J, IMU)*INDATA(3, IS + J)
                end if
                if (NSTOKES .eq. 4) then
                    SUMUV(4, M) = SUMUV(4, M) + CMU1(4, J, IMU)*INDATA(4, IS + J)
                end if
            end do
            SUMCS(:, 0) = SUMUV(:, 0)
            do M = 1, ME
                SUMCS(1, M) = SUMUV(1, M) + SUMUV(1, -M)
                SUMCS(1, -M) = SUMUV(1, -M) - SUMUV(1, M)
            end do
            if (NSTOKES .gt. 1) then
                do M = 1, ME
                    SUMCS(2, M) = SUMUV(2, M) + SUMUV(2, -M)
                    SUMCS(2, -M) = SUMUV(2, -M) - SUMUV(2, M)
                    SUMCS(3, -M) = SUMUV(3, M) - SUMUV(3, -M)
                    SUMCS(3, M) = SUMUV(3, M) + SUMUV(3, -M)
                end do
            end if
            if (NSTOKES .eq. 4) then
                do M = 1, ME
                    SUMCS(4, -M) = SUMUV(4, M) - SUMUV(4, -M)
                    SUMCS(4, M) = SUMUV(4, M) + SUMUV(4, -M)
                end do
            end if

            ! Then do the Fourier transform from m to phi for each mu
            if (FFTFLAG(IMU)) then
                do N = 1, NSTOKES
                    ! For sines and cosines use a real FFT, and rearrange terms
                    TMP(1) = SUMCS(N, 0)
                    MR = 2
                    do M = 1, ME
                        TMP(MR) = 0.5*SUMCS(N, M)
                        TMP(MR + 1) = -0.5*SUMCS(N, -M)
                        MR = MR + 2
                    end do
                    do M2 = MR, NPHI0
                        TMP(M2) = 0.0
                    end do
                    call DFFTB(NPHI0, TMP, WSAVE(1, IMU))
                    do K = 1, NPHI0
                        OUTDATA(N, K, I) = TMP(K)
                    end do
                end do

            else
                ! Else do a slow DFT
                do IPHI = 1, NPHI0
                    do N = 1, NSTOKES
                        SUM1 = 0.0
                        do M = MS, ME
                            SUM1 = SUM1 + CPHI1(M, IPHI, IMU)*SUMCS(N, M)
                        end do
                        OUTDATA(N, IPHI, I) = SUM1
                    end do
                end do
            end if
        end if
    end do
    deallocate (MOFJ, SUMUV, SUMCS, TMP)
    return
end subroutine shd_solv__SHtoDO

subroutine shd_solv__SHtoDO_(NPTS, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
    &SHPTR, INDATA, OUTDATA)
    !   Transforms the input data from spherical harmonic space to discrete
    ! ordinate space for one zenith angle (IMU) for the unpolarized case
    ! Successive transforms are done in zenith and then azimuth angle.
    ! The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(2*MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small primes
    ! (powers of 2 are best).  The spherical harmonic truncation of the
    ! input array is adaptive (with pointer SHPTR); the last term for
    ! each grid point must be the maximum m for the particular l term.
    ! The number of floating point operations for the DFT is 2*Nmu*Nphi0*Nm.
    ! The number of floating point operations for the zenith angle transform
    ! is 2*Nmu*Nlm.
    implicit none
    integer NPTS
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, SHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(*), OUTDATA(NPHI0MAX, NPTS)
    real CMU1(NLM, NMU), CPHI1(-16:16, 32, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, J, K, L, M, M2, N, IPHI, MS, ME, MR, IS, NSH
    real SUM1
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:), SUMCS(:), TMP(:)

    allocate (MOFJ(NLM), SUMUV(-MM:MM), SUMCS(-MM:MM), TMP(2*MM + 2))

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        IS = SHPTR(I)
        NSH = SHPTR(I + 1) - IS
        if (NSH .eq. 0) then
            do K = 1, NPHI0
                OUTDATA(K, I) = 0.0
            end do
        else
            ! Figure the max Fourier azimuthal mode we can do for this Nph0
            ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
            MS = -ME
            ! First do mu part of real generalized spherical harmonics
            !   transform by summing over l for each m.
            SUMUV(:) = 0.0
            do J = 1, NSH
                M = MOFJ(J)
                SUMUV(M) = SUMUV(M) + CMU1(J, IMU)*INDATA(IS + J)
            end do
            SUMCS(0) = SUMUV(0)
            do M = 1, ME
                SUMCS(M) = SUMUV(M) + SUMUV(-M)
                SUMCS(-M) = SUMUV(-M) - SUMUV(M)
            end do

            ! Then do the Fourier transform from m to phi for each mu
            if (FFTFLAG(IMU)) then
                ! For sines and cosines use a real FFT, and rearrange terms
                TMP(1) = SUMCS(0)
                MR = 2
                do M = 1, ME
                    TMP(MR) = 0.5*SUMCS(M)
                    TMP(MR + 1) = -0.5*SUMCS(-M)
                    MR = MR + 2
                end do
                do M2 = MR, NPHI0
                    TMP(M2) = 0.0
                end do
                call DFFTB(NPHI0, TMP, WSAVE(1, IMU))
                do K = 1, NPHI0
                    OUTDATA(K, I) = TMP(K)
                end do
            else
                ! Else do a slow DFT
                do IPHI = 1, NPHI0
                    SUM1 = 0.0
                    do M = MS, ME
                        SUM1 = SUM1 + CPHI1(M, IPHI, IMU)*SUMCS(M)
                    end do
                    OUTDATA(IPHI, I) = SUM1
                end do
            end if
        end if
    end do
    deallocate (MOFJ, SUMUV, SUMCS, TMP)
    return
end subroutine shd_solv__SHtoDO_

subroutine shd_solv__DOtoSH(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    &RSHPTR, INDATA, OUTDATA)

    !   Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See shd_solv__SHtoDO for more).

    implicit none
    integer NPTS, NSTOKES, NSTLEG
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, RSHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NSTOKES, NPHI0MAX, NPTS), OUTDATA(NSTOKES, *)
    real CMU2(NSTLEG, NMU, NLM), CPHI2(32, -16:16, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, IPHI, J, K, L, M, MS, ME, MR, N, IS, NSH
    real SUM1, DELPHI
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:, :), SUMCS(:, :), TMP(:)

    allocate (MOFJ(NLM), TMP(NPHI0))
    allocate (SUMUV(NSTOKES, -MM:MM), SUMCS(NSTOKES, -MM:MM))

    DELPHI = 2.0*acos(-1.0)/NPHI0

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    IS = 0
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        NSH = RSHPTR(I + 1) - RSHPTR(I)
        ! Figure the max Fourier azimuthal mode we can do for this Nphi0
        ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
        MS = -ME
        ! First do Fourier transform from phi to m for each mu
        if (FFTFLAG(IMU)) then
            do N = 1, NSTOKES
                do K = 1, NPHI0
                    TMP(K) = INDATA(N, K, I)
                end do
                ! For sines and cosines use a real FFT, and rearrange terms
                call DFFTF(NPHI0, TMP, WSAVE(1, IMU))
                SUMCS(N, 0) = TMP(1)*DELPHI
                MR = 2
                do M = 1, ME
                    SUMCS(N, M) = TMP(MR)*DELPHI
                    SUMCS(N, -M) = -TMP(MR + 1)*DELPHI
                    MR = MR + 2
                end do
            end do
        else
            ! Else do a slow DFT
            do M = MS, ME
                do N = 1, NSTOKES
                    SUM1 = 0.0
                    do IPHI = 1, NPHI0
                        SUM1 = SUM1 + CPHI2(IPHI, M, IMU)*INDATA(N, IPHI, I)
                    end do
                    SUMCS(N, M) = SUM1
                end do
            end do
        end if
        ! Then do mu part of real generalized spherical harmonics
        !   transform by summing over l for each m.
        SUMUV(:, ME + 1:MM) = 0.0
        SUMUV(:, -MM:MS - 1) = 0.0
        SUMUV(:, 0) = SUMCS(:, 0)
        do M = 1, ME
            SUMUV(1, M) = SUMCS(1, M) - SUMCS(1, -M)
            SUMUV(1, -M) = SUMCS(1, M) + SUMCS(1, -M)
        end do
        if (NSTOKES .gt. 1) then
            do M = 1, ME
                SUMUV(2, M) = SUMCS(2, M) - SUMCS(2, -M)
                SUMUV(2, -M) = SUMCS(2, M) + SUMCS(2, -M)
                SUMUV(3, M) = SUMCS(3, M) + SUMCS(3, -M)
                SUMUV(3, -M) = SUMCS(3, M) - SUMCS(3, -M)
            end do
        end if
        if (NSTOKES .eq. 4) then
            do M = 1, ME
                SUMUV(4, M) = SUMCS(4, M) + SUMCS(4, -M)
                SUMUV(4, -M) = SUMCS(4, M) - SUMCS(4, -M)
            end do
        end if
        do J = 1, NSH
            M = MOFJ(J)
            OUTDATA(1, IS + J) = OUTDATA(1, IS + J) + CMU2(1, IMU, J)*SUMUV(1, M)
            if (NSTOKES .gt. 1) then
                OUTDATA(2, IS + J) = OUTDATA(2, IS + J) + CMU2(2, IMU, J)*SUMUV(2, M)&
                    &+ CMU2(5, IMU, J)*SUMUV(3, M)
                OUTDATA(3, IS + J) = OUTDATA(3, IS + J) + CMU2(6, IMU, J)*SUMUV(2, M)&
                    &+ CMU2(3, IMU, J)*SUMUV(3, M)
            end if
            if (NSTOKES .eq. 4) then
                OUTDATA(4, IS + J) = OUTDATA(4, IS + J) + CMU2(4, IMU, J)*SUMUV(4, M)
            end if
        end do
        IS = IS + NSH
    end do
    deallocate (MOFJ, SUMUV, SUMCS, TMP)
    return
end subroutine shd_solv__DOtoSH

subroutine shd_solv__DOtoSH_(NPTS, ML, MM, NLM,&
    &NMU, NPHI0MAX, NPHI0,&
    &IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    &RSHPTR, INDATA, OUTDATA)

    !   Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See shd_solv__SHtoDO for more).

    implicit none
    integer NPTS
    integer ML, MM, NLM, NMU, NPHI0MAX, NPHI0
    integer IMU, RSHPTR(NPTS + 1)
    logical FFTFLAG(NMU)
    real INDATA(NPHI0MAX, NPTS), OUTDATA(*)
    real CMU2(NMU, NLM), CPHI2(32, -16:16, NMU)
    real WSAVE(3*NPHI0MAX + 15, NMU)
    integer I, IPHI, J, K, L, M, MS, ME, MR, N, IS, NSH
    real SUM1, DELPHI
    integer, allocatable :: MOFJ(:)
    real, allocatable :: SUMUV(:), SUMCS(:), TMP(:)

    allocate (MOFJ(NLM), SUMUV(-MM:MM), SUMCS(-MM:MM), TMP(NPHI0))

    DELPHI = 2.0*acos(-1.0)/NPHI0
    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = -min(L, MM), min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    IS = 0
    do I = 1, NPTS
        ! Determine the current SH truncation to use
        NSH = RSHPTR(I + 1) - RSHPTR(I)
        SUMCS(:) = 0.0
        ! Figure the max Fourier azimuthal mode we can do for this Nphi0
        ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
        MS = -ME
        ! First do Fourier transform from phi to m for each mu
        if (FFTFLAG(IMU)) then
            do K = 1, NPHI0
                TMP(K) = INDATA(K, I)
            end do
            ! For sines and cosines use a real FFT, and rearrange terms
            call DFFTF(NPHI0, TMP, WSAVE(1, IMU))
            SUMCS(0) = TMP(1)*DELPHI
            MR = 2
            do M = 1, ME
                SUMCS(M) = TMP(MR)*DELPHI
                SUMCS(-M) = -TMP(MR + 1)*DELPHI
                MR = MR + 2
            end do
        else
            ! Else do a slow DFT
            do M = MS, ME
                SUM1 = 0.0
                do IPHI = 1, NPHI0
                    SUM1 = SUM1 + CPHI2(IPHI, M, IMU)*INDATA(IPHI, I)
                end do
                SUMCS(M) = SUM1
            end do
        end if
        ! Then do mu part of real generalized spherical harmonics
        !   transform by summing over l for each m.
        SUMUV(0) = SUMCS(0)
        do M = 1, MM
            SUMUV(M) = SUMCS(M) - SUMCS(-M)
            SUMUV(-M) = SUMCS(M) + SUMCS(-M)
        end do
        do J = 1, NSH
            M = MOFJ(J)
            OUTDATA(IS + J) = OUTDATA(IS + J) + CMU2(IMU, J)*SUMUV(M)
        end do
        IS = IS + NSH
    end do
    deallocate (MOFJ, SUMUV, SUMCS, TMP)
    return
end subroutine shd_solv__DOtoSH_

subroutine shd_solv__splitGrid(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
    &DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
    &NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, NX, XGRID, NY, YGRID,&
    &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
    &TEMP, PLANCK, DIRFLUX, SRCTYPE,&
    &SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    &UNITS, WAVENO, WAVELEN, ADAPTCRIT, ADAPTIND)

    !   Splits the cells that have a cell dividing criterion greater than
    ! CURSPLITACC if DOSPLIT is true.  The current splitting criterion
    ! achieved is returned in SPLITCRIT.  If we are at the end of the
    ! grid point, grid cell, or spherical harmonic arrays (sizes given
    ! by MAXIG, MAXIC, MAXIV, MAXIDO) then the OUTOFMEM flag is returned true.

    use shd_grid
    use shd_ext_mpi

    implicit none
    integer MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX
    integer NPTS, NCELLS, BCFLAG, IPFLAG
    integer NX, NY, NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    real ADAPTIND(*)
    integer*2 CELLFLAGS(*), IPHASE(*)
    logical DELTAM, ACCELFLAG, DOSPLIT, OUTOFMEM
    real CURSPLITACC, SPLITCRIT
    real SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(NSTLEG, NLM)
    real WAVENO(2), WAVELEN
    real GRIDPOS(3, *), XGRID(*), YGRID(*), ADAPTCRIT(*)
    real EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    real RADIANCE(NSTOKES, *), SOURCE(NSTOKES, *)
    real TEMP(*), PLANCK(*), DIRFLUX(*)
    character(len=1) :: SRCTYPE, UNITS
    integer ICELL, ICELL1, IDIR, I, N, NEWPOINTS(3, 4), OLDNPTS
    integer MAXCELLS, MAXPTS, MAXWORK, MAXSH, MAXRAD
    logical OUTOFMEM0
    real ADAPT(3), FRAC

    if (DOSPLIT) then
        OLDNPTS = NPTS
        ! Do the cells in batchs, so the newly split cells can be tested
        !   for further splitting
        ICELL1 = 1
        do while (ICELL1 .le. NCELLS)
            ! For the current batch of cells, compute the adaptive grid
            !   splitting criterion for 3 directions for each cell.
            N = 0
            do ICELL = ICELL1, NCELLS
                if (TREEPTR(2, ICELL) .eq. 0) then
                    ! If this is a end node cell then compute adaptive criterion
                    call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
                        &SHPTR, NSTOKES, SOURCE, ICELL,&
                        &ADAPT, ADAPTCRIT(N + 1), IDIR)
                    ADAPTIND(N+1) = 4*ICELL + IDIR
                    N = N + 1
                end if
            end do
            ICELL1 = NCELLS + 1

            ! Sort the adaptive criterion in decreasing order
            call shd_comp__ssort(ADAPTCRIT, ADAPTIND, N, -2)

            ! Figure how many cells we can split, leaving some room
            !   for cells split during grid smoothing
            FRAC = 0.03
            MAXCELLS = MAXIC - FRAC*(MAXIC - NCELLS) - 2
            MAXPTS = MAXIG - FRAC*(MAXIG - NPTS) - 4
            MAXWORK = MAXIDO - FRAC*(MAXIDO - NPHI0MAX*NPTS) - 4*NPHI0MAX
            MAXSH = MAXIV - FRAC*(MAXIV - SHPTR(NPTS + 1)) - 4*NLM
            MAXRAD = MAXIV + MAXIG - FRAC*(MAXIV + MAXIG - RSHPTR(NPTS + 1)) - 4*NLM

            ! Go down the list from highest adaptive cell criterion on down
            OUTOFMEM0 = OUTOFMEM
            I = 1
            do while (ADAPTCRIT(I) .gt. CURSPLITACC .and. .not. OUTOFMEM0&
                &.and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    &NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    &SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    &RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM0 = .true.
                else
                    ICELL = int(ADAPTIND(I))/4
                    IDIR = iand(int(ADAPTIND(I)), 3)
                    ! Create the new grid points and make pointers
                    call shd_grid__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                    ! Interpolate the medium properties, radiance, and source
                    call shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
                        &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                        &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                        &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                        &UNITS, WAVENO, WAVELEN,&
                        &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                        &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                end if
                I = I + 1
            end do

            MAXCELLS = MAXIC - 2
            MAXPTS = MAXIG - 4
            MAXWORK = MAXIDO - 4*NPHI0MAX
            MAXSH = MAXIV - 4*NLM
            MAXRAD = MAXIV + MAXIG - 4*NLM

            ! Split some more cells to make the adaptive grid "smoother"
            do while (.not. OUTOFMEM .and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    &NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    &SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    &RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM = .true.
                else
                    ! See if cell should be split to make better grid:
                    !   IDIR=0 if no, IDIR=direction (1,2,3) otherwise
                    ICELL = int(ADAPTIND(I))/4
                    call shd_grid__testGridSmooth(GRIDPTR, GRIDPOS, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, SHPTR, ICELL, IDIR)
                    ! Set IDIR=0 to stop adaptive grid smoothing
                    ! IDIR = 0
                    if (IDIR .gt. 0) then
                        ! Create the new grid points and make pointers
                        call shd_grid__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                            &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                        ! Interpolate the medium properties, radiance, and source
                        call shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
                            &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                            &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                            &UNITS, WAVENO, WAVELEN,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                            &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                    end if
                end if
                I = I + 1
            end do

        end do
        if (OUTOFMEM0) then
            OUTOFMEM = .true.
            write (6, '(2A,I8,A,I8)')&
                &'Message [_splitGrid]: Further cell splitting requires larger array allocations: ',&
                &'Ncells = ', NCELLS, ', Npnts = ', NPTS
        end if

        ! If in parallel mode, calculate all the direct beam fluxes at once
        if (SRCTYPE .ne. 'T' .and.&
            &(btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) then
            call shd_ext__mpiMakeDirect(OLDNPTS + 1, NPTS, BCFLAG,&
                &IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
                &SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
                &NX, XGRID, NY, YGRID, DIRFLUX)
        end if
    end if

    ! Find the max adaptive cell criterion after the cell divisions
    SPLITCRIT = 0.0
    do ICELL = 1, NCELLS
        if (TREEPTR(2, ICELL) .eq. 0) then
            call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
                &SHPTR, NSTOKES, SOURCE, ICELL,&
                &ADAPT, ADAPTCRIT(1), IDIR)
            SPLITCRIT = max(SPLITCRIT, ADAPTCRIT(1))
        end if
    end do

    return

end subroutine shd_solv__splitGrid
