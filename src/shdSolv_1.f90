subroutine SOLVE_RTE(&
    & NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
    & ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE,&
    & NPHI0, MU, PHI, WTDO,&
    & MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG,&
    & BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD,&
    & SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    & SFCTYPE, GNDTEMP, GNDALBEDO,&
    & NXSFC, NYSFC, DELXSFC, DELYSFC,&
    & NSFCPAR, SFCPARMS, SFCGRIDPARMS,&
    & UNITS, WAVENO, WAVELEN,&
    & ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,&
    & SPLITACC, SHACC, XGRID, YGRID, ZGRID,&
    & TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,&
    & MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    & NPTS, GRIDPOS,&
    & NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,&
    & SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)

    ! Performs the SHDOM solution procedure.
    ! Output is returned in SOURCE, RADIANCE, FLUXES, DIRFLUX.
    ! The adaptive grid stucture (NCELLS,GRIDPTR,NEIGHPTR,TREEPTR,
    ! CELLFLAGS,NPTS,GRIDPOS) is input and modified for new grid
    ! cells and points on output.
    implicit none
    integer NSTOKES, NX, NY, NZ, NX1, NY1, NXSFC, NYSFC, NSFCPAR
    integer ML, MM, NLM, NMU, NPHI, NANG, NLEG, NSTLEG, NUMPHASE
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
    character SRCTYPE*1, UNITS*1, SFCTYPE*2

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

    ! Allocate local arrays
    ! NSTLEG=1 for NSTOKES=1, otherwise NSTLEG=6
    allocate (CMU1(NSTLEG, NLM, NMU), CMU2(NSTLEG, NLM, NMU))
    allocate (CPHI1(-16:16, 32, NMU), CPHI2(32, -16:16, NMU))
    allocate (YLMSUN(NSTLEG, NLM))
    allocate (FFTFLAG(NMU), WTMU(NMU), WPHISAVE(NMU*(3*NPHI + 15)))

    ! Set up some things before solution loop

    ! Transfer the medium properties to the internal grid and add
    call INTERP_GRID(NPTS, NSTLEG, NLEG, GRIDPOS,&
        & TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)

    ! If Delta-M then scale the extinction, albedo, and Legendre ter
    ! Put the Planck blackbody source in PLANCK.
    call PREPARE_PROP(ML, MM, NSTLEG, NLEG, NPTS, DELTAM, NUMPHASE,&
        & SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX,&
        & EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE)
    ! Get the maximum single scattering albedo over all processors
    call TOTAL_ALBEDO_MAX(ALBMAX)

    ! Compute the solar transmission in DIRFLUX.
    ! Precompute Ylm's for solar direction
    if (SRCTYPE .ne. 'T') then
        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            call MAKE_DIRECT_PAR(1, NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                & NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                & GRIDPOS, NX, XGRID, NY, YGRID, DIRFLUX)
        else
            call MAKE_DIRECT(NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                & NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                & GRIDPOS, DIRFLUX)
        end if
        call YLMALL(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG,&
            & YLMSUN)
    end if

    ! Make the discrete ordinates (angles)
    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
    ORDINATESET = 2
    call MAKE_ANGLE_SET(NMU, NPHI, ORDINATESET, NPHI0MAX,&
        & NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the Ylm transform coefficients
    call MAKE_SH_DO_COEF(NSTLEG, ML, MM, NLM, NMU, NPHI0,&
        & NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,&
        & CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

    FIXSH = .false.
    if (.not. INRADFLAG) then
        ! Initialize the radiance on the base grid using Eddington
        ! two-stream plane-parallel
        call INIT_RADIANCE(NSTOKES, NX1*NY1, NZ, NSTLEG, NLEG,&
            & RSHPTR, ZGRID,&
            & EXTINCT, ALBEDO, LEGEN, TEMP, NUMPHASE, IPHASE,&
            & SRCTYPE, SOLARFLUX, SOLARMU, GNDALBEDO, GNDTEMP, SKYRAD,&
            & UNITS, WAVENO, WAVELEN, RADIANCE)
        ! Interpolate the radiance on the non-base grid points from th
        !   base grid points
        NBPTS = NX1*NY1*NZ
        NBCELLS = (NZ - 1)*(NX + ibits(BCFLAG, 0, 1) - ibits(BCFLAG, 2, 1))&
            & *(NY + ibits(BCFLAG, 1, 1) - ibits(BCFLAG, 3, 1))
        call INTERP_RADIANCE(NSTOKES, NBPTS, NPTS, RSHPTR, RADIANCE,&
            & NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        ! Initialize the source function from the radiance field
        call COMPUTE_SOURCE(NSTOKES, ML, MM, NLM,&
            & NSTLEG, NLEG, NUMPHASE, NPTS,&
            & FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            & IPHASE, PLANCK, DIRFLUX, SHACC,&
            & MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            & .true., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)
    end if

    if (ACCELFLAG) then
        OSHPTR(1:NPTS + 1) = SHPTR(1:NPTS + 1)
        DELSOURCE(:, 1:OSHPTR(NPTS + 1)) = 0.0
    end if
    SOLCRIT = 1.0

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'
    ! Make the pointers to the grid points on the top and bottom bou
    call BOUNDARY_PNTS(NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
        & ZGRID(1), ZGRID(NZ), GRIDPOS,&
        & NTOPPTS, NBOTPTS, BCPTR)
    if (SFCTYPE(1:1) .eq. 'V') then
        ! If this is a variable surface then bilinearly interpolate
        !   the surface temperature and reflection parameters to the
        !   bottom grid points.
        call SURFACE_PARM_INTERP(NBOTPTS, BCPTR(1, 2), GRIDPOS,&
            & SRCTYPE, WAVENO, WAVELEN, UNITS,&
            & NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
            & SFCGRIDPARMS)
        call CHECK_UNIFORM_SFC(NXSFC, NYSFC, NSFCPAR, SFCPARMS,&
            & UNIFORM_SFC_BRDF)
        if (UNIFORM_SFC_BRDF) then
            if (NSTOKES*(NMU/2)*NPHI0MAX .gt.&
                & sqrt(1.0*huge(NMU)/kind(NMU))) then
                IERR = 1
            else
                allocate (SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX,&
                    & NSTOKES, NMU/2, NPHI0MAX), STAT=IERR)
            end if
            if (IERR .eq. 0) then
                call MAKE_SFC_BRDF_DO_MATRIX(SFCTYPE, WAVELEN,&
                    & NXSFC, NYSFC, NSFCPAR, SFCPARMS, NSTOKES,&
                    & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                    & SFC_BRDF_DO)
            else
                print '(A,A)', 'Not enough memory to allocate precomputed',&
                    & ' surface BRDF discrete ordinate matrix: '
                print *, '  Proceeding more slowly.'
                UNIFORM_SFC_BRDF = .false.
            end if
        end if
    end if

    ITER = 0
    if (MAXITER .le. ITER) return
    ! Starting values for the adaptive cell splitting controlling me
    ! The splitting accuracy used at a particular iteration is CURSP
    ! and it is started out high an lowered to the desired SPLITACC.
    ! The cell splitting starts when the solution criterion (SOLCRIT
    ! goes below STARTADAPTSOL and ends when the solution criterion
    ! goes below ENDADAPTSOL, though the current splitting accuracy
    ! is designed to get to the final SPLITACC over a smaller range
    ! (ADAPTRANGE).  The currect splitting accuracy is computed from
    ! a geometric average of the solution criterion (AVGSOLCRIT).
    ! The SPLIT_GRID routine is called and the splitting criterion r
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

    write (6, *) '! Iter Log(Sol)  SplitCrit  Npoints  Nsh(avg)'

    ! Main solution loop
    !   Iterate until the iterations exceed the limit, or the desire
    !   solution criterion is reached and the current splitting accu
    !   is down to the desired level.
    do while (ITER .lt. MAXITER .and. (SOLCRIT .gt. SOLACC&
        & .or. (SPLITCRIT .gt. SPLITACC .and.&
        & CURSPLITACC .gt. SPLITACC .and. .not. OUTOFMEM)))
        ITER = ITER + 1

        ! Don't do the adaptive grid cell splitting stuff unless wante
        if (SPLITACC .gt. 0.0) then
            ! If in the right solution range then try to split some cell
            !   interpolate the medium and source functions to the new p
            ! The current splitting accuracy is started high and then lo
            AVGSOLCRIT = sqrt(AVGSOLCRIT*SOLCRIT)
            DOSPLIT = SOLCRIT .le. STARTADAPTSOL .and.&
                & (SOLCRIT .gt. ENDADAPTSOL .or. CURSPLITACC .gt. SPLITACC)&
                & .and. .not. OUTOFMEM
            ! Make sure all processors are going to split cells if any wan
            call UNIFY_SPLITTING(DOSPLIT, STARTSPLITACC)
            BETA = log(STARTSPLITACC/SPLITACC)/log(ADAPTRANGE)
            CURSPLITACC = min(CURSPLITACC, max(SPLITACC,&
                & SPLITACC*(AVGSOLCRIT/(3.0*ENDADAPTSOL))**BETA))
            if (SOLCRIT .le. ENDADAPTSOL) CURSPLITACC = SPLITACC
            if (SPLITTESTING) then
                call SPLIT_GRID(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
                    & DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
                    & NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    & GRIDPOS, NX, XGRID, NY, YGRID,&
                    & NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                    & DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
                    & EXTINCT, ALBEDO, LEGEN, IPHASE,&
                    & RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
                    & TEMP, PLANCK, DIRFLUX, SRCTYPE,&
                    & SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                    & UNITS, WAVENO, WAVELEN, WORK, WORK2)
                if (SOLCRIT .gt. STARTADAPTSOL) STARTSPLITACC = SPLITCRIT
            end if
            ! Find the maximum splitting criterion over all processors
            call TOTAL_SPLITCRIT_MAX(SPLITCRIT)
            if (SOLCRIT .le. ENDADAPTSOL) SPLITTESTING = .false.
        end if

        ! Find the SH truncation to use for the radiance (RSHPTR)
        call RADIANCE_TRUNCATION(HIGHORDERRAD,&
            & NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE,&
            & NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
            & MAXIV + MAXIG, FIXSH, SHACC, RSHPTR)

        ! Integrate the source function along discrete ordinates
        !   to compute the radiance field.  The input source function
        !   and output radiance are in spherical harmonic expansions,
        !   and the SH/DO transforms are done during the loop over
        !   discrete ordinates.
        call PATH_INTEGRATION(NX, NY, NZ, NSTOKES, NSTLEG, ML, MM, NLM,&
            & NMU, NPHI0MAX, NPHI0, NANG, BCFLAG, IPFLAG,&
            & NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
            & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
            & SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN,&
            & UNITS, SFCTYPE, GNDTEMP, GNDALBEDO,&
            & NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
            & SFCGRIDPARMS, UNIFORM_SFC_BRDF, SFC_BRDF_DO,&
            & MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
            & MU, PHI, WTDO,&
            & CMU1, CPHI1, CMU2, CPHI2, WPHISAVE, FFTFLAG,&
            & DIRFLUX, FLUXES, EXTINCT,&
            & SHPTR, SOURCE, RSHPTR, RADIANCE,&
            & WORK, WORK1, WORK2)

        ! Compute the source function from the radiance field,
        !   do the adaptive spherical harmonics truncation, compute
        !   the solution criterion, and dot products for acceleration
        if (SOLCRIT .lt. ENDADAPTSOL .or. ITER .gt. 30) FIXSH = .true.
        call COMPUTE_SOURCE(NSTOKES, ML, MM, NLM,&
            & NSTLEG, NLEG, NUMPHASE, NPTS,&
            & FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            & IPHASE, PLANCK, DIRFLUX, SHACC,&
            & MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            & .false., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

        ! Calculate the acceleration parameter and solution criterion
        !   from all the processors
        call CALC_ACCEL_SOLCRIT(ACCELFLAG, DELJDOT, DELJOLD, DELJNEW,&
            & JNORM, ACCELPAR, SOLCRIT)

        ! Accelerate the convergence of the source function vector
        call ACCELERATE_SOLUTION(ACCELPAR, NPTS, SHPTR, OSHPTR,&
            & NSTOKES, SOURCE, DELSOURCE)

        ! If it is a not scattering medium then do only one iteration
        if (ALBMAX .lt. SOLACC) SOLCRIT = SOLACC

        ! Print out the log solution criterion, number of points,
        !   and average SH truncation.
        write (6, '(2X,I4,F8.3,1X,E10.3,1X,I8,1X,F8.2,1X,F6.3)')&
            & ITER, log10(max(SOLCRIT, 1.0e-20)), SPLITCRIT,&
            & NPTS, FLOAT(SHPTR(NPTS + 1))/NPTS,&
            & FLOAT(SHPTR(NPTS + 1))/(NPTS*NLM)
    end do

    write (6, '(1X,A,I6,A,F9.6)') '! Iterations: ', ITER,&
        & '     Final Criterion: ', SOLCRIT

    if (SFCTYPE(1:1) .eq. 'V' .and. UNIFORM_SFC_BRDF) then
        deallocate (SFC_BRDF_DO)
    end if
    deallocate (YLMSUN, FFTFLAG, WTMU, WPHISAVE)
    deallocate (CMU1, CMU2, CPHI1, CPHI2)

    ! Comment in for GLE graphics output of cell structure (see rou
    ! call VISUALIZE_CELLS('cellfinal.gle', YGRID(1), 1, IPFLAG, &
    !     & NX, XGRID(NX+1)-XGRID(1), ZGRID(NZ)-ZGRID(1), &
    !     & NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, &
    !     & CELLFLAGS, GRIDPOS)

    ! Comment in for output of cell splitting criterion for every c
    ! call OUTPUT_CELL_SPLIT ('cellsplit.dat',GRIDPTR, &
    !     & GRIDPOS, EXTINCT, SHPTR, SOURCE, NCELLS)

    return
end subroutine

subroutine CALC_SOURCE_PNT(&
    & NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
    & SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    & ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    & SOURCET)
    ! Calculates the polarized source function in real generalized
    ! spherical harmonics at one grid point from the radiance in
    ! real generalized spherical harmonics and thermal and/or solar
    ! source information.
    implicit none
    integer NSTOKES, NSTLEG, NLEG, NLM, LOFJ(NLM), IPH, NR
    character SRCTYPE*1
    real FLUX0, YLMSUN(NSTLEG, NLM), PLANCK, ALBEDO
    real LEGEN(NSTLEG, 0:NLEG, *)
    real RADIANCE(NSTOKES, NR)
    real SOURCET(NSTOKES, NLM)
    real C
    parameter(C=3.544907703)
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
                & + LEGEN(2, LOFJ(J), IPH)*RADIANCE(2, J)
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(3, LOFJ(J), IPH)*RADIANCE(3, J)
        end do
    end if
    if (NSTOKES .eq. 4) then
        do J = 1, NR
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(6, LOFJ(J), IPH)*RADIANCE(4, J)
            SOURCET(4, J) = SOURCET(4, J) - LEGEN(6, LOFJ(J), IPH)*RADIANCE(3, J)&
                & + LEGEN(4, LOFJ(J), IPH)*RADIANCE(4, J)
        end do
    end if
    ! Do the final multiplication by the single scattering albedo
    SOURCET(:, :) = ALBEDO*SOURCET(:, :)

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        do J = 1, NLM
            SOURCET(1, J) = SOURCET(1, J)&
                & + FLUX0*ALBEDO*LEGEN(1, LOFJ(J), IPH)*YLMSUN(1, J)
        end do
        if (NSTOKES .gt. 1) then
            do J = 5, NLM
                SOURCET(2, J) = SOURCET(2, J)&
                    & + FLUX0*ALBEDO*LEGEN(5, LOFJ(J), IPH)*YLMSUN(1, J)
            end do
        end if
    end if
    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1, 1) = SOURCET(1, 1) + C*PLANCK
    end if

    return
end subroutine

subroutine CALC_SOURCE_PNT_UNPOL(&
    & NLM, NLEG, LOFJ,&
    & SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    & ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    & SOURCET)
    ! Calculates the scalar source function in spherical harmonics at
    ! one grid point from the radiance in spherical harmonics and
    ! thermal and/or solar source information.
    implicit none
    integer NLEG, NLM, LOFJ(NLM), IPH, NR
    character SRCTYPE*1
    real FLUX0, YLMSUN(NLM), PLANCK, ALBEDO
    real LEGEN(0:NLEG, *)
    real RADIANCE(NR)
    real SOURCET(NLM)
    real C
    parameter(C=3.544907703)
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
end subroutine

subroutine COMPUTE_SOURCE(&
    & NSTOKES, ML, MM, NLM,&
    & NSTLEG, NLEG, NUMPHASE, NPTS, FIXSH,&
    & SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
    & IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,&
    & RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
    & FIRST, ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)
    ! Computes the source function (SOURCE) in spherical harmonic spac
    !   for all the grid points.  The thermal source and/or solar
    !   pseudo-source (in PLANCK or DIRFLUX) is added to the scattering so
    !   (computed from LEGEN and the spherical harmonic expansion in RADIA
    !   Also computes the adaptive spherical harmonic truncation for
    !   the source function, unless FIXSH is true.  The SH truncation is s
    !   so that the higher terms of the source function are smaller in
    !   absolute value than SHACC. The truncation is to the nearest l,
    !   so that a triangular (for M=L) truncation is preserved.  If there
    !   is little or no scattering then the SH truncation can be set to
    !   zero terms.  The source function SH truncation is output in terms
    !   of the pointer SHPTR (note that RSHPTR refers to the RADIANCE arra
    !   and OSHPTR to the DELSOURCE array).
    ! To save memory this routine also computes some things for the
    !   series acceleration: the dot products of the difference in source
    !   function between successive iterations and the new difference sour
    !   function vector (DELSOURCE).  To avoid using another array and
    !   because of the adaptive truncation pointer stuff, this requires
    !   two extra calculations of the source function.
    ! Computes the sum of squares of the new source function for calcu
    !   of the solution criterion (SOLCRIT, the RMS difference in succesiv
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
    character SRCTYPE*1
    integer IS, ISO, IR, I, IPH, J, JS, K, L, LS, M, MS, ME, NS, NR
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)
    real SRCMIN, C, SECMU0, D

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))

    SRCMIN = SHACC
    ! Make source function from scattering integral and real sources
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
                write (6, *) 'COMPUTE_SOURCE: NR>NLM 1', I, IR, RSHPTR(I + 1)
                stop
            end if
            if (NSTOKES .eq. 1) then
                call CALC_SOURCE_PNT_UNPOL(NLM, NLEG, LOFJ,&
                    & SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    & PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    & NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call CALC_SOURCE_PNT(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    & SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    & PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    & NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Compute the dot product of the new and old source differen
            !   and the solution criterion (don't use new points for Jno
            if (ACCELFLAG) then
                IS = SHPTR(I)
                ISO = OSHPTR(I)
                NS = min(SHPTR(I + 1) - IS, OSHPTR(I + 1) - ISO)
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJDOT = DELJDOT +&
                            & (SOURCET(K, J) - SOURCE(K, IS + J))*DELSOURCE(K, ISO + J)
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
                call CALC_SOURCE_PNT_UNPOL(NLM, NLEG, LOFJ,&
                    & SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    & PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    & NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call CALC_SOURCE_PNT(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    & SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    & PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    & NR, RADIANCE(1, IR + 1), SOURCET)
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
            write (6, *) 'COMPUTE_SOURCE: NR>NLM 3', I, IR, RSHPTR(I + 1)
            stop
        end if
        if (NSTOKES .eq. 1) then
            call CALC_SOURCE_PNT_UNPOL(NLM, NLEG, LOFJ,&
                & SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                & PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                & NR, RADIANCE(1, IR + 1), SOURCET)
        else
            call CALC_SOURCE_PNT(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                & SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                & PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                & NR, RADIANCE(1, IR + 1), SOURCET)
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
            write (6, *) 'COMPUTE_SOURCE: MAXIV exceeded', MAXIV
            write (6, *) 'Out of memory for more spherical harmonic terms.'
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
end subroutine

subroutine RADIANCE_TRUNCATION(HIGHORDERRAD,&
    & NSTOKES, ML, MM, NSTLEG, NLEG, NUMPHASE,&
    & NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
    & MAXIR, FIXSH, SHACC, RSHPTR)
    ! Computes the next radiance spherical harmonic series truncation
    ! to use based on the current source function truncation and the
    ! scattering properties.  The radiance truncation is the smaller
    ! of that based on a minimum SH term radiance and that needed for
    ! the computation of the source function.  The truncation is always
    ! made for a complete set in the L index.  The radiance SH truncatio
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
        !   We don't have radiances for any new points, so don't use t
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

        ! The radiance L truncation is the smaller of the radiance cut
        !   from above and a few more than the source truncation, but
        !   it is at least 1 (so low order terms are always available)
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
                & 'RADIANCE_TRUNCATION: Out of memory for more radiance terms.'
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
            write (6, *) 'RADIANCE_TRUNCATION: Really out of memory ',&
                & 'for more radiance terms.'
            write (6, *) 'Increase MAXIV.'
            stop
        end if
        RSHPTR(I + 1) = IR
    end do
    RSHPTR(NPTS + 2) = IR
    return
end subroutine

subroutine ACCELERATE_SOLUTION(ACCELPAR, NPTS, SHPTR, OSHPTR,&
    & NSTOKES, SOURCE, DELSOURCE)
    !   Accelerates the successive order of scattering series using the
    ! input acceleration parameter (ACCELPAR).
    implicit none
    integer NPTS, SHPTR(*), OSHPTR(*), NSTOKES
    real ACCELPAR, SOURCE(NSTOKES, *), DELSOURCE(NSTOKES, *)
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
                    & + ACCELPAR*DELSOURCE(:, ISD + J)
            end do
        end do
    end if

    return
end subroutine

subroutine PATH_INTEGRATION(NX, NY, NZ, NSTOKES, NSTLEG,&
    & ML, MM, NLM, NMU, NPHI0MAX, NPHI0, NANG,&
    & BCFLAG, IPFLAG,&
    & NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN,&
    & UNITS, SFCTYPE, GNDTEMP, GNDALBEDO,&
    & NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
    & SFCGRIDPARMS, UNIFORM_SFC_BRDF, SFC_BRDF_DO,&
    & MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    & MU, PHI, WTDO,&
    & CMU1, CPHI1, CMU2, CPHI2, WSAVE, FFTFLAG,&
    & DIRFLUX, FLUXES, EXTINCT,&
    & SHPTR, SOURCE, RSHPTR, RADIANCE,&
    & WORK, SWEEPORD, GRIDRAD)
    !   Performs the path integrations through the medium specified by
    ! the extinction (EXTINCT) and source function (SOURCE) in
    ! spherical harmonic space.  The source function is transformed to
    ! discrete ordinates one zenith angle at a time.  The integrations
    ! are done from the top boundary (with no reflection) down to the
    ! bottom boundary for the downwelling angles, and then up from the
    ! bottom to the top after computing the surface reflection/emission.
    ! After the path integrations produce the discrete angle radiance
    ! it is transformed and added to output radiance spherical harmonic
    ! expansion.  The hemispheric flux FLUXES is also computed.
    ! The discrete ordinates in MU, PHI must be a complete set with the
    ! downwelling (mu<0) angles first.
    !   The surface reflection is handled differently for Lambertian
    ! surfaces and more general surfaces specified by bidirectional
    ! reflection distribution functions (BRDF).  For the general BRDF,
    ! the bottom boundary radiance must be computed for each upwelling
    ! ordinate, while for the Lambertian surface the boundary radiances
    ! can be computed just once, since they are isotropic.  For Lamberti
    ! cases the BCRAD array holds only the isotropic radiance for the
    ! boundaries at the top (downwelling) and bottom (upwelling) boundar
    ! For the general BRDF, the BCRAD array in addition holds the
    ! downwelling radiances at all discrete ordinates for the bottom
    ! boundary points.
    implicit none
    integer NX, NY, NZ, NSTOKES, NSTLEG, ML, MM, NLM
    integer NMU, NPHI0MAX, NPHI0(NMU), NANG
    integer BCFLAG, IPFLAG
    integer NPTS, NCELLS, MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS
    integer NXSFC, NYSFC, NSFCPAR
    integer SHPTR(NPTS + 1), RSHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    integer SWEEPORD(NPTS, *)
    logical FFTFLAG(NMU), UNIFORM_SFC_BRDF
    real SOLARMU, SOLARAZ, SKYRAD
    real GNDTEMP, GNDALBEDO
    real DELXSFC, DELYSFC
    real SFCPARMS(NSFCPAR, *), SFCGRIDPARMS(NSFCPAR, MAXNBC)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    real WAVENO(2), WAVELEN
    real XGRID(NX), YGRID(NY), ZGRID(NZ), GRIDPOS(3, NPTS)
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real CMU1(*), CPHI1(*), CMU2(*), CPHI2(*), WSAVE(*)
    real BCRAD(NSTOKES, MAXBCRAD)
    real DIRFLUX(NPTS), FLUXES(2, NPTS), EXTINCT(NPTS)
    real GRIDRAD(NSTOKES, NPTS), WORK(NSTOKES, NPHI0MAX, NPTS)
    real SOURCE(NSTOKES, *), RADIANCE(NSTOKES, *)
    character SRCTYPE*1, UNITS*1, SFCTYPE*2

    integer I, I1, K, NR, IPHI, IMU, IBC, IANG, IUPDOWN, OLDNPTS
    integer joct, ipt, ip, iorder, ipcell, icorner
    logical LAMBERTIAN, BTEST, theflag
    real A, TRANSMIN
    save OLDNPTS

    ! Set the minimum transmission for the cell tracing (1 for singl
    TRANSMIN = 1.00

    ! Make the new grid cell/point sweeping order if need to
    if (NPTS .ne. OLDNPTS) then
        call SWEEPING_ORDER(NX, NY, NZ, NPTS, NCELLS,&
            & GRIDPTR, TREEPTR, CELLFLAGS,&
            & BCFLAG, IPFLAG, SWEEPORD, GRIDRAD)
    end if

    call FIND_BOUNDARY_POINTS(BCFLAG, IPFLAG, NPTS, SWEEPORD, GRIDPTR,&
        & GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'
    if (NPTS .ne. OLDNPTS) then
        ! Make the pointers to the grid points on the top and bottom bou
        call BOUNDARY_PNTS(NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
            & ZGRID(1), ZGRID(NZ), GRIDPOS,&
            & NTOPPTS, NBOTPTS, BCPTR)
        if (SFCTYPE(1:1) .eq. 'V') then
            ! If this is a variable surface then bilinearly interpolate
            !   the surface temperature and reflection parameters to the
            !   bottom grid points.
            call SURFACE_PARM_INTERP(NBOTPTS, BCPTR(1, 2), GRIDPOS,&
                & SRCTYPE, WAVENO, WAVELEN, UNITS,&
                & NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
                & SFCGRIDPARMS)
        end if
    end if
    ! Make the isotropic radiances for the top boundary
    call COMPUTE_TOP_RADIANCES(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        & UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))

    ! Zero the flux arrays and set temporary radiance to -1 for not
    do I = 1, NPTS
        FLUXES(1, I) = 0.0
        FLUXES(2, I) = 0.0
        GRIDRAD(1, I) = -1.0
    end do
    ! Zero the output spherical harmonic radiance
    NR = RSHPTR(NPTS + 1)
    do I = 1, NR
        RADIANCE(:, I) = 0.0
    end do

    ! print *, 'Path integration:'
    ! Loop over the zenith angles (downwelling angles must be first
    IANG = 1
    do IMU = 1, NMU

        if (IMU .eq. NMU/2 + 1) then
            ! For the first upwelling angle, make the bottom boundary radi
            !   for the Lambertian surfaces.  Compute the upwelling bottom
            !   radiances using the downwelling fluxes.
            if (SFCTYPE .eq. 'FL') then
                call FIXED_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1, 2),&
                    & DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
                    & WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
            else if (SFCTYPE .eq. 'VL') then
                call VARIABLE_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR(1, 2),&
                    & DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
                    & NSTOKES, BCRAD(1, 1 + NTOPPTS))
            end if
        end if
        if (NSTOKES .eq. 1) then
            call SH_TO_DO_UNPOL(NPTS, ML, MM, NLM,&
                & NMU, NPHI0MAX, NPHI0(IMU),&
                & IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
                & SHPTR, SOURCE, WORK)
        else
            call SH_TO_DO(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
                & NMU, NPHI0MAX, NPHI0(IMU),&
                & IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
                & SHPTR, SOURCE, WORK)
        end if

        ! Enforce the discrete ordinate I source function to be non-ne
        do I = 1, NPTS
            do IPHI = 1, NPHI0(IMU)
                WORK(1, IPHI, I) = max(0.0, WORK(1, IPHI, I))
            end do
        end do

        ! Loop over the azimuthal angles
        do IPHI = 1, NPHI0(IMU)

            if (MU(IMU) .lt. 0.0) then
                ! For downward ordinates, initialize the top boundary radi
                IUPDOWN = 1
                do IBC = 1, NTOPPTS
                    I = BCPTR(IBC, 1)
                    GRIDRAD(:, I) = BCRAD(:, IBC)
                end do
            else
                ! Initialize the bottom boundary radiances
                IUPDOWN = 2
                if (.not. LAMBERTIAN) then
                    ! If a Lambertian surface, use the radiances computed abov
                    !   otherwise, compute the radiance for this angle
                    !   by integrating over the stored downwelling radiances.
                    ! If the surface reflection is uniform use the precomputed
                    !   BRDF matrix to do the integration.
                    if (UNIFORM_SFC_BRDF) then
                        call UNIFORM_BRDF_SURFACE(NBOTPTS, BCPTR(1, 2),&
                            & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, IMU, IPHI,&
                            & SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                            & SFCTYPE, NSFCPAR, SFCGRIDPARMS, SFC_BRDF_DO,&
                            & NSTOKES, BCRAD(1, 1 + NTOPPTS))
                    else
                        call VARIABLE_BRDF_SURFACE(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                            & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                            & MU(IMU), PHI(IMU, IPHI),&
                            & SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                            & SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                            & NSTOKES, BCRAD(1, 1 + NTOPPTS))
                    end if
                end if
                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    GRIDRAD(:, I) = BCRAD(:, IBC + NTOPPTS)
                end do
            end if

            ! Do backward path integration from each grid point not done
            if (btest(IPFLAG, 1) .and. btest(IPFLAG, 0)) then
                call BACK_INT_GRID1D(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    & GRIDPOS, SWEEPORD, MU(IMU), PHI(IMU, IPHI), TRANSMIN,&
                    & EXTINCT, NSTOKES, WORK, IPHI, GRIDRAD)
            else if (btest(IPFLAG, 1) .and.&
                & .not. (btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) then
                call BACK_INT_GRID2D(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    & GRIDPOS, SWEEPORD, MU(IMU), PHI(IMU, IPHI), TRANSMIN,&
                    & EXTINCT, NSTOKES, WORK, IPHI, GRIDRAD)
            else
                if (NSTOKES .eq. 1) then
                    call BACK_INT_GRID3D_UNPOL(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                        & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                        & GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
                        & MU(IMU), PHI(IMU, IPHI), TRANSMIN, BCFLAG, IPFLAG,&
                        & EXTINCT, WORK, IPHI, GRIDRAD)
                else
                    call BACK_INT_GRID3D(NX, NY, NZ, NPHI0MAX, NPTS, NCELLS,&
                        & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                        & GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
                        & MU(IMU), PHI(IMU, IPHI), TRANSMIN, BCFLAG, IPFLAG,&
                        & EXTINCT, NSTOKES, WORK, IPHI, GRIDRAD)
                end if
            end if

            ! Transfer the radiance for the (IMU,IPHI) discrete ordinate
            !   GRIDRAD to the work array.
            do I = 1, NPTS
                if (GRIDRAD(1, I) .lt. -0.0001) then
                    write (6, *) 'PATH_INTEGRATION: Grid point has no value',&
                        & I, MU(IMU), PHI(IMU, IPHI), GRIDRAD(1, I), XGRID(1), YGRID(1),&
                        & GRIDPOS(1, I), GRIDPOS(2, I), GRIDPOS(3, I)
                    ! stop
                    GRIDRAD(1, I) = 0.0
                else
                    WORK(:, IPHI, I) = GRIDRAD(:, I)
                    GRIDRAD(1, I) = -1.0
                end if
            end do

            ! Do the sum for the hemispheric fluxes for this angle
            A = abs(MU(IMU))*WTDO(IMU, IPHI)
            do I = 1, NPTS
                FLUXES(IUPDOWN, I) = FLUXES(IUPDOWN, I) + A*WORK(1, IPHI, I)
            end do

            ! For a general BRDF surface save the downwelling radiance
            !   at the bottom boundary points.
            if (.not. LAMBERTIAN .and. IMU .le. NMU/2) then
                I1 = NTOPPTS + NBOTPTS*IANG
                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    BCRAD(:, IBC + I1) = WORK(:, IPHI, I)
                end do
            end if

            IANG = IANG + 1
        end do

        ! Transform the radiance to spherical harmonics for this zenit
        if (NSTOKES .eq. 1) then
            call DO_TO_SH_UNPOL(NPTS, ML, MM, NLM,&
                & NMU, NPHI0MAX, NPHI0(IMU),&
                & IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
                & RSHPTR, WORK, RADIANCE)
        else
            call DO_TO_SH(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
                & NMU, NPHI0MAX, NPHI0(IMU),&
                & IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
                & RSHPTR, WORK, RADIANCE)
        end if
    end do

    OLDNPTS = NPTS

    return
end subroutine

subroutine BOUNDARY_PNTS(NPTS, NANG, LAMBERTIAN, MAXNBC, MAXBCRAD,&
    & ZBOT, ZTOP, GRIDPOS,&
    & NTOPPTS, NBOTPTS, BCPTR)
    !   Returns the BCPTR pointer array for the NTOPPTS top boundary
    ! points (BCPTR(*,1)) and the NBOTPTS bottom boundary points (BCPTR(
    ! All the grid points are looked at to find the boundary ones.
    implicit none
    integer NPTS, NANG, MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS
    integer BCPTR(MAXNBC, 2)
    logical LAMBERTIAN
    real ZBOT, ZTOP, GRIDPOS(3, NPTS)
    integer I, IT, IB, NA

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
                write (6, *) 'BOUNDARY_PNTS: MAXNBC exceeded - top.'
                stop
            end if
            BCPTR(IT, 1) = I
        end if
        if (GRIDPOS(3, I) .le. ZBOT) then
            IB = IB + 1
            if (IB .gt. MAXNBC .or. IT + IB*NA .gt. MAXBCRAD) then
                write (6, *) 'BOUNDARY_PNTS: MAXNBC exceeded - bottom.'
                stop
            end if
            BCPTR(IB, 2) = I
        end if
    end do
    NTOPPTS = IT
    NBOTPTS = IB
    return
end subroutine

subroutine SURFACE_PARM_INTERP(NBOTPTS, BCPTR, GRIDPOS,&
    & SRCTYPE, WAVENO, WAVELEN, UNITS,&
    & NXSFC, NYSFC, DELXSFC, DELYSFC, NSFCPAR, SFCPARMS,&
    & SFCGRIDPARMS)
    !   Bilinearly interpolates the surface properties from the regular
    ! surface grid to the boundary points of the adaptive grid.
    ! The locations of the NBOTPTS grid points on the bottom boundary
    ! (pointed to by BCPTR) are in the GRIDPOS array.  The input regular
    ! surface info is specified by the number of grid lines (NXSFC,NYSFC
    ! the grid spacing (DELXSFC,DELYSFC), the number of parameters at
    ! each grid point (NSFCPAR), and the surface parameters (SFCPARMS).
    ! The output interpolated surface parameters at the NBOTPTS boundary
    ! points are in SFCGRIDPARMS.  The first parameter is input as
    ! temperature, but converted to Planck function for output.
    implicit none
    integer NBOTPTS, BCPTR(*), NXSFC, NYSFC, NSFCPAR
    real WAVENO(2), WAVELEN
    real DELXSFC, DELYSFC
    real GRIDPOS(3, *), SFCPARMS(NSFCPAR, NXSFC + 1, NYSFC + 1)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS)
    character SRCTYPE*1, UNITS*1
    integer I, IBC, IX, IY, J
    real RX, RY, U, V, TEMP, PLANCK

    ! Loop over all bottom points
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        RX = GRIDPOS(1, I)/DELXSFC
        RY = GRIDPOS(2, I)/DELYSFC
        IX = max(1, min(NXSFC, int(RX) + 1))
        IY = max(1, min(NYSFC, int(RY) + 1))
        U = max(0.0, min(1.0, RX - (IX - 1)))
        V = max(0.0, min(1.0, RY - (IY - 1)))
        do J = 1, NSFCPAR
            SFCGRIDPARMS(J, IBC) = (1 - U)*(1 - V)*SFCPARMS(J, IX, IY)&
                & + (1 - U)*V*SFCPARMS(J, IX, IY + 1)&
                & + U*(1 - V)*SFCPARMS(J, IX + 1, IY) + U*V*SFCPARMS(J, IX + 1, IY + 1)
        end do
        if (SRCTYPE .eq. 'S') then
            PLANCK = 0.0
        else
            TEMP = SFCGRIDPARMS(1, IBC)
            call PLANCK_FUNCTION(TEMP, UNITS, WAVENO, WAVELEN, PLANCK)
        end if
        SFCGRIDPARMS(1, IBC) = PLANCK
    end do
    return
end subroutine

subroutine CHECK_UNIFORM_SFC(NXSFC, NYSFC, NSFCPAR, SFCPARMS,&
    & UNIFORM_SFC_BRDF)
    !   Checks to see if the surface reflection parameters are uniform,
    ! and if so sets the UNIFORM_SFC_BRDF flag.  The temperature need
    ! not be uniform.
    implicit none
    integer NXSFC, NYSFC, NSFCPAR
    real SFCPARMS(NSFCPAR, NXSFC + 1, NYSFC + 1)
    logical UNIFORM_SFC_BRDF
    integer IX, IY, J

    UNIFORM_SFC_BRDF = .true.
    do IY = 1, NYSFC
        do IX = 1, NXSFC
            do J = 2, NSFCPAR
                if (SFCPARMS(J, IX, IY) .ne. SFCPARMS(J, 1, 1)) then
                    UNIFORM_SFC_BRDF = .false.
                    return
                end if
            end do
        end do
    end do
    return
end subroutine

subroutine COMPUTE_TOP_RADIANCES(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
    & UNITS, NTOPPTS, NSTOKES, BCRAD)
    !   Returns the beginning of the BCRAD array with the downwelling
    ! radiance for the NTOPPTS top boundary points.  Currently all
    ! points have the same unpolarized isotropic radiance.
    implicit none
    integer NTOPPTS, NSTOKES
    real SKYRAD, WAVENO(2), WAVELEN
    real BCRAD(NSTOKES, *)
    character SRCTYPE*1, UNITS*1
    integer IBC, K
    real SKYRAD2

    ! At top, boundary radiance is any isotropic sky radiation
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        SKYRAD2 = SKYRAD
    else
        call PLANCK_FUNCTION(SKYRAD, UNITS, WAVENO, WAVELEN, SKYRAD2)
    end if
    ! Loop over all points assigning the uniform radiance
    do IBC = 1, NTOPPTS
        BCRAD(1, IBC) = SKYRAD2
        BCRAD(2:, IBC) = 0.0
    end do
    return
end subroutine

subroutine FIXED_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR,&
    & DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
    & WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD)
    !   Returns the BCRAD array with the upwelling unpolarized radiance
    ! a Lambertian surface for the bottom points (pointed to by BCPTR).
    ! If there is reflection the downwelling hemispheric flux at the
    ! bottom points is gotten (from FLUXES). There may also be thermal
    ! emission (depending on temperature and emissivity) or reflected
    ! collimated solar radiation (DIRFLUX).
    implicit none
    integer NBOTPTS, BCPTR(*), NSTOKES
    real DIRFLUX(*), FLUXES(2, *)
    real GNDTEMP, GNDALBEDO, WAVENO(2), WAVELEN
    real BCRAD(NSTOKES, *)
    character SRCTYPE*1, UNITS*1
    integer I, IBC, K
    real ALB, GNDRAD

    ALB = GNDALBEDO/acos(-1.0)
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        call PLANCK_FUNCTION(GNDTEMP, UNITS, WAVENO, WAVELEN, GNDRAD)
        GNDRAD = GNDRAD*(1.0 - GNDALBEDO)
    end if

    ! Loop over all bottom points finding the upwelling radiance
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        if (SRCTYPE .eq. 'S') then
            BCRAD(1, IBC) = ALB*(DIRFLUX(I) + FLUXES(1, I))
        else if (SRCTYPE .eq. 'T') then
            BCRAD(1, IBC) = GNDRAD + ALB*FLUXES(1, I)
        else if (SRCTYPE .eq. 'B') then
            BCRAD(1, IBC) = ALB*(DIRFLUX(I) + FLUXES(1, I)) + GNDRAD
        end if
        BCRAD(2:, IBC) = 0.0
    end do
    return
end subroutine

subroutine VARIABLE_LAMBERTIAN_BOUNDARY(NBOTPTS, BCPTR,&
    & DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
    & NSTOKES, BCRAD)
    !   Returns the BCRAD array with the upwelling unpolarized radiance
    ! variable Lambertian surface for the bottom points (pointed to by B
    ! SFCGRIDPARMS contains the interpolated surface properties for the
    ! bottom boundary points; the first is Planck function and the secon
    ! albedo.  The reflection is computed from the downwelling hemispher
    ! flux (in FLUXES).  There may also be thermal emission (depending o
    ! temperature and emissivity) or reflected collimated solar radiatio
    ! (from DIRFLUX).
    implicit none
    integer NBOTPTS, BCPTR(*), NSFCPAR, NSTOKES
    real DIRFLUX(*), FLUXES(2, *)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS)
    real BCRAD(NSTOKES, *)
    character SRCTYPE*1
    integer I, IBC, K
    real OPI, ALB, GNDRAD

    OPI = 1.0/acos(-1.0)

    ! Loop over all bottom points
    do IBC = 1, NBOTPTS
        I = BCPTR(IBC)
        ALB = SFCGRIDPARMS(2, IBC)
        if (SRCTYPE .eq. 'S') then
            BCRAD(1, IBC) = OPI*ALB*(DIRFLUX(I) + FLUXES(1, I))
        else
            GNDRAD = SFCGRIDPARMS(1, IBC)*(1 - ALB)
            if (SRCTYPE .eq. 'T') then
                BCRAD(1, IBC) = GNDRAD + OPI*ALB*FLUXES(1, I)
            else if (SRCTYPE .eq. 'B') then
                BCRAD(1, IBC) = OPI*ALB*(DIRFLUX(I) + FLUXES(1, I)) + GNDRAD
            end if
        end if
        BCRAD(2:, IBC) = 0.0
    end do
    return
end subroutine

subroutine VARIABLE_BRDF_SURFACE(NBOTPTS, IBEG, IEND, BCPTR,&
    & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MU2, PHI2,&
    & SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
    & SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES, BCRAD)
    !   Computes the upwelling Stokes vector of the bottom boundary poin
    ! for one outgoing direction using the specified bidirectional
    ! reflectance distribution function.  The upwelling Stokes vector
    ! includes the reflection of the incident radiance, the thermal
    ! emission (emissivity*Planck function), and the reflected direct
    ! solar flux (if applicable).  The upwelling radiance vector is the
    ! integral over all incident directions of the BRDF times the
    ! downwelling radiance, so a discrete sum is done and the integral
    ! weights (WTDO) are included. The general BRDF function is called
    ! to compute the Stokes reflectance matrix for the particular BRDF
    ! type (SFCTYPE) with parameters (SFCGRIDPARMS) for the incident and
    ! outgoing directions.  The emissivity is computed implicitly from t
    ! integral of the BRDF.  The outgoing direction is specified with
    ! (MU2,PHI2), and the BRDF is computed for all incident directions
    ! (loop over JMU,JPHI).  The incident downwelling radiance vectors a
    ! input in BCRAD(*,*,2...) and the outgoing upwelling radiance vecto
    ! are output in BCRAD(*,*,1).
    implicit none
    integer NBOTPTS, IBEG, IEND, BCPTR(NBOTPTS)
    integer NMU, NPHI0MAX, NPHI0(NMU), NSFCPAR, NSTOKES
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real MU2, PHI2, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS), BCRAD(NSTOKES, NBOTPTS, *)
    character SRCTYPE*1, SFCTYPE*2
    integer JMU, JPHI, IBC, I, JANG, K, K1, K2
    real OPI, REFLECT(4, 4), W, SUM0, SUM1, RADSPEC(4)

    OPI = 1.0/acos(-1.0)

    do IBC = IBEG, IEND

        ! Initialize the upwelling boundary radiances to zero or to
        !   the reflection of direct unpolarized solar flux.
        BCRAD(:, IBC, 1) = 0.0
        if (SRCTYPE .ne. 'T') then
            I = BCPTR(IBC)
            call SURFACE_BRDF(SFCTYPE(2:2), SFCGRIDPARMS(2, IBC), WAVELEN,&
                & MU2, PHI2, SOLARMU, SOLARAZ, NSTOKES, REFLECT)
            BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                & + OPI*REFLECT(1:NSTOKES, 1)*DIRFLUX(I)
        end if

        ! Integrate over the incident discrete ordinate directions (JMU,
        JANG = 1
        do JMU = 1, NMU/2
            do JPHI = 1, NPHI0(JMU)
                call SURFACE_BRDF(SFCTYPE(2:2), SFCGRIDPARMS(2, IBC),&
                    & WAVELEN, MU2, PHI2, MU(JMU), PHI(JMU, JPHI),&
                    & NSTOKES, REFLECT)
                W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
                ! Add in the polarized reflection
                do K1 = 1, NSTOKES
                    BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                        & + W*REFLECT(1:NSTOKES, K1)*BCRAD(K1, IBC, JANG + 1)
                end do
                ! Add in the polarized thermal emission
                BCRAD(1, IBC, 1) = BCRAD(1, IBC, 1)&
                    & + W*(1 - REFLECT(1, 1))*SFCGRIDPARMS(1, IBC)
                BCRAD(2:, IBC, 1) = BCRAD(2:, IBC, 1)&
                    & - W*REFLECT(2:NSTOKES, 1)*SFCGRIDPARMS(1, IBC)
                JANG = JANG + 1
            end do
        end do

    end do
    return
end subroutine

subroutine MAKE_SFC_BRDF_DO_MATRIX(SFCTYPE, WAVELEN,&
    & NXSFC, NYSFC, NSFCPAR, SFCPARMS, NSTOKES,&
    & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    & SFC_BRDF_DO)
    !   Makes the precomputed BRDF surface reflection matrix for
    ! discrete ordinate incident and outgoing directions.  For each
    ! pair of discrete ordinates SURFACE_BRDF is called and the
    ! integration weights are included.
    implicit none
    integer NXSFC, NYSFC, NSFCPAR
    integer NSTOKES, NMU, NPHI0MAX, NPHI0(NMU)
    real SFCPARMS(NSFCPAR, NXSFC + 1, NYSFC + 1), WAVELEN
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    character SFCTYPE*2
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
                    call SURFACE_BRDF(SFCTYPE(2:2), SFCPARMS(2, 1, 1),&
                        & WAVELEN, MU(IMU), PHI(IMU, IPHI),&
                        & MU(JMU), PHI(JMU, JPHI), NSTOKES, REFLECT)
                    W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
                    SFC_BRDF_DO(:, IMU1, IPHI, :, JMU, JPHI)&
                        & = W*REFLECT(1:NSTOKES, 1:NSTOKES)
                end do
            end do
        end do
    end do
    return
end subroutine

subroutine UNIFORM_BRDF_SURFACE(NBOTPTS, BCPTR,&
    & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, IMU, IPHI,&
    & SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
    & SFCTYPE, NSFCPAR, SFCGRIDPARMS, SFC_BRDF_DO,&
    & NSTOKES, BCRAD)
    !   Computes the upwelling Stokes vector of the bottom boundary poin
    ! for one outgoing discrete ordinate using the precomputed bidirecti
    ! reflectance distribution matrix for the discrete ordinates.
    ! The upwelling Stokes vector includes the reflection of the inciden
    ! radiance, the thermal emission (emissivity*Planck function), and t
    ! reflected direct solar flux (if applicable).  The upwelling radian
    ! vector is the integral over all incident directions of the BRDF ti
    ! the downwelling radiance, so a discrete sum is done using the SFC_
    ! matrix, which includes the integral weights (WTDO).  The emissivit
    ! is computed implicitly from the integral of the BRDF.  The outgoin
    ! direction is specified with (IMU,IPHI), and the BRDF is computed f
    ! all incident directions (loop over JMU,JPHI).  The incident downwe
    ! radiance vectors are input in BCRAD(*,*,2...) and the outgoing
    ! upwelling radiance vectors are output in BCRAD(*,*,1).
    implicit none
    integer NBOTPTS, BCPTR(NBOTPTS)
    integer NMU, NPHI0MAX, NPHI0(NMU), IMU, IPHI, NSFCPAR, NSTOKES
    real WTDO(NMU, NPHI0MAX), MU(NMU), PHI(NMU, NPHI0MAX)
    real WAVELEN, SOLARMU, SOLARAZ, DIRFLUX(*)
    real SFCGRIDPARMS(NSFCPAR, NBOTPTS)
    real SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX, NSTOKES, NMU/2, NPHI0MAX)
    real BCRAD(NSTOKES, NBOTPTS, *)
    character SRCTYPE*1, SFCTYPE*2
    integer IMU1, JMU, JPHI, IBC, I, JANG, K, K1, K2
    real OPI, REFLECT(4, 4), W

    OPI = 1.0/acos(-1.0)
    IMU1 = IMU - NMU/2
    do IBC = 1, NBOTPTS

        ! Initialize the upwelling boundary radiances to zero or to
        !   the reflection of direct unpolarized solar flux.
        BCRAD(:, IBC, 1) = 0.0
        if (SRCTYPE .ne. 'T') then
            I = BCPTR(IBC)
            call SURFACE_BRDF(SFCTYPE(2:2), SFCGRIDPARMS(2, IBC), WAVELEN,&
                & MU(IMU), PHI(IMU, IPHI), SOLARMU, SOLARAZ,&
                & NSTOKES, REFLECT)
            BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1)&
                & + OPI*REFLECT(1:NSTOKES, 1)*DIRFLUX(I)
        end if

        ! Integrate over the incident discrete ordinate directions (JMU,
        JANG = 1
        do JMU = 1, NMU/2
            do JPHI = 1, NPHI0(JMU)
                ! Add in the polarized reflection
                do K1 = 1, NSTOKES
                    BCRAD(:, IBC, 1) = BCRAD(:, IBC, 1) + BCRAD(K1, IBC, JANG + 1)&
                        & *SFC_BRDF_DO(:, IMU1, IPHI, K1, JMU, JPHI)
                end do
                ! Add in the polarized thermal emission
                W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
                BCRAD(1, IBC, 1) = BCRAD(1, IBC, 1) + SFCGRIDPARMS(1, IBC)&
                    & *(W - SFC_BRDF_DO(1, IMU1, IPHI, 1, JMU, JPHI))
                BCRAD(2:, IBC, 1) = BCRAD(2:, IBC, 1) - SFCGRIDPARMS(1, IBC)&
                    & *SFC_BRDF_DO(2:, IMU1, IPHI, 1, JMU, JPHI)
                JANG = JANG + 1
            end do
        end do
    end do
    return
end subroutine

subroutine SH_TO_DO(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    & NMU, NPHI0MAX, NPHI0,&
    & IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
    & SHPTR, INDATA, OUTDATA)
    !   Transforms the input data from spherical harmonic space to discr
    ! ordinate space for one zenith angle (IMU).  Successive transforms
    ! done in zenith and then azimuth angle.  The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(2*MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small prime
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
            ! Figure the max Fourier azimuthal mode we can do for this Nph
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
                        & + CMU1(2, J, IMU)*INDATA(2, IS + J)&
                        & + CMU1(5, J, IMU)*INDATA(3, IS + J)
                    SUMUV(3, M) = SUMUV(3, M)&
                        & + CMU1(6, J, IMU)*INDATA(2, IS + J)&
                        & + CMU1(3, J, IMU)*INDATA(3, IS + J)
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
                    ! For sines and cosines use a real FFT, and rearrange term
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
                    call RFFTB(NPHI0, TMP, WSAVE(1, IMU))
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
end subroutine

subroutine SH_TO_DO_UNPOL(NPTS, ML, MM, NLM,&
    & NMU, NPHI0MAX, NPHI0,&
    & IMU, FFTFLAG, CMU1, CPHI1, WSAVE,&
    & SHPTR, INDATA, OUTDATA)
    !   Transforms the input data from spherical harmonic space to discr
    ! ordinate space for one zenith angle (IMU) for the unpolarized case
    ! Successive transforms are done in zenith and then azimuth angle.
    ! The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(2*MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small prime
    ! (powers of 2 are best).  The spherical harmonic truncation of the
    ! input array is adaptive (with pointer SHPTR); the last term for
    ! each grid point must be the maximum m for the particular l term.
    ! The number of floating point operations for the DFT is 2*Nmu*Nphi0
    ! The number of floating point operations for the zenith angle trans
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
            ! Figure the max Fourier azimuthal mode we can do for this Nph
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
                call RFFTB(NPHI0, TMP, WSAVE(1, IMU))
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
end subroutine

subroutine DO_TO_SH(NPTS, NSTOKES, NSTLEG, ML, MM, NLM,&
    & NMU, NPHI0MAX, NPHI0,&
    & IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    & RSHPTR, INDATA, OUTDATA)
    !   Transforms the input field from discrete ordinate space to spher
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See SH_TO_DO for more).
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
        ! Figure the max Fourier azimuthal mode we can do for this Nph
        ME = max(0, min((NPHI0/2) - 1, MOFJ(NSH)))
        MS = -ME
        ! First do Fourier transform from phi to m for each mu
        if (FFTFLAG(IMU)) then
            do N = 1, NSTOKES
                do K = 1, NPHI0
                    TMP(K) = INDATA(N, K, I)
                end do
                ! For sines and cosines use a real FFT, and rearrange term
                call RFFTF(NPHI0, TMP, WSAVE(1, IMU))
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
                    & + CMU2(5, IMU, J)*SUMUV(3, M)
                OUTDATA(3, IS + J) = OUTDATA(3, IS + J) + CMU2(6, IMU, J)*SUMUV(2, M)&
                    & + CMU2(3, IMU, J)*SUMUV(3, M)
            end if
            if (NSTOKES .eq. 4) then
                OUTDATA(4, IS + J) = OUTDATA(4, IS + J) + CMU2(4, IMU, J)*SUMUV(4, M)
            end if
        end do
        IS = IS + NSH
    end do
    deallocate (MOFJ, SUMUV, SUMCS, TMP)
    return
end subroutine

subroutine DO_TO_SH_UNPOL(NPTS, ML, MM, NLM,&
    & NMU, NPHI0MAX, NPHI0,&
    & IMU, FFTFLAG, CMU2, CPHI2, WSAVE,&
    & RSHPTR, INDATA, OUTDATA)
    !   Transforms the input field from discrete ordinate space to spher
    ! harmonic space for one zenith angle (IMU) by keeping a running
    ! sum in OUTDATA between successive calls.  (See SH_TO_DO for more).
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
            ! For sines and cosines use a real FFT, and rearrange term
            call RFFTF(NPHI0, TMP, WSAVE(1, IMU))
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
end subroutine

subroutine SWEEPING_ORDER(NX, NY, NZ, NPTS, NCELLS,&
    & GRIDPTR, TREEPTR, CELLFLAGS,&
    & BCFLAG, IPFLAG, SWEEPORD, GRIDRAD)
    !   Make the array of grid point sweeping order for each discrete
    ! ordinate octant (SWEEPORD).  Rather than pointing to the grid poin
    ! directly, the elements of SWEEPORD point to the cell (bits 3-30) a
    ! the corner of the cell (bits 0-2).  The first index is the grid
    ! point visiting order (1-NPTS), and the second index refers to the
    ! octant, but is not the octant number.  The octant index is transla
    ! to the octant number with IOCTORDER - this complication arises bec
    ! of the desire to use just 2 or 4 octants in 1D or 2D.
    implicit none
    integer NX, NY, NZ, NPTS, NCELLS, BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDRAD(NPTS)
    integer NOCT, IOCT, JOCT, NXC, NYC
    integer IPCELL, IPT, IORDER
    integer INDEXCORN, ICORNER, CORNDOG(8, 8), IOCTORDER(8)
    logical SWEEP_NEXT_CELL, BTEST
    data CORNDOG/8, 7, 6, 5, 4, 3, 2, 1, 7, 8, 5, 6, 3, 4, 1, 2, 6, 5, 8, 7, 2, 1, 4, 3,&
        & 5, 6, 7, 8, 1, 2, 3, 4, 4, 3, 2, 1, 8, 7, 6, 5, 3, 4, 1, 2, 7, 8, 5, 6,&
        & 2, 1, 4, 3, 6, 5, 8, 7, 1, 2, 3, 4, 5, 6, 7, 8/
    data IOCTORDER/1, 5, 2, 6, 3, 7, 4, 8/

    ! IOCT is the octant the ray come from or the discrete ordinate
    !   IOCT = 1 + BITX + 2*BITY + 4*BITZ

    if (btest(IPFLAG, 1) .and. btest(IPFLAG, 0)) then
        NOCT = 2
    else if (btest(IPFLAG, 1)) then
        NOCT = 4
    else
        NOCT = 8
    end if
    NXC = NX
    if (btest(BCFLAG, 0)) NXC = NX + 1
    if (btest(BCFLAG, 2) .and. .not. btest(IPFLAG, 0)) NXC = NX - 1
    NYC = NY
    if (btest(BCFLAG, 1)) NYC = NY + 1
    if (btest(BCFLAG, 3) .and. .not. btest(IPFLAG, 1)) NYC = NY - 1

    ! Loop over the number of octants needed; get the octant number
    do JOCT = 1, NOCT
        IOCT = IOCTORDER(JOCT)

        ! Set each grid point flag to indicate unvisited.
        do IPT = 1, NPTS
            GRIDRAD(IPT) = -1.0
        end do

        ! Sweep through all the grid points
        IORDER = 1
        IPCELL = 0
        INDEXCORN = 8
100     continue
        ! INDEXCORN is a logical ordering index to the 8 corners in a ce
        !   it is not the actual corner. ICORNER is the actual corner an
        !   depends on the ray direction (IOCT) and INDEXCORN.
        !   Scan the 8 points in a cell according to the ray direction:
        !   Do the ones furthest "upstream" first; loop until we find
        !   a point that is not already done.
        if (INDEXCORN .eq. 8) then
            ! If done with the current cell, then go to next
            INDEXCORN = 1
            if (.not. SWEEP_NEXT_CELL(BCFLAG, NXC, NYC, NZ, IOCT,&
                & IPCELL, TREEPTR, CELLFLAGS)) goto 900
        else
            INDEXCORN = INDEXCORN + 1
        end if
        ! Compute the corresponding point
        !   ICORNER = IAND(NOT(IEOR(INDEXCORN-1,IOCT-1)),7) + 1
        ICORNER = CORNDOG(INDEXCORN, IOCT)
        IPT = GRIDPTR(ICORNER, IPCELL)
        if (GRIDRAD(IPT) .ge. 0.0) goto 100

        GRIDRAD(IPT) = 1.0
        ! Fill in the next entry with the grid cell and corner of th
        SWEEPORD(IORDER, JOCT) = ior(ishft(IPCELL, 3), ICORNER - 1)
        IORDER = IORDER + 1
        goto 100
900     continue
    end do

    return
end subroutine

subroutine BACK_INT_GRID3D(NX, NY, NZ, NA, NPTS, NCELLS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
    & MU, PHI, TRANSMIN, BCFLAG, IPFLAG,&
    & EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)
    !   Sweeps through the spatial grid computing radiance vectors for t
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  The sweeping
    ! order through the grid points is in SWEEPORD. If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is tr
    ! back to a face containing known (valid) radiances.  As each cell i
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! at the end of the ray path is interpolated from the known grid poi
    ! radiances on the face.  The resulting radiances are put in GRIDRAD
    implicit none
    integer NX, NY, NZ, NA, NPTS, NCELLS, NSTOKES, KANG
    integer BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS), XGRID(*), YGRID(*), ZGRID(NZ)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITY, BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer I1, I2, I3, I4, IOPP, JFACE, KFACE, IC, IZ
    integer GRIDFACE(4, 6), OPPFACE(6)
    integer IORDER, ICORNER, JOCTORDER3(8), JOCTORDER2(8)
    logical VALIDRAD, VALIDFACE, IPINX, IPINY, BTEST
    double precision PI, CX, CY, CZ, CXINV, CYINV, CZINV
    double precision XE, YE, ZE
    double precision SO, SOX, SOY, SOZ, EPS, path
    double precision U, V, F1, F2, F3, F4
    double precision EXT, EXT0, EXT1
    double precision SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    double precision EXT0P, SRCEXT0P(NSTOKES)
    double precision TAU, TRANSCELL, ABSCELL, TRANSMIT
    double precision RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        & 1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/
    data JOCTORDER3/1, 3, 5, 7, 2, 4, 6, 8/, JOCTORDER2/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    PI = acos(-1.0)
    CX = sqrt(1.0 - MU**2)*cos(PHI + PI)
    CY = sqrt(1.0 - MU**2)*sin(PHI + PI)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-5) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negativ
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [BACK_INT_GRID]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    if (btest(IPFLAG, 1)) then
        JOCT = JOCTORDER2(IOCT)
    else
        JOCT = JOCTORDER3(IOCT)
    end if

    IZ = (NZ - 1)*(1 - BITZ) + 1
    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)

        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if ((GRIDPOS(3, IPT) .lt. ZGRID(IZ) .and. BITZ .eq. 0) .or.&
                &(GRIDPOS(3, IPT) .gt. ZGRID(IZ) .and. BITZ .eq. 1)) then
                IZ = IZ + 2*BITZ - 1
                call CALC_BOUNDARY_RADIANCES(BCFLAG, IPFLAG, JOCT, IZ,&
                    & NX, NY, NZ, XGRID, YGRID, ZGRID,&
                    & NA, NPTS, NCELLS, GRIDPTR,&
                    & NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,&
                    & MU, PHI, EXTINCT, NSTOKES, SOURCE,&
                    & KANG, GRIDRAD)
            end if
        end if

        if (GRIDRAD(1, IPT) .ge. 0.0) goto 500

        ! Start off at the grid point
        ICELL = IPCELL
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        RAD(:) = 0.0
        SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
        XE = GRIDPOS(1, IPT)
        YE = GRIDPOS(2, IPT)
        ZE = GRIDPOS(3, IPT)

        ! Loop until finding a face with known radiances
        VALIDRAD = .false.
        do while (.not. VALIDRAD)
            ! Make sure current cell is valid
            if (ICELL .le. 0) then
                write (6, *) 'BACK_INT_GRID: ICELL=0  ',&
                    & MU, PHI, IPT, INEXTCELL, xe, ye, ze,&
                    & gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                stop
            end if
            IPINX = btest(int(CELLFLAGS(ICELL)), 0)
            IPINY = btest(int(CELLFLAGS(ICELL)), 1)

            ! Find boundaries of the current cell
            ! Find the three possible intersection planes (X,Y,Z)
            !   from the coordinates of the opposite corner grid point
            IOPP = GRIDPTR(9 - IOCT, ICELL)
            ! Get the distances to the 3 planes and select the closest
            if (IPINX) then
                SOX = 1.0e20
            else
                SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
            end if
            if (IPINY) then
                SOY = 1.0e20
            else
                SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
            end if
            SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
            SO = min(SOX, SOY, SOZ)
            if (SO .lt. -EPS) then
                write (6, *) 'BACK_INT_GRID: SO<0  ',&
                    & MU, PHI, XE, YE, ZE, SOX, SOY, SOZ, SO, IPT, ICELL
                stop
            end if
            ! Compute the coordinates of the cell exitting location
            XE = XE + SO*CX
            YE = YE + SO*CY
            ZE = ZE + SO*CZ

            ! Get the intersection face number (i.e. neighptr index)
            if (SOX .le. SOZ .and. SOX .le. SOY) then
                IFACE = 2 - BITX
                JFACE = 1
            else if (SOY .le. SOZ) then
                IFACE = 4 - BITY
                JFACE = 2
            else
                IFACE = 6 - BITZ
                JFACE = 3
            end if
            ! Get the next cell to go to
            INEXTCELL = NEIGHPTR(IFACE, ICELL)
            if (INEXTCELL .lt. 0) then
                call NEXT_CELL(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
            end if

            ! Get the face grid pointers
            if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
                ! If going to same or larger face then use previous face
                KFACE = IFACE
                IC = ICELL
            else
                ! If going to smaller face then use next face (more accura
                KFACE = OPPFACE(IFACE)
                IC = INEXTCELL
            end if
            I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
            I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
            I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
            I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
            ! Compute the face interpolation factors
            if (JFACE .eq. 1) then
                U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                if (IPINY) then
                    V = 0.5
                else
                    V = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I2) - GRIDPOS(2, I1))
                end if
            else if (JFACE .eq. 2) then
                U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                if (IPINX) then
                    V = 0.5
                else
                    V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                end if
            else
                if (IPINY) then
                    U = 0.5
                else
                    U = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I3) - GRIDPOS(2, I1))
                end if
                if (IPINX) then
                    V = 0.5
                else
                    V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                end if
            end if
            if (u .lt. -0.0001 .or. v .lt. -0.0001 .or.&
                & u .gt. 1.0001 .or. v .gt. 1.0001) then
                print *, 'BACK_INT_GRID3D: u,v<0 or u,v>1: ',&
                    & mu, phi, xe, ye, ze, u, v, iface, jface, kface,&
                    & gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
            end if

            ! Get the location coordinate (does the boundary wrapping)
            if (INEXTCELL .gt. 0) then
                if (JFACE .eq. 1) then
                    XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
                else if (JFACE .eq. 2) then
                    YE = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
                else
                    ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                end if
            end if

            ! Interpolate extinction and source function at face interse
            F1 = (1 - U)*(1 - V)
            F2 = (1 - U)*V
            F3 = U*(1 - V)
            F4 = U*V
            EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)&
                & + F3*EXTINCT(I3) + F4*EXTINCT(I4)
            ! Correctly interpolate source using extinction*source
            SRCEXT0(:) = F1*SOURCE(:, KANG, I1)*EXTINCT(I1)&
                & + F2*SOURCE(:, KANG, I2)*EXTINCT(I2)&
                & + F3*SOURCE(:, KANG, I3)*EXTINCT(I3)&
                & + F4*SOURCE(:, KANG, I4)*EXTINCT(I4)
            ! Compute the cell radiance: integration of the source funct
            EXT = 0.5*(EXT0 + EXT1)
            TAU = EXT*SO
            if (TAU .ge. 0.5) then
                TRANSCELL = exp(-TAU)
                ABSCELL = 1.0 - TRANSCELL
            else
                ABSCELL = TAU*(1.0 - 0.5*TAU*&
                    & (1.0 - 0.33333333333*TAU*(1 - 0.25*TAU)))
                TRANSCELL = 1.0 - ABSCELL
            end if
            if (TAU .le. 2.0) then
                if (EXT .eq. 0.0) then
                    SRC(:) = 0.0
                else
                    ! Linear extinction, linear source*extinction, to first
                    SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + 0.08333333333&
                        & *(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                end if
            else
                ! Combined first order expansion and constant extinction f
                EXT0P = EXT0
                SRCEXT0P(:) = SRCEXT0(:)
                if (TAU .gt. 4.0) then
                    EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                    if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                end if
                SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                    & + (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                    & *(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
            end if
            SRC(1) = max(SRC(1), 0.0d0)

            ! Add in the cell radiance and update the transmission.
            ! If this is a boundary or the transmission is below the
            !   cutoff and we have a valid radiance then set the flag
            !   to stop the tracing and add in the interpolated face radia
            RAD(:) = RAD(:) + TRANSMIT*SRC(:)*ABSCELL
            TRANSMIT = TRANSMIT*TRANSCELL
            if (RAD(1) < -1.0e-5) then
                print '(A,3(1X,F6.4),8(1X,E12.5))',&
                    & 'BACK_INT_GRID3D: RAD<0 ',&
                    & xe, ye, ze, ext0, ext1, so, tau, transmit, src, abscell, rad(:)
                call ABORT_SHDOM_MPI('The end')
            end if
            VALIDFACE = (GRIDRAD(1, I1) .ge. -0.1 .and. GRIDRAD(1, I2) .ge. -0.1&
                & .and. GRIDRAD(1, I3) .ge. -0.1 .and. GRIDRAD(1, I4) .ge. -0.1)
            if (INEXTCELL .le. 0 .or.&
                & (TRANSMIT .le. TRANSMIN .and. VALIDFACE)) then
                if (VALIDFACE) then
                    VALIDRAD = .true.
                    RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)&
                        & + F3*GRIDRAD(:, I3) + F4*GRIDRAD(:, I4)
                    RAD(:) = RAD(:) + TRANSMIT*RAD0(:)
                else
                    print *, 'BACK_INT_GRID3D: INEXTCELL=0: ', xe, ye, ze, icell
                    print *, iz, joct, mu, phi
                    print '(i7,4(1x,f7.4))', i1, GRIDRAD(1, I1),&
                        & gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                    print '(i7,4(1x,f7.4))', i2, GRIDRAD(1, I2),&
                        & gridpos(1, i2), gridpos(2, i2), gridpos(3, i2)
                    print '(i7,4(1x,f7.4))', i3, GRIDRAD(1, I3),&
                        & gridpos(1, i3), gridpos(2, i3), gridpos(3, i3)
                    print '(i7,4(1x,f7.4))', i4, GRIDRAD(1, I4),&
                        & gridpos(1, i4), gridpos(2, i4), gridpos(3, i4)
                    print '(i7,8x,3(1x,f7.4))', ipt,&
                        & gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                    stop
                end if
            else
                EXT1 = EXT0
                SRCEXT1(:) = SRCEXT0(:)
                ICELL = INEXTCELL
            end if
        end do
        GRIDRAD(:, IPT) = RAD(:)

500     continue
    end do

    return
end subroutine

subroutine BACK_INT_GRID3D_UNPOL(NX, NY, NZ, NA, NPTS, NCELLS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & GRIDPOS, XGRID, YGRID, ZGRID, SWEEPORD,&
    & MU, PHI, TRANSMIN, BCFLAG, IPFLAG,&
    & EXTINCT, SOURCE, KANG, GRIDRAD)
    !   Sweeps through the spatial grid computing radiances for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  The sweeping
    ! order through the grid points is in SWEEPORD. If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is tr
    ! back to a face containing known (valid) radiances.  As each cell i
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! at the end of the ray path is interpolated from the known grid poi
    ! radiances on the face.  The resulting radiances are put in GRIDRAD
    integer NX, NY, NZ, NA, NPTS, NCELLS, KANG, BCFLAG, IPFLAG
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS), XGRID(*), YGRID(*), ZGRID(NZ)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NA, NPTS), GRIDRAD(NPTS)
    integer BITX, BITY, BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer I1, I2, I3, I4, IOPP, JFACE, KFACE, IC, IZ
    integer GRIDFACE(4, 6), OPPFACE(6)
    integer IORDER, ICORNER, JOCTORDER3(8), JOCTORDER2(8)
    logical VALIDRAD, VALIDFACE, IPINX, IPINY, BTEST
    double precision PI, CX, CY, CZ, CXINV, CYINV, CZINV
    double precision XE, YE, ZE
    double precision SO, SOX, SOY, SOZ, EPS, path
    double precision U, V, F1, F2, F3, F4
    double precision EXT, EXT0, EXT1, SRC, SRCEXT0, SRCEXT1
    double precision EXT0P, SRCEXT0P
    double precision TAU, TRANSCELL, ABSCELL, RAD, RAD0, TRANSMIT
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        & 1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/
    data JOCTORDER3/1, 3, 5, 7, 2, 4, 6, 8/, JOCTORDER2/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    PI = acos(-1.0)
    CX = sqrt(1.0 - MU**2)*cos(PHI + PI)
    CY = sqrt(1.0 - MU**2)*sin(PHI + PI)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-5) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negativ
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [BACK_INT_GRID]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    if (btest(IPFLAG, 1)) then
        JOCT = JOCTORDER2(IOCT)
    else
        JOCT = JOCTORDER3(IOCT)
    end if

    IZ = (NZ - 1)*(1 - BITZ) + 1
    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)

        if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
            ! For multiple processors get the subdomain boundary radiances
            !   for each Z slab
            if ((GRIDPOS(3, IPT) .lt. ZGRID(IZ) .and. BITZ .eq. 0) .or.&
                & (GRIDPOS(3, IPT) .gt. ZGRID(IZ) .and. BITZ .eq. 1)) then
                IZ = IZ + 2*BITZ - 1
                call CALC_BOUNDARY_RADIANCES(BCFLAG, IPFLAG, JOCT, IZ,&
                    & NX, NY, NZ, XGRID, YGRID, ZGRID,&
                    & NA, NPTS, NCELLS, GRIDPTR,&
                    & NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,&
                    & MU, PHI, EXTINCT, 1, SOURCE,&
                    & KANG, GRIDRAD)
            end if
        end if

        if (GRIDRAD(IPT) .ge. 0.0) goto 500

        ! Start off at the grid point
        ICELL = IPCELL
        RAD = 0.0
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        SRCEXT1 = EXT1*SOURCE(KANG, IPT)
        XE = GRIDPOS(1, IPT)
        YE = GRIDPOS(2, IPT)
        ZE = GRIDPOS(3, IPT)

        ! Loop until finding a face with known radiances
        VALIDRAD = .false.
        do while (.not. VALIDRAD)
            ! Make sure current cell is valid
            if (ICELL .le. 0) then
                write (6, *) 'BACK_INT_GRID: ICELL=0  ',&
                    & MU, PHI, IPT, INEXTCELL, xe, ye, ze,&
                    & gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                stop
            end if
            IPINX = btest(int(CELLFLAGS(ICELL)), 0)
            IPINY = btest(int(CELLFLAGS(ICELL)), 1)

            ! Find boundaries of the current cell
            ! Find the three possible intersection planes (X,Y,Z)
            !   from the coordinates of the opposite corner grid point
            IOPP = GRIDPTR(9 - IOCT, ICELL)
            ! Get the distances to the 3 planes and select the closest
            if (IPINX) then
                SOX = 1.0e20
            else
                SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
            end if
            if (IPINY) then
                SOY = 1.0e20
            else
                SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
            end if
            SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
            SO = min(SOX, SOY, SOZ)
            if (SO .lt. -EPS) then
                write (6, *) 'BACK_INT_GRID: SO<0  ',&
                    & MU, PHI, XE, YE, ZE, SOX, SOY, SOZ, SO, IPT, ICELL
                stop
            end if
            ! Compute the coordinates of the cell exitting location
            XE = XE + SO*CX
            YE = YE + SO*CY
            ZE = ZE + SO*CZ

            ! Get the intersection face number (i.e. neighptr index)
            if (SOX .le. SOZ .and. SOX .le. SOY) then
                IFACE = 2 - BITX
                JFACE = 1
            else if (SOY .le. SOZ) then
                IFACE = 4 - BITY
                JFACE = 2
            else
                IFACE = 6 - BITZ
                JFACE = 3
            end if
            ! Get the next cell to go to
            INEXTCELL = NEIGHPTR(IFACE, ICELL)
            if (INEXTCELL .lt. 0) then
                call NEXT_CELL(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
            end if

            ! Get the face grid pointers
            if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
                ! If going to same or larger face then use previous face
                KFACE = IFACE
                IC = ICELL
            else
                ! If going to smaller face then use next face (more accura
                KFACE = OPPFACE(IFACE)
                IC = INEXTCELL
            end if
            I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
            I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
            I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
            I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
            ! Compute the face interpolation factors
            if (JFACE .eq. 1) then
                U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                if (IPINY) then
                    V = 0.5
                else
                    V = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I2) - GRIDPOS(2, I1))
                end if
            else if (JFACE .eq. 2) then
                U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
                if (IPINX) then
                    V = 0.5
                else
                    V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                end if
            else
                if (IPINY) then
                    U = 0.5
                else
                    U = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I3) - GRIDPOS(2, I1))
                end if
                if (IPINX) then
                    V = 0.5
                else
                    V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                end if
            end if
            if (u .lt. -0.0001 .or. v .lt. -0.0001 .or.&
                & u .gt. 1.0001 .or. v .gt. 1.0001) then
                print *, 'BACK_INT_GRID3D: u,v<0 or u,v>1: ',&
                    & mu, phi, xe, ye, ze, u, v, iface, jface, kface,&
                    & gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
            end if

            ! Get the location coordinate (does the boundary wrapping)
            if (INEXTCELL .gt. 0) then
                if (JFACE .eq. 1) then
                    XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
                else if (JFACE .eq. 2) then
                    YE = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
                else
                    ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                end if
            end if

            ! Interpolate extinction and source function at face interse
            F1 = (1 - U)*(1 - V)
            F2 = (1 - U)*V
            F3 = U*(1 - V)
            F4 = U*V
            EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)&
                & + F3*EXTINCT(I3) + F4*EXTINCT(I4)
            ! Correctly interpolate source using extinction*source
            SRCEXT0 = (F1*SOURCE(KANG, I1)*EXTINCT(I1)&
                & + F2*SOURCE(KANG, I2)*EXTINCT(I2)&
                & + F3*SOURCE(KANG, I3)*EXTINCT(I3)&
                & + F4*SOURCE(KANG, I4)*EXTINCT(I4))
            ! Compute the cell radiance: integration of the source funct
            EXT = 0.5*(EXT0 + EXT1)
            TAU = EXT*SO
            if (TAU .ge. 0.5) then
                TRANSCELL = exp(-TAU)
                ABSCELL = 1.0 - TRANSCELL
            else
                ABSCELL = TAU*(1.0 - 0.5*TAU*&
                    & (1.0 - 0.33333333333*TAU*(1 - 0.25*TAU)))
                TRANSCELL = 1.0 - ABSCELL
            end if
            if (TAU .le. 2.0) then
                if (EXT .eq. 0.0) then
                    SRC = 0.0
                else
                    ! Linear extinction, linear source*extinction, to first
                    SRC = (0.5*(SRCEXT0 + SRCEXT1)&
                        & + 0.08333333333*(EXT0*SRCEXT1 - EXT1*SRCEXT0)*SO)/EXT
                end if
            else
                ! Combined first order expansion and constant extinction f
                ! SRC = 0.5/EXT *( SRCEXT0+SRCEXT1&
                !     & + (EXT0*SRCEXT1-EXT1*SRCEXT0)*SO&
                !     & *(1-2/TAU+2*TRANSCELL/ABSCELL)/TAU )
                EXT0P = EXT0
                SRCEXT0P = SRCEXT0
                if (TAU .gt. 4.0) then
                    EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                    if (EXT0 .gt. 0.0) SRCEXT0P = SRCEXT0*EXT0P/EXT0
                end if
                SRC = 1.0/(EXT0P + EXT1)*(SRCEXT0P + SRCEXT1&
                    & + (EXT0P*SRCEXT1 - EXT1*SRCEXT0P)*2.0/(EXT0P + EXT1)&
                    & *(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
            end if
            SRC = max(SRC, 0.0d0)

            ! Add in the cell radiance and update the transmission.
            ! If this is a boundary or the transmission is below the
            ! cutoff and we have a valid radiance then set the flag
            ! to stop the tracing and add in the interpolated face radia
            RAD = RAD + TRANSMIT*SRC*ABSCELL
            TRANSMIT = TRANSMIT*TRANSCELL
            if (RAD < -1.0e-5) then
                print '(A,3(1X,F6.4),8(1X,E12.5))',&
                    & 'BACK_INT_GRID3D: RAD<0 ',&
                    & xe, ye, ze, ext0, ext1, so, tau, transmit, src, abscell, rad
                call ABORT_SHDOM_MPI('The end')
            end if
            VALIDFACE = (GRIDRAD(I1) .ge. -0.1 .and. GRIDRAD(I2) .ge. -0.1&
                & .and. GRIDRAD(I3) .ge. -0.1 .and. GRIDRAD(I4) .ge. -0.1)
            if (INEXTCELL .le. 0 .or.&
                & (TRANSMIT .le. TRANSMIN .and. VALIDFACE)) then
                if (VALIDFACE) then
                    VALIDRAD = .true.
                    RAD0 = F1*GRIDRAD(I1) + F2*GRIDRAD(I2)&
                        & + F3*GRIDRAD(I3) + F4*GRIDRAD(I4)
                    RAD = RAD + TRANSMIT*RAD0
                else
                    print *, 'BACK_INT_GRID3D: INEXTCELL=0: ', xe, ye, ze, icell
                    print *, iz, joct, mu, phi
                    print '(i7,4(1x,f7.4))', i1, GRIDRAD(I1),&
                        & gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                    print '(i7,4(1x,f7.4))', i2, GRIDRAD(I2),&
                        & gridpos(1, i2), gridpos(2, i2), gridpos(3, i2)
                    print '(i7,4(1x,f7.4))', i3, GRIDRAD(I3),&
                        & gridpos(1, i3), gridpos(2, i3), gridpos(3, i3)
                    print '(i7,4(1x,f7.4))', i4, GRIDRAD(I4),&
                        & gridpos(1, i4), gridpos(2, i4), gridpos(3, i4)
                    print '(i7,8x,3(1x,f7.4))', ipt,&
                        & gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                    stop
                end if
            else
                EXT1 = EXT0
                SRCEXT1 = SRCEXT0
                ICELL = INEXTCELL
            end if
        end do

        GRIDRAD(IPT) = RAD

500     continue
    end do

    return
end subroutine

subroutine BACK_INT_GRID2D(NX, NY, NZ, NA, NPTS, NCELLS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & GRIDPOS, SWEEPORD, MU, PHI, TRANSMIN,&
    & EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)
    !   Sweeps through the spatial grid computing radiances for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is tr
    ! back to a face containing known (valid) radiances.  As each cell i
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! ad the end of the ray path is interpolated from the known grid poi
    ! radiances on the face.  The resulting radiances are put in GRIDRAD
    implicit none
    integer NX, NY, NZ, NA, NPTS, NCELLS, KANG, NSTOKES
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE
    integer I1, I2, IOPP, JFACE, KFACE, IC, IZ
    integer GRIDFACE(2, 6), OPPFACE(6)
    integer IORDER, ICORNER, JOCTORDER(8)
    logical VALIDRAD, VALIDFACE, IPINX, BTEST
    double precision PI, CX, CZ, CXINV, CZINV
    double precision XE, YE, ZE
    double precision SO, SOX, SOZ, EPS
    double precision U, F1, F2
    double precision EXT, EXT0, EXT1
    double precision SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    double precision EXT0P, SRCEXT0P(NSTOKES)
    double precision TAU, TRANSCELL, ABSCELL, TRANSMIT
    double precision RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/1, 5, 2, 6, 0, 0, 0, 0, 1, 2, 5, 6/
    data OPPFACE/2, 1, 4, 3, 6, 5/
    data JOCTORDER/1, 3, 1, 3, 2, 4, 2, 4/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    PI = acos(-1.0)
    CX = sqrt(1.0 - MU**2)*cos(PHI + PI)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negativ
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [BACK_INT_GRID]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate
    IOCT = 1 + BITX + 4*BITZ
    JOCT = JOCTORDER(IOCT)

    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)
        if (GRIDRAD(1, IPT) .ge. 0.0) goto 500

        ! Start off at the grid point
        ICELL = IPCELL
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        RAD(:) = 0.0
        SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
        XE = GRIDPOS(1, IPT)
        YE = GRIDPOS(2, IPT)
        ZE = GRIDPOS(3, IPT)

        ! Loop until finding a face with known radiances
        VALIDRAD = .false.
        do while (.not. VALIDRAD)
            ! Make sure current cell is valid
            if (ICELL .le. 0) then
                write (6, *) 'BACK_INT_GRID: ICELL=0  ',&
                    & MU, PHI, IPT, INEXTCELL
                stop
            end if
            IPINX = btest(int(CELLFLAGS(ICELL)), 0)

            ! Find boundaries of the current cell
            ! Find the three possible intersection planes (X,Y,Z)
            !   from the coordinates of the opposite corner grid point
            IOPP = GRIDPTR(9 - IOCT, ICELL)
            ! Get the distances to the 3 planes and select the closest
            if (IPINX) then
                SOX = 1.0e20
            else
                SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
            end if
            SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
            SO = min(SOX, SOZ)
            if (SO .lt. -EPS) then
                write (6, *) 'BACK_INT_GRID2D: SO<0  ',&
                    & MU, PHI, XE, ZE, SOX, SOZ, SO, IPT, ICELL
                stop
            end if
            ! Compute the coordinates of the cell exitting location
            XE = XE + SO*CX
            ZE = ZE + SO*CZ

            ! Get the intersection face number (i.e. neighptr index)
            if (SOX .le. SOZ) then
                IFACE = 2 - BITX
                JFACE = 1
            else
                IFACE = 6 - BITZ
                JFACE = 3
            end if
            ! Get the next cell to go to
            INEXTCELL = NEIGHPTR(IFACE, ICELL)
            if (INEXTCELL .lt. 0) then
                call NEXT_CELL(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
            end if

            ! Get the face grid pointers
            if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
                ! If going to same or larger face then use previous face
                KFACE = IFACE
                IC = ICELL
            else
                ! If going to smaller face then use next face (more accura
                KFACE = OPPFACE(IFACE)
                IC = INEXTCELL
            end if
            I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
            I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
            ! Compute the face interpolation factors
            if (JFACE .eq. 1) then
                U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I2) - GRIDPOS(3, I1))
            else
                if (IPINX) then
                    U = 0.5
                else
                    U = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
                end if
            end if

            ! Get the location coordinate (does the boundary wrapping)
            if (INEXTCELL .gt. 0) then
                if (JFACE .eq. 1) then
                    XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
                else
                    ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
                end if
            end if

            ! Interpolate extinction and source function at face interse
            F1 = 1 - U
            F2 = U
            EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2)
            ! Correctly interpolate source using extinction*source
            SRCEXT0(:) = (F1*SOURCE(:, KANG, I1)*EXTINCT(I1)&
                & + F2*SOURCE(:, KANG, I2)*EXTINCT(I2))
            ! Compute the cell radiance: integration of the source funct
            EXT = 0.5*(EXT0 + EXT1)
            TAU = EXT*SO
            if (TAU .ge. 0.5) then
                TRANSCELL = exp(-TAU)
                ABSCELL = 1.0 - TRANSCELL
            else
                ABSCELL = TAU*(1.0 - 0.5*TAU*&
                    & (1.0 - 0.33333333333*TAU*(1 - 0.25*TAU)))
                TRANSCELL = 1.0 - ABSCELL
            end if
            if (TAU .le. 2.0) then
                if (EXT .eq. 0.0) then
                    SRC(:) = 0.0
                else
                    ! Linear extinction, linear source*extinction, to first
                    SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + 0.08333333333&
                        & *(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                end if
            else
                ! Combined first order expansion and constant extinction f
                EXT0P = EXT0
                SRCEXT0P(:) = SRCEXT0(:)
                if (TAU .gt. 4.0) then
                    EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                    if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                end if
                SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                    & + (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                    & *(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
            end if
            SRC(1) = max(SRC(1), 0.0d0)

            ! Add in the cell radiance and update the transmission.
            ! If this is a boundary or the transmission is below the
            ! cutoff and we have a valid radiance then set the flag
            ! to stop the tracing and add in the interpolated face radia
            RAD(:) = RAD(:) + TRANSMIT*SRC(:)*ABSCELL
            TRANSMIT = TRANSMIT*TRANSCELL
            VALIDFACE = (GRIDRAD(1, I1) .ge. -0.1 .and. GRIDRAD(1, I2) .ge. -0.1)
            if (INEXTCELL .le. 0 .or.&
                & (TRANSMIT .le. TRANSMIN .and. VALIDFACE)) then
                if (VALIDFACE) then
                    VALIDRAD = .true.
                    RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)
                    RAD(:) = RAD(:) + TRANSMIT*RAD0(:)
                else
                    print *, 'BACK_INT_GRID2D: INEXTCELL=0: ', xe, ye, ze, icell
                    print *, iz, joct, mu, phi
                    print '(i7,4(1x,f7.4))', i1, GRIDRAD(1, I1),&
                        & gridpos(1, i1), gridpos(2, i1), gridpos(3, i1)
                    print '(i7,4(1x,f7.4))', i2, GRIDRAD(1, I2),&
                        & gridpos(1, i2), gridpos(2, i2), gridpos(3, i2)
                    print '(i7,8x,3(1x,f7.4))', ipt,&
                        & gridpos(1, ipt), gridpos(2, ipt), gridpos(3, ipt)
                    stop
                end if
            else
                EXT1 = EXT0
                SRCEXT1(:) = SRCEXT0(:)
                ICELL = INEXTCELL
            end if
        end do

        GRIDRAD(:, IPT) = RAD(:)
500     continue
    end do
    return
end subroutine

subroutine BACK_INT_GRID1D(NX, NY, NZ, NA, NPTS, NCELLS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & GRIDPOS, SWEEPORD, MU, PHI, TRANSMIN,&
    & EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD)
    !   Sweeps through the spatial grid computing radiances for the
    ! discrete ordinate direction (MU,PHI) by integrating the source
    ! function and extinction backward from each point.  If a point is
    ! already done (nonnegative in GRIDRAD) then it is skipped.
    ! For each grid point a ray opposite to the ordinate direction is tr
    ! back to a face containing known (valid) radiances.  As each cell i
    ! crossed (usually only one) the integration of the source function
    ! across the cell is done and added to the radiance.  The radiance
    ! ad the end of the ray path is interpolated from the known grid poi
    ! radiances on the face.  The resulting radiances are put in GRIDRAD
    implicit none
    integer NX, NY, NZ, NA, NPTS, NCELLS, KANG, NSTOKES
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SWEEPORD(NPTS, *)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS)
    real MU, PHI, TRANSMIN
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)
    integer BITZ, IOCT, JOCT
    integer IPCELL, ICELL, INEXTCELL, IPT, IFACE, I1, IOPP
    integer GRIDFACE(6)
    integer IORDER, ICORNER, JOCTORDER(8)
    logical VALIDRAD
    double precision PI, CZ, CZINV
    double precision XE, YE, ZE
    double precision SO, EPS
    double precision EXT, EXT0, EXT1
    double precision SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    double precision EXT0P, SRCEXT0P(NSTOKES)
    double precision TAU, TRANSCELL, ABSCELL, TRANSMIT
    double precision RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/0, 0, 0, 0, 1, 5/
    data JOCTORDER/1, 1, 1, 1, 2, 2, 2, 2/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    PI = acos(-1.0)
    CZ = -MU
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    ! BITc is 0 for positive direction trace back ray, 1 for negativ
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [BACK_INT_GRID3D]: Bad MU.'
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate
    IOCT = 1 + 4*BITZ
    JOCT = JOCTORDER(IOCT)

    ! Sweep through all the grid points
    do IORDER = 1, NPTS
        IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
        ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
        IPT = GRIDPTR(ICORNER, IPCELL)
        if (GRIDRAD(1, IPT) .ge. 0.0) goto 500

        ! Start off at the grid point
        ICELL = IPCELL
        TRANSMIT = 1.0
        EXT1 = EXTINCT(IPT)
        RAD(:) = 0.0
        SRCEXT1(:) = EXT1*SOURCE(:, KANG, IPT)
        XE = GRIDPOS(1, IPT)
        YE = GRIDPOS(2, IPT)
        ZE = GRIDPOS(3, IPT)

        ! Loop until finding a face with known radiances
        VALIDRAD = .false.
        do while (.not. VALIDRAD)
            ! Make sure current cell is valid
            if (ICELL .le. 0) then
                write (6, *) 'BACK_INT_GRID: ICELL=0  ',&
                    & MU, PHI, IPT, INEXTCELL
                stop
            end if
            ! Find boundaries of the current cell
            ! Find the three possible intersection planes (X,Y,Z)
            !   from the coordinates of the opposite corner grid point
            IOPP = GRIDPTR(9 - IOCT, ICELL)
            ! Get the distances to the 3 planes and select the closest
            SO = (GRIDPOS(3, IOPP) - ZE)*CZINV
            if (SO .lt. -EPS) then
                write (6, *) 'BACK_INT_GRID1D: SO<0  ',&
                    & MU, PHI, ZE, SO, IPT, ICELL
                stop
            end if
            ! Compute the coordinates of the cell exitting location
            ZE = ZE + SO*CZ

            ! Get the intersection face number (i.e. neighptr index)
            IFACE = 6 - BITZ
            ! Get the next cell to go to
            INEXTCELL = NEIGHPTR(IFACE, ICELL)

            ! Get the face grid pointers
            I1 = GRIDPTR(GRIDFACE(IFACE), ICELL)
            ! Get the location coordinate (does the boundary wrapping)
            if (INEXTCELL .gt. 0) then
                ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
            end if

            ! Get extinction and source function at face intersection
            EXT0 = EXTINCT(I1)
            SRCEXT0(:) = EXTINCT(I1)*SOURCE(:, KANG, I1)
            ! Compute the cell radiance: integration of the source funct
            EXT = 0.5*(EXT0 + EXT1)
            TAU = EXT*SO
            if (TAU .ge. 0.5) then
                TRANSCELL = exp(-TAU)
                ABSCELL = 1.0 - TRANSCELL
            else
                ABSCELL = TAU*(1.0 - 0.5*TAU*&
                    & (1.0 - 0.33333333333*TAU*(1 - 0.25*TAU)))
                TRANSCELL = 1.0 - ABSCELL
            end if
            if (TAU .le. 2.0) then
                if (EXT .eq. 0.0) then
                    SRC(:) = 0.0
                else
                    ! Linear extinction, linear source*extinction, to first
                    SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + 0.08333333333&
                        & *(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
                end if
            else
                ! Combined first order expansion and constant extinction f
                EXT0P = EXT0
                SRCEXT0P(:) = SRCEXT0(:)
                if (TAU .gt. 4.0) then
                    EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                    if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
                end if
                SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                    & + (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                    & *(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
            end if
            SRC(1) = max(SRC(1), 0.0d0)

            ! See if there are valid radiances for the intersection face
            ! If so, set flag, interpolate radiance, and add to
            !   radiance from source function in cell
            ! If not, only add in cell radiance and prepare for next cel
            if (GRIDRAD(1, I1) .ge. -0.1) then
                VALIDRAD = .true.
                RAD0(:) = GRIDRAD(:, I1)
                RAD(:) = RAD(:) + TRANSMIT*(RAD0(:)*TRANSCELL + SRC(:)*ABSCELL)
            else
                RAD(:) = RAD(:) + TRANSMIT*SRC(:)*ABSCELL
                SRCEXT1(:) = SRCEXT0(:)
                TRANSMIT = TRANSMIT*TRANSCELL
                EXT1 = EXT0
                ICELL = INEXTCELL
            end if
        end do

        GRIDRAD(:, IPT) = RAD(:)
500     continue
    end do
    return
end subroutine

subroutine NEXT_CELL(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXT)
    !   Gets the next cell (INEXT) to go to that is adjacent to this fac
    ! IFACE (1-6) has which face (-X,+X,-Y,+Y,-Z,+Z) of the current cell
    ! (ICELL) we are at, and JFACE (1-3) is the direction (X,Y,Z).
    ! The current location (XE,YE,ZE) on the face is used to find the ce
    ! if there are several adjacent to the current cell.
    implicit none
    integer IFACE, JFACE, ICELL, INEXT
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    double precision XE, YE, ZE
    integer IC, IC1, DIR

    ! Get the pointer to the neighbor cell
    INEXT = NEIGHPTR(IFACE, ICELL)
    if (INEXT .lt. 0) then
        ! If neighbor pointer is negative then we are going into a cel
        !   that has subcells, and the absolute value of the pointer is
        !   cell to start the tree search.
        IC = abs(INEXT)
        ! Go down the tree until there is no children
        do while (TREEPTR(2, IC) .gt. 0)
            ! Get the direction cell is split
            DIR = ibits(int(CELLFLAGS(IC)), 2, 2)
            if (DIR .eq. 0) stop 'Error [NEXT_CELL]: No split direction.'
            ! If split is in direction of face then choose the near chil
            !   which is the negative side cell if this is a positive fa
            IC1 = TREEPTR(2, IC)
            if (DIR .eq. JFACE) then
                IC = IC1 + 1 - mod(IFACE - 1, 2)
            else
                ! If split is not in direction of face then find which of tw
                !   has a face that contains the current location.
                ! Compare the appropriate coordinate (X,Y,Z) of the grid p
                !   on the split line with the current location.
                IC = IC1
                if (DIR .eq. 1) then
                    if (XE .gt. GRIDPOS(1, GRIDPTR(8, IC1))) IC = IC + 1
                else if (DIR .eq. 2) then
                    if (YE .gt. GRIDPOS(2, GRIDPTR(8, IC1))) IC = IC + 1
                else
                    if (ZE .gt. GRIDPOS(3, GRIDPTR(8, IC1))) IC = IC + 1
                end if
            end if
        end do
        INEXT = IC
    end if

    return
end subroutine

logical function SWEEP_BASE_CELL(BCFLAG, NXC, NYC, NZ, IOCT, ICELL)
    !   Returns the next base grid cell.  The octant the ray direction
    ! is coming from (IOCT) determines the order of sweeping.  ICELL is
    ! for the first call to set things up.  The function returns false i
    ! we are done sweeping the grid.  Do loops over IZ, IY, and IX are
    ! simulated.  For open boundary condition do the endpoint cells
    ! (in X,Y) first, so the boundary grid points will be done with
    ! independent pixel radiative transfer.
    implicit none
    integer BCFLAG, NXC, NYC, NZ, IOCT, ICELL
    logical BTEST
    integer IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ
    save IX, IY, IZ, SIX, SIY, SIZ, EIX, EIY, EIZ, DIX, DIY, DIZ

    SWEEP_BASE_CELL = .true.
    ! If first time here then set up the loops
    if (ICELL .eq. 0) then
        if (btest(IOCT - 1, 0)) then
            DIX = +1
            if (btest(BCFLAG, 0)) then
                SIX = NXC
                EIX = max(NXC - 1, 1)
            else
                SIX = 1
                EIX = NXC
            end if
        else
            DIX = -1
            if (btest(BCFLAG, 0)) then
                SIX = 1
                EIX = min(2, NXC)
            else
                SIX = NXC
                EIX = 1
            end if
        end if
        if (btest(IOCT - 1, 1)) then
            DIY = +1
            if (btest(BCFLAG, 1)) then
                SIY = NYC
                EIY = max(NYC - 1, 1)
            else
                SIY = 1
                EIY = NYC
            end if
        else
            DIY = -1
            if (btest(BCFLAG, 1)) then
                SIY = 1
                EIY = min(2, NYC)
            else
                SIY = NYC
                EIY = 1
            end if
        end if
        if (btest(IOCT - 1, 2)) then
            DIZ = +1
            SIZ = 1
            EIZ = NZ - 1
        else
            DIZ = -1
            SIZ = NZ - 1
            EIZ = 1
        end if
        IX = SIX
        IY = SIY
        IZ = SIZ
    else
        ! Sweep through grid, first in X, in Y, then in Z
        if (IX .eq. EIX) then
            if (IY .eq. EIY) then
                if (IZ .eq. EIZ) then
                    SWEEP_BASE_CELL = .false.
                    return
                else
                    IZ = IZ + DIZ
                end if
                IY = SIY
            else
                IY = mod(IY + DIY + NYC - 1, NYC) + 1
            end if
            IX = SIX
        else
            IX = mod(IX + DIX + NXC - 1, NXC) + 1
        end if
    end if

    ! Compute the base grid cell
    ICELL = IZ + (NZ - 1)*(IY - 1) + (NZ - 1)*NYC*(IX - 1)
    return
end function

logical function SWEEP_NEXT_CELL(BCFLAG, NXC, NYC, NZ, IOCT,&
    & ICELL, TREEPTR, CELLFLAGS)
    !   Returns the next grid cell given the last cell (ICELL).  The
    ! octant the ray direction is coming from (IOCT) determines the
    ! order of sweeping.  ICELL is 0 for the first call to set things up
    ! The function returns false if we are done sweeping the grid.
    ! Within a base cell the tree structure is traced using a stack unti
    ! we get to a new end cell.
    implicit none
    integer BCFLAG, NXC, NYC, NZ, IOCT, ICELL, TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    integer IC, IDIR, SP, STACK(50)
    logical DONE, SWEEP_BASE_CELL, BTEST
    save SP, STACK

    IC = ICELL

    ! Special case for first time here
    if (ICELL .eq. 0) then
        IC = 0
        SWEEP_NEXT_CELL = SWEEP_BASE_CELL(BCFLAG, NXC, NYC, NZ, IOCT, IC)
        SP = 0
    end if
    ! The scanning direction is the direction (pos or neg) for each
    !   of X,Y,Z that the ray is going in (indicated by IOCT).

    ! Go through the tree checking the cells:
    ! The stack holds the second child cells, which are the
    !   next cell to go to if at an end cell.
    DONE = .false.
    do while (.not. DONE)
        ! If this is an end cell then we are either at the desired cell
        !   or if this is the old cell then we need to pop the stack
        if (TREEPTR(2, IC) .eq. 0) then
            if (IC .ne. ICELL) then
                DONE = .true.
            else
                ! If no more to pop then we are at the base cell, so go to
                !   next base cell if there is one.
                if (SP .eq. 0) then
                    SWEEP_NEXT_CELL =&
                        & SWEEP_BASE_CELL(BCFLAG, NXC, NYC, NZ, IOCT, IC)
                    if (.not. SWEEP_NEXT_CELL) return
                else
                    ! Otherwise, pop the stack to get the next cell to go to.
                    IC = STACK(SP)
                    SP = SP - 1
                end if
            end if
        else
            ! Not an end cell, so this cell has children.
            ! Depending on the ordinate octant, decide which child to
            !   do first and put the second child on the stack.
            SP = SP + 1
            if (SP .gt. 50) stop 'Error [SWEEP_NEXT_CELL]: Stack exceeded.'
            IDIR = ibits(int(CELLFLAGS(IC)), 2, 2) - 1
            if (btest(IOCT - 1, IDIR)) then
                STACK(SP) = TREEPTR(2, IC) + 1
                IC = TREEPTR(2, IC)
            else
                STACK(SP) = TREEPTR(2, IC)
                IC = TREEPTR(2, IC) + 1
            end if
        end if
    end do

    ICELL = IC
    SWEEP_NEXT_CELL = .true.
    return
end function

subroutine SPLIT_GRID(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
    & DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
    & NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & GRIDPOS, NX, XGRID, NY, YGRID,&
    & NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    & DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
    & EXTINCT, ALBEDO, LEGEN, IPHASE,&
    & RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
    & TEMP, PLANCK, DIRFLUX, SRCTYPE,&
    & SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    & UNITS, WAVENO, WAVELEN, ADAPTCRIT, ADAPTIND)
    !   Splits the cells that have a cell dividing criterion greater tha
    ! CURSPLITACC if DOSPLIT is true.  The current splitting criterion
    ! achieved is returned in SPLITCRIT.  If we are at the end of the
    ! grid point, grid cell, or spherical harmonic arrays (sizes given
    ! by MAXIG, MAXIC, MAXIV, MAXIDO) then the OUTOFMEM flag is returned
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
    character SRCTYPE*1, UNITS*1
    integer ICELL, ICELL1, IDIR, I, N, NEWPOINTS(3, 4), OLDNPTS
    integer MAXCELLS, MAXPTS, MAXWORK, MAXSH, MAXRAD
    logical OUTOFMEM0
    real ADAPT(3), FRAC

    if (DOSPLIT) then
        OLDNPTS = NPTS
        ! Do the cells in batchs, so the newly split cells can be teste
        !   for further splitting
        ICELL1 = 1
        do while (ICELL1 .le. NCELLS)
            ! For the current batch of cells, compute the adaptive grid
            !   splitting criterion for 3 directions for each cell.
            N = 0
            do ICELL = ICELL1, NCELLS
                if (TREEPTR(2, ICELL) .eq. 0) then
                    ! If this is a end node cell then compute adaptive criteri
                    call CELL_SPLIT_TEST(GRIDPTR, GRIDPOS, EXTINCT,&
                        & SHPTR, NSTOKES, SOURCE, ICELL,&
                        & ADAPT, ADAPTCRIT(N + 1), IDIR)
                    ADAPTIND(N+1) = 4*ICELL + IDIR
                    N = N + 1
                end if
            end do
            ICELL1 = NCELLS + 1

            ! Sort the adaptive criterion in decreasing order
            call SSORT(ADAPTCRIT, ADAPTIND, N, -2)

            ! Figure how many cells we can split, leaving some room
            !   for cells split during grid smoothing
            FRAC = 0.03
            MAXCELLS = MAXIC - FRAC*(MAXIC - NCELLS) - 2
            MAXPTS = MAXIG - FRAC*(MAXIG - NPTS) - 4
            MAXWORK = MAXIDO - FRAC*(MAXIDO - NPHI0MAX*NPTS) - 4*NPHI0MAX
            MAXSH = MAXIV - FRAC*(MAXIV - SHPTR(NPTS + 1)) - 4*NLM
            MAXRAD = MAXIV + MAXIG - FRAC*(MAXIV + MAXIG - RSHPTR(NPTS + 1)) - 4*NLM

            ! Go down the list from highest adaptive cell criterion on d
            OUTOFMEM0 = OUTOFMEM
            I = 1
            do while (ADAPTCRIT(I) .gt. CURSPLITACC .and. .not. OUTOFMEM0&
                & .and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    & NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    & SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    & RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM0 = .true.
                else
                    ICELL = int(ADAPTIND(I))/4
                    IDIR = iand(int(ADAPTIND(I)), 3)
                    ! Create the new grid points and make pointers
                    call DIVIDE_CELL(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                        & TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                    ! Interpolate the medium properties, radiance, and sourc
                    call INTERPOLATE_POINT(NEWPOINTS, NPTS,&
                        & NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                        & DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                        & SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                        & UNITS, WAVENO, WAVELEN,&
                        & EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                        & RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
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
                    & NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    & SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    & RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM = .true.
                else
                    ! See if cell should be split to make better grid:
                    !   IDIR=0 if no, IDIR=direction (1,2,3) otherwise
                    ICELL = int(ADAPTIND(I))/4
                    call GRID_SMOOTH_TEST(GRIDPTR, GRIDPOS, NEIGHPTR,&
                        & TREEPTR, CELLFLAGS, SHPTR, ICELL, IDIR)
                    ! Set IDIR=0 to stop adaptive grid smoothing
                    ! IDIR = 0
                    if (IDIR .gt. 0) then
                        ! Create the new grid points and make pointers
                        call DIVIDE_CELL(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                            & TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                        ! Interpolate the medium properties, radiance, and sou
                        call INTERPOLATE_POINT(NEWPOINTS, NPTS,&
                            & NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            & DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                            & SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                            & UNITS, WAVENO, WAVELEN,&
                            & EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                            & RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                    end if
                end if
                I = I + 1
            end do

        end do
        if (OUTOFMEM0) then
            OUTOFMEM = .true.
            write (6, '(2A,I8,A,I8)')&
                & 'Further cell splitting requires larger array allocations:',&
                & '  Ncells=', NCELLS, '  Npnts=', NPTS
        end if

        ! If in parallel mode, calculate all the direct beam fluxes at on
        if (SRCTYPE .ne. 'T' .and.&
            & (btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) then
            call MAKE_DIRECT_PAR(OLDNPTS + 1, NPTS, BCFLAG,&
                & IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
                & SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
                & NX, XGRID, NY, YGRID, DIRFLUX)
        end if
    end if

    ! Find the max adaptive cell criterion after the cell divisions
    SPLITCRIT = 0.0
    do ICELL = 1, NCELLS
        if (TREEPTR(2, ICELL) .eq. 0) then
            call CELL_SPLIT_TEST(GRIDPTR, GRIDPOS, EXTINCT,&
                & SHPTR, NSTOKES, SOURCE, ICELL,&
                & ADAPT, ADAPTCRIT(1), IDIR)
            SPLITCRIT = max(SPLITCRIT, ADAPTCRIT(1))
        end if
    end do

    return
end subroutine

subroutine INTERPOLATE_POINT(NEWPOINTS, NPTS,&
    & NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    & DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
    & SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    & UNITS, WAVENO, WAVELEN,&
    & EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
    & RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
    !   Interpolates the medium properties, radiance, and source functio
    ! new grid points (specified in NEWPOINTS).  The medium properties a
    ! interpolated from the property grid, since interpolating from the
    ! internal grid may not faithfully represent the true medium.
    ! The delta-M scaling is done to the new grid point properties if ne
    ! The direct beam flux and the Planck function for the new points is
    ! recomputed (rather than interpolated).  The radiance is averaged
    ! from the two parent grid points, and used to compute the source
    ! function for the new points.
    implicit none
    integer NEWPOINTS(3, 4), NPTS
    integer NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer BCFLAG, IPFLAG
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    integer*2 IPHASE(NPTS)
    logical DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(NSTLEG, NLM)
    real WAVENO(2), WAVELEN
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real TEMP(NPTS), PLANCK(NPTS), DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *), SOURCE(NSTOKES, *)
    real GRIDPOS(3, *)
    character SRCTYPE*1, UNITS*1
    integer I, IP1, IP2, IP, IPH, J, K, L, M, ME, MS
    integer IS, NS, NS1, NS2
    integer IR, IR1, IR2, NR, NR1, NR2, SIDE
    logical VALIDBEAM
    real X, Y, Z, F, BB, C, D, SECMU0
    real RAD1(NSTOKES), RAD2(NSTOKES)
    real UNIFZLEV, XO, YO, ZO, DIRPATH
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))
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

    C = sqrt(4.0*acos(-1.0))
    if (SRCTYPE .ne. 'T') then
        SECMU0 = 1.0/abs(SOLARMU)
    end if

    ! Loop over the (potentially) four new grid points
    do I = 1, 4
        if (NEWPOINTS(3, I) .gt. 0) then
            ! Get the parent and new grid pointers
            IP1 = NEWPOINTS(1, I)
            IP2 = NEWPOINTS(2, I)
            IP = NEWPOINTS(3, I)

            ! Interpolate the medium properties from the property grid
            X = GRIDPOS(1, IP)
            Y = GRIDPOS(2, IP)
            Z = GRIDPOS(3, IP)
            call TRILIN_INTERP_PROP(X, Y, Z, .false., NSTLEG, NLEG,&
                & TEMP(IP), EXTINCT(IP), ALBEDO(IP),&
                & LEGEN(1, 0, IP), IPHASE(IP))
            ! Do the Delta-M scaling of extinction and albedo for this p
            if (DELTAM) then
                if (NUMPHASE .gt. 0) then
                    F = LEGEN(1, ML + 1, IPHASE(IP))
                else
                    F = LEGEN(1, ML + 1, IP)
                    do L = 0, ML
                        LEGEN(1, L, IP) = (LEGEN(1, L, IP) - F)/(1 - F)
                    end do
                end if
                EXTINCT(IP) = (1.0 - ALBEDO(IP)*F)*EXTINCT(IP)
                ALBEDO(IP) = (1.0 - F)*ALBEDO(IP)/(1.0 - ALBEDO(IP)*F)
            end if

            ! Compute the new Planck source function (if needed)
            if (SRCTYPE .ne. 'S') then
                call PLANCK_FUNCTION(TEMP(IP), UNITS, WAVENO, WAVELEN, BB)
                PLANCK(IP) = (1.0 - ALBEDO(IP))*BB
            end if

            ! Compute the new direct beam flux (if needed)
            if (SRCTYPE .ne. 'T') then
                if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
                    ! If parallel mode then calculate a temporary approximatio
                    !   to the direct flux (for source function calculation be
                    DIRFLUX(IP) = exp(0.5*(log(DIRFLUX(IP1)) + log(DIRFLUX(IP2))))
                else
                    ! Otherwise, calculate the exact direct beam from property
                    DIRPATH = 0.0
                    call DIRECT_BEAM_PROP(0, X, Y, Z, BCFLAG, IPFLAG, DELTAM,&
                        & ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                        & DIRFLUX(IP),&
                        & UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                end if
            end if

            ! Interpolate the radiance
            IR1 = RSHPTR(IP1)
            IR2 = RSHPTR(IP2)
            NR1 = RSHPTR(IP1 + 1) - RSHPTR(IP1)
            NR2 = RSHPTR(IP2 + 1) - RSHPTR(IP2)
            NR = max(NR1, NR2)
            IR = RSHPTR(IP)
            RSHPTR(IP + 1) = IR + NR
            do J = 1, NR
                if (J .le. NR1) then
                    RAD1(:) = RADIANCE(:, IR1 + J)
                else
                    RAD1(:) = 0.0
                end if
                if (J .le. NR2) then
                    RAD2(:) = RADIANCE(:, IR2 + J)
                else
                    RAD2(:) = 0.0
                end if
                RADIANCE(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
            end do

            ! Compute the source function for the new points
            NS1 = SHPTR(IP1 + 1) - SHPTR(IP1)
            NS2 = SHPTR(IP2 + 1) - SHPTR(IP2)
            NS = max(NS1, NS2)
            IS = SHPTR(IP)
            SHPTR(IP + 1) = IS + NS
            if (ACCELFLAG) OSHPTR(IP + 1) = OSHPTR(IP)
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(IP)
            else
                IPH = IP
            end if
            if (NSTOKES .eq. 1) then
                call CALC_SOURCE_PNT_UNPOL(NLM, NLEG, LOFJ,&
                    & SRCTYPE, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    & PLANCK(IP), ALBEDO(IP), IPH, LEGEN,&
                    & NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call CALC_SOURCE_PNT(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    & SRCTYPE, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    & PLANCK(IP), ALBEDO(IP), IPH, LEGEN,&
                    & NR, RADIANCE(1, IR + 1), SOURCET)
            end if
            do J = 1, NS
                SOURCE(:, IS + J) = SOURCET(:, J)
            end do
        end if
    end do
    deallocate (LOFJ, SOURCET)

    return
end subroutine

subroutine DIVIDE_CELL(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
    & TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
    !   Divides a single cell (ICELL) in the IDIR direction (1=Z,2=Y,3=X
    ! Makes the pointers to the two new cells (TREEPTR), set the flags
    ! (CELLFLAGS), creates the new grid points if they don't exist alrea
    ! and updates the neighbor pointers (NEIGHPTR) for the new cells and
    ! their neighboring cells.
    implicit none
    integer NPTS, NCELLS, ICELL, IDIR, NEWPOINTS(3, 4)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer NEWCELL, IFACE, I
    logical BTEST

    ! Only allow to divide an end cell (i.e. one not yet divided)
    if (TREEPTR(2, ICELL) .ne. 0) then
        stop 'Error [DIVIDE_CELL]: Cannot divide already split cell.'
    end if

    ! Do the tree pointers: parent points to negative child,
    !   children point back to parent.
    NEWCELL = NCELLS + 1
    NCELLS = NCELLS + 2
    TREEPTR(2, ICELL) = NEWCELL
    TREEPTR(1, NEWCELL) = ICELL
    TREEPTR(2, NEWCELL) = 0
    TREEPTR(1, NEWCELL + 1) = ICELL
    TREEPTR(2, NEWCELL + 1) = 0
    ! Do the flags: put the splitting direction in the parent and
    !   inherit the boundary flags for the children.
    CELLFLAGS(ICELL) = ior(int(CELLFLAGS(ICELL)), ishft(IDIR, 2))
    CELLFLAGS(NEWCELL) = 0
    CELLFLAGS(NEWCELL + 1) = 0
    do I = 0, 1
        if (btest(int(CELLFLAGS(ICELL)), I)) then
            CELLFLAGS(NEWCELL) = ibset(int(CELLFLAGS(NEWCELL)), I)
            CELLFLAGS(NEWCELL + 1) = ibset(int(CELLFLAGS(NEWCELL + 1)), I)
        end if
    end do

    ! Do the grid points:
    call NEW_GRID_POINTS(IDIR, ICELL, NEWCELL, NPTS,&
        & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
        & GRIDPOS, NEWPOINTS)

    ! Do the neighbor pointers:
    do IFACE = 1, 6
        ! Special case if the neighbor is the cell itself
        if (NEIGHPTR(IFACE, ICELL) .eq. ICELL) then
            NEIGHPTR(IFACE, NEWCELL) = NEWCELL
            NEIGHPTR(IFACE, NEWCELL + 1) = NEWCELL + 1
        ! For the faces in the split direction, the children either
        !   point to the other child or inherit the parent's neighbor.
        !   Also update the neighbor's neighbor pointer.
        else if (IFACE .eq. 2*IDIR) then
            NEIGHPTR(IFACE, NEWCELL) = NEWCELL + 1
            NEIGHPTR(IFACE, NEWCELL + 1) = NEIGHPTR(IFACE, ICELL)
        else if (IFACE .eq. 2*IDIR - 1) then
            NEIGHPTR(IFACE, NEWCELL + 1) = NEWCELL
            NEIGHPTR(IFACE, NEWCELL) = NEIGHPTR(IFACE, ICELL)
        else
            ! For the other four faces:
            !   The children inherit the parent's neighbor pointers
            NEIGHPTR(IFACE, NEWCELL) = NEIGHPTR(IFACE, ICELL)
            NEIGHPTR(IFACE, NEWCELL + 1) = NEIGHPTR(IFACE, ICELL)
        end if
        ! Then see if we can match faces further down the neighbor's tre
        !   also update the neighboring cells pointers
        call MATCH_NEIGHBOR_FACE(IFACE, NEWCELL, NEIGHPTR,&
            & CELLFLAGS, TREEPTR, GRIDPTR, GRIDPOS)
        call MATCH_NEIGHBOR_FACE(IFACE, NEWCELL + 1, NEIGHPTR,&
            & CELLFLAGS, TREEPTR, GRIDPTR, GRIDPOS)
    end do

    return
end subroutine

subroutine MATCH_NEIGHBOR_FACE(IFACE, IC, NEIGHPTR, CELLFLAGS,&
    & TREEPTR, GRIDPTR, GRIDPOS)
    !   Resolves the neighbor pointers for a new cell IC (from a subdivi
    ! for the IFACE face.  Based on the size of this new face and
    ! the corresponding faces of the neighboring cells (down its tree),
    ! the new cell's NEIGHPTR and those of the neighboring cells are upd
    implicit none
    integer IFACE, IC, NEIGHPTR(6, *), TREEPTR(2, *), GRIDPTR(8, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer IN, INN, JFACE, IC1, IC8, IN1, IN8, DIR, DIR1, DIR2
    logical DONE
    real POS(3)
    integer OPPFACE(6)
    data OPPFACE/2, 1, 4, 3, 6, 5/

    IN = abs(NEIGHPTR(IFACE, IC))
    ! If no neighbor then return
    if (IN .eq. 0) return
    JFACE = (IFACE + 1)/2
    ! Get the center of the new face
    IC1 = GRIDPTR(1, IC)
    IC8 = GRIDPTR(8, IC)
    DIR1 = mod(JFACE - 1 + 1, 3) + 1
    DIR2 = mod(JFACE - 1 + 2, 3) + 1
    POS(DIR1) = (GRIDPOS(DIR1, IC1) + GRIDPOS(DIR1, IC8))/2
    POS(DIR2) = (GRIDPOS(DIR2, IC1) + GRIDPOS(DIR2, IC8))/2
    ! Go down the neighbors cell tree
    DONE = .false.
    do while (.not. DONE .and. IC .ne. IN)
        ! Compare the coordinates of this neighbor cell to see
        !   if the new face is contained within the neighbor face.
        IN1 = GRIDPTR(1, IN)
        IN8 = GRIDPTR(8, IN)
        if (GRIDPOS(DIR1, IC1) .ge. GRIDPOS(DIR1, IN1) .and.&
            & GRIDPOS(DIR1, IC8) .le. GRIDPOS(DIR1, IN8) .and.&
            & GRIDPOS(DIR2, IC1) .ge. GRIDPOS(DIR2, IN1) .and.&
            & GRIDPOS(DIR2, IC8) .le. GRIDPOS(DIR2, IN8)) then
            ! If this is the bottom of the tree then we have a face matc
            if (TREEPTR(2, IN) .eq. 0) then
                NEIGHPTR(IFACE, IC) = IN
            else
                NEIGHPTR(IFACE, IC) = -IN
            end if
        end if
        ! Compare the size of the neighbor face with the new face.
        !   Update the neighbor's neighbor pointer to the new face
        !   if the neighbor face is smaller, otherwise make the
        !   pointer negative, because it no longer points to an end no
        if (GRIDPOS(DIR1, IN1) .ge. GRIDPOS(DIR1, IC1) .and.&
            & GRIDPOS(DIR1, IN8) .le. GRIDPOS(DIR1, IC8) .and.&
            & GRIDPOS(DIR2, IN1) .ge. GRIDPOS(DIR2, IC1) .and.&
            & GRIDPOS(DIR2, IN8) .le. GRIDPOS(DIR2, IC8)) then
            call INHERIT_NEIGHBOR(IN, OPPFACE(IFACE), IC,&
                & NEIGHPTR, TREEPTR, CELLFLAGS)
        else
            NEIGHPTR(OPPFACE(IFACE), IN) = -abs(NEIGHPTR(OPPFACE(IFACE), IN))
        end if

        ! Go down the neighbors cell tree until there is no children,
        !   finding the cell that contains the center of the new face.
        if (TREEPTR(2, IN) .eq. 0) then
            DONE = .true.
        else
            DIR = ibits(int(CELLFLAGS(IN)), 2, 2)
            INN = TREEPTR(2, IN)
            if (DIR .eq. JFACE) then
                IN = INN + 1 - mod(IFACE - 1, 2)
            else
                if (POS(DIR) .gt. GRIDPOS(DIR, GRIDPTR(8, INN))) then
                    IN = INN + 1
                else
                    IN = INN
                end if
            end if
        end if
    end do

    return
end subroutine

subroutine INHERIT_NEIGHBOR(ICELL, IFACE, IN,&
    & NEIGHPTR, TREEPTR, CELLFLAGS)
    !   Updates the ICELL cell's new neighbor cell (IN) for all
    ! children cells that border the face (IFACE).
    implicit none
    integer ICELL, IFACE, IN, NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    integer JFACE, IC, DIR, SP, STACK(50)
    logical DONE

    ! Go through the tree checking the cells:
    ! The stack holds the second child cells, which are the
    !   next cell to go to if at an end cell.
    JFACE = (IFACE + 1)/2
    IC = ICELL
    SP = 0
    DONE = .false.
    do while (.not. DONE)
        NEIGHPTR(IFACE, IC) = IN
        ! If this is an end cell then pop the stack
        if (TREEPTR(2, IC) .eq. 0) then
            ! If no more to pop then we are done
            if (SP .eq. 0) then
                DONE = .true.
            else
                ! Otherwise, pop the stack to get the next cell to go to
                IC = STACK(SP)
                SP = SP - 1
            end if
        else
            ! Not and end cell, so this cell has children.
            ! If split is in same direction as face then do only close c
            DIR = ibits(int(CELLFLAGS(IC)), 2, 2)
            if (DIR .eq. JFACE) then
                IC = TREEPTR(2, IC) + mod(IFACE - 1, 2)
            else
                ! Otherwise need to go to both children so put the second ch
                !   on the stack and set the pointer to go to the first chil
                SP = SP + 1
                if (SP .gt. 50) stop 'Error [INHERIT_NEIGHBOR]: Stack exceeded.'
                STACK(SP) = TREEPTR(2, IC) + 1
                IC = TREEPTR(2, IC)
            end if
        end if
    end do

    return
end subroutine

subroutine NEW_GRID_POINTS(IDIR, ICELL, NEWCELL, NPTS,&
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & GRIDPOS, NEWPOINTS)
    !   Sets up the pointers to the grid points for the two new cells.
    ! Creates the four new points if they don't already exist, changing
    ! NPTS and GRIDPOS accordingly.  Returns the NEWPOINTS(3,4) arrays
    ! that gives information about the new points: NEWPOINTS(3,) is the
    ! point address of the new points (or zero), and NEWPOINTS(1,) and
    ! NEWPOINTS(2,) are the "parent" grid points these should be
    ! interpolated from.
    implicit none
    integer IDIR, ICELL, NEWCELL, NPTS, NEWPOINTS(3, 4)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer I, K, I1, I2, IP1, IP2, IPT, IPMATCH
    integer IFACE1, IFACE2, ICELL2
    integer GRIDCORNER(2, 4, 3), FACEGRID(2, 4, 3)
    real XP, YP, ZP
    data GRIDCORNER/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        & 1, 5, 2, 6, 3, 7, 4, 8/
    data FACEGRID/3, 5, 4, 5, 3, 6, 4, 6, 1, 5, 2, 5, 1, 6, 2, 6,&
        & 1, 3, 2, 3, 1, 4, 2, 4/

    !   In general, when subdividing a cell four new points may be made.
    ! Each of these points belongs to two faces of the original cell.
    ! Use the MATCH_GRID_POINT routine to try to find matches in cells
    ! across the two faces that contain the point; also cross the face o
    ! one of these to get the third bounding cell.  If any of the approp
    ! gridpoints of these three cells has the same coordinates then the
    ! already exists, otherwise create the point.

    ! Temporarily set new cells' grid pointer to the parent cell's o
    do K = 1, 8
        GRIDPTR(K, NEWCELL) = GRIDPTR(K, ICELL)
        GRIDPTR(K, NEWCELL + 1) = GRIDPTR(K, ICELL)
    end do

    ! Loop over the four edges with potentially new points
    do I = 1, 4
        ! Get the grid corners and coordinates for this edge and split
        I1 = GRIDCORNER(1, I, IDIR)
        I2 = GRIDCORNER(2, I, IDIR)
        IP1 = GRIDPTR(I1, ICELL)
        IP2 = GRIDPTR(I2, ICELL)
        XP = (GRIDPOS(1, IP1) + GRIDPOS(1, IP2))/2
        YP = (GRIDPOS(2, IP1) + GRIDPOS(2, IP2))/2
        ZP = (GRIDPOS(3, IP1) + GRIDPOS(3, IP2))/2

        ! Look for gridpoints matching this new point in the three cell
        !   that border this edge
        IFACE1 = FACEGRID(1, I, IDIR)
        call MATCH_GRID_POINT(XP, YP, ZP, ICELL, IFACE1, GRIDPOS,&
            & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
        if (IPMATCH .eq. 0) then
            IFACE2 = FACEGRID(2, I, IDIR)
            call MATCH_GRID_POINT(XP, YP, ZP, ICELL, IFACE2, GRIDPOS,&
                & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
        end if
        if (IPMATCH .eq. 0) then
            ICELL2 = abs(NEIGHPTR(IFACE1, ICELL))
            if (ICELL2 .gt. 0) then
                call MATCH_GRID_POINT(XP, YP, ZP, ICELL2, IFACE2, GRIDPOS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, IPMATCH)
            end if
        end if

        ! Set some new cells' grid pointers to the new or matching one
        if (IPMATCH .eq. 0) then
            NPTS = NPTS + 1
            GRIDPTR(I2, NEWCELL) = NPTS
            GRIDPTR(I1, NEWCELL + 1) = NPTS
            GRIDPOS(1, NPTS) = XP
            GRIDPOS(2, NPTS) = YP
            GRIDPOS(3, NPTS) = ZP
            NEWPOINTS(1, I) = IP1
            NEWPOINTS(2, I) = IP2
            NEWPOINTS(3, I) = NPTS
        else
            GRIDPTR(I2, NEWCELL) = IPMATCH
            GRIDPTR(I1, NEWCELL + 1) = IPMATCH
            NEWPOINTS(3, I) = 0
        end if
    end do

    return
end subroutine

subroutine MATCH_GRID_POINT(XP, YP, ZP, ICELL, IFACE,&
    & GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    & IPMATCH)
    !   Determines whether the neighboring cells to ICELL adjacent to
    ! face IFACE (1-6: -X,+X,-Y,+Y,-Z,+Z) have gridpoints that match
    ! the point at XP,YP,ZP on face IFACE of cell ICELL.
    ! If the face of ICELL has several neighboring cells, then the tree
    ! is traced down to find those cells that might contain XP,YP,ZP.
    ! IPMATCH is returned 0 if there is no match, otherwise it has the
    ! number of the matching existing gridpoint.
    implicit none
    real XP, YP, ZP
    integer ICELL, IFACE
    real GRIDPOS(3, *)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*)
    integer IPMATCH
    integer MAXSTACK
    parameter(MAXSTACK=64)
    integer CELLSTACK(MAXSTACK), SP
    integer INEXT, IDIR, DIR, IC, IC1, KFACE, IPT, I
    integer OPPFACE(6), GRIDFACE(4, 6)
    data OPPFACE/2, 1, 4, 3, 6, 5/,&
        & GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8/

    IPMATCH = 0
    IDIR = (IFACE + 1)/2
    KFACE = OPPFACE(IFACE)

    ! Get the pointer to the neighbor cell
    IC = abs(NEIGHPTR(IFACE, ICELL))
    ! Don't process if the neighbor cell is a boundary
    if (IC .eq. 0) return

    ! Go down all relevant branches of the tree until there are no chil
    SP = 0
    do while (.true.)

        ! If at an end cell then check for a gridpoint match at the four
        !   corners of the cell's face
        do while (TREEPTR(2, IC) .eq. 0)
            do I = 1, 4
                IPT = GRIDPTR(GRIDFACE(I, KFACE), IC)
                if (XP .eq. GRIDPOS(1, IPT) .and. YP .eq. GRIDPOS(2, IPT)&
                    & .and. ZP .eq. GRIDPOS(3, IPT)) then
                    IPMATCH = IPT
                    return
                end if
            end do
            ! Return if there are no more branches on the stack, otherwise
            if (SP .eq. 0) return
            IC = CELLSTACK(SP)
            SP = SP - 1
        end do

        ! Get the direction current cell is split
        DIR = ibits(int(CELLFLAGS(IC)), 2, 2)
        if (DIR .eq. 0) stop 'Error [MATCH_GRID_POINT]: No split direction.'
        ! If split is in direction of face then choose the near child ce
        !   which is the negative side cell if this is a positive face.
        IC1 = TREEPTR(2, IC)
        if (DIR .eq. IDIR) then
            IC = IC1 + 1 - mod(IFACE - 1, 2)
        else
            !   If split is not in direction of face then find which of two c
            ! has a face that contains the current location by comparing th
            ! appropriate coordinate (X,Y,Z) of the grid point on the split
            ! with the current location.  If the split line matches a coord
            ! then we have to trace down both cells so put one on the stack
            IC = IC1
            if (DIR .eq. 1) then
                if (XP .eq. GRIDPOS(1, GRIDPTR(8, IC1))) then
                    SP = SP + 1
                    CELLSTACK(SP) = IC + 1
                else if (XP .gt. GRIDPOS(1, GRIDPTR(8, IC1))) then
                    IC = IC + 1
                end if
            else if (DIR .eq. 2) then
                if (YP .eq. GRIDPOS(2, GRIDPTR(8, IC1))) then
                    SP = SP + 1
                    CELLSTACK(SP) = IC + 1
                else if (YP .gt. GRIDPOS(2, GRIDPTR(8, IC1))) then
                    IC = IC + 1
                end if
            else
                if (ZP .eq. GRIDPOS(3, GRIDPTR(8, IC1))) then
                    SP = SP + 1
                    CELLSTACK(SP) = IC + 1
                else if (ZP .gt. GRIDPOS(3, GRIDPTR(8, IC1))) then
                    IC = IC + 1
                end if
            end if
            if (SP .ge. MAXSTACK)&
                & stop 'Error [MATCH_GRID_POINT]: MAXSTACK exceeded.'
        end if
    end do

end subroutine

subroutine CELL_SPLIT_TEST(GRIDPTR, GRIDPOS, EXTINCT,&
    & SHPTR, NSTOKES, SOURCE, ICELL,&
    & ADAPTCRIT, MAXADAPT, IDIR)
    !   Computes the adaptive cell splitting criterion for the X, Y, and
    ! directions for cell (ICELL).  The three criterion are returned in
    ! ADAPTCRIT, the max criterion is returned in MAXADAPT, and the
    ! direction of max criterion in IDIR.
    ! The adaptive cell criterion is computed by averaging the division
    ! criterion over the four pairs of grid points that cross the potent
    ! split direction. The cell division criterion is based on the
    ! integration of the source function across the cell (including the
    ! optical depth/transmission factor) that gives the delta radiance.
    ! The criterion is actually the part of the integral that depends on
    ! the difference in the extinction times source function across a ce
    ! The rms over all angles of the extinction/source function differen
    ! is computed in spherical harmonic space by summing the square of
    ! the difference over all terms.
    implicit none
    integer NSTOKES, ICELL, IDIR
    integer SHPTR(*), GRIDPTR(8, *)
    real GRIDPOS(3, *), EXTINCT(*), SOURCE(NSTOKES, *)
    real ADAPTCRIT(3), MAXADAPT
    integer ID, IE, IP1, IP2, IS1, IS2, NS1, NS2, NS
    integer J, NUM
    integer GRIDCORNER(2, 4, 3)
    real SPLIT, SUM1, JAY, TAU, EXT
    real C0
    data GRIDCORNER/1, 2, 3, 4, 5, 6, 7, 8, 1, 3, 2, 4, 5, 7, 6, 8,&
        & 1, 5, 2, 6, 3, 7, 4, 8/
    data C0/0.282095/

    MAXADAPT = -1.0
    ! Loop over the three possible directions to do the split (X,Y,Z
    do ID = 1, 3
        SUM1 = 0.0
        NUM = 0
        ! Loop over four edges that cross the ID direction
        do IE = 1, 4
            ! Get the grid points
            IP1 = GRIDPTR(GRIDCORNER(1, IE, ID), ICELL)
            IP2 = GRIDPTR(GRIDCORNER(2, IE, ID), ICELL)
            if (IP1 .ne. IP2) then
                NUM = NUM + 1
                IS1 = SHPTR(IP1)
                IS2 = SHPTR(IP2)
                NS1 = SHPTR(IP1 + 1) - IS1
                NS2 = SHPTR(IP2 + 1) - IS2
                JAY = 0.0
                NS = min(NS1, NS2)
                !   Compute the source function rms difference integrated
                ! over all angles, which is sum of squares of spherical
                ! harmonic term differences.
                do J = 1, NS
                    JAY = JAY + (EXTINCT(IP2)*SOURCE(1, IS2 + J)&
                        & - EXTINCT(IP1)*SOURCE(1, IS1 + J))**2
                end do
                do J = NS + 1, NS1
                    JAY = JAY + (EXTINCT(IP1)*SOURCE(1, IS1 + J))**2
                end do
                do J = NS + 1, NS2
                    JAY = JAY + (EXTINCT(IP2)*SOURCE(1, IS2 + J))**2
                end do
                EXT = 0.5*(EXTINCT(IP1) + EXTINCT(IP2))
                if (EXT .gt. 0.0) then
                    JAY = C0*sqrt(JAY)/EXT
                else
                    JAY = 0.0
                end if
                TAU = abs(EXT*(GRIDPOS(ID, IP2) - GRIDPOS(ID, IP1)))
                SPLIT = abs(JAY)*(1 - exp(-TAU))
                SUM1 = SUM1 + SPLIT
            end if
        end do

        if (NUM .gt. 0) then
            ADAPTCRIT(ID) = SUM1/NUM
        else
            ADAPTCRIT(ID) = 0.0
        end if
        if (ADAPTCRIT(ID) .gt. MAXADAPT) then
            MAXADAPT = ADAPTCRIT(ID)
            IDIR = ID
        end if
    end do

    return
end subroutine

subroutine GRID_SMOOTH_TEST(GRIDPTR, GRIDPOS, NEIGHPTR,&
    & TREEPTR, CELLFLAGS, SHPTR, ICELL, IDIR)
    !   Tests the grid cell ICELL for splitting in order to improve the
    ! adaptive grid smoothness.  The direction to split the cell is retu
    ! in IDIR (1-X, 2-Y, 3-Z, 0-don't).
    ! The criterion for splitting the cell is
    !   1) cell has finer grid cells on either side
    !   2) a neighboring cell has grid spacing more than 3 times finer
    implicit none
    integer ICELL, IDIR
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *), SHPTR(*)
    integer*2 CELLFLAGS(*)
    real GRIDPOS(3, *)
    integer ID, IFACE, IP, IP1, IP2, IN, I, J, ISUM, IBITS
    integer INCELL(2), DIRSPLIT(2)
    integer FACECORNER(4, 6), EDGECORNER(2, 3)
    double precision XE, YE, ZE
    real RATIO, SIZERATIO, CURGRIDSIZE(3), GRIDSIZE(2), GRIDSIZEINV
    data FACECORNER/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        & 1, 2, 3, 4, 5, 6, 7, 8/
    data EDGECORNER/1, 2, 1, 3, 1, 5/

    IDIR = 0
    SIZERATIO = 1.0
    ! Loop over the three possible directions to do the split (X,Y,Z
    do ID = 1, 3
        ! Get the grid points
        IP1 = GRIDPTR(EDGECORNER(1, ID), ICELL)
        IP2 = GRIDPTR(EDGECORNER(2, ID), ICELL)
        CURGRIDSIZE(ID) = 1.0e20
        if (IP1 .ne. IP2) then
            CURGRIDSIZE(ID) = abs(GRIDPOS(ID, IP1) - GRIDPOS(ID, IP2))
            GRIDSIZEINV = 1.0/CURGRIDSIZE(ID)
            ! Loop over two opposite faces
            do J = 1, 2
                DIRSPLIT(J) = 0
                IFACE = 2*(ID - 1) + J
                ! Get location of center of cell's face
                XE = 0.0
                YE = 0.0
                ZE = 0.0
                do I = 1, 4
                    IP = GRIDPTR(FACECORNER(I, IFACE), ICELL)
                    XE = XE + GRIDPOS(1, IP)*0.25
                    YE = YE + GRIDPOS(2, IP)*0.25
                    ZE = ZE + GRIDPOS(3, IP)*0.25
                end do
                ! Get the neighboring cell and its size in this direction
                call NEXT_CELL(XE, YE, ZE, IFACE, ID, ICELL, GRIDPOS,&
                    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INCELL(J))
                ! If the neighbor does not exist or is a base cell or
                !   is an independent pixel cell in X or Y then set
                !   its grid size the same as current cell (so don't split
                if (INCELL(J) .eq. 0) then
                    GRIDSIZE(J) = CURGRIDSIZE(ID)
                else
                    GRIDSIZE(J) = abs(GRIDPOS(ID, GRIDPTR(1, INCELL(J)))&
                        & - GRIDPOS(ID, GRIDPTR(8, INCELL(J))))

                    if (TREEPTR(1, INCELL(J)) .eq. 0) then
                        GRIDSIZE(J) = CURGRIDSIZE(ID)
                    end if
                    if (ID .le. 2 .and. btest(CELLFLAGS(INCELL(J)), ID - 1)) then
                        GRIDSIZE(J) = CURGRIDSIZE(ID)
                    end if
                end if
                ! If there are multiple neighboring cells then get
                !   direction the neighbor cell was split
                IN = NEIGHPTR(IFACE, ICELL)
                if (IN .lt. 0) then
                    DIRSPLIT(J) = ibits(int(CELLFLAGS(abs(IN))), 2, 2)
                end if
            end do

            ! Split cell this way if cells on opposite side are smaller
            !   this direction OR opposite neighbor cells are split in
            !   same other direction.
            if (GRIDSIZE(1)*GRIDSIZEINV .lt. 0.75 .and.&
                & GRIDSIZE(2)*GRIDSIZEINV .lt. 0.75) then
                IDIR = ID
            end if
            if (DIRSPLIT(1) .gt. 0 .and. DIRSPLIT(1) .ne. ID .and.&
                & DIRSPLIT(1) .eq. DIRSPLIT(2)) then
                IDIR = DIRSPLIT(1)
            end if

            ! Split cell this way if ratio of neighboring grid cells' si
            !   to this cells size (in ID direction) is less than 0.4.
            RATIO = min(GRIDSIZE(1)*GRIDSIZEINV, GRIDSIZE(2)*GRIDSIZEINV)
            if (RATIO .lt. 0.4 .and. RATIO .lt. SIZERATIO) then
                ! Check if there is any source function for this cell, bef
                !   deciding to split it (no point in dividing clear sky).
                ISUM = 0
                do I = 1, 8
                    IP = GRIDPTR(I, ICELL)
                    ISUM = ISUM + SHPTR(IP + 1) - SHPTR(IP)
                end do
                if (ISUM .gt. 0) then
                    IDIR = ID
                    SIZERATIO = RATIO
                end if
            end if

        end if
    end do

    return
end subroutine
