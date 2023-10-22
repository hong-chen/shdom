!*==SOLVE_RTE.spg  processed by SPAG 6.72Dc at 07:55 on 22 Oct 2023

subroutine SOLVE_RTE(Nstokes, Nx, Ny, Nz, Nx1, Ny1, Nang, Ml, Mm, Nlm, Nmu, &
                   & Nphi, Nleg, Nstleg, Numphase, Nphi0, Mu, Phi, Wtdo, &
                   & Maxiv, Maxic, Maxig, Maxido, Inradflag, Bcflag,   &
                   & Ipflag, Deltam, Srctype, Highorderrad, Solarflux,&
                   & Solarmu, Solaraz, Skyrad, Sfctype, Gndtemp,      &
                   & Gndalbedo, Nxsfc, Nysfc, Delxsfc, Delysfc,       &
                   & Nsfcpar, Sfcparms, Sfcgridparms, Units, Waveno,  &
                   & Wavelen, Accelflag, Solacc, Maxiter, Solcrit,    &
                   & Iter, Splitacc, Shacc, Xgrid, Ygrid, Zgrid, Temp,  &
                   & Planck, Extinct, Albedo, Legen, Iphase, Maxnbc,   &
                   & Maxbcrad, Ntoppts, Nbotpts, Bcptr, Bcrad, Npts,   &
                   & Gridpos, Ncells, Gridptr, Neighptr, Treeptr,     &
                   & Cellflags, Rshptr, Shptr, Oshptr, Work, Work1,    &
                   & Work2, Source, Delsource, Radiance, Fluxes,      &
                   & Dirflux)

!       Performs the SHDOM solution procedure.
!      Output is returned in SOURCE, RADIANCE, FLUXES, DIRFLUX.
!       The adaptive grid stucture (NCELLS,GRIDPTR,NEIGHPTR,TREEPTR,
!       CELLFLAGS,NPTS,GRIDPOS) is input and modified for new grid
!       cells and points on output.
  implicit none
!*--SOLVE_RTE25
  integer Nstokes, Nx, Ny, Nz, Nx1, Ny1, Nxsfc, Nysfc,      &
        & Nsfcpar
  integer Ml, Mm, Nlm, Nmu, Nphi, Nang, Nleg, Nstleg,       &
        & Numphase
  integer Nphi0(*), Maxiter, Iter, Bcflag, Ipflag
  integer Maxiv, Maxic, Maxig, Maxido
  integer Maxnbc, Maxbcrad, Ntoppts, Nbotpts, Bcptr(Maxnbc, 2)
  integer Npts, Ncells
  integer Rshptr(*), Shptr(*), Oshptr(*)
  integer Gridptr(8, *), Neighptr(6, *), Treeptr(2, *)
  integer*2 Cellflags(*), Iphase(*)
  integer Work1(*)
  logical Inradflag, Deltam, Accelflag, Highorderrad
  real Solarflux, Solarmu, Solaraz, Skyrad
  real Gndtemp, Gndalbedo
  real Delxsfc, Delysfc, Sfcparms(*), Sfcgridparms(*)
  real Waveno(2), Wavelen
  real Solacc, Solcrit, Splitacc, Shacc
  real Mu(*), Wtdo(*), Phi(*)
  real Temp(*), Planck(*)
  real Extinct(*), Albedo(*), Legen(Nstleg, 0:Nleg, *)
  real Xgrid(Nx1), Ygrid(Ny1), Zgrid(Nz)
  real Gridpos(3, *)
  real Bcrad(*)
  real Fluxes(2, *), Dirflux(*)
  real Source(Nstokes, *), Delsource(Nstokes, *)
  real Radiance(Nstokes, *)
  real Work(*), Work2(*)
  character Srctype*1, Units*1, Sfctype*2

  integer nphi0max, ordinateset, nbpts, nbcells
  integer ierr
  logical fixsh, splittesting, dosplit, outofmem, lambertian
  logical uniform_sfc_brdf
  logical, allocatable :: fftflag(:)
  real albmax
  real startadaptsol, endadaptsol, adaptrange, splitcrit
  real startsplitacc, cursplitacc, avgsolcrit, beta
  real deljdot, deljold, deljnew, jnorm, accelpar
  real, allocatable :: cmu1(:, :, :), cmu2(:, :, :)
  real, allocatable :: cphi1(:, :, :), cphi2(:, :, :)
  real, allocatable :: ylmsun(:, :), wtmu(:), wphisave(:)
  real, allocatable :: sfc_brdf_do(:, :, :, :, :, :)

!       Allocate local arrays
!         NSTLEG=1 for NSTOKES=1, otherwise NSTLEG=6
  allocate (cmu1(Nstleg, Nlm, Nmu), cmu2(Nstleg, Nlm, Nmu))
  allocate (cphi1(-16:16, 32, Nmu), cphi2(32, -16:16, Nmu))
  allocate (ylmsun(Nstleg, Nlm))
  allocate (fftflag(Nmu), wtmu(Nmu), wphisave(Nmu*(3*Nphi + 15)))

!       Set up some things before solution loop

!           Transfer the medium properties to the internal grid and add gas abs
  call INTERP_GRID(Npts, Nstleg, Nleg, Gridpos, Temp, Extinct, Albedo,    &
                 & Legen, Iphase)

!         If Delta-M then scale the extinction, albedo, and Legendre terms.
!         Put the Planck blackbody source in PLANCK.
  call PREPARE_PROP(Ml, Mm, Nstleg, Nleg, Npts, Deltam, Numphase, Srctype, &
                  & Units, Waveno, Wavelen, albmax, Extinct, Albedo,     &
                  & Legen, Temp, Planck, Iphase)
!       Get the maximum single scattering albedo over all processors
  call TOTAL_ALBEDO_MAX(albmax)

!           Compute the solar transmission in DIRFLUX.
!           Precompute Ylm's for solar direction
  if (Srctype .ne. 'T') then
    if (btest(Bcflag, 2) .or. btest(Bcflag, 3)) then
      call MAKE_DIRECT_PAR(1, Npts, Bcflag, Ipflag, Deltam, Ml, Nstleg, &
                         & Nleg, Solarflux, Solarmu, Solaraz, Gridpos,&
                         & Nx, Xgrid, Ny, Ygrid, Dirflux)
    else
      call MAKE_DIRECT(Npts, Bcflag, Ipflag, Deltam, Ml, Nstleg, Nleg,  &
                     & Solarflux, Solarmu, Solaraz, Gridpos, Dirflux)
    end if
    call YLMALL(.true., Solarmu, Solaraz, Ml, Mm, Nstleg, ylmsun)
  end if

!           Make the discrete ordinates (angles)
!             (set 2 is reduced gaussian, 3 is reduced double gauss)
  ordinateset = 2
  call MAKE_ANGLE_SET(Nmu, Nphi, ordinateset, nphi0max, Nphi0, Mu, Phi,   &
                    & wtmu, Wtdo, fftflag, Nang)

!           Make the Ylm transform coefficients
  call MAKE_SH_DO_COEF(Nstleg, Ml, Mm, Nlm, Nmu, Nphi0, nphi0max, Mu, Phi,  &
                     & wtmu, Wtdo, fftflag, cmu1, cmu2, cphi1, cphi2,     &
                     & wphisave)

  fixsh = .false.
  if (.not. Inradflag) then
!           Initialize the radiance on the base grid using Eddington
!             two-stream plane-parallel
    call INIT_RADIANCE(Nstokes, Nx1*Ny1, Nz, Nstleg, Nleg, Rshptr, Zgrid,&
                     & Extinct, Albedo, Legen, Temp, Numphase, Iphase,  &
                     & Srctype, Solarflux, Solarmu, Gndalbedo, Gndtemp,&
                     & Skyrad, Units, Waveno, Wavelen, Radiance)
!           Interpolate the radiance on the non-base grid points from the
!             base grid points
    nbpts = Nx1*Ny1*Nz
    nbcells = (Nz - 1)*(Nx + ibits(Bcflag, 0, 1) - ibits(Bcflag, 2, 1))      &
             & *(Ny + ibits(Bcflag, 1, 1) - ibits(Bcflag, 3, 1))
    call INTERP_RADIANCE(Nstokes, nbpts, Npts, Rshptr, Radiance,       &
                       & nbcells, Ncells, Treeptr, Gridptr, Gridpos)
!           Initialize the source function from the radiance field
    call COMPUTE_SOURCE(Nstokes, Ml, Mm, Nlm, Nstleg, Nleg, Numphase,    &
                      & Npts, fixsh, Srctype, Solarmu, ylmsun, Albedo,  &
                      & Legen, Iphase, Planck, Dirflux, Shacc, Maxiv,   &
                      & Rshptr, Radiance, Shptr, Source, Oshptr,       &
                      & Delsource, .true., Accelflag, deljdot, deljold,&
                      & deljnew, jnorm)
  end if

  if (Accelflag) then
    Oshptr(1:Npts + 1) = Shptr(1:Npts + 1)
    Delsource(:, 1:Oshptr(Npts + 1)) = 0.0
  end if
  Solcrit = 1.0

  lambertian = Sfctype(2:2) .eq. 'L'
!         Make the pointers to the grid points on the top and bottom boundaries
  call BOUNDARY_PNTS(Npts, Nang, lambertian, Maxnbc, Maxbcrad, Zgrid(1), &
                   & Zgrid(Nz), Gridpos, Ntoppts, Nbotpts, Bcptr)
  if (Sfctype(1:1) .eq. 'V') then
!             If this is a variable surface then bilinearly interpolate
!               the surface temperature and reflection parameters to the
!               bottom grid points.
    call SURFACE_PARM_INTERP(Nbotpts, Bcptr(1, 2), Gridpos, Srctype,   &
                           & Waveno, Wavelen, Units, Nxsfc, Nysfc,     &
                           & Delxsfc, Delysfc, Nsfcpar, Sfcparms,     &
                           & Sfcgridparms)
    call CHECK_UNIFORM_SFC(Nxsfc, Nysfc, Nsfcpar, Sfcparms,           &
                         & uniform_sfc_brdf)
    if (uniform_sfc_brdf) then
      if (Nstokes*(Nmu/2)                                        &
         & *nphi0max .gt. sqrt(1.0*huge(Nmu)/kind(Nmu))) then
        ierr = 1
      else
        allocate (sfc_brdf_do(Nstokes, Nmu/2, nphi0max, Nstokes,    &
                & Nmu/2, nphi0max), STAT=ierr)
      end if
      if (ierr .eq. 0) then
        call MAKE_SFC_BRDF_DO_MATRIX(Sfctype, Wavelen, Nxsfc, Nysfc,&
         & Nsfcpar, Sfcparms, Nstokes, Nmu, nphi0max, Nphi0, Mu, Phi,   &
         & Wtdo, sfc_brdf_do)
      else
        print '(A,A)',                                          &
             &'Not enough memory to allocate precomputed',      &
             &' surface BRDF discrete ordinate matrix: '
        print *, '  Proceeding more slowly.'
        uniform_sfc_brdf = .false.
      end if
    end if
  end if

  Iter = 0
  if (Maxiter .le. Iter) return
!         Starting values for the adaptive cell splitting controlling method.
!         The splitting accuracy used at a particular iteration is CURSPLITACC,
!         and it is started out high an lowered to the desired SPLITACC.
!         The cell splitting starts when the solution criterion (SOLCRIT)
!         goes below STARTADAPTSOL and ends when the solution criterion
!         goes below ENDADAPTSOL, though the current splitting accuracy
!         is designed to get to the final SPLITACC over a smaller range
!         (ADAPTRANGE).  The currect splitting accuracy is computed from
!         a geometric average of the solution criterion (AVGSOLCRIT).
!         The SPLIT_GRID routine is called and the splitting criterion returned
!         as long as SPLITTESTING is true.
  splittesting = .true.
  outofmem = .false.
  endadaptsol = 0.001
  startadaptsol = 0.1
  adaptrange = startadaptsol/(3.0*endadaptsol)
  cursplitacc = Splitacc*adaptrange
  startsplitacc = cursplitacc
  avgsolcrit = Solcrit
  splitcrit = 0.0

  write (6, *) '! Iter Log(Sol)  SplitCrit  Npoints  Nsh(avg)'

!         Main solution loop
!           Iterate until the iterations exceed the limit, or the desired
!           solution criterion is reached and the current splitting accuracy
!           is down to the desired level.
  do while (Iter .lt. Maxiter .and.                                  &
           & (Solcrit .gt. Solacc .or. (splitcrit .gt. Splitacc .and.   &
           & cursplitacc .gt. Splitacc .and. .not. outofmem)))
    Iter = Iter + 1

!           Don't do the adaptive grid cell splitting stuff unless wanted
    if (Splitacc .gt. 0.0) then
!             If in the right solution range then try to split some cells, and
!               interpolate the medium and source functions to the new points.
!             The current splitting accuracy is started high and then lowered.
      avgsolcrit = sqrt(avgsolcrit*Solcrit)
      dosplit = Solcrit .le. startadaptsol .and.                    &
              & (Solcrit .gt. endadaptsol .or.                      &
              & cursplitacc .gt. Splitacc) .and. .not. outofmem
!           Make sure all processors are going to split cells if any want to
      call UNIFY_SPLITTING(dosplit, startsplitacc)
      beta = log(startsplitacc/Splitacc)/log(adaptrange)
      cursplitacc = min(cursplitacc,                              &
                  & max(Splitacc, Splitacc*(avgsolcrit/            &
                  & (3.0*endadaptsol))**beta))
      if (Solcrit .le. endadaptsol) cursplitacc = Splitacc
      if (splittesting) then
        call SPLIT_GRID(Maxig, Maxic, Maxiv, Maxido, nphi0max,       &
                      & dosplit, outofmem, cursplitacc, splitcrit,  &
                      & Npts, Ncells, Gridptr, Neighptr, Treeptr,    &
                      & Cellflags, Gridpos, Nx, Xgrid, Ny, Ygrid,     &
                      & Nstokes, Ml, Mm, Nlm, Nstleg, Nleg, Numphase,  &
                      & Deltam, Bcflag, Ipflag, Accelflag, Extinct,  &
                      & Albedo, Legen, Iphase, Rshptr, Radiance,     &
                      & Shptr, Source, Oshptr, Temp, Planck, Dirflux, &
                      & Srctype, Solarflux, Solarmu, Solaraz, ylmsun,&
                      & Units, Waveno, Wavelen, Work, Work2)
        if (Solcrit .gt. startadaptsol) startsplitacc = splitcrit
      end if
!           Find the maximum splitting criterion over all processors
      call TOTAL_SPLITCRIT_MAX(splitcrit)
      if (Solcrit .le. endadaptsol) splittesting = .false.
    end if

!            Find the SH truncation to use for the radiance (RSHPTR)
    call RADIANCE_TRUNCATION(Highorderrad, Nstokes, Ml, Mm, Nstleg,    &
                           & Nleg, Numphase, Npts, Albedo, Legen,      &
                           & Iphase, Shptr, Radiance, Maxiv + Maxig,    &
                           & fixsh, Shacc, Rshptr)

!           Integrate the source function along discrete ordinates
!             to compute the radiance field.  The input source function
!             and output radiance are in spherical harmonic expansions,
!             and the SH/DO transforms are done during the loop over
!             discrete ordinates.
    call PATH_INTEGRATION(Nx, Ny, Nz, Nstokes, Nstleg, Ml, Mm, Nlm, Nmu,   &
                        & nphi0max, Nphi0, Nang, Bcflag, Ipflag, Npts,  &
                        & Ncells, Xgrid, Ygrid, Zgrid, Gridpos, Gridptr,&
                        & Neighptr, Treeptr, Cellflags, Srctype,      &
                        & Solarmu, Solaraz, Skyrad, Waveno, Wavelen,   &
                        & Units, Sfctype, Gndtemp, Gndalbedo, Nxsfc,   &
                        & Nysfc, Delxsfc, Delysfc, Nsfcpar, Sfcparms,  &
                        & Sfcgridparms, uniform_sfc_brdf,           &
                        & sfc_brdf_do, Maxnbc, Maxbcrad, Ntoppts,     &
                        & Nbotpts, Bcptr, Bcrad, Mu, Phi, Wtdo, cmu1,    &
                        & cphi1, cmu2, cphi2, wphisave, fftflag,       &
                        & Dirflux, Fluxes, Extinct, Shptr, Source,     &
                        & Rshptr, Radiance, Work, Work1, Work2)

!            Compute the source function from the radiance field,
!              do the adaptive spherical harmonics truncation, compute
!              the solution criterion, and dot products for acceleration.
    if (Solcrit .lt. endadaptsol .or. Iter .gt. 30) fixsh = .true.
    call COMPUTE_SOURCE(Nstokes, Ml, Mm, Nlm, Nstleg, Nleg, Numphase,    &
                      & Npts, fixsh, Srctype, Solarmu, ylmsun, Albedo,  &
                      & Legen, Iphase, Planck, Dirflux, Shacc, Maxiv,   &
                      & Rshptr, Radiance, Shptr, Source, Oshptr,       &
                      & Delsource, .false., Accelflag, deljdot,       &
                      & deljold, deljnew, jnorm)

!           Calculate the acceleration parameter and solution criterion
!            from all the processors
    call CALC_ACCEL_SOLCRIT(Accelflag, deljdot, deljold, deljnew,     &
                          & jnorm, accelpar, Solcrit)

!           Accelerate the convergence of the source function vector
    call ACCELERATE_SOLUTION(accelpar, Npts, Shptr, Oshptr, Nstokes,   &
                           & Source, Delsource)

!           If it is a not scattering medium then do only one iteration
    if (albmax .lt. Solacc) Solcrit = Solacc

!           Print out the log solution criterion, number of points,
!             and average SH truncation.
    write (6, '(2X,I4,F8.3,1X,E10.3,1X,I8,1X,F8.2,1X,F6.3)') Iter, &
         & log10(max(Solcrit, 1.0e-20)), splitcrit, Npts,        &
         & FLOAT(Shptr(Npts + 1))/Npts, FLOAT(Shptr(Npts + 1))        &
         & /(Npts*Nlm)
  end do

  write (6, '(1X,A,I6,A,F9.6)') '! Iterations: ', Iter,            &
                              &'     Final Criterion: ', Solcrit

  if (Sfctype(1:1) .eq. 'V' .and. uniform_sfc_brdf)                 &
     & deallocate (sfc_brdf_do)
  deallocate (ylmsun, fftflag, wtmu, wphisave)
  deallocate (cmu1, cmu2, cphi1, cphi2)

!          Comment in for GLE graphics output of cell structure (see routine)
!      CALL VISUALIZE_CELLS ('cellfinal.gle', YGRID(1), 1, IPFLAG,
!     .           NX, XGRID(NX+1)-XGRID(1), ZGRID(NZ)-ZGRID(1),
!     .           NCELLS, GRIDPTR, NEIGHPTR, TREEPTR,
!     .           CELLFLAGS, GRIDPOS)

!          Comment in for output of cell splitting criterion for every cell
!      CALL OUTPUT_CELL_SPLIT ('cellsplit.dat',GRIDPTR,
!     .                        GRIDPOS, EXTINCT, SHPTR, SOURCE, NCELLS)

end
!*==CALC_SOURCE_PNT.spg  processed by SPAG 6.72Dc at 07:55 on 22 Oct 2023

subroutine CALC_SOURCE_PNT(Nstokes, Nlm, Nstleg, Nleg, Lofj, Srctype,  &
                         & Flux0, Ylmsun, Planck, Albedo, Iph, Legen,  &
                         & Nr, Radiance, Sourcet)
!      Calculates the polarized source function in real generalized
!      spherical harmonics at one grid point from the radiance in
!      real generalized spherical harmonics and thermal and/or solar
!      source information.
  implicit none
!*--CALC_SOURCE_PNT342
  integer Nstokes, Nstleg, Nleg, Nlm, Lofj(Nlm), Iph, Nr
  character Srctype*1
  real Flux0, Ylmsun(Nstleg, Nlm), Planck, Albedo
  real Legen(Nstleg, 0:Nleg, *)
  real Radiance(Nstokes, Nr)
  real Sourcet(Nstokes, Nlm)
  real C
  parameter(C=3.544907703)
  integer j

  Sourcet(:, :) = 0.0
!       Calculate the generalized spherical harmonics source function
!       from the radiance vector in genSH.  Do different parts of the
!       matrix multiply depending on NSTOKES.
  do j = 1, Nr
    Sourcet(1, j) = Sourcet(1, j) + Legen(1, Lofj(j), Iph)             &
                   & *Radiance(1, j)
  end do
  if (Nstokes .gt. 1) then
    do j = 1, Nr
      Sourcet(1, j) = Sourcet(1, j) + Legen(5, Lofj(j), Iph)          &
                     & *Radiance(2, j)
    end do
    do j = 5, Nr
      Sourcet(2, j) = Sourcet(2, j) + Legen(5, Lofj(j), Iph)          &
                     & *Radiance(1, j) + Legen(2, Lofj(j), Iph)        &
                     & *Radiance(2, j)
      Sourcet(3, j) = Sourcet(3, j) + Legen(3, Lofj(j), Iph)          &
                     & *Radiance(3, j)
    end do
  end if
  if (Nstokes .eq. 4) then
    do j = 1, Nr
      Sourcet(3, j) = Sourcet(3, j) + Legen(6, Lofj(j), Iph)          &
                     & *Radiance(4, j)
      Sourcet(4, j) = Sourcet(4, j) - Legen(6, Lofj(j), Iph)          &
                     & *Radiance(3, j) + Legen(4, Lofj(j), Iph)        &
                     & *Radiance(4, j)
    end do
  end if
!       Do the final multiplication by the single scattering albedo
  Sourcet(:, :) = Albedo*Sourcet(:, :)

!       Calculate the generalized spherical harmonics source function
!       from the solar direct beam
  if (Srctype .eq. 'S' .or. Srctype .eq. 'B') then
    do j = 1, Nlm
      Sourcet(1, j) = Sourcet(1, j)                                 &
                   & + Flux0*Albedo*Legen(1, Lofj(j), Iph)          &
                   & *Ylmsun(1, j)
    end do
    if (Nstokes .gt. 1) then
      do j = 5, Nlm
        Sourcet(2, j) = Sourcet(2, j)                              &
                     & + Flux0*Albedo*Legen(5, Lofj(j), Iph)       &
                     & *Ylmsun(1, j)
      end do
    end if
  end if
!       Add in the isotropic unpolarized thermal emission
  if (Srctype .eq. 'T' .or. Srctype .eq. 'B') Sourcet(1, 1)            &
     & = Sourcet(1, 1) + C*Planck

end
