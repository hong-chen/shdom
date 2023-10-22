!*==SOLVE_RTE.spg  processed by SPAG 6.72Dc at 07:55 on 22 Oct 2023

      SUBROUTINE SOLVE_RTE(Nstokes,Nx,Ny,Nz,Nx1,Ny1,Nang,Ml,Mm,Nlm,Nmu, &
                         & Nphi,Nleg,Nstleg,Numphase,Nphi0,Mu,Phi,Wtdo, &
                         & Maxiv,Maxic,Maxig,Maxido,Inradflag,Bcflag,   &
                         & Ipflag,Deltam,Srctype,Highorderrad,Solarflux,&
                         & Solarmu,Solaraz,Skyrad,Sfctype,Gndtemp,      &
                         & Gndalbedo,Nxsfc,Nysfc,Delxsfc,Delysfc,       &
                         & Nsfcpar,Sfcparms,Sfcgridparms,Units,Waveno,  &
                         & Wavelen,Accelflag,Solacc,Maxiter,Solcrit,    &
                         & Iter,Splitacc,Shacc,Xgrid,Ygrid,Zgrid,Temp,  &
                         & Planck,Extinct,Albedo,Legen,Iphase,Maxnbc,   &
                         & Maxbcrad,Ntoppts,Nbotpts,Bcptr,Bcrad,Npts,   &
                         & Gridpos,Ncells,Gridptr,Neighptr,Treeptr,     &
                         & Cellflags,Rshptr,Shptr,Oshptr,Work,Work1,    &
                         & Work2,Source,Delsource,Radiance,Fluxes,      &
                         & Dirflux)
 
!       Performs the SHDOM solution procedure.
!      Output is returned in SOURCE, RADIANCE, FLUXES, DIRFLUX.
!       The adaptive grid stucture (NCELLS,GRIDPTR,NEIGHPTR,TREEPTR,
!       CELLFLAGS,NPTS,GRIDPOS) is input and modified for new grid
!       cells and points on output.
      IMPLICIT NONE
!*--SOLVE_RTE25
      INTEGER Nstokes , Nx , Ny , Nz , Nx1 , Ny1 , Nxsfc , Nysfc ,      &
            & Nsfcpar
      INTEGER Ml , Mm , Nlm , Nmu , Nphi , Nang , Nleg , Nstleg ,       &
            & Numphase
      INTEGER Nphi0(*) , Maxiter , Iter , Bcflag , Ipflag
      INTEGER Maxiv , Maxic , Maxig , Maxido
      INTEGER Maxnbc , Maxbcrad , Ntoppts , Nbotpts , Bcptr(Maxnbc,2)
      INTEGER Npts , Ncells
      INTEGER Rshptr(*) , Shptr(*) , Oshptr(*)
      INTEGER Gridptr(8,*) , Neighptr(6,*) , Treeptr(2,*)
      INTEGER*2 Cellflags(*) , Iphase(*)
      INTEGER Work1(*)
      LOGICAL Inradflag , Deltam , Accelflag , Highorderrad
      REAL Solarflux , Solarmu , Solaraz , Skyrad
      REAL Gndtemp , Gndalbedo
      REAL Delxsfc , Delysfc , Sfcparms(*) , Sfcgridparms(*)
      REAL Waveno(2) , Wavelen
      REAL Solacc , Solcrit , Splitacc , Shacc
      REAL Mu(*) , Wtdo(*) , Phi(*)
      REAL Temp(*) , Planck(*)
      REAL Extinct(*) , Albedo(*) , Legen(Nstleg,0:Nleg,*)
      REAL Xgrid(Nx1) , Ygrid(Ny1) , Zgrid(Nz)
      REAL Gridpos(3,*)
      REAL Bcrad(*)
      REAL Fluxes(2,*) , Dirflux(*)
      REAL Source(Nstokes,*) , Delsource(Nstokes,*)
      REAL Radiance(Nstokes,*)
      REAL Work(*) , Work2(*)
      CHARACTER Srctype*1 , Units*1 , Sfctype*2
 
      INTEGER nphi0max , ordinateset , nbpts , nbcells
      INTEGER ierr
      LOGICAL fixsh , splittesting , dosplit , outofmem , lambertian
      LOGICAL uniform_sfc_brdf
      LOGICAL , ALLOCATABLE :: fftflag(:)
      REAL albmax
      REAL startadaptsol , endadaptsol , adaptrange , splitcrit
      REAL startsplitacc , cursplitacc , avgsolcrit , beta
      REAL deljdot , deljold , deljnew , jnorm , accelpar
      REAL , ALLOCATABLE :: cmu1(:,:,:) , cmu2(:,:,:)
      REAL , ALLOCATABLE :: cphi1(:,:,:) , cphi2(:,:,:)
      REAL , ALLOCATABLE :: ylmsun(:,:) , wtmu(:) , wphisave(:)
      REAL , ALLOCATABLE :: sfc_brdf_do(:,:,:,:,:,:)
 
!       Allocate local arrays
!         NSTLEG=1 for NSTOKES=1, otherwise NSTLEG=6
      ALLOCATE (cmu1(Nstleg,Nlm,Nmu),cmu2(Nstleg,Nlm,Nmu))
      ALLOCATE (cphi1(-16:16,32,Nmu),cphi2(32,-16:16,Nmu))
      ALLOCATE (ylmsun(Nstleg,Nlm))
      ALLOCATE (fftflag(Nmu),wtmu(Nmu),wphisave(Nmu*(3*Nphi+15)))
 
 
!       Set up some things before solution loop
 
!           Transfer the medium properties to the internal grid and add gas abs
      CALL INTERP_GRID(Npts,Nstleg,Nleg,Gridpos,Temp,Extinct,Albedo,    &
                     & Legen,Iphase)
 
!         If Delta-M then scale the extinction, albedo, and Legendre terms.
!         Put the Planck blackbody source in PLANCK.
      CALL PREPARE_PROP(Ml,Mm,Nstleg,Nleg,Npts,Deltam,Numphase,Srctype, &
                      & Units,Waveno,Wavelen,albmax,Extinct,Albedo,     &
                      & Legen,Temp,Planck,Iphase)
!       Get the maximum single scattering albedo over all processors
      CALL TOTAL_ALBEDO_MAX(albmax)
 
!           Compute the solar transmission in DIRFLUX.
!           Precompute Ylm's for solar direction
      IF ( Srctype.NE.'T' ) THEN
         IF ( BTEST(Bcflag,2) .OR. BTEST(Bcflag,3) ) THEN
            CALL MAKE_DIRECT_PAR(1,Npts,Bcflag,Ipflag,Deltam,Ml,Nstleg, &
                               & Nleg,Solarflux,Solarmu,Solaraz,Gridpos,&
                               & Nx,Xgrid,Ny,Ygrid,Dirflux)
         ELSE
            CALL MAKE_DIRECT(Npts,Bcflag,Ipflag,Deltam,Ml,Nstleg,Nleg,  &
                           & Solarflux,Solarmu,Solaraz,Gridpos,Dirflux)
         ENDIF
         CALL YLMALL(.TRUE.,Solarmu,Solaraz,Ml,Mm,Nstleg,ylmsun)
      ENDIF
 
 
!           Make the discrete ordinates (angles)
!             (set 2 is reduced gaussian, 3 is reduced double gauss)
      ordinateset = 2
      CALL MAKE_ANGLE_SET(Nmu,Nphi,ordinateset,nphi0max,Nphi0,Mu,Phi,   &
                        & wtmu,Wtdo,fftflag,Nang)
 
 
!           Make the Ylm transform coefficients
      CALL MAKE_SH_DO_COEF(Nstleg,Ml,Mm,Nlm,Nmu,Nphi0,nphi0max,Mu,Phi,  &
                         & wtmu,Wtdo,fftflag,cmu1,cmu2,cphi1,cphi2,     &
                         & wphisave)
 
      fixsh = .FALSE.
      IF ( .NOT.Inradflag ) THEN
!           Initialize the radiance on the base grid using Eddington
!             two-stream plane-parallel
         CALL INIT_RADIANCE(Nstokes,Nx1*Ny1,Nz,Nstleg,Nleg,Rshptr,Zgrid,&
                          & Extinct,Albedo,Legen,Temp,Numphase,Iphase,  &
                          & Srctype,Solarflux,Solarmu,Gndalbedo,Gndtemp,&
                          & Skyrad,Units,Waveno,Wavelen,Radiance)
!           Interpolate the radiance on the non-base grid points from the
!             base grid points
         nbpts = Nx1*Ny1*Nz
         nbcells = (Nz-1)*(Nx+IBITS(Bcflag,0,1)-IBITS(Bcflag,2,1))      &
                 & *(Ny+IBITS(Bcflag,1,1)-IBITS(Bcflag,3,1))
         CALL INTERP_RADIANCE(Nstokes,nbpts,Npts,Rshptr,Radiance,       &
                            & nbcells,Ncells,Treeptr,Gridptr,Gridpos)
!           Initialize the source function from the radiance field
         CALL COMPUTE_SOURCE(Nstokes,Ml,Mm,Nlm,Nstleg,Nleg,Numphase,    &
                           & Npts,fixsh,Srctype,Solarmu,ylmsun,Albedo,  &
                           & Legen,Iphase,Planck,Dirflux,Shacc,Maxiv,   &
                           & Rshptr,Radiance,Shptr,Source,Oshptr,       &
                           & Delsource,.TRUE.,Accelflag,deljdot,deljold,&
                           & deljnew,jnorm)
      ENDIF
 
      IF ( Accelflag ) THEN
         Oshptr(1:Npts+1) = Shptr(1:Npts+1)
         Delsource(:,1:Oshptr(Npts+1)) = 0.0
      ENDIF
      Solcrit = 1.0
 
 
      lambertian = Sfctype(2:2).EQ.'L'
!         Make the pointers to the grid points on the top and bottom boundaries
      CALL BOUNDARY_PNTS(Npts,Nang,lambertian,Maxnbc,Maxbcrad,Zgrid(1), &
                       & Zgrid(Nz),Gridpos,Ntoppts,Nbotpts,Bcptr)
      IF ( Sfctype(1:1).EQ.'V' ) THEN
!             If this is a variable surface then bilinearly interpolate
!               the surface temperature and reflection parameters to the
!               bottom grid points.
         CALL SURFACE_PARM_INTERP(Nbotpts,Bcptr(1,2),Gridpos,Srctype,   &
                                & Waveno,Wavelen,Units,Nxsfc,Nysfc,     &
                                & Delxsfc,Delysfc,Nsfcpar,Sfcparms,     &
                                & Sfcgridparms)
         CALL CHECK_UNIFORM_SFC(Nxsfc,Nysfc,Nsfcpar,Sfcparms,           &
                              & uniform_sfc_brdf)
         IF ( uniform_sfc_brdf ) THEN
            IF ( Nstokes*(Nmu/2)                                        &
               & *nphi0max.GT.SQRT(1.0*HUGE(Nmu)/KIND(Nmu)) ) THEN
               ierr = 1
            ELSE
               ALLOCATE (sfc_brdf_do(Nstokes,Nmu/2,nphi0max,Nstokes,    &
                       & Nmu/2,nphi0max),STAT=ierr)
            ENDIF
            IF ( ierr.EQ.0 ) THEN
               CALL MAKE_SFC_BRDF_DO_MATRIX(Sfctype,Wavelen,Nxsfc,Nysfc,&
                & Nsfcpar,Sfcparms,Nstokes,Nmu,nphi0max,Nphi0,Mu,Phi,   &
                & Wtdo,sfc_brdf_do)
            ELSE
               PRINT '(A,A)' ,                                          &
                    &'Not enough memory to allocate precomputed' ,      &
                    &' surface BRDF discrete ordinate matrix: '
               PRINT * , '  Proceeding more slowly.'
               uniform_sfc_brdf = .FALSE.
            ENDIF
         ENDIF
      ENDIF
 
      Iter = 0
      IF ( Maxiter.LE.Iter ) RETURN
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
      splittesting = .TRUE.
      outofmem = .FALSE.
      endadaptsol = 0.001
      startadaptsol = 0.1
      adaptrange = startadaptsol/(3.0*endadaptsol)
      cursplitacc = Splitacc*adaptrange
      startsplitacc = cursplitacc
      avgsolcrit = Solcrit
      splitcrit = 0.0
 
      WRITE (6,*) '! Iter Log(Sol)  SplitCrit  Npoints  Nsh(avg)'
 
!         Main solution loop
!           Iterate until the iterations exceed the limit, or the desired
!           solution criterion is reached and the current splitting accuracy
!           is down to the desired level.
      DO WHILE ( Iter.LT.Maxiter .AND.                                  &
               & (Solcrit.GT.Solacc .OR. (splitcrit.GT.Splitacc .AND.   &
               & cursplitacc.GT.Splitacc .AND. .NOT.outofmem)) )
         Iter = Iter + 1
 
!           Don't do the adaptive grid cell splitting stuff unless wanted
         IF ( Splitacc.GT.0.0 ) THEN
!             If in the right solution range then try to split some cells, and
!               interpolate the medium and source functions to the new points.
!             The current splitting accuracy is started high and then lowered.
            avgsolcrit = SQRT(avgsolcrit*Solcrit)
            dosplit = Solcrit.LE.startadaptsol .AND.                    &
                    & (Solcrit.GT.endadaptsol .OR.                      &
                    & cursplitacc.GT.Splitacc) .AND. .NOT.outofmem
!           Make sure all processors are going to split cells if any want to
            CALL UNIFY_SPLITTING(dosplit,startsplitacc)
            beta = LOG(startsplitacc/Splitacc)/LOG(adaptrange)
            cursplitacc = MIN(cursplitacc,                              &
                        & MAX(Splitacc,Splitacc*(avgsolcrit/            &
                        & (3.0*endadaptsol))**beta))
            IF ( Solcrit.LE.endadaptsol ) cursplitacc = Splitacc
            IF ( splittesting ) THEN
               CALL SPLIT_GRID(Maxig,Maxic,Maxiv,Maxido,nphi0max,       &
                             & dosplit,outofmem,cursplitacc,splitcrit,  &
                             & Npts,Ncells,Gridptr,Neighptr,Treeptr,    &
                             & Cellflags,Gridpos,Nx,Xgrid,Ny,Ygrid,     &
                             & Nstokes,Ml,Mm,Nlm,Nstleg,Nleg,Numphase,  &
                             & Deltam,Bcflag,Ipflag,Accelflag,Extinct,  &
                             & Albedo,Legen,Iphase,Rshptr,Radiance,     &
                             & Shptr,Source,Oshptr,Temp,Planck,Dirflux, &
                             & Srctype,Solarflux,Solarmu,Solaraz,ylmsun,&
                             & Units,Waveno,Wavelen,Work,Work2)
               IF ( Solcrit.GT.startadaptsol ) startsplitacc = splitcrit
            ENDIF
!           Find the maximum splitting criterion over all processors
            CALL TOTAL_SPLITCRIT_MAX(splitcrit)
            IF ( Solcrit.LE.endadaptsol ) splittesting = .FALSE.
         ENDIF
 
!            Find the SH truncation to use for the radiance (RSHPTR)
         CALL RADIANCE_TRUNCATION(Highorderrad,Nstokes,Ml,Mm,Nstleg,    &
                                & Nleg,Numphase,Npts,Albedo,Legen,      &
                                & Iphase,Shptr,Radiance,Maxiv+Maxig,    &
                                & fixsh,Shacc,Rshptr)
 
!           Integrate the source function along discrete ordinates
!             to compute the radiance field.  The input source function
!             and output radiance are in spherical harmonic expansions,
!             and the SH/DO transforms are done during the loop over
!             discrete ordinates.
         CALL PATH_INTEGRATION(Nx,Ny,Nz,Nstokes,Nstleg,Ml,Mm,Nlm,Nmu,   &
                             & nphi0max,Nphi0,Nang,Bcflag,Ipflag,Npts,  &
                             & Ncells,Xgrid,Ygrid,Zgrid,Gridpos,Gridptr,&
                             & Neighptr,Treeptr,Cellflags,Srctype,      &
                             & Solarmu,Solaraz,Skyrad,Waveno,Wavelen,   &
                             & Units,Sfctype,Gndtemp,Gndalbedo,Nxsfc,   &
                             & Nysfc,Delxsfc,Delysfc,Nsfcpar,Sfcparms,  &
                             & Sfcgridparms,uniform_sfc_brdf,           &
                             & sfc_brdf_do,Maxnbc,Maxbcrad,Ntoppts,     &
                             & Nbotpts,Bcptr,Bcrad,Mu,Phi,Wtdo,cmu1,    &
                             & cphi1,cmu2,cphi2,wphisave,fftflag,       &
                             & Dirflux,Fluxes,Extinct,Shptr,Source,     &
                             & Rshptr,Radiance,Work,Work1,Work2)
 
!            Compute the source function from the radiance field,
!              do the adaptive spherical harmonics truncation, compute
!              the solution criterion, and dot products for acceleration.
         IF ( Solcrit.LT.endadaptsol .OR. Iter.GT.30 ) fixsh = .TRUE.
         CALL COMPUTE_SOURCE(Nstokes,Ml,Mm,Nlm,Nstleg,Nleg,Numphase,    &
                           & Npts,fixsh,Srctype,Solarmu,ylmsun,Albedo,  &
                           & Legen,Iphase,Planck,Dirflux,Shacc,Maxiv,   &
                           & Rshptr,Radiance,Shptr,Source,Oshptr,       &
                           & Delsource,.FALSE.,Accelflag,deljdot,       &
                           & deljold,deljnew,jnorm)
 
!           Calculate the acceleration parameter and solution criterion
!            from all the processors
         CALL CALC_ACCEL_SOLCRIT(Accelflag,deljdot,deljold,deljnew,     &
                               & jnorm,accelpar,Solcrit)
 
!           Accelerate the convergence of the source function vector
         CALL ACCELERATE_SOLUTION(accelpar,Npts,Shptr,Oshptr,Nstokes,   &
                                & Source,Delsource)
 
!           If it is a not scattering medium then do only one iteration
         IF ( albmax.LT.Solacc ) Solcrit = Solacc
 
!           Print out the log solution criterion, number of points,
!             and average SH truncation.
         WRITE (6,'(2X,I4,F8.3,1X,E10.3,1X,I8,1X,F8.2,1X,F6.3)') Iter , &
              & LOG10(MAX(Solcrit,1.0E-20)) , splitcrit , Npts ,        &
              & FLOAT(Shptr(Npts+1))/Npts , FLOAT(Shptr(Npts+1))        &
              & /(Npts*Nlm)
      ENDDO
 
      WRITE (6,'(1X,A,I6,A,F9.6)') '! Iterations: ' , Iter ,            &
                                  &'     Final Criterion: ' , Solcrit
 
      IF ( Sfctype(1:1).EQ.'V' .AND. uniform_sfc_brdf )                 &
         & DEALLOCATE (sfc_brdf_do)
      DEALLOCATE (ylmsun,fftflag,wtmu,wphisave)
      DEALLOCATE (cmu1,cmu2,cphi1,cphi2)
 
!          Comment in for GLE graphics output of cell structure (see routine)
!      CALL VISUALIZE_CELLS ('cellfinal.gle', YGRID(1), 1, IPFLAG,
!     .           NX, XGRID(NX+1)-XGRID(1), ZGRID(NZ)-ZGRID(1),
!     .           NCELLS, GRIDPTR, NEIGHPTR, TREEPTR,
!     .           CELLFLAGS, GRIDPOS)
 
!          Comment in for output of cell splitting criterion for every cell
!      CALL OUTPUT_CELL_SPLIT ('cellsplit.dat',GRIDPTR,
!     .                        GRIDPOS, EXTINCT, SHPTR, SOURCE, NCELLS)
 
      END
!*==CALC_SOURCE_PNT.spg  processed by SPAG 6.72Dc at 07:55 on 22 Oct 2023
 
 
 
 
      SUBROUTINE CALC_SOURCE_PNT(Nstokes,Nlm,Nstleg,Nleg,Lofj,Srctype,  &
                               & Flux0,Ylmsun,Planck,Albedo,Iph,Legen,  &
                               & Nr,Radiance,Sourcet)
!      Calculates the polarized source function in real generalized
!      spherical harmonics at one grid point from the radiance in
!      real generalized spherical harmonics and thermal and/or solar
!      source information.
      IMPLICIT NONE
!*--CALC_SOURCE_PNT342
      INTEGER Nstokes , Nstleg , Nleg , Nlm , Lofj(Nlm) , Iph , Nr
      CHARACTER Srctype*1
      REAL Flux0 , Ylmsun(Nstleg,Nlm) , Planck , Albedo
      REAL Legen(Nstleg,0:Nleg,*)
      REAL Radiance(Nstokes,Nr)
      REAL Sourcet(Nstokes,Nlm)
      REAL C
      PARAMETER (C=3.544907703)
      INTEGER j
 
      Sourcet(:,:) = 0.0
!       Calculate the generalized spherical harmonics source function
!       from the radiance vector in genSH.  Do different parts of the
!       matrix multiply depending on NSTOKES.
      DO j = 1 , Nr
         Sourcet(1,j) = Sourcet(1,j) + Legen(1,Lofj(j),Iph)             &
                      & *Radiance(1,j)
      ENDDO
      IF ( Nstokes.GT.1 ) THEN
         DO j = 1 , Nr
            Sourcet(1,j) = Sourcet(1,j) + Legen(5,Lofj(j),Iph)          &
                         & *Radiance(2,j)
         ENDDO
         DO j = 5 , Nr
            Sourcet(2,j) = Sourcet(2,j) + Legen(5,Lofj(j),Iph)          &
                         & *Radiance(1,j) + Legen(2,Lofj(j),Iph)        &
                         & *Radiance(2,j)
            Sourcet(3,j) = Sourcet(3,j) + Legen(3,Lofj(j),Iph)          &
                         & *Radiance(3,j)
         ENDDO
      ENDIF
      IF ( Nstokes.EQ.4 ) THEN
         DO j = 1 , Nr
            Sourcet(3,j) = Sourcet(3,j) + Legen(6,Lofj(j),Iph)          &
                         & *Radiance(4,j)
            Sourcet(4,j) = Sourcet(4,j) - Legen(6,Lofj(j),Iph)          &
                         & *Radiance(3,j) + Legen(4,Lofj(j),Iph)        &
                         & *Radiance(4,j)
         ENDDO
      ENDIF
!       Do the final multiplication by the single scattering albedo
      Sourcet(:,:) = Albedo*Sourcet(:,:)
 
!       Calculate the generalized spherical harmonics source function
!       from the solar direct beam
      IF ( Srctype.EQ.'S' .OR. Srctype.EQ.'B' ) THEN
         DO j = 1 , Nlm
            Sourcet(1,j) = Sourcet(1,j)                                 &
                         & + Flux0*Albedo*Legen(1,Lofj(j),Iph)          &
                         & *Ylmsun(1,j)
         ENDDO
         IF ( Nstokes.GT.1 ) THEN
            DO j = 5 , Nlm
               Sourcet(2,j) = Sourcet(2,j)                              &
                            & + Flux0*Albedo*Legen(5,Lofj(j),Iph)       &
                            & *Ylmsun(1,j)
            ENDDO
         ENDIF
      ENDIF
!       Add in the isotropic unpolarized thermal emission
      IF ( Srctype.EQ.'T' .OR. Srctype.EQ.'B' ) Sourcet(1,1)            &
         & = Sourcet(1,1) + C*Planck
 
      END
