MODULE SHDOMPPDA_FWD_TL_ADJ

! This module contains the single column forward, tangent linear, and  
! adjoint radiative transfer using SHDOMPPDA.  Single scattering 
! properties as a function of mass mean radius are input from an ascii 
! file table for each hydrometeor species and wavelength.  Molecular
! absorption is input through generic subroutine calls (FORWARD_ABSORPTION,  
! TANLIN_ABSORPTION, ADJOINT_ABSORPTION). 
! 
! SHDOMPPDA: Spherical harmonic discrete ordinate radiative transfer method 
! for plane-parallel atmospheres modified for data assimilation. 
! 
! 
!   Frank Evans     January 2006

  ! Dummy molecular absorption routines (replace with real ones you make)
USE molec_absorption_interface

IMPLICIT NONE

PRIVATE  ! All variables are private to this module
PUBLIC :: READ_SCATTERING_TABLES, FORWARD_COLUMN_RADIANCE, &
          TANLIN_COLUMN_RADIANCE, ADJOINT_COLUMN_RADIANCE, &
          DEALLOC_COLUMN_RADIANCE, DEALLOC_SCATTERING_TABLES


   ! Very small amount to add to hydrometeor mixing ratios (g/kg)
  REAL, PARAMETER :: MINMIXR = 1.0E-6
   ! MAXDELTAU  maximum scattering optical depth allowed in initial sublayers
  REAL, PARAMETER :: MAXDELTAU = 1.0
   ! SPLITTING_FACTOR  max number of internal grid points is this factor
   !                   times the number of input layers + optical depth
  REAL, PARAMETER :: SPLITTING_FACTOR = 10.0
   ! maximum number of SHDOMPP iterations
  INTEGER, PARAMETER :: MAXITER = 150
   ! Flag for SHDOMPP iteration printing
  LOGICAL, PARAMETER :: PRNT=.FALSE.

   ! Scattering table variables (allocated in READ_SCATTERING_TABLES):
  INTEGER :: NSCATTAB, NCHANTAB, MAXLEG, MAXNR
  INTEGER, ALLOCATABLE :: SPECIES2TABLE(:), NRTAB(:,:), NLEGTAB(:,:,:)
  REAL,    ALLOCATABLE :: WAVELEN1(:,:), WAVELEN2(:,:), PARDENS(:)
  REAL,    ALLOCATABLE :: RTAB(:,:,:)
  REAL,    ALLOCATABLE :: EXTINCTTAB(:,:,:), SSALBTAB(:,:,:), LEGENTAB(:,:,:,:)

   ! SHDOMPPDA and MOLEC_ABS_INTERFACE variables 
   !    (allocated in FORWARD_COLUMN_RADIANCE)
  INTEGER :: ML, MM, NLM, NPHI0MAX
  INTEGER :: NPTS, MAXIG, NMOREITER
  LOGICAL, ALLOCATABLE :: SPECIESFLAG(:,:)
  REAL,    ALLOCATABLE :: TEMPP(:)
  REAL,    ALLOCATABLE :: LWP(:,:), RMM(:,:), TAURAYL(:)
  REAL,    ALLOCATABLE :: TAUP(:), ALBEDOP(:), LEGENP(:,:), SSPHASE(:,:)
  INTEGER, ALLOCATABLE :: IXP(:,:)
  REAL,    ALLOCATABLE :: FRACTAUG(:), SHSOURCE(:,:)
  REAL,    ALLOCATABLE :: ACCELPAR(:)
  REAL,    ALLOCATABLE :: PRESLEV(:),PRES2(:), TEMP2(:), VAPMIXR2(:), O3MIXR2(:)
  REAL,    ALLOCATABLE :: MOLECTAU(:)
  REAL,    ALLOCATABLE :: RADIANCE(:)


CONTAINS


SUBROUTINE READ_SCATTERING_TABLES (NCOMP, NCHAN, SPECIESZRANGE, &
                                   SCATTABDIR, SCATTABFILES)
 ! Reads the tables of single scattering properties as a function of
 ! mean mass radius for the needed hydrometeor species and all the 
 ! NCHAN channels (over all platforms).  If the SPECIESZRANGE(2,I) is 
 ! negative for a component than those scattering table files are not read.
 ! The SPECIES2TABLE module array is made, which converts from hydrometeor 
 ! species number (1 to NCOMP) to the scattering table number (1 to NSCATTAB) 
 ! in the tabulated properties arrays.  The module arrays for the tabulated 
 ! scattering properties (NRTAB, RTAB, WAVELEN1, WAVELEN2,  PARDENS, 
 ! EXTINCTTAB, SSALBTAB, NLEGTAB, LEGENTAB) are allocated and filled in here.
 !
 ! Inputs:
 !  NCOMP          number of hydrometeor components or species
 !  NCHAN          number of wavelengths
 !  SPECIESZRANGE  height range (in km) for each hydrometeor species
 !  SCATTABDIR     scattering table directory for all files
 !  SCATTABFILES   scattering table files for all components and channels
 !
 ! Outputs:
 !  No output arguments.  Routine sets global variables in module.
 !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NCOMP, NCHAN
  REAL,    INTENT(IN) :: SPECIESZRANGE(2,NCOMP)
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABDIR, SCATTABFILES(NCOMP,NCHAN)
  INTEGER :: I, J, K, MAXNLEG
  CHARACTER(LEN=256) :: SCATFILE

   ! Make the SPECIES2TABLE index array, which converts from hydrometeor
   !   species number to scattering table number, using SPECIESZRANGE(1,I)
  ALLOCATE (SPECIES2TABLE(NCOMP))
  NSCATTAB = 0
  DO I = 1, NCOMP
    IF (SPECIESZRANGE(2,I) > 0.0) THEN
      NSCATTAB = NSCATTAB + 1
      SPECIES2TABLE(I) = NSCATTAB
    ELSE
      SPECIES2TABLE(I) = 0
    ENDIF
  ENDDO
  IF (NSCATTAB <= 0) THEN
    PRINT *, 'READ_SCATTERING_TABLES: No scattering tables selected due to SPECIESZRANGE'
    STOP
  ENDIF

  NCHANTAB = NCHAN
  ALLOCATE (NRTAB(NSCATTAB,NCHAN))
  ALLOCATE (WAVELEN1(NSCATTAB,NCHAN), WAVELEN2(NSCATTAB,NCHAN), PARDENS(NCOMP))
  MAXLEG = 1
  DO J = 1, NCHAN
    DO I = 1, NCOMP
      K = SPECIES2TABLE(I)
      IF (K > 0) THEN
        WRITE (SCATFILE,'(A,A)') TRIM(SCATTABDIR), TRIM(SCATTABFILES(I,J))
        CALL READ_SCAT_TABLE_SIZE (SCATFILE, WAVELEN1(K,J), WAVELEN2(K,J), &
                                   PARDENS(I), NRTAB(K,J), MAXNLEG)
        MAXLEG = MAX(MAXLEG,MAXNLEG)
        MAXNR = MAX(MAXNR,NRTAB(K,J))
        IF (K > 1 .AND. ABS((WAVELEN1(K,J)-WAVELEN1(1,J))/WAVELEN1(1,J)) > 0.01 &
           .OR. ABS((WAVELEN2(K,J)-WAVELEN2(1,J))/WAVELEN2(1,J)) > 0.01) THEN
          PRINT *, 'Warning: scattering table wavelengths do not agree for table',I,J
        ENDIF
      ENDIF
    ENDDO
  ENDDO
  ALLOCATE (RTAB(MAXNR,NSCATTAB,NCHAN))
  ALLOCATE (EXTINCTTAB(MAXNR,NSCATTAB,NCHAN), SSALBTAB(MAXNR,NSCATTAB,NCHAN))
  ALLOCATE (NLEGTAB(MAXNR,NSCATTAB,NCHAN), LEGENTAB(1:MAXLEG,MAXNR,NSCATTAB,NCHAN))

  ! Read in the optical properties from the scattering table files
  DO J = 1, NCHAN
    DO I = 1, NCOMP
      K = SPECIES2TABLE(I)
      IF (K > 0) THEN
        WRITE (SCATFILE,'(A,A)') TRIM(SCATTABDIR), TRIM(SCATTABFILES(I,J))
        CALL READ_SCAT_TABLE (SCATFILE, NRTAB(K,J), RTAB(:,K,J), &
                          EXTINCTTAB(:,K,J), SSALBTAB(:,K,J), &   
                          MAXLEG, NLEGTAB(:,K,J), LEGENTAB(:,:,K,J) )
      ENDIF
    ENDDO
  ENDDO

END SUBROUTINE READ_SCATTERING_TABLES 





SUBROUTINE FORWARD_COLUMN_RADIANCE (NZ, ZLEV, TEMP, PRES, VAPMIXR, O3MIXR, &
                                    TOPMOLECTAU, TOPTEMP, INDEXABS, &
                                    NCOMP, MIXR, CONC, SPECIESZRANGE, &
                                    SFCPRES, TEMPSKIN, SFCALBEDO, &
                                    ICHANSCAT, NOUT, MUOBS, PHIOBSREL, &
                                    SRCTYPE, SOLARFLUX, SOLARMU, WAVELEN, &
                                    NMU, NPHI, SPLITACC, SOLACC, &
                                    UNITS, RADOUT)
 ! Computes the upwelling radiance for the specified atmospheric state
 ! in the specified directions.
 ! Allocates many of the module global arrays (deallocate with 
 ! DEALLOC_COLUMN_RADIANCE).
 ! Model levels are treated as SHDOMPP layer boundary levels with the
 ! cloud mixing ratio linearly interpolated and integrated to get the
 ! hydrometeor mass in each layer.  The hydrometeor water paths are
 ! obtained from the mixing ratios and pressure profiles; ZLEV is
 ! used only for speciesZrange, not for optical properties.
 !
 ! Inputs:
 !  NZ         number of model levels
 !  ZLEV       model height levels (km) (from bottom to top of atmosphere)
 !  TEMP       model temperature profile (K)
 !  PRES       model pressure profile (mb)
 !  VAPMIXR    model water vapor mass mixing ratio (g/kg)
 !  O3MIXR     ozone mixing ratio (ppmv) at model levels
 !  TOPMOLECTAU equivalent top layer (above model top) molecular optical depth 
 !  TOPTEMP     equivalent top layer temperature (K)
 !  INDEXABS    absorption channel index number (<=0 for no absorption call)
 !  NCOMP      number of hydrometeor species input
 !  MIXR       model hydrometeor mass mixing ratio profiles (g/kg)
 !  CONC       model hydrometeor concentration profile for each species (#/kg) 
 !  SPECIESZRANGE  lower and upper range of layers to allow each species (km)
 !  SFCPRES    model surface pressure (mb)
 !  TEMPSKIN   model surface skin temperature (K)
 !  SFCALBEDO  surface Lambertian albedo
 !  ICHANSCAT  channel number for indexing into scattering tables
 !  NOUT       number of output radiance directions
 !  MUOBS      viewing cosine zenith angles
 !  PHIOBSREL  viewing azimuth angles relative to solar azimuth (deg; 0 is forward scatter)
 !  SRCTYPE    radiation source type: 'S'=solar, 'T'=thermal, 'B'=both
 !  SOLARFLUX  solar flux on a *horizontal* surface (W/m^2 um)
 !  SOLARMU    cosine solar zenith angle (must be < 0)
 !  WAVELEN    wavelength (microns) for thermal source
 !  NMU        number of discrete ordinate zenith angles in both hemispheres
 !  NPHI       number of discrete ordinate azimuth angles in 2\pi
 !  SPLITACC   adaptive layer splitting accuracy (0.03 to 0.0003 recommended)
 !  SOLACC     SHDOMPP solution accuracy (1E-6 to 1E-7 recommended for adjoint
 !               accuracy even though that is not needed for radiance accuracy)
 !  UNITS      output units: 'R' for reflectance (relative to SOLARFLUX), 
 !                           'T' for EBB brightness temperature
 !                           'W' for radiance output (W/(m^2 ster um))
 ! Outputs:
 !  RADOUT     reflectance or brightness temperature vector

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NZ, NCOMP, INDEXABS, ICHANSCAT, NOUT
  REAL,    INTENT(IN) :: ZLEV(NZ), TEMP(NZ), PRES(NZ), VAPMIXR(NZ), O3MIXR(NZ)
  REAL,    INTENT(IN) :: TOPMOLECTAU, TOPTEMP
  REAL,    INTENT(IN) :: MIXR(NCOMP,NZ), CONC(NCOMP,NZ)
  REAL,    INTENT(IN) :: SPECIESZRANGE(2,NCOMP)
  REAL,    INTENT(IN) :: SFCPRES, TEMPSKIN, SFCALBEDO
  INTEGER, INTENT(IN) :: NMU, NPHI
  REAL,    INTENT(IN) :: MUOBS(NOUT), PHIOBSREL(NOUT)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, WAVELEN, SPLITACC, SOLACC
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, UNITS
  REAL,    INTENT(OUT) :: RADOUT(NOUT)
  INTEGER :: NLAY, L, LR, I, ITER
  LOGICAL :: ACCELFLAG
  REAL    :: MIXAVG, CONCAVG, RMM_MID, DCONC
  REAL    :: MU0, WAVENO(2), SKYRAD, SFCTEMP, SFCPARMS(3)
  REAL    :: SOLCRIT
  CHARACTER(LEN=1) :: SFCTYPE

  NLAY = NZ-1+1   ! NLAY is number of SHDOMPP layers including added top layer
  ALLOCATE (TEMPP(NLAY+1))
  ALLOCATE (SPECIESFLAG(NCOMP,NLAY))
   ! Allocate the single column particle LWP and mean mass radius arrays
  ALLOCATE (LWP(NCOMP,NLAY), RMM(NCOMP,NLAY), MOLECTAU(NLAY), TAURAYL(NLAY))
   ! Allocate the molecular absorption profile arrays (these go from top down)
  ALLOCATE (PRESLEV(0:NZ-1), PRES2(NZ-1))
  ALLOCATE (TEMP2(NZ-1), VAPMIXR2(NZ-1), O3MIXR2(NZ-1))

   ! Reverse the temperature levels for SHDOMPP
  TEMPP(2:NLAY+1) = TEMP(NZ:1:-1)
   ! Set top level temp to get TOPTEMP average
  TEMPP(1) = 2*TOPTEMP - TEMPP(2)

   ! Call molecular absorption routine if channel index so indicates
  IF (INDEXABS > 0) THEN
     ! Fill in the layer center arrays that go from top down
    DO L = 1, NZ-1
      PRES2(L) = 0.5*(PRES(NZ-L)+PRES(NZ-L+1))
      TEMP2(L) = 0.5*(TEMP(NZ-L)+TEMP(NZ-L+1))
      VAPMIXR2(L) = 0.5*(VAPMIXR(NZ-L)+VAPMIXR(NZ-L+1))
      O3MIXR2(L) = 0.5*(O3MIXR(NZ-L)+O3MIXR(NZ-L+1))
    ENDDO
    PRESLEV(0:NZ-2) = PRES(NZ:2:-1)
    PRESLEV(NZ-1) = SFCPRES   ! Use surface pressure instead of model pres(1)
     ! Call interface to get molecular absorption layer optical depth
    CALL FORWARD_MOLEC_ABSORPTION (INDEXABS, NZ-1, PRESLEV, &
                                   PRES2, TEMP2, VAPMIXR2, O3MIXR2, &
                                   MOLECTAU(2:NLAY))
    MOLECTAU(1) = TOPMOLECTAU
  ELSE
    MOLECTAU(:) = 0.0
  ENDIF

   ! Convert from particle mass mixing ratio and number concentration
   !  to layer integrated mass content and mean mass radius.  Add in a
   !  very small amount of mixing ratio and concentration for clear layers.
  DO L = 2, NLAY
    LR = NLAY-L+1
    DO I = 1, NCOMP
      IF (SPECIES2TABLE(I) > 0) THEN
        MIXAVG = 0.5*(MIXR(I,LR)+MIXR(I,LR+1)) + MINMIXR
        RMM_MID = RTAB(NRTAB(I,1)/2+1,I,ICHANSCAT) + 0.1
        DCONC = MINMIXR/((4*3.14159*PARDENS(I)/3.0)*(RMM_MID/1.0E4)**3)
        CONCAVG = 0.5*(CONC(I,LR)+CONC(I,LR+1)) + DCONC
        LWP(I,L) = MIXAVG*(PRES(LR)-PRES(LR+1))*100/9.80
        RMM(I,L) = 1.0E4*( ((0.75/3.14159)/PARDENS(I)) &
                            *(MIXAVG/CONCAVG) )**(1.0/3)
      ELSE
        LWP(I,L) = 0.0
        RMM(I,L) = 0.0
      ENDIF
    ENDDO
  ENDDO
  LWP(:,1) = 0.0
  RMM(:,1) = 0.0


   ! Set the flags that decide where the hydrometeor species can be changed
  DO I = 1, NCOMP
    SPECIESFLAG(I,1) = .FALSE.
    DO L = 2, NLAY
      SPECIESFLAG(I,L) = .FALSE.
      LR = NLAY-L+1
      IF (ZLEV(LR)  <= SPECIESZRANGE(2,I) .AND. &
          ZLEV(LR+1) > SPECIESZRANGE(1,I))  SPECIESFLAG(I,L) = .TRUE.
    ENDDO
  ENDDO

  IF (ICHANSCAT < 1 .OR. ICHANSCAT > NCHANTAB) THEN
    WRITE (*,*) 'FORWARD_COLUMN_RADIANCE: Illegal ICHANSCAT: ', ICHANSCAT, NCHANTAB
    STOP
  ENDIF

   ! Calculate the spherical harmonic truncation from NMU and NPHI
  IF (NMU /= MAX(2,2*INT((NMU+1)/2)) ) STOP 'SOLVE_SHDOMPP: Bad NMU'
  IF (NPHI < 1) STOP 'SOLVE_SHDOMPP: NPHI must be greater than 0'
  ML = NMU-1
  MM = MAX(0,INT(NPHI/2)-1)
  NLM = (MM+1)*(ML+1) - (MM*(MM+1))/2
  NPHI0MAX = INT((NPHI+2)/2)
  IF (ANY(MUOBS(:) == 0.0) .OR. ANY(ABS(MUOBS(:)) > 1.0)) THEN
    WRITE (*,*) 'FORWARD_COLUMN_RADIANCE: Bad mu for radiance output'
    STOP
  ENDIF      
  IF (SRCTYPE /= "T" .AND. SOLARMU >= 0.0) &
      STOP 'FORWARD_COLUMN_RADIANCE: SOLARMU must be negative'


   ! Allocate optical properties arrays
  ALLOCATE (TAUP(NLAY), ALBEDOP(NLAY))
  ALLOCATE (LEGENP(1:ML+1,NLAY), SSPHASE(NLAY,NOUT))
  ALLOCATE (ACCELPAR(MAXITER))
  ACCELPAR(:) = 0.0

  IF (WAVELEN < 1.0) THEN
    TAURAYL(1) = PRES(NZ)
    TAURAYL(2:NLAY-1) = PRES(NZ-1:2:-1)-PRES(NZ:3:-1)
    TAURAYL(NLAY) = SFCPRES - PRES(2)
    TAURAYL(:) = (TAURAYL(:)/1013.25)* 0.0088*WAVELEN**(-4.15+0.2*WAVELEN)
  ELSE
    TAURAYL(1:NLAY) = 0.0
  ENDIF

   ! Calculate the optical properties from the particle properties and
   !  the scattering tables.
  CALL CALCULATE_OPTICAL_PROPERTIES (NLAY, NCOMP, ML, SPECIESFLAG, &
                                     LWP, RMM, MOLECTAU, TAURAYL, &
                                     SRCTYPE, SOLARMU, NOUT, MUOBS, PHIOBSREL, &
                                     TAUP, ALBEDOP, LEGENP, SSPHASE, &
                     NSCATTAB, SPECIES2TABLE, MAXNR, &
                     NRTAB(:,ICHANSCAT), RTAB(:,:,ICHANSCAT), &
                     EXTINCTTAB(:,:,ICHANSCAT), SSALBTAB(:,:,ICHANSCAT), &
                     MAXLEG, NLEGTAB(:,:,ICHANSCAT), LEGENTAB(:,:,:,ICHANSCAT))

!    print *, 'Lay Pres   Temp  Tau   Albedo  Legen1  SSphase'
!    DO L = 1, NLAY
!       print '(1X,I3,1X,F6.1,1X,F6.2,1X,F7.4,3(1X,F7.5))', &
!            L, PRES(NZ-L+1), TEMPP(L), TAUP(L), ALBEDOP(L), LEGENP(1,L), SSPHASE(L,1)
!    ENDDO
!    print '(1X,I3,1X,F6.1,1X,F6.2)', NLAY+1, SFCPRES, TEMPP(NLAY+1)
!    print *, ' '

   ! Guess the maximum size of the adaptive grid and allocate arrays
  MAXIG = NINT(SPLITTING_FACTOR*(NLAY+SUM(TAUP(:)*ALBEDOP(:))/MAXDELTAU))
  ALLOCATE (IXP(2,NLAY), FRACTAUG(MAXIG), SHSOURCE(NLM,MAXIG), RADIANCE(NOUT))

  SKYRAD = 0.0         ! could be Planck function for 2.7 K
  MU0 = -ABS(SOLARMU)
  SFCTEMP = TEMPSKIN
  SFCTYPE = 'L'   ! Lambertian surface 
  SFCPARMS(1) = SFCALBEDO
  ACCELFLAG = .TRUE.

  CALL SOLVE_ADAPTIVE_SHDOMPP (NLAY, TEMPP, TAUP, ALBEDOP, LEGENP, &
                               MAXIG, ML, MM, NLM, NMU, NPHI, NPHI0MAX, &
                               SRCTYPE, SOLARFLUX, MU0, &
                               WAVELEN, WAVENO, SKYRAD, &
                               SFCTEMP, SFCTYPE, SFCPARMS, &
                               MAXDELTAU, SPLITACC, MAXITER, ACCELFLAG, PRNT, &
                               NPTS, FRACTAUG, IXP, SHSOURCE, &
                               SOLCRIT, ITER)

  CALL SOLVE_FIXED_SHDOMPP1 (NLAY, TEMPP, &
                            TAUP, ALBEDOP, LEGENP, SSPHASE, &
                            NPTS, FRACTAUG, IXP, SHSOURCE, &
                            ML, MM, NLM, NMU, NPHI, NPHI0MAX, &
                            SRCTYPE, SOLARFLUX, MU0, &
                            WAVELEN, WAVENO, SKYRAD, &
                            SFCTEMP, SFCTYPE, SFCPARMS, &
                            SOLACC, MAXITER, ACCELFLAG, PRNT, &
                            SOLCRIT, ITER, NOUT, MUOBS, PHIOBSREL, &
                            RADIANCE, NMOREITER, ACCELPAR)
  NMOREITER = MAX(NMOREITER,2)

   ! Convert radiance to brightness temperature or reflectance if desired
  IF (UNITS == 'T') THEN
    RADOUT(:) = 1.4388E4/(WAVELEN*LOG( 1 + 1.1911E8/(WAVELEN**5*RADIANCE(:)) ))
  ELSE IF (UNITS == 'R') THEN
    RADOUT(:) = 3.14159*RADIANCE(:)/SOLARFLUX
  ELSE
    RADOUT(:) = RADIANCE(:)
  ENDIF

END SUBROUTINE FORWARD_COLUMN_RADIANCE






SUBROUTINE TANLIN_COLUMN_RADIANCE (NZ, ZLEV, TEMP, PRES, VAPMIXR, O3MIXR, &
                                    G_TEMP, G_VAPMIXR, G_O3MIXR, &
                                    INDEXABS, &
                                    NCOMP, MIXR, CONC, G_MIXR, G_CONC, &
                                    SFCPRES, TEMPSKIN, SFCALBEDO, &
                                    G_TEMPSKIN, G_SFCALBEDO, &
                                    ICHANSCAT, NOUT, MUOBS, PHIOBSREL, &
                                    SRCTYPE, SOLARFLUX, SOLARMU, WAVELEN, &
                                    NMU, NPHI, UNITS, G_RADOUT)
 ! Computes the tangent linear upwelling radiance for the specified 
 ! atmospheric state.
 ! Must be called after a call to FORWARD_COLUMN_RADIANCE with all the
 ! same parameters.
 !
 ! Inputs:
 !  NZ         number of model levels
 !  ZLEV      model height levels (km) (from bottom to top of atmosphere)
 !  TEMP       model temperature profile (K)
 !  PRES       model pressure profile (mb)
 !  VAPMIXR    model water vapor mass mixing ratio (g/kg)
 !  O3MIXR     ozone mixing ratio (ppmv) at model levels
 !  G_TEMP     perturbation temperature profile (K)
 !  G_VAPMIXR  perturbation water vapor mass mixing ratio (g/kg)
 !  G_O3MIXR   perturbation ozone mixing ratio (ppmv)
 !  INDEXABS    absorption channel index number (<=0 for no absorption call)
 !  NCOMP      number of hydrometeor species input (usually 7 for model)
 !  MIXR       model hydrometeor mass mixing ratio profiles (g/kg)
 !  CONC       model hydrometeor concentration profile for each species (#/kg) 
 !  G_MIXR     perturbation hydrometeor mass mixing ratio profiles (g/kg)
 !  G_CONC     perturbation hydrometeor concentration profile for each species (#/kg) 
 !  SFCPRES    model surface pressure (mb)
 !  TEMPSKIN   model surface skin temperature (K)
 !  SFCALBEDO  surface Lambertian albedo
 !  G_TEMPSKIN  perturbation surface skin temperature (K)
 !  G_SFCALBEDO perturbation surface Lambertian albedo
 !  ICHANSCAT  channel number for indexing into scattering tables
 !  NOUT       number of output radiance directions
 !  MUOBS      viewing cosine zenith angles
 !  PHIOBSREL  viewing azimuth angles relative to solar azimuth (deg; 0 is forward scatter)
 !  SRCTYPE    radiation source type: 'S'=solar, 'T'=thermal, 'B'=both
 !  SOLARFLUX  solar flux on a *horizontal* surface (W/m^2 um)
 !  SOLARMU    cosine solar zenith angle (must be < 0)
 !  WAVELEN    wavelength (microns) for thermal source
 !  NMU        number of discrete ordinate zenith angles in both hemispheres
 !  NPHI       number of discrete ordinate azimuth angles in 2\pi
 !  UNITS      output units: 'R' for reflectance (relative to SOLARFLUX), 
 !                           'T' for EBB brightness temperature
 !                           'W' for radiance output (W/(m^2 ster um))
 ! Outputs:
 !  G_RADOUT   tangent linear reflectance or brightness temperature vector

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NZ, NCOMP, INDEXABS, ICHANSCAT, NOUT
  REAL,    INTENT(IN) :: ZLEV(NZ), PRES(NZ)
  REAL,    INTENT(IN) :: TEMP(NZ), VAPMIXR(NZ), O3MIXR(NZ)
  REAL,    INTENT(IN) :: G_TEMP(NZ), G_VAPMIXR(NZ), G_O3MIXR(NZ)
  REAL,    INTENT(IN) :: MIXR(NCOMP,NZ), CONC(NCOMP,NZ)
  REAL,    INTENT(IN) :: G_MIXR(NCOMP,NZ), G_CONC(NCOMP,NZ)
  REAL,    INTENT(IN) :: SFCPRES, TEMPSKIN, SFCALBEDO
  REAL,    INTENT(IN) :: G_TEMPSKIN, G_SFCALBEDO
  INTEGER, INTENT(IN) :: NMU, NPHI
  REAL,    INTENT(IN) :: MUOBS(NOUT), PHIOBSREL(NOUT)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, WAVELEN
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, UNITS
  REAL,    INTENT(OUT) :: G_RADOUT(NOUT)
  INTEGER :: NLAY, L, LR, I, ITER
  LOGICAL :: ACCELFLAG
  REAL    :: MIXAVG, CONCAVG, RMM_MID, DCONC, G_MIXAVG, G_CONCAVG
  REAL    :: MU0, WAVENO(2), SKYRAD
  REAL    :: SFCTEMP, SFCPARMS(3), G_SFCTEMP, G_SFCPARMS(3)
  REAL    :: SOLCRIT
  CHARACTER(LEN=1) :: SFCTYPE
  REAL,    ALLOCATABLE :: G_TEMPP(:), G_LWP(:,:), G_RMM(:,:)
  REAL,    ALLOCATABLE :: G_MOLECTAU(:), G_TAURAYL(:)
  REAL,    ALLOCATABLE :: G_TEMP2(:), G_VAPMIXR2(:), G_O3MIXR2(:)
  REAL,    ALLOCATABLE :: G_TAUP(:), G_ALBEDOP(:), G_LEGENP(:,:), G_SSPHASE(:,:)
  REAL,    ALLOCATABLE :: RADIANCE2(:), G_RADIANCE(:)

  NLAY = NZ-1+1   ! NLAY is number of SHDOMPP layers including added top layer
   ! Allocate temporary arrays 
  ALLOCATE (G_TEMPP(NLAY+1), G_LWP(NCOMP,NLAY), G_RMM(NCOMP,NLAY))
  ALLOCATE (G_MOLECTAU(NLAY), G_TAURAYL(NLAY))
  ALLOCATE (G_TEMP2(NZ-1), G_VAPMIXR2(NZ-1), G_O3MIXR2(NZ-1))
  ALLOCATE (G_TAUP(NLAY), G_ALBEDOP(NLAY))
  ALLOCATE (G_LEGENP(1:ML+1,NLAY), G_SSPHASE(NLAY,NOUT))
  ALLOCATE (RADIANCE2(NOUT), G_RADIANCE(NOUT))

   ! Reverse the temperature levels for SHDOMPP
  G_TEMPP(2:NLAY+1) = G_TEMP(NZ:1:-1)
  G_TEMPP(1) = - G_TEMPP(2)

   ! Call molecular absorption routine if channel index so indicates
  IF (INDEXABS > 0) THEN
     ! Fill in the layer center arrays that go from top down
    DO L = 1, NZ-1
      G_TEMP2(L) = 0.5*(G_TEMP(NZ-L)+G_TEMP(NZ-L+1))
      G_VAPMIXR2(L) = 0.5*(G_VAPMIXR(NZ-L)+G_VAPMIXR(NZ-L+1))
      G_O3MIXR2(L) = 0.5*(G_O3MIXR(NZ-L)+G_O3MIXR(NZ-L+1))
    ENDDO
     ! Call interface to get molecular absorption layer optical depth
    CALL TANLIN_MOLEC_ABSORPTION (INDEXABS, NZ-1, PRESLEV, &
                                  PRES2, TEMP2, VAPMIXR2, O3MIXR2, &
                                  G_TEMP2, G_VAPMIXR2, G_O3MIXR2, &
                                  G_MOLECTAU(2:NLAY))
  ELSE
    G_MOLECTAU(:) = 0.0
  ENDIF

   ! Gradient of layer integrated mass content and mean mass radius
   !   from particle mass mixing ratio and number concentration.
  DO L = 2, NLAY
    LR = NLAY-L+1
    DO I = 1, NCOMP
      IF (SPECIES2TABLE(I) > 0) THEN
        MIXAVG = 0.5*(MIXR(I,LR)+MIXR(I,LR+1)) + MINMIXR
        RMM_MID = RTAB(NRTAB(I,1)/2+1,I,ICHANSCAT) + 0.1
        DCONC = MINMIXR/((4*3.14159*PARDENS(I)/3.0)*(RMM_MID/1.0E4)**3)
        CONCAVG = 0.5*(CONC(I,LR)+CONC(I,LR+1)) + DCONC
        G_MIXAVG = 0.5*(G_MIXR(I,LR)+G_MIXR(I,LR+1))
        G_CONCAVG = 0.5*(G_CONC(I,LR)+G_CONC(I,LR+1))
        G_LWP(I,L) = G_MIXAVG*(PRES(LR)-PRES(LR+1))*100/9.80
        G_RMM(I,L) = RMM(I,L)/(3.0*MIXAVG) * G_MIXAVG &
                   - RMM(I,L)/(3.0*CONCAVG) * G_CONCAVG
      ELSE
        G_LWP(I,L) = 0.0
        G_RMM(I,L) = 0.0
      ENDIF
    ENDDO
  ENDDO
  G_LWP(:,1) = 0.0
  G_RMM(:,1) = 0.0

   ! No Rayleigh scattering gradient since we assume no pressure perturbation
  G_TAURAYL(:) = 0.0

   ! Calculate the optical properties from the particle properties and
   !  the scattering tables.
  CALL G_CALCULATE_OPTICAL_PROPERTIES (NLAY, NCOMP, ML, SPECIESFLAG, &
                                     LWP, G_LWP, RMM, G_RMM, &
                                     MOLECTAU, G_MOLECTAU, TAURAYL, G_TAURAYL, &
                                     SRCTYPE, SOLARMU, NOUT, MUOBS, PHIOBSREL, &
                                     TAUP, G_TAUP, ALBEDOP, G_ALBEDOP, &
                                     LEGENP, G_LEGENP, SSPHASE, G_SSPHASE, &
                     NSCATTAB, SPECIES2TABLE, MAXNR, &
                     NRTAB(:,ICHANSCAT), RTAB(:,:,ICHANSCAT), &
                     EXTINCTTAB(:,:,ICHANSCAT), SSALBTAB(:,:,ICHANSCAT), &
                     MAXLEG, NLEGTAB(:,:,ICHANSCAT), LEGENTAB(:,:,:,ICHANSCAT))

!    PRINT *
!    PRINT *, '      Optical Properties and Their Perturbations'
!    PRINT *,'Layer   Tau     dTau     SSalb   dSSalb     Legen1  dLegen1   SSph$
!    DO I = 1, NLAY
!      WRITE (*,'(1X,I3,2X,F7.3,1X,F8.5,3(1X,2(1X,F8.6)))') &
!           I, TAUP(I), G_TAUP(I), ALBEDOP(I), G_ALBEDOP(I), &
!              LEGENP(1,I), G_LEGENP(1,I), SSPHASE(I,NOUT), G_SSPHASE(I,NOUT)
!    ENDDO
!    PRINT *

  SKYRAD = 0.0
  MU0 = -ABS(SOLARMU)
  SFCTEMP = TEMPSKIN
  G_SFCTEMP = G_TEMPSKIN
  SFCTYPE = 'L'   ! Lambertian surface 
  SFCPARMS(1) = SFCALBEDO
  G_SFCPARMS(1) = G_SFCALBEDO

   ! Call tangent linear for fixed grid part of SHDOMPP iterations
  CALL G_SOLVE_FIXED_SHDOMPP2 (NLAY, TEMPP, G_TEMPP, TAUP, G_TAUP, &
                              ALBEDOP, G_ALBEDOP, LEGENP, G_LEGENP, & 
                              SSPHASE, G_SSPHASE, &
                              NPTS, FRACTAUG, IXP, SHSOURCE, &
                              ML, MM, NLM, NMU, NPHI, NPHI0MAX, &
                              SRCTYPE, SOLARFLUX, MU0, &
                              WAVELEN, WAVENO, SKYRAD, &
                              SFCTEMP,G_SFCTEMP, SFCTYPE, SFCPARMS,G_SFCPARMS, &
                              NMOREITER, ACCELPAR, &
                              NOUT, MUOBS, PHIOBSREL, RADIANCE2, G_RADIANCE)

   ! Convert radiance gradient to brightness temperature or reflectance gradient
  IF (UNITS == 'T') THEN
    G_RADOUT(:) = G_RADIANCE(:)* 1.4388E4*1.1911E8 / &
                ( WAVELEN*RADIANCE(:)*(1.1911E8+WAVELEN**5*RADIANCE(:)) &
                  *LOG(1.0+1.1911E8/(WAVELEN**5*RADIANCE(:)))**2 )
  ELSE IF (UNITS == 'R') THEN
    G_RADOUT(:) = 3.14159*G_RADIANCE(:)/SOLARFLUX
  ELSE
    G_RADOUT(:) = G_RADIANCE(:)
  ENDIF

  DEALLOCATE (RADIANCE2, G_RADIANCE)
  DEALLOCATE (G_TEMP2, G_VAPMIXR2, G_O3MIXR2)
  DEALLOCATE (G_TAUP, G_ALBEDOP, G_LEGENP, G_SSPHASE)
  DEALLOCATE (G_TEMPP, G_MOLECTAU, G_TAURAYL, G_LWP, G_RMM)
END SUBROUTINE TANLIN_COLUMN_RADIANCE







SUBROUTINE ADJOINT_COLUMN_RADIANCE (NZ, ZLEV, TEMP, PRES, VAPMIXR, O3MIXR, &
                                    INDEXABS, &
                                    NCOMP, MIXR, CONC, &
                                    SFCPRES, TEMPSKIN, SFCALBEDO, &
                                    ICHANSCAT, NOUT, MUOBS, PHIOBSREL, &
                                    SRCTYPE, SOLARFLUX, SOLARMU, WAVELEN, &
                                    NMU, NPHI,  UNITS,  ADRADOUT, &
                                    ADTEMP, ADVAPMIXR, ADO3MIXR, &
                                    ADMIXR, ADCONC, ADTEMPSKIN, ADSFCALBEDO)
 ! Computes the single column radiative transfer adjoint. 
 ! Must be called after a call to FORWARD_COLUMN_RADIANCE with all the
 ! same parameters.
 !
 ! Inputs:
 !  NZ         number of model midpoint levels and SHDOMPP layers
 !  ZLEV      model height levels (km)
 !  TEMP       model temperature profile (K)
 !  PRES       model pressure profile (mb)
 !  VAPMIXR    model water vapor mass mixing ratio (g/kg)
 !  O3MIXR     ozone mixing ratio (ppmv) at model levels
 !  INDEXABS   absorption channel index number (<=0 for no call)
 !  NCOMP      number of hydrometeor species input (usually 7 for model)
 !  MIXR       model hydrometeor mass mixing ratio profiles (g/kg)
 !  CONC       model hydrometeor concentration profile for each species (#/kg) 
 !  SFCPRES    model surface pressure (mb)
 !  TEMPSKIN   model surface skin temperature (K)
 !  SFCALBEDO  surface Lambertian albedo
 !  ICHANSCAT  channel number for indexing into scattering tables
 !  NOUT       number of output radiance directions
 !  MUOBS      viewing cosine zenith angles
 !  PHIOBSREL  viewing azimuth angles relative to solar azimuth (deg; 0 is forward scatter)
 !  SRCTYPE    radiation source type: 'S'=solar, 'T'=thermal, 'B'=both
 !  SOLARFLUX  solar flux on a *horizontal* surface (W/m^2 um)
 !  SOLARMU    cosine solar zenith angle (forced to be < 0)
 !  WAVELEN    wavelength (microns) for thermal source
 !  NMU        number of discrete ordinate zenith angles in both hemispheres
 !  NPHI       number of discrete ordinate azimuth angles in 2\pi
 !  UNITS      output units: 'R' for reflectance (relative to SOLARFLUX), 
 !                           'T' for EBB brightness temperature
 !  ADRADOUT   reflectance or brightness temperature "forcing" vector
 !
 ! Adjoint Outputs:
 !  ADTEMP      adjoint temperature profile (model grid)
 !  ADVAPMIXR   adjoint water vapor profile (model grid)
 !  ADO3MIXR    adjoint ozone profile (model grid)
 !  ADMIXR      adjoint hydrometeor mass mixing ratio profiles (model grid)
 !  ADCONC      adjoint hydrometeor concentration profile for each species
 !  ADTEMPSKIN  adjoint skin temperature
 !  ADSFCALBEDO adjoint surface albedo

  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NZ, NCOMP, INDEXABS, ICHANSCAT, NOUT
  REAL,    INTENT(IN) :: ZLEV(NZ), TEMP(NZ), PRES(NZ), VAPMIXR(NZ), O3MIXR(NZ)
  REAL,    INTENT(IN) :: MIXR(NCOMP,NZ), CONC(NCOMP,NZ)
  REAL,    INTENT(IN) :: SFCPRES, TEMPSKIN, SFCALBEDO
  INTEGER, INTENT(IN) :: NMU, NPHI
  REAL,    INTENT(IN) :: MUOBS(NOUT), PHIOBSREL(NOUT)
  REAL,    INTENT(IN) :: SOLARFLUX, SOLARMU, WAVELEN
  CHARACTER(LEN=1), INTENT(IN) :: SRCTYPE, UNITS
  REAL,    INTENT(IN) :: ADRADOUT(NOUT)
  REAL,    INTENT(INOUT) :: ADTEMP(NZ), ADVAPMIXR(NZ), ADO3MIXR(NZ)
  REAL,    INTENT(INOUT) :: ADMIXR(NCOMP,NZ), ADCONC(NCOMP,NZ)
  REAL,    INTENT(INOUT) :: ADTEMPSKIN, ADSFCALBEDO
  INTEGER :: NLAY, L, LR, I
  REAL    :: MIXAVG, CONCAVG, ADMIXAVG, ADCONCAVG, RMM_MID, DCONC
  REAL    :: MU0, WAVENO(2), SKYRAD
  REAL    :: SFCTEMP, ADSFCTEMP, SFCPARMS(3), ADSFCPARMS(3)
  CHARACTER(LEN=1) :: SFCTYPE
  REAL,    ALLOCATABLE :: ADTEMPP(:), ADLWP(:,:), ADRMM(:,:)
  REAL,    ALLOCATABLE :: ADMOLECTAU(:), ADTAURAYL(:)
  REAL,    ALLOCATABLE :: ADTAUP(:), ADALBEDOP(:), ADLEGENP(:,:), ADSSPHASE(:,:)
  REAL,    ALLOCATABLE :: ADTEMP2(:), ADVAPMIXR2(:), ADO3MIXR2(:)
  REAL,    ALLOCATABLE :: ADRADIANCE(:)

  NLAY = NZ-1+1   ! NLAY is number of SHDOMPP layers including added top layer

   ! Allocate temporary arrays 
  ALLOCATE (ADTEMPP(NLAY+1), ADMOLECTAU(NLAY), ADTAURAYL(NLAY))
  ALLOCATE (ADLWP(NCOMP,NLAY), ADRMM(NCOMP,NLAY))
  ALLOCATE (ADTAUP(NLAY), ADALBEDOP(NLAY))
  ALLOCATE (ADLEGENP(1:ML+1,NLAY), ADSSPHASE(NLAY,NOUT))
  ALLOCATE (ADTEMP2(NZ-1), ADVAPMIXR2(NZ-1), ADO3MIXR2(NZ-1))
  ALLOCATE (ADRADIANCE(NOUT))
   ! Zero the adjoint arrays
  ADLWP(:,:) = 0.0
  ADRMM(:,:) = 0.0
  ADMOLECTAU(:) = 0.0
  ADTAURAYL(:) = 0.0
  ADTEMPP(:) = 0.0
  ADSFCTEMP = 0.0
  ADSFCPARMS(:) = 0.0
  ADTAUP(:) = 0.0
  ADALBEDOP(:) = 0.0
  ADLEGENP(:,:) = 0.0
  ADSSPHASE(:,:) = 0.0
  ADTEMP2(:) = 0.0
  ADVAPMIXR2(:) = 0.0
  ADO3MIXR2(:) = 0.0

   ! Convert adjoint radiance from brightness temperature or reflectance
  IF (UNITS == 'T') THEN
    ADRADIANCE(:) = ADRADOUT(:)* 1.4388E4*1.1911E8 / &
                ( WAVELEN*RADIANCE(:)*(1.1911E8+WAVELEN**5*RADIANCE(:)) &
                  *LOG(1.0+1.1911E8/(WAVELEN**5*RADIANCE(:)))**2 )
  ELSE IF (UNITS == 'R') THEN
    ADRADIANCE(:) = 3.14159*ADRADOUT(:)/SOLARFLUX
  ELSE
    ADRADIANCE(:) = ADRADOUT(:)
  ENDIF

  SKYRAD = 0.0
  MU0 = -ABS(SOLARMU)
  SFCTEMP = TEMPSKIN
  SFCTYPE = 'L'   ! Lambertian surface 
  SFCPARMS(1) = SFCALBEDO

    ! Call adjoint of fixed grid part of SHDOMPP iterations
  CALL ADSOLVE_FIXED_SHDOMPP2 (NLAY, TEMPP, ADTEMPP, TAUP, ADTAUP, &
                              ALBEDOP, ADALBEDOP, LEGENP, ADLEGENP, &
                              SSPHASE, ADSSPHASE, &
                              NPTS, FRACTAUG, IXP, SHSOURCE, &
                              ML, MM, NLM, NMU, NPHI, NPHI0MAX, &
                              SRCTYPE, SOLARFLUX, MU0, &
                              WAVELEN, WAVENO, SKYRAD, &
                              SFCTEMP,ADSFCTEMP, SFCTYPE, SFCPARMS,ADSFCPARMS, &
                              NMOREITER, ACCELPAR,  &
                              NOUT, MUOBS, PHIOBSREL, ADRADIANCE)
!    PRINT *, 'Lay Pres     AdTau   AdAlbedo  AdLegen1'
!    DO L = 1, NLAY
!       print '(1X,I3,1X,F6.1,1X,E10.3,1X,F8.3,1X,F8.3)', &
!            L, PRES(NZ-L+1), ADTAUP(L), ADALBEDOP(L), ADLEGENP(1,L)
!    ENDDO
!    print *, ' '
  ADTEMPSKIN = ADTEMPSKIN + ADSFCTEMP
  ADSFCALBEDO = ADSFCALBEDO + ADSFCPARMS(1)

    ! Call the adjoint optical properties routine
  CALL ADCALCULATE_OPTICAL_PROPERTIES (NLAY, NCOMP, ML, SPECIESFLAG, &
                       LWP, ADLWP, RMM, ADRMM, &
                       MOLECTAU, ADMOLECTAU, TAURAYL, ADTAURAYL, &
                       SRCTYPE, SOLARMU, NOUT, MUOBS, PHIOBSREL, &
                       TAUP, ADTAUP, ALBEDOP, ADALBEDOP, LEGENP, ADLEGENP, &
                       SSPHASE, ADSSPHASE, &
                   NSCATTAB, SPECIES2TABLE, MAXNR, &
                   NRTAB(:,ICHANSCAT), RTAB(:,:,ICHANSCAT), &
                   EXTINCTTAB(:,:,ICHANSCAT), SSALBTAB(:,:,ICHANSCAT), &
                   MAXLEG, NLEGTAB(:,:,ICHANSCAT), LEGENTAB(:,:,:,ICHANSCAT))

   ! Convert from adjoint integrated mass content and mean mass radius
   !  to adjoint particle mass mixing ratio and number concentration
  DO L = 2, NLAY
    LR = NLAY-L+1
    DO I = 1, NCOMP
      IF (SPECIES2TABLE(I) > 0) THEN
        MIXAVG = 0.5*(MIXR(I,LR)+MIXR(I,LR+1)) + MINMIXR
        RMM_MID = RTAB(NRTAB(I,1)/2+1,I,ICHANSCAT) + 0.1
        DCONC = MINMIXR/((4*3.14159*PARDENS(I)/3.0)*(RMM_MID/1.0E4)**3)
        CONCAVG = 0.5*(CONC(I,LR)+CONC(I,LR+1)) + DCONC
        ADMIXAVG = RMM(I,L)/(3.0*MIXAVG) *ADRMM(I,L)
        ADCONCAVG = -RMM(I,L)/(3.0*CONCAVG) *ADRMM(I,L)
        ADMIXAVG = ADMIXAVG + ADLWP(I,L)*(PRES(LR)-PRES(LR+1))*100/9.80
        ADMIXR(I,LR) = ADMIXR(I,LR) + 0.5*ADMIXAVG
        ADMIXR(I,LR+1) = ADMIXR(I,LR+1) + 0.5*ADMIXAVG
        ADCONC(I,LR) = ADCONC(I,LR) + 0.5*ADCONCAVG
        ADCONC(I,LR+1) = ADCONC(I,LR+1) + 0.5*ADCONCAVG
      ENDIF
    ENDDO
  ENDDO

   ! Call molecular absorption adjoint if channel index so indicates
  IF (INDEXABS > 0) THEN
    CALL ADJOINT_MOLEC_ABSORPTION (INDEXABS, NZ-1, PRESLEV, &
                                   PRES2, TEMP2, VAPMIXR2, O3MIXR2, &
                                   ADMOLECTAU(2:NLAY), &
                                   ADTEMP2, ADVAPMIXR2, ADO3MIXR2)
     ! Convert from top down arrays to bottom up arrays and 
     !   do adjoint of interpolation
    DO L = 1, NZ-1
      ADTEMP(NZ-L) = ADTEMP(NZ-L) + 0.5*ADTEMP2(L)
      ADTEMP(NZ-L+1) = ADTEMP(NZ-L+1) + 0.5*ADTEMP2(L)
      ADVAPMIXR(NZ-L) = ADVAPMIXR(NZ-L) + 0.5*ADVAPMIXR2(L)
      ADVAPMIXR(NZ-L+1) = ADVAPMIXR(NZ-L+1) + 0.5*ADVAPMIXR2(L)
      ADO3MIXR(NZ-L) = ADO3MIXR(NZ-L) + 0.5*ADO3MIXR2(L)
      ADO3MIXR(NZ-L+1) = ADO3MIXR(NZ-L+1) + 0.5*ADO3MIXR2(L)
    ENDDO
  ENDIF

   ! Reverse the temperature levels from SHDOMPP
  ADTEMP(1:NZ) = ADTEMP(1:NZ) + ADTEMPP(NLAY+1:2:-1)

  DEALLOCATE (ADTEMP2, ADVAPMIXR2, ADO3MIXR2)
  DEALLOCATE (ADTAUP, ADALBEDOP, ADLEGENP, ADSSPHASE)
  DEALLOCATE (ADTEMPP, ADMOLECTAU, ADTAURAYL, ADLWP, ADRMM)
END SUBROUTINE ADJOINT_COLUMN_RADIANCE





SUBROUTINE DEALLOC_COLUMN_RADIANCE (INDEXABS)
 ! Deallocates arrays allocated in FORWARD_COLUMN_RADIANCE
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: INDEXABS

  DEALLOCATE (IXP, FRACTAUG, SHSOURCE, ACCELPAR, RADIANCE)
  DEALLOCATE (TAUP, ALBEDOP, LEGENP, SSPHASE)
  DEALLOCATE (LWP, RMM, MOLECTAU, TAURAYL)
  DEALLOCATE (PRESLEV, PRES2, TEMP2, VAPMIXR2, O3MIXR2)
  DEALLOCATE (TEMPP, SPECIESFLAG)

  IF (INDEXABS > 0) THEN
    CALL dealloc_molec_abs_arrays
  ENDIF
END SUBROUTINE DEALLOC_COLUMN_RADIANCE



SUBROUTINE DEALLOC_SCATTERING_TABLES
 ! Deallocates arrays allocated in READ_SCATTERING_TABLES
  IMPLICIT NONE

  DEALLOCATE (RTAB, EXTINCTTAB, SSALBTAB, NLEGTAB, LEGENTAB)
  DEALLOCATE (SPECIES2TABLE, NRTAB, WAVELEN1, WAVELEN2, PARDENS)
END SUBROUTINE DEALLOC_SCATTERING_TABLES


END MODULE SHDOMPPDA_FWD_TL_ADJ

