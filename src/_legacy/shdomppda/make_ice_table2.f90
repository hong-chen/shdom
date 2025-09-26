PROGRAM MAKE_ICE_TABLE2

!   Reads the 0.25 degree gaussian forward peak smoothed version of 
! Ping Yang's 2005 ice crystal scattering database, and creates a
! scattering table for gamma size distributions of a sequence of mean
! mass radii.  Another input file specifies the mixing fraction for
! each of the six particle shapes as a function of the maximum diameter.
! The format of the mixing file:
!    one header line
!    Dmax  F1 F2 F3 F4 F5 F6    (Fn are the fractions for each shape)
!
! The scattering properties are averaged over the desired spectral range 
! using linear interpolation between wavelengths.  The phase functions in 
! the output scattering table are represented with Legendre series.
!
!  compile: pgf90 -fast -o make_ice_table2  make_ice_table2.f90
!
!    Frank Evans    University of Colorado     November 2005

  IMPLICIT NONE
  INTEGER :: NRTAB
  REAL    :: WAVELEN1, WAVELEN2, SRTAB, ERTAB, ALPHA
  CHARACTER(LEN=80) :: ICESCATDB, SHAPEMIXFILE, ICETABFILE
  INTEGER, PARAMETER :: NSHAPE=6   ! number of shapes in ice scattering db
  INTEGER, PARAMETER :: NSIZE=45   ! number of sizes in ice scattering db
  INTEGER, PARAMETER :: NANG=288   ! number of phase function angles in db
  INTEGER, PARAMETER :: NQUAD=2500 ! number of quad angles for Legendre conversion
  INTEGER, PARAMETER :: MAXLEG=NQUAD ! number of Legendre coeffs
  INTEGER :: I, J
  REAL    :: SCATTER, AREA, PI
  INTEGER, ALLOCATABLE :: NLEG(:)
  REAL, ALLOCATABLE :: DMAX(:), DIAMAREA(:,:), DIAMVOL(:,:), ND(:)
  REAL, ALLOCATABLE :: QEXT(:,:), QSCA(:,:)
  REAL, ALLOCATABLE :: RMM(:), REFF(:)
  REAL, ALLOCATABLE :: ANGLES(:), PHASEFUNC(:,:,:), PHASEDIST(:)
  REAL, ALLOCATABLE :: EXTINCT(:), SSALB(:), LEGCOEF(:,:)
  REAL(8), ALLOCATABLE :: QUADMU(:), QUADWTS(:)


  CALL USER_INPUT (WAVELEN1, WAVELEN2, ICESCATDB, SHAPEMIXFILE, &
                   NRTAB, SRTAB, ERTAB, ALPHA, ICETABFILE)


   ! Allocate all the arrays here
  ALLOCATE (QUADMU(NQUAD), QUADWTS(NQUAD))
  ALLOCATE (DMAX(NSIZE), DIAMAREA(NSIZE,NSHAPE), DIAMVOL(NSIZE,NSHAPE))
  ALLOCATE (QEXT(NSIZE,NSHAPE), QSCA(NSIZE,NSHAPE))
  ALLOCATE (ANGLES(NANG), PHASEFUNC(NANG,NSIZE,NSHAPE))
  ALLOCATE (ND(NSIZE), PHASEDIST(NANG))
  ALLOCATE (RMM(NRTAB), REFF(NRTAB))
  ALLOCATE (EXTINCT(NRTAB), SSALB(NRTAB))
  ALLOCATE (NLEG(NRTAB), LEGCOEF(0:MAXLEG,NRTAB))

   ! Make the Gauss-Legendre quadrature abscissas and weights
  CALL GAUSQUAD (NQUAD, QUADMU, QUADWTS)


   ! Read in the ice scattering properties averaged over the wavelength
   !   range for all shapes.  The diameters are in microns.
  CALL READ_ICE_SCAT_DB (ICESCATDB, WAVELEN1, WAVELEN2, &
                         NSHAPE, NSIZE, NANG, &
                         DMAX, DIAMAREA, DIAMVOL, &
                         QEXT, QSCA, ANGLES, PHASEFUNC)

   ! Average the properties over the shapes according to the shape mixture file
  CALL MIX_SHAPES (SHAPEMIXFILE, NSHAPE, NSIZE, NANG, &
                   DMAX, DIAMAREA, DIAMVOL, QEXT, QSCA, PHASEFUNC)

!  IF (SRTAB < 2.0*0.5*DIAMVOL(1) .OR. ERTAB > 0.5*0.5*DIAMVOL(NSIZE)) THEN
!    PRINT *, 'MAKE_ICE_TABLE: effective radius range outside possible range'
!    STOP
!  ENDIF

  PI = ACOS(-1.0)  
  ! Loop over the number of output tabulated mean mass radii
  DO I = 1, NRTAB
     ! Set tabulated mean mass radius
    IF (NRTAB <= 1) THEN
      RMM(I) = SRTAB
    ELSE
      RMM(I) = SRTAB* (ERTAB/SRTAB)** (FLOAT(I-1)/(NRTAB-1))
    ENDIF

    ! Calculate the discrete size number concentrations (ND), which vary
    !   according to a gamma distribution in volume equivalent diameter,
    !   that give the desired mean mass radius (RMM) and IWC (1 g/m^3).
    !   Returns the actual RMM and effective radius (REFF).
    CALL MAKE_SIZE_DIST (NSIZE, DIAMAREA(:,1), DIAMVOL(:,1), &
                         RMM(I), ALPHA, ND, REFF(I))

    ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    PHASEDIST(:) = 0.0
    DO J = 1, NSIZE
      AREA = PI*(0.5*DIAMAREA(J,1))**2 
      EXTINCT(I) = EXTINCT(I) + QEXT(J,1)*AREA*ND(J)
      SCATTER = SCATTER + QSCA(J,1)*AREA*ND(J)
      PHASEDIST(:) = PHASEDIST(:) + QSCA(J,1)*AREA*ND(J)*PHASEFUNC(:,J,1)
    ENDDO
    PHASEDIST(:) = PHASEDIST(:)/SCATTER
    SSALB(I) = SCATTER/EXTINCT(I)
    EXTINCT(I) = 0.001*EXTINCT(I)

    ! Convert phase function from angle to Legendre coefficients
    CALL CONVERT_LEGENDRE (NANG, ANGLES, PHASEDIST, NQUAD, QUADMU, QUADWTS, &
                           MAXLEG, NLEG(I), LEGCOEF(:,I))
  ENDDO  ! end of mean mass radius loop


  CALL WRITE_SCAT_TABLE (ICETABFILE, WAVELEN1, WAVELEN2, SHAPEMIXFILE, &
                         NRTAB, SRTAB, ERTAB, ALPHA, RMM, REFF, &
                         EXTINCT, SSALB, MAXLEG, NLEG, LEGCOEF)

  DEALLOCATE (RMM, REFF, DMAX,DIAMAREA, DIAMVOL, QEXT, QSCA, ANGLES, PHASEFUNC)
  DEALLOCATE (ND, PHASEDIST, EXTINCT, SSALB, NLEG, LEGCOEF)
END





SUBROUTINE USER_INPUT (WAVELEN1, WAVELEN2, ICESCATDB, SHAPEMIXFILE, &
                       NRTAB, SRTAB, ERTAB, ALPHA, ICETABFILE)
 ! Reads the input parameters from the standard input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NRTAB
  REAL,    INTENT(OUT) :: WAVELEN1, WAVELEN2, SRTAB, ERTAB, ALPHA
  CHARACTER(LEN=*), INTENT(OUT) :: ICESCATDB, SHAPEMIXFILE, ICETABFILE

  WRITE (*,*) 'Input ice scattering database name'
  READ (*,'(A)') ICESCATDB
    WRITE(*,'(1X,A70)') ICESCATDB

  WRITE (*,*) 'Wavelength range (micron)'
  READ (*,*) WAVELEN1, WAVELEN2
    WRITE (*,'(2(1X,F6.3))') WAVELEN1, WAVELEN2
  IF (WAVELEN1 > WAVELEN2) STOP 'USER_INPUT: wavelength1 must be <= wavelength2'

  WRITE (*,*) 'Input ice particle shape mixing file name'
  READ (*,'(A)') SHAPEMIXFILE
    WRITE(*,'(1X,A70)') SHAPEMIXFILE

  WRITE(*,*) 'Number, starting, and ending tabulated mean mass radius'
  READ(*,*) NRTAB, SRTAB, ERTAB
    WRITE (*,'(1X,I3,2(1X,F7.2))') NRTAB, SRTAB, ERTAB

  WRITE(*,*) 'Gamma size distribution shape parameter (alpha)'
  READ (*,*) ALPHA
    WRITE (*,'(1X,F6.3)') ALPHA

  WRITE (*,*) 'Output scattering table name'
  READ (*,'(A)') ICETABFILE
    WRITE(*,'(1X,A70)') ICETABFILE
END SUBROUTINE USER_INPUT





SUBROUTINE READ_ICE_SCAT_DB (ICESCATDB, WAVELEN1, WAVELEN2, &
                             NSHAPE, NSIZE, NANG, &
                             DMAX, DIAMAREA, DIAMVOL, &
                             QEXT, QSCA, ANGLES, PHASEFUNC)
 ! Reads the 0.25 degree Gaussian forward peak smoothed version of 
 ! Ping Yang's 2005 ice crystal scattering database.  The scattering
 ! properties are linearly interpolated and averaged over the desired 
 ! spectral range with uniform spectral weighting. 
 ! Inputs:
 !   ICESCATDB  database filename
 !   WAVELEN1, WAVELEN2  wavelength range (micron)
 !   NSHAPE     number of particle shapes (6)
 !   NSIZE      number of particle sizes (45)
 !   NANG       number of phase function angles (288)
 ! Outputs: 
 !   DMAX(:)       maximum diameter (microns)
 !   DIAMAREA(:,:) equivalent area spherical diameter (microns)
 !   DIAMVOL(:,;)  equivalent volume spherical diameter (microns)
 !   QEXT(:,:)     extinction efficiency (2 in geometric optics limit)
 !   QSCA(:,:)     scattering efficiency
 !   ANGLES(:)     angles of phase functions (degrees)
 !   PHASEFUNC(:,:,:) phase functions
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: ICESCATDB
  INTEGER, INTENT(IN)  :: NSIZE, NANG, NSHAPE
  REAL,    INTENT(IN)  :: WAVELEN1, WAVELEN2
  REAL,    INTENT(OUT) :: DMAX(NSIZE)
  REAL,    INTENT(OUT) :: DIAMAREA(NSIZE,NSHAPE), DIAMVOL(NSIZE,NSHAPE)
  REAL,    INTENT(OUT) :: QEXT(NSIZE,NSHAPE), QSCA(NSIZE,NSHAPE)
  REAL,    INTENT(OUT) :: ANGLES(NANG), PHASEFUNC(NANG,NSIZE,NSHAPE)
  INTEGER, PARAMETER :: NWAVE=65
  INTEGER :: I, J, K, M, N, ISH, D1, D2, D3
  REAL    :: DELWAVE, F1, F2, WAVE, EXT, SSALB, PHASE
  REAL, ALLOCATABLE :: WAVEDB(:), WTWAVE(:)
  CHARACTER(LEN=288*3) :: PHASESTRING

  QEXT(:,:) = 0.0
  QSCA(:,:) = 0.0
  PHASEFUNC(:,:,:) = 0.0

  OPEN (UNIT=1, FILE=ICESCATDB, STATUS='OLD')  
  ! Read useful parts of header: wavelengths and phase function angles
  DO I = 1, 9
    READ (1,*)
  ENDDO
  READ (1,*) N, (ANGLES(I), I=1,N)
  IF (N /= NANG) STOP 'READ_ICE_SCAT_DB: wrong ice scattering database'
  READ (1,*)
  ALLOCATE (WAVEDB(NWAVE), WTWAVE(NWAVE))
  READ (1,*) (WAVEDB(I), I=1,NWAVE)
  READ (1,*)

  IF (WAVELEN1 < WAVEDB(1) .OR. WAVELEN2 > WAVEDB(NWAVE)) &
    STOP 'READ_ICE_SCAT_DB: Wavelengths out of range of ice scattering database'

  ! Make the weights for each database wavelength
  DELWAVE = WAVELEN2-WAVELEN1
  WTWAVE(:) = 0.0
  DO K = 1, NWAVE-1
    F1 = (WAVELEN1-WAVEDB(K))/(WAVEDB(K+1)-WAVEDB(K))
    F2 = (WAVELEN2-WAVEDB(K))/(WAVEDB(K+1)-WAVEDB(K))
    IF (F1 > 1.0 .OR. F2 < 0.0) THEN
    ELSE IF (F1 >= 0.0 .AND. F2 <= 1.0) THEN
      WTWAVE(K) = 0.5*(2.0-F1-F2)
      WTWAVE(K+1) = 0.5*(F1+F2)
    ELSE IF (F1 >= 0.0 .AND. F1 <= 1.0 .AND. F2 > 1.0) THEN
      WTWAVE(K) = 0.5*(1-F1) *(WAVEDB(K+1)-WAVELEN1)/DELWAVE
      WTWAVE(K+1) = 0.5*(1+F1) *(WAVEDB(K+1)-WAVELEN1)/DELWAVE
    ELSE IF (F1 < 0.0 .AND. F2 > 1.0) THEN
      WTWAVE(K) = WTWAVE(K) + 0.5 *(WAVEDB(K+1)-WAVEDB(K))/DELWAVE
      WTWAVE(K+1) = 0.5 *(WAVEDB(K+1)-WAVEDB(K))/DELWAVE
    ELSE IF (F1 < 0.0 .AND. F2 >= 0.0 .AND. F2 <= 1.0) THEN
      WTWAVE(K) = WTWAVE(K) + 0.5*(2-F2) *(WAVELEN2-WAVEDB(K))/DELWAVE
      WTWAVE(K+1) = 0.5*F2 *(WAVELEN2-WAVEDB(K))/DELWAVE
    ENDIF
  ENDDO

  print *
  print *,' Wavelen  Weight  <- summing database wavelengths'
  do k = 1, nwave
    if (wtwave(k) > 0) then
      print '(1X,F7.3,2X,F7.4)', wavedb(k), wtwave(k)
    endif
  enddo

  DO I = 1, NSHAPE
    ! If this is a wavelength we need then read the data in and do the
    !   weighted average of single scattering properties.
    DO K = 1, NWAVE
      IF (WTWAVE(K) > 0.0) THEN
        DO J = 1, NSIZE
          READ (1,'(I2,1X,F8.4,3(1X,F8.3),2(1X,F8.6),10X,A)') &
              ISH, WAVE, DMAX(J), DIAMAREA(J,I), DIAMVOL(J,I), &
                   EXT, SSALB, PHASESTRING
          IF (WAVE /= WAVEDB(K)) &
            STOP 'READ_ICE_SCAT_DB: Error reading scattering database'
          IF (I /= ISH) &
            STOP 'READ_ICE_SCAT_DB: Error reading scattering database'
          QEXT(J,I) = QEXT(J,I) + WTWAVE(K)*EXT
          QSCA(J,I) = QSCA(J,I) + WTWAVE(K)*EXT*SSALB
          DO M = 1, NANG
            D1 = ICHAR(PHASESTRING(3*M-2:3*M-2))-32
            D2 = ICHAR(PHASESTRING(3*M-1:3*M-1))-32
            D3 = ICHAR(PHASESTRING(3*M:3*M))-32
            PHASE = 10**(10.D0*(D1/95.D0+D2/95.D0**2+D3/95.D0**3) - 4.D0)
            PHASEFUNC(M,J,I) = PHASEFUNC(M,J,I) + WTWAVE(K)*EXT*SSALB*PHASE
          ENDDO
        ENDDO
      ELSE  ! Otherwise skip this wavelength
        DO J = 1, NSIZE
          READ (1,*)
        ENDDO
      ENDIF
    ENDDO

    DO J = 1, NSIZE
      PHASEFUNC(:,J,I) = PHASEFUNC(:,J,I)/QSCA(J,I)
    ENDDO
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_ICE_SCAT_DB




SUBROUTINE MIX_SHAPES (SHAPEMIXFILE, NSHAPE, NSIZE, NANG, &
                       DMAX, DIAMAREA, DIAMVOL, QEXT, QSCA, PHASEFUNC)
 ! Reads the shape mixture file, which specifies the mixing fraction for
 ! each of the six shapes in the ice crystal scattering database as a
 ! function of maximum diameter. 
 ! The format of the mixing file (Dmax must increase):
 !    one header line
 !    Dmax  F1 F2 F3 F4 F5 F6    (Fn are the fractions for each shape)
 ! The mixing fractions are interpolated in Dmax between entries in the
 ! mixing table.  The mixing fractions of the smallest Dmax in the table
 ! are used for smaller Dmax database values, and fractions for the largest
 ! Dmax in the table are used for larger Dmax database values/
 ! The mixture properties are returned in the first shape index.  
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SHAPEMIXFILE
  INTEGER, INTENT(IN) :: NSHAPE, NSIZE, NANG
  REAL,    INTENT(IN) :: DMAX(NSIZE)
  REAL,    INTENT(INOUT) :: DIAMAREA(NSIZE,NSHAPE), DIAMVOL(NSIZE,NSHAPE)
  REAL,    INTENT(INOUT) :: QEXT(NSIZE,NSHAPE), QSCA(NSIZE,NSHAPE)
  REAL,    INTENT(INOUT) :: PHASEFUNC(NANG,NSIZE,NSHAPE)
  INTEGER :: N, I, J, M, IERR
  REAL    :: F, FMIX(6), TOTAREA, TOTSCA
  REAL, ALLOCATABLE :: DMIX(:), FRAC(:,:)

   ! Read the mixing file to determine the number of entries
  OPEN (UNIT=1, FILE=SHAPEMIXFILE, STATUS='OLD')
  READ (1,*)
  N = 0
  DO WHILE (IERR == 0)
    N = N + 1
    READ (1,*,IOSTAT=IERR) 
  ENDDO
  CLOSE (1)
  N = N - 1
  IF (N < 2) STOP 'MIX_SHAPES: Must be at least two entries in shape mixing file'

   ! Allocate the mixing fraction arrays and read in from the file
  ALLOCATE (DMIX(N), FRAC(6,N))
  OPEN (UNIT=1, FILE=SHAPEMIXFILE, STATUS='OLD')
  READ (1,*)
  DO I = 1, N
    READ (1,*) DMIX(I), FRAC(:,I)
  ENDDO
  CLOSE (1)

   ! Loop over the database maximum diameters
  DO J = 1, NSIZE
     ! Interpolate between Dmax entries in the shape mixing table to
     !   get the FMIX mixing fractions for this Dmax(j)
    I = 1
    DO WHILE (I<N .AND. DMIX(I) < DMAX(J))
      I = I + 1
    ENDDO
    I = MAX(I-1,1)
    F = (DMAX(J)-DMIX(I))/(DMIX(I+1)-DMIX(I))
    IF (F < 0.0) THEN
      FMIX(:) = FRAC(:,I)
    ELSE IF (F > 1.0) THEN
      FMIX(:) = FRAC(:,I+1)
    ELSE
      FMIX(:) = (1-F)*FRAC(:,I) + F*FRAC(:,I+1)
    ENDIF
    FMIX(:) = FMIX(:)/SUM(FMIX(:))  ! Force normalization of the mixing fracs

     ! Average extensive properties over the shapes
    TOTAREA = SUM(FMIX(:) * DIAMAREA(J,:)**2)
    TOTSCA = SUM(FMIX(:) * QSCA(J,:) * DIAMAREA(J,:)**2)
    DO M = 1, NANG
      PHASEFUNC(M,J,1) = &
          SUM(FMIX(:)*PHASEFUNC(M,J,:)*QSCA(J,:)*DIAMAREA(J,:)**2)/TOTSCA
    ENDDO
    QSCA(J,1) = TOTSCA/TOTAREA
    QEXT(J,1) = SUM(FMIX(:) * QEXT(J,:) * DIAMAREA(J,:)**2) / TOTAREA
    DIAMAREA(J,1) = SQRT(TOTAREA)
    DIAMVOL(J,1) = ( SUM(FMIX(:)*DIAMVOL(J,:)**3) )**(1.0/3.0)
  ENDDO
  DEALLOCATE (DMIX, FRAC)
END SUBROUTINE MIX_SHAPES



SUBROUTINE MAKE_SIZE_DIST (NSIZE, DIAMAREA, DIAMVOL, RMM, ALPHA, ND, REFF)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle sizes of a gamma size distribution (in volume
 ! spherical equivalent diameter DIAMVOL) with a distribution 
 ! mean mass radius  of RMM, gamma shape parameter ALPHA, and ice water 
 ! content of 1 g/m^3.  The effective radius (0.75*volume/projected area) 
 ! is also output in REFF.
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE
  REAL,    INTENT(IN)  :: DIAMAREA(NSIZE), DIAMVOL(NSIZE), ALPHA
  REAL,    INTENT(INOUT) :: RMM
  REAL,    INTENT(OUT) :: ND(NSIZE), REFF
  REAL, PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Rmm
  INTEGER :: I
  REAL    :: DEHI, DELO, DEMID, RM, RE


  ! Find Deff that gives Rmm above desired value
  DEHI = 2.0*RMM
  I = 0
  RM = 0.5*RMM
  DO WHILE (RM <= RMM .AND. I < 4)
    DEHI = 2.0*DEHI
    I = I + 1
    CALL DO_GAMMA_DIST (NSIZE, DEHI, ALPHA, DIAMAREA, DIAMVOL, ND, RE, RM)
  ENDDO
  IF (RM <= RMM) THEN
    PRINT *, 'MAKE_SIZE_DIST: mean mass radius cannot be achieved',RMM,RM
    STOP
  ENDIF

  ! Find Deff that gives Rmm below desired value
  DELO = 2.0*RMM
  I = 0
  RM = 2*RMM
  DO WHILE (RM >= RMM .AND. I < 4)
    DELO = 0.5*DELO
    I = I + 1
    CALL DO_GAMMA_DIST (NSIZE, DELO, ALPHA, DIAMAREA, DIAMVOL, ND, RE, RM)
  ENDDO
  IF (RM >= RMM) THEN
    PRINT *, 'MAKE_SIZE_DIST: mean mass radius cannot be achieved',RMM,RM
    STOP
  ENDIF

  ! Do bisection to get correct mean mass radius
  DO WHILE (ABS(RM-RMM) > TOL*RMM)
    DEMID = 0.5*(DELO+DEHI)
    CALL DO_GAMMA_DIST (NSIZE, DEMID, ALPHA, DIAMAREA, DIAMVOL, ND, RE, RM)
    IF (RM < RMM) THEN
      DELO = DEMID
    ELSE
      DEHI = DEMID
    ENDIF
  ENDDO
  RMM = RM
  REFF = RE
END SUBROUTINE MAKE_SIZE_DIST



SUBROUTINE DO_GAMMA_DIST (NSIZE, DEFF, ALPHA, DIAMAREA, DIAMVOL, ND, RE, RM)
 ! For the input effective diameter (DEFF) and ALPHA, returns the 
 ! number concentrations ND [cm^-3] and the calculated effective radius 
 ! RE [um] and mean mass radius RM [um] for a gamma size distribution 
 ! in DIAMVOL with total IWC of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: DEFF, ALPHA, DIAMAREA(NSIZE), DIAMVOL(NSIZE)
  REAL,    INTENT(OUT) :: ND(NSIZE), RE, RM
  REAL, PARAMETER :: DENSICE = 0.916
  INTEGER :: J
  REAL    :: PI, B, D, DELD, IWC, SUM2, SUM3, NT

  PI = ACOS(-1.0)
  B = (ALPHA+3)/DEFF
  IWC = 0.0
  SUM2 = 0.0
  SUM3 = 0.0
  DO J = 1, NSIZE
    D = DIAMVOL(J)
    DELD = SQRT(DIAMVOL(J)*DIAMVOL(MIN(NSIZE,J+1))) &
         - SQRT(DIAMVOL(J)*DIAMVOL(MAX(1,J-1)))
    ND(J) = D**ALPHA *EXP(-B*D) * DELD
    IWC = IWC + 1.0E-6*DENSICE*ND(J)*(PI/6)*DIAMVOL(J)**3
    SUM2 = SUM2 + ND(J)*DIAMAREA(J)**2
    SUM3 = SUM3 + ND(J)*DIAMVOL(J)**3
  ENDDO
  ND(:) = (1.0/IWC)*ND(:)
  RE = 0.5*SUM3/SUM2
  NT = SUM(ND(:))    ! total number concentration
  RM = 100* ( 3*1.0 / (4*PI*DENSICE*NT) ) **0.33333333  !  r_mm in um
END SUBROUTINE DO_GAMMA_DIST




SUBROUTINE CONVERT_LEGENDRE (NANG, ANGLES, PHASE, NQUAD, QUADMU, QUADWTS, &
                             MAXLEG, NLEG, LEGCOEF)
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NANG, NQUAD, MAXLEG
  REAL,    INTENT(IN)  :: ANGLES(NANG), PHASE(NANG)
  REAL(8), INTENT(IN)  :: QUADMU(NQUAD), QUADWTS(NQUAD)
  INTEGER, INTENT(OUT) :: NLEG
  REAL,    INTENT(OUT) :: LEGCOEF(0:MAXLEG)
  INTEGER :: N, L
  REAL(8) :: INTEG, PL, PL1, PL2
  REAL(8), ALLOCATABLE :: PHASEQUAD(:), LEGEN(:)

  ALLOCATE (PHASEQUAD(NQUAD), LEGEN(0:MAXLEG))

  ! Interpolate phase function to the gaussian quadrature points
  CALL SPLINE_INTERP_PHASE (NANG, ANGLES, PHASE, &
                            NQUAD, QUADMU, PHASEQUAD, QUADWTS, INTEG)

  ! Do final small normalization by multiplication
  PHASEQUAD(:) = (2.D0/INTEG)*PHASEQUAD(:)

  ! Compute the Legendre coefficients for the smoothed phase func
  LEGEN(:) = 0.0
  DO N = 1, NQUAD
    ! Use upward recurrence to find Legendre polynomials
    PL1 = 1.0
    PL = 1.0 
    DO L = 0, MAXLEG
      IF (L .GT. 0) PL = (2*L-1)*QUADMU(N)*PL1/L-(L-1)*PL2/L
      LEGEN(L)=LEGEN(L) + PL*QUADWTS(N)*PHASEQUAD(N)
      PL2 = PL1
      PL1 = PL 
    ENDDO
  ENDDO
  ! Find last significant Legendre coefficient
  DO L = 0, MAXLEG
    LEGCOEF(L) = 0.5*(2*L+1)*LEGEN(L)
    IF (LEGCOEF(L) .GT. 1.0E-5) THEN   
      NLEG = L
    ENDIF
  ENDDO  

  DEALLOCATE (PHASEQUAD, LEGEN)
END SUBROUTINE CONVERT_LEGENDRE



SUBROUTINE SPLINE_INTERP_PHASE (NANG, ANGLES, PHASE, &
                                NQUAD, QUADMU, PHASEQUAD, QUADWTS, INTEG)
 ! Interpolates the phase function from ANGLES (degrees) to 
 ! quadrature points (in cos theta).  Also computes integral.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NANG, NQUAD
  REAL,    INTENT(IN) :: ANGLES(NANG), PHASE(NANG)
  REAL(8), INTENT(IN) :: QUADMU(NQUAD), QUADWTS(NQUAD)
  REAL(8), INTENT(OUT) :: PHASEQUAD(NQUAD), INTEG
  INTEGER :: I, J
  REAL(8) :: RD, SEVAL
  REAL(8), ALLOCATABLE :: X(:), Y(:), B(:), C(:), D(:)

  ALLOCATE (X(NANG), Y(NANG), B(NANG), C(NANG), D(NANG))

  RD = ACOS(-1.0D0)/180.D0
  DO I = 1, NANG
     X(I) = COS(RD*ANGLES(NANG+1-I))
     Y(I) = PHASE(NANG+1-I)
  ENDDO
  CALL SPLINE (NANG, X, Y, B, C, D)

  INTEG = 0.0D0
  DO J = 1, NQUAD
    PHASEQUAD(J) = SEVAL (NANG, QUADMU(J), X, Y, B, C, D)
    INTEG = INTEG + QUADWTS(J)*PHASEQUAD(J)
  ENDDO
END SUBROUTINE SPLINE_INTERP_PHASE




SUBROUTINE WRITE_SCAT_TABLE (ICETABFILE, WAVELEN1, WAVELEN2, SHAPEMIXFILE, &
                             NRTAB, SRTAB, ERTAB, ALPHA, RMM, REFF, &
                             EXTINCT, SSALB, MAXLEG, NLEG, LEGCOEF)
 ! Writes the table of ice scattering properties as a function of radius.  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NRTAB, MAXLEG, NLEG(NRTAB)
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, SRTAB, ERTAB, ALPHA
  REAL,    INTENT(IN) :: RMM(NRTAB), REFF(NRTAB)
  REAL,    INTENT(IN) :: EXTINCT(NRTAB), SSALB(NRTAB)
  REAL,    INTENT(IN) :: LEGCOEF(0:MAXLEG,NRTAB)
  CHARACTER(LEN=*), INTENT(IN) :: ICETABFILE, SHAPEMIXFILE
  INTEGER :: I, J, L, NL

  OPEN (UNIT=3, FILE=ICETABFILE, STATUS='REPLACE')
  WRITE (3,'(A)') '! Ice scattering table vs. mean mass radius (IWC=1 g/m^3)'
  WRITE (3,'(2(1X,F6.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
  WRITE (3,'(1X,F5.3,A)') 0.916, '  particle density (g/cm^3) for ice'
  WRITE (3,'(A,A55)') '! Shape mixture file: ', SHAPEMIXFILE
  WRITE (3,'(1X,F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
  WRITE (3,'(1X,I3,2(1X,F8.3),A)') NRTAB, SRTAB, ERTAB, &
        '  number, starting, ending mean mass radius with log spacing'

  DO I = 1, NRTAB
    WRITE (3,'(2(1X,F8.4),1X,E12.5,1X,F8.6,1X,I6,A)') &
        RMM(I), REFF(I), EXTINCT(I), SSALB(I), NLEG(I), &
        '   Rmm  Reff  Ext  Alb  Nleg'
    WRITE (3,'(2X,201(1X,F10.5))') (LEGCOEF(L,I), L=0,MIN(NLEG(I),200))
    DO J = 200, NLEG(I)-1, 200
      WRITE (3,'(2X,200(1X,F10.5))') (LEGCOEF(J+L,I),L=1,MIN(200,NLEG(I)-J))
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE WRITE_SCAT_TABLE



SUBROUTINE GAUSQUAD (N, XA, WT)
 ! Generates the abscissas (X) and weights (W) for an N point
 ! Gauss-Legendre quadrature.  
  IMPLICIT NONE
  INTEGER :: N
  REAL(8) :: XA(N), WT(N)
  INTEGER :: K, I, J, L
  REAL(8) :: X, XP, PL, PL1, PL2, DPL

  K = (N+1)/2
  DO J = 1, K
    X = COS(3.141592654*(J-.25)/(N+.5))
    I = 0
    DO WHILE (I < 10)
      PL1 = 1
      PL = X
      DO L = 2, N
         PL2 = PL1
         PL1 = PL
         PL = ( (2*L-1)*X*PL1 - (L-1)*PL2 )/L
       ENDDO
       DPL = N*(X*PL-PL1)/(X*X-1)
       XP = X
       X = XP - PL/DPL
       I = I+1
       IF (ABS(X-XP) < 2*EPSILON(X)) EXIT
    ENDDO
    XA(J)     = -X
    XA(N-J+1) = X
    WT(J  )   = 2.0D0/((1.0D0-X*X)*DPL*DPL)
    WT(N-J+1) = WT(J)
  ENDDO
END SUBROUTINE GAUSQUAD



subroutine spline (n, x, y, b, c, d)
  implicit none
  integer n
  real(8) x(n), y(n), b(n), c(n), d(n)

!  the coefficients b(i), c(i), and d(i), i=1,2,...,n are computed
!  for a cubic interpolating spline
!
!    s(x) = y(i) + b(i)*(x-x(i)) + c(i)*(x-x(i))**2 + d(i)*(x-x(i))**3
!
!    for  x(i) .le. x .le. x(i+1)
!
!  input..
!
!    n = the number of data points or knots (n.ge.2)
!    x = the abscissas of the knots in strictly increasing order
!    y = the ordinates of the knots
!
!  output..
!
!    b, c, d  = arrays of spline coefficients as defined above.
!
!  using  p  to denote differentiation,
!
!    y(i) = s(x(i))
!    b(i) = sp(x(i))
!    c(i) = spp(x(i))/2
!    d(i) = sppp(x(i))/6  (derivative from the right)
!
!  the accompanying function subprogram  seval  can be used
!  to evaluate the spline.
!
!
      integer nm1, ib, i
      double precision t
!
      nm1 = n-1
      if ( n .lt. 2 ) return
      if ( n .lt. 3 ) go to 50
!
!  set up tridiagonal system
!
!  b = diagonal, d = offdiagonal, c = right hand side.
!
      d(1) = x(2) - x(1)
      c(2) = (y(2) - y(1))/d(1)
      do 10 i = 2, nm1
         d(i) = x(i+1) - x(i)
         b(i) = 2.*(d(i-1) + d(i))
         c(i+1) = (y(i+1) - y(i))/d(i)
         c(i) = c(i+1) - c(i)
   10 continue
!
!  end conditions.  third derivatives at  x(1)  and  x(n)
!  obtained from divided differences
!
      b(1) = -d(1)
      b(n) = -d(n-1)
      c(1) = 0.
      c(n) = 0.
      if ( n .eq. 3 ) go to 15
      c(1) = c(3)/(x(4)-x(2)) - c(2)/(x(3)-x(1))
      c(n) = c(n-1)/(x(n)-x(n-2)) - c(n-2)/(x(n-1)-x(n-3))
      c(1) = c(1)*d(1)**2/(x(4)-x(1))
      c(n) = -c(n)*d(n-1)**2/(x(n)-x(n-3))
!
!  forward elimination
!
   15 do 20 i = 2, n
         t = d(i-1)/b(i-1)
         b(i) = b(i) - t*d(i-1)
         c(i) = c(i) - t*c(i-1)
   20 continue
!
!  back substitution
!
      c(n) = c(n)/b(n)
      do 30 ib = 1, nm1
         i = n-ib
         c(i) = (c(i) - d(i)*c(i+1))/b(i)
   30 continue
!
!  c(i) is now the sigma(i) of the text
!
!  compute polynomial coefficients
!
      b(n) = (y(n) - y(nm1))/d(nm1) + d(nm1)*(c(nm1) + 2.*c(n))
      do 40 i = 1, nm1
         b(i) = (y(i+1) - y(i))/d(i) - d(i)*(c(i+1) + 2.*c(i))
         d(i) = (c(i+1) - c(i))/d(i)
         c(i) = 3.*c(i)
   40 continue
      c(n) = 3.*c(n)
      d(n) = d(n-1)
      return
!
   50 b(1) = (y(2)-y(1))/(x(2)-x(1))
      c(1) = 0.
      d(1) = 0.
      b(2) = b(1)
      c(2) = 0.
      d(2) = 0.
end subroutine spline



double precision function seval(n, u, x, y, b, c, d)
  implicit none
  integer n
  real(8)  u, x(n), y(n), b(n), c(n), d(n)
!  this subroutine evaluates the cubic spline function
!
!    seval = y(i) + b(i)*(u-x(i)) + c(i)*(u-x(i))**2 + d(i)*(u-x(i))**3
!
!    where  x(i) .lt. u .lt. x(i+1), using horner's rule
!
!  if  u .lt. x(1) then  i = 1  is used.
!  if  u .ge. x(n) then  i = n  is used.
!
!  input..
!
!    n = the number of data points
!    u = the abscissa at which the spline is to be evaluated
!    x,y = the arrays of data abscissas and ordinates
!    b,c,d = arrays of spline coefficients computed by spline
!
!  if  u  is not in the same interval as the previous call, then a
!  binary search is performed to determine the proper interval.
!
      integer i, j, k
      real*8 dx
      data i/1/
      if ( i .ge. n ) i = 1
      if ( u .lt. x(i) ) go to 10
      if ( u .le. x(i+1) ) go to 30
!
!  binary search
!
   10 i = 1
      j = n+1
   20 k = (i+j)/2
      if ( u .lt. x(k) ) j = k
      if ( u .ge. x(k) ) i = k
      if ( j .gt. i+1 ) go to 20
!
!  evaluate spline
!
   30 dx = u - x(i)
      seval = y(i) + dx*(b(i) + dx*(c(i) + dx*d(i)))
      return
end function seval


