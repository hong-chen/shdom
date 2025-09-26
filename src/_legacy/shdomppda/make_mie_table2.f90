PROGRAM MAKE_MIE_TABLE2
! 
! Does Mie computations to create a scattering table as a function of
! mean mass radius for gamma or lognormal size distributions of spherical
! particles.  The particles may be water or ice (in which case the 
! program provides the index of refraction depending on wavelength) or
! "aerosols" (in which case the index of refraction is user specified).
! For water or ice particles the scattering properties may be averaged 
! over the desired spectral range with Planck function weighting.  
! The phase functions in the output scattering table are represented 
! with Legendre series.
!
!  compile: pgf90 -fast -o make_mie_table2  make_mie_table2.f90 mieindsub.f
!    
!    Frank Evans    University of Colorado     March 2005
!
  IMPLICIT NONE
  INTEGER :: NRTAB
  REAL    :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
  REAL    :: SRTAB, ERTAB, ALPHA, MAXRADIUS
  COMPLEX :: RINDEX
  CHARACTER(LEN=1) :: PARTYPE, AVGFLAG, DISTFLAG
  CHARACTER(LEN=80) :: MIETABFILE
  INTEGER :: NSIZE, MAXLEG, I, J, L, NL
  REAL    :: PI, WAVELENCEN, XMAX, SCATTER
  INTEGER, ALLOCATABLE :: NLEG1(:), NLEG(:)
  REAL, ALLOCATABLE :: RADII(:), ND(:)
  REAL, ALLOCATABLE :: QEXT(:), QSCA(:)
  REAL, ALLOCATABLE :: RMM(:), REFF(:)
  REAL, ALLOCATABLE :: EXTINCT1(:), SCATTER1(:), LEGEN1(:,:)
  REAL, ALLOCATABLE :: EXTINCT(:), SSALB(:), LEGCOEF(:,:)


   ! Get the input parameters
  CALL USER_INPUT (WAVELEN1, WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                   AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, &
                   NRTAB, SRTAB, ERTAB, MAXRADIUS,  MIETABFILE)
  IF (PARTYPE == 'W') THEN
    PARDENS = 1.0
  ELSE IF (PARTYPE == 'I') THEN
    PARDENS = 0.916
  ENDIF
  
   ! Calculate the maximum size parameter and the max number of Legendre terms
  CALL GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
  PI = ACOS(-1.0)
  IF (AVGFLAG == 'A') THEN
    XMAX = 2*PI*MAXRADIUS/WAVELEN1
  ELSE
    XMAX = 2*PI*MAXRADIUS/WAVELENCEN
  ENDIF  
  MAXLEG = NINT(2*(XMAX + 4.0*XMAX**0.3334 + 2))

   ! Get the average index of refraction for water or ice
  IF (PARTYPE /= 'A') THEN
    CALL GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
  ENDIF

   ! Figure the number of radii there will be
  CALL GET_NSIZE (SRTAB, MAXRADIUS, WAVELENCEN, NSIZE)

   ! Allocate all the arrays here
  ALLOCATE (RADII(NSIZE), ND(NSIZE), NLEG1(NSIZE))
  ALLOCATE (EXTINCT1(NSIZE), SCATTER1(NSIZE), LEGEN1(0:MAXLEG,NSIZE))
  ALLOCATE (RMM(NRTAB), REFF(NRTAB))
  ALLOCATE (EXTINCT(NRTAB), SSALB(NRTAB))
  ALLOCATE (NLEG(NRTAB), LEGCOEF(0:MAXLEG,NRTAB))


   ! Make up the discrete particle radii to use
  CALL GET_SIZES (SRTAB, MAXRADIUS, WAVELENCEN, NSIZE, RADII)

   ! Do the Mie computations for each radius, which may involve several
   !   Mie calculation over the wavelength integration
  CALL COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, PARTYPE, &
                              WAVELENCEN, RINDEX, NSIZE, RADII, MAXLEG, &
                              EXTINCT1, SCATTER1, NLEG1, LEGEN1)


  ! Loop over the number of output tabulated radii
  DO I = 1, NRTAB
    ! Set tabulated radius
    IF (NRTAB <= 1) THEN
      RMM(I) = SRTAB
    ELSE
      RMM(I) = SRTAB* (ERTAB/SRTAB)** (FLOAT(I-1)/(NRTAB-1))
    ENDIF

    ! Calculate the discrete size number concentrations (ND), which vary
    !   according to a truncated gamma or lognormal distribution,
    !   that gives the desired mean mass radius (RMM) and LWC (1 g/m^3).
    !   Returns the actual RMM and effective radius (REFF).
    CALL MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, RMM(I), ALPHA, &
                         ND, REFF(I))


    ! Sum the scattering properties over the discrete size distribution
    EXTINCT(I) = 0.0
    SCATTER = 0.0
    LEGCOEF(:,I) = 0.0
    NL = 1
    DO J = 1, NSIZE
      EXTINCT(I) = EXTINCT(I) + ND(J)*EXTINCT1(J)
      SCATTER = SCATTER + ND(J)*SCATTER1(J)
      NL = MAX(NL,NLEG1(J))
      LEGCOEF(0:NL,I) = LEGCOEF(0:NL,I) + ND(J)*LEGEN1(0:NL,J)
    ENDDO
    DO L = 0, NL
      LEGCOEF(L,I) = LEGCOEF(L,I)/SCATTER
      IF (LEGCOEF(L,I) .GT. 0.5E-5) NLEG(I) = L
    ENDDO
    IF (ABS(LEGCOEF(0,I)-1.0) > 0.0001) THEN
      PRINT *,'Phase function not normalized for R=',RMM(I),LEGCOEF(0,I)
      STOP
    ENDIF 
    IF (EXTINCT(I) > 0.0) THEN
      SSALB(I) = SCATTER/EXTINCT(I)
    ENDIF
    EXTINCT(I) = 0.001*EXTINCT(I)

  ENDDO  ! end of radius loop


  CALL WRITE_MIE_TABLE (MIETABFILE, WAVELEN1, WAVELEN2, DELTAWAVE, &
                        PARTYPE, PARDENS, RINDEX, DISTFLAG, &
                        ALPHA, NRTAB, SRTAB, ERTAB, RMM, REFF, &
                        EXTINCT, SSALB, MAXLEG, NLEG, LEGCOEF)

  DEALLOCATE (RMM, REFF, EXTINCT, SSALB, NLEG, LEGCOEF)
  DEALLOCATE (ND, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
END




SUBROUTINE USER_INPUT (WAVELEN1,WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                       AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, &
                       NRTAB, SRTAB, ERTAB, MAXRADIUS,  MIETABFILE)
 ! Reads the input parameters from the standard input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NRTAB
  REAL,    INTENT(OUT) :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
  REAL,    INTENT(OUT) :: SRTAB, ERTAB, ALPHA, MAXRADIUS
  COMPLEX, INTENT(OUT) :: RINDEX
  CHARACTER(LEN=1), INTENT(OUT) :: PARTYPE, AVGFLAG, DISTFLAG
  CHARACTER(LEN=*), INTENT(OUT) :: MIETABFILE

  WRITE(*,*) 'Making Mie scattering tables for spherical particles'

  WRITE (*,*) 'Wavelength range (micron)'
  READ (*,*) WAVELEN1, WAVELEN2
    WRITE (*,'(2(1X,F6.3))') WAVELEN1, WAVELEN2
  IF (WAVELEN1 > WAVELEN2) STOP 'USER_INPUT: wavelength1 must be <= wavelength2'

  WRITE(*,*) 'Water, Ice, or Aerosol spherical particles (W,I,A)'
  READ(*,'(A1)') PARTYPE
    WRITE(*,'(1X,A1)') PARTYPE

  IF (PARTYPE == 'W' .OR. PARTYPE == 'I') THEN
    WRITE(*,*) 'Average Mie properties over wavelength or use Planck weighted'
    WRITE(*,*) '  index of refraction at center wavelength (A or C)'
    READ(*,'(A1)') AVGFLAG
    WRITE (*,'(1X,A1)') AVGFLAG   
    IF (AVGFLAG == 'A') THEN
      WRITE(*,*) 'Wavelength interval for averaging (micron)'
      READ(*,*) DELTAWAVE
      WRITE(*,'(1X,F6.3)') DELTAWAVE
    ELSE
      AVGFLAG = 'C'
      DELTAWAVE = -1.0
    ENDIF
  ELSE
    PARTYPE = 'A'
    WRITE (*,*) 'Aerosol complex index of refraction (negative imaginary part)'
    READ (*,*) RINDEX 
    WRITE (*,*) RINDEX
    WRITE (*,*) 'Aerosol particle bulk density (g/cm^3)'
    READ (*,*) PARDENS
    WRITE (*,'(1X,F5.3)') PARDENS
    AVGFLAG = 'C'
  ENDIF

  WRITE (*,*) 'Particle size distribution type: G = Gamma, L = Lognormal'
  READ (*,*) DISTFLAG
  WRITE (*,*) DISTFLAG
  IF (DISTFLAG == 'L') THEN
    WRITE (*,*) 'Log normal size distribution log standard deviation'
  ELSE
    DISTFLAG = 'G'
    WRITE(*,*) 'Gamma size distribution shape parameter (alpha)'
  ENDIF
  READ (*,*) ALPHA
  WRITE (*,'(1X,F6.3)') ALPHA

  WRITE(*,*) 'Number, starting, and ending tabulated radius (micron)'
  READ(*,*) NRTAB, SRTAB, ERTAB
    WRITE (*,'(1X,I3,2(1X,F7.2))') NRTAB, SRTAB, ERTAB

  WRITE(*,*) 'Maxium particle radius in size distribution (micron)'
  READ(*,*) MAXRADIUS
    WRITE (*,'(1X,F7.2)') MAXRADIUS

  WRITE (*,*) 'Output Mie scattering table name'
  READ (*,'(A)') MIETABFILE
    WRITE(*,'(1X,A70)') MIETABFILE
END SUBROUTINE USER_INPUT



SUBROUTINE GET_NSIZE (SRTAB, MAXRADIUS, WAVELEN, NSIZE)
 ! Calculates the number of radii for which the Mie computation will be run.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  REAL,    INTENT(IN)  :: SRTAB, MAXRADIUS, WAVELEN
  INTEGER, INTENT(OUT) :: NSIZE
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD

  TWOPI = 2.0*ACOS(-1.0)
  RADMIN = 0.02*SRTAB
  RAD = RADMIN
  NSIZE = 1
  DO WHILE (RAD < MAXRADIUS)
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = DELX*WAVELEN/TWOPI
    RAD = RAD + DELRAD
    NSIZE = NSIZE + 1
  ENDDO
END SUBROUTINE GET_NSIZE


SUBROUTINE GET_SIZES (SRTAB, MAXRADIUS, WAVELEN, NSIZE, RADII)
 ! Calculates the radii for which the Mie computation will be run and
 ! from which all the size distributions will be computed.
 ! The formula and spacing in size parameter can be changed to trade
 ! off size distribution integration accuracy vs. computer time.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: SRTAB, MAXRADIUS, WAVELEN
  REAL,    INTENT(OUT) :: RADII(NSIZE)
  INTEGER :: N
  REAL    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD

  TWOPI = 2.0*ACOS(-1.0)
  RAD = 0.02*SRTAB
  RADII(1) = RAD
  DO N = 2, NSIZE
    X = TWOPI*RAD/WAVELEN
    DELX = MAX(0.01,0.03*X**0.5)    ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
    DELRAD = DELX*WAVELEN/TWOPI
    RAD = RAD + DELRAD
    RADII(N) = RAD
  ENDDO
END SUBROUTINE GET_SIZES



SUBROUTINE GET_CENTER_WAVELEN (WAVELEN1, WAVELEN2, WAVELENCEN)
!  Returns the Planck weighted center wavelength averaged over the 
! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
! blackbody temperature of 5800 K is used for the Planck weighting
! if the average wavelength is less than 3 microns, no Planck weighting
! is done for an average wavelength between 3 and 5 microns, and a 
! blackbody temperature of 270 K is done for an average wavelength 
! greater than 5 microns.
  IMPLICIT NONE
  REAL, INTENT(IN)  :: WAVELEN1, WAVELEN2
  REAL, INTENT(OUT) :: WAVELENCEN
  REAL :: WAVECEN, DELWAVE, WAVE, BBTEMP, PLANCK, SUMP, SUMW

  IF (WAVELEN1 == WAVELEN2) THEN
    WAVELENCEN = WAVELEN1  
  ELSE
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF 
    DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
    SUMP = 0.0
    SUMW = 0.0
    WAVE = WAVELEN1   
    DO WHILE (WAVE <= WAVELEN2)
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      SUMW = SUMW + PLANCK*WAVE
      WAVE = WAVE + DELWAVE
    ENDDO
    WAVELENCEN = 0.001*NINT(1000*SUMW/SUMP)
  ENDIF
END SUBROUTINE GET_CENTER_WAVELEN




SUBROUTINE GET_REFRACT_INDEX (PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
 ! Returns the index of refraction for water or ice averaged over
 ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
 ! averaging is done at 0.05 micron intervals and is weighted by
 ! a Planck function at 5800 K temperature for central wavelengths
 ! less than 3 microns, a flat function between 3 and 5 microns, and
 ! 270 K Planck function beyond 5 microns.  
 ! The index of refraction is at -30 C for ice and +10 C for water
 ! (the temperature dependence is important in the microwave).
  IMPLICIT NONE
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE
  REAL, INTENT(IN) :: WAVELEN1, WAVELEN2
  COMPLEX, INTENT(OUT) :: RINDEX
  REAL :: WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
  REAL :: MRE, MIM, SUMP, SUMMR, SUMMI, A

  WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
  IF (WAVECEN < 3.0) THEN
    BBTEMP = 5800.0
  ELSE IF (WAVECEN > 5.0) THEN
    BBTEMP = 270.0
  ELSE
    BBTEMP = -1.0
    PLANCK = 1.0
  ENDIF 
  DELWAVE = MIN(WAVECEN/100.,0.1*ABS(WAVELEN2-WAVELEN1))
  DELWAVE = MAX(DELWAVE,WAVECEN*1.0E-5)
  SUMP = 0.0
  SUMMR = 0.0
  SUMMI = 0.0
  WAVE = WAVELEN1
  DO WHILE (WAVE <= WAVELEN2)
    IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
    SUMP = SUMP + PLANCK
    IF (PARTYPE == 'I') THEN
      CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
    ELSE
      CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
    ENDIF
    SUMMR = SUMMR + PLANCK*MRE
    SUMMI = SUMMI + PLANCK*MIM
    WAVE = WAVE + DELWAVE
  ENDDO
  MRE = SUMMR/SUMP
  MIM = SUMMI/SUMP
  RINDEX = CMPLX(MRE,-MIM)
END SUBROUTINE GET_REFRACT_INDEX


 

SUBROUTINE COMPUTE_MIE_ALL_SIZES (AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                  PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                  MAXLEG, EXTINCT1, SCATTER1, NLEG1, LEGEN1)
 ! Does a Mie computation for each particle radius in RADII and returns the
 ! optical properties in arrays EXTINCT1, SCATTER1, NLEG1, and LEGEN1.
 ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
 ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
 ! integration of the Mie properties over wavelength is performed for
 ! each radius.  For each wavelength, with spacing DELTAWAVE, the water 
 ! or ice (depending on PARTYPE) index of refraction is obtained and
 ! used in the Mie computation for that wavelength, and the Mie optical
 ! properties are averaged with Planck function weighting (blackbody
 ! temperature depends on wavelength).  The Legendre coefficients are
 ! returned with the product of the phase function times the scattering
 ! coefficient.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE, MAXLEG
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
  REAL,    INTENT(IN) :: RADII(NSIZE)
  COMPLEX, INTENT(IN) :: RINDEX
  CHARACTER(LEN=1), INTENT(IN) :: AVGFLAG, PARTYPE
  INTEGER, INTENT(OUT) :: NLEG1(NSIZE)
  REAL,    INTENT(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
  REAL,    INTENT(OUT) :: LEGEN1(0:MAXLEG,NSIZE)
  INTEGER :: I, NL
  REAL    :: WAVECEN, WAVE, BBTEMP, PLANCK, SUMP, A
  REAL    :: MRE, MIM, EXT, SCAT, LEG(0:MAXLEG)
  COMPLEX :: REFIND

  IF (AVGFLAG == 'C') THEN
     ! For using one central wavelength: just call Mie routine for each radius
    DO I = 1, NSIZE
      CALL MIE_ONE (WAVELENCEN, RINDEX, RADII(I), MAXLEG, &
                    EXTINCT1(I), SCATTER1(I), NLEG1(I), LEGEN1(0,I) )
    ENDDO

  ELSE
     ! For averaging over wavelength range:
    WAVECEN = 0.5*(WAVELEN1+WAVELEN2)
    IF (WAVECEN < 3.0) THEN
      BBTEMP = 5800.0
    ELSE IF (WAVECEN > 5.0) THEN
      BBTEMP = 270.0
    ELSE
      BBTEMP = -1.0
      PLANCK = 1.0
    ENDIF 
    EXTINCT1(:) = 0.0
    SCATTER1(:) = 0.0
    NLEG1(:) = 1
    LEGEN1(:,:) = 0.0
    SUMP = 0.0
    WAVE = WAVELEN1
    DO WHILE (WAVE <= WAVELEN2)   ! Loop over the wavelengths
      IF (BBTEMP > 0) PLANCK = (1.19E8/WAVE**5)/(EXP(1.439E4/(WAVE*BBTEMP))-1)
      SUMP = SUMP + PLANCK
      IF (PARTYPE == 'I') THEN   ! Get the index of refraction of water or ice
        CALL REFICE (0, WAVE, 243.0, MRE, MIM, A, A)
      ELSE
        CALL REFWAT (0, WAVE, 283.0, MRE, MIM, A, A)
      ENDIF
      REFIND = CMPLX(MRE,-MIM)
      DO I = 1, NSIZE
        CALL MIE_ONE (WAVE, REFIND, RADII(I), MAXLEG, EXT, SCAT, NL, LEG)
        EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
        SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
        NLEG1(I) = MAX(NLEG1(I),NL)
        LEGEN1(0:NL,I) = LEGEN1(0:NL,I) + PLANCK*LEG(0:NL)
      ENDDO
      WAVE = WAVE + DELTAWAVE
    ENDDO
    EXTINCT1(:) = EXTINCT1(:)/SUMP
    SCATTER1(:) = SCATTER1(:)/SUMP
    LEGEN1(:,:) = LEGEN1(:,:)/SUMP
  ENDIF
END SUBROUTINE COMPUTE_MIE_ALL_SIZES




SUBROUTINE MAKE_SIZE_DIST (DISTFLAG, PARDENS, NSIZE, RADII, RMM, ALPHA, &
                           ND, REFF)
 ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
 ! discrete particle radii (micron) of a gamma or lognormal size distribution
 ! with a mean mass radius of RMM (micron), gamma shape parameter or 
 ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
 ! Outputs the actual mean mass radius RMM and effective radius REFF (micron).
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NSIZE
  REAL,    INTENT(IN)  :: RADII(NSIZE), ALPHA, PARDENS
  REAL,    INTENT(INOUT)  :: RMM
  REAL,    INTENT(OUT) :: ND(NSIZE), REFF
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL, PARAMETER :: TOL=0.001  ! fractional tolerance in achieving Rmm
  INTEGER :: I
  REAL    :: REHI, RELO, REMID, RE, RM

  ! Find Reff that gives Rmm above desired value
  REHI = RMM
  I = 0
  RM = 0.5*RMM
  DO WHILE (RM <= RMM .AND. I < 4)
    REHI = 2*REHI
    I = I + 1
    CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, REHI, NSIZE,RADII, ND, RE, RM)
  ENDDO
  IF (RM <= RMM) THEN
    PRINT *, 'MAKE_SIZE_DIST: mean mass radius cannot be achieved',RMM,RM
    STOP
  ENDIF

  ! Find Reff that gives Rmm below desired value
  RELO = RMM
  I = 0
  RM = 2*RMM
  DO WHILE (RM >= RMM .AND. I < 4)
    RELO = 0.5*RELO
    I = I + 1
    CALL DO_SIZE_DIST (PARDENS,DISTFLAG, ALPHA, RELO, NSIZE,RADII, ND, RE, RM)
  ENDDO
  IF (RM >= RMM) THEN
    PRINT *, 'MAKE_SIZE_DIST: mean mass radius cannot be achieved',RMM, RM
    STOP
  ENDIF

  ! Do bisection to get correct mean mass radius
  DO WHILE (ABS(RM-RMM) > TOL*RMM)
    REMID = 0.5*(RELO+REHI)
    CALL DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REMID, NSIZE,RADII, ND, RE, RM)
    IF (RM < RMM) THEN
      RELO = REMID
    ELSE
      REHI = REMID
    ENDIF
  ENDDO
  RMM = RM
  REFF = RE
END SUBROUTINE MAKE_SIZE_DIST



SUBROUTINE DO_SIZE_DIST (PARDENS, DISTFLAG, ALPHA, REFF, NSIZE, RADII, &
                         ND, RE, RM)
 ! For the input effective radius (REFF) [um], returns the number 
 ! concentrations ND [cm^-3] and the calculated effective radius RE [um] and 
 ! mean mass radius RM [um] for a gamma or lognormal size distribution 
 ! with mass content of 1 g/m^3.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NSIZE
  REAL,    INTENT(IN) :: PARDENS, ALPHA, REFF, RADII(NSIZE)
  CHARACTER(LEN=1), INTENT(IN) :: DISTFLAG
  REAL,    INTENT(OUT) :: ND(NSIZE), RE, RM
  INTEGER :: J
  REAL    :: DENS, PI, A, B, LWC, R, DELR, SUM2, SUM3, NT, GAMMLN

  PI = ACOS(-1.0)
  IF (DISTFLAG == 'G') THEN
    B = (ALPHA+3)/REFF
    A = 1.E6/( (4*PI/3.)*PARDENS *B**(-ALPHA-4) *EXP(GAMMLN(ALPHA+4.)) )
  ELSE
    B = REFF*EXP(-2.5*ALPHA**2)
    A = 1.E6/( (4*PI/3.)*PARDENS *SQRT(2*PI)*ALPHA * B**3 *EXP(4.5*ALPHA**2) )
  ENDIF
  LWC = 0.0
  SUM2 = 0.0
  SUM3 = 0.0
  DO J = 1, NSIZE
    R = RADII(J)
    DELR = SQRT(RADII(J)*RADII(MIN(NSIZE,J+1))) &
         - SQRT(RADII(J)*RADII(MAX(1,J-1)))
    IF (DISTFLAG == 'G') THEN
      ND(J) = A* R**ALPHA *EXP(-B*R) *DELR
    ELSE
      ND(J) = (A/R)*EXP(-0.5*(LOG(R/B))**2/ALPHA**2) *DELR
    ENDIF
    LWC = LWC + 1.0E-6*PARDENS*ND(J)*(4*PI/3)*R**3
    SUM2 = SUM2 + ND(J)*R**2
    SUM3 = SUM3 + ND(J)*R**3
  ENDDO
  ND(:) = (1.0/LWC)*ND(:)
  RE = SUM3/SUM2
  NT = SUM(ND(:))    ! total number concentration
  RM = 100* ( 3*1.0 / (4*PI*PARDENS*NT) ) **0.33333333  !  r_mm in um
END SUBROUTINE DO_SIZE_DIST




SUBROUTINE WRITE_MIE_TABLE (MIETABFILE, WAVELEN1, WAVELEN2, DELTAWAVE, &
                            PARTYPE, PARDENS, RINDEX, DISTFLAG, &
                            ALPHA, NRTAB, SRTAB, ERTAB, RMM, REFF, &
                            EXTINCT, SSALB, MAXLEG, NLEG, LEGCOEF)
 ! Writes the table of Mie scattering properties as a function of radius.  
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NRTAB, MAXLEG, NLEG(NRTAB)
  COMPLEX, INTENT(IN) :: RINDEX
  REAL,    INTENT(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
  REAL,    INTENT(IN) :: PARDENS, SRTAB, ERTAB, ALPHA
  REAL,    INTENT(IN) :: RMM(NRTAB), REFF(NRTAB)
  REAL,    INTENT(IN) :: EXTINCT(NRTAB), SSALB(NRTAB)
  REAL,    INTENT(IN) :: LEGCOEF(0:MAXLEG,NRTAB)
  CHARACTER(LEN=1), INTENT(IN) :: PARTYPE, DISTFLAG
  CHARACTER(LEN=*), INTENT(IN) :: MIETABFILE
  INTEGER :: I, J, L, NL

  OPEN (UNIT=3, FILE=MIETABFILE, STATUS='REPLACE')
  WRITE (3,'(A)') '! Mie scattering table vs. mean mass radius (LWC=1 g/m^3)'
  IF (DELTAWAVE < 0.0) THEN
    WRITE (3,'(2(1X,F8.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
  ELSE
    WRITE (3,'(3(1X,F8.3),A)') WAVELEN1, WAVELEN2, DELTAWAVE, '  wavelength range and averaging step (micron)'
  ENDIF
  WRITE (3,'(1X,F5.3,2X,A1,A)') PARDENS, PARTYPE, '   particle density (g/cm^3) and type (Water, Ice, Aerosol)'
  WRITE (3,'(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
  IF (DISTFLAG == 'L') THEN
    WRITE (3,'(F7.5,A)') ALPHA, '  lognormal log standard deviation'
  ELSE
    WRITE (3,'(F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
  ENDIF
  WRITE (3,'(1X,I3,2(1X,F8.3),A)') NRTAB, SRTAB, ERTAB, &
        '  number, starting, ending mean mass radius with log spacing'

  DO I = 1, NRTAB
    WRITE (3,'(2(1X,F8.4),1X,E12.5,1X,F8.6,1X,I6,A)') &
        RMM(I), REFF(I), EXTINCT(I), SSALB(I), NLEG(I),&
        '   Rmm  Reff  Ext  Alb  Nleg'
    WRITE (3,'(2X,201(1X,F10.5))') (LEGCOEF(L,I), L=0,MIN(NLEG(I),200))
    DO J = 200, NLEG(I)-1, 200
      WRITE (3,'(2X,200(1X,F10.5))') (LEGCOEF(J+L,I),L=1,MIN(200,NLEG(I)-J))
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE WRITE_MIE_TABLE

