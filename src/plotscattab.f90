PROGRAM PLOTSCATTAB
! PLOTSCATTAB make a plotting file of the phase functions in a scattering 
! table file produced by make_mie_table.f90 or make_tmatrix_table.f90.  
! The scattering tables contains the Legendre polynomial expansion of the
! phase functions, or the Wigner d-function explansions of the six unique
! elements of the phase matrix, for several effective radii.  Selected 
! phase matrix elements for specified effective radii may be converted 
! to a function of angle for output.  The program can also plot the phase 
! functions or phase matrix elements in the beginning of an SHDOM tabulated 
! phase function property file.
!
!  compile: pgf90 -fast -o plotscattab  plotscattab.f90
! 
!    Frank Evans    University of Colorado     May 2003
!               Updated for polarized SHDOM    Feb 2014
  IMPLICIT NONE
  INTEGER :: NANGLE, NOUT, NPHASE, NPHASEPOL, MAXNCOEF
  INTEGER :: I, J, K, L
  LOGICAL :: POLTAB, FOUND
  REAL    :: PI, PL, PL1, PL2, SCAT, MU
  INTEGER, ALLOCATABLE :: NCOEF(:), PELEM(:)
  REAL,    ALLOCATABLE :: ANGLE(:), OUTREFF(:), REFF(:)
  REAL,    ALLOCATABLE :: WIGCOEF(:,:,:), PHASE(:,:)
  CHARACTER(LEN=1)  :: FILETYPE
  CHARACTER(LEN=72) :: INFILE, PLOTFILE
  CHARACTER(LEN=7)  :: ELEMNAMES(-6:6)

  ELEMNAMES(:) = (/ 'P34/P11', 'P12/P11', 'P44/P11', 'P33/P11', 'P22/P11', &
                    '  P11  ', '  P11  ', '  P11  ', &
                    '  P22  ', '  P33  ', '  P44  ', '  P12  ', '  P34  ' /)
  
  WRITE (*,*) 'Scattering table or SHDOM property file input (S or P)'
  READ (*,'(A1)') FILETYPE
    WRITE (*,'(1X,A1)') FILETYPE

  WRITE (*,*) 'Input file name'
  READ (*,'(A)') INFILE
    WRITE (*,'(A)') INFILE

  WRITE (*,'(1X,A)') 'Number of angles (<0 to input |n| angles)'
  READ (*,*) NANGLE
    WRITE (*,*) NANGLE
  ALLOCATE (ANGLE(ABS(NANGLE)))
  IF (NANGLE <= 0) THEN
    NANGLE = ABS(NANGLE)
    WRITE (*,*) 'Input the angles (degrees)'
    DO J = 1, NANGLE
      WRITE (*,'(1X,I3,A)') J, ' : '
      READ (*,*) ANGLE(J)
    ENDDO
  ELSE
    ANGLE(:) = (/ (180.0*FLOAT(J-1)/(NANGLE-1), J=1,NANGLE) /)
  ENDIF


  IF (FILETYPE == 'P') THEN
    WRITE (*,'(1X,A)') 'Number of phase functions to output'
    READ (*,*) NOUT
      WRITE (*,*) NOUT
    ALLOCATE (OUTREFF(NOUT))
    WRITE (*,*) 'Input the tabulated phase function indices'
    READ (*,*) OUTREFF(:)
      WRITE (*,'(20(1X,I4))') NINT(OUTREFF(:))
    WRITE (*,*) 'Phase matrix elements to output for each phase function index'
  ELSE
    WRITE (*,'(1X,A)') 'Number of effective radii/matrix elements to output'
    READ (*,*) NOUT
      WRITE (*,*) NOUT
    ALLOCATE (OUTREFF(NOUT))
    WRITE (*,*) 'Input the effective radii (micron)'
    READ (*,*) OUTREFF(:)
      WRITE (*,'(20(1X,F6.2))') OUTREFF(:)
    WRITE (*,*) 'Phase matrix elements to output for each effective radius'
  ENDIF
  ALLOCATE (PELEM(NOUT))
  WRITE (*,*) '(1=P11, 2=P22, 3=P33, 4=P44, 5=P12, 6=P34) (<0 for ratio to P11)'
  READ (*,*) PELEM(:)
    WRITE (*,'(20(1X,I2))') PELEM(:)

  WRITE (*,*) 'Plotting output file name'
  READ (*,'(A)') PLOTFILE
    WRITE (*,'(A)') PLOTFILE

  PELEM(:) = MIN(6,MAX(-6,PELEM(:)))
  WHERE (PELEM(:) == -1 .OR. PELEM(:) == 0)
    PELEM(:) = 1
  END WHERE

  IF (FILETYPE == 'P') THEN
    CALL READ_PROPERTY_SIZE (INFILE, NPHASEPOL, NPHASE, MAXNCOEF)
    ALLOCATE (REFF(NPHASE), NCOEF(NPHASE), WIGCOEF(NPHASEPOL,0:MAXNCOEF,NPHASE))
    CALL READ_PROPERTY_PHASEMAT (INFILE, NPHASEPOL, NPHASE, MAXNCOEF, &
                                 NCOEF, WIGCOEF)
    REFF(:) = (/ (FLOAT(J), J=1,NPHASE) /)
  ELSE
    CALL READ_SCAT_TABLE_SIZE (INFILE, POLTAB, NPHASEPOL, NPHASE, MAXNCOEF)
    ALLOCATE (REFF(NPHASE), NCOEF(NPHASE), WIGCOEF(NPHASEPOL,0:MAXNCOEF,NPHASE))
    CALL READ_SCAT_TABLE (INFILE, POLTAB, NPHASEPOL, NPHASE, REFF, &
                          MAXNCOEF, NCOEF, WIGCOEF)
  ENDIF

   ! Loop over each effective radius in scattering table (I), pulling out
   !   ones we want to output (K)
  ALLOCATE (PHASE(NANGLE,NOUT))
  PI = ACOS(-1.0)
  DO K = 1, NOUT
    FOUND = .FALSE.
    DO I = 1, NPHASE
      IF (REFF(I) == OUTREFF(K)) THEN
        FOUND = .TRUE.
        CALL TRANSFORM_WIGNERD_TO_PHASE (MAXNCOEF, NPHASEPOL, &
                                  PELEM(K), NCOEF(I), WIGCOEF(:,:,I), &
                                  NANGLE, ANGLE, PHASE(:,K))
      ENDIF
    ENDDO
    IF (.NOT. FOUND) THEN
      PRINT *, 'Phase function not found in input file:',OUTREFF(K)
      PHASE(:,K) = 1.0
    ENDIF
  ENDDO


   ! Output the phase functions
  OPEN (UNIT=1, FILE=PLOTFILE, STATUS='UNKNOWN')
  IF (FILETYPE == 'P') THEN
    WRITE (1,'(A,A40)') '!  Property file phase functions: ',INFILE
    WRITE (1,'(A)')'! Angle  cos(angle)  Phase functions with index numbers'
    WRITE (1,'(A,20(6X,I3,3X))') '!               ', (NINT(OUTREFF(K)), K=1,NOUT)
  ELSE
    WRITE (1,'(A,A40)') '!  Scattering table phase functions: ',INFILE
    WRITE (1,'(A)')'! Angle  cos(angle)  Phase functions for effective radii (um) / matrix element'
    WRITE (1,'(A,20(6X,F6.2))') '!               ', (OUTREFF(K), K=1,NOUT)
  ENDIF
  WRITE (1,'(A,20(5X,A7))') '!                ', (ELEMNAMES(PELEM(K)), K=1,NOUT)
  DO J = 1, NANGLE
      WRITE (1,'(1X,F7.2,1X,F9.6,20(E12.4))') &
          ANGLE(J), COS(ANGLE(J)*PI/180.0), (PHASE(J,K), K=1,NOUT)
  ENDDO
  CLOSE (1)

  DEALLOCATE (PHASE, REFF, NCOEF, WIGCOEF, OUTREFF, ANGLE, PELEM)
END





SUBROUTINE READ_SCAT_TABLE_SIZE (SCATTABFILE, POLTAB, NPHASEPOL, &
                                 NRETAB, MAXNCOEF)
 ! Reads the scattering table file to find the number of effective radii
 ! and the maximum order of the Wigner d-function series phase functions.
 ! Also returns whether the scattering table is polarized in POLTAB.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABFILE
  INTEGER, INTENT(OUT) :: NPHASEPOL, NRETAB, MAXNCOEF
  LOGICAL, INTENT(OUT) :: POLTAB
  INTEGER :: I, J, L, NCOEF, IERR
  REAL    :: REFF, EXT, SSALB, CHI

  OPEN (UNIT=3, FILE=SCATTABFILE, STATUS='OLD')
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*) NRETAB
   ! See if this is a polarized scattering table
  READ (3,*) REFF, EXT, SSALB, NCOEF
  POLTAB = .TRUE.
  DO I = 1, 6
    READ (3,*,IOSTAT=IERR) J, (CHI, L=0,NCOEF)
    IF (IERR>0) THEN
      POLTAB=.FALSE.
      EXIT
    ENDIF
    IF (J /= I) POLTAB=.FALSE.
  ENDDO
  REWIND (3)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*)
  READ (3,*) NRETAB
   ! Read the rest of the table to get the maximum Ncoef
  MAXNCOEF = 1
  DO I = 1, NRETAB
    READ (3,*) REFF, EXT, SSALB, NCOEF
    MAXNCOEF = MAX(MAXNCOEF,NCOEF)
    IF (POLTAB) THEN
      READ (3,*) J, (CHI, L=0,NCOEF)
      READ (3,*) J, (CHI, L=0,NCOEF)
      READ (3,*) J, (CHI, L=0,NCOEF)
      READ (3,*) J, (CHI, L=0,NCOEF)
      READ (3,*) J, (CHI, L=0,NCOEF)
      READ (3,*) J, (CHI, L=0,NCOEF)
    ELSE
      READ (3,*) (CHI, L=0,NCOEF)
    ENDIF
  ENDDO
  CLOSE (3)
  IF (POLTAB) THEN
    NPHASEPOL = 6
  ELSE
    NPHASEPOL = 1
  ENDIF
END SUBROUTINE READ_SCAT_TABLE_SIZE



SUBROUTINE READ_SCAT_TABLE (SCATTABFILE, POLTAB, NPHASEPOL, NRETAB,RETAB, &
                             MAXNCOEF, NCOEF, WIGCOEF)
 ! Reads the table of scattering properties as a function of effective 
 ! radius.  POLTAB=.true. means a polarized scattering table with Wigner 
 ! d-function coefficients for six phase matrix elements instead of one.
 ! The order of the six elements is the four diagonal ones (a1, a2, a3, a4)
 ! followed by the IQ and UV coefficients (b1, b2).
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABFILE
  INTEGER, INTENT(IN)  :: NPHASEPOL, NRETAB, MAXNCOEF
  LOGICAL, INTENT(IN)  :: POLTAB
  INTEGER, INTENT(OUT) :: NCOEF(NRETAB)
  REAL,    INTENT(OUT) :: RETAB(NRETAB)
  REAL,    INTENT(OUT) :: WIGCOEF(NPHASEPOL,0:MAXNCOEF,NRETAB)
  INTEGER :: I, J, L
  REAL    :: CHI, EXTINCT, SSALB

  OPEN (UNIT=1, FILE=SCATTABFILE, STATUS='OLD')
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*) 
  DO I = 1, NRETAB
    READ (1,*) RETAB(I), EXTINCT, SSALB, NCOEF(I)
    IF (POLTAB) THEN
      READ (1,*) J, (WIGCOEF(1,L,I), L=0,NCOEF(I))
      IF (NPHASEPOL > 1) THEN
        READ (1,*) J, (WIGCOEF(2,L,I), L=0,NCOEF(I))
        READ (1,*) J, (WIGCOEF(3,L,I), L=0,NCOEF(I))
        READ (1,*) J, (WIGCOEF(4,L,I), L=0,NCOEF(I))
        READ (1,*) J, (WIGCOEF(5,L,I), L=0,NCOEF(I))
        READ (1,*) J, (WIGCOEF(6,L,I), L=0,NCOEF(I))
      ELSE
        READ (1,*) J, (CHI, L=0,NCOEF(I))
        READ (1,*) J, (CHI, L=0,NCOEF(I))
        READ (1,*) J, (CHI, L=0,NCOEF(I))
        READ (1,*) J, (CHI, L=0,NCOEF(I))
        READ (1,*) J, (CHI, L=0,NCOEF(I))
      ENDIF
    ELSE
      READ (1,*) (WIGCOEF(1,L,I), L=0,NCOEF(I))
      IF (NPHASEPOL > 1) THEN
        WIGCOEF(2,:,I) = 0.0
        WIGCOEF(3,:,I) = 0.0
        WIGCOEF(4,:,I) = 0.0
        WIGCOEF(5,:,I) = 0.0
        WIGCOEF(6,:,I) = 0.0
      ENDIF
    ENDIF
    IF (ABS(WIGCOEF(1,0,I)-1.0) > 0.0001) THEN
      PRINT *, 'READ_SCAT_TABLE: Incorrect Legendre series; chi_0 is not 1'
      STOP
    ENDIF
    IF (I > 1 .AND. RETAB(I) <= RETAB(I-1)) THEN
      PRINT *,'READ_SCAT_TABLE: Effective radius not increasing in table:',SCATTABFILE
      STOP
    ENDIF
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_SCAT_TABLE



 
SUBROUTINE READ_PROPERTY_SIZE (PROPFILE, NPHASEPOL, NUMPHASE, MAXNCOEF)
 ! Reads the header of the tabulated phase function property file to 
 ! get the maximum array sizes needed for allocatable arrays.  
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER, INTENT(OUT) :: NPHASEPOL, NUMPHASE, MAXNCOEF
  INTEGER :: NPX, NPY, NPZ, NUML, I, J, K, L
  REAL    :: DELX, DELY, ZLEVELS, CHI
  CHARACTER(LEN=1) :: PROPTYPE
 
  ! Open the file, figure out the type, and get the grid size
  OPEN (UNIT=2, FILE=PROPFILE, STATUS='OLD')
  READ (2,'(A1)') PROPTYPE
  IF (PROPTYPE /= 'T' .AND. PROPTYPE /= 'P') THEN
    PRINT *, 'Must be a tabulated phase function property file.'
    STOP
  ENDIF
  READ (2,*) NPX, NPY, NPZ
  READ (2,*) DELX, DELY, (ZLEVELS, K=1,NPZ)
  READ (2,*) NUMPHASE
  IF (PROPTYPE .EQ. 'T') THEN
!     Property file type T is for tabulated phase function format
    NPHASEPOL = 1
    DO I = 1, NUMPHASE
      READ (2,*) NUML, (CHI, L=1,NUML)
      MAXNCOEF = MAX(NUML,MAXNCOEF)
    ENDDO
  ELSE IF (PROPTYPE .EQ. 'P') THEN
!    Property file type P is for polarized tabulated phase function format
    NPHASEPOL = 6
    DO I = 1, NUMPHASE 
      DO J = 1, 6
        READ (2,*) K, NUML, (CHI, L=0,NUML)
        IF (K /= J) THEN
          PRINT *, 'Error reading  polarized tabulated phase function format'
          STOP
        ENDIF
        MAXNCOEF = MAX(NUML,MAXNCOEF)   
      ENDDO
    ENDDO  
  ENDIF
  CLOSE (2)
END SUBROUTINE READ_PROPERTY_SIZE



SUBROUTINE READ_PROPERTY_PHASEMAT (PROPFILE, NPHASEPOL, NUMPHASE, MAXNCOEF,&
                                   NCOEF, WIGCOEF)
 ! Reads the Wigner d-function phase matrix coefficients from the header of
 ! a tabulated or polarized tabulated SHDOM optical property file.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER, INTENT(IN)  :: NPHASEPOL, NUMPHASE, MAXNCOEF
  INTEGER, INTENT(OUT) :: NCOEF(NUMPHASE)
  REAL,    INTENT(OUT) :: WIGCOEF(NPHASEPOL,0:MAXNCOEF,NUMPHASE)
  INTEGER :: NPX, NPY, NPZ, I, J, K, L
  REAL    :: DELX, DELY, ZLEVELS
  CHARACTER(LEN=1) :: PROPTYPE
 
  OPEN (UNIT=2, FILE=PROPFILE, STATUS='OLD')
  READ (2,'(A1)') PROPTYPE
  READ (2,*) NPX, NPY, NPZ
  READ (2,*) DELX, DELY, (ZLEVELS, K=1,NPZ)
  READ (2,*)
  DO I = 1, NUMPHASE
    IF (PROPTYPE .EQ. 'T') THEN
      READ (2,*) NCOEF(I), (WIGCOEF(1,L,I), L=1,NCOEF(I))
      WIGCOEF(1,0,I) = 1.0
    ELSE
      DO J = 1, 6
        READ (2,*) K, NCOEF(I), (WIGCOEF(K,L,I), L=0,NCOEF(I))
      ENDDO
    ENDIF  
  ENDDO
  CLOSE (2)
END SUBROUTINE READ_PROPERTY_PHASEMAT
 



SUBROUTINE TRANSFORM_WIGNERD_TO_PHASE (MAXNCOEF, NPHASEPOL, &
                                       PELEM, NCOEF, WIGCOEF, &
                                       NANGLE, ANGLE, PHASE)
 ! Transforms the phase matrix element (PELEM=1 to 6) from the Wigner 
 ! d-function based coefficients of the scattering matrix (WIGCOEF) to 
 ! a function of angle, PHASE(NANGLE) (at scattering angles in degrees
 ! in ANGLE(:)).  The order of the six elements in WIGCOEF is the four 
 ! diagonal ones (alpha1, alpha2, alpha3, alpha4) followed by the IQ and
 ! UV coefficients (beta1, beta2).  The phase matrix elements indexed by
 ! PELEM are P11, P22, P33, P44, P12, P34.  If PELEM<0 then the ABS(PELEM) 
 ! phase matrix element is normalized by the P11 element on output in PHASE.
 ! (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: MAXNCOEF, NPHASEPOL, PELEM, NCOEF, NANGLE
  REAL,    INTENT(IN) :: WIGCOEF(NPHASEPOL,0:MAXNCOEF), ANGLE(NANGLE)
  REAL,    INTENT(OUT) :: PHASE(NANGLE)
  INTEGER :: J, N
  DOUBLE PRECISION :: PI, X, XM, A1, A2, A3, A4, B1, B2, A2P3, A2M3
  DOUBLE PRECISION, ALLOCATABLE :: D00(:), D22P(:), D22M(:), D20(:)
  REAL, ALLOCATABLE :: NORM(:)

  IF (NPHASEPOL == 1 .AND. PELEM > 1) THEN
    PRINT '(A,I1,A)', 'Cannot do phase matrix element ',PELEM, &
       ' for unpolarized scattering table'
    STOP
  ENDIF

   ! Allocate arrays for Wigner functions D20, D22, D2-2 and D00
  ALLOCATE (D00(0:NCOEF), D22P(0:NCOEF), D22M(0:NCOEF), D20(0:NCOEF))
  ALLOCATE (NORM(NANGLE))

  PI = DACOS(-1.0D0)
  DO J = 1, NANGLE
    X = DCOS(ANGLE(J)*PI/180)
    XM = -X

     ! Calculate the desired scattering matrix element for this angle (X) 
     ! from the Wigner d-function coefficients
    SELECT CASE (ABS(PELEM))
    CASE (1)
      CALL WIGNERFCT (X, NCOEF, 0, 0, D00)
      A1 = 0.0D0
      DO N = 0, NCOEF
        A1 = A1 + WIGCOEF(1,N) * D00(N)
      ENDDO
      PHASE(J) = A1

    CASE (2:3)
      CALL WIGNERFCT (X, NCOEF, 2, 2, D22P)
      CALL WIGNERFCT (XM, NCOEF, 2, 2, D22M) ! multiply by (-1)**N
      A2P3 = 0.0D0
      A2M3 = 0.0D0
      DO N = 2, NCOEF
        A2P3 = A2P3 + (WIGCOEF(2,N) + WIGCOEF(3,N)) * D22P(N)
        A2M3 = A2M3 + (WIGCOEF(2,N) - WIGCOEF(3,N)) *(-1)**N * D22M(N)
      ENDDO
      A2 = 0.5D0 *(A2P3 + A2M3)
      A3 = 0.5D0 *(A2P3 - A2M3)
      IF (ABS(PELEM) == 2) THEN
        PHASE(J) = A2
      ELSE
        PHASE(J) = A3
      ENDIF

    CASE (4)
      CALL WIGNERFCT (X, NCOEF, 0, 0, D00)
      A4 = 0.0D0
      DO N = 0, NCOEF
        A4 = A4 + WIGCOEF(4,N) * D00(N)
      ENDDO
      PHASE(J) = A4

    CASE (5)
      CALL WIGNERFCT (X, NCOEF, 2, 0, D20)
      B1 = 0.0D0
      DO N = 2, NCOEF
        B1 = B1 - WIGCOEF(5,N) * D20(N)
      ENDDO
      PHASE(J) = B1

    CASE (6)
      CALL WIGNERFCT (X, NCOEF, 2, 0, D20)
      B2 = 0.0D0
      DO N = 2, NCOEF
        B2 = B2 - WIGCOEF(6,N) * D20(N)
      ENDDO
      PHASE(J) = B2

    CASE DEFAULT
      PRINT *, 'Illegal PELEM phase matrix element number'
      STOP
    END SELECT

     ! Calculate P11 function if normalization by it is desired
    IF (PELEM < 0) THEN
      CALL WIGNERFCT (X, NCOEF, 0, 0, D00)
      A1 = 0.0D0
      DO N = 0, NCOEF
        A1 = A1 + WIGCOEF(1,N) * D00(N)
      ENDDO
      NORM(J) = A1
    ENDIF
  ENDDO

   ! Do the P11 normalization if desired
  IF (PELEM < 0) THEN
    PHASE(:) = PHASE(:)/NORM(:)
  ENDIF
  DEALLOCATE (D00, D22P, D22M, D20, NORM)
END SUBROUTINE TRANSFORM_WIGNERD_TO_PHASE



SUBROUTINE WIGNERFCT ( X, NRANK, M, M1, DMM1 )
 ! The routine computes the DMM1N vector coefficients for
 ! M >= 0, M1 >= 0,  N = 0,...,NRANK, and -1 < X= COS(BETA) < 1
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NRANK, M, M1
  DOUBLE PRECISION, INTENT(IN) :: X
  DOUBLE PRECISION, INTENT(OUT) :: DMM1(0:NRANK)
  INTEGER :: N0, N
  DOUBLE PRECISION :: FACT1, FACT2, DMM1_N0

  N0 = MAX( M, M1 )
  DMM1 = 0.D0
  IF ( N0 .EQ. 0 ) THEN
    DMM1(0) = 1.D0
    DMM1(1) = X
    DO N = 1, NRANK - 1
      FACT1 = DBLE( 2 * N + 1 ) * X / DBLE( N + 1 )
      FACT2 = DBLE( N ) / DBLE( N + 1 )
      DMM1(N+1) = FACT1 * DMM1(N) - FACT2 * DMM1(N-1)
    END DO
  ELSE
    DMM1(N0) = DMM1_N0( X, M, M1 )
    DO N = N0, NRANK - 1
      FACT1 = DBLE( N * (N + 1) ) * X - DBLE( M * M1  )
      FACT1 = FACT1 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
      FACT1 = FACT1 / DSQRT( DBLE( (N + 1)**2 - M1**2 ) )
      FACT1 = FACT1 * DBLE( 2 * N + 1 ) / DBLE( N )
      FACT2 = DSQRT( DBLE( N**2 - M**2   ) )  &
                   * DSQRT( DBLE( N**2 - M1**2 ) )
      FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - M**2   ) )
      FACT2 = FACT2 / DSQRT( DBLE( (N + 1)**2 - M1**2 ) )
      FACT2 = FACT2 * DBLE( N + 1 ) / DBLE( N )

      DMM1(N+1) = FACT1 * DMM1(N) - FACT2 * DMM1(N-1)
    END DO
  END IF
END SUBROUTINE WIGNERFCT


FUNCTION DMM1_N0 ( X, M, M1 ) RESULT( DMM1N0 )
 ! The routine computes the Wigner functions for
 !   M >= 0, M1 >= 0, AND N0 = MAX(M,M1)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: M, M1
  DOUBLE PRECISION, INTENT(IN) :: X
  DOUBLE PRECISION :: DMM1N0
  INTEGER :: P, MAXM, MINM
  DOUBLE PRECISION :: PROD, FACT, CMM1

  IF ( M .EQ. M1 ) THEN
    FACT   = ( ( 1.D0 + X ) / 2.D0 )**M
    DMM1N0 = FACT
  ELSE
    IF ( M1 .GT. M ) THEN
      CMM1 = 1.D0
    ELSE
      CMM1 = ( - 1.D0 )**( M - M1 )
    END IF
    MAXM = MAX(M,M1)
    MINM = MIN(M,M1)
    PROD = 1.D0
    DO P = 1, MAXM - MINM
      FACT = DSQRT( DBLE( M + M1 + P ) / DBLE( P ) )
      PROD = PROD * FACT
    END DO
    FACT   = DSQRT( ( 1.D0 - X ) / 2.D0 )
    DMM1N0 = CMM1 * PROD * FACT**( MAXM - MINM )
    FACT   = DSQRT( ( 1.D0 + X ) / 2.D0 )
    DMM1N0 = DMM1N0 * FACT**( MAXM + MINM )
  END IF
END FUNCTION DMM1_N0

