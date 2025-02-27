! PROGRAM PROPGEN
!
!   Makes an SHDOM tabulated phase function optical properties file 
! from a particle properties file specifying the 3D distribution of mass 
! content and effective radius for several types of particles.  The 
! optical properties for these particle types are specified in scattering 
! table files.  The optical properties are linearly interpolated in
! effective radius within the scattering tables.  The extinction and
! single scattering albedo for the mixture specified at each grid point
! are calculated exactly.  The phase functions for the mixtures are
! approximate, however, because there is not enough memory to store
! a phase function for each grid point.  Instead, the closest phase
! function to the correct one is used.  If none of the phase functions
! in the scattering tables are within the user specified tolerances
! then the new mixture phase function is added to the list.  Tolerances
! are specified for the asymmetry parameter and the maximum fractional
! error in 90 phase function values from 2 to 180 degrees.
!  Particle property file format:
!   3        
!   Nx Ny Nz     [number of X, Y, Z grid points]
!   delX delY    [X and Y grid spacing in km]   
!   Z1 ... Zn    [heights of particle levels in km]
!   T1 ... Tn    [temperatures in Kelvin]
!   IX IY IZ  numcomp  Type1 Mass1 Reff1 ... TypeN MassN ReffN
! numcomp is the number of particle components, Type? are the type numbers
! of the components, Mass? are the mass contents [g/m^3] of the components, 
! and Reff? are the effective radii [microns] of the components.
! 
!  compile: pgf90 -fast -o propgen  propgen.f90
!
!    Frank Evans    University of Colorado     May 2003 
!                 Updated for polarization     July 2013


MODULE COMBINE_PROP

CONTAINS

SUBROUTINE USER_INPUT (NSCATTAB, SCATTABFILES, SCATNUMS, POLTABS, PARFILE, &
                       MAXNEWPHASE, ASYMTOL, FRACPHASETOL, &
                       RAYWAVELEN, RAYSFCPRES, NZO, ZOTHER, TEMPOTHER, &
                       POLOUT, PROPFILE)
 ! Reads the input parameters from the standard input
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NSCATTAB, NZO, MAXNEWPHASE
  REAL,    INTENT(OUT) :: RAYWAVELEN, RAYSFCPRES, ASYMTOL, FRACPHASETOL
  CHARACTER(LEN=*), INTENT(OUT) :: PARFILE, PROPFILE
  INTEGER, POINTER     :: SCATNUMS(:)
  LOGICAL, POINTER     :: POLTABS(:)
  REAL,    POINTER     :: ZOTHER(:), TEMPOTHER(:)
  CHARACTER(LEN=80), POINTER :: SCATTABFILES(:)
  CHARACTER(LEN=1) :: POLOUT
  INTEGER :: I

  WRITE(*,*) 'Number of input scattering tables'
  READ(*,*) NSCATTAB
    WRITE(*,*) NSCATTAB

  ALLOCATE (SCATTABFILES(NSCATTAB), SCATNUMS(NSCATTAB), POLTABS(NSCATTAB))
  IF (NSCATTAB > 0) THEN
    WRITE (*,*) 'Name of each scattering table'
    DO I = 1, NSCATTAB
      READ (*,'(A)') SCATTABFILES(I)
      WRITE(*,'(A72)') SCATTABFILES(I)
    ENDDO
    WRITE(*,*) 'Index or type of each scattering table'
    READ (*,*) SCATNUMS(1:NSCATTAB)
      WRITE (*,'(20(1X,I2))') SCATNUMS(1:NSCATTAB)
    WRITE(*,*) 'Polarization flag (T or F) for each scattering table'
    READ (*,*) POLTABS(1:NSCATTAB)
      WRITE (*,'(20(1X,L1))') POLTABS(1:NSCATTAB)
  ENDIF

  WRITE (*,*) 'Input particle properties filename (or NONE for Rayleigh only)'
  READ (*,'(A)') PARFILE
    WRITE (*,'(A72)') PARFILE

  WRITE(*,*) 'Maximum number of new phase functions created'
  READ (*,*) MAXNEWPHASE
    WRITE (*,'(1X,I4)') MAXNEWPHASE

  WRITE(*,*) 'Tolerance in asymmetry parameter for creating a new phase function'
  READ (*,*) ASYMTOL
    WRITE (*,'(1X,F6.4)') ASYMTOL

  WRITE(*,*) 'Maximum fractional error tolerance in value for creating a new phase function'
  READ (*,*) FRACPHASETOL
    WRITE (*,'(1X,F6.3)') FRACPHASETOL

  WRITE(*,*) 'Wavelength for Rayleigh molecular scattering (micron)'
  READ (*,*) RAYWAVELEN
    WRITE (*,*) RAYWAVELEN

  WRITE(*,*) 'Surface pressure for Rayleigh molecular scattering (mb or hPa)'
  READ (*,*) RAYSFCPRES
    WRITE (*,*) RAYSFCPRES

  WRITE (*,*) 'Number of extra Z levels'
  READ (*,*) NZO
    WRITE (*,'(1X,I3)') NZO
  ALLOCATE (ZOTHER(NZO), TEMPOTHER(NZO))
  DO I = 1, NZO
    WRITE (*,*) 'Height (km) and temperature (K) :', I
    READ (*,*) ZOTHER(I), TEMPOTHER(I)
    WRITE (*,'(1X,F7.3,1X,F6.2)') ZOTHER(I), TEMPOTHER(I)
  ENDDO

  WRITE (*,*) 'Output polarized (P) or unpolarized (U) SHDOM property file'
  READ (*,'(A)') POLOUT
  WRITE (*,*) POLOUT   

  WRITE(*,*) 'Output property file name'
  READ(*,'(A)') PROPFILE
    WRITE(*,'(A72)') PROPFILE
END SUBROUTINE USER_INPUT




SUBROUTINE READ_SCAT_TABLE_SIZE (SCATTABFILE, POLTAB, WAVELEN1, WAVELEN2, &
                                 NRETAB, MAXNCOEF)
 ! Reads the scattering table file to find the number of effective radii,
 ! the maximum order of the Wigner d-function series phase functions, and 
 ! also the wavelength range.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABFILE
  LOGICAL, INTENT(IN) :: POLTAB
  INTEGER, INTENT(OUT) :: NRETAB, MAXNCOEF
  REAL,    INTENT(OUT) :: WAVELEN1, WAVELEN2
  INTEGER :: I, J, L, NCOEF
  REAL    :: REFF, EXT, SSALB, CHI

  OPEN (UNIT=1, FILE=SCATTABFILE, STATUS='OLD')
  READ (1,*)
  READ (1,*) WAVELEN1, WAVELEN2
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*) NRETAB
  MAXNCOEF = 1
  DO I = 1, NRETAB
    READ (1,*) REFF, EXT, SSALB, NCOEF
    MAXNCOEF = MAX(MAXNCOEF,NCOEF)
    IF (POLTAB) THEN
      READ (1,*) J, (CHI, L=0,NCOEF)
      READ (1,*) J, (CHI, L=0,NCOEF)
      READ (1,*) J, (CHI, L=0,NCOEF)
      READ (1,*) J, (CHI, L=0,NCOEF)
      READ (1,*) J, (CHI, L=0,NCOEF)
      READ (1,*) J, (CHI, L=0,NCOEF)
    ELSE
      READ (1,*) (CHI, L=0,NCOEF)
    ENDIF
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_SCAT_TABLE_SIZE



SUBROUTINE READ_SCAT_TABLE (SCATTABFILE, POLTAB, NPHASEPOL, NRETAB,RETAB, &
                            EXTINCT, SSALB, MAXCOEF, NCOEF, WIGCOEF)
 ! Reads the table of scattering properties as a function of effective 
 ! radius.  POLTAB=.true. means a polarized scattering table with Wigner 
 ! d-function coefficients for six phase matrix elements instead of one.
 ! The order of the six elements is the four diagonal ones (a1, a2, a3, a4)
 ! followed by the IQ and UV coefficients (b1, b2).
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: SCATTABFILE
  INTEGER, INTENT(IN)  :: NPHASEPOL, NRETAB, MAXCOEF
  LOGICAL, INTENT(IN)  :: POLTAB
  INTEGER, INTENT(OUT) :: NCOEF(:)
  REAL,    INTENT(OUT) :: RETAB(NRETAB), EXTINCT(NRETAB), SSALB(NRETAB)
  REAL,    INTENT(OUT) :: WIGCOEF(NPHASEPOL,0:MAXCOEF,NRETAB)
  INTEGER :: I, J, L
  REAL    :: CHI

  OPEN (UNIT=1, FILE=SCATTABFILE, STATUS='OLD')
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*)
  READ (1,*) 
  DO I = 1, NRETAB
    READ (1,*) RETAB(I), EXTINCT(I), SSALB(I), NCOEF(I)
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



SUBROUTINE MAKE_TYPE_CONV_TABLE (NTYPES, PARTYPES, MAXTYPENUM, TYPECONV)
 ! Makes a table to convert from particle type number to the index in
 ! the arrays - basically the inverse of PARTYPES.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NTYPES, PARTYPES(NTYPES), MAXTYPENUM
  INTEGER, INTENT(OUT) :: TYPECONV(0:MAXTYPENUM)
  INTEGER :: I

  TYPECONV(:) = 0
  DO I = 1, NTYPES
    TYPECONV(PARTYPES(I)) = I
  ENDDO
END SUBROUTINE MAKE_TYPE_CONV_TABLE




SUBROUTINE READ_PARTICLE_FILE_SIZE (PARFILE, NX, NY, NZP)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PARFILE
  INTEGER, INTENT(OUT) :: NX, NY, NZP

  OPEN (UNIT=2, FILE=PARFILE, STATUS='OLD')
  READ (2,*)
  READ (2,*) NX, NY, NZP
  CLOSE (2)
END SUBROUTINE READ_PARTICLE_FILE_SIZE



SUBROUTINE READ_PARTICLE_FILE (PARFILE, NX, NY, NZP, NSCATTAB, PARTYPES, &
                               DELX, DELY, ZPAR, TEMPPAR, &
                               NCOMP, PTYPE, MASSCONT, REFF)
 ! Reads the particle physical properties file.  The file contains a header
 ! containing the array size (NX,NY,NZP), heights (ZPAR), and
 ! temperature profile (TEMPPAR).  Then each row has format:
 !    IX IY IZ  Ncomp  Ptype1 masscont1 Reff1 ... PtypeN masscontN ReffN
 ! where Ncomp is the number of components, Ptype is the particle type,
 ! masscont is the mass content (g/m^3), and Reff is the particle effective
 ! radius (micron).
 !  Also reads the old 2 parameter LWC files with LWC and Reff at each point.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PARFILE
  INTEGER, INTENT(IN) :: NX, NY, NZP, NSCATTAB, PARTYPES(NSCATTAB)
  REAL,    INTENT(OUT) :: DELX, DELY, ZPAR(NZP), TEMPPAR(NZP)
  INTEGER, INTENT(OUT) :: NCOMP(NZP,NY,NX), PTYPE(NSCATTAB,NZP,NY,NX)
  REAL,    INTENT(OUT) :: MASSCONT(NSCATTAB,NZP,NY,NX), REFF(NSCATTAB,NZP,NY,NX)
  INTEGER :: I, FILEKIND, IX, IY, IZ, NC, PT(NSCATTAB)
  REAL    :: MASS(NSCATTAB), RE(NSCATTAB)

   ! Initialize the output arrays (in case there are missing grid points)
  NCOMP(:,:,:) = 0
  PTYPE(:,:,:,:) = 0
  MASSCONT(:,:,:,:) = 0.0
  REFF(:,:,:,:) = 0.0

   ! Read in the header
  OPEN (UNIT=2, FILE=PARFILE, STATUS='OLD')
  READ (2,*) FILEKIND
  READ (2,*) ! NX, NY, NZP
  READ (2,*) DELX, DELY
  READ (2,*) (ZPAR(IZ), IZ=1, NZP)
  READ (2,*) (TEMPPAR(IZ), IZ=1, NZP)

  IF (FILEKIND == 2) THEN
    IF (NSCATTAB /= 1) THEN
      PRINT *,'READ_PARTICLE_FILE: Must have only one scattering table to use 2 parameter LWC file.'
      STOP
    ENDIF
     ! Read in the data 
    DO WHILE (.TRUE.)
      READ (2,*,END=190) IX, IY, IZ, MASS(1), RE(1)
      IF (IX >= 1 .AND. IX <= NX .AND. IY >= 1 .AND. IY <= NY .AND. &
          IZ >= 1 .AND. IZ <= NZP) THEN
        NCOMP(IZ,IY,IX) = 1
        PTYPE(1,IZ,IY,IX) = PARTYPES(1)
        MASSCONT(1,IZ,IY,IX) = MASS(1)
        REFF(1,IZ,IY,IX) = RE(1)
      ENDIF
    ENDDO

  ELSE IF (FILEKIND == 3) THEN
     ! Read in the data 
    DO WHILE (.TRUE.)
      READ (2,*,END=190) IX, IY, IZ, NC, &
                       (PT(I), MASS(I), RE(I), I=1,MIN(NC,NSCATTAB))
      IF (NC > NSCATTAB) THEN
        PRINT *, 'READ_PARTICLE_FILE: More particle components than scattering tables.'
        STOP
      ENDIF
      IF (IX >= 1 .AND. IX <= NX .AND. IY >= 1 .AND. IY <= NY .AND. &
          IZ >= 1 .AND. IZ <= NZP) THEN
        NCOMP(IZ,IY,IX) = NC
        PTYPE(1:NC,IZ,IY,IX) = PT(1:NC)
        MASSCONT(1:NC,IZ,IY,IX) = MASS(1:NC)
        REFF(1:NC,IZ,IY,IX) = RE(1:NC)
      ENDIF
    ENDDO

  ELSE
    PRINT *,'READ_PARTICLE_FILE: Must be type 2 or 3 particle properties file.'
    STOP
  ENDIF

190 CONTINUE
  CLOSE (2)
END SUBROUTINE READ_PARTICLE_FILE




SUBROUTINE ORGANIZE_LEVELS (NZP, ZPAR, TEMPPAR, NZO, ZOTHER, TEMPOTHER, &
                            NZT, ZLEVELS, TEMP, IPARLEV)
 ! Combines the particle and extra levels, returning the resulting
 ! heights (ZLEVELS) and temperature profile (TEMP), and the
 ! pointers to the particle levels (IPARLEV).
  IMPLICIT NONE
  INTEGER, INTENT(IN)  :: NZP, NZO, NZT
  INTEGER, INTENT(OUT) :: IPARLEV(NZT)
  REAL,    INTENT(IN)  :: ZPAR(NZP), TEMPPAR(NZP), ZOTHER(NZO), TEMPOTHER(NZO)
  REAL,    INTENT(OUT) :: ZLEVELS(NZT), TEMP(NZT)
  INTEGER :: J, K

  ! Put the particle and other Z levels (and temperatures) in one array
  ZLEVELS(1:NZP) = ZPAR(:)
  TEMP(1:NZP) = TEMPPAR(:)
  ZLEVELS(NZP+1:NZT) = ZOTHER(:)
  TEMP(NZP+1:NZT) = TEMPOTHER(:)

   ! Warn if extra levels are within the particle layer
  IF (NZP > 0) THEN
    DO J = 1, NZO
      IF (ZOTHER(J) >= ZPAR(1) .AND. ZOTHER(J)<= ZPAR(NZP)) THEN
        PRINT '(1X,A,F7.3)', 'Warning: other Z level inside range of particle levels:',ZOTHER(J)
        PRINT *, '         The other Z levels have only Rayleigh scattering.'
      ENDIF
    ENDDO
  ENDIF

  ! Sort the levels
  CALL SSORT (ZLEVELS, TEMP, NZT, 2)

  ! Check to see if two Z levels are same
  DO K = 1, NZT-1
    IF (ZLEVELS(K) == ZLEVELS(K+1)) THEN
      PRINT *, 'ORGANIZE_LEVELS: Same Z levels: ', ZLEVELS(K)
      STOP
    ENDIF
  ENDDO
  ! Make the particle level pointers
  IF (NZP == 0) THEN
    IPARLEV(:) = 0
  ELSE
    J = 1
    DO K = 1, NZT
      IPARLEV(K) = 0
      IF (ZLEVELS(K) == ZPAR(J)) THEN
        IPARLEV(K) = J
        J = MIN(NZP,J+1)
      ENDIF
    ENDDO 
  ENDIF
END SUBROUTINE ORGANIZE_LEVELS



SUBROUTINE SSORT (X, Y, N, KFLAG)
!   SSORT sorts array X and optionally makes the same interchanges in
!   array Y.  The array X may be sorted in increasing order or
!   decreasing order.  A slightly modified quicksort algorithm is used.
!
!   Description of Parameters
!      X - array of values to be sorted   (usually abscissas)
!      Y - array to be (optionally) carried along
!      N - number of values in array X to be sorted
!      KFLAG - control parameter
!            =  2  means sort X in increasing order and carry Y along.
!            =  1  means sort X in increasing order (ignoring Y)
!            = -1  means sort X in decreasing order (ignoring Y)
!            = -2  means sort X in decreasing order and carry Y along.
!
      IMPLICIT NONE
      INTEGER KFLAG, N
      REAL X(*), Y(*)
      REAL R, T, TT, TTY, TY
      INTEGER I, IJ, J, K, KK, L, M, NN
      INTEGER IL(21), IU(21)

      NN = N
      IF (NN .LT. 1) THEN
        STOP 'The number of values to be sorted is not positive.'
      ENDIF
      KK = ABS(KFLAG)
      IF (KK.NE.1 .AND. KK.NE.2) THEN
        STOP 'The sort control parameter, K, is not 2, 1, -1, or -2.'
      ENDIF

!       Alter array X to get decreasing order if needed
      IF (KFLAG .LE. -1) THEN
        DO 10 I=1,NN
            X(I) = -X(I)
   10    CONTINUE
      ENDIF

      IF (KK .EQ. 2) GO TO 100

!       Sort X only
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

!       Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)
!       If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
      ENDIF
      L = J
!       If last element of array is less than than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
!          If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
         ENDIF
      ENDIF
!       Find an element in the second half of the array which is smaller than T
   40 L = L-1
      IF (X(L) .GT. T) GO TO 40
!       Find an element in the first half of the array which is greater than T
   50 K = K+1
      IF (X(K) .LT. T) GO TO 50
!       Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         GO TO 40
      ENDIF
!       Save upper and lower subscripts of the array yet to be sorted
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
!       Begin again on another portion of the unsorted array
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
!       Sort X and carry Y along
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
  120 K = I
!       Select a central element of the array and save it in location T
      IJ = I + INT((J-I)*R)
      T = X(IJ)
      TY = Y(IJ)
!       If first element of array is greater than T, interchange with T
      IF (X(I) .GT. T) THEN
         X(IJ) = X(I)
         X(I) = T
         T = X(IJ)
         Y(IJ) = Y(I)
         Y(I) = TY
         TY = Y(IJ)
      ENDIF
      L = J
!       If last element of array is less than T, interchange with T
      IF (X(J) .LT. T) THEN
         X(IJ) = X(J)
         X(J) = T
         T = X(IJ)
         Y(IJ) = Y(J)
         Y(J) = TY
         TY = Y(IJ)
!          If first element of array is greater than T, interchange with T
         IF (X(I) .GT. T) THEN
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
         ENDIF
      ENDIF
!       Find an element in the second half of the array which is smaller than T
  130 L = L-1
      IF (X(L) .GT. T) GO TO 130
!       Find an element in the first half of the array which is greater than T
  140 K = K+1
      IF (X(K) .LT. T) GO TO 140
!       Interchange these elements
      IF (K .LE. L) THEN
         TT = X(L)
         X(L) = X(K)
         X(K) = TT
         TTY = Y(L)
         Y(L) = Y(K)
         Y(K) = TTY
         GO TO 130
      ENDIF
!       Save upper and lower subscripts of the array yet to be sorted
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

!       Begin again on another portion of the unsorted array
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
!       Clean up
  190 IF (KFLAG .LE. -1) THEN
         DO 200 I=1,NN
            X(I) = -X(I)
  200    CONTINUE
      ENDIF
END SUBROUTINE SSORT



SUBROUTINE RAYLEIGH_EXTINCT (NZT, ZLEVELS,TEMP, RAYSFCPRES, RAYLCOEF, EXTRAYL)
 ! Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
 ! from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
 ! a linear lapse rate between levels to compute the pressure at
 ! each level starting with the SFCPRES [mb].  The Rayleigh extinction is 
 ! proportional to air density, with the coefficient RAYLCOEF in [K/(mb km)].
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NZT
  REAL,    INTENT(IN) :: ZLEVELS(NZT), TEMP(NZT), RAYSFCPRES, RAYLCOEF
  REAL,    INTENT(OUT) :: EXTRAYL(NZT)
  INTEGER :: I
  REAL    :: PRES, LAPSE, TS, DZ
      
  ! Find surface pressure by integrating hydrostatic relation
  !  for a dry atmosphere up to surface height.
  PRES = RAYSFCPRES
  TS = TEMP(1)
  LAPSE = 6.5*0.001
  PRES = PRES*(TS/(TS+LAPSE*ZLEVELS(1)*1000.))**(9.8/(287.*LAPSE))

  ! Use layer mean temperature to compute fractional pressure change.
  DO I = 1, NZT-1
    EXTRAYL(I) = RAYLCOEF*PRES/TEMP(I)
    DZ = 1000.*(ZLEVELS(I+1)-ZLEVELS(I))
    LAPSE = (TEMP(I)-TEMP(I+1))/DZ
    IF (ABS(LAPSE) .GT. 0.00001) THEN
      PRES = PRES*(TEMP(I+1)/TEMP(I))**(9.8/(287.*LAPSE))
    ELSE
      PRES = PRES*EXP(-9.8*DZ/(287.*TEMP(I)))
    ENDIF
  ENDDO
  EXTRAYL(NZT) = RAYLCOEF*PRES/TEMP(NZT)
END SUBROUTINE RAYLEIGH_EXTINCT



SUBROUTINE RAYLEIGH_PHASE_FUNCTION (WAVELEN, NPHASEPOL, MAXCOEF, &
                                    NCOEF, RAYWIGCOEF)
 ! Returns the Wigner d-function coefficients for either the Rayleigh 
 ! scalar phase function or polarized phase matrix for molecular 
 ! scattering by air.  Includes the wavelength depolarization factor.  
 ! From Mishchenko's book "Multiple scattering of light by particles: 
 ! Radiative Transfer and Coherent Backscattering", Cambridge, 2006.  
 ! Thanks to Adrian Doicu.   WAVELEN is the wavelength in microns.  
 ! NPHASEPOL=1 for the scalar case, and 6 for the phase matrix case.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NPHASEPOL, MAXCOEF
  REAL,    INTENT(IN) :: WAVELEN
  INTEGER, INTENT(OUT) :: NCOEF
  REAL,    INTENT(OUT) :: RAYWIGCOEF(NPHASEPOL,0:MAXCOEF)
  DOUBLE PRECISION, PARAMETER :: AKING=1.0469541D0, &
                    BKING=3.2503153D-04, CKING=3.8622851D-05
  REAL :: FKING, DEPOL, DELTA, DELTAP

  NCOEF = 2
  RAYWIGCOEF(1:NPHASEPOL,:) = 0.0
  FKING = AKING + BKING/WAVELEN**2 + CKING/WAVELEN**4
  DEPOL = 6.D0*(FKING-1.D0) / (3.D0 + 7.D0*FKING)
  DELTA = (1.0 - DEPOL) / (1.0 + 0.5*DEPOL) 
  DELTAP = (1.0 - 2.0*DEPOL) / (1.0 - DEPOL) 
  RAYWIGCOEF(1,0) = 1.0  ; RAYWIGCOEF(1,2) = 0.5*DELTA
  IF (NPHASEPOL == 6) THEN
    RAYWIGCOEF(2,2) = 3.0*DELTA
    RAYWIGCOEF(4,1) = 1.5*DELTAP*DELTA
    RAYWIGCOEF(5,2) = SQRT(1.5)*DELTA
  ENDIF
END SUBROUTINE RAYLEIGH_PHASE_FUNCTION




SUBROUTINE PHASEFUNC_FROM_LEGENDRE (NANGLES, NCOEF, LEGCOEF, PHASE)
 ! Calculates the phase function at NANGLES equally spaced angles
 ! from the Legendres series expansion in LEGCOEF.  The scalar (I,I)
 ! element of the phase matrix is expanded in Wigner d^l_{0,0} functions,
 ! which are the same as Legendre polynomials in \cos\Theta.
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: NANGLES, NCOEF
  REAL,    INTENT(IN) :: LEGCOEF(0:NCOEF)
  REAL,    INTENT(OUT) :: PHASE(NANGLES)
  INTEGER :: J, L
  REAL(8) :: RD, MU, SUM, PL, PL1, PL2
  
  RD = ACOS(-1.0D0)/180.D0
  DO J = 1, NANGLES
    MU = COS(RD*J*180.D0/NANGLES)
    SUM = 0.0
      ! Use upward recurrence to find Legendre polynomials
    PL1 = 1.0
    PL = 1.0 
    DO L = 0, NCOEF
      IF (L .GT. 0) PL = (2*L-1)*MU*PL1/L-(L-1)*PL2/L
      SUM = SUM + LEGCOEF(L)*PL
      PL2 = PL1
      PL1 = PL 
    ENDDO
    PHASE(J) = SUM
  ENDDO
END SUBROUTINE PHASEFUNC_FROM_LEGENDRE



SUBROUTINE OUTPUT_PROPERTIES (PROPFILE, NPHASEPOL, NX, NY, NZT, &
                              DELX,DELY, ZLEVELS,TEMP, &
                              NPHASE, MAXCOEF, NCOEF, WIGCOEF, &
                              EXTINCT, SSALB, IPHASE)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: PROPFILE
  INTEGER, INTENT(IN) :: NPHASEPOL, NX, NY, NZT, NPHASE, MAXCOEF
  INTEGER, INTENT(IN) :: NCOEF(0:NPHASE), IPHASE(NZT,NY,NX)
  REAL,    INTENT(IN) :: DELX, DELY, ZLEVELS(NZT), TEMP(NZT)
  REAL,    INTENT(IN) :: WIGCOEF(NPHASEPOL,0:MAXCOEF,0:NPHASE)
  REAL,    INTENT(IN) :: EXTINCT(NZT,NY,NX), SSALB(NZT,NY,NX)
  INTEGER :: NPHASEOUT, I, IX, IY, IZ, J, K, L, NL
  INTEGER, ALLOCATABLE :: PHASECONV(:)

  ! Make a phase function index conversion array to output only the ones used.
  ALLOCATE (PHASECONV(0:NPHASE))
  PHASECONV(:) = 0
  NPHASEOUT = 0
  DO I = 0, NPHASE
    IF (ANY(IPHASE(:,:,:) == I)) THEN
      NPHASEOUT = NPHASEOUT + 1
      PHASECONV(I) = NPHASEOUT
    ENDIF
  ENDDO

  ! Open the property file and output the header including the
  !   tabulated phase functions
  OPEN (UNIT=3, FILE=PROPFILE, STATUS='REPLACE')
  IF (NPHASEPOL == 1) THEN
    WRITE (3,'(A)') 'Tabulated phase function property file'
  ELSE
    WRITE (3,'(A)') 'Polarized tabulated phase function property file'
  ENDIF
  WRITE (3,'(3(1X,I5))') NX, NY, NZT
  WRITE (3,'(2(1X,F7.4),202(1X,F8.4))') DELX, DELY, (ZLEVELS(I),I=1,NZT)
  WRITE (3,'(1X,I4)') NPHASEOUT
  DO I = 0, NPHASE
    IF (PHASECONV(I) > 0) THEN
      NL = NCOEF(I)
      DO J = 1, NPHASEPOL
        IF (NPHASEPOL == 1) THEN
          WRITE (3,'(1X,I4,5X,200(1X,F9.5))') NL,(WIGCOEF(J,L,I), L=1,MIN(NL,200))
        ELSE
          WRITE (3,'(1X,I1,1X,I4,5X,201(1X,F9.5))') J,NL,(WIGCOEF(J,L,I), L=0,MIN(NL,200))
        ENDIF
        DO K = 200, NL-1, 200
          WRITE (3,'(200(1X,F9.5))') (WIGCOEF(J,K+L,I), L=1,MIN(NL-K,200))
        ENDDO
      ENDDO
    ENDIF
  ENDDO  

  ! Write one line of optical properties for each grid point
  DO IX = 1, NX
    DO IY = 1, NY
      DO IZ = 1, NZT
          WRITE (3,'(3(1X,I4),1X,F5.1,1X,F9.4,1X,F8.6,1X,I3)') IX, IY, IZ, &
            TEMP(IZ), EXTINCT(IZ,IY,IX), SSALB(IZ,IY,IX), PHASECONV(IPHASE(IZ,IY,IX))
      ENDDO
    ENDDO
  ENDDO
  CLOSE (3)
END SUBROUTINE OUTPUT_PROPERTIES 


END MODULE COMBINE_PROP






PROGRAM PROPGEN
  USE COMBINE_PROP
  IMPLICIT NONE
  INTEGER :: NSCATTAB, MAXNEWPHASE, NZO
  REAL    :: ASYMTOL, FRACPHASETOL, RAYWAVELEN, RAYSFCPRES
  INTEGER, POINTER :: SCATNUMS(:)
  LOGICAL, POINTER :: POLTABS(:)
  REAL,    POINTER :: ZOTHER(:), TEMPOTHER(:)
  CHARACTER(LEN=80), POINTER :: SCATTABFILES(:)
  CHARACTER(LEN=80) :: PARFILE, PROPFILE
  CHARACTER(LEN=1) :: POLOUT
  INTEGER, PARAMETER :: NANGLES=90
  INTEGER :: NTYPES, NENTRIES, MAXNPHASE, NPHASE, NPHASEPOL
  INTEGER :: I, J, N1, N2, NRETAB, MAXCOEF, MAXNCOEF, MAXTYPENUM
  INTEGER :: NX, NY, NZP, NZT, IX, IY, IZ, IZP, NC, NL, IL, IM, IU, L, IPH
  REAL    :: RAYLCOEF
  REAL    :: DELX, DELY, EXT, SCAT, SCATTER, ASYM0, F, ERR, MINERR
  INTEGER, ALLOCATABLE :: NRE(:), NCOEF(:)
  INTEGER, ALLOCATABLE :: PARTYPES(:), TYPECONV(:)
  INTEGER, ALLOCATABLE :: NCOMP(:,:,:), PTYPE(:,:,:,:), IPARLEV(:)
  REAL, ALLOCATABLE    :: WAVELEN1(:), WAVELEN2(:)
  REAL, ALLOCATABLE    :: WIGCOEF(:,:), PHASE0(:)
  REAL, ALLOCATABLE    :: RETAB(:), EXTINCTTAB(:), SSALBTAB(:), WIGCOEFTAB(:,:,:)
  REAL, ALLOCATABLE    :: ZPAR(:), TEMPPAR(:), ZLEVELS(:), TEMP(:), EXTRAYL(:)
  REAL, ALLOCATABLE    :: ASYM(:), PHASE(:,:)
  REAL, ALLOCATABLE    :: MASSCONT(:,:,:,:), REFF(:,:,:,:)
  REAL, ALLOCATABLE    :: EXTINCT(:,:,:), SSALB(:,:,:)
  INTEGER, ALLOCATABLE :: IPHASE(:,:,:)

  ! Data structures:


  CALL USER_INPUT (NSCATTAB, SCATTABFILES, SCATNUMS, POLTABS, PARFILE, &
                   MAXNEWPHASE, ASYMTOL, FRACPHASETOL, &
                   RAYWAVELEN, RAYSFCPRES,  NZO, ZOTHER, TEMPOTHER, &
                   POLOUT, PROPFILE)
  MAXNEWPHASE = MAX(MAXNEWPHASE,2)
  ASYMTOL = MAX(ASYMTOL,0.002)
  FRACPHASETOL = MAX(FRACPHASETOL,0.01)

  IF (POLOUT == 'P' .OR. POLOUT == 'p') THEN
    NPHASEPOL = 6
  ELSE
    NPHASEPOL = 1
  ENDIF

  NTYPES = NSCATTAB
  ALLOCATE (WAVELEN1(NSCATTAB), WAVELEN2(NSCATTAB))
  ALLOCATE (NRE(0:NTYPES), PARTYPES(NTYPES))

  NRE(0) = 1
  MAXCOEF = 2
  DO I = 1, NSCATTAB
    CALL READ_SCAT_TABLE_SIZE (SCATTABFILES(I), POLTABS(I), &
                               WAVELEN1(I), WAVELEN2(I), NRETAB, MAXNCOEF)
    NRE(I) = NRE(I-1) + NRETAB
    MAXCOEF = MAX(MAXCOEF,MAXNCOEF)
    IF (I > 1 .AND. ABS((WAVELEN1(I)-WAVELEN1(1))/WAVELEN1(1)) > 0.01 &
        .OR. ABS((WAVELEN2(I)-WAVELEN2(1))/WAVELEN2(1)) > 0.01) THEN
      PRINT *, 'Warning: scattering table wavelengths do not agree for table',I
    ENDIF
  ENDDO
  DEALLOCATE (WAVELEN1, WAVELEN2)
  NENTRIES = NRE(NSCATTAB)-1
  MAXNPHASE = NENTRIES + MAXNEWPHASE

  ALLOCATE (RETAB(NENTRIES), EXTINCTTAB(NENTRIES), SSALBTAB(NENTRIES))
  ALLOCATE (NCOEF(0:MAXNPHASE), WIGCOEFTAB(NPHASEPOL,0:MAXCOEF,0:MAXNPHASE))

  WIGCOEFTAB(:,:,:) = 0.0
  ! Put the Rayleigh phase matrix in the zero index
  CALL RAYLEIGH_PHASE_FUNCTION (RAYWAVELEN, NPHASEPOL, MAXCOEF, &
                                NCOEF(0), WIGCOEFTAB(:,:,0))

  ! Read in the optical properties from the scattering table files
  DO I = 1, NSCATTAB
    N1 = NRE(I-1)  ;  N2 = NRE(I)-1
    CALL READ_SCAT_TABLE (SCATTABFILES(I), POLTABS(I), NPHASEPOL, &
                          N2-N1+1, RETAB(N1:N2), &
                          EXTINCTTAB(N1:N2), SSALBTAB(N1:N2), &
                          MAXCOEF, NCOEF(N1:N2), WIGCOEFTAB(:,0:,N1:N2) )
  ENDDO
  PARTYPES(:) = SCATNUMS(:)
  DEALLOCATE (SCATNUMS, SCATTABFILES)

  ! Make the conversion table from particle type to scat table index
  MAXTYPENUM = MAXVAL(PARTYPES(:))
  ALLOCATE (TYPECONV(0:MAXTYPENUM))
  CALL MAKE_TYPE_CONV_TABLE (NTYPES, PARTYPES, MAXTYPENUM, TYPECONV)


  ! Read in the particle properties file
  IF (PARFILE /= 'NONE' .AND. PARFILE /= 'none') THEN
    CALL READ_PARTICLE_FILE_SIZE (PARFILE, NX, NY, NZP)
    ALLOCATE (ZPAR(NZP), TEMPPAR(NZP))
    ALLOCATE (NCOMP(NZP,NY,NX), PTYPE(NSCATTAB,NZP,NY,NX))
    ALLOCATE (MASSCONT(NSCATTAB,NZP,NY,NX), REFF(NSCATTAB,NZP,NY,NX))
    CALL READ_PARTICLE_FILE (PARFILE, NX, NY, NZP, NSCATTAB, PARTYPES, &
                             DELX, DELY, ZPAR, TEMPPAR, &
                             NCOMP, PTYPE, MASSCONT, REFF)
  ELSE
    NX=1 ; NY=1 ; NZP=0
    DELX=1.0 ; DELY=1.0
  ENDIF

   ! Combine the particle levels and extra levels
  NZT = NZP+NZO
  ALLOCATE (ZLEVELS(NZT), TEMP(NZT), IPARLEV(NZT), EXTRAYL(NZT))
  CALL ORGANIZE_LEVELS (NZP, ZPAR, TEMPPAR, NZO, ZOTHER, TEMPOTHER, &
                        NZT, ZLEVELS, TEMP, IPARLEV)
 

   ! Use the parameterization of Bodhaine et al. (1999) eq 30 for tau_R at 
   ! sea level and convert to Rayleigh density coefficient:
   !   k = 0.03370*(p_sfc/1013.25)*tau_R  for k in K/(mb km) (from assumptions in RAYLEIGH_EXTINCT)
  RAYLCOEF = 0.03370*(RAYSFCPRES/1013.25)* 0.0021520 &
          *(1.0455996 - 341.29061/RAYWAVELEN**2 - 0.90230850*RAYWAVELEN**2) &
          / (1 + 0.0027059889/RAYWAVELEN**2 - 85.968563*RAYWAVELEN**2)
  print *, 'Rayleigh total atmosphere optical depth: ', RAYLCOEF/0.03370
    ! Compute the molecular Rayleigh scattering extinction profile
  CALL RAYLEIGH_EXTINCT (NZT, ZLEVELS, TEMP, RAYSFCPRES, RAYLCOEF, EXTRAYL)

        
  ALLOCATE (ASYM(0:MAXNPHASE), PHASE(NANGLES,0:MAXNPHASE))
  NPHASE = NENTRIES
   ! Make the asymmetry param and angle space low resolution phase functions
  DO I = 0, NPHASE
    ASYM(I) = WIGCOEFTAB(1,1,I)/3
    CALL PHASEFUNC_FROM_LEGENDRE (NANGLES, NCOEF(I), WIGCOEFTAB(1,:,I), PHASE(:,I))
  ENDDO

  ! Allocate the output optical property arrays
  ALLOCATE (EXTINCT(NZT,NY,NX), SSALB(NZT,NY,NX), IPHASE(NZT,NY,NX))

  ! Do the optical property mixing from the components specified at each
  !   grid point in the particle file.
  ALLOCATE (WIGCOEF(NPHASEPOL,0:MAXCOEF), PHASE0(NANGLES))
  DO IX = 1, NX
   DO IY = 1, NY
    DO IZ = 1, NZT
      IZP = IPARLEV(IZ)
      IF (IZP == 0) THEN
        EXTINCT(IZ,IY,IX) = EXTRAYL(IZ)  ! Only Rayleigh scattering if not
        SSALB(IZ,IY,IX) = 1.0            !   at a particle level
        IPHASE(IZ,IY,IX) = 0
      ELSE
        NC = NCOMP(IZP,IY,IX)
        ! Initialize optical property sums with Rayleigh scattering
        EXTINCT(IZ,IY,IX) = EXTRAYL(IZ)
        SCATTER = EXTRAYL(IZ)
        ASYM0 = 0.0
        PHASE0(:) = EXTRAYL(IZ)*PHASE(:,0)
        NL = 2
        WIGCOEF(:,:) = EXTRAYL(IZ)*WIGCOEFTAB(:,:,0)
        DO J = 1, NC
          L = TYPECONV(PTYPE(J,IZP,IY,IX))
          ! Do binary search to find effective radius entry in table
          IL = NRE(L-1)
          IU = NRE(L)-1
          DO WHILE (IU-IL > 1)
            IM = (IL+IU)/2
            IF (REFF(J,IZP,IY,IX) >= RETAB(IM)) THEN
              IL = IM
            ELSE
              IU = IM
            ENDIF
          ENDDO
          IF (IU > IL) THEN
            F = (REFF(J,IZP,IY,IX) - RETAB(IL)) / (RETAB(IU)-RETAB(IL))
          ELSE IF (IU == IL .AND. &
               ABS((REFF(J,IZP,IY,IX)-RETAB(IL))/RETAB(IL)) < 0.001) THEN
            F = 0.0
          ELSE
            F = -1.0
          ENDIF
          IF (F < 0.0 .OR. F > 1.0) THEN
            PRINT *, 'Warning: effective radius outside of table (IX,IY,IZ,type,Reff):'
            PRINT '(4(1X,I3),1X,F6.2)', IX, IY, IZP, PTYPE(J,IZP,IY,IX), REFF(J,IZP,IY,IX)
            F = MAX(MIN(F,1.0),0.0)
          ENDIF
          ! Interpolate optical properties in Reff for this component
          !   and sum over the components.
          EXT = MASSCONT(J,IZP,IY,IX)* &
                EXP((1-F)*LOG(EXTINCTTAB(IL)) + F*LOG(EXTINCTTAB(IU)))
          SCAT = EXT* ((1-F)*SSALBTAB(IL) + F*SSALBTAB(IU))
          ASYM0 = ASYM0 + SCAT* ((1-F)*ASYM(IL) + F*ASYM(IU))
          PHASE0(:) = PHASE0(:) + SCAT* ((1-F)*PHASE(:,IL) + F*PHASE(:,IU))
          NL = MAX(NL,NCOEF(IL),NCOEF(IU))
          WIGCOEF(:,0:NL) = WIGCOEF(:,0:NL) &
                      + SCAT* ((1-F)*WIGCOEFTAB(:,0:NL,IL) + F*WIGCOEFTAB(:,0:NL,IU))
          SCATTER = SCATTER + SCAT
          EXTINCT(IZ,IY,IX) = EXTINCT(IZ,IY,IX) + EXT
        ENDDO ! end of component loop

        ! Normalize to get the combined optical properties 
        IF (EXTINCT(IZ,IY,IX) > 0.0) THEN
          SSALB(IZ,IY,IX) = SCATTER/EXTINCT(IZ,IY,IX)
        ELSE
          SSALB(IZ,IY,IX) = 1.0
        ENDIF
        IF (SCATTER > 0.0) THEN
          ASYM0 = ASYM0/SCATTER
          PHASE0(:) = PHASE0(:)/SCATTER
          WIGCOEF(:,0:NL) = WIGCOEF(:,0:NL)/SCATTER

          ! Find the closest phase function in terms of the asymmetry parameter
          !   and max fractional function value error normalized by user tolerances
          IPH = 0    ! initialize with the Rayleigh phase function
          MINERR = ABS(ASYM0-ASYM(IPH))/ASYMTOL &
             + MAXVAL(ABS((PHASE0(:)-PHASE(:,IPH))/(0.001+PHASE(:,IPH))))/FRACPHASETOL
          DO I = 1, NPHASE         ! search all phase functions
            ERR = ABS(ASYM0-ASYM(I))/ASYMTOL &
              + MAXVAL(ABS((PHASE0(:)-PHASE(:,I))/MAX(0.001,PHASE(:,I))))/FRACPHASETOL
            IF (ERR < MINERR) THEN
              MINERR = ERR
              IPH = I
            ENDIF
          ENDDO
          ! If the closest phase function is within the tolerances then use it,
          !   otherwise add the current grid point phase function to the list.
          IF (ABS(ASYM0-ASYM(IPH)) < ASYMTOL .AND. FRACPHASETOL > &
             MAXVAL(ABS((PHASE0(:)-PHASE(:,IPH))/MAX(0.001,PHASE(:,IPH)))) ) THEN
            IPHASE(IZ,IY,IX) = IPH
          ELSE
            NPHASE = NPHASE + 1
            IF (NPHASE > MAXNPHASE) THEN
              PRINT *, 'Maximum number of new phase functions exceeded',NPHASE-NENTRIES
              STOP
            ENDIF
            NCOEF(NPHASE) = NL
            WIGCOEFTAB(:,:,NPHASE) = WIGCOEF(:,:)
            ASYM(NPHASE) = WIGCOEFTAB(1,1,NPHASE)/3
            CALL PHASEFUNC_FROM_LEGENDRE (NANGLES, NL, WIGCOEF(1,0:), &
                                          PHASE(:,NPHASE))
            IPHASE(IZ,IY,IX) = NPHASE
          ENDIF

        ELSE  ! if no scattering then just use the Rayleigh phase function
          IPHASE(IZ,IY,IX) = 0
        ENDIF
      ENDIF
    ENDDO  ! end of Z level loop
   ENDDO
  ENDDO

  PRINT *, 'Number of new phase functions added: ',NPHASE-NENTRIES

  CALL OUTPUT_PROPERTIES (PROPFILE, NPHASEPOL, NX, NY, NZT, &
                          DELX, DELY, ZLEVELS, TEMP, &
                          NPHASE, MAXCOEF, NCOEF, WIGCOEFTAB, &
                          EXTINCT, SSALB, IPHASE)

  DEALLOCATE (EXTINCT, SSALB, IPHASE, WIGCOEF, ASYM, PHASE0, PHASE)
  IF (ALLOCATED (MASSCONT)) THEN
    DEALLOCATE (ZPAR, TEMPPAR, NCOMP, PTYPE, MASSCONT, REFF)
  ENDIF
  DEALLOCATE (ZLEVELS, TEMP, IPARLEV)
  DEALLOCATE (NCOEF, WIGCOEFTAB, RETAB, EXTINCTTAB, SSALBTAB)
  DEALLOCATE (TYPECONV, PARTYPES, NRE)
  DEALLOCATE (POLTABS)

END



