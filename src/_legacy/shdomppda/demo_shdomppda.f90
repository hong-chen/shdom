PROGRAM demo_shdomppda
 !  
 ! Demonstration program for using the SHDOMPPDA module.
 ! See shdomppda_fwd_tl_adj.f90 for documentation on arguments of
 ! the subroutines.
 ! 
 ! Compile with "make" using the "makefile"

  use molec_absorption_interface  ! dummy molecular absorption routines
  use shdomppda_fwd_tl_adj

  implicit none
  integer, parameter :: maxcomp=7, maxnout=10
  character(len=75) :: instatefile, outstatefile, outradfile
  character(len=75) :: scattabdir, scattabfiles(maxcomp)
  integer :: nx, ny, nz
  integer :: ncomp, icomps(maxcomp)
  integer :: nout, nmu, nphi
  integer :: indexabs, ix, iy, i, ichan
  real    :: specieszrange(2,maxcomp)
  real    :: sfcpres, sfctemp, sfcalbedo, solarflux, solarmu, wavelen
  real    :: muobs(maxnout), phiobsrel(maxnout), splitacc, solacc
  real    :: topmolectau=0.0, toptemp=200.0
  real    :: pert_radout(maxnout), adradout(maxnout)
  real, allocatable :: zlev(:,:,:), temp(:,:,:), pres(:,:,:)
  real, allocatable :: vapmixr(:,:,:), o3mixr(:)
  real, allocatable :: mixr(:,:,:,:), conc(:,:,:,:)
  real, allocatable :: radout(:,:,:)
  real, allocatable :: adtemp(:,:,:), advapmixr(:,:,:), ado3mixr(:)
  real, allocatable :: admixr(:,:,:,:), adconc(:,:,:,:)
  real, allocatable :: adsfctemp(:,:), adsfcalbedo(:,:)
  character(len=1) :: srctype, units


  call user_input (instatefile, outradfile, outstatefile, &
                   ncomp, icomps, scattabfiles, specieszrange, &
                   srctype, solarflux, solarmu, nout, muobs, phiobsrel, &
                   sfctemp, sfcalbedo, wavelen, units, &
                   nmu, nphi, splitacc, solacc, &
                   pert_radout)
  if (nout > maxnout) stop 'maxnout exceeded'

   ! Dummy call to initialize molecular absorption module
  call initialize_molec_abs ()

   ! Read the SHDOMPPDA scattering tables for the channel
  scattabdir = ''
  call READ_SCATTERING_TABLES (ncomp, 1, specieszrange, &
                               scattabdir, scattabfiles)

   ! Read in the number of levels and allocate arrays
  call READ_STATE_SIZE (instatefile, nz, ny, nx)
  allocate (zlev(nz,ny,nx), temp(nz,ny,nx), pres(nz,ny,nx))
  allocate (vapmixr(nz,ny,nx), o3mixr(nz))
  allocate (mixr(ncomp,nz,ny,nx), conc(ncomp,nz,ny,nx))
  allocate (radout(nout,ny,nx))
  allocate (adtemp(nz,ny,nx), advapmixr(nz,ny,nx), ado3mixr(nz))
  allocate (admixr(ncomp,nz,ny,nx), adconc(ncomp,nz,ny,nx))
  allocate (adsfctemp(ny,nx), adsfcalbedo(ny,nx))

   ! Read in the single column of the state fields
  call READ_STATE_FILE (instatefile, ncomp, icomps, nz, ny, nx, &
                        zlev, temp, pres, vapmixr, mixr, conc)
  o3mixr(:) = 0.0
  indexabs = 0    ! no molecular absorption
  ichan = 1       ! only one wavelength in this problem

     ! Loop over columns in the state file
  do ix = 1, nx
    do iy = 1, ny

       ! Call the forward radiative transfer model
      sfcpres = pres(1,iy,ix)
      call FORWARD_COLUMN_RADIANCE (nz, zlev(:,iy,ix), temp(:,iy,ix), &
                       pres(:,iy,ix), vapmixr(:,iy,ix), o3mixr(:), &
                       topmolectau, toptemp, indexabs, &
                       ncomp, mixr(:,:,iy,ix), conc(:,:,iy,ix), specieszrange, &
                       sfcpres, sfctemp, sfcalbedo, &
                       ichan, nout, muobs, phiobsrel, &
                       srctype, solarflux, solarmu, wavelen, &
                       nmu, nphi, splitacc, solacc, &
                       units, radout(:,iy,ix))

       ! Compute the adjoint with a specified change in this single radiance
      adradout(1:nout) = pert_radout(1:nout)
      adtemp(:,iy,ix) = 0.0
      advapmixr(:,iy,ix) = 0.0
      ado3mixr(:) = 0.0
      admixr(:,:,iy,ix) = 0.0
      adconc(:,:,iy,ix) = 0.0
      adsfctemp(iy,ix) = 0.0
      adsfcalbedo(iy,ix) = 0.0
      call ADJOINT_COLUMN_RADIANCE (nz, zlev(:,iy,ix), temp(:,iy,ix), &
                      pres(:,iy,ix), vapmixr(:,iy,ix), o3mixr(:), indexabs, &
                      ncomp, mixr(:,:,iy,ix), conc(:,:,iy,ix), &
                      sfcpres, sfctemp, sfcalbedo, &
                      ichan, nout, muobs, phiobsrel, &
                      srctype, solarflux, solarmu, wavelen, &
                      nmu, nphi,  units,  adradout(1:nout), &
                      adtemp(:,iy,ix), advapmixr(:,iy,ix), ado3mixr(:), &
                      admixr(:,:,iy,ix), adconc(:,:,iy,ix), &
                      adsfctemp(iy,ix), adsfcalbedo(iy,ix))
      call dealloc_column_radiance (indexabs)

    enddo   ! End of loop over columns
  enddo


  call dealloc_scattering_tables

  call write_radiance_file (outradfile, nout, ny, nx, muobs, phiobsrel, radout)

  call write_state_file (outstatefile, ncomp, nz, ny, nx, &
                         zlev, temp, pres, advapmixr, admixr, adconc)

END





SUBROUTINE READ_STATE_SIZE (STATEFILE, NZ, NY, NX)
 ! Reads the array size state file.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: STATEFILE
  INTEGER, INTENT(OUT)  :: NZ, NY, NX

  ! Open the file
  OPEN (UNIT=1, FILE=STATEFILE, STATUS='OLD')
  READ (1,*) NX,NY,NZ
  CLOSE (1)
END SUBROUTINE READ_STATE_SIZE



SUBROUTINE READ_STATE_FILE (STATEFILE, NCOMP, ICOMPS, NZ, NY, NX, &
                            ZLEV, TEMP, PRES, VAPMIXR, MIXR, CONC)
 ! Reads the state file with the profiles of level heights, temperature, 
 ! and pressures, mixing ratios of water vapor, and mass mixing ratio 
 ! and number concentration for several hydrometeor species.   
 !  Format:
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: STATEFILE
  INTEGER, INTENT(IN)  :: NCOMP, ICOMPS(NCOMP), NZ, NY, NX
  REAL,    INTENT(OUT) :: ZLEV(NZ,NY,NX), TEMP(NZ,NY,NX), PRES(NZ,NY,NX)
  REAL,    INTENT(OUT) :: VAPMIXR(NZ,NY,NX)
  REAL,    INTENT(OUT) :: MIXR(NCOMP,NZ,NY,NX), CONC(NCOMP,NZ,NY,NX)
  INTEGER :: NX2,NY2,NZ2, NC, IX,IY,IZ, IX1,IY1,IZ1,  J
  REAL    :: Z, P, T, V, M(7), C(7)

  OPEN (UNIT=1, FILE=STATEFILE, STATUS='OLD')
  READ (1,*) NX2, NY2, NZ2, NC
  IF (NCOMP > NC) THEN
    PRINT *, 'Not enough hydrometeor species in state file:'
    PRINT *, STATEFILE
    STOP
  ENDIF
  READ (1,*) 
  READ (1,*) 
  READ (1,*) 
  DO IX1 = 1, NX
    DO IY1 = 1, NY
      DO IZ1 = 1, NZ
        READ (1,*) IX, IY, IZ, Z, P, T, V, (M(J), C(J), J=1, NC)
        ZLEV(IZ,IY,IX) = Z
        PRES(IZ,IY,IX) = P
        TEMP(IZ,IY,IX) = T 
        VAPMIXR(IZ,IY,IX) = V
        DO J = 1, NCOMP
          MIXR(J,IZ,IY,IX) = M(ICOMPS(J))
          CONC(J,IZ,IY,IX) = C(ICOMPS(J))
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  CLOSE (1)
END SUBROUTINE READ_STATE_FILE



SUBROUTINE WRITE_STATE_FILE (STATEFILE, NCOMP, NZ, NY, NX, &
                             ZLEV, TEMP, PRES, ADVAPMIXR, ADMIXR, ADCONC)
 ! Writes the state file with the profiles of level heights, temperature, 
 ! and pressures, and the layer adjoint water vapor mixing ratio,
 ! and adjoint mass mixing ratio and adjoint number concentration for 
 ! several hydrometeor species.
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: STATEFILE
  INTEGER, INTENT(IN) :: NCOMP, NZ, NY, NX
  REAL,    INTENT(IN) :: ZLEV(NZ,NY,NX)
  REAL,    INTENT(IN) :: TEMP(NZ,NY,NX), PRES(NZ,NY,NX)
  REAL,    INTENT(IN) :: ADVAPMIXR(NZ,NY,NX)
  REAL,    INTENT(IN) :: ADMIXR(NCOMP,NZ,NY,NX), ADCONC(NCOMP,NZ,NY,NX)
  INTEGER :: IX,IY,IZ, J, L

  IF (NCOMP /= 2) THEN
    PRINT *, 'WRITE_STATE_FILE: Warning. Output header assumes cloud water and pristine ice hydrometeors'
  ENDIF

   ! Open the file and write it out
  OPEN (UNIT=2, FILE=STATEFILE, STATUS='UNKNOWN')
  WRITE (2,'(4(1X,I3),A)') NX, NY, NZ, NCOMP, '  ! Nx Ny Nz Ncomp'
  WRITE (2,*) '!                                  Adjoint     Adjoint Cloud Water    Adjoint Pristine Ice'
  WRITE (2,*) '! IX  IY  IZ  Height  Pres Temp    VapMixR      MixRat    NumConc      MixRat    NumConc  '
  WRITE (2,*) '!              (km)   (mb)  (K)    (g/kg)      (g/kg)     (#/kg)      (g/kg)     (#/kg)   '
  DO IX = 1, NX
    DO IY = 1, NY
      DO IZ = 1, NZ
        WRITE (2,'(3(1X,I3),1X,F7.4,1X,F6.1,1X,F5.1,1X,ES10.3,2(1X,2(1X,ES10.3)))') &
          IX, IY, IZ, ZLEV(IZ,IY,IX), PRES(IZ,IY,IX), TEMP(IZ,IY,IX), &
          ADVAPMIXR(IZ,IY,IX), (ADMIXR(J,IZ,IY,IX), ADCONC(J,IZ,IY,IX), J=1, NCOMP)
      ENDDO
    ENDDO
  ENDDO
  CLOSE (2)
END SUBROUTINE WRITE_STATE_FILE




SUBROUTINE WRITE_RADIANCE_FILE (OUTRADFILE, NOUT, NY, NX, MUOUT, PHIOUT, RADOUT)
  IMPLICIT NONE
  CHARACTER(LEN=*), INTENT(IN) :: OUTRADFILE
  INTEGER, INTENT(IN)  :: NOUT, NY, NX
  REAL,    INTENT(IN) :: MUOUT(NOUT), PHIOUT(NOUT), RADOUT(NOUT,NY,NX)
  INTEGER :: IX, IY, J

  OPEN (UNIT=1, FILE=OUTRADFILE, STATUS='UNKNOWN')
  WRITE (1,'(A)') '! SHDOMPPDA_DEMO radiance output'
  WRITE (1,'(A,8(2X,F6.4,2X))') '!    mu:', (MUOUT(J), J=1, NOUT)
  WRITE (1,'(A,8(2X,F5.1,3X))') '! IX  IY', (PHIOUT(J), J=1, NOUT)
  DO IX = 1, NX
    DO IY = 1, NY
      WRITE (1,'(2(1X,I3),8(2X,F8.4))') IX, IY, (RADOUT(J,IY,IX), J=1, NOUT)
    ENDDO
  ENDDO
  CLOSE (1)
END SUBROUTINE WRITE_RADIANCE_FILE






SUBROUTINE USER_INPUT (INSTATEFILE, OUTRADFILE, OUTSTATEFILE, &
                       NCOMP, ICOMPS, SCATTABFILES, SPECIESZRANGE, &
                       SRCTYPE, SOLARFLUX, SOLARMU, NOUT, MUOBS, PHIOBSREL, &
                       SFCTEMP, SFCALBEDO, WAVELEN, UNITS, &
                       NMU, NPHI, SPLITACC, SOLACC, &
                       PERT_RADOUT)
!       Obtains the input parameters for the program by writing prompts
!     to and reading responses from the terminal.  Echos the inputs to
!     make a useful log file when running non-interactively.  
  IMPLICIT NONE
  INTEGER, INTENT(OUT) :: NCOMP, ICOMPS(*), NMU, NPHI, NOUT
  REAL,    INTENT(OUT) :: SOLARFLUX, SOLARMU, MUOBS(*), PHIOBSREL(*)
  REAL,    INTENT(OUT) :: SFCTEMP, SFCALBEDO, WAVELEN
  REAL,    INTENT(OUT) :: SPLITACC, SOLACC
  REAL,    INTENT(OUT) :: SPECIESZRANGE(2,*)
  REAL,    INTENT(OUT) :: PERT_RADOUT(*)
  CHARACTER(LEN=1),  INTENT(OUT) :: SRCTYPE, UNITS
  CHARACTER(LEN=*), INTENT(OUT) :: INSTATEFILE, OUTRADFILE, OUTSTATEFILE
  CHARACTER(LEN=*), INTENT(OUT) :: SCATTABFILES(*)
  INTEGER :: I

  WRITE (*,*)
  WRITE (*,'(1X,A)') 'Input state file name'
  READ (*,'(A)') INSTATEFILE
  WRITE (*,*) INSTATEFILE

  WRITE (*,'(1X,A)') 'Output radiance file name'
  READ (*,'(A)') OUTRADFILE
  WRITE (*,*) OUTRADFILE

  WRITE (*,'(1X,A)') 'Output adjoint state file name'
  READ (*,'(A)') OUTSTATEFILE
  WRITE (*,*) OUTSTATEFILE

  WRITE (*,'(1X,A)') 'Number of hydrometeor components'
  READ (*,*) NCOMP
    WRITE (*,*) NCOMP

  DO I = 1, NCOMP
    WRITE (*,'(1X,I2,A)') I,'th component position in input state file'
    READ (*,*) ICOMPS(I)
    WRITE (*,*) ICOMPS(I)
  ENDDO

  DO I = 1, NCOMP
    WRITE (*,'(1X,I2,A)') I,'th component scattering table file name'
    READ (*,'(A)') SCATTABFILES(I)
    WRITE (*,*) SCATTABFILES(I)
  ENDDO

  DO I = 1, NCOMP
    WRITE (*,'(1X,I2,A)') I,'th particle species height range (min, max; km)'
    READ (*,*) SPECIESZRANGE(1:2,I)
    WRITE (*,*) SPECIESZRANGE(1:2,I)
  ENDDO

  WRITE (*,'(1X,A)') 'Thermal, solar, or both source (T, S, B)'
  READ (*,'(A)') SRCTYPE
    WRITE (*,*) SRCTYPE
  IF (SRCTYPE .EQ. 'S' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (*,'(1X,A)') 'Solar flux and direction (F, mu0)'
    READ (*,*) SOLARFLUX, SOLARMU
    WRITE (*,*) SOLARFLUX, SOLARMU
    SOLARMU = -ABS(SOLARMU)
  ENDIF
  IF (SRCTYPE .EQ. 'T' .OR. SRCTYPE .EQ. 'B') THEN
    WRITE (*,'(1X,A)') 'Surface temperature'
    READ (*,*) SFCTEMP
      WRITE (*,*) SFCTEMP
  ENDIF

  WRITE (*,'(1X,A)') 'Surface albedo'
  READ (*,*) SFCALBEDO
    WRITE (*,*) SFCALBEDO

  WRITE (*,'(1X,A)') 'Wavelength (micron)'
  READ (*,*) WAVELEN
    WRITE (*,*) WAVELEN

  WRITE (*,'(1X,A)') 'Output units (R=reflectance, T=brightness temperature K)'
  READ (*,'(A)') UNITS
    WRITE (*,*) UNITS

  WRITE (*,'(1X,A)') 'Number of discrete ordinates in mu and phi'
  READ (*,*) NMU, NPHI
    WRITE (*,*) NMU, NPHI

  WRITE (*,'(1X,A)') 'Layer splitting accuracy'
  READ (*,*) SPLITACC
    WRITE (*,*) SPLITACC

  WRITE (*,'(1X,A,A)') 'Solution accuracy'
  READ (*,*) SOLACC
    WRITE (*,*) SOLACC

  WRITE (*,*) 'Number of output radiance directions'
  READ (*,*) NOUT
    WRITE (*,*) NOUT
    
  WRITE (*,*) 'Cosine zenith angles for output radiance'
  READ (*,*) MUOBS(1:NOUT)
    WRITE (*,*) MUOBS(1:NOUT)
    
  WRITE (*,*) 'Azimuth angles for output radiance (degrees)'
  READ (*,*) PHIOBSREL(1:NOUT)
    WRITE (*,*) PHIOBSREL(1:NOUT)
    
  WRITE (*,*) 'Radiance perturbations for each radiance output'
  READ (*,*) PERT_RADOUT(1:NOUT)
    WRITE (*,*) PERT_RADOUT(1:NOUT)

  WRITE (*,*)
END SUBROUTINE USER_INPUT




