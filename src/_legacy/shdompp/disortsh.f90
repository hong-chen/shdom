program DISORTSH
    ! Driver program that calls DISORT to run plane-parallel radiative transfer
    ! with an SHDOMPP property file input.

    implicit none
    ! DISORT array sizes
    integer :: MAXCLY, MAXMOM, MAXULV, MAXUMU, MAXPHI
    ! MAXCLY   :  Max. number of computational layers
    ! MAXMOM   :  Max. number of phase function moment coefficients
    ! MAXULV   :  Max. number of user levels
    ! MAXUMU   :  Max. number of user polar angles
    ! MAXPHI   :  Max. number of user azimuth angles
    !    DISORT variables - see below
    logical :: LAMBER, PLANK, ONLYFL, PRNT(5), USRANG, USRTAU, DELTAMPLUS, DO_PSEUDO_SPHERE
    integer :: IBCND, NLYR, NMOM, NUMU, NSTR, NPHI, NTAU
    real    :: FBEAM, UMU0, PHI0, FISOT, TEMIS, TTEMP, BTEMP, ALBEDO
    real    :: WVNMLO, WVNMHI, ACCUR
    real(kind=4),parameter :: EARTH_RADIUS = 6371.0
    real(kind=4),dimension(:,:,:),allocatable :: RHOQ, RHOU
    real(kind=4),dimension(:),allocatable     :: EMUST, BEMST
    real(kind=4),dimension(:,:),allocatable   :: RHO_ACCURATE
    real, allocatable :: DTAUC(:), SSALB(:), PMOM(:, :), TEMPER(:)
    real, allocatable :: UTAU(:), PHI(:), UMU(:)
    real, allocatable :: RFLDIR(:), RFLDN(:), FLUP(:), DFDT(:), UAVG(:)
    real, allocatable :: UU(:, :, :), ALBMED(:), TRNMED(:), H_LYR(:)
    character(LEN=127) ::  HEADER

    ! SHDOMPP variables
    ! Number of output files and inputs parameter per file
    integer, parameter   :: MAXPAR = 100 ! Max array size for user input arrays
    real                 :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)

    integer :: NMU, NUMOUT, MAXLEG, NZOUT, NMUOUT, NPHIOUT
    integer :: NRAD, IRAD, I, J, K, N, L, NCASES, ICASE, NG, NZCKD, IG
    integer, allocatable :: NLEG(:), LOUT(:)
    logical :: KDIST
    real    :: SOLARFLUX, SOLFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO(2)
    real    :: ZDIF, MINDIF
    real, allocatable :: HEIGHTS(:), TAUP(:), ALBEDOP(:)
    real, allocatable :: DELG(:), ZCKD(:), KABS(:, :)
    real, allocatable :: FLUXUP(:), FLUXDN(:), MEANRAD(:), RADOUT(:, :, :)
    character(LEN=1) :: SRCTYPE
    character(LEN=80) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE

    call USER_INPUT(PROPFILE, CKDFILE, KDIST, NMU, SRCTYPE, &
                    SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                    FLUXFILE, RADFILE, MAXPAR, &
                    NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT)

    DELTAMPLUS = .false.
    DO_PSEUDO_SPHERE = .false.

    NMU = max(2, 2*int((NMU + 1)/2))
    NSTR = NMU

    ! Find out how many radiance directions we will be generating
    allocate (UMU(NMUOUT), PHI(NPHIOUT), RADOUT(NZOUT, NMUOUT, NPHIOUT))
    NPHI = NPHIOUT
    MAXPHI = NPHIOUT
    NUMU = NMUOUT
    MAXUMU = NMUOUT
    UMU(1:NMUOUT) = MUOUT(1:NMUOUT)
    PHI(1:NPHIOUT) = PHIOUT(1:NPHIOUT)

    allocate( RHOQ(NSTR/2, 0:NSTR/2, 0:(NSTR-1)), RHOU(NUMU, 0:NSTR/2, 0:(NSTR-1)), &
        EMUST(NUMU), BEMST(NSTR/2), RHO_ACCURATE(NUMU, NPHI) )
    RHOQ = 0.0; RHOU = 0.0; EMUST = 0.0; BEMST = 0.0; RHO_ACCURATE = 0.0

    call READ_PROPERTY_SIZE(PROPFILE, NLYR, MAXLEG)
    MAXLEG = max(MAXLEG, NMU)
    allocate (HEIGHTS(NLYR + 1), TEMPER(NLYR + 1), TAUP(NLYR), ALBEDOP(NLYR))
    allocate (DTAUC(NLYR), SSALB(NLYR), NLEG(NLYR))
    allocate (PMOM(0:MAXLEG, NLYR))
    MAXCLY = NLYR
    MAXULV = NLYR + 1
    MAXMOM = MAXLEG
    NMOM = MAXLEG
    allocate (RFLDIR(MAXULV), RFLDN(MAXULV), FLUP(MAXULV), DFDT(MAXULV))
    allocate (UTAU(MAXULV), UAVG(MAXULV), UU(MAXUMU, MAXULV, MAXPHI))
    allocate (ALBMED(MAXUMU), TRNMED(MAXUMU))
    allocate (LOUT(MAXULV))
    allocate (FLUXUP(MAXULV), FLUXDN(MAXULV), MEANRAD(MAXULV))
    allocate (H_LYR(0:NLYR))

    ! Read in the properties of the medium
    call READ_PROPERTIES(PROPFILE, NLYR, MAXLEG, &
                         HEIGHTS, TEMPER, TAUP, ALBEDOP, NLEG, PMOM)

    H_LYR(:) = HEIGHTS(:)

    do L = 0, MAXLEG ! Convert Legendre coefficients to DISORT definition
        PMOM(L, :) = PMOM(L, :)/(2*L + 1)
    end do

    ! Have disort output radiances
    ONLYFL = .false.

    ! Radiant quantities are to be returned at boundary of every
    ! computational layer.
    USRTAU = .false.

    ! Radiant quantities are to be returned at user polar angles.
    USRANG = .true.

    ! General top and bottom boundary conditions
    IBCND = 0

    if (SRCTYPE == 'S' .or. SRCTYPE == 'B') then
        FISOT = SKYRAD ! Intensity of top-boundary isotropic illumination.
        ! Intensity of incident parallel beam at top boundary (relative units)
        FBEAM = SOLARFLUX/SOLARMU
    else
        TTEMP = SKYRAD ! Temperature of top boundary (K).
        TEMIS = 1.0 ! Emissivity of top boundary.
        FBEAM = 0.0
    end if

    ! Polar angle cosine of incident beam (positive).
    UMU0 = abs(SOLARMU)
    ! Azimuth angle of incident beam (0 to 360 degrees)
    PHI0 = 0.0

    ! If true, isotropically reflecting (Lambertian) bottom boundary
    LAMBER = .true.
    ! Bottom-boundary albedo
    ALBEDO = SFCALB

    if (SRCTYPE == 'T' .or. SRCTYPE == 'B') then
        ! Temperature of bottom boundary (K)  (bottom emissivity from ALBEDO)
        BTEMP = SFCTEMP
        PLANK = .true. ! Include thermal emission for thermal source
        ! Wavenumbers (inv cm) of spectral interval of interest
        WVNMLO = WAVENO(1)
        WVNMHI = WAVENO(2)
    else
        BTEMP = 0.0
        PLANK = .false.
    end if

    ! Find levels of the desired radiances for output (pick closest)
    do I = 1, NZOUT
        MINDIF = 1.0e10
        do L = 1, NLYR + 1
            ZDIF = abs(ZOUT(I) - HEIGHTS(L))
            if (ZDIF < MINDIF) then
                MINDIF = ZDIF
                LOUT(I) = L
            end if
        end do
        ZOUT(I) = HEIGHTS(LOUT(I))
    end do

    !              CONTROL FLAGS

    ! Convergence criterion for azimuthal (Fourier cosine) series.
    ACCUR = 0.0001

    ! Array of LOGICAL print flags causing the following prints:
    !         L        quantities printed
    !        --        ------------------
    !         1        input variables (except PMOM)
    !         2        fluxes
    !         3        intensities at user levels and angles
    !         4        planar transmissivity and planar albedo
    !                  as a function solar zenith angle ( IBCND = 1 )
    !         5        phase function moments PMOM for each layer
    !                  ( only if PRNT(1) = TRUE, and only for layers
    !                  with scattering )
    PRNT(1) = .false.
    PRNT(2) = .false.
    PRNT(3) = .false.
    PRNT(4) = .false.
    PRNT(5) = .false.

    ! A 127 character header for prints, embedded in a DISORT banner;
    ! setting HEADER = '' will eliminate both the banner and the header.
    HEADER = ' '

    ! If doing a k-distribution then get the band info from the CKD file
    if (KDIST) then
        call READ_CKD_SIZE(CKDFILE, WAVENO, NG, NZCKD)
        allocate (DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG))
        call READ_CKD(CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)
        FBEAM = SOLFLUX*SOLARFLUX
    else
        NG = 1
        NZCKD = 0
        allocate (DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG))
        DELG(1) = 1.0
    end if

    !  WRITE (*,*) 'Number of times to run DISORT routine for timing'
    !  READ (*,*) NCASES
    !      WRITE (*,*) NCASES
    !  DO ICASE = 1, NCASES

    ! Loop over the k-distribution g's from low to high absorption
    do IG = 1, NG
        ! Add in the layer gaseous absorption if doing a k-distribution
        call ADD_IN_KDIST(NLYR, HEIGHTS, TAUP, ALBEDOP, &
                          NZCKD, ZCKD, KABS(:, IG), DTAUC, SSALB)

        !call DISORT(NLYR, DTAUC, SSALB, NMOM, PMOM, TEMPER, &
        !            WVNMLO, WVNMHI, USRTAU, NTAU, UTAU, NSTR, USRANG, &
        !            NUMU, UMU, NPHI, PHI, IBCND, &
        !            FBEAM, UMU0, PHI0, FISOT, LAMBER, ALBEDO, &
        !            BTEMP, TTEMP, TEMIS, PLANK, ONLYFL, &
        !            ACCUR, PRNT, HEADER, &
        !            MAXCLY, MAXULV, MAXUMU, MAXPHI, MAXMOM, &
        !            RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU, ALBMED, TRNMED)

        CALL DISORT( NLYR, NMOM, NSTR, NUMU, NPHI, NTAU,           &
                     USRANG, USRTAU, IBCND, ONLYFL, PRNT,          &
                     PLANK, LAMBER, DELTAMPLUS, DO_PSEUDO_SPHERE,  &
                     DTAUC, SSALB, PMOM, TEMPER, WVNMLO, WVNMHI,   &
                     UTAU, UMU0, PHI0, UMU, PHI, FBEAM,            &
                     FISOT, ALBEDO, BTEMP, TTEMP, TEMIS,           &
                     EARTH_RADIUS, H_LYR,                          &
                     RHOQ, RHOU, RHO_ACCURATE, BEMST, EMUST,       &
                     ACCUR,  HEADER,                               &
                     RFLDIR, RFLDN, FLUP, DFDT, UAVG, UU,          &
                     ALBMED, TRNMED )

        if (KDIST) then
            call SUM_KDIST_OUTPUT(IG, DELG(IG), NLYR, FLUP, FLUXUP, &
                                  RFLDIR + RFLDN, FLUXDN, UAVG, MEANRAD, &
                                  NZOUT, NMUOUT, NPHIOUT, LOUT, UU, RADOUT)
        else
            ! Put the upwelling and downwelling fluxes in the output array
            FLUXUP(:) = FLUP(:)
            FLUXDN(:) = RFLDIR(:) + RFLDN(:)
            MEANRAD(:) = UAVG(:)
            ! Put the radiances at desired levels in the output array
            do J = 1, NMUOUT
                do K = 1, NPHIOUT
                    do I = 1, NZOUT
                        RADOUT(I, J, K) = UU(J, LOUT(I), K)
                    end do
                end do
            end do
        end if

    end do ! End of k-distribution loop

!  ENDDO

    call OUTPUT_RESULTS(NLYR, NSTR, NG, HEIGHTS, PROPFILE, CKDFILE, SRCTYPE, &
                        SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                        FLUXUP, FLUXDN, MEANRAD, NZOUT, ZOUT, &
                        NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                        'F', FLUXFILE)

    if (NZOUT*NMUOUT*NPHIOUT > 0) then
        call OUTPUT_RESULTS(NLYR, NSTR, NG, HEIGHTS, PROPFILE, CKDFILE, SRCTYPE, &
                            SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                            FLUXUP, FLUXDN, MEANRAD, NZOUT, ZOUT, &
                            NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                            'R', RADFILE)
    end if

    deallocate (RFLDIR, RFLDN, FLUP, DFDT, UAVG, UTAU, UU, ALBMED, TRNMED)
    deallocate (DTAUC, SSALB, NLEG, PMOM, TEMPER)
    deallocate (UMU, PHI, RADOUT, FLUXUP, FLUXDN, MEANRAD)
    deallocate (HEIGHTS, TAUP, ALBEDOP, DELG, ZCKD, KABS, H_LYR)
    deallocate (RHOQ, RHOU, EMUST, BEMST, RHO_ACCURATE)
end

subroutine ADD_IN_KDIST(NLAY, HEIGHTP, TAUP, ALBEDOP, &
                        NZCKD, ZCKD, GASABS, TAUPM, ALBEDOPM)
    ! Adds the gaseous absorption in with the particle optical properties.
    ! The gaseous absorption optical depth is found for each of the NLAY
    ! optical property layers by interpolated the absorption linearly
    ! between the levels in ZCKD.
    integer, intent(IN) :: NLAY, NZCKD
    real, intent(IN) :: HEIGHTP(NLAY + 1), TAUP(NLAY), ALBEDOP(NLAY)
    real, intent(IN) :: ZCKD(NZCKD), GASABS(NZCKD)
    real, intent(OUT) :: TAUPM(NLAY), ALBEDOPM(NLAY)
    integer :: LAY, I
    real    :: GASTAU, Z, Z1, Z2, F, ABS

    if (NZCKD <= 0) then
        TAUPM(1:NLAY) = TAUP(1:NLAY)
        ALBEDOPM(1:NLAY) = ALBEDOP(1:NLAY)
    else
        if (HEIGHTP(1) > ZCKD(1) .or. HEIGHTP(NLAY + 1) < ZCKD(NZCKD)) then
            write (*, *) 'Layer heights outside of CKD file height range.'
            stop
        end if
        I = 1
        do LAY = 1, NLAY
            ! First, use trapezoidal integration to get optical depth of
            ! gaseous absorption assuming linear interpolation.
            GASTAU = 0.0
            Z1 = HEIGHTP(LAY)
            Z2 = HEIGHTP(LAY + 1)
            do while (ZCKD(I) >= Z1)
                I = I + 1
            end do
            I = I - 1
            Z = Z1
            do while (ZCKD(I) > Z2)
                F = (Z - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
                ABS = (1 - F)*GASABS(I) + F*GASABS(I + 1)
                GASTAU = GASTAU + 0.5*(ABS + GASABS(I + 1))*(Z - ZCKD(I + 1))
                Z = ZCKD(I + 1)
                I = I + 1
            end do
            I = I - 1
            F = (Z2 - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
            ABS = (1 - F)*GASABS(I) + F*GASABS(I + 1)
            GASTAU = GASTAU - 0.5*(ABS + GASABS(I + 1))*(Z2 - ZCKD(I + 1))
            ! Then, add the gaseous optical depth to the particle properties
            TAUPM(LAY) = TAUP(LAY) + GASTAU
            if (TAUPM(LAY) > 0.0) then
                ALBEDOPM(LAY) = ALBEDOP(LAY)*TAUP(LAY)/TAUPM(LAY)
            else
                ALBEDOPM(LAY) = 0.0
            end if
            ! print *, 'add_kdist: ',LAY,GASTAU,TAUPM(LAY),ALBEDOPM(LAY)
        end do
    end if
end subroutine ADD_IN_KDIST

subroutine SUM_KDIST_OUTPUT(IG, WT, NLAY, FLUP, FLUXUP, &
                            FLDN, FLUXDN, UAVG, MEANRAD, &
                            NZOUT, NMUOUT, NPHIOUT, LOUT, UU, RADOUT)
    ! Sums the output arrays over the k-distribution.
    integer, intent(IN) :: IG, NLAY, NZOUT, NMUOUT, NPHIOUT, LOUT(NLAY + 1)
    real, intent(IN) :: WT, UU(NMUOUT, NLAY + 1, NPHIOUT)
    real, intent(IN) :: FLUP(NLAY + 1), FLDN(NLAY + 1), UAVG(NLAY + 1)
    real, intent(OUT) :: FLUXUP(NLAY + 1), FLUXDN(NLAY + 1), MEANRAD(NLAY + 1)
    real, intent(OUT) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)
    integer :: I, J, K

    if (IG == 1) then
        FLUXUP(:) = 0.0
        FLUXDN(:) = 0.0
        MEANRAD(:) = 0.0
        RADOUT(:, :, :) = 0.0
    end if

    FLUXUP(:) = FLUXUP(:) + WT*FLUP(:)
    FLUXDN(:) = FLUXDN(:) + WT*FLDN(:)
    MEANRAD(:) = MEANRAD(:) + WT*UAVG(:)

    ! Add in the radiances at desired levels in the output array
    do J = 1, NMUOUT
        do K = 1, NPHIOUT
            do I = 1, NZOUT
                RADOUT(I, J, K) = RADOUT(I, J, K) + WT*UU(J, LOUT(I), K)
            end do
        end do
    end do
end subroutine SUM_KDIST_OUTPUT

subroutine USER_INPUT(PROPFILE, CKDFILE, KDIST, NMU, &
                      SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                      SFCTEMP, SFCALB, WAVENO, &
                      FLUXFILE, RADFILE, MAXPAR, &
                      NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT)
!       Obtains the input parameters for the program by writing prompts
!     to and reading responses from the terminal.  Echos the inputs to
!     make a useful log file when running non-interactively.  See the overall
!     program documentation for the list of input parameters.
    implicit none
    integer, intent(IN)  :: MAXPAR
    integer, intent(OUT) :: NMU, NZOUT, NMUOUT, NPHIOUT
    logical, intent(OUT) :: KDIST
    real, intent(OUT) :: SOLARFLUX, SOLARMU
    real, intent(OUT) :: SFCTEMP, SFCALB, SKYRAD, WAVENO(2)
    real, intent(OUT) :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
    character(LEN=1), intent(OUT) :: SRCTYPE
    character(LEN=80), intent(OUT) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE

    write (*, *)
    write (*, '(1X,A)') 'Atmospheric optical properties file name'
    read (*, '(A)') PROPFILE
    write (*, *) PROPFILE

    write (*, '(1X,A)') 'Correlated k-distribution file name (or NONE)'
    read (*, '(A)') CKDFILE
    write (*, *) CKDFILE
    KDIST = CKDFILE(1:2) .ne. 'NO'

    write (*, '(1X,A)') 'Number of discrete ordinates in -1 < mu < 1'
    read (*, *) NMU
    write (*, *) NMU

    write (*, '(1X,A)') 'Thermal, solar, or both source (T, S, B)'
    read (*, '(A)') SRCTYPE
    write (*, *) SRCTYPE
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        write (*, '(1X,A)') 'Solar flux and direction (F, mu0)'
        read (*, *) SOLARFLUX, SOLARMU
        write (*, *) SOLARFLUX, SOLARMU
        write (*, '(1X,A)') 'Isotropic sky radiance'
        read (*, *) SKYRAD
        write (*, *) SKYRAD
        SOLARMU = abs(SOLARMU)
    end if
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        write (*, '(1X,A)') 'Surface temperature'
        read (*, *) SFCTEMP
        write (*, *) SFCTEMP
    end if
    if (SRCTYPE .eq. 'T') then
        write (*, '(1X,A)') 'Sky temperature'
        read (*, *) SKYRAD
        write (*, *) SKYRAD
    end if

    write (*, '(1X,A)') 'Lambertian surface albedo'
    read (*, *) SFCALB
    write (*, *) SFCALB

    if (KDIST .or. SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        write (*, '(1X,A)') 'Wavenumber range (cm^-1)'
        read (*, *) WAVENO(1:2)
        write (*, *) WAVENO(1:2)
    end if

    write (*, '(1X,A,A)') 'Flux output file name'
    read (*, '(A)') FLUXFILE
    write (*, '(A)') FLUXFILE

    write (*, '(1X,A,A)') 'Radiance output file name'
    read (*, '(A)') RADFILE
    write (*, '(A)') RADFILE

    write (*, *) 'Number of optical depth levels for output radiance'
    read (*, *) NZOUT
    write (*, *) NZOUT
    write (*, *) 'Height levels for output radiance'
    read (*, *) ZOUT(1:NZOUT)
    write (*, *) ZOUT(1:NZOUT)

    write (*, *) 'Number of cosine zenith angles for output radiance'
    read (*, *) NMUOUT
    write (*, *) NMUOUT
    write (*, *) 'Cosine zenith angles for output radiance'
    read (*, *) MUOUT(1:NMUOUT)
    write (*, *) MUOUT(1:NMUOUT)

    write (*, *) 'Number of azimuth angles for output radiance'
    read (*, *) NPHIOUT
    write (*, *) NPHIOUT
    write (*, *) 'Azimuth angles for output radiance (degrees)'
    read (*, *) PHIOUT(1:NPHIOUT)
    write (*, *) PHIOUT(1:NPHIOUT)

    write (*, *)
end subroutine USER_INPUT

subroutine READ_PROPERTY_SIZE(PROPFILE, NLAY, MAXLEG)
!       Reads parts of the property file to get the maximum array sizes
!     needed for allocatable arrays.
    implicit none
    integer, intent(OUT) :: NLAY, MAXLEG
    character(LEN=*), intent(IN) :: PROPFILE
    integer :: I, NLEG
    real    :: TAU, Z, TTOP, TBOT, SSALB

    ! Open the file, get the number of layers and the max number Legendre terms
    open (UNIT=1, FILE=PROPFILE, STATUS='OLD')
    read (1, *) NLAY
    read (1, *) Z, TTOP
    MAXLEG = 1
    do I = 1, NLAY
        read (1, *) Z, TBOT, TAU, SSALB, NLEG
        MAXLEG = max(MAXLEG, NLEG)
    end do
    close (1)
end subroutine READ_PROPERTY_SIZE

subroutine READ_PROPERTIES(PROPFILE, NLAY, MAXLEG, &
                           HEIGHTP, TEMPP, TAUP, ALBEDOP, NLEGP, LEGENP)
    ! Reads the plane-parallel medium properties from the file into the
    ! property arrays.  The optical properties are uniform within each layer.
    ! The heights must decrease (top of atmosphere downward).
    implicit none
    integer, intent(IN)  :: NLAY, MAXLEG
    integer, intent(OUT) :: NLEGP(NLAY)
    real, intent(OUT) :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
    real, intent(OUT) :: HEIGHTP(NLAY + 1), TEMPP(NLAY + 1)
    character(LEN=*), intent(IN) :: PROPFILE
    integer :: NL, I, L

    ! Open the file, get the number of layers and the max number Legendre terms
    open (UNIT=1, FILE=PROPFILE, STATUS='OLD')
    read (1, *) NL
    read (1, *) HEIGHTP(1), TEMPP(1)
    do I = 1, NLAY
        LEGENP(:, I) = 0.0
        read (1, *) HEIGHTP(I + 1), TEMPP(I + 1), TAUP(I), ALBEDOP(I), &
            NLEGP(I), (LEGENP(L, I), L=1, NLEGP(I))
        if (HEIGHTP(I + 1) > HEIGHTP(I)) then
            write (*, *) 'Heights must decrease in property file.'
            stop
        end if
        LEGENP(0, I) = 1.0
    end do
    close (1)
end subroutine READ_PROPERTIES

subroutine READ_CKD_SIZE(CKDFILE, WAVENO, NG, NZCKD)
    ! Find the size of the k-distribution arrays for this wavenumber band.
    integer, intent(OUT) :: NG, NZCKD
    real, intent(IN)  :: WAVENO(2)
    character(LEN=*), intent(IN) :: CKDFILE
    integer :: NB, IB, JB, KB, N
    real    :: WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=CKDFILE, STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, SF, N
        if (abs(WAVENUM1 - WAVENO(1)) < 1.0 .and. &
            abs(WAVENUM2 - WAVENO(2)) < 1.0) then
            backspace (1)
            NG = N
            read (1, *) JB, WAVENUM1, WAVENUM2, SF, N
        end if
    end do
    if (JB == 0) then
        write (*, *) ' Wavenumber range not found in CKD file : ', &
            WAVENO(1), WAVENO(2), CKDFILE
        stop
    end if
    read (1, *) NZCKD
    close (1)
end subroutine READ_CKD_SIZE

subroutine READ_CKD(CKDFILE, WAVENO, NG, NZCKD, SOLFLUX, DELG, ZCKD, KABS)
    ! Reads the information appropriate for one band from a
    ! correlated k-distribution file.   The number of "g"'s are input in NG,
    ! and the number of levels are input in NZCKD. The wavenumber range (cm^-1)
    ! in the file must match the desired in WAVENO.  The solar flux
    ! is returned in SOLFLUX, the weights or delta g's in DELG, the
    ! levels from the top down in ZCKD, and the absorption coefficients
    ! for each g and level in KABS.
    integer, intent(IN) :: NG, NZCKD
    real, intent(INOUT) :: WAVENO(2)
    real, intent(OUT) :: SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG)
    character(LEN=*), intent(IN) :: CKDFILE
    integer :: NB, IB, JB, KB, IG, I, J, N
    real    :: WAVENUM1, WAVENUM2, SF

    open (UNIT=1, FILE=CKDFILE, STATUS='OLD')
    read (1, *)
    read (1, *) NB
    read (1, *)
    ! Read band information until the right one is found
    JB = 0
    do IB = 1, NB
        read (1, *) KB, WAVENUM1, WAVENUM2, SF, N
        if (abs(WAVENUM1 - WAVENO(1)) < 1.0 .and. &
            abs(WAVENUM2 - WAVENO(2)) < 1.0) then
            backspace (1)
            read (1, *) JB, WAVENO(1), WAVENO(2), SOLFLUX, &
                N, (DELG(IG), IG=1, NG)
        end if
    end do
    ! Read the Z level heights
    read (1, *) N
    read (1, *)
    read (1, *)
    do I = 1, NZCKD
        read (1, *) ZCKD(I)
    end do
    ! Skip over the irrelevant bands and read the absorption coefficients
    read (1, *)
    do I = 1, (JB - 1)*NZCKD
        read (1, *)
    end do
    do I = 1, NZCKD
        read (1, *) IB, J, (KABS(I, IG), IG=1, NG)
        do IG = 1, NG
            if (KABS(I, IG) < 0.0) then
                write (*, '(1X,A,A,I3,A,I3,A,I2)') 'Negative absorption in CKD file. ', &
                    'Band: ', IB, '   Level: ', I, '   k: ', IG
            end if
        end do
    end do
    close (1)
end subroutine READ_CKD

subroutine OUTPUT_RESULTS(NLAY, NMU, NG, HEIGHTS, PROPFILE, CKDFILE, SRCTYPE, &
                          SOLARFLUX, SOLARMU, SKYRAD, SFCTEMP, SFCALB, WAVENO, &
                          FLUXUP, FLUXDN, MEANRAD, NZOUT, ZOUT, &
                          NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
                          OUTTYPE, OUTFILE)
!       Writes the desired type of output file from the output fields
!     with a header giving the input parameters.
!     There are two types (OUTTYPE) of output: 'R' - for radiance,
!     'F' - for hemispheric/actinic flux.
!     See the overall program documentation for output parameters.
    implicit none
    integer, intent(IN) :: NLAY, NMU, NG
    integer, intent(IN) :: NZOUT, NMUOUT, NPHIOUT
    real, intent(IN) :: SOLARFLUX, SOLARMU, SFCTEMP, SFCALB, SKYRAD
    real, intent(IN) :: WAVENO(2)
    real, intent(IN) :: HEIGHTS(NLAY + 1)
    real, intent(IN) :: FLUXUP(NLAY + 1), FLUXDN(NLAY + 1), MEANRAD(NLAY + 1)
    real, intent(IN) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
    real, intent(IN) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)
    character(LEN=1), intent(IN) :: SRCTYPE, OUTTYPE
    character(LEN=*), intent(IN) :: PROPFILE, CKDFILE, OUTFILE
    integer :: I, J, K, L, N
    real    :: FOURPI
    character(LEN=32) :: SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

    FOURPI = 4*acos(-1.0)
    if (SRCTYPE == 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE == 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE == 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    SFCNAME = 'LAMBERTIAN'

    if (OUTTYPE == 'R') then
        UNITSNAME = 'WATTS/(M^2 STER)'
    else
        UNITSNAME = 'WATTS/(M^2)'
    end if
    if (OUTTYPE == 'R') OUTNAME = 'RADIANCE'
    if (OUTTYPE == 'F') OUTNAME = 'FLUX'

    open (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
    write (2, '(A)') '! DISORT Plane-parallel Radiative Transfer Output'
    write (2, '(2(A,I3))') '!  N_STREAMS=', NMU, '    N_LAYERS=', NLAY
    write (2, '(A,A32)') '!  PROPERTY_FILE=', PROPFILE
    write (2, '(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE, &
        '   NUM_G=', NG
    write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    write (2, '(A,A20)') '!  SURFACE_TYPE=', SFCNAME
    write (2, '(A,F9.7)') '!  SURFACE_ALBEDO=', SFCALB
    write (2, '(A,F8.3,A,E12.5)') '!  SURFACE_TEMP=', SFCTEMP, '  SKY_RAD=', SKYRAD
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        write (2, '(A,E13.6,A,F10.7)') &
            '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU
    end if
    write (2, '(A,A18,A,2F10.3)') '!  UNITS=', UNITSNAME, &
        '   WAVENUMBER_RANGE(CM^-1)=', WAVENO(1:2)
    write (2, '(A,A20)') '!  OUTPUT_TYPE=', OUTNAME

    if (OUTTYPE .eq. 'R') then
        ! Radiance output
        write (2, '(3(A,I3))') '!   NZOUT=', NZOUT, &
            '  NMUOUT=', NMUOUT, '  NPHIOUT=', NPHIOUT
        write (2, '(A)') '!    Z       MU     PHI     RADIANCE'
        do I = 1, NZOUT
            do J = 1, NMUOUT
                do K = 1, NPHIOUT
                    write (2, '(1X,F7.3,1X,F8.5,1X,F6.2,1X,E12.5)') &
                        ZOUT(I), MUOUT(J), PHIOUT(K), RADOUT(I, J, K)
                end do
            end do
        end do

    else if (OUTTYPE .eq. 'F') then
        ! Hemispheric flux output
        write (2, '(A)') '!    Z         F_UP        F_DOWN       F_ACTINIC'
        do I = 1, NLAY + 1
            write (2, '(1X,F7.3,3(2X,E12.5))') HEIGHTS(I), &
                FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
        end do

    end if
    close (2)
end subroutine OUTPUT_RESULTS
