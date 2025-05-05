! PROGRAM propgen
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

module propgen_comb

contains

    subroutine USER_INPUT(NSCATTAB, SCATTABFILES, SCATNUMS, POLTABS, PARFILE, &
                          MAXNEWPHASE, ASYMTOL, FRACPHASETOL, &
                          RAYWAVELEN, RAYSFCPRES, NZO, ZOTHER, TEMPOTHER, &
                          POLOUT, PROPFILE)
        ! Reads the input parameters from the standard input
        implicit none
        integer, intent(out) :: NSCATTAB, NZO, MAXNEWPHASE
        real, intent(out) :: RAYWAVELEN, RAYSFCPRES, ASYMTOL, FRACPHASETOL
        character(len=256), intent(out) :: PARFILE, PROPFILE
        integer, pointer     :: SCATNUMS(:)
        logical, pointer     :: POLTABS(:)
        real, pointer     :: ZOTHER(:), TEMPOTHER(:)
        character(len=256), pointer :: SCATTABFILES(:)
        character(len=1) :: POLOUT
        integer :: I
        logical :: VERBOSE=.FALSE.

        if (VERBOSE) write (*, *) 'Number of input scattering tables'
        read (*, *) NSCATTAB
        if (VERBOSE) write (*, *) NSCATTAB

        allocate (SCATTABFILES(NSCATTAB), SCATNUMS(NSCATTAB), POLTABS(NSCATTAB))
        if (NSCATTAB > 0) then
            if (VERBOSE) write (*, *) 'Name of each scattering table'
            do I = 1, NSCATTAB
                read (*, '(A)') SCATTABFILES(I)
                if (VERBOSE) write (*, '(A256)') SCATTABFILES(I)
            end do
            if (VERBOSE) write (*, *) 'Index or type of each scattering table'
            read (*, *) SCATNUMS(1:NSCATTAB)
            if (VERBOSE) write (*, '(20(1X,I2))') SCATNUMS(1:NSCATTAB)
            if (VERBOSE) write (*, *) 'Polarization flag (T or F) for each scattering table'
            read (*, *) POLTABS(1:NSCATTAB)
            if (VERBOSE) write (*, '(20(1X,L1))') POLTABS(1:NSCATTAB)
        end if

        if (VERBOSE) write (*, *) 'Input particle properties filename (or NONE for Rayleigh only)'
        read (*, '(A)') PARFILE
        if (VERBOSE) write (*, '(A256)') PARFILE

        if (VERBOSE) write (*, *) 'Maximum number of new phase functions created'
        read (*, *) MAXNEWPHASE
        if (VERBOSE) write (*, '(1X,I4)') MAXNEWPHASE

        if (VERBOSE) write (*, *) 'Tolerance in asymmetry parameter for creating a new phase function'
        read (*, *) ASYMTOL
        if (VERBOSE) write (*, '(1X,F6.4)') ASYMTOL

        if (VERBOSE) write (*, *) 'Maximum fractional error tolerance in value for creating a new phase function'
        read (*, *) FRACPHASETOL
        if (VERBOSE) write (*, '(1X,F6.3)') FRACPHASETOL

        if (VERBOSE) write (*, *) 'Wavelength for Rayleigh molecular scattering (micron)'
        read (*, *) RAYWAVELEN
        if (VERBOSE) write (*, *) RAYWAVELEN

        if (VERBOSE) write (*, *) 'Surface pressure for Rayleigh molecular scattering (mb or hPa)'
        read (*, *) RAYSFCPRES
        if (VERBOSE) write (*, *) RAYSFCPRES

        if (VERBOSE) write (*, *) 'Number of extra Z levels'
        read (*, *) NZO
        if (VERBOSE) write (*, '(1X,I3)') NZO
        allocate (ZOTHER(NZO), TEMPOTHER(NZO))
        do I = 1, NZO
            if (VERBOSE) write (*, *) 'Height (km) and temperature (K) :', I
            read (*, *) ZOTHER(I), TEMPOTHER(I)
            if (VERBOSE) write (*, '(1X,F7.3,1X,F6.2)') ZOTHER(I), TEMPOTHER(I)
        end do

        if (VERBOSE) write (*, *) 'Output polarized (P) or unpolarized (U) SHDOM property file'
        read (*, '(A)') POLOUT
        if (VERBOSE) write (*, *) POLOUT

        if (VERBOSE) write (*, *) 'Output property file name'
        read (*, '(A)') PROPFILE
        if (VERBOSE) write (*, '(A256)') PROPFILE
    end subroutine USER_INPUT

    subroutine READ_SCAT_TABLE_SIZE(SCATTABFILE, POLTAB, WAVELEN1, WAVELEN2, &
                                    NRETAB, MAXNCOEF)
        ! Reads the scattering table file to find the number of effective radii,
        ! the maximum order of the Wigner d-function series phase functions, and
        ! also the wavelength range.
        implicit none
        character(len=256), intent(in) :: SCATTABFILE
        logical, intent(in) :: POLTAB
        integer, intent(out) :: NRETAB, MAXNCOEF
        real, intent(out) :: WAVELEN1, WAVELEN2
        integer :: I, J, L, NCOEF
        real    :: REFF, EXT, SSALB, CHI

        open (UNIT=1, FILE=trim(adjustl(SCATTABFILE)), STATUS='OLD')
        read (1, *)
        read (1, *) WAVELEN1, WAVELEN2
        read (1, *)
        read (1, *)
        read (1, *)
        read (1, *) NRETAB
        MAXNCOEF = 1
        do I = 1, NRETAB
            read (1, *) REFF, EXT, SSALB, NCOEF
            MAXNCOEF = max(MAXNCOEF, NCOEF)
            if (POLTAB) then
                read (1, *) J, (CHI, L=0, NCOEF)
                read (1, *) J, (CHI, L=0, NCOEF)
                read (1, *) J, (CHI, L=0, NCOEF)
                read (1, *) J, (CHI, L=0, NCOEF)
                read (1, *) J, (CHI, L=0, NCOEF)
                read (1, *) J, (CHI, L=0, NCOEF)
            else
                read (1, *) (CHI, L=0, NCOEF)
            end if
        end do
        close (1)
    end subroutine READ_SCAT_TABLE_SIZE

    subroutine READ_SCAT_TABLE(SCATTABFILE, POLTAB, NPHASEPOL, NRETAB, RETAB, &
                               EXTINCT, SSALB, MAXCOEF, NCOEF, WIGCOEF)
        ! Reads the table of scattering properties as a function of effective
        ! radius.  POLTAB=.true. means a polarized scattering table with Wigner
        ! d-function coefficients for six phase matrix elements instead of one.
        ! The order of the six elements is the four diagonal ones (a1, a2, a3, a4)
        ! followed by the IQ and UV coefficients (b1, b2).
        implicit none
        character(len=256), intent(in) :: SCATTABFILE
        integer, intent(in)  :: NPHASEPOL, NRETAB, MAXCOEF
        logical, intent(in)  :: POLTAB
        integer, intent(out) :: NCOEF(:)
        real, intent(out) :: RETAB(NRETAB), EXTINCT(NRETAB), SSALB(NRETAB)
        real, intent(out) :: WIGCOEF(NPHASEPOL, 0:MAXCOEF, NRETAB)
        integer :: I, J, L
        real    :: CHI

        open (UNIT=1, FILE=trim(adjustl(SCATTABFILE)), STATUS='OLD')
        read (1, *)
        read (1, *)
        read (1, *)
        read (1, *)
        read (1, *)
        read (1, *)
        do I = 1, NRETAB
            read (1, *) RETAB(I), EXTINCT(I), SSALB(I), NCOEF(I)
            if (POLTAB) then
                read (1, *) J, (WIGCOEF(1, L, I), L=0, NCOEF(I))
                if (NPHASEPOL > 1) then
                    read (1, *) J, (WIGCOEF(2, L, I), L=0, NCOEF(I))
                    read (1, *) J, (WIGCOEF(3, L, I), L=0, NCOEF(I))
                    read (1, *) J, (WIGCOEF(4, L, I), L=0, NCOEF(I))
                    read (1, *) J, (WIGCOEF(5, L, I), L=0, NCOEF(I))
                    read (1, *) J, (WIGCOEF(6, L, I), L=0, NCOEF(I))
                else
                    read (1, *) J, (CHI, L=0, NCOEF(I))
                    read (1, *) J, (CHI, L=0, NCOEF(I))
                    read (1, *) J, (CHI, L=0, NCOEF(I))
                    read (1, *) J, (CHI, L=0, NCOEF(I))
                    read (1, *) J, (CHI, L=0, NCOEF(I))
                end if
            else
                read (1, *) (WIGCOEF(1, L, I), L=0, NCOEF(I))
                if (NPHASEPOL > 1) then
                    WIGCOEF(2, :, I) = 0.0
                    WIGCOEF(3, :, I) = 0.0
                    WIGCOEF(4, :, I) = 0.0
                    WIGCOEF(5, :, I) = 0.0
                    WIGCOEF(6, :, I) = 0.0
                end if
            end if
            if (abs(WIGCOEF(1, 0, I) - 1.0) > 0.0001) then
                write (*, *) 'Error [READ_SCAT_TABLE]: Incorrect Legendre series; chi_0 is not 1.'
                stop
            end if
            if (I > 1 .and. RETAB(I) <= RETAB(I - 1)) then
                write (*, *) 'Error [READ_SCAT_TABLE]: Effective radius not increasing in table:', SCATTABFILE
                stop
            end if
        end do
        close (1)
    end subroutine READ_SCAT_TABLE

    subroutine MAKE_TYPE_CONV_TABLE(NTYPES, PARTYPES, MAXTYPENUM, TYPECONV)
        ! Makes a table to convert from particle type number to the index in
        ! the arrays - basically the inverse of PARTYPES.
        implicit none
        integer, intent(in) :: NTYPES, PARTYPES(NTYPES), MAXTYPENUM
        integer, intent(out) :: TYPECONV(0:MAXTYPENUM)
        integer :: I

        TYPECONV(:) = 0
        do I = 1, NTYPES
            TYPECONV(PARTYPES(I)) = I
        end do
    end subroutine MAKE_TYPE_CONV_TABLE

    subroutine READ_PARTICLE_FILE_SIZE(PARFILE, NX, NY, NZP)
        implicit none
        character(len=256), intent(in) :: PARFILE
        integer, intent(out) :: NX, NY, NZP

        open (UNIT=2, FILE=trim(adjustl(PARFILE)), STATUS='OLD')
        read (2, *)
        read (2, *) NX, NY, NZP
        close (2)
    end subroutine READ_PARTICLE_FILE_SIZE

    subroutine READ_PARTICLE_FILE(PARFILE, NX, NY, NZP, NSCATTAB, PARTYPES, &
                                  DELX, DELY, ZPAR, TEMPPAR, &
                                  NCOMP, PTYPE, MASSCONT, REFF, CKDFILE, ISEXT)
        ! Reads the particle physical properties file.  The file contains a header
        ! containing the array size (NX,NY,NZP), heights (ZPAR), and
        ! temperature profile (TEMPPAR).  Then each row has format:
        !    IX IY IZ  Ncomp  Ptype1 masscont1 Reff1 ... PtypeN masscontN ReffN
        ! where Ncomp is the number of components, Ptype is the particle type,
        ! masscont is the mass content (g/m^3), and Reff is the particle effective
        ! radius (micron).
        !  Also reads the old 2 parameter LWC files with LWC and Reff at each point.
        implicit none
        character(len=256), intent(in) :: PARFILE
        integer, intent(in) :: NX, NY, NZP, NSCATTAB, PARTYPES(NSCATTAB)
        real, intent(out) :: DELX, DELY, ZPAR(NZP), TEMPPAR(NZP)
        integer, intent(out) :: NCOMP(NZP, NY, NX), PTYPE(NSCATTAB, NZP, NY, NX)
        real, intent(out) :: MASSCONT(NSCATTAB, NZP, NY, NX), REFF(NSCATTAB, NZP, NY, NX)
        character(len=256), intent(out) :: CKDFILE
        logical, intent(out) :: ISEXT
        integer :: I, FILEKIND, IX, IY, IZ, NC, PT(NSCATTAB)
        real    :: MASS(NSCATTAB), RE(NSCATTAB)

        character(len=256) :: TMP_CHAR, OUTFILE_DATA, DIMS_CHAR
        character(len=20) :: DATA_POSTFIX
        integer :: NXX, NYY, NZZ, NPP

        ISEXT = .FALSE.
        CKDFILE = 'NONE'

        ! Initialize the output arrays (in case there are missing grid points)
        NCOMP(:, :, :) = 0
        PTYPE(:, :, :, :) = 0
        MASSCONT(:, :, :, :) = 0.0
        REFF(:, :, :, :) = 0.0

        ! Read in the header
        open (UNIT=2, FILE=trim(adjustl(PARFILE)), STATUS='OLD')
        read (2, *) FILEKIND
        read (2, *) ! NX, NY, NZP
        read (2, *) DELX, DELY
        read (2, *) (ZPAR(IZ), IZ=1, NZP)
        read (2, *) (TEMPPAR(IZ), IZ=1, NZP)

        if (FILEKIND == 2) then
            if (NSCATTAB /= 1) then
                write (*, *) 'Error [READ_PARTICLE_FILE]: Must have only one scattering table to use 2 parameter LWC file.'
                stop
            end if

            read (2, '(A)') TMP_CHAR ! skip one line

            if (TMP_CHAR(1:1) .eq. '!') then

                read (2, '(A)') DATA_POSTFIX
                read (2, '(A)') DIMS_CHAR

                read (DIMS_CHAR( 3: 8), '(I6)') NXX
                read (DIMS_CHAR(10:15), '(I6)') NYY
                read (DIMS_CHAR(17:22), '(I6)') NZZ
                read (DIMS_CHAR(24:29), '(I6)') NPP

                read (2, '(A)') CKDFILE
                CKDFILE = trim(adjustl(CKDFILE(2:)))

                OUTFILE_DATA = trim(adjustl(PARFILE))//trim(adjustl(DATA_POSTFIX(2:)))

                open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='OLD',&
                    &ACCESS='STREAM', FORM='UNFORMATTED')
                read (88) (((MASSCONT(1, IZ, IY, IX), IZ=1,NZZ), IY=1,NYY), IX=1,NXX)
                read (88) (((    REFF(1, IZ, IY, IX), IZ=1,NZZ), IY=1,NYY), IX=1,NXX)
                close (88)

                NCOMP = 1
                PTYPE(1, :, :, :) = PARTYPES(1)

                ISEXT = .True.

            else

                backspace(2)

                ! Read in the data
                do while (.true.)
                    read (2, *, end=190) IX, IY, IZ, MASS(1), RE(1)
                    if (IX >= 1 .and. IX <= NX .and. IY >= 1 .and. IY <= NY .and. &
                        IZ >= 1 .and. IZ <= NZP) then
                        NCOMP(IZ, IY, IX) = 1
                        PTYPE(1, IZ, IY, IX) = PARTYPES(1)
                        MASSCONT(1, IZ, IY, IX) = MASS(1)
                        REFF(1, IZ, IY, IX) = RE(1)
                    end if
                end do

            end if

        else if (FILEKIND == 3) then
            ! Read in the data
            do while (.true.)
                read (2, *, end=190) IX, IY, IZ, NC, &
                    (PT(I), MASS(I), RE(I), I=1, min(NC, NSCATTAB))
                if (NC > NSCATTAB) then
                    write (*, *) 'Error [READ_PARTICLE_FILE]: More particle components than scattering tables.'
                    stop
                end if
                if (IX >= 1 .and. IX <= NX .and. IY >= 1 .and. IY <= NY .and. &
                    IZ >= 1 .and. IZ <= NZP) then
                    NCOMP(IZ, IY, IX) = NC
                    PTYPE(1:NC, IZ, IY, IX) = PT(1:NC)
                    MASSCONT(1:NC, IZ, IY, IX) = MASS(1:NC)
                    REFF(1:NC, IZ, IY, IX) = RE(1:NC)
                end if
            end do

        else
            write (*, *) 'Error [READ_PARTICLE_FILE]: Must be type 2 or 3 particle properties file.'
            stop
        end if

190     continue
        close (2)
    end subroutine READ_PARTICLE_FILE

    subroutine ORGANIZE_LEVELS(NZP, ZPAR, TEMPPAR, NZO, ZOTHER, TEMPOTHER, &
                               NZT, ZLEVELS, TEMP, IPARLEV)
        ! Combines the particle and extra levels, returning the resulting
        ! heights (ZLEVELS) and temperature profile (TEMP), and the
        ! pointers to the particle levels (IPARLEV).
        implicit none
        integer, intent(in)  :: NZP, NZO, NZT
        integer, intent(out) :: IPARLEV(NZT)
        real, intent(in)  :: ZPAR(NZP), TEMPPAR(NZP), ZOTHER(NZO), TEMPOTHER(NZO)
        real, intent(out) :: ZLEVELS(NZT), TEMP(NZT)
        integer :: J, K

        ! Put the particle and other Z levels (and temperatures) in one array
        ZLEVELS(1:NZP) = ZPAR(:)
        TEMP(1:NZP) = TEMPPAR(:)
        ZLEVELS(NZP + 1:NZT) = ZOTHER(:)
        TEMP(NZP + 1:NZT) = TEMPOTHER(:)

        ! Warn if extra levels are within the particle layer
        if (NZP > 0) then
            do J = 1, NZO
                if (ZOTHER(J) >= ZPAR(1) .and. ZOTHER(J) <= ZPAR(NZP)) then
                    write (*, '(1X,A,F7.3)') 'Warning [ORGANIZE_LEVELS]: Other Z level inside range of particle levels:', ZOTHER(J)
                    write (*, *) '         The other Z levels have only Rayleigh scattering.'
                end if
            end do
        end if

        ! Sort the levels
        call SSORT(ZLEVELS, TEMP, NZT, 2)

        ! Check to see if two Z levels are same
        do K = 1, NZT - 1
            if (ZLEVELS(K) == ZLEVELS(K + 1)) then
                write (*, *) 'Error [ORGANIZE_LEVELS]: Same Z levels: ', ZLEVELS(K)
                stop
            end if
        end do
        ! Make the particle level pointers
        if (NZP == 0) then
            IPARLEV(:) = 0
        else
            J = 1
            do K = 1, NZT
                IPARLEV(K) = 0
                if (ZLEVELS(K) == ZPAR(J)) then
                    IPARLEV(K) = J
                    J = min(NZP, J + 1)
                end if
            end do
        end if
    end subroutine ORGANIZE_LEVELS

    subroutine SSORT(X, Y, N, KFLAG)
        ! SSORT sorts array X and optionally makes the same interchanges in
        ! array Y.  The array X may be sorted in increasing order or
        ! decreasing order.  A slightly modified quicksort algorithm is used.
        !
        ! Description of Parameters
        !   X - array of values to be sorted   (usually abscissas)
        !   Y - array to be (optionally) carried along
        !   N - number of values in array X to be sorted
        !   KFLAG - control parameter
        !         =  2  means sort X in increasing order and carry Y along.
        !         =  1  means sort X in increasing order (ignoring Y)
        !         = -1  means sort X in decreasing order (ignoring Y)
        !         = -2  means sort X in decreasing order and carry Y along.

        implicit none
        integer KFLAG, N
        real X(*), Y(*)
        real R, T, TT, TTY, TY
        integer I, IJ, J, K, KK, L, M, NN
        integer IL(21), IU(21)

        NN = N
        if (NN .lt. 1) then
            stop 'Error [SSORT]: The number of values to be sorted is not positive.'
        end if
        KK = abs(KFLAG)
        if (KK .ne. 1 .and. KK .ne. 2) then
            stop 'Error [SSORT]: The sort control parameter, K, is not 2, 1, -1, or -2.'
        end if

        ! Alter array X to get decreasing order if needed
        if (KFLAG .le. -1) then
            do 10 I = 1, NN
                X(I) = -X(I)
10          continue
        end if

        if (KK .eq. 2) GO TO 100

        ! Sort X only
        M = 1
        I = 1
        J = NN
        R = 0.375e0

20      if (I .eq. J) GO TO 60
        if (R .le. 0.5898437e0) then
            R = R + 3.90625e-2
        else
            R = R - 0.21875e0
        end if
30      K = I

        ! Select a central element of the array and save it in location T
        IJ = I + int((J - I)*R)
        T = X(IJ)
        ! If first element of array is greater than T, interchange with T
        if (X(I) .gt. T) then
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
        end if
        L = J
        ! If last element of array is less than than T, interchange with T
        if (X(J) .lt. T) then
            X(IJ) = X(J)
            X(J) = T
            T = X(IJ)
            ! If first element of array is greater than T, interchange with T
            if (X(I) .gt. T) then
                X(IJ) = X(I)
                X(I) = T
                T = X(IJ)
             end if
        end if
        ! Find an element in the second half of the array which is smaller than T
40      L = L - 1
        if (X(L) .gt. T) GO TO 40
        ! Find an element in the first half of the array which is greater than T
50      K = K + 1
        if (X(K) .lt. T) GO TO 50
        ! Interchange these elements
        if (K .le. L) then
            TT = X(L)
            X(L) = X(K)
            X(K) = TT
            GO TO 40
        end if
        ! Save upper and lower subscripts of the array yet to be sorted
        if (L - I .gt. J - K) then
            IL(M) = I
            IU(M) = L
            I = K
            M = M + 1
        else
            IL(M) = K
            IU(M) = J
            J = L
            M = M + 1
        end if
        GO TO 70
        ! Begin again on another portion of the unsorted array
60      M = M - 1
        if (M .eq. 0) GO TO 190
        I = IL(M)
        J = IU(M)
70      if (J - I .ge. 1) GO TO 30
        if (I .eq. 1) GO TO 20
        I = I - 1
80      I = I + 1
        if (I .eq. J) GO TO 60
        T = X(I + 1)
        if (X(I) .le. T) GO TO 80
        K = I
90      X(K + 1) = X(K)
        K = K - 1
        if (T .lt. X(K)) GO TO 90
        X(K + 1) = T
        GO TO 80
        ! Sort X and carry Y along
100     M = 1
        I = 1
        J = NN
        R = 0.375e0
110     if (I .eq. J) GO TO 150
        if (R .le. 0.5898437e0) then
            R = R + 3.90625e-2
        else
            R = R - 0.21875e0
        end if
120     K = I
        ! Select a central element of the array and save it in location T
        IJ = I + int((J - I)*R)
        T = X(IJ)
        TY = Y(IJ)
        ! If first element of array is greater than T, interchange with T
        if (X(I) .gt. T) then
            X(IJ) = X(I)
            X(I) = T
            T = X(IJ)
            Y(IJ) = Y(I)
            Y(I) = TY
            TY = Y(IJ)
        end if
        L = J
        ! If last element of array is less than T, interchange with T
        if (X(J) .lt. T) then
            X(IJ) = X(J)
            X(J) = T
            T = X(IJ)
            Y(IJ) = Y(J)
            Y(J) = TY
            TY = Y(IJ)
            ! If first element of array is greater than T, interchange with T
            if (X(I) .gt. T) then
                X(IJ) = X(I)
                X(I) = T
                T = X(IJ)
                Y(IJ) = Y(I)
                Y(I) = TY
                TY = Y(IJ)
            end if
        end if
        ! Find an element in the second half of the array which is smaller than T
130     L = L - 1
        if (X(L) .gt. T) GO TO 130
        ! Find an element in the first half of the array which is greater than T
140     K = K + 1
        if (X(K) .lt. T) GO TO 140
        ! Interchange these elements
        if (K .le. L) then
             TT = X(L)
             X(L) = X(K)
             X(K) = TT
             TTY = Y(L)
             Y(L) = Y(K)
             Y(K) = TTY
             GO TO 130
        end if
        ! Save upper and lower subscripts of the array yet to be sorted
        if (L - I .gt. J - K) then
            IL(M) = I
            IU(M) = L
            I = K
            M = M + 1
        else
            IL(M) = K
            IU(M) = J
            J = L
            M = M + 1
        end if
        GO TO 160

        ! Begin again on another portion of the unsorted array
150     M = M - 1
        if (M .eq. 0) GO TO 190
        I = IL(M)
        J = IU(M)
160     if (J - I .ge. 1) GO TO 120
        if (I .eq. 1) GO TO 110
        I = I - 1
170     I = I + 1
        if (I .eq. J) GO TO 150
        T = X(I + 1)
        TY = Y(I + 1)
        if (X(I) .le. T) GO TO 170
        K = I
180     X(K + 1) = X(K)
        Y(K + 1) = Y(K)
        K = K - 1
        if (T .lt. X(K)) GO TO 180
        X(K + 1) = T
        Y(K + 1) = TY
        GO TO 170
!       Clean up
190     if (KFLAG .le. -1) then
            do 200 I = 1, NN
                X(I) = -X(I)
200         continue
        end if
    end subroutine SSORT

    subroutine RAYLEIGH_EXTINCT(NZT, ZLEVELS, TEMP, RAYSFCPRES, RAYLCOEF, EXTRAYL)
        ! Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
        ! from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
        ! a linear lapse rate between levels to compute the pressure at
        ! each level starting with the SFCPRES [mb].  The Rayleigh extinction is
        ! proportional to air density, with the coefficient RAYLCOEF in [K/(mb km)].
        implicit none
        integer, intent(in) :: NZT
        real, intent(in) :: ZLEVELS(NZT), TEMP(NZT), RAYSFCPRES, RAYLCOEF
        real, intent(out) :: EXTRAYL(NZT)
        integer :: I
        real    :: PRES, LAPSE, TS, DZ

        ! Find surface pressure by integrating hydrostatic relation
        !  for a dry atmosphere up to surface height.
        PRES = RAYSFCPRES
        TS = TEMP(1)
        LAPSE = 6.5*0.001
        PRES = PRES*(TS/(TS + LAPSE*ZLEVELS(1)*1000.0))**(9.80665/(287.052874*LAPSE))

        ! Use layer mean temperature to compute fractional pressure change.
        do I = 1, NZT - 1
            EXTRAYL(I) = RAYLCOEF*PRES/TEMP(I)
            DZ = 1000.0*(ZLEVELS(I + 1) - ZLEVELS(I))
            LAPSE = (TEMP(I) - TEMP(I + 1))/DZ
            if (abs(LAPSE) .gt. 0.00001) then
                PRES = PRES*(TEMP(I + 1)/TEMP(I))**(9.80665/(287.052874*LAPSE))
            else
                PRES = PRES*exp(-9.80665*DZ/(287.052874*TEMP(I)))
            end if
        end do
        EXTRAYL(NZT) = RAYLCOEF*PRES/TEMP(NZT)

    end subroutine RAYLEIGH_EXTINCT

    subroutine RAYLEIGH_EXTINCT_FROM_CKD(CKDFILE, NZT, ZLEVELS, EXTRAYL)
        ! Computes the molecular Rayleigh extinction profile EXTRAYL [/km]
        ! from the temperature profile TEMP [K] at ZLEVELS [km].  Assumes
        ! a linear lapse rate between levels to compute the pressure at
        ! each level starting with the SFCPRES [mb].  The Rayleigh extinction is
        ! proportional to air density, with the coefficient RAYLCOEF in [K/(mb km)].
        implicit none
        character(len=*), intent(in) :: CKDFILE
        integer, intent(in) :: NZT
        real, intent(in) :: ZLEVELS(NZT)
        real, intent(out) :: EXTRAYL(NZT)
        integer :: I, IZ, IL, IU, IM, NZCKD
        real, allocatable :: ZCKD(:), EXTRAY(:)
        double precision :: W0

        ! read rayleigh extinction from CKDFILE
        !╭────────────────────────────────────────────────╮!
        open (UNIT=1, FILE=trim(adjustl(CKDFILE)), STATUS='OLD')
        read (1, *)
        read (1, *)
        read (1, *)
        read (1, *)
        read (1, *) NZCKD

        allocate(ZCKD(NZCKD), EXTRAY(NZCKD))

        read (1, *)
        read (1, *)
        do IZ = 1, NZCKD
            read (1, *) ZCKD(IZ), EXTRAY(IZ)
        end do

        close(1)
        !╰────────────────────────────────────────────────╯!

        ! interpolate
        !╭────────────────────────────────────────────────╮!
        EXTRAYL = 0.0

        do IZ = 1, NZT
            if (NZCKD .gt. 0) then
                IL = 1
                IU = NZCKD
                do while (IU - IL .gt. 1)
                    IM = (IU + IL)/2
                    if (ZLEVELS(IZ) .le. ZCKD(IM)) then
                        IL = IM
                    else
                        IU = IM
                    end if
                end do
                I = min(max(IL, 1), NZCKD - 1)
                W0 = (ZLEVELS(IZ) - ZCKD(I))/(ZCKD(I + 1) - ZCKD(I))
                W0 = min(max(W0, 0.0d0), 1.0d0)
                EXTRAYL(IZ) = (1.0 - W0)*EXTRAY(I) + W0*EXTRAY(I + 1)
            else
                EXTRAYL(IZ) = 0.0
            end if
        end do
        !╰────────────────────────────────────────────────╯!

        deallocate(ZCKD, EXTRAY)

    end subroutine RAYLEIGH_EXTINCT_FROM_CKD

    subroutine RAYLEIGH_PHASE_FUNCTION(WAVELEN, NPHASEPOL, MAXCOEF, &
                                       NCOEF, RAYWIGCOEF)
        ! Returns the Wigner d-function coefficients for either the Rayleigh
        ! scalar phase function or polarized phase matrix for molecular
        ! scattering by air.  Includes the wavelength depolarization factor.
        ! From Mishchenko's book "Multiple scattering of light by particles:
        ! Radiative Transfer and Coherent Backscattering", Cambridge, 2006.
        ! Thanks to Adrian Doicu.   WAVELEN is the wavelength in microns.
        ! NPHASEPOL=1 for the scalar case, and 6 for the phase matrix case.
        ! Reference:
        ! Section 5.1.1.3, Eq(5.13 - 5.23) in Stefano Corradini's PHD thesis:
        !   Aerosol extinction coefficient profile retrieval procedure for satellite
        !   measurements in multiple scattering atmosphere
        implicit none
        integer, intent(in) :: NPHASEPOL, MAXCOEF
        real, intent(in) :: WAVELEN
        integer, intent(out) :: NCOEF
        real, intent(out) :: RAYWIGCOEF(NPHASEPOL, 0:MAXCOEF)
        double precision, parameter :: AKING = 1.0469541d0, &
            BKING = 3.2503153d-04, CKING = 3.8622851d-05
        real :: FKING, DEPOL, DELTA, DELTAP

        NCOEF = 2
        RAYWIGCOEF(1:NPHASEPOL, :) = 0.0

        FKING = AKING + BKING/WAVELEN**2 + CKING/WAVELEN**4
        DEPOL = 6.d0*(FKING - 1.d0)/(3.d0 + 7.d0*FKING)

        DELTA  = (1.0 - DEPOL)/(1.0 + 0.5*DEPOL)
        DELTAP = (1.0 - 2.0*DEPOL)/(1.0 - DEPOL)

        RAYWIGCOEF(1, 0) = 1.0
        !RAYWIGCOEF(1, 2) = 0.5*DELTA ! original
        RAYWIGCOEF(1, 2) = 0.1 * DELTA ! correction, see Corradini 2004 Eq 5.20

        if (NPHASEPOL == 6) then
            RAYWIGCOEF(2, 2) = 3.0*DELTA
            RAYWIGCOEF(4, 1) = 1.5*DELTAP*DELTA
            RAYWIGCOEF(5, 2) = sqrt(1.5)*DELTA
        end if

    end subroutine RAYLEIGH_PHASE_FUNCTION

    subroutine PHASEFUNC_FROM_LEGENDRE(NANGLES, NCOEF, LEGCOEF, PHASE)
        ! Calculates the phase function at NANGLES equally spaced angles
        ! from the Legendres series expansion in LEGCOEF.  The scalar (I,I)
        ! element of the phase matrix is expanded in Wigner d^l_{0,0} functions,
        ! which are the same as Legendre polynomials in \cos\Theta.
        implicit none
        integer, intent(in) :: NANGLES, NCOEF
        real, intent(in) :: LEGCOEF(0:NCOEF)
        real, intent(out) :: PHASE(NANGLES)
        integer :: J, L
        real(8) :: RD, MU, SUM, PL, PL1, PL2

        RD = dacos(-1.0d0)/180.d0
        do J = 1, NANGLES
            MU = dcos(RD*J*180.d0/NANGLES)
            SUM = 0.0d0
            ! Use upward recurrence to find Legendre polynomials
            PL1 = 1.0d0
            PL = 1.0d0
            do L = 0, NCOEF
                if (L .gt. 0) PL = (2*L - 1)*MU*PL1/L - (L - 1)*PL2/L
                SUM = SUM + LEGCOEF(L)*PL
                PL2 = PL1
                PL1 = PL
            end do
            PHASE(J) = SUM
        end do
    end subroutine PHASEFUNC_FROM_LEGENDRE

    subroutine OUTPUT_PROPERTIES(PROPFILE, NPHASEPOL, NX, NY, NZT, &
                                 DELX, DELY, ZLEVELS, TEMP, &
                                 NPHASE, MAXCOEF, NCOEF, WIGCOEF, &
                                 EXTINCT, SSALB, IPHASE)
        implicit none
        character(len=*), intent(in) :: PROPFILE
        integer, intent(in) :: NPHASEPOL, NX, NY, NZT, NPHASE, MAXCOEF
        integer, intent(in) :: NCOEF(0:NPHASE), IPHASE(NZT, NY, NX)
        real, intent(in) :: DELX, DELY, ZLEVELS(NZT), TEMP(NZT)
        real, intent(in) :: WIGCOEF(NPHASEPOL, 0:MAXCOEF, 0:NPHASE)
        real, intent(in) :: EXTINCT(NZT, NY, NX), SSALB(NZT, NY, NX)
        integer :: NPHASEOUT, I, IX, IY, IZ, J, K, L, NL
        integer*2, allocatable :: PHASECONV(:)

        character(len=256) :: OUTFILE_DATA
        character(len=20)  :: DATA_POSTFIX

        DATA_POSTFIX = '.sHdOm-prp'
        OUTFILE_DATA = trim(adjustl(PROPFILE))//trim(adjustl(DATA_POSTFIX))

        ! Make a phase function index conversion array to output only the ones used.
        allocate (PHASECONV(0:NPHASE))
        PHASECONV(:) = 0
        NPHASEOUT = 0
        do I = 0, NPHASE
            if (any(IPHASE(:, :, :) == I)) then
                NPHASEOUT = NPHASEOUT + 1
                PHASECONV(I) = NPHASEOUT
            end if
        end do

        ! Open the property file and output the header including the
        !   tabulated phase functions
        open (UNIT=3, FILE=trim(adjustl(PROPFILE)), STATUS='REPLACE')
        if (NPHASEPOL == 1) then
            write (3, '(A)') 'Tabulated phase function property file'
        else
            write (3, '(A)') 'Polarized tabulated phase function property file'
        end if
        write (3, '(3(1X,I5))') NX, NY, NZT
        write (3, '(2(1X,F7.4),202(1X,F8.4))') DELX, DELY, (ZLEVELS(I), I=1, NZT)
        write (3, '(1X,I4)') NPHASEOUT
        do I = 0, NPHASE
            if (PHASECONV(I) > 0) then
                NL = NCOEF(I)
                do J = 1, NPHASEPOL
                    if (NPHASEPOL == 1) then
                        write (3, '(1X,I4,5X,200(1X,E15.6))') NL, (WIGCOEF(J, L, I), L=1, min(NL, 200))
                    else
                        write (3, '(1X,I1,1X,I4,5X,201(1X,E15.6))') J, NL, (WIGCOEF(J, L, I), L=0, min(NL, 200))
                    end if
                    do K = 200, NL - 1, 200
                        write (3, '(200(1X,E15.6))') (WIGCOEF(J, K + L, I), L=1, min(NL - K, 200))
                    end do
                end do
            end if
        end do

        ! description for binary data
        !╭────────────────────────────────────────────────╮!
        write (3, '(A)') '! The following provides information for interpreting binary data:'
        write (3, '(A)') '! '//trim(adjustl(DATA_POSTFIX))
        write (3, '(A, 4(I6,A))') '! ', NX, ',', NY, ',', NZT, ',', 4, ''
        !╰────────────────────────────────────────────────╯!

        ! Write one line of optical properties for each grid point (binary)
        !╭────────────────────────────────────────────────╮!
        open (UNIT=88, FILE=trim(adjustl(OUTFILE_DATA)), STATUS='REPLACE',&
            &form='UNFORMATTED', access='STREAM')
        write (88) (((                     TEMP(IZ), IZ=1,NZT), IY=1,NY), IX=1,NX)
        write (88) (((          EXTINCT(IZ, IY, IX), IZ=1,NZT), IY=1,NY), IX=1,NX)
        write (88) (((            SSALB(IZ, IY, IX), IZ=1,NZT), IY=1,NY), IX=1,NX)
        write (88) (((PHASECONV(IPHASE(IZ, IY, IX)), IZ=1,NZT), IY=1,NY), IX=1,NX)
        close (88)
        !╰────────────────────────────────────────────────╯!

        ! original
        !╭────────────────────────────────────────────────╮!
        !do IX = 1, NX
        !    do IY = 1, NY
        !        do IZ = 1, NZT
        !            write (3, '(3(1X,I4),1X,F5.1,1X,F9.4,1X,F8.6,1X,I3)') IX, IY, IZ, &
        !                TEMP(IZ), EXTINCT(IZ, IY, IX), SSALB(IZ, IY, IX), PHASECONV(IPHASE(IZ, IY, IX))
        !        end do
        !    end do
        !end do
        !╰────────────────────────────────────────────────╯!

        close (3)

    end subroutine OUTPUT_PROPERTIES

end module propgen_comb

program propgen

    use propgen_comb

    implicit none
    integer :: NSCATTAB, MAXNEWPHASE, NZO
    real    :: ASYMTOL, FRACPHASETOL, RAYWAVELEN, RAYSFCPRES
    integer, pointer :: SCATNUMS(:)
    logical, pointer :: POLTABS(:)
    real, pointer :: ZOTHER(:), TEMPOTHER(:)
    character(len=256), pointer :: SCATTABFILES(:)
    character(len=256) :: PARFILE, PROPFILE
    character(len=256) :: CKDFILE
    character(len=1) :: POLOUT
    integer, parameter :: NANGLES = 90
    integer :: NTYPES, NENTRIES, MAXNPHASE, NPHASE, NPHASEPOL
    integer :: I, J, N1, N2, NRETAB, MAXCOEF, MAXNCOEF, MAXTYPENUM
    integer :: NX, NY, NZP, NZT, IX, IY, IZ, IZP, NC, NL, IL, IM, IU, L, IPH
    real    :: RAYLCOEF
    real    :: DELX, DELY, EXT, SCAT, SCATTER, ASYM0, F, ERR, MINERR
    logical :: ISEXT
    integer, allocatable :: NRE(:), NCOEF(:)
    integer, allocatable :: PARTYPES(:), TYPECONV(:)
    integer, allocatable :: NCOMP(:, :, :), PTYPE(:, :, :, :), IPARLEV(:)
    real, allocatable    :: WAVELEN1(:), WAVELEN2(:)
    real, allocatable    :: WIGCOEF(:, :), PHASE0(:)
    real, allocatable    :: RETAB(:), EXTINCTTAB(:), SSALBTAB(:), WIGCOEFTAB(:, :, :)
    real, allocatable    :: ZPAR(:), TEMPPAR(:), ZLEVELS(:), TEMP(:), EXTRAYL(:)
    real, allocatable    :: ASYM(:), PHASE(:, :)
    real, allocatable    :: MASSCONT(:, :, :, :), REFF(:, :, :, :)
    real, allocatable    :: EXTINCT(:, :, :), SSALB(:, :, :)
    integer, allocatable :: IPHASE(:, :, :)

    logical :: VERBOSE=.FALSE.

    ! Data structures:

    call USER_INPUT(NSCATTAB, SCATTABFILES, SCATNUMS, POLTABS, PARFILE, &
                    MAXNEWPHASE, ASYMTOL, FRACPHASETOL, &
                    RAYWAVELEN, RAYSFCPRES, NZO, ZOTHER, TEMPOTHER, &
                    POLOUT, PROPFILE)
    MAXNEWPHASE = max(MAXNEWPHASE, 2)
    ASYMTOL = max(ASYMTOL, 0.002)
    FRACPHASETOL = max(FRACPHASETOL, 0.01)

    if (POLOUT == 'P' .or. POLOUT == 'p') then
        NPHASEPOL = 6
    else
        NPHASEPOL = 1
    end if

    NTYPES = NSCATTAB
    allocate (WAVELEN1(NSCATTAB), WAVELEN2(NSCATTAB))
    allocate (NRE(0:NTYPES), PARTYPES(NTYPES))

    NRE(0) = 1
    MAXCOEF = 2
    do I = 1, NSCATTAB
        call READ_SCAT_TABLE_SIZE(SCATTABFILES(I), POLTABS(I), &
                                  WAVELEN1(I), WAVELEN2(I), NRETAB, MAXNCOEF)
        NRE(I) = NRE(I - 1) + NRETAB
        MAXCOEF = max(MAXCOEF, MAXNCOEF)
        if (I > 1 .and. abs((WAVELEN1(I) - WAVELEN1(1))/WAVELEN1(1)) > 0.01 &
            .or. abs((WAVELEN2(I) - WAVELEN2(1))/WAVELEN2(1)) > 0.01) then
            write (*, *) 'Warning [PROPGEN]: Scattering table wavelengths do not agree for table', I
        end if
    end do
    deallocate (WAVELEN1, WAVELEN2)
    NENTRIES = NRE(NSCATTAB) - 1
    MAXNPHASE = NENTRIES + MAXNEWPHASE

    allocate (RETAB(NENTRIES), EXTINCTTAB(NENTRIES), SSALBTAB(NENTRIES))
    allocate (NCOEF(0:MAXNPHASE), WIGCOEFTAB(NPHASEPOL, 0:MAXCOEF, 0:MAXNPHASE))

    WIGCOEFTAB(:, :, :) = 0.0
    ! Put the Rayleigh phase matrix in the zero index
    call RAYLEIGH_PHASE_FUNCTION(RAYWAVELEN, NPHASEPOL, MAXCOEF, &
                                 NCOEF(0), WIGCOEFTAB(:, :, 0))

    ! Read in the optical properties from the scattering table files
    do I = 1, NSCATTAB
        N1 = NRE(I - 1); N2 = NRE(I) - 1
        call READ_SCAT_TABLE(SCATTABFILES(I), POLTABS(I), NPHASEPOL, &
                             N2 - N1 + 1, RETAB(N1:N2), &
                             EXTINCTTAB(N1:N2), SSALBTAB(N1:N2), &
                             MAXCOEF, NCOEF(N1:N2), WIGCOEFTAB(:, 0:, N1:N2))
    end do
    PARTYPES(:) = SCATNUMS(:)
    deallocate (SCATNUMS, SCATTABFILES)

    ! Make the conversion table from particle type to scat table index
    MAXTYPENUM = maxval(PARTYPES(:))
    allocate (TYPECONV(0:MAXTYPENUM))
    call MAKE_TYPE_CONV_TABLE(NTYPES, PARTYPES, MAXTYPENUM, TYPECONV)

    ! Read in the particle properties file
    if (PARFILE /= 'NONE' .and. PARFILE /= 'none') then
        call READ_PARTICLE_FILE_SIZE(PARFILE, NX, NY, NZP)
        allocate (ZPAR(NZP), TEMPPAR(NZP))
        allocate (NCOMP(NZP, NY, NX), PTYPE(NSCATTAB, NZP, NY, NX))
        allocate (MASSCONT(NSCATTAB, NZP, NY, NX), REFF(NSCATTAB, NZP, NY, NX))
        call READ_PARTICLE_FILE(PARFILE, NX, NY, NZP, NSCATTAB, PARTYPES, &
                                DELX, DELY, ZPAR, TEMPPAR, &
                                NCOMP, PTYPE, MASSCONT, REFF, CKDFILE, ISEXT)
    else
        NX = 1; NY = 1; NZP = 0
        DELX = 1.0; DELY = 1.0
        ISEXT = .FALSE.
        CKDFILE = 'NONE'
    end if

    ! Combine the particle levels and extra levels
    NZT = NZP + NZO
    allocate (ZLEVELS(NZT), TEMP(NZT), IPARLEV(NZT), EXTRAYL(NZT))
    call ORGANIZE_LEVELS(NZP, ZPAR, TEMPPAR, NZO, ZOTHER, TEMPOTHER, &
                         NZT, ZLEVELS, TEMP, IPARLEV)

    if (ISEXT) then

        call RAYLEIGH_EXTINCT_FROM_CKD(CKDFILE, NZT, ZLEVELS, EXTRAYL)

    else

        ! Use the parameterization of Bodhaine et al. (1999) eq 30 for tau_R at
        ! sea level and convert to Rayleigh density coefficient:
        !   k = 0.03370*(p_sfc/1013.25)*tau_R  for k in K/(mb km) (from assumptions in RAYLEIGH_EXTINCT)
        RAYLCOEF = 0.03370*(RAYSFCPRES/1013.25)*0.0021520 &
                   *(1.0455996 - 341.29061/RAYWAVELEN**2 - 0.90230850*RAYWAVELEN**2) &
                   /(1 + 0.0027059889/RAYWAVELEN**2 - 85.968563*RAYWAVELEN**2)

        if (VERBOSE) write (*, *) 'Rayleigh total atmosphere optical depth: ', RAYLCOEF/0.03370

        ! Compute the molecular Rayleigh scattering extinction profile
        call RAYLEIGH_EXTINCT(NZT, ZLEVELS, TEMP, RAYSFCPRES, RAYLCOEF, EXTRAYL)

    end if

    allocate (ASYM(0:MAXNPHASE), PHASE(NANGLES, 0:MAXNPHASE))
    NPHASE = NENTRIES
    ! Make the asymmetry param and angle space low resolution phase functions
    do I = 0, NPHASE
        ASYM(I) = WIGCOEFTAB(1, 1, I)/3
        call PHASEFUNC_FROM_LEGENDRE(NANGLES, NCOEF(I), WIGCOEFTAB(1, :, I), PHASE(:, I))
    end do

    ! Allocate the output optical property arrays
    allocate (EXTINCT(NZT, NY, NX), SSALB(NZT, NY, NX), IPHASE(NZT, NY, NX))

    ! Do the optical property mixing from the components specified at each
    !   grid point in the particle file.
    allocate (WIGCOEF(NPHASEPOL, 0:MAXCOEF), PHASE0(NANGLES))
    do IX = 1, NX
        do IY = 1, NY
            do IZ = 1, NZT
                IZP = IPARLEV(IZ)
                if (IZP == 0) then
                    EXTINCT(IZ, IY, IX) = EXTRAYL(IZ) ! Only Rayleigh scattering if not
                    SSALB(IZ, IY, IX) = 1.0 !   at a particle level
                    IPHASE(IZ, IY, IX) = 0
                else
                    NC = NCOMP(IZP, IY, IX)
                    ! Initialize optical property sums with Rayleigh scattering
                    EXTINCT(IZ, IY, IX) = EXTRAYL(IZ)
                    SCATTER = EXTRAYL(IZ)
                    ASYM0 = 0.0
                    PHASE0(:) = EXTRAYL(IZ)*PHASE(:, 0)
                    NL = 2
                    WIGCOEF(:, :) = EXTRAYL(IZ)*WIGCOEFTAB(:, :, 0)
                    do J = 1, NC
                        L = TYPECONV(PTYPE(J, IZP, IY, IX))
                        ! Do binary search to find effective radius entry in table
                        IL = NRE(L - 1)
                        IU = NRE(L) - 1
                        do while (IU - IL > 1)
                            IM = (IL + IU)/2
                            if (REFF(J, IZP, IY, IX) >= RETAB(IM)) then
                                IL = IM
                            else
                                IU = IM
                            end if
                        end do
                        if (IU > IL) then
                            F = (REFF(J, IZP, IY, IX) - RETAB(IL))/(RETAB(IU) - RETAB(IL))
                        else if (IU == IL .and. &
                                 abs((REFF(J, IZP, IY, IX) - RETAB(IL))/RETAB(IL)) < 0.001) then
                            F = 0.0
                        else
                            F = -1.0
                        end if
                        if (F < 0.0 .or. F > 1.0) then
                            write (*, *) 'Warning [PROPGEN]: Effective radius outside of table (IX,IY,IZ,type,Reff):'
                            write (*, '(4(1X,I3),1X,F6.2)') IX, IY, IZP, PTYPE(J, IZP, IY, IX), REFF(J, IZP, IY, IX)
                            F = max(min(F, 1.0), 0.0)
                        end if

                        if (ISEXT) then
                            EXT = MASSCONT(J, IZP, IY, IX)
                        else
                            ! Interpolate optical properties in Reff for this component
                            !   and sum over the components.
                            EXT = MASSCONT(J, IZP, IY, IX)* &
                                  exp((1 - F)*log(EXTINCTTAB(IL)) + F*log(EXTINCTTAB(IU)))
                        end if

                        SCAT = EXT*((1 - F)*SSALBTAB(IL) + F*SSALBTAB(IU))
                        ASYM0 = ASYM0 + SCAT*((1 - F)*ASYM(IL) + F*ASYM(IU))
                        PHASE0(:) = PHASE0(:) + SCAT*((1 - F)*PHASE(:, IL) + F*PHASE(:, IU))
                        NL = max(NL, NCOEF(IL), NCOEF(IU))
                        WIGCOEF(:, 0:NL) = WIGCOEF(:, 0:NL) &
                                        + SCAT*((1 - F)*WIGCOEFTAB(:, 0:NL, IL) + F*WIGCOEFTAB(:, 0:NL, IU))
                        SCATTER = SCATTER + SCAT
                        EXTINCT(IZ, IY, IX) = EXTINCT(IZ, IY, IX) + EXT

                    end do ! end of component loop

                    ! Normalize to get the combined optical properties
                    if (EXTINCT(IZ, IY, IX) > 0.0) then
                        SSALB(IZ, IY, IX) = SCATTER/EXTINCT(IZ, IY, IX)
                    else
                        SSALB(IZ, IY, IX) = 1.0
                    end if
                    if (SCATTER > 0.0) then
                        ASYM0 = ASYM0/SCATTER
                        PHASE0(:) = PHASE0(:)/SCATTER
                        WIGCOEF(:, 0:NL) = WIGCOEF(:, 0:NL)/SCATTER

                        ! Find the closest phase function in terms of the asymmetry parameter
                        !   and max fractional function value error normalized by user tolerances
                        IPH = 0 ! initialize with the Rayleigh phase function
                        MINERR = abs(ASYM0 - ASYM(IPH))/ASYMTOL &
                             + maxval(abs((PHASE0(:) - PHASE(:, IPH))/(0.001 + PHASE(:, IPH))))/FRACPHASETOL
                        do I = 1, NPHASE ! search all phase functions
                            ERR = abs(ASYM0 - ASYM(I))/ASYMTOL &
                               + maxval(abs((PHASE0(:) - PHASE(:, I))/max(0.001, PHASE(:, I))))/FRACPHASETOL
                            if (ERR < MINERR) then
                                MINERR = ERR
                                IPH = I
                            end if
                        end do
                        ! If the closest phase function is within the tolerances then use it,
                        !   otherwise add the current grid point phase function to the list.
                        if (abs(ASYM0 - ASYM(IPH)) < ASYMTOL .and. FRACPHASETOL > &
                            maxval(abs((PHASE0(:) - PHASE(:, IPH))/max(0.001, PHASE(:, IPH))))) then
                            IPHASE(IZ, IY, IX) = IPH
                        else
                            NPHASE = NPHASE + 1
                            if (NPHASE > MAXNPHASE) then
                                write (*, *) 'Error [PROPGEN]: Maximum number of new phase functions exceeded', NPHASE - NENTRIES
                                stop
                            end if
                            NCOEF(NPHASE) = NL
                            WIGCOEFTAB(:, :, NPHASE) = WIGCOEF(:, :)
                            ASYM(NPHASE) = WIGCOEFTAB(1, 1, NPHASE)/3
                            call PHASEFUNC_FROM_LEGENDRE(NANGLES, NL, WIGCOEF(1, 0:), &
                                                         PHASE(:, NPHASE))
                            IPHASE(IZ, IY, IX) = NPHASE
                        end if

                    else ! if no scattering then just use the Rayleigh phase function
                        IPHASE(IZ, IY, IX) = 0
                    end if
                end if
            end do ! end of Z level loop
        end do
    end do

    if (VERBOSE) write (*, *) 'Number of new phase functions added: ', NPHASE - NENTRIES

    call OUTPUT_PROPERTIES(PROPFILE, NPHASEPOL, NX, NY, NZT, &
                           DELX, DELY, ZLEVELS, TEMP, &
                           NPHASE, MAXCOEF, NCOEF, WIGCOEFTAB, &
                           EXTINCT, SSALB, IPHASE)

    deallocate (EXTINCT, SSALB, IPHASE, WIGCOEF, ASYM, PHASE0, PHASE)
    if (allocated(MASSCONT)) then
        deallocate (ZPAR, TEMPPAR, NCOMP, PTYPE, MASSCONT, REFF)
    end if
    deallocate (ZLEVELS, TEMP, IPARLEV)
    deallocate (NCOEF, WIGCOEFTAB, RETAB, EXTINCTTAB, SSALBTAB)
    deallocate (TYPECONV, PARTYPES, NRE)
    deallocate (POLTABS)

end
