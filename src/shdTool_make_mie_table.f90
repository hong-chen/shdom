program MAKE_MIE_TABLE
!
! Does Mie computations to create a scattering table as a function of
! effective radius for gamma, modified gamma, or lognormal size distributions
! of spherical particles.  The particles may be water or ice (in which case
! the program provides the index of refraction depending on wavelength) or
! "aerosols" (in which case the index of refraction is user specified).
! For water or ice particles the scattering properties may be averaged
! over the desired spectral range with Planck function weighting.
! The phase functions in the output scattering table are represented
! with Legendre series.  For polarized output, the six unique elements of
! the phase matrix are represented with Wigner d-function expansions
! (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
!
!  compile: pgf90 -fast -o make_mie_table  make_mie_table.f90
!                 indexwatice.f  miewig.f
!
!    Frank Evans    University of Colorado       May 2003
!  Modified for polarization                     March 2012
!  Modified for modified gamma distribution      September 2012
!  Modified for Wigner d-function coefficients   July 2012
    implicit none
    integer :: NRETAB
    logical :: POLTAB, LOGRE
    real    :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
    real    :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
    complex :: RINDEX
    character(LEN=1) :: PARTYPE, AVGFLAG, DISTFLAG
    character(LEN=80) :: MIETABFILE
    integer :: NSIZE, MAXRANK, I, J, L, NL
    real    :: PI, WAVELENCEN, XMAX, REFF, SCATTER
    integer, allocatable :: NRANK1(:), NRANK(:)
    real, allocatable :: RADII(:), ND(:)
    real, allocatable :: QEXT(:), QSCA(:)
    real, allocatable :: EXTINCT1(:), SCATTER1(:), WIGCOEF1(:, :, :)
    real, allocatable :: EXTINCT(:), SSALB(:), WIGCOEF(:, :, :)

    ! Get the input parameters
    call USER_INPUT(POLTAB, WAVELEN1, WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                    AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, GAMMA, &
                    NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS, MIETABFILE)
    if (PARTYPE == 'W') then
        PARDENS = 1.0
    else if (PARTYPE == 'I') then
        PARDENS = 0.916
    end if

    ! Calculate the maximum size parameter and the max number of Legendre terms
    call GET_CENTER_WAVELEN(WAVELEN1, WAVELEN2, WAVELENCEN)
    PI = acos(-1.0)
    if (AVGFLAG == 'A') then
        XMAX = 2*PI*MAXRADIUS/WAVELEN1
    else
        XMAX = 2*PI*MAXRADIUS/WAVELENCEN
    end if
    MAXRANK = nint(2*(XMAX + 4.0*XMAX**0.3334 + 2))

    ! Get the average index of refraction for water or ice
    if (PARTYPE /= 'A') then
        call GET_REFRACT_INDEX(PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
    end if

    ! Figure the number of radii there will be
    call GET_NSIZE(SRETAB, MAXRADIUS, WAVELENCEN, NSIZE)

    ! Allocate all the arrays here
    allocate (RADII(NSIZE), ND(NSIZE), NRANK1(NSIZE))
    allocate (EXTINCT1(NSIZE), SCATTER1(NSIZE), WIGCOEF1(6, 0:MAXRANK, NSIZE))
    allocate (EXTINCT(NRETAB), SSALB(NRETAB))
    allocate (NRANK(NRETAB), WIGCOEF(6, 0:MAXRANK, NRETAB))

    ! Make up the discrete particle radii to use
    call GET_SIZES(SRETAB, MAXRADIUS, WAVELENCEN, NSIZE, RADII)

    ! Do the Mie computations for each radius, which may involve several
    !   Mie calculation over the wavelength integration
    call COMPUTE_MIE_ALL_SIZES(AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, PARTYPE, &
                               WAVELENCEN, RINDEX, NSIZE, RADII, MAXRANK, &
                               EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)

    ! Loop over the number of output tabulated effective radii
    do I = 1, NRETAB
        ! Set tabulated effective radius
        if (NRETAB <= 1) then
            REFF = SRETAB
        else if (LOGRE) then
            REFF = exp((log(ERETAB) - log(SRETAB))*FLOAT(I - 1)/(NRETAB - 1) + log(SRETAB))
        else
            REFF = (ERETAB - SRETAB)*FLOAT(I - 1)/(NRETAB - 1) + SRETAB
        end if

        ! Calculate the discrete size number concentrations (ND), which vary
        ! according to a truncated gamma, modified gamma, or lognormal
        ! distribution that gives the desired effective radius (REFF) and LWC (1 g/m^3).
        call MAKE_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA, &
                            ND)

        ! Sum the scattering properties over the discrete size distribution
        EXTINCT(I) = 0.0
        SCATTER = 0.0
        WIGCOEF(:, :, I) = 0.0
        NL = 1
        do J = 1, NSIZE
            EXTINCT(I) = EXTINCT(I) + ND(J)*EXTINCT1(J)
            SCATTER = SCATTER + ND(J)*SCATTER1(J)
            NL = max(NL, NRANK1(J))
            WIGCOEF(:, 0:NL, I) = WIGCOEF(:, 0:NL, I) + ND(J)*WIGCOEF1(:, 0:NL, J)
        end do
        do L = 0, NL
            WIGCOEF(:, L, I) = WIGCOEF(:, L, I)/SCATTER
            if (WIGCOEF(1, L, I) .gt. 0.5e-5) NRANK(I) = L
        end do
        if (abs(WIGCOEF(1, 0, I) - 1.0) > 0.0001) then
            print *, 'Phase function not normalized for Reff=', REFF, WIGCOEF(1, 0, I)
            stop
        end if
        if (EXTINCT(I) > 0.0) then
            SSALB(I) = SCATTER/EXTINCT(I)
        end if
        EXTINCT(I) = 0.001*EXTINCT(I)

    end do ! end of effective radius loop

    call WRITE_MIE_TABLE(MIETABFILE, POLTAB, WAVELEN1, WAVELEN2, DELTAWAVE, &
                         PARTYPE, PARDENS, RINDEX, &
                         DISTFLAG, ALPHA, GAMMA, NRETAB, SRETAB, ERETAB, LOGRE, &
                         EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)

    deallocate (EXTINCT, SSALB, NRANK, WIGCOEF)
    deallocate (ND, EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)
end

subroutine USER_INPUT(POLTAB, WAVELEN1, WAVELEN2, PARTYPE, RINDEX, PARDENS, &
                      AVGFLAG, DELTAWAVE, DISTFLAG, ALPHA, GAMMA, &
                      NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS, MIETABFILE)
    ! Reads the input parameters from the standard input
    implicit none
    integer, intent(OUT) :: NRETAB
    logical, intent(OUT) :: POLTAB, LOGRE
    real, intent(OUT) :: WAVELEN1, WAVELEN2, DELTAWAVE, PARDENS
    real, intent(OUT) :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
    complex, intent(OUT) :: RINDEX
    character(LEN=1), intent(OUT) :: PARTYPE, AVGFLAG, DISTFLAG
    character(LEN=*), intent(OUT) :: MIETABFILE

    write (*, *) 'Making Mie scattering tables for spherical particles'

    write (*, *) 'Make polarized Mie table (T or F)'
    read (*, *) POLTAB
    write (*, '(L)') POLTAB

    write (*, *) 'Wavelength range (micron)'
    read (*, *) WAVELEN1, WAVELEN2
    write (*, '(2(1X,F9.3))') WAVELEN1, WAVELEN2
    if (WAVELEN1 > WAVELEN2) stop 'USER_INPUT: wavelength1 must be <= wavelength2'

    write (*, *) 'Water, Ice, or Aerosol spherical particles (W,I,A)'
    read (*, '(A1)') PARTYPE
    write (*, '(1X,A1)') PARTYPE

    if (PARTYPE == 'W' .or. PARTYPE == 'I') then
        write (*, *) 'Average Mie properties over wavelength or use Planck weighted'
        write (*, *) '  index of refraction at center wavelength (A or C)'
        read (*, '(A1)') AVGFLAG
        write (*, '(1X,A1)') AVGFLAG
        if (AVGFLAG == 'A') then
            write (*, *) 'Wavelength interval for averaging (micron)'
            read (*, *) DELTAWAVE
            write (*, '(1X,F6.3)') DELTAWAVE
        else
            AVGFLAG = 'C'
            DELTAWAVE = -1.0
        end if
    else
        PARTYPE = 'A'
        write (*, *) 'Aerosol complex index of refraction (negative imaginary part)'
        read (*, *) RINDEX
        write (*, *) RINDEX
        write (*, *) 'Aerosol particle bulk density (g/cm^3)'
        read (*, *) PARDENS
        write (*, '(1X,F5.3)') PARDENS
        AVGFLAG = 'C'
    end if

    write (*, *) 'Particle size distribution type: G = Gamma, M = modified gamma, L = Lognormal'
    read (*, *) DISTFLAG
    write (*, *) DISTFLAG
    if (DISTFLAG == 'L') then
        write (*, *) 'Log normal size distribution log standard deviation'
        read (*, *) ALPHA
        write (*, '(1X,F6.3)') ALPHA
    else if (DISTFLAG == 'G') then
        write (*, *) 'Gamma size distribution shape parameter (alpha)'
        read (*, *) ALPHA
        write (*, '(1X,F6.3)') ALPHA
    else if (DISTFLAG == 'M') then
        write (*, *) 'Modified gamma size distribution shape parameters (alpha & gamma)'
        read (*, *) ALPHA, GAMMA
        write (*, '(2(1X,F6.3))') ALPHA, GAMMA
    else
        write (*, *) 'Unrecognized size distribution type'
        stop
    end if

    write (*, *) 'Number, starting, and ending tabulated effective radius (micron)'
    read (*, *) NRETAB, SRETAB, ERETAB
    write (*, '(1X,I3,2(1X,F7.2))') NRETAB, SRETAB, ERETAB

    write (*, *) 'Log-spaced effective radius (T or F) (F for evenly spaced)'
    read (*, *) LOGRE
    write (*, '(L)') LOGRE

    write (*, *) 'Maximum particle radius in size distribution (micron)'
    read (*, *) MAXRADIUS
    write (*, '(1X,F7.2)') MAXRADIUS

    write (*, *) 'Output Mie scattering table name'
    read (*, *) MIETABFILE
    write (*, '(1X,A70)') MIETABFILE
end subroutine USER_INPUT

subroutine GET_NSIZE(SRETAB, MAXRADIUS, WAVELEN, NSIZE)
    ! Calculates the number of radii for which the Mie computation will be run.
    ! The formula and spacing in size parameter can be changed to trade
    ! off size distribution integration accuracy vs. computer time.
    implicit none
    real, intent(IN)  :: SRETAB, MAXRADIUS, WAVELEN
    integer, intent(OUT) :: NSIZE
    real    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD, DELRADMAX

    TWOPI = 2.0*acos(-1.0)
    RADMIN = 0.02*SRETAB
    RAD = RADMIN
    DELRADMAX = 0.002*MAXRADIUS
    NSIZE = 1
    do while (RAD < MAXRADIUS)
        X = TWOPI*RAD/WAVELEN
        DELX = max(0.01, 0.03*X**0.5) ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
        DELRAD = min(DELX*WAVELEN/TWOPI, DELRADMAX)
        RAD = RAD + DELRAD
        NSIZE = NSIZE + 1
    end do
end subroutine GET_NSIZE

subroutine GET_SIZES(SRETAB, MAXRADIUS, WAVELEN, NSIZE, RADII)
    ! Calculates the radii for which the Mie computation will be run and
    ! from which all the size distributions will be computed.
    ! The formula and spacing in size parameter can be changed to trade
    ! off size distribution integration accuracy vs. computer time.
    implicit none
    integer, intent(IN) :: NSIZE
    real, intent(IN) :: SRETAB, MAXRADIUS, WAVELEN
    real, intent(OUT) :: RADII(NSIZE)
    integer :: N
    real    :: TWOPI, RADMIN, RAD, X, DELX, DELRAD, DELRADMAX

    TWOPI = 2.0*acos(-1.0)
    RAD = 0.02*SRETAB
    DELRADMAX = 0.002*MAXRADIUS
    RADII(1) = RAD
    do N = 2, NSIZE
        X = TWOPI*RAD/WAVELEN
        DELX = max(0.01, 0.03*X**0.5) ! coarser spacing at large size parameters
!    DELX = 0.1                     ! One alternative method
        DELRAD = min(DELX*WAVELEN/TWOPI, DELRADMAX)
        RAD = RAD + DELRAD
        RADII(N) = RAD
    end do
end subroutine GET_SIZES

subroutine GET_CENTER_WAVELEN(WAVELEN1, WAVELEN2, WAVELENCEN)
!  Returns the Planck weighted center wavelength averaged over the
! wavelength interval (WAVELEN1 < WAVELEN2 [microns]).  A solar
! blackbody temperature of 5800 K is used for the Planck weighting
! if the average wavelength is less than 3 microns, no Planck weighting
! is done for an average wavelength between 3 and 5 microns, and a
! blackbody temperature of 270 K is done for an average wavelength
! greater than 5 microns.
    implicit none
    real, intent(IN)  :: WAVELEN1, WAVELEN2
    real, intent(OUT) :: WAVELENCEN
    real :: WAVECEN, DELWAVE, WAVE, BBTEMP, PLANCK, SUMP, SUMW

    if (WAVELEN1 == WAVELEN2) then
        WAVELENCEN = WAVELEN1
    else
        WAVECEN = 0.5*(WAVELEN1 + WAVELEN2)
        if (WAVECEN < 3.0) then
            BBTEMP = 5800.0
        else if (WAVECEN > 5.0) then
            BBTEMP = 270.0
        else
            BBTEMP = -1.0
            PLANCK = 1.0
        end if
        DELWAVE = min(WAVECEN/100., 0.1*abs(WAVELEN2 - WAVELEN1))
        SUMP = 0.0
        SUMW = 0.0
        WAVE = WAVELEN1
        do while (WAVE <= WAVELEN2)
            if (BBTEMP > 0) PLANCK = (1.19e8/WAVE**5)/(exp(1.439e4/(WAVE*BBTEMP)) - 1)
            SUMP = SUMP + PLANCK
            SUMW = SUMW + PLANCK*WAVE
            WAVE = WAVE + DELWAVE
        end do
        WAVELENCEN = 0.001*nint(1000*SUMW/SUMP)
    end if
end subroutine GET_CENTER_WAVELEN

subroutine GET_REFRACT_INDEX(PARTYPE, WAVELEN1, WAVELEN2, RINDEX)
    ! Returns the index of refraction for water or ice averaged over
    ! the wavelength interval (WAVELEN1 < WAVELEN2 [microns]).   The
    ! averaging is done at 0.05 micron intervals and is weighted by
    ! a Planck function at 5800 K temperature for central wavelengths
    ! less than 3 microns, a flat function between 3 and 5 microns, and
    ! 270 K Planck function beyond 5 microns.
    ! The index of refraction is at -30 C for ice and +10 C for water
    ! (the temperature dependence is important in the microwave).
    implicit none
    character(LEN=1), intent(IN) :: PARTYPE
    real, intent(IN) :: WAVELEN1, WAVELEN2
    complex, intent(OUT) :: RINDEX
    real :: WAVECEN, WAVECUT, DELWAVE, WAVE, BBTEMP, PLANCK
    real :: MRE, MIM, SUMP, SUMMR, SUMMI, A

    WAVECEN = 0.5*(WAVELEN1 + WAVELEN2)
    if (WAVECEN < 3.0) then
        BBTEMP = 5800.0
    else if (WAVECEN > 5.0) then
        BBTEMP = 270.0
    else
        BBTEMP = -1.0
        PLANCK = 1.0
    end if
    DELWAVE = min(WAVECEN/100., 0.1*abs(WAVELEN2 - WAVELEN1))
    DELWAVE = max(DELWAVE, WAVECEN*1.0e-5)
    SUMP = 0.0
    SUMMR = 0.0
    SUMMI = 0.0
    WAVE = WAVELEN1
    do while (WAVE <= WAVELEN2)
        if (BBTEMP > 0) PLANCK = (1.19e8/WAVE**5)/(exp(1.439e4/(WAVE*BBTEMP)) - 1)
        SUMP = SUMP + PLANCK
        if (PARTYPE == 'I') then
            call REFICE(0, WAVE, 243.0, MRE, MIM, A, A)
        else
            call REFWAT(0, WAVE, 283.0, MRE, MIM, A, A)
        end if
        SUMMR = SUMMR + PLANCK*MRE
        SUMMI = SUMMI + PLANCK*MIM
        WAVE = WAVE + DELWAVE
    end do
    MRE = SUMMR/SUMP
    MIM = SUMMI/SUMP
    RINDEX = cmplx(MRE, -MIM)
end subroutine GET_REFRACT_INDEX

subroutine COMPUTE_MIE_ALL_SIZES(AVGFLAG, WAVELEN1, WAVELEN2, DELTAWAVE, &
                                 PARTYPE, WAVELENCEN, RINDEX, NSIZE, RADII, &
                                 MAXRANK, EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)
    ! Does a Mie computation for each particle radius in RADII and returns the
    ! optical properties in arrays EXTINCT1, SCATTER1, NRANK1, and WIGCOEF1.
    ! For AVGFLAG='C' the computation is done at a single wavelength (WAVELENCEN),
    ! using the input index of refraction (RINDEX).  For AVGFLAG='A' an
    ! integration of the Mie properties over wavelength is performed for
    ! each radius.  For each wavelength, with spacing DELTAWAVE, the water
    ! or ice (depending on PARTYPE) index of refraction is obtained and
    ! used in the Mie computation for that wavelength, and the Mie optical
    ! properties are averaged with Planck function weighting (blackbody
    ! temperature depends on wavelength).  The Wigner d coefficients are
    ! returned with the product of the phase function times the scattering
    ! coefficient.
    implicit none
    integer, intent(IN) :: NSIZE, MAXRANK
    real, intent(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE, WAVELENCEN
    real, intent(IN) :: RADII(NSIZE)
    complex, intent(IN) :: RINDEX
    character(LEN=1), intent(IN) :: AVGFLAG, PARTYPE
    integer, intent(OUT) :: NRANK1(NSIZE)
    real, intent(OUT) :: EXTINCT1(NSIZE), SCATTER1(NSIZE)
    real, intent(OUT) :: WIGCOEF1(6, 0:MAXRANK, NSIZE)
    integer :: I, NL
    real    :: WAVECEN, WAVE, BBTEMP, PLANCK, SUMP, A
    real    :: MRE, MIM, EXT, SCAT, COEF(6, 0:MAXRANK)
    complex :: REFIND

    if (AVGFLAG == 'C') then
        ! For using one central wavelength: just call Mie routine for each radius
        do I = 1, NSIZE
            call MIE_ONE(WAVELENCEN, RINDEX, RADII(I), MAXRANK, &
                         EXTINCT1(I), SCATTER1(I), NRANK1(I), WIGCOEF1(1, 0, I))
        end do

    else
        ! For averaging over wavelength range:
        WAVECEN = 0.5*(WAVELEN1 + WAVELEN2)
        if (WAVECEN < 3.0) then
            BBTEMP = 5800.0
        else if (WAVECEN > 5.0) then
            BBTEMP = 270.0
        else
            BBTEMP = -1.0
            PLANCK = 1.0
        end if
        EXTINCT1(:) = 0.0
        SCATTER1(:) = 0.0
        NRANK1(:) = 1
        WIGCOEF1(:, :, :) = 0.0
        SUMP = 0.0
        WAVE = WAVELEN1
        do while (WAVE <= WAVELEN2) ! Loop over the wavelengths
            if (BBTEMP > 0) PLANCK = (1.19e8/WAVE**5)/(exp(1.439e4/(WAVE*BBTEMP)) - 1)
            SUMP = SUMP + PLANCK
            if (PARTYPE == 'I') then ! Get the index of refraction of water or ice
                call REFICE(0, WAVE, 243.0, MRE, MIM, A, A)
            else
                call REFWAT(0, WAVE, 283.0, MRE, MIM, A, A)
            end if
            REFIND = cmplx(MRE, -MIM)
            do I = 1, NSIZE
                call MIE_ONE(WAVE, REFIND, RADII(I), MAXRANK, EXT, SCAT, NL, COEF)
                EXTINCT1(I) = EXTINCT1(I) + PLANCK*EXT
                SCATTER1(I) = SCATTER1(I) + PLANCK*SCAT
                NRANK1(I) = max(NRANK1(I), NL)
                WIGCOEF1(:, 0:NL, I) = WIGCOEF1(:, 0:NL, I) + PLANCK*COEF(:, 0:NL)
            end do
            WAVE = WAVE + DELTAWAVE
        end do
        EXTINCT1(:) = EXTINCT1(:)/SUMP
        SCATTER1(:) = SCATTER1(:)/SUMP
        WIGCOEF1(:, :, :) = WIGCOEF1(:, :, :)/SUMP
    end if
end subroutine COMPUTE_MIE_ALL_SIZES

subroutine MAKE_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA, &
                          ND)
    ! Calculates the number concentrations (ND in cm^-3) for the NSIZE
    ! discrete particle radii (micron) of a gamma or lognormal size distribution
    ! with an effective radius of REFF (micron), gamma shape parameter or
    ! lognormal standard deviation of ALPHA, and mass content of 1 g/m^3.
    ! Also handles modified gamma distributions specified by ALPHA and GAMMA.
    implicit none
    integer, intent(IN)  :: NSIZE
    real, intent(IN)  :: RADII(NSIZE), REFF, ALPHA, GAMMA, PARDENS
    real, intent(OUT) :: ND(NSIZE)
    character(LEN=1), intent(IN) :: DISTFLAG
    real, parameter :: TOL = 0.001 ! fractional tolerance in achieving Reff
    integer :: I
    real    :: TRUERE, F, REHI, RELO, REMID

    ! See if the true effective radius is already close enough
    call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, REFF, &
                      NSIZE, RADII, ND, TRUERE)
    if (abs(TRUERE - REFF) < TOL*REFF) return
    F = REFF/TRUERE

    if (TRUERE < REFF) then
        ! Find Reff that gives true Reff above desired value
        RELO = REFF
        REHI = F*REFF
        I = 0
        TRUERE = REFF/F
        do while (TRUERE <= REFF .and. I < 8)
            REHI = F*REHI
            I = I + 1
            call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, REHI, &
                              NSIZE, RADII, ND, TRUERE)
        end do
        if (TRUERE <= REFF) then
            print *, 'MAKE_SIZE_DIST: effective radius cannot be achieved', REFF, TRUERE
            stop
        end if
    else
        ! Find Reff that gives true Reff below desired value
        REHI = REFF
        RELO = F*REFF
        I = 0
        TRUERE = REFF/F
        do while (TRUERE >= REFF .and. I < 8)
            RELO = F*RELO
            I = I + 1
            call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, RELO, &
                              NSIZE, RADII, ND, TRUERE)
        end do
        if (TRUERE >= REFF) then
            print *, 'MAKE_SIZE_DIST: effective radius cannot be achieved', REFF, TRUERE
            stop
        end if
    end if
    ! Do bisection to get correct effective radius
    do while (abs(TRUERE - REFF) > TOL*REFF)
        REMID = 0.5*(RELO + REHI)
        call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, REMID, &
                          NSIZE, RADII, ND, TRUERE)
        if (TRUERE < REFF) then
            RELO = REMID
        else
            REHI = REMID
        end if
    end do
end subroutine MAKE_SIZE_DIST

subroutine DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, RE, NSIZE, RADII, &
                        ND, TRUERE)
    ! For the input effective radius (RE) [um], returns the number concentrations
    ! ND [cm^-3] and the calculated effective radius TRUERE [um] for a
    ! gamma, modified gamma, or lognormal size distribution with mass content
    ! of 1 g/m^3.
    implicit none
    integer, intent(IN) :: NSIZE
    real, intent(IN) :: PARDENS, ALPHA, GAMMA, RE, RADII(NSIZE)
    character(LEN=1), intent(IN) :: DISTFLAG
    real, intent(OUT) :: ND(NSIZE), TRUERE
    integer :: J
    real    :: DENS, PI, A, B, LWC, R, DELR, SUM2, SUM3, GAMMLN

    PI = acos(-1.0)
    if (DISTFLAG == 'G') then
        B = (ALPHA + 3)/RE
        A = 1.e6/((4*PI/3.)*PARDENS*B**(-ALPHA - 4)*exp(GAMMLN(ALPHA + 4.)))
    else if (DISTFLAG == 'M') then
        B = (exp(GAMMLN((ALPHA + 4.)/GAMMA) - GAMMLN((ALPHA + 3.)/GAMMA))/RE)**GAMMA
        A = 1.e6*GAMMA*B**((ALPHA + 4)/GAMMA) &
            /((4*PI/3.)*PARDENS*exp(GAMMLN((ALPHA + 4.)/GAMMA)))
    else if (DISTFLAG == 'L') then
        B = RE*exp(-2.5*ALPHA**2)
        A = 1.e6/((4*PI/3.)*PARDENS*sqrt(2*PI)*ALPHA*B**3*exp(4.5*ALPHA**2))
    end if
    LWC = 0.0
    SUM2 = 0.0
    SUM3 = 0.0
    do J = 1, NSIZE
        R = RADII(J)
        DELR = sqrt(RADII(J)*RADII(min(NSIZE, J + 1))) &
               - sqrt(RADII(J)*RADII(max(1, J - 1)))
        if (DISTFLAG == 'G') then
            ND(J) = A*R**ALPHA*exp(-B*R)*DELR
        else if (DISTFLAG == 'M') then
            ND(J) = A*R**ALPHA*exp(-B*R**GAMMA)*DELR
        else if (DISTFLAG == 'L') then
            ND(J) = (A/R)*exp(-0.5*(log(R/B))**2/ALPHA**2)*DELR
        end if
        LWC = LWC + 1.0e-6*PARDENS*ND(J)*(4*PI/3)*R**3
        SUM2 = SUM2 + ND(J)*R**2
        SUM3 = SUM3 + ND(J)*R**3
    end do
    ND(:) = (1.0/LWC)*ND(:)
    TRUERE = SUM3/SUM2
end subroutine DO_SIZE_DIST

subroutine WRITE_MIE_TABLE(MIETABFILE, POLTAB, WAVELEN1, WAVELEN2, DELTAWAVE, &
                           PARTYPE, PARDENS, RINDEX, &
                           DISTFLAG, ALPHA, GAMMA, NRETAB, SRETAB, ERETAB, LOGRE, &
                           EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)
    ! Writes the table of Mie scattering properties as a function of
    ! effective radius.  There are two types of tables output:
    ! unpolarized with Legendre (Wigner d_l_00) series of a single phase
    ! matrix element (P11) or polarized with six Wigner d-function series
    ! of the Mie phase matrix elements.
    implicit none
    integer, intent(IN) :: NRETAB, MAXRANK, NRANK(NRETAB)
    logical, intent(IN) :: POLTAB, LOGRE
    complex, intent(IN) :: RINDEX
    real, intent(IN) :: WAVELEN1, WAVELEN2, DELTAWAVE
    real, intent(IN) :: PARDENS, SRETAB, ERETAB, ALPHA, GAMMA
    real, intent(IN) :: EXTINCT(NRETAB), SSALB(NRETAB)
    real, intent(IN) :: WIGCOEF(6, 0:MAXRANK, NRETAB)
    character(LEN=1), intent(IN) :: PARTYPE, DISTFLAG
    character(LEN=*), intent(IN) :: MIETABFILE
    integer :: I, J, JP, K, L, NL, NPHASEPOL
    real REFF

    open (UNIT=3, FILE=MIETABFILE, STATUS='REPLACE')
    if (POLTAB) then
        NPHASEPOL = 6
        write (3, '(A)') '! Polarized Mie scattering table vs. effective radius (LWC=1 g/m^3)'
    else
        NPHASEPOL = 1
        write (3, '(A)') '! Mie scattering table vs. effective radius (LWC=1 g/m^3)'
    end if
    if (DELTAWAVE < 0.0) then
        write (3, '(2(1X,F8.3),A)') WAVELEN1, WAVELEN2, '  wavelength range (micron)'
    else
        write (3, '(3(1X,F8.3),A)') WAVELEN1, WAVELEN2, DELTAWAVE, '  wavelength range and averaging step (micron)'
    end if
    write (3, '(1X,F5.3,2X,A1,A)') PARDENS, PARTYPE, '   particle density (g/cm^3) and type (Water, Ice, Aerosol)'
    write (3, '(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
    if (DISTFLAG == 'L') then
        write (3, '(F7.5,A)') ALPHA, '  lognormal log standard deviation'
    else if (DISTFLAG == 'G') then
        write (3, '(F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
    else if (DISTFLAG == 'M') then
        write (3, '(2(1X,F7.5),A)') ALPHA, GAMMA, ' modified gamma size distribution shape parameters'
    end if
    write (3, '(1X,I3,2(1X,F8.3),A)') NRETAB, SRETAB, ERETAB, &
        '  number, starting, ending effective radius'

    do I = 1, NRETAB
        if (NRETAB <= 1) then
            REFF = SRETAB
        else if (LOGRE) then
            REFF = exp((log(ERETAB) - log(SRETAB))*FLOAT(I - 1)/(NRETAB - 1) + log(SRETAB))
        else
            REFF = (ERETAB - SRETAB)*FLOAT(I - 1)/(NRETAB - 1) + SRETAB
        end if
        write (3, '(1X,F8.4,1X,E12.5,1X,F8.6,1X,I6,A)') &
            REFF, EXTINCT(I), SSALB(I), NRANK(I), '  Reff  Ext  Alb  Nrank'
        do J = 1, NPHASEPOL
            if (POLTAB) then
                write (3, '(I1,1X,201(1X,F10.5))') J, (WIGCOEF(J, L, I), L=0, min(NRANK(I), 200))
            else
                write (3, '(2X,201(1X,F10.5))') (WIGCOEF(J, L, I), L=0, min(NRANK(I), 200))
            end if
            do K = 200, NRANK(I) - 1, 200
                write (3, '(2X,200(1X,F10.5))') (WIGCOEF(J, K + L, I), L=1, min(200, NRANK(I) - K))
            end do
        end do
    end do
    close (3)
end subroutine WRITE_MIE_TABLE
