!╭────────────────────────────────────────────────────────────────────────────╮!
! @LICENSE_HEADER_START@
!
!   This file is part of SHDOMNG.
!
!   --
!   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
!   --
!
!   Copyright (C) 2025-
!   EaR³T Community: https://discord.gg/ntqsguwaWv
!
!   SHDOMNG is free software: you can redistribute it and/or modify
!   it under the terms of the GNU General Public License as published by
!   the Free Software Foundation, either version 3 of the License, or
!   (at your option) any later version.
!
!   SHDOMNG is distributed in the hope that it will be useful,
!   but WITHOUT ANY WARRANTY; without even the implied warranty of
!   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!   GNU General Public License for more details.
!
!   You should have received a copy of the GNU General Public License
!   along with SHDOMNG. If not, see <http://www.gnu.org/licenses/>.
!
! @LICENSE_HEADER_END@
!
! Acknowledgement:
!   SHDOM was originally developed by K. Franklin Evans and is available
!   at https://coloradolinux.com/shdom.
!
! History:
!   Created by:  K. Franklin Evans
!                evans@nit.colorado.edu
!
!   Modified by: Hong Chen
!                hong.chen-1@colorado.edu
!╰────────────────────────────────────────────────────────────────────────────╯!

program MAKE_TMATRIX_TABLE
    !
    ! Does T-matrix computations to create a scattering table as a function of
    ! equivalent volume effective radius for gamma, modified gamma, or lognormal
    ! size distributions of randomly oriented spheroidal or cylindrical particles
    ! with a specified axis ratio
    !   The particles may be water or ice (in which case the program provides
    ! the index of refraction depending on wavelength) or "aerosols" (in which
    ! case the index of refraction is user specified).
    !   The phase functions in the output scattering table are represented
    ! with Legendre series.  For polarized output, the six unique elements of
    ! the phase matrix are represented with Wigner d-function expansions
    ! (Doicu et al., 2013, JQSRT, http://dx.doi.org/10.1016/j.jqsrt.2012.12.009).
    !
    !  compile: pgf90 -fast -o make_tmatrix_table  make_tmatrix_table.f90
    !             indexwatice.f  tmatrixwig.f lpd.f  (uses tmd.par.f include file)
    !
    !    Frank Evans    University of Colorado       January 2014
    implicit none
    integer :: NRETAB
    logical :: POLTAB, CYLFLAG, LOGRE
    real    :: WAVELEN, PARDENS, AXRATIO
    real    :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
    complex :: RINDEX
    character(LEN=1) :: PARTYPE, DISTFLAG
    character(LEN=256) :: TABLEFILE
    integer :: NSIZE, MAXRANK, I, J, L, NL
    real    :: MRE, MIM, A, PI, XMAX, REFF, SCATTER, ALBEDO, ASYM
    integer, allocatable :: NRANK1(:), NRANK(:)
    real, allocatable :: RADII(:), ND(:)
    real, allocatable :: QEXT(:), QSCA(:)
    real, allocatable :: EXTINCT1(:), SCATTER1(:), WIGCOEF1(:, :, :)
    real, allocatable :: EXTINCT(:), SSALB(:), WIGCOEF(:, :, :)

    ! Get the input parameters
    call USER_INPUT(POLTAB, WAVELEN, PARTYPE, RINDEX, PARDENS,&
                  & CYLFLAG, AXRATIO, DISTFLAG, ALPHA, GAMMA,&
                  & NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS, TABLEFILE)
    if (PARTYPE == 'W') then
        PARDENS = 1.0
    else if (PARTYPE == 'I') then
        PARDENS = 0.916
    end if

    ! Calculate the maximum size parameter and the max number of Legendre terms
    PI = acos(-1.0)
    XMAX = 2*PI*MAXRADIUS/WAVELEN
    MAXRANK = nint(2*(XMAX + 4.0*XMAX**0.3334 + 2))

    ! Get the index of refraction for water or ice
    if (PARTYPE == 'W') then
        call REFWAT(0, WAVELEN, 283.0, MRE, MIM, A, A)
        RINDEX = cmplx(MRE, -MIM)
    end if
    if (PARTYPE == 'I') then
        call REFICE(0, WAVELEN, 243.0, MRE, MIM, A, A)
        RINDEX = cmplx(MRE, -MIM)
    end if

    ! Figure the number of radii there will be
    call GET_NSIZE(SRETAB, MAXRADIUS, WAVELEN, NSIZE)

    ! Allocate all the arrays here
    allocate(RADII(NSIZE), ND(NSIZE), NRANK1(NSIZE))
    allocate(EXTINCT1(NSIZE), SCATTER1(NSIZE), WIGCOEF1(6, 0:MAXRANK, NSIZE))
    allocate(EXTINCT(NRETAB), SSALB(NRETAB))
    allocate(NRANK(NRETAB), WIGCOEF(6, 0:MAXRANK, NRETAB))

    ! Make up the discrete particle radii to use
    call GET_SIZES(SRETAB, MAXRADIUS, WAVELEN, NSIZE, RADII)

    ! Do the T-matrix computations for each radius
    do I = 1, NSIZE
        call TMATRIXWIG(WAVELEN, RINDEX, RADII(I), CYLFLAG, AXRATIO, MAXRANK,&
                      & EXTINCT1(I), SCATTER1(I), ALBEDO, ASYM,&
                      & NRANK1(I), WIGCOEF1(1, 0, I))
        ! PRINT '(I4,1X,F7.3,1X,ES9.3,2(1X,F6.4),I4)',&
            ! & I, RADII(I), EXTINCT1(I), ALBEDO, ASYM, NRANK1(I)
    end do

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
        call MAKE_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA,&
                          & ND)

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
            print *, 'Error [MAKE_TMATRIX_TABLE]: Phase function not normalized for Reff=', REFF, WIGCOEF(1, 0, I)
            stop
        end if
        if (EXTINCT(I) > 0.0) then
            SSALB(I) = SCATTER/EXTINCT(I)
        end if
        EXTINCT(I) = 0.001*EXTINCT(I)

    end do ! end of effective radius loop

    call WRITE_SCAT_TABLE(TABLEFILE, POLTAB, WAVELEN,&
                        & PARTYPE, PARDENS, RINDEX, CYLFLAG, AXRATIO,&
                        & DISTFLAG, ALPHA, GAMMA, NRETAB, SRETAB, ERETAB, LOGRE,&
                        & EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)

    deallocate(EXTINCT, SSALB, NRANK, WIGCOEF)
    deallocate(ND, EXTINCT1, SCATTER1, NRANK1, WIGCOEF1)
end program

subroutine USER_INPUT(POLTAB, WAVELEN, PARTYPE, RINDEX, PARDENS,&
                    & CYLFLAG, AXRATIO, DISTFLAG, ALPHA, GAMMA,&
                    & NRETAB, SRETAB, ERETAB, LOGRE, MAXRADIUS, TABLEFILE)
    ! Reads the input parameters from the standard input
    implicit none
    integer, intent(OUT) :: NRETAB
    logical, intent(OUT) :: POLTAB, CYLFLAG, LOGRE
    real, intent(OUT) :: WAVELEN, PARDENS, AXRATIO
    real, intent(OUT) :: SRETAB, ERETAB, ALPHA, GAMMA, MAXRADIUS
    complex, intent(OUT) :: RINDEX
    character(LEN=1), intent(OUT) :: PARTYPE, DISTFLAG
    character(LEN=256), intent(OUT) :: TABLEFILE

    write (*, *) 'Making T-matrix scattering tables for spheroidal particles'

    write (*, *) 'Make polarized scattering table (T or F)'
    read (*, *) POLTAB
    write (*, '(L)') POLTAB

    write (*, *) 'Wavelength (micron)'
    read (*, *) WAVELEN
    write (*, '(1X,F9.3)') WAVELEN

    write (*, *) 'Water, Ice, or Aerosol spherical particles (W,I,A)'
    read (*, '(A1)') PARTYPE
    write (*, '(1X,A1)') PARTYPE

    if (PARTYPE /= 'W' .and. PARTYPE /= 'I') then
        PARTYPE = 'A'
        write (*, *) 'Aerosol complex index of refraction'
        read (*, *) RINDEX
        write (*, *) RINDEX
        write (*, *) 'Aerosol particle bulk density (g/cm^3)'
        read (*, *) PARDENS
        write (*, '(1X,F5.3)') PARDENS
    end if

    write (*, *) 'Cylindrical particle (T or F) (F for spheroidal)'
    read (*, *) CYLFLAG
    write (*, '(1X,L1)') CYLFLAG

    write (*, *) 'Particle axis ratio (>1 for oblate spheroid, <1 for prolate)'
    read (*, *) AXRATIO
    write (*, '(1X,F6.4)') AXRATIO

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
        write (*, *) 'Error [USER_INPUT]: Unrecognized size distribution type.'
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

    write (*, *) 'Output scattering table name'
    read (*, '(A)') TABLEFILE
    write (*, '(A)') TABLEFILE
end subroutine USER_INPUT

subroutine GET_NSIZE(SRETAB, MAXRADIUS, WAVELEN, NSIZE)
    ! Calculates the number of radii for which the T-matrix computation will be run.
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
        ! DELX = 0.1                     ! One alternative method
        DELRAD = min(DELX*WAVELEN/TWOPI, DELRADMAX)
        RAD = RAD + DELRAD
        NSIZE = NSIZE + 1
    end do
end subroutine GET_NSIZE

subroutine GET_SIZES(SRETAB, MAXRADIUS, WAVELEN, NSIZE, RADII)
    ! Calculates the radii for which the T-matrix computation will be run and
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
        ! DELX = 0.1                     ! One alternative method
        DELRAD = min(DELX*WAVELEN/TWOPI, DELRADMAX)
        RAD = RAD + DELRAD
        RADII(N) = RAD
    end do
end subroutine GET_SIZES

subroutine MAKE_SIZE_DIST(DISTFLAG, PARDENS, NSIZE, RADII, REFF, ALPHA, GAMMA,&
                        & ND)
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
    call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, REFF,&
                    & NSIZE, RADII, ND, TRUERE)
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
            call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, REHI,&
                            & NSIZE, RADII, ND, TRUERE)
        end do
        if (TRUERE <= REFF) then
            print *, 'Error [MAKE_SIZE_DIST]: Effective radius cannot be achieved', REFF, TRUERE
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
            call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, RELO,&
                            & NSIZE, RADII, ND, TRUERE)
        end do
        if (TRUERE >= REFF) then
            print *, 'Error [MAKE_SIZE_DIST]: Effective radius cannot be achieved', REFF, TRUERE
            stop
        end if
    end if
    ! Do bisection to get correct effective radius
    do while (abs(TRUERE - REFF) > TOL*REFF)
        REMID = 0.5*(RELO + REHI)
        call DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, REMID,&
                        & NSIZE, RADII, ND, TRUERE)
        if (TRUERE < REFF) then
            RELO = REMID
        else
            REHI = REMID
        end if
    end do
end subroutine MAKE_SIZE_DIST

subroutine DO_SIZE_DIST(PARDENS, DISTFLAG, ALPHA, GAMMA, RE, NSIZE, RADII,&
                      & ND, TRUERE)
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
        A = 1.e6*GAMMA*B**((ALPHA + 4)/GAMMA)&
            & /((4*PI/3.)*PARDENS*exp(GAMMLN((ALPHA + 4.)/GAMMA)))
    else if (DISTFLAG == 'L') then
        B = RE*exp(-2.5*ALPHA**2)
        A = 1.e6/((4*PI/3.)*PARDENS*sqrt(2*PI)*ALPHA*B**3*exp(4.5*ALPHA**2))
    end if
    LWC = 0.0
    SUM2 = 0.0
    SUM3 = 0.0
    do J = 1, NSIZE
        R = RADII(J)
        DELR = sqrt(RADII(J)*RADII(min(NSIZE, J + 1)))&
            & - sqrt(RADII(J)*RADII(max(1, J - 1)))
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

subroutine WRITE_SCAT_TABLE(TABLEFILE, POLTAB, WAVELEN,&
                          & PARTYPE, PARDENS, RINDEX, CYLFLAG, AXRATIO,&
                          & DISTFLAG, ALPHA, GAMMA, NRETAB, SRETAB, ERETAB, LOGRE,&
                          & EXTINCT, SSALB, MAXRANK, NRANK, WIGCOEF)
    ! Writes the table of scattering properties as a function of
    ! effective radius.  There are two types of tables output:
    ! unpolarized with Legendre (Wigner d_l_00) series of a single phase
    ! matrix element (P11) or polarized with six Wigner d-function series
    ! of the phase matrix elements.
    implicit none
    integer, intent(IN) :: NRETAB, MAXRANK, NRANK(NRETAB)
    logical, intent(IN) :: POLTAB, CYLFLAG, LOGRE
    complex, intent(IN) :: RINDEX
    real, intent(IN) :: WAVELEN, AXRATIO
    real, intent(IN) :: PARDENS, SRETAB, ERETAB, ALPHA, GAMMA
    real, intent(IN) :: EXTINCT(NRETAB), SSALB(NRETAB)
    real, intent(IN) :: WIGCOEF(6, 0:MAXRANK, NRETAB)
    character(LEN=1), intent(IN) :: PARTYPE, DISTFLAG
    character(LEN=256), intent(IN) :: TABLEFILE
    integer :: I, J, JP, K, L, NL, NPHASEPOL
    real REFF
    character(LEN=8) :: SHAPE

    if (CYLFLAG) then
        SHAPE = 'cylinder'
    else
        SHAPE = 'spheroid'
    end if
    open (UNIT=3, FILE=trim(adjustl(TABLEFILE)), STATUS='REPLACE')
    if (POLTAB) then
        NPHASEPOL = 6
        write (3, '(3A)') '! Polarized T-matrix ', trim(SHAPE), ' scattering table vs. effective radius (LWC=1 g/m^3)'
    else
        NPHASEPOL = 1
        write (3, '(3A)') '! T-matrix ', trim(SHAPE), ' scattering table vs. effective radius (LWC=1 g/m^3)'
    end if
    write (3, '(2(1X,F8.3),A)') WAVELEN, WAVELEN, '  wavelength range (micron)'
    write (3, '(1X,F5.3,1X,F6.4,2X,A1,A)') PARDENS, AXRATIO, PARTYPE,&
        & '   particle density (g/cm^3), axial ratio, and type (Water, Ice, Aerosol)'
    write (3, '(2(1X,E13.6),A)') RINDEX, '  particle index of refraction'
    if (DISTFLAG == 'L') then
        write (3, '(F7.5,A)') ALPHA, '  lognormal log standard deviation'
    else if (DISTFLAG == 'G') then
        write (3, '(F7.5,A)') ALPHA, ' gamma size distribution shape parameter'
    else if (DISTFLAG == 'M') then
        write (3, '(2(1X,F7.5),A)') ALPHA, GAMMA, ' modified gamma size distribution shape parameters'
    end if
    write (3, '(1X,I3,2(1X,F8.3),A)') NRETAB, SRETAB, ERETAB,&
        & '  number, starting, ending effective radius'

    do I = 1, NRETAB
        if (NRETAB <= 1) then
            REFF = SRETAB
        else if (LOGRE) then
            REFF = exp((log(ERETAB) - log(SRETAB))*FLOAT(I - 1)/(NRETAB - 1) + log(SRETAB))
        else
            REFF = (ERETAB - SRETAB)*FLOAT(I - 1)/(NRETAB - 1) + SRETAB
        end if
        write (3, '(1X,F8.4,1X,E12.5,1X,F8.6,1X,I6,A)')&
            & REFF, EXTINCT(I), SSALB(I), NRANK(I), '  Reff  Ext  Alb  Nrank'
        do J = 1, NPHASEPOL
            if (POLTAB) then
                write (3, '(I1,1X,201(1X,E15.6))') J, (WIGCOEF(J, L, I), L=0, min(NRANK(I), 200))
            else
                write (3, '(2X,201(1X,E15.6))') (WIGCOEF(J, L, I), L=0, min(NRANK(I), 200))
            end if
            do K = 200, NRANK(I) - 1, 200
                write (3, '(2X,200(1X,E15.6))') (WIGCOEF(J, K + L, I), L=1, min(200, NRANK(I) - K))
            end do
        end do
    end do
    close (3)
end subroutine WRITE_SCAT_TABLE
