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
!   Hong Chen, hong.chen-1@colorado.edu
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
! SHDOM was originally developed by Frank Evans and is available at https://coloradolinux.com/shdom.
!╰────────────────────────────────────────────────────────────────────────────╯!

program PLOTSCATTAB
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
    implicit none
    integer :: NANGLE, NOUT, NPHASE, NPHASEPOL, MAXNCOEF
    integer :: I, J, K, L
    logical :: POLTAB, FOUND
    real    :: PI, PL, PL1, PL2, SCAT, MU
    integer, allocatable :: NCOEF(:), PELEM(:)
    real, allocatable :: ANGLE(:), OUTREFF(:), REFF(:)
    real, allocatable :: WIGCOEF(:, :, :), PHASE(:, :)
    character(LEN=1)  :: FILETYPE
    character(LEN=256) :: INFILE, PLOTFILE
    character(LEN=7)  :: ELEMNAMES(-6:6)

    ELEMNAMES(:) = (/'P34/P11', 'P12/P11', 'P44/P11', 'P33/P11', 'P22/P11', &
                     '  P11  ', '  P11  ', '  P11  ', &
                     '  P22  ', '  P33  ', '  P44  ', '  P12  ', '  P34  '/)

    write (*, *) 'Scattering table or SHDOM property file input (S or P)'
    read (*, '(A1)') FILETYPE
    write (*, '(1X,A1)') FILETYPE

    write (*, *) 'Input file name'
    read (*, '(A)') INFILE
    write (*, '(A)') INFILE

    write (*, '(1X,A)') 'Number of angles (<0 to input |n| angles)'
    read (*, *) NANGLE
    write (*, *) NANGLE
    allocate (ANGLE(abs(NANGLE)))
    if (NANGLE <= 0) then
        NANGLE = abs(NANGLE)
        write (*, *) 'Input the angles (degrees)'
        do J = 1, NANGLE
            write (*, '(1X,I3,A)') J, ' : '
            read (*, *) ANGLE(J)
        end do
    else
        ANGLE(:) = (/(180.0*FLOAT(J - 1)/(NANGLE - 1), J=1, NANGLE)/)
    end if

    if (FILETYPE == 'P') then
        write (*, '(1X,A)') 'Number of phase functions to output'
        read (*, *) NOUT
        write (*, *) NOUT
        allocate (OUTREFF(NOUT))
        write (*, *) 'Input the tabulated phase function indices'
        read (*, *) OUTREFF(:)
        write (*, '(20(1X,I4))') nint(OUTREFF(:))
        write (*, *) 'Phase matrix elements to output for each phase function index'
    else
        write (*, '(1X,A)') 'Number of effective radii/matrix elements to output'
        read (*, *) NOUT
        write (*, *) NOUT
        allocate (OUTREFF(NOUT))
        write (*, *) 'Input the effective radii (micron)'
        read (*, *) OUTREFF(:)
        write (*, '(20(1X,F6.2))') OUTREFF(:)
        write (*, *) 'Phase matrix elements to output for each effective radius'
    end if
    allocate (PELEM(NOUT))
    write (*, *) '(1=P11, 2=P22, 3=P33, 4=P44, 5=P12, 6=P34) (<0 for ratio to P11)'
    read (*, *) PELEM(:)
    write (*, '(20(1X,I2))') PELEM(:)

    write (*, *) 'Plotting output file name'
    read (*, '(A)') PLOTFILE
    write (*, '(A)') PLOTFILE

    PELEM(:) = min(6, max(-6, PELEM(:)))
    where (PELEM(:) == -1 .or. PELEM(:) == 0)
        PELEM(:) = 1
    end where

    if (FILETYPE == 'P') then
        call READ_PROPERTY_SIZE(INFILE, NPHASEPOL, NPHASE, MAXNCOEF)
        allocate (REFF(NPHASE), NCOEF(NPHASE), WIGCOEF(NPHASEPOL, 0:MAXNCOEF, NPHASE))
        call READ_PROPERTY_PHASEMAT(INFILE, NPHASEPOL, NPHASE, MAXNCOEF, &
                                    NCOEF, WIGCOEF)
        REFF(:) = (/(FLOAT(J), J=1, NPHASE)/)
    else
        call READ_SCAT_TABLE_SIZE(INFILE, POLTAB, NPHASEPOL, NPHASE, MAXNCOEF)
        allocate (REFF(NPHASE), NCOEF(NPHASE), WIGCOEF(NPHASEPOL, 0:MAXNCOEF, NPHASE))
        call READ_SCAT_TABLE(INFILE, POLTAB, NPHASEPOL, NPHASE, REFF, &
                             MAXNCOEF, NCOEF, WIGCOEF)
    end if

    ! Loop over each effective radius in scattering table (I), pulling out
    !   ones we want to output (K)
    allocate (PHASE(NANGLE, NOUT))
    PI = acos(-1.0)
    do K = 1, NOUT
        FOUND = .false.
        do I = 1, NPHASE
            if (REFF(I) == OUTREFF(K)) then
                FOUND = .true.
                call TRANSFORM_WIGNERD_TO_PHASE(MAXNCOEF, NPHASEPOL, &
                                                PELEM(K), NCOEF(I), WIGCOEF(:, :, I), &
                                                NANGLE, ANGLE, PHASE(:, K))
            end if
        end do
        if (.not. FOUND) then
            print *, 'Phase function not found in input file:', OUTREFF(K)
            PHASE(:, K) = 1.0
        end if
    end do

    ! Output the phase functions
    open (UNIT=1, FILE=trim(adjustl(PLOTFILE)), STATUS='UNKNOWN')
    if (FILETYPE == 'P') then
        write (1, '(A,A40)') '!  Property file phase functions: ', INFILE
        write (1, '(A)') '! Angle  cos(angle)  Phase functions with index numbers'
        write (1, '(A,20(6X,I3,3X))') '!               ', (nint(OUTREFF(K)), K=1, NOUT)
    else
        write (1, '(A,A40)') '!  Scattering table phase functions: ', INFILE
        write (1, '(A)') '! Angle  cos(angle)  Phase functions for effective radii (um) / matrix element'
        write (1, '(A,20(6X,F6.2))') '!               ', (OUTREFF(K), K=1, NOUT)
    end if
    write (1, '(A,20(5X,A7))') '!                ', (ELEMNAMES(PELEM(K)), K=1, NOUT)
    do J = 1, NANGLE
        write (1, '(1X,F7.2,1X,F9.6,20(E12.4))') &
            ANGLE(J), cos(ANGLE(J)*PI/180.0), (PHASE(J, K), K=1, NOUT)
    end do
    close (1)

    deallocate (PHASE, REFF, NCOEF, WIGCOEF, OUTREFF, ANGLE, PELEM)
end program

subroutine READ_SCAT_TABLE_SIZE(SCATTABFILE, POLTAB, NPHASEPOL, &
                                NRETAB, MAXNCOEF)
    ! Reads the scattering table file to find the number of effective radii
    ! and the maximum order of the Wigner d-function series phase functions.
    ! Also returns whether the scattering table is polarized in POLTAB.
    implicit none
    character(LEN=256), intent(IN) :: SCATTABFILE
    integer, intent(OUT) :: NPHASEPOL, NRETAB, MAXNCOEF
    logical, intent(OUT) :: POLTAB
    integer :: I, J, L, NCOEF, IERR
    real    :: REFF, EXT, SSALB, CHI

    open (UNIT=3, FILE=trim(adjustl(SCATTABFILE)), STATUS='OLD')
    read (3, *)
    read (3, *)
    read (3, *)
    read (3, *)
    read (3, *)
    read (3, *) NRETAB
    ! See if this is a polarized scattering table
    read (3, *) REFF, EXT, SSALB, NCOEF
    POLTAB = .true.
    do I = 1, 6
        read (3, *, IOSTAT=IERR) J, (CHI, L=0, NCOEF)
        if (IERR > 0) then
            POLTAB = .false.
            exit
        end if
        if (J /= I) POLTAB = .false.
    end do
    rewind (3)
    read (3, *)
    read (3, *)
    read (3, *)
    read (3, *)
    read (3, *)
    read (3, *) NRETAB
    ! Read the rest of the table to get the maximum Ncoef
    MAXNCOEF = 1
    do I = 1, NRETAB
        read (3, *) REFF, EXT, SSALB, NCOEF
        MAXNCOEF = max(MAXNCOEF, NCOEF)
        if (POLTAB) then
            read (3, *) J, (CHI, L=0, NCOEF)
            read (3, *) J, (CHI, L=0, NCOEF)
            read (3, *) J, (CHI, L=0, NCOEF)
            read (3, *) J, (CHI, L=0, NCOEF)
            read (3, *) J, (CHI, L=0, NCOEF)
            read (3, *) J, (CHI, L=0, NCOEF)
        else
            read (3, *) (CHI, L=0, NCOEF)
        end if
    end do
    close (3)
    if (POLTAB) then
        NPHASEPOL = 6
    else
        NPHASEPOL = 1
    end if
end subroutine READ_SCAT_TABLE_SIZE

subroutine READ_SCAT_TABLE(SCATTABFILE, POLTAB, NPHASEPOL, NRETAB, RETAB, &
                           MAXNCOEF, NCOEF, WIGCOEF)
    ! Reads the table of scattering properties as a function of effective
    ! radius.  POLTAB=.true. means a polarized scattering table with Wigner
    ! d-function coefficients for six phase matrix elements instead of one.
    ! The order of the six elements is the four diagonal ones (a1, a2, a3, a4)
    ! followed by the IQ and UV coefficients (b1, b2).
    implicit none
    character(LEN=256), intent(IN) :: SCATTABFILE
    integer, intent(IN)  :: NPHASEPOL, NRETAB, MAXNCOEF
    logical, intent(IN)  :: POLTAB
    integer, intent(OUT) :: NCOEF(NRETAB)
    real, intent(OUT) :: RETAB(NRETAB)
    real, intent(OUT) :: WIGCOEF(NPHASEPOL, 0:MAXNCOEF, NRETAB)
    integer :: I, J, L
    real    :: CHI, EXTINCT, SSALB

    open (UNIT=1, FILE=trim(adjustl(SCATTABFILE)), STATUS='OLD')
    read (1, *)
    read (1, *)
    read (1, *)
    read (1, *)
    read (1, *)
    read (1, *)
    do I = 1, NRETAB
        read (1, *) RETAB(I), EXTINCT, SSALB, NCOEF(I)
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
            print *, 'Error [READ_SCAT_TABLE]: Incorrect Legendre series; chi_0 is not 1'
            stop
        end if
        if (I > 1 .and. RETAB(I) <= RETAB(I - 1)) then
            print *, 'Error [READ_SCAT_TABLE]: Effective radius not increasing in table:', SCATTABFILE
            stop
        end if
    end do
    close (1)
end subroutine READ_SCAT_TABLE

subroutine READ_PROPERTY_SIZE(PROPFILE, NPHASEPOL, NUMPHASE, MAXNCOEF)
    ! Reads the header of the tabulated phase function property file to
    ! get the maximum array sizes needed for allocatable arrays.
    implicit none
    character(LEN=256), intent(IN) :: PROPFILE
    integer, intent(OUT) :: NPHASEPOL, NUMPHASE, MAXNCOEF
    integer :: NPX, NPY, NPZ, NUML, I, J, K, L
    real    :: DELX, DELY, ZLEVELS, CHI
    character(LEN=1) :: PROPTYPE

    ! Open the file, figure out the type, and get the grid size
    open (UNIT=2, FILE=trim(adjustl(PROPFILE)), STATUS='OLD')
    read (2, '(A1)') PROPTYPE
    if (PROPTYPE /= 'T' .and. PROPTYPE /= 'P') then
        print *, 'Error [READ_PROPERTY_SIZE]: Must be a tabulated phase function property file.'
        stop
    end if
    read (2, *) NPX, NPY, NPZ
    read (2, *) DELX, DELY, (ZLEVELS, K=1, NPZ)
    read (2, *) NUMPHASE
    if (PROPTYPE .eq. 'T') then
        ! Property file type T is for tabulated phase function format
        NPHASEPOL = 1
        do I = 1, NUMPHASE
            read (2, *) NUML, (CHI, L=1, NUML)
            MAXNCOEF = max(NUML, MAXNCOEF)
        end do
    else if (PROPTYPE .eq. 'P') then
        ! Property file type P is for polarized tabulated phase function format
        NPHASEPOL = 6
        do I = 1, NUMPHASE
            do J = 1, 6
                read (2, *) K, NUML, (CHI, L=0, NUML)
                if (K /= J) then
                    print *, 'Error [READ_PROPERTY_SIZE]: Error reading  polarized tabulated phase function format.'
                    stop
                end if
                MAXNCOEF = max(NUML, MAXNCOEF)
            end do
        end do
    end if
    close (2)
end subroutine READ_PROPERTY_SIZE

subroutine READ_PROPERTY_PHASEMAT(PROPFILE, NPHASEPOL, NUMPHASE, MAXNCOEF, &
                                  NCOEF, WIGCOEF)
    ! Reads the Wigner d-function phase matrix coefficients from the header of
    ! a tabulated or polarized tabulated SHDOM optical property file.
    implicit none
    character(LEN=256), intent(IN) :: PROPFILE
    integer, intent(IN)  :: NPHASEPOL, NUMPHASE, MAXNCOEF
    integer, intent(OUT) :: NCOEF(NUMPHASE)
    real, intent(OUT) :: WIGCOEF(NPHASEPOL, 0:MAXNCOEF, NUMPHASE)
    integer :: NPX, NPY, NPZ, I, J, K, L
    real    :: DELX, DELY, ZLEVELS
    character(LEN=1) :: PROPTYPE

    open (UNIT=2, FILE=trim(adjustl(PROPFILE)), STATUS='OLD')
    read (2, '(A1)') PROPTYPE
    read (2, *) NPX, NPY, NPZ
    read (2, *) DELX, DELY, (ZLEVELS, K=1, NPZ)
    read (2, *)
    do I = 1, NUMPHASE
        if (PROPTYPE .eq. 'T') then
            read (2, *) NCOEF(I), (WIGCOEF(1, L, I), L=1, NCOEF(I))
            WIGCOEF(1, 0, I) = 1.0
        else
            do J = 1, 6
                read (2, *) K, NCOEF(I), (WIGCOEF(K, L, I), L=0, NCOEF(I))
            end do
        end if
    end do
    close (2)
end subroutine READ_PROPERTY_PHASEMAT

subroutine TRANSFORM_WIGNERD_TO_PHASE(MAXNCOEF, NPHASEPOL, &
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
    implicit none
    integer, intent(IN) :: MAXNCOEF, NPHASEPOL, PELEM, NCOEF, NANGLE
    real, intent(IN) :: WIGCOEF(NPHASEPOL, 0:MAXNCOEF), ANGLE(NANGLE)
    real, intent(OUT) :: PHASE(NANGLE)
    integer :: J, N
    double precision :: PI, X, XM, A1, A2, A3, A4, B1, B2, A2P3, A2M3
    double precision, allocatable :: D00(:), D22P(:), D22M(:), D20(:)
    real, allocatable :: NORM(:)

    if (NPHASEPOL == 1 .and. PELEM > 1) then
        print '(A,I1,A)', 'Error [TRANSFORM_WIGNERD_TO_PHASE]: Cannot do phase matrix element ', PELEM,&
            &' for unpolarized scattering table.'
        stop
    end if

    ! Allocate arrays for Wigner functions D20, D22, D2-2 and D00
    allocate (D00(0:NCOEF), D22P(0:NCOEF), D22M(0:NCOEF), D20(0:NCOEF))
    allocate (NORM(NANGLE))

    PI = DACOS(-1.0d0)
    do J = 1, NANGLE
        X = DCOS(ANGLE(J)*PI/180)
        XM = -X

        ! Calculate the desired scattering matrix element for this angle (X)
        ! from the Wigner d-function coefficients
        select case (abs(PELEM))
        case (1)
            call WIGNERFCT(X, NCOEF, 0, 0, D00)
            A1 = 0.0d0
            do N = 0, NCOEF
                A1 = A1 + WIGCOEF(1, N)*D00(N)
            end do
            PHASE(J) = A1

        case (2:3)
            call WIGNERFCT(X, NCOEF, 2, 2, D22P)
            call WIGNERFCT(XM, NCOEF, 2, 2, D22M) ! multiply by (-1)**N
            A2P3 = 0.0d0
            A2M3 = 0.0d0
            do N = 2, NCOEF
                A2P3 = A2P3 + (WIGCOEF(2, N) + WIGCOEF(3, N))*D22P(N)
                A2M3 = A2M3 + (WIGCOEF(2, N) - WIGCOEF(3, N))*(-1)**N*D22M(N)
            end do
            A2 = 0.5d0*(A2P3 + A2M3)
            A3 = 0.5d0*(A2P3 - A2M3)
            if (abs(PELEM) == 2) then
                PHASE(J) = A2
            else
                PHASE(J) = A3
            end if

        case (4)
            call WIGNERFCT(X, NCOEF, 0, 0, D00)
            A4 = 0.0d0
            do N = 0, NCOEF
                A4 = A4 + WIGCOEF(4, N)*D00(N)
            end do
            PHASE(J) = A4

        case (5)
            call WIGNERFCT(X, NCOEF, 2, 0, D20)
            B1 = 0.0d0
            do N = 2, NCOEF
                B1 = B1 - WIGCOEF(5, N)*D20(N)
            end do
            PHASE(J) = B1

        case (6)
            call WIGNERFCT(X, NCOEF, 2, 0, D20)
            B2 = 0.0d0
            do N = 2, NCOEF
                B2 = B2 - WIGCOEF(6, N)*D20(N)
            end do
            PHASE(J) = B2

        case DEFAULT
            print *, 'Error [TRANSFORM_WIGNERD_TO_PHASE]: Illegal PELEM phase matrix element number.'
            stop
        end select

        ! Calculate P11 function if normalization by it is desired
        if (PELEM < 0) then
            call WIGNERFCT(X, NCOEF, 0, 0, D00)
            A1 = 0.0d0
            do N = 0, NCOEF
                A1 = A1 + WIGCOEF(1, N)*D00(N)
            end do
            NORM(J) = A1
        end if
    end do

    ! Do the P11 normalization if desired
    if (PELEM < 0) then
        PHASE(:) = PHASE(:)/NORM(:)
    end if
    deallocate (D00, D22P, D22M, D20, NORM)
end subroutine TRANSFORM_WIGNERD_TO_PHASE

subroutine WIGNERFCT(X, NRANK, M, M1, DMM1)
    ! The routine computes the DMM1N vector coefficients for
    ! M >= 0, M1 >= 0,  N = 0,...,NRANK, and -1 < X= COS(BETA) < 1
    implicit none
    integer, intent(IN) :: NRANK, M, M1
    double precision, intent(IN) :: X
    double precision, intent(OUT) :: DMM1(0:NRANK)
    integer :: N0, N
    double precision :: FACT1, FACT2, DMM1_N0

    N0 = max(M, M1)
    DMM1 = 0.d0
    if (N0 .eq. 0) then
        DMM1(0) = 1.d0
        DMM1(1) = X
        do N = 1, NRANK - 1
            FACT1 = dble(2*N + 1)*X/dble(N + 1)
            FACT2 = dble(N)/dble(N + 1)
            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    else
        DMM1(N0) = DMM1_N0(X, M, M1)
        do N = N0, NRANK - 1
            FACT1 = dble(N*(N + 1))*X - dble(M*M1)
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M**2))
            FACT1 = FACT1/DSQRT(dble((N + 1)**2 - M1**2))
            FACT1 = FACT1*dble(2*N + 1)/dble(N)
            FACT2 = DSQRT(dble(N**2 - M**2)) &
                    *DSQRT(dble(N**2 - M1**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M**2))
            FACT2 = FACT2/DSQRT(dble((N + 1)**2 - M1**2))
            FACT2 = FACT2*dble(N + 1)/dble(N)

            DMM1(N + 1) = FACT1*DMM1(N) - FACT2*DMM1(N - 1)
        end do
    end if
end subroutine WIGNERFCT

function DMM1_N0(X, M, M1) result(DMM1N0)
    ! The routine computes the Wigner functions for
    !   M >= 0, M1 >= 0, AND N0 = MAX(M,M1)
    implicit none
    integer, intent(IN) :: M, M1
    double precision, intent(IN) :: X
    double precision :: DMM1N0
    integer :: P, MAXM, MINM
    double precision :: PROD, FACT, CMM1

    if (M .eq. M1) then
        FACT = ((1.d0 + X)/2.d0)**M
        DMM1N0 = FACT
    else
        if (M1 .gt. M) then
            CMM1 = 1.d0
        else
            CMM1 = (-1.d0)**(M - M1)
        end if
        MAXM = max(M, M1)
        MINM = min(M, M1)
        PROD = 1.d0
        do P = 1, MAXM - MINM
            FACT = DSQRT(dble(M + M1 + P)/dble(P))
            PROD = PROD*FACT
        end do
        FACT = DSQRT((1.d0 - X)/2.d0)
        DMM1N0 = CMM1*PROD*FACT**(MAXM - MINM)
        FACT = DSQRT((1.d0 + X)/2.d0)
        DMM1N0 = DMM1N0*FACT**(MAXM + MINM)
    end if
end function DMM1_N0
