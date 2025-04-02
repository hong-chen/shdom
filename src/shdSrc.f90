module shd_src

    use shd_grid_ind
    use shd_util
    use shd_pha

    implicit none

    public :: shd_src__computeOneSrc
    public :: shd_src__integrateSrc
    public :: shd_src__computeSrc
    public :: shd_src__calSrcPnt
    public :: shd_src__calSrcPnt_

    private :: shd_src__computeSolarStokesVect

contains

subroutine shd_src__computeOneSrc(NSTOKES,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NPTS, DELTAM, MU, PHI, SRCTYPE, SOLARMU, SOLARAZ,&
    &ALBEDO, LEGEN, IPHASE, DIRFLUX,&
    &SHPTR, SOURCE, SOURCE1)

    ! Computes the source function (SOURCE1) in the direction (MU,PHI)
    ! for the whole domain (NX,NZ).  The spherical harmonic source function
    ! series is input in SOURCE.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    integer NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer SHPTR(NPTS + 1)
    integer*2 IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ, MU, PHI
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS), DIRFLUX(NPTS)
    real SOURCE(NSTOKES, *), SOURCE1(NSTOKES, NPTS)
    character(len=1) :: SRCTYPE
    integer I, J, JT, K, L, M, N, IPH, IS, NS, ME, MS
    real OFOURPI, DA, F, A, A1, B1, SECMU0, COSSCAT
    real, allocatable :: YLMDIR(:, :), YLMSUN(:, :)
    double precision, allocatable :: SUNDIRLEG(:)
    real, allocatable :: UNSCLEGEN(:, :), SINGSCAT(:, :)
    logical TMS

    allocate (YLMDIR(NSTLEG, NLM), YLMSUN(NSTLEG, NLM))
    allocate (SUNDIRLEG(0:NLEG))
    allocate (UNSCLEGEN(2, 0:NLEG), SINGSCAT(NSTOKES, NUMPHASE))
    if (abs(MU) .gt. 1.0) stop 'Error [_computeOneSrc]: Bad mu.'

    TMS = DELTAM
    ! TMS = .FALSE.

    ! Precompute Ylm's for output direction and solar direction
    call shd_pha__computeYLM(.false., MU, PHI, ML, MM, NSTLEG, YLMDIR)
    if (SRCTYPE .ne. 'T') then
        OFOURPI = 1.0/(4.0*acos(-1.0))
        call shd_pha__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG,&
            &YLMSUN)
        SECMU0 = 1.0/abs(SOLARMU)
        COSSCAT = SOLARMU*MU + sqrt((1.0d0 - SOLARMU**2)*(1.0d0 - MU**2))&
            &*cos(SOLARAZ - PHI)
        if (NUMPHASE .eq. 0 .and. TMS) then
            ! If non-tabulate phase functions (unpolarized) then
            !   compute the Legendre polynomials for the scattering angle
            !   for the untruncated solar single scattering computation.
            call shd_pha__computeLegendre(dble(COSSCAT), NLEG, SUNDIRLEG)
            do L = 0, NLEG
                SUNDIRLEG(L) = SUNDIRLEG(L)*(2*L + 1)*OFOURPI
            end do
        end if
    end if

    ! If solar source and tabulated phase functions then precompute
    !   the single scattered solar Stokes vector to outgoing direction
    if (SRCTYPE .ne. 'T' .and. TMS .and. NUMPHASE .gt. 0) then
        do IPH = 1, NUMPHASE
            ! Unscaled the needed Legendre coefficients, but divide by 1-f
            !   because extinction is still scaled (for TMS method)
            F = LEGEN(1, ML + 1, IPH)
            do L = 0, NLEG
                if (L .le. ML) then
                    UNSCLEGEN(1, L) = LEGEN(1, L, IPH) + F/(1 - F)
                else
                    UNSCLEGEN(1, L) = LEGEN(1, L, IPH)/(1 - F)
                end if
                if (NSTLEG .gt. 1) then
                    if (L .le. ML) then
                        UNSCLEGEN(2, L) = LEGEN(5, L, IPH)
                    else
                        UNSCLEGEN(2, L) = LEGEN(5, L, IPH)/(1 - F)
                    end if
                end if
            end do
            ! Sum the Wigner series (for A1 and B1), rotate the plane of
            !   polarization to the outgoing direction, and return the Stoke
            call shd_src__computeSolarStokesVect(NLEG, UNSCLEGEN, COSSCAT,&
                &SOLARMU, MU, SOLARAZ - PHI, NSTOKES,&
                &SINGSCAT(1, IPH))
            if (SINGSCAT(1, IPH) .le. 0.0) then
                write (6, *) 'Error [_computeOneSrc]: negative source function',&
                    &' for tabulated phase function:', IPH, MU, PHI, SINGSCAT(1, IPH)
                stop
            end if
        end do
    end if

    ! Loop over all the grid points, computing the source function
    !   at the viewing angle from the spherical harmonic source func
    do I = 1, NPTS
        IS = SHPTR(I)
        NS = SHPTR(I + 1) - IS
        ! Sum over the real generalized spherical harmonic series
        !   of the source function
        SOURCE1(:, I) = 0.0
        do J = 1, NS
            SOURCE1(1, I) = SOURCE1(1, I) + SOURCE(1, IS + J)*YLMDIR(1, J)
        end do
        if (NSTOKES .gt. 1) then
            do J = 1, NS
                SOURCE1(2, I) = SOURCE1(2, I) + SOURCE(2, IS + J)*YLMDIR(2, J)&
                    &+ SOURCE(3, IS + J)*YLMDIR(5, J)
                SOURCE1(3, I) = SOURCE1(3, I) + SOURCE(2, IS + J)*YLMDIR(6, J)&
                    &+ SOURCE(3, IS + J)*YLMDIR(3, J)
            end do
        end if
        if (NSTOKES .eq. 4) then
            do J = 1, NS
                SOURCE1(4, I) = SOURCE1(4, I) + SOURCE(4, IS + J)*YLMDIR(4, J)
            end do
        end if

        ! Special case for solar source and Delta-M
        if (SRCTYPE .ne. 'T' .and. TMS) then
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            ! First subtract off the truncated single scattering
            DA = ALBEDO(I)*DIRFLUX(I)*SECMU0
            J = 1
            do L = 0, ML
                ME = min(L, MM)
                MS = -ME
                A1 = DA*LEGEN(1, L, IPH)
                B1 = DA*LEGEN(5, L, IPH)
                if (J .le. NS) then
                    JT = J
                    do M = MS, ME
                        SOURCE1(1, I) = SOURCE1(1, I) - A1*YLMDIR(1, J)*YLMSUN(1, J)
                        J = J + 1
                    end do
                    if (NSTOKES .gt. 1) then
                        J = JT
                        do M = MS, ME
                            SOURCE1(2, I) = SOURCE1(2, I) - B1*YLMDIR(2, J)*YLMSUN(1, J)
                            SOURCE1(3, I) = SOURCE1(3, I) - B1*YLMDIR(6, J)*YLMSUN(1, J)
                            J = J + 1
                        end do
                    end if
                end if
            end do

            ! Then add in the single scattering contribution for the
            !   original unscaled phase matrix.
            if (NUMPHASE .gt. 0) then
                SOURCE1(:, I) = SOURCE1(:, I) + DA*SINGSCAT(:, IPH)
            else if (NSTOKES .eq. 1) then
                ! If we have a phase function for each grid point (which
                !   should be the unpolarized case) then unscale the Legendre
                !   coefficients for L<=ML and sum over the Legendre series
                !   in cosine of the scattering angle.
                F = LEGEN(1, ML + 1, IPH)
                do L = 0, NLEG
                    if (L .le. ML) then
                        A = DA*(LEGEN(1, L, IPH) + F/(1 - F))
                    else
                        A = DA*LEGEN(1, L, IPH)/(1 - F)
                    end if
                    SOURCE1(1, I) = SOURCE1(1, I) + A*SUNDIRLEG(L)
                end do
            end if
        end if
    end do

    deallocate (YLMDIR, YLMSUN, SUNDIRLEG, UNSCLEGEN, SINGSCAT)

    return

end subroutine shd_src__computeOneSrc

subroutine shd_src__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
    &NX, NY, NZ, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &MU, PHI, GRIDRAD, EXTINCT, SOURCE,&
    &X0, Y0, Z0, XE, YE, ZE, SIDE,&
    &TRANSMIT, RADIANCE, VALIDRAD)

    ! Integrates the source function (SOURCE) through the extinction
    ! field (EXTINCT) backward in the direction opposite to (MU,PHI)
    ! to find the outgoing radiance (RADIANCE) at the point X0,Y0,Z0.
    ! The radiances at the opposite boundary are input in GRIDRAD.
    ! The transmission and radiance of the ray so far (TRANSMIT, RADIANCE)
    ! are input and returned after the integration along with the exitting
    ! ray location (XE,YE,ZE) and side of the domain (1=-X,2=+X,3=-Y,4=+Y,
    ! 5=-Z,6=+Z).

    integer BCFLAG, IPFLAG, NSTOKES, NX, NY, NZ, NPTS, NCELLS, SIDE
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    logical VALIDRAD
    real XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ), GRIDPOS(3, NPTS)
    real MU, PHI
    double precision X0, Y0, Z0, XE, YE, ZE, TRANSMIT, RADIANCE(NSTOKES)
    real EXTINCT(NPTS)
    real SOURCE(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    integer BITX, BITY, BITZ, IOCT, ICELL, INEXTCELL, IFACE
    integer I1, I2, I3, I4, IOPP, NTAU, IT, K
    logical IPINX, IPINY, OPENBCFACE, BTEST
    integer GRIDFACE(4, 6), OPPFACE(6), JFACE, KFACE, IC
    real EXT0(1), EXT1(1), EXTN(1)
    real XM, YM
    real SRCEXT0(NSTOKES), SRCEXT1(NSTOKES), RAD0(NSTOKES)
    double precision PI, CX, CY, CZ, CXINV, CYINV, CZINV
    double precision XN, YN, ZN, XI, YI, ZI
    double precision SO, SOX, SOY, SOZ, EPS
    double precision TAUTOL, TAUGRID, S, DELS
    double precision EXT, TAU, TRANSCELL, ABSCELL
    double precision SRC(NSTOKES), RAD(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8,&
        &1, 2, 3, 4, 5, 6, 7, 8/, OPPFACE/2, 1, 4, 3, 6, 5/

    ! TAUTOL is the maximum optical path for the subgrid intervals
    TAUTOL = 0.1
    EPS = 1.0e-4*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate dire
    PI = acos(-1.0)
    CX = sqrt(1.0 - MU**2)*cos(PHI + PI)
    CY = sqrt(1.0 - MU**2)*sin(PHI + PI)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-5) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup for the ray path direction
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        stop 'Error [_integrateSrc]: Bad MU.'
    end if
    IOCT = 1 + BITX + 2*BITY + 4*BITZ
    XM = 0.5*(XGRID(1) + XGRID(NX))
    YM = 0.5*(YGRID(1) + YGRID(NY))

    ! Start at the desired point, getting the extinction and source
    XE = X0
    YE = Y0
    ZE = Z0
    call shd_grid__findGridCell(NX, NY, NZ, XGRID, YGRID, ZGRID,&
        &NCELLS, TREEPTR, GRIDPTR, CELLFLAGS, GRIDPOS,&
        &BCFLAG, IPFLAG, XE, YE, ZE, ICELL)
    ! print '(A,5(1X,F8.5),1X,I6)', &
    !     &'_integrateSrc:', X0, Y0, Z0, MU, PHI, ICELL
    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
        &XE, YE, ZE, 1, EXTINCT, 1, 1, EXT1)
    call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
        &XE, YE, ZE, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT1)
    SRCEXT1(1) = max(0.0, SRCEXT1(1))

    ! Loop until finding a face with known radiances
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
        ! This cell is independent pixel if IP mode or open boundary
        !   conditions and ray is leaving domain (i.e. not entering)
        IPINX = btest(int(CELLFLAGS(ICELL)), 0) .and.&
            &.not. (btest(BCFLAG, 0) .and.&
            &((CX .gt. 0 .and. XE .lt. XM) .or. (CX .lt. 0 .and. XE .gt. XM)))
        IPINY = btest(int(CELLFLAGS(ICELL)), 1) .and.&
            &.not. (btest(BCFLAG, 1) .and.&
            &((CY .gt. 0 .and. YE .lt. YM) .or. (CY .lt. 0 .and. YE .gt. YM)))

        ! Find boundaries of the current cell
        ! Find the three possible intersection planes (X,Y,Z)
        ! from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9 - IOCT, ICELL)
        ! Get the distances to the 3 planes and select the closest
        !   (always need to deal with the cell that is wrapped)
        if (IPINX .or. CX .eq. 0.0) then
            SOX = 1.0e20
        else
            SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
        end if
        if (IPINY .or. CY .eq. 0.0) then
            SOY = 1.0e20
        else
            SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
        end if
        SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Error [_integrateSrc]: SO<0 ',&
                &MU, PHI, X0, Y0, Z0, XE, YE, ZE, SO, SOX, SOY, SOZ, XGRID(1), YGRID(1),&
                &GRIDPOS(1, IOPP), GRIDPOS(2, IOPP), GRIDPOS(3, IOPP)
            stop
        end if
        if (IPINX) then
            XN = XE
        else
            XN = XE + SO*CX
        end if
        if (IPINY) then
            YN = YE
        else
            YN = YE + SO*CY
        end if
        ZN = ZE + SO*CZ

        ! Find the optical path across the grid cell and figure how
        !   many subgrid intervals to use
        call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
            &XN, YN, ZN, 1, EXTINCT, 1, 1, EXTN)
        TAUGRID = SO*0.5*(EXT1(1) + EXTN(1))
        NTAU = max(1, 1 + int(TAUGRID/TAUTOL))
        DELS = SO/NTAU
        ! Loop over the subgrid cells
        do IT = 1, NTAU
            S = IT*DELS
            if (IPINX) then
                XI = XE
            else
                XI = XE + S*CX
            end if
            if (IPINY) then
                YI = YE
            else
                YI = YE + S*CY
            end if
            ZI = ZE + S*CZ
            ! Interpolate extinction and source function along path
            if (IT .eq. NTAU) then
                EXT0(1) = EXTN(1)
            else
                call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
                    &XI, YI, ZI, 1, EXTINCT, 1, 1, EXT0)
            end if
            call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS,&
                &XI, YI, ZI, NSTOKES, SOURCE, 1, NSTOKES, SRCEXT0)
            SRCEXT0(1) = max(0.0, SRCEXT0(1))

            ! Compute the subgrid radiance: integration of the source fun
            EXT = 0.5*(EXT0(1) + EXT1(1))
            if (EXT .ne. 0.0) then
                TAU = EXT*DELS
                ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - 0.33333333333*TAU))
                TRANSCELL = 1.0 - ABSCELL
                ! Linear extinction, linear source*extinction, to second
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:)) + 0.08333333333*&
                    &(EXT0(1)*SRCEXT1(:) - EXT1(1)*SRCEXT0(:))*DELS&
                    &*(1.0 - 0.05*(EXT1(1) - EXT0(1))*DELS))/EXT
            else
                ABSCELL = 0.0
                TRANSCELL = 1.0
                SRC(:) = 0.0
            end if

            RADIANCE(:) = RADIANCE(:) + TRANSMIT*SRC(:)*ABSCELL
            SRCEXT1(:) = SRCEXT0(:)
            EXT1(1) = EXT0(1)
            TRANSMIT = TRANSMIT*TRANSCELL
        ! End of sub grid cell loop
        end do

        ! Get the intersection face number (i.e. neighptr index)
        if (SOX .le. SOZ .and. SOX .le. SOY) then
            IFACE = 2 - BITX
            JFACE = 1
            OPENBCFACE = btest(int(CELLFLAGS(ICELL)), 0) .and. btest(BCFLAG, 0)
        else if (SOY .le. SOZ) then
            IFACE = 4 - BITY
            JFACE = 2
            OPENBCFACE = btest(int(CELLFLAGS(ICELL)), 1) .and. btest(BCFLAG, 1)
        else
            IFACE = 6 - BITZ
            JFACE = 3
            OPENBCFACE = .false.
        end if
        ! Get the next cell to go to
        INEXTCELL = NEIGHPTR(IFACE, ICELL)
        if (INEXTCELL .lt. 0) then
            call shd_grid__getNextCell(XN, YN, ZN, IFACE, JFACE, ICELL, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
        end if
        ! Get the four grid points on the intersection face for
        !   telling if at the boundary
        !   If going to same or larger face then use previous face
        if (NEIGHPTR(IFACE, ICELL) .ge. 0 .and. .not. OPENBCFACE) then
            KFACE = IFACE
            IC = ICELL
        else
            ! If going to smaller face then use next face (more accurate
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
        end if
        I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
        I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
        I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
        I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
        ! Get the location coordinate
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XN = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YN = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
            else
                ZN = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
            end if
        end if

        ! See if there are valid radiances for the intersection face:
        !   If so, set flag, interpolate radiance, and add in
        !   If not, only prepare for next cell
        if (TRANSMIT .lt. 1.0e-5) then
            VALIDRAD = .true.
        else if (GRIDRAD(1, I1) .ge. 0.0 .and. GRIDRAD(1, I2) .ge. 0.0 .and.&
            &GRIDRAD(1, I3) .ge. 0.0 .and. GRIDRAD(1, I4) .ge. 0.0) then
            VALIDRAD = .true.
            call shd_util__interpTriLin(ICELL, GRIDPTR, GRIDPOS, XI, YI, ZI,&
                &NSTOKES, GRIDRAD, 1, NSTOKES, RAD0)
            RADIANCE(:) = RADIANCE(:) + TRANSMIT*RAD0(:)
        else
            ICELL = INEXTCELL
        end if
        XE = XN
        YE = YN
        ZE = ZN

    end do
    SIDE = IFACE
    return
end subroutine shd_src__integrateSrc

subroutine shd_src__computeSrc(&
    &NSTOKES, ML, MM, NLM,&
    &NSTLEG, NLEG, NUMPHASE, NPTS, FIXSH,&
    &SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
    &IPHASE, PLANCK, DIRFLUX, SHACC, MAXIV,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
    &FIRST, ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM)

    ! Computes the source function (SOURCE) in spherical harmonic space
    !   for all the grid points.  The thermal source and/or solar
    !   pseudo-source (in PLANCK or DIRFLUX) is added to the scattering source
    !   (computed from LEGEN and the spherical harmonic expansion in RADIAIANCE).
    !   Also computes the adaptive spherical harmonic truncation for
    !   the source function, unless FIXSH is true.  The SH truncation is set
    !   so that the higher terms of the source function are smaller in
    !   absolute value than SHACC. The truncation is to the nearest l,
    !   so that a triangular (for M=L) truncation is preserved.  If there
    !   is little or no scattering then the SH truncation can be set to
    !   zero terms.  The source function SH truncation is output in terms
    !   of the pointer SHPTR (note that RSHPTR refers to the RADIANCE array
    !   and OSHPTR to the DELSOURCE array).
    ! To save memory this routine also computes some things for the
    !   series acceleration: the dot products of the difference in source
    !   function between successive iterations and the new difference source
    !   function vector (DELSOURCE).  To avoid using another array and
    !   because of the adaptive truncation pointer stuff, this requires
    !   two extra calculations of the source function.
    ! Computes the sum of squares of the new source function for calculation
    !   of the solution criterion (SOLCRIT, the RMS difference in succesive
    !   source function fields normalized by the RMS of the field).

    integer NSTOKES, NPTS, ML, MM, NLM, NSTLEG, NLEG
    integer NUMPHASE, MAXIV
    integer RSHPTR(NPTS + 1), SHPTR(NPTS + 1), OSHPTR(NPTS + 1)
    integer*2 IPHASE(NPTS)
    logical FIXSH, FIRST, ACCELFLAG
    real SHACC, YLMSUN(NSTLEG, NLM), SOLARMU
    real DELJDOT, DELJOLD, DELJNEW, JNORM
    real ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real PLANCK(NPTS), DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *)
    real SOURCE(NSTOKES, *), DELSOURCE(NSTOKES, *)
    character(len=1) :: SRCTYPE
    integer IS, ISO, IR, I, IPH, J, JS, K, L, LS, M, MS, ME, NS, NR
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)
    real SRCMIN, C, SECMU0, D

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))

    SRCMIN = SHACC
    ! Make source function from scattering integral and real sources:
    !   Solar source: DIRFLUX has F0*exp(-tau)
    !   Thermal source: PLANCK has (1-omega)*B(T)
    SECMU0 = 1.0/abs(SOLARMU)
    C = sqrt(4.0*acos(-1.0))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    ! Compute the dot products of the successive source differences
    if (.not. FIRST) then
        DELJDOT = 0.0
        DELJOLD = 0.0
        DELJNEW = 0.0
        JNORM = 0.0
        do I = 1, NPTS
            ! Compute the temporary source function for this point
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            IR = RSHPTR(I)
            NR = RSHPTR(I + 1) - IR
            if (NR .gt. NLM) then
                write (6, *) 'Error [_computeSrc]: NR>NLM 1', I, IR, RSHPTR(I + 1)
                stop
            end if
            if (NSTOKES .eq. 1) then
                call shd_src__calSrcPnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_src__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Compute the dot product of the new and old source difference,
            !   and the solution criterion (don't use new points for Jnorm).
            if (ACCELFLAG) then
                IS = SHPTR(I)
                ISO = OSHPTR(I)
                NS = min(SHPTR(I + 1) - IS, OSHPTR(I + 1) - ISO)
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJDOT = DELJDOT +&
                            &(SOURCET(K, J) - SOURCE(K, IS + J))*DELSOURCE(K, ISO + J)
                        DELJOLD = DELJOLD + DELSOURCE(K, ISO + J)**2
                        DELJNEW = DELJNEW + (SOURCET(K, J) - SOURCE(K, IS + J))**2
                        JNORM = JNORM + SOURCE(K, IS + J)**2
                    end do
                end do
            else
                IS = SHPTR(I)
                NS = SHPTR(I + 1) - IS
                do K = 1, NSTOKES
                    do J = 1, NS
                        DELJNEW = DELJNEW + (SOURCET(K, J) - SOURCE(K, IS + J))**2
                        JNORM = JNORM + SOURCE(K, IS + J)**2
                    end do
                end do
            end if
        end do
    end if

    ! Compute and store the new source function difference vector
    if (.not. FIRST .and. ACCELFLAG) then
        do I = 1, NPTS
            ! Compute the temporary source function for this point
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(I)
            else
                IPH = I
            end if
            IR = RSHPTR(I)
            NR = RSHPTR(I + 1) - IR
            if (NSTOKES .eq. 1) then
                call shd_src__calSrcPnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_src__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                    &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if

            ! Make the source difference vector
            IS = SHPTR(I)
            NS = SHPTR(I + 1) - IS
            OSHPTR(I) = IS
            do J = 1, NS
                DELSOURCE(:, IS + J) = SOURCET(:, J) - SOURCE(:, IS + J)
            end do
        end do
        OSHPTR(NPTS + 1) = SHPTR(NPTS + 1)
    end if

    ! Compute the source function in the temporary array.
    IS = 0
    do I = 1, NPTS
        if (NUMPHASE .gt. 0) then
            IPH = IPHASE(I)
        else
            IPH = I
        end if
        IR = RSHPTR(I)
        NR = RSHPTR(I + 1) - IR
        if (NR .gt. NLM) then
            write (6, *) 'Error [_computeSrc]: NR>NLM 3', I, IR, RSHPTR(I + 1)
            stop
        end if
        if (NSTOKES .eq. 1) then
            call shd_src__calSrcPnt_(NLM, NLEG, LOFJ,&
                &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                &NR, RADIANCE(1, IR + 1), SOURCET)
        else
            call shd_src__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                &SRCTYPE, DIRFLUX(I)*SECMU0, YLMSUN,&
                &PLANCK(I), ALBEDO(I), IPH, LEGEN,&
                &NR, RADIANCE(1, IR + 1), SOURCET)
        end if
        if (FIXSH) then
            NS = SHPTR(I + 1) - IS
        else
            ! Then do the adaptive spherical harmonics:
            !   Find how many source function terms to use, and how many
            !   terms wanted for the radiance.
            if (SRCTYPE .eq. 'S') then
                JS = 0
            else
                JS = 1
            end if
            do J = 1, NLM
                do K = 1, NSTOKES
                    if (abs(SOURCET(K, J)) .gt. SRCMIN) JS = J
                end do
            end do
            ! Use all the m's for the last l for the source function.
            if (JS .eq. 0) then
                NS = 0
            else
                LS = LOFJ(JS)
                if (LS .le. MM) then
                    NS = (LS*(LS + 1)) + LS + 1
                else
                    NS = (2*MM + 1)*LS - (MM*(1 + (MM - 1))) + MM + 1
                end if
            end if
            SHPTR(I) = IS
        end if
        if (IS + NS .gt. MAXIV) then
            write (6, *) 'Error [_computeSrc]: MAXIV exceeded.', MAXIV
            write (6, *) '  Out of memory for more spherical harmonic terms.'
            stop
        end if
        ! Transfer the new source function
        do J = 1, NS
            SOURCE(:, IS + J) = SOURCET(:, J)
        end do
        IS = IS + NS
    end do
    SHPTR(NPTS + 1) = IS
    deallocate (LOFJ, SOURCET)
    return
end subroutine shd_src__computeSrc

subroutine shd_src__calSrcPnt(&
    &NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    &SOURCET)

    ! Calculates the polarized source function in real generalized
    ! spherical harmonics at one grid point from the radiance in
    ! real generalized spherical harmonics and thermal and/or solar
    ! source information.

    integer, intent(in) :: NSTOKES, NSTLEG, NLEG, NLM, LOFJ(NLM), IPH, NR
    real, intent(in) :: FLUX0, YLMSUN(NSTLEG, NLM), PLANCK, ALBEDO
    real, intent(in) :: LEGEN(NSTLEG, 0:NLEG, *)
    real, intent(in) :: RADIANCE(NSTOKES, NR)
    character(len=1), intent(in) :: SRCTYPE

    real, intent(out) :: SOURCET(NSTOKES, NLM)

    real, parameter :: C = 3.544907703
    integer J

    SOURCET(:, :) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.  Do different parts of the
    ! matrix multiply depending on NSTOKES.
    do J = 1, NR
        SOURCET(1, J) = SOURCET(1, J) + LEGEN(1, LOFJ(J), IPH)*RADIANCE(1, J)
    end do

    if (NSTOKES .gt. 1) then

        do J = 1, NR
            SOURCET(1, J) = SOURCET(1, J) + LEGEN(5, LOFJ(J), IPH)*RADIANCE(2, J)
        end do

        do J = 5, NR
            SOURCET(2, J) = SOURCET(2, J) + LEGEN(5, LOFJ(J), IPH)*RADIANCE(1, J)&
                &+ LEGEN(2, LOFJ(J), IPH)*RADIANCE(2, J)
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(3, LOFJ(J), IPH)*RADIANCE(3, J)
        end do

    end if

    if (NSTOKES .eq. 4) then

        do J = 1, NR
            SOURCET(3, J) = SOURCET(3, J) + LEGEN(6, LOFJ(J), IPH)*RADIANCE(4, J)
            SOURCET(4, J) = SOURCET(4, J) - LEGEN(6, LOFJ(J), IPH)*RADIANCE(3, J)&
                &+ LEGEN(4, LOFJ(J), IPH)*RADIANCE(4, J)
        end do

    end if

    ! Do the final multiplication by the single scattering albedo
    SOURCET(:, :) = ALBEDO*SOURCET(:, :)

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        do J = 1, NLM
            SOURCET(1, J) = SOURCET(1, J)&
                &+ FLUX0*ALBEDO*LEGEN(1, LOFJ(J), IPH)*YLMSUN(1, J)
        end do

        if (NSTOKES .gt. 1) then
            do J = 5, NLM
                SOURCET(2, J) = SOURCET(2, J)&
                    &+ FLUX0*ALBEDO*LEGEN(5, LOFJ(J), IPH)*YLMSUN(1, J)
            end do
        end if

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1, 1) = SOURCET(1, 1) + C*PLANCK
    end if

    return

end subroutine shd_src__calSrcPnt

subroutine shd_src__calSrcPnt_(&
    &NLM, NLEG, LOFJ,&
    &SRCTYPE, FLUX0, YLMSUN, PLANCK,&
    &ALBEDO, IPH, LEGEN, NR, RADIANCE,&
    &SOURCET)

    ! Calculates the scalar source function in spherical harmonics at
    ! one grid point from the radiance in spherical harmonics and
    ! thermal and/or solar source information.

    integer, intent(in) :: NLEG, NLM, LOFJ(NLM), IPH, NR
    character(len=1), intent(in) :: SRCTYPE
    real, intent(in) :: FLUX0, YLMSUN(NLM), PLANCK, ALBEDO
    real, intent(in) :: LEGEN(0:NLEG, *)
    real, intent(in) :: RADIANCE(NR)

    real, intent(out) :: SOURCET(NLM)

    real, parameter :: C=3.544907703
    integer J

    SOURCET(:) = 0.0

    ! Calculate the generalized spherical harmonics source function
    ! from the solar direct beam
    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then

        do J = 1, NLM
            SOURCET(J) = FLUX0*ALBEDO*LEGEN(LOFJ(J), IPH)*YLMSUN(J)
        end do

    end if

    ! Add in the isotropic unpolarized thermal emission
    if (SRCTYPE .eq. 'T' .or. SRCTYPE .eq. 'B') then
        SOURCET(1) = SOURCET(1) + C*PLANCK
    end if

    ! Calculate the generalized spherical harmonics source function
    ! from the radiance vector in genSH.
    do J = 1, NR
        SOURCET(J) = SOURCET(J) + ALBEDO*LEGEN(LOFJ(J), IPH)*RADIANCE(J)
    end do

    return

end subroutine shd_src__calSrcPnt_

subroutine shd_src__computeSolarStokesVect(&
    &NLEG, SINGLEGEN, COSSCAT,&
    &SOLARMU, MU, DELPHI, NSTOKES,&
    &SCATVECT)

    ! Calculates the Stokes vector for a single scattering of unpolarized
    ! sunlight from the solar direction (SOLARMU) to the outgoing
    ! direction (MU) and with cosine scattering angle (COSSCAT).
    ! The A1 and B1 elements of the phase matrix are calculated from
    ! the Wigner function expansion coefficients (SINGLEGEN).
    ! The plane of polarization is rotated from the scattering plane to
    ! the outgoing meridional plane.  Also divides the output vector by 4*pi.

    integer NLEG, NSTOKES
    real SINGLEGEN(2, 0:NLEG), COSSCAT, SOLARMU, MU, DELPHI
    real SCATVECT(NSTOKES)
    integer I, J, K, L
    double precision OFOURPI, X, FCT, A1, B1
    double precision, allocatable :: DMM1(:)
    double precision SIN_SCAT, SIN_THETA1, SIN_THETA2
    double precision SINPHI, COSPHI, SIN2, COS2, SIN22, COS22
    parameter(OFOURPI=0.07957747151)

    allocate (DMM1(0:NLEG))

    ! Sum the first Wigner function series, which is actually a Legend
    X = dble(COSSCAT)
    call shd_pha__computeWigner(X, NLEG, 0, 0, DMM1)
    A1 = 0.0d0
    do L = 0, NLEG
        FCT = 2.0d0*L + 1.0d0
        A1 = A1 + FCT*dble(SINGLEGEN(1, L))*DMM1(L)
    end do
    SCATVECT(1) = A1*OFOURPI

    if (NSTOKES .gt. 1) then
        ! If doing polarization, sum the second Wigner function series
        call shd_pha__computeWigner(X, NLEG, 2, 0, DMM1)
        B1 = 0.0d0
        do L = 0, NLEG
            FCT = 2.0d0*L + 1.0d0
            B1 = B1 - FCT*dble(SINGLEGEN(2, L))*DMM1(L)
        end do

        ! Rotate the plane of polarization from the scattering plane
        !   to the outgoing meridional plane.
        SIN_SCAT = DSQRT(max(0.0d0, 1.d0 - COSSCAT**2))
        SIN_THETA1 = DSQRT(1.d0 - SOLARMU**2)
        SIN_THETA2 = DSQRT(1.d0 - MU**2)
        SINPHI = DSIN(dble(DELPHI))
        COSPHI = DCOS(dble(DELPHI))
        if (SIN_SCAT .eq. 0.0) then
            SIN2 = 0.0d0
            COS2 = -1.0d0
        else
            SIN2 = SIN_THETA1*SINPHI/SIN_SCAT
            COS2 = (SIN_THETA2*SOLARMU - SIN_THETA1*MU*COSPHI)/SIN_SCAT
        end if
        SIN22 = 2.0d0*SIN2*COS2
        COS22 = 1.0d0 - 2.0d0*SIN2**2
        SCATVECT(2) = B1*COS22*OFOURPI
        SCATVECT(3) = B1*SIN22*OFOURPI
    end if
    if (NSTOKES .eq. 4) then
        SCATVECT(4) = 0.0
    end if
    deallocate (DMM1)

    return

end subroutine shd_src__computeSolarStokesVect

end module shd_src
