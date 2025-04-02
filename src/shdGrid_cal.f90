subroutine shd_comp__interpGrid(&
    &NPTS, NSTLEG, NLEG, GRIDPOS,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)

    ! Calls shd_comp__triInterpProp to interpolate the input arrays from
    !   the property grid to each internal grid point.

    implicit none
    integer NPTS, NSTLEG, NLEG
    integer*2 IPHASE(*)
    real GRIDPOS(3, NPTS)
    real TEMP(*), EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    integer IP

    ! Initialize: transfer the tabulated phase functions
    call shd_comp__triInterpProp(0.0, 0.0, 0.0, .true., NSTLEG, NLEG,&
        &TEMP(1), EXTINCT(1), ALBEDO(1), LEGEN(1, 0, 1), IPHASE(1))

    ! Trilinearly interpolate from the property grid to the adaptive
    do IP = 1, NPTS
        call shd_comp__triInterpProp&
            &(GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &.false., NSTLEG, NLEG, TEMP(IP), EXTINCT(IP),&
            &ALBEDO(IP), LEGEN(1, 0, IP), IPHASE(IP))
    end do

    return
end subroutine shd_comp__interpGrid

subroutine shd_comp__triInterpProp(X, Y, Z, INIT, NSTLEG, NLEG,&
    &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE)
    !  Trilinearly interpolates the quantities on the input property
    ! grid at the single point (X,Y,Z) to get the output TEMP,EXTINCT,
    ! ALBEDO, and LEGEN or IPHASE.  Interpolation is done on the
    ! volume coefficients.  Also adds in the separate gaseous absorption.
    ! Divides the phase function Legendre coefficients LEGEN by 2*l+1.
    ! The phase function pointer IPHASE (for tabulated phase functions)
    ! is that of the maximum weighted scattering property grid point.
    ! If INIT=.TRUE. then transfers the tabulated phase functions.

    use globals_shdom

    implicit none
    integer NSTLEG, NLEG
    integer*2 IPHASE
    logical INIT
    real X, Y, Z, TEMP, EXTINCT, ALBEDO, LEGEN(NSTLEG, 0:NLEG, *)
    integer IX, IXP, IY, IYP, IZ, L, IL, IM, IU, J
    integer I1, I2, I3, I4, I5, I6, I7, I8, I
    double precision U, V, W, F1, F2, F3, F4, F5, F6, F7, F8, F
    double precision SCAT1, SCAT2, SCAT3, SCAT4, SCAT5, SCAT6, SCAT7, SCAT8
    double precision SCATTER, MAXSCAT, KG, EXTMIN, SCATMIN
    save EXTMIN, SCATMIN

    if (INIT) then
       ! If there are tabulated phase functions, then transfer them
        do I = 1, NUMPHASE_
            do L = 0, NLEG
                do J = 1, NSTLEG
                    LEGEN(J, L, I) = LEGENP_(J + NSTLEG*(L + (NLEG + 1)*(I - 1)))/(2*L + 1)
                end do
            end do
        end do
        EXTMIN = 1.0e-5/((ZLEVELS_(NPZ_) - ZLEVELS_(1))/NPZ_)
        SCATMIN = 0.1*EXTMIN
        return
    end if

    ! Find the grid location and compute the interpolation factors
    IL = 0
    IU = NPZ_
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZLEVELS_(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)
    W = dble(Z - ZLEVELS_(IZ))/(ZLEVELS_(IZ + 1) - ZLEVELS_(IZ))
    W = max(min(W, 1.0d0), 0.0d0)
    IX = int((X - XSTART_)/DELX_) + 1
    if (abs(X - XSTART_ - NPX_*DELX_) .lt. 0.01*DELX_) IX = NPX_
    if (IX .lt. 1 .or. IX .gt. NPX_) then
        write (6, *) 'TRILIN: Beyond X domain'
        stop
    end if
    IXP = mod(IX, NPX_) + 1
    U = dble(X - XSTART_ - DELX_*(IX - 1))/DELX_
    U = max(min(U, 1.0d0), 0.0d0)
    if (U .lt. 1.0d-5) U = 0.0d0
    if (U .gt. 1.0d0 - 1.0d-5) U = 1.0d0
    IY = int((Y - YSTART_)/DELY_) + 1
    if (abs(Y - YSTART_ - NPY_*DELY_) .lt. 0.01*DELY_) IY = NPY_
    if (IY .lt. 1 .or. IY .gt. NPY_) then
        write (6, *) 'TRILIN: Beyond Y domain'
        stop
    end if
    IYP = mod(IY, NPY_) + 1
    V = dble(Y - YSTART_ - DELY_*(IY - 1))/DELY_
    V = max(min(V, 1.0d0), 0.0d0)
    if (V .lt. 1.0d-5) V = 0.0d0
    if (V .gt. 1.0d0 - 1.0d-5) V = 1.0d0

    F1 = (1 - U)*(1 - V)*(1 - W)
    F2 = U*(1 - V)*(1 - W)
    F3 = (1 - U)*V*(1 - W)
    F4 = U*V*(1 - W)
    F5 = (1 - U)*(1 - V)*W
    F6 = U*(1 - V)*W
    F7 = (1 - U)*V*W
    F8 = U*V*W
    I1 = IZ + NPZ_*(IY - 1) + NPZ_*NPY_*(IX - 1)
    I2 = IZ + NPZ_*(IY - 1) + NPZ_*NPY_*(IXP - 1)
    I3 = IZ + NPZ_*(IYP - 1) + NPZ_*NPY_*(IX - 1)
    I4 = IZ + NPZ_*(IYP - 1) + NPZ_*NPY_*(IXP - 1)
    I5 = I1 + 1
    I6 = I2 + 1
    I7 = I3 + 1
    I8 = I4 + 1

    ! Trilinearly interpolate the temperature, extinction, scattering
    TEMP = F1*TEMPP_(I1) + F2*TEMPP_(I2) + F3*TEMPP_(I3) + F4*TEMPP_(I4)&
        &+ F5*TEMPP_(I5) + F6*TEMPP_(I6) + F7*TEMPP_(I7) + F8*TEMPP_(I8)
    EXTINCT = F1*EXTINCTP_(I1) + F2*EXTINCTP_(I2) + F3*EXTINCTP_(I3)&
        &+ F4*EXTINCTP_(I4) + F5*EXTINCTP_(I5) + F6*EXTINCTP_(I6)&
        &+ F7*EXTINCTP_(I7) + F8*EXTINCTP_(I8)
    SCAT1 = F1*EXTINCTP_(I1)*ALBEDOP_(I1)
    SCAT2 = F2*EXTINCTP_(I2)*ALBEDOP_(I2)
    SCAT3 = F3*EXTINCTP_(I3)*ALBEDOP_(I3)
    SCAT4 = F4*EXTINCTP_(I4)*ALBEDOP_(I4)
    SCAT5 = F5*EXTINCTP_(I5)*ALBEDOP_(I5)
    SCAT6 = F6*EXTINCTP_(I6)*ALBEDOP_(I6)
    SCAT7 = F7*EXTINCTP_(I7)*ALBEDOP_(I7)
    SCAT8 = F8*EXTINCTP_(I8)*ALBEDOP_(I8)
    SCATTER = SCAT1 + SCAT2 + SCAT3 + SCAT4 + SCAT5 + SCAT6 + SCAT7 + SCAT8
    if (EXTINCT .gt. EXTMIN) then
        ALBEDO = SCATTER/EXTINCT
    else
        ALBEDO = SCATTER/EXTMIN
    end if
    ! For tabulated phase functions pick the one we are on top of
    !   or the one with the most scattering weight.
    if (NUMPHASE_ .gt. 0) then
        MAXSCAT = -1.0
        if (SCAT1 .gt. MAXSCAT .or. abs(F1 - 1) .lt. 0.001) then
            MAXSCAT = SCAT1
            IPHASE = IPHASEP_(I1)
        end if
        if (SCAT2 .gt. MAXSCAT .or. abs(F2 - 1) .lt. 0.001) then
            MAXSCAT = SCAT2
            IPHASE = IPHASEP_(I2)
        end if
        if (SCAT3 .gt. MAXSCAT .or. abs(F3 - 1) .lt. 0.001) then
            MAXSCAT = SCAT3
            IPHASE = IPHASEP_(I3)
        end if
        if (SCAT4 .gt. MAXSCAT .or. abs(F4 - 1) .lt. 0.001) then
            MAXSCAT = SCAT4
            IPHASE = IPHASEP_(I4)
        end if
        if (SCAT5 .gt. MAXSCAT .or. abs(F5 - 1) .lt. 0.001) then
            MAXSCAT = SCAT5
            IPHASE = IPHASEP_(I5)
        end if
        if (SCAT6 .gt. MAXSCAT .or. abs(F6 - 1) .lt. 0.001) then
            MAXSCAT = SCAT6
            IPHASE = IPHASEP_(I6)
        end if
        if (SCAT7 .gt. MAXSCAT .or. abs(F7 - 1) .lt. 0.001) then
            MAXSCAT = SCAT7
            IPHASE = IPHASEP_(I7)
        end if
        if (SCAT8 .gt. MAXSCAT .or. abs(F8 - 1) .lt. 0.001) then
            MAXSCAT = SCAT8
            IPHASE = IPHASEP_(I8)
        end if
    else
        ! Standard property file format: unpolarized
        LEGEN(1, 0, 1) = 1.0
        do L = 1, NLEG
            LEGEN(1, L, 1) = SCAT1*LEGENP_(L + (NLEG + 1)*(I1 - 1))&
                &+ SCAT2*LEGENP_(L + (NLEG + 1)*(I2 - 1))&
                &+ SCAT3*LEGENP_(L + (NLEG + 1)*(I3 - 1))&
                &+ SCAT4*LEGENP_(L + (NLEG + 1)*(I4 - 1))&
                &+ SCAT5*LEGENP_(L + (NLEG + 1)*(I5 - 1))&
                &+ SCAT6*LEGENP_(L + (NLEG + 1)*(I6 - 1))&
                &+ SCAT7*LEGENP_(L + (NLEG + 1)*(I7 - 1))&
                &+ SCAT8*LEGENP_(L + (NLEG + 1)*(I8 - 1))
            if (SCATTER .gt. SCATMIN) then
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATTER
            else
                LEGEN(1, L, 1) = LEGEN(1, L, 1)/SCATMIN
            end if
            LEGEN(1, L, 1) = LEGEN(1, L, 1)/(2*L + 1)
        end do
    end if

    ! Add in the gaseous absorption to extinction and albedo
    if (NZCKD_ .gt. 0) then
        IL = 1
        IU = NZCKD_
        do while (IU - IL .gt. 1)
            IM = (IU + IL)/2
            if (Z .le. ZCKD_(IM)) then
                IL = IM
            else
                IU = IM
            end if
        end do
        I = min(max(IL, 1), NZCKD_ - 1)
        F = (Z - ZCKD_(I))/(ZCKD_(I + 1) - ZCKD_(I))
        F = min(max(F, 0.0d0), 1.0d0)
        KG = (1.0 - F)*GASABS_(I) + F*GASABS_(I + 1)
        if (EXTINCT + KG .gt. 0.0) then
            ALBEDO = ALBEDO*EXTINCT/(EXTINCT + KG)
        else
            ALBEDO = 0.0
        end if
        EXTINCT = EXTINCT + KG
    end if
    return
end subroutine shd_comp__triInterpProp

subroutine shd_grid__splitGrid(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
    &DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
    &NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, NX, XGRID, NY, YGRID,&
    &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
    &TEMP, PLANCK, DIRFLUX, SRCTYPE,&
    &SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    &UNITS, WAVENO, WAVELEN, ADAPTCRIT, ADAPTIND)

    !   Splits the cells that have a cell dividing criterion greater than
    ! CURSPLITACC if DOSPLIT is true.  The current splitting criterion
    ! achieved is returned in SPLITCRIT.  If we are at the end of the
    ! grid point, grid cell, or spherical harmonic arrays (sizes given
    ! by MAXIG, MAXIC, MAXIV, MAXIDO) then the OUTOFMEM flag is returned true.

    use shd_grid
    use shd_ext_mpi

    implicit none
    integer MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX
    integer NPTS, NCELLS, BCFLAG, IPFLAG
    integer NX, NY, NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    real ADAPTIND(*)
    integer*2 CELLFLAGS(*), IPHASE(*)
    logical DELTAM, ACCELFLAG, DOSPLIT, OUTOFMEM
    real CURSPLITACC, SPLITCRIT
    real SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(NSTLEG, NLM)
    real WAVENO(2), WAVELEN
    real GRIDPOS(3, *), XGRID(*), YGRID(*), ADAPTCRIT(*)
    real EXTINCT(*), ALBEDO(*), LEGEN(NSTLEG, 0:NLEG, *)
    real RADIANCE(NSTOKES, *), SOURCE(NSTOKES, *)
    real TEMP(*), PLANCK(*), DIRFLUX(*)
    character(len=1) :: SRCTYPE, UNITS
    integer ICELL, ICELL1, IDIR, I, N, NEWPOINTS(3, 4), OLDNPTS
    integer MAXCELLS, MAXPTS, MAXWORK, MAXSH, MAXRAD
    logical OUTOFMEM0
    real ADAPT(3), FRAC

    if (DOSPLIT) then
        OLDNPTS = NPTS
        ! Do the cells in batchs, so the newly split cells can be tested
        !   for further splitting
        ICELL1 = 1
        do while (ICELL1 .le. NCELLS)
            ! For the current batch of cells, compute the adaptive grid
            !   splitting criterion for 3 directions for each cell.
            N = 0
            do ICELL = ICELL1, NCELLS
                if (TREEPTR(2, ICELL) .eq. 0) then
                    ! If this is a end node cell then compute adaptive criterion
                    call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
                        &SHPTR, NSTOKES, SOURCE, ICELL,&
                        &ADAPT, ADAPTCRIT(N + 1), IDIR)
                    ADAPTIND(N+1) = 4*ICELL + IDIR
                    N = N + 1
                end if
            end do
            ICELL1 = NCELLS + 1

            ! Sort the adaptive criterion in decreasing order
            call shd_util__ssort(ADAPTCRIT, ADAPTIND, N, -2)

            ! Figure how many cells we can split, leaving some room
            !   for cells split during grid smoothing
            FRAC = 0.03
            MAXCELLS = MAXIC - FRAC*(MAXIC - NCELLS) - 2
            MAXPTS = MAXIG - FRAC*(MAXIG - NPTS) - 4
            MAXWORK = MAXIDO - FRAC*(MAXIDO - NPHI0MAX*NPTS) - 4*NPHI0MAX
            MAXSH = MAXIV - FRAC*(MAXIV - SHPTR(NPTS + 1)) - 4*NLM
            MAXRAD = MAXIV + MAXIG - FRAC*(MAXIV + MAXIG - RSHPTR(NPTS + 1)) - 4*NLM

            ! Go down the list from highest adaptive cell criterion on down
            OUTOFMEM0 = OUTOFMEM
            I = 1
            do while (ADAPTCRIT(I) .gt. CURSPLITACC .and. .not. OUTOFMEM0&
                &.and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    &NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    &SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    &RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM0 = .true.
                else
                    ICELL = int(ADAPTIND(I))/4
                    IDIR = iand(int(ADAPTIND(I)), 3)
                    ! Create the new grid points and make pointers
                    call shd_grid__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                    ! Interpolate the medium properties, radiance, and source
                    call shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
                        &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                        &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                        &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                        &UNITS, WAVENO, WAVELEN,&
                        &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                        &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                end if
                I = I + 1
            end do

            MAXCELLS = MAXIC - 2
            MAXPTS = MAXIG - 4
            MAXWORK = MAXIDO - 4*NPHI0MAX
            MAXSH = MAXIV - 4*NLM
            MAXRAD = MAXIV + MAXIG - 4*NLM

            ! Split some more cells to make the adaptive grid "smoother"
            do while (.not. OUTOFMEM .and. I .le. N)
                if (NCELLS .gt. MAXCELLS .or. NPTS .gt. MAXPTS .or.&
                    &NPTS*NPHI0MAX .gt. MAXWORK .or.&
                    &SHPTR(NPTS + 1) .gt. MAXSH .or.&
                    &RSHPTR(NPTS + 1) .gt. MAXRAD) then
                    OUTOFMEM = .true.
                else
                    ! See if cell should be split to make better grid:
                    !   IDIR=0 if no, IDIR=direction (1,2,3) otherwise
                    ICELL = int(ADAPTIND(I))/4
                    call shd_grid__testGridSmooth(GRIDPTR, GRIDPOS, NEIGHPTR,&
                        &TREEPTR, CELLFLAGS, SHPTR, ICELL, IDIR)
                    ! Set IDIR=0 to stop adaptive grid smoothing
                    ! IDIR = 0
                    if (IDIR .gt. 0) then
                        ! Create the new grid points and make pointers
                        call shd_grid__divideCell(NPTS, NCELLS, GRIDPTR, NEIGHPTR,&
                            &TREEPTR, CELLFLAGS, GRIDPOS, ICELL, IDIR, NEWPOINTS)
                        ! Interpolate the medium properties, radiance, and source
                        call shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
                            &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
                            &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
                            &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                            &UNITS, WAVENO, WAVELEN,&
                            &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
                            &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)
                    end if
                end if
                I = I + 1
            end do

        end do
        if (OUTOFMEM0) then
            OUTOFMEM = .true.
            write (6, '(2A,I8,A,I8)')&
                &'Message [_splitGrid]: Further cell splitting requires larger array allocations: ',&
                &'Ncells = ', NCELLS, ', Npnts = ', NPTS
        end if

        ! If in parallel mode, calculate all the direct beam fluxes at once
        if (SRCTYPE .ne. 'T' .and.&
            &(btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) then
            call shd_ext__mpiMakeDirect(OLDNPTS + 1, NPTS, BCFLAG,&
                &IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
                &SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
                &NX, XGRID, NY, YGRID, DIRFLUX)
        end if
    end if

    ! Find the max adaptive cell criterion after the cell divisions
    SPLITCRIT = 0.0
    do ICELL = 1, NCELLS
        if (TREEPTR(2, ICELL) .eq. 0) then
            call shd_grid__testCellSplit(GRIDPTR, GRIDPOS, EXTINCT,&
                &SHPTR, NSTOKES, SOURCE, ICELL,&
                &ADAPT, ADAPTCRIT(1), IDIR)
            SPLITCRIT = max(SPLITCRIT, ADAPTCRIT(1))
        end if
    end do

    return

end subroutine shd_grid__splitGrid

subroutine shd_grid__interpNewGridPnt(NEWPOINTS, NPTS,&
    &NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG, GRIDPOS,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
    &UNITS, WAVENO, WAVELEN,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP, PLANCK, DIRFLUX,&
    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR)

    !   Interpolates the medium properties, radiance, and source function for
    ! new grid points (specified in NEWPOINTS).  The medium properties are
    ! interpolated from the property grid, since interpolating from the
    ! internal grid may not faithfully represent the true medium.
    ! The delta-M scaling is done to the new grid point properties if needed.
    ! The direct beam flux and the Planck function for the new points is
    ! recomputed (rather than interpolated).  The radiance is averaged
    ! from the two parent grid points, and used to compute the source
    ! function for the new points.

    integer NEWPOINTS(3, 4), NPTS
    integer NSTOKES, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer BCFLAG, IPFLAG
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    integer*2 IPHASE(NPTS)
    logical DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN(NSTLEG, NLM)
    real WAVENO(2), WAVELEN
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real TEMP(NPTS), PLANCK(NPTS), DIRFLUX(NPTS)
    real RADIANCE(NSTOKES, *), SOURCE(NSTOKES, *)
    real GRIDPOS(3, *)
    character(len=1) :: SRCTYPE, UNITS
    integer I, IP1, IP2, IP, IPH, J, K, L, M, ME, MS
    integer IS, NS, NS1, NS2
    integer IR, IR1, IR2, NR, NR1, NR2, SIDE
    logical VALIDBEAM
    real X, Y, Z, F, BB, C, D, SECMU0
    real RAD1(NSTOKES), RAD2(NSTOKES)
    real UNIFZLEV, XO, YO, ZO, DIRPATH
    integer, allocatable :: LOFJ(:)
    real, allocatable :: SOURCET(:, :)

    allocate (LOFJ(NLM), SOURCET(NSTOKES, NLM))
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

    C = sqrt(4.0*acos(-1.0))
    if (SRCTYPE .ne. 'T') then
        SECMU0 = 1.0/abs(SOLARMU)
    end if

    ! Loop over the (potentially) four new grid points
    do I = 1, 4
        if (NEWPOINTS(3, I) .gt. 0) then
            ! Get the parent and new grid pointers
            IP1 = NEWPOINTS(1, I)
            IP2 = NEWPOINTS(2, I)
            IP = NEWPOINTS(3, I)

            ! Interpolate the medium properties from the property grid
            X = GRIDPOS(1, IP)
            Y = GRIDPOS(2, IP)
            Z = GRIDPOS(3, IP)
            call shd_comp__triInterpProp(X, Y, Z, .false., NSTLEG, NLEG,&
                &TEMP(IP), EXTINCT(IP), ALBEDO(IP),&
                &LEGEN(1, 0, IP), IPHASE(IP))
            ! Do the Delta-M scaling of extinction and albedo for this point
            if (DELTAM) then
                if (NUMPHASE .gt. 0) then
                    F = LEGEN(1, ML + 1, IPHASE(IP))
                else
                    F = LEGEN(1, ML + 1, IP)
                    do L = 0, ML
                        LEGEN(1, L, IP) = (LEGEN(1, L, IP) - F)/(1 - F)
                    end do
                end if
                EXTINCT(IP) = (1.0 - ALBEDO(IP)*F)*EXTINCT(IP)
                ALBEDO(IP) = (1.0 - F)*ALBEDO(IP)/(1.0 - ALBEDO(IP)*F)
            end if

            ! Compute the new Planck source function (if needed)
            if (SRCTYPE .ne. 'S') then
                call shd_util__computePlanck(TEMP(IP), UNITS, WAVENO, WAVELEN, BB)
                PLANCK(IP) = (1.0 - ALBEDO(IP))*BB
            end if

            ! Compute the new direct beam flux (if needed)
            if (SRCTYPE .ne. 'T') then
                if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) then
                    ! If parallel mode then calculate a temporary approximation
                    !   to the direct flux (for source function calculation below)
                    DIRFLUX(IP) = exp(0.5*(log(DIRFLUX(IP1)) + log(DIRFLUX(IP2))))
                else
                    ! Otherwise, calculate the exact direct beam from property grid
                    DIRPATH = 0.0
                    call shd_comp__computeDirectBeam(0, X, Y, Z, BCFLAG, IPFLAG, DELTAM,&
                        &ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ,&
                        &DIRFLUX(IP),&
                        &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                end if
            end if

            ! Interpolate the radiance
            IR1 = RSHPTR(IP1)
            IR2 = RSHPTR(IP2)
            NR1 = RSHPTR(IP1 + 1) - RSHPTR(IP1)
            NR2 = RSHPTR(IP2 + 1) - RSHPTR(IP2)
            NR = max(NR1, NR2)
            IR = RSHPTR(IP)
            RSHPTR(IP + 1) = IR + NR
            do J = 1, NR
                if (J .le. NR1) then
                    RAD1(:) = RADIANCE(:, IR1 + J)
                else
                    RAD1(:) = 0.0
                end if
                if (J .le. NR2) then
                    RAD2(:) = RADIANCE(:, IR2 + J)
                else
                    RAD2(:) = 0.0
                end if
                RADIANCE(:, IR + J) = 0.5*(RAD1(:) + RAD2(:))
            end do

            ! Compute the source function for the new points
            NS1 = SHPTR(IP1 + 1) - SHPTR(IP1)
            NS2 = SHPTR(IP2 + 1) - SHPTR(IP2)
            NS = max(NS1, NS2)
            IS = SHPTR(IP)
            SHPTR(IP + 1) = IS + NS
            if (ACCELFLAG) OSHPTR(IP + 1) = OSHPTR(IP)
            if (NUMPHASE .gt. 0) then
                IPH = IPHASE(IP)
            else
                IPH = IP
            end if
            if (NSTOKES .eq. 1) then
                call shd_solv__calSrcPnt_(NLM, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    &PLANCK(IP), ALBEDO(IP), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            else
                call shd_solv__calSrcPnt(NSTOKES, NLM, NSTLEG, NLEG, LOFJ,&
                    &SRCTYPE, DIRFLUX(IP)*SECMU0, YLMSUN,&
                    &PLANCK(IP), ALBEDO(IP), IPH, LEGEN,&
                    &NR, RADIANCE(1, IR + 1), SOURCET)
            end if
            do J = 1, NS
                SOURCE(:, IS + J) = SOURCET(:, J)
            end do
        end if
    end do
    deallocate (LOFJ, SOURCET)

    return

end subroutine shd_grid__interpNewGridPnt
