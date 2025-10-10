!╭────────────────────────────────────────────────────────────────────────────╮!
! @LICENSE_HEADER_START@
!
!   This file is part of SHDOMNG.
!
!   --
!   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
!   --
!
!   Copyright (C) 2025
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
!╰────────────────────────────────────────────────────────────────────────────╯!

!╭────────────────────────────────────────────────────────────────────────────╮!
! Additional notes:
!
! from original SHDOM:
! SHDOMPP: Spherical harmonic discrete ordinate radiative transfer method
!          for plane-parallel atmospheres.
!     See shdompp.doc for documentation.
!╰────────────────────────────────────────────────────────────────────────────╯!

program SHDOMPP

    use globals_const
    use globals_cfg
    use globals_prp
    use shd_init
    use shd_inp ! shd_io__
    use shd_out ! shd_io__

    implicit none

    !  MAXDELTAU  maximum scattering optical depth allowed in initial sublayers
    !  SPLITTING_FACTOR  max number of internal grid points is this factor
    !                    times the number of input layers + optical depth

    real, parameter :: MAXDELTAU = 1.0
    real, parameter :: SPLITTING_FACTOR = 10.0
    logical, parameter :: PRNT = .true. ! True for iteration printing

    integer, parameter :: MAXPAR = 100 ! Max array size for user input arrays
    integer, parameter :: MAXOUT = 20
    real :: OUTPARMS(MAXPAR, MAXOUT)
    character(len=1) :: OUTTYPES(MAXOUT)
    character(len=FNAMEMAXLEN_) :: OUTFILES(MAXOUT)

    real, allocatable :: SFCPARMS(:), SFCPARMS_TMP(:)

    integer :: NLAY, ILAY, IPHASE
    integer :: NZOUT, NMUOUT, NPHIOUT, FLUXPARM
    integer :: NCASES, ICASE
    real    :: SFCTEMP

    integer :: MAXPG, MAXIG, MAXIC, MAXIV, MAXIDO, MAXNZ, MAXNG
    integer :: MAXLEG, MAXPGL, MAXIGL
    integer :: MAXNBC, MAXBCRAD
    integer :: MAXSFCPTS, MAXSFCPARS, BIG_ARRAYS

    real :: MEMWORD, WANTMEM
    real :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT
    real, allocatable :: SOLFLUX(:)
    real, allocatable :: XGRID(:), YGRID(:), ZGRID(:), GRIDPOS(:, :)

    integer, allocatable :: GRIDPTR(:, :), NEIGHPTR(:, :), TREEPTR(:, :)
    integer*2, allocatable:: CELLFLAGS(:)
    integer, allocatable :: RSHPTR(:), SHPTR(:), OSHPTR(:)
    real, allocatable :: RADIANCE(:, :), SOURCE(:, :), DELSOURCE(:, :)
    real, allocatable :: FLUXES(:, :), DIRFLUX(:), DIRFLUX_TMP(:)
    integer, allocatable :: NPHI0(:)
    real, allocatable :: MU(:), WTDO(:, :), PHI(:, :)
    integer, allocatable :: BCPTR(:)
    real, allocatable :: BCRAD(:, :)
    real, allocatable :: SENPARMS(:)

    real, allocatable :: SUMFLUXES(:, :), SUMDIRFLUX(:), SUMFLUXDIV(:)
    real, allocatable :: SUMSHTERMS(:), SOURCE1OUT(:)
    real, allocatable :: SUMRADOUT(:), SUMVISOUT(:)
    real, allocatable :: ALLFLUXES(:, :, :, :)
    real, allocatable :: ALLFLUXDIV(:, :, :), ALLSHTERMS(:, :, :, :)
    real, allocatable :: XALLGRID(:), YALLGRID(:)
    real, allocatable :: TEMPPT(:), EXTINCTPT(:), ALBEDOPT(:), LEGENPT(:)
    integer*2, allocatable :: IPHASEPT(:)

    integer :: NXSFC, NYSFC, NSFCPAR, NTOPPTS, NBOTPTS
    integer :: NG, SIDE, IG, ITER, I
    integer :: ML, MM, NLM, NSTLEG, NLEG, NANG, NPHI0MAX, NX1, NY1
    integer :: NBPTS, NPTS, NSH, NCELLS, NBCELLS, OLDNPTS
    integer :: NRAD, IRAD, NXOUT, NYOUT, NVIS, IVIS, NSCAN, NPIX, NSHOUT
    integer :: NPTSTOT, NSHTOT
    logical :: INRADFLAG, NEWGRIDFLAG
    logical :: BASEOUT, HIGHORDERRAD, MASTERPROC, VALIDBEAM
    logical :: RADIANCE_OUT, VISUAL_OUT
    logical :: FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
    real :: DELXSFC, DELYSFC, MAXASYM, REDUCE_RATIO
    real :: UNIFZLEV, XO, YO, ZO, DIRPATH
    character(len=2) :: SFCTYPE
    integer :: NSENDATA, NSENPAR

    character(len=FNAMEMAXLEN_) :: RUNNAME, PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, SENFILE, OutFileNC
    integer :: NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    logical :: KDIST, DELTAM, ACCELFLAG
    character(len=1) :: GRIDTYPE, SRCTYPE, UNITS
    character(len=1) :: PROPTYPE
    real :: SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD
    real :: GNDTEMP, GNDALBEDO, WAVENO(2), WAVELEN
    real :: SOLACC, SOLCRIT, SPLITACC, SHACC
    integer :: TOTITER, MAXITER, NUMOUT
    real :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

    real :: SOLARFLUX

    real, allocatable :: ZOUT(:), MUOUT(:), PHIOUT(:)
    integer, allocatable :: NLEGP(:), IXP(:, :)
    real, allocatable :: TAUP(:), ALBEDOP(:), LEGENP(:, :)
    real, allocatable :: TAUPM(:), ALBEDOPM(:)
    real, allocatable :: HEIGHTP(:), TEMPP(:)
    real, allocatable :: DELG(:), ZCKD(:), KABS(:, :)
    real, allocatable :: HEIGHTG(:), TAUG(:)
    real, allocatable :: FLUXUP(:), FLUXDN(:), MEANRAD(:), RADOUT(:, :, :)
    real, allocatable :: FLUXUPM(:), FLUXDNM(:), MEANRADM(:), RADOUTM(:, :, :)

    integer :: count_start, count_rate, count_max
    integer :: count_end
    real :: elapsed_time
    real :: DELX0, DELY0
    logical :: mpiX, mpiY, mpiMode
    logical :: nmlMode = .true.

    ! namelist input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call shd_io__nmlInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
        INSAVEFILE, OUTSAVEFILE, SENFILE,&
        NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
        IPFLAG, KDIST, DELTAM, GRIDTYPE,&
        SRCTYPE, SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD,&
        GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
        ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
        MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
        OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! process input parameters
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! set calculation accuracy (1.0e-4 to 1.0e-8)
    RACC_ = max(REPS_, min(SOLACC*0.1, 1.0e-4))

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    NSTOKES = max(1, min(4, NSTOKES))
    if (NSTOKES .eq. 2) NSTOKES = 3

    ! X, Y, Z grids (force NX&NY to 1 and IPFLAG to 3)
    NX = 1; NY = 1; NZ = max(2, NZ)
    IPFLAG = 3

    ! Make ML and MM from NMU and NPHI
    !   NMU : number of discrete ordinate zenith angles in both hemispheres
    !   NPHI: number of discrete ordinate azimuth angles in 2*pi
    !   ML  : The maximum L of the spherical harmonic truncation
    !   MM  : The maximum M of the spherical harmonic truncation
    NMU  = max(2, 2*int((NMU + 1)/2)) ! force NMU an even number
    NPHI = max(1, NPHI)

    ! Make NPHI0MAX, ML, MM, and NLM from NMU and NPHI
    NPHI0MAX = int((NPHI + 2)/2)
    ML = NMU - 1
    MM = max(0, int(NPHI/2) - 1)
    NLM = (MM + 1)*(ML + 1) - (MM*(MM + 1))/2

    ! NLEG: The number of Legendre series phase function terms used in the solution procedure
    NLEG = ML
    if (DELTAM) NLEG = ML + 1

    ! NSTLEG is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (NSTOKES .eq. 1) then
        NSTLEG = 1
    else
        NSTLEG = 6
    end if

    NLEG = ML
    if (DELTAM) NLEG = ML + 1

    ! dummy variables for output
    NX1 = 1; NY1 = 1
    NBPTS = NX1*NY1*NZ


    allocate (MU(NMU), PHI(NPHI0MAX, NMU), NPHI0(NMU), WTDO(NMU, NPHI0MAX))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! read PROPFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    call shd_io__readPrpSize(PROPFILE, NSTLEG, NLEG, NPX_, NPY_, NPZ_, DELX_, DELY_,&
        &NUMPHASE_, MAXLEG, MAXPGL)

    MAXPG = NPZ_
    if (NUMPHASE_ .eq. 0) MAXPGL = MAXPG*(MAXLEG + 1)*NSTLEG
    MAXNZ = NPZ_

    allocate (ZLEVELS_(NPZ_), TEMPPT(MAXPG), EXTINCTPT(MAXPG), ALBEDOPT(MAXPG), LEGENPT(MAXPGL), IPHASEPT(MAXPG))

    call shd_io__readPrp(PROPFILE, NPX_, NPY_, NPZ_,&
        &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,&
        &PROPTYPE, DELX_, DELY_, ZLEVELS_, MAXASYM,&
        &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE_, IPHASEPT)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! process layer property
    ! reverse order (first value at top)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    NLAY = NPZ_ - 1
    allocate (HEIGHTP(NLAY+1), TEMPP(NLAY+1))
    HEIGHTP(1:NPZ_) = ZLEVELS_(NPZ_:1:-1)
    TEMPP(1:NPZ_) = TEMPPT(NPZ_:1:-1)

    allocate (TAUP(NLAY), ALBEDOP(NLAY))
    ! needs further investigation on the layer averaging issue
    TAUP(1:NLAY) = EXTINCTPT(NLAY+1:2:-1)*(HEIGHTP(1:NLAY)-HEIGHTP(2:NPZ_))
    ALBEDOP(1:NLAY) = ALBEDOPT(NLAY:1:-1)

    allocate (NLEGP(NLAY), LEGENP(0:MAXLEG, NLAY))
    do ILAY = 1, NLAY
        IPHASE = IPHASEPT(NLAY-ILAY+1)
        NLEGP(ILAY) = NLEG
        LEGENP(0:MAXLEG, ILAY) = LEGENPT(1+(MAXLEG+1)*(IPHASE-1):(MAXLEG+1)*IPHASE)
    end do

    FLUX_OUT = .false.
    RADIANCE_OUT = .false.
    NRAD = 0
    do I = 1, NUMOUT
        if (OUTTYPES(I) .eq. 'R') then
            RADIANCE_OUT = .true.
            NXOUT = 1
            NYOUT = 1
            NRAD = NRAD + NSTOKES*NXOUT*NYOUT*nint(OUTPARMS(6, I))
            IRAD = I

            NZOUT = 1
            allocate (ZOUT(NZOUT))
            ZOUT(1) = OUTPARMS(1, I)

        else if (OUTTYPES(I) .eq. 'F') then
            if (.not. FLUX_OUT) then
                allocate (SUMFLUXES(2, MAXIG), SUMDIRFLUX(MAXIG))
            end if
            FLUX_OUT = .true.
        end if
    end do

    allocate (SUMRADOUT(NRAD))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! allocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    allocate (TAUPM(NLAY), ALBEDOPM(NLAY))
    allocate (RADOUT(NZOUT, NMUOUT, NPHIOUT), RADOUTM(NZOUT, NMUOUT, NPHIOUT))

    ! Guess the maximum size of the adaptive grid and allocate arrays
    MAXIG = nint(SPLITTING_FACTOR*(NLAY + sum(TAUP(:)*ALBEDOP(:))/MAXDELTAU))
    allocate (IXP(2, NLAY), HEIGHTG(MAXIG), TAUG(MAXIG))
    allocate (FLUXUP(MAXIG), FLUXDN(MAXIG), MEANRAD(MAXIG))
    allocate (FLUXUPM(MAXIG), FLUXDNM(MAXIG), MEANRADM(MAXIG))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! read SFCFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (SFCFILE(1:4) .eq. 'NONE') then
        SFCTYPE = 'FL'
        SFCTEMP = GNDTEMP
        allocate (SFCPARMS(1))
        SFCPARMS(1) = GNDALBEDO
    else
        call shd_io__readSfcSize(SFCFILE, NXSFC, NYSFC, DELXSFC, DELYSFC, MAXSFCPTS, MAXSFCPARS)
        allocate (SFCPARMS_TMP(MAXSFCPARS*MAXSFCPTS), SFCPARMS(MAXSFCPARS-1))
        call shd_io__readSfc(SFCFILE, MAXSFCPTS, MAXSFCPARS,&
            &SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,&
            &NSFCPAR, SFCPARMS_TMP, GNDTEMP, GNDALBEDO)
        SFCTEMP = SFCPARMS_TMP(1)
        SFCPARMS(1:MAXSFCPARS-1) = SFCPARMS_TMP(2:MAXSFCPARS)
        deallocate (SFCPARMS_TMP)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! read CKDFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (KDIST) then
        call shd_io__readCkdSize(CKDFILE, WAVENO, NG, NZCKD_)

        allocate (ZCKD_(NZCKD_), GASABS_(NZCKD_), SOLFLUX(NG), DELG(NG), KABS(NZCKD_, NG))

        MAXNG = NG
        MAXNZ = NZCKD_

        call shd_io__readCkd(MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,&
            &NG, DELG, NZCKD_, ZCKD_, KABS)
    else
        NG = 1
        allocate (SOLFLUX(NG), DELG(NG))
        SOLFLUX(1) = 1.0
        DELG(1) = 1.0
        NZCKD_ = 0
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! read SENFILE
    !╭────────────────────────────────────────────────────────────────────────────╮!
    allocate (MUOUT(NRAD), PHIOUT(NRAD))
    if (SENFILE(1:4) .ne. 'NONE') then
        call shd_io__readSenSize(SENFILE, NSENDATA, NSENPAR)
        allocate (SENPARMS(NSENDATA*NSENPAR))
        call shd_io__readSen(SENFILE, NSENDATA, NSENPAR, SENPARMS)
        MUOUT(:) = SENPARMS(1:2*NRAD-1:2)
        PHIOUT(:) = SENPARMS(2:2*NRAD:2)
    else
        MUOUT(:) = OUTPARMS(7:2*NRAD+5:2, IRAD)
        PHIOUT(:) = OUTPARMS(8:2*NRAD+6:2, IRAD)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! grids
    !╭────────────────────────────────────────────────────────────────────────────╮!
    XSTART_ = 0.0
    YSTART_ = 0.0
    allocate (XGRID(1), YGRID(1), ZGRID(NZ), GRIDPOS(3, NZ))
    call shd_init__makeGrid(BCFLAG, GRIDTYPE,&
        &NPX_, NPY_, NPZ_, NX, NY, NZ,&
        &XSTART_, YSTART_, DELX_, DELY_, ZLEVELS_,&
        &XGRID, YGRID, ZGRID)

    NPTS = NZ
    NCELLS = NZ
    allocate (GRIDPTR(8, NZ), NEIGHPTR(6, NZ), TREEPTR(2, NZ), CELLFLAGS(NZ))
    call shd_init__makeCell(BCFLAG, IPFLAG,&
        &NX, NY, NZ, NX1, NY1, NPTS, NCELLS,&
        &XGRID, YGRID, ZGRID, GRIDPOS,&
        &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
    NBCELLS = NCELLS
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! picard iteration
    !╭────────────────────────────────────────────────────────────────────────────╮!
    TOTITER = 0

    ! Loop over the k-distribution g's from low to high absorption
    do IG = 1, NG
        SOLARFLUX = SOLFLUX(IG)*SOLARFLUX0*abs(SOLARMU)
        write (6, '(A,I2,A,I2,A)') 'Message [shdompp]: Starting ', IG, '/', NG, ' (ig/Ng) ...'

        ! Add in the layer gaseous absorption if doing a k-distribution
        call ADD_IN_KDIST(NLAY, HEIGHTP, TAUP, ALBEDOP, &
                          NZCKD_, ZCKD_, KABS(:, IG), TAUPM, ALBEDOPM)

        !call SOLVE_SHDOMPP(NSTOKES, NLAY, NANG,&
        !                   ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE_,&
        !                   NPHI0MAX, NPHI0, MU, PHI, WTDO,&
        !                   MAXLEG, HEIGHTP, TEMPP, TAUPM, ALBEDOPM, &
        !                   NLEGP, LEGENP, SRCTYPE, SOLARFLUX, SOLARMU, &
        !                   UNITS, WAVELEN, WAVENO, SKYRAD, SFCTYPE(2:2), SFCPARMS, SFCTEMP, &
        !                   MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
        !                   MAXITER, ACCELFLAG, PRNT, SOLCRIT, ITER, &
        !                   NPTS, HEIGHTG, TAUG, IXP, FLUXUPM, FLUXDNM, MEANRADM, &
        !                   NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
        !                   TAUPSC, ALBEDOSC, LEGENSC, TAUGSC, &
        !                   SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
        !stop

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in SOURCE (SHPTR), RADIANCE (RSHPTR), FLUXES, DIRFLUX.
        !call SOLVE_SHDOMPP(NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
        !    &ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE_,&
        !    &NPHI0, MU, PHI, WTDO,&
        !    &MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG,&
        !    &BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD,&
        !    &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
        !    &SFCTYPE, GNDTEMP, GNDALBEDO,&
        !    &NXSFC, NYSFC, DELXSFC, DELYSFC,&
        !    &NSFCPAR, SFCPARMS, SFCGRIDPARMS,&
        !    &UNITS, WAVENO, WAVELEN,&
        !    &ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,&
        !    &SPLITACC, SHACC, XGRID, YGRID, ZGRID,&
        !    &TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,&
        !    &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
        !    &NPTS, GRIDPOS,&
        !    &NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
        !    &RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,&
        !    &SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
        !INRADFLAG = .false.
        !TOTITER = TOTITER + ITER

        TOTITER = TOTITER + ITER

        ! post processing
        !╭────────────────────────────────────────────────────────────────────────────╮!
        !! Add the direct solar component to the downwelling flux and mean radiance
        !if (SRCTYPE /= 'T') then
        !    PI = acos(-1.0)
        !    do I = 1, NPTS
        !        DIRFLUX = SOLARFLUX*exp(-TAUGSC(I)/abs(SOLARMU))
        !        FLUXDN(I) = FLUXDN(I) + DIRFLUX
        !        MEANRAD(I) = MEANRAD(I) + DIRFLUX/(4*PI*abs(SOLARMU))
        !    end do
        !end if

        !! Make an unscaled optical depth grid for the output file
        !if (DELTAM) then
        !    call UNDELTA_SCALE(NLAY, TAUPSC, TAUP, IXP, NPTS, TAUGSC, TAUG)
        !else
        !    TAUG(:) = TAUGSC(:)
        !end if
        !call OUTPUT_HEIGHT_GRID(NLAY, TAUP, IXP, NPTS, TAUG, HEIGHTP, HEIGHTG)

        !if (NZOUT*NMUOUT*NPHIOUT > 0) then
        !    ! Compute the radiance output and put in array
        !    call COMPUTE_RADIANCE(NLAY, ALBEDOSC, LEGENSC, &
        !                          MAXLEG, NLEGP, LEGENP, DELTAM, &
        !                          NPTS, HEIGHTG, TAUG, TAUGSC, IXP, IXG, &
        !                          ML, MM, NLM, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
        !                          SRCTYPE, SOLARFLUX, SOLARMU, &
        !                          SKYRAD, UNITS, WAVELEN, WAVENO, SFCTEMP, &
        !                          SFCTYPE, SFCPARMS, SOURCE, RADIANCE, &
        !                          NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
        !                          RADOUT)
        !end if
        !╰────────────────────────────────────────────────────────────────────────────╯!

        !if (KDIST) then
        !    call SUM_KDIST_OUTPUT(IG, DELG(IG), NLAY, TAUP, IXP, NPTS, TAUG, &
        !                          HEIGHTP, HEIGHTG, FLUXUPM, FLUXUP, &
        !                          FLUXDNM, FLUXDN, MEANRADM, MEANRAD, &
        !                          NZOUT, NMUOUT, NPHIOUT, RADOUTM, RADOUT)
        !else
        !    FLUXUP(1:NPTS) = FLUXUPM(1:NPTS)
        !    FLUXDN(1:NPTS) = FLUXDNM(1:NPTS)
        !    MEANRAD(1:NPTS) = MEANRADM(1:NPTS)
        !    RADOUT(:, :, :) = RADOUTM(:, :, :)
        !end if

    end do
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! output data
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !call OUTPUT_RESULTS(NLAY, NPTS, NMU, NPHI, NG, &
    !                    PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
    !                    SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, &
    !                    SPLITACC, SOLACC, MAXITER, TOTITER, &
    !                    IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
    !                    NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
    !                    'F', FLUXFILE)

    !if (NZOUT*NMUOUT*NPHIOUT > 0) then
    !    call OUTPUT_RESULTS(NLAY, NPTS, NMU, NPHI, NG, &
    !                        PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
    !                        SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, &
    !                        SPLITACC, SOLACC, MAXITER, TOTITER, &
    !                        IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
    !                        NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, RADOUT, &
    !                        'R', RADFILE)
    !end if

    do I = 1, NUMOUT
        if ((OUTTYPES(I) .eq. 'F') .or. (OUTTYPES(I) .eq. 'R')) then
            if (SENFILE(1:4) .ne. 'NONE') then
                NSH = 0
                !print *, SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV
                !print *, IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT
                print *, 'haha'

                !call shd_io__outputResult(NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NPTS,&
                !    NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_, NMU, NPHI, NANG, NG,&
                !    PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                !    BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                !    SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                !    SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                !    SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                !    XGRID, YGRID, ZGRID, GRIDPOS,&
                !    TREEPTR, GRIDPTR, CELLFLAGS,&
                !    EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT, TEMPPT,&
                !    SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,&
                !    IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT,&
                !    OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I), SENPARMS)
            !else
            !    call shd_io__outputResult(NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NPTS,&
            !        NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_, NMU, NPHI, NANG, NG,&
            !        PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
            !        BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
            !        SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
            !        SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
            !        SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
            !        XGRID, YGRID, ZGRID, GRIDPOS,&
            !        TREEPTR, GRIDPTR, CELLFLAGS,&
            !        EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT, TEMPPT,&
            !        SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,&
            !        IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT,&
            !        OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I))
            end if
        else
            write (*, *) "Warning [shdompp]: shdompp does not support output type - ", OUTTYPES(I), ", skipping output ..."
        end if
    end do
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! deallocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !deallocate (MU, WTDO, PHI, NPHI0)
    !deallocate (TAUPSC, ALBEDOSC, LEGENSC, TAUGSC)
    !deallocate (SOURCE, RADIANCE, CONSTSOURCE, DELSOURCE, SRCINT)

    ! layer property
    deallocate (HEIGHTP, TEMPP, TAUP, ALBEDOP, NLEGP, LEGENP, TAUPM, ALBEDOPM)

    ! grid property
    deallocate (IXP, HEIGHTG, TAUG, FLUXUP, FLUXDN, MEANRAD, RADOUT)

    ! mean property
    deallocate (FLUXUPM, FLUXDNM, MEANRADM, RADOUTM)

    ! grids
    deallocate (XGRID, YGRID, ZGRID, GRIDPOS)
    deallocate (GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)

    ! CKD
    deallocate (SOLFLUX, DELG)
    if (KDIST) deallocate (ZCKD_, GASABS_, KABS)

    ! SFC
    if (allocated(SFCPARMS)) deallocate (SFCPARMS)

    ! SEN
    if (allocated(SENPARMS)) deallocate (SENPARMS)

    ! PRP
    deallocate (ZLEVELS_, TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT)
    !╰────────────────────────────────────────────────────────────────────────────╯!

end program SHDOMPP

!subroutine SOLVE_SHDOMPP(NSTOKES, NLAY, NANG,&
!                         ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE_,&
!                         NPHI0MAX, NPHI0, MU, PHI, WTDO,&
!                         MAXLEG, HEIGHTP, TEMPP, TAUP, ALBEDOP, &
!                         NLEGP, LEGENP, SRCTYPE, SOLARFLUX, SOLARMU, &
!                         UNITS, WAVELEN, WAVENO, SKYRAD, SFCTYPE, SFCPARMS, SFCTEMP, &
!                         MAXIG, MAXDELTAU, SPLITACC, SOLACC, &
!                         MAXITER, ACCELFLAG, PRNT, SOLCRIT, ITER, &
!                         NPTS, HEIGHTG, TAUG, IXP, FLUXUP, FLUXDN, MEANRAD, &
!                         NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
!                         TAUPSC, ALBEDOSC, LEGENSC, TAUGSC, &
!                         SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
!    ! Performs the SHDOMPP solution procedure to calculate plane-parallel
!    ! monochromatic, unpolarized radiative transfer.  The optical properties,
!    ! uniform in each layer, are input.  The output are the upwelling and
!    ! downwelling fluxes and mean radiance at the adaptive grid optical depth
!    ! levels (TAUG), and the radiances in the NMUOUT*NPHIOUT specified
!    ! directions at each of NZOUT radiance output levels.
!    !
!    ! Inputs:
!    !  NLAY      number of input property layers (TAUP, ALBEDOP, etc.)
!    !  MAXLEG    maximum Legendre series order for phase functions
!    !  HEIGHTP   heights (e.g. km) of input layer interfaces (NLAY+1 values)
!    !  TEMPP     temperatures (K) of input layer interfaces (NLAY+1 values)
!    !              temperature is linearly interpolated in tau within each layer
!    !  TAUP      optical depths of the input layers (each layer from top down)
!    !  ALBEDOP   single scattering albedo of the inputs layers
!    !  NLEGP     order of each layer Legendre series phase function
!    !  LEGENP    phase function Legendre series coefficients for each layer
!    !  SRCTYPE   'S' for solar source, 'T' for thermal source, 'B' for both
!    !  SOLARFLUX solar flux on a *horizontal* surface
!    !  SOLARMU   negative of cosine of solar zenith angle (< 0)
!    !              direct beam is assumed to be going in phi=0 azimuth
!    !  UNITS     units for radiation output
!    !  WAVELEN   wavelength (microns) for thermal emission
!    !  WAVENO(2) starting and ending wavenumber (cm^-1) for thermal (if WAVELEN=0)
!    !  SKYRAD    top of atmosphere downward sky radiance (or temp for SRCTYPE='T')
!    !  SFCTYPE   surface type (currently 'L' for Lambertian, 'R' for RPV
!    !  SFCPARMS  array of surface parameters (first is albedo for Lambertian)
!    !  SFCTEMP   surface temperature (K)
!    !  NMU       number of discrete ordinate zenith angles in both hemispheres
!    !  NPHI      number of discrete ordinate azimuth angles in 2\pi
!    !  MAXIG     maximum possible size of adaptive grid arrays
!    !  MAXDELTAU max sublayer scattering optical depth in adaptive grid
!    !  SPLITACC  adaptive layer splitting accuracy
!    !  SOLACC    desired solution accuracy (try 1E-4 or less)
!    !  MAXITER   maximum number of iterations allowed
!    !  ACCELFLAG logical flag: true to perform iteration acceleration
!    !  PRNT      logical flat: true to printer iteration information
!    !  NZOUT     number of radiance output height levels
!    !  ZOUT      height levels of radiance output (set to closest adaptive level)
!    !  NMUOUT    number of cosine zenith angles for radiance output
!    !  MUOUT     cosine zenith angles for radiance output
!    !  NPHIOUT   number of azimuth angles for radiance output
!    !  PHIOUT    azimuth angles (degrees) for radiance output
!    !
!    ! Outputs:
!    !  SOLCRIT   actual solution accuracy achieved
!    !  ITER      actual number of iterations
!    !  NPTS      number of adaptive grid points
!    !  HEIGHTG   heights of the adaptive grid levels (same units as HEIGHTP)
!    !  TAUG      sorted grid point optical depth array for the adaptive grid
!    !  IXP       index of starting and ending point in TAUG for each input layer
!    !  FLUXUP    upward hemispheric flux at each adaptive grid level
!    !  FLUXDN    downward hemispheric fluxes (includes direct solar) at each level
!    !  MEANRAD   mean radiance (includes direct solar) at each level
!    !  RADOUT    array of output radiances
!    !
!    !  TAUG(IXP(1,layer)) is the optical depth at the top of the input layer
!    !  TAUG(IXP(2,layer)) is the optical depth at the bottom of the input layer
!    !  FLUXUP(IXP(1,layer)) is the upwelling flux at the top of the layer
!    !  FLUXDN(IXP(2,layer)) is the downwelling flux at the bottom of the layer
!    !

!    use shd_shdo

!    implicit none
!    integer, intent(in) :: NSTOKES, NLAY, MAXLEG, NLEGP(NLAY), MAXIG
!    integer, intent(in) :: ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE_, NPHI0MAX
!    integer, intent(inout) :: NPHI0(*)
!    real, intent(inout) :: MU(*), WTDO(*), PHI(*)

!    integer, intent(in)  :: MAXITER
!    integer, intent(in)  :: NZOUT, NMUOUT, NPHIOUT
!    logical, intent(in)  :: ACCELFLAG, PRNT
!    real, intent(in)  :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
!    real, intent(in)  :: HEIGHTP(NLAY + 1), TEMPP(NLAY + 1)
!    real, intent(in)  :: SOLARFLUX, SOLARMU, WAVELEN, WAVENO(2), SKYRAD
!    real, intent(in)  :: SFCPARMS(*), SFCTEMP
!    real, intent(in)  :: MAXDELTAU, SPLITACC, SOLACC
!    real, intent(inout)  :: ZOUT(NZOUT)
!    character(len=1), intent(in) :: SRCTYPE, SFCTYPE, UNITS
!    integer, intent(out) :: ITER, NPTS, IXP(2, NLAY)
!    real, intent(out) :: SOLCRIT
!    real, intent(out) :: HEIGHTG(MAXIG), TAUG(MAXIG)
!    real, intent(out) :: FLUXUP(MAXIG), FLUXDN(MAXIG), MEANRAD(MAXIG)
!    real, intent(inout) :: RADOUT(*), MUOUT(*), PHIOUT(*)

!    integer :: ORDINATESET, OLDNPTS, PREVNPTS, LASTNPTS, NANG, I
!    logical, parameter :: DELTAM = .true.
!    logical :: FFTFLAG(NMU)
!    real, parameter :: MINDELTAUSPLIT = 1.0e-3
!    real    :: DELJDOT, DELJOLD, DELJNEW
!    real    :: SOLARAZ, PI, GNDALB, DIRFLUX
!    real    :: TRUESOLCRIT, CURSPLITACC
!    integer, allocatable :: IXG(:)
!    real, intent (inout) :: TAUPSC(*), ALBEDOSC(*), LEGENSC(*)
!    real, allocatable :: TAUGSC(:), RADIANCE(:, :, :), SOURCE(:, :)

!    real, allocatable :: WTMU(:), YLMSUN(:)
!    real, allocatable :: CMU1(:, :), CMU2(:, :), CPHI1(:, :, :), CPHI2(:, :, :)
!    real, allocatable :: WPHISAVE(:, :), SRCINT(:, :, :)
!    real, allocatable :: CONSTSOURCE(:, :), SHRADIANCE(:, :), DOSOURCE(:, :, :)
!    real, allocatable :: DELSOURCE(:, :)

!    allocate (TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY))
!    allocate (TAUGSC(MAXIG), IXG(MAXIG))

!    allocate (YLMSUN(NLM), WTMU(NMU), CMU1(NLM, NMU), CMU2(NMU, NLM))
!    allocate (CPHI1(0:16, 32, NMU), CPHI2(32, 0:16, NMU), WPHISAVE(3*NPHI0MAX + 15, NMU))

!    allocate (SOURCE(NLM, MAXIG), RADIANCE(NPHI0MAX, NMU, MAXIG))
!    allocate (CONSTSOURCE(NLM, MAXIG), DELSOURCE(NLM, MAXIG))
!    allocate (SRCINT(0:4, NMU, MAXIG))

!    ! Set up some things before solution loop

!    ! If Delta-M then scale the optical depth, albedo, and Legendre terms.
!    if (DELTAM) then
!        call DELTA_SCALE(NLAY, ML, TAUP, ALBEDOP, MAXLEG, NLEGP, LEGENP, &
!                         TAUPSC, ALBEDOSC, LEGENSC)
!    else
!        TAUPSC(:) = TAUP(:)
!        ALBEDOSC(:) = ALBEDOP(:)
!        LEGENSC(0:ML, :) = LEGENP(0:ML, :)
!    end if

!    call MAKE_TAU_GRID(NLAY, TAUPSC, ALBEDOSC, MAXDELTAU, MINDELTAUSPLIT, &
!                       NPTS, TAUGSC, IXP, IXG)

!    ! Precompute Ylm's for solar direction
!    if (SRCTYPE /= 'T') then
!        SOLARAZ = 0.0
!        call YLMALL(SOLARMU, SOLARAZ, ML, MM, 1, YLMSUN)
!    end if
!    ! Prepare the constant part of the level source function, which are
!    !   the Planck blackbody source and the solar attenuated single scattered
!    !   radiance, as a function of level.
!    CONSTSOURCE(1, 1:NPTS) = -1.0
!    call PREPARE_SOURCE(NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
!                        SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
!                        SOLARFLUX, SOLARMU, YLMSUN, UNITS, WAVELEN, WAVENO, &
!                        CONSTSOURCE)

!    ! Make the discrete ordinates (angles)
!    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
!    ORDINATESET = 3
!    call MAKE_ANGLE_SET(NMU, NPHI, 1, NPHI0MAX, ORDINATESET, &
!                        NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

!    ! Make the Ylm transform coefficients
!    !call MAKE_SH_DO_COEF(ML, MM, NLM, NMU, NPHI0, &
!    !                     NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
!    !                     CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)
!    call shd_shdo__makeSHDO(1, ML, MM, NLM, NMU, NPHI0, &
!                         NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
!                         CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

!    ! Initialize the radiance on the base grid using Eddington
!    !   two-stream plane-parallel
!    if (SFCTYPE == 'L') then
!        GNDALB = SFCPARMS(1)
!    else
!        GNDALB = 0.2
!    end if
!    allocate (SHRADIANCE(NLM, NPTS))
!    call INIT_RADIANCE(NLAY, ML, TAUPSC, ALBEDOSC, LEGENSC, TEMPP, &
!                       SRCTYPE, SOLARFLUX, SOLARMU, GNDALB, SFCTEMP, &
!                       SKYRAD, UNITS, WAVELEN, WAVENO, NPTS, TAUGSC, IXP, NLM, &
!                       SHRADIANCE)

!    ! Initialize the SH source function from the SH radiance field
!    OLDNPTS = 0
!    SOURCE(:, :) = 0.0
!    call COMPUTE_SOURCE(ML, MM, NLM, NLAY, NPTS, IXP, IXG, &
!                        ALBEDOSC, LEGENSC, &
!                        SHRADIANCE, CONSTSOURCE, SOURCE, &
!                        OLDNPTS, DELSOURCE, &
!                        DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
!    deallocate (SHRADIANCE)

!    SOLCRIT = 1.0
!    TRUESOLCRIT = 1.0
!    CURSPLITACC = 1.0
!    LASTNPTS = -1
!    OLDNPTS = NPTS
!    RADIANCE(:, :, :) = 0.0
!    DELSOURCE(:, :) = 0.0
!    ITER = 0
!    if (MAXITER <= ITER) return

!    if (PRNT) write (*, *) '! Iter Log(Sol)  Log(True)  Npoints'

!    ! Main solution loop
!    ! Iterate until the iterations exceed the limit, or the desired
!    ! solution criterion is reached
!    do while (ITER < MAXITER .and. (SOLCRIT > SOLACC .or. NPTS /= OLDNPTS))
!        ITER = ITER + 1

!        ! Transform the source from spherical harmonics to discrete ordinates
!        allocate (DOSOURCE(NPHI0MAX, NMU, MAXIG))
!        call SH_TO_DO(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
!                      FFTFLAG, CMU1, CPHI1, WPHISAVE, SOURCE, DOSOURCE)
!        ! Enforce the discrete ordinate source function to be non-negative
!        DOSOURCE(:, :, :) = max(0.0, DOSOURCE(:, :, :))

!        ! Keep on doing integrations as long as new levels are added
!        OLDNPTS = NPTS
!        PREVNPTS = 0
!        do while (NPTS > PREVNPTS)
!            ! Integrate the source function along discrete ordinates to
!            ! compute the radiance field.
!            call PATH_INTEGRATION(NLAY, NPTS, IXP, IXG, TAUGSC, &
!                                  NMU, NPHI0MAX, NPHI0, &
!                                  SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
!                                  SFCTYPE, SFCPARMS, SFCTEMP, UNITS, WAVELEN, WAVENO, &
!                                  MU, PHI, WTDO, SRCINT, LASTNPTS, &
!                                  FLUXUP, FLUXDN, MEANRAD, DOSOURCE, RADIANCE)
!            PREVNPTS = NPTS
!            ! Do the adaptive grid cell splitting stuff if desired
!            call SPLIT_CELLS(NLAY, NPTS, MAXIG, IXP, IXG, TAUGSC, &
!                             NMU, NPHI0MAX, NPHI0, MU, DOSOURCE, RADIANCE, &
!                             SPLITACC, TRUESOLCRIT, CURSPLITACC)
!            if (NPTS > PREVNPTS) then
!                CONSTSOURCE(1, PREVNPTS + 1:NPTS) = -1.0
!                call PREPARE_SOURCE(NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
!                                    SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
!                                    SOLARFLUX, SOLARMU, YLMSUN, UNITS, WAVELEN, WAVENO, &
!                                    CONSTSOURCE)
!            end if
!        end do
!        deallocate (DOSOURCE)

!        ! Transform the radiance from discrete ordinates to spherical harmonics
!        allocate (SHRADIANCE(NLM, NPTS))
!        call DO_TO_SH(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
!                      FFTFLAG, CMU2, CPHI2, WPHISAVE, RADIANCE, SHRADIANCE)

!        ! Calculate the source function from the radiance in spherical harmonics.
!        !   Also computes the solution criterion and the series acceleration
!        !   stuff (DELJDOT, DELJOLD, DELJNEW, DELSOURCE).
!        call COMPUTE_SOURCE(ML, MM, NLM, NLAY, NPTS, IXP, IXG, &
!                            ALBEDOSC, LEGENSC, &
!                            SHRADIANCE, CONSTSOURCE, SOURCE, &
!                            OLDNPTS, DELSOURCE, &
!                            DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
!        deallocate (SHRADIANCE)

!        ! Accelerate the convergence of the source function vector
!        if (ACCELFLAG) then
!            call ACCELERATE_SOLUTION(NLM, NPTS, OLDNPTS, SOURCE, DELSOURCE, &
!                                     DELJDOT, DELJOLD, DELJNEW)
!        end if
!        if (DELJNEW < DELJOLD) &
!            TRUESOLCRIT = SOLCRIT/(1.0 - sqrt(DELJNEW/DELJOLD))

!        ! Print out the iteration, log solution criterion, and number of points
!        if (PRNT) write (*, '(2X,I4,2F8.3,1X,I6)') ITER, &
!            log10(max(SOLCRIT, 1.0e-20)), log10(max(TRUESOLCRIT, 1.0e-20)), NPTS
!    end do

!    if (PRNT) write (*, '(1X,A,I6,A,F9.6)') &
!        '! Iterations: ', ITER, '     Final Criterion: ', SOLCRIT


!    deallocate (MU, PHI, NPHI0, WTDO, YLMSUN, WTMU)
!    deallocate (CPHI1, CPHI2, WPHISAVE, CMU1, CMU2)
!    return

!end subroutine SOLVE_SHDOMPP

subroutine DELTA_SCALE(NLAY, ML, TAUP, ALBEDOP, MAXLEG, NLEGP, LEGENP, &
                       TAUPSC, ALBEDOSC, LEGENSC)
    ! Delta scales the layer optical depth, albedo, and Legendre terms;
    ! only the 0 to ML LEGEN terms are scaled.
    implicit none
    integer, intent(in) :: NLAY, ML, MAXLEG, NLEGP(NLAY)
    real, intent(in) ::  TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
    real, intent(out) ::  TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    integer :: I, L
    real    :: F

    do I = 1, NLAY
        if (MAXLEG < ML + 1) then
            F = 0.0
            LEGENSC(:, I) = 0.0
        else
            F = LEGENP(ML + 1, I)/(2*(ML + 1) + 1)
        end if
        do L = 0, min(ML, MAXLEG)
            LEGENSC(L, I) = (2*L + 1)*(LEGENP(L, I)/(2*L + 1) - F)/(1 - F)
        end do
        TAUPSC(I) = (1.0 - ALBEDOP(I)*F)*TAUP(I)
        ALBEDOSC(I) = (1.0 - F)*ALBEDOP(I)/(1.0 - ALBEDOP(I)*F)
    end do
end subroutine DELTA_SCALE

subroutine UNDELTA_SCALE(NLAY, TAUPSC, TAUP, IXP, NPTS, TAUGSC, TAUG)
    ! Undoes the delta scaling to get the optical depth grid in the
    ! original scale.  Uses the scaled and unscaled input layer optical depths.
    implicit none
    integer, intent(in) :: NLAY, IXP(2, NLAY), NPTS
    real, intent(in) :: TAUPSC(NLAY), TAUP(NLAY), TAUGSC(NPTS)
    real, intent(out) :: TAUG(NPTS)
    integer :: LAY, I
    real    :: R

    TAUG(1) = 0.0
    do LAY = 1, NLAY
        if (TAUPSC(LAY) > 0.0) then
            R = TAUP(LAY)/TAUPSC(LAY)
        else
            R = 1.0
        end if
        do I = IXP(1, LAY), IXP(2, LAY) - 1
            TAUG(I + 1) = TAUG(I) + R*(TAUGSC(I + 1) - TAUGSC(I))
        end do
        if (LAY < NLAY) TAUG(IXP(1, LAY + 1)) = TAUG(IXP(2, LAY))
    end do
end subroutine UNDELTA_SCALE

subroutine OUTPUT_HEIGHT_GRID(NLAY, TAUP, IXP, NPTS, TAUG, HEIGHTP, HEIGHTG)
    ! Makes the output height grid (HEIGHTG) by interpolating between the
    ! levels in the input height grid (HEIGHTP) according to the optical depth
    ! grid (TAUG) assuming constant extinction in each input layer.
    ! The levels in HEIGHTG correspond to the optical depths in TAUG.
    implicit none
    integer, intent(in) :: NLAY, IXP(2, NLAY), NPTS
    real, intent(in) :: TAUP(NLAY), TAUG(NPTS), HEIGHTP(NLAY + 1)
    real, intent(out) :: HEIGHTG(NPTS)
    integer :: LAY, I
    real    :: E

    do LAY = 1, NLAY
        HEIGHTG(IXP(1, LAY)) = HEIGHTP(LAY)
        HEIGHTG(IXP(2, LAY)) = HEIGHTP(LAY + 1)
        if (TAUP(LAY) > 0.0) then
            E = (HEIGHTP(LAY) - HEIGHTP(LAY + 1))/TAUP(LAY)
            do I = IXP(1, LAY) + 1, IXP(2, LAY) - 1
                HEIGHTG(I) = HEIGHTP(LAY) - E*(TAUG(I) - TAUG(IXP(1, LAY)))
            end do
        else
            do I = IXP(1, LAY) + 1, IXP(2, LAY) - 1
                HEIGHTG(I) = HEIGHTP(LAY) + (HEIGHTP(LAY + 1) - HEIGHTP(LAY))* &
                             (I - IXP(1, LAY))/FLOAT(IXP(2, LAY) - IXP(1, LAY))
            end do
        end if
    end do
end subroutine OUTPUT_HEIGHT_GRID

subroutine MAKE_TAU_GRID(NLAY, TAUPSC, ALBEDOSC, MAXDELTAU, MINDELTAUSPLIT, &
                         NPTS, TAUGSC, IXP, IXG)
    ! Makes the initial optical depth grid, TAUGSC, of increasing values
    ! (representing decreasing altitudes in the atmosphere). TAUGSC(1) is zero.
    ! The "base" tau grid is the input property grid boundaries and even
    ! subdivision of each of these layers so that the sublayers are less
    ! than MAXDELTAU scattering optical depth.   All layers with scaled
    ! optical depth greater than MINDELTAUSPLIT are split into two sublayers.
    ! The IXP index array points to the top and bottom TAUGSC levels for each
    ! property layer. The IXG index array points to the source and
    ! radiance arrays for each TAUGSC level.
    implicit none
    integer, intent(in)  :: NLAY
    real, intent(in)  :: TAUPSC(NLAY), ALBEDOSC(NLAY)
    real, intent(in)  :: MAXDELTAU, MINDELTAUSPLIT
    real, intent(out) :: TAUGSC(*)
    integer, intent(out) :: NPTS, IXP(2, NLAY), IXG(*)
    integer :: I, L, N
    real    :: TAU

    NPTS = 0
    TAU = 0.0
    do L = 1, NLAY
        NPTS = NPTS + 1
        TAUGSC(NPTS) = TAU
        IXP(1, L) = NPTS
        if (TAUPSC(L) > MINDELTAUSPLIT) then
            N = max(1, int(ALBEDOSC(L)*TAUPSC(L)/MAXDELTAU))
        else
            N = 0
        end if
        do I = 1, N
            NPTS = NPTS + 1
            TAUGSC(NPTS) = TAUGSC(IXP(1, L)) + TAUPSC(L)*FLOAT(I)/(N + 1)
        end do
        NPTS = NPTS + 1
        IXP(2, L) = NPTS
        TAUGSC(NPTS) = TAUGSC(IXP(1, L)) + TAUPSC(L)
        TAU = TAUGSC(NPTS)
    end do
    IXG(1:NPTS) = (/(I, I=1, NPTS)/)
end subroutine MAKE_TAU_GRID

subroutine PREPARE_SOURCE(NLAY, TAUPSC, ALBEDOSC, TEMPP, LEGENSC, &
                          SRCTYPE, ML, MM, NLM, NPTS, TAUGSC, IXP, IXG, &
                          SOLARFLUX, SOLARMU, YLMSUN, UNITS, WAVELEN, WAVENO, &
                          CONSTSOURCE)
    ! Prepares the constant part of the source function for each level.
    ! This is either or both thermal emission or the single scattering
    ! solar pseudo-source.  Only those levels (points) with CONSTSOURCE(1,:)
    ! less than zero are changed, which allows this to be run for each time
    ! adaptive grid points are added.
    ! The thermal emission source is \sqrt{4\pi}(1-\omega)B(T), where the
    ! temperature T is assumed to vary linearly in optical depth.
    ! The solar source is
    !   S_0/\mu_0 e^{-\tau/\mu_0} Y_{lm}(\mu_0,\phi_0) {\omega \chi_l \over 2l+1}
    ! where S_0 is the solar flux on a horizontal plane, \chi_l are the Legendre
    ! coefficients, and \tau is the optical depth at the level.

    use shd_util

    implicit none
    integer, intent(in) :: NLAY, ML, MM, NLM, NPTS
    integer, intent(in) :: IXP(2, NLAY), IXG(NPTS)
    real, intent(in) :: SOLARFLUX, SOLARMU, YLMSUN(NLM), WAVELEN, WAVENO(2)
    real, intent(in) :: TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(in) :: TEMPP(NLAY + 1), TAUGSC(NPTS)
    real, intent(inout) :: CONSTSOURCE(NLM, NPTS)
    character(len=1), intent(in) :: SRCTYPE, UNITS
    integer :: I, J, K, L, LAY, ME, M, LOFJ(NLM)
    real    :: S0, C, TAU1, U, TEMP, BB, D, PLANCKFUNC

    S0 = SOLARFLUX/abs(SOLARMU)
    C = sqrt(4.0*acos(-1.0))
    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        do M = 0, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    do LAY = 1, NLAY
        TAU1 = TAUGSC(IXP(1, LAY))
        do I = IXP(1, LAY), IXP(2, LAY)
            K = IXG(I)
            if (CONSTSOURCE(1, K) < 0) then
                CONSTSOURCE(:, K) = 0.0
                if (SRCTYPE /= 'S') then
                    if (TAUPSC(LAY) > 0.0) then
                        U = (TAUGSC(I) - TAU1)/TAUPSC(LAY)
                    else
                        U = (I - IXP(1, LAY))/FLOAT(IXP(2, LAY) - IXP(1, LAY))
                    end if
                    TEMP = TEMPP(LAY) + U*(TEMPP(LAY + 1) - TEMPP(LAY))
                    call shd_util__computePlanck(TEMP, UNITS, WAVENO, WAVELEN, PLANCKFUNC)
                    CONSTSOURCE(1, K) = C*(1.0 - ALBEDOSC(LAY)) &
                                        *PLANCKFUNC
                end if
                if (SRCTYPE /= 'T') then
                    D = S0*exp(-TAUGSC(I)/abs(SOLARMU))
                    do J = 1, NLM
                        L = LOFJ(J)
                        CONSTSOURCE(J, K) = CONSTSOURCE(J, K) &
                                            + D*YLMSUN(J)*ALBEDOSC(LAY)*LEGENSC(L, LAY)/(2*L + 1)
                    end do
                end if
            end if
        end do
    end do
end subroutine PREPARE_SOURCE

subroutine MAKE_ANGLE_SET(NMU, NPHI, NCS, NPHI0MAX, ITYPE, &
                          NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)
    ! Make the set of angles for the discrete space representation.
    ! The number of mu's (cosine zenith angles) and maximum number of
    ! phi's (azimuth angles) is input.  The actual number of azimuthal
    ! angles for each mu is output (NPHI0).  There are three types of
    ! discrete ordinate sets: ITYPE=1 is a gaussian grid, 2 is a reduced
    ! gaussian grid, and 3 is a reduced double gaussian set.
    ! If NCS=1 then only about half the azimuthal angles (from 0 to pi)
    ! are used because the radiance is even in phi (cosine terms).
    ! The output is the NMU mu values, the NPHI0 phi values for each mu,
    ! and the integration weight for each ordinate. The first NMU/2 mu
    ! angles are the downwelling (mu<0) angles.  Also output are the
    ! the flags for doing an azimuthal FFT for each mu and the total
    ! number of angles (NANG).

    use shd_util

    implicit none
    integer, intent(in)  :: NMU, NPHI, NCS, NPHI0MAX, ITYPE
    integer, intent(out) :: NPHI0(NMU), NANG
    logical, intent(out) :: FFTFLAG(NMU)
    real, intent(out) :: MU(NMU), WTMU(NMU)
    real, intent(out) :: PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    integer :: J, K, MM
    real    :: DELPHI
    integer, parameter :: MAXNPHI = 1024
    integer :: GOODNFFT(MAXNPHI)

    call GOOD_FFT_LENGTHS(MAXNPHI, GOODNFFT)
    GOODNFFT(1:14) = (/(K, K=1, 14)/)

    MM = max(0, int(NPHI/2) - 1)
    NANG = 0
    if (ITYPE .ge. 1 .or. ITYPE .le. 3) then
        if (ITYPE .le. 2) then
            call shd_util__computeGausQuad(NMU, MU, WTMU)
        else if (ITYPE .eq. 3) then
            call shd_util__computeGausQuadDble(NMU, MU, WTMU)
        end if
        do J = 1, NMU
            if (NCS .eq. 1) then
                if (NPHI0MAX .ne. int((NPHI + 2)/2)) stop 'MAKE_ANGLE_SET: bad NPHI0MAX'
                if (ITYPE .eq. 1) then
                    NPHI0(J) = NPHI0MAX
                else
                    ! For the reduced gaussian set, make the smaller NPHI0 values
                    ! still be good for the FFT (but don't let NPHI0MAX be exceeded).
                    if (NPHI0MAX .gt. MAXNPHI) stop 'MAKE_ANGLE_SET: exceeded GOODNFFT'
                    NPHI0(J) = int(0.9 + 1 + (NPHI0MAX - 1)*sqrt(1 - MU(J)**2))
                    if (NPHI0(J) .gt. 1) NPHI0(J) = min(NPHI0MAX, GOODNFFT(NPHI0(J) - 1) + 1)
                end if
                ! Compute the azimuth angles and weights
                DELPHI = acos(-1.0)/max(1, NPHI0(J) - 1)
                do K = 1, NPHI0(J)
                    PHI(K, J) = (K - 1)*DELPHI
                    if ((K .eq. 1 .or. K .eq. NPHI0(J)) .and. NPHI0(J) .ne. 1) then
                        WTDO(J, K) = DELPHI*WTMU(J)
                    else
                        WTDO(J, K) = 2.0*DELPHI*WTMU(J)
                    end if
                end do
            else
                if (NPHI0MAX .ne. NPHI) stop 'MAKE_ANGLE_SET: bad NPHI0MAX'
                if (ITYPE .eq. 1) then
                    NPHI0(J) = NPHI0MAX
                else
                    if (NPHI0MAX .gt. MAXNPHI) stop 'MAKE_ANGLE_SET: exceeded GOODNFFT'
                    NPHI0(J) = int(0.9 + NPHI0MAX*sqrt(1 - MU(J)**2))
                    NPHI0(J) = min(NPHI0MAX, GOODNFFT(NPHI0(J)))
                end if
                DELPHI = 2.0*acos(-1.0)/NPHI0(J)
                do K = 1, NPHI0(J)
                    PHI(K, J) = (K - 1)*DELPHI
                    WTDO(J, K) = DELPHI*WTMU(J)
                end do
            end if
            NANG = NANG + NPHI0(J)
            FFTFLAG(J) = (NPHI0(J) .gt. 14) .or. (MM .gt. 16)
        end do

    else
        stop 'MAKE_ANGLE_SET: invalid discrete ordinate type'
    end if
end subroutine MAKE_ANGLE_SET

subroutine GOOD_FFT_LENGTHS(NMAX, GOODNFFT)
    ! Returns a list in GOODNFFT(I) that rounds up to the nearest I that
    ! has only divisors of 2 or 3.
    implicit none
    integer, intent(in)  :: NMAX
    integer, intent(out) :: GOODNFFT(NMAX)
    integer :: I, J, ILAST
    integer, parameter :: PRIMES(310) = &
                          (/2, 3, 5, 7, 11, 13, 17, 19, 23, 29, &
                            31, 37, 41, 43, 47, 53, 59, 61, 67, 71, &
                            73, 79, 83, 89, 97, 101, 103, 107, 109, 113, &
                            127, 131, 137, 139, 149, 151, 157, 163, 167, 173, &
                            179, 181, 191, 193, 197, 199, 211, 223, 227, 229, &
                            233, 239, 241, 251, 257, 263, 269, 271, 277, 281, &
                            283, 293, 307, 311, 313, 317, 331, 337, 347, 349, &
                            353, 359, 367, 373, 379, 383, 389, 397, 401, 409, &
                            419, 421, 431, 433, 439, 443, 449, 457, 461, 463, &
                            467, 479, 487, 491, 499, 503, 509, 521, 523, 541, &
                            547, 557, 563, 569, 571, 577, 587, 593, 599, 601, &
                            607, 613, 617, 619, 631, 641, 643, 647, 653, 659, &
                            661, 673, 677, 683, 691, 701, 709, 719, 727, 733, &
                            739, 743, 751, 757, 761, 769, 773, 787, 797, 809, &
                            811, 821, 823, 827, 829, 839, 853, 857, 859, 863, &
                            877, 881, 883, 887, 907, 911, 919, 929, 937, 941, &
                            947, 953, 967, 971, 977, 983, 991, 997, 1009, 1013, &
                            1019, 1021, 1031, 1033, 1039, 1049, 1051, 1061, 1063, 1069, &
                            1087, 1091, 1093, 1097, 1103, 1109, 1117, 1123, 1129, 1151, &
                            1153, 1163, 1171, 1181, 1187, 1193, 1201, 1213, 1217, 1223, &
                            1229, 1231, 1237, 1249, 1259, 1277, 1279, 1283, 1289, 1291, &
                            1297, 1301, 1303, 1307, 1319, 1321, 1327, 1361, 1367, 1373, &
                            1381, 1399, 1409, 1423, 1427, 1429, 1433, 1439, 1447, 1451, &
                            1453, 1459, 1471, 1481, 1483, 1487, 1489, 1493, 1499, 1511, &
                            1523, 1531, 1543, 1549, 1553, 1559, 1567, 1571, 1579, 1583, &
                            1597, 1601, 1607, 1609, 1613, 1619, 1621, 1627, 1637, 1657, &
                            1663, 1667, 1669, 1693, 1697, 1699, 1709, 1721, 1723, 1733, &
                            1741, 1747, 1753, 1759, 1777, 1783, 1787, 1789, 1801, 1811, &
                            1823, 1831, 1847, 1861, 1867, 1871, 1873, 1877, 1879, 1889, &
                            1901, 1907, 1913, 1931, 1933, 1949, 1951, 1973, 1979, 1987, &
                            1993, 1997, 1999, 2003, 2011, 2017, 2027, 2029, 2039, 2053/)

    GOODNFFT(:) = (/(I, I=1, NMAX)/)
    do I = 3, size(PRIMES)
        do J = PRIMES(I), NMAX, PRIMES(I)
            GOODNFFT(J) = 0
        end do
    end do
    ILAST = NMAX
    if (GOODNFFT(ILAST) == 0) stop 'GOOD_FFT_LENGTHS: NMAX not good'
    do I = NMAX, 1, -1
        if (GOODNFFT(I) == 0) then
            GOODNFFT(I) = ILAST
        else
            ILAST = I
        end if
    end do
end subroutine GOOD_FFT_LENGTHS

subroutine MAKE_SH_DO_COEF(ML, MM, NLM, NMU, NPHI0, &
                           NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG, &
                           CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)
    ! Makes the transformation coefficients for the spherical harmonic
    ! transform.  The main coefficients are output in four arrays: 1 is for
    ! the SH_TO_DO forward transform, 2 is for the DO_TO_SH back transform,
    ! which contains the discrete angle integration weights.
    ! The coefficients are divided up into the mu dependent set CMUn
    ! (function of l, m, mu_j), and the phi dependent set CPHIn
    ! (function of m, phi_k) for each mu_j.
    ! The FFTPACK phase coefficients for the FFT in azimuth are also
    ! output in WPHISAVE.
    implicit none
    integer, intent(in) :: ML, MM, NLM, NMU, NPHI0(NMU), NPHI0MAX
    logical, intent(in) :: FFTFLAG(NMU)
    real, intent(in) :: MU(NMU), WTMU(NMU)
    real, intent(in) :: PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(out) :: CMU1(NLM, NMU), CMU2(NMU, NLM)
    real, intent(out) :: CPHI1(0:16, 32, NMU), CPHI2(32, 0:16, NMU)
    real, intent(out) :: WPHISAVE(3*NPHI0MAX + 15, NMU)
    integer :: I, K, M, ier
    real    :: W

    ! Make the to and from associate Legendre coefficients
    do I = 1, NMU
        call YLMALL(MU(I), 0.0, ML, MM, -1, CMU1(1, I))
    end do
    do I = 1, NMU
        CMU2(I, 1:NLM) = CMU1(1:NLM, I)*WTMU(I)
    end do

    ! Make the to and from Fourier coefficients for each mu
    do I = 1, NMU
        if (.not. FFTFLAG(I)) then
            ! If not doing an FFT for this mu then make the DFT coefficient
            W = 1.0/WTMU(I)
            do K = 1, NPHI0(I)
                if (NPHI0(I) .gt. 32 .or. MM .gt. 16) &
                    stop 'MAKE_SH_DO_COEF: Fourier coefficients array exceeded'
                do M = 0, MM
                    CPHI1(M, K, I) = cos(M*PHI(K, I))
                    CPHI2(K, M, I) = CPHI1(M, K, I)*WTDO(I, K)*W
                end do
            end do
        else
            ! Precompute the phase factors for the FFTs
            call cost1i(NPHI0(I), WPHISAVE(1, I), 3*NPHI0MAX + 15, ier)
            if (ier /= 0) stop
        end if
    end do
end subroutine MAKE_SH_DO_COEF

subroutine SH_TO_DO(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                    FFTFLAG, CMU1, CPHI1, WSAVE, INDATA, OUTDATA)
    ! Transforms the input data from spherical harmonic space to discrete
    ! ordinate space for all points in the array.  Successive transforms are
    ! done in zenith and then azimuth angle.  The Fourier transform in
    ! azimuth is done using a matrix multiply DFT for smaller number of
    ! angles (NPHI0*(MM+1)<160) and with an FFT for larger number of
    ! angles.  The FFT is used for any number of angles above the limit,
    ! but is much more efficient for NPHI0 whose factors are small primes
    ! (powers of 2 are best).  The number of floating point operations for
    ! the DFT is 2*Nmu*Nphi0*Nm. The number of floating point operations for
    ! the zenith angle transform is 2*Nmu*Nlm.
    implicit none
    integer, intent(in) :: NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
    logical, intent(in) :: FFTFLAG(NMU)
    real, intent(in) :: INDATA(NLM, NPTS)
    real, intent(in) :: CMU1(NLM, NMU), CPHI1(0:16, 32, NMU)
    real, intent(in) :: WSAVE(3*NPHI0MAX + 15, NMU)
    real, intent(out) :: OUTDATA(NPHI0MAX, NMU, NPTS)
    integer :: I, J, L, M, ME, IMU, IPHI, ier
    integer :: MOFJ(NLM)
    real    :: S(0:MM), WORK(NPHI0MAX)

    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = 0, min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        do IMU = 1, NMU
            ! Figure the max Fourier azimuthal mode we can do for this Nphi0
            ME = max(0, min(NPHI0(IMU) - 2, MOFJ(NLM)))
            ! First do Legendre transform by summing over l for each m.
            S(:) = 0.0
            do J = 1, NLM
                M = MOFJ(J)
                S(M) = S(M) + CMU1(J, IMU)*INDATA(J, I)
            end do
            ! Then do the Fourier transform from m to phi for each mu
            if (FFTFLAG(IMU)) then
                ! For cosine modes only use a cosine FFT
                OUTDATA(1:ME + 1, IMU, I) = S(0:ME)
                OUTDATA(ME + 2:NPHI0(IMU), IMU, I) = 0.0
                call cost1b(NPHI0(IMU), 1, OUTDATA(:, IMU, I), NPHI0MAX, &
                            WSAVE(:, IMU), 3*NPHI0MAX + 15, WORK, NPHI0MAX, ier)
                if (ier /= 0) stop
            else
                ! Else do a slow DFT
                do IPHI = 1, NPHI0(IMU)
                    OUTDATA(IPHI, IMU, I) = sum(CPHI1(0:ME, IPHI, IMU)*S(0:ME))
                end do
            end if
        end do
    end do
end subroutine SH_TO_DO

subroutine DO_TO_SH(NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0, &
                    FFTFLAG, CMU2, CPHI2, WSAVE, INDATA, OUTDATA)
    ! Transforms the input field from discrete ordinate space to spherical
    ! harmonic space for all the points.  (See SH_TO_DO for more).
    ! The SELECTPNTS logical array indicates whether a point should be
    ! transformed or to use the previous points (for the input layer interfaces).
    implicit none
    integer, intent(in) :: NPTS, ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
    logical, intent(in) :: FFTFLAG(NMU)
    real, intent(in) :: INDATA(NPHI0MAX, NMU, NPTS)
    real, intent(in) :: CMU2(NMU, NLM), CPHI2(32, 0:16, NMU)
    real, intent(in) :: WSAVE(3*NPHI0MAX + 15, NMU)
    real, intent(out) :: OUTDATA(NLM, NPTS)
    integer :: I, IMU, IPHI, J, L, M, ME, N, ier
    real    :: PI, DELPHI
    integer :: MOFJ(NLM)
    real    :: S(0:NPHI0MAX - 1), WORK(NPHI0MAX)

    PI = acos(-1.0)
    ! Make the M for each J array
    J = 0
    do L = 0, ML
        do M = 0, min(L, MM)
            J = J + 1
            MOFJ(J) = M
        end do
    end do

    ! Do the transform for each grid point
    do I = 1, NPTS
        OUTDATA(:, I) = 0.0
        do IMU = 1, NMU
            N = NPHI0(IMU)
            ! Figure the max Fourier azimuthal mode we can do for this Nphi0
            ME = max(0, min(N - 2, MOFJ(NLM)))
            ! First do Fourier transform from phi to m for each mu
            if (FFTFLAG(IMU)) then
                ! For cosine modes only use a cosine FFT
                S(0:N - 1) = INDATA(1:N, IMU, I)
                call cost1f(N, 1, S(0:N - 1), NPHI0MAX, &
                            WSAVE(:, IMU), 3*NPHI0MAX + 15, WORK, NPHI0MAX, ier)
                if (ier /= 0) stop
                S(0:ME) = PI*S(0:ME)
                S(0) = 2.0*S(0)
                S(ME) = 2.0*S(ME)
            else
                ! Else do a slow DFT
                S(:) = 0.0
                do M = 0, ME
                    S(M) = sum(CPHI2(1:N, M, IMU)*INDATA(1:N, IMU, I))
                end do
            end if
            ! Then do Legendre transform by adding to output for each l, m
            do J = 1, NLM
                M = MOFJ(J)
                OUTDATA(J, I) = OUTDATA(J, I) + CMU2(IMU, J)*S(M)
            end do
        end do
    end do
end subroutine DO_TO_SH

subroutine COMPUTE_SOURCE(ML, MM, NLM, NLAY, NPTS, IXP, IXG, &
                          ALBEDOSC, LEGENSC, &
                          SHRADIANCE, CONSTSOURCE, SOURCE, &
                          OLDNPTS, DELSOURCE, &
                          DELJDOT, DELJOLD, DELJNEW, SOLCRIT)
    ! Computes the source function (SOURCE) in spherical harmonic space
    ! for all the grid levels in all the input layers.
    ! The thermal source and/or solar pseudo-sources in CONSTSOURCE is
    ! added to the scattering source (computed from LEGENP and the spherical
    ! harmonic expansion in SHRADIANCE).
    ! To save memory this routine also computes some things for the sequence
    ! acceleration: the dot products of the difference in source function
    ! between successive iterations and the new difference source function
    ! vector (DELSOURCE).
    ! Computes the new solution criterion (SOLCRIT, the RMS difference
    ! in succesive source function fields normalized by the RMS of the field
    implicit none
    integer, intent(in) :: ML, MM, NLM, NLAY, NPTS, OLDNPTS
    integer, intent(in) :: IXP(2, NLAY), IXG(NPTS)
    real, intent(in) :: ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(in) :: SHRADIANCE(NLM, NPTS), CONSTSOURCE(NLM, NPTS)
    real, intent(inout) :: SOURCE(NLM, NPTS), DELSOURCE(NLM, OLDNPTS)
    real, intent(out) :: DELJDOT, DELJOLD, DELJNEW, SOLCRIT
    integer :: LOFJ(NLM), I, J, K, L, LAY, M, ME
    real    :: SOURCET, NORM

    ! Make the l index as a function of SH term (J)
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        do M = 0, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    DELJDOT = 0.0
    DELJOLD = 0.0
    DELJNEW = 0.0
    NORM = 0.0
    do LAY = 1, NLAY
        do I = IXP(1, LAY), IXP(2, LAY)
            K = IXG(I)
            do J = 1, NLM
                L = LOFJ(J)
                SOURCET = CONSTSOURCE(J, K) &
                          + SHRADIANCE(J, K)*ALBEDOSC(LAY)*LEGENSC(L, LAY)/(2*L + 1)
                if (K <= OLDNPTS) then
                    DELJDOT = DELJDOT + (SOURCET - SOURCE(J, K))*DELSOURCE(J, K)
                    DELJOLD = DELJOLD + DELSOURCE(J, K)**2
                    DELSOURCE(J, K) = SOURCET - SOURCE(J, K)
                    DELJNEW = DELJNEW + (SOURCET - SOURCE(J, K))**2
                    NORM = NORM + SOURCET**2
                end if
                SOURCE(J, K) = SOURCET
            end do
        end do
    end do
    if (NORM > 0.0) then
        SOLCRIT = sqrt(DELJNEW/NORM)
    else
        SOLCRIT = 1.0
    end if
end subroutine COMPUTE_SOURCE

subroutine ACCELERATE_SOLUTION(NLM, NPTS, OLDNPTS, SOURCE, DELSOURCE, &
                               DELJDOT, DELJOLD, DELJNEW)
    ! Accelerates the successive order of scattering series.  Uses
    ! information about the ratio in lengths and angle between successive
    ! source function differences to decide how far along the current
    ! source function difference to extrapolate.
    implicit none
    integer, intent(in) :: NLM, NPTS, OLDNPTS
    real, intent(in) :: DELSOURCE(NLM, OLDNPTS)
    real, intent(inout) :: SOURCE(NLM, NPTS)
    real, intent(in) :: DELJDOT, DELJOLD, DELJNEW
    integer :: I, J, IS, ISD, NS, NSD, NSC
    real, save :: A, R, THETA

    ! Accelerate if the grids are the same, didn't last time, and
    !   things are converging.
    if (NPTS == OLDNPTS .and. A == 0.0 .and. DELJNEW < DELJOLD) then
        ! Compute the acceleration extrapolation factor and apply it.
        R = sqrt(DELJNEW/DELJOLD)
        THETA = acos(DELJDOT/sqrt(DELJOLD*DELJNEW))
        A = (1 - R*cos(THETA) + R**(1 + 0.5*3.14159/THETA)) &
            /(1 + R**2 - 2*R*cos(THETA)) - 1.0
        A = max(0.0, A)
        ! WRITE (*,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
        SOURCE(:, :) = SOURCE(:, :) + A*DELSOURCE(:, :)
    else
        A = 0.0
    end if
end subroutine ACCELERATE_SOLUTION

subroutine PATH_INTEGRATION(NLAY, NPTS, IXP, IXG, TAUGSC, &
                            NMU, NPHI0MAX, NPHI0, &
                            SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                            SFCTYPE, SFCPARMS, SFCTEMP, UNITS, WAVELEN, WAVENO, &
                            MU, PHI, WTDO, SRCINT, LASTNPTS, &
                            FLUXUP, FLUXDN, MEANRAD, DOSOURCE, RADIANCE)
    ! Performs the path integrations through the medium specified by
    ! the optical depth (TAUGSC) and source function (DOSOURCE) in
    ! discrete ordinates.  The integrations are done from the top boundary
    ! (with no reflection) down to the bottom boundary for the downwelling
    ! angles, and then up from the bottom to the top after computing the
    ! surface reflection/emission.   The hemispheric flux (FLUXES) and mean
    ! radiance (MEANRAD) are also computed.
    ! The discrete ordinates in MU, PHI must be a complete set with the
    ! downwelling (mu<0) angles first and downwelling and upwelling matching
    ! [MU(J) = -MU(NMU/2+J)].
    ! The surface reflection is handled differently for Lambertian surfaces
    ! and more general surfaces specified by bidirectional reflection
    ! distribution functions (BRDF).  For the general BRDF, the bottom
    ! boundary radiance must be computed for each upwelling ordinate, while
    ! for the Lambertian surface the boundary radiances can be computed just
    ! once, since they are isotropic.

    use shd_util

    implicit none
    integer, intent(in) :: NLAY, NPTS, IXP(2, NLAY), IXG(NPTS)
    integer, intent(in) :: NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(inout) :: LASTNPTS
    real, intent(in) :: SOLARFLUX, SOLARMU, SKYRAD
    real, intent(in) :: SFCPARMS(*), SFCTEMP, WAVELEN, WAVENO(2)
    real, intent(in) :: MU(NMU), PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(inout) :: SRCINT(0:4, NMU, NPTS)
    real, intent(in) :: TAUGSC(NPTS), DOSOURCE(NPHI0MAX, NMU, NPTS)
    real, intent(out) :: FLUXUP(NPTS), FLUXDN(NPTS), MEANRAD(NPTS)
    real, intent(out) :: RADIANCE(NPHI0MAX, NMU, NPTS)
    character(len=1), intent(in) :: SRCTYPE, SFCTYPE, UNITS

    integer :: LAY, I, I1, I2, J, K, IPHI, IMU, M, N
    real    :: PI, TEMP, SFCPLANCK, ALB
    real    :: DIRFLUXSFC, TOPRAD, BOTRAD
    double precision, allocatable :: A(:, :), B(:, :), C(:, :), D(:, :), E(:)
    double precision, allocatable :: T(:), T1(:), T3(:), T4(:), TR(:)

    PI = acos(-1.0)

    M = NMU/2

    if (LASTNPTS /= NPTS) then
        ! Calculate and store the source function integration coefficients
        allocate (A(4, M), B(4, M), C(4, M), D(4, M), E(M), T(M), T1(M), T3(M), T4(M), TR(M))
        do LAY = 1, NLAY
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then ! special case for layer with just two points
                T(:) = (TAUGSC(I1 + 1) - TAUGSC(I1))/abs(MU(1:M))
                TR(:) = exp(-T(:)) ! sublayer transmission along ordinate
                if (T(1) < 1.0e-4) then
                    D(1, :) = 0.5*T(:) - (1.0d0/3)*T(:)**2
                    D(2, :) = 0.5*T(:) - (1.0d0/6)*T(:)**2
                else
                    D(1, :) = (1 - TR(:)*(1 + T(:)))/T(:)
                    D(2, :) = (T(:) - (1 - TR(:)))/T(:)
                end if
                SRCINT(0, 1:M, I1) = TR(:)
                SRCINT(0, M + 1:NMU, I1) = TR(:)
                SRCINT(1, 1:M, I1) = D(2, :)
                SRCINT(2, 1:M, I1) = D(1, :)
                SRCINT(1, M + 1:NMU, I1) = D(1, :)
                SRCINT(2, M + 1:NMU, I1) = D(2, :)
            else
                do I = I1, I2
                    A(:, :) = 0.0; B(:, :) = 0.0; C(:, :) = 0.0; D(:, :) = 0.0
                    T3(:) = (TAUGSC(I + 1) - TAUGSC(I))/abs(MU(1:M))
                    if (I == I1) then
                        T4(:) = (TAUGSC(I + 2) - TAUGSC(I))/abs(MU(1:M))
                        E = T3*T4; B(2, :) = 1/E; C(2, :) = -(T3 + T4)/E; D(2, :) = 1.0
                        E = T3*(T3 - T4); B(3, :) = 1/E; C(3, :) = -T4/E; D(3, :) = 0.0
                        E = T4*(T4 - T3); B(4, :) = 1/E; C(4, :) = -T3/E; D(4, :) = 0.0
                    else if (I < I2) then
                        T1(:) = (TAUGSC(I - 1) - TAUGSC(I))/abs(MU(1:M))
                        T4(:) = (TAUGSC(I + 2) - TAUGSC(I))/abs(MU(1:M))
                        E = T1*(T1 - T3)*(T1 - T4); A(1, :) = 1/E
                        B(1, :) = -(T3 + T4)/E; C(1, :) = T3*T4/E; D(1, :) = 0.0
                        E = -T1*T3*T4; A(2, :) = 1/E
                        B(2, :) = -(T1 + T3 + T4)/E; C(2, :) = (T1*T3 + T1*T4 + T3*T4)/E; D(2, :) = 1.0
                        E = T3*(T3 - T1)*(T3 - T4); A(3, :) = 1/E
                        B(3, :) = -(T1 + T4)/E; C(3, :) = T1*T4/E; D(3, :) = 0.0
                        E = T4*(T4 - T1)*(T4 - T3); A(4, :) = 1/E
                        B(4, :) = -(T1 + T3)/E; C(4, :) = T1*T3/E; D(4, :) = 0.0
                    else
                        T1(:) = (TAUGSC(I - 1) - TAUGSC(I))/abs(MU(1:M))
                        E = T1*(T1 - T3); B(1, :) = 1/E; C(1, :) = -T3/E; D(1, :) = 0.0
                        E = T1*T3; B(2, :) = 1/E; C(2, :) = -(T1 + T3)/E; D(2, :) = 1.0
                        E = T3*(T3 - T1); B(3, :) = 1/E; C(3, :) = -T1/E; D(3, :) = 0.0
                    end if
                    T(:) = T3(:)
                    TR(:) = exp(-T(:))
                    SRCINT(0, 1:M, I) = TR(:)
                    SRCINT(0, M + 1:NMU, I) = TR(:)
                    do J = 1, 4
                        SRCINT(J, 1:M, I) = A(J, :)*(6*TR + T**3 - 3*T**2 + 6*T - 6) &
                                            + B(J, :)*(-2*TR + T**2 - 2*T + 2) &
                                            + C(J, :)*(TR + T - 1) + D(J, :)*(1 - TR)
                        SRCINT(J, M + 1:NMU, I) = A(J, :)*(6 - TR*(T**3 + 3*T**2 + 6*T + 6)) &
                                                  + B(J, :)*(2 - TR*(T**2 + 2*T + 2)) &
                                                  + C(J, :)*(1 - TR*(T + 1)) + D(J, :)*(1 - TR)
                    end do
                end do
            end if
        end do
        deallocate (A, B, C, D, E, T, T1, T3, T4, TR)
    end if

    ! Zero the fluxes and mean radiance arrays
    FLUXUP(:) = 0.0; FLUXDN(:) = 0.0; MEANRAD(:) = 0.0

    ! Make the isotropic radiances for the top boundary
    if (SRCTYPE .eq. 'T') then
        TEMP = SKYRAD
        call shd_util__computePlanck(TEMP, UNITS, WAVENO, WAVELEN, TOPRAD)
    else
        TOPRAD = SKYRAD
    end if
    RADIANCE(:, 1:NMU/2, IXG(IXP(1, 1))) = TOPRAD

    ! First integrate the downwelling zenith angles (these angles must be first)
    ! For downwelling ordinates integrate source function from the top down
    !   Integrations are performed over the vector of azimuths
    do IMU = 1, NMU/2
        N = NPHI0(IMU)
        do LAY = 1, NLAY
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADIANCE(1:N, IMU, IXG(I1 + 1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(1, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1))
            else
                RADIANCE(1:N, IMU, IXG(I1 + 1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                                  + SRCINT(3, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1)) &
                                                  + SRCINT(4, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 2))
                do I = I1 + 1, I2 - 1
                    RADIANCE(1:N, IMU, IXG(I + 1)) = SRCINT(0, IMU, I)*RADIANCE(1:N, IMU, IXG(I)) &
                                                     + SRCINT(1, IMU, I)*DOSOURCE(1:N, IMU, IXG(I - 1)) &
                                                     + SRCINT(2, IMU, I)*DOSOURCE(1:N, IMU, IXG(I)) &
                                                     + SRCINT(3, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 1)) &
                                                     + SRCINT(4, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 2))
                end do
                RADIANCE(1:N, IMU, IXG(I2 + 1)) = SRCINT(0, IMU, I2)*RADIANCE(1:N, IMU, IXG(I2)) &
                                                  + SRCINT(1, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 - 1)) &
                                                  + SRCINT(2, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2)) &
                                                  + SRCINT(3, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 + 1))
            end if
            if (LAY < NLAY) RADIANCE(:, IMU, IXG(IXP(1, LAY + 1))) = RADIANCE(:, IMU, IXG(IXP(2, LAY)))
        end do

        ! Integrate over the DO radiances to get the downwelling fluxes
        do I = 1, NPTS
            do K = 1, NPHI0(IMU)
                FLUXDN(I) = FLUXDN(I) + abs(MU(IMU))*WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))
                MEANRAD(I) = MEANRAD(I) + WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))/(4*PI)
            end do
        end do
    end do

    ! Now do all the surface reflection stuff
    if (SRCTYPE /= 'T') then
        DIRFLUXSFC = SOLARFLUX*exp(-TAUGSC(IXP(2, NLAY))/abs(SOLARMU))
    end if
    SFCPLANCK = 0.0
    if (SRCTYPE /= 'S') then
        call shd_util__computePlanck(SFCTEMP, UNITS, WAVENO, WAVELEN, SFCPLANCK)
    end if

    ! For a Lambertian surface calculate the isotropic surface radiance
    !   from the thermal emission or direct reflection and the reflected
    !   downwelling flux at the surface.
    if (SFCTYPE == 'L') then
        ALB = SFCPARMS(1)
        if (SRCTYPE == 'T') then
            BOTRAD = (1.0 - ALB)*SFCPLANCK
        else if (SRCTYPE == 'S') then
            BOTRAD = ALB*DIRFLUXSFC/PI
        else if (SRCTYPE .eq. 'B') then
            BOTRAD = (1.0 - ALB)*SFCPLANCK + ALB*DIRFLUXSFC/PI
        end if
        BOTRAD = BOTRAD + FLUXDN(IXP(2, NLAY))*ALB/PI
        RADIANCE(:, NMU/2 + 1:NMU, IXG(IXP(2, NLAY))) = BOTRAD
    else
        ! If not a Lambertian surface, compute the radiance for each of the
        !   upwelling ordinates by integrating the BRDF over the stored
        !   downwelling radiances.
        do IMU = NMU/2 + 1, NMU
            do IPHI = 1, NPHI0(IMU)
                call VARIABLE_BRDF_SURFACE(NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                           MU(IMU), PHI(IPHI, IMU), WAVELEN, &
                                           SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                           SFCTYPE, SFCPARMS, SFCPLANCK, &
                                           RADIANCE(:, :, IXG(IXP(2, NLAY))), BOTRAD)
                RADIANCE(IPHI, IMU, IXG(IXP(2, NLAY))) = BOTRAD
            end do
        end do
    end if

    ! Then integrate the upwelling zenith angles
    ! For upwelling ordinates integrate source function from the bottom up
    do IMU = NMU/2 + 1, NMU
        N = NPHI0(IMU)
        do LAY = NLAY, 1, -1
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADIANCE(1:N, IMU, IXG(I1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1 + 1)) &
                                              + SRCINT(1, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                              + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1))
            else
                RADIANCE(1:N, IMU, IXG(I2)) = SRCINT(0, IMU, I2)*RADIANCE(1:N, IMU, IXG(I2 + 1)) &
                                              + SRCINT(1, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 - 1)) &
                                              + SRCINT(2, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2)) &
                                              + SRCINT(3, IMU, I2)*DOSOURCE(1:N, IMU, IXG(I2 + 1))
                do I = I2 - 1, I1 + 1, -1
                    RADIANCE(1:N, IMU, IXG(I)) = SRCINT(0, IMU, I)*RADIANCE(1:N, IMU, IXG(I + 1)) &
                                                 + SRCINT(1, IMU, I)*DOSOURCE(1:N, IMU, IXG(I - 1)) &
                                                 + SRCINT(2, IMU, I)*DOSOURCE(1:N, IMU, IXG(I)) &
                                                 + SRCINT(3, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 1)) &
                                                 + SRCINT(4, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 2))
                end do
                RADIANCE(1:N, IMU, IXG(I1)) = SRCINT(0, IMU, I1)*RADIANCE(1:N, IMU, IXG(I1 + 1)) &
                                              + SRCINT(2, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1)) &
                                              + SRCINT(3, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 1)) &
                                              + SRCINT(4, IMU, I1)*DOSOURCE(1:N, IMU, IXG(I1 + 2))
            end if
            if (LAY > 1) RADIANCE(:, IMU, IXG(IXP(2, LAY - 1))) = RADIANCE(:, IMU, IXG(IXP(1, LAY)))
        end do

        ! Integrate over the DO radiances to get the upwelling fluxes
        do I = 1, NPTS
            do K = 1, NPHI0(IMU)
                FLUXUP(I) = FLUXUP(I) + abs(MU(IMU))*WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))
                MEANRAD(I) = MEANRAD(I) + WTDO(IMU, K)*RADIANCE(K, IMU, IXG(I))/(4*PI)
            end do
        end do
    end do

    LASTNPTS = NPTS
end subroutine PATH_INTEGRATION

subroutine SPLIT_CELLS(NLAY, NPTS, MAXIG, IXP, IXG, TAUGSC, &
                       NMU, NPHI0MAX, NPHI0, MU, DOSOURCE, RADIANCE, &
                       SPLITACC, TRUESOLCRIT, CURSPLITACC)
    ! Determines the current cell splitting accuracy (CURSPLITACC) from the
    ! solution accuracy (SOLCRIT) and the desired final cell splitting accuracy
    ! (SPLITACC).  Calculates the cell splitting criterion for each sublayer,
    ! The cell splitting criterion for a sublayer is the rms difference over
    ! outgoing ordinates in radiance between the original cubic source function
    ! interpolation integration (in RADIANCE) and a linear interpolation
    ! integration, normalized by the rms of RADIANCE over all ordinates and
    ! levels.  Splits cells with criterion greater than the current splitting
    ! accuracy.  Interpolates the discrete ordinate source function to the new
    ! levels.  The number of levels (NPTS), the adaptive index arrays (IXP, IXG),
    ! and the optical depth grid array (TAUGSC) are also updated.
    ! The boundary condition information is not needed for this routine because
    ! it has already been computed by PATH_INTEGRATION and stored in RADIANCE.
    ! SPLITCRITRATIO is the factor needed to divide the cell splitting criterion
    ! to convert it to a typical radiance accuracy of the split sublayer with
    ! the cubic interpolation integration.

    use shd_util

    implicit none
    integer, intent(in) :: NLAY, MAXIG
    integer, intent(inout) :: NPTS, IXP(2, NLAY), IXG(*)
    integer, intent(in) :: NMU, NPHI0MAX, NPHI0(NMU)
    real, intent(in) :: MU(NMU)
    real, intent(in) :: SPLITACC, TRUESOLCRIT
    real, intent(in) :: RADIANCE(NPHI0MAX, NMU, NPTS)
    real, intent(inout) :: CURSPLITACC
    real, intent(inout) :: TAUGSC(*), DOSOURCE(NPHI0MAX, NMU, *)
    integer :: LAY, I, IMU, IMU2, N
    integer :: IXPN(2, NLAY), NEWNPTS, I1, I2, I3, I4
    real    :: SPLITCRIT, NORM, T1, T2, T3, T4
    double precision  :: TAU(NMU/2), TR(NMU/2)
    double precision, allocatable :: SRCINT1(:, :, :), DIFFRAD(:, :)

    NORM = sqrt(sum(RADIANCE(:, :, 1:NPTS)**2)/(NPTS*sum(NPHI0(:))))
    ! Leave immediately if splitting is not desired
    if (SPLITACC <= 0.0 .or. NORM == 0.0) return

    allocate (SRCINT1(0:2, NMU/2, NPTS), DIFFRAD(NPHI0MAX, NMU))
    DIFFRAD(:, :) = 0.0

    CURSPLITACC = max(SPLITACC, min(CURSPLITACC, TRUESOLCRIT))

    ! Calculate and store the linear source function integration coefficients
    do LAY = 1, NLAY
        do I = IXP(1, LAY), IXP(2, LAY) - 1
            TAU(:) = (TAUGSC(I + 1) - TAUGSC(I))/abs(MU(1:NMU/2))
            TR(:) = exp(-TAU(:)) ! sublayer transmission along ordinate
            SRCINT1(0, :, I) = TR(:)
            if (TAU(1) < 1.0e-4) then
                SRCINT1(1, :, I) = 0.5*TAU(:) - (1.0d0/3)*TAU(:)**2
                SRCINT1(2, :, I) = 0.5*TAU(:) - (1.0d0/6)*TAU(:)**2
            else
                SRCINT1(1, :, I) = (1 - TR(:)*(1 + TAU(:)))/TAU(:)
                SRCINT1(2, :, I) = (TAU(:) - (1 - TR(:)))/TAU(:)
            end if
        end do
    end do

    ! For each sublayer integrate the downwelling zenith angles from top to
    ! bottom and the upwelling angles from bottom to top using linear
    ! source function interpolation and difference from the previous
    ! cubic/quadratic interpolation results.
    NEWNPTS = NPTS
    IXPN(:, :) = IXP(:, :)
    do LAY = 1, NLAY
        if (IXP(1, LAY) < IXP(2, LAY) - 1) then
            do I = IXP(1, LAY), IXP(2, LAY) - 1
                do IMU = 1, NMU/2
                    IMU2 = IMU + NMU/2
                    N = NPHI0(IMU)
                    DIFFRAD(1:N, IMU) = -RADIANCE(1:N, IMU, IXG(I + 1)) &
                                        + SRCINT1(0, IMU, I)*RADIANCE(1:N, IMU, IXG(I)) &
                                        + SRCINT1(1, IMU, I)*DOSOURCE(1:N, IMU, IXG(I)) &
                                        + SRCINT1(2, IMU, I)*DOSOURCE(1:N, IMU, IXG(I + 1))
                    DIFFRAD(1:N, IMU2) = -RADIANCE(1:N, IMU2, IXG(I)) &
                                         + SRCINT1(0, IMU, I)*RADIANCE(1:N, IMU2, IXG(I + 1)) &
                                         + SRCINT1(1, IMU, I)*DOSOURCE(1:N, IMU2, IXG(I + 1)) &
                                         + SRCINT1(2, IMU, I)*DOSOURCE(1:N, IMU2, IXG(I))
                end do
                SPLITCRIT = sqrt(sum(DIFFRAD(:, :)**2)/sum(NPHI0(:)))/NORM
                ! If the splitting criterion is large enough then add a level in
                !   the center of the current sublayer.  Put the new level at the
                !   end of the radiance array and update the pointers.
                if (SPLITCRIT > CURSPLITACC) then
                    if (NEWNPTS >= MAXIG) then
                        write (*, *)
                        write (*, *) 'Adaptive grid exceeded MAXIG=', MAXIG
                        write (*, *) 'Increase SPLITTING_FACTOR'
                        stop
                    end if
                    NEWNPTS = NEWNPTS + 1
                    TAUGSC(NEWNPTS) = 0.5*(TAUGSC(I) + TAUGSC(I + 1))
                    ! print *, newnpts,taug(newnpts), splitcrit
                    IXG(NEWNPTS) = NEWNPTS
                    IXPN(2, LAY:NLAY) = IXPN(2, LAY:NLAY) + 1
                    IXPN(1, LAY + 1:NLAY) = IXPN(1, LAY + 1:NLAY) + 1
                    I1 = max(I - 1, IXP(1, LAY)); I2 = I
                    I3 = I + 1; I4 = min(I + 2, IXP(2, LAY))
                    call INTERP_SOURCE(NPHI0MAX, NMU, TAUGSC(I1), DOSOURCE(:, :, IXG(I1)), &
                                       TAUGSC(I2), DOSOURCE(:, :, IXG(I2)), &
                                       TAUGSC(I3), DOSOURCE(:, :, IXG(I3)), &
                                       TAUGSC(I4), DOSOURCE(:, :, IXG(I4)), &
                                       TAUGSC(NEWNPTS), DOSOURCE(:, :, IXG(NEWNPTS)))
                end if
            end do
        end if
    end do

    if (NEWNPTS > NPTS) then
        call shd_util__ssort(TAUGSC, IXG, NEWNPTS, 2)
        NPTS = NEWNPTS
        IXP(:, :) = IXPN(:, :)
        ! Make sure that the IXG pointers are in the correct order
        do I = 1, NPTS - 1
            if (TAUGSC(I) == TAUGSC(I + 1) .and. IXG(I) > IXG(I + 1)) then
                I1 = IXG(I)
                IXG(I) = IXG(I + 1)
                IXG(I + 1) = I1
            end if
        end do
    end if
    deallocate (SRCINT1, DIFFRAD)
end subroutine SPLIT_CELLS

subroutine INTERP_SOURCE(NPH, NMU, T1, S1, T2, S2, T3, S3, T4, S4, T, S)
    ! Does either quadratic or cubic polynomial interpolation of the discrete
    ! ordinate source function as a function of optical depth.  The source
    ! function (S1,S2,S3,S4) is specified at three or four optical depths
    ! (T1,T2,T3,T4), interpolated to T, and output in S.
    implicit none
    integer, intent(in) :: NPH, NMU
    real, intent(in) :: T1, T2, T3, T4, T
    real, intent(in) :: S1(NPH, NMU), S2(NPH, NMU), S3(NPH, NMU), S4(NPH, NMU)
    real, intent(out) :: S(NPH, NMU)

    if (abs(T1 - T2) < 4*spacing(T1)) then
        S = (T - T3)*(T - T4)*S2/((T2 - T3)*(T2 - T4)) &
            + (T - T2)*(T - T4)*S3/((T3 - T2)*(T3 - T4)) &
            + (T - T2)*(T - T3)*S4/((T4 - T2)*(T4 - T3))
    else if (abs(T3 - T4) < 4*spacing(T3)) then
        S = (T - T2)*(T - T3)*S1/((T1 - T2)*(T1 - T3)) &
            + (T - T1)*(T - T3)*S2/((T2 - T1)*(T2 - T3)) &
            + (T - T1)*(T - T2)*S3/((T3 - T1)*(T3 - T2))
    else
        S = (T - T2)*(T - T3)*(T - T4)*S1/((T1 - T2)*(T1 - T3)*(T1 - T4)) &
            + (T - T1)*(T - T3)*(T - T4)*S2/((T2 - T1)*(T2 - T3)*(T2 - T4)) &
            + (T - T1)*(T - T2)*(T - T4)*S3/((T3 - T1)*(T3 - T2)*(T3 - T4)) &
            + (T - T1)*(T - T2)*(T - T3)*S4/((T4 - T1)*(T4 - T2)*(T4 - T3))
    end if
end subroutine INTERP_SOURCE

subroutine VARIABLE_BRDF_SURFACE(NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                 MU2, PHI2, WAVELEN, SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                 SFCTYPE, SFCPARMS, SFCPLANCK, &
                                 SFCRADIANCE, BOTRAD)
    ! Computes the upwelling radiance at the bottom boundary for one
    ! outgoing direction using the specified bidirectional reflectance
    ! distribution function.  The upwelling radiance includes the reflection of
    ! the incident radiance, the thermal emission (emissivity*Planck function),
    ! and the reflected direct solar flux (if applicable).  The upwelling
    ! radiance is the integral over all incident directions of the BRDF times
    ! the downwelling radiance, so a discrete sum is done and the integral
    ! weights (WTDO) are included.  The general BRDF function is called
    ! to compute the reflectance for the particular BRDF type (SFCTYPE)
    ! with parameters (SFCPARMS) for the incident and outgoing directions.
    ! The emissivity is computed implicitly from the integral of the BRDF.
    ! The outgoing direction is specified with (MU2,PHI2), and the BRDF is
    ! computed for all incident directions (loop over JMU,JPHI).  The incident
    ! downwelling radiances are input in SFCRADIANCE(:,:) and the outgoing
    ! upwelling radiance is output in BOTRAD.

    use shd_sfc

    implicit none
    integer, intent(in) :: NMU, NPHI0MAX, NPHI0(NMU)
    real, intent(in) :: MU(NMU), PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(in) :: MU2, PHI2, WAVELEN, SOLARMU, DIRFLUXSFC, SFCPLANCK
    real, intent(in) :: SFCPARMS(*), SFCRADIANCE(NPHI0MAX, NMU)
    real, intent(out):: BOTRAD
    character(len=1), intent(in) :: SRCTYPE, SFCTYPE
    integer :: JMU, JPHI
    real    :: OPI, SOLARAZ, REFLECT, W, REFLECT_TMP(4, 4)

    OPI = 1.0/acos(-1.0)
    SOLARAZ = 0.0

    ! Initialize the upwelling boundary radiances to zero (for thermal)
    !   or to the reflected direct solar flux.
    if (SRCTYPE == 'T') then
        BOTRAD = 0.0
    else
        call shd_sfc__brdfSurface(SFCTYPE, SFCPARMS, WAVELEN, MU2, PHI2, SOLARMU, SOLARAZ, 1, REFLECT_TMP)
        REFLECT = REFLECT_TMP(1, 1)
        BOTRAD = OPI*REFLECT*DIRFLUXSFC
    end if

    ! Integrate over the incident discrete ordinate directions (JMU,JPHI)
    do JMU = 1, NMU/2
        do JPHI = 1, NPHI0(JMU)
            call shd_sfc__brdfSurface(SFCTYPE, SFCPARMS, WAVELEN, MU2, PHI2, MU(JMU), PHI(JPHI, JMU), 1, REFLECT_TMP)
            REFLECT = REFLECT_TMP(1, 1)
            W = OPI*abs(MU(JMU))*WTDO(JMU, JPHI)
            BOTRAD = BOTRAD + W*REFLECT*SFCRADIANCE(JPHI, JMU) &
                     + W*(1 - REFLECT)*SFCPLANCK
        end do
    end do

    return

end subroutine VARIABLE_BRDF_SURFACE

subroutine COMPUTE_RADIANCE(NLAY, ALBEDOSC, LEGENSC, &
                            MAXLEG, NLEGP, LEGENP, DELTAM, &
                            NPTS, HEIGHTG, TAUG, TAUGSC, IXP, IXG, &
                            ML, MM, NLM, NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                            SRCTYPE, SOLARFLUX, SOLARMU, &
                            SKYRAD, UNITS, WAVELEN, WAVENO, SFCTEMP, &
                            SFCTYPE, SFCPARMS, SOURCE, RADIANCE, &
                            NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                            RADOUT)

    ! Computes the radiances for the specified tau's, mu's, and phi's.

    use shd_util

    implicit none
    integer, intent(in) :: NLAY, MAXLEG, NLEGP(NLAY)
    integer, intent(in) :: NPTS, IXP(2, NLAY), IXG(NPTS)
    integer, intent(in) :: ML, MM, NLM, NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(in) :: NZOUT, NMUOUT, NPHIOUT
    logical, intent(in) :: DELTAM
    character(len=1), intent(in) :: SRCTYPE, SFCTYPE, UNITS
    real, intent(in) :: ALBEDOSC(NLAY), LEGENSC(0:Ml, NLAY)
    real, intent(in) :: LEGENP(0:MAXLEG, NLAY)
    real, intent(in) :: HEIGHTG(NPTS), TAUG(NPTS), TAUGSC(NPTS)
    real, intent(in) :: MU(NMU), PHI(NPHI0MAX, NMU), WTDO(NMU, NPHI0MAX)
    real, intent(in) :: SOLARFLUX, SOLARMU, SKYRAD
    real, intent(in) :: SFCPARMS(*), SFCTEMP, WAVELEN, WAVENO(2)
    real, intent(in) :: SOURCE(NLM, NPTS), RADIANCE(NPHI0MAX, NMU, NPTS)
    real, intent(inout) :: ZOUT(NZOUT)
    real, intent(inout) :: RADOUT(*), MUOUT(*), PHIOUT(*)

    integer :: I, J, K, L, IRAD
    real    :: ZDIF, RADDEG, PHIRAD
    real    :: DIRFLUXSFC, SFCPLANCK, BNDRAD
    integer, allocatable :: LZOUT(:)
    real, allocatable :: DIRFLUX(:), SOURCE1(:), RADVEC(:)

    RADDEG = acos(-1.0)/180.0
    allocate (DIRFLUX(NPTS), SOURCE1(NPTS), RADVEC(NPTS), LZOUT(NZOUT))

    ! Find the closest level in height to the output heights
    do I = 1, NZOUT
        ZDIF = 1.0e10
        do L = 1, NPTS
            if (abs(HEIGHTG(L) - ZOUT(I)) < ZDIF) then
                LZOUT(I) = L
                ZDIF = abs(HEIGHTG(L) - ZOUT(I))
            end if
        end do
        ZOUT(I) = HEIGHTG(LZOUT(I))
    end do

    do L = 1, NPTS
        ! Make the perpendicular direct beam flux for all the levels
        DIRFLUX(IXG(L)) = SOLARFLUX/abs(SOLARMU)*exp(-TAUGSC(L)/abs(SOLARMU))
    end do

    if (SFCTYPE /= 'L') then
        if (SRCTYPE /= 'T') then
            DIRFLUXSFC = SOLARFLUX*exp(-TAUGSC(IXP(2, NLAY))/abs(SOLARMU))
            SFCPLANCK = 0.0
        else
            call shd_util__computePlanck(SFCTEMP, UNITS, WAVENO, WAVELEN, SFCPLANCK)
        end if
    end if

    ! Loop over the output radiance directions
    do IRAD = 1, NPTS

        PHIRAD = PHIOUT(IRAD)*RADDEG
        if (MUOUT(IRAD) == 0.0 .or. abs(MUOUT(IRAD)) > 1.0) then
            write (*, *) 'COMPUTE_RADIANCE: Bad mu for radiance', MUOUT(IRAD)
            RADOUT(IRAD) = -999.0
            cycle
        end if

        ! Compute the source function throughout grid for this angle
        call COMPUTE_ONE_SOURCE(NLAY, ALBEDOSC, LEGENSC, &
                                MAXLEG, NLEGP, LEGENP, &
                                NPTS, IXP, IXG, ML, MM, NLM, &
                                MUOUT(IRAD), PHIRAD, SRCTYPE, DELTAM, &
                                SOLARMU, DIRFLUX, SOURCE, SOURCE1)

        ! Get boundary radiances: either top or bottom
        !  Isotropic top boundary or Lambertian bottom boundary can use
        !  the previously computed boundary radiances in RADIANCE,
        !  otherwise, compute the radiance for this angle by integrating
        !  over the stored downwelling radiances.
        if (MUOUT(IRAD) < 0.0) then
            BNDRAD = RADIANCE(1, 1, IXG(IXP(1, 1)))
        else
            if (SFCTYPE == 'L') then
                BNDRAD = RADIANCE(1, NMU/2 + 1, IXG(IXP(2, NLAY)))
            else
                call VARIABLE_BRDF_SURFACE(NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
                                           MUOUT(IRAD), PHIRAD, WAVELEN, &
                                           SRCTYPE, SOLARMU, DIRFLUXSFC, &
                                           SFCTYPE, SFCPARMS, SFCPLANCK, &
                                           RADIANCE(:, :, IXG(IXP(2, NLAY))), BNDRAD)
            end if
        end if

        ! Integrate the source function for this direction to get the radiances
        !   at all the levels (RADVEC), then pull out the desired one.
        call INTEGRATE_SOURCE(NLAY, NPTS, IXP, IXG, TAUGSC, &
                              MUOUT(IRAD), BNDRAD, SOURCE1, RADVEC)
        do I = 1, NZOUT
            RADOUT(IRAD) = RADVEC(LZOUT(I))
        end do
    end do

    deallocate (RADVEC, SOURCE1, DIRFLUX, LZOUT)

    return

end subroutine COMPUTE_RADIANCE

subroutine COMPUTE_ONE_SOURCE(NLAY, ALBEDOSC, LEGENSC, &
                              MAXLEG, NLEGP, LEGENP, &
                              NPTS, IXP, IXG, ML, MM, NLM, &
                              MUO, PHIO, SRCTYPE, DELTAM, &
                              SOLARMU, DIRFLUX, SOURCE, SOURCE1)
    ! Computes the source function (SOURCE1) in the direction (MUO,PHIO)
    ! for all levels (NPTS).  The spherical harmonic source function
    ! series is input in SOURCE.  DIRFLUX is the attenuated direct beam
    ! perpendicular flux.
    ! For a solar source if delta-M then use Nakajima and Tanaka TMS
    ! procedure, replacing delta-M single scattering with single scattering
    ! for unscaled untruncated phase function.

    use shd_scat

    implicit none
    integer, intent(in) :: NLAY, MAXLEG, NLEGP(NLAY), NPTS
    integer, intent(in) :: ML, MM, NLM, IXP(2, NLAY), IXG(NPTS)
    logical, intent(in) :: DELTAM
    character(len=1), intent(in) :: SRCTYPE
    real, intent(in) :: ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(in) :: LEGENP(0:MAXLEG, NLAY)
    real, intent(in) :: MUO, PHIO, SOLARMU, DIRFLUX(NPTS)
    real, intent(in) :: SOURCE(NLM, NPTS)
    real, intent(out) :: SOURCE1(NPTS)
    integer :: I, J, K, L, LAY, M, ME
    double precision  :: F, A, COSSCAT
    real, allocatable :: YLMDIR(:), YLMSUN(:)
    double precision, allocatable :: SUNDIRLEG(:), SINGSCAT(:)

    allocate (YLMDIR(NLM), YLMSUN(NLM), SUNDIRLEG(0:MAXLEG), SINGSCAT(NLAY))

    ! Precompute Ylm's for output direction
    call YLMALL(MUO, PHIO, ML, MM, 1, YLMDIR)

    ! Special case for solar source and Delta-M
    if (SRCTYPE /= 'T' .and. DELTAM) then
        ! Compute the Legendre polynomials for the scattering angle
        !   for the untruncated solar single scattering computation.
        COSSCAT = SOLARMU*MUO + sqrt((1.0d0 - SOLARMU**2)*(1.0d0 - MUO**2))*cos(PHIO)
        call shd_scat__computeLegendre(COSSCAT, MAXLEG, SUNDIRLEG)
        SUNDIRLEG(:) = SUNDIRLEG(:)/(4.0d0*acos(-1.0d0))
        call YLMALL(SOLARMU, 0.0, ML, MM, 1, YLMSUN)

        ! Precompute the single scattering contribution for the original
        ! unscaled phase function (except for the direct beam transmission).
        ! Need to put the inverse of the tau scaling in the source function
        ! because tau is still scaled.
        do LAY = 1, NLAY
            SINGSCAT(LAY) = 0.0
            if (ML + 1 > NLEGP(LAY)) then
                F = 0
            else
                F = LEGENP(ML + 1, LAY)/(2*(ML + 1) + 1)
            end if
            do L = 0, NLEGP(LAY)
                SINGSCAT(LAY) = SINGSCAT(LAY) &
                                + SUNDIRLEG(L)*ALBEDOSC(LAY)*LEGENP(L, LAY)/(1 - F)
            end do
            if (SINGSCAT(LAY) < 0.0) then
                write (*, *) 'COMPUTE_ONE_SOURCE: negative source function for layer: ', LAY
                stop
            end if
        end do

        ! Loop over all the grid points, computing the source function
        !   at the viewing angle from the spherical harmonic source function.
        do LAY = 1, NLAY
            do I = IXP(1, LAY), IXP(2, LAY)
                K = IXG(I)
                ! Sum over the spherical harmonic series of the source function
                SOURCE1(K) = sum(SOURCE(:, K)*YLMDIR(:))
                ! First subtract off the truncated single scattering
                J = 1
                do L = 0, ML
                    ME = min(L, MM)
                    A = DIRFLUX(K)*ALBEDOSC(LAY)*LEGENSC(L, LAY)/(2*L + 1)
                    do M = 0, ME
                        SOURCE1(K) = SOURCE1(K) - A*YLMDIR(J)*YLMSUN(J)
                        J = J + 1
                    end do
                end do
                ! Then add in the single scattering contribution
                SOURCE1(K) = SOURCE1(K) + DIRFLUX(K)*SINGSCAT(LAY)
            end do
        end do

    else ! for thermal sources or no delta-M it is very simple:
        ! Sum over the spherical harmonic series of the source function
        do K = 1, NPTS
            SOURCE1(K) = sum(SOURCE(:, K)*YLMDIR(:))
        end do
    end if

    deallocate (YLMDIR, YLMSUN, SUNDIRLEG, SINGSCAT)
end subroutine COMPUTE_ONE_SOURCE

subroutine INTEGRATE_SOURCE(NLAY, NPTS, IXP, IXG, TAUGSC, &
                            MUO, BNDRAD, SOURCE1, RADVEC)
    ! Integrates the radiative transfer equation through the medium specified
    ! by the optical depth (TAUGSC) and source function (SOURCE1) in direction
    ! cosine MUO.  The radiance at the approprate boundary (top or bottom
    ! depending on the sign of MUO) is input (BNDRAD).  The radiance for
    ! the direction the source function was made at all the adaptive grid
    ! levels is output (RADVEC with indexing corresponding to TAUGSC).
    ! Uses the same quadratic/cubic polynomial source function interpolation
    ! as in PATH_INTEGRATION.
    implicit none
    integer, intent(in) :: NLAY, NPTS, IXP(2, NLAY), IXG(NPTS)
    real, intent(in) :: MUO, BNDRAD
    real, intent(in) :: TAUGSC(NPTS), SOURCE1(NPTS)
    real, intent(out) :: RADVEC(NPTS)
    integer :: LAY, I, I1, I2
    real, allocatable :: SRCINT(:, :)
    double precision :: A(4), B(4), C(4), D(4), E, T, T1, T3, T4, TR

    ! Calculate and store the source function integration coefficients
    allocate (SRCINT(0:4, NPTS))
    do LAY = 1, NLAY
        I1 = IXP(1, LAY)
        I2 = IXP(2, LAY) - 1
        if (I1 == I2) then ! special case for layer with just two points
            T = (TAUGSC(I1 + 1) - TAUGSC(I1))/abs(MUO)
            TR = exp(-T) ! sublayer transmission along ordinate
            if (T < 1.0e-4) then
                D(1) = 0.5*T - (1.0d0/3)*T**2
                D(2) = 0.5*T - (1.0d0/6)*T**2
            else
                D(1) = (1 - TR*(1 + T))/T
                D(2) = (T - (1 - TR))/T
            end if
            SRCINT(0, I1) = TR
            if (MUO < 0.0) then
                SRCINT(1, I1) = D(2)
                SRCINT(2, I1) = D(1)
            else
                SRCINT(1, I1) = D(1)
                SRCINT(2, I1) = D(2)
            end if
        else
            do I = I1, I2
                A(:) = 0.0; B(:) = 0.0; C(:) = 0.0; D(:) = 0.0
                T3 = (TAUGSC(I + 1) - TAUGSC(I))/abs(MUO)
                if (I == I1) then
                    T4 = (TAUGSC(I + 2) - TAUGSC(I))/abs(MUO)
                    E = T3*T4; B(2) = 1/E; C(2) = -(T3 + T4)/E; D(2) = 1.0
                    E = T3*(T3 - T4); B(3) = 1/E; C(3) = -T4/E; D(3) = 0.0
                    E = T4*(T4 - T3); B(4) = 1/E; C(4) = -T3/E; D(4) = 0.0
                else if (I < I2) then
                    T1 = (TAUGSC(I - 1) - TAUGSC(I))/abs(MUO)
                    T4 = (TAUGSC(I + 2) - TAUGSC(I))/abs(MUO)
                    E = T1*(T1 - T3)*(T1 - T4); A(1) = 1/E
                    B(1) = -(T3 + T4)/E; C(1) = T3*T4/E; D(1) = 0.0
                    E = -T1*T3*T4; A(2) = 1/E
                    B(2) = -(T1 + T3 + T4)/E; C(2) = (T1*T3 + T1*T4 + T3*T4)/E; D(2) = 1.0
                    E = T3*(T3 - T1)*(T3 - T4); A(3) = 1/E
                    B(3) = -(T1 + T4)/E; C(3) = T1*T4/E; D(3) = 0.0
                    E = T4*(T4 - T1)*(T4 - T3); A(4) = 1/E
                    B(4) = -(T1 + T3)/E; C(4) = T1*T3/E; D(4) = 0.0
                else
                    T1 = (TAUGSC(I - 1) - TAUGSC(I))/abs(MUO)
                    E = T1*(T1 - T3); B(1) = 1/E; C(1) = -T3/E; D(1) = 0.0
                    E = T1*T3; B(2) = 1/E; C(2) = -(T1 + T3)/E; D(2) = 1.0
                    E = T3*(T3 - T1); B(3) = 1/E; C(3) = -T1/E; D(3) = 0.0
                end if
                T = T3
                TR = exp(-T)
                SRCINT(0, I) = TR
                if (MUO < 0.0) then
                    SRCINT(1:4, I) = A(:)*(6*TR + T**3 - 3*T**2 + 6*T - 6) + D(:)*(1 - TR) &
                                     + B(:)*(-2*TR + T**2 - 2*T + 2) + C(:)*(TR + T - 1)
                else
                    SRCINT(1:4, I) = A(:)*(6 - TR*(T**3 + 3*T**2 + 6*T + 6)) + D(:)*(1 - TR) &
                                     + B(:)*(2 - TR*(T**2 + 2*T + 2)) + C(:)*(1 - TR*(T + 1))
                end if
            end do
        end if
    end do

    if (MUO < 0.0) then
        ! For downwelling direction integrate source function from the top down
        RADVEC(IXP(1, 1)) = BNDRAD
        do LAY = 1, NLAY
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADVEC(I1 + 1) = SRCINT(0, I1)*RADVEC(I1) &
                                 + SRCINT(1, I1)*SOURCE1(IXG(I1)) + SRCINT(2, I1)*SOURCE1(IXG(I1 + 1))
            else
                RADVEC(I1 + 1) = SRCINT(0, I1)*RADVEC(I1) + SRCINT(2, I1)*SOURCE1(IXG(I1)) &
                                 + SRCINT(3, I1)*SOURCE1(IXG(I1 + 1)) + SRCINT(4, I1)*SOURCE1(IXG(I1 + 2))
                do I = I1 + 1, I2 - 1
                    RADVEC(I + 1) = SRCINT(0, I)*RADVEC(I) + SRCINT(1, I)*SOURCE1(IXG(I - 1)) &
                                    + SRCINT(2, I)*SOURCE1(IXG(I)) + SRCINT(3, I)*SOURCE1(IXG(I + 1)) &
                                    + SRCINT(4, I)*SOURCE1(IXG(I + 2))
                end do
                RADVEC(I2 + 1) = SRCINT(0, I2)*RADVEC(I2) + SRCINT(1, I2)*SOURCE1(IXG(I2 - 1)) &
                                 + SRCINT(2, I2)*SOURCE1(IXG(I2)) + SRCINT(3, I2)*SOURCE1(IXG(I2 + 1))
            end if
            if (LAY < NLAY) RADVEC(IXP(1, LAY + 1)) = RADVEC(IXP(2, LAY))
        end do
    else
        ! For upwelling direction integrate source function from the bottom up
        RADVEC(IXP(2, NLAY)) = BNDRAD
        do LAY = NLAY, 1, -1
            I1 = IXP(1, LAY)
            I2 = IXP(2, LAY) - 1
            if (I1 == I2) then
                RADVEC(I1) = SRCINT(0, I1)*RADVEC(I1 + 1) &
                             + SRCINT(1, I1)*SOURCE1(IXG(I1)) + SRCINT(2, I1)*SOURCE1(IXG(I1 + 1))
            else
                RADVEC(I2) = SRCINT(0, I2)*RADVEC(I2 + 1) + SRCINT(1, I2)*SOURCE1(IXG(I2 - 1)) &
                             + SRCINT(2, I2)*SOURCE1(IXG(I2)) + SRCINT(3, I2)*SOURCE1(IXG(I2 + 1))
                do I = I2 - 1, I1 + 1, -1
                    RADVEC(I) = SRCINT(0, I)*RADVEC(I + 1) + SRCINT(1, I)*SOURCE1(IXG(I - 1)) &
                                + SRCINT(2, I)*SOURCE1(IXG(I)) + SRCINT(3, I)*SOURCE1(IXG(I + 1)) &
                                + SRCINT(4, I)*SOURCE1(IXG(I + 2))
                end do
                RADVEC(I1) = SRCINT(0, I1)*RADVEC(I1 + 1) + SRCINT(2, I1)*SOURCE1(IXG(I1)) &
                             + SRCINT(3, I1)*SOURCE1(IXG(I1 + 1)) + SRCINT(4, I1)*SOURCE1(IXG(I1 + 2))
            end if
            if (LAY > 1) RADVEC(IXP(2, LAY - 1)) = RADVEC(IXP(1, LAY))
        end do
    end if
    deallocate (SRCINT)
end subroutine INTEGRATE_SOURCE

subroutine INIT_RADIANCE(NLAY, ML, TAUPSC, ALBEDOSC, LEGENSC, TEMPP, &
                         SRCTYPE, SOLARFLUX, SOLARMU, GNDALB, SFCTEMP, &
                         SKYRAD, UNITS, WAVELEN, WAVENO, NPTS, TAUGSC, IXP, NLM, &
                         SHRADIANCE)
    ! Initializes the spherical harmonic radiance field by solving the
    ! plane-parallel Eddington system to give the first two terms of the
    ! spherical harmonic series at the NPTS sublayer interfaces.

    use shd_init

    implicit none
    integer, intent(in) :: NLAY, ML, NPTS, IXP(2, NLAY), NLM
    real, intent(in) :: TAUPSC(NLAY), ALBEDOSC(NLAY), LEGENSC(0:ML, NLAY)
    real, intent(in) :: TEMPP(NLAY + 1)
    real, intent(in) :: SOLARFLUX, SOLARMU, GNDALB, SFCTEMP
    real, intent(in) :: SKYRAD, WAVELEN, WAVENO(2), TAUGSC(NPTS)
    character(len=1), intent(in) :: SRCTYPE, UNITS
    real, intent(out):: SHRADIANCE(NLM, NPTS)

    integer :: NLAYER, LAY, I, L
    logical, parameter :: DELTAM = .false.
    real    :: PI, C0, C1, GNDEMIS, TAU1, U
    real, allocatable :: OPTDEPTHS(:), ALBEDOS(:), ASYMMETRIES(:)
    real, allocatable :: TEMPS(:), FLUXES(:, :)
    character(len=1) :: SRCTYPE2

    NLAYER = sum(IXP(2, :) - IXP(1, :))
    allocate (OPTDEPTHS(NLAYER), ALBEDOS(NLAYER), ASYMMETRIES(NLAYER))
    allocate (TEMPS(NLAYER + 1), FLUXES(3, NLAYER + 1))

    SRCTYPE2 = SRCTYPE
    if (SRCTYPE == 'B') SRCTYPE2 = 'T'
    PI = acos(-1.0)
    C0 = sqrt(1.0/PI)
    C1 = sqrt(3.0/(4*PI))

    ! Make layer properties for the Eddington routine
    L = 0
    do LAY = 1, NLAY
        TAU1 = TAUGSC(IXP(1, LAY))
        do I = IXP(1, LAY), IXP(2, LAY) - 1
            L = L + 1
            OPTDEPTHS(L) = TAUGSC(I + 1) - TAUGSC(I)
            ALBEDOS(L) = ALBEDOSC(LAY)
            ASYMMETRIES(L) = LEGENSC(1, LAY)/3.0
            if (TAUPSC(LAY) > 0.0) then
                U = (TAUGSC(I) - TAU1)/TAUPSC(LAY)
            else
                U = (I - IXP(1, LAY))/FLOAT(IXP(2, LAY) - IXP(1, LAY))
            end if
            TEMPS(L) = TEMPP(LAY) + U*(TEMPP(LAY + 1) - TEMPP(LAY))
        end do
    end do
    if (L /= NLAYER) stop 'INIT_RADIANCE: initial grid error'
    TEMPS(NLAYER + 1) = TEMPP(NLAY + 1)

    ! Call the Eddington flux routine
    GNDEMIS = 1.0 - GNDALB
    call shd_init__makeFlux(NLAYER, OPTDEPTHS, ALBEDOS, ASYMMETRIES, &
                TEMPS, DELTAM, SRCTYPE2, SOLARFLUX, SOLARMU, &
                SFCTEMP, GNDEMIS, SKYRAD, UNITS, WAVENO, WAVELEN, FLUXES)

    ! Convert fluxes to first two moments of spherical harmonics
    SHRADIANCE(:, :) = 0.0
    L = 0
    do LAY = 1, NLAY
        do I = IXP(1, LAY), IXP(2, LAY)
            L = L + 1
            SHRADIANCE(1, I) = C0*(FLUXES(1, L) + FLUXES(2, L))
            SHRADIANCE(2, I) = C1*(FLUXES(1, L) - FLUXES(2, L))
        end do
        L = L - 1
    end do

    deallocate (OPTDEPTHS, ALBEDOS, ASYMMETRIES, TEMPS, FLUXES)

    return

end subroutine INIT_RADIANCE

subroutine YLMALL(MU, PHI, ML, MM, NCS, P)
    ! This subroutine computes a set of normalized spherical harmonic
    ! functions, P(J), for a particular direction mu,phi.
    ! ML is the maximum meridional mode, MM is the maximum azimuthal mode,
    ! and NCS is the azimuthal mode flag (|NCS|=1 for cosine only, |NCS|=2 for
    ! sines and cosines).  Returns normalized associated Legendre functions
    ! only if NCS<0. The set is returned for triangular truncation:
    ! J = NCS*(L*(L+1))/2 + M+1  for L<=MM
    ! J = (NCS*MM+1)*L-MM*(2+NCS*(MM-1))/2 + M+1  for L>MM
    implicit none
    integer ML, MM, NCS
    real MU, PHI, P(*)
    integer J, L, M, C
    double precision X, Y, A, PMM, PL, PL1, PL2, PHI8

    C = abs(NCS)
    if (C .ne. 1 .and. C .ne. 2) stop 'YLMALL: bad NCS'
    if (MM .gt. ML) stop 'YLMALL: MM greater than LM'
    if (MU*MU .gt. 1.) stop 'YLMALL: |MU| larger than 1'
    X = dble(MU)
    Y = sqrt(1.0d0 - X*X)
    ! Use the stable upward recursion on l, starting from P^m_m
    ! Put in the spherical harmonic normalization as it goes
    PMM = 1.0d0/sqrt(4.0d0*acos(-1.0d0))
    do M = 0, MM
        if (M .gt. 0) PMM = -PMM*Y*sqrt((2*M + 1.0d0)/(2.0d0*M))
        J = C*(M*(M + 1))/2 + M + 1
        P(J) = PMM
        PL2 = PMM
        if (M .lt. ML) then
            if (M + 1 .le. MM) J = C*((M + 1)*(M + 2))/2 + M + 1
            if (M + 1 .gt. MM) J = (C*MM + 1)*(M + 1) - (MM*(2 + C*(MM - 1)))/2 + M + 1
            PL1 = sqrt(2*M + 3.0d0)*X*PMM
            P(J) = PL1
        end if
        do L = M + 1, ML - 1
            if (L + 1 .le. MM) J = C*((L + 1)*(L + 2))/2 + M + 1
            if (L + 1 .gt. MM) J = (C*MM + 1)*(L + 1) - (MM*(2 + C*(MM - 1)))/2 + M + 1
            A = 1.0d0/((L + M + 1.d0)*(L - M + 1.d0))
            PL = sqrt((2*L + 1)*A*(2*L + 3))*X*PL1 &
                 - sqrt((2*L + 3)*A*(L + M)*(L - M)/(2*L - 1.))*PL2
            P(J) = PL
            PL2 = PL1
            PL1 = PL
        end do
        if (M .eq. 0) PMM = PMM*sqrt(2.0d0)
    end do
    ! If there are M<0 terms then fill them in
    if (C .eq. 2) then
        do L = 0, ML
            do M = 1, min(L, MM)
                if (L .le. MM) J = L*(L + 1) + M + 1
                if (L .gt. MM) J = MM*(2*L - MM) + L + M + 1
                P(J - 2*M) = P(J)
            end do
        end do
    end if
    ! Put in the azimuthal dependence
    PHI8 = dble(PHI)
    if (NCS .gt. 0) then
        do M = (1 - NCS)*MM, MM
            if (M .lt. 0) then
                A = sin(-M*PHI8)
            else if (M .gt. 0) then
                A = cos(M*PHI8)
            else
                A = 1.0d0
            end if
            do L = abs(M), ML
                if (L .le. MM) then
                    J = C*(L*(L + 1))/2 + M + 1
                else
                    J = (C*MM + 1)*L - (MM*(2 + C*(MM - 1)))/2 + M + 1
                end if
                P(J) = A*P(J)
            end do
        end do
    end if
end subroutine YLMALL

subroutine ADD_IN_KDIST(NLAY, HEIGHTP, TAUP, ALBEDOP, &
                        NZCKD, ZCKD, GASABS, TAUPM, ALBEDOPM)
    ! Adds the gaseous absorption in with the particle optical properties.
    ! The gaseous absorption optical depth is found for each of the NLAY
    ! optical property layers by interpolated the absorption linearly
    ! between the levels in ZCKD.
    integer, intent(in) :: NLAY, NZCKD
    real, intent(in) :: HEIGHTP(NLAY + 1), TAUP(NLAY), ALBEDOP(NLAY)
    real, intent(in) :: ZCKD(NZCKD), GASABS(NZCKD)
    real, intent(out) :: TAUPM(NLAY), ALBEDOPM(NLAY)
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
            !print *, 'add_kdist: ',LAY,GASTAU,TAUPM(LAY),ALBEDOPM(LAY)
        end do
    end if

    return
end subroutine ADD_IN_KDIST

subroutine SUM_KDIST_OUTPUT(IG, WT, NLAY, TAUP, IXP, NPTS, TAUG, &
                            HEIGHTP, HEIGHTG, FLUXUPM, FLUXUP, &
                            FLUXDNM, FLUXDN, MEANRADM, MEANRAD, &
                            NZOUT, NMUOUT, NPHIOUT, RADOUTM, RADOUT)
    ! Sums the output arrays over the k-distribution.  Only sums the output
    ! at the input layer interfaces because the adaptive grid changes for
    ! each k.  Also fixes up the tau grid and height grid accordingly.
    integer, intent(in) :: IG, NLAY, NZOUT, NMUOUT, NPHIOUT
    integer, intent(inout) :: NPTS, IXP(2, NLAY)
    real, intent(in) :: WT, TAUP(NLAY), HEIGHTP(NLAY + 1)
    real, intent(inout) :: TAUG(*), HEIGHTG(*)
    real, intent(in) :: FLUXUPM(*), FLUXDNM(*), MEANRADM(*)
    real, intent(in) :: RADOUTM(NZOUT, NMUOUT, NPHIOUT)
    real, intent(out) :: FLUXUP(*), FLUXDN(*), MEANRAD(*)
    real, intent(out) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)

    if (IG == 1) then
        FLUXUP(1:NLAY + 1) = 0.0
        FLUXDN(1:NLAY + 1) = 0.0
        MEANRAD(1:NLAY + 1) = 0.0
        RADOUT(:, :, :) = 0.0
    end if

    do L = 1, NLAY
        FLUXUP(L) = FLUXUP(L) + WT*FLUXUPM(IXP(1, L))
        FLUXDN(L) = FLUXDN(L) + WT*FLUXDNM(IXP(1, L))
        MEANRAD(L) = MEANRAD(L) + WT*MEANRADM(IXP(1, L))
    end do
    FLUXUP(NLAY + 1) = FLUXUP(NLAY + 1) + WT*FLUXUPM(IXP(2, NLAY))
    FLUXDN(NLAY + 1) = FLUXDN(NLAY + 1) + WT*FLUXDNM(IXP(2, NLAY))
    MEANRAD(NLAY + 1) = MEANRAD(NLAY + 1) + WT*MEANRADM(IXP(2, NLAY))
    RADOUT(:, :, :) = RADOUT(:, :, :) + WT*RADOUTM(:, :, :)

    HEIGHTG(1:NLAY + 1) = HEIGHTP(1:NLAY + 1)
    TAUG(1) = 0.0
    do L = 1, NLAY
        TAUG(L + 1) = TAUG(L) + TAUP(L)
        IXP(1, L) = L
        IXP(2, L) = L + 1
    end do
    NPTS = NLAY + 1
end subroutine SUM_KDIST_OUTPUT

subroutine USER_INPUT(PROPFILE, CKDFILE, KDIST, NMU, NPHI, &
                      SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD, &
                      SFCTEMP, SFCTYPE, SFCPARMS, WAVELEN, WAVENO, &
                      ACCELFLAG, SOLACC, MAXITER, SPLITACC, &
                      FLUXFILE, FLUXPARM, RADFILE, MAXPAR, NZOUT, ZOUT, &
                      NMUOUT, MUOUT, NPHIOUT, PHIOUT)
    !   Obtains the input parameters for the program by writing prompts
    ! to and reading responses from the terminal.  Echos the inputs to
    ! make a useful log file when running non-interactively.  See the overall
    ! program documentation for the list of input parameters.
    implicit none
    integer, intent(in)  :: MAXPAR
    integer, intent(out) :: NMU, NPHI
    integer, intent(out) :: MAXITER, FLUXPARM, NZOUT, NMUOUT, NPHIOUT
    logical, intent(out) :: KDIST, ACCELFLAG
    real, intent(out) :: SOLARFLUX, SOLARMU
    real, intent(out) :: SFCTEMP, SKYRAD, WAVELEN, WAVENO(2)
    real, intent(out) :: SOLACC, SPLITACC
    real, intent(out) :: SFCPARMS(MAXPAR)
    real, intent(out) :: ZOUT(MAXPAR), MUOUT(MAXPAR), PHIOUT(MAXPAR)
    character(len=1), intent(out) :: SRCTYPE, SFCTYPE
    character(len=80), intent(out) :: PROPFILE, CKDFILE, FLUXFILE, RADFILE

    write (*, *)
    write (*, '(1X,A)') 'Atmospheric optical properties file name'
    read (*, '(A)') PROPFILE
    write (*, *) PROPFILE

    write (*, '(1X,A)') 'Correlated k-distribution file name (or NONE)'
    read (*, '(A)') CKDFILE
    write (*, *) CKDFILE
    KDIST = CKDFILE(1:2) .ne. 'NO'

    write (*, '(1X,A)') 'Number of discrete ordinates in mu and phi'
    read (*, *) NMU, NPHI
    write (*, *) NMU, NPHI

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
        SOLARMU = -abs(SOLARMU)
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

    write (*, '(1X,A)') 'Surface type (L-Lambertian, R-RPV vegetation BRDF)'
    read (*, '(A1)') SFCTYPE
    write (*, *) SFCTYPE
    if (SFCTYPE == 'L') then
        write (*, '(1X,A)') 'Surface albedo'
        read (*, *) SFCPARMS(1)
        write (*, *) SFCPARMS(1)
    else if (SFCTYPE == 'R') then
        write (*, '(1X,A)') 'RPV parameters (rho0, k, Theta)'
        read (*, *) SFCPARMS(1:3)
        write (*, *) SFCPARMS(1:3)
    else
        write (*, *) 'Unknown surface type'
        stop
    end if

    write (*, '(1X,A)') 'Wavelength (micron) (<=0 for wavenumber range)'
    read (*, *) WAVELEN
    write (*, *) WAVELEN
    if (KDIST .and. WAVELEN > 0.0) then
        write (*, *) 'Must use wavenumber range for k-distribution'
        stop
    end if
    if (WAVELEN <= 0.0) then
        write (*, '(1X,A)') 'Wavenumber range (cm^-1)'
        read (*, *) WAVENO(1:2)
        write (*, *) WAVENO(1:2)
    end if

    write (*, '(1X,A)') 'Layer splitting accuracy'
    read (*, *) SPLITACC
    write (*, *) SPLITACC

    write (*, '(1X,A,A)') 'Solution accuracy'
    read (*, *) SOLACC
    write (*, *) SOLACC

    write (*, '(1X,A,A)') 'Acceleration flag (T or F) and max iterations'
    read (*, *) ACCELFLAG, MAXITER
    write (*, *) ACCELFLAG, MAXITER

    write (*, '(1X,A,A)') 'Flux output file name'
    read (*, '(A)') FLUXFILE
    write (*, '(A)') FLUXFILE
    write (*, '(1X,A)') 'Fluxes output at input layer boundaries or adaptive levels (1,2)'
    read (*, *) FLUXPARM
    write (*, *) FLUXPARM

    write (*, '(1X,A,A)') 'Radiance output file name'
    read (*, '(A)') RADFILE
    write (*, '(A)') RADFILE

    write (*, *) 'Number of heights for output radiance'
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
    !   Reads parts of the property file to get the maximum array sizes
    ! needed for allocatable arrays.
    implicit none
    integer, intent(out) :: NLAY, MAXLEG
    character(len=*), intent(in) :: PROPFILE
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
    integer, intent(in)  :: NLAY, MAXLEG
    integer, intent(out) :: NLEGP(NLAY)
    real, intent(out) :: TAUP(NLAY), ALBEDOP(NLAY), LEGENP(0:MAXLEG, NLAY)
    real, intent(out) :: HEIGHTP(NLAY + 1), TEMPP(NLAY + 1)
    character(len=*), intent(in) :: PROPFILE
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

subroutine OUTPUT_RESULTS(NLAY, NPTS, NMU, NPHI, NG, &
                          PROPFILE, CKDFILE, SRCTYPE, SOLARFLUX, SOLARMU, &
                          SKYRAD, SFCTYPE, SFCTEMP, SFCPARMS, WAVELEN, WAVENO, &
                          SPLITACC, SOLACC, MAXITER, ITER, &
                          IXP, HEIGHTG, FLUXUP, FLUXDN, MEANRAD, FLUXPARM, &
                          NZOUT, ZOUT, NMUOUT, MUOUT, NPHIOUT, PHIOUT, &
                          RADOUT, OUTTYPE, OUTFILE)
    !   Writes the desired type of output file from the output fields
    ! with a header giving the input parameters.
    ! There are two types (OUTTYPE) of output: 'R' - for radiance,
    ! 'F' - for hemispheric/actinic flux.
    ! See the overall program documentation for output parameters.
    implicit none
    integer, intent(in) :: NLAY, NPTS, NMU, NPHI, NG
    integer, intent(in) :: NZOUT, NMUOUT, NPHIOUT, FLUXPARM
    integer, intent(in) :: MAXITER, ITER, IXP(2, NLAY)
    real, intent(in) :: SOLARFLUX, SOLARMU, SFCTEMP, SFCPARMS(*), SKYRAD
    real, intent(in) :: SOLACC, SPLITACC, WAVELEN, WAVENO(2)
    real, intent(in) :: HEIGHTG(NPTS)
    real, intent(in) :: FLUXUP(NPTS), FLUXDN(NPTS), MEANRAD(NPTS)
    real, intent(in) :: ZOUT(NZOUT), MUOUT(NMUOUT), PHIOUT(NPHIOUT)
    real, intent(in) :: RADOUT(NZOUT, NMUOUT, NPHIOUT)
    character(len=1), intent(in) :: SRCTYPE, SFCTYPE, OUTTYPE
    character(len=*), intent(in) :: PROPFILE, CKDFILE, OUTFILE
    integer :: I, J, K, L, N
    real    :: FOURPI
    character(len=32) :: SOURCENAME, UNITSNAME, OUTNAME, SFCNAME

    FOURPI = 4*acos(-1.0)
    if (SRCTYPE == 'S') then
        SOURCENAME = 'SOLAR'
    else if (SRCTYPE == 'T') then
        SOURCENAME = 'THERMAL'
    else if (SRCTYPE == 'B') then
        SOURCENAME = 'SOLAR/THERMAL'
    end if

    if (SFCTYPE == 'L') then
        SFCNAME = 'LAMBERTIAN'
    else if (SFCTYPE == 'R') then
        SFCNAME = 'RPV-BRDF'
    end if

    if (WAVELEN > 0.0) then
        if (OUTTYPE == 'R') then
            UNITSNAME = 'WATTS/(M^2 MICRON STER)'
        else
            UNITSNAME = 'WATTS/(M^2 MICRON)'
        end if
    else
        if (OUTTYPE == 'R') then
            UNITSNAME = 'WATTS/(M^2 STER)'
        else
            UNITSNAME = 'WATTS/(M^2)'
        end if
    end if
    if (OUTTYPE == 'R') OUTNAME = 'RADIANCE'
    if (OUTTYPE == 'F') OUTNAME = 'FLUX'

    open (UNIT=2, FILE=OUTFILE, STATUS='UNKNOWN')
    write (2, '(A)') '! SHDOM Plane-parallel Radiative Transfer Output'
    write (2, '(2(A,I3))') '!  NMU=', NMU, '  NPHI=', NPHI
    write (2, '(2(A,I4))') '!  NLAYER(input)=', NLAY, '    NPTS=', NPTS
    write (2, '(A,A32)') '!  PROPERTY_FILE=', PROPFILE
    write (2, '(A,A32,A,I2)') '!  CORRELATED_K-DIST_FILE=', CKDFILE, &
        '   NUM_G=', NG
    write (2, '(A,A14)') '!  SOURCE_TYPE=', SOURCENAME
    write (2, '(A,A20)') '!  SURFACE_TYPE=', SFCNAME
    if (SFCTYPE == 'L') then
        write (2, '(A,F9.7)') '!  SURFACE_ALBEDO=', SFCPARMS(1)
    else if (SFCTYPE == 'R') then
        write (2, '(3(A,F8.5))') '!  SURFACE_RHO0=', SFCPARMS(1), &
            '   SURFACE_K=', SFCPARMS(2), '   SURFACE_THETA=', SFCPARMS(3)
    end if
    write (2, '(A,F8.3,A,E12.5)') &
        '!  SURFACE_TEMP=', SFCTEMP, '  SKY_RAD=', SKYRAD

    if (SRCTYPE .eq. 'S' .or. SRCTYPE .eq. 'B') then
        write (2, '(A,E13.6,A,F10.7)') &
            '!  SOLAR_FLUX=', SOLARFLUX, '   SOLAR_MU=', SOLARMU
    end if
    if (WAVELEN > 0.0) then
        write (2, '(A,A24,A,F10.2)') '!  UNITS=', UNITSNAME, &
            '   WAVELENGTH(MICRON)=', WAVELEN
    else
        write (2, '(A,A18,A,2F10.3)') '!  UNITS=', UNITSNAME, &
            '   WAVENUMBER_RANGE(CM^-1)=', WAVENO(1:2)
    end if
    write (2, '(A,E10.3)') '!  SPLITTING_ACCURACY=', SPLITACC
    write (2, '(A,E10.3)') '!  SOLUTION_ACCURACY=', SOLACC
    write (2, '(2(A,I4))') '!  MAXIMUM_ITERATIONS=', MAXITER, &
        '   NUMBER_ITERATIONS=', ITER
    write (2, '(A,A20)') '!  OUTPUT_TYPE=', OUTNAME

    if (OUTTYPE .eq. 'R') then
        ! Radiance output
        write (2, '(3(A,I3))') '!   NZOUT=', NZOUT, &
            '  NMUOUT=', NMUOUT, '  NPHIOUT=', NPHIOUT
        write (2, '(A)') '!    Z       MU     PHI     RADIANCE'
        do I = 1, NZOUT
            do J = 1, NMUOUT
                do K = 1, NPHIOUT
                    write (2, '(1X,F7.3,1X,F8.5,1X,F6.2,1X,E13.6)') &
                        ZOUT(I), MUOUT(J), PHIOUT(K), RADOUT(I, J, K)
                end do
            end do
        end do

    else if (OUTTYPE .eq. 'F') then
        ! Hemispheric flux output
        !   Format 1 is fluxes at every property grid point
        !   Format 2 is fluxes at every grid point
        write (2, '(A)') '!    Z         F_UP        F_DOWN       F_ACTINIC'
        if (FLUXPARM == 1) then
            do L = 1, NLAY
                I = IXP(1, L)
                write (2, '(1X,F7.3,3(2X,E12.5))') &
                    HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
            end do
            I = IXP(2, NLAY)
            write (2, '(1X,F7.3,3(2X,E12.5))') &
                HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
        else
            do L = 1, NLAY
                do I = IXP(1, L), IXP(2, L) - 1
                    write (2, '(1X,F7.3,3(2X,E12.5))') &
                        HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
                end do
            end do
            I = IXP(2, NLAY)
            write (2, '(1X,F7.3,3(2X,E12.5))') &
                HEIGHTG(I), FLUXUP(I), FLUXDN(I), FOURPI*MEANRAD(I)
        end if

    end if
    close (2)
end subroutine OUTPUT_RESULTS
