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
!╰────────────────────────────────────────────────────────────────────────────╯!

program shdom

    use globals_const
    use globals_cfg
    use globals_prp

    use shd_init
    use shd_src
    use shd_shdo
    use shd_solv
    use shd_rad
    use shd_inp ! shd_io__
    use shd_out ! shd_io__
    use shd_mpi
    use shd_util

    implicit none

    integer, parameter :: R__ = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__ = selected_int_kind(9)  ! default precision for integer

    ! Number of output files and inputs parameter per file
    integer(kind=I__), parameter :: MAXOUT = 20, MAXPAR = 100
    real(kind=R__)   :: OUTPARMS(MAXPAR, MAXOUT)
    character(len=1) :: OUTTYPES(MAXOUT)
    character(len=FNAMEMAXLEN_) :: OUTFILES(MAXOUT)

    integer(kind=2)  , allocatable :: IPHASE(:)
    integer(kind=2)  , allocatable :: IPHASEPT(:)
    integer(kind=2)  , allocatable :: CELLFLAGS(:)
    integer(kind=I__), allocatable :: GRIDPTR(:, :), NEIGHPTR(:, :), TREEPTR(:, :)
    integer(kind=I__), allocatable :: RSHPTR(:), SHPTR(:), OSHPTR(:)
    integer(kind=I__), allocatable :: WORK1(:)
    integer(kind=I__), allocatable :: NPHI0(:)
    integer(kind=I__), allocatable :: BCPTR(:)

    real(kind=R__), allocatable :: SOLFLUX(:), DELG(:), KABS(:, :)
    real(kind=R__), allocatable :: XGRID(:), YGRID(:), ZGRID(:)
    real(kind=R__), allocatable :: TEMP(:), PLANCK(:), EXTINCT(:), ALBEDO(:)
    real(kind=R__), allocatable :: LEGEN(:)
    real(kind=R__), allocatable :: GRIDPOS(:, :)
    real(kind=R__), allocatable :: RADIANCE(:, :), SOURCE(:, :), DELSOURCE(:, :)
    real(kind=R__), allocatable :: WORK(:), WORK2(:)
    real(kind=R__), allocatable :: FLUXES(:, :), DIRFLUX(:), DIRFLUX_TMP(:)
    real(kind=R__), allocatable :: MU(:), WTDO(:), PHI(:)
    real(kind=R__), allocatable :: BCRAD(:, :)
    real(kind=R__), allocatable :: SFCPARMS(:)
    real(kind=R__), allocatable :: SFCGRIDPARMS(:)
    real(kind=R__), allocatable :: SENPARMS(:)
    real(kind=R__), allocatable :: SUMFLUXES(:, :), SUMDIRFLUX(:), SUMFLUXDIV(:)
    real(kind=R__), allocatable :: SUMSHTERMS(:), SOURCE1OUT(:)
    real(kind=R__), allocatable :: SUMRADOUT(:), SUMVISOUT(:)
    real(kind=R__), allocatable :: TEMPPT(:), EXTINCTPT(:), ALBEDOPT(:)
    real(kind=R__), allocatable :: LEGENPT(:)
    real(kind=R__), allocatable :: ALLFLUXES(:, :, :, :)
    real(kind=R__), allocatable :: ALLFLUXDIV(:, :, :), ALLSHTERMS(:, :, :, :)
    real(kind=R__), allocatable :: XALLGRID(:), YALLGRID(:)

    integer(kind=I__) :: MAXPG, MAXIG, MAXIC, MAXIV, MAXIDO, MAXNZ, MAXNG
    integer(kind=I__) :: MAXLEG, MAXPGL, MAXIGL
    integer(kind=I__) :: MAXNBC, MAXBCRAD
    integer(kind=I__) :: MAXSFCPTS, MAXSFCPARS, BIG_ARRAYS

    real(kind=R__) :: MEMWORD, WANTMEM, NMAXIG
    real(kind=R__) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real(kind=R__) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real(kind=R__) :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer(kind=I__) :: NUMPROC, NSTOKES, NX, NY, NZ, NX1, NY1
    integer(kind=I__) :: NXT, NYT, NPXT, NPYT, MAXPGT, MAXPGLT
    integer(kind=I__) :: NXSFC, NYSFC, NSFCPAR, NTOPPTS, NBOTPTS
    integer(kind=I__) :: NSENDATA, NSENPAR
    integer(kind=I__) :: NG, SIDE
    integer(kind=I__) :: ML, MM, NLM, NSTLEG, NLEG, NANG, NMU, NPHI, NPHI0MAX
    integer(kind=I__) :: NBPTS, NPTS, NCELLS, NSH, NBCELLS, OLDNPTS
    integer(kind=I__) :: MAXITER, ITER, TOTITER, IG, I, NUMOUT, BCFLAG, IPFLAG
    integer(kind=I__) :: NRAD, IRAD, NXOUT, NYOUT, NVIS, IVIS, NSCAN, NPIX, NSHOUT
    integer(kind=I__) :: NCELLSTOT, NPTSTOT, NSHTOT
    integer(kind=I__) :: COUNT_START, COUNT_RATE, COUNT_MAX

    real(kind=R__) :: SOLARFLUX0, SOLARFLUXCKD, SOLARFLUX, SOLARMU, SOLARAZ
    real(kind=R__) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real(kind=R__) :: SOLACC, SOLCRIT, SPLITACC, SHACC
    real(kind=R__) :: DELXSFC, DELYSFC, MAXASYM, REDUCE_RATIO
    real(kind=R__) :: MUOUT, PHIOUT
    real(kind=R__) :: UNIFZLEV, XO, YO, ZO, DIRPATH
    real(kind=R__) :: DELX0, DELY0

    character(len=1) :: SRCTYPE, UNITS
    character(len=1) :: GRIDTYPE, PROPTYPE
    character(len=2) :: SFCTYPE
    character(len=FNAMEMAXLEN_) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
    character(len=FNAMEMAXLEN_) :: SENFILE, OutFileNC, RUNNAME
    character(len=FNAMEMAXLEN_) :: MESSAGE

    logical :: KDIST, INRADFLAG, NEWGRIDFLAG, DELTAM, ACCELFLAG
    logical :: BASEOUT, MASTERPROC, VALIDBEAM
    logical :: RADIANCE_OUT, VISUAL_OUT
    logical :: FLUX_OUT, FLUXDIV_OUT, SH_OUT, SOURCE_OUT
    logical :: MPIX, MPIY, MPIMODE
    logical :: NMLMODE = .true.


    ! start timing
    call system_clock(COUNT_START, COUNT_RATE, COUNT_MAX)


    ! Start up MPI for multiple processors if available, and
    !   return flag for master process (procnum=0)
    call shd_mpi__start(MASTERPROC, NUMPROC)


    ! Read input from namelist or user input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! In MPI, the master processor performs the main computation,
    !   and the results are then broadcast to the child processors.
    if (MASTERPROC) then
        if (NMLMODE) then
            call shd_io__nmlInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                INSAVEFILE, OUTSAVEFILE, SENFILE,&
                NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
                IPFLAG, KDIST, DELTAM, GRIDTYPE,&
                SRCTYPE, SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD,&
                GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
                ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
                MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
                OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, VERBOSE_)
        else
            call shd_io__userInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
                INSAVEFILE, OUTSAVEFILE, SENFILE,&
                NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
                IPFLAG, KDIST, DELTAM, GRIDTYPE,&
                SRCTYPE, SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD,&
                GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
                ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
                MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES,&
                OutFileNC, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, VERBOSE_)
        end if
    end if

    call shd_mpi__bcastInput(RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
        &INSAVEFILE, OUTSAVEFILE, SENFILE,&
        &NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
        &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
        &SRCTYPE, SOLARFLUX0, SOLARMU, SOLARAZ, SKYRAD,&
        &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
        &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
        &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! account for initialization
    MEMWORD = 2000.0 + MAXOUT * (MAXPAR + 0.25 + (FNAMEMAXLEN_/4.0))


    ! Post-process parameters from input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! set calculation accuracy (1.0e-4 to 1.0e-8)
    RACC_ = max(REPS_, min(SOLACC*0.1, 1.0e-4))

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    NSTOKES = max(1, min(4, NSTOKES))
    if (NSTOKES .eq. 2) then
        NSTOKES = 3
        if (VERBOSE_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NSTOKES cannot be 2, reset NSTOKES to 3.'
    end if

    ! NSTLEG is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (NSTOKES .eq. 1) then
        NSTLEG = 1
    else
        NSTLEG = 6
    end if

    ! Make ML and MM from NMU and NPHI (for discrete ordinate)
    !   NMU : number of discrete ordinate zenith angles in both hemispheres
    !   NPHI: number of discrete ordinate azimuth angles in 2*pi
    NMU  = max(2, 2*int((NMU + 1)/2)) ! force NMU an even number
    NPHI = max(1, NPHI)

    ! ML  : The maximum L of the spherical harmonic truncation
    ! MM  : The maximum M of the spherical harmonic truncation
    ML = NMU - 1
    MM = max(0, int(NPHI/2) - 1)

    ! Compute NLM, NPHI0MAX.
    !   NLM: The number of spherical harmonic terms
    !   NPHI0MAX: The maximum number of azimuth angles actually used
    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    NPHI0MAX = NPHI

    ! NLEG: The number of Legendre series phase function terms used in the solution procedure
    NLEG = ML
    if (DELTAM) NLEG = ML + 1

    ! for thermal, reset solar azimuth angle
    if (SRCTYPE .eq. 'T') SOLARAZ = 0.0d0

    ! angle grids for DO
    allocate (MU(NMU), WTDO(NMU*NPHI), PHI(NMU*NPHI), NPHI0(NMU))
    MEMWORD = MEMWORD + 2.0*NMU*(1.0 + NPHI)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read PROPFILE (property file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (MASTERPROC) then
        ! Find out the size of the input property file grid
        call shd_io__readPrpSize(PROPFILE, NSTLEG, NLEG, NPX_, NPY_, NPZ_, DELX_, DELY_,&
            &NUMPHASE_, MAXLEG, MAXPGL)
    end if
    call shd_mpi__bcastPrpSize(NPX_, NPY_, NPZ_, DELX_, DELY_,&
        &NUMPHASE_, MAXLEG, MAXPGL)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read SFCFILE (surface file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (SFCFILE(1:4) .eq. 'NONE') then
        SFCTYPE = 'FL'
        NXSFC = 1
        NYSFC = 1
        DELXSFC = DELX_
        DELYSFC = DELY_
    else
        if (MASTERPROC) then
            ! find out the size of the surface file grid
            call shd_io__readSfcSize(SFCFILE, NXSFC, NYSFC, DELXSFC, DELYSFC, MAXSFCPTS, MAXSFCPARS)
        end if
        ! Broadcast the surface array sizes
        call shd_mpi__bcastSfcSize(NXSFC, NYSFC, DELXSFC, DELYSFC, MAXSFCPTS, MAXSFCPARS)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Process grids
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! X, Y, Z grids
    NX = max(1, NX); NY = max(1, NY); NZ = max(2, NZ)

    ! If single plane or column then force IP mode
    !   if   3D IPFLAG=0, bit representation '0000'
    !   if NX=1 IPFLAG=1, bit representation '0001'
    !   if NY=1 IPFLAG=2, bit representation '0010'
    !   if  IPA IPFLAG=3, bit representation '0011'
    if (NX .eq. 1) then
        IPFLAG = ior(IPFLAG, 1) ! '0000'(3D,   0) + '0001'(NX=1, 1) = '0001' (NX=1, 1)
        if (VERBOSE_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NX=1, reset IPFLAG to 1.'
    end if
    if (NY .eq. 1) then
        IPFLAG = ior(IPFLAG, 2) ! '0001'(NX=1, 1) + '0010'(NY=1, 2) = '0011' (IPA,  3)
        if (VERBOSE_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NY=1, reset IPFLAG to 2 (or 3).'
    end if

    ! Set the full domain grid sizes (variables end in T)
    NXT  = NX  ; NYT  = NY
    NPXT = NPX_; NPYT = NPY_

    if (NPXT .eq. 1) then
        DELX_ = NXSFC*DELXSFC ! overwrite DELX_ if only 1 grid in X direction
        if (VERBOSE_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NPXT=1, reset DELX_ to match surface domain.'
    end if
    if (NPYT .eq. 1) then
        DELY_ = NYSFC*DELYSFC ! overwrite DELY_ if only 1 grid in Y direction
        if (VERBOSE_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NPYT=1, reset DELY_ to match surface domain.'
    end if

    XDOMAIN_ = NPXT*DELX_
    YDOMAIN_ = NPYT*DELY_

    ! Choose the option to compute higher order radiance SH terms
    HIGHORDERRAD_ = .false.
    if ((NPXT .eq. 1) .and. (NPYT .eq. 1) .and.&
        &(NXSFC .le. 32) .and. (NXSFC .le. 32)) then
        HIGHORDERRAD_ = .true.
        if (VERBOSE_ .ge. 4) write (6, '(A)') 'Warning [shdom]: 1D calculation, set HIGHORDERRAD_ to .TRUE.'
    end if

    NSHOUT = 4
    if (HIGHORDERRAD_) NSHOUT = 5

    ! x, y, z grids
    allocate (XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ))
    MEMWORD = MEMWORD + (NX+1) + (NY+1) + (NZ)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! If MPI, divide domain into mini-domains
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and NX,NY to the subdomain size, and output
    !   the X and Y starting position of this subdomain. Set bits
    !   2 and 3 in BCFLAG to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    call shd_mpi__map(BCFLAG, NPX_, NPY_, NX, NY, NZ, DELX_, DELY_, PROPFILE,&
        &XSTART_, YSTART_, RUNNAME)

    ! BCFLAG bit-representation '...01xx' MPIX
    if (btest(BCFLAG, 2)) then
        MPIX = .true.
    else
        MPIX = .false.
    end if

    ! BCFLAG bit-representation '...1xxx' MPIY
    if (btest(BCFLAG, 3)) then
        MPIY = .true.
    else
        MPIY = .false.
    end if

    if (MPIX .or. MPIY) then
        MPIMODE = .true.
    else
        MPIMODE = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Base grids setup after MPI map
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Choose type of initialization/output for k-distributions:
    !   True means initiatize with and output only the base grid.
    BASEOUT = .false.

    ! Set up base grid point actual size (NX1xNY1xNZ)
    !   btest(IPFLAG, 0) checks whether IPFLAG=1 or IPFLAG=0 ('...0001' [IPA in X] or '...0000' [3D in X])
    !   btest(BCFLAG, 0) checks whether BCFLAG=1 or BCFLAG=0 ('...0001' [open in X] or '...0000' [periodic in X])
    NX1 = NX + 1
    if (btest(IPFLAG, 0) .or. btest(BCFLAG, 0) .or. MPIX) NX1 = NX

    ! btest(IPFLAG, 1) checks whether IPFLAG=2 [or 3] or IPFLAG=0 [or 1]
    !   ('...0010' [or '...0011', IPA in Y] or '...0000' [or '...0001', 3D in Y])
    ! btest(BCFLAG, 1) checks whether BCFLAG=2 [or 3] or BCFLAG=0 [or 1]
    !   ('...0010' [or '...0011', open in Y] or '...0000' [or '...0001', periodic in Y])
    NY1 = NY + 1
    if (btest(IPFLAG, 1) .or. btest(BCFLAG, 1) .or. MPIY) NY1 = NY

    ! Number of base points
    NBPTS = NX1*NY1*NZ

    ! Calculate the number of base grid cells depending on the BCFLAG
    NBCELLS = (NZ - 1)*(NX + ibits(BCFLAG, 0, 1) - ibits(BCFLAG, 2, 1))&
        *(NY + ibits(BCFLAG, 1, 1) - ibits(BCFLAG, 3, 1))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read PRPFILE and distribute data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Determine the size of the property grid arrays so can allocate
    MAXNZ = NPZ_
    MAXPG = NPX_*NPY_*NPZ_
    if (NUMPHASE_ .eq. 0) MAXPGL = MAXPG*(MAXLEG + 1)*NSTLEG

    allocate (ZLEVELS_(NPZ_))
    allocate (TEMPP_(MAXPG), EXTINCTP_(MAXPG), EXTDIRP_(MAXPG), ALBEDOP_(MAXPG))
    allocate (IPHASEP_(MAXPG), LEGENP_(MAXPGL))
    MEMWORD = MEMWORD + NPZ_ + 6.0*MAXPG + MAXPGL

    if (MPIMODE) then

        ! Read the properties of the medium to the master processor and
        !   distribute to all the processors
        if (MASTERPROC) then

            MAXPGT = NPXT*NPYT*NPZ_
            MAXPGLT = MAXPGL
            if (NUMPHASE_ .eq. 0) MAXPGLT = MAXPGT*(MAXLEG + 1)*NSTLEG

            allocate (TEMPPT(MAXPGT), EXTINCTPT(MAXPGT), ALBEDOPT(MAXPGT))
            allocate (IPHASEPT(MAXPGT), LEGENPT(MAXPGLT))

            ! Only use DELX and DELY from initial read, DELX0 and DELY0 are discarded
            call shd_io__readPrp(PROPFILE, NPXT, NPYT, NPZ_,&
                &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPGT, MAXPGLT, DELTAM,&
                &PROPTYPE, DELX0, DELY0, ZLEVELS_, MAXASYM,&
                &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE_, IPHASEPT)

            if (NUMPHASE_ .eq. 0) then
                stop 'Error [shdom]: Standard property file (with phase function Legendre series at every point)&
                    & is not allowed with multiple processors.'
            end if

            ! Check values of medium properties
            call shd_io__checkPrp(NPXT, NPYT, NPZ_, NSTLEG, NLEG,&
                &DELX_, DELY_, ZLEVELS_,&
                &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, NUMPHASE_, IPHASEPT)

        end if

        ! distribute file grid data to working grids (different CPU processes different grids)
        call shd_mpi__distPrpParam(NPXT, NPYT, NPX_, NPY_, NPZ_,&
            &NSTLEG, NLEG, NUMPHASE_,&
            &ZLEVELS_, MAXASYM,&
            &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, IPHASEP_)

        if (MASTERPROC) deallocate (TEMPPT, EXTINCTPT, ALBEDOPT, IPHASEPT, LEGENPT)

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readPrp(PROPFILE, NPX_, NPY_, NPZ_,&
            &NSTLEG, MAXLEG, NLEG, MAXNZ, MAXPG, MAXPGL, DELTAM,&
            &PROPTYPE, DELX_, DELY_, ZLEVELS_, MAXASYM,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, NUMPHASE_, IPHASEP_)

        ! Check values of medium properties
        call shd_io__checkPrp(NPX_, NPY_, NPZ_, NSTLEG, NLEG,&
            &DELX_, DELY_, ZLEVELS_,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, NUMPHASE_, IPHASEP_)

    end if

    if (MASTERPROC) then
        ! Check input parameters
        call shd_io__checkInput(NMU, NPHI, DELTAM, MAXASYM,&
            &GRIDTYPE, SRCTYPE, UNITS, SOLARFLUX0, SOLARMU,&
            &WAVENO, WAVELEN, GNDTEMP, GNDALBEDO,&
            &SPLITACC, SHACC)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read SFCFILE and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in the surface properties
    if (SFCFILE(1:4) .eq. 'NONE') then
        SFCTYPE = 'FL'
    else
        allocate (SFCPARMS(MAXSFCPARS*MAXSFCPTS))
        MEMWORD = MEMWORD + MAXSFCPARS*MAXSFCPTS

        if (MASTERPROC) then
            call shd_io__readSfc(SFCFILE, MAXSFCPTS, MAXSFCPARS,&
                &SFCTYPE, NXSFC, NYSFC, DELXSFC, DELYSFC,&
                &NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
        end if
        ! Broadcast the surface array data
        call shd_mpi__bcastSfcParam(SFCTYPE, NXSFC, NYSFC,&
            &NSFCPAR, SFCPARMS, GNDTEMP, GNDALBEDO)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read CKDFILE and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (KDIST) then
        if (MASTERPROC) then
            call shd_io__readCkdSize(CKDFILE, WAVENO, NG, NZCKD_)
        end if
        ! Broadcast the k-distribution array sizes
        call shd_mpi__bcastCkdSize(NG, NZCKD_)
        allocate (ZCKD_(NZCKD_), GASABS_(NZCKD_), KABS(NZCKD_, NG))
        allocate (SOLFLUX(NG), DELG(NG))
        MEMWORD = MEMWORD + NZCKD_*(2.0+NG) + 2.0*NG

        if (MASTERPROC) then
            MAXNG = NG; MAXNZ = NZCKD_
            call shd_io__readCkd(MAXNG, MAXNZ, CKDFILE, WAVENO, SOLFLUX,&
                &NG, DELG, NZCKD_, ZCKD_, KABS)
        end if
        call shd_mpi__bcastCkdParam(SOLFLUX, NG, DELG, NZCKD_, ZCKD_, KABS)
    else
        NG = 1
        allocate (SOLFLUX(NG), DELG(NG))
        MEMWORD = MEMWORD + 2.0*NG

        SOLFLUX(1) = 1.0
        DELG(1) = 1.0
        NZCKD_ = 0
        BASEOUT = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read SENFILE and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in user-specified sensor properties if given
    !   this is a binary angle input implemented for 'R' and 'V3'
    if (SENFILE(1:4) .ne. 'NONE') then
        if (MASTERPROC) then
            call shd_io__readSenSize(SENFILE, NSENDATA, NSENPAR)
        end if
        call shd_mpi__bcastSenSize(NSENDATA, NSENPAR)
        allocate (SENPARMS(NSENDATA*NSENPAR))
        MEMWORD = MEMWORD + NSENDATA*NSENPAR

        if (MASTERPROC) then
            call shd_io__readSen(SENFILE, NSENDATA, NSENPAR,&
                &SENPARMS)
        end if
        call shd_mpi__bcastSenParam(NSENDATA, NSENPAR, SENPARMS)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Check outputs
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Figure out what type of output we will be generating
    NRAD = 0
    NVIS = 0
    RADIANCE_OUT = .false.
    VISUAL_OUT   = .false.
    SOURCE_OUT   = .false.
    FLUX_OUT     = .false.
    FLUXDIV_OUT  = .false.
    SH_OUT       = .false.
    do I = 1, NUMOUT
        if (OUTTYPES(I) .eq. 'R') then
            RADIANCE_OUT = .true.
            if (OUTPARMS(2, I) .le. 0.0) then
                NXOUT = 1
            else
                if (MPIX) then
                    NXOUT = max(1, nint(XDOMAIN_/OUTPARMS(2, I)))
                else
                    NXOUT = max(1, nint((XDOMAIN_ + 2*OUTPARMS(4, I))/OUTPARMS(2, I)))
                end if
            end if
            if (OUTPARMS(3, I) .le. 0.0) then
                NYOUT = 1
            else
                if (MPIY) then
                    NYOUT = max(1, nint(YDOMAIN_/OUTPARMS(3, I)))
                else
                    NYOUT = max(1, nint((YDOMAIN_ + 2*OUTPARMS(5, I))/OUTPARMS(3, I)))
                end if
            end if
            NRAD = NRAD + NSTOKES*NXOUT*NYOUT*nint(OUTPARMS(6, I))
        else if (OUTTYPES(I) .eq. 'V') then
            VISUAL_OUT = .true.
            if (nint(OUTPARMS(1, I)) .eq. 1) then
                NVIS = NVIS + NSTOKES*nint(OUTPARMS(10, I))*nint(OUTPARMS(11, I))
            else if (nint(OUTPARMS(1, I)) .eq. 2) then
                NSCAN = 1 + sqrt((OUTPARMS(4, I) - OUTPARMS(7, I))**2&
                    + (OUTPARMS(5, I) - OUTPARMS(8, I))**2&
                    + (OUTPARMS(6, I) - OUTPARMS(9, I))**2)/OUTPARMS(10, I)
                NPIX = 1 + abs(OUTPARMS(12, I) - OUTPARMS(11, I))/OUTPARMS(13, I)
                NVIS = NVIS + NSTOKES*NSCAN*NPIX
            ! user-specified sensor
            else if (nint(OUTPARMS(1, I)) .eq. 3) then
                NVIS = NVIS + NSTOKES*nint(OUTPARMS(2, I))
            end if
        else if (OUTTYPES(I) .eq. 'F') then
            FLUX_OUT = .true.
        else if (OUTTYPES(I) .eq. 'H') then
            FLUXDIV_OUT = .true.
        else if (OUTTYPES(I) .eq. 'S') then
            SH_OUT = .true.
        else if (OUTTYPES(I) .eq. 'J') then
            SOURCE_OUT = .true.
        end if
    end do
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Notify user that certain outputs are not supported with multiple processors
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (MPIMODE .and. MASTERPROC) then
        do I = 1, NUMOUT
            if (OUTTYPES(I) .eq. 'F') then
                if (nint(OUTPARMS(1, I)) .eq. 2 .or. nint(OUTPARMS(1, I)) .eq. 5) then
                    if (VERBOSE_ .ge. 4) write (6, '(A,I1,A,A)') 'Warning [shdom]: F output format ', nint(OUTPARMS(1, I)),&
                        &' not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
                end if
            else if (OUTTYPES(I) .eq. 'H') then
                if (nint(OUTPARMS(1, I)) .eq. 3) then
                    if (VERBOSE_ .ge. 4) write (6, '(A,I1,A,A)') 'Warning [shdom]: H output format ', nint(OUTPARMS(1, I)),&
                        &' not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
                end if
            else if (OUTTYPES(I) .eq. 'S') then
                if (nint(OUTPARMS(1, I)) .eq. 2) then
                    if (VERBOSE_ .ge. 4) write (6, '(A,I1,A,A)') 'Warning [shdom]: S output format ', nint(OUTPARMS(1, I)),&
                        &' not supported for multiple processors: ', trim(adjustl(OUTFILES(I)))
                end if
            else if (OUTTYPES(I) .eq. 'J') then
                if (VERBOSE_ .ge. 4) write (6, '(A,A)') 'Warning [shdom]: J output format not supported for multiple processors: ',&
                    &trim(adjustl(OUTFILES(I)))
            else if (OUTTYPES(I) .eq. 'M') then
                if (VERBOSE_ .ge. 4) write (6, '(A,A)') 'Warning [shdom]: M output format not supported for multiple processors: ',&
                    &trim(adjustl(OUTFILES(I)))
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Memory check
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (MPIMODE) then
        ! Read the memory parameters from a file, if it exists, and send
        !   the parameters to all processors. Otherwise use the input ones.
        call shd_mpi__readBcastMemParam(MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
            &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, RUNNAME)
    end if

    ! Guess maximum number of grid points, cells, SH vector size needed
    !   but don't let MAX_TOTAL_MB be exceeded
    MAXCPUMEM_MB_ = 0.9*TOTSYSMEM_MB_/real(NUMPROC)
    if (MAX_TOTAL_MB .gt. MAXCPUMEM_MB_) then
        MAX_TOTAL_MB = MAXCPUMEM_MB_
        if (VERBOSE_ .ge. 5) write (6, '(A,F0.2," MB")') 'Message [shdom]: MAX_TOTAL_MB reduced to: ', MAX_TOTAL_MB
    else
        if (VERBOSE_ .ge. 5) write (6, '(A,F0.2," MB")') 'Message [shdom]: MAX_TOTAL_MB is set to: ', MAX_TOTAL_MB
    end if

    if ((SPLITACC .le. 1.0d-12) .and. (INSAVEFILE(1:4) .eq. 'NONE')) then
        ADAPT_GRID_FACTOR = 1.0
    end if

    BIG_ARRAYS = 3
    if (.not. ACCELFLAG) BIG_ARRAYS = 2

    ! NMAXIG
    !╭──────────────────────────────────────────────────────────────╮!
    NMAXIG = 0.0

    NMAXIG = NMAXIG + 4.5 ! TEMP, PLANCK, EXTINCT, ALBEDO, IPHASE[0.5, integer(kind=2)]

    if (NUMPHASE_ .gt. 0) then
        MEMWORD = MEMWORD + NUMPHASE_*NSTLEG*(MAXLEG + 1) ! offset phase table
    else
        NMAXIG = NMAXIG + NSTLEG*(MAXLEG + 1)
    end if

    if (SFCFILE(1:4) .ne. 'NONE') then
        if (SFCTYPE(2:2) .eq. 'L') then
            NMAXIG = real(3)/real(NZ) * (MAXSFCPARS + 2.0 + NSTOKES*2.0)
        else
            NMAXIG = real(3)/real(NZ) * (MAXSFCPARS + 2.0 + NSTOKES*(2.0+ NMU*NPHI0MAX/2))
        end if
    end if

    if (FLUX_OUT)    NMAXIG = NMAXIG + 3.0
    if (FLUXDIV_OUT) NMAXIG = NMAXIG + 1.0
    if (SH_OUT)      NMAXIG = NMAXIG + NSHOUT
    if (SOURCE_OUT)  NMAXIG = NMAXIG + NSTOKES

    NMAXIG = NMAXIG + 16.5*CELL_TO_POINT_RATIO ! GRIDPTR[8.0], NEIGHPTR[6.0], TREEPTR[2.0], CELLFLAGS[0.5]

    NMAXIG = NMAXIG + 6.0   ! GRIDPOS[3.0], RSHPTR[1.0], SHPTR[1.0], OSHPTR[1.0]
    MEMWORD = MEMWORD + 4.0 ! offset RSHPTR[+2], SHPTR[+1], OSHPTR[+1]

    if (ACCELFLAG) then
        NMAXIG = NMAXIG + NSTOKES*(1.0+3.0*NUM_SH_TERM_FACTOR*NLM) ! RADIANCE, SOURCE, DELSOURCE
    else
        NMAXIG = NMAXIG + NSTOKES*(1.0+2.0*NUM_SH_TERM_FACTOR*NLM) ! RADIANCE, SOURCE
    end if

    NMAXIG = NMAXIG + 3.0 ! FLUXES[2.0], DIRFLUX[1.0]

    MEMWORD = MEMWORD + (NRAD+NVIS) ! offset SUMRADOUT[NRAD] and SUMVISOUT[NVIS]

    NMAXIG = NMAXIG + (NSTOKES*(max(NPHI0MAX, 4)+1.0) + 8.0) !WORK, WORK1, WORK2

    MEMWORD = MEMWORD + NRAD + NVIS ! WORK2 (over-count but memory safe)
    !╰──────────────────────────────────────────────────────────────╯!

    ! MEMWORD is now more accurate than before
    ! before:
    !   MEMWORD = NMU*(2 + 2*NPHI + 2*NLM + 2*33*32)&            ! SHDO allocation
    !       + 4.5*MAXPG + MAXPGL + NUMPHASE_*NSTLEG*(MAXLEG + 1) ! phase allocation
    ! now: update MEMWORD after every allocation step
    !
    ! NMAXIG is now more accurate than before
    ! before:
    !   NMAXIG = 28.0 + 16.5*CELL_TO_POINT_RATIO&
    !       + NPHI0MAX*NSTOKES + NUM_SH_TERM_FACTOR*NSTOKES*NLM*BIG_ARRAYS
    ! now: update NMAXIG from every allocation step

    WANTMEM = ADAPT_GRID_FACTOR*NBPTS*NMAXIG

    if (HIGHORDERRAD_) then
        REDUCE_RATIO = min(2.0, ((MAX_TOTAL_MB*1024**2)/4 - MEMWORD)/WANTMEM)
    else
        REDUCE_RATIO = min(1.5, ((MAX_TOTAL_MB*1024**2)/4 - MEMWORD)/WANTMEM)
    end if

    ADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR*REDUCE_RATIO
    if (ADAPT_GRID_FACTOR .lt. 0.999999999d0) then
        stop 'Error [shdom]: MAX_TOTAL_MB memory limit exceeded with just base grid.'
    end if

    if (REDUCE_RATIO .lt. 1.0) then
        if (VERBOSE_ .ge. 5) write (6, '(A,F8.2)') 'Message [shdom]: ADAPT_GRID_FACTOR reduced to:', ADAPT_GRID_FACTOR
    else
        if (VERBOSE_ .ge. 5) write (6, '(A,F8.2)') 'Message [shdom]: ADAPT_GRID_FACTOR increased to:', ADAPT_GRID_FACTOR
    end if

    WANTMEM = REDUCE_RATIO*WANTMEM
    if (WANTMEM .gt. huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Number of words of memory exceeds max integer size:', WANTMEM
        stop
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Retrieve max numbers for array allocation after memory check
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! max number of internal grid points
    MAXIG = nint(ADAPT_GRID_FACTOR*NBPTS)

    ! max number of grid cells, between 1 and 2 times grid points
    MAXIC = nint(CELL_TO_POINT_RATIO*MAXIG)

    ! max words for source function and radiance arrays
    MAXIV = nint(NUM_SH_TERM_FACTOR*NLM*MAXIG)

    ! max words for temporary discrete ordinate azimuthal array
    MAXIDO = MAXIG*max(NPHI0MAX, 4)

    if (4.0*(MAXIV + MAXIG)*NSTOKES .gt. huge(MAXIV)) then
        write (6, '(A,I12)') 'Error [shdom]: Size of big SH arrays (MAXIV) probably exceeds max integer number of bytes:', MAXIV
        stop
    end if

    if (4.0*8.0*MAXIC .gt. huge(MAXIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Size of GRIDPTR array (8*MAXIC) probably exceeds max integer number of bytes:', 8*MAXIC
        stop
    end if

    MAXNBC = MAXIG*3/NZ
    if (SFCTYPE(2:2) .eq. 'L') then ! Lambertian surface
        MAXBCRAD = 2*MAXNBC
    else
        MAXBCRAD = (2 + NMU*NPHI0MAX/2)*MAXNBC
    end if

    ! max number of legendre coefficients
    if (NUMPHASE_ .gt. 0) then
        MAXIGL = NUMPHASE_*NSTLEG*(MAXLEG + 1)
    else
        MAXIGL = MAXIG*NSTLEG*(MAXLEG + 1)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Allocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! temperature and planck radiance (maybe we don't need this for shortwave)
    allocate (TEMP(MAXIG), PLANCK(MAXIG))

    ! extinction, single scattering albedo, and phase function indices
    allocate (EXTINCT(MAXIG), ALBEDO(MAXIG), IPHASE(MAXIG))

    ! legendre coefficients
    allocate (LEGEN(MAXIGL))

    ! boundary points pointers and radiances
    !   BCPTR(:, 1): pointers to grid points that are on top
    !   BCPTR(:, 2): pointers to grid points that are on bottom
    !   BCRAD(:[NSTOKES], 1:NTOPPTS): isotropic top radiances
    !   BCRAD(:[NSTOKES], NTOPPTS+1:NTOPPTS+NBOTPTS): upwelling bottom radiances
    !   BCRAD(:[NSTOKES], NTOPPTS+NBOTPTS+1:NTOPPTS+NBOTPTS+NBOTPTS*NANG/2):
    !         downwelling bottom radiances at discrete ordinates
    if (SFCFILE(1:4) .ne. 'NONE') allocate (SFCGRIDPARMS(MAXSFCPARS*MAXNBC))
    allocate (BCPTR(2*MAXNBC), BCRAD(NSTOKES, MAXBCRAD))

    ! output
    if (FLUX_OUT)    allocate (SUMFLUXES(2, MAXIG), SUMDIRFLUX(MAXIG))
    if (FLUXDIV_OUT) allocate (SUMFLUXDIV(MAXIG))
    if (SH_OUT)      allocate (SUMSHTERMS(MAXIG*NSHOUT))
    if (SOURCE_OUT)  allocate (SOURCE1OUT(MAXIG*NSTOKES))

    ! pointer to 8 grid points belonging to a cell
    ! pointer to 6 neighboring cells of a cell
    ! pointer to parent and children cells in tree structure
    allocate (GRIDPTR(8, MAXIC), NEIGHPTR(6, MAXIC), TREEPTR(2, MAXIC))

    ! a collection of binary flags for the cells
    allocate (CELLFLAGS(MAXIC))

    ! x,y,z coordinate values of grid point locations
    allocate (GRIDPOS(3, MAXIG))

    ! SH related array
    allocate (RSHPTR(MAXIG + 2), SHPTR(MAXIG + 1), OSHPTR(MAXIG + 1))

    ! radiance and source
    allocate (RADIANCE(NSTOKES, MAXIV + MAXIG), SOURCE(NSTOKES, MAXIV))

    ! delta source
    if (ACCELFLAG) allocate (DELSOURCE(NSTOKES, MAXIV))

    ! fluxes
    allocate (FLUXES(2, MAXIG), DIRFLUX(MAXIG))

    ! output radiances to be summed over g
    allocate (SUMRADOUT(NRAD), SUMVISOUT(NVIS))

    ! WORK: Work array for path integration; radiance field, Nazimuth*Npts
    allocate (WORK(NSTOKES*MAXIDO))

    ! WORK1: Work array for path integration; grid point sweeping order
    allocate (WORK1(8*MAXIG))

    ! WORK2: Work array for path integration; single angle radiance
    allocate (WORK2(max(NVIS, NRAD, NSTOKES*MAXIG)))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read in previous run from binary file if desired
    call shd_io__restoreState(INSAVEFILE, NSTOKES, NX, NY, NZ, ML, MM, NLM,&
        &INRADFLAG, NEWGRIDFLAG, XGRID, YGRID, ZGRID,&
        &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
        &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)


    ! Check if we can use pre-calculated "state", if not, generate new grids
    !   NEWGRIDFLAG is from shd_io__restoreState call
    if (NEWGRIDFLAG) then
        ! Make the internal base grid lines
        call shd_init__makeGrid(BCFLAG, GRIDTYPE, NPX_, NPY_, NPZ_, NX, NY, NZ,&
            &XSTART_, YSTART_, DELX_, DELY_, ZLEVELS_,&
            &XGRID, YGRID, ZGRID)
    end if


    MESSAGE = 'Message [shdom]: Total time used for preparation: '
    call shd_util__ticktock(MESSAGE, COUNT_START)
    if (MASTERPROC) then
        if (VERBOSE_ .ge. 5) write (6, '(A)') trim(adjustl(MESSAGE))
    end if


    ! solve radiative transfer equation
    !   1. Use Picard iteration to solve radiance and source
    !   2. Compute actual radiance (and/or other radiative quantities specified by user)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (VERBOSE_ .ge. 5) write (6, '(A,"<",A,">")') 'Message [shdom]: Starting solution iterations for SHDOM run: ', trim(adjustl(RUNNAME))

    TOTITER = 0
    OLDNPTS = 0
    SOLARFLUX = 0.0d0

    ! Loop over the k-distribution g's from high to low absorption
    CKD_LOOP: do IG = NG, 1, -1

        if (VERBOSE_ .ge. 5) write (6, '(A,I2,A,I2,A)') 'Message [shdom]: Starting ', NG-IG+1, '/', NG, ' (ig/Ng) ...'

        ! solar flux for given g of CKD
        SOLARFLUXCKD = SOLFLUX(IG)*SOLARFLUX0*abs(SOLARMU)

        ! weighted average solar flux
        SOLARFLUX = SOLARFLUX + SOLFLUX(IG)*DELG(IG)*SOLARFLUX0*abs(SOLARMU)

        ! If doing k-distribution then get the gas absorption profile.
        if (KDIST) then
            do I = 1, NZCKD_
                GASABS_(I) = KABS(I, IG)
            end do
        end if

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (BASEOUT .or. ((IG .eq. NG) .and. NEWGRIDFLAG)) then
            call shd_init__makeCell(BCFLAG, IPFLAG,&
                &NX, NY, NZ, NX1, NY1, NPTS, NCELLS,&
                &XGRID, YGRID, ZGRID, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
            NBCELLS = NCELLS
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in SOURCE (SHPTR), RADIANCE (RSHPTR), FLUXES, DIRFLUX.
        call shd_solv__solveRte(NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
            &ML, MM, NLM, NMU, NPHI, NLEG, NSTLEG, NUMPHASE_,&
            &NPHI0, MU, PHI, WTDO,&
            &MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG,&
            &BCFLAG, IPFLAG, DELTAM, SRCTYPE, HIGHORDERRAD_,&
            &SOLARFLUXCKD, SOLARMU, SOLARAZ, SKYRAD,&
            &SFCTYPE, GNDTEMP, GNDALBEDO,&
            &NXSFC, NYSFC, DELXSFC, DELYSFC,&
            &NSFCPAR, SFCPARMS, SFCGRIDPARMS,&
            &UNITS, WAVENO, WAVELEN,&
            &ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,&
            &SPLITACC, SHACC, XGRID, YGRID, ZGRID,&
            &TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,&
            &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
            &NPTS, GRIDPOS,&
            &NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
            &RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,&
            &SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX)
        INRADFLAG = .false.
        TOTITER = TOTITER + ITER

        ! Compute the output and put in arrays. First the radiance.
        if (RADIANCE_OUT) then
            IRAD = 0
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'R') then
                    if (MPIMODE) then
                        if (SENFILE(1:4) .ne. 'NONE') then
                            call shd_mpi__computeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, XDOMAIN_, YDOMAIN_, IPFLAG,&
                                &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, WORK, WORK(1 + NPTS*NSTOKES),&
                                &OUTPARMS(1, I), IRAD, WORK2, SENPARMS)
                        else
                            call shd_mpi__computeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, XDOMAIN_, YDOMAIN_, IPFLAG,&
                                &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, WORK, WORK(1 + NPTS*NSTOKES),&
                                &OUTPARMS(1, I), IRAD, WORK2)
                        end if
                    else
                        if (SENFILE(1:4) .ne. 'NONE') then
                            call shd_rad__computeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, WORK(1), WORK(1 + NPTS*NSTOKES),&
                                &OUTPARMS(1, I), IRAD, WORK2, SENPARMS)
                        else
                            call shd_rad__computeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, WORK(1), WORK(1 + NPTS*NSTOKES),&
                                &OUTPARMS(1, I), IRAD, WORK2)
                        end if
                    end if
                end if
            end do
            call shd_util__sumData(NG - IG, DELG(IG), NRAD, WORK2, SUMRADOUT)
        end if

        if (VISUAL_OUT) then
            IVIS = 0
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'V') then
                    if (MPIMODE) then
                        if (SENFILE(1:4) .ne. 'NONE') then
                            call shd_mpi__computeCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, XDOMAIN_, YDOMAIN_, IPFLAG,&
                                &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2, SENPARMS)
                        else
                            call shd_mpi__computeCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, XDOMAIN_, YDOMAIN_, IPFLAG,&
                                &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2)
                        end if
                    else
                        if (SENFILE(1:4) .ne. 'NONE') then
                            call shd_rad__computeCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2, SENPARMS)
                        else
                            call shd_rad__computeCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                                &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_,&
                                &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                                &BCFLAG, IPFLAG, SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
                                &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                                &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                                &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
                                &XGRID, YGRID, ZGRID, GRIDPOS,&
                                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                                &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
                                &SHPTR, SOURCE, OUTPARMS(1, I), IVIS, WORK2)
                        end if
                    end if
                end if
            end do
            call shd_util__sumData(NG-IG, DELG(IG), NVIS, WORK2, SUMVISOUT)
        end if

        if (SOURCE_OUT) then
            ! Compute source function throughout grid for this one angle
            do I = 1, NUMOUT
                if (OUTTYPES(I) .eq. 'J') then
                    MUOUT = OUTPARMS(2, I)
                    PHIOUT = OUTPARMS(3, I)*acos(-1.0)/180.0
                    call shd_src__computeSrcOneDir(NSTOKES, ML, MM, NLM,&
                        &NSTLEG, NLEG, NUMPHASE_,&
                        &NPTS, DELTAM, MUOUT, PHIOUT,&
                        &SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE,&
                        &DIRFLUX, SHPTR, SOURCE, SOURCE1OUT)
                end if
            end do
        end if

        if (BASEOUT) then
            ! If outputting only the base grid then throw out the other points
            call shd_init__makeCell(BCFLAG, IPFLAG,&
                &NX, NY, NZ, NX1, NY1,&
                &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS)
        else
            ! For the output from the grid point arrays (Flux, Net flux, SH)
            !   we may need to interpolate the running k-distribution sum
            !   to the newly created grid points.
            call shd_util__interpData(OLDNPTS, NPTS, FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
                SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
                NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        end if

        if (FLUX_OUT) then

            ! unscale delta-m fluxes
            if ((SRCTYPE .ne. 'T') .and. DELTAM) then

                allocate(DIRFLUX_TMP(NPTS))

                if (MPIMODE) then
                    call shd_mpi__makeDirect(1, NPTS, BCFLAG, IPFLAG, .false., ML,&
                        &NSTLEG, NLEG, SOLARFLUXCKD, SOLARMU, SOLARAZ,&
                        &GRIDPOS, NX, XGRID, NY, YGRID, DIRFLUX_TMP)
                else
                    call shd_rad__makeDirect(NPTS, BCFLAG, IPFLAG, .false., ML,&
                        &NSTLEG, NLEG, SOLARFLUXCKD, SOLARMU, SOLARAZ,&
                        &GRIDPOS, DIRFLUX_TMP)
                end if

                FLUXES(1, 1:NPTS) = FLUXES(1, 1:NPTS) + DIRFLUX(1:NPTS) - DIRFLUX_TMP(1:NPTS)
                DIRFLUX(1:NPTS) = DIRFLUX_TMP(1:NPTS)

                deallocate(DIRFLUX_TMP)

            end if

            call shd_util__sumData(NG - IG, DELG(IG), 2*NPTS, FLUXES, SUMFLUXES)
            call shd_util__sumData(NG - IG, DELG(IG), NPTS, DIRFLUX, SUMDIRFLUX)
        end if

        if (FLUXDIV_OUT) then
            ! Compute the net flux divergence at every point
            call shd_rad__computeNetFluxDiv(NSTOKES, NPTS, RSHPTR, SRCTYPE, SOLARMU,&
                &EXTINCT, ALBEDO, PLANCK, DIRFLUX, RADIANCE, WORK)
            call shd_util__sumData(NG - IG, DELG(IG), NPTS, WORK, SUMFLUXDIV)
        end if

        if (SH_OUT) then
            call shd_shdo__computeSH(NSHOUT, NSTOKES, NPTS, SRCTYPE, SOLARMU, SOLARAZ,&
                &DIRFLUX, RSHPTR, ML, MM, RADIANCE, WORK)
            call shd_util__sumData(NG - IG, DELG(IG), NSHOUT*NPTS, WORK, SUMSHTERMS)
        end if

        OLDNPTS = NPTS

    end do CKD_LOOP
    !╰────────────────────────────────────────────────────────────────────────────╯!


    MESSAGE = 'Message [shdom]: Total time used for runs: '
    call shd_util__ticktock(MESSAGE, COUNT_START)
    if (MASTERPROC) then
        if (VERBOSE_ .ge. 5) write (6, '(A)') trim(adjustl(MESSAGE))
    end if


    ! Output radiance and/or flux (and/or other radiative quantities)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    NSH = SHPTR(NPTS + 1)
    IRAD = 0
    IVIS = 0
    if (MPIMODE) then

        if (MASTERPROC) then
            if (FLUX_OUT) allocate (ALLFLUXES(3, NZ, NYT, NXT))
            if (FLUXDIV_OUT) allocate (ALLFLUXDIV(NZ, NYT, NXT))
            if (SH_OUT) allocate (ALLSHTERMS(NSHOUT, NZ, NYT, NXT))
            allocate (XALLGRID(NXT), YALLGRID(NYT))
        end if

        call shd_mpi__gatherOutput(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
            &IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,&
            &SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
            &XDOMAIN_, YDOMAIN_, XALLGRID, YALLGRID,&
            &ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,&
            &NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)

        if (MASTERPROC) then

            do I = 1, NUMOUT
                if (OUTTYPES(I) .ne. 'V' .and. OUTFILES(I)(1:4) .ne. 'NONE') then
                    if (SENFILE(1:4) .ne. 'NONE') then
                        call shd_io__outputResultMpi(NSTOKES, NXT, NYT, NZ,&
                            &NPTSTOT, NCELLSTOT, NSHTOT,&
                            &ML, MM, NLM, NMU, NPHI, NANG, NG,&
                            &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                            &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                            &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                            &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                            &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                            &XDOMAIN_, YDOMAIN_, XALLGRID, YALLGRID, ZGRID,&
                            &ALLFLUXES, ALLFLUXDIV,&
                            &NSHOUT, ALLSHTERMS, IRAD, SUMRADOUT,&
                            &OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I), SENPARMS)
                    else
                        call shd_io__outputResultMpi(NSTOKES, NXT, NYT, NZ,&
                            &NPTSTOT, NCELLSTOT, NSHTOT,&
                            &ML, MM, NLM, NMU, NPHI, NANG, NG,&
                            &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                            &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                            &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                            &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                            &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                            &XDOMAIN_, YDOMAIN_, XALLGRID, YALLGRID, ZGRID,&
                            &ALLFLUXES, ALLFLUXDIV,&
                            &NSHOUT, ALLSHTERMS, IRAD, SUMRADOUT,&
                            &OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I))
                    end if
                end if
            end do

            if (FLUX_OUT) deallocate (ALLFLUXES)
            if (FLUXDIV_OUT) deallocate (ALLFLUXDIV)
            if (SH_OUT) deallocate (ALLSHTERMS)
            deallocate (XALLGRID, YALLGRID)
        end if

    else

        do I = 1, NUMOUT
            if ((OUTTYPES(I) .ne. 'V') .and. (OUTFILES(I)(1:4) .ne. 'NONE')) then
                if (SENFILE(1:4) .ne. 'NONE') then
                    call shd_io__outputResult(NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NCELLS,&
                        NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_, NMU, NPHI, NANG, NG,&
                        PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                        BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                        SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                        SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                        SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                        XGRID, YGRID, ZGRID, GRIDPOS,&
                        TREEPTR, GRIDPTR, CELLFLAGS,&
                        EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,&
                        SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,&
                        IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT,&
                        OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I), SENPARMS)
                else
                    call shd_io__outputResult(NSTOKES, NX, NY, NX1, NY1, NZ, NBPTS, NPTS, NCELLS,&
                        NSH, ML, MM, NLM, NSTLEG, NLEG, NUMPHASE_, NMU, NPHI, NANG, NG,&
                        PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                        BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                        SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                        SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                        SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                        XGRID, YGRID, ZGRID, GRIDPOS,&
                        TREEPTR, GRIDPTR, CELLFLAGS,&
                        EXTINCT, ALBEDO, LEGEN, IPHASE, TEMP,&
                        SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,&
                        IRAD, SUMRADOUT, NSHOUT, SUMSHTERMS, SOURCE1OUT,&
                        OUTTYPES(I), OUTPARMS(1, I), OUTFILES(I))
                end if
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Output camera radiance if there are any
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (MASTERPROC) then
        do I = 1, NUMOUT
            if ((OUTTYPES(I) .eq. 'V')) then
                if (SENFILE(1:4) .ne. 'NONE') then
                    call shd_io__outputCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                        &NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
                        &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                        &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                        &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                        &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                        &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                        &IVIS, SUMVISOUT, OUTPARMS(1, I), OUTFILES(I), SENPARMS)
                else
                    call shd_io__outputCamRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
                        &NSH, ML, MM, NLM, NMU, NPHI, NANG, NG,&
                        &PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE,&
                        &BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE,&
                        &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
                        &SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS,&
                        &SPLITACC, SHACC, SOLACC, MAXITER, TOTITER,&
                        &IVIS, SUMVISOUT, OUTPARMS(1, I), OUTFILES(I))
                end if
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! save state
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (MAXITER .gt. 0) then
        call shd_io__saveState(OUTSAVEFILE, NSTOKES, NX, NY, NZ, ML, MM, NLM,&
            &WORK, XGRID, YGRID, ZGRID,&
            &NPTS, NCELLS, GRIDPOS, GRIDPTR, NEIGHPTR, TREEPTR,&
            &CELLFLAGS, FLUXES, SHPTR, SOURCE, RSHPTR, RADIANCE)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    MESSAGE = 'Message [shdom]: Total time used for output: '
    call shd_util__ticktock(MESSAGE, COUNT_START)
    if (MASTERPROC) then
        if (VERBOSE_ .ge. 5) write (6, '(A)') trim(adjustl(MESSAGE))
    end if


    ! Calculate the actual memory parameters and output them to the log file
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! MAXMB is now more accurate than before
    ! before:
    !  MAXMB_OUT = 4*(NMU*(2 + 2*NPHI + 2*NLM + 2*33*32)&
    !      + 4.5*MAXPG + MAXPGL + NUMPHASE_*NSTLEG*(MAXLEG + 1)&
    !      + 16.5*NCELLS + NPTS*(28 + NPHI0MAX*NSTOKES)&
    !      + NSTOKES*NSH*BIG_ARRAYS)/1024**2
    ! now:
    !   use NMAXIG
    MAXMB_OUT = 4*(MEMWORD + 16.5*NCELLS&
        + NPTS*(NMAXIG-16.5*CELL_TO_POINT_RATIO))/1024**2

    ADAPTGRIDFAC_OUT = float(NPTS)/NBPTS
    SHTERM_FAC_OUT = float(NSH)/(NLM*NPTS)
    CELL_POINT_OUT = float(NCELLS)/NPTS
    if (VERBOSE_ .ge. 5) write (6, '(A)')'Message [shdom]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    if (VERBOSE_ .ge. 5) write (6, '(26X,F9.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Deallocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    deallocate (ZLEVELS_, TEMPP_, EXTINCTP_, ALBEDOP_, EXTDIRP_, IPHASEP_, LEGENP_)
    if (FLUX_OUT) deallocate (SUMFLUXES, SUMDIRFLUX)
    if (FLUXDIV_OUT) deallocate (SUMFLUXDIV)
    if (SH_OUT) deallocate (SUMSHTERMS)
    if (SOURCE_OUT) deallocate (SOURCE1OUT)
    if (KDIST) deallocate (ZCKD_, GASABS_, KABS)
    if (allocated(SFCPARMS)) deallocate (SFCPARMS, SFCGRIDPARMS)
    if (allocated(SENPARMS)) deallocate (SENPARMS)
    deallocate (DELG, SOLFLUX)
    deallocate (BCPTR, BCRAD)

    deallocate (MU, WTDO, PHI, NPHI0)
    deallocate (XGRID, YGRID, ZGRID)
    deallocate (TEMP, PLANCK)
    deallocate (EXTINCT, ALBEDO, IPHASE)
    deallocate (LEGEN)
    deallocate (GRIDPTR, NEIGHPTR, TREEPTR)
    deallocate (CELLFLAGS)
    deallocate (GRIDPOS)
    deallocate (RSHPTR, SHPTR, OSHPTR)
    deallocate (RADIANCE, SOURCE)
    if (ACCELFLAG) deallocate (DELSOURCE)
    deallocate (FLUXES, DIRFLUX)
    deallocate (SUMRADOUT, SUMVISOUT)
    deallocate (WORK, WORK1, WORK2)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    call shd_mpi__end(NPTS, GRIDPOS, NPX_, NPY_, XSTART_, YSTART_, XDOMAIN_, YDOMAIN_,&
        &NPXT, NPYT, PROPFILE, RUNNAME)

end program shdom
