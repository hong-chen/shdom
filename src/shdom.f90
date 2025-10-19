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

program shdom

    use globals_const, only: REPS_
    use globals_cfg,   only: verbose_, fnameMaxLen_, RAcc_, &
                             TOTSYSMEM_MB_, MAXCPUMEM_MB_, &
                             isNmlInp_, isHighOrder_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, isMpi_, &
                             is1DX_, is1DYx_, is1D_, is3Ddir1Ddif_, &
                             is1D1_, isPol_, isPol4_, &
                             isSolar_, isThermal_


    use globals_prp,   only: NPX_, NPY_, NPZ_, NUMPHASE_, &
                             delX_, delY_, XSTART_, YSTART_, XDOMAIN_, YDOMAIN_, &
                             ZLEVELS_, TEMPP_, &
                             EXTINCTP_, EXTDIRP_, ALBEDOP_, &
                             LEGENP_, IPHASEP_, &
                             NzCkd_, zCkd_, GASABS_

    use shd_init, only: shd_init__makeGrid, &
                        shd_init__makeCell

    use shd_src,  only: shd_src__computeSrcOneDir

    use shd_shdo, only: shd_shdo__computeSH

    use shd_solv, only: shd_solv__solveRte

    use shd_rad,  only: shd_rad__computeRad, &
                        shd_rad__computeCamRad, &
                        shd_rad__makeDirect, &
                        shd_rad__computeNetFluxDiv

    use shd_inp,  only: shd_io__nmlInp,      shd_io__userInp, shd_io__checkInp, &
                        shd_io__readPrpSize, shd_io__readPrp, shd_io__checkPrp, &
                        shd_io__readSfcSize, shd_io__readSfc, &
                        shd_io__readCkdSize, shd_io__readCkd, &
                        shd_io__readSenSize, shd_io__readSen, &
                        shd_io__restoreState

    use shd_out,  only: shd_io__outputResult, &
                        shd_io__outputResultMpi, &
                        shd_io__outputCamRad, &
                        shd_io__saveState

    use shd_mpi,  only: shd_mpi__start, shd_mpi__map, shd_mpi__end, &
                        shd_mpi__bcastInp, &
                        shd_mpi__bcastMem, &
                        shd_mpi__bcastPrpSize, shd_mpi__bcastPrp, &
                        shd_mpi__bcastSfcSize, shd_mpi__bcastSfc, &
                        shd_mpi__bcastCkdSize, shd_mpi__bcastCkd, &
                        shd_mpi__bcastSenSize, shd_mpi__bcastSen, &
                        shd_mpi__computeRad, &
                        shd_mpi__computeCamRad, &
                        shd_mpi__makeDirect, &
                        shd_mpi__gatherOut

    use shd_util, only: shd_util__sumData, &
                        shd_util__interpData, &
                        shd_util__tickTock

    implicit none

    integer, parameter :: R__ = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__ = selected_int_kind(9)  ! default precision for integer

    ! Output
    !╭────────────────────────────────────────────────────────────────────────────╮!
    integer(kind=I__), parameter :: maxNpar_out=100, maxNset_out=20
    real(kind=R__)   :: parms_Out(maxNpar_out, maxNset_out)
    character(len=1) :: type_Out(maxNset_out)
    character(len=fnameMaxLen_) :: file_Out(maxNset_out)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=2)  , allocatable :: iPhase(:)
    integer(kind=2)  , allocatable :: iPhasePT(:)
    integer(kind=I__), allocatable :: RSHptr(:), SHptr(:), OSHptr(:)
    integer(kind=I__), allocatable :: Nphi0(:)
    integer(kind=I__), allocatable :: work1(:)

    real(kind=R__), allocatable :: temp(:), planck(:), extinct(:), ssAlb(:)
    real(kind=R__), allocatable :: legen(:)
    real(kind=R__), allocatable :: radiance(:, :), source(:, :), delSource(:, :)
    real(kind=R__), allocatable :: work(:), work2(:)
    real(kind=R__), allocatable :: fluxes(:, :), dirFlux(:), dirFlux_tmp(:)
    real(kind=R__), allocatable :: mu(:), wtDO(:), phi(:)
    real(kind=R__), allocatable :: sumFluxes(:, :), sumDirFlux(:), sumFluxDiv(:)
    real(kind=R__), allocatable :: sumSHterms(:), source1out(:)
    real(kind=R__), allocatable :: sumRadOut(:), sumCamOut(:)
    real(kind=R__), allocatable :: tempPT(:), extinctPT(:), ssAlbPT(:)
    real(kind=R__), allocatable :: legenPT(:)
    real(kind=R__), allocatable :: allFluxes(:, :, :, :)
    real(kind=R__), allocatable :: allFluxDiv(:, :, :), allSHterms(:, :, :, :)
    real(kind=R__), allocatable :: xAllGrid(:), yAllGrid(:)

    ! SH and DO
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! PRP
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=fnameMaxLen_) :: file_prp
    character(len=1) :: type_prp
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SFC
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=fnameMaxLen_) :: file_sfc
    character(len=2)  :: type_sfc
    integer(kind=I__) :: Nx_sfc, Ny_sfc, Npar_sfc
    integer(kind=I__) :: maxPts_sfc, maxPar_sfc
    real(kind=R__)    :: dx_sfc, dy_sfc
    real(kind=R__), allocatable :: parms_Sfc(:), &
                                   parmsGrid_Sfc(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! CKD
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=fnameMaxLen_) :: file_ckd
    integer(kind=I__) :: Ng_ckd, ig_ckd
    real(kind=R__)    :: solar0_ckd
    real(kind=R__), allocatable :: solar_Ckd(:), &
                                   kabs_Ckd(:, :), &
                                   delg_Ckd(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SEN
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=fnameMaxLen_) :: file_sen
    integer(kind=I__) :: Npts_sen, Npar_sen
    real(kind=R__), allocatable :: parms_Sen(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Bndry
    !╭────────────────────────────────────────────────────────────────────────────╮!
    integer(kind=I__) :: NptsTop_bndry, NptsBot_bndry
    integer(kind=I__) :: maxNpts_bndry, maxNrad_bndry
    integer(kind=I__), allocatable :: ptr_Bndry(:)
    real(kind=R__), allocatable    :: rad_Bndry(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! grids and cells
    !╭────────────────────────────────────────────────────────────────────────────╮!
    real(kind=R__), allocatable    :: x_Grid(:), y_Grid(:), z_Grid(:)
    real(kind=R__), allocatable    :: xyz_Grid(:, :)

    integer(kind=2)  , allocatable :: flags_Cell(:)
    integer(kind=I__), allocatable :: ptr8Grid_Cell(:, :), &
                                      ptr6Cell_Cell(:, :), &
                                      ptr2Tree_Cell(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    integer(kind=I__) :: maxPG, maxIG, maxIC, maxIV, maxIDO, maxNz, maxNg
    integer(kind=I__) :: maxLeg, maxPGL, maxIGL, NbigArr

    real(kind=R__) :: memUsed, memNeed, scaleRatio, NmaxIG
    real(kind=R__) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real(kind=R__) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real(kind=R__) :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer(kind=I__) :: Nproc, Nstokes, Nx, Ny, Nz, Nx1, Ny1
    integer(kind=I__) :: NxT, NyT, NpxT, NpyT, maxPGT, maxPGLT
    integer(kind=I__) :: maxNl, maxNm, Nlm, Nstleg, Nleg, Nang, Nmu, Nphi, Nphi0max
    integer(kind=I__) :: NptsB, Npts, Ncells, Nsh, NcellsB, oldNpts
    integer(kind=I__) :: maxIter, iter, totIter, iz, iout, Nout, bcFlag, ipFlag
    integer(kind=I__) :: Nrad, irad, NxOut, NyOut, Nlos, ilos, Nscan, Npix, NshOut
    integer(kind=I__) :: NcellsTot, NptsTot, NshTot
    integer(kind=I__) :: countStart, countRate, countMax

    real(kind=R__) :: solar0, solar, mu_sol, az_sol
    real(kind=R__) :: gndTemp, gndAlbedo, skyRad, waveNum(2), waveLen
    real(kind=R__) :: solAcc, solCrit, splitAcc, shAcc, maxAsym
    real(kind=R__) :: mu_out, phi_out
    real(kind=R__) :: delX0, delY0

    character(len=1) :: srcType, units
    character(len=1) :: gridType
    character(len=fnameMaxLen_) :: inSaveFile, outSaveFile
    character(len=fnameMaxLen_) :: outNcFile, runName
    character(len=fnameMaxLen_) :: message

    logical :: isCkd, INRADFLAG, NEWGRIDFLAG, isDeltaM, isAccel
    logical :: isBaseOut, isMasterProc
    logical :: isRadOut, isCamOut, isFluxOut, isHROut, isSHOut, isSrcOut


    ! start timing
    call system_clock(countStart, countRate, countMax)


    ! Start up MPI for multiple processors if available, and
    !   return flag for master process (procnum=0)
    call shd_mpi__start(isMasterProc, Nproc)


    ! Read input from namelist or user input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! In MPI, the master processor performs the main computation,
    !   and the results are then broadcast to the child processors.
    if (isMasterProc) then
        if (isNmlInp_) then
            call shd_io__nmlInp(runName, file_prp, file_sfc, file_ckd,&
                inSaveFile, outSaveFile, file_sen,&
                Nstokes, Nx, Ny, Nz, Nmu, Nphi, bcFlag,&
                ipFlag, isCkd, isDeltaM, gridType,&
                srcType, solar0, mu_sol, az_sol, skyRad,&
                gndTemp, gndAlbedo, units, waveNum, waveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                outNcFile, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        else
            call shd_io__userInp(runName, file_prp, file_sfc, file_ckd,&
                inSaveFile, outSaveFile, file_sen,&
                Nstokes, Nx, Ny, Nz, Nmu, Nphi, bcFlag,&
                ipFlag, isCkd, isDeltaM, gridType,&
                srcType, solar0, mu_sol, az_sol, skyRad,&
                gndTemp, gndAlbedo, units, waveNum, waveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                outNcFile, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        end if
    end if

    call shd_mpi__bcastInp(runName, file_prp, file_sfc, file_ckd,&
        &inSaveFile, outSaveFile, file_sen,&
        &Nstokes, Nx, Ny, Nz, Nmu, Nphi, bcFlag,&
        &ipFlag, isCkd, isDeltaM, gridType,&
        &srcType, solar0, mu_sol, az_sol, skyRad,&
        &gndTemp, gndAlbedo, units, waveNum, waveLen,&
        &isAccel, solAcc, maxIter, splitAcc, shAcc,&
        &maxNset_out, maxNpar_out, Nout, type_Out, parms_Out,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Memory used for initialization (rough estimate, 1 WORD = 4 bytes)
    memUsed = 2000.0 + maxNset_out * (maxNpar_out + 0.25 + (fnameMaxLen_/4.0))


    ! Post-process parameters from input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! set calculation accuracy (1.0e-4 to 1.0e-8)
    RAcc_ = max(REPS_, min(solAcc*0.1, 1.0e-4))

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    Nstokes = max(1, min(4, Nstokes))
    if (Nstokes .eq. 2) then
        Nstokes = 3
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: Nstokes cannot be 2, reset Nstokes to 3.'
    end if

    ! Nstleg is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (Nstokes .eq. 1) then
        Nstleg = 1
        isPol_ = .false.
    else
        Nstleg = 6
        isPol_ = .true.
    end if

    if (Nstokes .eq. 4) then
        isPol4_ = .true.
    else
        isPol4_ = .false.
    end if

    ! Make maxNl and maxNm from Nmu and Nphi (for discrete ordinate)
    !   Nmu : number of discrete ordinate zenith angles in both hemispheres
    !   Nphi: number of discrete ordinate azimuth angles in 2*pi
    Nmu  = max(2, 2*int((Nmu + 1)/2)) ! force Nmu an even number
    Nphi = max(1, Nphi)

    ! maxNl : The maximum L of the spherical harmonic truncation
    ! maxNm : The maximum M of the spherical harmonic truncation
    maxNl = Nmu - 1
    maxNm = max(0, int(Nphi/2) - 1)

    ! Compute Nlm, Nphi0max.
    !   Nlm: The number of spherical harmonic terms
    !   Nphi0max: The maximum number of azimuth angles actually used
    Nlm = (2*maxNm + 1)*(maxNl + 1) - maxNm*(maxNm + 1)
    Nphi0max = Nphi

    ! Nleg: The number of Legendre series phase function terms used in the solution procedure
    Nleg = maxNl
    if (isDeltaM) Nleg = maxNl + 1

    if (srcType .eq. 'B') then
        isSolar_ = .true.
        isThermal_ = .true.
    else
        if (srcType .eq. 'S') then
            isSolar_ = .true.
            isThermal_ = .false.
        else if (srcType .eq. 'T') then
            isSolar_ = .false.
            isThermal_ = .true.
        else
            stop 'Error [shdom]: Unknown srcType.'
        end if
    end if

    ! for thermal, reset solar azimuth angle
    if (srcType .eq. 'T') az_sol= 0.0d0

    ! angle grids for DO
    allocate (mu(Nmu), wtDO(Nmu*Nphi), phi(Nmu*Nphi), Nphi0(Nmu))

    memUsed = memUsed + 2.0*Nmu*(1.0 + Nphi)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read file_prp (property file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        ! Find out the size of the input property file grid
        call shd_io__readPrpSize(file_prp, Nstleg, Nleg, NPX_, NPY_, NPZ_, delX_, delY_,&
            &NUMPHASE_, maxLeg, maxPGL)
    end if
    call shd_mpi__bcastPrpSize(NPX_, NPY_, NPZ_, delX_, delY_,&
        &NUMPHASE_, maxLeg, maxPGL)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read file_sfc (surface file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (file_sfc(1:4) .eq. 'NONE') then
        type_sfc = 'FL'
        Nx_sfc = 1
        Ny_sfc = 1
        dx_sfc = delX_
        dy_sfc = delY_
    else
        if (isMasterProc) then
            ! find out the size of the surface file grid
            call shd_io__readSfcSize(file_sfc, Nx_sfc, Ny_sfc, dx_sfc, dy_sfc, maxPts_sfc, maxPar_sfc)
        end if
        ! Broadcast the surface array sizes
        call shd_mpi__bcastSfcSize(Nx_sfc, Ny_sfc, dx_sfc, dy_sfc, maxPts_sfc, maxPar_sfc)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Process grids
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! X, Y, Z grids
    Nx = max(1, Nx); Ny = max(1, Ny); Nz = max(2, Nz)

    ! If single plane or column then force IP mode
    !   if   3D ipFlag=0, bit representation '0000'
    !   if Nx=1 ipFlag=1, bit representation '0001'
    !   if Ny=1 ipFlag=2, bit representation '0010'
    !   if  IPA ipFlag=3, bit representation '0011'
    if (Nx .eq. 1) then
        ipFlag = ior(ipFlag, 1) ! '0000'(3D,   0) + '0001'(Nx=1, 1) = '0001' (Nx=1, 1)
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: Nx=1, reset ipFlag to 1.'
    end if
    if (Ny .eq. 1) then
        ipFlag = ior(ipFlag, 2) ! '0001'(Nx=1, 1) + '0010'(Ny=1, 2) = '0011' (IPA,  3)
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: Ny=1, reset ipFlag to 2 (or 3).'
    end if

    if ((Nx .eq. 1) .and. (Ny .eq. 1)) then
        is1D1_ = .true.
    else
        is1D1_ = .false.
    end if

    ! Set the full domain grid sizes (variables end in T)
    NxT  = Nx  ; NyT  = Ny
    NpxT = NPX_; NpyT = NPY_

    if (NpxT .eq. 1) then
        delX_ = Nx_sfc*dx_sfc! overwrite delX_ if only 1 grid in X direction
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NpxT=1, reset delX_ to match surface domain.'
    end if
    if (NpyT .eq. 1) then
        delY_ = Ny_sfc*dy_sfc! overwrite delY_ if only 1 grid in Y direction
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: NpyT=1, reset delY_ to match surface domain.'
    end if

    XDOMAIN_ = NpxT*delX_
    YDOMAIN_ = NpyT*delY_

    ! Choose the option to compute higher order radiance SH terms
    if ((NpxT .eq. 1) .and. (NpyT .eq. 1) .and.&
        &(Nx_sfc.le. 32) .and. (Nx_sfc.le. 32)) then
        isHighOrder_ = .true.
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: 1D calculation, set isHighOrder_ to .TRUE.'
    else
        isHighOrder_ = .false.
    end if

    NshOut = 4
    if (isHighOrder_) NshOut = 5

    ! x, y, z grids
    allocate (x_Grid(Nx + 1), y_Grid(Ny + 1), z_Grid(Nz))
    memUsed = memUsed + (Nx+1) + (Ny+1) + (Nz)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! If MPI, divide domain into mini-domains
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and Nx,Ny to the subdomain size, and output
    !   the X and Y starting position of this subdomain. Set bits
    !   2 and 3 in bcFlag to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    call shd_mpi__map(bcFlag, NPX_, NPY_, Nx, Ny, Nz, delX_, delY_, file_prp,&
        &XSTART_, YSTART_, runName)

    ! bcFlag bit-representation '...0001' isOpenX_
    if (btest(bcFlag, 0)) then
        isOpenX_ = .true.
    else
        isOpenX_ = .false.
    end if

    ! bcFlag bit-representation '...001x' isOpenYx_
    if (btest(bcFlag, 1)) then
        isOpenYx_ = .true.
    else
        isOpenYx_ = .false.
    end if

    ! bcFlag bit-representation '...01xx' isMpiX_
    if (btest(bcFlag, 2)) then
        isMpiX_ = .true.
    else
        isMpiX_ = .false.
    end if

    ! bcFlag bit-representation '...1xxx' isMpiY_
    if (btest(bcFlag, 3)) then
        isMpiY_ = .true.
    else
        isMpiY_ = .false.
    end if

    if (isMpiX_ .or. isMpiY_) then
        isMpi_ = .true.
    else
        isMpi_ = .false.
    end if

    ! btest(ipFlag, 0) checks whether ipFlag=1 or ipFlag=0 ('...0001' [IPA in X] or '...0000' [3D in X])
    if (btest(ipFlag, 0)) then
        is1DX_ = .true.
    else
        is1DX_ = .false.
    end if

    ! btest(ipFlag, 1) checks whether ipFlag=2 [or 3] or ipFlag=0 [or 1]
    !   ('...0010' [or '...0011', IPA in Y] or '...0000' [or '...0001', 3D in Y])
    if (btest(ipFlag, 1)) then
        is1DYx_ = .true.
    else
        is1DYx_ = .false.
    end if

    ! btest(ipFlag, 1) checks whether ipFlag=7 [or 8]
    !   ('...0110' [or '...0111'], 3D direct beam but IP diffuse radiative transfer)
    if (btest(ipFlag, 2)) then
        is3Ddir1Ddif_ = .true.
    else
        is3Ddir1Ddif_ = .false.
    end if

    if (is1DX_ .and. is1DYx_) then
        is1D_ = .true.
    else
        is1D_ = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Base grids setup after MPI map
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Choose type of initialization/output for k-distributions:
    !   True means initiatize with and output only the base grid.
    isBaseOut = .false.

    ! Set up base grid point actual size (Nx1xNY1xNZ)
    !   btest(bcFlag, 0) checks whether bcFlag=1 or bcFlag=0 ('...0001' [open in X] or '...0000' [periodic in X])
    Nx1 = Nx + 1
    if (is1DX_ .or. isOpenX_ .or. isMpiX_) Nx1 = Nx

    ! btest(bcFlag, 1) checks whether bcFlag=2 [or 3] or bcFlag=0 [or 1]
    !   ('...0010' [or '...0011', open in Y] or '...0000' [or '...0001', periodic in Y])
    Ny1 = Ny + 1
    if (is1DYx_ .or. isOpenYx_ .or. isMpiY_) Ny1 = Ny

    ! Number of base points
    NptsB = Nx1*Ny1*Nz

    ! Calculate the number of base grid cells depending on the bcFlag
    NcellsB = (Nz - 1)*(Nx + ibits(bcFlag, 0, 1) - ibits(bcFlag, 2, 1))&
        *(Ny + ibits(bcFlag, 1, 1) - ibits(bcFlag, 3, 1))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read file_prp and distribute data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Determine the size of the property grid arrays so can allocate
    maxNz = NPZ_
    maxPG = NPX_*NPY_*NPZ_
    if (NUMPHASE_ .eq. 0) maxPGL = maxPG*(maxLeg + 1)*Nstleg

    allocate (ZLEVELS_(NPZ_))
    allocate (TEMPP_(maxPG), EXTINCTP_(maxPG), EXTDIRP_(maxPG), ALBEDOP_(maxPG))
    allocate (IPHASEP_(maxPG), LEGENP_(maxPGL))
    memUsed = memUsed + NPZ_ + 6.0*maxPG + maxPGL

    if (isMpi_) then

        ! Read the properties of the medium to the master processor and
        !   distribute to all the processors
        if (isMasterProc) then

            maxPGT = NpxT*NpyT*NPZ_
            maxPGLT = maxPGL
            if (NUMPHASE_ .eq. 0) maxPGLT = maxPGT*(maxLeg + 1)*Nstleg

            allocate (tempPT(maxPGT), extinctPT(maxPGT), ssAlbPT(maxPGT))
            allocate (iPhasePT(maxPGT), legenPT(maxPGLT))

            ! Only use delX and delY from initial read, delX0 and delY0 are discarded
            call shd_io__readPrp(file_prp, NpxT, NpyT, NPZ_,&
                &Nstleg, maxLeg, Nleg, maxNz, maxPGT, maxPGLT, isDeltaM,&
                &type_prp, delX0, delY0, ZLEVELS_, maxAsym,&
                &tempPT, extinctPT, ssAlbPT, legenPT, NUMPHASE_, iPhasePT)

            if (NUMPHASE_ .eq. 0) then
                stop 'Error [shdom]: Standard property file (with phase function Legendre series at every point)&
                    & is not allowed with multiple processors.'
            end if

            ! Check values of medium properties
            call shd_io__checkPrp(NpxT, NpyT, NPZ_, Nstleg, Nleg,&
                &delX_, delY_, ZLEVELS_,&
                &tempPT, extinctPT, ssAlbPT, legenPT, NUMPHASE_, iPhasePT)

        end if

        ! distribute file grid data to working grids (different CPU processes different grids)
        call shd_mpi__bcastPrp(NpxT, NpyT, NPX_, NPY_, NPZ_,&
            &Nstleg, Nleg, NUMPHASE_,&
            &ZLEVELS_, maxAsym,&
            &tempPT, extinctPT, ssAlbPT, legenPT, iPhasePT,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, IPHASEP_)

        if (isMasterProc) deallocate (tempPT, extinctPT, ssAlbPT, iPhasePT, legenPT)

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readPrp(file_prp, NPX_, NPY_, NPZ_,&
            &Nstleg, maxLeg, Nleg, maxNz, maxPG, maxPGL, isDeltaM,&
            &type_prp, delX_, delY_, ZLEVELS_, maxAsym,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, NUMPHASE_, IPHASEP_)

        ! Check values of medium properties
        call shd_io__checkPrp(NPX_, NPY_, NPZ_, Nstleg, Nleg,&
            &delX_, delY_, ZLEVELS_,&
            &TEMPP_, EXTINCTP_, ALBEDOP_, LEGENP_, NUMPHASE_, IPHASEP_)

    end if

    if (isMasterProc) then
        ! Check input parameters
        call shd_io__checkInp(Nmu, Nphi, isDeltaM, maxAsym,&
            &gridType, srcType, units, solar0, mu_sol,&
            &waveNum, waveLen, gndTemp, gndAlbedo,&
            &splitAcc, shAcc)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read file_sfc and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in the surface properties
    if (file_sfc(1:4) .eq. 'NONE') then
        type_sfc = 'FL'
    else
        allocate (parms_Sfc(maxPar_sfc*maxPts_sfc))
        memUsed = memUsed + maxPar_sfc*maxPts_sfc

        if (isMasterProc) then
            call shd_io__readSfc(file_sfc, maxPts_sfc, maxPar_sfc,&
                &type_sfc, Nx_sfc, Ny_sfc, dx_sfc, dy_sfc,&
                &Npar_sfc, parms_Sfc, gndTemp, gndAlbedo)
        end if
        ! Broadcast the surface array data
        call shd_mpi__bcastSfc(type_sfc, Nx_sfc, Ny_sfc,&
            &Npar_sfc, parms_Sfc, gndTemp, gndAlbedo)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read file_ckd and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (isCkd) then
        if (isMasterProc) then
            call shd_io__readCkdSize(file_ckd, waveNum, Ng_ckd, NzCkd_)
        end if
        ! Broadcast the k-distribution array sizes
        call shd_mpi__bcastCkdSize(Ng_ckd, NzCkd_)
        allocate (zCkd_(NzCkd_), GASABS_(NzCkd_), kabs_Ckd(NzCkd_, Ng_ckd))
        allocate (solar_Ckd(Ng_ckd), delg_Ckd(Ng_ckd))
        memUsed = memUsed + NzCkd_*(2.0+Ng_ckd) + 2.0*Ng_ckd

        if (isMasterProc) then
            maxNg = Ng_ckd; maxNz = NzCkd_
            call shd_io__readCkd(maxNg, maxNz, file_ckd, waveNum, solar_Ckd,&
                &Ng_ckd, delg_Ckd, NzCkd_, zCkd_, kabs_Ckd)
        end if
        call shd_mpi__bcastCkd(solar_Ckd, Ng_ckd, delg_Ckd, NzCkd_, zCkd_, kabs_Ckd)
    else
        Ng_ckd = 1
        allocate (solar_Ckd(Ng_ckd), delg_Ckd(Ng_ckd))
        memUsed = memUsed + 2.0*Ng_ckd

        solar_Ckd(1) = 1.0
        delg_Ckd(1) = 1.0
        NzCkd_ = 0
        isBaseOut = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read file_sen and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in user-specified sensor properties if given
    !   this is a binary angle input implemented for 'R' and 'V3'
    if (file_sen(1:4) .ne. 'NONE') then
        if (isMasterProc) then
            call shd_io__readSenSize(file_sen, Npts_sen, Npar_sen)
        end if
        call shd_mpi__bcastSenSize(Npts_sen, Npar_sen)
        allocate (parms_Sen(Npts_sen*Npar_sen))
        memUsed = memUsed + Npts_sen*Npar_sen

        if (isMasterProc) then
            call shd_io__readSen(file_sen, Npts_sen, Npar_sen,&
                &parms_Sen)
        end if
        call shd_mpi__bcastSen(Npts_sen, Npar_sen, parms_Sen)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Check outputs
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Figure out what type of output we will be generating
    Nrad = 0
    Nlos = 0
    isRadOut  = .false.
    isCamOut  = .false.
    isSrcOut  = .false.
    isFluxOut = .false.
    isHROut  = .false.
    isSHOut  = .false.
    do iout = 1, Nout
        if (type_Out(iout) .eq. 'R') then
            isRadOut = .true.
            if (parms_Out(2, iout) .le. 0.0) then
                NxOut = 1
            else
                if (isMpiX_) then
                    NxOut = max(1, nint(XDOMAIN_/parms_Out(2, iout)))
                else
                    NxOut = max(1, nint((XDOMAIN_ + 2*parms_Out(4, iout))/parms_Out(2, iout)))
                end if
            end if
            if (parms_Out(3, iout) .le. 0.0) then
                NyOut = 1
            else
                if (isMpiY_) then
                    NyOut = max(1, nint(YDOMAIN_/parms_Out(3, iout)))
                else
                    NyOut = max(1, nint((YDOMAIN_ + 2*parms_Out(5, iout))/parms_Out(3, iout)))
                end if
            end if
            Nrad = Nrad + Nstokes*NxOut*NyOut*nint(parms_Out(6, iout))
        else if (type_Out(iout) .eq. 'V') then
            isCamOut = .true.
            if (nint(parms_Out(1, iout)) .eq. 1) then
                Nlos = Nlos + Nstokes*nint(parms_Out(10, iout))*nint(parms_Out(11, iout))
            else if (nint(parms_Out(1, iout)) .eq. 2) then
                Nscan = 1 + sqrt((parms_Out(4, iout) - parms_Out(7, iout))**2&
                    + (parms_Out(5, iout) - parms_Out(8, iout))**2&
                    + (parms_Out(6, iout) - parms_Out(9, iout))**2)/parms_Out(10, iout)
                Npix = 1 + abs(parms_Out(12, iout) - parms_Out(11, iout))/parms_Out(13, iout)
                Nlos = Nlos + Nstokes*Nscan*Npix
            ! user-specified sensor
            else if (nint(parms_Out(1, iout)) .eq. 3) then
                Nlos = Nlos + Nstokes*nint(parms_Out(2, iout))
            end if
        else if (type_Out(iout) .eq. 'F') then
            isFluxOut = .true.
        else if (type_Out(iout) .eq. 'H') then
            isHROut = .true.
        else if (type_Out(iout) .eq. 'S') then
            isSHOut = .true.
        else if (type_Out(iout) .eq. 'J') then
            isSrcOut = .true.
        end if
    end do
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Notify user that certain outputs are not supported with multiple processors
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMpi_ .and. isMasterProc) then
        do iout = 1, Nout
            if (type_Out(iout) .eq. 'F') then
                if (nint(parms_Out(1, iout)) .eq. 2 .or. nint(parms_Out(1, iout)) .eq. 5) then
                    if (verbose_ .ge. 4) write (6, '(A,I1,A,A)') 'Warning [shdom]: F output format ', nint(parms_Out(1, iout)),&
                        &' not supported for multiple processors: ', trim(adjustl(file_Out(iout)))
                end if
            else if (type_Out(iout) .eq. 'H') then
                if (nint(parms_Out(1, iout)) .eq. 3) then
                    if (verbose_ .ge. 4) write (6, '(A,I1,A,A)') 'Warning [shdom]: H output format ', nint(parms_Out(1, iout)),&
                        &' not supported for multiple processors: ', trim(adjustl(file_Out(iout)))
                end if
            else if (type_Out(iout) .eq. 'S') then
                if (nint(parms_Out(1, iout)) .eq. 2) then
                    if (verbose_ .ge. 4) write (6, '(A,I1,A,A)') 'Warning [shdom]: S output format ', nint(parms_Out(1, iout)),&
                        &' not supported for multiple processors: ', trim(adjustl(file_Out(iout)))
                end if
            else if (type_Out(iout) .eq. 'J') then
                if (verbose_ .ge. 4) write (6, '(A,A)') 'Warning [shdom]: J output format not supported for multiple processors: ',&
                    &trim(adjustl(file_Out(iout)))
            else if (type_Out(iout) .eq. 'M') then
                if (verbose_ .ge. 4) write (6, '(A,A)') 'Warning [shdom]: M output format not supported for multiple processors: ',&
                    &trim(adjustl(file_Out(iout)))
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Memory check
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMpi_) then
        ! Read the memory parameters from a file, if it exists, and send
        !   the parameters to all processors. Otherwise use the input ones.
        call shd_mpi__bcastMem(MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
            &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, runName)
    end if

    ! Guess maximum number of grid points, cells, SH vector size needed
    !   but don't let MAX_TOTAL_MB be exceeded
    MAXCPUMEM_MB_ = 0.9*TOTSYSMEM_MB_/real(Nproc)
    if (MAX_TOTAL_MB .gt. MAXCPUMEM_MB_) then
        MAX_TOTAL_MB = MAXCPUMEM_MB_
        if (verbose_ .ge. 5) write (6, '(A,F0.2," MB")') 'Message [shdom]: MAX_TOTAL_MB reduced to: ', MAX_TOTAL_MB
    else
        if (verbose_ .ge. 5) write (6, '(A,F0.2," MB")') 'Message [shdom]: MAX_TOTAL_MB is set to: ', MAX_TOTAL_MB
    end if

    if ((splitAcc .le. 1.0d-12) .and. (inSaveFile(1:4) .eq. 'NONE')) then
        ADAPT_GRID_FACTOR = 1.0
    end if

    NbigArr = 3
    if (.not. isAccel) NbigArr = 2

    ! NmaxIG
    !╭──────────────────────────────────────────────────────────────╮!
    NmaxIG = 0.0

    NmaxIG = NmaxIG + 4.5 ! temp, planck, extinct, ssAlb, iPhase[0.5, integer(kind=2)]

    if (NUMPHASE_ .gt. 0) then
        memUsed = memUsed + NUMPHASE_*Nstleg*(maxLeg + 1) ! offset phase table
    else
        NmaxIG = NmaxIG + Nstleg*(maxLeg + 1)
    end if

    if (file_sfc(1:4) .ne. 'NONE') then
        if (type_sfc(2:2) .eq. 'L') then
            NmaxIG = real(3)/real(Nz) * (maxPar_sfc + 2.0 + Nstokes*2.0)
        else
            NmaxIG = real(3)/real(Nz) * (maxPar_sfc + 2.0 + Nstokes*(2.0+ Nmu*Nphi0max/2))
        end if
    end if

    if (isFluxOut)  NmaxIG = NmaxIG + 3.0
    if (isHROut)    NmaxIG = NmaxIG + 1.0
    if (isSHOut)    NmaxIG = NmaxIG + NshOut
    if (isSrcOut)   NmaxIG = NmaxIG + Nstokes

    NmaxIG = NmaxIG + 16.5*CELL_TO_POINT_RATIO ! ptr8Grid_Cell[8.0], ptr6Cell_Cell[6.0], ptr2Tree_Cell[2.0], flags_Cell[0.5]

    NmaxIG = NmaxIG + 6.0   ! xyz_Grid[3.0], RSHptr[1.0], SHptr[1.0], OSHptr[1.0]
    memUsed = memUsed + 4.0 ! offset RSHptr[+2], SHptr[+1], OSHptr[+1]

    if (isAccel) then
        NmaxIG = NmaxIG + Nstokes*(1.0+3.0*NUM_SH_TERM_FACTOR*Nlm) ! radiance, source, delSource
    else
        NmaxIG = NmaxIG + Nstokes*(1.0+2.0*NUM_SH_TERM_FACTOR*Nlm) ! radiance, source
    end if

    NmaxIG = NmaxIG + 3.0 ! fluxes[2.0], dirFlux[1.0]

    memUsed = memUsed + (Nrad+Nlos) ! offset sumRadOut[Nrad] and sumCamOut[Nlos]

    NmaxIG = NmaxIG + (Nstokes*(max(Nphi0max, 4)+1.0) + 8.0) !work, work1, work2

    memUsed = memUsed + Nrad + Nlos ! work2 (over-count but memory safe)
    !╰──────────────────────────────────────────────────────────────╯!

    ! memUsed is now more accurate than before
    ! before:
    !   memUsed = Nmu*(2 + 2*Nphi + 2*Nlm + 2*33*32)&            ! SHDO allocation
    !       + 4.5*maxPG + maxPGL + NUMPHASE_*Nstleg*(maxLeg + 1) ! phase allocation
    ! now: update memUsed after every allocation step
    !
    ! NmaxIG is now more accurate than before
    ! before:
    !   NmaxIG = 28.0 + 16.5*CELL_TO_POINT_RATIO&
    !       + Nphi0max*Nstokes + NUM_SH_TERM_FACTOR*Nstokes*Nlm*NbigArr
    ! now: update NmaxIG from every allocation step

    memNeed = ADAPT_GRID_FACTOR*NptsB*NmaxIG

    if (isHighOrder_) then
        scaleRatio= min(2.0, ((MAX_TOTAL_MB*1024**2)/4 - memUsed)/memNeed)
    else
        scaleRatio= min(1.5, ((MAX_TOTAL_MB*1024**2)/4 - memUsed)/memNeed)
    end if

    ADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR*scaleRatio
    if (ADAPT_GRID_FACTOR .lt. 0.999999999d0) then
        stop 'Error [shdom]: MAX_TOTAL_MB memory limit exceeded with just base grid.'
    end if

    if (scaleRatio .lt. 1.0) then
        if (verbose_ .ge. 5) write (6, '(A,F8.2)') 'Message [shdom]: ADAPT_GRID_FACTOR reduced to:', ADAPT_GRID_FACTOR
    else
        if (verbose_ .ge. 5) write (6, '(A,F8.2)') 'Message [shdom]: ADAPT_GRID_FACTOR increased to:', ADAPT_GRID_FACTOR
    end if

    memNeed = scaleRatio*memNeed
    if (memNeed .gt. huge(maxIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Number of words of memory exceeds max integer size:', memNeed
        stop
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Retrieve max numbers for array allocation after memory check
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! max number of internal grid points
    maxIG = nint(ADAPT_GRID_FACTOR*NptsB)

    ! max number of grid cells, between 1 and 2 times grid points
    maxIC = nint(CELL_TO_POINT_RATIO*maxIG)

    ! max words for source function and radiance arrays
    maxIV = nint(NUM_SH_TERM_FACTOR*Nlm*maxIG)

    ! max words for temporary discrete ordinate azimuthal array
    maxIDO = maxIG*max(Nphi0max, 4)

    if (4.0*(maxIV + maxIG)*Nstokes .gt. huge(maxIV)) then
        write (6, '(A,I12)') 'Error [shdom]: Size of big SH arrays (maxIV) probably exceeds max integer number of bytes:', maxIV
        stop
    end if

    if (4.0*8.0*maxIC .gt. huge(maxIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Size of ptr8Grid_Cell array (8*maxIC) probably exceeds max integer number of bytes:', 8*maxIC
        stop
    end if

    maxNpts_bndry = maxIG*3/Nz
    if (type_sfc(2:2) .eq. 'L') then ! Lambertian surface
        maxNrad_bndry = 2 * maxNpts_bndry
    else
        maxNrad_bndry = (2 + Nmu*Nphi0max/2)*maxNpts_bndry
    end if

    ! max number of legendre coefficients
    if (NUMPHASE_ .gt. 0) then
        maxIGL = NUMPHASE_*Nstleg*(maxLeg + 1)
    else
        maxIGL = maxIG*Nstleg*(maxLeg + 1)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Allocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! temperature and planck radiance (maybe we don't need this for shortwave)
    allocate (temp(maxIG), planck(maxIG))

    ! extinction, single scattering albedo, and phase function indices
    allocate (extinct(maxIG), ssAlb(maxIG), iPhase(maxIG))

    ! legendre coefficients
    allocate (legen(maxIGL))

    ! boundary points pointers and radiances
    !   ptr_Bndry(:, 1): pointers to grid points that are on top
    !   ptr_Bndry(:, 2): pointers to grid points that are on bottom
    !   rad_Bndry(:[Nstokes], 1:NptsTop_bndry): isotropic top radiances
    !   rad_Bndry(:[Nstokes], NptsTop_bndry+1:NptsTop_bndry+NptsBot_bndry): upwelling bottom radiances
    !   rad_Bndry(:[Nstokes], NptsTop_bndry+NptsBot_bndry+1:NptsTop_bndry+NptsBot_bndry+NptsBot_bndry*Nang/2):
    !         downwelling bottom radiances at discrete ordinates
    if (file_sfc(1:4) .ne. 'NONE') allocate (parmsGrid_Sfc(maxPar_sfc*maxNpts_bndry))
    allocate (ptr_Bndry(2*maxNpts_bndry), rad_Bndry(Nstokes, maxNrad_bndry))

    ! output
    if (isFluxOut) allocate (sumFluxes(2, maxIG), sumDirFlux(maxIG))
    if (isHROut)   allocate (sumFluxDiv(maxIG))
    if (isSHOut)   allocate (sumSHterms(maxIG*NshOut))
    if (isSrcOut)  allocate (source1out(maxIG*Nstokes))

    ! pointer to 8 grid points belonging to a cell
    ! pointer to 6 neighboring cells of a cell
    ! pointer to parent and children cells in tree structure
    allocate (ptr8Grid_Cell(8, maxIC), ptr6Cell_Cell(6, maxIC), ptr2Tree_Cell(2, maxIC))

    ! a collection of binary flags for the cells
    allocate (flags_Cell(maxIC))

    ! x,y,z coordinate values of grid point locations
    allocate (xyz_Grid(3, maxIG))

    ! SH related array
    allocate (RSHptr(maxIG + 2), SHptr(maxIG + 1), OSHptr(maxIG + 1))

    ! radiance and source
    allocate (radiance(Nstokes, maxIV + maxIG), source(Nstokes, maxIV))

    ! delta source
    if (isAccel) allocate (delSource(Nstokes, maxIV))

    ! fluxes
    allocate (fluxes(2, maxIG), dirFlux(maxIG))

    ! output radiances to be summed over g
    allocate (sumRadOut(Nrad), sumCamOut(Nlos))

    ! work: work array for path integration; radiance field, Nazimuth*Npts
    allocate (work(Nstokes*maxIDO))

    ! work1: work array for path integration; grid point sweeping order
    allocate (work1(8*maxIG))

    ! work2: work array for path integration; single angle radiance
    allocate (work2(max(Nlos, Nrad, Nstokes*maxIG)))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read in previous run from binary file if desired
    call shd_io__restoreState(inSaveFile, Nstokes, Nx, Ny, Nz, maxNl, maxNm, Nlm,&
        &INRADFLAG, NEWGRIDFLAG, x_Grid, y_Grid, z_Grid,&
        &Npts, Ncells, xyz_Grid, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell,&
        &flags_Cell, fluxes, SHptr, source, RSHptr, radiance)


    ! Check if we can use pre-calculated "state", if not, generate new grids
    !   NEWGRIDFLAG is from shd_io__restoreState call
    if (NEWGRIDFLAG) then
        ! Make the internal base grid lines
        call shd_init__makeGrid(bcFlag, gridType, NPX_, NPY_, NPZ_, Nx, Ny, Nz,&
            &XSTART_, YSTART_, delX_, delY_, ZLEVELS_,&
            &x_Grid, y_Grid, z_Grid)
    end if


    message = 'Message [shdom]: Total time used for preparation: '
    call shd_util__tickTock(message, countStart)
    if (isMasterProc) then
        if (verbose_ .ge. 5) write (6, '(A)') trim(adjustl(message))
    end if


    ! solve radiative transfer equation
    !   1. Use Picard iteration to solve radiance and source
    !   2. Compute actual radiance (and/or other radiative quantities specified by user)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (verbose_ .ge. 5) write (6, '(A,"<",A,">")') 'Message [shdom]: Starting solution iterations for SHDOM run: ', trim(adjustl(runName))

    totIter = 0
    oldNpts = 0
    solar = 0.0d0

    ! Loop over the k-distribution g's from high to low absorption
    CKD_LOOP: do ig_ckd = Ng_ckd, 1, -1

        if (verbose_ .ge. 5) write (6, '(A,I2,A,I2,A)') 'Message [shdom]: Starting ', Ng_ckd-ig_ckd+1, '/', Ng_ckd, ' (ig/Ng) ...'

        ! solar flux for given g of CKD
        solar0_ckd = solar_Ckd(ig_ckd)*solar0*abs(mu_sol)

        ! weighted average solar flux
        solar = solar + solar_Ckd(ig_ckd)*delg_Ckd(ig_ckd)*solar0*abs(mu_sol)

        ! If doing k-distribution then get the gas absorption profile.
        if (isCkd) then
            do iz = 1, NzCkd_
                GASABS_(iz) = kabs_Ckd(iz, ig_ckd)
            end do
        end if

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (isBaseOut .or. ((ig_ckd.eq. Ng_ckd) .and. NEWGRIDFLAG)) then
            call shd_init__makeCell(bcFlag, ipFlag,&
                &Nx, Ny, Nz, Nx1, Ny1, Npts, Ncells,&
                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell)
            NcellsB = Ncells
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in source (SHptr), radiance (RSHptr), fluxes, dirFlux.
        call shd_solv__solveRte(Nstokes, Nx, Ny, Nz, Nx1, Ny1, Nang,&
            &maxNl, maxNm, Nlm, Nmu, Nphi, Nleg, Nstleg, NUMPHASE_,&
            &Nphi0, mu, phi, wtDO,&
            &maxIV, maxIC, maxIG, maxIDO, INRADFLAG,&
            &bcFlag, ipFlag, isDeltaM, srcType, isHighOrder_,&
            &solar0_ckd, mu_sol, az_sol, skyRad,&
            &type_sfc, gndTemp, gndAlbedo,&
            &Nx_sfc, Ny_sfc, dx_sfc, dy_sfc,&
            &Npar_sfc, parms_Sfc, parmsGrid_Sfc,&
            &units, waveNum, waveLen,&
            &isAccel, solAcc, maxIter, solCrit, iter,&
            &splitAcc, shAcc, x_Grid, y_Grid, z_Grid,&
            &temp, planck, extinct, ssAlb, legen, iPhase,&
            &maxNpts_bndry, maxNrad_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
            &Npts, xyz_Grid,&
            &Ncells, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
            &RSHptr, SHptr, OSHptr, work, work1, work2,&
            &source, delSource, radiance, fluxes, dirFlux)
        INRADFLAG = .false.
        totIter = totIter + iter

        ! Compute the output and put in arrays. First the radiance.
        if (isRadOut) then
            irad = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'R') then
                    if (isMpi_) then
                        if (file_sen(1:4) .ne. 'NONE') then
                            call shd_mpi__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0mAX, Nphi0, mu, phi, wtDO,&
                                &bcFlag, XDOMAIN_, YDOMAIN_, ipFlag,&
                                &srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, work, work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2, parms_Sen)
                        else
                            call shd_mpi__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, XDOMAIN_, YDOMAIN_, ipFlag,&
                                &srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, work, work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2)
                        end if
                    else
                        if (file_sen(1:4) .ne. 'NONE') then
                            call shd_rad__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, ipFlag, srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, work(1), work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2, parms_Sen)
                        else
                            call shd_rad__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, ipFlag, srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, work(1), work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2)
                        end if
                    end if
                end if
            end do
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Nrad, work2, sumRadOut)
        end if

        if (isCamOut) then
            ilos = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'V') then
                    if (isMpi_) then
                        if (file_sen(1:4) .ne. 'NONE') then
                            call shd_mpi__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, XDOMAIN_, YDOMAIN_, ipFlag,&
                                &srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, parms_Out(1, iout), ilos, work2, parms_Sen)
                        else
                            call shd_mpi__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, XDOMAIN_, YDOMAIN_, ipFlag,&
                                &srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, parms_Out(1, iout), ilos, work2)
                        end if
                    else
                        if (file_sen(1:4) .ne. 'NONE') then
                            call shd_rad__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, ipFlag, srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, parms_Out(1, iout), ilos, work2, parms_Sen)
                        else
                            call shd_rad__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                                &maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_,&
                                &Nmu, Nphi0max, Nphi0, mu, phi, wtDO,&
                                &bcFlag, ipFlag, srcType, isDeltaM, mu_sol, az_sol,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bndry, NptsTop_bndry, NptsBot_bndry, ptr_Bndry, rad_Bndry,&
                                &gndTemp, gndAlbedo, skyRad, waveNum, waveLen, units,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell,&
                                &extinct, ssAlb, legen, iPhase, dirFlux, fluxes,&
                                &SHptr, source, parms_Out(1, iout), ilos, work2)
                        end if
                    end if
                end if
            end do
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Nlos, work2, sumCamOut)
        end if

        if (isSrcOut) then
            ! Compute source function throughout grid for this one angle
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'J') then
                    mu_out  = parms_Out(2, iout)
                    phi_out = parms_Out(3, iout)*acos(-1.0)/180.0
                    call shd_src__computeSrcOneDir(Nstokes, maxNl, maxNm, Nlm,&
                        &Nstleg, Nleg, NUMPHASE_,&
                        &Npts, isDeltaM, mu_out, phi_out,&
                        &srcType, mu_sol, az_sol, ssAlb, legen, iPhase,&
                        &dirFlux, SHptr, source, source1out)
                end if
            end do
        end if

        if (isBaseOut) then
            ! If outputting only the base grid then throw out the other points
            call shd_init__makeCell(bcFlag, ipFlag,&
                &Nx, Ny, Nz, Nx1, Ny1,&
                &Npts, Ncells, x_Grid, y_Grid, z_Grid, xyz_Grid,&
                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flags_Cell)
        else
            ! For the output from the grid point arrays (Flux, Net flux, SH)
            !   we may need to interpolate the running k-distribution sum
            !   to the newly created grid points.
            call shd_util__interpData(oldNpts, Npts, isFluxOut, isHROut, isSHOut,&
                sumFluxes, sumDirFlux, sumFluxDiv, sumSHterms,&
                NcellsB, Ncells, ptr2Tree_Cell, ptr8Grid_Cell, xyz_Grid)
        end if

        if (isFluxOut) then

            ! unscale delta-m fluxes
            if ((srcType .ne. 'T') .and. isDeltaM) then

                allocate(dirFlux_tmp(Npts))

                if (isMpi_) then
                    call shd_mpi__makeDirect(1, Npts, bcFlag, ipFlag, .false., maxNl,&
                        &Nstleg, Nleg, solar0_ckd, mu_sol, az_sol,&
                        &xyz_Grid, Nx, x_Grid, Ny, y_Grid, dirFlux_tmp)
                else
                    call shd_rad__makeDirect(Npts, bcFlag, ipFlag, .false., maxNl,&
                        &Nstleg, Nleg, solar0_ckd, mu_sol, az_sol,&
                        &xyz_Grid, dirFlux_tmp)
                end if

                fluxes(1, 1:Npts) = fluxes(1, 1:Npts) + dirFlux(1:Npts) - dirFlux_tmp(1:Npts)
                dirFlux(1:Npts) = dirFlux_tmp(1:Npts)

                deallocate(dirFlux_tmp)

            end if

            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), 2*Npts, fluxes, sumFluxes)
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Npts, dirFlux, sumDirFlux)
        end if

        if (isHROut) then
            ! Compute the net flux divergence at every point
            call shd_rad__computeNetFluxDiv(Nstokes, Npts, RSHptr, srcType, mu_sol,&
                &extinct, ssAlb, planck, dirFlux, radiance, work)
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Npts, work, sumFluxDiv)
        end if

        if (isSHOut) then
            call shd_shdo__computeSH(NshOut, Nstokes, Npts, srcType, mu_sol, az_sol,&
                &dirFlux, RSHptr, maxNl, maxNm, radiance, work)
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), NshOut*Npts, work, sumSHterms)
        end if

        oldNpts = Npts

    end do CKD_LOOP
    !╰────────────────────────────────────────────────────────────────────────────╯!


    message = 'Message [shdom]: Total time used for runs: '
    call shd_util__tickTock(message, countStart)
    if (isMasterProc) then
        if (verbose_ .ge. 5) write (6, '(A)') trim(adjustl(message))
    end if


    ! Output radiance and/or flux (and/or other radiative quantities)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    Nsh = SHptr(Npts + 1)
    irad = 0
    ilos = 0
    if (isMpi_) then

        if (isMasterProc) then
            if (isFluxOut) allocate (allFluxes(3, Nz, NyT, NxT))
            if (isHROut)   allocate (allFluxDiv(Nz, NyT, NxT))
            if (isSHOut)   allocate (allSHterms(NshOut, Nz, NyT, NxT))
            allocate (xAllGrid(NxT), yAllGrid(NyT))
        end if

        call shd_mpi__gatherOut(isFluxOut, isHROut, isSHOut,&
            &ipFlag, bcFlag, NptsB, NxT, NyT, Nz, NshOut,&
            &sumFluxes, sumDirFlux, sumFluxDiv, sumSHterms,&
            &XDOMAIN_, YDOMAIN_, xAllGrid, yAllGrid,&
            &allFluxes, allFluxDiv, allSHterms,&
            &Ncells, Npts, Nsh, NcellsTot, NptsTot, NshTot)

        if (isMasterProc) then

            do iout = 1, Nout
                if (type_Out(iout) .ne. 'V' .and. file_Out(iout)(1:4) .ne. 'NONE') then
                    if (file_sen(1:4) .ne. 'NONE') then
                        call shd_io__outputResultMpi(Nstokes, NxT, NyT, Nz,&
                            &NptsTot, NcellsTot, NshTot,&
                            &maxNl, maxNm, Nlm, Nmu, Nphi, Nang, Ng_ckd,&
                            &file_prp, file_sfc, file_ckd, inSaveFile, outSaveFile,&
                            &bcFlag, ipFlag, isDeltaM, gridType, srcType,&
                            &solar, mu_sol, az_sol, skyRad,&
                            &type_sfc, gndTemp, gndAlbedo, waveNum, waveLen, units,&
                            &splitAcc, shAcc, solAcc, maxIter, totIter,&
                            &XDOMAIN_, YDOMAIN_, xAllGrid, yAllGrid, z_Grid,&
                            &allFluxes, allFluxDiv,&
                            &NshOut, allSHterms, irad, sumRadOut,&
                            &type_Out(iout), parms_Out(1, iout), file_Out(iout), parms_Sen)
                    else
                        call shd_io__outputResultMpi(Nstokes, NxT, NyT, Nz,&
                            &NptsTot, NcellsTot, NshTot,&
                            &maxNl, maxNm, Nlm, Nmu, Nphi, Nang, Ng_ckd,&
                            &file_prp, file_sfc, file_ckd, inSaveFile, outSaveFile,&
                            &bcFlag, ipFlag, isDeltaM, gridType, srcType,&
                            &solar, mu_sol, az_sol, skyRad,&
                            &type_sfc, gndTemp, gndAlbedo, waveNum, waveLen, units,&
                            &splitAcc, shAcc, solAcc, maxIter, totIter,&
                            &XDOMAIN_, YDOMAIN_, xAllGrid, yAllGrid, z_Grid,&
                            &allFluxes, allFluxDiv,&
                            &NshOut, allSHterms, irad, sumRadOut,&
                            &type_Out(iout), parms_Out(1, iout), file_Out(iout))
                    end if
                end if
            end do

            if (isFluxOut) deallocate (allFluxes)
            if (isHROut)   deallocate (allFluxDiv)
            if (isSHOut)   deallocate (allSHterms)
            deallocate (xAllGrid, yAllGrid)
        end if

    else

        do iout = 1, Nout
            if ((type_Out(iout) .ne. 'V') .and. (file_Out(iout)(1:4) .ne. 'NONE')) then
                if (file_sen(1:4) .ne. 'NONE') then
                    call shd_io__outputResult(Nstokes, Nx, Ny, Nx1, Ny1, Nz, NptsB, Npts, Ncells,&
                        Nsh, maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_, Nmu, Nphi, Nang, Ng_ckd,&
                        file_prp, file_sfc, file_ckd, inSaveFile, outSaveFile,&
                        bcFlag, ipFlag, isDeltaM, gridType, srcType,&
                        solar, mu_sol, az_sol, skyRad,&
                        type_sfc, gndTemp, gndAlbedo, waveNum, waveLen, units,&
                        splitAcc, shAcc, solAcc, maxIter, totIter,&
                        x_Grid, y_Grid, z_Grid, xyz_Grid,&
                        ptr2Tree_Cell, ptr8Grid_Cell, flags_Cell,&
                        extinct, ssAlb, legen, iPhase, temp,&
                        sumFluxes, sumDirFlux, sumFluxDiv,&
                        irad, sumRadOut, NshOut, sumSHterms, source1out,&
                        type_Out(iout), parms_Out(1, iout), file_Out(iout), parms_Sen)
                else
                    call shd_io__outputResult(Nstokes, Nx, Ny, Nx1, Ny1, Nz, NptsB, Npts, Ncells,&
                        Nsh, maxNl, maxNm, Nlm, Nstleg, Nleg, NUMPHASE_, Nmu, Nphi, Nang, Ng_ckd,&
                        file_prp, file_sfc, file_ckd, inSaveFile, outSaveFile,&
                        bcFlag, ipFlag, isDeltaM, gridType, srcType,&
                        solar, mu_sol, az_sol, skyRad,&
                        type_sfc, gndTemp, gndAlbedo, waveNum, waveLen, units,&
                        splitAcc, shAcc, solAcc, maxIter, totIter,&
                        x_Grid, y_Grid, z_Grid, xyz_Grid,&
                        ptr2Tree_Cell, ptr8Grid_Cell, flags_Cell,&
                        extinct, ssAlb, legen, iPhase, temp,&
                        sumFluxes, sumDirFlux, sumFluxDiv,&
                        irad, sumRadOut, NshOut, sumSHterms, source1out,&
                        type_Out(iout), parms_Out(1, iout), file_Out(iout))
                end if
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Output camera radiance if there are any
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        do iout = 1, Nout
            if ((type_Out(iout) .eq. 'V')) then
                if (file_sen(1:4) .ne. 'NONE') then
                    call shd_io__outputCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                        &Nsh, maxNl, maxNm, Nlm, Nmu, Nphi, Nang, Ng_ckd,&
                        &file_prp, file_sfc, file_ckd, inSaveFile, outSaveFile,&
                        &bcFlag, ipFlag, isDeltaM, gridType, srcType,&
                        &solar, mu_sol, az_sol, skyRad,&
                        &type_sfc, gndTemp, gndAlbedo, waveNum, waveLen, units,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &ilos, sumCamOut, parms_Out(1, iout), file_Out(iout), parms_Sen)
                else
                    call shd_io__outputCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncells,&
                        &Nsh, maxNl, maxNm, Nlm, Nmu, Nphi, Nang, Ng_ckd,&
                        &file_prp, file_sfc, file_ckd, inSaveFile, outSaveFile,&
                        &bcFlag, ipFlag, isDeltaM, gridType, srcType,&
                        &solar, mu_sol, az_sol, skyRad,&
                        &type_sfc, gndTemp, gndAlbedo, waveNum, waveLen, units,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &ilos, sumCamOut, parms_Out(1, iout), file_Out(iout))
                end if
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! save state
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (maxIter .gt. 0) then
        call shd_io__saveState(outSaveFile, Nstokes, Nx, Ny, Nz, maxNl, maxNm, Nlm,&
            &work, x_Grid, y_Grid, z_Grid,&
            &Npts, Ncells, xyz_Grid, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell,&
            &flags_Cell, fluxes, SHptr, source, RSHptr, radiance)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    message = 'Message [shdom]: Total time used for output: '
    call shd_util__tickTock(message, countStart)
    if (isMasterProc) then
        if (verbose_ .ge. 5) write (6, '(A)') trim(adjustl(message))
    end if


    ! Calculate the actual memory parameters and output them to the log file
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! MAXMB is now more accurate than before
    ! before:
    !  MAXMB_OUT = 4*(Nmu*(2 + 2*Nphi + 2*Nlm + 2*33*32)&
    !      + 4.5*maxPG + maxPGL + NUMPHASE_*Nstleg*(maxLeg + 1)&
    !      + 16.5*Ncells + Npts*(28 + Nphi0max*Nstokes)&
    !      + Nstokes*Nsh*NbigArr)/1024**2
    ! now:
    !   use NmaxIG
    MAXMB_OUT = 4*(memUsed + 16.5*Ncells&
        + Npts*(NmaxIG-16.5*CELL_TO_POINT_RATIO))/1024**2

    ADAPTGRIDFAC_OUT = float(Npts)/NptsB
    SHTERM_FAC_OUT = float(Nsh)/(Nlm*Npts)
    CELL_POINT_OUT = float(Ncells)/Npts
    if (verbose_ .ge. 5) write (6, '(A)')'Message [shdom]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    if (verbose_ .ge. 5) write (6, '(26X,F9.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Deallocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    deallocate (ZLEVELS_, TEMPP_, EXTINCTP_, ALBEDOP_, EXTDIRP_, IPHASEP_, LEGENP_)
    if (isFluxOut) deallocate (sumFluxes, sumDirFlux)
    if (isHROut)   deallocate (sumFluxDiv)
    if (isSHOut)   deallocate (sumSHterms)
    if (isSrcOut)  deallocate (source1out)
    if (isCkd) deallocate (zCkd_, GASABS_, kabs_Ckd)
    if (allocated(parms_Sfc)) deallocate (parms_Sfc, parmsGrid_Sfc)
    if (allocated(parms_Sen)) deallocate (parms_Sen)
    deallocate (delg_Ckd, solar_Ckd)
    deallocate (ptr_Bndry, rad_Bndry)

    deallocate (mu, wtDO, phi, Nphi0)
    deallocate (x_Grid, y_Grid, z_Grid)
    deallocate (temp, planck)
    deallocate (extinct, ssAlb, iPhase)
    deallocate (legen)
    deallocate (ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell)
    deallocate (flags_Cell)
    deallocate (xyz_Grid)
    deallocate (RSHptr, SHptr, OSHptr)
    deallocate (radiance, source)
    if (isAccel) deallocate (delSource)
    deallocate (fluxes, dirFlux)
    deallocate (sumRadOut, sumCamOut)
    deallocate (work, work1, work2)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    call shd_mpi__end(Npts, xyz_Grid, NpX_, NPY_, XSTART_, YSTART_, XDOMAIN_, YDOMAIN_,&
        &NpxT, NpyT, file_prp, runName)

end program shdom
