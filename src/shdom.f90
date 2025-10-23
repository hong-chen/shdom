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
    use globals_cfg,   only: verbose_, charMaxLen_, accSolv_, &
                             memTotSysMB_, memPerCpuMB_, &
                             isNmlInp_, isHighOrder_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, isMpi_, &
                             is1DX_, is1DYx_, is1D_, is3Ddir1Ddif_, &
                             is1D1_, isPol_, isPol4_, &
                             isSolar_, isThermal_

    use globals_prp,   only: Nx_prp_, Ny_prp_, Nz_prp_, Npha_prp_, &
                             dx_prp_, dy_prp_, xStart_prp_, yStart_prp_, &
                             xDomain_prp_, yDomain_prp_, &
                             z_Prp_, temp_Prp_, &
                             ext_Prp_, ssa_Prp_, &
                             extDir_Prp_, &
                             leg_Prp_, ipha_Prp_

    use globals_ckd,   only: Nz_ckd_, z_Ckd_, kabs0_Ckd_

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
    real(kind=R__)   :: mu_out, phi_out
    real(kind=R__)   :: parms_Out(maxNpar_out, maxNset_out)
    character(len=1) :: type_Out(maxNset_out)
    character(len=charMaxLen_)  :: file_Out(maxNset_out)
    real(kind=R__), allocatable :: sumFluxTot_Out(:, :), sumFluxDir_Out(:), sumFluxDiv_Out(:)
    real(kind=R__), allocatable :: sumSH_Out(:), src1_Out(:)
    real(kind=R__), allocatable :: sumRad_Out(:), sumCam_Out(:)
    logical :: isBaseOut, isRadOut, isCamOut, isFluxOut, isHROut, isSHOut, isSrcOut

    character(len=charMaxLen_) :: fileOut_nc
    !╰────────────────────────────────────────────────────────────────────────────╯!


    real(kind=R__), allocatable  :: planck(:)
    real(kind=R__), allocatable  :: temp(:), ext(:), ssa(:)
    real(kind=R__), allocatable  :: leg(:)
    integer(kind=2), allocatable :: ipha(:)

    real(kind=R__), allocatable  :: temp_PT(:), ext_PT(:), ssa_PT(:)
    real(kind=R__), allocatable  :: leg_PT(:)
    integer(kind=2), allocatable :: ipha_PT(:)

    real(kind=R__), allocatable    :: work(:)
    integer(kind=I__), allocatable :: work1(:)
    real(kind=R__), allocatable    :: work2(:)

    real(kind=R__), allocatable :: fluxTot(:, :), fluxDir(:)
    real(kind=R__), allocatable :: fluxDirTmp(:)

    real(kind=R__), allocatable :: fluxTot_all(:, :, :, :), fluxDiv_all(:, :, :)
    real(kind=R__), allocatable :: SH_all(:, :, :, :)
    real(kind=R__), allocatable :: xGrid_all(:), yGrid_all(:)

    ! SH and DO
    !╭────────────────────────────────────────────────────────────────────────────╮!
    integer(kind=I__) :: maxNl_sh, maxNm_sh, Nlm_sh
    integer(kind=I__) :: Nmu_do, Nphi_do, maxNphi0_do, Nang_do
    integer(kind=I__), allocatable :: Nphi0_Do(:)
    real(kind=R__), allocatable    :: mu_Do(:), phi_Do(:), wt_Do(:)
    real(kind=R__), allocatable    :: rad_Sh(:, :), src_Sh(:, :), srcDel_Sh(:, :)
    integer(kind=I__), allocatable :: ptrRad_Sh(:), ptrSrc_Sh(:), ptrSrcDel_Sh(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SAV
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_sav, fileOut_sav
    logical :: isNewRad, isNewGrid
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! PRP
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_prp
    character(len=1) :: type_prp
    integer(kind=I__) :: NlegSet_prp, Nleg_prp
    real(kind=R__)    :: dxTmp_prp, dyTmp_prp

    ! <globals_prp> contains the following:
    !
    ! integer, save :: Nx_prp_, Ny_prp_, Nz_prp_, Npha_prp_
    ! real, save    :: dx_prp_, dy_prp_, xStart_prp_, yStart_prp_
    ! real, save    :: xDomain_prp_, yDomain_prp_
    ! real, save, allocatable :: z_Prp_(:)
    ! real, save, allocatable :: temp_Prp_(:), ext_Prp_(:), ssa_Prp_(:)
    ! real, save, allocatable :: extDir_Prp_(:)
    ! real, save, allocatable :: leg_Prp_(:)
    ! integer(2), save, allocatable :: ipha_Prp_(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SFC
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_sfc
    character(len=2)  :: type_sfc
    integer(kind=I__) :: Nx_sfc, Ny_sfc, Npar_sfc
    integer(kind=I__) :: maxNpts_sfc, maxNpar_sfc
    real(kind=R__)    :: dx_sfc, dy_sfc
    real(kind=R__)    :: temp0_sfc, alb0_sfc
    real(kind=R__), allocatable :: parms_Sfc(:), &
                                   parmsGrid_Sfc(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! CKD
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_ckd
    integer(kind=I__) :: Ng_ckd, ig_ckd
    real(kind=R__)    :: solar0_ckd
    real(kind=R__), allocatable :: solar_Ckd(:), &
                                   kabs_Ckd(:, :), &
                                   delg_Ckd(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SEN
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_sen
    integer(kind=I__) :: Npts_sen, Npar_sen
    real(kind=R__), allocatable :: parms_Sen(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Boundary (borrowed BC [boundary condition] as postfix)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    integer(kind=I__) :: NptsTop_bc, NptsBot_bc
    integer(kind=I__) :: maxNpts_bc, maxNrad_bc
    integer(kind=I__), allocatable :: ptr_Bc(:)
    real(kind=R__), allocatable    :: rad_Bc(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! grids and cells
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=1) :: type_grid
    real(kind=R__), allocatable    :: x_Grid(:), y_Grid(:), z_Grid(:)
    real(kind=R__), allocatable    :: xyz_Grid(:, :)

    integer(kind=2), allocatable   :: flag_Cell(:)
    integer(kind=I__), allocatable :: ptr8Grid_Cell(:, :), &
                                      ptr6Cell_Cell(:, :), &
                                      ptr2Tree_Cell(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(kind=R__) :: NmaxIG
    real(kind=R__) :: memUsed, memNeed, scaleFactor

    integer(kind=I__) :: maxG_P, maxIG, maxIC, maxIV, maxIDO, maxNz, maxNg
    integer(kind=I__) :: maxLeg, maxGL_P, maxIGL, NbigArr

    real(kind=R__) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real(kind=R__) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real(kind=R__) :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer(kind=I__) :: Nstokes, Nx, Ny, Nz, Nx1, Ny1
    integer(kind=I__) :: Nx_T, Ny_T, Nx_PT, Ny_PT, maxG_PT, maxGL_PT

    integer(kind=I__) :: Npts_base, Npts, Ncell, Nsh, Ncell_base, oldNpts
    integer(kind=I__) :: maxIter, iter, totIter, iz, iout, Nout, flag_bc, ipFlag
    integer(kind=I__) :: Nrad, irad, Nx_out, Ny_out, Nlos, ilos, Nscan, Npix, Nsh_out
    integer(kind=I__) :: Ncell_tot, Npts_tot, Nsh_tot

    character(len=1) :: type_src, units_src
    real(kind=R__) :: solar0, solar, mu_solar, az_solar

    real(kind=R__) :: skyRad, waveNum(2), waveLen
    real(kind=R__) :: solAcc, solCrit, splitAcc, shAcc, maxAsym

    integer(kind=I__) :: Nproc
    logical :: isMasterProc
    character(len=charMaxLen_) :: runName

    integer(kind=I__) :: countStart, countRate, countMax
    character(len=charMaxLen_) :: message

    logical :: isCkd, isDeltaM, isAccel


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
            call shd_io__nmlInp(runName, fileInp_prp, fileInp_sfc, fileInp_ckd,&
                fileInp_sav, fileOut_sav, fileInp_sen,&
                Nstokes, Nx, Ny, Nz, Nmu_do, Nphi_do, flag_bc,&
                ipFlag, isCkd, isDeltaM, type_grid,&
                type_src, solar0, mu_solar, az_solar, skyRad,&
                temp0_sfc, alb0_sfc, units_src, waveNum, waveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        else
            call shd_io__userInp(runName, fileInp_prp, fileInp_sfc, fileInp_ckd,&
                fileInp_sav, fileOut_sav, fileInp_sen,&
                Nstokes, Nx, Ny, Nz, Nmu_do, Nphi_do, flag_bc,&
                ipFlag, isCkd, isDeltaM, type_grid,&
                type_src, solar0, mu_solar, az_solar, skyRad,&
                temp0_sfc, alb0_sfc, units_src, waveNum, waveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        end if
    end if

    call shd_mpi__bcastInp(runName, fileInp_prp, fileInp_sfc, fileInp_ckd,&
        &fileInp_sav, fileOut_sav, fileInp_sen,&
        &Nstokes, Nx, Ny, Nz, Nmu_do, Nphi_do, flag_bc,&
        &ipFlag, isCkd, isDeltaM, type_grid,&
        &type_src, solar0, mu_solar, az_solar, skyRad,&
        &temp0_sfc, alb0_sfc, units_src, waveNum, waveLen,&
        &isAccel, solAcc, maxIter, splitAcc, shAcc,&
        &maxNset_out, maxNpar_out, Nout, type_Out, parms_Out,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Memory used for initialization (rough estimate, 1 WORD = 4 bytes)
    memUsed = 2000.0 + maxNset_out * (maxNpar_out + 0.25 + (charMaxLen_/4.0))


    ! Post-process parameters from input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! set RT solving accuracy (1.0e-4 to 1.0e-8)
    accSolv_ = max(REPS_, min(solAcc*0.1, 1.0e-4))

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    Nstokes = max(1, min(4, Nstokes))
    if (Nstokes .eq. 2) then
        Nstokes = 3
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: Nstokes cannot be 2, reset Nstokes to 3.'
    end if

    ! NlegSet_prp is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (Nstokes .eq. 1) then
        NlegSet_prp = 1
        isPol_ = .false.
    else
        NlegSet_prp = 6
        isPol_ = .true.
    end if

    if (Nstokes .eq. 4) then
        isPol4_ = .true.
    else
        isPol4_ = .false.
    end if

    ! Make maxNl_sh and maxNm_sh from Nmu_do and Nphi_do (for discrete ordinate)
    !   Nmu_do: number of discrete ordinate zenith angles in both hemispheres
    !   Nphi_do: number of discrete ordinate azimuth angles in 2*pi
    Nmu_do = max(2, 2*int((Nmu_do + 1)/2)) ! force Nmu_do an even number
    Nphi_do = max(1, Nphi_do)

    ! maxNl_sh: The maximum L of the spherical harmonic truncation
    ! maxNm_sh: The maximum M of the spherical harmonic truncation
    maxNl_sh = Nmu_do - 1
    maxNm_sh = max(0, int(Nphi_do/2) - 1)

    ! Compute Nlm_sh, maxNphi0_do.
    !   Nlm_sh: The number of spherical harmonic terms
    !   maxNphi0_do: The maximum number of azimuth angles actually used
    Nlm_sh = (2*maxNm_sh + 1)*(maxNl_sh + 1) - maxNm_sh*(maxNm_sh + 1)
    maxNphi0_do = Nphi_do

    ! Nleg_prp: The number of Legendre series phase function terms used in the solution procedure
    Nleg_prp = maxNl_sh
    if (isDeltaM) Nleg_prp = maxNl_sh + 1

    if (type_src .eq. 'B') then
        isSolar_ = .true.
        isThermal_ = .true.
    else
        if (type_src .eq. 'S') then
            isSolar_ = .true.
            isThermal_ = .false.
        else if (type_src.eq. 'T') then
            isSolar_ = .false.
            isThermal_ = .true.
        else
            stop 'Error [shdom]: Unknown type_src.'
        end if
    end if

    ! for thermal, reset solar azimuth angle
    if (type_src.eq. 'T') az_solar = 0.0d0

    ! angle grids for DO
    allocate (mu_Do(Nmu_do), phi_Do(Nmu_do*Nphi_do), Nphi0_Do(Nmu_do), wt_Do(Nmu_do*Nphi_do))

    memUsed = memUsed + 2.0*Nmu_do*(1.0 + Nphi_do)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read fileInp_prp (property file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        ! Find out the size of the input property file grid
        call shd_io__readPrpSize(fileInp_prp, NlegSet_prp, Nleg_prp, Nx_prp_, Ny_prp_, Nz_prp_, dx_prp_, dy_prp_,&
            &Npha_prp_, maxLeg, maxGL_P)
    end if
    call shd_mpi__bcastPrpSize(Nx_prp_, Ny_prp_, Nz_prp_, dx_prp_, dy_prp_,&
        &Npha_prp_, maxLeg, maxGL_P)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read fileInp_sfc (surface file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (fileInp_sfc(1:4) .eq. 'NONE') then
        type_sfc = 'FL'
        Nx_sfc = 1
        Ny_sfc = 1
        dx_sfc = dx_prp_
        dy_sfc = dy_prp_
    else
        if (isMasterProc) then
            ! find out the size of the surface file grid
            call shd_io__readSfcSize(fileInp_sfc, Nx_sfc, Ny_sfc, dx_sfc, dy_sfc, maxNpts_sfc, maxNpar_sfc)
        end if
        ! Broadcast the surface array sizes
        call shd_mpi__bcastSfcSize(Nx_sfc, Ny_sfc, dx_sfc, dy_sfc, maxNpts_sfc, maxNpar_sfc)
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
    Nx_T = Nx  ; Ny_T = Ny
    Nx_PT = Nx_prp_; Ny_PT = Ny_prp_

    if (Nx_PT .eq. 1) then
        dx_prp_ = Nx_sfc*dx_sfc! overwrite dx_prp_ if only 1 grid in X direction
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: Nx_PT=1, reset dx_prp_ to match surface domain.'
    end if
    if (Ny_PT .eq. 1) then
        dy_prp_ = Ny_sfc*dy_sfc! overwrite dy_prp_ if only 1 grid in Y direction
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: Ny_PT=1, reset dy_prp_ to match surface domain.'
    end if

    xDomain_prp_ = Nx_PT*dx_prp_
    yDomain_prp_ = Ny_PT*dy_prp_

    ! Choose the option to compute higher order radiance SH terms
    if ((Nx_PT .eq. 1) .and. (Ny_PT .eq. 1) .and.&
        &(Nx_sfc.le. 32) .and. (Nx_sfc.le. 32)) then
        isHighOrder_ = .true.
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: 1D calculation, set isHighOrder_ to .TRUE.'
    else
        isHighOrder_ = .false.
    end if

    Nsh_out = 4
    if (isHighOrder_) Nsh_out = 5

    ! x, y, z grids
    allocate (x_Grid(Nx + 1), y_Grid(Ny + 1), z_Grid(Nz))
    memUsed = memUsed + (Nx+1) + (Ny+1) + (Nz)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! If MPI, divide domain into mini-domains
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and Nx,Ny to the subdomain size, and output
    !   the X and Y starting position of this subdomain. Set bits
    !   2 and 3 in flag_bc to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    call shd_mpi__map(flag_bc, Nx_prp_, Ny_prp_, Nx, Ny, Nz, dx_prp_, dy_prp_, fileInp_prp,&
        &xStart_prp_, yStart_prp_, runName)

    ! flag_bc bit-representation '...0001' isOpenX_
    if (btest(flag_bc, 0)) then
        isOpenX_ = .true.
    else
        isOpenX_ = .false.
    end if

    ! flag_bc bit-representation '...001x' isOpenYx_
    if (btest(flag_bc, 1)) then
        isOpenYx_ = .true.
    else
        isOpenYx_ = .false.
    end if

    ! flag_bc bit-representation '...01xx' isMpiX_
    if (btest(flag_bc, 2)) then
        isMpiX_ = .true.
    else
        isMpiX_ = .false.
    end if

    ! flag_bc bit-representation '...1xxx' isMpiY_
    if (btest(flag_bc, 3)) then
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
    !   btest(flag_bc, 0) checks whether flag_bc=1 or flag_bc=0 ('...0001' [open in X] or '...0000' [periodic in X])
    Nx1 = Nx + 1
    if (is1DX_ .or. isOpenX_ .or. isMpiX_) Nx1 = Nx

    ! btest(flag_bc, 1) checks whether flag_bc=2 [or 3] or flag_bc=0 [or 1]
    !   ('...0010' [or '...0011', open in Y] or '...0000' [or '...0001', periodic in Y])
    Ny1 = Ny + 1
    if (is1DYx_ .or. isOpenYx_ .or. isMpiY_) Ny1 = Ny

    ! Number of base points
    Npts_base = Nx1*Ny1*Nz

    ! Calculate the number of base grid cells depending on the boundary condition flag
    Ncell_base = (Nz - 1)*(Nx + ibits(flag_bc, 0, 1) - ibits(flag_bc, 2, 1))&
        *(Ny + ibits(flag_bc, 1, 1) - ibits(flag_bc, 3, 1))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read fileInp_prp and distribute data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Determine the size of the property grid arrays so can allocate
    maxNz  = Nz_prp_
    maxG_P = Nx_prp_*Ny_prp_*Nz_prp_
    if (Npha_prp_ .eq. 0) maxGL_P = maxG_P*(maxLeg + 1)*NlegSet_prp

    allocate (z_Prp_(Nz_prp_))
    allocate (temp_Prp_(maxG_P), ext_Prp_(maxG_P), ssa_Prp_(maxG_P))
    allocate (extDir_Prp_(maxG_P), ipha_Prp_(maxG_P), leg_Prp_(maxGL_P))
    memUsed = memUsed + Nz_prp_+ 6.0*maxG_P + maxGL_P

    if (isMpi_) then

        ! Read the properties of the medium to the master processor and
        !   distribute to all the processors
        if (isMasterProc) then

            maxG_PT = Nx_PT*Ny_PT*Nz_prp_
            maxGL_PT = maxGL_P
            if (Npha_prp_ .eq. 0) maxGL_PT = maxG_PT*(maxLeg + 1)*NlegSet_prp

            allocate (temp_PT(maxG_PT), ext_PT(maxG_PT), ssa_PT(maxG_PT))
            allocate (ipha_PT(maxG_PT), leg_PT(maxGL_PT))

            ! Only use delX and delY from initial read, dxTmp_prp and dyTmp_prp are discarded
            call shd_io__readPrp(fileInp_prp, Nx_PT, Ny_PT, Nz_prp_,&
                &NlegSet_prp, maxLeg, Nleg_prp, maxNz, maxG_PT, maxGL_PT, isDeltaM,&
                &type_prp, dxTmp_prp, dyTmp_prp, z_Prp_, maxAsym,&
                &temp_PT, ext_PT, ssa_PT, leg_PT, Npha_prp_, ipha_PT)

            if (Npha_prp_ .eq. 0) then
                stop 'Error [shdom]: Standard property file (with phase function Legendre series at every point)&
                    & is not allowed with multiple processors.'
            end if

            ! Check values of medium properties
            call shd_io__checkPrp(Nx_PT, Ny_PT, Nz_prp_, NlegSet_prp, Nleg_prp,&
                &dx_prp_, dy_prp_, z_Prp_,&
                &temp_PT, ext_PT, ssa_PT, leg_PT, Npha_prp_, ipha_PT)

        end if

        ! distribute file grid data to working grids (different CPU processes different grids)
        call shd_mpi__bcastPrp(Nx_PT, Ny_PT, Nx_prp_, Ny_prp_, Nz_prp_,&
            &NlegSet_prp, Nleg_prp, Npha_prp_,&
            &z_Prp_, maxAsym,&
            &temp_PT, ext_PT, ssa_PT, leg_PT, ipha_PT,&
            &temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, ipha_Prp_)

        if (isMasterProc) deallocate (temp_PT, ext_PT, ssa_PT, ipha_PT, leg_PT)

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readPrp(fileInp_prp, Nx_prp_, Ny_prp_, Nz_prp_,&
            &NlegSet_prp, maxLeg, Nleg_prp, maxNz, maxG_P, maxGL_P, isDeltaM,&
            &type_prp, dx_prp_, dy_prp_, z_Prp_, maxAsym,&
            &temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, Npha_prp_, ipha_Prp_)

        ! Check values of medium properties
        call shd_io__checkPrp(Nx_prp_, Ny_prp_, Nz_prp_, NlegSet_prp, Nleg_prp,&
            &dx_prp_, dy_prp_, z_Prp_,&
            &temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, Npha_prp_, ipha_Prp_)

    end if

    if (isMasterProc) then
        ! Check input parameters
        call shd_io__checkInp(Nmu_do, Nphi_do, isDeltaM, maxAsym,&
            &type_grid, type_src, units_src, solar0, mu_solar,&
            &waveNum, waveLen, temp0_sfc, alb0_sfc,&
            &splitAcc, shAcc)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read fileInp_sfc and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in the surface properties
    if (fileInp_sfc(1:4) .ne. 'NONE') then
        allocate (parms_Sfc(maxNpar_sfc*maxNpts_sfc))
        memUsed = memUsed + maxNpar_sfc*maxNpts_sfc

        if (isMasterProc) then
            call shd_io__readSfc(fileInp_sfc, maxNpts_sfc, maxNpar_sfc,&
                &type_sfc, Nx_sfc, Ny_sfc, dx_sfc, dy_sfc,&
                &Npar_sfc, parms_Sfc, temp0_sfc, alb0_sfc)
        end if
        ! Broadcast the surface array data
        call shd_mpi__bcastSfc(type_sfc, Nx_sfc, Ny_sfc,&
            &Npar_sfc, parms_Sfc, temp0_sfc, alb0_sfc)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read fileInp_ckd and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (isCkd) then
        if (isMasterProc) then
            call shd_io__readCkdSize(fileInp_ckd, waveNum, Ng_ckd, Nz_ckd_)
        end if
        ! Broadcast the k-distribution array sizes
        call shd_mpi__bcastCkdSize(Ng_ckd, Nz_ckd_)
        allocate (z_Ckd_(Nz_ckd_), kabs0_Ckd_(Nz_ckd_), kabs_Ckd(Nz_ckd_, Ng_ckd))
        allocate (solar_Ckd(Ng_ckd), delg_Ckd(Ng_ckd))
        memUsed = memUsed + Nz_ckd_*(2.0+Ng_ckd) + 2.0*Ng_ckd

        if (isMasterProc) then
            maxNg = Ng_ckd; maxNz = Nz_ckd_
            call shd_io__readCkd(maxNg, maxNz, fileInp_ckd, waveNum, solar_Ckd,&
                &Ng_ckd, delg_Ckd, Nz_ckd_, z_Ckd_, kabs_Ckd)
        end if
        call shd_mpi__bcastCkd(solar_Ckd, Ng_ckd, delg_Ckd, Nz_ckd_, z_Ckd_, kabs_Ckd)
    else
        Ng_ckd = 1
        allocate (solar_Ckd(Ng_ckd), delg_Ckd(Ng_ckd))
        memUsed = memUsed + 2.0*Ng_ckd

        solar_Ckd(1) = 1.0
        delg_Ckd(1) = 1.0
        Nz_ckd_ = 0
        isBaseOut = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read fileInp_sen and broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in user-specified sensor properties if given
    !   this is a binary angle input implemented for 'R' and 'V3'
    if (fileInp_sen(1:4) .ne. 'NONE') then
        if (isMasterProc) then
            call shd_io__readSenSize(fileInp_sen, Npts_sen, Npar_sen)
        end if
        call shd_mpi__bcastSenSize(Npts_sen, Npar_sen)
        allocate (parms_Sen(Npts_sen*Npar_sen))
        memUsed = memUsed + Npts_sen*Npar_sen

        if (isMasterProc) then
            call shd_io__readSen(fileInp_sen, Npts_sen, Npar_sen,&
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
                Nx_out = 1
            else
                if (isMpiX_) then
                    Nx_out = max(1, nint(xDomain_prp_/parms_Out(2, iout)))
                else
                    Nx_out = max(1, nint((xDomain_prp_+ 2*parms_Out(4, iout))/parms_Out(2, iout)))
                end if
            end if
            if (parms_Out(3, iout) .le. 0.0) then
                Ny_out = 1
            else
                if (isMpiY_) then
                    Ny_out = max(1, nint(yDomain_prp_/parms_Out(3, iout)))
                else
                    Ny_out = max(1, nint((yDomain_prp_+ 2*parms_Out(5, iout))/parms_Out(3, iout)))
                end if
            end if
            Nrad = Nrad + Nstokes*Nx_out*Ny_out*nint(parms_Out(6, iout))
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
    memPerCpuMB_ = 0.9*memTotSysMB_/real(Nproc)
    if (MAX_TOTAL_MB .gt. memPerCpuMB_) then
        MAX_TOTAL_MB = memPerCpuMB_
        if (verbose_ .ge. 5) write (6, '(A,F0.2," MB")') 'Message [shdom]: MAX_TOTAL_MB reduced to: ', MAX_TOTAL_MB
    else
        if (verbose_ .ge. 5) write (6, '(A,F0.2," MB")') 'Message [shdom]: MAX_TOTAL_MB is set to: ', MAX_TOTAL_MB
    end if

    if ((splitAcc .le. 1.0d-12) .and. (fileInp_sav(1:4) .eq. 'NONE')) then
        ADAPT_GRID_FACTOR = 1.0
    end if

    NbigArr = 3
    if (.not. isAccel) NbigArr = 2

    ! NmaxIG
    !╭──────────────────────────────────────────────────────────────╮!
    NmaxIG = 0.0

    NmaxIG = NmaxIG + 4.5 ! temp, planck, ext, ssa, ipha[0.5, integer(kind=2)]

    if (Npha_prp_ .gt. 0) then
        memUsed = memUsed + Npha_prp_*NlegSet_prp*(maxLeg + 1) ! offset phase table
    else
        NmaxIG = NmaxIG + NlegSet_prp*(maxLeg + 1)
    end if

    if (fileInp_sfc(1:4) .ne. 'NONE') then
        if (type_sfc(2:2) .eq. 'L') then
            NmaxIG = real(3)/real(Nz) * (maxNpar_sfc + 2.0 + Nstokes*2.0)
        else
            NmaxIG = real(3)/real(Nz) * (maxNpar_sfc + 2.0 + Nstokes*(2.0+ Nmu_do*maxNphi0_do/2))
        end if
    end if

    if (isFluxOut)  NmaxIG = NmaxIG + 3.0
    if (isHROut)    NmaxIG = NmaxIG + 1.0
    if (isSHOut)    NmaxIG = NmaxIG + Nsh_out
    if (isSrcOut)   NmaxIG = NmaxIG + Nstokes

    NmaxIG = NmaxIG + 16.5*CELL_TO_POINT_RATIO ! ptr8Grid_Cell[8.0], ptr6Cell_Cell[6.0], ptr2Tree_Cell[2.0], flag_Cell[0.5]

    NmaxIG = NmaxIG + 6.0   ! xyz_Grid[3.0], ptrRad_Sh[1.0], ptrSrc_Sh[1.0], ptrSrcDel_Sh[1.0]
    memUsed = memUsed + 4.0 ! offset ptrRad_Sh[+2], ptrSrc_Sh[+1], ptrSrcDel_Sh[+1]

    if (isAccel) then
        NmaxIG = NmaxIG + Nstokes*(1.0+3.0*NUM_SH_TERM_FACTOR*Nlm_sh) ! radiance, source, delta source
    else
        NmaxIG = NmaxIG + Nstokes*(1.0+2.0*NUM_SH_TERM_FACTOR*Nlm_sh) ! radiance, source
    end if

    NmaxIG = NmaxIG + 3.0 ! fluxTot[2.0], fluxDir[1.0]

    memUsed = memUsed + (Nrad+Nlos) ! offset sumRad_Out[Nrad] and sumCam_Out[Nlos]

    NmaxIG = NmaxIG + (Nstokes*(max(maxNphi0_do, 4)+1.0) + 8.0) !work, work1, work2

    memUsed = memUsed + Nrad + Nlos ! work2 (over-count but memory safe)
    !╰──────────────────────────────────────────────────────────────╯!

    ! memUsed is now more accurate than before
    ! before:
    !   memUsed = Nmu_do*(2 + 2*Nphi_do+ 2*Nlm_sh + 2*33*32)&            ! SHDO allocation
    !       + 4.5*maxG_P + maxGL_P+ Npha_prp_*NlegSet_prp*(maxLeg + 1) ! phase allocation
    ! now: update memUsed after every allocation step
    !
    ! NmaxIG is now more accurate than before
    ! before:
    !   NmaxIG = 28.0 + 16.5*CELL_TO_POINT_RATIO&
    !       + maxNphi0_do*Nstokes + NUM_SH_TERM_FACTOR*Nstokes*Nlm_sh*NbigArr
    ! now: update NmaxIG from every allocation step

    memNeed = ADAPT_GRID_FACTOR*Npts_base*NmaxIG

    if (isHighOrder_) then
        scaleFactor = min(2.0, ((MAX_TOTAL_MB*1024**2)/4 - memUsed)/memNeed)
    else
        scaleFactor = min(1.5, ((MAX_TOTAL_MB*1024**2)/4 - memUsed)/memNeed)
    end if

    ADAPT_GRID_FACTOR = ADAPT_GRID_FACTOR*scaleFactor
    if (ADAPT_GRID_FACTOR .lt. 0.999999999d0) then
        stop 'Error [shdom]: MAX_TOTAL_MB memory limit exceeded with just base grid.'
    end if

    if (scaleFactor .lt. 1.0) then
        if (verbose_ .ge. 5) write (6, '(A,F8.2)') 'Message [shdom]: ADAPT_GRID_FACTOR reduced to:', ADAPT_GRID_FACTOR
    else
        if (verbose_ .ge. 5) write (6, '(A,F8.2)') 'Message [shdom]: ADAPT_GRID_FACTOR increased to:', ADAPT_GRID_FACTOR
    end if

    memNeed = scaleFactor*memNeed
    if (memNeed .gt. huge(maxIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Number of words of memory exceeds max integer size:', memNeed
        stop
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Retrieve max numbers for array allocation after memory check
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! max number of internal grid points
    maxIG = nint(ADAPT_GRID_FACTOR*Npts_base)

    ! max number of grid cells, between 1 and 2 times grid points
    maxIC = nint(CELL_TO_POINT_RATIO*maxIG)

    ! max words for source function and radiance arrays
    maxIV = nint(NUM_SH_TERM_FACTOR*Nlm_sh*maxIG)

    ! max words for temporary discrete ordinate azimuthal array
    maxIDO = maxIG*max(maxNphi0_do, 4)

    if (4.0*(maxIV + maxIG)*Nstokes .gt. huge(maxIV)) then
        write (6, '(A,I12)') 'Error [shdom]: Size of big SH arrays (maxIV) probably exceeds max integer number of bytes:', maxIV
        stop
    end if

    if (4.0*8.0*maxIC .gt. huge(maxIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Size of ptr8Grid_Cell array (8*maxIC) probably exceeds max integer number of bytes:', 8*maxIC
        stop
    end if

    maxNpts_bc = maxIG*3/Nz
    if (type_sfc(2:2) .eq. 'L') then ! Lambertian surface
        maxNrad_bc = 2 * maxNpts_bc
    else
        maxNrad_bc = (2 + Nmu_do*maxNphi0_do/2)*maxNpts_bc
    end if

    ! max number of legendre coefficients
    if (Npha_prp_ .gt. 0) then
        maxIGL = Npha_prp_*NlegSet_prp*(maxLeg + 1)
    else
        maxIGL = maxIG*NlegSet_prp*(maxLeg + 1)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Allocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! temperature and planck radiance (maybe we don't need this for shortwave)
    allocate (temp(maxIG), planck(maxIG))

    ! extinction, single scattering albedo, and phase function indices
    allocate (ext(maxIG), ssa(maxIG), ipha(maxIG))

    ! legendre coefficients
    allocate (leg(maxIGL))

    ! boundary points pointers and radiances
    !   ptr_Bc(:, 1): pointers to grid points that are on top
    !   ptr_Bc(:, 2): pointers to grid points that are on bottom
    !   rad_Bc(:[Nstokes], 1:NptsTop_bc): isotropic top radiances
    !   rad_Bc(:[Nstokes], NptsTop_bc+1:NptsTop_bc+NptsBot_bc): upwelling bottom radiances
    !   rad_Bc(:[Nstokes], NptsTop_bc+NptsBot_bc+1:NptsTop_bc+NptsBot_bc+NptsBot_bc*Nang_do/2):
    !         downwelling bottom radiances at discrete ordinates
    if (fileInp_sfc(1:4) .ne. 'NONE') allocate (parmsGrid_Sfc(maxNpar_sfc*maxNpts_bc))
    allocate (ptr_Bc(2*maxNpts_bc), rad_Bc(Nstokes, maxNrad_bc))

    ! output
    if (isFluxOut) allocate (sumFluxTot_Out(2, maxIG), sumFluxDir_Out(maxIG))
    if (isHROut)   allocate (sumFluxDiv_Out(maxIG))
    if (isSHOut)   allocate (sumSH_Out(maxIG*Nsh_out))
    if (isSrcOut)  allocate (src1_Out(maxIG*Nstokes))

    ! pointer to 8 grid points belonging to a cell
    ! pointer to 6 neighboring cells of a cell
    ! pointer to parent and children cells in tree structure
    allocate (ptr8Grid_Cell(8, maxIC), ptr6Cell_Cell(6, maxIC), ptr2Tree_Cell(2, maxIC))

    ! a collection of binary flags for the cells
    allocate (flag_Cell(maxIC))

    ! x,y,z coordinate values of grid point locations
    allocate (xyz_Grid(3, maxIG))

    ! SH related array
    allocate (ptrRad_Sh(maxIG + 2), ptrSrc_Sh(maxIG + 1), ptrSrcDel_Sh(maxIG + 1))

    ! radiance and source
    allocate (rad_Sh(Nstokes, maxIV + maxIG), src_Sh(Nstokes, maxIV))

    ! delta source
    if (isAccel) allocate (srcDel_Sh(Nstokes, maxIV))

    ! fluxes
    allocate (fluxTot(2, maxIG), fluxDir(maxIG))

    ! output radiances to be summed over g
    allocate (sumRad_Out(Nrad), sumCam_Out(Nlos))

    ! work: work array for path integration; radiance field, Nazimuth*Npts
    allocate (work(Nstokes*maxIDO))

    ! work1: work array for path integration; grid point sweeping order
    allocate (work1(8*maxIG))

    ! work2: work array for path integration; single angle radiance
    allocate (work2(max(Nlos, Nrad, Nstokes*maxIG)))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read in previous run from binary file if desired
    call shd_io__restoreState(fileInp_sav, Nstokes, Nx, Ny, Nz, maxNl_sh, maxNm_sh, Nlm_sh,&
        &isNewRad, isNewGrid, x_Grid, y_Grid, z_Grid,&
        &Npts, Ncell, xyz_Grid, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell,&
        &flag_Cell, fluxTot, ptrSrc_Sh, src_Sh, ptrRad_Sh, rad_Sh)


    ! Check if we can use pre-calculated "state", if not, generate new grids
    !   isNewGrid is from shd_io__restoreState call
    if (isNewGrid) then
        ! Make the internal base grid lines
        call shd_init__makeGrid(flag_bc, type_grid, Nx_prp_, Ny_prp_, Nz_prp_, Nx, Ny, Nz,&
            &xStart_prp_, yStart_prp_, dx_prp_, dy_prp_, z_Prp_,&
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
        solar0_ckd = solar_Ckd(ig_ckd)*solar0*abs(mu_solar)

        ! weighted average solar flux
        solar = solar + solar_Ckd(ig_ckd)*delg_Ckd(ig_ckd)*solar0*abs(mu_solar)

        ! If doing k-distribution then get the gas absorption profile.
        if (isCkd) then
            do iz = 1, Nz_ckd_
                kabs0_Ckd_(iz) = kabs_Ckd(iz, ig_ckd)
            end do
        end if

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (isBaseOut .or. ((ig_ckd.eq. Ng_ckd) .and. isNewGrid)) then
            call shd_init__makeCell(flag_bc, ipFlag,&
                &Nx, Ny, Nz, Nx1, Ny1, Npts, Ncell,&
                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell)
            Ncell_base = Ncell
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in source (ptrSrc_Sh), radiance (ptrRad_Sh), fluxes, direct Flux.
        call shd_solv__solveRte(Nstokes, Nx, Ny, Nz, Nx1, Ny1, Nang_do,&
            &maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nleg_prp, NlegSet_prp, Npha_prp_,&
            &Nphi0_Do, mu_Do, phi_Do, wt_Do,&
            &maxIV, maxIC, maxIG, maxIDO, isNewRad,&
            &flag_bc, ipFlag, isDeltaM, type_src, isHighOrder_,&
            &solar0_ckd, mu_solar, az_solar, skyRad,&
            &type_sfc, temp0_sfc, alb0_sfc,&
            &Nx_sfc, Ny_sfc, dx_sfc, dy_sfc,&
            &Npar_sfc, parms_Sfc, parmsGrid_Sfc,&
            &units_src, waveNum, waveLen,&
            &isAccel, solAcc, maxIter, solCrit, iter,&
            &splitAcc, shAcc, x_Grid, y_Grid, z_Grid,&
            &temp, planck, ext, ssa, leg, ipha,&
            &maxNpts_bc, maxNrad_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
            &Npts, xyz_Grid,&
            &Ncell, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
            &ptrRad_Sh, ptrSrc_Sh, ptrSrcDel_Sh, work, work1, work2,&
            &src_Sh, srcDel_Sh, rad_Sh, fluxTot, fluxDir)
        isNewRad = .false.
        totIter = totIter + iter

        ! Compute the output and put in arrays. First the radiance.
        if (isRadOut) then
            irad = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'R') then
                    if (isMpi_) then
                        if (fileInp_sen(1:4) .ne. 'NONE') then
                            call shd_mpi__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, xDomain_prp_, yDomain_prp_, ipFlag,&
                                &type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, work, work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2, parms_Sen)
                        else
                            call shd_mpi__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, xDomain_prp_, yDomain_prp_, ipFlag,&
                                &type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, work, work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2)
                        end if
                    else
                        if (fileInp_sen(1:4) .ne. 'NONE') then
                            call shd_rad__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, ipFlag, type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, work(1), work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2, parms_Sen)
                        else
                            call shd_rad__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, ipFlag, type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, work(1), work(1 + Npts*Nstokes),&
                                &parms_Out(1, iout), irad, work2)
                        end if
                    end if
                end if
            end do
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Nrad, work2, sumRad_Out)
        end if

        if (isCamOut) then
            ilos = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'V') then
                    if (isMpi_) then
                        if (fileInp_sen(1:4) .ne. 'NONE') then
                            call shd_mpi__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, xDomain_prp_, yDomain_prp_, ipFlag,&
                                &type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, parms_Out(1, iout), ilos, work2, parms_Sen)
                        else
                            call shd_mpi__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, xDomain_prp_, yDomain_prp_, ipFlag,&
                                &type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, parms_Out(1, iout), ilos, work2)
                        end if
                    else
                        if (fileInp_sen(1:4) .ne. 'NONE') then
                            call shd_rad__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, ipFlag, type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, parms_Out(1, iout), ilos, work2, parms_Sen)
                        else
                            call shd_rad__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                                &maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_,&
                                &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                                &flag_bc, ipFlag, type_src, isDeltaM, mu_solar, az_solar,&
                                &type_sfc, Npar_sfc, parmsGrid_Sfc,&
                                &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                                &temp0_sfc, alb0_sfc, skyRad, waveNum, waveLen, units_src,&
                                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                                &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                                &ptrSrc_Sh, src_Sh, parms_Out(1, iout), ilos, work2)
                        end if
                    end if
                end if
            end do
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Nlos, work2, sumCam_Out)
        end if

        if (isSrcOut) then
            ! Compute source function throughout grid for this one angle
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'J') then
                    mu_out  = parms_Out(2, iout)
                    phi_out = parms_Out(3, iout)*acos(-1.0)/180.0
                    call shd_src__computeSrcOneDir(Nstokes, maxNl_sh, maxNm_sh, Nlm_sh,&
                        &NlegSet_prp, Nleg_prp, Npha_prp_,&
                        &Npts, isDeltaM, mu_out, phi_out,&
                        &type_src, mu_solar, az_solar, ssa, leg, ipha,&
                        &fluxDir, ptrSrc_Sh, src_Sh, src1_Out)
                end if
            end do
        end if

        if (isBaseOut) then
            ! If outputting only the base grid then throw out the other points
            call shd_init__makeCell(flag_bc, ipFlag,&
                &Nx, Ny, Nz, Nx1, Ny1,&
                &Npts, Ncell, x_Grid, y_Grid, z_Grid, xyz_Grid,&
                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell)
        else
            ! For the output from the grid point arrays (Flux, Net flux, SH)
            !   we may need to interpolate the running k-distribution sum
            !   to the newly created grid points.
            call shd_util__interpData(oldNpts, Npts, isFluxOut, isHROut, isSHOut,&
                sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out, sumSH_Out,&
                Ncell_base, Ncell, ptr2Tree_Cell, ptr8Grid_Cell, xyz_Grid)
        end if

        if (isFluxOut) then

            ! unscale delta-m fluxes
            if ((type_src.ne. 'T') .and. isDeltaM) then

                allocate(fluxDirTmp(Npts))

                if (isMpi_) then
                    call shd_mpi__makeDirect(1, Npts, flag_bc, ipFlag, .false., maxNl_sh,&
                        &NlegSet_prp, Nleg_prp, solar0_ckd, mu_solar, az_solar,&
                        &xyz_Grid, Nx, x_Grid, Ny, y_Grid, fluxDirTmp)
                else
                    call shd_rad__makeDirect(Npts, flag_bc, ipFlag, .false., maxNl_sh,&
                        &NlegSet_prp, Nleg_prp, solar0_ckd, mu_solar, az_solar,&
                        &xyz_Grid, fluxDirTmp)
                end if

                fluxTot(1, 1:Npts) = fluxTot(1, 1:Npts) + fluxDir(1:Npts) - fluxDirTmp(1:Npts)
                fluxDir(1:Npts) = fluxDirTmp(1:Npts)

                deallocate(fluxDirTmp)

            end if

            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), 2*Npts, fluxTot, sumFluxTot_Out)
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Npts, fluxDir, sumFluxDir_Out)
        end if

        if (isHROut) then
            ! Compute the net flux divergence at every point
            call shd_rad__computeNetFluxDiv(Nstokes, Npts, ptrRad_Sh, type_src, mu_solar,&
                &ext, ssa, planck, fluxDir, rad_Sh, work)
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Npts, work, sumFluxDiv_Out)
        end if

        if (isSHOut) then
            call shd_shdo__computeSH(Nsh_out, Nstokes, Npts, type_src, mu_solar, az_solar,&
                &fluxDir, ptrRad_Sh, maxNl_sh, maxNm_sh, rad_Sh, work)
            call shd_util__sumData(Ng_ckd-ig_ckd, delg_Ckd(ig_ckd), Nsh_out*Npts, work, sumSH_Out)
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
    Nsh = ptrSrc_Sh(Npts + 1)
    irad = 0
    ilos = 0
    if (isMpi_) then

        if (isMasterProc) then
            if (isFluxOut) allocate (fluxTot_all(3, Nz, Ny_T, Nx_T))
            if (isHROut)   allocate (fluxDiv_all(Nz, Ny_T, Nx_T))
            if (isSHOut)   allocate (SH_all(Nsh_out, Nz, Ny_T, Nx_T))
            allocate (xGrid_all(Nx_T), yGrid_all(Ny_T))
        end if

        call shd_mpi__gatherOut(isFluxOut, isHROut, isSHOut,&
            &ipFlag, flag_bc, Npts_base, Nx_T, Ny_T, Nz, Nsh_out,&
            &sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out, sumSH_Out,&
            &xDomain_prp_, yDomain_prp_, xGrid_all, yGrid_all,&
            &fluxTot_all, fluxDiv_all, SH_all,&
            &Ncell, Npts, Nsh, Ncell_tot, Npts_tot, Nsh_tot)

        if (isMasterProc) then

            do iout = 1, Nout
                if (type_Out(iout) .ne. 'V' .and. file_Out(iout)(1:4) .ne. 'NONE') then
                    if (fileInp_sen(1:4) .ne. 'NONE') then
                        call shd_io__outputResultMpi(Nstokes, Nx_T, Ny_T, Nz,&
                            &Npts_tot, Ncell_tot, Nsh_tot,&
                            &maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nang_do, Ng_ckd,&
                            &fileInp_prp, fileInp_sfc, fileInp_ckd, fileInp_sav, fileOut_sav,&
                            &flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                            &solar, mu_solar, az_solar, skyRad,&
                            &type_sfc, temp0_sfc, alb0_sfc, waveNum, waveLen, units_src,&
                            &splitAcc, shAcc, solAcc, maxIter, totIter,&
                            &xDomain_prp_, yDomain_prp_, xGrid_all, yGrid_all, z_Grid,&
                            &fluxTot_all, fluxDiv_all,&
                            &Nsh_out, SH_all, irad, sumRad_Out,&
                            &type_Out(iout), parms_Out(1, iout), file_Out(iout), parms_Sen)
                    else
                        call shd_io__outputResultMpi(Nstokes, Nx_T, Ny_T, Nz,&
                            &Npts_tot, Ncell_tot, Nsh_tot,&
                            &maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nang_do, Ng_ckd,&
                            &fileInp_prp, fileInp_sfc, fileInp_ckd, fileInp_sav, fileOut_sav,&
                            &flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                            &solar, mu_solar, az_solar, skyRad,&
                            &type_sfc, temp0_sfc, alb0_sfc, waveNum, waveLen, units_src,&
                            &splitAcc, shAcc, solAcc, maxIter, totIter,&
                            &xDomain_prp_, yDomain_prp_, xGrid_all, yGrid_all, z_Grid,&
                            &fluxTot_all, fluxDiv_all,&
                            &Nsh_out, SH_all, irad, sumRad_Out,&
                            &type_Out(iout), parms_Out(1, iout), file_Out(iout))
                    end if
                end if
            end do

            if (isFluxOut) deallocate (fluxTot_all)
            if (isHROut)   deallocate (fluxDiv_all)
            if (isSHOut)   deallocate (SH_all)
            deallocate (xGrid_all, yGrid_all)
        end if

    else

        do iout = 1, Nout
            if ((type_Out(iout) .ne. 'V') .and. (file_Out(iout)(1:4) .ne. 'NONE')) then
                if (fileInp_sen(1:4) .ne. 'NONE') then
                    call shd_io__outputResult(Nstokes, Nx, Ny, Nx1, Ny1, Nz, Npts_base, Npts, Ncell,&
                        Nsh, maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_, Nmu_do, Nphi_do, Nang_do, Ng_ckd,&
                        fileInp_prp, fileInp_sfc, fileInp_ckd, fileInp_sav, fileOut_sav,&
                        flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                        solar, mu_solar, az_solar, skyRad,&
                        type_sfc, temp0_sfc, alb0_sfc, waveNum, waveLen, units_src,&
                        splitAcc, shAcc, solAcc, maxIter, totIter,&
                        x_Grid, y_Grid, z_Grid, xyz_Grid,&
                        ptr2Tree_Cell, ptr8Grid_Cell, flag_Cell,&
                        ext, ssa, leg, ipha, temp,&
                        sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out,&
                        irad, sumRad_Out, Nsh_out, sumSH_Out, src1_Out,&
                        type_Out(iout), parms_Out(1, iout), file_Out(iout), parms_Sen)
                else
                    call shd_io__outputResult(Nstokes, Nx, Ny, Nx1, Ny1, Nz, Npts_base, Npts, Ncell,&
                        Nsh, maxNl_sh, maxNm_sh, Nlm_sh, NlegSet_prp, Nleg_prp, Npha_prp_, Nmu_do, Nphi_do, Nang_do, Ng_ckd,&
                        fileInp_prp, fileInp_sfc, fileInp_ckd, fileInp_sav, fileOut_sav,&
                        flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                        solar, mu_solar, az_solar, skyRad,&
                        type_sfc, temp0_sfc, alb0_sfc, waveNum, waveLen, units_src,&
                        splitAcc, shAcc, solAcc, maxIter, totIter,&
                        x_Grid, y_Grid, z_Grid, xyz_Grid,&
                        ptr2Tree_Cell, ptr8Grid_Cell, flag_Cell,&
                        ext, ssa, leg, ipha, temp,&
                        sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out,&
                        irad, sumRad_Out, Nsh_out, sumSH_Out, src1_Out,&
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
                if (fileInp_sen(1:4) .ne. 'NONE') then
                    call shd_io__outputCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                        &Nsh, maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nang_do, Ng_ckd,&
                        &fileInp_prp, fileInp_sfc, fileInp_ckd, fileInp_sav, fileOut_sav,&
                        &flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                        &solar, mu_solar, az_solar, skyRad,&
                        &type_sfc, temp0_sfc, alb0_sfc, waveNum, waveLen, units_src,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &ilos, sumCam_Out, parms_Out(1, iout), file_Out(iout), parms_Sen)
                else
                    call shd_io__outputCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                        &Nsh, maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nang_do, Ng_ckd,&
                        &fileInp_prp, fileInp_sfc, fileInp_ckd, fileInp_sav, fileOut_sav,&
                        &flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                        &solar, mu_solar, az_solar, skyRad,&
                        &type_sfc, temp0_sfc, alb0_sfc, waveNum, waveLen, units_src,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &ilos, sumCam_Out, parms_Out(1, iout), file_Out(iout))
                end if
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! save state
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (maxIter .gt. 0) then
        call shd_io__saveState(fileOut_sav, Nstokes, Nx, Ny, Nz, maxNl_sh, maxNm_sh, Nlm_sh,&
            &work, x_Grid, y_Grid, z_Grid,&
            &Npts, Ncell, xyz_Grid, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell,&
            &flag_Cell, fluxTot, ptrSrc_Sh, src_Sh, ptrRad_Sh, rad_Sh)
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
    !  MAXMB_OUT = 4*(Nmu_do*(2 + 2*Nphi_do + 2*Nlm_sh + 2*33*32)&
    !      + 4.5*maxG_P + maxGL_P + Npha_prp_*NlegSet_prp*(maxLeg + 1)&
    !      + 16.5*Ncell + Npts*(28 + maxNphi0_do*Nstokes)&
    !      + Nstokes*Nsh*NbigArr)/1024**2
    ! now:
    !   use NmaxIG
    MAXMB_OUT = 4*(memUsed + 16.5*Ncell&
        + Npts*(NmaxIG-16.5*CELL_TO_POINT_RATIO))/1024**2

    ADAPTGRIDFAC_OUT = float(Npts)/Npts_base
    SHTERM_FAC_OUT = float(Nsh)/(Nlm_sh*Npts)
    CELL_POINT_OUT = float(Ncell)/Npts
    if (verbose_ .ge. 5) write (6, '(A)')'Message [shdom]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    if (verbose_ .ge. 5) write (6, '(26X,F9.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Deallocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    deallocate (z_Prp_, temp_Prp_, ext_Prp_, ssa_Prp_, ipha_Prp_, leg_Prp_, extDir_Prp_)
    if (isFluxOut) deallocate (sumFluxTot_Out, sumFluxDir_Out)
    if (isHROut)   deallocate (sumFluxDiv_Out)
    if (isSHOut)   deallocate (sumSH_Out)
    if (isSrcOut)  deallocate (src1_Out)
    if (isCkd) deallocate (z_Ckd_, kabs0_Ckd_, kabs_Ckd)
    if (allocated(parms_Sfc)) deallocate (parms_Sfc, parmsGrid_Sfc)
    if (allocated(parms_Sen)) deallocate (parms_Sen)
    deallocate (delg_Ckd, solar_Ckd)
    deallocate (ptr_Bc, rad_Bc)

    deallocate (mu_Do, phi_Do, Nphi0_Do, wt_Do)
    deallocate (x_Grid, y_Grid, z_Grid)
    deallocate (temp, planck)
    deallocate (ext, ssa, ipha)
    deallocate (leg)
    deallocate (ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell)
    deallocate (flag_Cell)
    deallocate (xyz_Grid)
    deallocate (ptrRad_Sh, ptrSrc_Sh, ptrSrcDel_Sh)
    deallocate (rad_Sh, src_Sh)
    if (isAccel) deallocate (srcDel_Sh)
    deallocate (fluxTot, fluxDir)
    deallocate (sumRad_Out, sumCam_Out)
    deallocate (work, work1, work2)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    call shd_mpi__end(Npts, xyz_Grid, Nx_prp_, Ny_prp_, xStart_prp_, yStart_prp_, xDomain_prp_, yDomain_prp_,&
        &Nx_PT, Ny_PT, fileInp_prp, runName)

end program shdom
