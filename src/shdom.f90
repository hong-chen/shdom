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
                             isSolar_, isThermal_, &
                             shd_type_prp, shd_type_sfc, &
                             shd_type_ckd, shd_type_sen

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
    type(shd_type_prp) :: prp_obj, prpAll_obj

    real(kind=R__)   :: dxTmp_prp, dyTmp_prp

    real(kind=R__), allocatable, target :: dataRaw_z_Prp(:)
    real(kind=R__), allocatable, target :: dataRaw_temp_Prp(:), dataRaw_ext_Prp(:), dataRaw_ssa_Prp(:)
    real(kind=R__), allocatable, target :: dataRaw_extDir_Prp(:)
    real(kind=R__), allocatable, target :: dataRaw_leg_Prp(:)
    integer(kind=2), allocatable, target :: dataRaw_ipha_Prp(:)

    real(kind=R__), allocatable, target  :: dataRaw_temp_PrpAll(:), dataRaw_ext_PrpAll(:), dataRaw_ssa_PrpAll(:)
    real(kind=R__), allocatable, target  :: dataRaw_leg_PrpAll(:)
    integer(kind=2), allocatable, target :: dataRaw_ipha_PrpAll(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! sfc_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_sfc) :: sfc_obj
    real(kind=R__), allocatable, target :: dataRaw_brdf_Sfc(:)
    real(kind=R__), allocatable, target :: dataGrid_brdf_Sfc(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! ckd_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_ckd) :: ckd_obj
    integer(kind=I__) :: ig_ckd
    real(kind=R__)    :: solar0_ckd
    real(kind=R__), allocatable, target :: dataRaw_solar_Ckd(:), &
                                           dataRaw_delg_Ckd(:), &
                                           dataRaw_z_Ckd(:), &
                                           dataRaw_kabs_Ckd(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! sen_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_sen) :: sen_obj
    real(kind=R__), allocatable, target :: dataRaw_view_Sen(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! BC (borrowed from Boundary Condition)
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

    integer(kind=I__) :: maxIG, maxIC, maxIV, maxIDO, maxNg
    integer(kind=I__) :: maxIGL, NbigArr

    real(kind=R__) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real(kind=R__) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real(kind=R__) :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer(kind=I__) :: Nstokes, Nx, Ny, Nz, Nx1, Ny1
    integer(kind=I__) :: Nx_T, Ny_T

    integer(kind=I__) :: Npts_base, Npts, Ncell, Nsh, Ncell_base, oldNpts
    integer(kind=I__) :: maxIter, iter, totIter, iz, iout, Nout, flag_bc, ipFlag
    integer(kind=I__) :: Nrad, irad, Nx_out, Ny_out, Nlos, ilos, Nscan, Npix, Nsh_out
    integer(kind=I__) :: Ncell_tot, Npts_tot, Nsh_tot

    character(len=1) :: type_src, units_src
    real(kind=R__) :: solar0, solar, mu_solar, az_solar

    real(kind=R__) :: skyRad, waveNum(2), waveLen
    real(kind=R__) :: solAcc, solCrit, splitAcc, shAcc

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
            call shd_io__nmlInp(runName, prpAll_obj%fname, sfc_obj%fname, ckd_obj%fname,&
                fileInp_sav, fileOut_sav, sen_obj%fname,&
                Nstokes, Nx, Ny, Nz, Nmu_do, Nphi_do, flag_bc,&
                ipFlag, isCkd, isDeltaM, type_grid,&
                type_src, solar0, mu_solar, az_solar, skyRad,&
                sfc_obj%avgTemp, sfc_obj%avgAlb, units_src, waveNum, waveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        else
            call shd_io__userInp(runName, prpAll_obj%fname, sfc_obj%fname, ckd_obj%fname,&
                fileInp_sav, fileOut_sav, sen_obj%fname,&
                Nstokes, Nx, Ny, Nz, Nmu_do, Nphi_do, flag_bc,&
                ipFlag, isCkd, isDeltaM, type_grid,&
                type_src, solar0, mu_solar, az_solar, skyRad,&
                sfc_obj%avgTemp, sfc_obj%avgAlb, units_src, waveNum, waveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        end if
    end if

    call shd_mpi__bcastInp(runName,&
        &fileInp_sav, fileOut_sav,&
        &Nstokes, Nx, Ny, Nz, Nmu_do, Nphi_do, flag_bc,&
        &ipFlag, isCkd, isDeltaM, type_grid,&
        &type_src, solar0, mu_solar, az_solar, skyRad,&
        &units_src, waveNum, waveLen,&
        &isAccel, solAcc, maxIter, splitAcc, shAcc,&
        &maxNset_out, maxNpar_out, Nout, type_Out, parms_Out,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, prpAll_obj, sfc_obj, ckd_obj, sen_obj)
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

    ! prpAll_obj%NlegSet is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (Nstokes .eq. 1) then
        prpAll_obj%NlegSet = 1
        isPol_ = .false.
    else
        prpAll_obj%NlegSet = 6
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
    !   maxNphi0_do: The maximum number of azimuth angles actually used
    Nmu_do = max(2, 2*int((Nmu_do + 1)/2)) ! force Nmu_do an even number
    Nphi_do = max(1, Nphi_do)
    maxNphi0_do = Nphi_do

    ! maxNl_sh: The maximum L of the spherical harmonic truncation
    ! maxNm_sh: The maximum M of the spherical harmonic truncation
    maxNl_sh = Nmu_do - 1
    maxNm_sh = max(0, int(Nphi_do/2) - 1)

    ! Nlm_sh: The number of spherical harmonic terms
    Nlm_sh = (2*maxNm_sh + 1)*(maxNl_sh + 1) - maxNm_sh*(maxNm_sh + 1)

    ! prpAll_obj%Nleg: The number of Legendre series phase function terms used in the solution procedure
    prpAll_obj%Nleg = maxNl_sh
    if (isDeltaM) prpAll_obj%Nleg = maxNl_sh + 1

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


    ! Read prpAll_obj%fname (property file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        ! Find out the size of the input property file grid
        call shd_io__readPrpSize(prpAll_obj)
    end if

    call shd_mpi__bcastPrpSize(prpAll_obj)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sfc_obj%fname (surface file) to get general info, broadcast sfc_obj if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (sfc_obj%fname(1:4) .eq. 'NONE') then
        sfc_obj%fname = 'NONE'
        sfc_obj%type  = 'FL'
        sfc_obj%Nx = 1
        sfc_obj%Ny = 1
        sfc_obj%dx = prpAll_obj%dx
        sfc_obj%dy = prpAll_obj%dy
        sfc_obj%maxNpts = 1
        sfc_obj%maxNpar = 2
        sfc_obj%isUniform    = .true.
        sfc_obj%isLambertian = .true.
    else
        if (isMasterProc) then
            ! find out the size of the surface file grid
            call shd_io__readSfcSize(sfc_obj)
        end if
        ! Broadcast the surface object
        call shd_mpi__bcastSfcSize(sfc_obj)
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

    if (prpAll_obj%Nx .eq. 1) then
        prpAll_obj%dx = sfc_obj%Nx * sfc_obj%dx ! overwrite prpAll_obj%dx if only 1 grid in X direction
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: prpAll_obj%Nx=1, reset prpAll_obj%dx to match surface domain.'
    end if

    if (prpAll_obj%Ny .eq. 1) then
        prpAll_obj%dy = sfc_obj%Ny * sfc_obj%dy ! overwrite prpAll_obj%dy if only 1 grid in Y direction
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: prpAll_obj%Ny=1, reset prpAll_obj%dy to match surface domain.'
    end if

    prpAll_obj%xDomain = prpAll_obj%Nx * prpAll_obj%dx
    prpAll_obj%yDomain = prpAll_obj%Ny * prpAll_obj%dy

    ! Choose the option to compute higher order radiance SH terms
    if  ((prpAll_obj%Nx .eq. 1) .and. (prpAll_obj%Ny .eq. 1) .and.&
        &(sfc_obj%Nx .le. 32) .and. (sfc_obj%Ny .le. 32)) then
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
    prp_obj = prpAll_obj
    call shd_mpi__map(flag_bc, Nx, Ny, Nz, runName, prp_obj)

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

    ! Read prp_obj%fname and distribute property object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Determine the size of the property grid arrays so can allocate
    prp_obj%maxNpts = prp_obj%Nx * prp_obj%Ny * prp_obj%Nz
    if (prp_obj%Npha .eq. 0) prp_obj%maxNptsLeg = prp_obj%maxNpts * (prp_obj%maxNleg + 1) * prp_obj%NlegSet

    allocate (dataRaw_z_Prp(prp_obj%Nz))
    allocate (dataRaw_temp_Prp(prp_obj%maxNpts), &
              dataRaw_ext_Prp( prp_obj%maxNpts), &
              dataRaw_ssa_Prp( prp_obj%maxNpts))
    allocate (dataRaw_extDir_Prp(prp_obj%maxNpts), &
              dataRaw_ipha_Prp(prp_obj%maxNpts), &
              dataRaw_leg_Prp( prp_obj%maxNptsLeg))
    memUsed = memUsed + prp_obj%Nz + 6.0*prp_obj%maxNpts + prp_obj%maxNptsLeg

    prp_obj%dataRaw_z      => dataRaw_z_Prp
    prp_obj%dataRaw_temp   => dataRaw_temp_Prp
    prp_obj%dataRaw_ext    => dataRaw_ext_Prp
    prp_obj%dataRaw_extDir => dataRaw_extDir_Prp
    prp_obj%dataRaw_ssa    => dataRaw_ssa_Prp
    prp_obj%dataRaw_ipha   => dataRaw_ipha_Prp
    prp_obj%dataRaw_leg    => dataRaw_leg_Prp
    prp_obj%isAllocated = .true.

    if (isMpi_) then

        ! Read the properties of the medium to the master processor and
        !   distribute to all the processors
        if (isMasterProc) then

            allocate (dataRaw_temp_PrpAll(prpAll_obj%maxNpts), &
                      dataRaw_ext_PrpAll( prpAll_obj%maxNpts), &
                      dataRaw_ssa_PrpAll( prpAll_obj%maxNpts))
            allocate (dataRaw_ipha_PrpAll(prpAll_obj%maxNpts))
                      !dataRaw_leg_PrpAll( prpAll_obj%maxNptsLeg))

            prpAll_obj%dataRaw_z    => dataRaw_z_Prp
            prpAll_obj%dataRaw_leg  => dataRaw_leg_Prp
            prpAll_obj%dataRaw_temp => dataRaw_temp_PrpAll
            prpAll_obj%dataRaw_ext  => dataRaw_ext_PrpAll
            prpAll_obj%dataRaw_ssa  => dataRaw_ssa_PrpAll
            prpAll_obj%dataRaw_ipha => dataRaw_ipha_PrpAll
            prpAll_obj%isAllocated = .true.

            call shd_io__readPrp(isDeltaM, prpAll_obj)

            if (prpAll_obj%Npha .eq. 0) then
                stop 'Error [shdom]: Standard property file (with phase function Legendre series at every point)&
                    & is not allowed with multiple processors.'
            end if

            ! Check values of medium properties
            call shd_io__checkPrp(prpAll_obj)

        end if

        ! distribute file grid data to working grids (different CPU processes different grids)
        call shd_mpi__bcastPrp(prpAll_obj, prp_obj)

        if (isMasterProc) then
            deallocate(dataRaw_temp_PrpAll, &
                       dataRaw_ext_PrpAll, &
                       dataRaw_ssa_PrpAll, &
                       dataRaw_ipha_PrpAll)
                       !dataRaw_leg_PrpAll)
            if (prpAll_obj%isAllocated) then
                nullify(prpAll_obj%dataRaw_z, &
                        prpAll_obj%dataRaw_temp, &
                        prpAll_obj%dataRaw_ext, &
                        prpAll_obj%dataRaw_ssa, &
                        prpAll_obj%dataRaw_ipha, &
                        prpAll_obj%dataRaw_leg)
                prpAll_obj%isAllocated = .false.
            end if
        end if

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readPrp(isDeltaM, prp_obj)

        ! Check values of medium properties
        call shd_io__checkPrp(prp_obj)

    end if

    if (isMasterProc) then
        ! Check input parameters
        call shd_io__checkInp(Nmu_do, Nphi_do, isDeltaM, prp_obj%maxAsym,&
            &type_grid, type_src, units_src, solar0, mu_solar,&
            &waveNum, waveLen, sfc_obj%avgTemp, sfc_obj%avgAlb,&
            &splitAcc, shAcc)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sfc_obj%fname and broadcast the surface object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in the surface properties
    if (sfc_obj%fname(1:4) .ne. 'NONE') then

        allocate (dataRaw_brdf_Sfc(sfc_obj%maxNpar * sfc_obj%maxNpts))
        sfc_obj%dataRaw_brdf => dataRaw_brdf_Sfc
        sfc_obj%isAllocated = .true.

        memUsed = memUsed + sfc_obj%maxNpar * sfc_obj%maxNpts

        if (isMasterProc) then
            call shd_io__readSfc(sfc_obj)
        end if

        ! Broadcast the surface object
        call shd_mpi__bcastSfc(sfc_obj)

    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read ckd_obj%fname and broadcast correlated-k object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (isCkd) then

        if (isMasterProc) then
            call shd_io__readCkdSize(waveNum, ckd_obj)
        end if

        ! Broadcast the k-distribution object
        call shd_mpi__bcastCkdSize(ckd_obj)

        allocate (dataRaw_solar_Ckd(ckd_obj%Ng), dataRaw_delg_Ckd(ckd_obj%Ng), dataRaw_z_Ckd(ckd_obj%Nz))
        allocate (dataRaw_kabs_Ckd(ckd_obj%Nz, ckd_obj%Ng))

        memUsed = memUsed + ckd_obj%Nz*(2.0+ckd_obj%Ng) + 2.0*ckd_obj%Ng

        ckd_obj%dataRaw_solar => dataRaw_solar_Ckd
        ckd_obj%dataRaw_delg  => dataRaw_delg_Ckd
        ckd_obj%dataRaw_z     => dataRaw_z_Ckd
        ckd_obj%dataRaw_kabs  => dataRaw_kabs_Ckd
        ckd_obj%isAllocated = .true.

        if (isMasterProc) then
            ckd_obj%maxNg = ckd_obj%Ng; ckd_obj%maxNz = ckd_obj%Nz
            call shd_io__readCkd(waveNum, ckd_obj)
        end if

        call shd_mpi__bcastCkd(ckd_obj)

    else
        ckd_obj%Ng = 1

        allocate (dataRaw_solar_Ckd(ckd_obj%Ng), dataRaw_delg_Ckd(ckd_obj%Ng))
        memUsed = memUsed + 2.0*ckd_obj%Ng

        dataRaw_solar_Ckd(1) = 1.0
        dataRaw_delg_Ckd(1) = 1.0

        ckd_obj%Nz = 0
        isBaseOut = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sen_obj%fname and broadcast sensor object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in user-specified sensor properties if given
    !   this is a binary angle input implemented for 'R' and 'V3'
    if (sen_obj%fname(1:4) .ne. 'NONE') then

        if (isMasterProc) then
            call shd_io__readSenSize(sen_obj)
        end if

        call shd_mpi__bcastSenSize(sen_obj)

        allocate (dataRaw_view_Sen(sen_obj%Npar*sen_obj%Npts))

        memUsed = memUsed + sen_obj%Npar*sen_obj%Npts

        sen_obj%dataRaw_view => dataRaw_view_Sen
        sen_obj%isAllocated = .true.

        if (isMasterProc) then
            call shd_io__readSen(sen_obj)
        end if

        call shd_mpi__bcastSen(sen_obj)
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
    isHROut   = .false.
    isSHOut   = .false.
    do iout = 1, Nout
        if (type_Out(iout) .eq. 'R') then
            isRadOut = .true.
            if (parms_Out(2, iout) .le. 0.0) then
                Nx_out = 1
            else
                if (isMpiX_) then
                    Nx_out = max(1, nint(prp_obj%xDomain/parms_Out(2, iout)))
                else
                    Nx_out = max(1, nint((prp_obj%xDomain + 2*parms_Out(4, iout))/parms_Out(2, iout)))
                end if
            end if
            if (parms_Out(3, iout) .le. 0.0) then
                Ny_out = 1
            else
                if (isMpiY_) then
                    Ny_out = max(1, nint(prp_obj%yDomain/parms_Out(3, iout)))
                else
                    Ny_out = max(1, nint((prp_obj%yDomain+ 2*parms_Out(5, iout))/parms_Out(3, iout)))
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

    if (prp_obj%Npha .gt. 0) then
        memUsed = memUsed + prp_obj%Npha * prp_obj%NlegSet * (prp_obj%maxNleg + 1) ! offset phase table
    else
        NmaxIG = NmaxIG + prp_obj%NlegSet*(prp_obj%maxNleg + 1)
    end if

    if (sfc_obj%fname(1:4) .ne. 'NONE') then
        if (sfc_obj%type(2:2) .eq. 'L') then
            NmaxIG = real(3)/real(Nz) * (sfc_obj%maxNpar + 2.0 + Nstokes*2.0)
        else
            NmaxIG = real(3)/real(Nz) * (sfc_obj%maxNpar + 2.0 + Nstokes*(2.0+ Nmu_do*maxNphi0_do/2))
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
    !       + 4.5*prp_obj%maxNpts + prp_obj%maxNptsLeg + prp_obj%Npha*prp_obj%NlegSet*(prp_obj%maxNleg + 1) ! phase allocation
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
    if (sfc_obj%isLambertian) then
        maxNrad_bc = 2 * maxNpts_bc
    else
        maxNrad_bc = (2 + Nmu_do*maxNphi0_do/2)*maxNpts_bc
    end if

    ! max number of legendre coefficients
    if (prp_obj%Npha .gt. 0) then
        maxIGL = prp_obj%Npha * prp_obj%NlegSet * (prp_obj%maxNleg + 1)
    else
        maxIGL = maxIG*prp_obj%NlegSet*(prp_obj%maxNleg + 1)
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
    if (sfc_obj%fname(1:4) .ne. 'NONE') then
        allocate (dataGrid_brdf_Sfc(sfc_obj%maxNpar*maxNpts_bc))
        sfc_obj%dataGrid_brdf => dataGrid_brdf_Sfc
        sfc_obj%isAllocated = .true.
    end if

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
        call shd_init__makeGrid(flag_bc, type_grid, prp_obj%Nx, prp_obj%Ny, prp_obj%Nz, Nx, Ny, Nz,&
            &prp_obj%xStart, prp_obj%yStart, prp_obj%dx, prp_obj%dy, prp_obj%dataRaw_z,&
            &x_Grid, y_Grid, z_Grid)
    end if


    message = 'Message [shdom]: Total time used for preparation: '
    call shd_util__tickTock(message, countStart)
    if (isMasterProc) then
        if (verbose_ .ge. 5) write (6, '(A)') trim(adjustl(message))
    end if

    !print *, 'haha1', prp_obj%Nx, prp_obj%Ny, prp_obj%Nz, prp_obj%Npha, prp_obj%Nleg, prp_obj%NlegSet
    !print *, 'haha2', prp_obj%maxNz, prp_obj%maxNpts, prp_obj%maxNleg, prp_obj%maxNptsLeg, prp_obj%maxAsym
    !print *, 'haha3', prp_obj%dx, prp_obj%dy, prp_obj%xStart, prp_obj%yStart, prp_obj%xDomain, prp_obj%yDomain
    !stop

    ! solve radiative transfer equation
    !   1. Use Picard iteration to solve radiance and source
    !   2. Compute actual radiance (and/or other radiative quantities specified by user)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (verbose_ .ge. 5) write (6, '(A,"<",A,">")') 'Message [shdom]: Starting solution iterations for SHDOM run: ', trim(adjustl(runName))

    totIter = 0
    oldNpts = 0
    solar = 0.0d0

    ! Loop over the k-distribution g's from high to low absorption
    CKD_LOOP: do ig_ckd = ckd_obj%Ng, 1, -1

        if (verbose_ .ge. 5) write (6, '(A,I2,A,I2,A)') 'Message [shdom]: Starting ', ckd_obj%Ng-ig_ckd+1, '/', ckd_obj%Ng, ' (ig/Ng) ...'

        ! solar flux for given g of CKD
        solar0_ckd = dataRaw_solar_Ckd(ig_ckd)*solar0*abs(mu_solar)

        ! weighted average solar flux
        solar = solar + dataRaw_solar_Ckd(ig_ckd)*dataRaw_delg_Ckd(ig_ckd)*solar0*abs(mu_solar)

        ! If doing k-distribution then get the gas absorption profile.
        ckd_obj%ig = ig_ckd

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (isBaseOut .or. ((ig_ckd.eq. ckd_obj%Ng) .and. isNewGrid)) then
            call shd_init__makeCell(flag_bc, ipFlag,&
                &Nx, Ny, Nz, Nx1, Ny1, Npts, Ncell,&
                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell)
            Ncell_base = Ncell
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in source (ptrSrc_Sh), radiance (ptrRad_Sh), fluxes, direct Flux.
        call shd_solv__solveRte(Nstokes, Nx, Ny, Nz, Nx1, Ny1, Nang_do,&
            &maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, prp_obj%Nleg, prp_obj%NlegSet, prp_obj%Npha,&
            &Nphi0_Do, mu_Do, phi_Do, wt_Do,&
            &maxIV, maxIC, maxIG, maxIDO, isNewRad,&
            &flag_bc, ipFlag, isDeltaM, type_src, isHighOrder_,&
            &solar0_ckd, mu_solar, az_solar, skyRad,&
            &units_src, waveNum, waveLen,&
            &isAccel, solAcc, maxIter, solCrit, iter,&
            &splitAcc, shAcc, x_Grid, y_Grid, z_Grid,&
            &temp, planck, ext, ssa, leg, ipha,&
            &maxNpts_bc, maxNrad_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
            &Npts, xyz_Grid,&
            &Ncell, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
            &ptrRad_Sh, ptrSrc_Sh, ptrSrcDel_Sh, work, work1, work2,&
            &src_Sh, srcDel_Sh, rad_Sh, fluxTot, fluxDir, sfc_obj, ckd_obj)
        isNewRad = .false.
        totIter = totIter + iter

        ! Compute the output and put in arrays. First the radiance.
        if (isRadOut) then
            irad = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'R') then
                    if (isMpi_) then
                        call shd_mpi__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                            &maxNl_sh, maxNm_sh, Nlm_sh, prp_obj%NlegSet, prp_obj%Nleg, prp_obj%Npha,&
                            &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                            &flag_bc, prp_obj%xDomain, prp_obj%yDomain, ipFlag,&
                            &type_src, isDeltaM, mu_solar, az_solar,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, waveNum, waveLen, units_src,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &ptrSrc_Sh, src_Sh, work, work(1 + Npts*Nstokes),&
                            &parms_Out(1, iout), irad, work2, sfc_obj, sen_obj)
                    else
                        call shd_rad__computeRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                            &maxNl_sh, maxNm_sh, Nlm_sh, prp_obj%NlegSet, prp_obj%Nleg, prp_obj%Npha,&
                            &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                            &flag_bc, ipFlag, type_src, isDeltaM, mu_solar, az_solar,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, waveNum, waveLen, units_src,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &ptrSrc_Sh, src_Sh, work(1), work(1 + Npts*Nstokes),&
                            &parms_Out(1, iout), irad, work2, sfc_obj, sen_obj)
                    end if
                end if
            end do
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, dataRaw_delg_Ckd(ig_ckd), Nrad, work2, sumRad_Out)
        end if

        if (isCamOut) then
            ilos = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'V') then
                    if (isMpi_) then
                        call shd_mpi__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                            &maxNl_sh, maxNm_sh, Nlm_sh, prp_obj%NlegSet, prp_obj%Nleg, prp_obj%Npha,&
                            &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                            &flag_bc, prp_obj%xDomain, prp_obj%yDomain, ipFlag,&
                            &type_src, isDeltaM, mu_solar, az_solar,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, waveNum, waveLen, units_src,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &ptrSrc_Sh, src_Sh, parms_Out(1, iout), ilos, work2, sfc_obj, sen_obj)
                    else
                        call shd_rad__computeCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                            &maxNl_sh, maxNm_sh, Nlm_sh, prp_obj%NlegSet, prp_obj%Nleg, prp_obj%Npha,&
                            &Nmu_do, maxNphi0_do, Nphi0_Do, mu_Do, phi_Do, wt_Do,&
                            &flag_bc, ipFlag, type_src, isDeltaM, mu_solar, az_solar,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, waveNum, waveLen, units_src,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &ptrSrc_Sh, src_Sh, parms_Out(1, iout), ilos, work2, sfc_obj, sen_obj)
                    end if
                end if
            end do
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, dataRaw_delg_Ckd(ig_ckd), Nlos, work2, sumCam_Out)
        end if

        if (isSrcOut) then
            ! Compute source function throughout grid for this one angle
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'J') then
                    mu_out  = parms_Out(2, iout)
                    phi_out = parms_Out(3, iout)*acos(-1.0)/180.0
                    call shd_src__computeSrcOneDir(Nstokes, maxNl_sh, maxNm_sh, Nlm_sh,&
                        &prp_obj%NlegSet, prp_obj%Nleg, prp_obj%Npha,&
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
                        &prp_obj%NlegSet, prp_obj%Nleg, solar0_ckd, mu_solar, az_solar,&
                        &xyz_Grid, Nx, x_Grid, Ny, y_Grid, fluxDirTmp, ckd_obj)
                else
                    call shd_rad__makeDirect(Npts, flag_bc, ipFlag, .false., maxNl_sh,&
                        &prp_obj%NlegSet, prp_obj%Nleg, solar0_ckd, mu_solar, az_solar,&
                        &xyz_Grid, fluxDirTmp, ckd_obj)
                end if

                fluxTot(1, 1:Npts) = fluxTot(1, 1:Npts) + fluxDir(1:Npts) - fluxDirTmp(1:Npts)
                fluxDir(1:Npts) = fluxDirTmp(1:Npts)

                deallocate(fluxDirTmp)

            end if

            call shd_util__sumData(ckd_obj%Ng-ig_ckd, dataRaw_delg_Ckd(ig_ckd), 2*Npts, fluxTot, sumFluxTot_Out)
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, dataRaw_delg_Ckd(ig_ckd), Npts, fluxDir, sumFluxDir_Out)
        end if

        if (isHROut) then
            ! Compute the net flux divergence at every point
            call shd_rad__computeNetFluxDiv(Nstokes, Npts, ptrRad_Sh, type_src, mu_solar,&
                &ext, ssa, planck, fluxDir, rad_Sh, work)
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, dataRaw_delg_Ckd(ig_ckd), Npts, work, sumFluxDiv_Out)
        end if

        if (isSHOut) then
            call shd_shdo__computeSH(Nsh_out, Nstokes, Npts, type_src, mu_solar, az_solar,&
                &fluxDir, ptrRad_Sh, maxNl_sh, maxNm_sh, rad_Sh, work)
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, dataRaw_delg_Ckd(ig_ckd), Nsh_out*Npts, work, sumSH_Out)
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
            &prp_obj%xDomain, prp_obj%yDomain, xGrid_all, yGrid_all,&
            &fluxTot_all, fluxDiv_all, SH_all,&
            &Ncell, Npts, Nsh, Ncell_tot, Npts_tot, Nsh_tot)

        if (isMasterProc) then

            do iout = 1, Nout
                if (type_Out(iout) .ne. 'V' .and. file_Out(iout)(1:4) .ne. 'NONE') then
                    call shd_io__outputResultMpi(Nstokes, Nx_T, Ny_T, Nz,&
                        &Npts_tot, Ncell_tot, Nsh_tot,&
                        &maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nang_do,&
                        &prp_obj%fname, fileInp_sav, fileOut_sav,&
                        &flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                        &solar, mu_solar, az_solar, skyRad,&
                        &waveNum, waveLen, units_src,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &prp_obj%xDomain, prp_obj%yDomain, xGrid_all, yGrid_all, z_Grid,&
                        &fluxTot_all, fluxDiv_all,&
                        &Nsh_out, SH_all, irad, sumRad_Out,&
                        &type_Out(iout), parms_Out(1, iout), file_Out(iout), sfc_obj, ckd_obj, sen_obj)
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
                call shd_io__outputResult(Nstokes, Nx, Ny, Nx1, Ny1, Nz, Npts_base, Npts, Ncell,&
                    Nsh, maxNl_sh, maxNm_sh, Nlm_sh, prp_obj%NlegSet, prp_obj%Nleg, prp_obj%Npha, Nmu_do, Nphi_do, Nang_do,&
                    prp_obj%fname, fileInp_sav, fileOut_sav,&
                    flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                    solar, mu_solar, az_solar, skyRad,&
                    waveNum, waveLen, units_src,&
                    splitAcc, shAcc, solAcc, maxIter, totIter,&
                    x_Grid, y_Grid, z_Grid, xyz_Grid,&
                    ptr2Tree_Cell, ptr8Grid_Cell, flag_Cell,&
                    ext, ssa, leg, ipha, temp,&
                    sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out,&
                    irad, sumRad_Out, Nsh_out, sumSH_Out, src1_Out,&
                    type_Out(iout), parms_Out(1, iout), file_Out(iout), sfc_obj, ckd_obj, sen_obj)
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Output camera radiance if there are any
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        do iout = 1, Nout
            if ((type_Out(iout) .eq. 'V')) then
                call shd_io__outputCamRad(Nstokes, Nx, Ny, Nz, Npts, Ncell,&
                    &Nsh, maxNl_sh, maxNm_sh, Nlm_sh, Nmu_do, Nphi_do, Nang_do,&
                    &prp_obj%fname, fileInp_sav, fileOut_sav,&
                    &flag_bc, ipFlag, isDeltaM, type_grid, type_src,&
                    &solar, mu_solar, az_solar, skyRad,&
                    &waveNum, waveLen, units_src,&
                    &splitAcc, shAcc, solAcc, maxIter, totIter,&
                    &ilos, sumCam_Out, parms_Out(1, iout), file_Out(iout), sfc_obj, ckd_obj, sen_obj)
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
    !      + 4.5*prp_obj%maxNpts + prp_obj%maxNptsLeg + prp_obj%Npha*prp_obj%NlegSet*(prp_obj%maxNleg + 1)&
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
    deallocate (dataRaw_z_Prp, dataRaw_temp_Prp, dataRaw_ext_Prp, dataRaw_ssa_Prp)
    deallocate (dataRaw_ipha_Prp, dataRaw_leg_Prp, dataRaw_extDir_Prp)
    if (prp_obj%isAllocated) then
        nullify(prp_obj%dataRaw_z, &
                prp_obj%dataRaw_temp, &
                prp_obj%dataRaw_ext, &
                prp_obj%dataRaw_ssa, &
                prp_obj%dataRaw_ipha, &
                prp_obj%dataRaw_leg)
        prp_obj%isAllocated = .false.
    end if

    if (isFluxOut) deallocate (sumFluxTot_Out, sumFluxDir_Out)
    if (isHROut)   deallocate (sumFluxDiv_Out)
    if (isSHOut)   deallocate (sumSH_Out)
    if (isSrcOut)  deallocate (src1_Out)

    if (allocated(dataRaw_kabs_Ckd)) then
        deallocate (dataRaw_solar_Ckd, dataRaw_delg_Ckd, dataRaw_z_Ckd, dataRaw_kabs_Ckd)
        if (ckd_obj%isAllocated) then
            nullify(ckd_obj%dataRaw_solar, ckd_obj%dataRaw_delg, ckd_obj%dataRaw_z, ckd_obj%dataRaw_kabs)
            ckd_obj%isAllocated = .false.
        end if
    end if

    if (allocated(dataRaw_brdf_Sfc)) then
        deallocate(dataRaw_brdf_Sfc, dataGrid_brdf_Sfc)
        if (sfc_obj%isAllocated) then
            nullify(sfc_obj%dataRaw_brdf, sfc_obj%dataGrid_brdf)
            sfc_obj%isAllocated = .false.
        end if
    end if

    if (allocated(dataRaw_view_Sen)) then
        deallocate (dataRaw_view_Sen)
        if (sen_obj%isAllocated) then
            nullify(sen_obj%dataRaw_view)
            sen_obj%isAllocated = .false.
        end if
    end if

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


    call shd_mpi__end(Npts, xyz_Grid, prp_obj%Nx, prp_obj%Ny, prp_obj%xStart, prp_obj%yStart, prp_obj%xDomain, prp_obj%yDomain,&
        &prpAll_obj%Nx, prpAll_obj%Ny, prp_obj%fname, runName)

end program shdom
