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
    use globals_type,  only: shd_type_prp, &
                             shd_type_sfc, &
                             shd_type_ckd, &
                             shd_type_sen, &
                             shd_type_shdo, &
                             shd_type_rtm, &
                             shd_type_grid, &
                             shd_type_src

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


    ! rtm_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_rtm) :: rtmAll_obj, rtm_obj

    real(kind=R__), allocatable  :: planck(:)
    real(kind=R__), allocatable  :: temp(:), ext(:), ssa(:)
    real(kind=R__), allocatable  :: leg(:)
    integer(kind=2), allocatable :: ipha(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! grid_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_grid) :: grid_obj
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! src_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_src) :: src_obj
    !╰────────────────────────────────────────────────────────────────────────────╯!


    real(kind=R__), allocatable    :: work(:)
    integer(kind=I__), allocatable :: work1(:)
    real(kind=R__), allocatable    :: work2(:)

    real(kind=R__), allocatable :: fluxTot(:, :), fluxDir(:)
    real(kind=R__), allocatable :: fluxDirTmp(:)

    real(kind=R__), allocatable :: fluxTot_all(:, :, :, :), fluxDiv_all(:, :, :)
    real(kind=R__), allocatable :: SH_all(:, :, :, :)
    real(kind=R__), allocatable :: xGrid_all(:), yGrid_all(:)

    ! shdo_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_shdo) :: shdo_obj
    integer(kind=I__), allocatable, target :: Nphi0_Shdo(:)
    real(kind=R__),    allocatable, target :: mu_Shdo(:), &
                                              phi_Shdo(:), &
                                              wgt_Shdo(:)
    real(kind=R__),    allocatable, target :: radI_Shdo(:, :), &
                                              srcJ_Shdo(:, :)
    integer(kind=I__), allocatable, target :: iradI_Shdo(:), &
                                              isrcJ_Shdo(:)

    real(kind=R__),    allocatable :: dsrc_Shdo(:, :)
    integer(kind=I__), allocatable :: idsrc_Shdo(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SAV
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_sav, fileOut_sav
    logical :: isNewRad, isNewGrid
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! prp_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_prp) :: prpAll_obj, prp_obj

    real(kind=R__),  allocatable, target :: z_Prp(:)
    real(kind=R__),  allocatable, target :: temp_Prp(:), &
                                            ext_Prp(:), &
                                            ssa_Prp(:)
    real(kind=R__),  allocatable, target :: extDir_Prp(:)
    real(kind=R__),  allocatable, target :: leg_Prp(:)
    integer(kind=2), allocatable, target :: ipha_Prp(:)

    real(kind=R__),  allocatable, target :: temp_PrpAll(:), &
                                            ext_PrpAll(:), &
                                            ssa_PrpAll(:)
    real(kind=R__),  allocatable, target :: leg_PrpAll(:)
    integer(kind=2), allocatable, target :: ipha_PrpAll(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! sfc_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_sfc) :: sfc_obj
    real(kind=R__), allocatable, target :: brdf_Sfc(:)
    real(kind=R__), allocatable, target :: brdf_grid_Sfc(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! ckd_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_ckd) :: ckd_obj
    integer(kind=I__) :: ig_ckd
    real(kind=R__)    :: solar0_ckd
    real(kind=R__), allocatable, target :: solar_Ckd(:), &
                                           wgt_Ckd(:), &
                                           z_Ckd(:), &
                                           kabs_Ckd(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! sen_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_sen) :: sen_obj
    real(kind=R__), allocatable, target :: view_Sen(:)
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

    integer(kind=I__) :: Npts_base, Npts, Ncell, Nsh, Ncell_base, oldNpts
    integer(kind=I__) :: maxIter, iter, totIter, iz, iout, Nout
    integer(kind=I__) :: Nrad, irad, Nx_out, Ny_out, Nlos, ilos, Nscan, Npix, Nsh_out
    integer(kind=I__) :: Ncell_tot, Npts_tot, Nsh_tot

    real(kind=R__) :: skyRad
    real(kind=R__) :: solAcc, solCrit, splitAcc, shAcc

    integer(kind=I__) :: Nproc
    logical :: isMasterProc

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
            call shd_io__nmlInp(rtmAll_obj%runName, prpAll_obj%fname, sfc_obj%fname, ckd_obj%fname,&
                fileInp_sav, fileOut_sav, sen_obj%fname,&
                shdo_obj%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtmAll_obj%Nz, shdo_obj%Nmu, shdo_obj%Nphi, rtmAll_obj%flag_bc,&
                rtmAll_obj%flag_ip, isCkd, isDeltaM, grid_obj%type,&
                src_obj%type, src_obj%avgSolar0, src_obj%avgMu, src_obj%avgAz, skyRad,&
                sfc_obj%avgTemp, sfc_obj%avgAlb, src_obj%unit, src_obj%waveNum, src_obj%avgWaveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        else
            call shd_io__userInp(rtmAll_obj%runName, prpAll_obj%fname, sfc_obj%fname, ckd_obj%fname,&
                fileInp_sav, fileOut_sav, sen_obj%fname,&
                shdo_obj%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtmAll_obj%Nz, shdo_obj%Nmu, shdo_obj%Nphi, rtmAll_obj%flag_bc,&
                rtmAll_obj%flag_ip, isCkd, isDeltaM, grid_obj%type,&
                src_obj%type, src_obj%avgSolar0, src_obj%avgMu, src_obj%avgAz, skyRad,&
                sfc_obj%avgTemp, sfc_obj%avgAlb, src_obj%unit, src_obj%waveNum, src_obj%avgWaveLen,&
                isAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        end if
    end if

    call shd_mpi__bcastInp(rtmAll_obj%runName,&
        &fileInp_sav, fileOut_sav,&
        &shdo_obj%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtmAll_obj%Nz, shdo_obj%Nmu, shdo_obj%Nphi, rtmAll_obj%flag_bc,&
        &rtmAll_obj%flag_ip, isCkd, isDeltaM, grid_obj%type,&
        &src_obj%type, src_obj%avgSolar0, src_obj%avgMu, src_obj%avgAz, skyRad,&
        &src_obj%unit, src_obj%waveNum, src_obj%avgWaveLen,&
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
    shdo_obj%Nstokes = max(1, min(4, shdo_obj%Nstokes))
    if (shdo_obj%Nstokes .eq. 2) then
        shdo_obj%Nstokes = 3
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: shdo_obj%Nstokes cannot be 2, reset shdo_obj%Nstokes to 3.'
    end if

    ! prpAll_obj%NlegSet is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (shdo_obj%Nstokes .eq. 1) then
        prpAll_obj%NlegSet = 1
        shdo_obj%NlegSet = 1
        isPol_ = .false.
    else
        prpAll_obj%NlegSet = 6
        shdo_obj%NlegSet = 6
        isPol_ = .true.
    end if

    if (shdo_obj%Nstokes .eq. 4) then
        isPol4_ = .true.
    else
        isPol4_ = .false.
    end if

    ! Make shdo_obj%maxNl and shdo_obj%maxNm from shdo_obj%Nmu and shdo_obj%Nphi (for discrete ordinate)
    !   shdo_obj%Nmu : number of discrete ordinate zenith angles in both hemispheres
    !   shdo_obj%Nphi: number of discrete ordinate azimuth angles in 2*pi
    !   shdo_obj%maxNphi0: The maximum number of azimuth angles actually used
    shdo_obj%Nmu  = max(2, 2*int((shdo_obj%Nmu + 1)/2)) ! force shdo_obj%Nmu an even number
    shdo_obj%Nphi = max(1, shdo_obj%Nphi)
    shdo_obj%maxNphi0 = shdo_obj%Nphi

    ! shdo_obj%maxNl: The maximum L of the spherical harmonic truncation
    ! shdo_obj%maxNm: The maximum M of the spherical harmonic truncation
    shdo_obj%maxNl = shdo_obj%Nmu - 1
    shdo_obj%maxNm = max(0, int(shdo_obj%Nphi/2) - 1)

    ! shdo_obj%Nlm: The number of spherical harmonic terms
    shdo_obj%Nlm = (2*shdo_obj%maxNm + 1)*(shdo_obj%maxNl + 1) - shdo_obj%maxNm*(shdo_obj%maxNm + 1)

    ! prpAll_obj%Nleg: The number of Legendre series phase function terms used in the solution procedure
    prpAll_obj%Nleg = shdo_obj%maxNl
    if (isDeltaM) prpAll_obj%Nleg = shdo_obj%maxNl + 1

    if (src_obj%type .eq. 'B') then
        isSolar_ = .true.
        isThermal_ = .true.
    else
        if (src_obj%type .eq. 'S') then
            isSolar_ = .true.
            isThermal_ = .false.
        else if (src_obj%type .eq. 'T') then
            isSolar_ = .false.
            isThermal_ = .true.
        else
            stop 'Error [shdom]: Unknown src_obj%type.'
        end if
    end if

    ! for thermal, reset solar azimuth angle
    if (src_obj%type .eq. 'T') src_obj%avgAz = 0.0d0

    ! angle grids for DO
    allocate (Nphi0_Shdo(shdo_obj%Nmu), &
              mu_Shdo(shdo_obj%Nmu), &
              phi_Shdo(shdo_obj%Nmu*shdo_obj%Nphi), &
              wgt_Shdo(shdo_obj%Nmu*shdo_obj%Nphi)); &
        memUsed = memUsed + 2.0*shdo_obj%Nmu*(1.0 + shdo_obj%Nphi)

    shdo_obj%Nphi0 => Nphi0_Shdo
    shdo_obj%mu    => mu_Shdo
    shdo_obj%phi   => phi_Shdo
    shdo_obj%wgt   => wgt_Shdo
    shdo_obj%isAllocated = .true.
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
    rtmAll_obj%Nx = max(1, rtmAll_obj%Nx)
    rtmAll_obj%Ny = max(1, rtmAll_obj%Ny)
    rtmAll_obj%Nz = max(2, rtmAll_obj%Nz)

    ! If single plane or column then force IP mode
    !   if              3D rtmAll_obj%flag_ip=0, bit representation '0000'
    !   if rtmAll_obj%Ny=1 rtmAll_obj%flag_ip=1, bit representation '0001'
    !   if rtmAll_obj%Ny=1 rtmAll_obj%flag_ip=2, bit representation '0010'
    !   if             IPA rtmAll_obj%flag_ip=3, bit representation '0011'
    if (rtmAll_obj%Nx .eq. 1) then
        ! '0000'(3D,   0) + '0001'(rtmAll_obj%Nx=1, 1) = '0001' (rtmAll_obj%Nx=1, 1)
        rtmAll_obj%flag_ip = ior(rtmAll_obj%flag_ip, 1)
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: rtmAll_obj%Nx=1, reset rtmAll_obj%flag_ip to 1.'
    end if

    if (rtmAll_obj%Ny .eq. 1) then
        ! '0001'(rtmAll_obj%Nx=1, 1) + '0010'(rtmAll_obj%Ny=1, 2) = '0011' (IPA,  3)
        rtmAll_obj%flag_ip = ior(rtmAll_obj%flag_ip, 2)
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: rtmAll_obj%Ny=1, reset rtmAll_obj%flag_ip to 2 (or 3).'
    end if

    if ((rtmAll_obj%Nx .eq. 1) .and. (rtmAll_obj%Ny .eq. 1)) then
        is1D1_ = .true.
    else
        is1D1_ = .false.
    end if

    !==================

    ! Overwrite prpAll_obj%dx if only 1 grid in X direction
    if (prpAll_obj%Nx .eq. 1) then
        prpAll_obj%dx = sfc_obj%Nx * sfc_obj%dx
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: prpAll_obj%Nx=1, reset prpAll_obj%dx to match surface domain.'
    end if

    ! overwrite prpAll_obj%dy if only 1 grid in Y direction
    if (prpAll_obj%Ny .eq. 1) then
        prpAll_obj%dy = sfc_obj%Ny * sfc_obj%dy
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: prpAll_obj%Ny=1, reset prpAll_obj%dy to match surface domain.'
    end if

    prpAll_obj%xDomain = prpAll_obj%Nx * prpAll_obj%dx
    prpAll_obj%yDomain = prpAll_obj%Ny * prpAll_obj%dy

    ! Choose the option to compute higher order radiance SH terms
    if  ((prpAll_obj%Nx .eq. 1) .and. (prpAll_obj%Ny .eq. 1) .and.&
        &(sfc_obj%Nx .le. 32) .and. (sfc_obj%Ny .le. 32)) then
        isHighOrder_ = .true.
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: 1D calculation, set isHighOrder_ to .TRUE.'
    else
        isHighOrder_ = .false.
    end if

    Nsh_out = 4
    if (isHighOrder_) Nsh_out = 5
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! If MPI, divide domain into mini-domains
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and Nx,Ny to the subdomain size, and output
    !   the X and Y starting position of this subdomain. Set bits
    !   2 and 3 in rtmAll_obj%flag_bc to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    rtm_obj = rtmAll_obj; prp_obj = prpAll_obj
    call shd_mpi__map(rtmAll_obj, prpAll_obj, rtm_obj, prp_obj)

    ! rtm_obj%flag_bc bit-representation '...0001' isOpenX_
    if (btest(rtm_obj%flag_bc, 0)) then
        isOpenX_ = .true.
    else
        isOpenX_ = .false.
    end if

    ! rtm_obj%flag_bc bit-representation '...001x' isOpenYx_
    if (btest(rtm_obj%flag_bc, 1)) then
        isOpenYx_ = .true.
    else
        isOpenYx_ = .false.
    end if

    ! rtm_obj%flag_bc bit-representation '...01xx' isMpiX_
    if (btest(rtm_obj%flag_bc, 2)) then
        isMpiX_ = .true.
    else
        isMpiX_ = .false.
    end if

    ! rtm_obj%flag_bc bit-representation '...1xxx' isMpiY_
    if (btest(rtm_obj%flag_bc, 3)) then
        isMpiY_ = .true.
    else
        isMpiY_ = .false.
    end if

    if (isMpiX_ .or. isMpiY_) then
        isMpi_ = .true.
    else
        isMpi_ = .false.
    end if

    ! btest(rtm_obj%flag_ip, 0) checks whether rtm_obj%flag_ip=1 or rtm_obj%flag_ip=0 ('...0001' [IPA in X] or '...0000' [3D in X])
    if (btest(rtm_obj%flag_ip, 0)) then
        is1DX_ = .true.
    else
        is1DX_ = .false.
    end if

    ! btest(rtm_obj%flag_ip, 1) checks whether rtm_obj%flag_ip=2 [or 3] or rtm_obj%flag_ip=0 [or 1]
    !   ('...0010' [or '...0011', IPA in Y] or '...0000' [or '...0001', 3D in Y])
    if (btest(rtm_obj%flag_ip, 1)) then
        is1DYx_ = .true.
    else
        is1DYx_ = .false.
    end if

    ! btest(rtm_obj%flag_ip, 1) checks whether rtm_obj%flag_ip=7 [or 8]
    !   ('...0110' [or '...0111'], 3D direct beam but IP diffuse radiative transfer)
    if (btest(rtm_obj%flag_ip, 2)) then
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
    !   btest(rtm_obj%flag_bc, 0) checks whether rtm_obj%flag_bc=1 or rtm_obj%flag_bc=0 ('...0001' [open in X] or '...0000' [periodic in X])
    grid_obj%Nx = rtm_obj%Nx + 1
    if (is1DX_ .or. isOpenX_ .or. isMpiX_) grid_obj%Nx = rtm_obj%Nx

    ! btest(rtm_obj%flag_bc, 1) checks whether rtm_obj%flag_bc=2 [or 3] or rtm_obj%flag_bc=0 [or 1]
    !   ('...0010' [or '...0011', open in Y] or '...0000' [or '...0001', periodic in Y])
    grid_obj%Ny = rtm_obj%Ny + 1
    if (is1DYx_ .or. isOpenYx_ .or. isMpiY_) grid_obj%Ny = rtm_obj%Ny

    ! Number of base points
    Npts_base = grid_obj%Nx*grid_obj%Ny*rtm_obj%Nz

    grid_obj%Nz = rtm_obj%Nz

    ! Calculate the number of base grid cells depending on the boundary condition flag
    Ncell_base = (rtm_obj%Nz - 1)*(rtm_obj%Nx + ibits(rtm_obj%flag_bc, 0, 1) - ibits(rtm_obj%flag_bc, 2, 1))&
                                 *(rtm_obj%Ny + ibits(rtm_obj%flag_bc, 1, 1) - ibits(rtm_obj%flag_bc, 3, 1))
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Read prp_obj%fname and distribute property object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Determine the size of the property grid arrays so can allocate
    prp_obj%maxNpts = prp_obj%Nx * prp_obj%Ny * prp_obj%Nz
    if (prp_obj%Npha .eq. 0) prp_obj%maxNptsLeg = prp_obj%maxNpts * (prp_obj%maxNleg + 1) * prp_obj%NlegSet

    allocate (z_Prp(prp_obj%Nz))
    allocate (temp_Prp(prp_obj%maxNpts), &
              ext_Prp( prp_obj%maxNpts), &
              ssa_Prp( prp_obj%maxNpts))
    allocate (extDir_Prp(prp_obj%maxNpts), &
              ipha_Prp(prp_obj%maxNpts), &
              leg_Prp( prp_obj%maxNptsLeg)); &
        memUsed = memUsed + prp_obj%Nz + 6.0*prp_obj%maxNpts + prp_obj%maxNptsLeg

    prp_obj%z      => z_Prp
    prp_obj%temp   => temp_Prp
    prp_obj%ext    => ext_Prp
    prp_obj%extDir => extDir_Prp
    prp_obj%ssa    => ssa_Prp
    prp_obj%ipha   => ipha_Prp
    prp_obj%leg    => leg_Prp
    prp_obj%isAllocated = .true.

    if (isMpi_) then

        ! Read the properties of the medium to the master processor and
        !   distribute to all the processors
        if (isMasterProc) then

            allocate (temp_PrpAll(prpAll_obj%maxNpts), &
                      ext_PrpAll( prpAll_obj%maxNpts), &
                      ssa_PrpAll( prpAll_obj%maxNpts))
            allocate (ipha_PrpAll(prpAll_obj%maxNpts))

            prpAll_obj%z    => z_Prp
            prpAll_obj%leg  => leg_Prp
            prpAll_obj%temp => temp_PrpAll
            prpAll_obj%ext  => ext_PrpAll
            prpAll_obj%ssa  => ssa_PrpAll
            prpAll_obj%ipha => ipha_PrpAll
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
            deallocate(temp_PrpAll, &
                       ext_PrpAll, &
                       ssa_PrpAll, &
                       ipha_PrpAll)
            if (prpAll_obj%isAllocated) then
                nullify(prpAll_obj%z, &
                        prpAll_obj%temp, &
                        prpAll_obj%ext, &
                        prpAll_obj%ssa, &
                        prpAll_obj%ipha, &
                        prpAll_obj%leg)
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
        call shd_io__checkInp(shdo_obj%Nmu, shdo_obj%Nphi, isDeltaM,&
            &grid_obj%type, src_obj%type, src_obj%unit, src_obj%avgSolar0, src_obj%avgMu,&
            &src_obj%waveNum, src_obj%avgWaveLen,&
            &splitAcc, shAcc, prp_obj, sfc_obj)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sfc_obj%fname and broadcast the surface object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in the surface properties
    if (sfc_obj%fname(1:4) .ne. 'NONE') then

        allocate (brdf_Sfc(sfc_obj%maxNpar * sfc_obj%maxNpts)); &
            memUsed = memUsed + sfc_obj%maxNpar * sfc_obj%maxNpts

        sfc_obj%brdf => brdf_Sfc
        sfc_obj%isAllocated = .true.


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
            call shd_io__readCkdSize(src_obj%waveNum, ckd_obj)
        end if

        ! Broadcast the k-distribution object
        call shd_mpi__bcastCkdSize(ckd_obj)

        allocate (solar_Ckd(ckd_obj%Ng), wgt_Ckd(ckd_obj%Ng), z_Ckd(ckd_obj%Nz))
        allocate (kabs_Ckd(ckd_obj%Nz, ckd_obj%Ng)); &
            memUsed = memUsed + ckd_obj%Nz*(2.0+ckd_obj%Ng) + 2.0*ckd_obj%Ng

        ckd_obj%solar => solar_Ckd
        ckd_obj%wgt  => wgt_Ckd
        ckd_obj%z     => z_Ckd
        ckd_obj%kabs  => kabs_Ckd
        ckd_obj%isAllocated = .true.

        if (isMasterProc) then
            ckd_obj%maxNg = ckd_obj%Ng; ckd_obj%maxNz = ckd_obj%Nz
            call shd_io__readCkd(src_obj%waveNum, ckd_obj)
        end if

        call shd_mpi__bcastCkd(ckd_obj)

    else
        ckd_obj%Ng = 1
        ckd_obj%Nz = 0

        allocate (solar_Ckd(ckd_obj%Ng), wgt_Ckd(ckd_obj%Ng)); &
            memUsed = memUsed + 2.0*ckd_obj%Ng

        solar_Ckd(1) = 1.0
        wgt_Ckd(1) = 1.0

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

        allocate (view_Sen(sen_obj%Npar*sen_obj%Npts)); &
            memUsed = memUsed + sen_obj%Npar*sen_obj%Npts

        sen_obj%view => view_Sen
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
            Nrad = Nrad + shdo_obj%Nstokes*Nx_out*Ny_out*nint(parms_Out(6, iout))
        else if (type_Out(iout) .eq. 'V') then
            isCamOut = .true.
            if (nint(parms_Out(1, iout)) .eq. 1) then
                Nlos = Nlos + shdo_obj%Nstokes*nint(parms_Out(10, iout))*nint(parms_Out(11, iout))
            else if (nint(parms_Out(1, iout)) .eq. 2) then
                Nscan = 1 + sqrt((parms_Out(4, iout) - parms_Out(7, iout))**2&
                    + (parms_Out(5, iout) - parms_Out(8, iout))**2&
                    + (parms_Out(6, iout) - parms_Out(9, iout))**2)/parms_Out(10, iout)
                Npix = 1 + abs(parms_Out(12, iout) - parms_Out(11, iout))/parms_Out(13, iout)
                Nlos = Nlos + shdo_obj%Nstokes*Nscan*Npix
            ! user-specified sensor
            else if (nint(parms_Out(1, iout)) .eq. 3) then
                Nlos = Nlos + shdo_obj%Nstokes*nint(parms_Out(2, iout))
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
            &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, rtm_obj%runName)
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

    ! x, y, z grid
    memUsed = memUsed + (rtm_obj%Nx+1) + (rtm_obj%Ny+1) + (rtm_obj%Nz)

    if (prp_obj%Npha .gt. 0) then
        memUsed = memUsed + prp_obj%Npha * prp_obj%NlegSet * (prp_obj%maxNleg + 1) ! offset phase table
    else
        NmaxIG = NmaxIG + prp_obj%NlegSet*(prp_obj%maxNleg + 1)
    end if

    if (sfc_obj%fname(1:4) .ne. 'NONE') then
        if (sfc_obj%type(2:2) .eq. 'L') then
            NmaxIG = real(3)/real(rtm_obj%Nz) * (sfc_obj%maxNpar + 2.0 + shdo_obj%Nstokes*2.0)
        else
            NmaxIG = real(3)/real(rtm_obj%Nz) * (sfc_obj%maxNpar + 2.0 + shdo_obj%Nstokes*(2.0+ shdo_obj%Nmu*shdo_obj%maxNphi0/2))
        end if
    end if

    if (isFluxOut)  NmaxIG = NmaxIG + 3.0
    if (isHROut)    NmaxIG = NmaxIG + 1.0
    if (isSHOut)    NmaxIG = NmaxIG + Nsh_out
    if (isSrcOut)   NmaxIG = NmaxIG + shdo_obj%Nstokes

    NmaxIG = NmaxIG + 16.5*CELL_TO_POINT_RATIO ! ptr8Grid_Cell[8.0], ptr6Cell_Cell[6.0], ptr2Tree_Cell[2.0], flag_Cell[0.5]

    NmaxIG = NmaxIG + 6.0   ! xyz_Grid[3.0], shdo_obj%iradI[1.0], shdo_obj%isrcJ[1.0], idsrc_Shdo[1.0]
    memUsed = memUsed + 4.0 ! offset shdo_obj%iradI[+2], shdo_obj%isrcJ[+1], idsrc_Shdo[+1]

    if (isAccel) then
        NmaxIG = NmaxIG + shdo_obj%Nstokes*(1.0+3.0*NUM_SH_TERM_FACTOR*shdo_obj%Nlm) ! radiance, source, delta source
    else
        NmaxIG = NmaxIG + shdo_obj%Nstokes*(1.0+2.0*NUM_SH_TERM_FACTOR*shdo_obj%Nlm) ! radiance, source
    end if

    NmaxIG = NmaxIG + 3.0 ! fluxTot[2.0], fluxDir[1.0]

    memUsed = memUsed + (Nrad+Nlos) ! offset sumRad_Out[Nrad] and sumCam_Out[Nlos]

    NmaxIG = NmaxIG + (shdo_obj%Nstokes*(max(shdo_obj%maxNphi0, 4)+1.0) + 8.0) !work, work1, work2

    memUsed = memUsed + Nrad + Nlos ! work2 (over-count but memory safe)
    !╰──────────────────────────────────────────────────────────────╯!

    ! memUsed is now more accurate than before
    ! before:
    !   memUsed = shdo_obj%Nmu*(2 + 2*shdo_obj%Nphi + 2*shdo_obj%Nlm + 2*33*32)&            ! SHDO allocation
    !       + 4.5*prp_obj%maxNpts + prp_obj%maxNptsLeg + prp_obj%Npha*prp_obj%NlegSet*(prp_obj%maxNleg + 1) ! phase allocation
    ! now: update memUsed after every allocation step
    !
    ! NmaxIG is now more accurate than before
    ! before:
    !   NmaxIG = 28.0 + 16.5*CELL_TO_POINT_RATIO&
    !       + shdo_obj%maxNphi0*shdo_obj%Nstokes + NUM_SH_TERM_FACTOR*shdo_obj%Nstokes*shdo_obj%Nlm*NbigArr
    ! now: update NmaxIG from every allocation step

    memNeed = ADAPT_GRID_FACTOR*Npts_base*NmaxIG

    if (isHighOrder_) then
        scaleFactor = min(1.5, ((MAX_TOTAL_MB*1024**2)/4 - memUsed)/memNeed)
    else
        scaleFactor = min(1.0, ((MAX_TOTAL_MB*1024**2)/4 - memUsed)/memNeed)
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
    maxIV = nint(NUM_SH_TERM_FACTOR*shdo_obj%Nlm*maxIG)

    ! max words for temporary discrete ordinate azimuthal array
    maxIDO = maxIG*max(shdo_obj%maxNphi0, 4)

    if (4.0*(maxIV + maxIG)*shdo_obj%Nstokes .gt. huge(maxIV)) then
        write (6, '(A,I12)') 'Error [shdom]: Size of big SH arrays (maxIV) probably exceeds max integer number of bytes:', maxIV
        stop
    end if

    if (4.0*8.0*maxIC .gt. huge(maxIV)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Size of ptr8Grid_Cell array (8*maxIC) probably exceeds max integer number of bytes:', 8*maxIC
        stop
    end if

    maxNpts_bc = maxIG*3/rtm_obj%Nz
    if (sfc_obj%isLambertian) then
        maxNrad_bc = 2 * maxNpts_bc
    else
        maxNrad_bc = (2 + shdo_obj%Nmu*shdo_obj%maxNphi0/2)*maxNpts_bc
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
    ! x, y, z grids ???
    allocate (x_Grid(rtm_obj%Nx + 1), y_Grid(rtm_obj%Ny + 1), z_Grid(rtm_obj%Nz))

    ! temperature and planck radiance (maybe we don't need this for shortwave)
    allocate (temp(maxIG), planck(maxIG))

    ! extinction, single scattering albedo, and phase function indices
    allocate (ext(maxIG), ssa(maxIG), ipha(maxIG))

    ! legendre coefficients
    allocate (leg(maxIGL))

    ! boundary points pointers and radiances
    !   ptr_Bc(:, 1): pointers to grid points that are on top
    !   ptr_Bc(:, 2): pointers to grid points that are on bottom
    !   rad_Bc(:[shdo_obj%Nstokes], 1:NptsTop_bc): isotropic top radiances
    !   rad_Bc(:[shdo_obj%Nstokes], NptsTop_bc+1:NptsTop_bc+NptsBot_bc): upwelling bottom radiances
    !   rad_Bc(:[shdo_obj%Nstokes], NptsTop_bc+NptsBot_bc+1:NptsTop_bc+NptsBot_bc+NptsBot_bc*shdo_obj%Nang/2):
    !         downwelling bottom radiances at discrete ordinates
    if (sfc_obj%fname(1:4) .ne. 'NONE') then
        allocate (brdf_grid_Sfc(sfc_obj%maxNpar*maxNpts_bc))
        sfc_obj%brdf_grid => brdf_grid_Sfc
        sfc_obj%isAllocated = .true.
    end if

    allocate (ptr_Bc(2*maxNpts_bc), rad_Bc(shdo_obj%Nstokes, maxNrad_bc))

    ! output
    if (isFluxOut) allocate (sumFluxTot_Out(2, maxIG), sumFluxDir_Out(maxIG))
    if (isHROut)   allocate (sumFluxDiv_Out(maxIG))
    if (isSHOut)   allocate (sumSH_Out(maxIG*Nsh_out))
    if (isSrcOut)  allocate (src1_Out(maxIG*shdo_obj%Nstokes))

    ! pointer to 8 grid points belonging to a cell
    ! pointer to 6 neighboring cells of a cell
    ! pointer to parent and children cells in tree structure
    allocate (ptr8Grid_Cell(8, maxIC), ptr6Cell_Cell(6, maxIC), ptr2Tree_Cell(2, maxIC))

    ! a collection of binary flags for the cells
    allocate (flag_Cell(maxIC))

    ! x,y,z coordinate values of grid point locations
    allocate (xyz_Grid(3, maxIG))

    ! SH related array
    allocate (iradI_Shdo(maxIG + 2), isrcJ_Shdo(maxIG + 1))
    ! radiance and source
    allocate (radI_Shdo(shdo_obj%Nstokes, maxIV + maxIG), srcJ_Shdo(shdo_obj%Nstokes, maxIV))
    shdo_obj%iradI => iradI_Shdo
    shdo_obj%isrcJ => isrcJ_Shdo
    shdo_obj%radI  => radI_Shdo
    shdo_obj%srcJ  => srcJ_Shdo

    ! delta source
    if (isAccel) then
        allocate (idsrc_Shdo(maxIG + 1), dsrc_Shdo(shdo_obj%Nstokes, maxIV))
    end if

    ! fluxes
    allocate (fluxTot(2, maxIG), fluxDir(maxIG))

    ! output radiances to be summed over g
    allocate (sumRad_Out(Nrad), sumCam_Out(Nlos))

    ! work: work array for path integration; radiance field, Nazimuth*Npts
    allocate (work(shdo_obj%Nstokes*maxIDO))

    ! work1: work array for path integration; grid point sweeping order
    allocate (work1(8*maxIG))

    ! work2: work array for path integration; single angle radiance
    allocate (work2(max(Nlos, Nrad, shdo_obj%Nstokes*maxIG)))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read in previous run from binary file if desired
    call shd_io__restoreState(fileInp_sav, shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
        &isNewRad, isNewGrid, x_Grid, y_Grid, z_Grid,&
        &Npts, Ncell, xyz_Grid, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell,&
        &flag_Cell, fluxTot, shdo_obj%isrcJ, shdo_obj%srcJ, shdo_obj%iradI, shdo_obj%radI)


    ! Check if we can use pre-calculated "state", if not, generate new grids
    !   isNewGrid is from shd_io__restoreState call
    if (isNewGrid) then
        ! Make the internal base grid lines
        call shd_init__makeGrid(rtm_obj%flag_bc, grid_obj%type, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz,&
            &x_Grid, y_Grid, z_Grid, prp_obj)
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
    if (verbose_ .ge. 5) &
        write (6, '(A,"<",A,">")') 'Message [shdom]: Starting solution iterations for SHDOM run: ', &
            trim(adjustl(rtm_obj%runName))

    totIter = 0
    oldNpts = 0
    src_obj%avgSolar = 0.0d0

    ! Loop over the k-distribution g's from high to low absorption
    CKD_LOOP: do ig_ckd = ckd_obj%Ng, 1, -1

        if (verbose_ .ge. 5) write (6, '(A,I2,A,I2,A)') 'Message [shdom]: Starting ', ckd_obj%Ng-ig_ckd+1, '/', ckd_obj%Ng, ' (ig/Ng) ...'

        ! solar flux for given g of CKD
        solar0_ckd = solar_Ckd(ig_ckd)*src_obj%avgSolar0*abs(src_obj%avgMu)

        ! weighted average solar flux
        src_obj%avgSolar = src_obj%avgSolar + solar_Ckd(ig_ckd)*wgt_Ckd(ig_ckd)*src_obj%avgSolar0*abs(src_obj%avgMu)

        ! If doing k-distribution then get the gas absorption profile.
        ckd_obj%ig = ig_ckd

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (isBaseOut .or. ((ig_ckd.eq. ckd_obj%Ng) .and. isNewGrid)) then
            call shd_init__makeCell(rtm_obj%flag_bc, rtm_obj%flag_ip,&
                &rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, grid_obj%Nx, grid_obj%Ny, Npts, Ncell,&
                &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell)
            Ncell_base = Ncell
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in source (shdo_obj%srcJ), radiance (shdo_obj%radI), fluxes, direct Flux.
        call shd_solv__solveRte(&
            &maxIV, maxIC, maxIG, maxIDO, isNewRad,&
            &isDeltaM, isHighOrder_,&
            &solar0_ckd, skyRad,&
            &isAccel, solAcc, maxIter, solCrit, iter,&
            &splitAcc, shAcc, x_Grid, y_Grid, z_Grid,&
            &temp, planck, ext, ssa, leg, ipha,&
            &maxNpts_bc, maxNrad_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
            &Npts, xyz_Grid,&
            &Ncell, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
            &idsrc_Shdo, work, work1, work2,&
            &dsrc_Shdo, fluxTot, fluxDir,&
            &rtm_obj, grid_obj, shdo_obj, src_obj,&
            &prp_obj, sfc_obj, ckd_obj)
        isNewRad = .false.
        totIter = totIter + iter

        ! Compute the output and put in arrays. First the radiance.
        if (isRadOut) then
            irad = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'R') then
                    if (isMpi_) then
                        call shd_mpi__computeRad(shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, Npts, Ncell,&
                            &shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
                            &shdo_obj%Nmu, shdo_obj%maxNphi0, shdo_obj%Nphi0, shdo_obj%mu, shdo_obj%phi, shdo_obj%wgt,&
                            &rtm_obj%flag_bc, rtm_obj%flag_ip,&
                            &src_obj%type, isDeltaM, src_obj%avgMu, src_obj%avgAz,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &shdo_obj%isrcJ, shdo_obj%srcJ, work, work(1 + Npts*shdo_obj%Nstokes),&
                            &parms_Out(1, iout), irad, work2, prp_obj, sfc_obj, sen_obj)
                    else
                        call shd_rad__computeRad(shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, Npts, Ncell,&
                            &shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
                            &shdo_obj%Nmu, shdo_obj%maxNphi0, shdo_obj%Nphi0, shdo_obj%mu, shdo_obj%phi, shdo_obj%wgt,&
                            &rtm_obj%flag_bc, rtm_obj%flag_ip, src_obj%type, isDeltaM, src_obj%avgMu, src_obj%avgAz,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &shdo_obj%isrcJ, shdo_obj%srcJ, work(1), work(1 + Npts*shdo_obj%Nstokes),&
                            &parms_Out(1, iout), irad, work2, prp_obj, sfc_obj, sen_obj)
                    end if
                end if
            end do
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, wgt_Ckd(ig_ckd), Nrad, work2, sumRad_Out)
        end if

        if (isCamOut) then
            ilos = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'V') then
                    if (isMpi_) then
                        call shd_mpi__computeCamRad(shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, Npts, Ncell,&
                            &shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
                            &shdo_obj%Nmu, shdo_obj%maxNphi0, shdo_obj%Nphi0, shdo_obj%mu, shdo_obj%phi, shdo_obj%wgt,&
                            &rtm_obj%flag_bc, rtm_obj%flag_ip,&
                            &src_obj%type, isDeltaM, src_obj%avgMu, src_obj%avgAz,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &shdo_obj%isrcJ, shdo_obj%srcJ, parms_Out(1, iout), ilos, work2,&
                            &prp_obj, sfc_obj, sen_obj)
                    else
                        call shd_rad__computeCamRad(shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, Npts, Ncell,&
                            &shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
                            &shdo_obj%Nmu, shdo_obj%maxNphi0, shdo_obj%Nphi0, shdo_obj%mu, shdo_obj%phi, shdo_obj%wgt,&
                            &rtm_obj%flag_bc, rtm_obj%flag_ip, src_obj%type, isDeltaM, src_obj%avgMu, src_obj%avgAz,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &skyRad, src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                            &x_Grid, y_Grid, z_Grid, xyz_Grid,&
                            &ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell, flag_Cell,&
                            &ext, ssa, leg, ipha, fluxDir, fluxTot,&
                            &shdo_obj%isrcJ, shdo_obj%srcJ, parms_Out(1, iout), ilos, work2,&
                            &prp_obj, sfc_obj, sen_obj)
                    end if
                end if
            end do
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, wgt_Ckd(ig_ckd), Nlos, work2, sumCam_Out)
        end if

        if (isSrcOut) then
            ! Compute source function throughout grid for this one angle
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'J') then
                    mu_out  = parms_Out(2, iout)
                    phi_out = parms_Out(3, iout)*acos(-1.0)/180.0
                    call shd_src__computeSrcOneDir(shdo_obj%Nstokes, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
                        &Npts, isDeltaM, mu_out, phi_out,&
                        &src_obj%type, src_obj%avgMu, src_obj%avgAz, ssa, leg, ipha,&
                        &fluxDir, shdo_obj%isrcJ, shdo_obj%srcJ, src1_Out, prp_obj)
                end if
            end do
        end if

        if (isBaseOut) then
            ! If outputting only the base grid then throw out the other points
            call shd_init__makeCell(rtm_obj%flag_bc, rtm_obj%flag_ip,&
                &rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, grid_obj%Nx, grid_obj%Ny,&
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
            if ((src_obj%type.ne. 'T') .and. isDeltaM) then

                allocate(fluxDirTmp(Npts))

                if (isMpi_) then
                    call shd_mpi__makeDirect(1, Npts, rtm_obj%flag_bc, rtm_obj%flag_ip, .false., shdo_obj%maxNl,&
                        &solar0_ckd, src_obj%avgMu, src_obj%avgAz,&
                        &xyz_Grid, rtm_obj%Nx, x_Grid, rtm_obj%Ny, y_Grid, fluxDirTmp, prp_obj, ckd_obj)
                else
                    call shd_rad__makeDirect(Npts, rtm_obj%flag_bc, rtm_obj%flag_ip, .false., shdo_obj%maxNl,&
                        &solar0_ckd, src_obj%avgMu, src_obj%avgAz,&
                        &xyz_Grid, fluxDirTmp, prp_obj, ckd_obj)
                end if

                fluxTot(1, 1:Npts) = fluxTot(1, 1:Npts) + fluxDir(1:Npts) - fluxDirTmp(1:Npts)
                fluxDir(1:Npts) = fluxDirTmp(1:Npts)

                deallocate(fluxDirTmp)

            end if

            call shd_util__sumData(ckd_obj%Ng-ig_ckd, wgt_Ckd(ig_ckd), 2*Npts, fluxTot, sumFluxTot_Out)
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, wgt_Ckd(ig_ckd), Npts, fluxDir, sumFluxDir_Out)
        end if

        if (isHROut) then
            ! Compute the net flux divergence at every point
            call shd_rad__computeNetFluxDiv(shdo_obj%Nstokes, Npts, shdo_obj%iradI, src_obj%type, src_obj%avgMu,&
                &ext, ssa, planck, fluxDir, shdo_obj%radI, work)
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, wgt_Ckd(ig_ckd), Npts, work, sumFluxDiv_Out)
        end if

        if (isSHOut) then
            call shd_shdo__computeSH(Nsh_out, shdo_obj%Nstokes, Npts, src_obj%type, src_obj%avgMu, src_obj%avgAz,&
                &fluxDir, shdo_obj%iradI, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%radI, work)
            call shd_util__sumData(ckd_obj%Ng-ig_ckd, wgt_Ckd(ig_ckd), Nsh_out*Npts, work, sumSH_Out)
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
    Nsh = shdo_obj%isrcJ(Npts + 1)
    irad = 0
    ilos = 0
    if (isMpi_) then

        if (isMasterProc) then
            if (isFluxOut) allocate (fluxTot_all(3, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))
            if (isHROut)   allocate (fluxDiv_all(rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))
            if (isSHOut)   allocate (SH_all(Nsh_out, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))
            allocate (xGrid_all(rtmAll_obj%Nx), yGrid_all(rtmAll_obj%Ny))
        end if

        call shd_mpi__gatherOut(isFluxOut, isHROut, isSHOut,&
            &rtm_obj%flag_ip, rtm_obj%flag_bc, Npts_base, rtmAll_obj%Nx, rtmAll_obj%Ny, rtm_obj%Nz, Nsh_out,&
            &sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out, sumSH_Out,&
            &xGrid_all, yGrid_all,&
            &fluxTot_all, fluxDiv_all, SH_all,&
            &Ncell, Npts, Nsh, Ncell_tot, Npts_tot, Nsh_tot, prp_obj)

        if (isMasterProc) then

            do iout = 1, Nout
                if (type_Out(iout) .ne. 'V' .and. file_Out(iout)(1:4) .ne. 'NONE') then
                    call shd_io__outputResultMpi(shdo_obj%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtm_obj%Nz,&
                        &Npts_tot, Ncell_tot, Nsh_tot,&
                        &shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm, shdo_obj%Nmu, shdo_obj%Nphi, shdo_obj%Nang,&
                        &fileInp_sav, fileOut_sav,&
                        &rtm_obj%flag_bc, rtm_obj%flag_ip, isDeltaM, grid_obj%type, src_obj%type,&
                        &src_obj%avgSolar, src_obj%avgMu, src_obj%avgAz, skyRad,&
                        &src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &xGrid_all, yGrid_all, z_Grid,&
                        &fluxTot_all, fluxDiv_all,&
                        &Nsh_out, SH_all, irad, sumRad_Out,&
                        &type_Out(iout), parms_Out(1, iout), file_Out(iout),&
                        &prp_obj, sfc_obj, ckd_obj, sen_obj)
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
                call shd_io__outputResult(shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, grid_obj%Nx, grid_obj%Ny, rtm_obj%Nz, Npts_base, Npts, Ncell,&
                    Nsh, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm, shdo_obj%Nmu, shdo_obj%Nphi, shdo_obj%Nang,&
                    fileInp_sav, fileOut_sav,&
                    rtm_obj%flag_bc, rtm_obj%flag_ip, isDeltaM, grid_obj%type, src_obj%type,&
                    src_obj%avgSolar, src_obj%avgMu, src_obj%avgAz, skyRad,&
                    src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                    splitAcc, shAcc, solAcc, maxIter, totIter,&
                    x_Grid, y_Grid, z_Grid, xyz_Grid,&
                    ptr2Tree_Cell, ptr8Grid_Cell, flag_Cell,&
                    ext, ssa, leg, ipha, temp,&
                    sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out,&
                    irad, sumRad_Out, Nsh_out, sumSH_Out, src1_Out,&
                    type_Out(iout), parms_Out(1, iout), file_Out(iout),&
                    prp_obj, sfc_obj, ckd_obj, sen_obj)
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Output camera radiance if there are any
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        do iout = 1, Nout
            if ((type_Out(iout) .eq. 'V')) then
                call shd_io__outputCamRad(shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, Npts, Ncell,&
                    &Nsh, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm, shdo_obj%Nmu, shdo_obj%Nphi, shdo_obj%Nang,&
                    &fileInp_sav, fileOut_sav,&
                    &rtm_obj%flag_bc, rtm_obj%flag_ip, isDeltaM, grid_obj%type, src_obj%type,&
                    &src_obj%avgSolar, src_obj%avgMu, src_obj%avgAz, skyRad,&
                    &src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit,&
                    &splitAcc, shAcc, solAcc, maxIter, totIter,&
                    &ilos, sumCam_Out, parms_Out(1, iout), file_Out(iout),&
                    &prp_obj, sfc_obj, ckd_obj, sen_obj)
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! save state
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (maxIter .gt. 0) then
        call shd_io__saveState(fileOut_sav, shdo_obj%Nstokes, rtm_obj%Nx, rtm_obj%Ny, rtm_obj%Nz, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm,&
            &work, x_Grid, y_Grid, z_Grid,&
            &Npts, Ncell, xyz_Grid, ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell,&
            &flag_Cell, fluxTot, shdo_obj%isrcJ, shdo_obj%srcJ, shdo_obj%iradI, shdo_obj%radI)
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
    !  MAXMB_OUT = 4*(shdo_obj%Nmu*(2 + 2*shdo_obj%Nphi + 2*shdo_obj%Nlm + 2*33*32)&
    !      + 4.5*prp_obj%maxNpts + prp_obj%maxNptsLeg + prp_obj%Npha*prp_obj%NlegSet*(prp_obj%maxNleg + 1)&
    !      + 16.5*Ncell + Npts*(28 + shdo_obj%maxNphi0*shdo_obj%Nstokes)&
    !      + shdo_obj%Nstokes*Nsh*NbigArr)/1024**2
    ! now:
    !   use NmaxIG
    MAXMB_OUT = 4*(memUsed + 16.5*Ncell&
        + Npts*(NmaxIG-16.5*CELL_TO_POINT_RATIO))/1024**2

    ADAPTGRIDFAC_OUT = float(Npts)/Npts_base
    SHTERM_FAC_OUT = float(Nsh)/(shdo_obj%Nlm*Npts)
    CELL_POINT_OUT = float(Ncell)/Npts
    if (verbose_ .ge. 5) write (6, '(A)')'Message [shdom]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    if (verbose_ .ge. 5) write (6, '(26X,F9.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Deallocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    deallocate (z_Prp, temp_Prp, ext_Prp, ssa_Prp)
    deallocate (ipha_Prp, leg_Prp, extDir_Prp)
    if (prp_obj%isAllocated) then
        nullify(prp_obj%z, &
                prp_obj%temp, &
                prp_obj%ext, &
                prp_obj%ssa, &
                prp_obj%ipha, &
                prp_obj%leg)
        prp_obj%isAllocated = .false.
    end if

    if (isFluxOut) deallocate (sumFluxTot_Out, sumFluxDir_Out)
    if (isHROut)   deallocate (sumFluxDiv_Out)
    if (isSHOut)   deallocate (sumSH_Out)
    if (isSrcOut)  deallocate (src1_Out)

    if (allocated(kabs_Ckd)) then
        deallocate (solar_Ckd, wgt_Ckd, z_Ckd, kabs_Ckd)
        if (ckd_obj%isAllocated) then
            nullify(ckd_obj%solar, ckd_obj%wgt, ckd_obj%z, ckd_obj%kabs)
            ckd_obj%isAllocated = .false.
        end if
    end if

    if (allocated(brdf_Sfc)) then
        deallocate(brdf_Sfc, brdf_grid_Sfc)
        if (sfc_obj%isAllocated) then
            nullify(sfc_obj%brdf, sfc_obj%brdf_grid)
            sfc_obj%isAllocated = .false.
        end if
    end if

    if (allocated(view_Sen)) then
        deallocate (view_Sen)
        if (sen_obj%isAllocated) then
            nullify(sen_obj%view)
            sen_obj%isAllocated = .false.
        end if
    end if

    deallocate (ptr_Bc, rad_Bc)

    deallocate (mu_Shdo, phi_Shdo, wgt_Shdo, Nphi0_Shdo)
    deallocate (radI_Shdo, iradI_Shdo)
    deallocate (srcJ_Shdo, isrcJ_Shdo)
    if (isAccel) deallocate (dsrc_Shdo, idsrc_Shdo)
    if (shdo_obj%isAllocated) then
        nullify(shdo_obj%Nphi0, &
                shdo_obj%iradI, &
                shdo_obj%isrcJ, &
                shdo_obj%radI, &
                shdo_obj%srcJ, &
                shdo_obj%mu, &
                shdo_obj%phi, &
                shdo_obj%wgt)
    end if

    deallocate (x_Grid, y_Grid, z_Grid)
    deallocate (temp, planck)
    deallocate (ext, ssa, ipha)
    deallocate (leg)
    deallocate (ptr8Grid_Cell, ptr6Cell_Cell, ptr2Tree_Cell)
    deallocate (flag_Cell)
    deallocate (xyz_Grid)
    deallocate (fluxTot, fluxDir)
    deallocate (sumRad_Out, sumCam_Out)
    deallocate (work, work1, work2)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    call shd_mpi__end(Npts, xyz_Grid, rtm_obj%runName, prpAll_obj, prp_obj)

end program shdom
