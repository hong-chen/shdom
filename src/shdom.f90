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
                             doNmlInp_, doHighOrder_, &
                             isOpenX_, isOpenYx_, isMpiX_, isMpiY_, isMpi_, &
                             is1DX_, is1DYx_, is1D_, is3Ddir1Ddif_, &
                             is1D1_, isPol_, isPol4_, &
                             isSolar_, isThermal_
    use globals_type,  only: shd_type_prp, &
                             shd_type_sfc, &
                             shd_type_rtm, &
                             shd_type_grid

    use globals_obj,   only: shdo_obj_, iradI_Shdo_, radI_Shdo_, isrcJ_Shdo_, srcJ_Shdo_, &
                             src_obj_, &
                             prp_obj_, temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, ipha_Prp_, &
                             sfc_obj_, brdf_Sfc_, brdf_Rtm_, &
                             ckd_obj_, solar_Ckd_, wgt_Ckd_, kabs_Ckd_, &
                             sen_obj_, view_Sen_

    use shd_init, only: shd_init__makeGrid, &
                        shd_init__makeCell

    use shd_src,  only: shd_src__computeSrcOneDir

    use shd_shdo, only: shd_shdo__computeSH

    use shd_solv, only: shd_solv__solveRte

    use shd_rad,  only: shd_rad__computeRad, &
                        shd_rad__computeCamRad, &
                        shd_rad__makeDirect, &
                        shd_rad__computeFluxDiv

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

    real(kind=R__), allocatable, target :: x_Grid(:), y_Grid(:), z_Grid(:)
    real(kind=R__), allocatable, target :: x_Grid_all(:), y_Grid_all(:)

    real(kind=R__),  allocatable, target :: temp_Rtm(:), &
                                          planck_Rtm(:), &
                                             ext_Rtm(:), &
                                             ssa_Rtm(:), &
                                             leg_Rtm(:)
    integer(kind=2), allocatable, target :: ipha_Rtm(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! grid_obj
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_grid) :: grid_obj

    real(kind=R__), allocatable, target :: xyz_Grid(:, :)

    integer(kind=2),   allocatable, target :: flag_Cell(:)
    integer(kind=I__), allocatable, target :: i8Grid_Cell(:, :), &
                                              i6Cell_Cell(:, :), &
                                              i2Tree_Cell(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    real(kind=R__),    allocatable :: work(:)
    integer(kind=I__), allocatable :: work1(:)
    real(kind=R__),    allocatable :: work2(:)

    real(kind=R__), allocatable :: fluxTot(:, :), fluxDir(:)
    real(kind=R__), allocatable :: fluxDirTmp(:)

    real(kind=R__), allocatable :: fluxTot_all(:, :, :, :), fluxDiv_all(:, :, :)
    real(kind=R__), allocatable :: SH_all(:, :, :, :)

    ! shdo_obj_
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! SAV
    !╭────────────────────────────────────────────────────────────────────────────╮!
    character(len=charMaxLen_) :: fileInp_sav, fileOut_sav
    logical :: doRadInit, doGridInit
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! prp_obj_
    !╭────────────────────────────────────────────────────────────────────────────╮!
    type(shd_type_prp) :: prpAll_obj
    real(kind=R__),  allocatable :: temp_PrpAll(:), &
                                     ext_PrpAll(:), &
                                     ssa_PrpAll(:), &
                                     leg_PrpAll(:, :, :)
    integer(kind=2), allocatable :: ipha_PrpAll(:)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! sfc_obj_
    !╭────────────────────────────────────────────────────────────────────────────╮!
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! CKD
    !╭────────────────────────────────────────────────────────────────────────────╮!
    integer(kind=I__) :: ig_ckd
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! BC (borrowed from Boundary Condition)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    integer(kind=I__) :: NptsTop_bc, NptsBot_bc
    integer(kind=I__) :: maxNpts_bc, maxNrad_bc
    integer(kind=I__), allocatable :: ptr_Bc(:)
    real(kind=R__), allocatable    :: rad_Bc(:, :)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    real(kind=R__) :: NmaxNpts
    real(kind=R__) :: memUsed, memNeed, scaleFactor

    integer(kind=I__) :: NbigArr

    real(kind=R__) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real(kind=R__) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    real(kind=R__) :: MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT

    integer(kind=I__) :: oldNpts
    integer(kind=I__) :: maxIter, iter, totIter, iz, iout, Nout
    integer(kind=I__) :: Nrad, irad, Nx_out, Ny_out, Nlos, ilos, Nscan, Npix, Nsh_out
    integer(kind=I__) :: Ncell_tot, Npts_tot, Nsh_tot

    real(kind=R__) :: solAcc, solCrit, splitAcc, shAcc

    integer(kind=I__) :: Nproc
    logical :: isMasterProc

    integer(kind=I__) :: countStart, countRate, countMax
    character(len=charMaxLen_) :: message

    logical :: doCkd, doDeltaM, doAccel


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
        if (doNmlInp_) then
            call shd_io__nmlInp(rtmAll_obj%runName, prpAll_obj%fname, sfc_obj_%fname, ckd_obj_%fname,&
                fileInp_sav, fileOut_sav, sen_obj_%fname,&
                shdo_obj_%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtmAll_obj%Nz, shdo_obj_%Nmu, shdo_obj_%Nphi, rtmAll_obj%flag_bc,&
                rtmAll_obj%flag_ip, doCkd, doDeltaM, rtmAll_obj%type,&
                src_obj_%type, src_obj_%inpSolar, src_obj_%avgMu, src_obj_%avgAz, src_obj_%skyRad,&
                sfc_obj_%avgTemp, sfc_obj_%avgAlb, src_obj_%unit, src_obj_%waveNum, src_obj_%avgWaveLen,&
                doAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        else
            call shd_io__userInp(rtmAll_obj%runName, prpAll_obj%fname, sfc_obj_%fname, ckd_obj_%fname,&
                fileInp_sav, fileOut_sav, sen_obj_%fname,&
                shdo_obj_%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtmAll_obj%Nz, shdo_obj_%Nmu, shdo_obj_%Nphi, rtmAll_obj%flag_bc,&
                rtmAll_obj%flag_ip, doCkd, doDeltaM, rtmAll_obj%type,&
                src_obj_%type, src_obj_%inpSolar, src_obj_%avgMu, src_obj_%avgAz, src_obj_%skyRad,&
                sfc_obj_%avgTemp, sfc_obj_%avgAlb, src_obj_%unit, src_obj_%waveNum, src_obj_%avgWaveLen,&
                doAccel, solAcc, maxIter, splitAcc, shAcc,&
                maxNset_out, maxNpar_out, Nout, type_Out, parms_Out, file_Out,&
                fileOut_nc, MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
                NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, verbose_)
        end if
    end if

    call shd_mpi__bcastInp(rtmAll_obj%runName,&
        &fileInp_sav, fileOut_sav,&
        &shdo_obj_%Nstokes, rtmAll_obj%Nx, rtmAll_obj%Ny, rtmAll_obj%Nz, shdo_obj_%Nmu, shdo_obj_%Nphi, rtmAll_obj%flag_bc,&
        &rtmAll_obj%flag_ip, doCkd, doDeltaM, rtmAll_obj%type,&
        &src_obj_%type, src_obj_%inpSolar, src_obj_%avgMu, src_obj_%avgAz, src_obj_%skyRad,&
        &src_obj_%unit, src_obj_%waveNum, src_obj_%avgWaveLen,&
        &doAccel, solAcc, maxIter, splitAcc, shAcc,&
        &maxNset_out, maxNpar_out, Nout, type_Out, parms_Out,&
        &MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, prpAll_obj, sfc_obj_, ckd_obj_, sen_obj_)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Memory used for initialization (rough estimate, 1 WORD = 4 bytes)
    memUsed = 2000.0 + maxNset_out * (maxNpar_out + 0.25 + (charMaxLen_/4.0))

    ! Post-process parameters from input
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! set RT solving accuracy (1.0e-4 to 1.0e-8)
    accSolv_ = max(REPS_, min(solAcc*0.1, 1.0e-4))

    ! The number of Stokes parameters (and components of the real
    !   generalized spherical harmonics vectors) can be 1, 3, or 4.
    shdo_obj_%Nstokes = max(1, min(4, shdo_obj_%Nstokes))
    if (shdo_obj_%Nstokes .eq. 2) then
        shdo_obj_%Nstokes = 3
        if (verbose_ .ge. 4) write (6, '(A)') 'Warning [shdom]: shdo_obj_%Nstokes cannot be 2, reset shdo_obj_%Nstokes to 3.'
    end if

    ! prpAll_obj%NlegSet is the number of phase function elements used (1 or 6)
    !   Also used for the number of unique elements of the generalized
    !   spherical harmonics matrices.
    if (shdo_obj_%Nstokes .eq. 1) then
        prpAll_obj%NlegSet = 1
        shdo_obj_%NlegSet = 1
        isPol_ = .false.
    else
        prpAll_obj%NlegSet = 6
        shdo_obj_%NlegSet = 6
        isPol_ = .true.
    end if

    if (shdo_obj_%Nstokes .eq. 4) then
        isPol4_ = .true.
    else
        isPol4_ = .false.
    end if

    ! Make shdo_obj_%maxNl and shdo_obj_%maxNm from shdo_obj_%Nmu and shdo_obj_%Nphi (for discrete ordinate)
    !   shdo_obj_%Nmu : number of discrete ordinate zenith angles in both hemispheres
    !   shdo_obj_%Nphi: number of discrete ordinate azimuth angles in 2*pi
    !   shdo_obj_%maxNphi0: The maximum number of azimuth angles actually used
    shdo_obj_%Nmu  = max(2, 2*int((shdo_obj_%Nmu + 1)/2)) ! force shdo_obj_%Nmu an even number
    shdo_obj_%Nphi = max(1, shdo_obj_%Nphi)
    shdo_obj_%maxNphi0 = shdo_obj_%Nphi

    ! shdo_obj_%maxNl: The maximum L of the spherical harmonic truncation
    ! shdo_obj_%maxNm: The maximum M of the spherical harmonic truncation
    shdo_obj_%maxNl = shdo_obj_%Nmu - 1
    shdo_obj_%maxNm = max(0, int(shdo_obj_%Nphi/2) - 1)

    ! shdo_obj_%Nlm: The number of spherical harmonic terms
    shdo_obj_%Nlm = (2*shdo_obj_%maxNm + 1)*(shdo_obj_%maxNl + 1) - shdo_obj_%maxNm*(shdo_obj_%maxNm + 1)

    ! prpAll_obj%Nleg: The number of Legendre series phase function terms used in the solution procedure
    prpAll_obj%Nleg = shdo_obj_%maxNl
    if (doDeltaM) prpAll_obj%Nleg = shdo_obj_%maxNl + 1

    if (src_obj_%type .eq. 'B') then
        isSolar_ = .true.
        isThermal_ = .true.
    else
        if (src_obj_%type .eq. 'S') then
            isSolar_ = .true.
            isThermal_ = .false.
        else if (src_obj_%type .eq. 'T') then
            isSolar_ = .false.
            isThermal_ = .true.
        else
            stop 'Error [shdom]: Unknown src_obj_%type.'
        end if
    end if

    ! for thermal, reset solar azimuth angle
    if (src_obj_%type .eq. 'T') src_obj_%avgAz = 0.0d0

    ! angle grids for DO
    allocate(shdo_obj_%Nphi0(shdo_obj_%Nmu), &
             shdo_obj_%mu(shdo_obj_%Nmu), &
             shdo_obj_%phi(shdo_obj_%Nmu, shdo_obj_%Nphi), &
             shdo_obj_%wgt(shdo_obj_%Nmu, shdo_obj_%Nphi))
    shdo_obj_%isAllocated = .true.

    memUsed = memUsed + 2.0*shdo_obj_%Nmu*(1.0 + shdo_obj_%Nphi)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read prpAll_obj%fname (property file) to get data size info, broadcast data if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        ! Find out the size of the input property file grid
        call shd_io__readPrpSize(prpAll_obj)
    end if

    call shd_mpi__bcastPrpSize(prpAll_obj)
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sfc_obj_%fname (surface file) to get general info, broadcast sfc_obj_ if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (sfc_obj_%fname(1:4) .eq. 'NONE') then
        sfc_obj_%fname = 'NONE'
        sfc_obj_%type  = 'FL'
        sfc_obj_%Nx = 1
        sfc_obj_%Ny = 1
        sfc_obj_%dx = prpAll_obj%dx
        sfc_obj_%dy = prpAll_obj%dy
        sfc_obj_%maxNpts = 1
        sfc_obj_%maxNpar = 2
        sfc_obj_%isUniform    = .true.
        sfc_obj_%isLambertian = .true.
    else
        if (isMasterProc) then
            ! find out the size of the surface file grid
            call shd_io__readSfcSize(sfc_obj_)
        end if
        ! Broadcast the surface object
        call shd_mpi__bcastSfcSize(sfc_obj_)
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

    ! Overwrite prpAll_obj%dx if only 1 grid in X direction
    if (prpAll_obj%Nx .eq. 1) then
        prpAll_obj%dx = sfc_obj_%Nx * sfc_obj_%dx
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: prpAll_obj%Nx=1, reset prpAll_obj%dx to match surface domain.'
    end if

    ! overwrite prpAll_obj%dy if only 1 grid in Y direction
    if (prpAll_obj%Ny .eq. 1) then
        prpAll_obj%dy = sfc_obj_%Ny * sfc_obj_%dy
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: prpAll_obj%Ny=1, reset prpAll_obj%dy to match surface domain.'
    end if

    prpAll_obj%xDomain = prpAll_obj%Nx * prpAll_obj%dx
    prpAll_obj%yDomain = prpAll_obj%Ny * prpAll_obj%dy

    ! Choose the option to compute higher order radiance SH terms
    if  ((prpAll_obj%Nx .eq. 1) .and. (prpAll_obj%Ny .eq. 1) .and.&
        &(sfc_obj_%Nx .le. 32) .and. (sfc_obj_%Ny .le. 32)) then
        doHighOrder_ = .true.
        if (verbose_ .ge. 4) write (6, '(A)') &
            'Warning [shdom]: 1D calculation, set doHighOrder_ to .TRUE.'
    else
        doHighOrder_ = .false.
    end if

    Nsh_out = 4
    if (doHighOrder_) Nsh_out = 5
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! If MPI, divide domain into mini-domains
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Map the full SHDOM domain onto the multiple processors:
    !   change NPX,NPY and Nx,Ny to the subdomain size, and output
    !   the X and Y starting position of this subdomain. Set bits
    !   2 and 3 in rtmAll_obj%flag_bc to indicate multiple processors in X and Y.
    !   Also open the slave processor log files.
    rtm_obj = rtmAll_obj; prp_obj_ = prpAll_obj
    call shd_mpi__map(rtmAll_obj, prpAll_obj, rtm_obj, prp_obj_)

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

    grid_obj%Nz = rtm_obj%Nz

    ! Number of base points
    grid_obj%Npts0 = grid_obj%Nx*grid_obj%Ny*rtm_obj%Nz

    ! Calculate the number of base grid cells depending on the boundary condition flag
    grid_obj%Ncell0 = (rtm_obj%Nz - 1)*(rtm_obj%Nx + ibits(rtm_obj%flag_bc, 0, 1) - ibits(rtm_obj%flag_bc, 2, 1))&
                                      *(rtm_obj%Ny + ibits(rtm_obj%flag_bc, 1, 1) - ibits(rtm_obj%flag_bc, 3, 1))
    !╰────────────────────────────────────────────────────────────────────────────╯!

    ! Read prp_obj_%fname and distribute property object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Determine the size of the property grid arrays so can allocate
    prp_obj_%maxNpts = prp_obj_%Nx * prp_obj_%Ny * prp_obj_%Nz
    if (prp_obj_%Npha .eq. 0) &
        prp_obj_%maxNptsLeg = prp_obj_%maxNpts * (prp_obj_%maxNleg + 1) * prp_obj_%NlegSet

    allocate(prp_obj_%z(prp_obj_%Nz))
    prp_obj_%isAllocated = .true.

    allocate(temp_Prp_(prp_obj_%maxNpts),  &
              ext_Prp_(prp_obj_%maxNpts),  &
              ssa_Prp_(prp_obj_%maxNpts),  &
              leg_Prp_(1:prp_obj_%NlegSet, &
                       0:prp_obj_%maxNleg, &
                       1:prp_obj_%maxNptsLeg/(prp_obj_%maxNleg+1)/prp_obj_%NlegSet), &
             ipha_Prp_(prp_obj_%maxNpts))

    ! note: in addition to xxx_Prp_, here memUsed accounts for extDir_Prp
    !       that will be used in shd_rad__computeDirect (shdRad.f90)
    memUsed = memUsed + prp_obj_%Nz + 6.0*prp_obj_%maxNpts + prp_obj_%maxNptsLeg

    if (isMpi_) then

        ! Read the properties of the medium to the master processor and
        !   distribute to all the processors
        if (isMasterProc) then

            allocate(prpAll_obj%z(prpAll_obj%Nz))
            prpAll_obj%isAllocated = .true.

            allocate(temp_PrpAll(prpAll_obj%maxNpts),  &
                      ext_PrpAll(prpAll_obj%maxNpts),  &
                      ssa_PrpAll(prpAll_obj%maxNpts),  &
                      leg_PrpAll(prpAll_obj%NlegSet,   &
                                 0:prpAll_obj%maxNleg, &
                                 prpAll_obj%maxNptsLeg/(prpAll_obj%maxNleg+1)/prpAll_obj%NlegSet), &
                     ipha_PrpAll(prpAll_obj%maxNpts))

            call shd_io__readPrp(doDeltaM, &
                                 prpAll_obj, temp_PrpAll, ext_PrpAll, ssa_PrpAll, leg_PrpAll, ipha_PrpAll)

            if (prpAll_obj%Npha .eq. 0) then
                stop 'Error [shdom]: Standard property file (with phase function Legendre series at every point)&
                    & is not allowed with multiple processors.'
            end if

            ! Check values of medium properties
            call shd_io__checkPrp(prpAll_obj, temp_PrpAll, ext_PrpAll, ssa_PrpAll, leg_PrpAll, ipha_PrpAll)

        end if

        ! distribute file grid data to working grids (different CPU processes different grids)
        call shd_mpi__bcastPrp(prpAll_obj,  temp_PrpAll, ext_PrpAll, ssa_PrpAll, leg_PrpAll, ipha_PrpAll, &
                                  prp_obj_, temp_Prp_,   ext_Prp_,   ssa_Prp_,    leg_Prp_,    ipha_Prp_)

        if (isMasterProc) then
            deallocate(temp_PrpAll, &
                        ext_PrpAll,  &
                        ssa_PrpAll,  &
                        leg_PrpAll,  &
                       ipha_PrpAll)
            if (prpAll_obj%isAllocated) then
                deallocate(prpAll_obj%z)
                prpAll_obj%isAllocated = .false.
            end if
        end if

    else

        ! Read the properties of the medium to the single processor
        call shd_io__readPrp(doDeltaM, &
                             prp_obj_, temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, ipha_Prp_)

        ! Check values of medium properties
        call shd_io__checkPrp(prp_obj_, temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, ipha_Prp_)

    end if

    if (isMasterProc) then
        ! Check input parameters
        call shd_io__checkInp(doDeltaM, splitAcc, shAcc,&
            &rtm_obj, shdo_obj_, src_obj_, prp_obj_, sfc_obj_)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sfc_obj_%fname and broadcast the surface object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in the surface properties
    if (sfc_obj_%fname(1:4) .ne. 'NONE') then

        allocate(brdf_Sfc_(sfc_obj_%maxNpar, sfc_obj_%maxNpts))

        memUsed = memUsed + sfc_obj_%maxNpar * sfc_obj_%maxNpts

        if (isMasterProc) then
            call shd_io__readSfc(sfc_obj_, brdf_Sfc_)
        end if

        ! Broadcast the surface object
        call shd_mpi__bcastSfc(sfc_obj_, brdf_Sfc_)

    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read ckd_obj_%fname and broadcast correlated-k object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! If doing a k-distribution then get the band info from the CKD file
    if (doCkd) then

        if (isMasterProc) then
            call shd_io__readCkdSize(src_obj_%waveNum, &
                                     ckd_obj_)
        end if

        ! Broadcast the k-distribution object
        call shd_mpi__bcastCkdSize(ckd_obj_)

        allocate(ckd_obj_%z(ckd_obj_%Nz))
        ckd_obj_%isAllocated = .true.

        allocate(solar_Ckd_(ckd_obj_%Ng), &
                   wgt_Ckd_(ckd_obj_%Ng),   &
                  kabs_Ckd_(ckd_obj_%Nz, ckd_obj_%Ng))

        memUsed = memUsed + ckd_obj_%Nz*(2.0+ckd_obj_%Ng) + 2.0*ckd_obj_%Ng

        if (isMasterProc) then
            ckd_obj_%maxNz = ckd_obj_%Nz
            ckd_obj_%maxNg = ckd_obj_%Ng
            call shd_io__readCkd(src_obj_%waveNum, &
                                 ckd_obj_, solar_Ckd_, wgt_Ckd_, kabs_Ckd_)
        end if

        call shd_mpi__bcastCkd(ckd_obj_, solar_Ckd_, wgt_Ckd_, kabs_Ckd_)

    else
        ckd_obj_%Ng = 1
        ckd_obj_%Nz = 0

        allocate(solar_Ckd_(ckd_obj_%Ng), wgt_Ckd_(ckd_obj_%Ng))

        memUsed = memUsed + 2.0*ckd_obj_%Ng

        solar_Ckd_(1) = 1.0
        wgt_Ckd_(1) = 1.0

        isBaseOut = .false.
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read sen_obj_%fname and broadcast sensor object if MPI
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Read in user-specified sensor properties if given
    !   this is a binary angle input implemented for 'R' and 'V3'
    if (sen_obj_%fname(1:4) .ne. 'NONE') then

        if (isMasterProc) then
            call shd_io__readSenSize(sen_obj_)
        end if

        call shd_mpi__bcastSenSize(sen_obj_)

        allocate(view_Sen_(sen_obj_%Npar, sen_obj_%Npts))

        memUsed = memUsed + sen_obj_%Npar*sen_obj_%Npts

        if (isMasterProc) then
            call shd_io__readSen(sen_obj_, view_Sen_)
        end if

        call shd_mpi__bcastSen(sen_obj_, view_Sen_)
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
                    Nx_out = max(1, nint(prp_obj_%xDomain/parms_Out(2, iout)))
                else
                    Nx_out = max(1, nint((prp_obj_%xDomain + 2*parms_Out(4, iout))/parms_Out(2, iout)))
                end if
            end if
            if (parms_Out(3, iout) .le. 0.0) then
                Ny_out = 1
            else
                if (isMpiY_) then
                    Ny_out = max(1, nint(prp_obj_%yDomain/parms_Out(3, iout)))
                else
                    Ny_out = max(1, nint((prp_obj_%yDomain+ 2*parms_Out(5, iout))/parms_Out(3, iout)))
                end if
            end if
            Nrad = Nrad + shdo_obj_%Nstokes*Nx_out*Ny_out*nint(parms_Out(6, iout))
        else if (type_Out(iout) .eq. 'V') then
            isCamOut = .true.
            if (nint(parms_Out(1, iout)) .eq. 1) then
                Nlos = Nlos + shdo_obj_%Nstokes*nint(parms_Out(10, iout))*nint(parms_Out(11, iout))
            else if (nint(parms_Out(1, iout)) .eq. 2) then
                Nscan = 1 + sqrt((parms_Out(4, iout) - parms_Out(7, iout))**2&
                    + (parms_Out(5, iout) - parms_Out(8, iout))**2&
                    + (parms_Out(6, iout) - parms_Out(9, iout))**2)/parms_Out(10, iout)
                Npix = 1 + abs(parms_Out(12, iout) - parms_Out(11, iout))/parms_Out(13, iout)
                Nlos = Nlos + shdo_obj_%Nstokes*Nscan*Npix
            ! user-specified sensor
            else if (nint(parms_Out(1, iout)) .eq. 3) then
                Nlos = Nlos + shdo_obj_%Nstokes*nint(parms_Out(2, iout))
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
    if (.not. doAccel) NbigArr = 2

    ! NmaxNpts
    !╭──────────────────────────────────────────────────────────────╮!
    NmaxNpts = 0.0

    NmaxNpts = NmaxNpts + 4.5 ! temp_Rtm, planck_Rtm, ext_Rtm, ssa_Rtm, ipha_Rtm[0.5, integer(kind=2)]

    ! x, y, z grid
    memUsed = memUsed + (rtm_obj%Nx+1) + (rtm_obj%Ny+1) + (rtm_obj%Nz)

    if (prp_obj_%Npha .gt. 0) then
        memUsed = memUsed + prp_obj_%Npha * prp_obj_%NlegSet * (prp_obj_%maxNleg + 1) ! offset phase table
    else
        NmaxNpts = NmaxNpts + prp_obj_%NlegSet*(prp_obj_%maxNleg + 1)
    end if

    if (sfc_obj_%fname(1:4) .ne. 'NONE') then
        if (sfc_obj_%type(2:2) .eq. 'L') then
            NmaxNpts = 3.0/real(rtm_obj%Nz) * (sfc_obj_%maxNpar + 2.0 + shdo_obj_%Nstokes*2.0)
        else
            NmaxNpts = 3.0/real(rtm_obj%Nz) * (sfc_obj_%maxNpar + 2.0 + shdo_obj_%Nstokes*(2.0+ shdo_obj_%Nmu*shdo_obj_%maxNphi0/2))
        end if
    end if

    if (isFluxOut)  NmaxNpts = NmaxNpts + 3.0
    if (isHROut)    NmaxNpts = NmaxNpts + 1.0
    if (isSHOut)    NmaxNpts = NmaxNpts + Nsh_out
    if (isSrcOut)   NmaxNpts = NmaxNpts + shdo_obj_%Nstokes

    NmaxNpts = NmaxNpts + 16.5*CELL_TO_POINT_RATIO ! i8Grid_Cell[8.0], i6Cell_Cell[6.0], i2Tree_Cell[2.0], flag_Cell[0.5]

    ! note: in addition to xxx_Shdo_, here memUsed accounts for isrcD_Shdo and srcD_Shdo
    !       that will be used in shd_solv__solveRte(shdSolv.f90)
    NmaxNpts = NmaxNpts + 6.0   ! xyz_Grid[3.0], iradI_Shdo_[1.0], isrcJ_Shdo_[1.0], isrcD_Shdo[1.0]
    memUsed = memUsed + 4.0 ! offset iradI_Shdo_[+2], isrcJ_Shdo_[+1], isrcD_Shdo[+1]

    if (doAccel) then
        NmaxNpts = NmaxNpts + shdo_obj_%Nstokes*(1.0+3.0*NUM_SH_TERM_FACTOR*shdo_obj_%Nlm) ! radiance, source, delta source
    else
        NmaxNpts = NmaxNpts + shdo_obj_%Nstokes*(1.0+2.0*NUM_SH_TERM_FACTOR*shdo_obj_%Nlm) ! radiance, source
    end if

    NmaxNpts = NmaxNpts + 3.0 ! fluxTot[2.0], fluxDir[1.0]

    memUsed = memUsed + (Nrad+Nlos) ! offset sumRad_Out[Nrad] and sumCam_Out[Nlos]

    NmaxNpts = NmaxNpts + (shdo_obj_%Nstokes*(max(shdo_obj_%maxNphi0, 4)+1.0) + 8.0) !work, work1, work2

    memUsed = memUsed + Nrad + Nlos ! work2 (over-count but memory safe)
    !╰──────────────────────────────────────────────────────────────╯!

    ! memUsed is now more accurate than before
    ! before:
    !   memUsed = shdo_obj_%Nmu*(2 + 2*shdo_obj_%Nphi + 2*shdo_obj_%Nlm + 2*33*32)&            ! SHDO allocation
    !       + 4.5*prp_obj_%maxNpts + prp_obj_%maxNptsLeg + prp_obj_%Npha*prp_obj_%NlegSet*(prp_obj_%maxNleg + 1) ! phase allocation
    ! now: update memUsed after every allocation step
    !
    ! NmaxNpts is now more accurate than before
    ! before:
    !   NmaxNpts = 28.0 + 16.5*CELL_TO_POINT_RATIO&
    !       + shdo_obj_%maxNphi0*shdo_obj_%Nstokes + NUM_SH_TERM_FACTOR*shdo_obj_%Nstokes*shdo_obj_%Nlm*NbigArr
    ! now: update NmaxNpts from every allocation step

    memNeed = ADAPT_GRID_FACTOR*grid_obj%Npts0*NmaxNpts

    if (doHighOrder_) then
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
    if (memNeed .gt. huge(shdo_obj_%maxNptsSH)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Number of words of memory exceeds max integer size:', memNeed
        stop
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Retrieve max numbers for array allocation after memory check
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! max number of internal grid points
    grid_obj%maxNpts = nint(ADAPT_GRID_FACTOR*grid_obj%Npts0)

    ! max number of grid cells, between 1 and 2 times grid points
    grid_obj%maxNcell = nint(CELL_TO_POINT_RATIO*grid_obj%maxNpts)

    ! max words for source function and radiance arrays
    shdo_obj_%maxNptsSH = nint(NUM_SH_TERM_FACTOR*shdo_obj_%Nlm*grid_obj%maxNpts)

    ! max words for temporary discrete ordinate azimuthal array
    shdo_obj_%maxNptsDO = grid_obj%maxNpts*max(shdo_obj_%maxNphi0, 4)

    if (4.0*(shdo_obj_%maxNptsSH + grid_obj%maxNpts)*shdo_obj_%Nstokes .gt. huge(shdo_obj_%maxNptsSH)) then
        write (6, '(A,I12)') 'Error [shdom]: Size of big SH arrays (shdo_obj_%maxNptsSH) probably exceeds max integer number of bytes:', shdo_obj_%maxNptsSH
        stop
    end if

    if (4.0*8.0*grid_obj%maxNcell .gt. huge(shdo_obj_%maxNptsSH)) then
        write (6, '(A,F8.2)') 'Error [shdom]: Size of i8Grid_Cell array (8*grid_obj%maxNcell) probably exceeds max integer number of bytes:', 8*grid_obj%maxNcell
        stop
    end if

    maxNpts_bc = grid_obj%maxNpts*3/rtm_obj%Nz
    if (sfc_obj_%isLambertian) then
        maxNrad_bc = 2 * maxNpts_bc
    else
        maxNrad_bc = (2 + shdo_obj_%Nmu*shdo_obj_%maxNphi0/2)*maxNpts_bc
    end if

    ! max number of legendre coefficients
    if (prp_obj_%Npha .gt. 0) then
        rtm_obj%maxNptsLeg = prp_obj_%Npha * prp_obj_%NlegSet * (prp_obj_%maxNleg + 1)
    else
        rtm_obj%maxNptsLeg = grid_obj%maxNpts*prp_obj_%NlegSet*(prp_obj_%maxNleg + 1)
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Allocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! rtm
    !╭────────────────────────────────────────────────╮!
    ! x, y, z grids
    allocate(x_Grid(rtm_obj%Nx+1), y_Grid(rtm_obj%Ny+1), z_Grid(grid_obj%Nz))

    ! temperature and planck radiance (maybe we don't need this for shortwave)
    allocate(temp_Rtm(grid_obj%maxNpts), planck_Rtm(grid_obj%maxNpts))

    ! extinction, single scattering albedo, and phase function indices
    allocate(ext_Rtm(grid_obj%maxNpts), ssa_Rtm(grid_obj%maxNpts), ipha_Rtm(grid_obj%maxNpts))

    ! legendre coefficients
    allocate(leg_Rtm(rtm_obj%maxNptsLeg))

    rtm_obj%x => x_Grid
    rtm_obj%y => y_Grid
    rtm_obj%z => z_Grid

    rtm_obj%temp   => temp_Rtm
    rtm_obj%planck => planck_Rtm
    rtm_obj%ext    => ext_Rtm
    rtm_obj%ssa    => ssa_Rtm
    rtm_obj%ipha   => ipha_Rtm

    rtm_obj%leg(1:prp_obj_%NlegSet, &
                0:prp_obj_%maxNleg, &
                1:rtm_obj%maxNptsLeg/(prp_obj_%maxNleg+1)/prp_obj_%NlegSet) => leg_Rtm

    rtm_obj%isAllocated = .true.
    !╰────────────────────────────────────────────────╯!

    ! grid
    !╭────────────────────────────────────────────────╮!
    ! x,y,z coordinate values of grid point locations
    allocate(xyz_Grid(3, grid_obj%maxNpts))

    ! a collection of binary flags for the cells
    allocate(flag_Cell(grid_obj%maxNcell))

    ! pointer to 8 grid points belonging to a cell
    ! pointer to 6 neighboring cells of a cell
    ! pointer to parent and children cells in tree structure
    allocate(i8Grid_Cell(8, grid_obj%maxNcell), i6Cell_Cell(6, grid_obj%maxNcell), i2Tree_Cell(2, grid_obj%maxNcell))

    grid_obj%xyz    => xyz_Grid
    grid_obj%flag   => flag_Cell
    grid_obj%i8Grid => i8Grid_Cell
    grid_obj%i6Cell => i6Cell_Cell
    grid_obj%i2Tree => i2Tree_Cell
    grid_obj%isAllocated = .true.
    !╰────────────────────────────────────────────────╯!

    ! boundary points pointers and radiances
    !   ptr_Bc(:, 1): pointers to grid points that are on top
    !   ptr_Bc(:, 2): pointers to grid points that are on bottom
    !   rad_Bc(:[shdo_obj_%Nstokes], 1:NptsTop_bc): isotropic top radiances
    !   rad_Bc(:[shdo_obj_%Nstokes], NptsTop_bc+1:NptsTop_bc+NptsBot_bc): upwelling bottom radiances
    !   rad_Bc(:[shdo_obj_%Nstokes], NptsTop_bc+NptsBot_bc+1:NptsTop_bc+NptsBot_bc+NptsBot_bc*shdo_obj_%Nang/2):
    !         downwelling bottom radiances at discrete ordinates
    if (sfc_obj_%fname(1:4) .ne. 'NONE') then
        allocate(brdf_Rtm_(sfc_obj_%maxNpar, maxNpts_bc))
    end if

    allocate(ptr_Bc(2*maxNpts_bc), rad_Bc(shdo_obj_%Nstokes, maxNrad_bc))

    ! output
    if (isFluxOut) allocate(sumFluxTot_Out(2, grid_obj%maxNpts), &
                            sumFluxDir_Out(grid_obj%maxNpts))
    if (isHROut)   allocate(sumFluxDiv_Out(grid_obj%maxNpts))
    if (isSHOut)   allocate(sumSH_Out(grid_obj%maxNpts*Nsh_out))
    if (isSrcOut)  allocate(src1_Out(grid_obj%maxNpts*shdo_obj_%Nstokes))

    ! shdo
    !╭────────────────────────────────────────────────╮!
    ! SH related array
    allocate(iradI_Shdo_(grid_obj%maxNpts+2), &
              radI_Shdo_(shdo_obj_%Nstokes, shdo_obj_%maxNptsSH + grid_obj%maxNpts))

    ! radiance and source
    allocate(isrcJ_Shdo_(grid_obj%maxNpts+1), &
              srcJ_Shdo_(shdo_obj_%Nstokes, shdo_obj_%maxNptsSH))
    !╰────────────────────────────────────────────────╯!

    ! fluxes
    allocate(fluxTot(2, grid_obj%maxNpts), fluxDir(grid_obj%maxNpts))

    ! output radiances to be summed over g
    allocate(sumRad_Out(Nrad), sumCam_Out(Nlos))

    ! work: work array for path integration; radiance field, Nazimuth*grid_obj%Npts
    allocate(work(shdo_obj_%Nstokes*shdo_obj_%maxNptsDO))

    ! work1: work array for path integration; grid point sweeping order
    allocate(work1(8*grid_obj%maxNpts))

    ! work2: work array for path integration; single angle radiance
    allocate(work2(max(Nlos, Nrad, shdo_obj_%Nstokes*grid_obj%maxNpts)))
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Read in previous run from binary file if desired
    call shd_io__restoreState(fileInp_sav, &
        doRadInit, doGridInit, fluxTot, &
        rtm_obj, grid_OBJ, &
        shdo_obj_, iradI_Shdo_, radI_Shdo_, isrcJ_Shdo_, srcJ_Shdo_)


    ! Check if we can use pre-calculated "state", if not, generate new grids
    !   doGridInit is from shd_io__restoreState call
    if (doGridInit) then
        ! Make the internal base grid lines
        ! note: we might not need <rtm_obj> after all
        call shd_init__makeGrid(rtm_OBJ, prp_OBJ_)
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
    src_obj_%avgSolar = 0.0d0

    ! Loop over the k-distribution g's from high to low absorption
    CKD_LOOP: do ig_ckd = ckd_obj_%Ng, 1, -1

        if (verbose_ .ge. 5) write (6, '(A,I2,A,I2,A)') 'Message [shdom]: Starting ', ckd_obj_%Ng-ig_ckd+1, '/', ckd_obj_%Ng, ' (ig/Ng) ...'

        ! solar flux for given g of CKD
        src_obj_%avgSolar0 = solar_Ckd_(ig_ckd)*src_obj_%inpSolar*abs(src_obj_%avgMu)

        ! weighted average solar flux
        src_obj_%avgSolar = src_obj_%avgSolar + solar_Ckd_(ig_ckd)*wgt_Ckd_(ig_ckd)*src_obj_%inpSolar*abs(src_obj_%avgMu)

        ! If doing k-distribution then get the gas absorption profile.
        ckd_obj_%ig = ig_ckd

        ! If this is the first time through or we always want a base grid,
        !   then make the base grid point positions and grid cell structure.
        if (isBaseOut .or. ((ig_ckd .eq. ckd_obj_%Ng) .and. doGridInit)) then
            ! note: we might not need <rtm_obj> after all
            call shd_init__makeCell(rtm_obj, grid_OBJ)
            grid_obj%Ncell0 = grid_obj%Ncell
        end if

        ! Find the monochromatic SHDOM solution: adaptive grid output
        !   in source (srcJ_Shdo_), radiance (radI_Shdo_), fluxes, direct Flux.
        call shd_solv__solveRte(&
            &doRadInit,&
            &doDeltaM, doHighOrder_,&
            &doAccel, solAcc, maxIter, solCrit, iter,&
            &splitAcc, shAcc,&
            &maxNpts_bc, maxNrad_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
            &work, work1, work2,&
            &fluxTot, fluxDir,&
            &rtm_obj, grid_obj)

        doRadInit = .false.
        totIter = totIter + iter

        ! Compute the output and put in arrays. First the radiance.
        if (isRadOut) then
            irad = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'R') then
                    if (isMpi_) then
                        call shd_mpi__computeRad(&
                            &doDeltaM,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &fluxDir, fluxTot,&
                            &work, work(1 + grid_obj%Npts*shdo_obj_%Nstokes),&
                            &parms_Out(1, iout), irad, work2,&
                            &rtm_obj, grid_obj)
                    else
                        call shd_rad__computeRad(&
                            &doDeltaM,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &fluxDir, fluxTot,&
                            &work(1), work(1 + grid_obj%Npts*shdo_obj_%Nstokes),&
                            &parms_Out(1, iout), irad, work2,&
                            &rtm_obj, grid_obj)
                    end if
                end if
            end do
            call shd_util__sumData(ckd_obj_%Ng-ig_ckd, wgt_Ckd_(ig_ckd), Nrad, work2, sumRad_Out)
        end if

        if (isCamOut) then
            ilos = 0
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'V') then
                    if (isMpi_) then
                        call shd_mpi__computeCamRad(&
                            &doDeltaM,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &fluxDir, fluxTot,&
                            &parms_Out(1, iout), ilos, work2,&
                            &rtm_obj, grid_obj)
                    else
                        call shd_rad__computeCamRad(&
                            &doDeltaM,&
                            &maxNpts_bc, NptsTop_bc, NptsBot_bc, ptr_Bc, rad_Bc,&
                            &fluxDir, fluxTot,&
                            &parms_Out(1, iout), ilos, work2,&
                            &rtm_obj, grid_obj)
                    end if
                end if
            end do
            call shd_util__sumData(ckd_obj_%Ng-ig_ckd, wgt_Ckd_(ig_ckd), Nlos, work2, sumCam_Out)
        end if

        if (isSrcOut) then
            ! Compute source function throughout grid for this one angle
            do iout = 1, Nout
                if (type_Out(iout) .eq. 'J') then
                    mu_out  = parms_Out(2, iout)
                    phi_out = parms_Out(3, iout)*acos(-1.0)/180.0
                    call shd_src__computeSrcOneDir(doDeltaM,&
                        &mu_out, phi_out,&
                        &fluxDir, src1_Out,&
                        &rtm_obj, grid_obj, isrcJ_Shdo_, srcJ_Shdo_)
                end if
            end do
        end if

        if (isBaseOut) then
            ! If outputting only the base grid then throw out the other points
            call shd_init__makeCell(rtm_obj, grid_OBJ)
        else
            ! For the output from the grid point arrays (Flux, Net flux, SH)
            !   we may need to interpolate the running k-distribution sum
            !   to the newly created grid points.
            call shd_util__interpData(oldNpts, isFluxOut, isHROut, isSHOut,&
                sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out, sumSH_Out,&
                grid_obj)
        end if

        if (isFluxOut) then

            ! unscale delta-m fluxes
            if ((src_obj_%type.ne. 'T') .and. doDeltaM) then

                allocate(fluxDirTmp(grid_obj%Npts))

                if (isMpi_) then
                    call shd_mpi__makeDirect(1, .false., shdo_obj_%maxNl,&
                        &src_obj_%avgSolar0, src_obj_%avgMu, src_obj_%avgAz,&
                        &fluxDirTmp,&
                        &rtm_obj, grid_obj)
                else
                    call shd_rad__makeDirect(grid_obj%Npts, .false., shdo_obj_%maxNl,&
                        &src_obj_%avgSolar0, src_obj_%avgMu, src_obj_%avgAz,&
                        &fluxDirTmp,&
                        &rtm_obj, grid_obj)
                end if

                fluxTot(1, 1:grid_obj%Npts) = fluxTot(1, 1:grid_obj%Npts) + fluxDir(1:grid_obj%Npts) - fluxDirTmp(1:grid_obj%Npts)
                fluxDir(1:grid_obj%Npts) = fluxDirTmp(1:grid_obj%Npts)

                deallocate(fluxDirTmp)

            end if

            call shd_util__sumData(ckd_obj_%Ng-ig_ckd, wgt_Ckd_(ig_ckd), 2*grid_obj%Npts, fluxTot, sumFluxTot_Out)
            call shd_util__sumData(ckd_obj_%Ng-ig_ckd, wgt_Ckd_(ig_ckd), grid_obj%Npts, fluxDir, sumFluxDir_Out)
        end if

        if (isHROut) then
            ! Compute the net flux divergence at every point
            call shd_rad__computeFluxDiv(shdo_obj_%Nstokes, grid_obj%Npts, iradI_Shdo_, src_obj_%type, src_obj_%avgMu,&
                &rtm_obj%ext, rtm_obj%ssa, rtm_obj%planck, fluxDir, radI_Shdo_, work)
            call shd_util__sumData(ckd_obj_%Ng-ig_ckd, wgt_Ckd_(ig_ckd), grid_obj%Npts, work, sumFluxDiv_Out)
        end if

        if (isSHOut) then
            call shd_shdo__computeSH(Nsh_out, shdo_obj_%Nstokes, grid_obj%Npts, src_obj_%type, src_obj_%avgMu, src_obj_%avgAz,&
                &fluxDir, iradI_Shdo_, shdo_obj_%maxNl, shdo_obj_%maxNm, radI_Shdo_, work)
            call shd_util__sumData(ckd_obj_%Ng-ig_ckd, wgt_Ckd_(ig_ckd), Nsh_out*grid_obj%Npts, work, sumSH_Out)
        end if

        oldNpts = grid_obj%Npts

    end do CKD_LOOP
    !╰────────────────────────────────────────────────────────────────────────────╯!


    message = 'Message [shdom]: Total time used for runs: '
    call shd_util__tickTock(message, countStart)
    if (isMasterProc) then
        if (verbose_ .ge. 5) write (6, '(A)') trim(adjustl(message))
    end if

    ! Output radiance and/or flux (and/or other radiative quantities)
    !╭────────────────────────────────────────────────────────────────────────────╮!
    shdo_obj_%Nsh = isrcJ_Shdo_(grid_obj%Npts + 1)
    irad = 0
    ilos = 0
    if (isMpi_) then

        if (isMasterProc) then
            if (isFluxOut) allocate(fluxTot_all(3, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))
            if (isHROut)   allocate(fluxDiv_all(rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))
            if (isSHOut)   allocate(SH_all(Nsh_out, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx))
            allocate(x_Grid_all(rtmAll_obj%Nx), y_Grid_all(rtmAll_obj%Ny))
            rtmAll_obj%x => x_Grid_all
            rtmAll_obj%y => y_Grid_all
            rtmAll_obj%z => z_Grid
            rtmAll_obj%isAllocated = .true.
        end if

        call shd_mpi__gatherOut(isFluxOut, isHROut, isSHOut,&
            &Nsh_out,&
            &sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out, sumSH_Out,&
            &fluxTot_all, fluxDiv_all, SH_all,&
            &shdo_obj_%Nsh, Ncell_tot, Npts_tot, Nsh_tot,&
            &rtmAll_obj, grid_obj)

        if (isMasterProc) then

            do iout = 1, Nout
                if (type_Out(iout) .ne. 'V' .and. file_Out(iout)(1:4) .ne. 'NONE') then
                    call shd_io__outputResultMpi(&
                        &Npts_tot, Ncell_tot, Nsh_tot,&
                        &fileInp_sav, fileOut_sav,&
                        &doDeltaM,&
                        &splitAcc, shAcc, solAcc, maxIter, totIter,&
                        &fluxTot_all, fluxDiv_all,&
                        &Nsh_out, SH_all, irad, sumRad_Out,&
                        &type_Out(iout), parms_Out(1, iout), file_Out(iout),&
                        &rtmAll_obj, grid_obj)
                end if
            end do

            if (isFluxOut) deallocate(fluxTot_all)
            if (isHROut)   deallocate(fluxDiv_all)
            if (isSHOut)   deallocate(SH_all)

            deallocate(x_Grid_all, y_Grid_all)
            if (rtmAll_obj%isAllocated) then
                nullify(rtmAll_obj%x, rtmAll_obj%y, rtmAll_obj%z)
                rtmAll_obj%isAllocated = .false.
            end if
        end if

    else

        do iout = 1, Nout
            if ((type_Out(iout) .ne. 'V') .and. (file_Out(iout)(1:4) .ne. 'NONE')) then
                call shd_io__outputResult(&
                    fileInp_sav, fileOut_sav,&
                    doDeltaM,&
                    splitAcc, shAcc, solAcc, maxIter, totIter,&
                    sumFluxTot_Out, sumFluxDir_Out, sumFluxDiv_Out,&
                    irad, sumRad_Out, Nsh_out, sumSH_Out, src1_Out,&
                    type_Out(iout), parms_Out(1, iout), file_Out(iout),&
                    rtm_obj, grid_obj)
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Output camera radiance if there are any
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (isMasterProc) then
        do iout = 1, Nout
            if ((type_Out(iout) .eq. 'V')) then
                call shd_io__outputCamRad(&
                    &fileInp_sav, fileOut_sav,&
                    &doDeltaM,&
                    &splitAcc, shAcc, solAcc, maxIter, totIter,&
                    &ilos, sumCam_Out, parms_Out(1, iout), file_Out(iout),&
                    &rtm_obj, grid_obj)
            end if
        end do
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! save state
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (maxIter .gt. 0) then
        call shd_io__saveState(fileOut_sav, work, fluxTot, rtm_obj, grid_obj)
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
    !  MAXMB_OUT = 4*(shdo_obj_%Nmu*(2 + 2*shdo_obj_%Nphi + 2*shdo_obj_%Nlm + 2*33*32)&
    !      + 4.5*prp_obj_%maxNpts + prp_obj_%maxNptsLeg + prp_obj_%Npha*prp_obj_%NlegSet*(prp_obj_%maxNleg + 1)&
    !      + 16.5*grid_obj%Ncell + grid_obj%Npts*(28 + shdo_obj_%maxNphi0*shdo_obj_%Nstokes)&
    !      + shdo_obj_%Nstokes*shdo_obj_%Nsh*NbigArr)/1024**2
    ! now:
    !   use NmaxNpts
    MAXMB_OUT = 4*(memUsed + 16.5*grid_obj%Ncell&
        + grid_obj%Npts*(NmaxNpts-16.5*CELL_TO_POINT_RATIO))/1024**2

    ADAPTGRIDFAC_OUT = float(grid_obj%Npts)/grid_obj%Npts0
    SHTERM_FAC_OUT = float(shdo_obj_%Nsh)/(shdo_obj_%Nlm*grid_obj%Npts)
    CELL_POINT_OUT = float(grid_obj%Ncell)/grid_obj%Npts
    if (verbose_ .ge. 5) write (6, '(A)')'Message [shdom]: Actual MAX_TOTAL_MB, ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO'
    if (verbose_ .ge. 5) write (6, '(26X,F9.2,6X,F7.3,14X,F5.3,15X,F5.3)') MAXMB_OUT, ADAPTGRIDFAC_OUT, SHTERM_FAC_OUT, CELL_POINT_OUT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    call shd_mpi__end(rtm_obj%runName, grid_obj, prpAll_obj)


    ! Deallocate array
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! rtm_obj
    !╭────────────────────────────────────────────────╮!
    deallocate(x_Grid, y_Grid, z_Grid)
    deallocate(ptr_Bc, rad_Bc)
    deallocate(temp_Rtm, planck_Rtm)
    deallocate(ext_Rtm, ssa_Rtm, ipha_Rtm)
    deallocate(leg_Rtm)
    deallocate(fluxTot, fluxDir)
    deallocate(sumRad_Out, sumCam_Out)
    deallocate(work, work1, work2)
    if (rtm_obj%isAllocated) then
        nullify(rtm_obj%x, &
                rtm_obj%y, &
                rtm_obj%z, &
                rtm_obj%temp,   &
                rtm_obj%planck, &
                rtm_obj%ext, &
                rtm_obj%ssa, &
                rtm_obj%leg, &
                rtm_obj%ipha)
    end if

    !╰────────────────────────────────────────────────╯!

    ! grid_obj
    !╭────────────────────────────────────────────────╮!
    deallocate(xyz_Grid)
    deallocate(flag_Cell, i8Grid_Cell, i6Cell_Cell, i2Tree_Cell)
    if (grid_obj%isAllocated) then
        nullify(grid_obj%xyz, &
                grid_obj%flag, &
                grid_obj%i8Grid, &
                grid_obj%i6Cell, &
                grid_obj%i2Tree)
    end if
    !╰────────────────────────────────────────────────╯!

    ! prp_obj_
    !╭────────────────────────────────────────────────╮!
    deallocate(temp_Prp_, ext_Prp_, ssa_Prp_, leg_Prp_, ipha_Prp_)
    if (prp_obj_%isAllocated) then
        deallocate(prp_obj_%z)
        prp_obj_%isAllocated = .false.
    end if
    !╰────────────────────────────────────────────────╯!

    ! sfc_obj_
    !╭────────────────────────────────────────────────╮!
    if (allocated(brdf_Sfc_)) then
        deallocate(brdf_Sfc_, brdf_Rtm_)
    end if
    !╰────────────────────────────────────────────────╯!

    ! ckd_obj_
    !╭────────────────────────────────────────────────╮!
    if (allocated(kabs_Ckd_)) then
        deallocate(solar_Ckd_, wgt_Ckd_, kabs_Ckd_)
        if (ckd_obj_%isAllocated) then
            deallocate(ckd_obj_%z)
            ckd_obj_%isAllocated = .false.
        end if
    else
        deallocate(solar_Ckd_, wgt_Ckd_)
    end if
    !╰────────────────────────────────────────────────╯!

    ! sen_obj_
    !╭────────────────────────────────────────────────╮!
    if (allocated(view_Sen_)) then
        deallocate(view_Sen_)
    end if
    !╰────────────────────────────────────────────────╯!

    ! shdo_obj_
    !╭────────────────────────────────────────────────╮!
    deallocate(radI_Shdo_, iradI_Shdo_)
    deallocate(srcJ_Shdo_, isrcJ_Shdo_)
    if (shdo_obj_%isAllocated) then
        deallocate(shdo_obj_%Nphi0, &
                   shdo_obj_%mu,    &
                   shdo_obj_%phi,   &
                   shdo_obj_%wgt)
    end if
    !╰────────────────────────────────────────────────╯!

    ! out_obj
    !╭────────────────────────────────────────────────╮!
    if (isFluxOut) &
        deallocate(sumFluxTot_Out, sumFluxDir_Out)
    if (isHROut)   &
        deallocate(sumFluxDiv_Out)
    if (isSHOut)   &
        deallocate(sumSH_Out)
    if (isSrcOut)  &
        deallocate(src1_Out)
    !╰────────────────────────────────────────────────╯!
    !╰────────────────────────────────────────────────────────────────────────────╯!

end program shdom
