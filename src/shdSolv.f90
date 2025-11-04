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

module shd_solv

    use globals_cfg, only:  verbose_, isMpi_, &
                            isSolar_
    use globals_type, only: shd_type_rtm,  &
                            shd_type_grid

    use globals_obj, only: shdo_obj, iradI_Shdo, radI_Shdo, isrcJ_Shdo, srcJ_Shdo, &
                            src_obj, &
                            prp_obj, &
                            sfc_obj, brdf_Sfc, brdf_Rtm

    use shd_init, only: shd_init__makePrp, &
                        shd_init__makeRad, &
                        shd_init__checkSfcU

    use shd_src,  only: shd_src__computeSrc

    use shd_rad,  only: shd_rad__makeDirect, &
                        shd_rad__interpRad

    use shd_sfc,  only: shd_sfc__interpSfc

    use shd_shdo, only: shd_shdo__computeYlm, &
                        shd_shdo__makeAngle, &
                        shd_shdo__makeSHDO, &
                        shd_shdo__makeSfcDO

    use shd_grid, only: shd_grid__interpGrid, &
                        shd_grid__findBndryPnt, &
                        shd_grid__splitGrid

    use shd_intg, only: shd_intg__integratePath

    use shd_mpi,  only: shd_mpi__findMaxSsa, &
                        shd_mpi__makeDirect, &
                        shd_mpi__unifySplit, &
                        shd_mpi__findTotSplitCritMax, &
                        shd_mpi__calAccelSolCrit

    implicit none
    private

    integer, parameter :: R__  = selected_real_kind(6) ! default precision for float
    integer, parameter :: I__  = selected_int_kind(9)  ! default precision for integer

    public  :: shd_solv__solveRte

    private :: shd_solv__accelSol
    private :: shd_solv__truncRad

contains

subroutine shd_solv__solveRte(&
    &doRadInit,&
    &doDeltaM, doHighOrder,&
    &ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,&
    &SPLITACC, SHACC,&
    &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &isrcD_Shdo, WORK, WORK1, WORK2,&
    &srcD_Shdo, FLUXES, DIRFLUX,&
    &rtm_obj, grid_obj)

    ! Performs the SHDOM solution procedure.
    ! Output is returned in srcJ_Shdo, radI_Shdo, FLUXES, DIRFLUX.
    ! The adaptive grid stucture (grid_obj%Ncell,grid_obj%i8Grid,grid_obj%i6Cell,grid_obj%i2Tree,
    ! grid_obj%flag,grid_obj%Npts,grid_obj%xyz) is input and modified for new grid
    ! cells and points on output.

    type(shd_type_rtm),  intent(inout) :: rtm_obj
    type(shd_type_grid), intent(inout) :: grid_obj

    integer MAXITER, ITER
    integer MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR(MAXNBC, 2)
    integer isrcD_Shdo(*)
    integer WORK1(*)
    logical doRadInit, doDeltaM, ACCELFLAG, doHighOrder
    real SOLACC, SOLCRIT, SPLITACC, SHACC
    real BCRAD(*)
    real FLUXES(2, *), DIRFLUX(*)
    real srcD_Shdo(shdo_obj%Nstokes, *)
    real WORK(*), WORK2(*)

    integer I, isetDO
    integer MAXNANGBND, NANGBND(4, 2), ierr
    logical doFixSH, SPLITTESTING, doSplit, isOutMem
    logical, allocatable :: flagFFT_Shdo(:)
    real maxSsa
    real STARTADAPTSOL, ENDADAPTSOL, ADAPTRANGE, SPLITCRIT
    real STARTSPLITACC, CURSPLITACC, AVGSOLCRIT, BETA
    real DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR

    real, allocatable :: CMU1(:, :, :), CMU2(:, :, :)
    real, allocatable :: CPHI1(:, :, :), CPHI2(:, :, :)
    real, allocatable :: YlmSolar_Shdo(:, :), wgtMu_Shdo(:), wPhiFFT_Shdo(:)
    real, allocatable :: brdf_Shdo(:, :, :, :, :, :)

    integer :: COUNT_START, COUNT_RATE, COUNT_MAX
    integer :: COUNT_END
    real :: ELAPSED_TIME

    ! Allocate local arrays
    ! prp_obj%NlegSet=1 for shdo_obj%Nstokes=1, otherwise prp_obj%NlegSet=6
    allocate(CMU1(prp_obj%NlegSet, shdo_obj%Nlm, shdo_obj%Nmu), CMU2(prp_obj%NlegSet, shdo_obj%Nlm, shdo_obj%Nmu))
    allocate(CPHI1(-16:16, 32, shdo_obj%Nmu), CPHI2(32, -16:16, shdo_obj%Nmu))
    allocate(YlmSolar_Shdo(prp_obj%NlegSet, shdo_obj%Nlm))
    allocate(flagFFT_Shdo(shdo_obj%Nmu), wgtMu_Shdo(shdo_obj%Nmu), wPhiFFT_Shdo(shdo_obj%Nmu*(3*shdo_obj%Nphi+ 15)))

    ! Set up some things before solution loop

    ! Transfer the medium properties to the internal grid and add gas absorption
    call shd_grid__interpGrid(rtm_obj, grid_obj)

    ! If Delta-M then scale the extinction, albedo, and Legendre terms
    ! Put the Planck blackbody source in rtm_obj%planck.
    call shd_init__makePrp(doDeltaM, maxSsa,&
        &rtm_obj, grid_obj)

    ! Get the maximum single scattering albedo over all processors
    call shd_mpi__findMaxSsa(maxSsa)

    ! Compute the solar transmission in DIRFLUX.
    ! Precompute Ylm's for solar direction
    if (isSolar_) then
        if (isMpi_) then
            call shd_mpi__makeDirect(1, doDeltaM, shdo_obj%maxNl,&
                &src_obj%avgSolar0, src_obj%avgMu, src_obj%avgAz,&
                &DIRFLUX,&
                &rtm_obj, grid_obj)
        else
            call shd_rad__makeDirect(grid_obj%Npts, doDeltaM, shdo_obj%maxNl,&
                &src_obj%avgSolar0, src_obj%avgMu, src_obj%avgAz,&
                &DIRFLUX,&
                &rtm_obj, grid_obj)
        end if
        call shd_shdo__computeYlm(.true., src_obj%avgMu, src_obj%avgAz, shdo_obj%maxNl, shdo_obj%maxNm, prp_obj%NlegSet,&
            &YlmSolar_Shdo)
    end if

    ! Make the discrete ordinates (angles)
    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
    isetDO = 3
    call shd_shdo__makeAngle(shdo_obj%Nmu, shdo_obj%Nphi, isetDO, shdo_obj%maxNphi0,&
        &shdo_obj%Nphi0, shdo_obj%mu, shdo_obj%phi, wgtMu_Shdo, shdo_obj%wgt, flagFFT_Shdo, shdo_obj%Nang)

    ! Make the Ylm transform coefficients
    call shd_shdo__makeSHDO(prp_obj%NlegSet, shdo_obj%maxNl, shdo_obj%maxNm, shdo_obj%Nlm, shdo_obj%Nmu, shdo_obj%Nphi0,&
        &shdo_obj%maxNphi0, shdo_obj%mu, shdo_obj%phi, wgtMu_Shdo, shdo_obj%wgt, flagFFT_Shdo,&
        &CMU1, CMU2, CPHI1, CPHI2, wPhiFFT_Shdo)

    doFixSH = .false.
    if (.not. doRadInit) then
        ! Initialize the radiance on the base grid using Eddington
        ! two-stream plane-parallel
        call shd_init__makeRad(shdo_obj%Nstokes, grid_obj%Nx*grid_obj%Ny, grid_obj%Nz,&
            &iradI_Shdo, rtm_obj%z,&
            &rtm_obj%ext, rtm_obj%ssa, rtm_obj%leg, rtm_obj%temp, rtm_obj%ipha,&
            &src_obj%type, src_obj%avgSolar0, src_obj%avgMu, src_obj%skyRad,&
            &src_obj%unit, src_obj%waveNum, src_obj%avgWaveLen, radI_Shdo)
        ! Interpolate the radiance on the non-base grid points from the
        !   base grid points
        call shd_rad__interpRad(shdo_obj%Nstokes, grid_obj%Npts0, grid_obj%Npts, iradI_Shdo, radI_Shdo,&
            &grid_obj%Ncell0, grid_obj%Ncell, grid_obj%i2Tree, grid_obj%i8Grid, grid_obj%xyz)
        ! Initialize the source function from the radiance field
        call shd_src__computeSrc(&
            &doFixSH, YlmSolar_Shdo,&
            &DIRFLUX, SHACC,&
            &isrcD_Shdo, srcD_Shdo,&
            &.true., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM,&
            &rtm_obj, grid_obj)
    end if

    if (ACCELFLAG) then
        isrcD_Shdo(1:grid_obj%Npts + 1) = isrcJ_Shdo(1:grid_obj%Npts + 1)
        srcD_Shdo(:, 1:isrcD_Shdo(grid_obj%Npts + 1)) = 0.0
    end if
    SOLCRIT = 1.0


    ! Surface initialization
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Make the pointers to the grid points on the top and bottom boundaries
    call shd_grid__findBndryPnt(shdo_obj%Nang, MAXNBC, MAXBCRAD,&
        &rtm_obj%z(1), rtm_obj%z(grid_obj%Nz),&
        &NTOPPTS, NBOTPTS, BCPTR, grid_obj)
    if (sfc_obj%type(1:1) .eq. 'V') then
        ! If this is a variable surface then bilinearly interpolate
        !   the surface temperature and reflection parameters to the
        !   bottom grid points.
        call shd_sfc__interpSfc(NBOTPTS, BCPTR(1, 2), grid_obj%xyz, &
            src_obj%type, src_obj%waveNum, src_obj%avgWaveLen, src_obj%unit, &
            brdf_Rtm)
        call shd_init__checkSfcU()
        if (sfc_obj%isUniform) then
            if (shdo_obj%Nstokes*(shdo_obj%Nmu/2)*shdo_obj%maxNphi0 .gt.&
                &sqrt(1.0*huge(shdo_obj%Nmu)/kind(shdo_obj%Nmu))) then
                ierr = 1
            else
                allocate(brdf_Shdo(shdo_obj%Nstokes, shdo_obj%Nmu/2, shdo_obj%maxNphi0,&
                    &shdo_obj%Nstokes, shdo_obj%Nmu/2, shdo_obj%maxNphi0), stat=ierr)
            end if
            if (ierr .eq. 0) then
                call shd_shdo__makeSfcDO(sfc_obj%type, src_obj%avgWaveLen,&
                    &sfc_obj%Nx, sfc_obj%Ny, sfc_obj%Npar, brdf_Sfc, shdo_obj%Nstokes,&
                    &shdo_obj%Nmu, shdo_obj%maxNphi0, shdo_obj%Nphi0, shdo_obj%mu, shdo_obj%phi, shdo_obj%wgt,&
                    &brdf_Shdo)
            else
                print '(A,A)', 'Not enough memory to allocate precomputed',&
                    &' surface BRDF discrete ordinate matrix: '
                print *, '  Proceeding more slowly.'
                sfc_obj%isUniform = .false.
            end if
        end if
    end if
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Picard iteration
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ITER = 0

    if (MAXITER .le. ITER) return
    ! Starting values for the adaptive cell splitting controlling method.
    ! The splitting accuracy used at a particular iteration is CURSPLITACC,
    ! and it is started out high an lowered to the desired SPLITACC.
    ! The cell splitting starts when the solution criterion (SOLCRIT)
    ! goes below STARTADAPTSOL and ends when the solution criterion
    ! goes below ENDADAPTSOL, though the current splitting accuracy
    ! is designed to get to the final SPLITACC over a smaller range
    ! (ADAPTRANGE).  The currect splitting accuracy is computed from
    ! a geometric average of the solution criterion (AVGSOLCRIT).
    ! The shd_grid__splitGrid routine is called and the splitting criterion returned
    ! as long as SPLITTESTING is true.
    SPLITTESTING = .true.
    isOutMem = .false.
    ENDADAPTSOL = 0.001
    STARTADAPTSOL = 0.1
    ADAPTRANGE = STARTADAPTSOL/(3.0*ENDADAPTSOL)
    CURSPLITACC = SPLITACC*ADAPTRANGE
    STARTSPLITACC = CURSPLITACC
    AVGSOLCRIT = SOLCRIT
    SPLITCRIT = 0.0

    if (verbose_ .ge. 5) then
        write (6, '(A)') 'Message [_solveRte]:'
        write (6, '(A)') '+ Iter  Log(Sol)  SplitCrit  Npoints   Ncells  Nsh(avg)'
    end if

    ! Main solution loop
    !   Iterate until the iterations exceed the limit, or the desired
    !   solution criterion is reached and the current splitting accuracy
    !   is down to the desired level.
    PICARD_ITER_LOOP: do while ((ITER .lt. MAXITER) .and.&
                               &((SOLCRIT .gt. SOLACC) .or.&
                                &((SPLITCRIT .gt. SPLITACC) .and.&
                                 &(CURSPLITACC .gt. SPLITACC) .and.&
                                 &(.not. isOutMem))))

        call system_clock(COUNT_START, COUNT_RATE, COUNT_MAX)

        ITER = ITER + 1

        ! Don't do the adaptive grid cell splitting stuff unless wanted
        if (SPLITACC .gt. 0.0) then
            ! If in the right solution range then try to split some cells, and
            !   interpolate the medium and source functions to the new points.
            ! The current splitting accuracy is started high and then lowered.
            AVGSOLCRIT = sqrt(AVGSOLCRIT*SOLCRIT)
            doSplit = (SOLCRIT .le. STARTADAPTSOL) .and.&
                     &((SOLCRIT .gt. ENDADAPTSOL) .or.&
                      &(CURSPLITACC .gt. SPLITACC)) .and.&
                     &(.not. isOutMem)
            ! Make sure all processors are going to split cells if any want to
            call shd_mpi__unifySplit(doSplit, STARTSPLITACC)
            BETA = log(STARTSPLITACC/SPLITACC)/log(ADAPTRANGE)
            CURSPLITACC = min(CURSPLITACC, max(SPLITACC,&
                &SPLITACC*(AVGSOLCRIT/(3.0*ENDADAPTSOL))**BETA))
            if (SOLCRIT .le. ENDADAPTSOL) CURSPLITACC = SPLITACC
            if (SPLITTESTING) then
                call shd_grid__splitGrid(&
                    &doSplit, isOutMem, CURSPLITACC, SPLITCRIT,&
                    &doDeltaM, ACCELFLAG,&
                    &isrcD_Shdo,&
                    &DIRFLUX,&
                    &src_obj%avgSolar0, YlmSolar_Shdo,&
                    &WORK, WORK2,&
                    &rtm_obj, grid_obj)
                if (SOLCRIT .gt. STARTADAPTSOL) STARTSPLITACC = SPLITCRIT
            end if
            ! Find the maximum splitting criterion over all processors
            call shd_mpi__findTotSplitCritMax(SPLITCRIT)
            if (SOLCRIT .le. ENDADAPTSOL) SPLITTESTING = .false.
        else
            doSplit = .false.
        end if

        ! Find the SH truncation to use for the radiance (iradI_Shdo)
        call shd_solv__truncRad(doHighOrder, doFixSH, SHACC,&
            &rtm_obj, grid_obj)

        ! srcSH to srcDO, then srcDO to radDO, then radDO to radSH
        !╭────────────────────────────────────────────────╮!
        ! Integrate the source function along discrete ordinates
        !   to compute the radiance field.  The input source function
        !   and output radiance are in spherical harmonic expansions,
        !   and the SH/DO transforms are done during the loop over
        !   discrete ordinates.
        call shd_intg__integratePath(&
            &brdf_Shdo,&
            &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
            &CMU1, CPHI1, CMU2, CPHI2, wPhiFFT_Shdo, flagFFT_Shdo,&
            &DIRFLUX, FLUXES,&
            &WORK, WORK1, WORK2,&
            &rtm_obj, grid_obj)
        !╰────────────────────────────────────────────────╯!

        ! Compute the source function from the radiance field,
        !   do the adaptive spherical harmonics truncation, compute
        !   the solution criterion, and dot products for acceleration
        if (((SOLCRIT .lt. ENDADAPTSOL) .or. (ITER .gt. 30)) .and.&
            &((rtm_obj%Nx .gt. 36) .or. (rtm_obj%Ny .gt. 36))) doFixSH = .true.

        ! radSH to srcSH
        !╭────────────────────────────────────────────────╮!
        call shd_src__computeSrc(&
            &doFixSH, YlmSolar_Shdo,&
            &DIRFLUX, SHACC,&
            &isrcD_Shdo, srcD_Shdo,&
            &.false., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM,&
            &rtm_obj, grid_obj)
        !╰────────────────────────────────────────────────╯!

        ! Calculate the acceleration parameter and solution criterion
        !   from all the processors
        call shd_mpi__calAccelSolCrit(ACCELFLAG, DELJDOT, DELJOLD, DELJNEW,&
            &JNORM, ACCELPAR, SOLCRIT)

        ! Accelerate the convergence of the source function vector
        call shd_solv__accelSol(ACCELPAR, grid_obj%Npts, isrcD_Shdo,&
            &srcD_Shdo)

        ! If it is a not scattering medium then do only one iteration
        if (maxSsa .lt. SOLACC) SOLCRIT = SOLACC

        call system_clock(COUNT_END, COUNT_RATE, COUNT_MAX)
        ELAPSED_TIME = (real(COUNT_END) - real(COUNT_START)) / real(COUNT_RATE)

        ! Print out the log solution criterion, number of points,
        !   and average SH truncation.
        if (doSplit) then
            if (verbose_ .ge. 5) write (6, '("S",1X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                &SPLITCRIT,&
                &grid_obj%Npts,&
                &grid_obj%Ncell,&
                &float(isrcJ_Shdo(grid_obj%Npts + 1))/grid_obj%Npts,&
                &ELAPSED_TIME
                !&float(isrcJ_Shdo(grid_obj%Npts + 1))/(grid_obj%Npts*shdo_obj%Nlm)
        else
            if (isOutMem) then
                if (verbose_ .ge. 5) write (6, '("M",1X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                    &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                    &SPLITCRIT,&
                    &grid_obj%Npts,&
                    &grid_obj%Ncell,&
                    &float(isrcJ_Shdo(grid_obj%Npts + 1))/grid_obj%Npts,&
                    &ELAPSED_TIME
                    !&float(isrcJ_Shdo(grid_obj%Npts + 1))/(grid_obj%Npts*shdo_obj%Nlm)
            else
                if (verbose_ .ge. 5) write (6, '(2X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                    &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                    &SPLITCRIT,&
                    &grid_obj%Npts,&
                    &grid_obj%Ncell,&
                    &float(isrcJ_Shdo(grid_obj%Npts + 1))/grid_obj%Npts,&
                    &ELAPSED_TIME
                    !&float(isrcJ_Shdo(grid_obj%Npts + 1))/(grid_obj%Npts*shdo_obj%Nlm)
            end if
        end if
    end do PICARD_ITER_LOOP

    if (verbose_ .ge. 5) write (6, '(A,I4,A,F9.6)') '- Iterations: ', ITER,&
        &', Final Criterion: ', SOLCRIT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Deallocate arrays
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if ((sfc_obj%type(1:1) .eq. 'V') .and. sfc_obj%isUniform) then
        deallocate(brdf_Shdo)
    end if
    deallocate(YlmSolar_Shdo, flagFFT_Shdo, wgtMu_Shdo, wPhiFFT_Shdo)
    deallocate(CMU1, CMU2, CPHI1, CPHI2)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    return

end subroutine shd_solv__solveRte

subroutine shd_solv__accelSol(&
    &ACCELPAR, NPTS, isrcD_Shdo,&
    &srcD_Shdo)

    !   Accelerates the successive order of scattering series using the
    ! input acceleration parameter (ACCELPAR).

    integer, intent(in) :: NPTS, isrcD_Shdo(*)
    real, intent(in)    :: ACCELPAR, srcD_Shdo(shdo_obj%Nstokes, *)

    integer :: I, J, IS, ISD, K, NS, NSD, NSC

    if (ACCELPAR .gt. 0.0) then
        do I = 1, NPTS
            IS = isrcJ_Shdo(I)
            NS = isrcJ_Shdo(I + 1) - IS
            ISD = isrcD_Shdo(I)
            NSD = isrcD_Shdo(I + 1) - ISD
            NSC = min(NS, NSD)
            do J = 1, NSC
                srcJ_Shdo(:, IS + J) = srcJ_Shdo(:, IS + J)&
                    &+ ACCELPAR*srcD_Shdo(:, ISD + J)
            end do
        end do
    end if

    return
end subroutine shd_solv__accelSol

subroutine shd_solv__truncRad(doHighOrder, doFixSH, SHACC,&
    &rtm_obj, grid_obj)

    ! Computes the next radiance spherical harmonic series truncation
    ! to use based on the current source function truncation and the
    ! scattering properties.  The radiance truncation is the smaller
    ! of that based on a minimum SH term radiance and that needed for
    ! the computation of the source function.  The truncation is always
    ! made for a complete set in the L index.  The radiance SH truncation
    ! is set to a few higher L index than the source function, so that
    ! the source function truncation may grow (since initialized with
    ! a low order truncation).  The minimum radiance truncation is
    ! L=1 so that the mean radiance and net flux terms may be computed,
    ! unless doHighOrder is on, in which case all terms are kept.

    type(shd_type_rtm),  intent(in)    :: rtm_obj
    type(shd_type_grid), intent(in)    :: grid_obj

    integer MAXIR
    logical doHighOrder, doFixSH
    real SHACC
    integer IRO, IR, NRO, NR, I, J, K, L, LR, LS, M, MS, ME, NS, NLM
    integer, allocatable :: LOFJ(:)
    logical NOTEND
    real RADMIN

    MAXIR = shdo_obj%maxNptsSH + grid_obj%maxNpts

    if (doFixSH) goto 190

    RADMIN = SHACC
    NOTEND = .true.
    ! Make the l index as a function of SH term (J)
    NLM = (2*shdo_obj%maxNm + 1)*(shdo_obj%maxNl + 1) - shdo_obj%maxNm*(shdo_obj%maxNm + 1)
    allocate(LOFJ(NLM))
    J = 0
    do L = 0, shdo_obj%maxNl
        ME = min(L, shdo_obj%maxNm)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    IRO = 0
    IR = 0
    iradI_Shdo(1) = 0
    do I = 1, grid_obj%Npts
        ! Find how many terms wanted for the radiance based on the
        !   radiance threshold (RADMIN) (use the mean radiance as an
        !   upper limit to the SH term radiance).
        !   We don't have radiances for any new points, so don't use them
        NRO = iradI_Shdo(I + 1) - IRO
        if (NRO .eq. 0) NOTEND = .false.
        if (NOTEND) then
            if (prp_obj%Npha .gt. 0) then
                K = rtm_obj%ipha(I)
            else
                K = I
            end if
            LR = 1
            do L = 1, shdo_obj%maxNl
                if (rtm_obj%ssa(I)*rtm_obj%leg(1, L, K)*radI_Shdo(1, IRO + 1) .gt. RADMIN) LR = L
            end do
        else
            LR = shdo_obj%maxNl
        end if
        IRO = IRO + NRO

        ! The radiance L truncation is the smaller of the radiance cutoff
        !   from above and a few more than the source truncation, but
        !   it is at least 1 (so low order terms are always available).
        !   For higher order radiance option, use all terms.
        !   But if out of memory, use the source function truncation.
        NS = max(1, isrcJ_Shdo(I + 1) - isrcJ_Shdo(I))
        LS = LOFJ(NS)
        LR = min(LR, LS + shdo_obj%maxNl/8 + 2)
        if (doHighOrder) LR = shdo_obj%maxNl
        if (LR .le. shdo_obj%maxNm) then
            NR = (LR*(LR + 1)) + LR + 1
        else
            NR = (2*shdo_obj%maxNm + 1)*LR - (shdo_obj%maxNm*(1 + (shdo_obj%maxNm - 1))) + shdo_obj%maxNm + 1
        end if
        IR = IR + NR
        iradI_Shdo(I + 1) = IR
        if (IR .gt. MAXIR) then
            if (verbose_ .ge. 4) write (6, *)&
                &'Warning [_truncRad]: Out of memory for more radiance terms.'
            goto 190
        end if
    end do
    iradI_Shdo(grid_obj%Npts + 2) = IR
    deallocate(LOFJ)
    return

    ! If out of memory or the SH truncation is fixed then
    ! just set the radiance truncation to that of the source.
190 continue

    iradI_Shdo(1) = 0
    IR = 0
    do I = 1, grid_obj%Npts
        NR = max(4, isrcJ_Shdo(I + 1) - isrcJ_Shdo(I))
        if (doHighOrder) then
            if (shdo_obj%maxNl .le. shdo_obj%maxNm) then
                NR = (shdo_obj%maxNl*(shdo_obj%maxNl + 1)) + shdo_obj%maxNl + 1
            else
                NR = (2*shdo_obj%maxNm + 1)*shdo_obj%maxNl - (shdo_obj%maxNm*(1 + (shdo_obj%maxNm - 1))) + shdo_obj%maxNm + 1
            end if
        end if
        IR = IR + NR
        if (IR .gt. MAXIR) then
            write (6, *) 'Error [_truncRad]: Really out of memory ',&
                &'for more radiance terms.'
            write (6, *) '  Increase shdo_obj%maxNptsSH.'
            stop
        end if
        iradI_Shdo(I + 1) = IR
    end do
    iradI_Shdo(grid_obj%Npts + 2) = IR

    return

end subroutine shd_solv__truncRad

end module shd_solv
