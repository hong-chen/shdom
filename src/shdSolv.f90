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
    use globals_type, only: shd_type_prp, shd_type_sfc, shd_type_ckd

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

    use shd_mpi,  only: shd_mpi__findAlbMax, &
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
    &NSTOKES, NX, NY, NZ, NX1, NY1, NANG,&
    &ML, MM, NLM, NMU, NPHI,&
    &NPHI0, MU, PHI, WTDO,&
    &MAXIV, MAXIC, MAXIG, MAXIDO, INRADFLAG,&
    &BCFLAG, IPFLAG, DELTAM, SRCTYPE, isHighOrder,&
    &SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &UNITS, WAVENO, WAVELEN,&
    &ACCELFLAG, SOLACC, MAXITER, SOLCRIT, ITER,&
    &SPLITACC, SHACC, XGRID, YGRID, ZGRID,&
    &TEMP, PLANCK, EXTINCT, ALBEDO, LEGEN, IPHASE,&
    &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &NPTS, GRIDPOS,&
    &NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &RSHPTR, SHPTR, OSHPTR, WORK, WORK1, WORK2,&
    &SOURCE, DELSOURCE, RADIANCE, FLUXES, DIRFLUX,&
    &prp_obj, sfc_obj, ckd_obj)

    ! Performs the SHDOM solution procedure.
    ! Output is returned in SOURCE, RADIANCE, FLUXES, DIRFLUX.
    ! The adaptive grid stucture (NCELLS,GRIDPTR,NEIGHPTR,TREEPTR,
    ! CELLFLAGS,NPTS,GRIDPOS) is input and modified for new grid
    ! cells and points on output.

    type(shd_type_prp), intent(in)    :: prp_obj
    type(shd_type_sfc), intent(inout) :: sfc_obj
    type(shd_type_ckd), intent(in)    :: ckd_obj

    integer, intent(in) :: NSTOKES, NX, NY, NZ, NX1, NY1
    integer :: ML, MM, NLM, NMU, NPHI, NANG
    integer NPHI0(*), MAXITER, ITER, BCFLAG, IPFLAG
    integer MAXIV, MAXIC, MAXIG, MAXIDO
    integer MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR(MAXNBC, 2)
    integer NPTS, NCELLS
    integer RSHPTR(*), SHPTR(*), OSHPTR(*)
    integer GRIDPTR(8, *), NEIGHPTR(6, *), TREEPTR(2, *)
    integer*2 CELLFLAGS(*), IPHASE(*)
    integer WORK1(*)
    logical INRADFLAG, DELTAM, ACCELFLAG, isHighOrder
    real SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD
    real WAVENO(2), WAVELEN
    real SOLACC, SOLCRIT, SPLITACC, SHACC
    real MU(*), WTDO(*), PHI(*)
    real TEMP(*), PLANCK(*)
    real EXTINCT(*), ALBEDO(*), LEGEN(prp_obj%NlegSet, 0:prp_obj%Nleg, *)
    real XGRID(NX1), YGRID(NY1), ZGRID(NZ)
    real GRIDPOS(3, *)
    real BCRAD(*)
    real FLUXES(2, *), DIRFLUX(*)
    real SOURCE(NSTOKES, *), DELSOURCE(NSTOKES, *)
    real RADIANCE(NSTOKES, *)
    real WORK(*), WORK2(*)
    character(len=1) :: SRCTYPE, UNITS

    integer NPHI0MAX, I, ORDINATESET, NBPTS, NBCELLS
    integer MAXNANGBND, NANGBND(4, 2), IERR
    logical FIXSH, SPLITTESTING, DOSPLIT, OUTOFMEM
    logical, allocatable :: FFTFLAG(:)
    real ALBMAX
    real STARTADAPTSOL, ENDADAPTSOL, ADAPTRANGE, SPLITCRIT
    real STARTSPLITACC, CURSPLITACC, AVGSOLCRIT, BETA
    real DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR
    real, allocatable :: CMU1(:, :, :), CMU2(:, :, :)
    real, allocatable :: CPHI1(:, :, :), CPHI2(:, :, :)
    real, allocatable :: YLMSUN(:, :), WTMU(:), WPHISAVE(:)
    real, allocatable :: SFC_BRDF_DO(:, :, :, :, :, :)

    integer :: COUNT_START, COUNT_RATE, COUNT_MAX
    integer :: COUNT_END
    real :: ELAPSED_TIME

    ! Allocate local arrays
    ! prp_obj%NlegSet=1 for NSTOKES=1, otherwise prp_obj%NlegSet=6
    allocate (CMU1(prp_obj%NlegSet, NLM, NMU), CMU2(prp_obj%NlegSet, NLM, NMU))
    allocate (CPHI1(-16:16, 32, NMU), CPHI2(32, -16:16, NMU))
    allocate (YLMSUN(prp_obj%NlegSet, NLM))
    allocate (FFTFLAG(NMU), WTMU(NMU), WPHISAVE(NMU*(3*NPHI + 15)))

    ! Set up some things before solution loop

    ! Transfer the medium properties to the internal grid and add gas absorption
    call shd_grid__interpGrid(NPTS, GRIDPOS,&
        &TEMP, EXTINCT, ALBEDO, LEGEN, IPHASE, prp_obj, ckd_obj)

    ! If Delta-M then scale the extinction, albedo, and Legendre terms
    ! Put the Planck blackbody source in PLANCK.
    call shd_init__makePrp(ML, MM, NPTS, DELTAM,&
        &SRCTYPE, UNITS, WAVENO, WAVELEN, ALBMAX,&
        &EXTINCT, ALBEDO, LEGEN, TEMP, PLANCK, IPHASE, prp_obj)

    ! Get the maximum single scattering albedo over all processors
    call shd_mpi__findAlbMax(ALBMAX)

    ! Compute the solar transmission in DIRFLUX.
    ! Precompute Ylm's for solar direction
    if (isSolar_) then
        if (isMpi_) then
            call shd_mpi__makeDirect(1, NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                &SOLARFLUX, SOLARMU, SOLARAZ,&
                &GRIDPOS, NX, XGRID, NY, YGRID, DIRFLUX, prp_obj, ckd_obj)
        else
            call shd_rad__makeDirect(NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
                &SOLARFLUX, SOLARMU, SOLARAZ,&
                &GRIDPOS, DIRFLUX, prp_obj, ckd_obj)
        end if
        call shd_shdo__computeYlm(.true., SOLARMU, SOLARAZ, ML, MM, prp_obj%NlegSet,&
            &YLMSUN)
    end if

    ! Make the discrete ordinates (angles)
    !   (set 2 is reduced gaussian, 3 is reduced double gauss)
    ORDINATESET = 3
    call shd_shdo__makeAngle(NMU, NPHI, ORDINATESET, NPHI0MAX,&
        &NPHI0, MU, PHI, WTMU, WTDO, FFTFLAG, NANG)

    ! Make the Ylm transform coefficients
    call shd_shdo__makeSHDO(prp_obj%NlegSet, ML, MM, NLM, NMU, NPHI0,&
        &NPHI0MAX, MU, PHI, WTMU, WTDO, FFTFLAG,&
        &CMU1, CMU2, CPHI1, CPHI2, WPHISAVE)

    FIXSH = .false.
    if (.not. INRADFLAG) then
        ! Initialize the radiance on the base grid using Eddington
        ! two-stream plane-parallel
        call shd_init__makeRad(NSTOKES, NX1*NY1, NZ,&
            &RSHPTR, ZGRID,&
            &EXTINCT, ALBEDO, LEGEN, TEMP, IPHASE,&
            &SRCTYPE, SOLARFLUX, SOLARMU, SKYRAD,&
            &UNITS, WAVENO, WAVELEN, RADIANCE, prp_obj, sfc_obj)
        ! Interpolate the radiance on the non-base grid points from the
        !   base grid points
        NBPTS = NX1*NY1*NZ
        NBCELLS = (NZ - 1)*(NX + ibits(BCFLAG, 0, 1) - ibits(BCFLAG, 2, 1))&
            &*(NY + ibits(BCFLAG, 1, 1) - ibits(BCFLAG, 3, 1))
        call shd_rad__interpRad(NSTOKES, NBPTS, NPTS, RSHPTR, RADIANCE,&
            &NBCELLS, NCELLS, TREEPTR, GRIDPTR, GRIDPOS)
        ! Initialize the source function from the radiance field
        call shd_src__computeSrc(NSTOKES, ML, MM, NLM,&
            &NPTS,&
            &FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            &IPHASE, PLANCK, DIRFLUX, SHACC,&
            &MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            &.true., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM, prp_obj)
    end if

    if (ACCELFLAG) then
        OSHPTR(1:NPTS + 1) = SHPTR(1:NPTS + 1)
        DELSOURCE(:, 1:OSHPTR(NPTS + 1)) = 0.0
    end if
    SOLCRIT = 1.0


    ! Surface initialization
    !╭────────────────────────────────────────────────────────────────────────────╮!
    ! Make the pointers to the grid points on the top and bottom boundaries
    call shd_grid__findBndryPnt(NPTS, NANG, MAXNBC, MAXBCRAD,&
        &ZGRID(1), ZGRID(NZ), GRIDPOS,&
        &NTOPPTS, NBOTPTS, BCPTR, sfc_obj)
    if (sfc_obj%type(1:1) .eq. 'V') then
        ! If this is a variable surface then bilinearly interpolate
        !   the surface temperature and reflection parameters to the
        !   bottom grid points.
        call shd_sfc__interpSfc(NBOTPTS, BCPTR(1, 2), GRIDPOS,&
            &SRCTYPE, WAVENO, WAVELEN, UNITS,&
            &prp_obj, sfc_obj)
        call shd_init__checkSfcU(sfc_obj)
        if (sfc_obj%isUniform) then
            if (NSTOKES*(NMU/2)*NPHI0MAX .gt.&
                &sqrt(1.0*huge(NMU)/kind(NMU))) then
                IERR = 1
            else
                allocate (SFC_BRDF_DO(NSTOKES, NMU/2, NPHI0MAX,&
                    &NSTOKES, NMU/2, NPHI0MAX), STAT=IERR)
            end if
            if (IERR .eq. 0) then
                call shd_shdo__makeSfcDO(sfc_obj%type, WAVELEN,&
                    &sfc_obj%Nx, sfc_obj%Ny, sfc_obj%Npar, sfc_obj%brdf, NSTOKES,&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                    &SFC_BRDF_DO)
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
    OUTOFMEM = .false.
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
                                 &(.not. OUTOFMEM))))

        call system_clock(COUNT_START, COUNT_RATE, COUNT_MAX)

        ITER = ITER + 1

        ! Don't do the adaptive grid cell splitting stuff unless wanted
        if (SPLITACC .gt. 0.0) then
            ! If in the right solution range then try to split some cells, and
            !   interpolate the medium and source functions to the new points.
            ! The current splitting accuracy is started high and then lowered.
            AVGSOLCRIT = sqrt(AVGSOLCRIT*SOLCRIT)
            DOSPLIT = (SOLCRIT .le. STARTADAPTSOL) .and.&
                     &((SOLCRIT .gt. ENDADAPTSOL) .or.&
                      &(CURSPLITACC .gt. SPLITACC)) .and.&
                     &(.not. OUTOFMEM)
            ! Make sure all processors are going to split cells if any want to
            call shd_mpi__unifySplit(DOSPLIT, STARTSPLITACC)
            BETA = log(STARTSPLITACC/SPLITACC)/log(ADAPTRANGE)
            CURSPLITACC = min(CURSPLITACC, max(SPLITACC,&
                &SPLITACC*(AVGSOLCRIT/(3.0*ENDADAPTSOL))**BETA))
            if (SOLCRIT .le. ENDADAPTSOL) CURSPLITACC = SPLITACC
            if (SPLITTESTING) then
                call shd_grid__splitGrid(MAXIG, MAXIC, MAXIV, MAXIDO, NPHI0MAX,&
                    &DOSPLIT, OUTOFMEM, CURSPLITACC, SPLITCRIT,&
                    &NPTS, NCELLS, GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &GRIDPOS, NX, XGRID, NY, YGRID,&
                    &NSTOKES, ML, MM, NLM,&
                    &DELTAM, BCFLAG, IPFLAG, ACCELFLAG,&
                    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
                    &RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR,&
                    &TEMP, PLANCK, DIRFLUX, SRCTYPE,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, YLMSUN,&
                    &UNITS, WAVENO, WAVELEN, WORK, WORK2, prp_obj, ckd_obj)
                if (SOLCRIT .gt. STARTADAPTSOL) STARTSPLITACC = SPLITCRIT
            end if
            ! Find the maximum splitting criterion over all processors
            call shd_mpi__findTotSplitCritMax(SPLITCRIT)
            if (SOLCRIT .le. ENDADAPTSOL) SPLITTESTING = .false.
        else
            DOSPLIT = .false.
        end if

        ! Find the SH truncation to use for the radiance (RSHPTR)
        call shd_solv__truncRad(isHighOrder,&
            &NSTOKES, ML, MM,&
            &NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
            &MAXIV + MAXIG, FIXSH, SHACC, RSHPTR, prp_obj)

        ! Integrate the source function along discrete ordinates
        !   to compute the radiance field.  The input source function
        !   and output radiance are in spherical harmonic expansions,
        !   and the SH/DO transforms are done during the loop over
        !   discrete ordinates.
        call shd_intg__integratePath(NX, NY, NZ, NSTOKES, ML, MM, NLM,&
            &NMU, NPHI0MAX, NPHI0, NANG, BCFLAG, IPFLAG,&
            &NPTS, NCELLS, XGRID, YGRID, ZGRID, GRIDPOS,&
            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
            &SRCTYPE, SOLARMU, SOLARAZ, SKYRAD, WAVENO, WAVELEN,&
            &UNITS, SFC_BRDF_DO,&
            &MAXNBC, MAXBCRAD, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
            &MU, PHI, WTDO,&
            &CMU1, CPHI1, CMU2, CPHI2, WPHISAVE, FFTFLAG,&
            &DIRFLUX, FLUXES, EXTINCT,&
            &SHPTR, SOURCE, RSHPTR, RADIANCE,&
            &WORK, WORK1, WORK2, prp_obj, sfc_obj)

        ! Compute the source function from the radiance field,
        !   do the adaptive spherical harmonics truncation, compute
        !   the solution criterion, and dot products for acceleration
        if (((SOLCRIT .lt. ENDADAPTSOL) .or. (ITER .gt. 30)) .and.&
            &((NX .gt. 36) .or. (NY .gt. 36))) FIXSH = .true.
        call shd_src__computeSrc(NSTOKES, ML, MM, NLM,&
            &NPTS,&
            &FIXSH, SRCTYPE, SOLARMU, YLMSUN, ALBEDO, LEGEN,&
            &IPHASE, PLANCK, DIRFLUX, SHACC,&
            &MAXIV, RSHPTR, RADIANCE, SHPTR, SOURCE, OSHPTR, DELSOURCE,&
            &.false., ACCELFLAG, DELJDOT, DELJOLD, DELJNEW, JNORM, prp_obj)

        ! Calculate the acceleration parameter and solution criterion
        !   from all the processors
        call shd_mpi__calAccelSolCrit(ACCELFLAG, DELJDOT, DELJOLD, DELJNEW,&
            &JNORM, ACCELPAR, SOLCRIT)

        ! Accelerate the convergence of the source function vector
        call shd_solv__accelSol(ACCELPAR, NPTS, SHPTR, OSHPTR,&
            &NSTOKES, SOURCE, DELSOURCE)

        ! If it is a not scattering medium then do only one iteration
        if (ALBMAX .lt. SOLACC) SOLCRIT = SOLACC

        call system_clock(COUNT_END, COUNT_RATE, COUNT_MAX)
        ELAPSED_TIME = (real(COUNT_END) - real(COUNT_START)) / real(COUNT_RATE)

        ! Print out the log solution criterion, number of points,
        !   and average SH truncation.
        if (DOSPLIT) then
            if (verbose_ .ge. 5) write (6, '("S",1X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                &SPLITCRIT,&
                &NPTS,&
                &NCELLS,&
                &float(SHPTR(NPTS + 1))/NPTS,&
                &ELAPSED_TIME
                !&float(SHPTR(NPTS + 1))/(NPTS*NLM)
        else
            if (OUTOFMEM) then
                if (verbose_ .ge. 5) write (6, '("M",1X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                    &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                    &SPLITCRIT,&
                    &NPTS,&
                    &NCELLS,&
                    &float(SHPTR(NPTS + 1))/NPTS,&
                    &ELAPSED_TIME
                    !&float(SHPTR(NPTS + 1))/(NPTS*NLM)
            else
                if (verbose_ .ge. 5) write (6, '(2X,I4,2X,F7.3,2X,E10.3,1X,I8,1X,I8,1X,F8.2,2X,F6.2," sec")')&
                    &ITER, log10(max(SOLCRIT, 1.0e-20)),&
                    &SPLITCRIT,&
                    &NPTS,&
                    &NCELLS,&
                    &float(SHPTR(NPTS + 1))/NPTS,&
                    &ELAPSED_TIME
                    !&float(SHPTR(NPTS + 1))/(NPTS*NLM)
            end if
        end if
    end do PICARD_ITER_LOOP

    if (verbose_ .ge. 5) write (6, '(A,I4,A,F9.6)') '- Iterations: ', ITER,&
        &', Final Criterion: ', SOLCRIT
    !╰────────────────────────────────────────────────────────────────────────────╯!


    ! Deallocate arrays
    !╭────────────────────────────────────────────────────────────────────────────╮!
    if (sfc_obj%type(1:1) .eq. 'V' .and. sfc_obj%isUniform) then
        deallocate (SFC_BRDF_DO)
    end if
    deallocate (YLMSUN, FFTFLAG, WTMU, WPHISAVE)
    deallocate (CMU1, CMU2, CPHI1, CPHI2)
    !╰────────────────────────────────────────────────────────────────────────────╯!

    return

end subroutine shd_solv__solveRte

subroutine shd_solv__accelSol(&
    &ACCELPAR, NPTS, SHPTR, OSHPTR,&
    &NSTOKES, SOURCE, DELSOURCE)

    !   Accelerates the successive order of scattering series using the
    ! input acceleration parameter (ACCELPAR).

    integer, intent(in) :: NPTS, SHPTR(*), OSHPTR(*), NSTOKES
    real, intent(in)    :: ACCELPAR, DELSOURCE(NSTOKES, *)

    real, intent(inout) :: SOURCE(NSTOKES, *)

    integer :: I, J, IS, ISD, K, NS, NSD, NSC

    if (ACCELPAR .gt. 0.0) then
        do I = 1, NPTS
            IS = SHPTR(I)
            NS = SHPTR(I + 1) - IS
            ISD = OSHPTR(I)
            NSD = OSHPTR(I + 1) - ISD
            NSC = min(NS, NSD)
            do J = 1, NSC
                SOURCE(:, IS + J) = SOURCE(:, IS + J)&
                    &+ ACCELPAR*DELSOURCE(:, ISD + J)
            end do
        end do
    end if

    return
end subroutine shd_solv__accelSol

subroutine shd_solv__truncRad(isHighOrder,&
    &NSTOKES, ML, MM,&
    &NPTS, ALBEDO, LEGEN, IPHASE, SHPTR, RADIANCE,&
    &MAXIR, FIXSH, SHACC, RSHPTR, prp_obj)

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
    ! unless isHighOrder is on, in which case all terms are kept.

    type(shd_type_prp), intent(in) :: prp_obj

    integer NPTS, NSTOKES, ML, MM, MAXIR
    integer SHPTR(NPTS + 1), RSHPTR(NPTS + 2)
    integer*2 IPHASE(NPTS)
    logical isHighOrder, FIXSH
    real ALBEDO(NPTS), LEGEN(prp_obj%NlegSet, 0:prp_obj%Nleg, *), RADIANCE(NSTOKES, *)
    real SHACC
    integer IRO, IR, NRO, NR, I, J, K, L, LR, LS, M, MS, ME, NS, NLM
    integer, allocatable :: LOFJ(:)
    logical NOTEND
    real RADMIN

    if (FIXSH) goto 190

    RADMIN = SHACC
    NOTEND = .true.
    ! Make the l index as a function of SH term (J)
    NLM = (2*MM + 1)*(ML + 1) - MM*(MM + 1)
    allocate (LOFJ(NLM))
    J = 0
    do L = 0, ML
        ME = min(L, MM)
        MS = -ME
        do M = MS, ME
            J = J + 1
            LOFJ(J) = L
        end do
    end do

    IRO = 0
    IR = 0
    RSHPTR(1) = 0
    do I = 1, NPTS
        ! Find how many terms wanted for the radiance based on the
        !   radiance threshold (RADMIN) (use the mean radiance as an
        !   upper limit to the SH term radiance).
        !   We don't have radiances for any new points, so don't use them
        NRO = RSHPTR(I + 1) - IRO
        if (NRO .eq. 0) NOTEND = .false.
        if (NOTEND) then
            if (prp_obj%Npha .gt. 0) then
                K = IPHASE(I)
            else
                K = I
            end if
            LR = 1
            do L = 1, ML
                if (ALBEDO(I)*LEGEN(1, L, K)*RADIANCE(1, IRO + 1) .gt. RADMIN) LR = L
            end do
        else
            LR = ML
        end if
        IRO = IRO + NRO

        ! The radiance L truncation is the smaller of the radiance cutoff
        !   from above and a few more than the source truncation, but
        !   it is at least 1 (so low order terms are always available).
        !   For higher order radiance option, use all terms.
        !   But if out of memory, use the source function truncation.
        NS = max(1, SHPTR(I + 1) - SHPTR(I))
        LS = LOFJ(NS)
        LR = min(LR, LS + ML/8 + 2)
        if (isHighOrder) LR = ML
        if (LR .le. MM) then
            NR = (LR*(LR + 1)) + LR + 1
        else
            NR = (2*MM + 1)*LR - (MM*(1 + (MM - 1))) + MM + 1
        end if
        IR = IR + NR
        RSHPTR(I + 1) = IR
        if (IR .gt. MAXIR) then
            if (verbose_ .ge. 4) write (6, *)&
                &'Warning [_truncRad]: Out of memory for more radiance terms.'
            goto 190
        end if
    end do
    RSHPTR(NPTS + 2) = IR
    deallocate (LOFJ)
    return

    ! If out of memory or the SH truncation is fixed then
    ! just set the radiance truncation to that of the source.
190 continue

    RSHPTR(1) = 0
    IR = 0
    do I = 1, NPTS
        NR = max(4, SHPTR(I + 1) - SHPTR(I))
        if (isHighOrder) then
            if (ML .le. MM) then
                NR = (ML*(ML + 1)) + ML + 1
            else
                NR = (2*MM + 1)*ML - (MM*(1 + (MM - 1))) + MM + 1
            end if
        end if
        IR = IR + NR
        if (IR .gt. MAXIR) then
            write (6, *) 'Error [_truncRad]: Really out of memory ',&
                &'for more radiance terms.'
            write (6, *) '  Increase MAXIV.'
            stop
        end if
        RSHPTR(I + 1) = IR
    end do
    RSHPTR(NPTS + 2) = IR

    return

end subroutine shd_solv__truncRad

end module shd_solv
