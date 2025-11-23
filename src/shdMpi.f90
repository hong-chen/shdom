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

module shd_mpi

    use globals_const, only: RD_, PI_, DPI_, PIH_, FRAC13_, FRAC112_
    use globals_cfg,   only: verbose_, charMaxLen_, accSolv_, &
                             isMpiX_, isMpiY_, isMpi_, &
                             is1DX_, is1DYx_, is1D_, &
                             isPol_, isPol4_

    use globals_type,  only: shd_type_rtm,  &
                             shd_type_cell, &
                             shd_type_prp,  &
                             shd_type_sfc,  &
                             shd_type_ckd,  &
                             shd_type_sen

    use globals_obj,   only:  rtm_obj_, xyz_Rtm_, &
                                        ext_Rtm_, &
                                        leg_Rtm_, &
                             cell_obj_, i8Grid_Cell_, &
                                        i6Cell_Cell_, &
                                        i2Tree_Cell_, &
                                          flag_Cell_, &
                             shdo_obj_, isrcJ_Shdo_, &
                                         srcJ_Shdo_, &
                              src_obj_, &
                              prp_obj_, &
                              sfc_obj_, &
                              sen_obj_, view_Sen_

#if HAS_MPI
    use mpi_f08,  only: MPI_comm_world, MPI_comm, MPI_request, &
                        MPI_status, MPI_success, MPI_errors_return, &
                        MPI_integer, MPI_real, MPI_logical, MPI_character, &
                        MPI_integer2, MPI_double_precision, &
                        MPI_max, MPI_sum, MPI_lor, MPI_land, &
                        mpi_isend, mpi_irecv, mpi_reduce, mpi_allreduce, &
                        mpi_bcast, mpi_wtime, &
                        mpi_comm_set_errhandler

    use shd_init, only: shd_init__makeTopRad, &
                        shd_init__makeBotRadUL, &
                        shd_init__makeBotRadVL, &
                        shd_init__makeBotRadV

    use shd_util, only: shd_util__ssort, &
                        shd_util__findNextCell

    use shd_src,  only: shd_src__computeSrcOneDir, &
                        shd_src__integrateSrc, &
                        shd_src__integrateSrcOneRay

    use shd_rad,  only: shd_rad__computeDirect

    use shd_scat, only: shd_scat__computePhase

    use shd_shdo, only: shd_shdo__computeYlm
#endif

    implicit none
    private

    integer, parameter :: I__ = selected_int_kind(9)
    integer, parameter :: R__ = selected_real_kind(6)

    public :: shd_mpi__start
    public :: shd_mpi__map
    public :: shd_mpi__bcastInp
    public :: shd_mpi__bcastMem
    public :: shd_mpi__bcastPrpSize
    public :: shd_mpi__bcastPrp
    public :: shd_mpi__bcastSfcSize
    public :: shd_mpi__bcastSfc
    public :: shd_mpi__bcastCkdSize
    public :: shd_mpi__bcastCkd
    public :: shd_mpi__bcastSenSize
    public :: shd_mpi__bcastSen
    public :: shd_mpi__gatherOut
    public :: shd_mpi__end
    public :: shd_mpi__abort
    public :: shd_mpi__sumCpuTime

    public :: shd_mpi__makeDirect
    public :: shd_mpi__computeRad
    public :: shd_mpi__computeCamRad
    public :: shd_mpi__findMaxSsa
    public :: shd_mpi__findMaxSplitCrit
    public :: shd_mpi__findBndryPts
    public :: shd_mpi__calBndryRad
    public :: shd_mpi__calAccelSolCrit
    public :: shd_mpi__unifySplit

#if HAS_MPI
    private :: shd_mpi__optimizeProcDomain
    private :: shd_mpi__ran
    private :: shd_mpi__integrateRay
    private :: loadObject

    integer :: Nproc_, myProc_, iProcNeigh_(4), ixyPrp_(2, 2), ixySfc_(2, 2), ixyOut_(2, 2)
    logical, parameter :: loadBalance_ = .false.
    real(RD_) :: startTime_
    integer, allocatable :: BxyPtr_(:, :, :), iBxyPtr_(:, :, :)
    real, allocatable :: BxyPos_(:, :, :)
    type(MPI_comm) :: comm2D_
#endif

contains

subroutine shd_mpi__start(isMasterProc, Nproc)

    ! Initializes the MPI system, and gets the number of processors in use
    ! and the current processor.  If this is processor 0 then isMasterProc
    ! is returned true.  Also starts the logging files "shdom???.log"
    ! for the standard output from the non-master processors.

    logical, intent(out) :: isMasterProc
    integer, intent(out) :: Nproc

#if HAS_MPI
! note: this if clause has <else>
    integer :: ierr

    call mpi_init(ierr)
    if (ierr .ne. MPI_success) then
        write (6, *) 'Error [mpi__start]: Terminating ...'
        stop
    end if

    call mpi_comm_size(MPI_comm_world, Nproc_, ierr)
    call mpi_comm_rank(MPI_comm_world, myProc_, ierr)

    call mpi_comm_set_errhandler(MPI_comm_world, MPI_errors_return)

    startTime_ = mpi_wtime()

    if (myProc_ .eq. 0) then
        isMasterProc = .true.
        if (verbose_ .gt. 5) write (6, '(A)') 'Message [mpi__start]: Starting MPI version of SHDOM ...'
    else
        isMasterProc = .false.
    end if

    Nproc = Nproc_
#else
    isMasterProc = .true.
    Nproc = 1
#endif

    return

end subroutine shd_mpi__start

subroutine shd_mpi__map(rtmAll_obj, prpAll_obj, sfcAll_obj, &
                        rtm_obj,    prp_obj,    sfc_obj)

    ! Sets up the mapping of processors for the SHDOM domain (creating the
    ! 2D communicator, comm2D_).  Determines the subdomain property grid
    ! (prp_obj%Nx,prp_obj%Ny) and internal base grid (rtm_obj%Nx,rtm_obj%Ny) sizes, and outputs the
    ! starting position of this subdomain (prp_obj%xStart,prp_obj%yStart).  Bits 2 and 3
    ! of rtm_obj%flag_bc are set to indicate multiple processors in the X and Y
    ! directions, respectively.  Open boundary conditions (bits 0 and 1 in
    ! rtm_obj%flag_bc) are overridden by the multiple processor flags (i.e. multiple
    ! processor requires a periodic full domain).
    !   The neighboring processors are stored in iProcNeigh_(:), with the
    ! following indexing scheme: 1=-X, 2=+X, 3-=-Y, 4=+Y.   The ixyPrp_ array
    ! (module global variable) is set with the starting and ending property
    ! grid indices of this subdomain (ixyPrp_(:,1) is for X, and ixyPrp_(:,2)
    ! is for Y).  The ixyOut_ array is set with the starting index and number
    ! of grid points of the internal base grid for this subdomain (for the
    ! output arrays).

    type(shd_type_rtm), intent(in) :: rtmAll_obj
    type(shd_type_prp), intent(in) :: prpAll_obj
    type(shd_type_sfc), intent(in) :: sfcAll_obj

    type(shd_type_rtm), intent(inout) :: rtm_obj
    type(shd_type_prp), intent(inout) :: prp_obj
    type(shd_type_sfc), intent(inout) :: sfc_obj

#if HAS_MPI
! note: this if clause has <else>
    integer :: i, ierr, dims(2), coords(2)
    logical :: periodic(2) = (/.true., .true./), reorder = .false.
    real    :: xEnd_prp, yEnd_prp, xEnd_sfc, yEnd_sfc
    character(len=4) :: procStr
    character(len=charMaxLen_) :: logFile

    !if ((Nproc_ .eq. 1) .or. ((prpAll_obj%Nx .eq. 1) .and. (prpAll_obj%Ny .eq. 1))) then
    if (Nproc_ .eq. 1) then
        prp_obj%xStart = 0.0; prp_obj%yStart = 0.0
        return
    end if

    if (myProc_ .gt. 0) then
        write (procStr, '(I4.4)') myProc_
        logFile = trim(adjustl(rtm_obj%runName))//'.sHdOmNG-log.iproc-'//procStr
        open (UNIT=6, FILE=trim(adjustl(logFile)), STATUS='UNKNOWN')
    end if

    ! Partition 2D Cartesian topology:
    dims(:) = 0
    call mpi_dims_create(Nproc_, 2, dims, ierr)

    ! Create 2D Cartesian mapping of processors:
    call mpi_cart_create(MPI_comm_world, 2, dims, periodic, reorder,&
        &comm2D_, ierr)

    if (dims(1) .gt. 1) then
        rtm_obj%flag_bc = ibset(rtm_obj%flag_bc, 2)
        isMpiX_ = .true.
    else
        isMpiX_ = .false.
    end if

    if (dims(2) .gt. 1) then
        rtm_obj%flag_bc = ibset(rtm_obj%flag_bc, 3)
        isMpiY_ = .true.
    else
        isMpiY_ = .false.
    end if

    ! Multi-processor BC takes precedence over open BCs
    if (isMpiX_ .or. isMpiY_) then
        rtm_obj%flag_bc = iand(rtm_obj%flag_bc, 12)
        isMpi_ = .true.
    else
        isMpi_ = .false.
    end if

    ! Find where in the 2D grid of processors this is is at:
    call mpi_cart_coords(comm2D_, myProc_, 2, coords, ierr)

    ! Make the optimal (or default) property grid divisions (in ixyPrp_)
    call shd_mpi__optimizeProcDomain(dims, coords, rtm_obj%runName, prpAll_obj, prp_obj)

    ! Make the new property and base grid sizes and the ixyOut_ arrays
    if (dims(1) .gt. 1) then

        if (prpAll_obj%Nx .gt. 1) then

            prp_obj%Nx = ixyPrp_(2, 1) - ixyPrp_(1, 1) + 1
            rtm_obj%Nx = nint(rtmAll_obj%Nx*float(prp_obj%Nx - 1)/prpAll_obj%Nx) + 1
            if (abs(prpAll_obj%Nx/float(rtmAll_obj%Nx) - (prp_obj%Nx - 1)/float(rtm_obj%Nx - 1)) .gt. 0.0001) then
                call shd_mpi__abort('Error [mpi__map]: Even X base grid cannot be preserved&
                    & for this rtm_obj%Nx and number of processors')
            end if
            ixyOut_(1, 1) = nint(rtmAll_obj%Nx*float(ixyPrp_(1, 1) - 1)/prpAll_obj%Nx) + 1
            ixyOut_(2, 1) = nint(rtmAll_obj%Nx*float(prp_obj%Nx    - 1)/prpAll_obj%Nx)
            prp_obj%xStart = (ixyPrp_(1, 1) - 1)*prp_obj%dx
            xEnd_prp = (ixyPrp_(2, 1) - 1)*prp_obj%dx

        else

            prp_obj%Nx = 1

            ixyOut_(1, 1) = nint(real(rtmAll_obj%Nx)/real(dims(1)))*coords(1) + 1
            ixyOut_(2, 1) = min(nint(real(rtmAll_obj%Nx)/real(dims(1)))*(coords(1)+1), rtmAll_obj%Nx)

            rtm_obj%Nx = ixyOut_(2, 1) - ixyOut_(1, 1) + 2

            prp_obj%xStart = real(prp_obj%dx)/real(rtmAll_obj%Nx) * real((ixyOut_(1, 1)-1))
            xEnd_prp = real(prp_obj%dx)/real(rtmAll_obj%Nx) * real((ixyOut_(2, 1)))

            prp_obj%dx = xEnd_prp - prp_obj%xStart


        end if

    else

        ixyOut_(1, 1) = 1
        ixyOut_(2, 1) =  rtmAll_obj%Nx
        prp_obj%xStart = 0.0
        xEnd_prp = prpAll_obj%Nx*prp_obj%dx

    end if

    if (dims(2) .gt. 1) then

        if (prpAll_obj%Ny .gt. 1) then

            prp_obj%Ny = ixyPrp_(2, 2) - ixyPrp_(1, 2) + 1
            rtm_obj%Ny = nint(rtmAll_obj%Ny*float(prp_obj%Ny - 1)/prpAll_obj%Ny) + 1
            if (abs(prpAll_obj%Ny/float(rtmAll_obj%Ny) - (prp_obj%Ny - 1)/float(rtm_obj%Ny - 1)) .gt. 0.0001) then
                call shd_mpi__abort('Error [mpi__map]: Even Y base grid cannot be preserved&
                    & for this rtm_obj%Ny and number of processors')
            end if
            ixyOut_(1, 2) = nint(rtmAll_obj%Ny*float(ixyPrp_(1, 2) - 1)/prpAll_obj%Ny) + 1
            ixyOut_(2, 2) = nint(rtmAll_obj%Ny*float(prp_obj%Ny - 1)/prpAll_obj%Ny)
            prp_obj%yStart = (ixyPrp_(1, 2) - 1)*prp_obj%dy
            yEnd_prp = (ixyPrp_(2, 2) - 1)*prp_obj%dy

        else

            prp_obj%Ny = 1

            ixyOut_(1, 2) =     nint(real(rtmAll_obj%Ny)/real(dims(2)))*coords(2) + 1
            ixyOut_(2, 2) = min(nint(real(rtmAll_obj%Ny)/real(dims(2)))*(coords(2)+1), rtmAll_obj%Ny)

            rtm_obj%Ny = ixyOut_(2, 2) - ixyOut_(1, 2) + 2

            prp_obj%yStart = real(prp_obj%dy)/real(rtmAll_obj%Ny) * real((ixyOut_(1, 2)-1))
            yEnd_prp = real(prp_obj%dy)/real(rtmAll_obj%Ny) * real((ixyOut_(2, 2)))

            prp_obj%dy = yEnd_prp - prp_obj%yStart

        end if

    else

        ixyOut_(1, 2) = 1
        ixyOut_(2, 2) = rtmAll_obj%Ny
        prp_obj%yStart = 0.0
        yEnd_prp = prpAll_obj%Ny*prp_obj%dy

    end if

    if (sfcAll_obj%Nx .gt. 1) then
        sfc_obj%xStart = max(0.0, (int(prp_obj%xStart/sfc_obj%dx)-3)*sfc_obj%dx)
        xEnd_sfc       = min((sfcAll_obj%Nx+1)*sfcAll_obj%dx, (int(xEnd_prp/sfc_obj%dx)+3)*sfc_obj%dx)
        sfc_obj%Nx     = max(1, int((xEnd_sfc-sfc_obj%xStart)/sfc_obj%dx))

        ixySfc_(1, 1) = nint(sfc_obj%xStart/sfc_obj%dx) + 1
        ixySfc_(2, 1) = min(sfcAll_obj%Nx+1, ixySfc_(1, 1) + sfc_obj%Nx)

        sfc_obj%xStart = sfc_obj%dx * (ixySfc_(1, 1) - 1)
        sfc_obj%Nx     = ixySfc_(2, 1) - ixySfc_(1, 1) + 1
    else
        sfc_obj%xStart = 0.0
        sfc_obj%Nx = 1

        ixySfc_(1, 1) = 1
        ixySfc_(2, 1) = 1
    end if

    if (sfcAll_obj%Ny .gt. 1) then
        sfc_obj%yStart = max(0.0, (int(prp_obj%yStart/sfc_obj%dy)-3)*sfc_obj%dy)
        yEnd_sfc       = min((sfcAll_obj%Ny+1)*sfcAll_obj%dy, (int(yEnd_prp/sfc_obj%dy)+3)*sfc_obj%dy)
        sfc_obj%Ny     = max(1, int((yEnd_sfc-sfc_obj%yStart)/sfc_obj%dy))

        ixySfc_(1, 2) = nint(sfc_obj%yStart/sfc_obj%dy) + 1
        ixySfc_(2, 2) = min(sfcAll_obj%Ny+1, ixySfc_(1, 2) + sfc_obj%Ny)

        sfc_obj%yStart = sfc_obj%dy * (ixySfc_(1, 2) - 1)
        sfc_obj%Ny     = ixySfc_(2, 2) - ixySfc_(1, 2) + 1
    else
        sfc_obj%yStart = 0.0
        sfc_obj%Ny = 1

        ixySfc_(1, 2) = 1
        ixySfc_(2, 2) = 1
    end if

    sfc_obj%maxNpts = sfc_obj%Nx * sfc_obj%Ny

    ! Get the processors (rank numbers) of the neighbors
    call mpi_cart_shift(comm2D_, 0, 1, iProcNeigh_(1), iProcNeigh_(2), ierr)
    call mpi_cart_shift(comm2D_, 1, 1, iProcNeigh_(3), iProcNeigh_(4), ierr)

    if (verbose_ .ge. 5) write (6, '(A)') 'Message [mpi__map]:'
    if (verbose_ .ge. 5) write (6, '(2X,A,2(I0,2X))')    'Nproc in X and Y: ', dims(1:2)
    if (verbose_ .ge. 5) write (6, '(2X,A,I0)')          '  This processor: ', myProc_
    if (verbose_ .ge. 5) write (6, '(2X,A,4(I0,2X))')    '       Neighbors: ', iProcNeigh_(1:4)
    if (verbose_ .ge. 5) write (6, '(2X,A,2(F0.2," km",2X))') '         X range: ', prp_obj%xStart, xEnd_prp
    if (verbose_ .ge. 5) write (6, '(2X,A,2(F0.2," km",2X))') '         Y range: ', prp_obj%yStart, yEnd_prp
#else
    prp_obj%xStart = 0.0; prp_obj%yStart = 0.0
#endif

    return

end subroutine shd_mpi__map

subroutine shd_mpi__bcastInp(&
    &runName,&
    &fileInp_sav, fileOut_sav,&
    &Nstokes, Nx, Ny, Nz, NMU, NPHI, BCFLAG,&
    &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &UNITS, WAVENO, WAVELEN,&
    &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
    &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, parms_Out,&
    &MAX_TOTAL_MB, SPLITTING_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO,&
    &prp_obj, sfc_obj, ckd_obj, sen_obj)

    ! Broadcasts all the user input parameters (except unneeded file names)
    ! from the master processor to the slave processors.

    type(shd_type_prp), intent(inout) :: prp_obj
    type(shd_type_sfc), intent(inout) :: sfc_obj
    type(shd_type_ckd), intent(inout) :: ckd_obj
    type(shd_type_sen), intent(inout) :: sen_obj

    integer, intent(inout) :: Nstokes, Nx, Ny, Nz, NMU, NPHI, BCFLAG, IPFLAG
    integer, intent(in)    :: MAXOUT, MAXPAR
    integer, intent(inout) :: MAXITER, NUMOUT
    logical, intent(inout) :: KDIST, DELTAM, ACCELFLAG
    real, intent(inout) :: SOLARFLUX, SOLARMU, SOLARAZ
    real, intent(inout) :: SKYRAD, WAVENO(2), WAVELEN
    real, intent(inout) :: SOLACC, SPLITACC, SHACC, parms_Out(MAXPAR, MAXOUT)
    real, intent(inout) :: MAX_TOTAL_MB, SPLITTING_FACTOR
    real, intent(inout) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(len=1), intent(inout) :: SRCTYPE, GRIDTYPE, UNITS, OUTTYPES(*)
    character(len=charMaxLen_), intent(inout) :: runName
    character(len=charMaxLen_), intent(inout) :: fileInp_sav, fileOut_sav

#if HAS_MPI
    integer :: ierr, intbuf(10)
    logical :: logbuf(3)
    real    :: realbuf(16)
    character(len=1) :: char1buf(3)
    character(len=charMaxLen_) :: char256buf(7)
    character(len=4) :: procStr

    if (Nproc_ .le. 1) return

    intbuf(:) = (/Nstokes, Nx, Ny, Nz, NMU, NPHI, BCFLAG, IPFLAG, MAXITER, NUMOUT/)
    call mpi_bcast(intbuf, 10, MPI_integer, 0, MPI_comm_world, ierr)
    Nstokes = intbuf(1); Nx = intbuf(2); Ny = intbuf(3); Nz = intbuf(4)
    NMU = intbuf(5); NPHI = intbuf(6); BCFLAG = intbuf(7); IPFLAG = intbuf(8)
    MAXITER = intbuf(9); NUMOUT = intbuf(10)

    logbuf(:) = (/KDIST, DELTAM, ACCELFLAG/)
    call mpi_bcast(logbuf, 3, MPI_logical, 0, MPI_comm_world, ierr)
    KDIST = logbuf(1); DELTAM = logbuf(2); ACCELFLAG = logbuf(3)

    realbuf(:) = (/SOLARFLUX, SOLARMU, SOLARAZ, sfc_obj%avgTemp, sfc_obj%avgAlb, SKYRAD,&
        &WAVENO(1:2), WAVELEN, SOLACC, SPLITACC, SHACC,&
        &MAX_TOTAL_MB, SPLITTING_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO/)
    call mpi_bcast(realbuf, 16, MPI_real, 0, MPI_comm_world, ierr)
    SOLARFLUX = realbuf(1); SOLARMU = realbuf(2); SOLARAZ = realbuf(3)
    sfc_obj%avgTemp = realbuf(4); sfc_obj%avgAlb = realbuf(5); SKYRAD = realbuf(6)
    WAVENO(1:2) = realbuf(7:8); WAVELEN = realbuf(9)
    SOLACC = realbuf(10); SPLITACC = realbuf(11); SHACC = realbuf(12)
    MAX_TOTAL_MB = realbuf(13); SPLITTING_FACTOR = realbuf(14)
    NUM_SH_TERM_FACTOR = realbuf(15); CELL_TO_POINT_RATIO = realbuf(16)

    char1buf(:) = (/SRCTYPE, GRIDTYPE, UNITS/)
    call mpi_bcast(char1buf, 3, MPI_character, 0, MPI_comm_world, ierr)
    SRCTYPE = char1buf(1); GRIDTYPE = char1buf(2); UNITS = char1buf(3)

    char256buf(:) = (/prp_obj%fname, sfc_obj%fname, ckd_obj%fname, fileInp_sav, fileOut_sav, sen_obj%fname, runName/)
    call mpi_bcast(char256buf, 7*256, MPI_character, 0, MPI_comm_world, ierr)
    prp_obj%fname = char256buf(1); sfc_obj%fname = char256buf(2); ckd_obj%fname = char256buf(3)
    fileInp_sav = char256buf(4); fileOut_sav = char256buf(5); sen_obj%fname = char256buf(6); runName = char256buf(7)

    call mpi_bcast(parms_Out, MAXPAR*NUMOUT, MPI_real, 0, MPI_comm_world, ierr)
    call mpi_bcast(OUTTYPES, NUMOUT*1, MPI_character, 0, MPI_comm_world, ierr)

    write (procStr, '(I4.4)') myProc_
    if (fileInp_sav(1:4) .ne. 'NONE') fileInp_sav = trim(fileInp_sav)//'.iproc-'//procStr
    if (fileOut_sav(1:4) .ne. 'NONE') fileOut_sav = trim(fileOut_sav)//'.iproc-'//procStr
#endif

    return

end subroutine shd_mpi__bcastInp

subroutine shd_mpi__bcastPrpSize(prp_obj)

    ! Broadcasts the property file size information from the master processor
    ! to the slave processors.

    type(shd_type_prp), intent(inout) :: prp_obj

#if HAS_MPI
    integer :: ierr, intbuf(7)
    real    :: realbuf(2)

    if (Nproc_ .le. 1) return

    intbuf(:) = (/prp_obj%Nx, prp_obj%Ny, prp_obj%Nz, prp_obj%Npha, prp_obj%maxNz, prp_obj%maxNleg, prp_obj%maxNptsLeg/)
    call mpi_bcast(intbuf, 7, MPI_integer, 0, MPI_comm_world, ierr)
    prp_obj%Nx      = intbuf(1); prp_obj%Ny    = intbuf(2); prp_obj%Nz = intbuf(3)
    prp_obj%Npha    = intbuf(4); prp_obj%maxNz = intbuf(5)
    prp_obj%maxNleg = intbuf(6); prp_obj%maxNptsLeg = intbuf(7)

    realbuf(:) = (/prp_obj%dx, prp_obj%dy/)
    call mpi_bcast(realbuf, 2, MPI_real, 0, MPI_comm_world, ierr)
    prp_obj%dx = realbuf(1); prp_obj%dy = realbuf(2)
#endif

    return

end subroutine shd_mpi__bcastPrpSize

subroutine shd_mpi__bcastPrp(prpAll_obj, temp_PrpAll, ext_PrpAll, ssa_PrpAll, leg_PrpAll, ipha_PrpAll, &
                                prp_obj, temp_Prp,    ext_Prp,    ssa_Prp,    leg_Prp,    ipha_Prp)

    ! Sends portions of the optical properties on the property grid to each
    ! of the processors.  Deals with the duplicating the grid points on the
    ! boundaries between processors.

    type(shd_type_prp), intent(inout) :: prpAll_obj
    real(kind=R__),  intent(in) :: temp_PrpAll(prpAll_obj%Nz, prpAll_obj%Ny, prpAll_obj%Nx), &
                                   ext_PrpAll( prpAll_obj%Nz, prpAll_obj%Ny, prpAll_obj%Nx), &
                                   ssa_PrpAll( prpAll_obj%Nz, prpAll_obj%Ny, prpAll_obj%Nx), &
                                   leg_PrpAll( prpAll_obj%maxNptsLeg)
    integer(kind=2), intent(in) :: ipha_PrpAll(prpAll_obj%Nz, prpAll_obj%Ny, prpAll_obj%Nx)

    type(shd_type_prp), intent(inout) :: prp_obj
    real(kind=R__),  intent(out) :: temp_Prp(prp_obj%maxNpts), &
                                    ext_Prp( prp_obj%maxNpts), &
                                    ssa_Prp( prp_obj%maxNpts), &
                                    leg_Prp( prp_obj%maxNptsLeg)
    integer(kind=2), intent(out) :: ipha_Prp(prp_obj%maxNpts)

#if HAS_MPI
    integer :: Sx, Ex, Nx, Sy, Ey, Ny, N, maxN, maxNx, maxNy
    integer :: iproc, irq, ierr, Nrq
    integer, allocatable :: ixyAll(:, :, :)
    real, allocatable :: rbuf(:, :, :), tempbuf(:, :), extbuf(:, :), ssabuf(:, :)
    integer(kind=2), allocatable :: ibuf(:, :, :), iphabuf(:, :)

    type(MPI_request), allocatable :: requests(:)
    type(MPI_status),  allocatable :: status(:)

    if (prp_obj%Npha .gt. 0) then
        Nrq = 8
    else
        Nrq = 7
    end if

    ! Get the ixyPrp_ arrays from all processors to tell where they are in property grid
    allocate(ixyAll(2, 2, 0 : Nproc_-1))
    call mpi_gather(ixyPrp_, 4, MPI_integer, ixyAll, 4, MPI_integer, 0, comm2D_, ierr)

    allocate(requests(Nrq*Nproc_ + Nrq), status(Nrq*Nproc_ + Nrq))

    irq = 0
    if (myProc_ .eq. 0) then

        ! Allocate the send buffers for each processor
        MAXN = 1; maxNx = 1; maxNy = 1
        do iproc = 0, Nproc_ - 1
            Nx = ixyAll(2, 1, iproc) - ixyAll(1, 1, iproc) + 1; maxNx = max(maxNx, Nx)
            Ny = ixyAll(2, 2, iproc) - ixyAll(1, 2, iproc) + 1; maxNy = max(maxNy, Ny)
            MAXN = max(MAXN, prpAll_obj%Nz*Ny*Nx)
        end do

        allocate(rbuf(prpAll_obj%Nz, maxNy, maxNx), ibuf(prpAll_obj%Nz, maxNy, maxNx))
        allocate(tempbuf(MAXN, 0:Nproc_ - 1),  extbuf(MAXN, 0:Nproc_ - 1))
        allocate( ssabuf(MAXN, 0:Nproc_ - 1), iphabuf(MAXN, 0:Nproc_ - 1))
        ! The master processor reformats the four property arrays (to deal with
        ! the wrap around at the full domain edge) and sends the arrays with MPI.
        do iproc = 0, Nproc_ - 1
            Sx = ixyAll(1, 1, iproc); Ex = ixyAll(2, 1, iproc)
            Sy = ixyAll(1, 2, iproc); Ey = ixyAll(2, 2, iproc)
            Nx = Ex - Sx + 1; Ny = Ey - Sy + 1; N = Nx*Ny*prpAll_obj%Nz
            rbuf(:, 1:Ny - 1, 1:Nx - 1) = temp_PrpAll(:, Sy:Ey - 1, Sx:Ex - 1)
            rbuf(:, 1:Ny - 1,       Nx) = temp_PrpAll(:, Sy:Ey - 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            rbuf(:,       Ny, 1:Nx - 1) = temp_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, Sx:Ex - 1)
            rbuf(:,       Ny,       Nx) = temp_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            tempbuf(1:N, iproc) = reshape(rbuf(1:prpAll_obj%Nz, 1:Ny, 1:Nx), (/N/))
            call mpi_isend(tempbuf(1:N, iproc), N, MPI_real, iproc, 1, comm2D_, requests(irq + 1), ierr)

            rbuf(:, 1:Ny - 1, 1:Nx - 1) = ext_PrpAll(:, Sy:Ey - 1, Sx:Ex - 1)
            rbuf(:, 1:Ny - 1,       Nx) = ext_PrpAll(:, Sy:Ey - 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            rbuf(:,       Ny, 1:Nx - 1) = ext_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, Sx:Ex - 1)
            rbuf(:,       Ny,       Nx) = ext_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            extbuf(1:N, iproc) = reshape(rbuf(1:prpAll_obj%Nz, 1:Ny, 1:Nx), (/N/))
            call mpi_isend(extbuf(1:N, iproc), N, MPI_real, iproc, 2, comm2D_, requests(irq + 2), ierr)

            rbuf(:, 1:Ny - 1, 1:Nx - 1) = ssa_PrpAll(:, Sy:Ey - 1, Sx:Ex - 1)
            rbuf(:, 1:Ny - 1,       Nx) = ssa_PrpAll(:, Sy:Ey - 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            rbuf(:,       Ny, 1:Nx - 1) = ssa_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, Sx:Ex - 1)
            rbuf(:,       Ny,       Nx) = ssa_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            ssabuf(1:N, iproc) = reshape(rbuf(1:prpAll_obj%Nz, 1:Ny, 1:Nx), (/N/))
            call mpi_isend(ssabuf(1:N, iproc), N, MPI_real, iproc, 3, comm2D_, requests(irq + 3), ierr)

            ibuf(:, 1:Ny - 1, 1:Nx - 1) = ipha_PrpAll(:, Sy:Ey - 1, Sx:Ex - 1)
            ibuf(:, 1:Ny - 1,       Nx) = ipha_PrpAll(:, Sy:Ey - 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            ibuf(:,       Ny, 1:Nx - 1) = ipha_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, Sx:Ex - 1)
            ibuf(:,       Ny,       Nx) = ipha_PrpAll(:, mod(Ey - 1, prpAll_obj%Ny) + 1, mod(Ex - 1, prpAll_obj%Nx) + 1)
            iphabuf(1:N, iproc) = reshape(ibuf(1:prpAll_obj%Nz, 1:Ny, 1:Nx), (/N/))
            call mpi_isend(iphabuf(1:N, iproc), N, MPI_integer2, iproc, 4, comm2D_, requests(irq + 4), ierr)

            call mpi_isend(prpAll_obj%Nleg   , 1, MPI_integer, iproc, 5, comm2D_, requests(irq + 5), ierr)
            call mpi_isend(prpAll_obj%maxAsym, 1, MPI_real   , iproc, 6, comm2D_, requests(irq + 6), ierr)

            call mpi_isend(prpAll_obj%z, prpAll_obj%Nz, MPI_real, iproc, 7, comm2D_, requests(irq + 7), ierr)

            if (prpAll_obj%Npha .gt. 0) then
                call mpi_isend(leg_PrpAll, prpAll_obj%maxNptsLeg, MPI_real, iproc, 8, comm2D_, requests(irq + 8), ierr)
            end if

            irq = irq + Nrq

        end do


    end if

    N = prp_obj%Nx*prp_obj%Ny*prp_obj%Nz
    call mpi_irecv(temp_Prp, N, MPI_real    , 0, 1, comm2D_, requests(irq+1), ierr)
    call mpi_irecv(ext_Prp , N, MPI_real    , 0, 2, comm2D_, requests(irq+2), ierr)
    call mpi_irecv(ssa_Prp , N, MPI_real    , 0, 3, comm2D_, requests(irq+3), ierr)
    call mpi_irecv(ipha_Prp, N, MPI_integer2, 0, 4, comm2D_, requests(irq+4), ierr)

    call mpi_irecv(prp_obj%Nleg   , 1, MPI_integer, 0, 5, comm2D_, requests(irq+5), ierr)
    call mpi_irecv(prp_obj%maxAsym, 1, MPI_real   , 0, 6, comm2D_, requests(irq+6), ierr)

    call mpi_irecv(prp_obj%z, prp_obj%Nz, MPI_real, 0, 7, comm2D_, requests(irq+7), ierr)

    if (prp_obj%Npha .gt. 0) then
        call mpi_irecv(leg_Prp, prp_obj%maxNptsLeg, MPI_real, 0, 8, comm2D_, requests(irq+8), ierr)
    end if

    irq = irq + Nrq
    call mpi_waitall(irq, requests, status, ierr)

    if (myProc_ .eq. 0) deallocate(tempbuf, extbuf, ssabuf, iphabuf, rbuf, ibuf)
    deallocate(requests, status, ixyAll)
#endif

    return

end subroutine shd_mpi__bcastPrp

subroutine shd_mpi__bcastSfcSize(sfc_obj)

    ! Broadcasts the surface array sizes from the master process to the others.

    type(shd_type_sfc), intent(inout) :: sfc_obj

#if HAS_MPI
    integer :: ierr, intbuf(4)
    real :: realbuf(2)

    if (Nproc_ .le. 1) return

    intbuf(1:4) = (/sfc_obj%Nx, sfc_obj%Ny, sfc_obj%maxNpts, sfc_obj%maxNpar/)
    call mpi_bcast(intbuf, 4, MPI_integer, 0, comm2D_, ierr)
    sfc_obj%Nx = intbuf(1); sfc_obj%Ny = intbuf(2); sfc_obj%maxNpts = intbuf(3); sfc_obj%maxNpar = intbuf(4)

    realbuf(1:2) = (/sfc_obj%dx, sfc_obj%dy/)
    call mpi_bcast(realbuf, 2, MPI_real, 0, comm2D_, ierr)
    sfc_obj%dx = realbuf(1); sfc_obj%dy = realbuf(2)
#endif

    return

end subroutine shd_mpi__bcastSfcSize

subroutine shd_mpi__bcastSfc(sfcAll_obj, brdf_SfcAll, &
                             sfc_obj,    brdf_Sfc)

    ! Broadcasts the surface file parameters from the master process to the others.

    type(shd_type_sfc), intent(inout) :: sfcAll_obj
    real(kind=R__),  intent(in) :: brdf_SfcAll(sfcAll_obj%maxNpar, sfcAll_obj%Ny+1, sfcAll_obj%Nx+1)

    type(shd_type_sfc), intent(inout) :: sfc_obj
    real(kind=R__), intent(out) :: brdf_Sfc(sfc_obj%maxNpar*sfc_obj%maxNpts)

#if HAS_MPI
    integer :: Sx, Ex, Nx, Sy, Ey, Ny, N, maxN, maxNx, maxNy
    integer :: iproc, irq, ierr, Nrq
    integer, allocatable :: ixyAll(:, :, :)
    real, allocatable :: rbuf(:, :, :), brdfbuf(:, :)

    type(MPI_request), allocatable :: requests(:)
    type(MPI_status),  allocatable :: status(:)

    Nrq = 5

    ! Get the ixySfc_ arrays from all processors to tell where they are in property grid
    allocate(ixyAll(2, 2, 0 : Nproc_-1))
    call mpi_gather(ixySfc_, 4, MPI_integer, ixyAll, 4, MPI_integer, 0, comm2D_, ierr)

    allocate(requests(Nrq*Nproc_ + Nrq), status(Nrq*Nproc_ + Nrq))

    irq = 0
    if (myProc_ .eq. 0) then

        ! Allocate the send buffers for each processor
        maxN = 1; maxNx = 1; maxNy = 1
        do iproc = 0, Nproc_ - 1
            Nx = ixyAll(2, 1, iproc) - ixyAll(1, 1, iproc) + 1; maxNx = max(maxNx, Nx)
            Ny = ixyAll(2, 2, iproc) - ixyAll(1, 2, iproc) + 1; maxNy = max(maxNy, Ny)
            maxN = max(maxN, sfcAll_obj%maxNpar*Ny*Nx)
        end do

        allocate(rbuf(sfcAll_obj%maxNpar, maxNy, maxNx))
        allocate(brdfbuf(maxN, 0:Nproc_ - 1))

        ! The master processor reformats the four property arrays (to deal with
        ! the wrap around at the full domain edge) and sends the arrays with MPI.
        do iproc = 0, Nproc_ - 1

            Sx = ixyAll(1, 1, iproc); Ex = ixyAll(2, 1, iproc)
            Sy = ixyAll(1, 2, iproc); Ey = ixyAll(2, 2, iproc)
            Nx = Ex - Sx + 1; Ny = Ey - Sy + 1; N = Nx*Ny*sfcAll_obj%maxNpar

            write (6, *) iproc, Nx, Sx, Ex
            write (6, *) iproc, Ny, Sy, Ey
            write (6, *) Nx*Ny, sfc_obj%maxNpts
            write (6, *)

            rbuf(:, 1:Ny, 1:Nx) = brdf_SfcAll(:, Sy:Ey, Sx:Ex)

            brdfbuf(1:N, iproc) = reshape(rbuf(1:sfcAll_obj%maxNpar, 1:Ny, 1:Nx), (/N/))

            call mpi_isend(brdfbuf(1:N, iproc), N, MPI_real,      iproc, 1, comm2D_, requests(irq+1), ierr)
            call mpi_isend(sfcAll_obj%type,     2, MPI_character, iproc, 2, comm2D_, requests(irq+2), ierr)
            call mpi_isend(sfcAll_obj%Npar,     1, MPI_integer,   iproc, 3, comm2D_, requests(irq+3), ierr)
            call mpi_isend(sfcAll_obj%avgTemp,  1, MPI_real,      iproc, 4, comm2D_, requests(irq+4), ierr)
            call mpi_isend(sfcAll_obj%avgAlb,   1, MPI_real,      iproc, 5, comm2D_, requests(irq+5), ierr)

            irq = irq + Nrq

        end do

    end if

    N = sfc_obj%maxNpar*sfc_obj%maxNpts
    call mpi_irecv(brdf_Sfc,        N, MPI_real,      0, 1, comm2D_, requests(irq+1), ierr)
    call mpi_irecv(sfc_obj%type,    2, MPI_character, 0, 2, comm2D_, requests(irq+2), ierr)
    call mpi_irecv(sfc_obj%Npar,    1, MPI_integer,   0, 3, comm2D_, requests(irq+3), ierr)
    call mpi_irecv(sfc_obj%avgTemp, 1, MPI_real,      0, 4, comm2D_, requests(irq+4), ierr)
    call mpi_irecv(sfc_obj%avgAlb,  1, MPI_real,      0, 5, comm2D_, requests(irq+5), ierr)

    irq = irq + Nrq
    call mpi_waitall(irq, requests, status, ierr)

    if (myProc_ .eq. 0) deallocate(brdfbuf, rbuf)
    deallocate(requests, status, ixyAll)
#endif

    return

end subroutine shd_mpi__bcastSfc

subroutine shd_mpi__bcastCkdSize(ckd_obj)

    ! Broadcasts the k-distribution array sizes from the master process to the others.

    type(shd_type_ckd), intent(inout) :: ckd_obj

#if HAS_MPI
    integer :: ierr, intbuf(2)

    if (Nproc_ .le. 1) return

    intbuf(1:2) = (/ckd_obj%Ng, ckd_obj%Nz/)
    call mpi_bcast(intbuf, 2, MPI_integer, 0, comm2D_, ierr)
    ckd_obj%Ng = intbuf(1); ckd_obj%Nz = intbuf(2)
#endif

    return

end subroutine shd_mpi__bcastCkdSize

subroutine shd_mpi__bcastCkd(ckd_obj, solar_Ckd, wgt_Ckd, kabs_Ckd)

    ! Broadcasts the k-distribution file parameters from the master process
    ! to the other processes.

    type(shd_type_ckd), intent(inout) :: ckd_obj
    real(kind=R__), intent(inout) :: solar_Ckd(ckd_obj%Ng), &
                                       wgt_Ckd(ckd_obj%Ng), &
                                      kabs_Ckd(ckd_obj%Nz, ckd_obj%Ng)

#if HAS_MPI
    integer :: ierr
    real, allocatable :: realbuf(:)

    if (Nproc_ .le. 1) return

    allocate(realbuf(2*ckd_obj%Ng + ckd_obj%Nz))

    realbuf(1              : ckd_obj%Ng             ) = solar_Ckd(:)
    realbuf(ckd_obj%Ng+1   : 2*ckd_obj%Ng           ) = wgt_Ckd(:)
    realbuf(2*ckd_obj%Ng+1 : 2*ckd_obj%Ng+ckd_obj%Nz) = ckd_obj%z(:)

    call mpi_bcast(realbuf, 2*ckd_obj%Ng+ckd_obj%Nz, MPI_real, 0, comm2D_, ierr)

    solar_Ckd(:) = realbuf(             1 :  ckd_obj%Ng)
    wgt_Ckd(:)   = realbuf(  ckd_obj%Ng+1 : 2*ckd_obj%Ng)
    ckd_obj%z(:) = realbuf(2*ckd_obj%Ng+1 : 2*ckd_obj%Ng+ckd_obj%Nz)

    deallocate(realbuf)

    call mpi_bcast(kabs_Ckd, ckd_obj%Nz*ckd_obj%Ng, MPI_real, 0, comm2D_, ierr)
#endif

    return

end subroutine shd_mpi__bcastCkd

subroutine shd_mpi__bcastSenSize(sen_obj)

    ! Broadcasts the user-specified sensor data array sizes from
    ! the master process to the others.

    type(shd_type_sen), intent(inout) :: sen_obj

#if HAS_MPI
    integer :: ierr, intbuf(2)

    if (Nproc_ .le. 1) return

    intbuf(1:2) = (/sen_obj%Npts, sen_obj%Npar/)
    call mpi_bcast(intbuf, 2, MPI_integer, 0, comm2D_, ierr)
    sen_obj%Npts = intbuf(1); sen_obj%Npar = intbuf(2)
#endif

    return

end subroutine shd_mpi__bcastSenSize

subroutine shd_mpi__bcastSen(sen_obj, view_Sen)

    ! Broadcasts the user-specified sensor data array from
    ! the master process to the others.

    type(shd_type_sen), intent(inout) :: sen_obj
    real(kind=R__),     intent(inout) :: view_Sen(sen_obj%Npar, sen_obj%Npts)

#if HAS_MPI
    integer :: ierr, Ndata

    if (Nproc_ .le. 1) return

    Ndata = sen_obj%Npts*sen_obj%Npar
    call mpi_bcast(view_Sen, Ndata, MPI_real, 0, comm2D_, ierr)
#endif

    return

end subroutine shd_mpi__bcastSen

subroutine shd_mpi__bcastMem(MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, runName)

    ! The master process reads the runName//"_mem_params.inp" file (if there
    ! is one) and broadcasts the memory parameters to all the other processors.

    real, intent(inout) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real, intent(inout) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(len=*), intent(in) :: runName

#if HAS_MPI
    integer :: Nproc, ierr, i, j
    real, allocatable :: MemParam(:, :)

    allocate(MemParam(4, 0:Nproc_ - 1))
    MemParam(1, :) = MAX_TOTAL_MB; MemParam(2, :) = ADAPT_GRID_FACTOR
    MemParam(3, :) = NUM_SH_TERM_FACTOR; MemParam(4, :) = CELL_TO_POINT_RATIO

    if (myProc_ .eq. 0) then
        ! See if we can read the memory parameter file and number of processors matches
        open (UNIT=11, FILE=trim(runName)//'_mem_params.inp', STATUS='OLD', IOSTAT=ierr)
        if (ierr .eq. 0) then
            read (11, *)
            read (11, *) Nproc
            if (Nproc .ne. Nproc_) then
                write (6, '(A,A,A)') 'Memory parameter file ', trim(runName),&
                    &'_mem_params.inp is for a different number of processors.'
            else
                do i = 0, Nproc - 1
                    read (11, *) MemParam(:, i)
                end do
            end if
            close (11)
        end if
    end if

    call mpi_bcast(MemParam, 4*Nproc_, MPI_real, 0, comm2D_, ierr)
    MAX_TOTAL_MB = MemParam(1, myProc_)
    ADAPT_GRID_FACTOR = MemParam(2, myProc_)
    NUM_SH_TERM_FACTOR = MemParam(3, myProc_)
    CELL_TO_POINT_RATIO = MemParam(4, myProc_)
    deallocate(MemParam)
    if (verbose_ .ge. 5) write (6, '(A)') 'Message [mpi__bcastMemParam]:'
    if (verbose_ .ge. 5) write (6, '(2X,A,F0.2,"MB",3(2X,F0.3))') 'Memory parameters used: ', MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
#endif

    return

end subroutine shd_mpi__bcastMem

real function shd_mpi__sumCpuTime(cpuTimes)

    ! Sums the CPU time used across all processors

    real :: cpuTimes
    integer :: ierr

    shd_mpi__sumCpuTime = cpuTimes

#if HAS_MPI
    if (Nproc_ .gt. 1) then
        call mpi_reduce(cpuTimes, shd_mpi__sumCpuTime, 1, MPI_real, MPI_sum, 0, comm2D_, ierr)
    end if
#endif

    return

end function shd_mpi__sumCpuTime

subroutine shd_mpi__gatherOut(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
    &NSHOUT,&
    &SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
    &ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,&
    &NSH, NCELLSTOT, NPTSTOT, NSHTOT,&
    &rtmAll_obj, rtm_obj, cell_obj)

    ! Gathers the base grid output data in SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,
    ! and SUMSHTERMS from all the processors and assembles in the master
    ! processor output arrays ALLFLUXES, ALLFLUXDIV, ALLSHTERMS.  Gathers the
    ! ixyOut_ array from each processor to determine where each subdomain is.
    ! Also gathers the full base grid position arrays (rtmAll_obj%x/rtmAll_obj%y) and
    ! sums cell_obj%Npts,rtm_obj%Npts,NSH over all processors.

    type(shd_type_rtm),  intent(inout) :: rtmAll_obj, rtm_obj
    type(shd_type_cell), intent(in)    :: cell_obj

    logical, intent(in) :: FLUX_OUT, FLUXDIV_OUT, SH_OUT
    integer, intent(in) :: NSHOUT
    real, intent(in) :: SUMFLUXES(2, rtm_obj%Npts0), SUMDIRFLUX(rtm_obj%Npts0)
    real, intent(in) :: SUMFLUXDIV(rtm_obj%Npts0), SUMSHTERMS(NSHOUT, rtm_obj%Npts0)
    real, intent(out) :: ALLFLUXES(3, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx), ALLFLUXDIV(rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx)
    real, intent(out) :: ALLSHTERMS(NSHOUT, rtmAll_obj%Nz, rtmAll_obj%Ny, rtmAll_obj%Nx)
    integer, intent(in) :: NSH
    integer, intent(out) :: NCELLSTOT, NPTSTOT, NSHTOT

#if HAS_MPI
    integer :: IX, IY, Nx, Ny, NX1, NY1, maxN, I, send3int(3), recv3int(3)
    integer :: iproc, ierr, irq
    integer, allocatable :: ixyAll(:, :, :)
    real, allocatable :: sendbuf(:, :), recvbuf(:, :, :), buf(:, :, :, :)

    type(MPI_request), allocatable :: requests(:)
    type(MPI_status), allocatable :: status(:)

    allocate(ixyAll(2, 2, 0:Nproc_ - 1))
    call mpi_gather(ixyOut_, 4, MPI_integer, ixyAll, 4, MPI_integer, 0, comm2D_, ierr)

    allocate(requests(Nproc_), status(Nproc_))

    if (FLUX_OUT) then
        irq = 0
        if (myProc_ .gt. 0) then
            allocate(sendbuf(3, rtm_obj%Npts0))
            sendbuf(1:2, :) = SUMFLUXES(:, :); sendbuf(3, :) = SUMDIRFLUX(:)
            irq = irq + 1
            call mpi_isend(sendbuf, 3*rtm_obj%Npts0, MPI_real, 0, 1, comm2D_,&
                &requests(irq), ierr)
        else
            maxN = (maxval(ixyAll(2, 1, :)) + 1)*(maxval(ixyAll(2, 2, :)) + 1)*rtmAll_obj%Nz
            allocate(recvbuf(3, maxN, 1:Nproc_ - 1))
            do iproc = 1, Nproc_ - 1
                irq = irq + 1
                call mpi_irecv(recvbuf(:, :, iproc), 3*maxN, MPI_real, iproc, 1,&
                    &comm2D_, requests(irq), ierr)
            end do
        end if
        call mpi_waitall(irq, requests, status, ierr)

        if (myProc_ .eq. 0) then
            IX = ixyAll(1, 1, 0); Nx = ixyAll(2, 1, 0)
            IY = ixyAll(1, 2, 0); Ny = ixyAll(2, 2, 0)
            NX1 = Nx + 1; NY1 = Ny + 1
            if (is1DX_ .and. .not. isMpiX_) NX1 = Nx
            if (is1DYx_ .and. .not. isMpiY_) NY1 = Ny
            allocate(buf(3, rtmAll_obj%Nz, NY1, NX1))
            buf(1:2, :, :, :) = reshape(SUMFLUXES(1:2, :), (/2, rtmAll_obj%Nz, NY1, NX1/))
            buf(3, :, :, :) = reshape(SUMDIRFLUX(:), (/rtmAll_obj%Nz, NY1, NX1/))
            ALLFLUXES(1:3, 1:rtmAll_obj%Nz, IY:IY + Ny - 1, IX:IX + Nx - 1) = buf(:, 1:rtmAll_obj%Nz, 1:Ny, 1:Nx)
            deallocate(buf)
            do iproc = 1, Nproc_ - 1
                IX = ixyAll(1, 1, iproc); Nx = ixyAll(2, 1, iproc)
                IY = ixyAll(1, 2, iproc); Ny = ixyAll(2, 2, iproc)
                NX1 = Nx + 1; NY1 = Ny + 1
                if (is1DX_ .and. .not. isMpiX_) NX1 = Nx
                if (is1DYx_ .and. .not. isMpiY_) NY1 = Ny
                allocate(buf(3, rtmAll_obj%Nz, NY1, NX1))
                buf(:, :, :, :) = reshape(recvbuf(:, 1:rtmAll_obj%Nz*NY1*NX1, iproc), (/3, rtmAll_obj%Nz, NY1, NX1/))
                ALLFLUXES(1:3, 1:rtmAll_obj%Nz, IY:IY + Ny - 1, IX:IX + Nx - 1) = buf(:, 1:rtmAll_obj%Nz, 1:Ny, 1:Nx)
                deallocate(buf)
            end do
            deallocate(recvbuf)
        else
            deallocate(sendbuf)
        end if
    end if

    if (FLUXDIV_OUT) then
        irq = 0
        if (myProc_ .gt. 0) then
            irq = irq + 1
            call mpi_isend(SUMFLUXDIV, rtm_obj%Npts0, MPI_real, 0, 2, comm2D_,&
                &requests(irq), ierr)
        else
            maxN = (maxval(ixyAll(2, 1, :)) + 1)*(maxval(ixyAll(2, 2, :)) + 1)*rtmAll_obj%Nz
            allocate(recvbuf(1, maxN, 1:Nproc_ - 1))
            do iproc = 1, Nproc_ - 1
                irq = irq + 1
                call mpi_irecv(recvbuf(:, :, iproc), maxN, MPI_real, iproc, 2,&
                    &comm2D_, requests(irq), ierr)
            end do
        end if
        call mpi_waitall(irq, requests, status, ierr)

        if (myProc_ .eq. 0) then
            IX = ixyAll(1, 1, 0); Nx = ixyAll(2, 1, 0)
            IY = ixyAll(1, 2, 0); Ny = ixyAll(2, 2, 0)
            NX1 = Nx + 1; NY1 = Ny + 1
            if (is1DX_ .and. .not. isMpiX_) NX1 = Nx
            if (is1DYx_ .and. .not. isMpiY_) NY1 = Ny
            allocate(buf(1, rtmAll_obj%Nz, NY1, NX1))
            buf(1, :, :, :) = reshape(SUMFLUXDIV(:), (/rtmAll_obj%Nz, NY1, NX1/))
            ALLFLUXDIV(1:rtmAll_obj%Nz, IY:IY + Ny - 1, IX:IX + Nx - 1) = buf(1, 1:rtmAll_obj%Nz, 1:Ny, 1:Nx)
            deallocate(buf)
            do iproc = 1, Nproc_ - 1
                IX = ixyAll(1, 1, iproc); Nx = ixyAll(2, 1, iproc)
                IY = ixyAll(1, 2, iproc); Ny = ixyAll(2, 2, iproc)
                NX1 = Nx + 1; NY1 = Ny + 1
                if (is1DX_ .and. .not. isMpiX_) NX1 = Nx
                if (is1DYx_ .and. .not. isMpiY_) NY1 = Ny
                allocate(buf(1, rtmAll_obj%Nz, NY1, NX1))
                buf(:, :, :, :) = reshape(recvbuf(:, 1:rtmAll_obj%Nz*NY1*NX1, iproc), (/1, rtmAll_obj%Nz, NY1, NX1/))
                ALLFLUXDIV(1:rtmAll_obj%Nz, IY:IY + Ny - 1, IX:IX + Nx - 1) = buf(1, 1:rtmAll_obj%Nz, 1:Ny, 1:Nx)
                deallocate(buf)
            end do
            deallocate(recvbuf)
        end if
    end if

    if (SH_OUT) then
        irq = 0
        if (myProc_ .gt. 0) then
            irq = irq + 1
            call mpi_isend(SUMSHTERMS, NSHOUT*rtm_obj%Npts0, MPI_real, 0, 3, comm2D_,&
                &requests(irq), ierr)
        else
            maxN = (maxval(ixyAll(2, 1, :)) + 1)*(maxval(ixyAll(2, 2, :)) + 1)*rtmAll_obj%Nz
            allocate(recvbuf(NSHOUT, maxN, 1:Nproc_ - 1))
            do iproc = 1, Nproc_ - 1
                irq = irq + 1
                call mpi_irecv(recvbuf(:, :, iproc), NSHOUT*maxN, MPI_real, iproc, 3,&
                    &comm2D_, requests(irq), ierr)
            end do
        end if
        call mpi_waitall(irq, requests, status, ierr)

        if (myProc_ .eq. 0) then
            IX = ixyAll(1, 1, 0); Nx = ixyAll(2, 1, 0)
            IY = ixyAll(1, 2, 0); Ny = ixyAll(2, 2, 0)
            NX1 = Nx + 1; NY1 = Ny + 1
            if (is1DX_ .and. .not. isMpiX_) NX1 = Nx
            if (is1DYx_ .and. .not. isMpiY_) NY1 = Ny
            allocate(buf(NSHOUT, rtmAll_obj%Nz, NY1, NX1))
            buf(:, :, :, :) = reshape(SUMSHTERMS(:, :), (/NSHOUT, rtmAll_obj%Nz, NY1, NX1/))
            ALLSHTERMS(:, 1:rtmAll_obj%Nz, IY:IY + Ny - 1, IX:IX + Nx - 1) = buf(:, 1:rtmAll_obj%Nz, 1:Ny, 1:Nx)
            deallocate(buf)
            do iproc = 1, Nproc_ - 1
                IX = ixyAll(1, 1, iproc); Nx = ixyAll(2, 1, iproc)
                IY = ixyAll(1, 2, iproc); Ny = ixyAll(2, 2, iproc)
                NX1 = Nx + 1; NY1 = Ny + 1
                if (is1DX_ .and. .not. isMpiX_) NX1 = Nx
                if (is1DYx_ .and. .not. isMpiY_) NY1 = Ny
                allocate(buf(NSHOUT, rtmAll_obj%Nz, NY1, NX1))
                buf(:, :, :, :) = reshape(recvbuf(:, 1:rtmAll_obj%Nz*NY1*NX1, iproc), (/NSHOUT, rtmAll_obj%Nz, NY1, NX1/))
                ALLSHTERMS(:, 1:rtmAll_obj%Nz, IY:IY + Ny - 1, IX:IX + Nx - 1) = buf(:, 1:rtmAll_obj%Nz, 1:Ny, 1:Nx)
                deallocate(buf)
            end do
            deallocate(recvbuf)
        end if
    end if

    ! Make the X and Y grids for the full domain
    if (myProc_ .eq. 0) then
        rtmAll_obj%x(1:rtmAll_obj%Nx) = (/(prp_obj_%xDomain*(I - 1)/rtmAll_obj%Nx, I=1, rtmAll_obj%Nx)/)
        rtmAll_obj%y(1:rtmAll_obj%Ny) = (/(prp_obj_%yDomain*(I - 1)/rtmAll_obj%Ny, I=1, rtmAll_obj%Ny)/)
    end if

    ! Get the sum of the cell_obj%Npts, rtm_obj%Npts, and NSH over the processors
    send3int = (/cell_obj%Npts, rtm_obj%Npts, NSH/)
    call mpi_reduce(send3int, recv3int, 3, MPI_integer, MPI_sum, 0, comm2D_, ierr)
    NCELLSTOT = recv3int(1); NPTSTOT = recv3int(2); NSHTOT = recv3int(3)

    deallocate(ixyAll, requests, status)
#endif

    return

end subroutine shd_mpi__gatherOut

subroutine shd_mpi__findMaxSsa(maxSsa)

    real, intent(inout) :: maxSsa

#if HAS_MPI
    integer :: ierr
    real    :: sendrbuf

    if (Nproc_ .gt. 1) then
        ! Call mpi_allreduce to find the max maxSsa over all processors
        sendrbuf = maxSsa
        call mpi_allreduce(sendrbuf, maxSsa, 1, MPI_real, MPI_max, comm2D_, ierr)
        if (ierr .ne. MPI_success) call shd_mpi__abort('Error [mpi__findSsaMax]: mpi_allreduce error.')
    end if
#endif

    return

end subroutine shd_mpi__findMaxSsa

subroutine shd_mpi__unifySplit(DOSPLIT, STARTSPLITACC)

    logical, intent(inout) :: DOSPLIT
    real, intent(inout) :: STARTSPLITACC

#if HAS_MPI
    integer :: ierr(2)
    logical :: sendlbuf
    real    :: sendrbuf

    if (Nproc_ .gt. 1) then
        ! Call mpi_allreduce to
        sendrbuf = STARTSPLITACC
        call mpi_allreduce(sendrbuf, STARTSPLITACC, 1, MPI_real, MPI_max, comm2D_, ierr(1))
        sendlbuf = DOSPLIT
        call mpi_allreduce(sendlbuf, DOSPLIT, 1, MPI_logical, MPI_lor, comm2D_, ierr(2))
        if (any(ierr(:) .ne. MPI_success)) call shd_mpi__abort('Error [mpi__unifySplit]: mpi_allreduce error.')
    end if
#endif

    return

end subroutine shd_mpi__unifySplit

subroutine shd_mpi__findMaxSplitCrit(SPLITCRIT)

    real, intent(inout) :: SPLITCRIT

#if HAS_MPI
    integer :: ierr
    real    :: sendrbuf

    if (Nproc_ .gt. 1) then
        ! Call mpi_allreduce to find the max SPLITCRIT over all processors
        sendrbuf = SPLITCRIT
        call mpi_allreduce(sendrbuf, SPLITCRIT, 1, MPI_real, MPI_max, comm2D_, ierr)
        if (ierr .ne. MPI_success) call shd_mpi__abort('Error [mpi__findTotSplitCritMax]: mpi_allreduce error.')
    end if
#endif

    return

end subroutine shd_mpi__findMaxSplitCrit

subroutine shd_mpi__makeDirect(SPT,&
    &DELTAM, ML,&
    &SOLARFLUX, SOLARMU, SOLARAZ,&
    &DIRFLUX)

    ! Makes the direct beam solar flux for the rtm_obj_%Npts gridpoints started with
    ! index SPT having X/Y/Z positions in xyz_Rtm_.
    ! DIRFLUX is set to F*exp(-tau_sun).

    integer, intent(in) :: SPT, ML
    logical, intent(in) :: DELTAM
    real, intent(in) :: SOLARFLUX, SOLARMU, SOLARAZ
    real, intent(out) :: DIRFLUX(rtm_obj_%Npts)

#if HAS_MPI
    logical :: VALIDBEAM
    real    :: UNIFZLEV, UZLEV, DIRPATH
    integer :: MAXNBEAM, Lupx, Lupy, Ldownx, Ldowny, Nbeam
    integer :: SIDE, NINVALID(0:4)
    integer :: IP, I, J, L
    integer :: ierr, irq, iproc
    integer :: intsendbuf(2), intrecvbuf(2)
    logical :: VALIDRAD, AllDone
    real    :: XO, YO, ZO, X, Y, Z, XMAX, YMAX, EPS, DIRJUNK
    integer, allocatable :: Ndone(:)
    integer, allocatable :: BndDirIntInfo(:, :, :)
    real, allocatable    :: BndDirRealInfo(:, :, :)
    integer, allocatable :: DirPtrOut(:, :), DirPtrIn(:)
    real, allocatable    :: DirAllOut(:, :), DirAllIn(:)

    type(MPI_request) :: requests(12)
    type(MPI_status) :: status(12)

    MAXNBEAM = 2*rtm_obj_%Npts
    allocate(BndDirRealInfo(4, MAXNBEAM, 4), BndDirIntInfo(2, MAXNBEAM, 4))
    allocate(DirAllOut(MAXNBEAM, 0:Nproc_ - 1), DirPtrOut(MAXNBEAM, 0:Nproc_ - 1))
    allocate(DirAllIn(MAXNBEAM), DirPtrIn(MAXNBEAM), Ndone(0:Nproc_ - 1))

    EPS = 1.0e-5*max(rtm_obj_%x(2) - rtm_obj_%x(1), rtm_obj_%y(2) - rtm_obj_%y(1))

    if (SPT .eq. 1) then
        ! Initialize the direct beam calculation
        call shd_rad__computeDirect(1, 0.0, 0.0, 0.0, rtm_obj_%flag_bc, rtm_obj_%flag_ip,&
            &DELTAM, ML, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
        ! Find the max lowest uniform Z level over all processors and then
        !   tell shd_rad__computeDirect about it
        UZLEV = UNIFZLEV
        call mpi_allreduce(UZLEV, UNIFZLEV, 1, MPI_real, MPI_max, comm2D_, ierr)
        call shd_rad__computeDirect(2, 0.0, 0.0, 0.0, rtm_obj_%flag_bc, rtm_obj_%flag_ip,&
            &DELTAM, ML, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end if

    Lupx = 0; Lupy = 0; Ldownx = 0; Ldowny = 0
    if (isMpiX_) then
        if (cos(SOLARAZ) .gt. 1.0e-5) then
            Ldownx = 2; Lupx = 1
        else
            Ldownx = 1; Lupx = 2
        end if
        XMAX = rtm_obj_%x(rtm_obj_%Nx)
    else
        XMAX = rtm_obj_%x(rtm_obj_%Nx + 1)
    end if
    if (isMpiY_) then
        if (sin(SOLARAZ) .gt. 1.0e-5) then
            Ldowny = 4; Lupy = 3
        else
            Ldowny = 3; Lupy = 4
        end if
        YMAX = rtm_obj_%y(rtm_obj_%Ny)
    else
        YMAX = rtm_obj_%y(rtm_obj_%Ny + 1)
    end if


    ! Do the initial direct beam tracing starting at the base grid points
    NINVALID(:) = 0
    do IP = SPT, rtm_obj_%Npts
        DIRPATH = 0.0
        call shd_rad__computeDirect(0, xyz_Rtm_(1, IP), xyz_Rtm_(2, IP), xyz_Rtm_(3, IP),&
            &rtm_obj_%flag_bc, rtm_obj_%flag_ip, DELTAM, ML,&
            &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(IP),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
        ! If we don't have a good direct beam then save the information
        !   needed to pass to two upstream processors
        if (.not. VALIDBEAM) then
            if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) call shd_mpi__abort('Error [mpi__makeDirect]: Bad SIDE 1.')
            NINVALID(SIDE) = NINVALID(SIDE) + 1
            if (NINVALID(SIDE) .gt. MAXNBEAM) call shd_mpi__abort('Error [mpi__makeDirect]: NINVALID(SIDE) exceeds MAXNBEAM 1.')
            BndDirRealInfo(:, NINVALID(SIDE), SIDE) = (/XO, YO, ZO, DIRPATH/)
            BndDirIntInfo(:, NINVALID(SIDE), SIDE) = (/myProc_, IP/)
        end if
    end do

    ! Cycle over the passing of info between processors and doing more
    !  partial direct beam integrations until all rays have been finished
    Ndone(:) = 0
    AllDone = .false.
    do while (.not. AllDone)
        ! Do the send and receives of the number of incomplete rays and
        !   the real and integer info needed about the radiance integrations
        irq = 0
        if (isMpiX_) then
            call mpi_isend(NINVALID(Lupx), 1, MPI_integer,&
                &iProcNeigh_(Lupx), 1, comm2D_, requests(irq + 1), ierr)
            call mpi_isend(BndDirRealInfo(:, :, Lupx), 4*NINVALID(Lupx), MPI_real,&
                &iProcNeigh_(Lupx), 2, comm2D_, requests(irq + 2), ierr)
            call mpi_isend(BndDirIntInfo(:, :, Lupx), 2*NINVALID(Lupx), MPI_integer,&
                &iProcNeigh_(Lupx), 3, comm2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        if (isMpiY_) then
            call mpi_isend(NINVALID(Lupy), 1, MPI_integer,&
                &iProcNeigh_(Lupy), 4, comm2D_, requests(irq + 1), ierr)
            call mpi_isend(BndDirRealInfo(:, :, Lupy), 4*NINVALID(Lupy), MPI_real,&
                &iProcNeigh_(Lupy), 5, comm2D_, requests(irq + 2), ierr)
            call mpi_isend(BndDirIntInfo(:, :, Lupy), 2*NINVALID(Lupy), MPI_integer,&
                &iProcNeigh_(Lupy), 6, comm2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        if (isMpiX_) then
            call mpi_irecv(NINVALID(Ldownx), 1, MPI_integer,&
                &iProcNeigh_(Ldownx), 1, comm2D_, requests(irq + 1), ierr)
            call mpi_irecv(BndDirRealInfo(:, :, Ldownx), 4*MAXNBEAM, MPI_real,&
                &iProcNeigh_(Ldownx), 2, comm2D_, requests(irq + 2), ierr)
            call mpi_irecv(BndDirIntInfo(:, :, Ldownx), 2*MAXNBEAM, MPI_integer,&
                &iProcNeigh_(Ldownx), 3, comm2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        if (isMpiY_) then
            call mpi_irecv(NINVALID(Ldowny), 1, MPI_integer,&
                &iProcNeigh_(Ldowny), 4, comm2D_, requests(irq + 1), ierr)
            call mpi_irecv(BndDirRealInfo(:, :, Ldowny), 4*MAXNBEAM, MPI_real,&
                &iProcNeigh_(Ldowny), 5, comm2D_, requests(irq + 2), ierr)
            call mpi_irecv(BndDirIntInfo(:, :, Ldowny), 2*MAXNBEAM, MPI_integer,&
                &iProcNeigh_(Ldowny), 6, comm2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        call mpi_waitall(irq, requests, status, ierr)

        NINVALID(Lupx) = 0; NINVALID(Lupy) = 0
        if (isMpiX_) then
            ! Continue the backwards ray integrations from the X boundary
            do J = 1, NINVALID(Ldownx)
                if (Ldownx .eq. 1) BndDirRealInfo(1, J, Ldownx) = rtm_obj_%x(1)
                if (Ldownx .eq. 2) BndDirRealInfo(1, J, Ldownx) = rtm_obj_%x(rtm_obj_%Nx)
                X = BndDirRealInfo(1, J, Ldownx); Y = BndDirRealInfo(2, J, Ldownx)
                Z = BndDirRealInfo(3, J, Ldownx); DIRPATH = BndDirRealInfo(4, J, Ldownx)
                if (Y .lt. (rtm_obj_%y(1) - EPS)) Y = YMAX; if (Y .gt. (YMAX + EPS)) Y = rtm_obj_%y(1)
                call shd_rad__computeDirect(0, X, Y, Z, rtm_obj_%flag_bc, rtm_obj_%flag_ip,&
                    &DELTAM, ML,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, DIRJUNK,&
                    &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                ! If we got a valid direct beam then store it in DirAllOut otherwise
                ! save the information to continue passing to neighboring processors
                if (VALIDBEAM) then
                    iproc = BndDirIntInfo(1, J, Ldownx)
                    Ndone(iproc) = Ndone(iproc) + 1
                    if (Ndone(iproc) .gt. MAXNBEAM)&
                        &call shd_mpi__abort('Error [mpi__makeDirect]: Ndone(iproc) exceeds MAXNBEAM x.')
                    DirAllOut(Ndone(iproc), iproc) = DIRPATH
                    DirPtrOut(Ndone(iproc), iproc) = BndDirIntInfo(2, J, Ldownx)
                else
                    if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) call shd_mpi__abort('Error [mpi__makeDirect]: Bad SIDE 2x.')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) .gt. MAXNBEAM)&
                        &call shd_mpi__abort('Error [mpi__makeDirect]: NINVALID(SIDE) exceeds MAXNBEAM 2x.')
                    BndDirRealInfo(:, NINVALID(SIDE), SIDE) = (/XO, YO, ZO, DIRPATH/)
                    BndDirIntInfo(:, NINVALID(SIDE), SIDE) = BndDirIntInfo(:, J, Ldownx)
                end if
            end do
        end if

        if (isMpiY_) then
            ! Continue the backwards ray integrations for the Y boundary
            do J = 1, NINVALID(Ldowny)
                if (Ldowny .eq. 3) BndDirRealInfo(2, J, Ldowny) = rtm_obj_%y(1)
                if (Ldowny .eq. 4) BndDirRealInfo(2, J, Ldowny) = rtm_obj_%y(rtm_obj_%Ny)
                X = BndDirRealInfo(1, J, Ldowny); Y = BndDirRealInfo(2, J, Ldowny)
                Z = BndDirRealInfo(3, J, Ldowny); DIRPATH = BndDirRealInfo(4, J, Ldowny)
                if (X .lt. (rtm_obj_%x(1) - EPS)) X = XMAX; if (X .gt. (XMAX + EPS)) X = rtm_obj_%x(1)
                call shd_rad__computeDirect(0, X, Y, Z, rtm_obj_%flag_bc, rtm_obj_%flag_ip,&
                    &DELTAM, ML,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, DIRJUNK,&
                    &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                ! If we got a valid direct beam then store it in DirAllOut otherwise
                ! save the information to continue passing to neighboring processors
                if (VALIDBEAM) then
                    iproc = BndDirIntInfo(1, J, Ldowny)
                    Ndone(iproc) = Ndone(iproc) + 1
                    if (Ndone(iproc) .gt. MAXNBEAM)&
                        &call shd_mpi__abort('Error [mpi__makeDirect]: Ndone(iproc) exceeds MAXNBEAM y.')
                    DirAllOut(Ndone(iproc), iproc) = DIRPATH
                    DirPtrOut(Ndone(iproc), iproc) = BndDirIntInfo(2, J, Ldowny)
                else
                    if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) call shd_mpi__abort('Error [mpi__makeDirect]: Bad SIDE 2y.')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) .gt. MAXNBEAM)&
                        &call shd_mpi__abort('Error [mpi__makeDirect]: NINVALID(SIDE) exceeds MAXNBEAM 2x.')
                    BndDirRealInfo(:, NINVALID(SIDE), SIDE) = (/XO, YO, ZO, DIRPATH/)
                    BndDirIntInfo(:, NINVALID(SIDE), SIDE) = BndDirIntInfo(:, J, Ldowny)
                end if
            end do
        end if

        ! See if there are any more invalid radiance rays on all the processors
        intsendbuf = (/NINVALID(Lupx), NINVALID(Lupy)/)
        call mpi_allreduce(intsendbuf, intrecvbuf, 2, MPI_integer, MPI_max,&
            &comm2D_, ierr)
        AllDone = all(intrecvbuf(1:2) .eq. 0)
        if ((maxval(intrecvbuf(1:2)) .gt. MAXNBEAM) .or. (maxval(Ndone(:)) .gt. MAXNBEAM)) then
            call shd_mpi__abort('Error [mpi__makeDirect]: MAXNBEAM exceeded 2.')
        end if
    end do

    ! Exchange the finished direct beam paths with all of the other processors
    do iproc = 0, Nproc_ - 1
        if (iproc .ne. myProc_) then
            call mpi_isend(Ndone(iproc), 1, MPI_integer,&
                &iproc, 7, comm2D_, requests(1), ierr)
            call mpi_isend(DirAllOut(:, iproc), Ndone(iproc), MPI_real,&
                &iproc, 8, comm2D_, requests(2), ierr)
            call mpi_isend(DirPtrOut(:, iproc), Ndone(iproc), MPI_integer,&
                &iproc, 9, comm2D_, requests(3), ierr)
            call mpi_irecv(Nbeam, 1, MPI_integer,&
                &iproc, 7, comm2D_, requests(4), ierr)
            call mpi_irecv(DirAllIn, MAXNBEAM, MPI_real,&
                &iproc, 8, comm2D_, requests(5), ierr)
            call mpi_irecv(DirPtrIn, MAXNBEAM, MPI_integer,&
                &iproc, 9, comm2D_, requests(6), ierr)
            call mpi_waitall(6, requests, status, ierr)
            ! Put the direct beam flux in the correct grid point
            do I = 1, Nbeam
                IP = DirPtrIn(I)
                DIRFLUX(IP) = SOLARFLUX*exp(-DirAllIn(I))
            end do
        else
            ! If we have the direct paths for this processor then avoid MPI and
            !   put the direct beam flux in the grid point
            do I = 1, Ndone(myProc_)
                IP = DirPtrOut(I, myProc_)
                DIRFLUX(IP) = SOLARFLUX*exp(-DirAllOut(I, myProc_))
            end do
        end if
    end do

    deallocate(DirAllOut, DirPtrOut, DirAllIn, DirPtrIn, Ndone)
    deallocate(BndDirRealInfo, BndDirIntInfo)
#endif

    return

end subroutine shd_mpi__makeDirect

subroutine shd_mpi__findBndryPts(BCFLAG, IPFLAG, SWEEPORD,&
    &Nx, Ny, Nz)

    ! Makes the pointers to and positions of the horizontal boundary grid
    ! points of the processors neighboring this one.
    ! For each of the potentially eight octants there are
    ! two sets of boundary pointers; 1 is for the current subdomain,
    ! and 2 is for the neighboring subdomain.
    ! Fills in the BxyPtr_(:,IBC,JOCT) and BxyPos_(:,IBC,JOCT) arrays
    ! (module variables), where IBC is a boundary point index and JOCT is
    ! a sort of octant number. BxyPtr_ holds (iproc,IPT,SIDE) for the boundary
    ! grid points from grid point IPT in processor iproc that are on the SIDE
    ! boundary of this subdomain (SIDE: 1=-X, 2=+X, 3=-Y, 4=+Y).  BxyPos_ holds
    ! the gridpoint positions (X,Y,Z) for these boundary points.

    integer, intent(in) :: BCFLAG, IPFLAG, SWEEPORD(rtm_obj_%Npts, *)
    integer, intent(in) :: Nx, Ny, Nz

#if HAS_MPI
    integer :: tag, ierr, irq
    integer :: MAXNBXY, NOCT, JOCT, IORDER, IPCELL, ICORNER, IPT, SIDE
    integer :: L, LI, N(4), N1(8, 4), N2(8, 4), NS, NBX, NBY
    integer :: MAXBC, IBC, BITZ, SZ, EZ, DZ, LX, LY, IX, IY, IZ, J
    integer :: INOUT_ARR(4) = (/2, 1, 4, 3/)
    integer, save :: OLDNPTS = 0
    logical :: SAMEGRID
    real    :: XEDGE(8), YEDGE(8), XYB(4)
    integer, allocatable :: BNDPTR(:, :, :, :), INDX(:)
    real, allocatable :: BNDPOS(:, :, :, :, :), ZARRAY(:)

    type(MPI_request) :: requests(24)
    type(MPI_status) :: status(24)

    if (.not. (isMpi_)) return
    if (is1D_) return
    call mpi_allreduce(rtm_obj_%Npts .eq. OLDNPTS, SAMEGRID, 1, MPI_logical, MPI_land, comm2D_, ierr)
    OLDNPTS = rtm_obj_%Npts
    if (SAMEGRID) return

    ! IOCT=1+BITX+2*BITY+4*BITZ ; JOCT=JOCTORDER(IOCT) = (/ 1,3,5,7,2,4,6,8 /)
    ! CX<0: BITX=1, JOCT=3,7,4,8      CY<0: BITY=1, JOCT=5,7,6,8
    XEDGE(:) = (/rtm_obj_%x(Nx), rtm_obj_%x(Nx), rtm_obj_%x(1), rtm_obj_%x(1), rtm_obj_%x(Nx), rtm_obj_%x(Nx), rtm_obj_%x(1), rtm_obj_%x(1)/)
    YEDGE(:) = (/rtm_obj_%y(Ny), rtm_obj_%y(Ny), rtm_obj_%y(Ny), rtm_obj_%y(Ny), rtm_obj_%y(1), rtm_obj_%y(1), rtm_obj_%y(1), rtm_obj_%y(1)/)
    XYB(1) = rtm_obj_%x(1); XYB(2) = rtm_obj_%x(Nx)
    XYB(3) = rtm_obj_%y(1); XYB(4) = rtm_obj_%y(Ny)

    if (is1D_) then
        NOCT = 2
    else if (is1DYx_) then
        NOCT = 4
    else
        NOCT = 8
    end if
    ! Allocate grid point pointer and position arrays (second from last
    !   index is 1 for local, 2 for neighbor; last index is boundary L)
    if (allocated(iBxyPtr_)) then
        MAXNBXY = nint(3.0*maxval(iBxyPtr_(2, :, :)))
    else
        MAXNBXY = 3.0*(Nz - 2)*max(Nx, Ny)
    end if
    allocate(BNDPTR(NOCT, MAXNBXY, 2, 4), BNDPOS(3, NOCT, MAXNBXY, 2, 4))

    N1(:, :) = 0
    do JOCT = 1, NOCT
        BITZ = 1 - ibits(JOCT, 0, 1)
        N(:) = 0
        do IORDER = 1, rtm_obj_%Npts
            IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
            ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
            IPT = i8Grid_Cell_(ICORNER, IPCELL)

            if  ((xyz_Rtm_(3, IPT) .lt. rtm_obj_%z(Nz) .and. BITZ .eq. 0) .or.&
                &(xyz_Rtm_(3, IPT) .gt. rtm_obj_%z(1) .and. BITZ .eq. 1)) then
                if (isMpiX_ .and. xyz_Rtm_(1, IPT) .eq. XEDGE(JOCT)&
                    &.and. .not. is1DX_) then
                    SIDE = 2 - ibits(JOCT - 1, 1, 1)
                    N(SIDE) = N(SIDE) + 1
                    if (N(SIDE) .gt. MAXNBXY) then
                        call shd_mpi__abort('Error [mpi__findBndryPts]: MAXNBXY exceeded by N(SIDE) 1x.')
                    end if
                    BNDPTR(JOCT, N(SIDE), 1, SIDE) = IPT
                    BNDPOS(:, JOCT, N(SIDE), 1, SIDE) = xyz_Rtm_(:, IPT)
                end if
                if (isMpiY_ .and. xyz_Rtm_(2, IPT) .eq. YEDGE(JOCT)&
                    &.and. .not. is1DYx_) then
                    SIDE = 4 - ibits(JOCT - 1, 2, 1)
                    N(SIDE) = N(SIDE) + 1
                    if (N(SIDE) .gt. MAXNBXY) then
                        call shd_mpi__abort('Error [mpi__findBndryPts]: MAXNBXY exceeded by N(SIDE) 1y.')
                    end if
                    BNDPTR(JOCT, N(SIDE), 1, SIDE) = IPT
                    BNDPOS(:, JOCT, N(SIDE), 1, SIDE) = xyz_Rtm_(:, IPT)
                end if
            end if
        end do
        N1(JOCT, :) = N(:)
    end do

    NS = maxval(N1(:, :))
    if (NS .gt. MAXNBXY) call shd_mpi__abort('Error [mpi__findBndryPts]: MAXNBXY might be exceeded.')
    N2(:, :) = 0
    irq = 0
    do L = 1, 4
        if ((isMpiX_ .and. (L .eq. 1 .or. L .eq. 2)) .or.&
            &(isMpiY_ .and. (L .eq. 3 .or. L .eq. 4))) then
            ! Exchange the lists of boundary points with neighbors
            LI = INOUT_ARR(L)
            tag = L
            call mpi_isend(N1(:, L), 8, MPI_integer,&
                &iProcNeigh_(L), tag, comm2D_, requests(irq + 1), ierr)
            call mpi_irecv(N2(:, LI), 8, MPI_integer,&
                &iProcNeigh_(LI), tag, comm2D_, requests(irq + 2), ierr)
            tag = 1000 + L
            call mpi_isend(BNDPTR(:, :, 1, L), NOCT*NS, MPI_integer,&
                &iProcNeigh_(L), tag, comm2D_, requests(irq + 3), ierr)
            call mpi_irecv(BNDPTR(:, :, 2, LI), NOCT*MAXNBXY, MPI_integer,&
                &iProcNeigh_(LI), tag, comm2D_, requests(irq + 4), ierr)
            tag = 2000 + L
            call mpi_isend(BNDPOS(:, :, :, 1, L), 3*NOCT*NS, MPI_real,&
                &iProcNeigh_(L), tag, comm2D_, requests(irq + 5), ierr)
            call mpi_irecv(BNDPOS(:, :, :, 2, LI), 3*NOCT*MAXNBXY, MPI_real,&
                &iProcNeigh_(LI), tag, comm2D_, requests(irq + 6), ierr)
            irq = irq + 6
        end if
    end do
    ! Wait for the transfers among all the neighbors to resolve
    call mpi_waitall(irq, requests, status, ierr)
    NS = maxval(N2(:, :))
    if (NS .gt. MAXNBXY) call shd_mpi__abort('Error [mpi__findBndryPts]: MAXNBXY was exceeded.')

    if (allocated(BxyPtr_)) deallocate(BxyPtr_, BxyPos_, iBxyPtr_)
    MAXBC = maxval((/N2(:, 1) + N2(:, 3), N2(:, 1) + N2(:, 4), N2(:, 2) + N2(:, 3), N2(:, 2) + N2(:, 4)/))
    allocate(BxyPtr_(3, MAXBC, NOCT), BxyPos_(3, MAXBC, NOCT), iBxyPtr_(2, Nz, NOCT))

    ! Put the boundary points from the neighboring subdomains in the
    ! BxyPtr_ and BxyPos_ arrays.  Merge in X and Y boundary points one Z slab
    ! at a time, and set pointers (iBxyPtr_) to start and end of each slab.
    !  DO IBC = iBxyPtr_(1,IZ,JOCT), iBxyPtr_(2,IZ,JOCT)
    !    X=BxyPos_(1,IBC,JOCT) ; Y=BxyPos_(2,IBC,JOCT) ; Z=BxyPos_(3,IBC,JOCT)
    !    iproc=BxyPtr_(1,IBC,JOCT) ; IPT=BxyPtr_(2,IBC,JOCT) ; SIDE=BxyPtr_(2,IBC,JOCT)

    do JOCT = 1, NOCT
        BITZ = 1 - ibits(JOCT, 0, 1)
        DZ = 2*BITZ - 1
        LX = 1 + ibits(JOCT - 1, 1, 1)
        LY = 3 + ibits(JOCT - 1, 2, 1)
        ! Combine the indices to the X and Y boundary points and sort on Z
        NBX = 0; if (isMpiX_) NBX = N2(JOCT, LX)
        NBY = 0; if (isMpiY_) NBY = N2(JOCT, LY)
        allocate(ZARRAY(NBX + NBY), INDX(NBX + NBY))
        ZARRAY(1:NBX) = BNDPOS(3, JOCT, 1:NBX, 2, LX)
        INDX(1:NBX) = (/(J, J=1, NBX)/)
        ZARRAY(NBX + 1:NBX + NBY) = BNDPOS(3, JOCT, 1:NBY, 2, LY)
        INDX(NBX + 1:NBX + NBY) = (/(-J, J=1, NBY)/)
        call shd_util__ssort(ZARRAY, INDX, NBX + NBY, 2*DZ)

        ! For each Z slab use Z sorted index to transfer the X or Y boundary points
        iBxyPtr_(:, :, JOCT) = 0
        SZ = (Nz - 1)*(1 - BITZ) + 1 + DZ
        EZ = (Nz - 1)*(BITZ) + 1
        IBC = 0; J = 1
        do IZ = SZ, EZ, DZ
            iBxyPtr_(1, IZ, JOCT) = IBC + 1
            ! Store the boundary point info for this Z slab
            do while ((J .le. (NBX + NBY)) .and.&
                &(((ZARRAY(J) .ge. rtm_obj_%z(IZ)) .and. (BITZ .eq. 0)) .or.&
                 &((ZARRAY(J) .le. rtm_obj_%z(IZ)) .and. (BITZ .eq. 1))))
                IBC = IBC + 1
                if (INDX(J) .gt. 0) then
                    IX = INDX(J)
                    BxyPtr_(:, IBC, JOCT) = (/iProcNeigh_(LX), BNDPTR(JOCT, IX, 2, LX), LX/)
                    BxyPos_(:, IBC, JOCT) = BNDPOS(:, JOCT, IX, 2, LX)
                    BxyPos_(int((LX + 1)/2), IBC, JOCT) = XYB(LX)
                else
                    IY = abs(INDX(J))
                    BxyPtr_(:, IBC, JOCT) = (/iProcNeigh_(LY), BNDPTR(JOCT, IY, 2, LY), LY/)
                    BxyPos_(:, IBC, JOCT) = BNDPOS(:, JOCT, IY, 2, LY)
                    BxyPos_(int((LY + 1)/2), IBC, JOCT) = XYB(LY)
                end if
                J = J + 1
            end do
            iBxyPtr_(2, IZ, JOCT) = IBC
        end do
        deallocate(ZARRAY, INDX)
    end do

    deallocate(BNDPTR, BNDPOS)
#endif

    return

end subroutine shd_mpi__findBndryPts

subroutine shd_mpi__calBndryRad(BCFLAG, IPFLAG, JOCT, IZ,&
    &Nx, Ny, Nz,&
    &NA,&
    &MU, PHI, EXTINCT, Nstokes, SOURCE,&
    &KANG, GRIDRAD)

    ! Calculates the discrete ordinate boundary radiances for the KANG discrete
    ! ordinate (going in direction MU/PHI) for the grid points in the IZ slab.
    ! shd_mpi__findBndryPts has already been called to set up the BxyPtr_ and
    ! BxyPos_ arrays for the grid points on the boundaries, and the iBxyPtr_
    ! array which has the starting and ending points for the IZ slab.  The
    ! grid points are done in the sweeping order, and so the octant (JOCT)
    ! is used to index into the iBxyPtr_, BxyPtr_, and BxyPos_ arrays.  The
    ! grid points specified in BxyPtr_/BxyPos_ are on the boundaries of this
    ! processor, and need to be traced back along the discrete ordinate
    ! first within this subdomain.  The radiative transfer integrations
    ! with the extinction (EXTINCT) and source function (SOURCE) are done
    ! with the shd_mpi__integrateRay routine.  If the top or bottom boundary is
    ! hit when tracing back then the radiance is saved, otherwise the
    ! required information (X,Y,Z, transmission, partial radiance, original
    ! processor, and starting grid point) are saved to be sent to the
    ! appropriate neighboring processor.  These data are received from the
    ! neighboring processors and the radiative transfer integrations are
    ! continued, looping until all the rays have been finished.  Finally,
    ! the finished radiance and original grid points are exchanged with
    ! all the other processors.

    integer, intent(in) :: BCFLAG, IPFLAG, JOCT, IZ
    integer, intent(in) :: Nstokes, Nx, Ny, Nz, NA, KANG
    real, intent(in) :: MU, PHI
    real, intent(in) :: EXTINCT(rtm_obj_%Npts), SOURCE(Nstokes, NA, rtm_obj_%Npts)
    real, intent(inout) :: GRIDRAD(Nstokes, rtm_obj_%Npts)

#if HAS_MPI
    integer :: MAXNBND, NptsBC, Lupx, Lupy, Ldownx, Ldowny
    integer :: SIDE, NINVALID(0:4)
    integer :: IBC, IPT, I, J, L
    integer :: ierr(12), irq, iproc
    integer :: intsendbuf(2), intrecvbuf(2)
    logical :: VALIDRAD, AllDone
    real    :: X, Y, Z, RADIANCE(Nstokes), TRANSMIT, xo, yo, zo
    integer, allocatable :: Ndone(:)
    integer, allocatable :: BndRadIntInfo(:, :, :)
    real, allocatable    :: BndRadRealInfo(:, :, :)
    integer, allocatable :: RadPtrOut(:, :), RadPtrIn(:)
    real, allocatable    :: RadAllOut(:, :, :), RadAllIn(:, :)

    type(MPI_request) :: requests(12)
    type(MPI_status) :: status(12)

    MAXNBND = Nproc_*maxval(iBxyPtr_(2, :, :) - iBxyPtr_(1, :, :) + 1)
    allocate(BndRadRealInfo(4 + Nstokes, MAXNBND, 4), BndRadIntInfo(2, MAXNBND, 4))
    allocate(RadAllOut(Nstokes, MAXNBND, 0:Nproc_ - 1), RadPtrOut(MAXNBND, 0:Nproc_ - 1))
    allocate(RadAllIn(Nstokes, MAXNBND), RadPtrIn(MAXNBND), Ndone(0:Nproc_ - 1))

    Lupx = 0; Lupy = 0; Ldownx = 0; Ldowny = 0
    if (isMpiX_) then
        if (cos(PHI) .gt. 1.0e-5) then
            Ldownx = 2; Lupx = 1
        else
            Ldownx = 1; Lupx = 2
        end if
    end if
    if (isMpiY_) then
        if (sin(PHI) .gt. 1.0e-5) then
            Ldowny = 4; Lupy = 3
        else
            Ldowny = 3; Lupy = 4
        end if
    end if

    ! Do the initial ray tracing
    Ndone(:) = 0
    NINVALID(:) = 0
    ! Loop over the downstream boundary points from the neighbor subdomains
    do IBC = iBxyPtr_(1, IZ, JOCT), iBxyPtr_(2, IZ, JOCT)
        ! Call shd_mpi__integrateRay for this
        X = BxyPos_(1, IBC, JOCT); Y = BxyPos_(2, IBC, JOCT); Z = BxyPos_(3, IBC, JOCT)
        iproc = BxyPtr_(1, IBC, JOCT); IPT = BxyPtr_(2, IBC, JOCT); SIDE = BxyPtr_(3, IBC, JOCT)
        RADIANCE(:) = 0.0; TRANSMIT = 1.0
        xo = x; yo = y; zo = z
        call shd_mpi__integrateRay(BCFLAG, IPFLAG, Nx, Ny, Nz, NA,&
            &EXTINCT, Nstokes, SOURCE, KANG, GRIDRAD,&
            &MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

        ! If we got a radiance then store it otherwise save
        !  the information needed to pass to other processors
        if (VALIDRAD) then
            Ndone(iproc) = Ndone(iproc) + 1
            RadAllOut(:, Ndone(iproc), iproc) = RADIANCE(:)
            RadPtrOut(Ndone(iproc), iproc) = IPT
        else
            if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) then
                write (6, *) 'Bad SIDE 1:', myProc_, mu, phi, SIDE, x, y, z, xo, yo, zo
                call shd_mpi__abort('Error [mpi__calBndryRad]')
            end if
            NINVALID(SIDE) = NINVALID(SIDE) + 1
            BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/X, Y, Z, TRANSMIT, RADIANCE(:)/)
            BndRadIntInfo(:, NINVALID(SIDE), SIDE) = (/iproc, IPT/)
        end if
    end do

    ! Cycle over the passing of info between processors and doing more
    !  partial radiance integrations until all rays have been finished
    AllDone = .false.
    do while (.not. AllDone)
        ! Do the send and receives of the number of incomplete rays and
        !   the real and integer info needed about the radiance integrations
        irq = 0
        if (isMpiX_) then
            call mpi_isend(NINVALID(Lupx), 1, MPI_integer,&
                &iProcNeigh_(Lupx), 1, comm2D_, requests(irq + 1), ierr(irq + 1))
            call mpi_isend(BndRadRealInfo(:, :, Lupx), (4 + Nstokes)*NINVALID(Lupx), MPI_real,&
                &iProcNeigh_(Lupx), 2, comm2D_, requests(irq + 2), ierr(irq + 2))
            call mpi_isend(BndRadIntInfo(:, :, Lupx), 2*NINVALID(Lupx), MPI_integer,&
                &iProcNeigh_(Lupx), 3, comm2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (isMpiY_) then
            call mpi_isend(NINVALID(Lupy), 1, MPI_integer,&
                &iProcNeigh_(Lupy), 4, comm2D_, requests(irq + 1), ierr(irq + 1))
            call mpi_isend(BndRadRealInfo(:, :, Lupy), (4 + Nstokes)*NINVALID(Lupy), MPI_real,&
                &iProcNeigh_(Lupy), 5, comm2D_, requests(irq + 2), ierr(irq + 2))
            call mpi_isend(BndRadIntInfo(:, :, Lupy), 2*NINVALID(Lupy), MPI_integer,&
                &iProcNeigh_(Lupy), 6, comm2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (isMpiX_) then
            call mpi_irecv(NINVALID(Ldownx), 1, MPI_integer,&
                &iProcNeigh_(Ldownx), 1, comm2D_, requests(irq + 1), ierr(irq + 1))
            call mpi_irecv(BndRadRealInfo(:, :, Ldownx), (4 + Nstokes)*MAXNBND, MPI_real,&
                &iProcNeigh_(Ldownx), 2, comm2D_, requests(irq + 2), ierr(irq + 2))
            call mpi_irecv(BndRadIntInfo(:, :, Ldownx), 2*MAXNBND, MPI_integer,&
                &iProcNeigh_(Ldownx), 3, comm2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (isMpiY_) then
            call mpi_irecv(NINVALID(Ldowny), 1, MPI_integer,&
                &iProcNeigh_(Ldowny), 4, comm2D_, requests(irq + 1), ierr(irq + 1))
            call mpi_irecv(BndRadRealInfo(:, :, Ldowny), (4 + Nstokes)*MAXNBND, MPI_real,&
                &iProcNeigh_(Ldowny), 5, comm2D_, requests(irq + 2), ierr(irq + 2))
            call mpi_irecv(BndRadIntInfo(:, :, Ldowny), 2*MAXNBND, MPI_integer,&
                &iProcNeigh_(Ldowny), 6, comm2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (any(ierr(1:irq) .ne. MPI_success)) then
            write (6, *) 'MPI error: ', ierr(1:irq)
            call shd_mpi__abort('Error [mpi__calBndryRad]: MPI error 1.')
        end if
        call mpi_waitall(irq, requests, status, ierr(1))
        if (ierr(1) .ne. MPI_success) then
            write (6, *) 'mpi_waitall error: ', ierr(1)
            call shd_mpi__abort('Error [mpi__calBndryRad]: MPI error 2.')
        end if

        NINVALID(Lupx) = 0; NINVALID(Lupy) = 0
        if (isMpiX_) then
            ! Continue the backwards ray integrations for the X boundary
            do J = 1, NINVALID(Ldownx)
                if (Ldownx .eq. 1) BndRadRealInfo(1, J, Ldownx) = rtm_obj_%x(1)
                if (Ldownx .eq. 2) BndRadRealInfo(1, J, Ldownx) = rtm_obj_%x(Nx)
                X = BndRadRealInfo(1, J, Ldownx); Y = BndRadRealInfo(2, J, Ldownx)
                Z = BndRadRealInfo(3, J, Ldownx)
                SIDE = Ldownx; TRANSMIT = BndRadRealInfo(4, J, Ldownx)
                RADIANCE(:) = BndRadRealInfo(5:4 + Nstokes, J, Ldownx)
                call shd_mpi__integrateRay(BCFLAG, IPFLAG, Nx, Ny, Nz, NA,&
                    &EXTINCT, Nstokes, SOURCE, KANG, GRIDRAD,&
                    &MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

                ! If we got a radiance then store it in RadAllOut (if for another
                !   processor) otherwise save the information needed to continue
                !   passing to neighboring processors
                if (VALIDRAD) then
                    iproc = BndRadIntInfo(1, J, Ldownx)
                    Ndone(iproc) = Ndone(iproc) + 1
                    RadAllOut(:, Ndone(iproc), iproc) = RADIANCE(:)
                    RadPtrOut(Ndone(iproc), iproc) = BndRadIntInfo(2, J, Ldownx)
                else
                    if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) then
                        write (6, *) 'Bad SIDE 2x:', myProc_, mu, phi, SIDE, x, y, z
                        call shd_mpi__abort('Error [mpi__calBndryRad]: Bad SIDE 2x, see above.')
                    end if
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/X, Y, Z, TRANSMIT, RADIANCE(:)/)
                    BndRadIntInfo(:, NINVALID(SIDE), SIDE) = BndRadIntInfo(:, J, Ldownx)
                end if
            end do
        end if

        if (isMpiY_) then
            ! Continue the backwards ray integrations for the Y boundary
            do J = 1, NINVALID(Ldowny)
                if (Ldowny .eq. 3) BndRadRealInfo(2, J, Ldowny) = rtm_obj_%y(1)
                if (Ldowny .eq. 4) BndRadRealInfo(2, J, Ldowny) = rtm_obj_%y(Ny)
                X = BndRadRealInfo(1, J, Ldowny); Y = BndRadRealInfo(2, J, Ldowny)
                Z = BndRadRealInfo(3, J, Ldowny)
                SIDE = Ldowny; TRANSMIT = BndRadRealInfo(4, J, Ldowny)
                RADIANCE(:) = BndRadRealInfo(5:4 + Nstokes, J, Ldowny)
                call shd_mpi__integrateRay(BCFLAG, IPFLAG, Nx, Ny, Nz, NA,&
                    &EXTINCT, Nstokes, SOURCE, KANG, GRIDRAD,&
                    &MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

                ! If we got a radiance then store it in RadAllOut (if for another
                !   processor) otherwise save the information needed to continue
                !   passing to neighboring processors
                if (VALIDRAD) then
                    iproc = BndRadIntInfo(1, J, Ldowny)
                    Ndone(iproc) = Ndone(iproc) + 1
                    RadAllOut(:, Ndone(iproc), iproc) = RADIANCE(:)
                    RadPtrOut(Ndone(iproc), iproc) = BndRadIntInfo(2, J, Ldowny)
                else
                    if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) then
                        write (6, *) 'Bad SIDE 2y:', myProc_, mu, phi, SIDE, x, y, z
                        call shd_mpi__abort('Error [mpi__calBndryRad]: Bad SIDE 2y, see above.')
                    end if
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/X, Y, Z, TRANSMIT, RADIANCE(:)/)
                    BndRadIntInfo(:, NINVALID(SIDE), SIDE) = BndRadIntInfo(:, J, Ldowny)
                end if
            end do
        end if

        ! See if there are any more invalid radiance rays on all the processors
        intsendbuf = (/NINVALID(Lupx), NINVALID(Lupy)/)
        call mpi_allreduce(intsendbuf, intrecvbuf, 2, MPI_integer, MPI_max,&
            &comm2D_, ierr(1))
        AllDone = all(intrecvbuf(1:2) .eq. 0)
    end do

    ! Exchange the finished boundary radiances with all of the other processors
    do iproc = 0, Nproc_ - 1
        if (iproc .ne. myProc_) then
            call mpi_isend(Ndone(iproc), 1, MPI_integer,&
                &iproc, 7, comm2D_, requests(1), ierr(1))
            call mpi_isend(RadAllOut(:, :, iproc), Nstokes*Ndone(iproc), MPI_real,&
                &iproc, 8, comm2D_, requests(2), ierr(2))
            call mpi_isend(RadPtrOut(:, iproc), Ndone(iproc), MPI_integer,&
                &iproc, 9, comm2D_, requests(3), ierr(3))
            call mpi_irecv(NptsBC, 1, MPI_integer,&
                &iproc, 7, comm2D_, requests(4), ierr(4))
            call mpi_irecv(RadAllIn, Nstokes*MAXNBND, MPI_real,&
                &iproc, 8, comm2D_, requests(5), ierr(5))
            call mpi_irecv(RadPtrIn, MAXNBND, MPI_integer,&
                &iproc, 9, comm2D_, requests(6), ierr(6))
            call mpi_waitall(6, requests, status, ierr(7))
            if (any(ierr(1:7) .ne. MPI_success)) then
                write (6, *) 'MPI error: ', ierr(1:7)
                call shd_mpi__abort('Error [mpi__calBndryRad]: MPI error 3.')
            end if
            do I = 1, NptsBC
                IPT = RadPtrIn(I)
                GRIDRAD(:, IPT) = RadAllIn(:, I)
            end do
        else
            do I = 1, Ndone(myProc_)
                IPT = RadPtrOut(I, myProc_)
                GRIDRAD(:, IPT) = RadAllOut(:, I, myProc_)
            end do
        end if
    end do

    deallocate(RadAllOut, RadPtrOut, RadAllIn, RadPtrIn, Ndone)
    deallocate(BndRadRealInfo, BndRadIntInfo)
#endif

    return

end subroutine shd_mpi__calBndryRad

subroutine shd_mpi__computeRad(&
    &DELTAM,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &DIRFLUX, FLUXES,&
    &src1_Shdo, GRIDRAD,&
    &parms_Out, NRAD, RADOUT)

    ! Computes the radiances for the locations and directions specified
    ! in parms_Out, putting the results in RADOUT and increasing NRAD by
    ! the number of radiances output.  Does the radiative transfer integrations
    ! across multiple processors.

    integer, intent(inout) :: NRAD
    integer, intent(in) :: MAXNBC, NTOPPTS, NBOTPTS
    integer, intent(in) :: BCPTR(MAXNBC, 2)
    logical, intent(in) :: DELTAM
    real, intent(in) :: BCRAD(shdo_obj_%Nstokes, *)
    real, intent(in) :: DIRFLUX(rtm_obj_%Npts), FLUXES(2, rtm_obj_%Npts)
    real, intent(inout) :: src1_Shdo(shdo_obj_%Nstokes, rtm_obj_%Npts), GRIDRAD(shdo_obj_%Nstokes, rtm_obj_%Npts)
    real, intent(in) :: parms_Out(*)
    real, intent(out) :: RADOUT(shdo_obj_%Nstokes, *)

#if HAS_MPI
    integer :: IBC, JX, JY, K
    integer :: NANGOUT, Nx_out, Ny_out
    logical :: VALIDRAD
    real(RD_) :: X0, Y0, Z0, XE, YE, ZE, X, Y, Z, TRANSMIT, RADIANCE(shdo_obj_%Nstokes)
    real    :: mu_out, phi_out, phiDeg_out
    real    :: STARTX, STARTY, EPS = 1.0e-5

    integer :: MAXNRAD, NDOWN, Lupx, Lupy, Ldownx, Ldowny
    integer :: SIDE, NINVALID(0:4)
    integer :: I, J, L, Ndone, intsendbuf(2), intrecvbuf(2)
    integer :: ierr, irq, iproc
    logical :: AllDone
    integer, allocatable :: NradRecv(:)
    integer, allocatable :: IndexSend(:), IndexRecv(:, :), BndRadIntInfo(:, :)
    real, allocatable    :: RadSend(:, :), RadRecv(:, :, :), BndRadRealInfo(:, :, :)

    type(MPI_request), allocatable :: requests(:)
    type(MPI_status), allocatable :: status(:)

    ! added for fixing the "domain shift" issue
    !╭──────────────────────────────────────────────────────────────╮!
    real(RD_) :: DS_Z0, DS_VZA_RAD, DS_XOFFSET, DS_YOFFSET, TOT_XOFFSET, TOT_YOFFSET
    !╰──────────────────────────────────────────────────────────────╯!

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(src_obj_%type, src_obj_%skyRad, src_obj_%waveNum, src_obj_%avgWaveLen,&
        &src_obj_%unit, NTOPPTS, shdo_obj_%Nstokes, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces.
    ! Compute the upwelling bottom radiances using the downwelling fluxes.
    if (sfc_obj_%type .eq. 'FL') then
        call shd_init__makeBotRadUL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, src_obj_%type,&
            &src_obj_%waveNum, src_obj_%avgWaveLen, src_obj_%unit, shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
    else if (sfc_obj_%type .eq. 'VL') then
        call shd_init__makeBotRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, src_obj_%type,&
            &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
    end if

    ! Setup the regularly spaced radiance locations
    STARTX = 0.0
    if (parms_Out(2) .le. 0.0) then
        Nx_out = 1
    else
        Nx_out = max(1, nint(prp_obj_%xDomain/parms_Out(2)))
    end if
    STARTY = 0.0
    if (parms_Out(3) .le. 0.0) then
        Ny_out = 1
    else
        Ny_out = max(1, nint(prp_obj_%yDomain/parms_Out(3)))
    end if
    Z0 = min(max(parms_Out(1), rtm_obj_%z(1)), rtm_obj_%z(rtm_obj_%Nz))
    NANGOUT = nint(parms_Out(6))

    MAXNRAD = nint(4.0*(Nx_out*Ny_out)/Nproc_)
    ! MAXNRAD = Nx_out*Ny_out
    allocate(RadSend(shdo_obj_%Nstokes, MAXNRAD), IndexSend(MAXNRAD))
    if (myProc_ .eq. 0) allocate(NradRecv(Nproc_ - 1),&
        &RadRecv(shdo_obj_%Nstokes, MAXNRAD, Nproc_ - 1), IndexRecv(MAXNRAD, Nproc_ - 1))
    allocate(BndRadRealInfo(4 + shdo_obj_%Nstokes, MAXNRAD, 4), BndRadIntInfo(MAXNRAD, 4))
    allocate(requests(max(12, 3*Nproc_)), status(max(12, 3*Nproc_)))

    ! Loop over the radiance directions
    ANGOUT_LOOP: do K = 1, NANGOUT
        if (sen_obj_%fname(1:4) .ne. 'NONE') then
            mu_out     = view_Sen_(1, K)
            phiDeg_out = view_Sen_(2, K)
        else
            mu_out = parms_Out(2*K + 5)
            phiDeg_out = parms_Out(2*K + 6)
        endif
        phi_out = phiDeg_out*acos(-1.0)/180.0

        ! added for fixing the "domain shift" issue (only for 3D)
        ! this is the first step correction:
        !   single-pixel-wise mini adjustment, addressing
        !   interpolation effect
        ! the second step correction is done in the data output
        !   <shd_io__radShiftCorr> in shdOut.f90
        !╭────────────────────────────────────────────────╮!
        if (rtm_obj_%flag_ip .eq. 0) then
            DS_Z0 = Z0 - rtm_obj_%z(1)
            DS_VZA_RAD = acos(mu_out)
            DS_XOFFSET = DS_Z0*dtan(DS_VZA_RAD)*cos(phi_out)
            DS_YOFFSET = DS_Z0*dtan(DS_VZA_RAD)*sin(phi_out)
            TOT_XOFFSET = DS_XOFFSET
            TOT_YOFFSET = DS_YOFFSET
        else
            TOT_XOFFSET = 0.0
            TOT_YOFFSET = 0.0
        endif

        ! to avoid vertical stripes
        ! this extra pixel shift is considered in the data output (second step correction)
        STARTX = TOT_XOFFSET - parms_Out(2)*nint(TOT_XOFFSET/parms_Out(2))
        if (STARTX .lt. 0.0d0) then
            STARTX = STARTX + parms_Out(2)
        end if

        ! to avoid horizontal stripes
        ! this extra pixel shift is considered in the data output (second step correction)
        STARTY = TOT_YOFFSET - parms_Out(3)*nint(TOT_YOFFSET/parms_Out(3))
        if (STARTY .lt. 0.0d0) then
            STARTY = STARTY + parms_Out(3)
        end if
        !╰────────────────────────────────────────────────╯!

        if (mu_out .eq. 0.0) then
            if (verbose_ .ge. 4) write (6, *) 'Warning [mpi__computeRad]: Bad mu_out for radiance', K, mu_out, phiDeg_out
            mu_out = accSolv_
            write (6, *) '                           New mu_out =', mu_out
        end if

        if (abs(mu_out) .gt. 1.0) then
            if (verbose_ .ge. 4) write (6, *) 'Warning [mpi__computeRad]: Bad mu_out for radiance', K, mu_out, phiDeg_out
            if (mu_out .gt.  1.0) mu_out =  1.0
            if (mu_out .lt. -1.0) mu_out = -1.0
            write (6, *) '                           New mu_out =', mu_out
        end if

        ! Compute the source function throughout grid for this angle
        call shd_src__computeSrcOneDir(&
            &DELTAM, mu_out, phi_out,&
            &DIRFLUX, src1_Shdo,&
            &isrcJ_Shdo_, srcJ_Shdo_)

        ! Set the radiance field to -1, so we can determine valid radiances
        !   Also set the source array to the extinction times the source
        !   function, so the grid interpolation is correct.
        do I = 1, rtm_obj_%Npts
            GRIDRAD(1, I) = -1.0
            src1_Shdo(:, I) = src1_Shdo(:, I)*ext_Rtm_(I)
        end do
        ! Get boundary radiances: either top or bottom
        ! Isotropic top boundary or Lambertian bottom boundary can use
        !   the previously computed boundary radiances in BCRAD,
        !   otherwise, compute the radiance for this angle
        !   by integrating over the stored downwelling radiances.
        if (mu_out .lt. 0.0) then
            do IBC = 1, NTOPPTS
                I = BCPTR(IBC, 1)
                GRIDRAD(:, I) = BCRAD(:, IBC)
            end do
        else
            if (.not. sfc_obj_%isLambertian) then
                call shd_init__makeBotRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                    &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt, mu_out, phi_out,&
                    &src_obj_%type, src_obj_%avgWaveLen, src_obj_%avgMu, src_obj_%avgAz, DIRFLUX,&
                    &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
            end if
            do IBC = 1, NBOTPTS
                I = BCPTR(IBC, 2)
                GRIDRAD(:, I) = BCRAD(:, NTOPPTS + IBC)
            end do
        end if

        ! Set up the upstream/downstream boundaries for this direction
        Lupx = 0; Lupy = 0; Ldownx = 0; Ldowny = 0
        if (isMpiX_) then
            if ((cos(phi_out) .gt. 1.0e-5)) then
                Ldownx = 2; Lupx = 1
            else
                Ldownx = 1; Lupx = 2
            end if
        end if
        if (isMpiY_) then
            !if ((sin(phi_out) .gt. 1.0e-5) .and. (prp_obj_%Ny .gt. 1)) then
            if ((sin(phi_out) .gt. 1.0e-5)) then
                Ldowny = 4; Lupy = 3
            else
                Ldowny = 3; Lupy = 4
            end if
        end if

        NINVALID(:) = 0
        Ndone = 0

        ! Loop over all the radiance starting locations and integrate
        !   backwards for those in this subdomain.
        Y0 = STARTY
        do JY = 1, Ny_out
            X0 = STARTX
            do JX = 1, Nx_out
                NRAD = NRAD + 1

                if (X0 .ge. (rtm_obj_%x(1) - EPS) .and. X0 .le. (rtm_obj_%x(rtm_obj_%Nx) + EPS) .and.&
                    Y0 .ge. (rtm_obj_%y(1) - EPS) .and. Y0 .le. (rtm_obj_%y(rtm_obj_%Ny) + EPS)) then
                    TRANSMIT = 1.0d0; RADIANCE(:) = 0.0d0
                    call shd_src__integrateSrc(shdo_obj_%Nstokes,&
                        &mu_out, phi_out, GRIDRAD, src1_Shdo,&
                        &X0, Y0, Z0, XE, YE, ZE, SIDE,&
                        &TRANSMIT, RADIANCE, VALIDRAD)

                    ! If we got a radiance then store it in RadSend otherwise save
                    !  the information needed to pass to other processors
                    if (VALIDRAD) then
                        Ndone = Ndone + 1
                        if (Ndone .gt. MAXNRAD)&
                            &call shd_mpi__abort('Error [mpi__computeRad]: MAXNRAD exceeded by Ndone 1.')
                        RadSend(:, Ndone) = RADIANCE(:)
                        IndexSend(Ndone) = NRAD
                    else
                        if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) then
                            write (6, *) 'Bad SIDE 1:', myProc_, mu_out, phi_out, x0, y0, z0, SIDE, xe, ye, ze
                            call shd_mpi__abort('Error [mpi__computeRad]: Bad SIDE 1, see above.')
                        end if
                        NINVALID(SIDE) = NINVALID(SIDE) + 1
                        if (NINVALID(SIDE) .gt. MAXNRAD)&
                            call shd_mpi__abort('Error [mpi__computeRad]: MAXNRAD exceeded by NINVALID(SIDE) 1.')

                        BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, TRANSMIT, RADIANCE(:)/)
                        BndRadIntInfo(NINVALID(SIDE), SIDE) = NRAD
                    end if
                end if
                X0 = X0 + parms_Out(2)
            end do
            Y0 = Y0 + parms_Out(3)
        end do

        ! Cycle over the passing of info between processors and doing more
        !  partial radiance integrations until all rays have been finished
        AllDone = .false.
        RAYMPI_LOOP: do while (.not. AllDone)
            ! Do the send and receives of the number of incomplete rays and
            !   the real and integer info needed about the radiance integrations
            irq = 0
            if (isMpiX_) then
                call mpi_isend(NINVALID(Lupx), 1, MPI_integer,&
                    &iProcNeigh_(Lupx), 6*K + 1, comm2D_, requests(irq + 1), ierr)
                call mpi_isend(BndRadRealInfo(:, :, Lupx), (4 + shdo_obj_%Nstokes)*NINVALID(Lupx), MPI_double_precision,&
                    &iProcNeigh_(Lupx), 6*K + 2, comm2D_, requests(irq + 2), ierr)
                call mpi_isend(BndRadIntInfo(:, Lupx), NINVALID(Lupx), MPI_integer,&
                    &iProcNeigh_(Lupx), 6*K + 3, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
            if (isMpiY_) then
                call mpi_isend(NINVALID(Lupy), 1, MPI_integer,&
                    &iProcNeigh_(Lupy), 6*K + 4, comm2D_, requests(irq + 1), ierr)
                call mpi_isend(BndRadRealInfo(:, :, Lupy), (4 + shdo_obj_%Nstokes)*NINVALID(Lupy), MPI_double_precision,&
                    &iProcNeigh_(Lupy), 6*K + 5, comm2D_, requests(irq + 2), ierr)
                call mpi_isend(BndRadIntInfo(:, Lupy), NINVALID(Lupy), MPI_integer,&
                    &iProcNeigh_(Lupy), 6*K + 6, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
            if (isMpiX_) then
                call mpi_irecv(NINVALID(Ldownx), 1, MPI_integer,&
                    &iProcNeigh_(Ldownx), 6*K + 1, comm2D_, requests(irq + 1), ierr)
                call mpi_irecv(BndRadRealInfo(:, :, Ldownx), (4 + shdo_obj_%Nstokes)*MAXNRAD, MPI_double_precision,&
                    &iProcNeigh_(Ldownx), 6*K + 2, comm2D_, requests(irq + 2), ierr)
                call mpi_irecv(BndRadIntInfo(:, Ldownx), MAXNRAD, MPI_integer,&
                    &iProcNeigh_(Ldownx), 6*K + 3, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
            if (isMpiY_) then
                call mpi_irecv(NINVALID(Ldowny), 1, MPI_integer,&
                    &iProcNeigh_(Ldowny), 6*K + 4, comm2D_, requests(irq + 1), ierr)
                call mpi_irecv(BndRadRealInfo(:, :, Ldowny), (4 + shdo_obj_%Nstokes)*MAXNRAD, MPI_double_precision,&
                    &iProcNeigh_(Ldowny), 6*K + 5, comm2D_, requests(irq + 2), ierr)
                call mpi_irecv(BndRadIntInfo(:, Ldowny), MAXNRAD, MPI_integer,&
                    &iProcNeigh_(Ldowny), 6*K + 6, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
            call mpi_waitall(irq, requests, status, ierr)

            NINVALID(Lupx) = 0; NINVALID(Lupy) = 0
            if (isMpiX_) then
                ! Continue the backwards ray integrations for the X boundary
                do J = 1, NINVALID(Ldownx)
                    if (Ldownx .eq. 1) BndRadRealInfo(1, J, Ldownx) = rtm_obj_%x(1)
                    if (Ldownx .eq. 2) BndRadRealInfo(1, J, Ldownx) = rtm_obj_%x(rtm_obj_%Nx)
                    X = BndRadRealInfo(1, J, Ldownx); Y = BndRadRealInfo(2, J, Ldownx)
                    Z = BndRadRealInfo(3, J, Ldownx); TRANSMIT = BndRadRealInfo(4, J, Ldownx)
                    RADIANCE(:) = BndRadRealInfo(5:4 + shdo_obj_%Nstokes, J, Ldownx)
                    call shd_src__integrateSrc(shdo_obj_%Nstokes,&
                        &mu_out, phi_out, GRIDRAD, src1_Shdo,&
                        &X, Y, Z, XE, YE, ZE, SIDE,&
                        &TRANSMIT, RADIANCE, VALIDRAD)

                    ! If we got a radiance then store it in RadSend, otherwise
                    ! save information needed to pass to neighboring processors
                    if (VALIDRAD) then
                        Ndone = Ndone + 1
                        if (Ndone .gt. MAXNRAD)&
                            &call shd_mpi__abort('Error [mpi__computeRad]: MAXNRAD exceeded by Ndone 2.')
                        RadSend(:, Ndone) = RADIANCE(:)
                        IndexSend(Ndone) = BndRadIntInfo(J, Ldownx)
                    else
                        if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) then
                            write (6, *) 'Bad SIDE 2x', myProc_, mu_out, phi_out, x, y, z, SIDE, xe, ye, ze
                            call shd_mpi__abort('Error [mpi__computeRad]: Bad SIDE 2x, see above.')
                        end if
                        NINVALID(SIDE) = NINVALID(SIDE) + 1
                        if (NINVALID(SIDE) .gt. MAXNRAD)&
                            &call shd_mpi__abort('Error [mpi__computeRad]: MAXNRAD exceeded by NINVALID(SIDE) 2.')
                        BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, TRANSMIT, RADIANCE(:)/)
                        BndRadIntInfo(NINVALID(SIDE), SIDE) = BndRadIntInfo(J, Ldownx)
                    end if
                end do
            end if

            if (isMpiY_) then
                ! Continue the backwards ray integrations for the Y boundary
                do J = 1, NINVALID(Ldowny)
                    if (Ldowny .eq. 3) BndRadRealInfo(2, J, Ldowny) = rtm_obj_%y(1)
                    if (Ldowny .eq. 4) BndRadRealInfo(2, J, Ldowny) = rtm_obj_%y(rtm_obj_%Ny)
                    X = BndRadRealInfo(1, J, Ldowny); Y = BndRadRealInfo(2, J, Ldowny)
                    Z = BndRadRealInfo(3, J, Ldowny); TRANSMIT = BndRadRealInfo(4, J, Ldowny)
                    RADIANCE(:) = BndRadRealInfo(5:4 + shdo_obj_%Nstokes, J, Ldowny)
                    call shd_src__integrateSrc(shdo_obj_%Nstokes,&
                        &mu_out, phi_out, GRIDRAD, src1_Shdo,&
                        &X, Y, Z, XE, YE, ZE, SIDE,&
                        &TRANSMIT, RADIANCE, VALIDRAD)

                    ! If we got a radiance then store it in RadSend, otherwise
                    ! save information needed to pass to neighboring processors
                    if (VALIDRAD) then
                        Ndone = Ndone + 1
                        if (Ndone .gt. MAXNRAD)&
                            &call shd_mpi__abort('Error [mpi__computeRad]: MAXNRAD exceeded by Ndone 3.')
                        RadSend(:, Ndone) = RADIANCE(:)
                        IndexSend(Ndone) = BndRadIntInfo(J, Ldowny)
                    else
                        if ((SIDE .ne. Lupx) .and. (SIDE .ne. Lupy)) then
                            write (6, *) 'Bad SIDE 2y', myProc_, mu_out, phi_out, x, y, z, SIDE, xe, ye, ze
                            call shd_mpi__abort('Error [mpi__computeRad]: Bad SIDE 2y, see above.')
                        end if
                        NINVALID(SIDE) = NINVALID(SIDE) + 1
                        if (NINVALID(SIDE) .gt. MAXNRAD)&
                            &call shd_mpi__abort('Error [mpi__computeRad]: MAXNRAD exceeded by NINVALID(SIDE) 3.')
                        BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, TRANSMIT, RADIANCE(:)/)
                        BndRadIntInfo(NINVALID(SIDE), SIDE) = BndRadIntInfo(J, Ldowny)
                    end if
                end do
            end if

            ! See if there are any more invalid radiance rays on all the processors
            intsendbuf = (/NINVALID(Lupx), NINVALID(Lupy)/)
            call mpi_allreduce(intsendbuf, intrecvbuf, 2, MPI_integer, MPI_max,&
                &comm2D_, ierr)
            AllDone = all(intrecvbuf(1:2) .eq. 0)
        end do RAYMPI_LOOP

        if (myProc_ .gt. 0) then
            ! Send the finished radiances to the master processor
            call mpi_isend(Ndone, 1, MPI_integer, 0, 6*K + 7, comm2D_, requests(1), ierr)
            call mpi_isend(RadSend, shdo_obj_%Nstokes*Ndone, MPI_real, 0, 6*K + 8, comm2D_, requests(2), ierr)
            call mpi_isend(IndexSend, Ndone, MPI_integer, 0, 6*K + 9, comm2D_, requests(3), ierr)
            call mpi_waitall(3, requests, status, ierr)
        else
            ! Put the radiances finished by the master processor in the list
            do i = 1, Ndone
                RADOUT(:, IndexSend(i)) = RadSend(:, i)
            end do
            ! Receive the radiances from all the other processors
            irq = 0
            do iproc = 1, Nproc_ - 1
                call mpi_irecv(NradRecv(iproc), 1, MPI_integer, iproc,&
                    &6*K + 7, comm2D_, requests(irq + 1), ierr)
                call mpi_irecv(RadRecv(:, :, iproc), shdo_obj_%Nstokes*MAXNRAD, MPI_real, iproc,&
                    &6*K + 8, comm2D_, requests(irq + 2), ierr)
                call mpi_irecv(IndexRecv(:, iproc), MAXNRAD, MPI_integer, iproc,&
                    &6*K + 9, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end do
            call mpi_waitall(irq, requests, status, ierr)
            do iproc = 1, Nproc_ - 1
                do i = 1, NradRecv(iproc)
                    RADOUT(:, IndexRecv(i, iproc)) = RadRecv(:, i, iproc)
                end do
            end do
        end if

    end do ANGOUT_LOOP

    deallocate(requests, status, RadSend, IndexSend)
    if (myProc_ .eq. 0) deallocate(NradRecv, RadRecv, IndexRecv)
    deallocate(BndRadRealInfo, BndRadIntInfo)
#endif

    return

end subroutine shd_mpi__computeRad

subroutine shd_mpi__computeCamRad(&
    &DELTAM,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &DIRFLUX, FLUXES,&
    &parms_Out, IVIS, VISOUT)

    ! Computes Stokes radiances (output in VISOUT) for the visualization
    !   modes: 1) camera mode,
    !          2) cross track scanning,
    !          3) user-specified sensor mode.

    integer, intent(inout) :: IVIS
    integer, intent(in) :: MAXNBC, NTOPPTS, NBOTPTS
    integer, intent(in) :: BCPTR(MAXNBC, 2)
    logical, intent(in) :: DELTAM
    real, intent(in) :: BCRAD(shdo_obj_%Nstokes, *)
    real, intent(in) :: DIRFLUX(rtm_obj_%Npts), FLUXES(2, rtm_obj_%Npts)
    real, intent(in) :: parms_Out(*)
    real, intent(out) :: VISOUT(shdo_obj_%Nstokes, *)

#if HAS_MPI
    integer :: NSCATANGLE, NSTPHASE, I, J, L, N
    integer :: NL, NS, LINE, SAMP
    logical :: CAMERA_MODE, SENSOR_MODE, VALIDRAD
    real, parameter  :: EPS = 1.0e-5
    real(RD_) :: XE, YE, ZE, TRANSMIT, VISRAD(shdo_obj_%Nstokes), RADIANCE(shdo_obj_%Nstokes)
    real(RD_) :: X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X, Y, Z
    real(RD_) :: THETA0, THETA1, PHIR, PHI0, COSSCAT, MU2, PHI2
    real(RD_) :: COSTH, SINTH, SINTH0, SINTH1
    real(RD_) :: U, V, UP, VP, COSDELPHI, ROTANG, DEGRAD, R
    real(RD_) :: DIST, D, RX, RY, RZ, SCANANG
    integer, parameter :: MAXSCATANG = 721
    real, allocatable :: YLMSUN(:, :), PHASETAB(:, :, :)
    real              :: MEAN, STD1, STD2
    real, allocatable :: AOLP(:)

    ! Domain SIDEs: 1=-X, 2=+X, 3=-Y, 4=+Y
    integer :: MAXNRAD
    integer, parameter :: RECVSIDE(4) = (/2, 1, 4, 3/)
    integer :: SIDE, IBND, NINVALID(4), NINCOMING(4)
    integer :: Ndone, intrecvbuf(4)
    integer :: ierr, irq, iproc, K
    logical :: AllDone
    integer, allocatable :: NradRecv(:)
    integer, allocatable :: IndexSend(:), IndexRecv(:, :)
    integer, allocatable :: BndSendRadIntInfo(:, :), BndRecvRadIntInfo(:, :)
    real, allocatable    :: RadSend(:, :), RadRecv(:, :, :)
    real, allocatable    :: BndSendRadRealInfo(:, :, :), BndRecvRadRealInfo(:, :, :)

    type(MPI_request), allocatable :: requests(:)
    type(MPI_status), allocatable :: status(:)

    allocate(YLMSUN(prp_obj_%NlegSet, shdo_obj_%Nlm))

    if (src_obj_%doSolar) then
        call shd_shdo__computeYlm(.true., src_obj_%avgMu, src_obj_%avgAz, shdo_obj_%maxNl, shdo_obj_%maxNm, prp_obj_%NlegSet, YLMSUN)
        if (DELTAM .and. prp_obj_%Npha .gt. 0) then
            NSCATANGLE = max(36, min(MAXSCATANG, 2*prp_obj_%Nleg))
            NSTPHASE = min(prp_obj_%NlegSet, 2)
            allocate(PHASETAB(NSTPHASE, prp_obj_%Npha, NSCATANGLE))
            call shd_scat__computePhase(NSCATANGLE, prp_obj_%Npha, NSTPHASE,&
                &shdo_obj_%Nstokes, shdo_obj_%maxNl, prp_obj_%NlegSet, prp_obj_%Nleg, leg_Rtm_, PHASETAB)
        end if
    end if

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(src_obj_%type, src_obj_%skyRad, src_obj_%waveNum, src_obj_%avgWaveLen,&
        &src_obj_%unit, NTOPPTS, shdo_obj_%Nstokes, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces.
    ! Compute the upwelling bottom radiances using the downwelling fluxes.
    if (sfc_obj_%type .eq. 'FL') then
        call shd_init__makeBotRadUL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, src_obj_%type,&
            &src_obj_%waveNum, src_obj_%avgWaveLen, src_obj_%unit, shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
    else if (sfc_obj_%type .eq. 'VL') then
        call shd_init__makeBotRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, src_obj_%type,&
            &shdo_obj_%Nstokes, BCRAD(1, 1 + NTOPPTS))
    end if

    CAMERA_MODE = nint(parms_Out(1)) .eq. 1
    SENSOR_MODE = nint(parms_Out(1)) .eq. 3
    if (CAMERA_MODE) then
        NL = nint(parms_Out(10))
        NS = nint(parms_Out(11))
    else if (SENSOR_MODE) then
        NL = 1
        NS = parms_Out(2)
    else
        NL = 1 + int(sqrt((parms_Out(4) - parms_Out(7))**2&
            &+ (parms_Out(5) - parms_Out(8))**2&
            &+ (parms_Out(6) - parms_Out(9))**2)/parms_Out(10))
        NS = 1 + int(abs(parms_Out(12) - parms_Out(11))/parms_Out(13))
    end if

    DEGRAD = DPI_/180.0d0

    ! Allocate sending, receiving, etc. arrays for MPI
    MAXNRAD = NL*NS
    allocate(RadSend(shdo_obj_%Nstokes, MAXNRAD), IndexSend(MAXNRAD))
    if (myProc_ .eq. 0) allocate(NradRecv(Nproc_ - 1),&
        &RadRecv(shdo_obj_%Nstokes, MAXNRAD, Nproc_ - 1), IndexRecv(MAXNRAD, Nproc_ - 1))
    allocate(BndSendRadRealInfo(6 + shdo_obj_%Nstokes, MAXNRAD, 4), BndSendRadIntInfo(MAXNRAD, 4))
    allocate(BndRecvRadRealInfo(6 + shdo_obj_%Nstokes, MAXNRAD, 4), BndRecvRadIntInfo(MAXNRAD, 4))
    allocate(requests(max(24, 3*Nproc_)), status(max(24, 3*Nproc_)))

    NINVALID(:) = 0
    Ndone = 0

    ! Loop over pixels in image
    do LINE = 1, NL
        do SAMP = 1, NS
            IVIS = IVIS + 1

            if (CAMERA_MODE) then

                ! Camera mode:
                !   1, bytes, scale, X,Y,Z, theta, phi, rotang, NL, NS, delline, delsamp
                ! Use spherical trig to find the pixel direction (MURAY,PHIRAY)
                ! from the camera center (THETA0,PHI0) and the relative pixel angles (U,V).
                UP = (SAMP - NS/2 - 1)*parms_Out(13)*DEGRAD
                VP = (LINE - NL/2 - 1)*parms_Out(12)*DEGRAD
                ROTANG = parms_Out(9)*DEGRAD
                if (ROTANG .eq. 0.0) then
                    U = UP
                    V = VP
                else
                    U = cos(ROTANG)*UP - sin(ROTANG)*VP
                    V = sin(ROTANG)*UP + cos(ROTANG)*VP
                end if
                THETA0 = DEGRAD*parms_Out(7)
                PHI0 = DEGRAD*parms_Out(8)
                THETA1 = THETA0 + V
                if (V .eq. 0.0) then
                    COSTH = cos(U)*cos(THETA0)
                else
                    COSTH = cos(U)*(sin(THETA1)*cos(V) - sin(THETA0))/sin(V)
                    COSTH = min(+1.0d0, max(-1.0d0, COSTH))
                end if
                SINTH = sqrt(1 - COSTH**2)
                SINTH0 = sin(THETA0)
                SINTH1 = sin(THETA1)
                if (abs(SINTH) .lt. accSolv_) then
                    PHIR = 0.0
                else
                    if ((abs(SINTH0) .lt. accSolv_) .and. (abs(SINTH1) .le. accSolv_)) then
                        COSDELPHI = 0.0d0
                    else if (abs(SINTH1) .gt. accSolv_) then
                        COSDELPHI = (cos(U) - COSTH*cos(THETA1))/(SINTH1*SINTH)
                    else if (abs(SINTH0) .gt. accSolv_) then
                        COSDELPHI = (cos(U)*cos(V) - COSTH*cos(THETA0))/(SINTH0*SINTH)
                    end if
                    COSDELPHI = min(+1.0d0, max(-1.0d0, COSDELPHI))
                    if (U .ge. 0.0) then
                        PHIR = PHI0 - acos(COSDELPHI)
                    else
                        PHIR = PHI0 + acos(COSDELPHI)
                    end if
                end if

                X0 = parms_Out(4)
                Y0 = parms_Out(5)
                Z0 = parms_Out(6)

            else if (SENSOR_MODE) then

                if (sen_obj_%fname(1:4) .ne. 'NONE') then
                    X0    =  view_Sen_(1, SAMP)
                    Y0    =  view_Sen_(2, SAMP)
                    Z0    =  view_Sen_(3, SAMP)
                    COSTH =  view_Sen_(4, SAMP)
                    PHIR  =  view_Sen_(5, SAMP)
                else
                    stop 'Error [mpi__computeCamRad]: Need sen_obj_ to proceed.'
                end if

            else

                ! Cross track scanning in the vertical plane:
                !  2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, delscan
                ! Start and end of scan angle (scan1, scan2) are +/- relative
                ! to nadir, and positive is on right side.
                X1 = parms_Out(4)
                Y1 = parms_Out(5)
                Z1 = parms_Out(6)
                X2 = parms_Out(7)
                Y2 = parms_Out(8)
                Z2 = parms_Out(9)
                DIST = sqrt((X2 - X1)**2 + (Y2 - Y1)**2 + (Z2 - Z1)**2)
                RX = (X2 - X1)/DIST
                RY = (Y2 - Y1)/DIST
                RZ = (Z2 - Z1)/DIST
                D = (NL - LINE)*parms_Out(10)
                X0 = X1 + D*RX
                Y0 = Y1 + D*RY
                Z0 = Z1 + D*RZ
                SCANANG = DEGRAD*(parms_Out(11) + (SAMP - 1)*parms_Out(13))
                COSTH = cos(DPI_ - abs(SCANANG))
                if (SCANANG .gt. 0.0) then
                    PHIR = atan2(-RX, RY)
                else
                    PHIR = atan2(RX, -RY)
                end if

            end if

            !print *, IVIS, X0, Y0, Z0, 180.0/DPI_*COSTH, 180.0/DPI_*PHIR

            ! Extrapolate ray to domain top if above
            if (Z0 .gt. rtm_obj_%z(rtm_obj_%Nz)) then
                if (COSTH .ge. 0.0) then
                    Ndone = Ndone + 1
                    RadSend(:, Ndone) = 0.0
                    IndexSend(Ndone) = IVIS
                    cycle
                end if
                R = (rtm_obj_%z(rtm_obj_%Nz) - Z0)/COSTH
                X0 = X0 + R*sqrt(1 - COSTH**2)*cos(PHIR)
                Y0 = Y0 + R*sqrt(1 - COSTH**2)*sin(PHIR)
                Z0 = rtm_obj_%z(rtm_obj_%Nz)
            else if (Z0 .lt. rtm_obj_%z(1)) then
                write (6, *) 'Error [mpi__visualizeRad]: Level below domain.'
                stop
            end if

            ! Relocate the starting point in the overall periodic domain
            ! (we know we have a periodic domain if doing MPI)
            if (X0 .ge. 0.0) then
                X0 = mod(X0, dble(prp_obj_%xDomain))
            else
                X0 = mod(X0 + prp_obj_%xDomain*int(1.0 - X0/prp_obj_%xDomain), dble(prp_obj_%xDomain))
            end if
            if (Y0 .ge. 0.0) then
                Y0 = mod(Y0, dble(prp_obj_%yDomain))
            else
                Y0 = mod(Y0 + prp_obj_%yDomain*int(1.0 - Y0/prp_obj_%yDomain), dble(prp_obj_%yDomain))
            end if

            ! Start the radiance calculation for those rays starting in this subdomain
            if ((X0 .ge. rtm_obj_%x(1)) .and. (X0 .le. rtm_obj_%x(rtm_obj_%Nx)).and.&
                (Y0 .ge. rtm_obj_%y(1)) .and. (Y0 .le. rtm_obj_%y(rtm_obj_%Ny))) then

                ! MU2,PHI2 is radiance travel direction (opposite the camera viewing direction)
                MU2 = -COSTH
                PHI2 = PHIR + DPI_

                ! Integrate the extinction and source function along this ray
                ! to calculate the Stokes radiance vector for this pixel
                TRANSMIT = 1.0d0; RADIANCE(:) = 0.0d0
                call shd_src__integrateSrcOneRay(shdo_obj_%Nstokes,&
                    &NSTPHASE, NSCATANGLE, PHASETAB,&
                    &shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                    &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt,&
                    &DELTAM, src_obj_%type, src_obj_%avgWaveLen, src_obj_%avgMu, src_obj_%avgAz,&
                    &DIRFLUX, isrcJ_Shdo_, srcJ_Shdo_, YLMSUN,&
                    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                    &MU2, PHI2, X0, Y0, Z0,&
                    &XE, YE, ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

                ! If we got a radiance then store it in RadSend otherwise save
                !  the information needed to pass to other processors
                if (VALIDRAD) then
                    Ndone = Ndone + 1
                    if (Ndone .gt. MAXNRAD)&
                        &call shd_mpi__abort('Error [mpi__visualizeRad]: MAXNRAD exceeded by Ndone 1')
                    RadSend(:, Ndone) = RADIANCE(:)
                    IndexSend(Ndone) = IVIS
                else
                    if (SIDE .lt. 1 .or. SIDE .gt. 4) call shd_mpi__abort('Error [mpi__visualizeRad]: invalid SIDE 1')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) .gt. MAXNRAD)&
                        &call shd_mpi__abort('Error [mpi__visualizeRad]: MAXNRAD exceeded by NINVALID(SIDE) 1')
                    BndSendRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, MU2, PHI2, TRANSMIT, RADIANCE(:)/)
                    BndSendRadIntInfo(NINVALID(SIDE), SIDE) = IVIS
                end if
            end if
        end do
    end do

    ! Cycle over the passing of info between processors and doing more
    !  partial radiance integrations until all rays have been finished
    AllDone = .false.
    do while (.not. AllDone)
        irq = 0
        ! Do the sends for all four sides of the subdomain.  The number of
        ! incomplete rays and the real and integer info needed about the
        ! radiance integrations is sent.
        do SIDE = 1, 4
            if ((isMpiX_ .and. (SIDE .eq. 1 .or. SIDE .eq. 2)) .or.&
                (isMpiY_ .and. (SIDE .eq. 3 .or. SIDE .eq. 4))) then
                K = 3*SIDE
                call mpi_isend(NINVALID(SIDE), 1, MPI_integer,&
                    &iProcNeigh_(SIDE), K + 1, comm2D_, requests(irq + 1), ierr)
                call mpi_isend(BndSendRadRealInfo(:, :, SIDE), (6 + shdo_obj_%Nstokes)*NINVALID(SIDE), MPI_double_precision,&
                    &iProcNeigh_(SIDE), K + 2, comm2D_, requests(irq + 2), ierr)
                call mpi_isend(BndSendRadIntInfo(:, SIDE), NINVALID(SIDE), MPI_integer,&
                    &iProcNeigh_(SIDE), K + 3, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
        end do
        ! Do the receives for all four sides of the subdomain.
        do SIDE = 1, 4
            if ((isMpiX_ .and. (SIDE .eq. 1 .or. SIDE .eq. 2)) .or.&
                (isMpiY_ .and. (SIDE .eq. 3 .or. SIDE .eq. 4))) then
                K = 3*RECVSIDE(SIDE)
                call mpi_irecv(NINCOMING(SIDE), 1, MPI_integer,&
                    &iProcNeigh_(SIDE), K + 1, comm2D_, requests(irq + 1), ierr)
                call mpi_irecv(BndRecvRadRealInfo(:, :, SIDE), (6 + shdo_obj_%Nstokes)*MAXNRAD, MPI_double_precision,&
                    &iProcNeigh_(SIDE), K + 2, comm2D_, requests(irq + 2), ierr)
                call mpi_irecv(BndRecvRadIntInfo(:, SIDE), MAXNRAD, MPI_integer,&
                    &iProcNeigh_(SIDE), K + 3, comm2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
        end do
        call mpi_waitall(irq, requests, status, ierr)

        ! Continue the backwards ray integrations for the four boundaries
        NINVALID(:) = 0
        do IBND = 1, 4
            do J = 1, NINCOMING(IBND)
                if (IBND .eq. 1) BndRecvRadRealInfo(1, J, IBND) = rtm_obj_%x(1)
                if (IBND .eq. 2) BndRecvRadRealInfo(1, J, IBND) = rtm_obj_%x(rtm_obj_%Nx)
                if (IBND .eq. 3) BndRecvRadRealInfo(2, J, IBND) = rtm_obj_%y(1)
                if (IBND .eq. 4) BndRecvRadRealInfo(2, J, IBND) = rtm_obj_%y(rtm_obj_%Ny)
                X = BndRecvRadRealInfo(1, J, IBND)
                Y = BndRecvRadRealInfo(2, J, IBND)
                Z = BndRecvRadRealInfo(3, J, IBND)
                MU2 = BndRecvRadRealInfo(4, J, IBND)
                PHI2 = BndRecvRadRealInfo(5, J, IBND)
                TRANSMIT = BndRecvRadRealInfo(6, J, IBND)
                RADIANCE(:) = BndRecvRadRealInfo(7 : 6+shdo_obj_%Nstokes, J, IBND)
                call shd_src__integrateSrcOneRay(shdo_obj_%Nstokes,&
                    &NSTPHASE, NSCATANGLE, PHASETAB,&
                    &shdo_obj_%maxNl, shdo_obj_%maxNm, shdo_obj_%Nlm,&
                    &shdo_obj_%Nmu, shdo_obj_%maxNphiA, shdo_obj_%NphiA, shdo_obj_%mu, shdo_obj_%phi, shdo_obj_%wgt,&
                    &DELTAM, src_obj_%type, src_obj_%avgWaveLen, src_obj_%avgMu, src_obj_%avgAz,&
                    &DIRFLUX, isrcJ_Shdo_, srcJ_Shdo_, YLMSUN,&
                    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                    &MU2, PHI2, X, Y, Z,&
                    &XE, YE, ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

                ! If we got a valid radiance then store it in RadSend, otherwise
                ! save information needed to pass to neighboring processors
                if (VALIDRAD) then
                    Ndone = Ndone + 1
                    if (Ndone .gt. MAXNRAD)&
                        &call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 2')
                    RadSend(:, Ndone) = RADIANCE(:)
                    IndexSend(Ndone) = BndRecvRadIntInfo(J, IBND)
                else
                    if (SIDE .lt. 1 .or. SIDE .gt. 4) call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: invalid SIDE 2')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) .gt. MAXNRAD)&
                        &call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 2')
                    BndSendRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, MU2, PHI2, TRANSMIT, RADIANCE(:)/)
                    BndSendRadIntInfo(NINVALID(SIDE), SIDE) = BndRecvRadIntInfo(J, IBND)
                end if
            end do
        end do

        ! See if there are any more invalid radiance rays on all the processors
        call mpi_allreduce(NINVALID, intrecvbuf, 4, MPI_integer, MPI_max,&
            &comm2D_, ierr)
        AllDone = all(intrecvbuf(1:4) .eq. 0)
    end do
    !print *,'Done with ray tracing.', Ndone

    if (myProc_ .gt. 0) then
        ! Send the finished radiances to the master processor
        call mpi_isend(Ndone, 1, MPI_integer, 0, 21, comm2D_, requests(1), ierr)
        call mpi_isend(RadSend, shdo_obj_%Nstokes*Ndone, MPI_real, 0, 22, comm2D_, requests(2), ierr)
        call mpi_isend(IndexSend, Ndone, MPI_integer, 0, 23, comm2D_, requests(3), ierr)
        call mpi_waitall(3, requests, status, ierr)
    else
        ! Put the radiances finished by the master processor in the list
        do i = 1, Ndone
            VISOUT(:, IndexSend(i)) = RadSend(:, i)
        end do
        ! Receive the radiances from all the other processors
        irq = 0
        do iproc = 1, Nproc_ - 1
            call mpi_irecv(NradRecv(iproc), 1, MPI_integer, iproc,&
                &21, comm2D_, requests(irq + 1), ierr)
            call mpi_irecv(RadRecv(:, :, iproc), shdo_obj_%Nstokes*MAXNRAD, MPI_real, iproc,&
                &22, comm2D_, requests(irq + 2), ierr)
            call mpi_irecv(IndexRecv(:, iproc), MAXNRAD, MPI_integer, iproc,&
                &23, comm2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end do
        call mpi_waitall(irq, requests, status, ierr)
        do iproc = 1, Nproc_ - 1
            do i = 1, NradRecv(iproc)
                VISOUT(:, IndexRecv(i, iproc)) = RadRecv(:, i, iproc)
            end do
        end do
    end if

    deallocate(requests, status, RadSend, IndexSend)
    if (myProc_ .eq. 0) deallocate(NradRecv, RadRecv, IndexRecv)
    deallocate(BndSendRadRealInfo, BndSendRadIntInfo)
    deallocate(BndRecvRadRealInfo, BndRecvRadIntInfo)
    deallocate(YLMSUN)
    if (allocated(PHASETAB)) deallocate(PHASETAB)

    ! On the master processor convert VISOUT(2:,*) from Stokes parameters
    ! to DOLP, AOLP, DOCP
    if (myProc_ .eq. 0) then
        if (isPol_) then
            N = NL*NS
            do I = IVIS - N + 1, IVIS
                VISRAD(:) = VISOUT(:, I)
                if (VISRAD(1) .gt. 0.0) then
                    ! Output degree (0 to 1) and angle (-180 to 180) of linear polarization
                    VISOUT(2, I) = sqrt(VISRAD(2)**2 + VISRAD(3)**2)/VISRAD(1)
                    VISOUT(3, I) = (180.0/DPI_)*0.5*atan2(VISRAD(3), VISRAD(2))
                    if (isPol4_) then
                        ! Output degree of circular polarization (-1 to 1)
                        VISOUT(4, I) = VISRAD(4)/VISRAD(1)
                    end if
                end if
            end do
            ! Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
            allocate(AOLP(N))
            AOLP(:) = VISOUT(3, IVIS - N + 1:IVIS)
            MEAN = sum(AOLP(:))/N
            STD1 = sqrt(sum((AOLP(:) - MEAN)**2)/N)
            where (AOLP(:) .lt. 0.0)
                AOLP(:) = AOLP(:) + 180.0
            end where
            MEAN = sum(AOLP(:))/N
            STD2 = sqrt(sum((AOLP(:) - MEAN)**2)/N)
            if (STD2 .lt. STD1) then
                VISOUT(3, IVIS - N + 1:IVIS) = AOLP(:)
            end if
            deallocate(AOLP)
        end if
    end if
#endif

    return

end subroutine shd_mpi__computeCamRad

subroutine shd_mpi__calAccelSolCrit(DOACCEL, DELJDOT, DELJOLD, DELJNEW, JNORM,&
    &ACCELPAR, SOLCRIT)

    ! Does an mpi_allreduce to sum the delta source function vector dot
    ! products over all the processors, and from these calculates the
    ! full-domain-wide acceleration parameter and solution criterion.

    logical, intent(in) :: DOACCEL
    real, intent(inout) :: DELJDOT, DELJOLD, DELJNEW, JNORM
    real, intent(out) :: ACCELPAR, SOLCRIT
    real    ::  R, THETA
    real, save :: A = 0.0

#if HAS_MPI
    integer :: ierr
    real :: sendbuf(4), recvbuf(4)

    if (Nproc_ .gt. 1) then
        ! Call mpi_allreduce to sum up the four dot products
        sendbuf(:) = (/DELJDOT, DELJOLD, DELJNEW, JNORM/)
        call mpi_allreduce(sendbuf, recvbuf, 4, MPI_real, MPI_sum, comm2D_, ierr)
        if (ierr .ne. MPI_success) call shd_mpi__abort('Error [mpi__calAccelSolCrit]: mpi_allreduce error.')
        DELJDOT = recvbuf(1); DELJOLD = recvbuf(2)
        DELJNEW = recvbuf(3); JNORM = recvbuf(4)
    end if
#endif

    ! Accelerate if desired, didn't last time, and things are converging.
    if (DOACCEL .and. A .eq. 0.0 .and. DELJNEW .lt. DELJOLD) then
        ! Compute the acceleration extrapolation factor and apply it.
        R = sqrt(DELJNEW/DELJOLD)
        THETA = acos(DELJDOT/sqrt(DELJOLD*DELJNEW))
        A = (1 - R*cos(THETA) + R**(1.0 + PIH_/THETA))&
            /(1 + R**2 - 2*R*cos(THETA)) - 1.0
        A = min(10.0, max(0.0, A))
        ! WRITE (6,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
    else
        A = 0.0
    end if
    ACCELPAR = A

    if (JNORM .gt. 0.0) then
        SOLCRIT = sqrt(DELJNEW/JNORM)
    else if (DELJNEW .eq. 0.0) then
        SOLCRIT = 0.0
    end if

    return

end subroutine shd_mpi__calAccelSolCrit

subroutine shd_mpi__end(runName, rtm_obj, prpAll_obj)

    ! Shuts down MPI at the end, but mainly finds the number of gridpoints
    ! around each property grid point column, sends that array to the master
    ! process, and the master process write out the runName//"_load_balance.out"
    ! file.

    type(shd_type_rtm), intent(in) :: rtm_obj
    type(shd_type_prp), intent(in) :: prpAll_obj

    character(len=*), intent(in) :: runName

#if HAS_MPI
    integer :: I, IX, IY, maxN, Sx, Sy, Ex, Ey, Nx, Ny
    integer :: iproc, ierr, irq
    integer, allocatable :: NUMGRID(:, :), NUMGRIDALL(:, :), NORM(:, :), ixyAll(:, :, :)
    integer, allocatable :: recvbuf(:, :), buf(:, :)
    real :: dx_prp, dy_prp
    real :: ElapsedTime
    type(MPI_request), allocatable :: requests(:)
    type(MPI_status), allocatable :: status(:)

    dx_prp = real(prp_obj_%xDomain)/real(prpAll_obj%Nx)
    dy_prp = real(prp_obj_%yDomain)/real(prpAll_obj%Ny)

    if (loadBalance_) then
        ! Calculate the number of grid points around each property grid column
        allocate(NUMGRID(prp_obj_%Nx, prp_obj_%Ny))
        NUMGRID(:, :) = 0
        do I = 1, rtm_obj%Npts
            IX = nint((xyz_Rtm_(1, I) - prp_obj_%xStart)/dx_prp) + 1
            IY = nint((xyz_Rtm_(2, I) - prp_obj_%yStart)/dy_prp) + 1
            if ((IX .le. prp_obj_%Nx) .and. (IY .le. prp_obj_%Ny)) then
                NUMGRID(IX, IY) = NUMGRID(IX, IY) + 1
            end if
        end do

        if (Nproc_ .gt. 1) then
            ! Have the master process collect the property grid locations of all
            allocate(ixyAll(2, 2, 0:Nproc_ - 1))
            call mpi_gather(ixyPrp_, 4, MPI_integer, ixyAll, 4, MPI_integer, 0, comm2D_, ierr)

            ! Send the NUMGRIDs to the master processor
            allocate(requests(Nproc_), status(Nproc_))
            irq = 0
            if (myProc_ .gt. 0) then
                irq = irq + 1
                call mpi_isend(NUMGRID, prp_obj_%Nx*prp_obj_%Ny, MPI_integer, 0, 777, comm2D_,&
                    &requests(irq), ierr)
            else
                maxN = (maxval(ixyAll(2, 1, :) - ixyAll(1, 1, :)) + 1)*(maxval(ixyAll(2, 2, :) - ixyAll(1, 2, :)) + 1)
                allocate(recvbuf(maxN, 1:Nproc_ - 1))
                do iproc = 1, Nproc_ - 1
                    irq = irq + 1
                    call mpi_irecv(recvbuf(:, iproc), maxN, MPI_integer, iproc, 777,&
                        &comm2D_, requests(irq), ierr)
                end do
            end if
            call mpi_waitall(irq, requests, status, ierr)

            ! Put the NUMGRIDS into the full domain array on the master processor.
            !   Use a normalization counter (NORM) to deal with grid points that
            !   get contribution from two processors.
            if (myProc_ .eq. 0) then
                allocate(NUMGRIDALL(prpAll_obj%Nx, prpAll_obj%Ny), NORM(prpAll_obj%Nx, prpAll_obj%Ny))
                NUMGRIDALL(:, :) = 0; NORM(:, :) = 0
                Sx = ixyAll(1, 1, 0); Ex = ixyAll(2, 1, 0); Nx = Ex - Sx + 1
                Sy = ixyAll(1, 2, 0); Ey = ixyAll(2, 2, 0); Ny = Ey - Sy + 1
                NUMGRIDALL(Sx:Ex, Sy:Ey) = NUMGRID(:, :)
                NORM(Sx:Ex, Sy:Ey) = reshape((/1/), (/Nx, Ny/), (/1/))
                do iproc = 1, Nproc_ - 1
                    Sx = ixyAll(1, 1, iproc); Ex = ixyAll(2, 1, iproc); Nx = Ex - Sx + 1
                    Sy = ixyAll(1, 2, iproc); Ey = ixyAll(2, 2, iproc); Ny = Ey - Sy + 1
                    allocate(buf(Nx, Ny))
                    buf(:, :) = reshape(recvbuf(1:Nx*Ny, iproc), (/Nx, Ny/))
                    NUMGRIDALL(Sx:Ex - 1, Sy:Ey - 1) = NUMGRIDALL(Sx:Ex - 1, Sy:Ey - 1) + buf(1:Nx - 1, 1:Ny - 1)
                    NORM(Sx:Ex - 1, Sy:Ey - 1) = NORM(Sx:Ex - 1, Sy:Ey - 1) + reshape((/1/), (/Nx - 1, Ny - 1/), (/1/))
                    IX = mod(Ex - 1, prpAll_obj%Nx) + 1; IY = mod(Ey - 1, prpAll_obj%Ny) + 1
                    NUMGRIDALL(IX, Sy:Ey - 1) = NUMGRIDALL(IX, Sy:Ey - 1) + buf(Nx, 1:Ny - 1)
                    NORM(IX, Sy:Ey - 1) = NORM(IX, Sy:Ey - 1) + reshape((/1/), (/Ny - 1/), (/1/))
                    NUMGRIDALL(Sx:Ex - 1, IY) = NUMGRIDALL(Sx:Ex - 1, IY) + buf(1:Nx - 1, Ny)
                    NORM(Sx:Ex - 1, IY) = NORM(Sx:Ex - 1, IY) + reshape((/1/), (/Nx - 1/), (/1/))
                    NUMGRIDALL(IX, IY) = NUMGRIDALL(IX, IY) + buf(Nx, Ny)
                    NORM(IX, IY) = NORM(IX, IY) + 1

                    ! checked
                    deallocate(buf)
                end do

                ! checked
                deallocate(recvbuf)
                NUMGRIDALL(:, :) = NUMGRIDALL(:, :)/NORM(:, :)

                ! checked
                deallocate(NORM)
            end if

            ! checked
            deallocate(requests, status, ixyAll)

        else
            allocate(NUMGRIDALL(prpAll_obj%Nx, prpAll_obj%Ny))
            NUMGRIDALL(1:prp_obj_%Nx, 1:prp_obj_%Ny) = NUMGRID(:, :)
        end if

        ! checked
        deallocate(NUMGRID)

        if (myProc_ .eq. 0) then
            ! Write out the array to the load balancing file
            open (UNIT=10, FILE=trim(runName)//'_load_balance.out', STATUS='UNKNOWN')
            write (10, '(A)') '! SHDOM grid point load balancing file for property file:'
            write (10, '(A)') prp_obj_%fname
            write (10, '(2(1X,I4),A)') prpAll_obj%Nx, prpAll_obj%Ny, ' ! property grid size (Nx,Ny)'
            do IX = 1, prpAll_obj%Nx
                do IY = 1, prpAll_obj%Ny
                    write (10, '(I4,1X,I4,1X,I5)') IX, IY, NUMGRIDALL(IX, IY)
                end do
            end do
            close (10)

            ! checked
            deallocate(NUMGRIDALL)
        end if
    end if

    ElapsedTime = mpi_wtime() - startTime_

    if (ElapsedTime .lt. 60.0) then
        if (verbose_ .ge. 5) write (6, '(A,F0.2,A)') 'Message [mpi__end]: Total time used: ', ElapsedTime, ' sec'
    else
        if (verbose_ .ge. 5) write (6, '(A,F0.2,A)') 'Message [mpi__end]: Total time used: ', ElapsedTime/60.0, ' min'
    end if

    if (verbose_ .ge. 5) write (6, '(A,"<",A,">")') 'Message [mpi__end]: Finished with SHDOM run: ', trim(adjustl(runName))

    ! Shut down MPI
    if (allocated(BxyPtr_)) deallocate(BxyPtr_)
    if (allocated(iBxyPtr_)) deallocate(iBxyPtr_)
    if (allocated(BxyPos_)) deallocate(BxyPos_)
    call mpi_finalize(ierr)

#endif

    return

end subroutine shd_mpi__end

subroutine shd_mpi__abort(ERRSTR)

    character(len=*), intent(in) :: ERRSTR
    integer :: ierr, errorcode = 7

    write (6, *) trim(adjustl(ERRSTR))

#if HAS_MPI
! note: this if clause has <else>
    call MPI_ABORT(comm2D_, errorcode, ierr)
#else
    stop
#endif

end subroutine shd_mpi__abort

#if HAS_MPI
subroutine shd_mpi__optimizeProcDomain(ProcDims, ProcCoord,&
    &runName, prpAll_obj, prp_obj)

    ! Calculates the location of this processor (subdomain) in the property
    ! grid.  The starting and ending property grid point indices for X and Y
    ! are returned in ixyPrp_.  Each processor sends its location in the
    ! 2D array of subdomains to the master process and the master process
    ! reads the runName//"_load_balance.inp" file (if there is one) and does
    ! the optimization.  If the "shdom_proc_times.inp" file of relative
    ! processor times exists then the optimization uses this information on
    ! heterogeneous processors to try to give more work to faster processors.

    type(shd_type_prp), intent(in) :: prpAll_obj, prp_obj

    integer, intent(in) :: ProcDims(2), ProcCoord(2)
    character(len=*), intent(in) :: runName
    integer :: NPXTI, NPYTI
    integer :: IX, IY, JX, JY, Nxd, Nyd, N, I, ierr
    integer :: BndLines(ProcDims(1) + ProcDims(2) - 2)
    integer :: IDX(ProcDims(1) + 1), IDY(ProcDims(2) + 1)
    logical :: DoLoadBalancing
    integer, parameter :: NCYCLE = 30 ! number of annealing cooling cycles
    integer, parameter :: NiterFac = 100 ! number of iterations per length of optimization vector (BndLines)
    real, parameter :: Temp0 = 0.1 ! initial temperature (units of change in objective function "energy")
    real, parameter :: Tempfac = 4.0/NCYCLE ! exponential cooling rate
    integer :: icycle, ITER, NITER
    integer :: TrialLines(ProcDims(1) + ProcDims(2) - 2), CurLines(ProcDims(1) + ProcDims(2) - 2)
    integer :: iseed = -3
    real    :: Temp, DelLines, u, shd_mpi__ran
    real    :: loadObject, Energy, NewEnergy, LowestEnergy, DeltaE, Boltz
    integer, allocatable :: ProcCoordAll(:, :), Ngrid(:, :)
    real, allocatable :: wtproc(:)
    character(len=charMaxLen_) :: fileInp_prp0

    ! Initialize the subdomain boundary lines to evenly spaced
    Nxd = ProcDims(1); Nyd = ProcDims(2); N = Nxd + Nyd - 2
    BndLines(   1:Nxd - 1) = (/(nint(I*float(prpAll_obj%Nx)/Nxd) + 1, I=1, Nxd - 1)/)
    BndLines( Nxd:      N) = (/(nint(I*float(prpAll_obj%Ny)/Nyd) + 1, I=1, Nyd - 1)/)

    if (loadBalance_) then
        DoLoadBalancing = .false.
        if (myProc_ .eq. 0) then
            ! See if we can read the load balancing file and header matches
            open (UNIT=10, FILE=trim(runName)//'_load_balance.inp', STATUS='OLD', IOSTAT=ierr)
            if (ierr .eq. 0) then
                read (10, *)
                read (10, '(A)') fileInp_prp0
                read (10, *) NPXTI, NPYTI
                if ((NPXTI .ne. prpAll_obj%Nx) .or. (NPYTI .ne. prpAll_obj%Ny) .or. (fileInp_prp0 .ne. prp_obj%fname)) then
                    write (6, '(A,A,A)') 'Load balancing file ', trim(runName),&
                        &'_load_balance.inp is not compatible with property file:'
                    write (6, *) prp_obj%fname
                    write (6, *)
                else
                    DoLoadBalancing = .true.
                end if
                close (10)
            end if
        end if

        ! Have the master process collect the locations of all
        allocate(ProcCoordAll(2, 0:Nproc_ - 1))
        call mpi_gather(ProcCoord, 2, MPI_integer, ProcCoordAll, 2, MPI_integer, 0, comm2D_, ierr)

        if (myProc_ .eq. 0 .and. DoLoadBalancing) then
            ! Read the load balancing file to get array of number of grid points
            allocate(Ngrid(prpAll_obj%Nx + 1, prpAll_obj%Ny + 1))
            open (UNIT=10, FILE=trim(runName)//'_load_balance.inp', STATUS='OLD')
            read (10, *)
            read (10, *)
            read (10, *)
            do IX = 1, prpAll_obj%Nx
                do IY = 1, prpAll_obj%Ny
                    read (10, *) JX, JY, Ngrid(IX, IY)
                end do
            end do
            close (10)
            Ngrid(prpAll_obj%Nx + 1, :) = Ngrid(1, :)
            Ngrid(:, prpAll_obj%Ny + 1) = Ngrid(:, 1)
            Ngrid(prpAll_obj%Nx + 1, prpAll_obj%Ny+ 1) = Ngrid(1, 1)

            ! Read the relative processor times if the "shdom_proc_times.inp" file exists
            allocate(wtproc(Nproc_))
            open (UNIT=11, FILE='shdom_proc_times.inp', STATUS='OLD', IOSTAT=ierr)
            if (ierr .eq. 0) then
                do i = 1, Nproc_
                    read (11, *) wtproc(i)
                end do
                close (11)
                wtproc(:) = wtproc(:)*Nproc_/sum(wtproc(:))
            else
                wtproc(:) = 1.0
            end if

            ! Do the simulated annealing to minimize the objective function
            u = shd_mpi__ran(iseed)
            NITER = NiterFac*N
            Energy = loadObject(BndLines, Nxd, Nyd, ProcCoordAll, wtproc, Ngrid, prpAll_obj)
            LowestEnergy = Energy
            CurLines(:) = BndLines(:)
            ! print *, 'Annealing cycles: Cycle Temp Energy  BndLines'
            ! print '(I3,2(1X,F9.6),10(1X,I3))', 0,Temp0,Energy, BndLines(:)

            ! Loop through the cooling cycles
            do icycle = 1, NCYCLE
                Temp = Temp0*exp(-Tempfac*(icycle - 1)) ! lower the temp
                ! Adjust the maximum distance a boundary line is perturbed with the cycle number
                u = float(icycle - 1)/(NCYCLE - 1)
                DelLines = (1 - u)*0.5*max(prpAll_obj%Nx/Nxd, prpAll_obj%Ny/Nyd) + u*1.0

                ! Do NITER iterations per temp cycle
                do ITER = 1, NITER
                    ! Perturb the boundary lines
                    TrialLines(:) = CurLines(:)
                    do while (all(TrialLines(:) .eq. CurLines(:)))
                        do I = 1, N
                            TrialLines(I) = nint(CurLines(I) + (2*shd_mpi__ran(iseed) - 1.0)*DelLines)
                        end do
                    end do
                    ! Calculate the energy of these new boundary lines
                    NewEnergy = loadObject(TrialLines, Nxd, Nyd, ProcCoordAll, wtproc, Ngrid, prpAll_obj)
                    DeltaE = NewEnergy - Energy
                    if (DeltaE .lt. 0.0) then ! If lower energy accept change
                        CurLines(:) = TrialLines(:)
                        Energy = NewEnergy
                    else !  otherwise sometimes accept it anyway
                        Boltz = exp(-DeltaE/Temp)
                        if (shd_mpi__ran(iseed) .lt. Boltz) then
                            CurLines(:) = TrialLines(:)
                            Energy = NewEnergy
                        end if
                    end if
                    if (Energy .lt. LowestEnergy) then
                        BndLines(:) = CurLines(:)
                        LowestEnergy = Energy
                    end if
                end do
                ! print '(I3,2(1X,F9.6),10(1X,I3))', CYCLE,Temp,LowestEnergy, BndLines(:)
            end do ! end of temperature cycles

            deallocate(Ngrid)
        end if

        call mpi_bcast(BndLines, N, MPI_integer, 0, comm2D_, ierr)
        deallocate(ProcCoordAll)
        if (myProc_ .eq. 0 .and. DoLoadBalancing) then
            write (6, '(A)') 'Optimum load balanced boundary lines in the property grid:'
            write (6, '(A,16I4)') '  X lines:', BndLines(1:Nxd - 1)
            write (6, '(A,16I4)') '  Y lines:', BndLines(Nxd:N)
        end if
    end if

    ! Make the subdomain position array for this processor from BndLines
    IDX(1) = 1; IDX(Nxd + 1) = prpAll_obj%Nx + 1; if (Nxd .eq. 1) IDX(Nxd + 1) = prpAll_obj%Nx
    IDY(1) = 1; IDY(Nyd + 1) = prpAll_obj%Ny + 1; if (Nyd .eq. 1) IDY(Nyd + 1) = prpAll_obj%Ny
    IDX(2:Nxd) = BndLines(1:Nxd - 1); IDY(2:Nyd) = BndLines(Nxd:N)
    ixyPrp_(1:2, 1) = IDX(ProcCoord(1) + 1:ProcCoord(1) + 2)
    ixyPrp_(1:2, 2) = IDY(ProcCoord(2) + 1:ProcCoord(2) + 2)

end subroutine shd_mpi__optimizeProcDomain

real function shd_mpi__ran(idum)

    ! Returns uniform random deviate between 0.0 and 1.0 (exclusive of endpoints).
    ! Call ran with idum negative to initialize; thereafter, do not alter idum.

    integer, parameter :: K4B = selected_int_kind(9)
    integer(K4B), intent(inout) :: idum
    integer(K4B), parameter :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
    real, save :: am
    integer(K4B), save :: ix = -1, iy = -1, k

    if ((idum .le. 0) .or. (iy .lt. 0)) then
        am = nearest(1.0, -1.0)/IM
        iy = ior(ieor(888889999, abs(idum)), 1)
        ix = ieor(777755555, abs(idum))
        idum = abs(idum) + 1
    end if
    ix = ieor(ix, ishft(ix, 13))
    ix = ieor(ix, ishft(ix, -17))
    ix = ieor(ix, ishft(ix, 5))
    k = iy/IQ
    iy = IA*(iy - k*IQ) - IR*k
    if (iy .lt. 0) iy = iy + IM
    shd_mpi__ran = am*ior(iand(IM, ieor(ix, iy)), 1)

end function shd_mpi__ran

subroutine shd_mpi__integrateRay(BCFLAG, IPFLAG, Nx, Ny, Nz, NA,&
    &EXTINCT, Nstokes, SOURCE, KANG,&
    &GRIDRAD, MU, PHI,&
    &X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

    ! Integrates the radiative transfer equation tracing a ray starting
    ! at X,Y,Z on boundary (face) SIDE to get the radiance in direction MU,PHI.
    ! The source function and extinction are integrated backwards from the
    ! starting point.  As each cell is crossed the integration of the source
    ! function across the cell is done and added to the radiance.  If there
    ! are known radiances at the end of the ray path in GRIDRAD, then the
    ! radiance there is interpolated from the four grid points on the face
    ! and VALIDRAD is returned true. Otherwise VALIDRAD is returned false.
    ! TRANSMIT is returned with the transmission from the ray starting point
    ! to the ending point, and RADIANCE has either the total Stokes radiance or
    ! the accumulated radiance until the boundary (TRANSMIT and RADIANCE must be
    ! properly initialized on input).  The ending location X,Y,Z and the
    ! subdomain SIDE (1=-X,2=+X,3=-Y,4=+Y) are also returned.

    integer, intent(in) :: BCFLAG, IPFLAG, Nx, Ny, Nz, NA
    integer, intent(in) :: Nstokes, KANG
    real, intent(in) :: EXTINCT(rtm_obj_%Npts)
    real, intent(in) :: SOURCE(Nstokes, NA, rtm_obj_%Npts), GRIDRAD(Nstokes, rtm_obj_%Npts)
    real, intent(in) :: MU, PHI
    real, intent(inout) :: X, Y, Z
    integer, intent(inout) :: SIDE
    real, intent(inout) :: TRANSMIT, RADIANCE(Nstokes)
    logical, intent(out) :: VALIDRAD

    integer :: BITX, BITY, BITZ, IOCT, NXC, NYC, BTEST
    integer :: IX, IY, IZ, IL, IM, IU, IPTR, DIR, ICELL, INEXTCELL
    integer :: I1, I2, I3, I4, IOPP, IFACE, JFACE, KFACE, IC
    integer :: GRIDFACE(4, 6), OPPFACE(6)
    logical :: IPINX, IPINY
    real(RD_) :: CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(Nstokes), SRCEXT0(Nstokes), SRCEXT1(Nstokes)
    real(RD_) :: EXT0P, SRCEXT0P(Nstokes)
    real(RD_) :: TAU, TRANSCELL, ABSCELL, TRANS, RAD(Nstokes), RAD0(Nstokes)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8/
    data OPPFACE/2, 1, 4, 3, 6, 5/

    EPS = 1.0e-3*(xyz_Rtm_(3, i8Grid_Cell_(8, 1)) - xyz_Rtm_(3, i8Grid_Cell_(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    CX = sqrt(1.0 - MU**2)*cos(PHI + DPI_)
    CY = sqrt(1.0 - MU**2)*sin(PHI + DPI_)
    CZ = -MU
    if (abs(CX) .gt. 1.0e-5) then
        CXINV = 1.0d0/CX
    else
        CX = 0.0
        CXINV = 1.0e6
    end if
    if (abs(CY) .gt. 1.0e-5) then
        CYINV = 1.0d0/CY
    else
        CY = 0.0
        CYINV = 1.0e6
    end if
    CZINV = 1.0d0/CZ

    ! Setup the sweeping direction and the gridpoint location
    !   BITc is 0 for positive direction trace back ray, 1 for negative.
    if (CX .lt. 0.0) then
        BITX = 1
    else
        BITX = 0
    end if
    if (CY .lt. 0.0) then
        BITY = 1
    else
        BITY = 0
    end if
    if (CZ .lt. -1.0e-3) then
        BITZ = 1
    else if (CZ .gt. 1.0e-3) then
        BITZ = 0
    else
        call shd_mpi__abort('Error [mpi__integrateRay]: Bad MU.')
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate goes to.
    IOCT = 1 + BITX + 2*BITY + 4*BITZ

    ! Find the base grid cell that the starting point is in; first in Z
    IL = 0
    IU = Nz
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. rtm_obj_%z(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)

    ! Then get the base grid cell in X and Y assuming evenly spaced grid
    NXC = Nx
    if (isMpiX_ .and. .not. is1DX_) NXC = Nx - 1
    NYC = Ny
    if (isMpiY_ .and. .not. is1DYx_) NYC = Ny - 1
    if (is1DX_) then
        IX = 1 + nint((Nx - 1)*(X - rtm_obj_%x(1))/(rtm_obj_%x(Nx) - rtm_obj_%x(1)))
    else
        IX = 1 +  int((Nx - 1)*(X - rtm_obj_%x(1))/(rtm_obj_%x(Nx) - rtm_obj_%x(1)))
    end if
    if (is1DYx_) then
        IY = 1 + nint((Ny - 1)*(Y - rtm_obj_%y(1))/(rtm_obj_%y(Ny) - rtm_obj_%y(1)))
    else
        IY = 1 +  int((Ny - 1)*(Y - rtm_obj_%y(1))/(rtm_obj_%y(Ny) - rtm_obj_%y(1)))
    end if
    if ((ix .lt. 1) .or. (iy .lt. 1) .or. (ix .gt. (Nx + 1)) .or. (iy .gt. (Ny + 1))) then
        write (6, *) 'Bad ix and/or iy', x, y, z, ix, iy, iz, rtm_obj_%x(1), rtm_obj_%x(Nx), rtm_obj_%y(1), rtm_obj_%y(ny), myProc_
        write (6, *) 'Bad ix and/or iy', side, transmit
        call shd_mpi__abort('Error [mpi_integrateRay]: Bad ix and/or iy, see above.')
    end if
    IX = min(NXC, IX)
    IY = min(NYC, IY)

    ! Make the base grid cell pointer
    ICELL = IZ + (Nz - 1)*(IY - 1) + (Nz - 1)*NYC*(IX - 1)

    ! Trace down the tree: point to the positive child cell, get the grid
    ! point that is on the negative (X,Y,Z) side, use the flags to find
    ! which dimension (X,Y,Z) the parent cell was split in, then compare
    ! the grid point with the test point (X,Y,Z) to find which child cell
    ! (pos or neg) the test point is in.
    do while (i2Tree_Cell_(2, ICELL) .gt. 0)
        DIR = ibits(int(flag_Cell_(ICELL)), 2, 2)
        IC = i2Tree_Cell_(2, ICELL) + 1
        IPTR = i8Grid_Cell_(1, IC)
        if (DIR .eq. 1) then
            if (X .lt. xyz_Rtm_(1, IPTR)) IC = IC - 1
        else if (DIR .eq. 2) then
            if (Y .lt. xyz_Rtm_(2, IPTR)) IC = IC - 1
        else if (DIR .eq. 3) then
            if (Z .lt. xyz_Rtm_(3, IPTR)) IC = IC - 1
        end if
        ICELL = IC
    end do

    IFACE = SIDE
    ! Get the four gridpoints for this face
    I1 = i8Grid_Cell_(GRIDFACE(1, IFACE), ICELL)
    I2 = i8Grid_Cell_(GRIDFACE(2, IFACE), ICELL)
    I3 = i8Grid_Cell_(GRIDFACE(3, IFACE), ICELL)
    I4 = i8Grid_Cell_(GRIDFACE(4, IFACE), ICELL)
    IPINX = btest(int(flag_Cell_(ICELL)), 0)
    IPINY = btest(int(flag_Cell_(ICELL)), 1)
    ! Compute the face interpolation factors
    if (IFACE .le. 2) then
        U = (Z - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
        if (IPINY) then
            V = 0.5
        else
            V = (Y - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I2) - xyz_Rtm_(2, I1))
        end if
    else
        U = (Z - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
        if (IPINX) then
            V = 0.5
        else
            V = (X - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
        end if
    end if

    ! Bilinearly interpolate the starting extinction and source function
    F1 = (1 - U)*(1 - V); F2 = (1 - U)*V; F3 = U*(1 - V); F4 = U*V
    EXT1 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + F3*EXTINCT(I3) + F4*EXTINCT(I4)
    ! Correctly interpolate source using extinction*source
    SRCEXT1(:) = F1*SOURCE(:, KANG, I1)*EXTINCT(I1) + F2*SOURCE(:, KANG, I2)*EXTINCT(I2)&
        &+ F3*SOURCE(:, KANG, I3)*EXTINCT(I3) + F4*SOURCE(:, KANG, I4)*EXTINCT(I4)
    EXT1 = max(0.0d0, EXT1)
    SRCEXT1(1) = max(0.0d0, SRCEXT1(1))

    ! Loop until finding a face with known radiances or reaching the boundary
    XE = X; YE = Y; ZE = Z
    TRANS = TRANSMIT
    RAD(:) = RADIANCE(:)
    VALIDRAD = .false.
    do while (.not. VALIDRAD .and. ICELL .gt. 0)
        IPINX = btest(int(flag_Cell_(ICELL)), 0)
        IPINY = btest(int(flag_Cell_(ICELL)), 1)

        ! Find boundaries of the current cell
        !   Find the three possible intersection planes (X,Y,Z)
        !   from the coordinates of the opposite corner grid point
        IOPP = i8Grid_Cell_(9 - IOCT, ICELL)
        ! Get the distances to the 3 planes and select the closest
        if (IPINX) then
            SOX = 1.0e20
        else
            SOX = (xyz_Rtm_(1, IOPP) - XE)*CXINV
        end if
        if (IPINY) then
            SOY = 1.0e20
        else
            SOY = (xyz_Rtm_(2, IOPP) - YE)*CYINV
        end if
        SOZ = (xyz_Rtm_(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'SO<0: ', MU, PHI, X, Y, Z, SOX, SOY, SOZ, SO, ICELL, TRANS
            call shd_mpi__abort('Error [mpi__integrateRay]: SO<0, see above.')
        end if
        ! Compute the coordinates of the cell exitting location
        if (.not. IPINX) XE = XE + SO*CX
        if (.not. IPINY) YE = YE + SO*CY
        ZE = ZE + SO*CZ

        ! Get the intersection face number (i.e. ptr6Cell_Cellindex)
        if (SOX .le. SOZ .and. SOX .le. SOY) then
            IFACE = 2 - BITX
            JFACE = 1
        else if (SOY .le. SOZ) then
            IFACE = 4 - BITY
            JFACE = 2
        else
            IFACE = 6 - BITZ
            JFACE = 3
        end if
        ! Get the next cell to go to
        INEXTCELL = i6Cell_Cell_(IFACE, ICELL)
        if (INEXTCELL .lt. 0) then
            call shd_util__findNextCell(XE, YE, ZE, IFACE, JFACE, ICELL,&
                &INEXTCELL)
        end if

        ! Get the face grid pointers
        if (i6Cell_Cell_(IFACE, ICELL) .ge. 0) then
            ! If going to same or larger face then use previous face
            KFACE = IFACE
            IC = ICELL
        else
            ! If going to smaller face then use next face (more accurate)
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
        end if
        I1 = i8Grid_Cell_(GRIDFACE(1, KFACE), IC)
        I2 = i8Grid_Cell_(GRIDFACE(2, KFACE), IC)
        I3 = i8Grid_Cell_(GRIDFACE(3, KFACE), IC)
        I4 = i8Grid_Cell_(GRIDFACE(4, KFACE), IC)
        ! Compute the face interpolation factors
        if (JFACE .eq. 1) then
            U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
            if (IPINY) then
                V = 0.5
            else
                V = (YE - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I2) - xyz_Rtm_(2, I1))
            end if
        else if (JFACE .eq. 2) then
            U = (ZE - xyz_Rtm_(3, I1))/(xyz_Rtm_(3, I3) - xyz_Rtm_(3, I1))
            if (IPINX) then
                V = 0.5
            else
                V = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
            end if
        else
            if (IPINY) then
                U = 0.5
            else
                U = (YE - xyz_Rtm_(2, I1))/(xyz_Rtm_(2, I3) - xyz_Rtm_(2, I1))
            end if
            if (IPINX) then
                V = 0.5
            else
                V = (XE - xyz_Rtm_(1, I1))/(xyz_Rtm_(1, I2) - xyz_Rtm_(1, I1))
            end if
        end if
        if ((u .lt. -1.0e-4) .or. (v .lt. -1.0e-4) .or. (u .gt. 1.0001) .or. (v .gt. 1.0001)) then
            write (6, *) 'u,v<0 and/or u,v>1: ', mu, phi, x, y, z, xe, ye, ze, u, v, ipinx, ipiny
            call shd_mpi__abort('Error [mpi__integrateRay]: u,v<0 and/or u,v>1, see above.')
        end if

        ! Get the location coordinate (does the boundary wrapping)
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XE = xyz_Rtm_(1, i8Grid_Cell_(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YE = xyz_Rtm_(2, i8Grid_Cell_(IOCT, INEXTCELL))
            else
                ZE = xyz_Rtm_(3, i8Grid_Cell_(IOCT, INEXTCELL))
            end if
        end if

        ! Interpolate extinction and source function at face intersection
        F1 = (1 - U)*(1 - V); F2 = (1 - U)*V; F3 = U*(1 - V); F4 = U*V
        EXT0 = F1*EXTINCT(I1) + F2*EXTINCT(I2) + F3*EXTINCT(I3) + F4*EXTINCT(I4)
        ! Correctly interpolate source using extinction*source
        SRCEXT0(:) = F1*SOURCE(:, KANG, I1)*EXTINCT(I1) + F2*SOURCE(:, KANG, I2)*EXTINCT(I2)&
            &+ F3*SOURCE(:, KANG, I3)*EXTINCT(I3) + F4*SOURCE(:, KANG, I4)*EXTINCT(I4)
        EXT0 = max(0.0d0, EXT0)
        SRCEXT0(1) = max(0.0d0, SRCEXT0(1))
        ! Compute the cell radiance: integration of the source function
        EXT = 0.5*(EXT0 + EXT1)
        TAU = EXT*SO
        if (TAU .ge. 0.5) then
            TRANSCELL = exp(-TAU)
            ABSCELL = 1.0 - TRANSCELL
        else
            ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - FRAC13_*TAU*(1 - 0.25*TAU)))
            TRANSCELL = 1.0 - ABSCELL
        end if
        if (TAU .le. 2.0) then
            if (EXT .eq. 0.0) then
                SRC(:) = 0.0
            else
                ! Linear extinction, linear source*extinction, to first order
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:))&
                    &+ FRAC112_*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
            end if
        else
            ! Combined first order expansion and constant extinction formula
            EXT0P = EXT0
            SRCEXT0P(:) = SRCEXT0(:)
            if (TAU .gt. 4.0) then
                EXT0P = EXT1 + (EXT0 - EXT1)*4.0/TAU
                if (EXT0 .gt. 0.0) SRCEXT0P(:) = SRCEXT0(:)*EXT0P/EXT0
            end if
            SRC(:) = 1.0/(EXT0P + EXT1)*(SRCEXT0P(:) + SRCEXT1(:)&
                &+ (EXT0P*SRCEXT1(:) - EXT1*SRCEXT0P(:))*2.0/(EXT0P + EXT1)&
                &*(1 - 2/TAU + 2*TRANSCELL/ABSCELL))
        end if
        SRC(1) = max(SRC(1), 0.0d0)

        ! Add in the cell radiance and update the transmission.
        !  If this is a boundary or the transmission is below the
        !  cutoff and we have a valid radiance then set the flag
        !  to stop the tracing and add in the interpolated face radiance.
        RAD(:) = RAD(:) + TRANS*SRC(:)*ABSCELL
        TRANS = TRANS*TRANSCELL
        if (RAD(1) .lt. -1.0e-5) then
            write (6, '(3(1X,F5.3),8(1X,E12.5))')&
                &xe, ye, ze, ext0, ext1, so, tau, trans, src, abscell, rad(:)
        end if
        if (GRIDRAD(1, I1) .ge. 0.0 .and. GRIDRAD(1, I2) .ge. 0.0&
            .and. GRIDRAD(1, I3) .ge. 0.0 .and. GRIDRAD(1, I4) .ge. 0.0) then
            VALIDRAD = .true.
            RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)&
                &+ F3*GRIDRAD(:, I3) + F4*GRIDRAD(:, I4)
            RAD(:) = RAD(:) + TRANS*RAD0(:)
        !else if (TRANS .lt. 1.0e-5) then
        else if (TRANS .lt. accSolv_) then
            VALIDRAD = .true.
        else
            EXT1 = EXT0
            SRCEXT1(:) = SRCEXT0(:)
            ICELL = INEXTCELL
        end if
    end do
    X = XE; Y = YE; Z = ZE
    RADIANCE(:) = RAD(:)
    TRANSMIT = TRANS
    SIDE = IFACE

end subroutine shd_mpi__integrateRay

real function loadObject(BndLines, Nxd, Nyd, ProcCoordAll, wtproc,&
    &Ngrid, prpAll_obj)

    ! Calculates the objective function to minimize in doing the processor
    ! load balancing optimization.  The boundary line indices in X and Y to
    ! the property grid, which is the optimization control vector is input
    ! in BndLines.  The objective function is the normalized sum of the
    ! squares of the difference between the number of grid points (Ngrid)
    ! of each domain and the mean number, weighted by the relative processor
    ! times.

    type(shd_type_prp), intent(in) :: prpAll_obj

    integer, intent(in) :: Nxd, Nyd, BndLines(Nxd + Nyd - 2)
    integer, intent(in) :: ProcCoordAll(2, 0:Nxd*Nyd - 1)
    real, intent(in) :: wtproc(0:Nxd*Nyd)
    integer, intent(in) :: Ngrid(prpAll_obj%Nx + 1, prpAll_obj%Ny + 1)
    integer :: N, I, J, IX, IY, IDX(Nxd + 1), IDY(Nyd + 1)
    real    :: OBJA
    real, save :: sumN = -1.0, a

    if (sumN .lt. 0) then
        sumN = sum(Ngrid(:, :)); a = float(Nxd*Nyd)/sumN
    end if
    N = Nxd + Nyd - 2
    IDX(1) = 1; IDX(Nxd + 1) = prpAll_obj%Nx + 1; if (Nxd .eq. 1) IDX(Nxd + 1) = prpAll_obj%Nx
    IDY(1) = 1; IDY(Nyd + 1) = prpAll_obj%Ny + 1; if (Nyd .eq. 1) IDY(Nyd + 1) = prpAll_obj%Ny

    OBJA = 0.0
    do J = 0, Nxd*Nyd - 1
        IDX(2:Nxd) = BndLines(1:Nxd - 1); IDY(2:Nyd) = BndLines(Nxd:N)
        IX = ProcCoordAll(1, J) + 1; IY = ProcCoordAll(2, J) + 1
        OBJA = OBJA + (1.0/Nxd*Nyd)*&
            &(a*wtproc(J)*sum(Ngrid(IDX(IX):IDX(IX + 1), IDY(IY):IDY(IY + 1))) - 1.0)**2
    end do
    loadObject = OBJA

end function loadObject
#endif

end module shd_mpi
