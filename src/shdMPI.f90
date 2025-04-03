! All of the MPI routines for the parallel version of SHDOM are
! in this file.  If you do not have MPI, compile with shdom_nompi.f90
! instead.

module shd_mpi

    use globals
    use mpi_f08

    use shd_init
    use shd_util
    use shd_grid
    use shd_src
    use shd_comp

    implicit none
    private

    integer :: NUMPROC_, MYPROC_, IPROCNEIGH_(4), IXYPRP_(2, 2), IXYOUT_(2, 2)
    logical, parameter :: LOADBALANCE_ = .false.
    real(RD_) :: STARTTIME_
    integer, allocatable :: BXYPTR_(:, :, :), IBXYPTR_(:, :, :)
    real, allocatable :: BXYPOS_(:, :, :)
    type(MPI_Comm) :: COMM2D_

    public :: shd_mpi__start
    public :: shd_mpi__map
    public :: shd_mpi__bcastInput
    public :: shd_mpi__bcastPropSize
    public :: shd_mpi__scatProp
    public :: shd_mpi__readBcastMemParam
    public :: shd_mpi__bcastSurfSize
    public :: shd_mpi__bcastSurfParam
    public :: shd_mpi__bcastKdistSize
    public :: shd_mpi__bcastKdistParam
    public :: shd_mpi__gatherOutput
    public :: shd_mpi__end
    public :: shd_mpi__abort

    public :: shd_mpi__makeDirect
    public :: shd_mpi__computeRad
    public :: shd_mpi__visualizeRad
    public :: shd_mpi__findAlbMax
    public :: shd_mpi__findTotSplitCritMax
    public :: shd_mpi__findBoundaryPnt
    public :: shd_mpi__calBoundaryRad
    public :: shd_mpi__calAccelSolCrit
    public :: shd_mpi__unifySplit
    public :: shd_mpi__sumCpuTime

contains

subroutine shd_mpi__start(MASTERPROC)

    ! Initializes the MPI system, and gets the number of processors in use
    ! and the current processor.  If this is processor 0 then MASTERPROC
    ! is returned true.  Also starts the logging files "shdom???.log"
    ! for the standard output from the non-master processors.

    logical, intent(out) :: MASTERPROC
    integer :: ierr

    call MPI_INIT(ierr)
    if (ierr .ne. MPI_SUCCESS) then
        write (6, *) 'Error [START_MPI]: Terminating ...'
        stop
    end if

    call MPI_COMM_SIZE(MPI_COMM_WORLD, NUMPROC_, ierr)
    call MPI_COMM_RANK(MPI_COMM_WORLD, MYPROC_, ierr)

    call MPI_Errhandler_set(MPI_COMM_WORLD, MPI_ERRORS_RETURN, ierr)
    !call MPI_Errhandler_set (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL, ierr)

    STARTTIME_ = MPI_WTIME()

    if (MYPROC_ .eq. 0) then
        MASTERPROC = .true.
        write (6, '(A)') 'Message [START_MPI]: Starting MPI version of SHDOM ...'
    else
        MASTERPROC = .false.
    end if
end subroutine shd_mpi__start

subroutine shd_mpi__map(BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY, PROPFILE,&
    &XSTART, YSTART, RUNNAME)

    ! Sets up the mapping of processors for the SHDOM domain (creating the
    ! 2D communicator, COMM2D_).  Determines the subdomain property grid
    ! (NPX,NPY) and internal base grid (NX,NY) sizes, and outputs the
    ! starting position of this subdomain (XSTART,YSTART).  Bits 2 and 3
    ! of BCFLAG are set to indicate multiple processors in the X and Y
    ! directions, respectively.  Open boundary conditions (bits 0 and 1 in
    ! BCFLAG) are overridden by the multiple processor flags (i.e. multiple
    ! processor requires a periodic full domain).
    !   The neighboring processors are stored in iprocneigh_(:), with the
    ! following indexing scheme: 1=-X, 2=+X, 3-=-Y, 4=+Y.   The IXYPRP_ array
    ! (module global variable) is set with the starting and ending property
    ! grid indices of this subdomain (IXYPRP_(:,1) is for X, and IXYPRP_(:,2)
    ! is for Y).  The IXYOUT_ array is set with the starting index and number
    ! of grid points of the internal base grid for this subdomain (for the
    ! output arrays).

    integer, intent(inout) :: BCFLAG, NPX, NPY, NX, NY, NZ
    character(len=*), intent(in) :: PROPFILE, RUNNAME
    real, intent(inout) :: DELX, DELY
    real, intent(out) :: XSTART, YSTART
    integer :: NPXT, NPYT, NXT, NYT
    integer :: i, ierr, dims(2), coords(2)
    logical :: periodic(2) = (/.true., .true./), reorder = .false.
    real    :: XEND, YEND
    character(len=4) :: procstr
    character(len=FNAMEMAXLEN_) :: logfile

    NPXT = NPX; NPYT = NPY; NXT = NX; NYT = NY

    !if ((NUMPROC_ .eq. 1) .or. ((NPXT .eq. 1) .and. (NPYT .eq. 1))) then
    if (NUMPROC_ .eq. 1) then
        XSTART = 0.0; YSTART = 0.0
        return
    end if

    if (MYPROC_ > 0) then
        write (procstr, '(I4.4)') MYPROC_
        logfile = trim(adjustl(RUNNAME))//'.sHdOm-log.iproc-'//procstr
        open (UNIT=6, FILE=trim(adjustl(logfile)), STATUS='UNKNOWN')
    end if

    ! Partition 2D Cartesian topology:
    dims(:) = 0
    call MPI_DIMS_CREATE(NUMPROC_, 2, dims, ierr)

    ! Create 2D Cartesian mapping of processors:
    call MPI_CART_CREATE(MPI_COMM_WORLD, 2, dims, periodic, reorder,&
        &COMM2D_, ierr)

    if (dims(1) > 1) BCFLAG = ibset(BCFLAG, 2)
    if (dims(2) > 1) BCFLAG = ibset(BCFLAG, 3)
    ! Multi-processor BC takes precedence over open BCs
    if (btest(BCFLAG, 2) .or. btest(BCFLAG, 3)) BCFLAG = iand(BCFLAG, 12)

    ! Find where in the 2D grid of processors this is is at:
    call MPI_CART_COORDS(COMM2D_, MYPROC_, 2, coords, ierr)

    ! Make the optimal (or default) property grid divisions (in IXYPRP_)
    call shd_mpi__optimizeProcDomain(dims, coords, NPXT, NPYT, PROPFILE, RUNNAME)

    ! Make the new property and base grid sizes and the IXYOUT_ arrays
    if (dims(1) > 1) then
        if (NPXT .gt. 1) then

            NPX = IXYPRP_(2, 1) - IXYPRP_(1, 1) + 1
            NX = nint(NXT*FLOAT(NPX - 1)/NPXT) + 1
            if (abs(NPXT/FLOAT(NXT) - (NPX - 1)/FLOAT(NX - 1)) > 0.0001) then
                call shd_mpi__abort('MAP_SHDOM_MPI: Even X base grid cannot be preserved&
                    &for this NX and number of processors')
            end if
            IXYOUT_(1, 1) = nint(NXT*FLOAT(IXYPRP_(1, 1) - 1)/NPXT) + 1
            IXYOUT_(2, 1) = nint(NXT*FLOAT(NPX - 1)/NPXT)
            XSTART = (IXYPRP_(1, 1) - 1)*DELX
            XEND = (IXYPRP_(2, 1) - 1)*DELX

        else

            NPX = 1

            IXYOUT_(1, 1) = nint(real(NXT)/real(dims(1)))*coords(1) + 1
            IXYOUT_(2, 1) = min(nint(real(NXT)/real(dims(1)))*(coords(1)+1), NXT)

            NX = IXYOUT_(2, 1) - IXYOUT_(1, 1) + 2

            XSTART = real(DELX)/real(NXT) * real((IXYOUT_(1, 1)-1))
            XEND   = real(DELX)/real(NXT) * real((IXYOUT_(2, 1)))

            DELX = XEND - XSTART


        end if
    else
        IXYOUT_(1, 1) = 1
        IXYOUT_(2, 1) = NXT
        XSTART = 0.0
        XEND = NPXT*DELX
    end if

    if (dims(2) > 1) then

        if (NPYT .gt. 1) then

            NPY = IXYPRP_(2, 2) - IXYPRP_(1, 2) + 1
            NY = nint(NYT*FLOAT(NPY - 1)/NPYT) + 1
            if (abs(NPYT/FLOAT(NYT) - (NPY - 1)/FLOAT(NY - 1)) > 0.0001) then
                call shd_mpi__abort('MAP_SHDOM_MPI: Even Y base grid cannot be preserved for this NY and number of processors')
            end if
            IXYOUT_(1, 2) = nint(NYT*FLOAT(IXYPRP_(1, 2) - 1)/NPYT) + 1
            IXYOUT_(2, 2) = nint(NYT*FLOAT(NPY - 1)/NPYT)
            YSTART = (IXYPRP_(1, 2) - 1)*DELY
            YEND = (IXYPRP_(2, 2) - 1)*DELY

        else

            NPY = 1

            IXYOUT_(1, 2) = nint(real(NYT)/real(dims(2)))*coords(2) + 1
            IXYOUT_(2, 2) = min(nint(real(NYT)/real(dims(2)))*(coords(2)+1), NYT)

            NY = IXYOUT_(2, 2) - IXYOUT_(1, 2) + 2

            YSTART = real(DELY)/real(NYT) * real((IXYOUT_(1, 2)-1))
            YEND   = real(DELY)/real(NYT) * real((IXYOUT_(2, 2)))

            DELY = YEND - YSTART

        end if
    else
        IXYOUT_(1, 2) = 1
        IXYOUT_(2, 2) = NYT
        YSTART = 0.0
        YEND = NPYT*DELY
    end if

    ! Get the processors (rank numbers) of the neighbors
    call MPI_CART_SHIFT(COMM2D_, 0, 1, IPROCNEIGH_(1), IPROCNEIGH_(2), ierr)
    call MPI_CART_SHIFT(COMM2D_, 1, 1, IPROCNEIGH_(3), IPROCNEIGH_(4), ierr)

    write (6, '(A)') 'Message [_mpiMap]:'
    write (6, '(2X,A,2(1X,I2))') 'Number of processors in X and Y:', dims(1:2)
    write (6, '(2X,A,I3)') 'This processor: ', MYPROC_
    write (6, '(2X,A,4(1X,I3))') 'Neighbors:', IPROCNEIGH_(1:4)
    write (6, '(2X,2(A,2F7.3))') 'Processor X range: ', XSTART, XEND, '   Y range: ', YSTART, YEND

end subroutine shd_mpi__map

subroutine shd_mpi__optimizeProcDomain(ProcDims, ProcCoord, NPXT, NPYT,&
    &PROPFILE, RUNNAME)

    ! Calculates the location of this processor (subdomain) in the property
    ! grid.  The starting and ending property grid point indices for X and Y
    ! are returned in IXYPRP_.  Each processor sends its location in the
    ! 2D array of subdomains to the master process and the master process
    ! reads the RUNNAME//"_load_balance.inp" file (if there is one) and does
    ! the optimization.  If the "shdom_proc_times.inp" file of relative
    ! processor times exists then the optimization uses this information on
    ! heterogeneous processors to try to give more work to faster processors.

    integer, intent(in) :: ProcDims(2), ProcCoord(2), NPXT, NPYT
    character(len=*), intent(in) :: PROPFILE, RUNNAME
    integer :: NPXTI, NPYTI
    integer :: IX, IY, JX, JY, Nxd, Nyd, N, I, ierr
    integer :: BndLines(ProcDims(1) + ProcDims(2) - 2)
    integer :: IDX(ProcDims(1) + 1), IDY(ProcDims(2) + 1)
    logical :: DoLoadBalancing
    integer, parameter :: NCYCLE = 30 ! number of annealing cooling cycles
    integer, parameter :: NiterFac = 100 ! number of iterations per length of optimization vector (BndLines)
    real, parameter :: Temp0 = 0.1 ! initial temperature (units of change in objective function "energy")
    real, parameter :: Tempfac = 4.0/NCYCLE ! exponential cooling rate
    integer :: cycle, ITER, NITER
    integer :: TrialLines(ProcDims(1) + ProcDims(2) - 2), CurLines(ProcDims(1) + ProcDims(2) - 2)
    integer :: iseed = -3
    real    :: Temp, DelLines, u, shd_mpi__ran
    real    :: loadObject, Energy, NewEnergy, LowestEnergy, DeltaE, Boltz
    integer, allocatable :: ProcCoordAll(:, :), Ngrid(:, :)
    real, allocatable :: wtproc(:)
    character(len=FNAMEMAXLEN_) :: PROPFILEI

    ! Initialize the subdomain boundary lines to evenly spaced
    Nxd = ProcDims(1); Nyd = ProcDims(2); N = Nxd + Nyd - 2
    BndLines(1:Nxd - 1) = (/(nint(I*FLOAT(NPXT)/Nxd) + 1, I=1, Nxd - 1)/)
    BndLines(Nxd:N) = (/(nint(I*FLOAT(NPYT)/Nyd) + 1, I=1, Nyd - 1)/)

    if (LOADBALANCE_) then
        DoLoadBalancing = .false.
        if (MYPROC_ .eq. 0) then
            ! See if we can read the load balancing file and header matches
            open (UNIT=10, FILE=trim(RUNNAME)//'_load_balance.inp', STATUS='OLD', IOSTAT=ierr)
            if (ierr .eq. 0) then
                read (10, *)
                read (10, '(A)') PROPFILEI
                read (10, *) NPXTI, NPYTI
                if (NPXTI .ne. NPXT .or. NPYTI .ne. NPYT .or. PROPFILEI .ne. PROPFILE) then
                    write (6, '(A,A,A)') 'Load balancing file ', trim(RUNNAME),&
                        &'_load_balance.inp is not compatible with property file:'
                    write (6, *) PROPFILE
                    write (6, *)
                else
                    DoLoadBalancing = .true.
                end if
                close (10)
            end if
        end if

        ! Have the master process collect the locations of all
        allocate (ProcCoordAll(2, 0:NUMPROC_ - 1))
        call MPI_GATHER(ProcCoord, 2, MPI_INTEGER, ProcCoordAll, 2, MPI_INTEGER, 0, COMM2D_, ierr)

        if (MYPROC_ .eq. 0 .and. DoLoadBalancing) then
            ! Read the load balancing file to get array of number of grid points
            allocate (Ngrid(NPXT + 1, NPYT + 1))
            open (UNIT=10, FILE=trim(RUNNAME)//'_load_balance.inp', STATUS='OLD')
            read (10, *)
            read (10, *)
            read (10, *)
            do IX = 1, NPXT
                do IY = 1, NPYT
                    read (10, *) JX, JY, Ngrid(IX, IY)
                end do
            end do
            close (10)
            Ngrid(NPXT + 1, :) = Ngrid(1, :)
            Ngrid(:, NPYT + 1) = Ngrid(:, 1)
            Ngrid(NPXT + 1, NPYT + 1) = Ngrid(1, 1)

            ! Read the relative processor times if the "shdom_proc_times.inp" file exists
            allocate (wtproc(NUMPROC_))
            open (UNIT=11, FILE='shdom_proc_times.inp', STATUS='OLD', IOSTAT=ierr)
            if (ierr .eq. 0) then
                do i = 1, NUMPROC_
                    read (11, *) wtproc(i)
                end do
                close (11)
                wtproc(:) = wtproc(:)*NUMPROC_/sum(wtproc(:))
            else
                wtproc(:) = 1.0
            end if

            ! Do the simulated annealing to minimize the objective function
            u = shd_mpi__ran(iseed)
            NITER = NiterFac*N
            Energy = loadObject(BndLines, Nxd, Nyd, ProcCoordAll, wtproc, NPXT, NPYT, Ngrid)
            LowestEnergy = Energy
            CurLines(:) = BndLines(:)
            ! print *, 'Annealing cycles: Cycle Temp Energy  BndLines'
            ! print '(I3,2(1X,F9.6),10(1X,I3))', 0,Temp0,Energy, BndLines(:)

            ! Loop through the cooling cycles
            do cycle = 1, NCYCLE
                Temp = Temp0*exp(-Tempfac*(cycle - 1)) ! lower the temp
                ! Adjust the maximum distance a boundary line is perturbed with the cycle number
                u = FLOAT(cycle - 1)/(NCYCLE - 1)
                DelLines = (1 - u)*0.5*max(NPXT/Nxd, NPYT/Nyd) + u*1.0

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
                    NewEnergy = loadObject(TrialLines, Nxd, Nyd, ProcCoordAll, wtproc, NPXT, NPYT, Ngrid)
                    DeltaE = NewEnergy - Energy
                    if (DeltaE < 0.0) then ! If lower energy accept change
                        CurLines(:) = TrialLines(:)
                        Energy = NewEnergy
                    else !  otherwise sometimes accept it anyway
                        Boltz = exp(-DeltaE/Temp)
                        if (shd_mpi__ran(iseed) < Boltz) then
                            CurLines(:) = TrialLines(:)
                            Energy = NewEnergy
                        end if
                    end if
                    if (Energy < LowestEnergy) then
                        BndLines(:) = CurLines(:)
                        LowestEnergy = Energy
                    end if
                end do
                ! print '(I3,2(1X,F9.6),10(1X,I3))', CYCLE,Temp,LowestEnergy, BndLines(:)
            end do ! end of temperature cycles

            deallocate (Ngrid)
        end if

        call MPI_BCAST(BndLines, N, MPI_INTEGER, 0, COMM2D_, ierr)
        deallocate (ProcCoordAll)
        if (MYPROC_ .eq. 0 .and. DoLoadBalancing) then
            write (6, '(A)') 'Optimum load balanced boundary lines in the property grid:'
            write (6, '(A,16I4)') '  X lines:', BndLines(1:Nxd - 1)
            write (6, '(A,16I4)') '  Y lines:', BndLines(Nxd:N)
        end if
    end if

    ! Make the subdomain position array for this processor from BndLines
    IDX(1) = 1; IDX(Nxd + 1) = NPXT + 1; if (Nxd .eq. 1) IDX(Nxd + 1) = NPXT
    IDY(1) = 1; IDY(Nyd + 1) = NPYT + 1; if (Nyd .eq. 1) IDY(Nyd + 1) = NPYT
    IDX(2:Nxd) = BndLines(1:Nxd - 1); IDY(2:Nyd) = BndLines(Nxd:N)
    IXYPRP_(1:2, 1) = IDX(ProcCoord(1) + 1:ProcCoord(1) + 2)
    IXYPRP_(1:2, 2) = IDY(ProcCoord(2) + 1:ProcCoord(2) + 2)

end subroutine shd_mpi__optimizeProcDomain

real function loadObject(BndLines, Nxd, Nyd, ProcCoordAll, wtproc,&
    &NPXT, NPYT, Ngrid)

    ! Calculates the objective function to minimize in doing the processor
    ! load balancing optimization.  The boundary line indices in X and Y to
    ! the property grid, which is the optimization control vector is input
    ! in BndLines.  The objective function is the normalized sum of the
    ! squares of the difference between the number of grid points (Ngrid)
    ! of each domain and the mean number, weighted by the relative processor
    ! times.

    integer, intent(IN) :: Nxd, Nyd, BndLines(Nxd + Nyd - 2)
    integer, intent(IN) :: ProcCoordAll(2, 0:Nxd*Nyd - 1)
    real, intent(IN) :: wtproc(0:Nxd*Nyd)
    integer, intent(IN) :: NPXT, NPYT, Ngrid(NPXT + 1, NPYT + 1)
    integer :: N, I, J, IX, IY, IDX(Nxd + 1), IDY(Nyd + 1)
    real    :: OBJA
    real, save :: sumN = -1.0, a

    if (sumN < 0) then
        sumN = sum(Ngrid(:, :)); a = FLOAT(Nxd*Nyd)/sumN
    end if
    N = Nxd + Nyd - 2
    IDX(1) = 1; IDX(Nxd + 1) = NPXT + 1; if (Nxd .eq. 1) IDX(Nxd + 1) = NPXT
    IDY(1) = 1; IDY(Nyd + 1) = NPYT + 1; if (Nyd .eq. 1) IDY(Nyd + 1) = NPYT

    OBJA = 0.0
    do J = 0, Nxd*Nyd - 1
        IDX(2:Nxd) = BndLines(1:Nxd - 1); IDY(2:Nyd) = BndLines(Nxd:N)
        IX = ProcCoordAll(1, J) + 1; IY = ProcCoordAll(2, J) + 1
        OBJA = OBJA + (1.0/Nxd*Nyd)*&
            &(a*wtproc(J)*sum(Ngrid(IDX(IX):IDX(IX + 1), IDY(IY):IDY(IY + 1))) - 1.0)**2
    end do
    loadObject = OBJA

end function loadObject

real function shd_mpi__ran(idum)

    ! Returns uniform random deviate between 0.0 and 1.0 (exclusive of endpoints).
    ! Call ran with idum negative to initialize; thereafter, do not alter idum.

    integer, parameter :: K4B = selected_int_kind(9)
    integer(K4B), intent(INOUT) :: idum
    integer(K4B), parameter :: IA = 16807, IM = 2147483647, IQ = 127773, IR = 2836
    real, save :: am
    integer(K4B), save :: ix = -1, iy = -1, k

    if (idum <= 0 .or. iy < 0) then
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
    if (iy < 0) iy = iy + IM
    shd_mpi__ran = am*ior(iand(IM, ieor(ix, iy)), 1)

end function shd_mpi__ran

subroutine shd_mpi__bcastInput(&
    &RUNNAME, PROPFILE, SFCFILE, CKDFILE,&
    &INSAVEFILE, OUTSAVEFILE,&
    &NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
    &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
    &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
    &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS,&
    &MAX_TOTAL_MB, SPLITTING_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

    ! Broadcasts all the user input parameters (except unneeded file names)
    ! from the master processor to the slave processors.

    integer, intent(INOUT) :: NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    integer, intent(IN)    :: MAXOUT, MAXPAR
    integer, intent(INOUT) :: MAXITER, NUMOUT
    logical, intent(INOUT) :: KDIST, DELTAM, ACCELFLAG
    real, intent(INOUT) :: SOLARFLUX, SOLARMU, SOLARAZ
    real, intent(INOUT) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real, intent(INOUT) :: SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR, MAXOUT)
    real, intent(INOUT) :: MAX_TOTAL_MB, SPLITTING_FACTOR
    real, intent(INOUT) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(LEN=1), intent(INOUT) :: SRCTYPE, GRIDTYPE, UNITS, OUTTYPES(*)
    character(LEN=FNAMEMAXLEN_), intent(INOUT) :: RUNNAME, PROPFILE, SFCFILE, CKDFILE
    character(LEN=FNAMEMAXLEN_), intent(INOUT) :: INSAVEFILE, OUTSAVEFILE

    integer :: ierr, intbuf(10)
    logical :: logbuf(3)
    real    :: realbuf(16)
    character(LEN=1) :: char1buf(3)
    character(LEN=FNAMEMAXLEN_) :: char256buf(6)
    character(LEN=4) :: procstr

    if (NUMPROC_ <= 1) return

    intbuf(:) = (/NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG, MAXITER, NUMOUT/)
    call MPI_BCAST(intbuf, 10, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    NSTOKES = intbuf(1); NX = intbuf(2); NY = intbuf(3); NZ = intbuf(4)
    NMU = intbuf(5); NPHI = intbuf(6); BCFLAG = intbuf(7); IPFLAG = intbuf(8)
    MAXITER = intbuf(9); NUMOUT = intbuf(10)

    logbuf(:) = (/KDIST, DELTAM, ACCELFLAG/)
    call MPI_BCAST(logbuf, 3, MPI_LOGICAL, 0, MPI_COMM_WORLD, ierr)
    KDIST = logbuf(1); DELTAM = logbuf(2); ACCELFLAG = logbuf(3)

    realbuf(:) = (/SOLARFLUX, SOLARMU, SOLARAZ, GNDTEMP, GNDALBEDO, SKYRAD,&
        &WAVENO(1:2), WAVELEN, SOLACC, SPLITACC, SHACC,&
        &MAX_TOTAL_MB, SPLITTING_FACTOR, NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO/)
    call MPI_BCAST(realbuf, 16, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    SOLARFLUX = realbuf(1); SOLARMU = realbuf(2); SOLARAZ = realbuf(3)
    GNDTEMP = realbuf(4); GNDALBEDO = realbuf(5); SKYRAD = realbuf(6)
    WAVENO(1:2) = realbuf(7:8); WAVELEN = realbuf(9)
    SOLACC = realbuf(10); SPLITACC = realbuf(11); SHACC = realbuf(12)
    MAX_TOTAL_MB = realbuf(13); SPLITTING_FACTOR = realbuf(14)
    NUM_SH_TERM_FACTOR = realbuf(15); CELL_TO_POINT_RATIO = realbuf(16)

    char1buf(:) = (/SRCTYPE, GRIDTYPE, UNITS/)
    call MPI_BCAST(char1buf, 3, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    SRCTYPE = char1buf(1); GRIDTYPE = char1buf(2); UNITS = char1buf(3)

    char256buf(:) = (/PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, RUNNAME/)
    call MPI_BCAST(char256buf, 6*256, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)
    PROPFILE = char256buf(1); SFCFILE = char256buf(2); CKDFILE = char256buf(3)
    INSAVEFILE = char256buf(4); OUTSAVEFILE = char256buf(5); RUNNAME = char256buf(6)

    call MPI_BCAST(OUTPARMS, MAXPAR*NUMOUT, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    call MPI_BCAST(OUTTYPES, NUMOUT*1, MPI_CHARACTER, 0, MPI_COMM_WORLD, ierr)

    write (procstr, '(I4.4)') MYPROC_
    if (INSAVEFILE(1:4) .ne. 'NONE') INSAVEFILE = trim(INSAVEFILE)//'.iproc-'//procstr
    if (OUTSAVEFILE(1:4) .ne. 'NONE') OUTSAVEFILE = trim(OUTSAVEFILE)//'.iproc-'//procstr
end subroutine shd_mpi__bcastInput

subroutine shd_mpi__bcastPropSize(NPX, NPY, NPZ, DELX, DELY,&
    &NUMPHASE, MAXLEG, MAXPGL)

    ! Broadcasts the property file size information from the master processor
    ! to the slave processors.

    integer, intent(INOUT) :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
    real, intent(INOUT) :: DELX, DELY

    integer :: ierr, intbuf(6)
    real    :: realbuf(2)

    if (NUMPROC_ <= 1) return

    intbuf(:) = (/NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL/)
    call MPI_BCAST(intbuf, 6, MPI_INTEGER, 0, MPI_COMM_WORLD, ierr)
    NPX = intbuf(1); NPY = intbuf(2); NPZ = intbuf(3)
    NUMPHASE = intbuf(4); MAXLEG = intbuf(5); MAXPGL = intbuf(6)

    realbuf(:) = (/DELX, DELY/)
    call MPI_BCAST(realbuf, 2, MPI_REAL, 0, MPI_COMM_WORLD, ierr)
    DELX = realbuf(1); DELY = realbuf(2)
end subroutine shd_mpi__bcastPropSize

subroutine shd_mpi__scatProp(NPXT, NPYT, NPX, NPY, NPZ,&
    &NSTLEG, NLEG, NUMPHASE,&
    &ZLEVELS, MAXASYM,&
    &TEMPPT, EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,&
    &TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)

    ! Sends portions of the optical properties on the property grid to each
    ! of the processors.  Deals with the duplicating the grid points on the
    ! boundaries between processors.

    integer, intent(IN) :: NPXT, NPYT, NPX, NPY, NPZ, NSTLEG, NUMPHASE
    integer, intent(INOUT) :: NLEG
    real, intent(INOUT) :: ZLEVELS(NPZ), MAXASYM
    real, intent(IN) :: TEMPPT(NPZ, NPYT, NPXT), EXTINCTPT(NPZ, NPYT, NPXT)
    real, intent(IN) :: ALBEDOPT(NPZ, NPYT, NPXT), LEGENPT(*)
    integer*2, intent(IN) :: IPHASEPT(NPZ, NPYT, NPXT)
    real, intent(OUT) :: TEMPP(NPZ, NPY, NPX), EXTINCTP(NPZ, NPY, NPX)
    real, intent(OUT) :: ALBEDOP(NPZ, NPY, NPX), LEGENP(*)
    integer*2, intent(OUT) :: IPHASEP(NPZ, NPY, NPX)

    integer :: SX, EX, NX, SY, EY, NY, N, MAXN, MaxNX, MaxNY
    integer :: iproc, irq, ierr
    integer, allocatable :: IXYALL(:, :, :)
    real, allocatable :: realbuf(:), rbuf(:, :, :)
    real, allocatable :: tempbuf(:, :), extbuf(:, :), albbuf(:, :)
    integer*2, allocatable :: ibuf(:, :, :), iphbuf(:, :)

    type(MPI_Request), allocatable :: requests(:)
    type(MPI_Status), allocatable :: status(:)

    ! Broadcast the small stuff that is going to all processors
    allocate (realbuf(NPZ + 1))
    realbuf(1:NPZ + 1) = (/ZLEVELS(1:NPZ), MAXASYM/)
    call MPI_BCAST(realbuf, NPZ + 1, MPI_REAL, 0, COMM2D_, ierr)
    ZLEVELS(1:NPZ) = realbuf(1:NPZ); MAXASYM = realbuf(NPZ + 1)

    ! If we have a tabulated phase function, then broadcast to all processors
    if (NUMPHASE > 0) then
        if (MYPROC_ .eq. 0) LEGENP(1:NUMPHASE*(NLEG + 1)*NSTLEG) =&
            &LEGENPT(1:NUMPHASE*(NLEG + 1)*NSTLEG)
        call MPI_BCAST(NLEG, 1, MPI_INTEGER, 0, COMM2D_, ierr)
        call MPI_BCAST(LEGENP, NSTLEG*(NLEG + 1)*NUMPHASE, MPI_REAL, 0, COMM2D_, ierr)
    end if

    ! Get the IXYPRP_ arrays from all processors to tell where they are in property grid
    allocate (IXYALL(2, 2, 0:NUMPROC_ - 1))
    call MPI_GATHER(IXYPRP_, 4, MPI_INTEGER, IXYALL, 4, MPI_INTEGER, 0, COMM2D_, ierr)

    allocate (requests(4*NUMPROC_ + 4), status(4*NUMPROC_ + 4))

    irq = 0
    if (MYPROC_ .eq. 0) then
        ! Allocate the send buffers for each processor
        MAXN = 1; MaxNX = 1; MaxNY = 1
        do iproc = 0, NUMPROC_ - 1
            NX = IXYALL(2, 1, iproc) - IXYALL(1, 1, iproc) + 1; MaxNX = max(MaxNX, NX)
            NY = IXYALL(2, 2, iproc) - IXYALL(1, 2, iproc) + 1; MaxNY = max(MaxNY, NY)
            MAXN = max(MAXN, NPZ*NY*NX)
        end do
        allocate (rbuf(NPZ, MaxNY, MaxNX), ibuf(NPZ, MaxNY, MaxNX))
        allocate (tempbuf(MAXN, 0:NUMPROC_ - 1), extbuf(MAXN, 0:NUMPROC_ - 1))
        allocate (albbuf(MAXN, 0:NUMPROC_ - 1), iphbuf(MAXN, 0:NUMPROC_ - 1))
        ! The master processor reformats the four property arrays (to deal with
        ! the wrap around at the full domain edge) and sends the arrays with MPI.
        do iproc = 0, NUMPROC_ - 1
            SX = IXYALL(1, 1, iproc); EX = IXYALL(2, 1, iproc)
            SY = IXYALL(1, 2, iproc); EY = IXYALL(2, 2, iproc)
            NX = EX - SX + 1; NY = EY - SY + 1; N = NX*NY*NPZ
            rbuf(:, 1:NY - 1, 1:NX - 1) = TEMPPT(:, SY:EY - 1, SX:EX - 1)
            rbuf(:, 1:NY - 1, NX) = TEMPPT(:, SY:EY - 1, mod(EX - 1, NPXT) + 1)
            rbuf(:, NY, 1:NX - 1) = TEMPPT(:, mod(EY - 1, NPYT) + 1, SX:EX - 1)
            rbuf(:, NY, NX) = TEMPPT(:, mod(EY - 1, NPYT) + 1, mod(EX - 1, NPXT) + 1)
            tempbuf(1:N, iproc) = reshape(rbuf(1:NPZ, 1:NY, 1:NX), (/N/))
            call MPI_ISEND(tempbuf(1:N, iproc), N, MPI_REAL, iproc, 1, COMM2D_, requests(irq + 1), ierr)

            rbuf(:, 1:NY - 1, 1:NX - 1) = EXTINCTPT(:, SY:EY - 1, SX:EX - 1)
            rbuf(:, 1:NY - 1, NX) = EXTINCTPT(:, SY:EY - 1, mod(EX - 1, NPXT) + 1)
            rbuf(:, NY, 1:NX - 1) = EXTINCTPT(:, mod(EY - 1, NPYT) + 1, SX:EX - 1)
            rbuf(:, NY, NX) = EXTINCTPT(:, mod(EY - 1, NPYT) + 1, mod(EX - 1, NPXT) + 1)
            extbuf(1:N, iproc) = reshape(rbuf(1:NPZ, 1:NY, 1:NX), (/N/))
            call MPI_ISEND(extbuf(1:N, iproc), N, MPI_REAL, iproc, 2, COMM2D_, requests(irq + 2), ierr)

            rbuf(:, 1:NY - 1, 1:NX - 1) = ALBEDOPT(:, SY:EY - 1, SX:EX - 1)
            rbuf(:, 1:NY - 1, NX) = ALBEDOPT(:, SY:EY - 1, mod(EX - 1, NPXT) + 1)
            rbuf(:, NY, 1:NX - 1) = ALBEDOPT(:, mod(EY - 1, NPYT) + 1, SX:EX - 1)
            rbuf(:, NY, NX) = ALBEDOPT(:, mod(EY - 1, NPYT) + 1, mod(EX - 1, NPXT) + 1)
            albbuf(1:N, iproc) = reshape(rbuf(1:NPZ, 1:NY, 1:NX), (/N/))
            call MPI_ISEND(albbuf(1:N, iproc), N, MPI_REAL, iproc, 3, COMM2D_, requests(irq + 3), ierr)

            ibuf(:, 1:NY - 1, 1:NX - 1) = IPHASEPT(:, SY:EY - 1, SX:EX - 1)
            ibuf(:, 1:NY - 1, NX) = IPHASEPT(:, SY:EY - 1, mod(EX - 1, NPXT) + 1)
            ibuf(:, NY, 1:NX - 1) = IPHASEPT(:, mod(EY - 1, NPYT) + 1, SX:EX - 1)
            ibuf(:, NY, NX) = IPHASEPT(:, mod(EY - 1, NPYT) + 1, mod(EX - 1, NPXT) + 1)
            iphbuf(1:N, iproc) = reshape(ibuf(1:NPZ, 1:NY, 1:NX), (/N/))
            call MPI_ISEND(iphbuf(1:N, iproc), N, MPI_INTEGER2, iproc, 4, COMM2D_, requests(irq + 4), ierr)
            irq = irq + 4
        end do
    end if

    N = NPX*NPY*NPZ
    call MPI_IRECV(TEMPP   , N, MPI_REAL    , 0, 1, COMM2D_, requests(irq + 1), ierr)
    call MPI_IRECV(EXTINCTP, N, MPI_REAL    , 0, 2, COMM2D_, requests(irq + 2), ierr)
    call MPI_IRECV(ALBEDOP , N, MPI_REAL    , 0, 3, COMM2D_, requests(irq + 3), ierr)
    call MPI_IRECV(IPHASEP , N, MPI_INTEGER2, 0, 4, COMM2D_, requests(irq + 4), ierr)
    irq = irq + 4
    call MPI_WAITALL(irq, requests, status, ierr)

    if (MYPROC_ .eq. 0) deallocate (tempbuf, extbuf, albbuf, iphbuf, rbuf, ibuf)
    deallocate (requests, status, IXYALL, realbuf)

end subroutine shd_mpi__scatProp

subroutine shd_mpi__bcastSurfSize(MAXSFCPTS, MAXSFCPARS)

    ! Broadcasts the surface array sizes from the master process to the others.

    integer, intent(INOUT) :: MAXSFCPTS, MAXSFCPARS
    integer :: ierr, intbuf(2)

    if (NUMPROC_ <= 1) return

    intbuf(1:2) = (/MAXSFCPTS, MAXSFCPARS/)
    call MPI_BCAST(intbuf, 2, MPI_INTEGER, 0, COMM2D_, ierr)
    MAXSFCPTS = intbuf(1); MAXSFCPARS = intbuf(2)

end subroutine shd_mpi__bcastSurfSize

subroutine shd_mpi__bcastSurfParam(SFCTYPE, NXSFC, NYSFC, NSFCPAR,&
    &DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)

    ! Broadcasts the surface file parameters from the master process to the others.

    character(LEN=2), intent(INOUT) :: SFCTYPE
    integer, intent(INOUT) :: NXSFC, NYSFC, NSFCPAR
    real, intent(INOUT) :: DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO
    integer :: ierr, intbuf(3), n
    real    :: realbuf(4)

    if (NUMPROC_ <= 1) return

    call MPI_BCAST(SFCTYPE, 2, MPI_CHARACTER, 0, COMM2D_, ierr)

    intbuf(1:3) = (/NXSFC, NYSFC, NSFCPAR/)
    call MPI_BCAST(intbuf, 3, MPI_INTEGER, 0, COMM2D_, ierr)
    NXSFC = intbuf(1); NYSFC = intbuf(2); NSFCPAR = intbuf(3)

    realbuf(1:4) = (/DELXSFC, DELYSFC, GNDTEMP, GNDALBEDO/)
    call MPI_BCAST(realbuf, 4, MPI_REAL, 0, COMM2D_, ierr)
    DELXSFC = realbuf(1); DELYSFC = realbuf(2)
    GNDTEMP = realbuf(3); GNDALBEDO = realbuf(4)

    n = NSFCPAR*(NXSFC + 1)*(NYSFC + 1)
    call MPI_BCAST(SFCPARMS, n, MPI_REAL, 0, COMM2D_, ierr)

end subroutine shd_mpi__bcastSurfParam

subroutine shd_mpi__bcastKdistSize(NG, NZCKD)

    ! Broadcasts the k-distribution array sizes from the master process to the others.

    integer, intent(INOUT) :: NG, NZCKD
    integer :: ierr, intbuf(2)

    if (NUMPROC_ <= 1) return

    intbuf(1:2) = (/NG, NZCKD/)
    call MPI_BCAST(intbuf, 2, MPI_INTEGER, 0, COMM2D_, ierr)
    NG = intbuf(1); NZCKD = intbuf(2)

end subroutine shd_mpi__bcastKdistSize

subroutine shd_mpi__bcastKdistParam(SOLFLUX, NG, DELG, NZCKD, ZCKD, KABS)

    ! Broadcasts the k-distribution file parameters from the master process
    ! to the other processes.

    integer, intent(IN) :: NG, NZCKD
    real, intent(INOUT) :: SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG)
    integer :: ierr
    real, allocatable :: realbuf(:)

    if (NUMPROC_ <= 1) return

    allocate (realbuf(NG + NZCKD + 1))
    realbuf(1) = SOLFLUX; realbuf(2:NG + 1) = DELG(:)
    realbuf(NG + 2:NG + NZCKD + 1) = ZCKD(:)
    call MPI_BCAST(realbuf, NG + NZCKD + 1, MPI_REAL, 0, COMM2D_, ierr)
    SOLFLUX = realbuf(1); DELG(:) = realbuf(2:NG + 1)
    ZCKD(:) = realbuf(NG + 2:NG + NZCKD + 1)

    call MPI_BCAST(KABS, NZCKD*NG, MPI_REAL, 0, COMM2D_, ierr)

end subroutine shd_mpi__bcastKdistParam

subroutine shd_mpi__readBcastMemParam(MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO, RUNNAME)

    ! The master process reads the RUNNAME//"_mem_params.inp" file (if there
    ! is one) and broadcasts the memory parameters to all the other processors.

    real, intent(inout) :: MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real, intent(inout) :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(len=*), intent(in) :: RUNNAME
    integer :: Nproc, ierr, i, j
    real, allocatable :: MemParam(:, :)

    allocate (MemParam(4, 0:NUMPROC_ - 1))
    MemParam(1, :) = MAX_TOTAL_MB; MemParam(2, :) = ADAPT_GRID_FACTOR
    MemParam(3, :) = NUM_SH_TERM_FACTOR; MemParam(4, :) = CELL_TO_POINT_RATIO

    if (MYPROC_ .eq. 0) then
        ! See if we can read the memory parameter file and number of processors matches
        open (UNIT=11, FILE=trim(RUNNAME)//'_mem_params.inp', STATUS='OLD', IOSTAT=ierr)
        if (ierr .eq. 0) then
            read (11, *)
            read (11, *) Nproc
            if (Nproc .ne. NUMPROC_) then
                write (6, '(A,A,A)') 'Memory parameter file ', trim(RUNNAME),&
                    &'_mem_params.inp is for a different number of processors.'
            else
                do i = 0, Nproc - 1
                    read (11, *) MemParam(:, i)
                end do
            end if
            close (11)
        end if
    end if

    call MPI_BCAST(MemParam, 4*NUMPROC_, MPI_REAL, 0, COMM2D_, ierr)
    MAX_TOTAL_MB = MemParam(1, MYPROC_)
    ADAPT_GRID_FACTOR = MemParam(2, MYPROC_)
    NUM_SH_TERM_FACTOR = MemParam(3, MYPROC_)
    CELL_TO_POINT_RATIO = MemParam(4, MYPROC_)
    deallocate (MemParam)
    write (6, '(A)') 'Message [READ_BROADCAST_MEM_PARAMS]:'
    write (6, '(2X,A,F7.2,1X,F7.3,2(1X,F5.3))') 'Memory parameters used: ', MAX_TOTAL_MB, ADAPT_GRID_FACTOR,&
        &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

end subroutine shd_mpi__readBcastMemParam

real function shd_mpi__sumCpuTime(cpuTimes)

    ! Sums the CPU time used across all processors

    real cpuTimes
    integer ierr

    if (NUMPROC_ > 1) then
        call MPI_REDUCE(cpuTimes, shd_mpi__sumCpuTime, 1, MPI_REAL, MPI_SUM, 0, COMM2D_, ierr)
    else
        shd_mpi__sumCpuTime = cpuTimes
    end if

end function shd_mpi__sumCpuTime

subroutine shd_mpi__gatherOutput(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
    &IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,&
    &SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
    &XDOMAIN, YDOMAIN, XALLGRID, YALLGRID,&
    &ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,&
    &NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)

    ! Gathers the base grid output data in SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV,
    ! and SUMSHTERMS from all the processors and assembles in the master
    ! processor output arrays ALLFLUXES, ALLFLUXDIV, ALLSHTERMS.  Gathers the
    ! IXYOUT_ array from each processor to determine where each subdomain is.
    ! Also gathers the full base grid position arrays (XALLGRID/YALLGRID) and
    ! sums NCELLS,NPTS,NSH over all processors.

    logical, intent(in) :: FLUX_OUT, FLUXDIV_OUT, SH_OUT
    integer, intent(in) :: IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT
    real, intent(in) :: XDOMAIN, YDOMAIN
    real, intent(in) :: SUMFLUXES(2, NBPTS), SUMDIRFLUX(NBPTS)
    real, intent(in) :: SUMFLUXDIV(NBPTS), SUMSHTERMS(NSHOUT, NBPTS)
    real, intent(out) :: XALLGRID(NXT), YALLGRID(NYT)
    real, intent(out) :: ALLFLUXES(3, NZ, NYT, NXT), ALLFLUXDIV(NZ, NYT, NXT)
    real, intent(out) :: ALLSHTERMS(NSHOUT, NZ, NYT, NXT)
    integer, intent(in) :: NCELLS, NPTS, NSH
    integer, intent(out) :: NCELLSTOT, NPTSTOT, NSHTOT

    integer :: IX, IY, NX, NY, NX1, NY1, MAXN, I, send3int(3), recv3int(3)
    integer :: iproc, ierr, irq
    integer, allocatable :: IXYALL(:, :, :)
    real, allocatable :: sendbuf(:, :), recvbuf(:, :, :), buf(:, :, :, :)

    type(MPI_Request), allocatable :: requests(:)
    type(MPI_Status), allocatable :: status(:)

    allocate (IXYALL(2, 2, 0:NUMPROC_ - 1))
    call MPI_GATHER(IXYOUT_, 4, MPI_INTEGER, IXYALL, 4, MPI_INTEGER, 0, COMM2D_, ierr)

    allocate (requests(NUMPROC_), status(NUMPROC_))

    if (FLUX_OUT) then
        irq = 0
        if (MYPROC_ > 0) then
            allocate (sendbuf(3, NBPTS))
            sendbuf(1:2, :) = SUMFLUXES(:, :); sendbuf(3, :) = SUMDIRFLUX(:)
            irq = irq + 1
            call MPI_ISEND(sendbuf, 3*NBPTS, MPI_REAL, 0, 1, COMM2D_,&
                &requests(irq), ierr)
        else
            MAXN = (maxval(IXYALL(2, 1, :)) + 1)*(maxval(IXYALL(2, 2, :)) + 1)*NZ
            allocate (recvbuf(3, MAXN, 1:NUMPROC_ - 1))
            do iproc = 1, NUMPROC_ - 1
                irq = irq + 1
                call MPI_IRECV(recvbuf(:, :, iproc), 3*MAXN, MPI_REAL, iproc, 1,&
                    &COMM2D_, requests(irq), ierr)
            end do
        end if
        call MPI_WAITALL(irq, requests, status, ierr)

        if (MYPROC_ .eq. 0) then
            IX = IXYALL(1, 1, 0); NX = IXYALL(2, 1, 0)
            IY = IXYALL(1, 2, 0); NY = IXYALL(2, 2, 0)
            NX1 = NX + 1; NY1 = NY + 1
            if (btest(IPFLAG, 0) .and. .not. btest(BCFLAG, 2)) NX1 = NX
            if (btest(IPFLAG, 1) .and. .not. btest(BCFLAG, 3)) NY1 = NY
            allocate (buf(3, NZ, NY1, NX1))
            buf(1:2, :, :, :) = reshape(SUMFLUXES(1:2, :), (/2, NZ, NY1, NX1/))
            buf(3, :, :, :) = reshape(SUMDIRFLUX(:), (/NZ, NY1, NX1/))
            ALLFLUXES(1:3, 1:NZ, IY:IY + NY - 1, IX:IX + NX - 1) = buf(:, 1:NZ, 1:NY, 1:NX)
            deallocate (buf)
            do iproc = 1, NUMPROC_ - 1
                IX = IXYALL(1, 1, iproc); NX = IXYALL(2, 1, iproc)
                IY = IXYALL(1, 2, iproc); NY = IXYALL(2, 2, iproc)
                NX1 = NX + 1; NY1 = NY + 1
                if (btest(IPFLAG, 0) .and. .not. btest(BCFLAG, 2)) NX1 = NX
                if (btest(IPFLAG, 1) .and. .not. btest(BCFLAG, 3)) NY1 = NY
                allocate (buf(3, NZ, NY1, NX1))
                buf(:, :, :, :) = reshape(recvbuf(:, 1:NZ*NY1*NX1, iproc), (/3, NZ, NY1, NX1/))
                ALLFLUXES(1:3, 1:NZ, IY:IY + NY - 1, IX:IX + NX - 1) = buf(:, 1:NZ, 1:NY, 1:NX)
                deallocate (buf)
            end do
            deallocate (recvbuf)
        else
            deallocate (sendbuf)
        end if
    end if

    if (FLUXDIV_OUT) then
        irq = 0
        if (MYPROC_ > 0) then
            irq = irq + 1
            call MPI_ISEND(SUMFLUXDIV, NBPTS, MPI_REAL, 0, 2, COMM2D_,&
                &requests(irq), ierr)
        else
            MAXN = (maxval(IXYALL(2, 1, :)) + 1)*(maxval(IXYALL(2, 2, :)) + 1)*NZ
            allocate (recvbuf(1, MAXN, 1:NUMPROC_ - 1))
            do iproc = 1, NUMPROC_ - 1
                irq = irq + 1
                call MPI_IRECV(recvbuf(:, :, iproc), MAXN, MPI_REAL, iproc, 2,&
                    &COMM2D_, requests(irq), ierr)
            end do
        end if
        call MPI_WAITALL(irq, requests, status, ierr)

        if (MYPROC_ .eq. 0) then
            IX = IXYALL(1, 1, 0); NX = IXYALL(2, 1, 0)
            IY = IXYALL(1, 2, 0); NY = IXYALL(2, 2, 0)
            NX1 = NX + 1; NY1 = NY + 1
            if (btest(IPFLAG, 0) .and. .not. btest(BCFLAG, 2)) NX1 = NX
            if (btest(IPFLAG, 1) .and. .not. btest(BCFLAG, 3)) NY1 = NY
            allocate (buf(1, NZ, NY1, NX1))
            buf(1, :, :, :) = reshape(SUMFLUXDIV(:), (/NZ, NY1, NX1/))
            ALLFLUXDIV(1:NZ, IY:IY + NY - 1, IX:IX + NX - 1) = buf(1, 1:NZ, 1:NY, 1:NX)
            deallocate (buf)
            do iproc = 1, NUMPROC_ - 1
                IX = IXYALL(1, 1, iproc); NX = IXYALL(2, 1, iproc)
                IY = IXYALL(1, 2, iproc); NY = IXYALL(2, 2, iproc)
                NX1 = NX + 1; NY1 = NY + 1
                if (btest(IPFLAG, 0) .and. .not. btest(BCFLAG, 2)) NX1 = NX
                if (btest(IPFLAG, 1) .and. .not. btest(BCFLAG, 3)) NY1 = NY
                allocate (buf(1, NZ, NY1, NX1))
                buf(:, :, :, :) = reshape(recvbuf(:, 1:NZ*NY1*NX1, iproc), (/1, NZ, NY1, NX1/))
                ALLFLUXDIV(1:NZ, IY:IY + NY - 1, IX:IX + NX - 1) = buf(1, 1:NZ, 1:NY, 1:NX)
                deallocate (buf)
            end do
            deallocate (recvbuf)
        end if
    end if

    if (SH_OUT) then
        irq = 0
        if (MYPROC_ > 0) then
            irq = irq + 1
            call MPI_ISEND(SUMSHTERMS, NSHOUT*NBPTS, MPI_REAL, 0, 3, COMM2D_,&
                &requests(irq), ierr)
        else
            MAXN = (maxval(IXYALL(2, 1, :)) + 1)*(maxval(IXYALL(2, 2, :)) + 1)*NZ
            allocate (recvbuf(NSHOUT, MAXN, 1:NUMPROC_ - 1))
            do iproc = 1, NUMPROC_ - 1
                irq = irq + 1
                call MPI_IRECV(recvbuf(:, :, iproc), NSHOUT*MAXN, MPI_REAL, iproc, 3,&
                    &COMM2D_, requests(irq), ierr)
            end do
        end if
        call MPI_WAITALL(irq, requests, status, ierr)

        if (MYPROC_ .eq. 0) then
            IX = IXYALL(1, 1, 0); NX = IXYALL(2, 1, 0)
            IY = IXYALL(1, 2, 0); NY = IXYALL(2, 2, 0)
            NX1 = NX + 1; NY1 = NY + 1
            if (btest(IPFLAG, 0) .and. .not. btest(BCFLAG, 2)) NX1 = NX
            if (btest(IPFLAG, 1) .and. .not. btest(BCFLAG, 3)) NY1 = NY
            allocate (buf(NSHOUT, NZ, NY1, NX1))
            buf(:, :, :, :) = reshape(SUMSHTERMS(:, :), (/NSHOUT, NZ, NY1, NX1/))
            ALLSHTERMS(:, 1:NZ, IY:IY + NY - 1, IX:IX + NX - 1) = buf(:, 1:NZ, 1:NY, 1:NX)
            deallocate (buf)
            do iproc = 1, NUMPROC_ - 1
                IX = IXYALL(1, 1, iproc); NX = IXYALL(2, 1, iproc)
                IY = IXYALL(1, 2, iproc); NY = IXYALL(2, 2, iproc)
                NX1 = NX + 1; NY1 = NY + 1
                if (btest(IPFLAG, 0) .and. .not. btest(BCFLAG, 2)) NX1 = NX
                if (btest(IPFLAG, 1) .and. .not. btest(BCFLAG, 3)) NY1 = NY
                allocate (buf(NSHOUT, NZ, NY1, NX1))
                buf(:, :, :, :) = reshape(recvbuf(:, 1:NZ*NY1*NX1, iproc), (/NSHOUT, NZ, NY1, NX1/))
                ALLSHTERMS(:, 1:NZ, IY:IY + NY - 1, IX:IX + NX - 1) = buf(:, 1:NZ, 1:NY, 1:NX)
                deallocate (buf)
            end do
            deallocate (recvbuf)
        end if
    end if

    ! Make the X and Y grids for the full domain
    if (MYPROC_ .eq. 0) then
        XALLGRID(1:NXT) = (/(XDOMAIN*(I - 1)/NXT, I=1, NXT)/)
        YALLGRID(1:NYT) = (/(YDOMAIN*(I - 1)/NYT, I=1, NYT)/)
    end if

    ! Get the sum of the NCELLS, NPTS, and NSH over the processors
    send3int = (/NCELLS, NPTS, NSH/)
    call MPI_REDUCE(send3int, recv3int, 3, MPI_INTEGER, MPI_SUM, 0, COMM2D_, ierr)
    NCELLSTOT = recv3int(1); NPTSTOT = recv3int(2); NSHTOT = recv3int(3)

    deallocate (IXYALL, requests, status)

end subroutine shd_mpi__gatherOutput

subroutine shd_mpi__findAlbMax(ALBMAX)

    real, intent(INOUT) :: ALBMAX
    integer :: ierr
    real    :: sendrbuf

    if (NUMPROC_ > 1) then
        ! Call MPI_ALLREDUCE to find the max ALBMAX over all processors
        sendrbuf = ALBMAX
        call MPI_ALLREDUCE(sendrbuf, ALBMAX, 1, MPI_REAL, MPI_MAX, COMM2D_, ierr)
        if (ierr .ne. MPI_SUCCESS) call shd_mpi__abort('TOTAL_ALBEDO_MAX: MPI_ALLREDUCE error')
    end if

end subroutine shd_mpi__findAlbMax

subroutine shd_mpi__unifySplit(DOSPLIT, STARTSPLITACC)

    logical, intent(INOUT) :: DOSPLIT
    real, intent(INOUT) :: STARTSPLITACC
    integer :: ierr(2)
    logical :: sendlbuf
    real    :: sendrbuf

    if (NUMPROC_ > 1) then
        ! Call MPI_ALLREDUCE to
        sendrbuf = STARTSPLITACC
        call MPI_ALLREDUCE(sendrbuf, STARTSPLITACC, 1, MPI_REAL, MPI_MAX, COMM2D_, ierr(1))
        sendlbuf = DOSPLIT
        call MPI_ALLREDUCE(sendlbuf, DOSPLIT, 1, MPI_LOGICAL, MPI_LOR, COMM2D_, ierr(2))
        if (any(ierr(:) .ne. MPI_SUCCESS)) call shd_mpi__abort('UNIFY_SPLITTING: MPI_ALLREDUCE error')
    end if

end subroutine shd_mpi__unifySplit

subroutine shd_mpi__findTotSplitCritMax(SPLITCRIT)

    real, intent(INOUT) :: SPLITCRIT
    integer :: ierr
    real    :: sendrbuf

    if (NUMPROC_ > 1) then
        ! Call MPI_ALLREDUCE to find the max SPLITCRIT over all processors
        sendrbuf = SPLITCRIT
        call MPI_ALLREDUCE(sendrbuf, SPLITCRIT, 1, MPI_REAL, MPI_MAX, COMM2D_, ierr)
        if (ierr .ne. MPI_SUCCESS) call shd_mpi__abort('TOTAL_SPLITCRIT_MAX: MPI_ALLREDUCE error')
    end if

end subroutine shd_mpi__findTotSplitCritMax

subroutine shd_mpi__makeDirect(SPT, NPTS, BCFLAG, IPFLAG,&
    &DELTAM, ML, NSTLEG, NLEG,&
    &SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
    &NX, XGRID, NY, YGRID, DIRFLUX)

    ! Makes the direct beam solar flux for the NPTS gridpoints started with
    ! index SPT having X/Y/Z positions in GRIDPOS.
    ! DIRFLUX is set to F*exp(-tau_sun).

    integer, intent(IN) :: SPT, NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NX, NY
    logical, intent(IN) :: DELTAM
    real, intent(IN) :: SOLARFLUX, SOLARMU, SOLARAZ
    real, intent(IN) :: GRIDPOS(3, NPTS), XGRID(NX + 1), YGRID(NY + 1)
    real, intent(OUT) :: DIRFLUX(NPTS)

    logical :: VALIDBEAM
    real    :: UNIFZLEV, UZLEV, DIRPATH
    integer :: MAXNBEAM, Lupx, Lupy, Ldownx, Ldowny, Nbeam
    integer :: SIDE, NINVALID(0:4)
    integer :: IP, I, J, L
    integer :: ierr, irq, iproc
    integer :: intsendbuf(2), intrecvbuf(2)
    logical :: VALIDRAD, BTEST, AllDone
    real    :: XO, YO, ZO, X, Y, Z, XMAX, YMAX, EPS, DIRJUNK
    integer, allocatable :: Ndone(:)
    integer, allocatable :: BndDirIntInfo(:, :, :)
    real, allocatable    :: BndDirRealInfo(:, :, :)
    integer, allocatable :: DirPtrOut(:, :), DirPtrIn(:)
    real, allocatable    :: DirAllOut(:, :), DirAllIn(:)

    type(MPI_Request) :: requests(12)
    type(MPI_Status) :: status(12)

    MAXNBEAM = 2*NPTS
    allocate (BndDirRealInfo(4, MAXNBEAM, 4), BndDirIntInfo(2, MAXNBEAM, 4))
    allocate (DirAllOut(MAXNBEAM, 0:NUMPROC_ - 1), DirPtrOut(MAXNBEAM, 0:NUMPROC_ - 1))
    allocate (DirAllIn(MAXNBEAM), DirPtrIn(MAXNBEAM), Ndone(0:NUMPROC_ - 1))

    EPS = 1.0e-5*max(XGRID(2) - XGRID(1), YGRID(2) - YGRID(1))

    if (SPT .eq. 1) then
        ! Initialize the direct beam calculation
        call shd_comp__computeDirect(1, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
            &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
        ! Find the max lowest uniform Z level over all processors and then
        !   tell shd_comp__computeDirect about it
        UZLEV = UNIFZLEV
        call MPI_ALLREDUCE(UZLEV, UNIFZLEV, 1, MPI_REAL, MPI_MAX, COMM2D_, ierr)
        call shd_comp__computeDirect(2, 0.0, 0.0, 0.0, BCFLAG, IPFLAG,&
            &DELTAM, ML, NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(1),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
    end if

    Lupx = 0; Lupy = 0; Ldownx = 0; Ldowny = 0
    if (btest(BCFLAG, 2)) then
        if (cos(SOLARAZ) .gt. 1.0e-5) then
            Ldownx = 2; Lupx = 1
        else
            Ldownx = 1; Lupx = 2
        end if
        XMAX = XGRID(NX)
    else
        XMAX = XGRID(NX + 1)
    end if
    if (btest(BCFLAG, 3)) then
        if (sin(SOLARAZ) .gt. 1.0e-5) then
            Ldowny = 4; Lupy = 3
        else
            Ldowny = 3; Lupy = 4
        end if
        YMAX = YGRID(NY)
    else
        YMAX = YGRID(NY + 1)
    end if


    ! Do the initial direct beam tracing starting at the base grid points
    NINVALID(:) = 0
    do IP = SPT, NPTS
        DIRPATH = 0.0
        call shd_comp__computeDirect(0, GRIDPOS(1, IP), GRIDPOS(2, IP), GRIDPOS(3, IP),&
            &BCFLAG, IPFLAG, DELTAM, ML, NSTLEG, NLEG,&
            &SOLARFLUX, SOLARMU, SOLARAZ, DIRFLUX(IP),&
            &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
        ! If we don't have a good direct beam then save the information
        !   needed to pass to two upstream processors
        if (.not. VALIDBEAM) then
            if (SIDE .ne. Lupx .and. SIDE .ne. Lupy) call shd_mpi__abort('MAKE_DIRECT_PAR: Bad SIDE 1')
            NINVALID(SIDE) = NINVALID(SIDE) + 1
            if (NINVALID(SIDE) > MAXNBEAM) call shd_mpi__abort('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 1')
            BndDirRealInfo(:, NINVALID(SIDE), SIDE) = (/XO, YO, ZO, DIRPATH/)
            BndDirIntInfo(:, NINVALID(SIDE), SIDE) = (/MYPROC_, IP/)
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
        if (btest(BCFLAG, 2)) then
            call MPI_ISEND(NINVALID(Lupx), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Lupx), 1, COMM2D_, requests(irq + 1), ierr)
            call MPI_ISEND(BndDirRealInfo(:, :, Lupx), 4*NINVALID(Lupx), MPI_REAL,&
                &IPROCNEIGH_(Lupx), 2, COMM2D_, requests(irq + 2), ierr)
            call MPI_ISEND(BndDirIntInfo(:, :, Lupx), 2*NINVALID(Lupx), MPI_INTEGER,&
                &IPROCNEIGH_(Lupx), 3, COMM2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        if (btest(BCFLAG, 3)) then
            call MPI_ISEND(NINVALID(Lupy), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Lupy), 4, COMM2D_, requests(irq + 1), ierr)
            call MPI_ISEND(BndDirRealInfo(:, :, Lupy), 4*NINVALID(Lupy), MPI_REAL,&
                &IPROCNEIGH_(Lupy), 5, COMM2D_, requests(irq + 2), ierr)
            call MPI_ISEND(BndDirIntInfo(:, :, Lupy), 2*NINVALID(Lupy), MPI_INTEGER,&
                &IPROCNEIGH_(Lupy), 6, COMM2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        if (btest(BCFLAG, 2)) then
            call MPI_IRECV(NINVALID(Ldownx), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Ldownx), 1, COMM2D_, requests(irq + 1), ierr)
            call MPI_IRECV(BndDirRealInfo(:, :, Ldownx), 4*MAXNBEAM, MPI_REAL,&
                &IPROCNEIGH_(Ldownx), 2, COMM2D_, requests(irq + 2), ierr)
            call MPI_IRECV(BndDirIntInfo(:, :, Ldownx), 2*MAXNBEAM, MPI_INTEGER,&
                &IPROCNEIGH_(Ldownx), 3, COMM2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        if (btest(BCFLAG, 3)) then
            call MPI_IRECV(NINVALID(Ldowny), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Ldowny), 4, COMM2D_, requests(irq + 1), ierr)
            call MPI_IRECV(BndDirRealInfo(:, :, Ldowny), 4*MAXNBEAM, MPI_REAL,&
                &IPROCNEIGH_(Ldowny), 5, COMM2D_, requests(irq + 2), ierr)
            call MPI_IRECV(BndDirIntInfo(:, :, Ldowny), 2*MAXNBEAM, MPI_INTEGER,&
                &IPROCNEIGH_(Ldowny), 6, COMM2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end if
        call MPI_WAITALL(irq, requests, status, ierr)

        NINVALID(Lupx) = 0; NINVALID(Lupy) = 0
        if (btest(BCFLAG, 2)) then
            ! Continue the backwards ray integrations from the X boundary
            do J = 1, NINVALID(Ldownx)
                if (Ldownx .eq. 1) BndDirRealInfo(1, J, Ldownx) = XGRID(1)
                if (Ldownx .eq. 2) BndDirRealInfo(1, J, Ldownx) = XGRID(NX)
                X = BndDirRealInfo(1, J, Ldownx); Y = BndDirRealInfo(2, J, Ldownx)
                Z = BndDirRealInfo(3, J, Ldownx); DIRPATH = BndDirRealInfo(4, J, Ldownx)
                if (Y < YGRID(1) - EPS) Y = YMAX; if (Y > YMAX + EPS) Y = YGRID(1)
                call shd_comp__computeDirect(0, X, Y, Z, BCFLAG, IPFLAG,&
                    &DELTAM, ML, NSTLEG, NLEG,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, DIRJUNK,&
                    &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                ! If we got a valid direct beam then store it in DirAllOut otherwise
                ! save the information to continue passing to neighboring processors
                if (VALIDBEAM) then
                    iproc = BndDirIntInfo(1, J, Ldownx)
                    Ndone(iproc) = Ndone(iproc) + 1
                    if (Ndone(iproc) > MAXNBEAM)&
                        &call shd_mpi__abort('MAKE_DIRECT_PAR: Ndone(iproc) exceeds MAXNBEAM x')
                    DirAllOut(Ndone(iproc), iproc) = DIRPATH
                    DirPtrOut(Ndone(iproc), iproc) = BndDirIntInfo(2, J, Ldownx)
                else
                    if (SIDE .ne. Lupx .and. SIDE .ne. Lupy) call shd_mpi__abort('MAKE_DIRECT_PAR: Bad SIDE 2x')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) > MAXNBEAM)&
                        &call shd_mpi__abort('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 2x')
                    BndDirRealInfo(:, NINVALID(SIDE), SIDE) = (/XO, YO, ZO, DIRPATH/)
                    BndDirIntInfo(:, NINVALID(SIDE), SIDE) = BndDirIntInfo(:, J, Ldownx)
                end if
            end do
        end if

        if (btest(BCFLAG, 3)) then
            ! Continue the backwards ray integrations for the Y boundary
            do J = 1, NINVALID(Ldowny)
                if (Ldowny .eq. 3) BndDirRealInfo(2, J, Ldowny) = YGRID(1)
                if (Ldowny .eq. 4) BndDirRealInfo(2, J, Ldowny) = YGRID(NY)
                X = BndDirRealInfo(1, J, Ldowny); Y = BndDirRealInfo(2, J, Ldowny)
                Z = BndDirRealInfo(3, J, Ldowny); DIRPATH = BndDirRealInfo(4, J, Ldowny)
                if (X < XGRID(1) - EPS) X = XMAX; if (X > XMAX + EPS) X = XGRID(1)
                call shd_comp__computeDirect(0, X, Y, Z, BCFLAG, IPFLAG,&
                    &DELTAM, ML, NSTLEG, NLEG,&
                    &SOLARFLUX, SOLARMU, SOLARAZ, DIRJUNK,&
                    &UNIFZLEV, XO, YO, ZO, DIRPATH, SIDE, VALIDBEAM)
                ! If we got a valid direct beam then store it in DirAllOut otherwise
                ! save the information to continue passing to neighboring processors
                if (VALIDBEAM) then
                    iproc = BndDirIntInfo(1, J, Ldowny)
                    Ndone(iproc) = Ndone(iproc) + 1
                    if (Ndone(iproc) > MAXNBEAM)&
                        &call shd_mpi__abort('MAKE_DIRECT_PAR: Ndone(iproc) exceeds MAXNBEAM y')
                    DirAllOut(Ndone(iproc), iproc) = DIRPATH
                    DirPtrOut(Ndone(iproc), iproc) = BndDirIntInfo(2, J, Ldowny)
                else
                    if (SIDE .ne. Lupx .and. SIDE .ne. Lupy) call shd_mpi__abort('MAKE_DIRECT_PAR: Bad SIDE 2y')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) > MAXNBEAM)&
                        &call shd_mpi__abort('MAKE_DIRECT_PAR: NINVALID(SIDE) exceeds MAXNBEAM 2x')
                    BndDirRealInfo(:, NINVALID(SIDE), SIDE) = (/XO, YO, ZO, DIRPATH/)
                    BndDirIntInfo(:, NINVALID(SIDE), SIDE) = BndDirIntInfo(:, J, Ldowny)
                end if
            end do
        end if

        ! See if there are any more invalid radiance rays on all the processors
        intsendbuf = (/NINVALID(Lupx), NINVALID(Lupy)/)
        call MPI_ALLREDUCE(intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX,&
            &COMM2D_, ierr)
        AllDone = all(intrecvbuf(1:2) .eq. 0)
        if (maxval(intrecvbuf(1:2)) > MAXNBEAM .or. maxval(Ndone(:)) > MAXNBEAM) then
            call shd_mpi__abort('MAKE_DIRECT_PAR: MAXNBEAM exceeded 2')
        end if
    end do

    ! Exchange the finished direct beam paths with all of the other processors
    do iproc = 0, NUMPROC_ - 1
        if (iproc .ne. MYPROC_) then
            call MPI_ISEND(Ndone(iproc), 1, MPI_INTEGER,&
                &iproc, 7, COMM2D_, requests(1), ierr)
            call MPI_ISEND(DirAllOut(:, iproc), Ndone(iproc), MPI_REAL,&
                &iproc, 8, COMM2D_, requests(2), ierr)
            call MPI_ISEND(DirPtrOut(:, iproc), Ndone(iproc), MPI_INTEGER,&
                &iproc, 9, COMM2D_, requests(3), ierr)
            call MPI_IRECV(Nbeam, 1, MPI_INTEGER,&
                &iproc, 7, COMM2D_, requests(4), ierr)
            call MPI_IRECV(DirAllIn, MAXNBEAM, MPI_REAL,&
                &iproc, 8, COMM2D_, requests(5), ierr)
            call MPI_IRECV(DirPtrIn, MAXNBEAM, MPI_INTEGER,&
                &iproc, 9, COMM2D_, requests(6), ierr)
            call MPI_WAITALL(6, requests, status, ierr)
            ! Put the direct beam flux in the correct grid point
            do I = 1, Nbeam
                IP = DirPtrIn(I)
                DIRFLUX(IP) = SOLARFLUX*exp(-DirAllIn(I))
            end do
        else
            ! If we have the direct paths for this processor then avoid MPI and
            !   put the direct beam flux in the grid point
            do I = 1, Ndone(MYPROC_)
                IP = DirPtrOut(I, MYPROC_)
                DIRFLUX(IP) = SOLARFLUX*exp(-DirAllOut(I, MYPROC_))
            end do
        end if
    end do

    deallocate (DirAllOut, DirPtrOut, DirAllIn, DirPtrIn, Ndone)
    deallocate (BndDirRealInfo, BndDirIntInfo)

end subroutine shd_mpi__makeDirect

subroutine shd_mpi__findBoundaryPnt(BCFLAG, IPFLAG, NPTS, SWEEPORD, GRIDPTR,&
    &GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)

    ! Makes the pointers to and positions of the horizontal boundary grid
    ! points of the processors neighboring this one.
    ! For each of the potentially eight octants there are
    ! two sets of boundary pointers; 1 is for the current subdomain,
    ! and 2 is for the neighboring subdomain.
    ! Fills in the BXYPTR_(:,IBC,JOCT) and BXYPOS_(:,IBC,JOCT) arrays
    ! (module variables), where IBC is a boundary point index and JOCT is
    ! a sort of octant number. BXYPTR_ holds (iproc,IPT,SIDE) for the boundary
    ! grid points from grid point IPT in processor iproc that are on the SIDE
    ! boundary of this subdomain (SIDE: 1=-X, 2=+X, 3=-Y, 4=+Y).  BXYPOS_ holds
    ! the gridpoint positions (X,Y,Z) for these boundary points.

    integer, intent(IN) :: BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS, *), GRIDPTR(8, *)
    integer, intent(IN) :: NX, NY, NZ
    real, intent(IN) :: GRIDPOS(3, NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)

    integer :: tag, ierr, irq
    integer :: MAXNBXY, NOCT, JOCT, IORDER, IPCELL, ICORNER, IPT, SIDE
    integer :: L, LI, N(4), N1(8, 4), N2(8, 4), NS, NBX, NBY
    integer :: MAXBC, IBC, BITZ, SZ, EZ, DZ, LX, LY, IX, IY, IZ, J
    integer :: INOUT(4) = (/2, 1, 4, 3/)
    integer, save :: OLDNPTS = 0
    logical :: SAMEGRID
    real    :: XEDGE(8), YEDGE(8), XYB(4)
    integer, allocatable :: BNDPTR(:, :, :, :), INDX(:)
    real, allocatable :: BNDPOS(:, :, :, :, :), ZARRAY(:)

    type(MPI_Request) :: requests(24)
    type(MPI_Status) :: status(24)

    if (.not. (btest(BCFLAG, 2) .or. btest(BCFLAG, 3))) return
    if (btest(IPFLAG, 0) .and. btest(IPFLAG, 1)) return
    call MPI_ALLREDUCE(NPTS .eq. OLDNPTS, SAMEGRID, 1, MPI_LOGICAL, MPI_LAND, COMM2D_, ierr)
    OLDNPTS = NPTS
    if (SAMEGRID) return

    ! IOCT=1+BITX+2*BITY+4*BITZ ; JOCT=JOCTORDER(IOCT) = (/ 1,3,5,7,2,4,6,8 /)
    ! CX<0: BITX=1, JOCT=3,7,4,8      CY<0: BITY=1, JOCT=5,7,6,8
    XEDGE(:) = (/XGRID(NX), XGRID(NX), XGRID(1), XGRID(1), XGRID(NX), XGRID(NX), XGRID(1), XGRID(1)/)
    YEDGE(:) = (/YGRID(NY), YGRID(NY), YGRID(NY), YGRID(NY), YGRID(1), YGRID(1), YGRID(1), YGRID(1)/)
    XYB(1) = XGRID(1); XYB(2) = XGRID(NX)
    XYB(3) = YGRID(1); XYB(4) = YGRID(NY)

    if (btest(IPFLAG, 1) .and. btest(IPFLAG, 0)) then
        NOCT = 2
    else if (btest(IPFLAG, 1)) then
        NOCT = 4
    else
        NOCT = 8
    end if
    ! Allocate grid point pointer and position arrays (second from last
    !   index is 1 for local, 2 for neighbor; last index is boundary L)
    if (allocated(IBXYPTR_)) then
        MAXNBXY = nint(3.0*maxval(IBXYPTR_(2, :, :)))
    else
        MAXNBXY = 3.0*(NZ - 2)*max(NX, NY)
    end if
    allocate (BNDPTR(NOCT, MAXNBXY, 2, 4), BNDPOS(3, NOCT, MAXNBXY, 2, 4))

    N1(:, :) = 0
    do JOCT = 1, NOCT
        BITZ = 1 - ibits(JOCT, 0, 1)
        N(:) = 0
        do IORDER = 1, NPTS
            IPCELL = ishft(SWEEPORD(IORDER, JOCT), -3)
            ICORNER = ibits(SWEEPORD(IORDER, JOCT), 0, 3) + 1
            IPT = GRIDPTR(ICORNER, IPCELL)

            if ((GRIDPOS(3, IPT) < ZGRID(NZ) .and. BITZ .eq. 0) .or.&
                &(GRIDPOS(3, IPT) > ZGRID(1) .and. BITZ .eq. 1)) then
                if (btest(BCFLAG, 2) .and. GRIDPOS(1, IPT) .eq. XEDGE(JOCT)&
                    &.and. .not. btest(IPFLAG, 0)) then
                    SIDE = 2 - ibits(JOCT - 1, 1, 1)
                    N(SIDE) = N(SIDE) + 1
                    if (N(SIDE) > MAXNBXY) then
                        call shd_mpi__abort('FIND_BOUNDARY_POINTS: MAXNBXY exceeded by N(SIDE) 1x')
                    end if
                    BNDPTR(JOCT, N(SIDE), 1, SIDE) = IPT
                    BNDPOS(:, JOCT, N(SIDE), 1, SIDE) = GRIDPOS(:, IPT)
                end if
                if (btest(BCFLAG, 3) .and. GRIDPOS(2, IPT) .eq. YEDGE(JOCT)&
                    &.and. .not. btest(IPFLAG, 1)) then
                    SIDE = 4 - ibits(JOCT - 1, 2, 1)
                    N(SIDE) = N(SIDE) + 1
                    if (N(SIDE) > MAXNBXY) then
                        call shd_mpi__abort('FIND_BOUNDARY_POINTS: MAXNBXY exceeded by N(SIDE) 1y')
                    end if
                    BNDPTR(JOCT, N(SIDE), 1, SIDE) = IPT
                    BNDPOS(:, JOCT, N(SIDE), 1, SIDE) = GRIDPOS(:, IPT)
                end if
            end if
        end do
        N1(JOCT, :) = N(:)
    end do

    NS = maxval(N1(:, :))
    if (NS > MAXNBXY) call shd_mpi__abort('FIND_BOUNDARY_POINTS: MAXNBXY might be exceeded')
    N2(:, :) = 0
    irq = 0
    do L = 1, 4
        if ((btest(BCFLAG, 2) .and. (L .eq. 1 .or. L .eq. 2)) .or.&
            &(btest(BCFLAG, 3) .and. (L .eq. 3 .or. L .eq. 4))) then
            ! Exchange the lists of boundary points with neighbors
            LI = INOUT(L)
            tag = L
            call MPI_ISEND(N1(:, L), 8, MPI_INTEGER,&
                &IPROCNEIGH_(L), tag, COMM2D_, requests(irq + 1), ierr)
            call MPI_IRECV(N2(:, LI), 8, MPI_INTEGER,&
                &IPROCNEIGH_(LI), tag, COMM2D_, requests(irq + 2), ierr)
            tag = 1000 + L
            call MPI_ISEND(BNDPTR(:, :, 1, L), NOCT*NS, MPI_INTEGER,&
                &IPROCNEIGH_(L), tag, COMM2D_, requests(irq + 3), ierr)
            call MPI_IRECV(BNDPTR(:, :, 2, LI), NOCT*MAXNBXY, MPI_INTEGER,&
                &IPROCNEIGH_(LI), tag, COMM2D_, requests(irq + 4), ierr)
            tag = 2000 + L
            call MPI_ISEND(BNDPOS(:, :, :, 1, L), 3*NOCT*NS, MPI_REAL,&
                &IPROCNEIGH_(L), tag, COMM2D_, requests(irq + 5), ierr)
            call MPI_IRECV(BNDPOS(:, :, :, 2, LI), 3*NOCT*MAXNBXY, MPI_REAL,&
                &IPROCNEIGH_(LI), tag, COMM2D_, requests(irq + 6), ierr)
            irq = irq + 6
        end if
    end do
    ! Wait for the transfers among all the neighbors to resolve
    call MPI_WAITALL(irq, requests, status, ierr)
    NS = maxval(N2(:, :))
    if (NS > MAXNBXY) call shd_mpi__abort('FIND_BOUNDARY_POINTS: MAXNBXY was exceeded')

    if (allocated(BXYPTR_)) deallocate (BXYPTR_, BXYPOS_, IBXYPTR_)
    MAXBC = maxval((/N2(:, 1) + N2(:, 3), N2(:, 1) + N2(:, 4), N2(:, 2) + N2(:, 3), N2(:, 2) + N2(:, 4)/))
    allocate (BXYPTR_(3, MAXBC, NOCT), BXYPOS_(3, MAXBC, NOCT), IBXYPTR_(2, NZ, NOCT))

    ! Put the boundary points from the neighboring subdomains in the
    ! BXYPTR_ and BXYPOS_ arrays.  Merge in X and Y boundary points one Z slab
    ! at a time, and set pointers (IBXYPTR_) to start and end of each slab.
    !  DO IBC = IBXYPTR_(1,IZ,JOCT), IBXYPTR_(2,IZ,JOCT)
    !    X=BXYPOS_(1,IBC,JOCT) ; Y=BXYPOS_(2,IBC,JOCT) ; Z=BXYPOS_(3,IBC,JOCT)
    !    iproc=BXYPTR_(1,IBC,JOCT) ; IPT=BXYPTR_(2,IBC,JOCT) ; SIDE=BXYPTR_(2,IBC,JOCT)

    do JOCT = 1, NOCT
        BITZ = 1 - ibits(JOCT, 0, 1)
        DZ = 2*BITZ - 1
        LX = 1 + ibits(JOCT - 1, 1, 1)
        LY = 3 + ibits(JOCT - 1, 2, 1)
        ! Combine the indices to the X and Y boundary points and sort on Z
        NBX = 0; if (btest(BCFLAG, 2)) NBX = N2(JOCT, LX)
        NBY = 0; if (btest(BCFLAG, 3)) NBY = N2(JOCT, LY)
        allocate (ZARRAY(NBX + NBY), INDX(NBX + NBY))
        ZARRAY(1:NBX) = BNDPOS(3, JOCT, 1:NBX, 2, LX)
        INDX(1:NBX) = (/(J, J=1, NBX)/)
        ZARRAY(NBX + 1:NBX + NBY) = BNDPOS(3, JOCT, 1:NBY, 2, LY)
        INDX(NBX + 1:NBX + NBY) = (/(-J, J=1, NBY)/)
        call shd_util__ssort_i(ZARRAY, INDX, NBX + NBY, 2*DZ)

        ! For each Z slab use Z sorted index to transfer the X or Y boundary points
        IBXYPTR_(:, :, JOCT) = 0
        SZ = (NZ - 1)*(1 - BITZ) + 1 + DZ
        EZ = (NZ - 1)*(BITZ) + 1
        IBC = 0; J = 1
        do IZ = SZ, EZ, DZ
            IBXYPTR_(1, IZ, JOCT) = IBC + 1
            ! Store the boundary point info for this Z slab
            do while (J <= NBX + NBY .and.&
                &((ZARRAY(J) >= ZGRID(IZ) .and. BITZ .eq. 0) .or.&
                &(ZARRAY(J) <= ZGRID(IZ) .and. BITZ .eq. 1)))
                IBC = IBC + 1
                if (INDX(J) > 0) then
                    IX = INDX(J)
                    BXYPTR_(:, IBC, JOCT) = (/IPROCNEIGH_(LX), BNDPTR(JOCT, IX, 2, LX), LX/)
                    BXYPOS_(:, IBC, JOCT) = BNDPOS(:, JOCT, IX, 2, LX)
                    BXYPOS_(int((LX + 1)/2), IBC, JOCT) = XYB(LX)
                else
                    IY = abs(INDX(J))
                    BXYPTR_(:, IBC, JOCT) = (/IPROCNEIGH_(LY), BNDPTR(JOCT, IY, 2, LY), LY/)
                    BXYPOS_(:, IBC, JOCT) = BNDPOS(:, JOCT, IY, 2, LY)
                    BXYPOS_(int((LY + 1)/2), IBC, JOCT) = XYB(LY)
                end if
                J = J + 1
            end do
            IBXYPTR_(2, IZ, JOCT) = IBC
        end do
        deallocate (ZARRAY, INDX)
    end do

    deallocate (BNDPTR, BNDPOS)

end subroutine shd_mpi__findBoundaryPnt

subroutine shd_mpi__calBoundaryRad(BCFLAG, IPFLAG, JOCT, IZ,&
    &NX, NY, NZ, XGRID, YGRID, ZGRID,&
    &NA, NPTS, NCELLS, GRIDPTR,&
    &NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,&
    &MU, PHI, EXTINCT, NSTOKES, SOURCE,&
    &KANG, GRIDRAD)

    ! Calculates the discrete ordinate boundary radiances for the KANG discrete
    ! ordinate (going in direction MU/PHI) for the grid points in the IZ slab.
    ! FIND_BOUNDARY_POINTS has already been called to set up the BXYPTR_ and
    ! BXYPOS_ arrays for the grid points on the boundaries, and the IBXYPTR_
    ! array which has the starting and ending points for the IZ slab.  The
    ! grid points are done in the sweeping order, and so the octant (JOCT)
    ! is used to index into the IBXYPTR_, BXYPTR_, and BXYPOS_ arrays.  The
    ! grid points specified in BXYPTR_/BXYPOS_ are on the boundaries of this
    ! processor, and need to be traced back along the discrete ordinate
    ! first within this subdomain.  The radiative transfer integrations
    ! with the extinction (EXTINCT) and source function (SOURCE) are done
    ! with the INTEGRATE_RAY routine.  If the top or bottom boundary is
    ! hit when tracing back then the radiance is saved, otherwise the
    ! required information (X,Y,Z, transmission, partial radiance, original
    ! processor, and starting grid point) are saved to be sent to the
    ! appropriate neighboring processor.  These data are received from the
    ! neighboring processors and the radiative transfer integrations are
    ! continued, looping until all the rays have been finished.  Finally,
    ! the finished radiance and original grid points are exchanged with
    ! all the other processors.

    integer, intent(IN) :: BCFLAG, IPFLAG, JOCT, IZ
    integer, intent(IN) :: NSTOKES, NX, NY, NZ, NA, NPTS, NCELLS, KANG
    real, intent(IN) :: XGRID(NX), YGRID(NY), ZGRID(NZ)
    integer, intent(IN) :: GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS)
    integer, intent(IN) :: TREEPTR(2, NCELLS)
    integer*2, intent(IN) :: CELLFLAGS(NCELLS)
    real, intent(IN) :: GRIDPOS(3, NPTS)
    real, intent(IN) :: MU, PHI
    real, intent(IN) :: EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real, intent(INOUT) :: GRIDRAD(NSTOKES, NPTS)

    integer :: MAXNBND, Nbcpnts, Lupx, Lupy, Ldownx, Ldowny
    integer :: SIDE, NINVALID(0:4)
    integer :: IBC, IPT, I, J, L
    integer :: ierr(12), irq, iproc
    integer :: intsendbuf(2), intrecvbuf(2)
    logical :: VALIDRAD, BTEST, AllDone
    real    :: X, Y, Z, RADIANCE(NSTOKES), TRANSMIT, xo, yo, zo
    integer, allocatable :: Ndone(:)
    integer, allocatable :: BndRadIntInfo(:, :, :)
    real, allocatable    :: BndRadRealInfo(:, :, :)
    integer, allocatable :: RadPtrOut(:, :), RadPtrIn(:)
    real, allocatable    :: RadAllOut(:, :, :), RadAllIn(:, :)

    type(MPI_Request) :: requests(12)
    type(MPI_Status) :: status(12)

    MAXNBND = NUMPROC_*maxval(IBXYPTR_(2, :, :) - IBXYPTR_(1, :, :) + 1)
    allocate (BndRadRealInfo(4 + NSTOKES, MAXNBND, 4), BndRadIntInfo(2, MAXNBND, 4))
    allocate (RadAllOut(NSTOKES, MAXNBND, 0:NUMPROC_ - 1), RadPtrOut(MAXNBND, 0:NUMPROC_ - 1))
    allocate (RadAllIn(NSTOKES, MAXNBND), RadPtrIn(MAXNBND), Ndone(0:NUMPROC_ - 1))

    Lupx = 0; Lupy = 0; Ldownx = 0; Ldowny = 0
    if (btest(BCFLAG, 2)) then
        if (cos(PHI) .gt. 1.0e-5) then
            Ldownx = 2; Lupx = 1
        else
            Ldownx = 1; Lupx = 2
        end if
    end if
    if (btest(BCFLAG, 3)) then
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
    do IBC = IBXYPTR_(1, IZ, JOCT), IBXYPTR_(2, IZ, JOCT)
        ! Call INTEGRATE_RAY for this
        X = BXYPOS_(1, IBC, JOCT); Y = BXYPOS_(2, IBC, JOCT); Z = BXYPOS_(3, IBC, JOCT)
        iproc = BXYPTR_(1, IBC, JOCT); IPT = BXYPTR_(2, IBC, JOCT); SIDE = BXYPTR_(3, IBC, JOCT)
        RADIANCE(:) = 0.0; TRANSMIT = 1.0
        xo = x; yo = y; zo = z
        call shd_mpi__integrateRay(BCFLAG, IPFLAG, NX, NY, NZ, NA,&
            &XGRID, YGRID, ZGRID, NPTS, NCELLS,&
            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
            &GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD,&
            &MU, PHI, X, Y, Z, SIDE, TRANSMIT, RADIANCE, VALIDRAD)
        ! If we got a radiance then store it otherwise save
        !  the information needed to pass to other processors
        if (VALIDRAD) then
            Ndone(iproc) = Ndone(iproc) + 1
            RadAllOut(:, Ndone(iproc), iproc) = RADIANCE(:)
            RadPtrOut(Ndone(iproc), iproc) = IPT
        else
            if (SIDE .ne. Lupx .and. SIDE .ne. Lupy) then
                write (6, *) 'Bad SIDE 1:', MYPROC_, mu, phi, SIDE, x, y, z, xo, yo, zo
                call shd_mpi__abort('CALC_BOUNDARY_RADIANCES')
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
        if (btest(BCFLAG, 2)) then
            call MPI_ISEND(NINVALID(Lupx), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Lupx), 1, COMM2D_, requests(irq + 1), ierr(irq + 1))
            call MPI_ISEND(BndRadRealInfo(:, :, Lupx), (4 + NSTOKES)*NINVALID(Lupx), MPI_REAL,&
                &IPROCNEIGH_(Lupx), 2, COMM2D_, requests(irq + 2), ierr(irq + 2))
            call MPI_ISEND(BndRadIntInfo(:, :, Lupx), 2*NINVALID(Lupx), MPI_INTEGER,&
                &IPROCNEIGH_(Lupx), 3, COMM2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (btest(BCFLAG, 3)) then
            call MPI_ISEND(NINVALID(Lupy), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Lupy), 4, COMM2D_, requests(irq + 1), ierr(irq + 1))
            call MPI_ISEND(BndRadRealInfo(:, :, Lupy), (4 + NSTOKES)*NINVALID(Lupy), MPI_REAL,&
                &IPROCNEIGH_(Lupy), 5, COMM2D_, requests(irq + 2), ierr(irq + 2))
            call MPI_ISEND(BndRadIntInfo(:, :, Lupy), 2*NINVALID(Lupy), MPI_INTEGER,&
                &IPROCNEIGH_(Lupy), 6, COMM2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (btest(BCFLAG, 2)) then
            call MPI_IRECV(NINVALID(Ldownx), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Ldownx), 1, COMM2D_, requests(irq + 1), ierr(irq + 1))
            call MPI_IRECV(BndRadRealInfo(:, :, Ldownx), (4 + NSTOKES)*MAXNBND, MPI_REAL,&
                &IPROCNEIGH_(Ldownx), 2, COMM2D_, requests(irq + 2), ierr(irq + 2))
            call MPI_IRECV(BndRadIntInfo(:, :, Ldownx), 2*MAXNBND, MPI_INTEGER,&
                &IPROCNEIGH_(Ldownx), 3, COMM2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (btest(BCFLAG, 3)) then
            call MPI_IRECV(NINVALID(Ldowny), 1, MPI_INTEGER,&
                &IPROCNEIGH_(Ldowny), 4, COMM2D_, requests(irq + 1), ierr(irq + 1))
            call MPI_IRECV(BndRadRealInfo(:, :, Ldowny), (4 + NSTOKES)*MAXNBND, MPI_REAL,&
                &IPROCNEIGH_(Ldowny), 5, COMM2D_, requests(irq + 2), ierr(irq + 2))
            call MPI_IRECV(BndRadIntInfo(:, :, Ldowny), 2*MAXNBND, MPI_INTEGER,&
                &IPROCNEIGH_(Ldowny), 6, COMM2D_, requests(irq + 3), ierr(irq + 3))
            irq = irq + 3
        end if
        if (any(ierr(1:irq) .ne. MPI_SUCCESS)) then
            write (6, *) 'MPI error: ', ierr(1:irq)
            call shd_mpi__abort('CALC_BOUNDARY_RADIANCES: MPI error 1')
        end if
        call MPI_WAITALL(irq, requests, status, ierr(1))
        if (ierr(1) .ne. MPI_SUCCESS) then
            write (6, *) 'MPI_WAITALL error: ', ierr(1)
            call shd_mpi__abort('CALC_BOUNDARY_RADIANCES: MPI error 2')
        end if

        NINVALID(Lupx) = 0; NINVALID(Lupy) = 0
        if (btest(BCFLAG, 2)) then
            ! Continue the backwards ray integrations for the X boundary
            do J = 1, NINVALID(Ldownx)
                if (Ldownx .eq. 1) BndRadRealInfo(1, J, Ldownx) = XGRID(1)
                if (Ldownx .eq. 2) BndRadRealInfo(1, J, Ldownx) = XGRID(NX)
                X = BndRadRealInfo(1, J, Ldownx); Y = BndRadRealInfo(2, J, Ldownx)
                Z = BndRadRealInfo(3, J, Ldownx)
                SIDE = Ldownx; TRANSMIT = BndRadRealInfo(4, J, Ldownx)
                RADIANCE(:) = BndRadRealInfo(5:4 + NSTOKES, J, Ldownx)
                call shd_mpi__integrateRay(BCFLAG, IPFLAG, NX, NY, NZ, NA,&
                    &XGRID, YGRID, ZGRID, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD,&
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
                    if (SIDE .ne. Lupx .and. SIDE .ne. Lupy) then
                        write (6, *) 'Bad SIDE 2x:', MYPROC_, mu, phi, SIDE, x, y, z
                        call shd_mpi__abort('CALC_BOUNDARY_RADIANCES')
                    end if
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/X, Y, Z, TRANSMIT, RADIANCE(:)/)
                    BndRadIntInfo(:, NINVALID(SIDE), SIDE) = BndRadIntInfo(:, J, Ldownx)
                end if
            end do
        end if

        if (btest(BCFLAG, 3)) then
            ! Continue the backwards ray integrations for the Y boundary
            do J = 1, NINVALID(Ldowny)
                if (Ldowny .eq. 3) BndRadRealInfo(2, J, Ldowny) = YGRID(1)
                if (Ldowny .eq. 4) BndRadRealInfo(2, J, Ldowny) = YGRID(NY)
                X = BndRadRealInfo(1, J, Ldowny); Y = BndRadRealInfo(2, J, Ldowny)
                Z = BndRadRealInfo(3, J, Ldowny)
                SIDE = Ldowny; TRANSMIT = BndRadRealInfo(4, J, Ldowny)
                RADIANCE(:) = BndRadRealInfo(5:4 + NSTOKES, J, Ldowny)
                call shd_mpi__integrateRay(BCFLAG, IPFLAG, NX, NY, NZ, NA,&
                    &XGRID, YGRID, ZGRID, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG, GRIDRAD,&
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
                    if (SIDE .ne. Lupx .and. SIDE .ne. Lupy) then
                        write (6, *) 'Bad SIDE 2y:', MYPROC_, mu, phi, SIDE, x, y, z
                        call shd_mpi__abort('CALC_BOUNDARY_RADIANCES')
                    end if
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/X, Y, Z, TRANSMIT, RADIANCE(:)/)
                    BndRadIntInfo(:, NINVALID(SIDE), SIDE) = BndRadIntInfo(:, J, Ldowny)
                end if
            end do
        end if

        ! See if there are any more invalid radiance rays on all the processors
        intsendbuf = (/NINVALID(Lupx), NINVALID(Lupy)/)
        call MPI_ALLREDUCE(intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX,&
            &COMM2D_, ierr(1))
        AllDone = all(intrecvbuf(1:2) .eq. 0)
    end do

    ! Exchange the finished boundary radiances with all of the other processors
    do iproc = 0, NUMPROC_ - 1
        if (iproc .ne. MYPROC_) then
            call MPI_ISEND(Ndone(iproc), 1, MPI_INTEGER,&
                &iproc, 7, COMM2D_, requests(1), ierr(1))
            call MPI_ISEND(RadAllOut(:, :, iproc), NSTOKES*Ndone(iproc), MPI_REAL,&
                &iproc, 8, COMM2D_, requests(2), ierr(2))
            call MPI_ISEND(RadPtrOut(:, iproc), Ndone(iproc), MPI_INTEGER,&
                &iproc, 9, COMM2D_, requests(3), ierr(3))
            call MPI_IRECV(Nbcpnts, 1, MPI_INTEGER,&
                &iproc, 7, COMM2D_, requests(4), ierr(4))
            call MPI_IRECV(RadAllIn, NSTOKES*MAXNBND, MPI_REAL,&
                &iproc, 8, COMM2D_, requests(5), ierr(5))
            call MPI_IRECV(RadPtrIn, MAXNBND, MPI_INTEGER,&
                &iproc, 9, COMM2D_, requests(6), ierr(6))
            call MPI_WAITALL(6, requests, status, ierr(7))
            if (any(ierr(1:7) .ne. MPI_SUCCESS)) then
                write (6, *) 'MPI error: ', ierr(1:7)
                call shd_mpi__abort('CALC_BOUNDARY_RADIANCES: MPI error 3')
            end if
            do I = 1, Nbcpnts
                IPT = RadPtrIn(I)
                GRIDRAD(:, IPT) = RadAllIn(:, I)
            end do
        else
            do I = 1, Ndone(MYPROC_)
                IPT = RadPtrOut(I, MYPROC_)
                GRIDRAD(:, IPT) = RadAllOut(:, I, MYPROC_)
            end do
        end if
    end do

    deallocate (RadAllOut, RadPtrOut, RadAllIn, RadPtrIn, Ndone)
    deallocate (BndRadRealInfo, BndRadIntInfo)
end subroutine shd_mpi__calBoundaryRad

subroutine shd_mpi__integrateRay(BCFLAG, IPFLAG, NX, NY, NZ, NA,&
    &XGRID, YGRID, ZGRID, NPTS, NCELLS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &GRIDPOS, EXTINCT, NSTOKES, SOURCE, KANG,&
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

    integer, intent(IN) :: BCFLAG, IPFLAG, NX, NY, NZ, NA, NPTS, NCELLS
    integer, intent(IN) :: NSTOKES, KANG
    real, intent(IN) :: XGRID(NX), YGRID(NY), ZGRID(NZ)
    integer, intent(IN) :: GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS)
    integer, intent(IN) :: TREEPTR(2, NCELLS)
    integer*2, intent(IN) ::  CELLFLAGS(NCELLS)
    real, intent(IN) :: GRIDPOS(3, NPTS)
    real, intent(IN) :: EXTINCT(NPTS)
    real, intent(IN) :: SOURCE(NSTOKES, NA, NPTS), GRIDRAD(NSTOKES, NPTS)
    real, intent(IN) :: MU, PHI
    real, intent(INOUT) :: X, Y, Z
    integer, intent(INOUT) :: SIDE
    real, intent(INOUT) :: TRANSMIT, RADIANCE(NSTOKES)
    logical, intent(OUT) :: VALIDRAD

    integer :: BITX, BITY, BITZ, IOCT, NXC, NYC, BTEST
    integer :: IX, IY, IZ, IL, IM, IU, IPTR, DIR, ICELL, INEXTCELL
    integer :: I1, I2, I3, I4, IOPP, IFACE, JFACE, KFACE, IC
    integer :: GRIDFACE(4, 6), OPPFACE(6)
    logical :: IPINX, IPINY
    real(RD_) :: PI, CX, CY, CZ, CXINV, CYINV, CZINV
    real(RD_) :: XE, YE, ZE
    real(RD_) :: SO, SOX, SOY, SOZ, EPS
    real(RD_) :: U, V, F1, F2, F3, F4
    real(RD_) :: EXT, EXT0, EXT1
    real(RD_) :: SRC(NSTOKES), SRCEXT0(NSTOKES), SRCEXT1(NSTOKES)
    real(RD_) :: EXT0P, SRCEXT0P(NSTOKES)
    real(RD_) :: TAU, TRANSCELL, ABSCELL, TRANS, RAD(NSTOKES), RAD0(NSTOKES)
    data GRIDFACE/1, 3, 5, 7, 2, 4, 6, 8, 1, 2, 5, 6, 3, 4, 7, 8, 1, 2, 3, 4, 5, 6, 7, 8/
    data OPPFACE/2, 1, 4, 3, 6, 5/

    EPS = 1.0e-3*(GRIDPOS(3, GRIDPTR(8, 1)) - GRIDPOS(3, GRIDPTR(1, 1)))

    ! Make the ray direction (opposite to the discrete ordinate direction)
    PI = acos(-1.0)
    CX = sqrt(1.0 - MU**2)*cos(PHI + PI)
    CY = sqrt(1.0 - MU**2)*sin(PHI + PI)
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
        call shd_mpi__abort('Message [_mpiIntegrateRay]: Bad MU.')
    end if
    ! IOCT is the octant the ray come from or the discrete ordinate goes to.
    IOCT = 1 + BITX + 2*BITY + 4*BITZ

    ! Find the base grid cell that the starting point is in; first in Z
    IL = 0
    IU = NZ
    do while (IU - IL .gt. 1)
        IM = (IU + IL)/2
        if (Z .ge. ZGRID(IM)) then
            IL = IM
        else
            IU = IM
        end if
    end do
    IZ = max(IL, 1)

    ! Then get the base grid cell in X and Y assuming evenly spaced grid
    NXC = NX
    if (btest(BCFLAG, 2) .and. .not. btest(IPFLAG, 0)) NXC = NX - 1
    NYC = NY
    if (btest(BCFLAG, 3) .and. .not. btest(IPFLAG, 1)) NYC = NY - 1
    if (btest(IPFLAG, 0)) then
        IX = 1 + nint((NX - 1)*(X - XGRID(1))/(XGRID(NX) - XGRID(1)))
    else
        IX = 1 + int((NX - 1)*(X - XGRID(1))/(XGRID(NX) - XGRID(1)))
    end if
    if (btest(IPFLAG, 1)) then
        IY = 1 + nint((NY - 1)*(Y - YGRID(1))/(YGRID(NY) - YGRID(1)))
    else
        IY = 1 + int((NY - 1)*(Y - YGRID(1))/(YGRID(NY) - YGRID(1)))
    end if
    if (ix < 1 .or. iy < 1 .or. ix > NX + 1 .or. iy > NY + 1) then
        write (6, *) 'Bad ix,iy', x, y, z, ix, iy, iz, xgrid(1), xgrid(nx), ygrid(1), ygrid(ny), MYPROC_
        write (6, *) 'Bad ix,iy', side, transmit
        call shd_mpi__abort('Message [_mpiIntegrateRay]:')
    end if
    IX = min(NXC, IX)
    IY = min(NYC, IY)

    ! Make the base grid cell pointer
    ICELL = IZ + (NZ - 1)*(IY - 1) + (NZ - 1)*NYC*(IX - 1)

    ! Trace down the tree: point to the positive child cell, get the grid
    ! point that is on the negative (X,Y,Z) side, use the flags to find
    ! which dimension (X,Y,Z) the parent cell was split in, then compare
    ! the grid point with the test point (X,Y,Z) to find which child cell
    ! (pos or neg) the test point is in.
    do while (TREEPTR(2, ICELL) .gt. 0)
        DIR = ibits(int(CELLFLAGS(ICELL)), 2, 2)
        IC = TREEPTR(2, ICELL) + 1
        IPTR = GRIDPTR(1, IC)
        if (DIR .eq. 1) then
            if (X .lt. GRIDPOS(1, IPTR)) IC = IC - 1
        else if (DIR .eq. 2) then
            if (Y .lt. GRIDPOS(2, IPTR)) IC = IC - 1
        else if (DIR .eq. 3) then
            if (Z .lt. GRIDPOS(3, IPTR)) IC = IC - 1
        end if
        ICELL = IC
    end do

    IFACE = SIDE
    ! Get the four gridpoints for this face
    I1 = GRIDPTR(GRIDFACE(1, IFACE), ICELL)
    I2 = GRIDPTR(GRIDFACE(2, IFACE), ICELL)
    I3 = GRIDPTR(GRIDFACE(3, IFACE), ICELL)
    I4 = GRIDPTR(GRIDFACE(4, IFACE), ICELL)
    IPINX = btest(int(CELLFLAGS(ICELL)), 0)
    IPINY = btest(int(CELLFLAGS(ICELL)), 1)
    ! Compute the face interpolation factors
    if (IFACE .le. 2) then
        U = (Z - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
        if (IPINY) then
            V = 0.5
        else
            V = (Y - GRIDPOS(2, I1))/(GRIDPOS(2, I2) - GRIDPOS(2, I1))
        end if
    else
        U = (Z - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
        if (IPINX) then
            V = 0.5
        else
            V = (X - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
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
        IPINX = btest(int(CELLFLAGS(ICELL)), 0)
        IPINY = btest(int(CELLFLAGS(ICELL)), 1)

        ! Find boundaries of the current cell
        !   Find the three possible intersection planes (X,Y,Z)
        !   from the coordinates of the opposite corner grid point
        IOPP = GRIDPTR(9 - IOCT, ICELL)
        ! Get the distances to the 3 planes and select the closest
        if (IPINX) then
            SOX = 1.0e20
        else
            SOX = (GRIDPOS(1, IOPP) - XE)*CXINV
        end if
        if (IPINY) then
            SOY = 1.0e20
        else
            SOY = (GRIDPOS(2, IOPP) - YE)*CYINV
        end if
        SOZ = (GRIDPOS(3, IOPP) - ZE)*CZINV
        SO = min(SOX, SOY, SOZ)
        if (SO .lt. -EPS) then
            write (6, *) 'Message [_mpiIntegrateRay]: SO<0  ', MU, PHI, X, Y, Z, SOX, SOY, SOZ, SO, ICELL, TRANS
            call shd_mpi__abort('Message [_mpiIntegrateRay]: SO<0.')
        end if
        ! Compute the coordinates of the cell exitting location
        if (.not. IPINX) XE = XE + SO*CX
        if (.not. IPINY) YE = YE + SO*CY
        ZE = ZE + SO*CZ

        ! Get the intersection face number (i.e. neighptr index)
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
        INEXTCELL = NEIGHPTR(IFACE, ICELL)
        if (INEXTCELL .lt. 0) then
            call shd_grid__getNextCell(XE, YE, ZE, IFACE, JFACE, ICELL, GRIDPOS,&
                &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, INEXTCELL)
        end if

        ! Get the face grid pointers
        if (NEIGHPTR(IFACE, ICELL) .ge. 0) then
            ! If going to same or larger face then use previous face
            KFACE = IFACE
            IC = ICELL
        else
            ! If going to smaller face then use next face (more accurate)
            KFACE = OPPFACE(IFACE)
            IC = INEXTCELL
        end if
        I1 = GRIDPTR(GRIDFACE(1, KFACE), IC)
        I2 = GRIDPTR(GRIDFACE(2, KFACE), IC)
        I3 = GRIDPTR(GRIDFACE(3, KFACE), IC)
        I4 = GRIDPTR(GRIDFACE(4, KFACE), IC)
        ! Compute the face interpolation factors
        if (JFACE .eq. 1) then
            U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
            if (IPINY) then
                V = 0.5
            else
                V = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I2) - GRIDPOS(2, I1))
            end if
        else if (JFACE .eq. 2) then
            U = (ZE - GRIDPOS(3, I1))/(GRIDPOS(3, I3) - GRIDPOS(3, I1))
            if (IPINX) then
                V = 0.5
            else
                V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
            end if
        else
            if (IPINY) then
                U = 0.5
            else
                U = (YE - GRIDPOS(2, I1))/(GRIDPOS(2, I3) - GRIDPOS(2, I1))
            end if
            if (IPINX) then
                V = 0.5
            else
                V = (XE - GRIDPOS(1, I1))/(GRIDPOS(1, I2) - GRIDPOS(1, I1))
            end if
        end if
        if (u < -1.0e-4 .or. v < -1.0e-4 .or. u > 1.0001 .or. v > 1.0001) then
            write (6, *) 'u,v<0 or u,v>1: ', mu, phi, x, y, z, xe, ye, ze, u, v, ipinx, ipiny
            call shd_mpi__abort('Message [_mpiIntegrateRay]:')
        end if

        ! Get the location coordinate (does the boundary wrapping)
        if (INEXTCELL .gt. 0) then
            if (JFACE .eq. 1) then
                XE = GRIDPOS(1, GRIDPTR(IOCT, INEXTCELL))
            else if (JFACE .eq. 2) then
                YE = GRIDPOS(2, GRIDPTR(IOCT, INEXTCELL))
            else
                ZE = GRIDPOS(3, GRIDPTR(IOCT, INEXTCELL))
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
            ABSCELL = TAU*(1.0 - 0.5*TAU*(1.0 - 0.33333333333*TAU*(1 - 0.25*TAU)))
            TRANSCELL = 1.0 - ABSCELL
        end if
        if (TAU .le. 2.0) then
            if (EXT .eq. 0.0) then
                SRC(:) = 0.0
            else
                ! Linear extinction, linear source*extinction, to first order
                SRC(:) = (0.5*(SRCEXT0(:) + SRCEXT1(:))&
                    &+ 0.08333333333*(EXT0*SRCEXT1(:) - EXT1*SRCEXT0(:))*SO)/EXT
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
        if (RAD(1) < -1.0e-5) then
            write (6, '(3(1X,F5.3),8(1X,E12.5))')&
                &xe, ye, ze, ext0, ext1, so, tau, trans, src, abscell, rad(:)
        end if
        if (GRIDRAD(1, I1) .ge. 0.0 .and. GRIDRAD(1, I2) .ge. 0.0&
            .and. GRIDRAD(1, I3) .ge. 0.0 .and. GRIDRAD(1, I4) .ge. 0.0) then
            VALIDRAD = .true.
            RAD0(:) = F1*GRIDRAD(:, I1) + F2*GRIDRAD(:, I2)&
                &+ F3*GRIDRAD(:, I3) + F4*GRIDRAD(:, I4)
            RAD(:) = RAD(:) + TRANS*RAD0(:)
        else if (TRANS .lt. 1.0e-5) then
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

subroutine shd_mpi__computeRad(&
    &NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &BCFLAG, XDOMAIN, YDOMAIN, IPFLAG,&
    &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
    &SHPTR, SOURCE, SOURCE1, GRIDRAD,&
    &OUTPARMS, NRAD, RADOUT)

    ! Computes the radiances for the locations and directions specified
    ! in OUTPARMS, putting the results in RADOUT and increasing NRAD by
    ! the number of radiances output.  Does the radiative transfer integrations
    ! across multiple processors.

    integer, intent(IN) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer, intent(IN) :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer, intent(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(INOUT) :: NRAD
    integer, intent(IN) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
    integer, intent(IN) :: GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer, intent(IN) :: SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2, intent(IN) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical, intent(IN) :: DELTAM
    real, intent(IN) :: SOLARMU, SOLARAZ
    real, intent(IN) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real, intent(IN) :: MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real, intent(IN) :: XDOMAIN, YDOMAIN, XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real, intent(IN) :: GRIDPOS(3, NPTS)
    real, intent(IN) :: SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
    real, intent(IN) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real, intent(IN) :: DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(NSTOKES, *)
    real, intent(INOUT) :: SOURCE1(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    real, intent(IN) :: OUTPARMS(*)
    real, intent(OUT) :: RADOUT(NSTOKES, *)
    character, intent(IN) ::  SRCTYPE*1, SFCTYPE*2, UNITS*1

    integer :: IBC, JX, JY, K
    integer :: NANGOUT, NXOUT, NYOUT
    logical :: LAMBERTIAN, VALIDRAD
    real(RD_) :: X0, Y0, Z0, XE, YE, ZE, X, Y, Z, TRANSMIT, RADIANCE(NSTOKES)
    real    :: MUOUT, PHIOUT, PHID
    real    :: STARTX, STARTY, EPS = 1.0e-5

    integer :: MAXNRAD, NDOWN, Lupx, Lupy, Ldownx, Ldowny
    integer :: SIDE, NINVALID(0:4)
    integer :: I, J, L, Ndone, intsendbuf(2), intrecvbuf(2)
    integer :: ierr, irq, iproc
    logical :: AllDone
    integer, allocatable :: NradRecv(:)
    integer, allocatable :: IndexSend(:), IndexRecv(:, :), BndRadIntInfo(:, :)
    real, allocatable    :: RadSend(:, :), RadRecv(:, :, :), BndRadRealInfo(:, :, :)

    type(MPI_Request), allocatable :: requests(:)
    type(MPI_Status), allocatable :: status(:)

    ! added for fixing the "domain shift" issue
    !╭──────────────────────────────────────────────────────────────╮!
    real :: DS_Z0, DS_VZA_RAD, DS_XOFFSET, DS_YOFFSET, TOT_XOFFSET, TOT_YOFFSET
    !╰──────────────────────────────────────────────────────────────╯!

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces.
    ! Compute the upwelling bottom radiances using the downwelling fluxes.
    if (SFCTYPE .eq. 'FL') then
        call shd_init__makeSfcRadFL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
    else if (SFCTYPE .eq. 'VL') then
        call shd_init__makeSfcRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
            &NSTOKES, BCRAD(1, 1 + NTOPPTS))
    end if

    ! Setup the regularly spaced radiance locations
    STARTX = 0.0
    if (OUTPARMS(2) .le. 0.0) then
        NXOUT = 1
    else
        NXOUT = max(1, nint(XDOMAIN/OUTPARMS(2)))
    end if
    STARTY = 0.0
    if (OUTPARMS(3) .le. 0.0) then
        NYOUT = 1
    else
        NYOUT = max(1, nint(YDOMAIN/OUTPARMS(3)))
    end if
    Z0 = min(max(OUTPARMS(1), ZGRID(1)), ZGRID(NZ))
    NANGOUT = nint(OUTPARMS(6))

    LAMBERTIAN = SFCTYPE(2:2) .eq. 'L'

    MAXNRAD = nint(4.0*(NXOUT*NYOUT)/NUMPROC_)
    ! MAXNRAD = NXOUT*NYOUT
    allocate (RadSend(NSTOKES, MAXNRAD), IndexSend(MAXNRAD))
    if (MYPROC_ .eq. 0) allocate (NradRecv(NUMPROC_ - 1),&
        &RadRecv(NSTOKES, MAXNRAD, NUMPROC_ - 1), IndexRecv(MAXNRAD, NUMPROC_ - 1))
    allocate (BndRadRealInfo(4 + NSTOKES, MAXNRAD, 4), BndRadIntInfo(MAXNRAD, 4))
    allocate (requests(max(12, 3*NUMPROC_)), status(max(12, 3*NUMPROC_)))

    ! Loop over the radiance directions
    do K = 1, NANGOUT
        MUOUT = OUTPARMS(2*K + 5)
        PHID = OUTPARMS(2*K + 6)
        PHIOUT = PHID*acos(-1.0)/180.0

        ! added for fixing the "domain shift" issue (only for 3D)
        ! this is a first step correction:
        !   single-pixel-wise mini adjustment, adding
        !   interpolation effect
        ! second step correction is done in the data output
        !   <shd_io__radShiftCorr> in shdIO_out.f90
        !╭────────────────────────────────────────────────╮!
        if (IPFLAG .eq. 0) then
            DS_Z0 = Z0-ZGRID(1)
            DS_VZA_RAD = acos(MUOUT)
            DS_XOFFSET = DS_Z0*tan(DS_VZA_RAD)*cos(PHIOUT)
            DS_YOFFSET = DS_Z0*tan(DS_VZA_RAD)*sin(PHIOUT)
            TOT_XOFFSET = DS_XOFFSET
            TOT_YOFFSET = DS_YOFFSET
        else
            TOT_XOFFSET = 0.0
            TOT_YOFFSET = 0.0
        endif

        ! to avoid vertical stripes
        ! this extra pixel shift is considered in the data output (second step correction)
        STARTX = TOT_XOFFSET - OUTPARMS(2)*nint(TOT_XOFFSET/OUTPARMS(2))
        if (STARTX .lt. 0.0) then
            STARTX = STARTX + OUTPARMS(2)
        end if

        ! to avoid horizontal stripes
        ! this extra pixel shift is considered in the data output (second step correction)
        STARTY = TOT_YOFFSET - OUTPARMS(3)*nint(TOT_YOFFSET/OUTPARMS(3))
        if (STARTY .lt. 0.0) then
            STARTY = STARTY + OUTPARMS(3)
        end if
        !╰────────────────────────────────────────────────╯!

        if (MUOUT .eq. 0.0 .or. abs(MUOUT) .gt. 1.0) then
            write (6, *) 'Warning [COMPUTE_RADIANCE_PAR]: Bad mu for radiance', MUOUT
        else

            ! Compute the source function throughout grid for this angle
            call shd_src__computeSrcOneDir(NSTOKES, ML, MM, NLM, NSTLEG,&
                &NLEG, NUMPHASE, NPTS, DELTAM, MUOUT, PHIOUT,&
                &SRCTYPE, SOLARMU, SOLARAZ, ALBEDO, LEGEN, IPHASE,&
                &DIRFLUX, SHPTR, SOURCE, SOURCE1)

            ! Set the radiance field to -1, so we can determine valid radiances
            !   Also set the source array to the extinction times the source
            !   function, so the grid interpolation is correct.
            do I = 1, NPTS
                GRIDRAD(1, I) = -1.0
                SOURCE1(:, I) = SOURCE1(:, I)*EXTINCT(I)
            end do
            ! Get boundary radiances: either top or bottom
            ! Isotropic top boundary or Lambertian bottom boundary can use
            !   the previously computed boundary radiances in BCRAD,
            !   otherwise, compute the radiance for this angle
            !   by integrating over the stored downwelling radiances.
            if (MUOUT .lt. 0.0) then
                do IBC = 1, NTOPPTS
                    I = BCPTR(IBC, 1)
                    GRIDRAD(:, I) = BCRAD(:, IBC)
                end do
            else
                if (.not. LAMBERTIAN) then
                    call shd_init__makeSfcRadV(NBOTPTS, 1, NBOTPTS, BCPTR(1, 2),&
                        &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, MUOUT, PHIOUT,&
                        &SRCTYPE, WAVELEN, SOLARMU, SOLARAZ, DIRFLUX,&
                        &SFCTYPE, NSFCPAR, SFCGRIDPARMS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
                end if
                do IBC = 1, NBOTPTS
                    I = BCPTR(IBC, 2)
                    GRIDRAD(:, I) = BCRAD(:, NTOPPTS + IBC)
                end do
            end if

            ! Set up the upstream/downstream boundaries for this direction
            Lupx = 0; Lupy = 0; Ldownx = 0; Ldowny = 0
            if (btest(BCFLAG, 2)) then
                if ((cos(PHIOUT) .gt. 1.0e-5)) then
                    Ldownx = 2; Lupx = 1
                else
                    Ldownx = 1; Lupx = 2
                end if
            end if
            if (btest(BCFLAG, 3)) then
                !if ((sin(PHIOUT) .gt. 1.0e-5) .and. (NPY .gt. 1)) then
                if ((sin(PHIOUT) .gt. 1.0e-5)) then
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
            do JY = 1, NYOUT
                X0 = STARTX
                do JX = 1, NXOUT
                    NRAD = NRAD + 1

                    if (X0 .ge. (XGRID(1) - EPS) .and. X0 .le. (XGRID(NX) + EPS) .and.&
                        Y0 .ge. (YGRID(1) - EPS) .and. Y0 .le. (YGRID(NY) + EPS)) then
                        TRANSMIT = 1.0d0; RADIANCE(:) = 0.0d0
                        call shd_src__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
                            &NX, NY, NZ, NPTS, NCELLS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1,&
                            &X0, Y0, Z0, XE, YE, ZE, SIDE,&
                            &TRANSMIT, RADIANCE, VALIDRAD)

                        ! If we got a radiance then store it in RadSend otherwise save
                        !  the information needed to pass to other processors
                        if (VALIDRAD) then
                            Ndone = Ndone + 1
                            if (Ndone > MAXNRAD)&
                                &call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: MAXNRAD exceeded by Ndone 1.')
                            RadSend(:, Ndone) = RADIANCE(:)
                            IndexSend(Ndone) = NRAD
                        else
                            if ((SIDE .ne. Lupx .and. SIDE .ne. Lupy)) then
                                write (6, *) 'Bad SIDE 1:', MYPROC_, muout, phiout, x0, y0, z0, SIDE, xe, ye, ze
                                call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: Bad SIDE 1')
                            end if
                            NINVALID(SIDE) = NINVALID(SIDE) + 1
                            if (NINVALID(SIDE) > MAXNRAD)&
                                call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: MAXNRAD exceeded by NINVALID(SIDE) 1.')

                            BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, TRANSMIT, RADIANCE(:)/)
                            BndRadIntInfo(NINVALID(SIDE), SIDE) = NRAD
                        end if
                    end if
                    X0 = X0 + OUTPARMS(2)
                end do
                Y0 = Y0 + OUTPARMS(3)
            end do

            ! Cycle over the passing of info between processors and doing more
            !  partial radiance integrations until all rays have been finished
            AllDone = .false.
            do while (.not. AllDone)
                ! Do the send and receives of the number of incomplete rays and
                !   the real and integer info needed about the radiance integrations
                irq = 0
                if (btest(BCFLAG, 2)) then
                    call MPI_ISEND(NINVALID(Lupx), 1, MPI_INTEGER,&
                        &IPROCNEIGH_(Lupx), 6*K + 1, COMM2D_, requests(irq + 1), ierr)
                    call MPI_ISEND(BndRadRealInfo(:, :, Lupx), (4 + NSTOKES)*NINVALID(Lupx), MPI_DOUBLE_PRECISION,&
                        &IPROCNEIGH_(Lupx), 6*K + 2, COMM2D_, requests(irq + 2), ierr)
                    call MPI_ISEND(BndRadIntInfo(:, Lupx), NINVALID(Lupx), MPI_INTEGER,&
                        &IPROCNEIGH_(Lupx), 6*K + 3, COMM2D_, requests(irq + 3), ierr)
                    irq = irq + 3
                end if
                if (btest(BCFLAG, 3)) then
                    call MPI_ISEND(NINVALID(Lupy), 1, MPI_INTEGER,&
                        &IPROCNEIGH_(Lupy), 6*K + 4, COMM2D_, requests(irq + 1), ierr)
                    call MPI_ISEND(BndRadRealInfo(:, :, Lupy), (4 + NSTOKES)*NINVALID(Lupy), MPI_DOUBLE_PRECISION,&
                        &IPROCNEIGH_(Lupy), 6*K + 5, COMM2D_, requests(irq + 2), ierr)
                    call MPI_ISEND(BndRadIntInfo(:, Lupy), NINVALID(Lupy), MPI_INTEGER,&
                        &IPROCNEIGH_(Lupy), 6*K + 6, COMM2D_, requests(irq + 3), ierr)
                    irq = irq + 3
                end if
                if (btest(BCFLAG, 2)) then
                    call MPI_IRECV(NINVALID(Ldownx), 1, MPI_INTEGER,&
                        &IPROCNEIGH_(Ldownx), 6*K + 1, COMM2D_, requests(irq + 1), ierr)
                    call MPI_IRECV(BndRadRealInfo(:, :, Ldownx), (4 + NSTOKES)*MAXNRAD, MPI_DOUBLE_PRECISION,&
                        &IPROCNEIGH_(Ldownx), 6*K + 2, COMM2D_, requests(irq + 2), ierr)
                    call MPI_IRECV(BndRadIntInfo(:, Ldownx), MAXNRAD, MPI_INTEGER,&
                        &IPROCNEIGH_(Ldownx), 6*K + 3, COMM2D_, requests(irq + 3), ierr)
                    irq = irq + 3
                end if
                if (btest(BCFLAG, 3)) then
                    call MPI_IRECV(NINVALID(Ldowny), 1, MPI_INTEGER,&
                        &IPROCNEIGH_(Ldowny), 6*K + 4, COMM2D_, requests(irq + 1), ierr)
                    call MPI_IRECV(BndRadRealInfo(:, :, Ldowny), (4 + NSTOKES)*MAXNRAD, MPI_DOUBLE_PRECISION,&
                        &IPROCNEIGH_(Ldowny), 6*K + 5, COMM2D_, requests(irq + 2), ierr)
                    call MPI_IRECV(BndRadIntInfo(:, Ldowny), MAXNRAD, MPI_INTEGER,&
                        &IPROCNEIGH_(Ldowny), 6*K + 6, COMM2D_, requests(irq + 3), ierr)
                    irq = irq + 3
                end if
                call MPI_WAITALL(irq, requests, status, ierr)

                NINVALID(Lupx) = 0; NINVALID(Lupy) = 0
                if (btest(BCFLAG, 2)) then
                    ! Continue the backwards ray integrations for the X boundary
                    do J = 1, NINVALID(Ldownx)
                        if (Ldownx .eq. 1) BndRadRealInfo(1, J, Ldownx) = XGRID(1)
                        if (Ldownx .eq. 2) BndRadRealInfo(1, J, Ldownx) = XGRID(NX)
                        X = BndRadRealInfo(1, J, Ldownx); Y = BndRadRealInfo(2, J, Ldownx)
                        Z = BndRadRealInfo(3, J, Ldownx); TRANSMIT = BndRadRealInfo(4, J, Ldownx)
                        RADIANCE(:) = BndRadRealInfo(5:4 + NSTOKES, J, Ldownx)
                        call shd_src__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
                            &NX, NY, NZ, NPTS, NCELLS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1,&
                            &X, Y, Z, XE, YE, ZE, SIDE,&
                            &TRANSMIT, RADIANCE, VALIDRAD)
                        ! If we got a radiance then store it in RadSend, otherwise
                        ! save information needed to pass to neighboring processors
                        if (VALIDRAD) then
                            Ndone = Ndone + 1
                            if (Ndone > MAXNRAD)&
                                &call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: MAXNRAD exceeded by Ndone 2.')
                            RadSend(:, Ndone) = RADIANCE(:)
                            IndexSend(Ndone) = BndRadIntInfo(J, Ldownx)
                        else
                            if ((SIDE .ne. Lupx .and. SIDE .ne. Lupy)) then
                                write (6, *) 'Bad SIDE 2x', MYPROC_, muout, phiout, x, y, z, SIDE, xe, ye, ze
                                call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: Bad SIDE 2x')
                            end if
                            NINVALID(SIDE) = NINVALID(SIDE) + 1
                            if (NINVALID(SIDE) > MAXNRAD)&
                                &call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: MAXNRAD exceeded by NINVALID(SIDE) 2.')
                            BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, TRANSMIT, RADIANCE(:)/)
                            BndRadIntInfo(NINVALID(SIDE), SIDE) = BndRadIntInfo(J, Ldownx)
                        end if
                    end do
                end if

                if (btest(BCFLAG, 3)) then
                    ! Continue the backwards ray integrations for the Y boundary
                    do J = 1, NINVALID(Ldowny)
                        if (Ldowny .eq. 3) BndRadRealInfo(2, J, Ldowny) = YGRID(1)
                        if (Ldowny .eq. 4) BndRadRealInfo(2, J, Ldowny) = YGRID(NY)
                        X = BndRadRealInfo(1, J, Ldowny); Y = BndRadRealInfo(2, J, Ldowny)
                        Z = BndRadRealInfo(3, J, Ldowny); TRANSMIT = BndRadRealInfo(4, J, Ldowny)
                        RADIANCE(:) = BndRadRealInfo(5:4 + NSTOKES, J, Ldowny)
                        call shd_src__integrateSrc(BCFLAG, IPFLAG, NSTOKES,&
                            &NX, NY, NZ, NPTS, NCELLS,&
                            &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                            &XGRID, YGRID, ZGRID, GRIDPOS,&
                            &MUOUT, PHIOUT, GRIDRAD, EXTINCT, SOURCE1,&
                            &X, Y, Z, XE, YE, ZE, SIDE,&
                            &TRANSMIT, RADIANCE, VALIDRAD)
                        ! If we got a radiance then store it in RadSend, otherwise
                        ! save information needed to pass to neighboring processors
                        if (VALIDRAD) then
                            Ndone = Ndone + 1
                            if (Ndone > MAXNRAD)&
                                &call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: MAXNRAD exceeded by Ndone 3.')
                            RadSend(:, Ndone) = RADIANCE(:)
                            IndexSend(Ndone) = BndRadIntInfo(J, Ldowny)
                        else
                            if ((SIDE .ne. Lupx .and. SIDE .ne. Lupy)) then
                                write (6, *) 'Bad SIDE 2y', MYPROC_, muout, phiout, x, y, z, SIDE, xe, ye, ze
                                call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: Bad SIDE 2y')
                            end if
                            NINVALID(SIDE) = NINVALID(SIDE) + 1
                            if (NINVALID(SIDE) > MAXNRAD)&
                                &call shd_mpi__abort('Error [COMPUTE_RADIANCE_PAR]: MAXNRAD exceeded by NINVALID(SIDE) 3.')
                            BndRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, TRANSMIT, RADIANCE(:)/)
                            BndRadIntInfo(NINVALID(SIDE), SIDE) = BndRadIntInfo(J, Ldowny)
                        end if
                    end do
                end if

                ! See if there are any more invalid radiance rays on all the processors
                intsendbuf = (/NINVALID(Lupx), NINVALID(Lupy)/)
                call MPI_ALLREDUCE(intsendbuf, intrecvbuf, 2, MPI_INTEGER, MPI_MAX,&
                    &COMM2D_, ierr)
                AllDone = all(intrecvbuf(1:2) .eq. 0)
            end do

            if (MYPROC_ > 0) then
                ! Send the finished radiances to the master processor
                call MPI_ISEND(Ndone, 1, MPI_INTEGER, 0, 6*K + 7, COMM2D_, requests(1), ierr)
                call MPI_ISEND(RadSend, NSTOKES*Ndone, MPI_REAL, 0, 6*K + 8, COMM2D_, requests(2), ierr)
                call MPI_ISEND(IndexSend, Ndone, MPI_INTEGER, 0, 6*K + 9, COMM2D_, requests(3), ierr)
                call MPI_WAITALL(3, requests, status, ierr)
            else
                ! Put the radiances finished by the master processor in the list
                do i = 1, Ndone
                    RADOUT(:, IndexSend(i)) = RadSend(:, i)
                end do
                ! Receive the radiances from all the other processors
                irq = 0
                do iproc = 1, NUMPROC_ - 1
                    call MPI_IRECV(NradRecv(iproc), 1, MPI_INTEGER, iproc,&
                        &6*K + 7, COMM2D_, requests(irq + 1), ierr)
                    call MPI_IRECV(RadRecv(:, :, iproc), NSTOKES*MAXNRAD, MPI_REAL, iproc,&
                        &6*K + 8, COMM2D_, requests(irq + 2), ierr)
                    call MPI_IRECV(IndexRecv(:, iproc), MAXNRAD, MPI_INTEGER, iproc,&
                        &6*K + 9, COMM2D_, requests(irq + 3), ierr)
                    irq = irq + 3
                end do
                call MPI_WAITALL(irq, requests, status, ierr)
                do iproc = 1, NUMPROC_ - 1
                    do i = 1, NradRecv(iproc)
                        RADOUT(:, IndexRecv(i, iproc)) = RadRecv(:, i, iproc)
                    end do
                end do
            end if

        end if
    end do

    deallocate (requests, status, RadSend, IndexSend)
    if (MYPROC_ .eq. 0) deallocate (NradRecv, RadRecv, IndexRecv)
    deallocate (BndRadRealInfo, BndRadIntInfo)
end subroutine shd_mpi__computeRad

subroutine shd_mpi__visualizeRad(NSTOKES, NX, NY, NZ, NPTS, NCELLS,&
    &ML, MM, NLM, NSTLEG, NLEG, NUMPHASE,&
    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
    &BCFLAG, XDOMAIN, YDOMAIN, IPFLAG,&
    &SRCTYPE, DELTAM, SOLARMU, SOLARAZ,&
    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
    &GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS,&
    &XGRID, YGRID, ZGRID, GRIDPOS,&
    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
    &EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES,&
    &SHPTR, SOURCE, OUTPARMS, IVIS, VISOUT)

    ! Computes Stokes radiances (output in VISOUT) for the visualization
    ! modes: 1) camera mode, and 2) cross track scanning.

    integer, intent(IN) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer, intent(IN) :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer, intent(IN) :: NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(INOUT) :: IVIS
    integer, intent(IN) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
    integer, intent(IN) :: GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer, intent(IN) :: SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2, intent(IN) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical, intent(IN) :: DELTAM
    real, intent(IN) :: SOLARMU, SOLARAZ
    real, intent(IN) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real, intent(IN) :: MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real, intent(IN) :: XDOMAIN, YDOMAIN, XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real, intent(IN) :: GRIDPOS(3, NPTS)
    real, intent(IN) :: SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
    real, intent(IN) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real, intent(IN) :: DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(NSTOKES, *)
    real, intent(IN) :: OUTPARMS(*)
    real, intent(OUT) :: VISOUT(NSTOKES, *)
    character, intent(IN) :: SRCTYPE*1, SFCTYPE*2, UNITS*1

    integer :: NSCATANGLE, NSTPHASE, I, J, L, N
    integer :: NL, NS, LINE, SAMP
    logical :: CAMERA_MODE, VALIDRAD
    real, parameter  :: EPS = 1.0e-5
    real(RD_) :: XE, YE, ZE, TRANSMIT, VISRAD(NSTOKES), RADIANCE(NSTOKES)
    real(RD_) :: X0, Y0, Z0, X1, Y1, Z1, X2, Y2, Z2, X, Y, Z
    real(RD_) :: THETA0, THETA1, PHIR, PHI0, COSSCAT, MU2, PHI2
    real(RD_) :: COSTH, SINTH, SINTH0, SINTH1
    real(RD_) :: U, V, UP, VP, COSDELPHI, PI, ROTANG, DEGRAD, R
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

    type(MPI_Request), allocatable :: requests(:)
    type(MPI_Status), allocatable :: status(:)

    allocate (YLMSUN(NSTLEG, NLM))

    if (SRCTYPE .ne. 'T') then
        call shd_pha__computeYLM(.true., SOLARMU, SOLARAZ, ML, MM, NSTLEG, YLMSUN)
        if (DELTAM .and. NUMPHASE .gt. 0) then
            NSCATANGLE = max(36, min(MAXSCATANG, 2*NLEG))
            NSTPHASE = min(NSTLEG, 2)
            allocate (PHASETAB(NSTPHASE, NUMPHASE, NSCATANGLE))
            call shd_pha__computePhase(NSCATANGLE, NUMPHASE, NSTPHASE,&
                &NSTOKES, ML, NSTLEG, NLEG, LEGEN, PHASETAB)
        end if
    end if

    ! Make the isotropic radiances for the top boundary
    call shd_init__makeTopRad(SRCTYPE, SKYRAD, WAVENO, WAVELEN,&
        &UNITS, NTOPPTS, NSTOKES, BCRAD(1, 1))
    ! Make the bottom boundary radiances for the Lambertian surfaces.
    ! Compute the upwelling bottom radiances using the downwelling fluxes.
    if (SFCTYPE .eq. 'FL') then
        call shd_init__makeSfcRadFL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, GNDTEMP, GNDALBEDO,&
            &WAVENO, WAVELEN, UNITS, NSTOKES, BCRAD(1, 1 + NTOPPTS))
    else if (SFCTYPE .eq. 'VL') then
        call shd_init__makeSfcRadVL(NBOTPTS, BCPTR(1, 2),&
            &DIRFLUX, FLUXES, SRCTYPE, NSFCPAR, SFCGRIDPARMS,&
            &NSTOKES, BCRAD(1, 1 + NTOPPTS))
    end if

    CAMERA_MODE = nint(OUTPARMS(1)) .eq. 1
    if (CAMERA_MODE) then
        NL = nint(OUTPARMS(10))
        NS = nint(OUTPARMS(11))
    else
        NL = 1 + int(sqrt((OUTPARMS(4) - OUTPARMS(7))**2&
            &+ (OUTPARMS(5) - OUTPARMS(8))**2&
            &+ (OUTPARMS(6) - OUTPARMS(9))**2)/OUTPARMS(10))
        NS = 1 + int(abs(OUTPARMS(12) - OUTPARMS(11))/OUTPARMS(13))
    end if
    PI = acos(-1.0d0)
    DEGRAD = PI/180.

    ! Allocate sending, receiving, etc. arrays for MPI
    MAXNRAD = NL*NS
    allocate (RadSend(NSTOKES, MAXNRAD), IndexSend(MAXNRAD))
    if (MYPROC_ .eq. 0) allocate (NradRecv(NUMPROC_ - 1),&
        &RadRecv(NSTOKES, MAXNRAD, NUMPROC_ - 1), IndexRecv(MAXNRAD, NUMPROC_ - 1))
    allocate (BndSendRadRealInfo(6 + NSTOKES, MAXNRAD, 4), BndSendRadIntInfo(MAXNRAD, 4))
    allocate (BndRecvRadRealInfo(6 + NSTOKES, MAXNRAD, 4), BndRecvRadIntInfo(MAXNRAD, 4))
    allocate (requests(max(24, 3*NUMPROC_)), status(max(24, 3*NUMPROC_)))

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
                UP = (SAMP - NS/2 - 1)*OUTPARMS(13)*DEGRAD
                VP = (LINE - NL/2 - 1)*OUTPARMS(12)*DEGRAD
                ROTANG = OUTPARMS(9)*DEGRAD
                if (ROTANG .eq. 0.0) then
                    U = UP
                    V = VP
                else
                    U = cos(ROTANG)*UP - sin(ROTANG)*VP
                    V = sin(ROTANG)*UP + cos(ROTANG)*VP
                end if
                THETA0 = DEGRAD*OUTPARMS(7)
                PHI0 = DEGRAD*OUTPARMS(8)
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
                if (abs(SINTH) .lt. 1.0e-6) then
                    PHIR = 0.0
                else
                    if (abs(SINTH0) .lt. 1e-6 .and. abs(SINTH1) .le. 1e-6) then
                        COSDELPHI = 0.0d0
                    else if (abs(SINTH1) .gt. 1.0e-6) then
                        COSDELPHI = (cos(U) - COSTH*cos(THETA1))/(SINTH1*SINTH)
                    else if (abs(SINTH0) .gt. 1.0e-6) then
                        COSDELPHI = (cos(U)*cos(V) - COSTH*cos(THETA0))/(SINTH0*SINTH)
                    end if
                    COSDELPHI = min(+1.0d0, max(-1.0d0, COSDELPHI))
                    if (U .ge. 0.0) then
                        PHIR = PHI0 - acos(COSDELPHI)
                    else
                        PHIR = PHI0 + acos(COSDELPHI)
                    end if
                end if

                X0 = OUTPARMS(4)
                Y0 = OUTPARMS(5)
                Z0 = OUTPARMS(6)

            else

                ! Cross track scanning in the vertical plane:
                !  2, bytes, scale, X1,Y1,Z1, X2,Y2,Z2, spacing, scan1, scan2, delscan
                ! Start and end of scan angle (scan1, scan2) are +/- relative
                ! to nadir, and positive is on right side.
                X1 = OUTPARMS(4)
                Y1 = OUTPARMS(5)
                Z1 = OUTPARMS(6)
                X2 = OUTPARMS(7)
                Y2 = OUTPARMS(8)
                Z2 = OUTPARMS(9)
                DIST = sqrt((X2 - X1)**2 + (Y2 - Y1)**2 + (Z2 - Z1)**2)
                RX = (X2 - X1)/DIST
                RY = (Y2 - Y1)/DIST
                RZ = (Z2 - Z1)/DIST
                D = (NL - LINE)*OUTPARMS(10)
                X0 = X1 + D*RX
                Y0 = Y1 + D*RY
                Z0 = Z1 + D*RZ
                SCANANG = DEGRAD*(OUTPARMS(11) + (SAMP - 1)*OUTPARMS(13))
                COSTH = cos(PI - abs(SCANANG))
                if (SCANANG .gt. 0.0) then
                    PHIR = atan2(-RX, RY)
                else
                    PHIR = atan2(RX, -RY)
                end if
            end if

            ! Extrapolate ray to domain top if above
            if (Z0 .gt. ZGRID(NZ)) then
                if (COSTH .ge. 0.0) then
                    Ndone = Ndone + 1
                    RadSend(:, Ndone) = 0.0
                    IndexSend(Ndone) = IVIS
                    cycle
                end if
                R = (ZGRID(NZ) - Z0)/COSTH
                X0 = X0 + R*sqrt(1 - COSTH**2)*cos(PHIR)
                Y0 = Y0 + R*sqrt(1 - COSTH**2)*sin(PHIR)
                Z0 = ZGRID(NZ)
            else if (Z0 .lt. ZGRID(1)) then
                write (6, *) 'Error [VISUALIZE_RADIANCE]: Level below domain.'
                stop
            end if

            ! Relocate the starting point in the overall periodic domain
            ! (we know we have a periodic domain if doing MPI)
            if (X0 >= 0.0) then
                X0 = mod(X0, dble(XDOMAIN))
            else
                X0 = mod(X0 + XDOMAIN*int(1.0 - X0/XDOMAIN), dble(XDOMAIN))
            end if
            if (Y0 >= 0.0) then
                Y0 = mod(Y0, dble(YDOMAIN))
            else
                Y0 = mod(Y0 + YDOMAIN*int(1.0 - Y0/YDOMAIN), dble(YDOMAIN))
            end if

            ! Start the radiance calculation for those rays starting in this subdomain
            if (X0 .ge. XGRID(1) .and. X0 .le. XGRID(NX) .and.&
                Y0 .ge. YGRID(1) .and. Y0 .le. YGRID(NY)) then

                ! MU2,PHI2 is radiance travel direction (opposite the camera viewing direction)
                MU2 = -COSTH
                PHI2 = PHIR + PI

                ! Integrate the extinction and source function along this ray
                ! to calculate the Stokes radiance vector for this pixel
                TRANSMIT = 1.0d0; RADIANCE(:) = 0.0d0
                call shd_src__integrateSrcOneRay(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
                    &NSTPHASE, NSCATANGLE, PHASETAB,&
                    &NX, NY, NZ, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &XGRID, YGRID, ZGRID, GRIDPOS,&
                    &ML, MM, NLM, NLEG, NUMPHASE,&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                    &DELTAM, SRCTYPE, WAVELEN, SOLARMU, SOLARAZ,&
                    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
                    &DIRFLUX, SHPTR, SOURCE, YLMSUN,&
                    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                    &MU2, PHI2, X0, Y0, Z0,&
                    &XE, YE, ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

                ! If we got a radiance then store it in RadSend otherwise save
                !  the information needed to pass to other processors
                if (VALIDRAD) then
                    Ndone = Ndone + 1
                    if (Ndone > MAXNRAD)&
                        &call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 1')
                    RadSend(:, Ndone) = RADIANCE(:)
                    IndexSend(Ndone) = IVIS
                else
                    if (SIDE < 1 .or. SIDE > 4) call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: invalid SIDE 1')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) > MAXNRAD)&
                        &call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 1')
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
            if ((btest(BCFLAG, 2) .and. (SIDE .eq. 1 .or. SIDE .eq. 2)) .or.&
                (btest(BCFLAG, 3) .and. (SIDE .eq. 3 .or. SIDE .eq. 4))) then
                K = 3*SIDE
                call MPI_ISEND(NINVALID(SIDE), 1, MPI_INTEGER,&
                    &IPROCNEIGH_(SIDE), K + 1, COMM2D_, requests(irq + 1), ierr)
                call MPI_ISEND(BndSendRadRealInfo(:, :, SIDE), (6 + NSTOKES)*NINVALID(SIDE), MPI_DOUBLE_PRECISION,&
                    &IPROCNEIGH_(SIDE), K + 2, COMM2D_, requests(irq + 2), ierr)
                call MPI_ISEND(BndSendRadIntInfo(:, SIDE), NINVALID(SIDE), MPI_INTEGER,&
                    &IPROCNEIGH_(SIDE), K + 3, COMM2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
        end do
        ! Do the receives for all four sides of the subdomain.
        do SIDE = 1, 4
            if ((btest(BCFLAG, 2) .and. (SIDE .eq. 1 .or. SIDE .eq. 2)) .or.&
                (btest(BCFLAG, 3) .and. (SIDE .eq. 3 .or. SIDE .eq. 4))) then
                K = 3*RECVSIDE(SIDE)
                call MPI_IRECV(NINCOMING(SIDE), 1, MPI_INTEGER,&
                    &IPROCNEIGH_(SIDE), K + 1, COMM2D_, requests(irq + 1), ierr)
                call MPI_IRECV(BndRecvRadRealInfo(:, :, SIDE), (6 + NSTOKES)*MAXNRAD, MPI_DOUBLE_PRECISION,&
                    &IPROCNEIGH_(SIDE), K + 2, COMM2D_, requests(irq + 2), ierr)
                call MPI_IRECV(BndRecvRadIntInfo(:, SIDE), MAXNRAD, MPI_INTEGER,&
                    &IPROCNEIGH_(SIDE), K + 3, COMM2D_, requests(irq + 3), ierr)
                irq = irq + 3
            end if
        end do
        call MPI_WAITALL(irq, requests, status, ierr)

        ! Continue the backwards ray integrations for the four boundaries
        NINVALID(:) = 0
        do IBND = 1, 4
            do J = 1, NINCOMING(IBND)
                if (IBND .eq. 1) BndRecvRadRealInfo(1, J, IBND) = XGRID(1)
                if (IBND .eq. 2) BndRecvRadRealInfo(1, J, IBND) = XGRID(NX)
                if (IBND .eq. 3) BndRecvRadRealInfo(2, J, IBND) = YGRID(1)
                if (IBND .eq. 4) BndRecvRadRealInfo(2, J, IBND) = YGRID(NY)
                X = BndRecvRadRealInfo(1, J, IBND)
                Y = BndRecvRadRealInfo(2, J, IBND)
                Z = BndRecvRadRealInfo(3, J, IBND)
                MU2 = BndRecvRadRealInfo(4, J, IBND)
                PHI2 = BndRecvRadRealInfo(5, J, IBND)
                TRANSMIT = BndRecvRadRealInfo(6, J, IBND)
                RADIANCE(:) = BndRecvRadRealInfo(7:6 + NSTOKES, J, IBND)
                call shd_src__integrateSrcOneRay(BCFLAG, IPFLAG, NSTOKES, NSTLEG,&
                    &NSTPHASE, NSCATANGLE, PHASETAB,&
                    &NX, NY, NZ, NPTS, NCELLS,&
                    &GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS,&
                    &XGRID, YGRID, ZGRID, GRIDPOS,&
                    &ML, MM, NLM, NLEG, NUMPHASE,&
                    &NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO,&
                    &DELTAM, SRCTYPE, WAVELEN, SOLARMU, SOLARAZ,&
                    &EXTINCT, ALBEDO, LEGEN, IPHASE,&
                    &DIRFLUX, SHPTR, SOURCE, YLMSUN,&
                    &MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD,&
                    &SFCTYPE, NSFCPAR, SFCGRIDPARMS,&
                    &MU2, PHI2, X, Y, Z,&
                    &XE, YE, ZE, SIDE, TRANSMIT, RADIANCE, VALIDRAD)

                ! If we got a valid radiance then store it in RadSend, otherwise
                ! save information needed to pass to neighboring processors
                if (VALIDRAD) then
                    Ndone = Ndone + 1
                    if (Ndone > MAXNRAD)&
                        &call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by Ndone 2')
                    RadSend(:, Ndone) = RADIANCE(:)
                    IndexSend(Ndone) = BndRecvRadIntInfo(J, IBND)
                else
                    if (SIDE < 1 .or. SIDE > 4) call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: invalid SIDE 2')
                    NINVALID(SIDE) = NINVALID(SIDE) + 1
                    if (NINVALID(SIDE) > MAXNRAD)&
                        &call shd_mpi__abort('VISUALIZE_RADIANCE_PAR: MAXNRAD exceeded by NINVALID(SIDE) 2')
                    BndSendRadRealInfo(:, NINVALID(SIDE), SIDE) = (/XE, YE, ZE, MU2, PHI2, TRANSMIT, RADIANCE(:)/)
                    BndSendRadIntInfo(NINVALID(SIDE), SIDE) = BndRecvRadIntInfo(J, IBND)
                end if
            end do
        end do

        ! See if there are any more invalid radiance rays on all the processors
        call MPI_ALLREDUCE(NINVALID, intrecvbuf, 4, MPI_INTEGER, MPI_MAX,&
            &COMM2D_, ierr)
        AllDone = all(intrecvbuf(1:4) .eq. 0)
    end do
    !print *,'Done with ray tracing.', Ndone

    if (MYPROC_ > 0) then
        ! Send the finished radiances to the master processor
        call MPI_ISEND(Ndone, 1, MPI_INTEGER, 0, 21, COMM2D_, requests(1), ierr)
        call MPI_ISEND(RadSend, NSTOKES*Ndone, MPI_REAL, 0, 22, COMM2D_, requests(2), ierr)
        call MPI_ISEND(IndexSend, Ndone, MPI_INTEGER, 0, 23, COMM2D_, requests(3), ierr)
        call MPI_WAITALL(3, requests, status, ierr)
    else
        ! Put the radiances finished by the master processor in the list
        do i = 1, Ndone
            VISOUT(:, IndexSend(i)) = RadSend(:, i)
        end do
        ! Receive the radiances from all the other processors
        irq = 0
        do iproc = 1, NUMPROC_ - 1
            call MPI_IRECV(NradRecv(iproc), 1, MPI_INTEGER, iproc,&
                &21, COMM2D_, requests(irq + 1), ierr)
            call MPI_IRECV(RadRecv(:, :, iproc), NSTOKES*MAXNRAD, MPI_REAL, iproc,&
                &22, COMM2D_, requests(irq + 2), ierr)
            call MPI_IRECV(IndexRecv(:, iproc), MAXNRAD, MPI_INTEGER, iproc,&
                &23, COMM2D_, requests(irq + 3), ierr)
            irq = irq + 3
        end do
        call MPI_WAITALL(irq, requests, status, ierr)
        do iproc = 1, NUMPROC_ - 1
            do i = 1, NradRecv(iproc)
                VISOUT(:, IndexRecv(i, iproc)) = RadRecv(:, i, iproc)
            end do
        end do
    end if

    deallocate (requests, status, RadSend, IndexSend)
    if (MYPROC_ .eq. 0) deallocate (NradRecv, RadRecv, IndexRecv)
    deallocate (BndSendRadRealInfo, BndSendRadIntInfo)
    deallocate (BndRecvRadRealInfo, BndRecvRadIntInfo)
    deallocate (YLMSUN)
    if (allocated(PHASETAB)) deallocate (PHASETAB)

    ! On the master processor convert VISOUT(2:,*) from Stokes parameters
    ! to DOLP, AOLP, DOCP
    if (MYPROC_ .eq. 0) then
        if (NSTOKES .gt. 1) then
            N = NL*NS
            do I = IVIS - N + 1, IVIS
                VISRAD(:) = VISOUT(:, I)
                if (VISRAD(1) > 0.0) then
                    if (NSTOKES .gt. 1) then
                        ! Output degree (0 to 1) and angle (-180 to 180) of linear polarization
                        VISOUT(2, I) = sqrt(VISRAD(2)**2 + VISRAD(3)**2)/VISRAD(1)
                        VISOUT(3, I) = (180/PI)*0.5*atan2(VISRAD(3), VISRAD(2))
                    end if
                    if (NSTOKES .eq. 4) then
                        ! Output degree of circular polarization (-1 to 1)
                        VISOUT(4, I) = VISRAD(4)/VISRAD(1)
                    end if
                end if
            end do
            ! Choose the best range for the angle of linear polarization (-90 to 90 or 0 to 180)
            allocate (AOLP(N))
            AOLP(:) = VISOUT(3, IVIS - N + 1:IVIS)
            MEAN = sum(AOLP(:))/N
            STD1 = sqrt(sum((AOLP(:) - MEAN)**2)/N)
            where (AOLP(:) < 0.0)
                AOLP(:) = AOLP(:) + 180.0
            end where
            MEAN = sum(AOLP(:))/N
            STD2 = sqrt(sum((AOLP(:) - MEAN)**2)/N)
            if (STD2 < STD1) then
                VISOUT(3, IVIS - N + 1:IVIS) = AOLP(:)
            end if
            deallocate (AOLP)
        end if
    end if
end subroutine shd_mpi__visualizeRad

subroutine shd_mpi__calAccelSolCrit(DOACCEL, DELJDOT, DELJOLD, DELJNEW, JNORM,&
    &ACCELPAR, SOLCRIT)

    ! Does an MPI_ALLREDUCE to sum the delta source function vector dot
    ! products over all the processors, and from these calculates the
    ! full-domain-wide acceleration parameter and solution criterion.

    logical, intent(IN) :: DOACCEL
    real, intent(INOUT) :: DELJDOT, DELJOLD, DELJNEW, JNORM
    real, intent(OUT) :: ACCELPAR, SOLCRIT
    integer :: ierr
    real    ::  R, THETA, sendbuf(4), recvbuf(4)
    real, save :: A = 0.0

    if (NUMPROC_ > 1) then
        ! Call MPI_ALLREDUCE to sum up the four dot products
        sendbuf(:) = (/DELJDOT, DELJOLD, DELJNEW, JNORM/)
        call MPI_ALLREDUCE(sendbuf, recvbuf, 4, MPI_REAL, MPI_SUM, COMM2D_, ierr)
        if (ierr .ne. MPI_SUCCESS) call shd_mpi__abort('CALC_ACCEL_SOLCRIT MPI_ALLREDUCE error')
        DELJDOT = recvbuf(1); DELJOLD = recvbuf(2)
        DELJNEW = recvbuf(3); JNORM = recvbuf(4)
    end if

    ! Accelerate if desired, didn't last time, and things are converging.
    if (DOACCEL .and. A .eq. 0.0 .and. DELJNEW .lt. DELJOLD) then
        ! Compute the acceleration extrapolation factor and apply it.
        R = sqrt(DELJNEW/DELJOLD)
        THETA = acos(DELJDOT/sqrt(DELJOLD*DELJNEW))
        A = (1 - R*cos(THETA) + R**(1 + 0.5*3.14159/THETA))&
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
end subroutine shd_mpi__calAccelSolCrit

subroutine shd_mpi__end(NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART, XDOMAIN, YDOMAIN,&
    &NPXT, NPYT, PROPFILE, RUNNAME)

    ! Shuts down MPI at the end, but mainly finds the number of gridpoints
    ! around each property grid point column, sends that array to the master
    ! process, and the master process write out the RUNNAME//"_load_balance.out"
    ! file.

    integer, intent(in) :: NPTS, NPX, NPY, NPXT, NPYT
    real, intent(in) :: GRIDPOS(3, NPTS), XSTART, YSTART, XDOMAIN, YDOMAIN
    character(len=*), intent(in) :: PROPFILE, RUNNAME
    integer :: I, IX, IY, MAXN, SX, SY, EX, EY, NX, NY
    integer :: iproc, ierr, irq
    integer, allocatable :: NUMGRID(:, :), NUMGRIDALL(:, :), NORM(:, :), IXYALL(:, :, :)
    integer, allocatable :: recvbuf(:, :), buf(:, :)
    real :: DELX, DELY
    real :: ElapsedTime
    type(MPI_Request), allocatable :: requests(:)
    type(MPI_Status), allocatable :: status(:)

    DELX = real(XDOMAIN)/real(NPXT)
    DELY = real(YDOMAIN)/real(NPYT)

    if (LOADBALANCE_) then
        ! Calculate the number of grid points around each property grid column
        allocate (NUMGRID(NPX, NPY))
        NUMGRID(:, :) = 0
        do I = 1, NPTS
            IX = nint((GRIDPOS(1, I) - XSTART)/DELX) + 1
            IY = nint((GRIDPOS(2, I) - YSTART)/DELY) + 1
            if (IX <= NPX .and. IY <= NPY) then
                NUMGRID(IX, IY) = NUMGRID(IX, IY) + 1
            end if
        end do

        if (NUMPROC_ > 1) then
            ! Have the master process collect the property grid locations of all
            allocate (IXYALL(2, 2, 0:NUMPROC_ - 1))
            call MPI_GATHER(IXYPRP_, 4, MPI_INTEGER, IXYALL, 4, MPI_INTEGER, 0, COMM2D_, ierr)

            ! Send the NUMGRIDs to the master processor
            allocate (requests(NUMPROC_), status(NUMPROC_))
            irq = 0
            if (MYPROC_ > 0) then
                irq = irq + 1
                call MPI_ISEND(NUMGRID, NPX*NPY, MPI_INTEGER, 0, 777, COMM2D_,&
                    &requests(irq), ierr)
            else
                MAXN = (maxval(IXYALL(2, 1, :) - IXYALL(1, 1, :)) + 1)*(maxval(IXYALL(2, 2, :) - IXYALL(1, 2, :)) + 1)
                allocate (recvbuf(MAXN, 1:NUMPROC_ - 1))
                do iproc = 1, NUMPROC_ - 1
                    irq = irq + 1
                    call MPI_IRECV(recvbuf(:, iproc), MAXN, MPI_INTEGER, iproc, 777,&
                        &COMM2D_, requests(irq), ierr)
                end do
            end if
            call MPI_WAITALL(irq, requests, status, ierr)

            ! Put the NUMGRIDS into the full domain array on the master processor.
            !   Use a normalization counter (NORM) to deal with grid points that
            !   get contribution from two processors.
            if (MYPROC_ .eq. 0) then
                allocate (NUMGRIDALL(NPXT, NPYT), NORM(NPXT, NPYT))
                NUMGRIDALL(:, :) = 0; NORM(:, :) = 0
                SX = IXYALL(1, 1, 0); EX = IXYALL(2, 1, 0); NX = EX - SX + 1
                SY = IXYALL(1, 2, 0); EY = IXYALL(2, 2, 0); NY = EY - SY + 1
                NUMGRIDALL(SX:EX, SY:EY) = NUMGRID(:, :)
                NORM(SX:EX, SY:EY) = reshape((/1/), (/NX, NY/), (/1/))
                do iproc = 1, NUMPROC_ - 1
                    SX = IXYALL(1, 1, iproc); EX = IXYALL(2, 1, iproc); NX = EX - SX + 1
                    SY = IXYALL(1, 2, iproc); EY = IXYALL(2, 2, iproc); NY = EY - SY + 1
                    allocate (buf(NX, NY))
                    buf(:, :) = reshape(recvbuf(1:NX*NY, iproc), (/NX, NY/))
                    NUMGRIDALL(SX:EX - 1, SY:EY - 1) = NUMGRIDALL(SX:EX - 1, SY:EY - 1) + buf(1:NX - 1, 1:NY - 1)
                    NORM(SX:EX - 1, SY:EY - 1) = NORM(SX:EX - 1, SY:EY - 1) + reshape((/1/), (/NX - 1, NY - 1/), (/1/))
                    IX = mod(EX - 1, NPXT) + 1; IY = mod(EY - 1, NPYT) + 1
                    NUMGRIDALL(IX, SY:EY - 1) = NUMGRIDALL(IX, SY:EY - 1) + buf(NX, 1:NY - 1)
                    NORM(IX, SY:EY - 1) = NORM(IX, SY:EY - 1) + reshape((/1/), (/NY - 1/), (/1/))
                    NUMGRIDALL(SX:EX - 1, IY) = NUMGRIDALL(SX:EX - 1, IY) + buf(1:NX - 1, NY)
                    NORM(SX:EX - 1, IY) = NORM(SX:EX - 1, IY) + reshape((/1/), (/NX - 1/), (/1/))
                    NUMGRIDALL(IX, IY) = NUMGRIDALL(IX, IY) + buf(NX, NY)
                    NORM(IX, IY) = NORM(IX, IY) + 1
                    deallocate (buf)
                end do
                deallocate (recvbuf)
                NUMGRIDALL(:, :) = NUMGRIDALL(:, :)/NORM(:, :)
                deallocate (NORM)
            end if
            deallocate (requests, status, IXYALL)

        else
            allocate (NUMGRIDALL(NPXT, NPYT))
            NUMGRIDALL(1:NPX, 1:NPY) = NUMGRID(:, :)
        end if
        deallocate (NUMGRID)

        if (MYPROC_ .eq. 0) then
            ! Write out the array to the load balancing file
            open (UNIT=10, FILE=trim(RUNNAME)//'_load_balance.out', STATUS='UNKNOWN')
            write (10, '(A)') '! SHDOM grid point load balancing file for property file:'
            write (10, '(A)') PROPFILE
            write (10, '(2(1X,I4),A)') NPXT, NPYT, ' ! property grid size (Nx,Ny)'
            do IX = 1, NPXT
                do IY = 1, NPYT
                    write (10, '(I4,1X,I4,1X,I5)') IX, IY, NUMGRIDALL(IX, IY)
                end do
            end do
            close (10)
            deallocate (NUMGRIDALL)
        end if
    end if

    ElapsedTime = MPI_WTIME() - STARTTIME_

    if (ElapsedTime .lt. 60.0) then
        write (6, '(A,F6.2)') 'Message [_mpiEnd]: Total elapsed time (sec):', ElapsedTime
    else
        write (6, '(A,F6.2)') 'Message [_mpiEnd]: Total elapsed time (min):', ElapsedTime/60.0
    end if

    write (6, '(A,"<",A,">")') 'Message [_mpiEnd]: Finished with SHDOM run: ', trim(adjustl(RUNNAME))

    ! Shut down MPI
    if (allocated(BXYPTR_)) deallocate (BXYPTR_, IBXYPTR_, BXYPOS_)
    call MPI_FINALIZE(ierr)

end subroutine shd_mpi__end

subroutine shd_mpi__abort(ERRSTR)

    character(len=*), intent(in) :: ERRSTR
    integer :: ierr, errorcode = 7

    write (6, *) trim(adjustl(ERRSTR))
    call MPI_ABORT(COMM2D_, errorcode, ierr)

end subroutine shd_mpi__abort

end module shd_mpi
