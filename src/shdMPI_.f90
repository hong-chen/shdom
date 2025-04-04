module shd_mpi

    implicit none
    private

    public :: shd_mpi__start
    public :: shd_mpi__map
    public :: shd_mpi__bcastInput
    public :: shd_mpi__bcastPrpSize
    public :: shd_mpi__scatPrp
    public :: shd_mpi__readBcastMemParam
    public :: shd_mpi__bcastSfcSize
    public :: shd_mpi__bcastSfcParam
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
    ! Dummy routines to allow parallel SHDOM to compile without MPI.
    logical, intent(out) :: MASTERPROC

    MASTERPROC = .true.
    return
end subroutine shd_mpi__start

subroutine shd_mpi__map(&
    &BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY,&
    &PROPFILE, XSTART, YSTART, RUNNAME)

    integer BCFLAG, NPX, NPY, NX, NY, NZ
    character(len=*), intent(in) :: PROPFILE, RUNNAME
    real, intent(in) :: DELX, DELY
    real, intent(out) :: XSTART, YSTART

    XSTART = 0.0
    YSTART = 0.0
    return
end subroutine shd_mpi__map

subroutine shd_mpi__bcastInput(RUNNAME,&
    &PROPFILE, SFCFILE, CKDFILE,&
    &INSAVEFILE, OUTSAVEFILE,&
    &NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG,&
    &IPFLAG, KDIST, DELTAM, GRIDTYPE,&
    &SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD,&
    &GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN,&
    &ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC,&
    &MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS,&
    &MAX_TOTAL_MB, SPLITTING_FACTOR,&
    &NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO)

    integer NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    integer MAXOUT, MAXPAR, MAXITER, NUMOUT
    logical KDIST, DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR, MAXOUT)
    character SRCTYPE*1, GRIDTYPE*1, UNITS*1, OUTTYPES(*)*1
    character RUNNAME*256, PROPFILE*256, SFCFILE*256, CKDFILE*256
    character INSAVEFILE*256, OUTSAVEFILE*256
    real :: MAX_TOTAL_MB, SPLITTING_FACTOR
    real :: NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO

    return
end subroutine shd_mpi__bcastInput

subroutine shd_mpi__bcastPrpSize(NPX, NPY, NPZ, DELX, DELY,&
    &NUMPHASE, MAXLEG, MAXPGL)

    integer NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
    real DELX, DELY

    return
end subroutine shd_mpi__bcastPrpSize

subroutine shd_mpi__scatPrp(NPXT, NPYT, NPX, NPY, NPZ, NSTLEG,&
    &NLEG, NUMPHASE, ZLEVELS, MAXASYM, TEMPPT,&
    &EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT,&
    &TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)

    integer NPXT, NPYT, NPX, NPY, NPZ, NUMPHASE, NLEG, NSTLEG
    real ZLEVELS(NPZ), MAXASYM
    real TEMPPT(NPZ, NPYT, NPXT), EXTINCTPT(NPZ, NPYT, NPXT)
    real ALBEDOPT(NPZ, NPYT, NPXT), LEGENPT(*)
    integer*2 IPHASEPT(NPZ, NPYT, NPXT)
    real TEMPP(NPZ, NPY, NPX), EXTINCTP(NPZ, NPY, NPX)
    real ALBEDOP(NPZ, NPY, NPX), LEGENP(*)
    integer*2 IPHASEP(NPZ, NPY, NPX)

    return
end subroutine shd_mpi__scatPrp

subroutine shd_mpi__readBcastMemParam(MAX_TOTAL_MB,&
    &ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR,&
    &CELL_TO_POINT_RATIO, RUNNAME)

    real MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(*) RUNNAME

    return
end subroutine shd_mpi__readBcastMemParam

subroutine shd_mpi__bcastSfcSize(MAXSFCPTS, MAXSFCPARS)

    integer MAXSFCPTS, MAXSFCPARS

    return
end subroutine shd_mpi__bcastSfcSize

subroutine shd_mpi__bcastSfcParam(SFCTYPE, NXSFC, NYSFC, NSFCPAR,&
    &DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)

    character SFCTYPE*2
    integer NXSFC, NYSFC, NSFCPAR
    real DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO

    return
end subroutine shd_mpi__bcastSfcParam

subroutine shd_mpi__bcastKdistSize(NG, NZCKD)

    integer NG, NZCKD

    return
end subroutine shd_mpi__bcastKdistSize

subroutine shd_mpi__bcastKdistParam(SOLFLUX, NG, DELG,&
    &NZCKD, ZCKD, KABS)

    integer NG, NZCKD
    real SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG)

    return
end subroutine shd_mpi__bcastKdistParam

subroutine shd_mpi__gatherOutput(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
    &IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT,&
    &SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS,&
    &XDOMAIN, YDOMAIN, XALLGRID, YALLGRID,&
    &ALLFLUXES, ALLFLUXDIV, ALLSHTERMS,&
    &NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)

    logical FLUX_OUT, FLUXDIV_OUT, SH_OUT
    integer IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPX, NPY
    real, intent(in) :: XDOMAIN, YDOMAIN
    real SUMFLUXES(2, NBPTS), SUMDIRFLUX(NBPTS)
    real SUMFLUXDIV(NBPTS), SUMSHTERMS(NSHOUT, NBPTS)
    real XALLGRID(NXT), YALLGRID(NYT)
    real ALLFLUXES(3, NZ, NYT, NXT), ALLFLUXDIV(NZ, NYT, NXT)
    real ALLSHTERMS(NSHOUT, NZ, NYT, NXT)
    integer NCELLS, NPTS, NSH
    integer NCELLSTOT, NPTSTOT, NSHTOT

    return
end subroutine shd_mpi__gatherOutput

real function shd_mpi__sumCpuTime(cpuTimes)

    real cpuTimes
    shd_mpi__sumCpuTime = cpuTimes
    return
end function shd_mpi__sumCpuTime

subroutine shd_mpi__findAlbMax(ALBMAX)

    real ALBMAX
    return
end subroutine shd_mpi__findAlbMax

subroutine shd_mpi__unifySplit(DOSPLIT, STARTSPLITACC)

    logical DOSPLIT
    real STARTSPLITACC
    return
end subroutine shd_mpi__unifySplit

subroutine shd_mpi__findTotSplitCritMax(SPLITCRIT)

    real SPLITCRIT
    return
end subroutine shd_mpi__findTotSplitCritMax

subroutine shd_mpi__makeDirect(&
    &SPT, NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
    &NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
    &NX, XGRID, NY, YGRID, DIRFLUX)

    integer SPT, NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NX, NY
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), XGRID(NX + 1), YGRID(NY + 1)
    real DIRFLUX(NPTS)

    return
end subroutine shd_mpi__makeDirect

subroutine shd_mpi__findBoundaryPnt(&
    &BCFLAG, IPFLAG, NPTS, SWEEPORD,&
    &GRIDPTR, GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)

    integer BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS, *), GRIDPTR(8, *)
    integer NX, NY, NZ
    real GRIDPOS(3, NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)

    return
end subroutine shd_mpi__findBoundaryPnt

subroutine shd_mpi__calBoundaryRad(&
    &BCFLAG, IPFLAG, JOCT, IZ,&
    &NX, NY, NZ, XGRID, YGRID, ZGRID,&
    &NA, NPTS, NCELLS, GRIDPTR,&
    &NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS,&
    &MU, PHI, EXTINCT, NSTOKES, SOURCE,&
    &KANG, GRIDRAD)

    integer BCFLAG, IPFLAG, JOCT, IZ
    integer NX, NY, NZ, NA, NPTS, NCELLS, NSTOKES, KANG
    real XGRID(NX), YGRID(NY), ZGRID(NZ)
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS), MU, PHI
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)

    return
end subroutine shd_mpi__calBoundaryRad

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

    integer, intent(in) :: NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer, intent(in) :: ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer, intent(in) :: NMU, NPHI0MAX, NPHI0(NMU)
    integer, intent(inout) :: NRAD
    integer, intent(in) :: MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
    integer, intent(in) :: GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer, intent(in) :: SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2, intent(in) :: CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical, intent(in) :: DELTAM
    real, intent(in) :: SOLARMU, SOLARAZ
    real, intent(in) :: GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real, intent(in) :: MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real, intent(in) :: XDOMAIN, YDOMAIN, XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real, intent(in) :: GRIDPOS(3, NPTS)
    real, intent(in) :: SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
    real, intent(in) :: EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real, intent(in) :: DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(NSTOKES, *)
    real, intent(inout) :: SOURCE1(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    real, intent(in) :: OUTPARMS(*)
    real, intent(out) :: RADOUT(NSTOKES, *)
    character, intent(in) ::  SRCTYPE*1, SFCTYPE*2, UNITS*1

    return
end subroutine shd_mpi__computeRad

subroutine shd_mpi__visualizeRad(&
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
    &SHPTR, SOURCE, OUTPARMS, IVIS, VISOUT)

    integer NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer ML, MM, NLM, NSTLEG, NLEG, NUMPHASE
    integer NMU, NPHI0MAX, NPHI0(NMU), IVIS
    integer MAXNBC, NTOPPTS, NBOTPTS, NSFCPAR
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer SHPTR(NPTS + 1), BCPTR(MAXNBC, 2)
    integer*2 CELLFLAGS(NCELLS), IPHASE(NPTS)
    logical DELTAM
    real SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real MU(NMU), PHI(NMU, NPHI0MAX), WTDO(NMU, NPHI0MAX)
    real XDOMAIN, YDOMAIN, XGRID(NX + 1), YGRID(NY + 1), ZGRID(NZ)
    real GRIDPOS(3, NPTS)
    real SFCGRIDPARMS(*), BCRAD(NSTOKES, *)
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, *)
    real DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(NSTOKES, *)
    real OUTPARMS(*), VISOUT(NSTOKES, *)
    character SRCTYPE*1, SFCTYPE*2, UNITS*1

    return
end subroutine shd_mpi__visualizeRad

subroutine shd_mpi__calAccelSolCrit(&
    &DOACCEL, DELJDOT, DELJOLD, DELJNEW,&
    &JNORM, ACCELPAR, SOLCRIT)
    ! Calculates the acceleration parameter and solution criterion from
    ! the delta source function vector dot  products.

    logical DOACCEL
    real DELJDOT, DELJOLD, DELJNEW, JNORM, ACCELPAR, SOLCRIT
    real R, THETA, A
    save A
    data A/0.0/

    ! Accelerate if desired, didn't last time, and things are convergi
    if (DOACCEL .and. A .eq. 0.0 .and. DELJNEW .lt. DELJOLD) then
        ! Compute the acceleration extrapolation factor and apply it.
        R = sqrt(DELJNEW/DELJOLD)
        THETA = acos(DELJDOT/sqrt(DELJOLD*DELJNEW))
        A = (1 - R*cos(THETA) + R**(1 + 0.5*3.14159/THETA))&
            &/(1 + R**2 - 2*R*cos(THETA)) - 1.0
        A = min(10.0, max(0.0, A))
        ! WRITE (*,'(1X,A,3(1X,F7.3))') '! Acceleration: ', A,R,THETA
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

subroutine shd_mpi__end(NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART,&
    &XDOMAIN, YDOMAIN, NPXT, NPYT, PROPFILE, RUNNAME)

    integer NPTS, NPX, NPY, NPXT, NPYT
    real GRIDPOS(3, NPTS), XSTART, YSTART, XDOMAIN, YDOMAIN
    character(len=*), intent(in) :: PROPFILE, RUNNAME

    return
end subroutine shd_mpi__end

subroutine shd_mpi__abort(ERRSTR)

    character(*) ERRSTR

    write (6, *) ERRSTR
    stop
end subroutine shd_mpi__abort

end module shd_mpi
