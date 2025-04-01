module shd_ext_mpi

    implicit none
    private

    public :: shd_ext__mpiStart
    public :: shd_ext__mpiMap
    public :: shd_ext__mpiBcastInput
    public :: shd_ext__mpiBcastPropSize
    public :: shd_ext__mpiScatProp
    public :: shd_ext__mpiReadBcastMemParam
    public :: shd_ext__mpiBcastSurfSize
    public :: shd_ext__mpiBcastSurfParam
    public :: shd_ext__mpiBcastKdistSize
    public :: shd_ext__mpiBcastKdistParam
    public :: shd_ext__mpiGatherOutput
    public :: shd_ext__mpiMakeDirect
    public :: shd_ext__mpiComputeRad
    public :: shd_ext__mpiVisualizeRad
    public :: shd_ext__mpiEnd
    public :: shd_ext__mpiAbort

    public :: shd_ext__mpiFindAlbMax
    public :: shd_ext__mpiFindTotSplitCritMax
    public :: shd_ext__mpiFindBoundaryPnt
    public :: shd_ext__mpiCalBoundaryRad
    public :: shd_ext__mpiCalAccelSolCrit
    public :: shd_ext__mpiUnifySplit
    public :: shd_ext__mpiSumCpuTime

contains

subroutine shd_ext__mpiStart(MASTERPROC)
    ! Dummy routines to allow parallel SHDOM to compile without MPI.
    logical, intent(out) :: MASTERPROC

    MASTERPROC = .true.
    return
end subroutine shd_ext__mpiStart

subroutine shd_ext__mpiMap(&
    &BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY,&
    &PROPFILE, XSTART, YSTART, RUNNAME)

    integer BCFLAG, NPX, NPY, NX, NY, NZ
    character(len=*), intent(in) :: PROPFILE, RUNNAME
    real, intent(in) :: DELX, DELY
    real, intent(out) :: XSTART, YSTART

    XSTART = 0.0
    YSTART = 0.0
    return
end subroutine shd_ext__mpiMap

subroutine shd_ext__mpiBcastInput(RUNNAME,&
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
end subroutine shd_ext__mpiBcastInput

subroutine shd_ext__mpiBcastPropSize(NPX, NPY, NPZ, DELX, DELY,&
    &NUMPHASE, MAXLEG, MAXPGL)

    integer NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
    real DELX, DELY

    return
end subroutine shd_ext__mpiBcastPropSize

subroutine shd_ext__mpiScatProp(NPXT, NPYT, NPX, NPY, NPZ, NSTLEG,&
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
end subroutine shd_ext__mpiScatProp

subroutine shd_ext__mpiReadBcastMemParam(MAX_TOTAL_MB,&
    &ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR,&
    &CELL_TO_POINT_RATIO, RUNNAME)

    real MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(*) RUNNAME

    return
end subroutine shd_ext__mpiReadBcastMemParam

subroutine shd_ext__mpiBcastSurfSize(MAXSFCPTS, MAXSFCPARS)

    integer MAXSFCPTS, MAXSFCPARS

    return
end subroutine shd_ext__mpiBcastSurfSize

subroutine shd_ext__mpiBcastSurfParam(SFCTYPE, NXSFC, NYSFC, NSFCPAR,&
    &DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)

    character SFCTYPE*2
    integer NXSFC, NYSFC, NSFCPAR
    real DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO

    return
end subroutine shd_ext__mpiBcastSurfParam

subroutine shd_ext__mpiBcastKdistSize(NG, NZCKD)

    integer NG, NZCKD

    return
end subroutine shd_ext__mpiBcastKdistSize

subroutine shd_ext__mpiBcastKdistParam(SOLFLUX, NG, DELG,&
    &NZCKD, ZCKD, KABS)

    integer NG, NZCKD
    real SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG)

    return
end subroutine shd_ext__mpiBcastKdistParam

subroutine shd_ext__mpiGatherOutput(FLUX_OUT, FLUXDIV_OUT, SH_OUT,&
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
end subroutine shd_ext__mpiGatherOutput

real function shd_ext__mpiSumCpuTime(cpuTimes)

    real cpuTimes
    shd_ext__mpiSumCpuTime = cpuTimes
    return
end function shd_ext__mpiSumCpuTime

subroutine shd_ext__mpiFindAlbMax(ALBMAX)

    real ALBMAX
    return
end subroutine shd_ext__mpiFindAlbMax

subroutine shd_ext__mpiUnifySplit(DOSPLIT, STARTSPLITACC)

    logical DOSPLIT
    real STARTSPLITACC
    return
end subroutine shd_ext__mpiUnifySplit

subroutine shd_ext__mpiFindTotSplitCritMax(SPLITCRIT)

    real SPLITCRIT
    return
end subroutine shd_ext__mpiFindTotSplitCritMax

subroutine shd_ext__mpiMakeDirect(&
    &SPT, NPTS, BCFLAG, IPFLAG, DELTAM, ML,&
    &NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS,&
    &NX, XGRID, NY, YGRID, DIRFLUX)

    integer SPT, NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NX, NY
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), XGRID(NX + 1), YGRID(NY + 1)
    real DIRFLUX(NPTS)

    return
end subroutine shd_ext__mpiMakeDirect

subroutine shd_ext__mpiFindBoundaryPnt(&
    &BCFLAG, IPFLAG, NPTS, SWEEPORD,&
    &GRIDPTR, GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)

    integer BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS, *), GRIDPTR(8, *)
    integer NX, NY, NZ
    real GRIDPOS(3, NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)

    return
end subroutine shd_ext__mpiFindBoundaryPnt

subroutine shd_ext__mpiCalBoundaryRad(&
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
end subroutine shd_ext__mpiCalBoundaryRad

subroutine shd_ext__mpiComputeRad(&
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

    return
end subroutine shd_ext__mpiComputeRad

subroutine shd_ext__mpiVisualizeRad(&
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
end subroutine shd_ext__mpiVisualizeRad

subroutine shd_ext__mpiCalAccelSolCrit(&
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
end subroutine shd_ext__mpiCalAccelSolCrit

subroutine shd_ext__mpiEnd(NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART,&
    &XDOMAIN, YDOMAIN, NPXT, NPYT, PROPFILE, RUNNAME)

    integer NPTS, NPX, NPY, NPXT, NPYT
    real GRIDPOS(3, NPTS), XSTART, YSTART, XDOMAIN, YDOMAIN
    character(len=*), intent(in) :: PROPFILE, RUNNAME

    return
end subroutine shd_ext__mpiEnd

subroutine shd_ext__mpiAbort(ERRSTR)

    character(*) ERRSTR

    write (6, *) ERRSTR
    stop
end subroutine shd_ext__mpiAbort

end module shd_ext_mpi
