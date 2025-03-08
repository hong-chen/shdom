subroutine START_MPI(MASTERPROC)
    ! Dummy routines to allow parallel SHDOM to compile without MPI.
    implicit none
    logical MASTERPROC

    MASTERPROC = .true.
    return
end subroutine

subroutine MAP_SHDOM_MPI( &
    & BCFLAG, NPX, NPY, NX, NY, NZ, DELX, DELY, &
    & PROPFILE, XSTART, YSTART)
    implicit none
    integer BCFLAG, NPX, NPY, NX, NY, NZ
    real DELX, DELY
    character PROPFILE*256
    real XSTART, YSTART

    XSTART = 0.0
    YSTART = 0.0
    return
end subroutine

subroutine BROADCAST_USER_INPUT(RUNNAME, &
    & PROPFILE, SFCFILE, CKDFILE, &
    & INSAVEFILE, OUTSAVEFILE, &
    & NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, &
    & IPFLAG, KDIST, DELTAM, GRIDTYPE, &
    & SRCTYPE, SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
    & GNDTEMP, GNDALBEDO, UNITS, WAVENO, WAVELEN, &
    & ACCELFLAG, SOLACC, MAXITER, SPLITACC, SHACC, &
    & MAXOUT, MAXPAR, NUMOUT, OUTTYPES, OUTPARMS, OUTFILES)
    implicit none
    integer NSTOKES, NX, NY, NZ, NMU, NPHI, BCFLAG, IPFLAG
    integer MAXOUT, MAXPAR, MAXITER, NUMOUT
    logical KDIST, DELTAM, ACCELFLAG
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GNDTEMP, GNDALBEDO, SKYRAD, WAVENO(2), WAVELEN
    real SOLACC, SPLITACC, SHACC, OUTPARMS(MAXPAR, MAXOUT)
    character SRCTYPE*1, GRIDTYPE*1, UNITS*1, OUTTYPES(*)*1
    character RUNNAME*256, PROPFILE*256, SFCFILE*256, CKDFILE*256
    character OUTFILES(*)*256, INSAVEFILE*256, OUTSAVEFILE*256

    return
end subroutine

subroutine BROADCAST_PROPERTY_SIZE(NPX, NPY, NPZ, DELX, DELY, &
    & NUMPHASE, MAXLEG, MAXPGL)
    implicit none
    integer NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
    real DELX, DELY

    return
end subroutine

subroutine SCATTER_PROPERTIES(NPXT, NPYT, NPX, NPY, NPZ, NSTLEG, &
    & NLEG, NUMPHASE, ZLEVELS, MAXASYM, TEMPPT, &
    & EXTINCTPT, ALBEDOPT, LEGENPT, IPHASEPT, &
    & TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)
    implicit none
    integer NPXT, NPYT, NPX, NPY, NPZ, NUMPHASE, NLEG, NSTLEG
    real ZLEVELS(NPZ), MAXASYM
    real TEMPPT(NPZ, NPYT, NPXT), EXTINCTPT(NPZ, NPYT, NPXT)
    real ALBEDOPT(NPZ, NPYT, NPXT), LEGENPT(*)
    integer*2 IPHASEPT(NPZ, NPYT, NPXT)
    real TEMPP(NPZ, NPY, NPX), EXTINCTP(NPZ, NPY, NPX)
    real ALBEDOP(NPZ, NPY, NPX), LEGENP(*)
    integer*2 IPHASEP(NPZ, NPY, NPX)

    return
end subroutine

subroutine READ_BROADCAST_MEM_PARAMS(MAX_TOTAL_MB, &
    & ADAPT_GRID_FACTOR, NUM_SH_TERM_FACTOR, &
    & CELL_TO_POINT_RATIO, RUNNAME)
    implicit none
    real MAX_TOTAL_MB, ADAPT_GRID_FACTOR
    real NUM_SH_TERM_FACTOR, CELL_TO_POINT_RATIO
    character(*) RUNNAME

    return
end subroutine

subroutine BROADCAST_SURFACE_SIZE(MAXSFCPTS, MAXSFCPARS)
    implicit none
    integer MAXSFCPTS, MAXSFCPARS

    return
end subroutine

subroutine BROADCAST_SURFACE_PARMS(SFCTYPE, NXSFC, NYSFC, NSFCPAR, &
    & DELXSFC, DELYSFC, SFCPARMS, GNDTEMP, GNDALBEDO)
    implicit none
    character SFCTYPE*2
    integer NXSFC, NYSFC, NSFCPAR
    real DELXSFC, DELYSFC, SFCPARMS(*), GNDTEMP, GNDALBEDO

    return
end subroutine

subroutine BROADCAST_KDIST_SIZE(NG, NZCKD)
    implicit none
    integer NG, NZCKD

    return
end subroutine

subroutine BROADCAST_KDIST_PARMS(SOLFLUX, NG, DELG, &
    & NZCKD, ZCKD, KABS)
    implicit none
    integer NG, NZCKD
    real SOLFLUX, DELG(NG), ZCKD(NZCKD), KABS(NZCKD, NG)

    return
end subroutine

subroutine GATHER_OUTPUT(FLUX_OUT, FLUXDIV_OUT, SH_OUT, &
    & IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, &
    & SUMFLUXES, SUMDIRFLUX, SUMFLUXDIV, SUMSHTERMS, &
    & NPX, NPY, DELX, DELY, XALLGRID, YALLGRID, &
    & ALLFLUXES, ALLFLUXDIV, ALLSHTERMS, &
    & NCELLS, NPTS, NSH, NCELLSTOT, NPTSTOT, NSHTOT)
    implicit none
    logical FLUX_OUT, FLUXDIV_OUT, SH_OUT
    integer IPFLAG, BCFLAG, NBPTS, NXT, NYT, NZ, NSHOUT, NPX, NPY
    real DELX, DELY
    real SUMFLUXES(2, NBPTS), SUMDIRFLUX(NBPTS)
    real SUMFLUXDIV(NBPTS), SUMSHTERMS(NSHOUT, NBPTS)
    real XALLGRID(NXT), YALLGRID(NYT)
    real ALLFLUXES(3, NZ, NYT, NXT), ALLFLUXDIV(NZ, NYT, NXT)
    real ALLSHTERMS(NSHOUT, NZ, NYT, NXT)
    integer NCELLS, NPTS, NSH
    integer NCELLSTOT, NPTSTOT, NSHTOT

    return
end subroutine

real function SUM_CPU_TIME(cpuTimes)
    implicit none
    real cpuTimes
    SUM_CPU_TIME = cpuTimes
    return
end function

subroutine TOTAL_ALBEDO_MAX(ALBMAX)
    implicit none
    real ALBMAX
    return
end subroutine

subroutine UNIFY_SPLITTING(DOSPLIT, STARTSPLITACC)
    implicit none
    logical DOSPLIT
    real STARTSPLITACC
    return
end subroutine

subroutine TOTAL_SPLITCRIT_MAX(SPLITCRIT)
    implicit none
    real SPLITCRIT
    return
end subroutine

subroutine MAKE_DIRECT_PAR( &
    & SPT, NPTS, BCFLAG, IPFLAG, DELTAM, ML, &
    & NSTLEG, NLEG, SOLARFLUX, SOLARMU, SOLARAZ, GRIDPOS, &
    & NX, XGRID, NY, YGRID, DIRFLUX)
    implicit none
    integer SPT, NPTS, BCFLAG, IPFLAG, ML, NSTLEG, NLEG, NX, NY
    logical DELTAM
    real SOLARFLUX, SOLARMU, SOLARAZ
    real GRIDPOS(3, NPTS), XGRID(NX + 1), YGRID(NY + 1)
    real DIRFLUX(NPTS)

    return
end subroutine

subroutine FIND_BOUNDARY_POINTS( &
    & BCFLAG, IPFLAG, NPTS, SWEEPORD, &
    & GRIDPTR, GRIDPOS, NX, NY, NZ, XGRID, YGRID, ZGRID)
    implicit none
    integer BCFLAG, IPFLAG, NPTS, SWEEPORD(NPTS, *), GRIDPTR(8, *)
    integer NX, NY, NZ
    real GRIDPOS(3, NPTS), XGRID(NX), YGRID(NY), ZGRID(NZ)

    return
end subroutine

subroutine CALC_BOUNDARY_RADIANCES( &
    & BCFLAG, IPFLAG, JOCT, IZ, &
    & NX, NY, NZ, XGRID, YGRID, ZGRID, &
    & NA, NPTS, NCELLS, GRIDPTR, &
    & NEIGHPTR, TREEPTR, CELLFLAGS, GRIDPOS, &
    & MU, PHI, EXTINCT, NSTOKES, SOURCE, &
    & KANG, GRIDRAD)
    implicit none
    integer BCFLAG, IPFLAG, JOCT, IZ
    integer NX, NY, NZ, NA, NPTS, NCELLS, NSTOKES, KANG
    real XGRID(NX), YGRID(NY), ZGRID(NZ)
    integer GRIDPTR(8, NCELLS), NEIGHPTR(6, NCELLS), TREEPTR(2, NCELLS)
    integer*2 CELLFLAGS(NCELLS)
    real GRIDPOS(3, NPTS), MU, PHI
    real EXTINCT(NPTS), SOURCE(NSTOKES, NA, NPTS)
    real GRIDRAD(NSTOKES, NPTS)

    return
end subroutine

subroutine COMPUTE_RADIANCE_PAR( &
    & NSTOKES, NX, NY, NZ, NPTS, NCELLS, &
    & ML, MM, NSTLEG, NLEG, NUMPHASE, &
    & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
    & BCFLAG, XDOMAIN, YDOMAIN, IPFLAG, &
    & SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
    & SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
    & MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
    & GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
    & XGRID, YGRID, ZGRID, GRIDPOS, &
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
    & EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
    & SHPTR, SOURCE, SOURCE1, GRIDRAD, &
    & OUTPARMS, NRAD, RADOUT)
    implicit none
    integer NSTOKES, NX, NY, NZ, BCFLAG, IPFLAG, NPTS, NCELLS
    integer ML, MM, NSTLEG, NLEG, NUMPHASE
    integer NMU, NPHI0MAX, NPHI0(NMU), NRAD
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
    real EXTINCT(NPTS), ALBEDO(NPTS), LEGEN(NSTLEG, 0:NLEG, NPTS)
    real DIRFLUX(NPTS), FLUXES(2, NPTS), SOURCE(*)
    real SOURCE1(NSTOKES, NPTS), GRIDRAD(NSTOKES, NPTS)
    real OUTPARMS(*), RADOUT(NSTOKES, *)
    character SRCTYPE*1, SFCTYPE*2, UNITS*1

    return
end subroutine

subroutine VISUALIZE_RADIANCE_PAR( &
    & NSTOKES, NX, NY, NZ, NPTS, NCELLS, &
    & ML, MM, NLM, NSTLEG, NLEG, NUMPHASE, &
    & NMU, NPHI0MAX, NPHI0, MU, PHI, WTDO, &
    & BCFLAG, XDOMAIN, YDOMAIN, IPFLAG, &
    & SRCTYPE, DELTAM, SOLARMU, SOLARAZ, &
    & SFCTYPE, NSFCPAR, SFCGRIDPARMS, &
    & MAXNBC, NTOPPTS, NBOTPTS, BCPTR, BCRAD, &
    & GNDTEMP, GNDALBEDO, SKYRAD, WAVENO, WAVELEN, UNITS, &
    & XGRID, YGRID, ZGRID, GRIDPOS, &
    & GRIDPTR, NEIGHPTR, TREEPTR, CELLFLAGS, &
    & EXTINCT, ALBEDO, LEGEN, IPHASE, DIRFLUX, FLUXES, &
    & SHPTR, SOURCE, OUTPARMS, IVIS, VISOUT)
    implicit none
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
end subroutine

subroutine CALC_ACCEL_SOLCRIT( &
    & DOACCEL, DELJDOT, DELJOLD, DELJNEW, &
    & JNORM, ACCELPAR, SOLCRIT)
    ! Calculates the acceleration parameter and solution criterion from
    ! the delta source function vector dot  products.
    implicit none
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
        A = (1 - R*cos(THETA) + R**(1 + 0.5*3.14159/THETA)) &
            & /(1 + R**2 - 2*R*cos(THETA)) - 1.0
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
end subroutine

subroutine END_SHDOM_MPI(NPTS, GRIDPOS, NPX, NPY, XSTART, YSTART, &
    & DELX, DELY, NPXT, NPYT, PROPFILE)
    implicit none
    integer NPTS, NPX, NPY, NPXT, NPYT
    real GRIDPOS(3, NPTS), XSTART, YSTART, DELX, DELY
    character PROPFILE*256

    return
end subroutine

subroutine ABORT_SHDOM_MPI(ERRSTR)
    implicit none
    character(*) ERRSTR

    write (6, *) ERRSTR
    stop
end subroutine
