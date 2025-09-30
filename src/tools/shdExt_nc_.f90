!╭────────────────────────────────────────────────────────────────────────────╮!
! @LICENSE_HEADER_START@
!
!   This file is part of SHDOMNG.
!
!   --
!   SHDOMNG: Spherical Harmonic Discrete Ordinate Method - Next Generation
!
!
!   Copyright (C) 2025-
!   Hong Chen, hong.chen-1@colorado.edu
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
! SHDOM was originally developed by Frank Evans and is available at https://coloradolinux.com/shdom.
!╰────────────────────────────────────────────────────────────────────────────╯!

module shdom_nc
    !
    ! Dummy SHDOM netcdf interface for linking without netcdf.
    !
    implicit none

contains

    subroutine read_property_size_netcdf(PROPFILE, NSTLEG, NLEG, NPX, NPY, NPZ, &
                                         NUMPHASE, MAXLEG, MAXPGL, DELX, DELY)
        character(len=256) :: propfile
        integer          :: NSTLEG, NLEG
        integer          :: NPX, NPY, NPZ, NUMPHASE, MAXLEG, MAXPGL
        real             :: DELX, DELY

        print *, 'Cannot read the sizes from the netcdf property file because SHDOM was compiled without the netcdf library.'
        stop
    end subroutine read_property_size_netcdf

    subroutine read_properties_netcdf(PROPFILE, NPX, NPY, NPZ, NUMPHASE, &
                                      NSTLEG, MAXLEG, MAXPG, MAXPGL, DELTAM, NLEG, &
                                      PROPTYPE, ZLEVELS, MAXASYM, &
                                      TEMPP, EXTINCTP, ALBEDOP, LEGENP, IPHASEP)
        character(len=*) :: propfile
        integer          :: NPX, NPY, NPZ, NUMPHASE
        integer          :: MAXPG, MAXPGL, NSTLEG, MAXLEG
        logical          :: DELTAM
        integer          :: NLEG
        character        :: PROPTYPE
        real             :: zLevels(:)
        real             :: MAXASYM
        real             :: TEMPP(:), EXTINCTP(:), ALBEDOP(:), LEGENP(:)
        integer*2        :: IPHASEP(:)

        print *, 'Cannot read the netcdf property file because SHDOM was compiled without the netcdf library.'
        stop
    end subroutine read_properties_netcdf

    subroutine output_results_netcdf_par(NSTOKES, NX, NY, NZ, NPTS, NCELLS, &
                                         NSH, ML, MM, NLM, NMU, NPHI, NANG, NG, &
                                         PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE, &
                                         BCFLAG, IPFLAG, DELTAM, GRIDTYPE, SRCTYPE, &
                                         SOLARFLUX, SOLARMU, SOLARAZ, SKYRAD, &
                                         SFCTYPE, GNDTEMP, GNDALBEDO, WAVENO, WAVELEN, UNITS, &
                                         SPLITACC, SHACC, SOLACC, MAXITER, TOTITER, CPUTIME, &
                                         XDOMAIN, YDOMAIN, XGRID, YGRID, ZGRID, &
                                         FLUXES, FLUXDIV, NSHOUT, SHTERMS, IRAD, RADOUT, &
                                         NUMOUT, OUTTYPES, OUTPARMS, OutFileNC)
        integer :: NSTOKES, NX, NY, NZ, NPTS, NCELLS
        integer :: NSH, ML, MM, NLM, NMU, NPHI, NANG, NG
        integer :: BCFLAG, IPFLAG, MAXITER, TOTITER, IRAD, NSHOUT, NUMOUT
        logical :: deltaM
        real    :: SOLARFLUX, SOLARMU, SOLARAZ, GNDTEMP, GNDALBEDO, SKYRAD
        real    :: SOLACC, SPLITACC, SHACC, XDOMAIN, YDOMAIN, WAVELEN, CPUTIME
        character(len=*) :: SRCTYPE, SFCTYPE, UNITS, GRIDTYPE
        character(len=*) :: PROPFILE, SFCFILE, CKDFILE, INSAVEFILE, OUTSAVEFILE
        real    :: XGRID(:), YGRID(:), ZGRID(:), WAVENO(:)
        real    :: fluxes(:, :, :, :), fluxdiv(:, :, :), SHTERMS(:, :, :, :), radout(:)
        character(len=1) :: OUTTYPES(:)
        real    :: OUTPARMS(:, :)
        character(len=*) :: OutFileNC

        print *, 'Not outputting a netcdf file because SHDOM was compiled without netcdf library.'
    end subroutine output_results_netcdf_par

end module shdom_nc
