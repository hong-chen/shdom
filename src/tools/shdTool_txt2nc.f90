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

program PropertyToNetcdf
  use shdom_nc
  implicit none

  character(len=256) :: asciiFile, netcdfFile
  integer              :: numArgs
  ! -------------------------------------
  numArgs = command_argument_count()
  if(numArgs < 1) stop "Error [PropertyToNetcdf]: Need to provide an input file name."
  call get_command_argument(1, asciiFile)
  if(numArgs < 2) then
    netcdfFile = trim(asciiFile) // ".nc"
  else
    call get_command_argument(2, netcdfFile)
  end if

  print *, "Message [PropertyToNetcdf]: Converting ascii file ", trim(asciiFile), " to netcdf file ", trim(netcdfFile)
  call convert_Prp_to_netcdf(asciiFile, netcdfFile)

end program PropertyToNetcdf
