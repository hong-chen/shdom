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
