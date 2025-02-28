program PropertyToNetcdf
  use shdom_netcdf
  implicit none 
  
  interface 
    subroutine getarg(pos, value)
      integer,            intent(in)  :: pos
      character(len = *), intent(out) :: value
    end subroutine getarg
    
    function iargc()
      integer :: iargc
    end function iargc
  end interface

  character(len = 128) :: asciiFile, netcdfFile
  integer              :: numArgs
  ! -------------------------------------
  numArgs = iargc() 
  if(numArgs < 1) stop "Need to provide an input file name." 
  call getarg(1, asciiFile) 
  if(numArgs < 2) then 
    netcdfFile = trim(asciiFile) // ".nc"
  else
    call getarg(2, netcdfFile)
  end if 
  
  print *, "Converting ascii file ", trim(asciiFile), " to netcdf file ", trim(netcdfFile)
  call convert_Prp_to_netcdf(asciiFile, netcdfFile)
  
end program PropertyToNetcdf