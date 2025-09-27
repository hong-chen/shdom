#!/bin/csh
# Script to run SHDOMPPDA forward and adjoint demonstration program.
# demo_shdomppda reads an ascii file containing the RAMS height, pressure, 
# temperature, water vapor mixing ratio, and mixing ratio (g/kg) and 
# number concentration (#/kg) for a number of hydrometer species 
# (e.g. cloud water and pristine ice) for each grid point.  
# demo_shdomppda also reads the scattering tables produced by 
# make_mie_table2 and make_ice_table2.  It outputs a file of radiances 
# for each column and a file of the adjoint with respect to cloud 
# species mixing ratio and number concentration.


 # Set variable to 1 to run that section:
set SolarAdjointDemo=1
set ThermalAdjointDemo=0


 # Compile the "put" command
if (!(-e put))  cc -o put  put.c


if ($SolarAdjointDemo) then
  # Runs demo_shdomppda for solar RT on a RAMS slice with 
  #   cloud droplets and pristine ice

  set wavelen=0.63
  set Nmu=8
  set instatefile=rams_twomom2_y360.dat
  set outradfile=twomom2_rad${wavelen}_nmu${Nmu}.dat
  set outstatefile=twomom2_adjoint${wavelen}_nmu${Nmu}.dat

  set ncomp=2
  set icomps=(1 3)     # components to use in input state file
  set scatfiles=(goes_ch1_cloud.scat goes_ch1_ice.scat)
  set speciesZrange=(0.0 4.5  3.0 8.0)
  set muobs=(0.707 0.866 0.966  0.707 0.866 0.966)
  set phiobs=(150    150  150   30    30    30)   
  set radpertvect=(0.01 0.01 0.01  0.01 0.01 0.01)
  set solacc=1.0E-6
  set splitacc=3.0E-3
  @ Nphi = 2 * $Nmu

   # Set some parameters for the radiative transfer
  set solarflux=1.0      # solar flux on a horizontal surface
  set solarmu=0.707      # cosine solar zenith angles
  set sfcalb=0.20        # Lambertian surface albedo
  set units=R            # reflectance units

   # Run demo_shdomppda
  put $instatefile $outradfile $outstatefile  \
      $ncomp $icomps  $scatfiles $speciesZrange \
      S "$solarflux $solarmu" $sfcalb $wavelen $units \
      "$Nmu $Nphi"  $splitacc $solacc \
      $#muobs "$muobs" "$phiobs"  "$radpertvect"  \
    | /usr/bin/time bin/demo_shdomppda
endif



if ($ThermalAdjointDemo) then
  # Runs demo_shdomppda for solar thermal on a RAMS slice with 
  # cloud droplets and pristine ice

  set wavelen=10.73
  set Nmu=4
  set instatefile=rams_twomom2_y360.dat
  set outradfile=twomom2_rad${wavelen}_nmu${Nmu}.dat
  set outstatefile=twomom2_adjoint${wavelen}_nmu${Nmu}.dat

  set ncomp=2
  set icomps=(1 3)     # components to use in input state file
  set scatfiles=(goes_ch4_cloud.scat goes_ch4_ice.scat)
  set speciesZrange=(0.0 4.5  3.0 8.0)
  set muobs=(0.500 0.707 0.866 0.966 1.00)
  set phiobs=(0     0     0     0     0  )
  set radpertvect=(1.0 1.0 1.0 1.0 1.0)
  set solacc=1.0E-6
  set splitacc=2.0E-3

   # Set some parameters for the radiative transfer
  set sfctemp=290.0      # surface temperature
  set sfcalb=0.02        # surface albedo
  set units=T            # brightness temperature units (K)
  set Nphi=1

   # Run demo_shdomppda
  put $instatefile $outradfile $outstatefile  \
      $ncomp $icomps  $scatfiles $speciesZrange \
      T $sfctemp $sfcalb $wavelen $units \
      "$Nmu $Nphi"  $splitacc $solacc \
      $#muobs "$muobs" "$phiobs"  "$radpertvect" \
    | bin/demo_shdomppda
endif

