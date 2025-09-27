#!/bin/csh
#   Creates the hydrometeor single scattering tables for GOES imager channels.
# Loops over the channels and a specified subset of three hydrometeor species.
# The single scattering properties are either from Mie theory (for spheres)
# or from Ping Yang's database (for combinations of his 6 shapes).
# The single scattering properties are integrated over gamma particle 
# size distributions.  Each scattering table has a range of log spaced
# mean mass radius for one specified gamma distribution width parameter (alpha).
# Currently the scattering tables are prepared for the GOES channel
# central wavelengths, though it might be better to integrate across 
# the approximation channel bandwidth.
# The ice particle shape mixture files are now set only to use
# Bryan Baum's mixture recipe:
#    http://www.ssec.wisc.edu/~baum/Cirrus/MODIS_V2models.html
# The maximum for this mixture is mean mass radius is about 750 um.
#   Requires: make_mie_table2 and make_ice_table2 to be compiled,
# and the input files Baum_MODIS_shape_mixture.dat and yang2005_ice_scatter.db.

  # Hydrometeor species to compute:
  #   species = (cloud  ice  snow)
set ispecies = (1 2 3)

set scattabledir="."

  # GOES Imager central wavelegths
set cwavelen=(0.6288 3.9225 6.7685 10.7247 12.0105)


set Nchan=$#cwavelen


set speciesname=(cloud  ice    snow)
set alpha   =   ( 7      1      1)

  # starting and ending mean mass radius in tables (micron)
set Srtab   =   ( 1.0     5     25) 
set Ertab   =   (25.0   400    750)
  # number of mean mass radius in table
set Nrtab   =   ( 29     39     30)

  # scattering code: 0 = make_mie_table2, 1 = make_ice_table2
set code    =  (   0     1      1)

  # maximum radius in Mie calculation (micron)
set maxradius = (  50  1000   2000)

  #  particle substance in Mie calculation (W for water, I for ice)
set partype =   (   W     I      I )

  # shape mixture file for make_ice_table2
set speciesname=(cloud  ice   snow)
set shapemixfile = (none  Baum_MODIS_shape_mixture.dat \
                    Baum_MODIS_shape_mixture.dat)




set k=1
while ($k <= $Nchan)   # Loop over channels (k)
  set chnum=$k

  set j=1
  while ($j <= $#ispecies)  # Loop over hydrometeors
    set i = $ispecies[$j]   # i is the hydrometeor number

    set scatfile = "${scattabledir}/goes_ch${chnum}_$speciesname[$i].scat"

    set wavelen1 = $cwavelen[$k]
    set wavelen2 = $cwavelen[$k]

    if ($code[$i] == 0) then
       #  Makes the Mie scattering table for spherical water droplets for
       # gamma size distributions with a range of mean mass radius. 
      set distflag=G        # G=gamma, L=lognormal size distribution
      set avgflag=C         # A for averaging over wavelength range
                            #  C for using central wavelength
       # set deltawave=0.003   # wavelength step for averaging (microns)


      if ($avgflag == "A") then
        put "$wavelen1 $wavelen2" $partype[$i] $avgflag $deltawave $distflag  \
          $alpha[$i] "$Nrtab[$i] $Srtab[$i] $Ertab[$i]" $maxradius[$i] $scatfile \
         | bin/make_mie_table2
      else
        put "$wavelen1 $wavelen2" $partype[$i] $avgflag  $distflag  \
          $alpha[$i] "$Nrtab[$i] $Srtab[$i] $Ertab[$i]" $maxradius[$i] $scatfile \
         | bin/make_mie_table2
      endif
    endif


    if ($code[$i] == 1) then
      #  Makes scattering tables for several ice particle shapes for 
      # gamma size distributions with a range of effective radius. 
      #    Shape index:
      put  yang2005_ice_scatter.db  "$wavelen1 $wavelen2" $shapemixfile[$i] \
        "$Nrtab[$i] $Srtab[$i] $Ertab[$i]" $alpha[$i] $scatfile \
       | bin/make_ice_table2
    endif

    echo; echo

    @ j++
  end


  @ k++
end
