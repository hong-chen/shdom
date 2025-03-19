#!/bin/csh
# Script for demonstrating visualization mode output of SHDOM.
# Uses the stratocumulus cloud from phase 2 of the I3RC.
# Warning: the SHDOM runs are CPU intensive (5.5 minutes for the SHDOM
# solution and 33 minutes for the 100 frame camera visualization on
# a 3.0 GHz Intel Xeon processor).
# The PGM images may be viewed with ImageMagick "display" and many other
# image viewers.  The PGM images from the visualization may be converted
# to other formats with ImageMagick "convert".

set AWK=awk

set fdir_data="data"
set fdir_data_tmp="tmp-data"
mkdir $fdir_data_tmp

# Set each of these to 1 in turn to run the
# sections of the script.
set MakeMieTable=1
set MakePropertyFile=1
set RunSHDOMsolve=0
set SingleFrames=0
set RunCameraVisual=0

set basefile="modis"
set prpfile="${fdir_data_tmp}/shdom-prp_${basefile}.txt"
set savefile=${basefile}a.bin
set outbase=${basefile}a

set wavelen=0.65
set wavename=065

if ($MakeMieTable) then
  #  Makes the Mie scattering table for spherical water droplets for
  # gamma size distributions with a range of effective radius.
  set polarizedtable=F             # F=unpolarized output, T=polarized output

  set outfile="${fdir_data_tmp}/shdom-sca-${wavename}_${basefile}.txt"

  set partype="W"  # W for water

  set distflag=G   # G=gamma, L=lognormal size distribution

  set alpha=7 i    # gamma dist shape parameter

  set Nretab=50    # number of effective radius in table

  set Sretab=0.5

  set Eretab=25    # starting, ending effective radius (micron)

  set logre=F      # evenly spaced r_e

  set maxradius=75

  set avgflag=C    #  C for using central wavelength

  put $polarizedtable "$wavelen $wavelen" $partype $avgflag $distflag $alpha \
      "$Nretab $Sretab $Eretab" $logre $maxradius "$outfile" | make_mie_table

endif



if ($MakePropertyFile) then

  set lwcfile="${fdir_data}/shdom-lwc_${basefile}.txt"

  set scattable="${fdir_data_tmp}/shdom-sca-${wavename}_${basefile}.txt"

  set poltabs = (F)

  set maxnewphase=1000

  set asymtol=0.01

  set fracphasetol=0.1

  set sfcpres=1013.0

  set Nzother=2

  set ZTother=(0.0 286.3  20.0 180.0)

  set polout=U

  put 1 $scattable 1 "$poltabs" "$lwcfile" \
      $maxnewphase $asymtol $fracphasetol $wavelen $sfcpres \
      $Nzother $ZTother $polout $prpfile  | propgen

endif

set Nx=`$AWK 'BEGIN {getline; getline; print $1;}' $prpfile`
set Ny=`$AWK 'BEGIN {getline; getline; print $2;}' $prpfile`
set Nz=`$AWK 'BEGIN {getline; getline; print $3;}' $prpfile`

echo $Nx
echo $Ny
echo $Nz
exit
