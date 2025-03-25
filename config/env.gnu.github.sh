##/bin/bash -xe

# config for Netcdf
export NC_INCLUDE=`nf-config --fflags|cut -d " " -f1 | cut -d "I" -f 2`
export NC_LIBS=`nf-config --flibs`

# config for HDF4 
export H4_INCLUDE=""
export H4_LIBS=`h4fc -show TESTSRC | awk -F"TESTSRC " '{print $2}'`


# config for HDF5
export H5_INCLUDE=`h5fc -show TESTSRC | awk -F"TESTSRC" '{print $1}' | awk -F"-I" '{print $2}'|cut -d " " -f 1`
export H5_LIBS=`h5fc -show TESTSRC | awk -F"TESTSRC " '{print $2}'`


echo "====================================="
echo "          configurations:"
echo 
echo "NC_INCLUDE=$NC_INCLUDE"
echo "NC_LIBS=$NC_LIBS"

echo "H4_INCLUDE=$H4_INCLUDE"
echo "H4_LIBS=$H4_LIBS"

echo "H5_INCLUDE=$H5_INCLUDE"
echo "H5_LIBS=$H5_LIBS"
