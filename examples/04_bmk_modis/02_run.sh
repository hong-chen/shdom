#!/bin/bash

#export OMP_NUM_THREADS=8

mkdir -p tmp-data

fname_nml="data/shdom-inp_modis.nml"
if [ "$1" == "mpi" ]; then
    mpirun -np 4 shdom $fname_nml
else
    shdom $fname_nml
fi
