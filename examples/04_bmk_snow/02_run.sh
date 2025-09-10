#!/bin/bash

#export OMP_NUM_THREADS=8

mkdir -p tmp-data

fname_nml="data/shdom-inp_snow_save.nml"
if [ "$1" == "mpi" ]; then
    mpirun -np 1 shdom $fname_nml
else
    shdom $fname_nml
fi

fname_nml="data/shdom-inp_snow_restore.nml"
if [ "$1" == "mpi" ]; then
    mpirun -np 1 shdom $fname_nml
else
    shdom $fname_nml
fi
