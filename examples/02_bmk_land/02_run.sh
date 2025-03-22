#!/bin/bash

export OMP_NUM_THREADS=8

mkdir -p tmp-data
fname_nml="data/shdom-inp_land_save.nml"
shdom $fname_nml
fname_nml="data/shdom-inp_land_restore.nml"
shdom $fname_nml

#mpirun -np 4 shdom data/shdom-inp_land.nml
