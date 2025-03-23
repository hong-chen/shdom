#!/bin/bash

export OMP_NUM_THREADS=8

mkdir -p tmp-data

#fname_nml="data/shdom-inp_les_save.nml"
#shdom $fname_nml
#mpirun -np 4 shdom $fname_nml

fname_nml="data/shdom-inp_les_rad.nml"
shdom $fname_nml
#mpirun -np 4 shdom $fname_nml

#fname_nml="data/shdom-inp_les_flux.nml"
#shdom $fname_nml
#mpirun -np 4 shdom $fname_nml

