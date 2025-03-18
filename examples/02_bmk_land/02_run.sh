#!/bin/bash

mkdir -p tmp-data
fname_nml="data/shdom-inp_land.nml"
#shdom $fname_nml
mpirun -np 4 shdom data/shdom-inp_land.nml
