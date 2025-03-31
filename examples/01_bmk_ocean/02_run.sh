#!/bin/bash

mkdir -p tmp-data

fname_nml="data/shdom-inp_ocean.nml"
#shdom $fname_nml
mpirun -np 4 shdom $fname_nml
