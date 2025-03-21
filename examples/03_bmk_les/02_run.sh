#!/bin/bash

mkdir -p tmp-data

fname_nml="data/shdom-inp_les_rad.nml"
#fname_nml="data/shdom-inp_les_flux.nml"

shdom $fname_nml
#mpirun -np 4 shdom $fname_nml
