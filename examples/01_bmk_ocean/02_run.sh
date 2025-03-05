#!/bin/bash

mkdir -p tmp-data
fname_nml="data/shdom-inp_ocean.nml"
put $fname_nml | shdom
