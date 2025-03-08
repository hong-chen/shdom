#!/bin/bash

mkdir -p tmp-data
fname_nml="data/shdom-inp_land.nml"
put $fname_nml | shdom
