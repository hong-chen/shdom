#!/bin/bash

mkdir -p tmp-data
fname_nml="data/shdom-inp_les.nml"
put $fname_nml | shdom
