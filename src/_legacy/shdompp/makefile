# Makefile for SHDOMPP in distribution

FFLAGS= -O2 -fallow-argument-mismatch
FC= gfortran

CURRENT_DIR = $(shell pwd)
BIN_DIR     = $(abspath $(CURRENT_DIR)/bin)

all: bin shdompp_new shdompp disortsh ppmieprp put

bin:
	mkdir -p $(BIN_DIR)

shdompp_new: shdompp_new.f90 fftpack5.1cos.o
	$(FC) $(FFLAGS) shdompp_new.f90 fftpack5.1cos.o -o $(BIN_DIR)/shdompp_new

shdompp: shdompp.f90 fftpack5.1cos.o
	$(FC) $(FFLAGS) shdompp.f90 fftpack5.1cos.o -o $(BIN_DIR)/shdompp

fftpack5.1cos.o: fftpack5.1cos.f90
	$(FC) $(FFLAGS) -c fftpack5.1cos.f90

disortsh: disortsh.f90 disort4.o
	$(FC) $(FFLAGS) disortsh.f90 disort4.o  -o $(BIN_DIR)/disortsh

ppmieprp: ppmieprp.f90
	$(FC) $(FFLAGS) ppmieprp.f90 -o $(BIN_DIR)/ppmieprp

put:
	$(FC) $(FFLAGS) put.f90 -o $(BIN_DIR)/put

.f.o :; $(FC) -c $(FFLAGS) $*.f

clean:
	-/bin/rm -rf *.mod *.o *~ bin
