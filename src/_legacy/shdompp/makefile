# Makefile for SHDOMPP in distribution

FFLAGS= -O2 -fallow-argument-mismatch
FC= gfortran


all: shdompp disortsh ppmieprp

shdompp: shdompp.f90 fftpack5.1cos.o
	$(FC) $(FFLAGS) shdompp.f90 fftpack5.1cos.o -o shdompp

fftpack5.1cos.o: fftpack5.1cos.f90
	$(FC) $(FFLAGS) -c fftpack5.1cos.f90


disortsh: disortsh.f90 disort2.o
	$(FC) $(FFLAGS) disortsh.f90 disort2.o  -o disortsh

ppmieprp: ppmieprp.f90 
	$(FC) $(FFLAGS) ppmieprp.f90 -o ppmieprp


.f.o :; $(FC) -c $(FFLAGS) $*.f

clean:
	@ rm -f *.o *.mod *~ shdompp disortsh ppmieprp
