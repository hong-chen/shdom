# Makefile for SHDOMPPDA distribution

CompileF90 = gfortran -c -Ofast -fallow-argument-mismatch
LinkF90    = gfortran -Ofast

all: demo_shdomppda make_mie_table2 make_ice_table2


Objs =  molec_absorption_interface.o shdomppda_fwd_tl_adj.o \
        shdomppda.o shdomppda_ftls.o shdomppda_ads.o

demo_shdomppda: demo_shdomppda.f90  $(Objs)
	$(LinkF90)  -o $@ $< $(Objs)

molec_absorption_interface.o: molec_absorption_interface.f90
	$(CompileF90)  $<

shdomppda_fwd_tl_adj.o: shdomppda_fwd_tl_adj.f90
	$(CompileF90)  $<

shdomppda.o: shdomppda.f90
	$(CompileF90)  $<

shdomppda_ftls.o: shdomppda_ftls.f90
	$(CompileF90)  $<

shdomppda_ads.o: shdomppda_ads.f90
	$(CompileF90)  $<



make_mie_table2: make_mie_table2.f90 mieindsub.o
	$(LinkF90) -o $@ $< mieindsub.o
        
mieindsub.o: mieindsub.f
	$(CompileF90) $<


make_ice_table2: make_ice_table2.f90
	$(LinkF90) -o $@ $< 

clean::
	-/bin/rm  *.mod *.o  demo_shdomppda make_mie_table2 make_ice_table2

