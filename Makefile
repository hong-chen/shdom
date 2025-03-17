# This file is referenced from https://github.com/cd10kfsu/F90GIO
all:
	mkdir build
	. config/env.sh; cd build; cmake ..; make; make test
