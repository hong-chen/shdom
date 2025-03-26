# This file is referenced from https://github.com/cd10kfsu/F90GIO
CURRENT_DIR = $(shell pwd)
BUILD_DIR = $(abspath $(CURRENT_DIR)/build)

all:
	mkdir build
	cd build; cmake ..; make VERBOSE=1
