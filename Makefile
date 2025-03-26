# This file is referenced from https://github.com/cd10kfsu/F90GIO
CURRENT_DIR = $(shell pwd)
BUILD_DIR = $(abspath $(CURRENT_DIR)/build)

all:
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
	cd $(BUILD_DIR); \
		cmake $(CURRENT_DIR); \
		make

omp:
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
	cd $(BUILD_DIR); \
		cmake $(CURRENT_DIR) \
			-DUSE_OMP=on; \
		make

mpi:
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
	cd $(BUILD_DIR); \
		cmake $(CURRENT_DIR) \
			-DUSE_MPI=on; \
		make VERBOSE=1

nc:
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
	cd $(BUILD_DIR); \
		cmake $(CURRENT_DIR) \
			-DUSE_NETCDF=on \
			-DNC_INC="/opt/homebrew/Cellar/netcdf-fortran/4.6.1_1/include" \
			-DNC_LIBS="-L/opt/homebrew/Cellar/netcdf-fortran/4.6.1_1/lib -lnetcdff"; \
		make VERBOSE=1

clean:
	rm -rf $(BUILD_DIR)
