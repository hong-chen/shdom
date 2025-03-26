# This file is referenced from https://github.com/cd10kfsu/F90GIO
CURRENT_DIR = $(shell pwd)
BUILD_DIR = $(abspath $(CURRENT_DIR)/build)

all:
	rm -rf $(BUILD_DIR)
	mkdir $(BUILD_DIR)
	cd $(BUILD_DIR); cmake $(CURRENT_DIR); make VERBOSE=1

clean:
	rm -rf $(BUILD_DIR)
