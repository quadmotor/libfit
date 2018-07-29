

TOP_DIR = $(shell pwd)
SRC_DIR = $(TOP_DIR)/src
BUILD_DIR = $(TOP_DIR)/build
EXAMPLES_DIR = $(TOP_DIR)/examples
INCLUDE_DIR = $(TOP_DIR)/include
LIB_DIR = $(TOP_DIR)/lib

FC = gfortran
FFLAGS = -O2 -Ofast -g -I$(INCLUDE_DIR)

OBJECTS = $(BUILD_DIR)/lsqr_data.o \
	$(BUILD_DIR)/lsqrblas.o \
	$(BUILD_DIR)/lsqr.o \
	$(BUILD_DIR)/block.o \
	$(BUILD_DIR)/grid.o \
	$(BUILD_DIR)/param.o \
	$(BUILD_DIR)/bezier_helper.o \
	$(BUILD_DIR)/bezier.o \
	$(BUILD_DIR)/numpy.o \
	$(BUILD_DIR)/plot3d.o \

include $(EXAMPLES_DIR)/Makefile

.phony: all
all: $(BUILD_DIR)/libFit.a $(EXAMPLES_DIR)/examples

$(LIB_DIR)/libFit.a: $(OBJECTS)
	ar rc $@ $^
	
$(BUILD_DIR)/%.o: $(SRC_DIR)/%.f90
	$(FC) $(FFLAGS) -c $^ -o $@
	-@mv *.mod $(INCLUDE_DIR) 2>/dev/null || true

clean: $(EXAMPLES_DIR)/clean
	-rm $(BUILD_DIR)/* 2>/dev/null || true
	-rm $(INCLUDE_DIR)/* 2>/dev/null || true
	-rm $(LIB_DIR)/* 2>/dev/null || true
