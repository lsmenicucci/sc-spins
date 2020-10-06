PYTHON = python3.7
SRC_FILES = $(wildcard src/*.f90)
EXT_SUFFIX := $(shell python3.7-config --extension-suffix)
LIB_FILES = ${SRC_FILES:src/%.f90=%}

%$(EXT_SUFFIX): src/%.f90
	$(PYTHON) -m numpy.f2py -c --f90exec=mpif90 -m $* $<

default: $(LIB_FILES)$(EXT_SUFFIX)
