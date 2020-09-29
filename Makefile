PYTHON = python3.7
SRC_FILES = $(wildcard src/*.f90)
LIB_FILES = ${SRC_FILES:src/%.f90=build/%.so}

%.so: %.f90
	$(PYTHON) -m numpy.f2py -c $^

build: 
	$(PYTHON) -m numpy.f2py -c src/spintronics.f90 -m spintronics