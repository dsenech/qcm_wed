COMPILER	= c++
LINKER	= $(COMPILER)
PYTHON_LIBDIR:=$(shell python3 -c 'from distutils import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')
LINK	= -L/usr/local/lib/ -lcuba -lomp -framework Accelerate
LOCAL_INCLUDE 	= $(shell python3-config --cflag) -w -I$(shell python3 -c 'import numpy; print(numpy.get_include())')/numpy
OPTIONS	= -std=c++14 -fPIC -Xpreprocessor -fopenmp 
LDFLAGS 	= -shared  -fPIC -undefined dynamic_lookup
DEP	=
EXEC	= $(HOME)/lib/qcm.so

#the common core of the makefiles
include ../qcm_object_list.txt
include ../makefile_nodep.mk
