COMPILER	= c++
LINKER	= $(COMPILER)
PYTHON_LIBDIR:=$(shell python -c 'from distutils import sysconfig; print(sysconfig.get_config_var("LIBDIR"))')

# LINK	= $(shell python3-config --ldflags) -L/usr/local/lib/ -lcuba -L/opt/homebrew/Caskroom/miniforge/base/lib -lomp -lpython3.9 -lopenblas
LINK	= $(shell python3-config --ldflags) -L/usr/local/lib/ -lcuba -lomp -L$(PYTHON_LIBDIR) -lpython3.9 -lopenblas
# LINK	= $(shell python3-config --ldflags) -L/usr/local/lib/ -lcuba -lomp -L$(PYTHON_LIBDIR) -lpython3.9 -framework Accelerate

# LOCAL_INCLUDE 	= $(shell python3-config --cflag) -I$(shell python3 -c 'import numpy; print(numpy.get_include())')/numpy
LOCAL_INCLUDE 	= $(shell python3-config --cflag) -w -I$(shell python3 -c 'import numpy; print(numpy.get_include())')/numpy

OPTIONS	= -std=c++14 -fPIC -Xpreprocessor -fopenmp 

LDFLAGS 	= -shared  -fPIC

DEP	= -MM

EXEC	= $(HOME)/lib/qcm.so

#the common core of the makefiles
include ../qcm_object_list.txt
include ../makefile.mk


