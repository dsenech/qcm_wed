#Your favorite compiler goes here.
COMPILER	= c++

#linker; Should be the compiler in most cases.
#A notable exception in making a MPI app using clang to compile, then the linker should be mpicxx
LINKER	= $(COMPILER)

#dynamic link library will vary quite a bit from one platform to the other
PYPREFIX = $(shell python3-config --prefix)
LINK	= -L/usr/local/lib/ -lcuba -framework Accelerate $(shell python3-config --ldflags) -lpython3.9 -lomp
# LINK	= -L/usr/local/lib/ -lcuba -framework Accelerate $(shell python3-config --ldflags) -lpython3.9

#include path, should be left empty for most platform
LOCAL_INCLUDE 	= -I/usr/local/include/  -I${MKLROOT}/include \
	$(shell python3-config --includes)  \
	-I$(shell python3 -c 'import numpy; print(numpy.get_include())')/numpy	

#options and macro for the compilation. Do not put optimisation
#or debug option in there as it will interfere with the debug task and all task defined
# in the common parts of the makefiles
OPTIONS	= -std=c++11 -fPIC -DMKL_ILP64 -Xpreprocessor -fopenmp
# OPTIONS	= -std=c++11 -fPIC -DMKL_ILP64

#flags and search path for the linker.
LDFLAGS 	= -shared  -fPIC

#option for generating dependency files. "-MM" is the correct option for clang, gcc and
# intel's compiler.
DEP	= -MM

#the resulting executable
EXEC	= $(HOME)/lib/qcm.so

#the common core of the makefiles
include ../qcm_object_list.txt
include ../makefile.mk

