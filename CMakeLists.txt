cmake_minimum_required(VERSION 3.15)

# set the project name and version
project(qcm LANGUAGES CXX VERSION 1.0 DESCRIPTION "Quantum Cluster Methods shared library")

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_POSITION_INDEPENDENT_CODE ON)

message(STATUS "Host processor : ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "Architecture: ${ARCHITECTURE}")

find_package(Python3 3.7 COMPONENTS Interpreter Development NumPy REQUIRED)
message(STATUS "python executable: ${Python3_EXECUTABLE}")
find_package(PythonLibs REQUIRED)
message(STATUS "python include dir: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "python libs: ${PYTHON_LIBRARIES}")
set(NUMPY_INCLUDE_DIR ${Python3_NumPy_INCLUDE_DIRS})

if(SKBUILD)
message(STATUS "this is built with SciKit-build")
else()
execute_process (COMMAND ${Python3_EXECUTABLE} -c "import numpy; print(numpy.get_include())" OUTPUT_VARIABLE NUMPY_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
endif()

message(STATUS "numpy dir : ${NUMPY_INCLUDE_DIR}")

include(FindLAPACK)

find_library(CUBA_LIBRARY
    NAMES libcuba.a
    HINTS "/usr/local/lib" "$ENV{HOME}/lib"
)

add_library(qcm SHARED
src_qcm/CPT.cpp
src_qcm/Chern.cpp
src_qcm/QCM.cpp
src_qcm/average.cpp
src_qcm/basis3D.cpp
src_qcm/lattice3D.cpp
src_qcm/lattice_model.cpp
src_qcm/lattice_model_instance.cpp
src_qcm/lattice_operator.cpp
src_qcm/parameter_set.cpp
src_qcm/profile.cpp
src_ed/ED_basis.cpp
src_ed/Heisenberg_operator.cpp
src_ed/Hund_operator.cpp
src_ed/Lanczos.cpp
src_ed/binary_state.cpp
src_ed/continued_fraction.cpp
src_ed/continued_fraction_set.cpp
src_ed/destruction_operator.cpp
src_ed/model.cpp
src_ed/model_instance.cpp
src_ed/model_instance_base.cpp
src_ed/qcm_ED.cpp
src_ed/sector.cpp
src_ed/symmetry_group.cpp
src_python/common_Py.cpp
src_python/qcm_lib.cpp
src_util/ED_global_parameter.cpp
src_util/global_parameter.cpp
src_util/integrate.cpp
src_util/matrix.cpp
src_util/parser.cpp
src_util/vector_num.cpp
)

set_target_properties(qcm PROPERTIES VERSION ${PROJECT_VERSION} PREFIX "" SUFFIX ".so")

target_link_libraries(qcm ${PYTHON_LIBRARIES} ${CUBA_LIBRARY} ${LAPACK_LIBRARIES})

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
  target_link_libraries(qcm OpenMP::OpenMP_CXX)
endif()

target_include_directories(qcm PUBLIC src_qcm src_ed src_python src_util ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR}/numpy)

# install(TARGETS qcm DESTINATION ${PYTHON_INCLUDE_DIRS}/site-packages)
install(TARGETS qcm DESTINATION pyqcm)

