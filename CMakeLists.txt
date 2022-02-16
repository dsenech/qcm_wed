cmake_minimum_required(VERSION 3.22)


# specify the C++ standard
# set(CMAKE_C_COMPILER_WORKS 1)
# set(CMAKE_CXX_COMPILER_WORKS 1)
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)

if(APPLE)
    set(CMAKE_THREAD_LIBS_INIT "-lpthread")
    set(CMAKE_HAVE_THREADS_LIBRARY 1)
    set(CMAKE_USE_WIN32_THREADS_INIT 0)
    set(CMAKE_USE_PTHREADS_INIT 1)
    set(THREADS_PREFER_PTHREAD_FLAG ON)
    set(LAPACK_LIBRARIES "-framework accelerate")
else()
    find_package(LAPACK REQUIRED)
endif()

# if(NOT CMAKE_BUILD_TYPE)
#   set(CMAKE_BUILD_TYPE Release)
# endif()

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

message(STATUS "Host processor : ${CMAKE_HOST_SYSTEM_PROCESSOR}")
message(STATUS "Architecture: ${ARCHITECTURE}")

# set the project name and version
project(qcm LANGUAGES CXX VERSION 1.0 DESCRIPTION "Quantum Cluster Methods shared library")

find_package(Python3 3.7 COMPONENTS Interpreter NumPy REQUIRED)
find_package(PythonLibs REQUIRED)
message(STATUS "python include dir: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "numpy include dir: ${Python3_NumPy_INCLUDE_DIRS}")

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

target_include_directories(qcm PUBLIC src_qcm src_ed src_python src_util ${PYTHON_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS}/numpy)

if(SKBUILD)
  message(STATUS "The project is built using scikit-build")
  install(TARGETS qcm DESTINATION lib/python3.9/site-packages)
else()
  install(TARGETS qcm
  DESTINATION $ENV{HOME}/lib
)
endif()

# message(STATUS "output file = ${OUTPUT_NAME}")
# message(STATUS "source directory = ${CMAKE_SOURCE_DIR}")

# file(COPY ${CMAKE_SOURCE_DIR}/pyqcm DESTINATION $ENV{HOME}/lib/)
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/pyqcm $ENV{HOME}/lib/pyqcm SYMBOLIC)


#configure_file(qcm.pc.in qcm.pc @ONLY)


