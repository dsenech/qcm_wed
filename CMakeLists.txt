cmake_minimum_required(VERSION 3.10)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED True)


if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

#set(CMAKE_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")


# set the project name and version
project(qcm VERSION 1.0 DESCRIPTION "Quantum Cluster Methods shared library")
message(STATUS "Host processor : ${CMAKE_SYSTEM_PROCESSOR}")

find_package(Python3 3.7 COMPONENTS Interpreter NumPy REQUIRED)
find_package(PythonLibs REQUIRED)
message(STATUS "python include dir: ${PYTHON_INCLUDE_DIRS}")
message(STATUS "numpy include dir: ${Python3_NumPy_INCLUDE_DIRS}")

find_library(CUBA_LIBRARY
    NAMES libcuba.a
    HINTS "/usr/local/lib" "$ENV{HOME}/lib"
)

#set(BLA_VENDOR Intel10_64lp)
find_package(LAPACK)

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

target_include_directories(qcm PUBLIC src_qcm src_ed src_python src_util PRIVATE ${PYTHON_INCLUDE_DIRS} ${Python3_NumPy_INCLUDE_DIRS}/numpy)

 
# set(CMAKE_INSTALL_LIBDIR $ENV{HOME}/lib)
install(TARGETS qcm
  DESTINATION $ENV{HOME}/lib
  RENAME qcm.so
)
message(STATUS "output file = ${OUTPUT_NAME}")
message(STATUS "source directory = ${CMAKE_SOURCE_DIR}")

# file(COPY ${CMAKE_SOURCE_DIR}/pyqcm DESTINATION $ENV{HOME}/lib/)
file(CREATE_LINK ${CMAKE_SOURCE_DIR}/pyqcm $ENV{HOME}/lib/pyqcm SYMBOLIC)
# install(FILES ${OUTPUT_NAME}
#   DIRECTORY $ENV{HOME}/lib
#   RENAME qcm.so
# )

# add_custom_command(OUTPUT $ENV{HOME}/lib/qcm.so
#   COMMAND ln -s ${OUTPUT_NAME} $ENV{HOME}/lib/qcm.so
#   COMMENT "Generating symbolic link in $ENV{HOME}/lib/"
# )

#configure_file(qcm.pc.in qcm.pc @ONLY)


