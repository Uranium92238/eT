# Fortran submodule suport
cmake_minimum_required(VERSION 3.7)

# set project name and languanges that are involved
project(eT Fortran CXX)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    /usr/local/share/cmake/Modules
    )

# set options with standard values 
option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"      	                    ON)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    ON)
option(ENABLE_THREADED_MKL   "Enable OpenMP parallelization"                    ON)

message(${LIBINT_LIB})

# check for LIBINT_LIB environment variable 
# otherwise, set to standard path 
if(DEFINED ENV{LIBINT_LIB})
  set(LIBINT_LIB $ENV{LIBINT_LIB})
else()
  set(LIBINT_LIB /usr/local/libint)
endif()

message(${LIBINT_LIB})

# check for EIGEN_LIB environment variable 
# otherwise, set to standard path 
if(DEFINED ENV{EIGEN_LIB})
  set(EIGEN_LIB $ENV{EIGEN_LIB})
else()
  set(EIGEN_LIB /usr/local/include/eigen3)
endif()

# Including the source files
include(eT_sources)
include(Fortran_flags)
include(eT_CXX_sources)
include(CXX_flags)
include(ConfigTesting)

# math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

set(libint_libraries "${LIBINT_LIB}/2.4.2/lib/libint2.a")

set(eT_fortran_sources
        ${eT_fortran_sources}
        )

set(eT_CXX_sources
        ${eT_CXX_sources}
        )

add_executable(eT ${eT_fortran_sources} ${eT_CXX_sources})

target_link_libraries(eT ${lapackblas_libraries} ${libint_libraries})

set_target_properties(eT PROPERTIES LINKER_LANGUAGE Fortran)

file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})

if(EXISTS ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest)
   file(COPY ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest DESTINATION ${CMAKE_BINARY_DIR}/tests)
endif()

configure_file(${CMAKE_SOURCE_DIR}/eT_launch.in ${CMAKE_SOURCE_DIR}/eT_launch)
file(COPY ${CMAKE_SOURCE_DIR}/eT_launch 
     DESTINATION ${CMAKE_BINARY_DIR} 
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ)
file(REMOVE  ${CMAKE_SOURCE_DIR}/eT_launch)

