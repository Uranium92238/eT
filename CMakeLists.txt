#
#
#   eT - a coupled cluster program
#   Copyright (C) 2016-2019 the authors of eT
#
#   eT is free software: you can redistribute it and/or modify
#   it under the terms of the GNU General Public License as published by
#   the Free Software Foundation, either version 3 of the License, or
#   (at your option) any later version.
#
#   eT is distributed in the hope that it will be useful,
#   but WITHOUT ANY WARRANTY; without even the implied warranty of
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#   GNU General Public License for more details.
#
#   You should have received a copy of the GNU General Public License
#   along with this program. If not, see <https://www.gnu.org/licenses/>.
#
#
# Set minimum to guarantee Fortran submodule suport
cmake_minimum_required(VERSION 3.7)

# Project name and programming languanges
project(eT Fortran CXX)

# Paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    /usr/local/share/cmake/Modules
    )

# Options and their standard values 
option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"       			              ON)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    ON)
option(ENABLE_THREADED_MKL   "Enable OpenMP parallelization"                    ON)

# Check for LIBINT_LIB environment variable 
if(DEFINED ENV{LIBINT_LIB})
  set(LIBINT_LIB $ENV{LIBINT_LIB})
else()
  set(LIBINT_LIB /usr/local/libint)
endif()

# Check for EIGEN_LIB environment variable 
if(DEFINED ENV{EIGEN_LIB})
  set(EIGEN_LIB $ENV{EIGEN_LIB})
else()
  set(EIGEN_LIB /usr/local/include/eigen3)
endif()

# Include source files
include(eT_sources)
include(Fortran_flags)
include(eT_CXX_sources)
include(CXX_flags)
include(ConfigTesting)

# Math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

# Set lapack and blas libraries 
if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

# Set Libint library 
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

