#
#
#   eT - a coupled cluster program
#   Copyright (C) 2016-2022 the authors of eT
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
project(eT Fortran CXX C)

# Paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    ${CMAKE_SOURCE_DIR}/cmake/findlibint
    ${CMAKE_SOURCE_DIR}/cmake/findPCMSolver
    /usr/local/share/cmake/Modules
    )

# Set options with standard values

option(ENABLE_64BIT_INTEGERS  "Enable 64-bit integers"                                   ON)
option(ENABLE_AUTO_BLAS       "Enable CMake to autodetect BLAS"                          ON)
option(ENABLE_AUTO_LAPACK     "Enable CMake to autodetect LAPACK"                        ON)
option(ENABLE_OMP             "Enable OpenMP parallelization"                            ON)
option(ENABLE_THREADED_MKL    "Enable OpenMP parallelization in MKL"                     ON)
option(ENABLE_RUNTIME_CHECKS  "Enable all runtime checks (for debug purposes)"           OFF)
option(ENABLE_FORCED_BATCHING "Enable batching in every batch loop (for debug purposes)" OFF)
option(ENABLE_PCMSOLVER       "Enable PCMSolver library"                                 OFF)
option(ENABLE_PFUNIT          "Enable unit tests with pFUnit"                            OFF)

# Find packages needed for libint
find_package(Libint2 REQUIRED)
find_package(Eigen3)
find_package(Boost REQUIRED)

include_directories(${Boost_INCLUDE_DIRS})
include_directories(${EIGEN3_INCLUDE_DIR})
include_directories(${Libint2_H_DIR})
include_directories(${Libint2_INCLUDE_DIR})

message(STATUS "Boost location: ${Boost_INCLUDE_DIRS}")

# Including the source files
include(eT_sources)
include(Fortran_flags)
include(eT_CXX_sources)
include(CXX_flags)
include(ConfigTesting)

if(NOT ENABLE_64BIT_INTEGERS)
   add_definitions(-DHAS_32BIT_INTEGERS)
   message(STATUS "32bit integers enabled")
endif()

# Math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

# Set lapack and blas libraries
if(LAPACK_FOUND AND BLAS_FOUND)
    set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

#Handle PCMSolver
if(ENABLE_PCMSOLVER)
   add_definitions(-DHAS_PCMSOLVER)
   message(STATUS "Polarizable Continuum Model enabled")
   find_package(PCMSolver REQUIRED)
   include_directories(${PCMSolver_INCLUDE_DIR})
   set(eT_fortran_sources
   ${PCMSolver_INCLUDE_DIR}/PCMSolver/pcmsolver.f90
   ${eT_fortran_sources}
   )
endif()

# Make library for compilation information
add_library(get_compilation_info ${CMAKE_BINARY_DIR}/get_compilation_info.F90)

# Make eT library
include(eT_library)

# Make eT executable from eT_program.F90 and link the executable
# to eT library source as well as other libraries
add_executable(eT ${CMAKE_SOURCE_DIR}/src/eT_program.F90)

target_link_libraries(eT ${lapackblas_libraries} ${Libint2_LIBRARY}
                         ${PCMSolver_LIBRARY} ${EXTRA_LINKER_FLAGS}
                         get_compilation_info eT_library)

set_target_properties(eT PROPERTIES LINKER_LANGUAGE Fortran)

# Set compiler definition if forced batching is enabled
if(ENABLE_FORCED_BATCHING)
    message(STATUS "Forced batching enabled")
    target_compile_definitions(eT PUBLIC _FORCED_BATCHING)
    target_compile_definitions(eT_library PUBLIC _FORCED_BATCHING)
endif()

if(ENABLE_RUNTIME_CHECKS)
    message(STATUS "Runtime checks enabled")
endif()

file(COPY ${CMAKE_SOURCE_DIR}/tests DESTINATION ${CMAKE_BINARY_DIR})

if(EXISTS ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest)
   file(COPY ${CMAKE_SOURCE_DIR}/submodules/runtest/runtest DESTINATION ${CMAKE_BINARY_DIR}/tests)
endif()

#Configure eT_launch.in, set permissions toexecutable and copy into build
configure_file(${CMAKE_SOURCE_DIR}/eT_launch.in ${CMAKE_SOURCE_DIR}/eT_launch.py)
file(COPY ${CMAKE_SOURCE_DIR}/eT_launch.py
     DESTINATION ${CMAKE_BINARY_DIR}
     FILE_PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ WORLD_READ)
file(REMOVE  ${CMAKE_SOURCE_DIR}/eT_launch.py)

# Configure get_complilation_info library
include(CompilationInfo)
configure_file(${CMAKE_SOURCE_DIR}/cmake/get_compilation_info.in ${CMAKE_BINARY_DIR}/get_compilation_info.F90)

if(ENABLE_PFUNIT)

  # pFUnit - Fortran unit testing package

  cmake_minimum_required(VERSION 3.12) # Due to GFTL (used by pFUnit),
                                       # which sets <prefix>_ROOT CMake variables,
                                       # assuming >= 3.12 policy (CMP0074)

  find_package(PFUNIT)

  if(PFUNIT_FOUND)
    enable_testing()
    add_definitions(-DINCLUDE_INTENTIONALLY_BROKEN)
    include(UnitTestseT)
    message(STATUS "Unit testing requested, and pFUnit successfully located!")
  else()
    message(FATAL_ERROR "Unit testing requested, but could not find pFUnit!")
  endif()

endif()
