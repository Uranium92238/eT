# Try to find PCMSolver
#
# If successful, will define the required variables
# PCMSolver_FOUND - true if PCMSolver was found
# PCMSolver_INCLUDE_DIR - Include dir PCMSolver
# PCMSolver_LIBRARY - Location of PCMSolver Library
#
# Copyright (c) 2019 Marco Scavino and Tommaso Giovannini
# Distributed under the GNU Lesser General Public License.

#Where to look for PCMSolver

set(_PCMSolver_NORMAL_SEARCH
  "/usr/local"
  "/usr/local/PCMSolver"
  "/usr/opt/local/PCMSolver")

set(PCMSolver_ENV $ENV{PCMSolver_ROOT})

if(PCMSolver_ENV)
  set(PCMSolver_ROOT ${PCMSolver_ENV})
endif()

if(PCMSolver_ROOT)
   message("-- PCMSolver will be searched for based on PCMSolver_ROOT=" ${PCMSolver_ROOT})
else()
   message("-- PCMSolver will be searched for based on default path ${_PCMSolver_NORMAL_SEARCH}")
endif()

# Use the find_package macro to search in some possible paths

find_package(PCMSolver
 CONFIG
 PATHS ${PCMSolver_ROOT} ${_PCMSolver_NORMAL_SEARCH}
 NO_DEFAULT_PATH
 COMPONENTS exe
 QUIET)

set(PCMSolver_INCLUDE_DIR ${PCMSolver_INCLUDE_DIR}    CACHE STRING "PCMSolver include path")
set(PCMSolver_PYMOD ${PCMSolver_PYMOD}                CACHE STRING "PCMSolver python module path")
set(PCMSolver_EXECUTABLE ${PCMSolver_EXECUTABLE}      CACHE STRING "PCMSolver binary path")

# Let the module handle the variables
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PCMSolver
   FOUND_VAR PCMSolver_FOUND
   REQUIRED_VARS PCMSolver_INCLUDE_DIR
                 PCMSolver_LIBRARY
                 PCMSolver_PYMOD
                 PCMSolver_EXECUTABLE)
