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

#Eirik F. Kjønstad, Nov 2019: assume .dylib if mac, .so otherwise

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
 QUIET)

if(PCMSolver_LIBRARY)
   # If found PCMSolver, get the library path
   get_filename_component(PCMSolver_DIR ${PCMSolver_LIBRARY} DIRECTORY)

   set(PCMSolver_LIBRARY ${PCMSolver_LIBRARY})
   set(PCMSolver_INCLUDE_DIR ${PCMSolver_INCLUDE_DIRS})

else()

   #prepare for external project

   include(GNUInstallDirs)
   include(ExternalProject)

   set(PCMSolver_library ${CMAKE_CURRENT_SOURCE_DIR}/external/PCMSolver)

   set(PCMSolver_install ${PCMSolver_library}/install)

   message("-- Setting up the external project for PCMSolver in: ${PCMSolver_library}")

   #External project to construct
   ExternalProject_add(project_PCMSolver
      # git repository url
      GIT_REPOSITORY https://github.com/eirik-kjonstad/pcmsolver.git
      GIT_TAG patched_release_1.2
      # root directory of the project
      PREFIX ${PCMSolver_library}
      # source directory in which git clones the repository
      SOURCE_DIR ${PCMSolver_library}/repository
      # arguments pass to cmake to configure the library project
      # CMAKE_INSTALL_PREFIX set for the project the install directory
      CMAKE_ARGS -DCMAKE_INSTALL_PREFIX:PATH=${PCMSolver_install}
                 -DCMAKE_C_COMPILER=${CMAKE_C_COMPILER}
                 -DCMAKE_CXX_COMPILER=${CMAKE_CXX_COMPILER}
                 -DCMAKE_Fortran_COMPILER=${CMAKE_Fortran_COMPILER}
      # build directory
      BINARY_DIR ${PCMSolver_library}/build
      #  temporary directory for the project
      STAMP_DIR ${PCMSolver_library}/stamp
      TMP_DIR ${PCMSolver_library}/tmp
      # command used to compile the code in the build directory
      BUILD_COMMAND $(MAKE)
      # install command
      INSTALL_COMMAND $(MAKE) install
      UPDATE_COMMAND ""
   )

   # add the two libraries
   add_library(PCMSolver SHARED IMPORTED)

   # link the IMPORTED library to its files
   if(${CMAKE_SYSTEM_NAME} STREQUAL "Darwin")
     set_property(TARGET PCMSolver
      PROPERTY IMPORTED_LOCATION ${PCMSolver_install}/${CMAKE_INSTALL_LIBDIR}/libpcm.dylib)
   else()
     set_property(TARGET PCMSolver
      PROPERTY IMPORTED_LOCATION ${PCMSolver_install}/${CMAKE_INSTALL_LIBDIR}/libpcm.so)
   endif()

   # mandatory. set the library dependencies
   add_dependencies(PCMSolver project_PCMSolver)

   #set DFT variables
   set(PCMSolver_INCLUDE_DIR ${PCMSolver_install}/include)
   set(PCMSolver_LIBRARY PCMSolver)

endif()

# Let the module handle the variables
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(PCMSolver
   FOUND_VAR PCMSolver_FOUND
   REQUIRED_VARS PCMSolver_INCLUDE_DIR
                 PCMSolver_LIBRARY)

