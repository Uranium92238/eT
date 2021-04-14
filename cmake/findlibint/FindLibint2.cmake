#  Try to find Libint2
#
#  If successful, will define the required variables
#  Libint2_FOUND - true if Libint2 was found
#  Libint2_INCLUDE_DIR - Containing libint2.h
#  Libint2_H_DIR - Directory with Libint2 header files
#  Libint2_LIBRARY - Location of libint2.a
#
#  And the extra variable
#  Libint2_BASIS_DIR
#
#  Copyright (c) 2019 Rolf Heilemann Myhre
#  Distributed under the GNU Lesser General Public License.


#Where to look for Libint2
set(_LIBINT_NORMAL_SEARCH
    /usr/local/libint
    /usr/opt/local/libint)

#Save location if set in environment
set(LIBINT2_ENV $ENV{LIBINT2_ROOT})

set(LIBINT2_GLOB_DIR)

# Must use glob instead of find_path because path depends on version
# First see if provided as command line argument
if(LIBINT2_ROOT)
   message(STATUS "Libint2 will be searched for based on LIBINT2_ROOT=" ${LIBINT2_ROOT})
   file(GLOB_RECURSE _GLOBFILE ${LIBINT2_ROOT}/*/libint2.h)
   if(_GLOBFILE)
      get_filename_component(_GLOBDIR ${_GLOBFILE} PATH)
      list(APPEND LIBINT2_GLOB_DIR ${_GLOBDIR})
   endif()

# Else look in environment path
elseif(LIBINT2_ENV)
   message(STATUS "Libint2 will be searched for based on LIBINT2_ROOT=" ${LIBINT2_ENV})
   file(GLOB_RECURSE _GLOBFILE ${LIBINT2_ENV}/*/libint2.h)
   if(_GLOBFILE)
      get_filename_component(_GLOBDIR ${_GLOBFILE} PATH)
      list(APPEND LIBINT2_GLOB_DIR ${_GLOBDIR})
   endif()

# Else look in standard paths
else()
   message(STATUS "Libint2 will be searched for based on default path ${_LIBINT_NORMAL_SEARCH}")
   foreach(search ${_LIBINT_NORMAL_SEARCH})
      file(GLOB_RECURSE _GLOBFILE ${search}/*/libint2.h)
      if(_GLOBFILE)
         list(REVERSE _GLOBFILE)
         foreach(hit ${_GLOBFILE})
            get_filename_component(_GLOBDIR ${hit} PATH)
            list(APPEND LIBINT2_GLOB_DIR ${_GLOBDIR})
         endforeach()
      endif()
   endforeach()
endif()


# Store the actual variables
if(LIBINT2_GLOB_DIR)
   find_path(Libint2_INCLUDE_DIR
             NAMES libint2.h
             PATHS ${LIBINT2_GLOB_DIR} NO_DEFAULT_PATH)

   find_path(Libint2_H_DIR
             NAMES config.h
             PATHS ${Libint2_INCLUDE_DIR} NO_DEFAULT_PATH
             PATH_SUFFIXES libint2)

   get_filename_component(_BASE_PATH ${Libint2_INCLUDE_DIR} PATH)
   find_library(Libint2_LIBRARY
                NAMES libint2.a
                PATHS ${_BASE_PATH} NO_DEFAULT_PATH
                PATH_SUFFIXES lib lib/.libs)

   #Try to find the Libint data path
   file(GLOB_RECURSE _BASIS_DIR ${_BASE_PATH}/*/3-21g.g94)
   if(_BASIS_DIR)
      get_filename_component(Libint2_BASIS_DIR ${_BASIS_DIR} PATH)
   endif()

endif()

# Let the module handle the variables
include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Libint2
      FOUND_VAR Libint2_FOUND
      REQUIRED_VARS Libint2_INCLUDE_DIR
                    Libint2_H_DIR
                    Libint2_LIBRARY)

