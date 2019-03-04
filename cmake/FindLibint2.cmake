#
#

#Normal search
set(_LIBINT_NORMAL_SEARCH
    /usr/local/libint 
    /usr/opt/local/libint)

set(LIBINT2_ENV $ENV{LIBINT2_ROOT})

set(LIBINT2_GLOB_DIR)
if(LIBINT2_ROOT)
   message("-- Libint2 will be searched for based on LIBINT2_ROOT=" ${LIBINT2_ROOT})
   file(GLOB_RECURSE _GLOBFILE ${LIBINT2_ROOT}/*/libint2.h)
   if(_GLOBFILE)
      get_filename_component(_GLOBDIR ${_GLOBFILE} PATH)
      list(APPEND LIBINT2_GLOB_DIR ${_GLOBDIR})
   endif()
elseif(LIBINT2_ENV)
   message("-- Libint2 will be searched for based on LIBINT2_ROOT=" ${LIBINT2_ENV})
   file(GLOB_RECURSE _GLOBFILE ${LIBINT2_ENV}/*/libint2.h)
   if(_GLOBFILE)
      get_filename_component(_GLOBDIR ${_GLOBFILE} PATH)
      list(APPEND LIBINT2_GLOB_DIR ${_GLOBDIR})
   endif()
else()
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
                
   file(GLOB_RECURSE _BASIS_DIR ${_BASE_PATH}/share/libint/*/basis/3-21g.g94)
   get_filename_component(Libint2_BASIS_DIR ${_BASIS_DIR} PATH)

endif()


include(FindPackageHandleStandardArgs)
FIND_PACKAGE_HANDLE_STANDARD_ARGS(Libint2 
      FOUND_VAR Libint2_FOUND
      REQUIRED_VARS Libint2_INCLUDE_DIR
                    Libint2_H_DIR
                    Libint2_LIBRARY
                    Libint2_BASIS_DIR)

