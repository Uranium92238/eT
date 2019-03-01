#
#CMake is awful, but everything else is worse
#

set(LIBINT_ROOT "xxxxxx")

message(${LIBINT_ROOT})
if(LIBINT_ROOT)
   set(_LIBINT_ROOT_SEARCH ${LIBINT_ROOT} NO_DEFAULT_PATH)
endif()

#Normal search
set(_LIBINT_NORMAL_SEARCH
    /usr/local/libint 
    /usr/opt/local/libint 
    /home/rolf/progs/libint-2.4.2/
    $ENV{LIBINT_LIB}) 

message("${_LIBINT_NORMAL_SEARCH}")
message("${_LIBINT_ROOT_SEARCH}")

find_path(BLARG
          NAMES libint2.h
          PATHS ${_LIBINT_NORMAL_SEARCH}
          PATH_SUFFIXES include)

message(${BLARG})

foreach(search ${_LIBINT_NORMAL_SEARCH})
   find_path(BLARGI
             NAMES libint2.h
             PATHS ${search}
             PATH_SUFFIXES 2.4.2/include)
endforeach()

message(${BLARGI})

