# 
#   :: Tests for C++ compiler and setting of compiler flags accordingly
# 
#   1. Intel icpc 
# 
if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
# 
#   Set standard flags 
#  
    set(CMAKE_CXX_FLAGS "-std=c++11 -xHost -O3")
    if(DEVELOPMENT_CODE)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
    else()
        # suppress warnings in exported code
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -w")
    endif()
    set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
# 
#   Flags needed for Libint package
# 
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -fexceptions -I${Libint2_INCLUDE_DIR} -I${Libint2_H_DIR} -I${EIGEN3_INCLUDE_DIR}"
        )
# 
#   MKL flags 
# 
    if(DEFINED MKL_FLAG)
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${MKL_FLAG}")
    endif()
# 
#   Enable openmp if requested (default) 
# 
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -qopenmp -parallel"
            )
    endif()
endif ()
# 
#   2. GCC g++
# 
if(CMAKE_CXX_COMPILER_ID MATCHES GNU)
# 
#   Set standard flags 
# 
    set(CMAKE_CXX_FLAGS "-O3 -std=c++11")
# 
#   Testing processor 32-bit or 64-bit
# 
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -m32"
            )
    endif()
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -m64"
            )
    endif()
# 
#   Flags needed for Libint package
# 
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -fexceptions -I${Libint2_INCLUDE_DIR} -I${Libint2_H_DIR} -I${EIGEN3_INCLUDE_DIR}"
        )
# 
#   Enable openmp if requested (default) 
#  
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -fopenmp"
            )
    endif()
endif()

# Needed on Dirac: -I/usr/local/include/c++/8.2.0 -gcc-version=8.2
if(DEFINED EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()
