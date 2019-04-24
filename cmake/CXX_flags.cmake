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
#   C++ flags for icpc and g++
#
#   1. Intel icpc 
# 
if (CMAKE_CXX_COMPILER_ID MATCHES Intel)
# 
#   Set standard flags 
#  
    set(CMAKE_CXX_FLAGS "-std=c++11 -xHost -O3")
# 
    set (CMAKE_CXX_LINK_FLAGS "${CMAKE_CXX_LINK_FLAGS} -shared-intel")
# 
#   Flags needed for Libint package
# 
    set(CMAKE_CXX_FLAGS
        "${CMAKE_CXX_FLAGS} -fexceptions"
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
        "${CMAKE_CXX_FLAGS} -fexceptions"
        )
# 
#   Enable openmp if requested (default) 
#  
    if(ENABLE_OMP)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -fopenmp"
            )
    endif()
# 
#   Enable profiling 
# 
    if(ENABLE_PROFILING)
        set(CMAKE_CXX_FLAGS
            "${CMAKE_CXX_FLAGS} -pg"
            )
    endif()
endif()

# Needed on Dirac: -I/usr/local/include/c++/8.2.0 -gcc-version=8.2
if(DEFINED EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif()
