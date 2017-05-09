#
# Testing for Fortran compiler, so far only gfortran (GNU)
#  
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    add_definitions(-DVAR_GFORTRAN)
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -ffloat-store -fcray-pointer")
#
#   Testing processor 32-bit or 64-bit
#
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "i386")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m32"
            )
    endif()
#
    if(${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -m64"
            )
    endif()
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -fbacktrace -fcray-pointer -Wuninitialized")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ffast-math -funroll-loops -ftree-vectorize")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
#
#   If ENABLE_OMP is set to ON
#
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fopenmp"
            )
    endif()
#
#   If ENABLE_64BIT_INTEGERS ON
#
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
endif()