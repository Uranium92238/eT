#
# Testing for Fortran compiler and setting compiler flags
#
## GNU ##  
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
    add_definitions(-DVAR_GFORTRAN)
    set(CMAKE_Fortran_FLAGS         "-DVAR_GFORTRAN -ffloat-store -fcray-pointer -finit-local-zero")
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
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg -fbounds-check -o0")
#
#   If ENABLE_64BIT_INTEGERS ON
#
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
endif()
#
## Intel iFort ##
#
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
    add_definitions(-DVAR_IFORT)
    set(CMAKE_Fortran_FLAGS         "-fpp -assume byterecl -DVAR_IFORT")
    set(CMAKE_Fortran_FLAGS_DEBUG   "-O0 -g -traceback")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O3 -ip -diag-disable 8290 -diag-disable 8291")
    set(CMAKE_Fortran_FLAGS_PROFILE "${CMAKE_Fortran_FLAGS_RELEASE} -g -pg")
#
#   If ENABLE_64BIT_INTEGERS ON
#
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    endif()
#
    if(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
        message("--Switch off warnings due to incompatibility XCode 4 and Intel 11 on OsX 10.6")
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -Qoption,ld,-w"
            )
    endif()
    set(reorder_definitions " --nocollapse ${reorder_definitions}")
endif()
