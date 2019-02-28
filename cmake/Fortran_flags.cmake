# 
#   :: Tests for Fortran compiler and setting of compiler flags accordingly
# 
#   1. Intel ifort
# 
if(CMAKE_Fortran_COMPILER_ID MATCHES Intel)
# 
#   Set standard flags 
# 
    set(CMAKE_Fortran_FLAGS "-fpp -O3 -warn all -xHost")
# 
#   Enable 64 bit flag if requested (default)
# 
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i8"
            )
    else()
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -i4"
            )        
    endif()
# 
#   Enable openmp if requested (default) 
# 
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -qopenmp -parallel"
            )
    endif()
# 
endif()
# 
#   2. GNU gfortran
# 
if(CMAKE_Fortran_COMPILER_ID MATCHES GNU)
# 
#   Set standard flags 
# 
    set(CMAKE_Fortran_FLAGS "-std=f2008 -O3 -Wall -march=native")
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
# 
#   Enable 64 bit flag if requested (default)
# 
    if(ENABLE_64BIT_INTEGERS)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fdefault-integer-8"
            )
    endif()
# 
#   Enable openmp if requested (default) 
# 
    if(ENABLE_OMP)
        set(CMAKE_Fortran_FLAGS
            "${CMAKE_Fortran_FLAGS} -fopenmp"
            )
    endif()
# 
endif()
