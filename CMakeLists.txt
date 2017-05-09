cmake_minimum_required(VERSION 3.0)

# set project name and languanges that are involved
project(eT Fortran)

# To link with fortran not c++
set(CMAKE_CXX_LINKER_PREFERENCE_PROPAGATES FALSE)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    /usr/local/share/cmake/Modules
    )

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"       ON)

# Including the source files
include(eT_sources)
include(Fortran_flags)
include(FindBLAS)
include(FindLAPACK)

find_package(BLAS)
find_package(LAPACK)

if(LAPACK_FOUND AND BLAS_FOUND)
   set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

set(WORK_MEM_WORDS "64000000")
add_definitions(-DINSTALL_WRKMEM=${WORK_MEM_WORDS})

set(eT_fortran_sources
        ${eT_fortran_sources}
        )

add_executable(eT ${eT_fortran_sources})

target_link_libraries(eT ${lapackblas_libraries})