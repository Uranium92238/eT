cmake_minimum_required(VERSION 3.6)

# set project name and languanges that are involved
project(eT Fortran CXX)

# these are paths that CMake will search for cmake
# module files that end with .cmake
set(CMAKE_MODULE_PATH
    ${CMAKE_SOURCE_DIR}/cmake
    ${CMAKE_SOURCE_DIR}/cmake/math
    /usr/local/share/cmake/Modules
    )

option(ENABLE_64BIT_INTEGERS "Enable 64-bit integers"       			OFF)
option(ENABLE_AUTO_BLAS      "Enable CMake to autodetect BLAS"                  ON)
option(ENABLE_AUTO_LAPACK    "Enable CMake to autodetect LAPACK"                ON)
option(ENABLE_OMP            "Enable OpenMP parallelization"                    ON)

# Including the source files
include(eT_sources)
include(Fortran_flags)
include(eT_CXX_sources)
include(CXX_flags)

# math detection
set(BLAS_LANG "Fortran")
set(LAPACK_LANG "Fortran")
set(MKL_COMPILER_BINDINGS ${CMAKE_Fortran_COMPILER_ID})
include(ConfigMath)

if(LAPACK_FOUND AND BLAS_FOUND)
  set(lapackblas_libraries ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES})
endif()

set(libint_libraries "/usr/local/libint/2.4.2/lib/libint2.a")

set(eT_fortran_sources
        ${eT_fortran_sources}
        )

set(eT_CXX_sources
        ${eT_CXX_sources}
        )

add_executable(eT ${eT_fortran_sources} ${eT_CXX_sources})

target_link_libraries(eT ${lapackblas_libraries} ${libint_libraries})

set_target_properties(eT PROPERTIES LINKER_LANGUAGE Fortran)
