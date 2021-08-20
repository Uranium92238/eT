add_library(eT_library ${eT_fortran_sources} ${eT_CXX_sources})

# Supress ranlib no symbols warning for modules with no subroutines
if(APPLE)
    SET(CMAKE_C_ARCHIVE_CREATE   "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_CXX_ARCHIVE_CREATE "<CMAKE_AR> Scr <TARGET> <LINK_FLAGS> <OBJECTS>")
    SET(CMAKE_C_ARCHIVE_FINISH   "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
    SET(CMAKE_CXX_ARCHIVE_FINISH "<CMAKE_RANLIB> -no_warning_for_no_symbols -c <TARGET>")
endif(APPLE)

# Link library to Lapack, Libint, and PCMSolver 
target_link_libraries(eT_library ${lapackblas_libraries} ${Libint2_LIBRARY} 
          ${PCMSolver_LIBRARY} ${EXTRA_LINKER_FLAGS} get_compilation_info)

target_include_directories(eT_library PUBLIC ${CMAKE_CURRENT_BINARY_DIR})
