#
#
#   eT - a coupled cluster program
#   Copyright (C) 2016-2022 the authors of eT
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
#  Based on the cmake/TestseT.cmake file of the public
#   Dalton program (LGPL v2.1)
#
#  Copied and modified for eT by Rolf H. Myhre, Feb 2019
#

macro(add_eT_unittest _name _labels)

   set (filename unit_${_name})

   add_pfunit_ctest (${filename}
      TEST_SOURCES ${CMAKE_BINARY_DIR}/${filename}.pf
      LINK_LIBRARIES eT_library ${lapackblas_libraries} ${Libint2_LIBRARY} ${PCMSolver_LIBRARY} ${EXTRA_LINKER_FLAGS} get_compilation_info
   )

   # add "unit-test" as a label for all unit tests
   string(CONCAT labels "${_labels}" "; unit-test")

   set_tests_properties(${filename} PROPERTIES LABELS "${labels}")

endmacro()

add_eT_unittest(warning_suppressor "dummy; warning-suppression")
#
add_eT_unittest(angular_momentum "angular-momentum")
add_eT_unittest(cart_angular_momentum_7 "angular-momentum")
#
add_eT_unittest(tools "tools")
add_eT_unittest(memory "memory")
