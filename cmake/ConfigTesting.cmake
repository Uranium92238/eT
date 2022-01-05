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
#   Based on the cmake/ConfigTesting.cmake file of the public 
#   Dalton program (LGPL v2.1)
#
#   Copied and modified for eT by Rolf H. Myhre, Feb 2019
#
# set cdash buildname
set(BUILDNAME
    "${CMAKE_SYSTEM_NAME}-${CMAKE_HOST_SYSTEM_PROCESSOR}-${CMAKE_Fortran_COMPILER_ID}-${CMAKE_CXX_COMPILER_ID}-${BLAS_TYPE}-${CMAKE_BUILD_TYPE}"
    CACHE STRING
    "Name of build on the dashboard"
    )

# set ctest own timeout
if(ENABLE_LARGE_TEST)
   set(DART_TESTING_TIMEOUT
      "18000" # 5 hours 
      CACHE STRING
      "Set timeout in seconds for every single test"
      )
else()
   set(DART_TESTING_TIMEOUT
      "1200"  # 20 minutes
      CACHE STRING
      "Set timeout in seconds for every single test"
      )
endif()

include(TestseT)
include(CTest)
enable_testing()
