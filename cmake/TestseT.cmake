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
# 	Based on the cmake/DaltonTests.cmake file of the public 
#   Dalton program (LGPL v2.1)
#
# 	Copied and modified for eT by Rolf H. Myhre, Feb 2019
#
macro(add_eT_runtest _name _labels)
    add_test(
        ${_name}
        python3 ${PROJECT_BINARY_DIR}/tests/${_name}/test --binary-dir=${PROJECT_BINARY_DIR}  --work-dir=${PROJECT_BINARY_DIR}/tests/${_name} --verbose)
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

# All tests here should contain the label "eT"

# Add a keyword for the length of the test: 
# 
# 	short < 30 seconds
# 	medium > 30 seconds < 120 seconds
# 	long > 120 seconds < 200 seconds
# 	verylong > 200 seconds
# 
# NEVER comment out tests

add_eT_runtest(hf_energy                       "eT;short;hf;sad")
add_eT_runtest(hf_scf_energy                   "eT;short;hf;sad")
add_eT_runtest(mp2_energy                      "eT;short;hf;mp2")
add_eT_runtest(cc2_gs_energy                   "eT;short;cc2;gs")
add_eT_runtest(cc2_gs_energy_batched           "eT;short;cc2;gs;cholesky")
add_eT_runtest(cc2_right_es_energies           "eT;short;cc2;es")
add_eT_runtest(cc2_right_es_energies_diis      "eT;short;cc2;es")
add_eT_runtest(cc2_right_cvs_es_energies       "eT;short;cc2;es;cvs")
add_eT_runtest(cc2_left_es_energies            "eT;short;cc2;es")
add_eT_runtest(cc2_lowmem_gs_energy            "eT;short;lowmem-cc2;gs")
add_eT_runtest(cc2_lowmem_right_es_energies    "eT;short;lowmem-cc2;es")
add_eT_runtest(ccsd_gs_energy                  "eT;short;ccsd;gs")
add_eT_runtest(ccsd_right_es_energies          "eT;short;ccsd;es")
add_eT_runtest(ccsd_right_es_replace_red_space "eT;ccsd;es;red_space")
add_eT_runtest(ccsd_left_es_energies           "eT;short;ccsd;es")
add_eT_runtest(ccsd_right_es_diis              "eT;short;ccsd;es;right;diis")
add_eT_runtest(ccsd_left_es_diis               "eT;short;ccsd;es;left;diis")
add_eT_runtest(ccsd_right_cvs_es_energies      "eT;short;ccsd;es;cvs")
add_eT_runtest(ccsd_dipole                     "eT;short;ccsd;gs;dipole")
add_eT_runtest(ccsd_quadrupole                 "eT;short;ccsd;gs;quadrupole")
add_eT_runtest(cc3_gs_energy                   "eT;short;cc3;gs")
add_eT_runtest(cc3_right_es_energies           "eT;short;cc3;es")
