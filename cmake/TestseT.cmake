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
add_eT_runtest(hf_geoopt                       "eT;short;hf;bfgs;gradient")
add_eT_runtest(hf_scf_energy                   "eT;short;hf;sad")
add_eT_runtest(uhf_energy                      "eT;short;uhf")
add_eT_runtest(hf_energy_qmmmnopol             "eT;short;hf;sad;qmmm;mmnopol")
add_eT_runtest(hf_energy_qmfq                  "eT;short;hf;sad;qmmm;fq")
add_eT_runtest(hf_scf_energy_qmmmnopol         "eT;short;hf;sad;qmmm;mmnopol")
add_eT_runtest(hf_scf_energy_qmfq              "eT;short;hf;sad;qmmm;fq")
add_eT_runtest(uhf_energy_qmmmnopol            "eT;short;uhf;sad;qmmm;mmnopol")
add_eT_runtest(uhf_energy_qmfq                 "eT;short;uhf;sad;qmmm;fq")
add_eT_runtest(mp2_energy                      "eT;short;hf;mp2")
add_eT_runtest(mp2_energy_qmmmnopol            "eT;short;hf;mp2;qmmm;mmnopol")
add_eT_runtest(mp2_energy_qmfq                 "eT;short;hf;mp2;qmmm;fq")
add_eT_runtest(ccs_frozen_core                 "eT;short;ccs;gs;es;frozen core")
add_eT_runtest(ccs_oscillator_strength_eom     "eT;short;ccs;es;eom")
add_eT_runtest(cc2_gs_energy                   "eT;short;cc2;gs")
add_eT_runtest(cc2_gs_energy_batched_cd        "eT;short;cc2;gs;pcd")
add_eT_runtest(cc2_dipole                      "eT;short;cc2;gs;dipole")
add_eT_runtest(cc2_frozen_core                 "eT;short;cc2;gs;es;frozen core")
add_eT_runtest(cc2_right_es_energies           "eT;short;cc2;es;right")
add_eT_runtest(cc2_right_es_energies_diis      "eT;short;cc2;es;right;diis")
add_eT_runtest(cc2_right_cvs_es_energies       "eT;short;cc2;es;right;cvs")
add_eT_runtest(cc2_right_cvs_es_energies_diis  "eT;short;cc2;es;right;cvs;diis")
add_eT_runtest(cc2_right_ip_energies           "eT;short;cc2;es;right")
add_eT_runtest(cc2_right_ip_diis               "eT;short;cc2;es;right;diis")
add_eT_runtest(cc2_left_es_energies            "eT;short;cc2;es;left")
add_eT_runtest(cc2_left_cvs_es_energies        "eT;short;cc2;es;left;cvs")
add_eT_runtest(cc2_left_ip_energies            "eT;short;cc2;es;left")
add_eT_runtest(cc2_left_ip_diis                "eT;short;cc2;es;left;diis")
add_eT_runtest(cc2_oscillator_strength_eom     "eT;short;cc2;eom")
add_eT_runtest(cc2_lowmem_gs_energy            "eT;short;lowmem-cc2;gs")
add_eT_runtest(cc2_lowmem_frozen_core          "eT;short;lowmem-cc2;es;gs;frozen core")
add_eT_runtest(cc2_lowmem_right_es_energies    "eT;short;lowmem-cc2;es;right")
add_eT_runtest(cc2_lowmem_left_es_energies     "eT;short;lowmem-cc2;es")
add_eT_runtest(ccsd_gs_energy                  "eT;short;ccsd;gs")
add_eT_runtest(ccsd_gs_energy_qmmmnopol        "eT;short;ccsd;gs;qmmm;mmnopol")
add_eT_runtest(ccsd_gs_energy_qmfq             "eT;short;ccsd;gs;qmmm;fq")
add_eT_runtest(ccsd_gs_energy_1c_cd            "eT;short;ccsd;gs;1c-cd")
add_eT_runtest(ccsd_gs_energy_NR               "eT;short;ccsd;gs")
add_eT_runtest(ccsd_dipole                     "eT;short;ccsd;gs;dipole")
add_eT_runtest(ccsd_quadrupole                 "eT;short;ccsd;gs;quadrupole")
add_eT_runtest(ccsd_frozen_core                "eT;short;ccsd;es;gs;frozen core")
add_eT_runtest(ccsd_right_es_energies          "eT;short;ccsd;es;right")
add_eT_runtest(ccsd_right_es_replace_red_space "eT;short;ccsd;es;right;red_space")
add_eT_runtest(ccsd_right_es_diis              "eT;short;ccsd;es;right;diis")
add_eT_runtest(ccsd_right_cvs_es_energies      "eT;short;ccsd;es;right;cvs")
add_eT_runtest(ccsd_right_cvs_es_energies_diis "eT;short;ccsd;es;right;cvs;diis")
add_eT_runtest(ccsd_right_ip_energies          "eT;short;ccsd;ip;right")
add_eT_runtest(ccsd_right_ip_diis              "eT;short;ccsd;ip;right;diis")
add_eT_runtest(ccsd_left_ip_energies           "eT;short;ccsd;ip;left")
add_eT_runtest(ccsd_left_ip_diis               "eT;short;ccsd;ip;left;diis")
add_eT_runtest(ccsd_left_es_energies           "eT;short;ccsd;es;left")
add_eT_runtest(ccsd_left_es_diis               "eT;short;ccsd;es;left;diis")
add_eT_runtest(ccsd_oscillator_strength_eom    "eT;short;ccsd;es;eom")
add_eT_runtest(cc3_gs_energy                   "eT;short;cc3;gs")
add_eT_runtest(cc3_dipole                      "eT;short;cc3;gs;dipole")
add_eT_runtest(cc3_frozen_core                 "eT;short;cc3;gs;es;frozen core")
add_eT_runtest(cc3_right_es_energies           "eT;short;cc3;es;right;diis")
add_eT_runtest(cc3_right_cvs_es_energies       "eT;short;cc3;es;right;cvs;diis")
add_eT_runtest(cc3_left_es_energies            "eT;short;cc3;es;left;diis")
add_eT_runtest(cc3_left_cvs_es_energies        "eT;short;cc3;es;left;cvs;diis")
add_eT_runtest(mlcc2_cnto_full_right_es         "eT;mlcc2;cnto")
add_eT_runtest(mlcc2_pao_full_right_es          "eT;mlcc2;pao")
add_eT_runtest(mlcc2_pao_right_es               "eT;mlcc2;pao")
add_eT_runtest(mlcc2_pao_left_es                "eT;mlcc2;pao")
add_eT_runtest(mlcc2_nto_full_right_es          "eT;mlcc2;nto")
add_eT_runtest(mlcc2_nto_right_es               "eT;mlcc2;nto")
add_eT_runtest(mlcc2_nto_left_es                "eT;mlcc2;nto")
add_eT_runtest(mlcc2_cholesky_full_right_es     "eT;mlcc2;cholesky")
add_eT_runtest(mlcc2_cholesky_right_es          "eT;mlcc2;cholesky")
add_eT_runtest(mlcc2_cholesky_full_left_es      "eT;mlcc2;cholesky")
add_eT_runtest(mlcc2_cholesky_full_right_cvs_es "eT;mlcc2;cholesky;cvs")
add_eT_runtest(mlcc2_cholesky_full_left_cvs_es  "eT;mlcc2;cholesky;cvs")
add_eT_runtest(mlcc2_cholesky_right_cvs_es 		"eT;mlcc2;cholesky;cvs")
add_eT_runtest(mlcc2_cnto_right_es              "eT;mlcc2;cnto")
add_eT_runtest(mlcc2_cnto_left_es               "eT;mlcc2;cnto")
