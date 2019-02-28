
# Inspired by Dalton's test macro written by Radovan Bast 
# This file contains all tests of the eT program 


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
add_eT_runtest(cc2_multipliers                 "eT;short;cc2;gs;multipliers")
add_eT_runtest(cc2_right_es_energies           "eT;short;cc2;es")
add_eT_runtest(cc2_right_es_energies_diis      "eT;short;cc2;es")
add_eT_runtest(cc2_right_cvs_es_energies       "eT;short;cc2;es;cvs")
add_eT_runtest(cc2_left_es_energies            "eT;short;cc2;es")
add_eT_runtest(cc2_lowmem_gs_energy            "eT;short;lowmem-cc2;gs")
add_eT_runtest(cc2_lowmem_right_es_energies    "eT;short;lowmem-cc2;es")
add_eT_runtest(ccsd_gs_energy                  "eT;short;ccsd;gs")
add_eT_runtest(ccsd_right_es_energies          "eT;short;ccsd;es")
add_eT_runtest(ccsd_left_es_energies           "eT;short;ccsd;es")
add_eT_runtest(ccsd_right_cvs_es_energies      "eT;short;ccsd;es;cvs")
add_eT_runtest(ccsd_multipliers                "eT;short;ccsd;gs;multipliers")
add_eT_runtest(ccsd_multipliers_diis           "eT;short;ccsd;gs;multipliers")
add_eT_runtest(cc3_gs_energy                   "eT;short;cc3;gs")
