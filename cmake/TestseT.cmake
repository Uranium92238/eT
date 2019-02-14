macro(add_eT_runtest _name _labels)
    add_test(
        ${_name}
        python3 ${PROJECT_BINARY_DIR}/tests/${_name}/test --binary-dir=${PROJECT_BINARY_DIR}  --work-dir=${PROJECT_BINARY_DIR}/tests/${_name} --verbose)
    if(NOT "${_labels}" STREQUAL "")
        set_tests_properties(${_name} PROPERTIES LABELS "${_labels}")
    endif()
endmacro()

#shamelessly stolen from Dalton

# all tests here should contain the label "eT"

# all tests with label short should run in less than 30 seconds sequential on a 2GHz CPU
# short < 30 seconds
# medium > 30 seconds < 120 seconds
# long > 120 seconds < 200 seconds
# verylong > 200 seconds

# NEVER comment out tests, this will bite you in future in a terrible way
# also to minimize conflict release/master use the existing variables to distinguish
# tests which are run only on master

# "long" tests should placed apart on purpose to make it less likely that they
# are run at the same time (then they can be significantly slower)

add_eT_runtest(hf_energy                       "eT;hf;sad")
add_eT_runtest(hf_scf_energy                   "eT;hf;sad")
add_eT_runtest(mp2_energy                      "eT;hf;mp2")
add_eT_runtest(ccsd_gs_energy                  "eT;ccsd;gs")
add_eT_runtest(ccsd_right_es_energies          "eT;ccsd;es")
add_eT_runtest(ccsd_left_es_energies           "eT;ccsd;es")
add_eT_runtest(ccsd_right_cvs_es_energies      "eT;ccsd;es;cvs")
add_eT_runtest(cc3_gs_energy                   "eT;cc3;gs")
