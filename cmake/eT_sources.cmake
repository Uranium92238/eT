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
#
set(eT_fortran_sources
   src/eT_program.F90
   src/solvers/hf/abstract_hf_solver_class.F90
   src/solvers/hf/scf_diis_hf_class.F90
   src/solvers/hf/scf_hf_class.F90
   src/solvers/cc/diis_cc_gs_class.F90
   src/solvers/cc/diis_cc_es_class.F90
   src/solvers/cc/diis_cc_multipliers_class.F90
   src/solvers/cc/davidson_cc_multipliers_class.F90
   src/solvers/cc/davidson_cc_es_class.F90
   src/solvers/cc/davidson_cc_ip_class.F90
   src/solvers/cc/davidson_cvs_cc_es_class.F90
   src/solvers/cholesky/eri_cd_class.F90
   src/engines/hf_engine_class.F90
   src/engines/abstract_engine_class.F90
   src/engines/gs_engine_class.F90
   src/engines/es_engine_class.F90
   src/engines/multipliers_engine_class.F90
   src/io/disk_manager_class.F90
   src/io/file_class.F90
   src/io/io_utilities.F90
   src/io/io_eT_program.F90
   src/integrals/ao_integral_tool_class.F90
   src/integrals/mo_integral_tool_class.F90
   src/solver_tools/cholesky_array_list_class.F90
   src/memory/batching_index_class.F90
   src/memory/memory_manager_class.F90
   src/various/kinds.F90
   src/various/parameters.F90
   src/tools/interval_class.F90
   src/tools/index_invert.F90
   src/tools/reordering.F90
   src/tools/timings_class.F90
   src/tools/array_utilities.F90
   src/tools/array_analysis.F90
   src/tools/linked_list/array_list_class.F90
   src/tools/linked_list/array_node_class.F90
   src/wavefunctions/wavefunction_class.F90
   src/wavefunctions/hf_class.F90
   src/wavefunctions/uhf_class.F90
   src/wavefunctions/mlhf_class.F90
   src/wavefunctions/ccs_class.F90
   src/wavefunctions/cc2_class.F90
   src/wavefunctions/lowmem_cc2_class.F90
   src/wavefunctions/ccsd_class.F90
   src/wavefunctions/cc3_class.F90
   src/wavefunctions/mp2_class.F90
   src/molecular_system/atomic_class.F90
   src/molecular_system/molecular_system_class.F90
   src/molecular_system/periodic_table.F90
   src/molecular_system/shell_class.F90
   src/molecular_system/basis_set_info.F90
   src/libint/libint_initialization.F90
   src/molecular_system/active_atoms_info_class.F90
   src/solver_tools/diis_tool_class.F90
   src/solver_tools/davidson_tool_class.F90
   src/solver_tools/eigen_davidson_tool_class.F90
   src/solver_tools/linear_davidson_tool_class.F90
   src/submodules/ccsd/omega_ccsd.F90
   src/submodules/ccsd/jacobian_ccsd.F90
   src/submodules/ccsd/jacobian_transpose_ccsd.F90
   src/submodules/cc2/omega_cc2.F90
   src/submodules/cc2/jacobian_cc2.F90
   src/submodules/cc2/jacobian_transpose_cc2.F90
   src/submodules/lowmem_cc2/omega_lowmem_cc2.F90
   src/submodules/lowmem_cc2/jacobian_lowmem_cc2.F90
   src/submodules/cc3/omega_cc3.F90
   src/submodules/cc3/jacobian_cc3.F90
   src/submodules/cc3/jacobian_transpose_cc3.F90
)
