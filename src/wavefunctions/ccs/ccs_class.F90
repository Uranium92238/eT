!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module ccs_class
!
!!
!!    Coupled cluster singles (ccs) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
!
   use mo_integral_tool_class, only : mo_integral_tool
!
   use reordering
!
   use sequential_file_class, only : sequential_file
   use string_utilities, only : convert_to_uppercase
   use array_utilities, only : zero_array_complex
   use array_utilities, only : get_l2_norm, copy_and_scale, copy_and_scale_complex, copy_, our_zdotu
   use array_utilities, only : get_abs_max_w_index, get_n_lowest, get_n_highest
   use array_utilities, only: quicksort_with_index_descending, are_vectors_parallel
   use index_invert, only : invert_compound_index, invert_packed_index
   use batching_index_class, only : batching_index
   use timings_class, only : timings
   use file_storer_class, only: file_storer
!
   implicit none
!
   type, extends(wavefunction) :: ccs
!
      integer :: n_gs_amplitudes
      integer :: n_es_amplitudes
      integer :: n_t1
      integer :: n_singlet_states
      integer :: n_bath_orbitals
!
      integer :: n_core_MOs
      logical :: cvs 
!
      real(dp), dimension(:), allocatable :: left_excitation_energies
      real(dp), dimension(:), allocatable :: right_excitation_energies
!
      logical :: bath_orbital
!
      logical :: need_g_abcd
!
      type(sequential_file) :: t_file, tbar_file
      type(sequential_file) :: excitation_energies_file
      type(sequential_file) :: restart_file
!
      type(file_storer), allocatable :: r_files
      type(file_storer), allocatable :: l_files
!
      type(mo_integral_tool) :: integrals
!
      real(dp)    :: hf_energy
      complex(dp) :: hf_energy_complex
!
      real(dp),    dimension(:,:), allocatable :: t1
      complex(dp), dimension(:,:), allocatable :: t1_complex
!
      real(dp),    dimension(:,:), allocatable :: t1bar
      complex(dp), dimension(:,:), allocatable :: t1bar_complex
!
      real(dp),    dimension(:,:), allocatable :: fock_ij
      complex(dp), dimension(:,:), allocatable :: fock_ij_complex
!
      real(dp),    dimension(:,:), allocatable :: fock_ia
      complex(dp), dimension(:,:), allocatable :: fock_ia_complex
!
      real(dp),    dimension(:,:), allocatable :: fock_ai
      complex(dp), dimension(:,:), allocatable :: fock_ai_complex
!
      real(dp),    dimension(:,:), allocatable :: fock_ab
      complex(dp), dimension(:,:), allocatable :: fock_ab_complex
!
      real(dp),    dimension(:,:), allocatable :: density
      complex(dp), dimension(:,:), allocatable :: density_complex

      real(dp),    dimension(:,:), allocatable :: left_transition_density
      real(dp),    dimension(:,:), allocatable :: right_transition_density
!
      integer, dimension(:), allocatable :: core_MOs
!
   contains
!
!     Initialization/destruction procedures
!
      procedure :: initialize_amplitudes                         => initialize_amplitudes_ccs
      procedure :: initialize_amplitudes_complex                 => initialize_amplitudes_ccs_complex
!
      procedure :: destruct_amplitudes                           => destruct_amplitudes_ccs
      procedure :: destruct_amplitudes_complex                   => destruct_amplitudes_ccs_complex
!
      procedure :: initialize_multipliers                        => initialize_multipliers_ccs
      procedure :: initialize_multipliers_complex                => initialize_multipliers_ccs_complex
!
      procedure :: destruct_multipliers                          => destruct_multipliers_ccs
      procedure :: destruct_multipliers_complex                  => destruct_multipliers_ccs_complex
!
      procedure :: initialize_fock                               => initialize_fock_ccs
      procedure :: initialize_fock_complex                       => initialize_fock_ccs_complex
!
      procedure :: initialize_fock_ij                            => initialize_fock_ij_ccs
      procedure :: initialize_fock_ij_complex                    => initialize_fock_ij_ccs_complex
!
      procedure :: destruct_fock_ij                              => destruct_fock_ij_ccs
      procedure :: destruct_fock_ij_complex                      => destruct_fock_ij_ccs_complex
!
      procedure :: initialize_fock_ia                            => initialize_fock_ia_ccs
      procedure :: initialize_fock_ia_complex                    => initialize_fock_ia_ccs_complex
!
      procedure :: destruct_fock_ia                              => destruct_fock_ia_ccs
      procedure :: destruct_fock_ia_complex                      => destruct_fock_ia_ccs_complex
!
      procedure :: initialize_fock_ai                            => initialize_fock_ai_ccs
      procedure :: initialize_fock_ai_complex                    => initialize_fock_ai_ccs_complex
!
      procedure :: destruct_fock_ai                              => destruct_fock_ai_ccs
      procedure :: destruct_fock_ai_complex                      => destruct_fock_ai_ccs_complex
!
      procedure :: initialize_fock_ab                            => initialize_fock_ab_ccs
      procedure :: initialize_fock_ab_complex                    => initialize_fock_ab_ccs_complex
!
      procedure :: destruct_fock_ab                              => destruct_fock_ab_ccs
      procedure :: destruct_fock_ab_complex                      => destruct_fock_ab_ccs_complex
!
      procedure :: initialize_t1                                 => initialize_t1_ccs
      procedure :: initialize_t1_complex                         => initialize_t1_ccs_complex
!
      procedure :: destruct_t1                                   => destruct_t1_ccs
      procedure :: destruct_t1_complex                           => destruct_t1_ccs_complex
!
      procedure :: initialize_t1bar                              => initialize_t1bar_ccs
      procedure :: initialize_t1bar_complex                      => initialize_t1bar_ccs_complex
!
      procedure :: destruct_t1bar                                => destruct_t1bar_ccs
      procedure :: destruct_t1bar_complex                        => destruct_t1bar_ccs_complex
!
      procedure :: initialize_gs_density                         => initialize_gs_density_ccs
      procedure :: initialize_gs_density_complex                 => initialize_gs_density_ccs_complex
!
      procedure :: destruct_gs_density                           => destruct_gs_density_ccs
      procedure :: destruct_gs_density_complex                   => destruct_gs_density_ccs_complex
!
      procedure :: initialize_transition_densities               => initialize_transition_densities_ccs
!
      procedure :: destruct_transition_densities                 => destruct_transition_densities_ccs
!
      procedure :: initialize_right_excitation_energies          => initialize_right_excitation_energies_ccs
!
      procedure :: destruct_right_excitation_energies            => destruct_right_excitation_energies_ccs
!
      procedure :: initialize_left_excitation_energies           => initialize_left_excitation_energies_ccs
!
      procedure :: destruct_left_excitation_energies             => destruct_left_excitation_energies_ccs
!
      procedure :: initialize_core_MOs                           => initialize_core_MOs_ccs 
      procedure :: destruct_core_MOs                             => destruct_core_MOs_ccs
!
!     File handling procedures
!
      procedure :: read_hf                                       => read_hf_ccs
      procedure :: initialize_files                              => initialize_files_ccs
      procedure :: initialize_cc_files                           => initialize_cc_files_ccs
      procedure :: initialize_ground_state_files                 => initialize_ground_state_files_ccs
      procedure :: initialize_excited_state_files                => initialize_excited_state_files_ccs
!
      procedure :: read_settings                                 => read_settings_ccs
      procedure :: read_cvs_settings                             => read_cvs_settings_ccs
!
      procedure :: read_singles_vector                           => read_singles_vector_ccs
      procedure :: save_amplitudes                               => save_amplitudes_ccs
      procedure :: read_amplitudes                               => read_amplitudes_ccs
      procedure :: save_multipliers                              => save_multipliers_ccs
      procedure :: read_multipliers                              => read_multipliers_ccs
      procedure :: save_excited_state                            => save_excited_state_ccs
      procedure :: read_excited_state                            => read_excited_state_ccs
      procedure :: save_excitation_energies                      => save_excitation_energies_ccs
      procedure :: read_excitation_energies                      => read_excitation_energies_ccs
      procedure :: read_frozen_orbital_terms                     => read_frozen_orbital_terms_ccs
!
      procedure :: write_cc_restart                              => write_cc_restart_ccs
!
      procedure :: save_tbar_intermediates                       => save_tbar_intermediates_ccs
!
      procedure :: get_n_excited_states_on_file                  => get_n_excited_states_on_file_ccs
      procedure :: get_n_excitation_energies_on_file             => get_n_excitation_energies_on_file_ccs
!
!     Set/get procedures
!
      procedure :: set_amplitudes                                => set_amplitudes_ccs
      procedure :: set_amplitudes_complex                        => set_amplitudes_ccs_complex
!
      procedure :: get_amplitudes                                => get_amplitudes_ccs
      procedure :: get_amplitudes_complex                        => get_amplitudes_ccs_complex
!
      procedure :: set_multipliers                               => set_multipliers_ccs
      procedure :: set_multipliers_complex                       => set_multipliers_ccs_complex
!
      procedure :: get_multipliers                               => get_multipliers_ccs
      procedure :: get_multipliers_complex                       => get_multipliers_ccs_complex
!
      procedure :: set_fock                                      => set_fock_ccs
      procedure :: set_fock_complex                              => set_fock_ccs_complex
!
      procedure :: get_gs_orbital_differences                    => get_gs_orbital_differences_ccs
      procedure :: get_es_orbital_differences                    => get_gs_orbital_differences_ccs
!
!     Procedures related to the Fock matrix
!
      procedure :: construct_fock                                => construct_fock_ccs
      procedure :: construct_fock_complex                        => construct_fock_ccs_complex
!
      procedure :: add_frozen_core_fock_term                     => add_frozen_core_fock_term_ccs
      procedure :: add_frozen_core_fock_term_complex             => add_frozen_core_fock_term_ccs_complex

      procedure :: add_frozen_hf_fock_term                       => add_frozen_hf_fock_term_ccs
      procedure :: add_frozen_hf_fock_term_complex               => add_frozen_hf_fock_term_ccs_complex

      procedure :: add_molecular_mechanics_fock_term             => add_molecular_mechanics_fock_term_ccs
      procedure :: add_molecular_mechanics_fock_term_complex     => add_molecular_mechanics_fock_term_ccs_complex
!
      procedure :: add_t1_fock_length_dipole_term                => add_t1_fock_length_dipole_term_ccs
      procedure :: add_t1_fock_length_dipole_term_complex        => add_t1_fock_length_dipole_term_ccs_complex
!
      procedure :: add_pcm_fock_contribution                     => add_pcm_fock_contribution_ccs
      procedure :: add_pcm_fock_contribution_complex             => add_pcm_fock_contribution_ccs_complex
!
!     Procedures related to the omega vector
!
      procedure :: construct_omega                               => construct_omega_ccs
      procedure :: construct_omega_complex                       => construct_omega_ccs_complex
!
      procedure :: omega_ccs_a1                                  => omega_ccs_a1_ccs
      procedure :: omega_ccs_a1_complex                          => omega_ccs_a1_ccs_complex
!
      procedure :: construct_Jacobian_transform                  => construct_Jacobian_transform_ccs
!
!     Procedures related to the multiplier equation vector
!
      procedure :: construct_multiplier_equation                 => construct_multiplier_equation_ccs
      procedure :: construct_multiplier_equation_complex         => construct_multiplier_equation_ccs_complex
!
      procedure :: construct_eta                                 => construct_eta_ccs
      procedure :: construct_eta_complex                         => construct_eta_ccs_complex
!
!     Procedures related to the Jacobian transformation
!
      procedure :: jacobian_transformation                       => jacobian_transformation_ccs
      procedure :: jacobian_ccs_a1                               => jacobian_ccs_a1_ccs
      procedure :: jacobian_ccs_b1                               => jacobian_ccs_b1_ccs
!
!     Procedures related to the Jacobian transpose transformation
!
      procedure :: jacobian_transpose_transformation             => jacobian_transpose_transformation_ccs
      procedure :: jacobian_transpose_transformation_complex     => jacobian_transpose_transformation_ccs_complex
!
      procedure :: jacobian_transpose_ccs_a1                     => jacobian_transpose_ccs_a1_ccs
      procedure :: jacobian_transpose_ccs_a1_complex             => jacobian_transpose_ccs_a1_ccs_complex
!
      procedure :: jacobian_transpose_ccs_b1                     => jacobian_transpose_ccs_b1_ccs
      procedure :: jacobian_transpose_ccs_b1_complex             => jacobian_transpose_ccs_b1_ccs_complex
!
!     Procedures related to properties
!
      procedure :: calculate_expectation_value                   => calculate_expectation_value_ccs
      procedure :: calculate_expectation_value_complex           => calculate_expectation_value_ccs_complex
!
      procedure :: calculate_energy                              => calculate_energy_ccs
      procedure :: calculate_energy_complex                      => calculate_energy_ccs_complex
!
      procedure :: calculate_energy_omega_term                   => calculate_energy_omega_term_ccs
      procedure :: calculate_energy_omega_term_complex           => calculate_energy_omega_term_ccs_complex
!
      procedure :: calculate_energy_length_dipole_term           => calculate_energy_length_dipole_term_ccs
      procedure :: calculate_energy_length_dipole_term_complex   => calculate_energy_length_dipole_term_ccs_complex
!
      procedure :: construct_gs_density                          => construct_gs_density_ccs
      procedure :: construct_gs_density_complex                  => construct_gs_density_ccs_complex
!
      procedure :: construct_right_transition_density            => construct_right_transition_density_ccs
      procedure :: construct_left_transition_density             => construct_left_transition_density_ccs
!
      procedure :: density_ccs_ref_ref_oo                        => density_ccs_ref_ref_oo_ccs
      procedure :: density_ccs_ref_ref_oo_complex                => density_ccs_ref_ref_oo_ccs_complex
!
      procedure :: density_ccs_mu_ref_vo                         => density_ccs_mu_ref_vo_ccs
      procedure :: density_ccs_mu_ref_vo_complex                 => density_ccs_mu_ref_vo_ccs_complex
!
      procedure :: density_ccs_ref_mu_ov                         => density_ccs_ref_mu_ov_ccs
      procedure :: density_ccs_mu_nu_oo                          => density_ccs_mu_nu_oo_ccs
      procedure :: density_ccs_mu_nu_vv                          => density_ccs_mu_nu_vv_ccs
!
      procedure :: density_mu_mu_oo                              => density_mu_mu_oo_ccs
      procedure :: density_mu_ref                                => density_mu_ref_ccs
!
      procedure :: construct_etaX                                => construct_etaX_ccs
      procedure :: construct_eom_etaX                            => construct_eom_etaX_ccs
      procedure :: etaX_ccs_a1                                   => etaX_ccs_a1_ccs
      procedure :: etaX_ccs_b1                                   => etaX_ccs_b1_ccs
      procedure :: construct_csiX                                => construct_csiX_ccs
      procedure :: csiX_ccs_a1                                   => csiX_ccs_a1_ccs
      procedure :: etaX_eom_a                                    => etaX_eom_a_ccs
      procedure :: calculate_transition_strength                 => calculate_transition_strength_ccs
!
!     Routines related to post-processing excited states
!
      procedure :: biorthonormalize_L_and_R                      => biorthonormalize_L_and_R_ccs
      procedure :: L_R_overlap                                   => L_R_overlap_ccs
      procedure :: check_for_degeneracies                        => check_for_degeneracies_ccs
!
!     One-electron interals
!
      procedure :: t1_transform                                  => t1_transform_ccs
      procedure :: t1_transform_complex                          => t1_transform_ccs_complex
!
      procedure :: ao_to_t1_transformation                       => ao_to_t1_transformation_ccs
      procedure :: ao_to_t1_transformation_complex               => ao_to_t1_transformation_ccs_complex
!
      procedure :: construct_h                                   => construct_h_ccs
      procedure :: construct_h_complex                           => construct_h_ccs_complex
!
      procedure :: construct_mu                                  => construct_mu_ccs
      procedure :: construct_mu_complex                          => construct_mu_ccs_complex
!
      procedure :: construct_q                                   => construct_q_ccs
      procedure :: construct_q_complex                           => construct_q_ccs_complex
!
!     Two-electron integrals
!
      procedure :: t1_transform_4                                => t1_transform_4_ccs
      procedure :: t1_transform_4_complex                        => t1_transform_4_ccs_complex
!
      procedure :: get_ovov                                      => get_ovov_ccs
      procedure :: get_ovov_complex                              => get_ovov_ccs_complex
!
      procedure :: get_vovo                                      => get_vovo_ccs
      procedure :: get_vovo_complex                              => get_vovo_ccs_complex
!
      procedure :: get_vvoo                                      => get_vvoo_ccs
      procedure :: get_vvoo_complex                              => get_vvoo_ccs_complex
!
      procedure :: get_voov                                      => get_voov_ccs
      procedure :: get_voov_complex                              => get_voov_ccs_complex
!
      procedure :: get_ovvo                                      => get_ovvo_ccs
      procedure :: get_ovvo_complex                              => get_ovvo_ccs_complex
!
      procedure :: get_oovv                                      => get_oovv_ccs
      procedure :: get_oovv_complex                              => get_oovv_ccs_complex
!
      procedure :: get_oooo                                      => get_oooo_ccs
      procedure :: get_oooo_complex                              => get_oooo_ccs_complex
!
      procedure :: get_vvvv                                      => get_vvvv_ccs
      procedure :: get_vvvv_complex                              => get_vvvv_ccs_complex
!
      procedure :: get_ooov                                      => get_ooov_ccs
      procedure :: get_ooov_complex                              => get_ooov_ccs_complex
!
      procedure :: get_oovo                                      => get_oovo_ccs
      procedure :: get_oovo_complex                              => get_oovo_ccs_complex
!
      procedure :: get_ovoo                                      => get_ovoo_ccs
      procedure :: get_ovoo_complex                              => get_ovoo_ccs_complex
!
      procedure :: get_vooo                                      => get_vooo_ccs
      procedure :: get_vooo_complex                              => get_vooo_ccs_complex
!
      procedure :: get_vvvo                                      => get_vvvo_ccs
      procedure :: get_vvvo_complex                              => get_vvvo_ccs_complex
!
      procedure :: get_vvov                                      => get_vvov_ccs
      procedure :: get_vvov_complex                              => get_vvov_ccs_complex
!
      procedure :: get_vovv                                      => get_vovv_ccs
      procedure :: get_vovv_complex                              => get_vovv_ccs_complex
!
      procedure :: get_ovvv                                      => get_ovvv_ccs
      procedure :: get_ovvv_complex                              => get_ovvv_ccs_complex

!
!     Preparation procedures
!
      procedure :: prepare_for_jacobian                          => prepare_for_jacobian_ccs
!
      procedure :: prepare_for_jacobian_transpose                => prepare_for_jacobian_transpose_ccs
      procedure :: prepare_for_jacobian_transpose_complex        => prepare_for_jacobian_transpose_ccs_complex
!
      procedure :: prepare_for_multiplier_equation               => prepare_for_multiplier_equation_ccs
      procedure :: prepare_for_multiplier_equation_complex       => prepare_for_multiplier_equation_ccs_complex
!
      procedure :: prepare_for_density                           => prepare_for_density_ccs
      procedure :: prepare_for_density_complex                   => prepare_for_density_ccs_complex
!
      procedure :: approximate_double_excitation_vectors         => approximate_double_excitation_vectors_ccs
!
!     Frozen core
!
      procedure :: construct_t1_fock_fc_term                     => construct_t1_fock_fc_term_ccs
      procedure :: construct_t1_fock_fc_term_complex             => construct_t1_fock_fc_term_ccs_complex
!
      procedure :: construct_t1_fock_frozen_hf_term              => construct_t1_fock_frozen_hf_term_ccs
      procedure :: construct_t1_fock_frozen_hf_term_complex      => construct_t1_fock_frozen_hf_term_ccs_complex
!
!     MO preparations
!
      procedure :: mo_preparations                               => mo_preparations_ccs
!   
!     Debug 
!
      procedure :: omega_for_jacobian_debug                      => omega_for_jacobian_debug_ccs
      procedure :: amplitudes_for_jacobian_debug                 => amplitudes_for_jacobian_debug_ccs
      procedure :: normalization_for_jacobian_debug              => normalization_for_jacobian_debug_ccs
      procedure :: numerical_test_jacobian                       => numerical_test_jacobian_ccs
!
!     MM and PCM reading fock matrices
!
      procedure :: read_mm_fock_contributions                  => read_mm_fock_contributions_ccs
      procedure :: read_pcm_fock_contributions                 => read_pcm_fock_contributions_ccs
!
!     Core-valence separation procedures
!
      procedure :: get_cvs_projector                             => get_cvs_projector_ccs
      procedure :: set_cvs_start_indices                         => set_cvs_start_indices_ccs
!
!     Other procedures
!
      procedure :: cleanup                                       => cleanup_ccs
      procedure :: general_cc_preparations                       => general_cc_preparations_ccs
!
      procedure :: is_restart_safe                               => is_restart_safe_ccs
!
      procedure :: set_initial_amplitudes_guess                  => set_initial_amplitudes_guess_ccs
      procedure :: set_ip_start_indices                          => set_ip_start_indices_ccs
      procedure :: get_ip_projector                              => get_ip_projector_ccs
      procedure :: get_t1_diagnostic                             => get_t1_diagnostic_ccs
!
      procedure :: form_newton_raphson_t_estimate                => form_newton_raphson_t_estimate_ccs
!
      procedure :: print_dominant_x_amplitudes                   => print_dominant_x_amplitudes_ccs
      procedure :: print_dominant_amplitudes                     => print_dominant_amplitudes_ccs
      procedure :: print_dominant_x1                             => print_dominant_x1_ccs
!
      procedure :: make_bath_orbital                             => make_bath_orbital_ccs
!
      procedure :: construct_molecular_gradient                  => construct_molecular_gradient_ccs
!
      procedure :: get_g_pqrs_required                           => get_g_pqrs_required_ccs
!
!     Procedures related to time dependency
!
      procedure :: make_complex                                  => make_complex_ccs
      procedure :: make_ccs_complex                              => make_ccs_complex_ccs
!
      procedure :: construct_complex_time_derivative             => construct_complex_time_derivative_ccs
!
      procedure :: construct_complex_time_derivative_amplitudes  => construct_complex_time_derivative_amplitudes_ccs
      procedure :: construct_complex_time_derivative_multipliers => construct_complex_time_derivative_multipliers_ccs
!

   end type ccs
!
!
   interface 
!
      include "file_handling_ccs_interface.F90"
      include "initialize_destruct_ccs_interface.F90"
      include "set_get_ccs_interface.F90"
      include "omega_ccs_interface.F90"
      include "multiplier_equation_ccs_interface.F90"
      include "jacobian_ccs_interface.F90"
      include "jacobian_transpose_ccs_interface.F90"
      include "zop_ccs_interface.F90"
      include "fop_ccs_interface.F90"
      include "oei_ccs_interface.F90"
      include "tei_ccs_interface.F90"
      include "t1_ccs_interface.F90"
      include "fock_ccs_interface.F90"
      include "debug_jacobian_ccs_interface.F90"
!
      include "complex_ccs_interface.F90"
!
      include "autogenerated_complex_files/fock_ccs_interface_complex.F90"
      include "autogenerated_complex_files/initialize_destruct_ccs_interface_complex.F90"
      include "autogenerated_complex_files/jacobian_transpose_ccs_interface_complex.F90"
      include "autogenerated_complex_files/multiplier_equation_ccs_interface_complex.F90"
      include "autogenerated_complex_files/omega_ccs_interface_complex.F90"
      include "autogenerated_complex_files/set_get_ccs_interface_complex.F90"
      include "autogenerated_complex_files/t1_ccs_interface_complex.F90"
      include "autogenerated_complex_files/tei_ccs_interface_complex.F90"
      include "autogenerated_complex_files/zop_ccs_interface_complex.F90"
!
   end interface
!
!
   interface ccs 
!
      procedure :: new_ccs 
!
   end interface ccs
!
!
contains
!
!
   function new_ccs(system) result(wf)
!!
!!    New CCS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(ccs) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%name_ = 'ccs'
!
      call wf%general_cc_preparations(system)
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
      wf%need_g_abcd     = .false.
!
      call wf%initialize_fock()
!
   end function new_ccs
!
!
   subroutine general_cc_preparations_ccs(wf, system)
!!
!!    General CC preparations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%system   => system
!
!     Initialize CC files
!
      call wf%initialize_files()
!
!     Logicals for special methods
!
      wf%bath_orbital = .false.
      wf%frozen_core = .false.
      wf%frozen_hf_mos = .false.
      wf%cvs = .false.
      wf%need_g_abcd = .false.
!
!     Read CC settings from eT.inp (from cc section)
!
      call wf%read_settings()
!
!     Read necessary information from HF
!
      call wf%read_hf()
!
!     Handle changes in the number of MOs as a result of 
!     special methods
!
      if (wf%bath_orbital) call wf%make_bath_orbital()
!
      if (wf%frozen_core .or. wf%frozen_hf_mos) call wf%read_frozen_orbital_terms()
!
!     Read MM or PCM file if present
!
      if (wf%system%mm_calculation)  call wf%read_mm_fock_contributions()
      if (wf%system%pcm_calculation) call wf%read_pcm_fock_contributions()
!
!     print orbital space info for cc
      call output%printf(' - Number of orbitals for coupled cluster calculation', fs='(/t3,a)', pl='minimal')
      call output%printf('Number of occupied orbitals:    (i12)',ints=[wf%n_o], fs='(/t6,a)', pl='minimal')
      call output%printf('Number of virtual orbitals:     (i12)',ints=[wf%n_v], fs='(t6,a)', pl='minimal')
      call output%printf('Number of molecular orbitals:   (i12)',ints=[wf%n_mo], fs='(t6,a)', pl='minimal')
      call output%printf('Number of atomic orbitals:      (i12)',ints=[wf%n_ao], fs='(t6,a)', pl='minimal')
!
   end subroutine general_cc_preparations_ccs
!
!
   subroutine cleanup_ccs(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad and
!!    Alexander C. Paul , 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%destruct_amplitudes()
      call wf%destruct_multipliers()
      call wf%destruct_right_excitation_energies()
      call wf%destruct_left_excitation_energies()
      call wf%destruct_mo_fock_fc_term()
      call wf%destruct_mo_fock_frozen_hf_term()
!
      call wf%destruct_mm_matrices()
!
      call wf%destruct_pcm_matrices()
!
      call output%printf('- Cleaning up ' // trim(wf%name_) // ' wavefunction', pl='v', fs='(/t3,a)')
!
   end subroutine cleanup_ccs
!
!
   subroutine read_settings_ccs(wf)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad, Aug 2019
!!
!!    Reads the cc-section of the input file 
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call wf%read_frozen_orbitals_settings()
!
      if (input%requested_section('cc')) then
!
         if (input%requested_keyword_in_section('bath orbital','cc')) then 
!
            wf%bath_orbital = .true.
            wf%n_bath_orbitals = 1 ! May want to expand the number of bath orbitals later on
!
         endif
!
      endif
!
   end subroutine read_settings_ccs
!
!
   subroutine read_hf_ccs(wf)
!!
!!    Read HF file
!!    Written by Rolf H. Myhre, May 2019
!!    Short routine to read the HF information from disk
!!
      implicit none
!
      class(ccs) :: wf
!
      type(sequential_file) :: orbital_information_file
      type(sequential_file) :: CC_orbitals_file, CC_orbital_energies_file
!
      orbital_information_file = sequential_file('orbital_information')
      call orbital_information_file%open_('read', 'rewind')
!
      call orbital_information_file%read_(wf%n_o)     
      call orbital_information_file%read_(wf%n_v)         
      call orbital_information_file%read_(wf%n_ao)     
      call orbital_information_file%read_(wf%n_mo)     
      call orbital_information_file%read_(wf%hf_energy)    
!
      call orbital_information_file%close_()
!
!     Set orbital coefficients and energies
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      CC_orbitals_file = sequential_file('cc_orbital_coefficients')
      call CC_orbitals_file%open_('read', 'rewind')
!
      call CC_orbitals_file%read_(wf%orbital_coefficients, wf%n_ao*wf%n_mo)
!
      call CC_orbitals_file%close_('keep')
!
      CC_orbital_energies_file = sequential_file('cc_orbital_energies')
      call CC_orbital_energies_file%open_('read', 'rewind')
!
      call CC_orbital_energies_file%read_(wf%orbital_energies, wf%n_mo)
!
      call CC_orbital_energies_file%close_('keep')
!
   end subroutine read_hf_ccs
!
!
   subroutine write_cc_restart_ccs(wf)
!!
!!    Write CC restart file
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
!     Write information to restart file 
!
      call wf%restart_file%open_('write', 'rewind')
!
      call wf%restart_file%write_(wf%n_o)
      call wf%restart_file%write_(wf%n_v)
      call wf%restart_file%write_(wf%n_gs_amplitudes)
      call wf%restart_file%write_(wf%n_es_amplitudes)
!
      call wf%restart_file%close_()
!
   end subroutine write_cc_restart_ccs
!
!
   subroutine construct_molecular_gradient_ccs(wf, E_qk)
!!
!!    Construct molecular gradient 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(3, wf%system%n_atoms), intent(inout) :: E_qk 
!
      E_qk = zero 
!
      call output%error_msg('Molecular gradient not implemented for ' // trim(wf%name_))
!
   end subroutine construct_molecular_gradient_ccs
!
!
   subroutine is_restart_safe_ccs(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      character(len=*), intent(in) :: task 
!
      integer :: n_o, n_v, n_gs_amplitudes, n_es_amplitudes
!
      call wf%restart_file%open_('read', 'rewind')
!
      call wf%restart_file%read_(n_o)
      call wf%restart_file%read_(n_v)
      call wf%restart_file%read_(n_gs_amplitudes)
      call wf%restart_file%read_(n_es_amplitudes)
!
      call wf%restart_file%close_()
!
      if (n_o .ne. wf%n_o) call output%error_msg('attempted to restart from inconsistent number ' // &
                                                   'of occupied orbitals.')
!
      if (n_v .ne. wf%n_v) call output%error_msg('attempted to restart from inconsistent number ' // &
                                                   'of virtual orbitals.')
!
      if (trim(task) == 'ground state') then 
!
         if (n_gs_amplitudes .ne. wf%n_gs_amplitudes) &
            call output%error_msg('attempted to restart from inconsistent number ' // &
                                    'of ground state amplitudes.')    
!
      elseif (trim(task) == 'excited state') then    
!
         if (n_es_amplitudes .ne. wf%n_es_amplitudes) &
            call output%error_msg('attempted to restart from inconsistent number ' // &
                                    'of excited state amplitudes.')     
!
      else
!
         call output%error_msg('attempted to restart, but the task was not recognized: ' // task)
!
      endif   
!
   end subroutine is_restart_safe_ccs
!
!
   subroutine set_initial_amplitudes_guess_ccs(wf)
!!
!!    Set initial amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call zero_array(wf%t1, wf%n_t1)
!
   end subroutine set_initial_amplitudes_guess_ccs
!
!
   subroutine get_gs_orbital_differences_ccs(wf, orbital_differences, N)
!!
!!    Get orbital differences
!!    Written by Sarai D. Folkestad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in) :: N
      real(dp), dimension(N), intent(inout) :: orbital_differences
!
      integer :: a, i, ai
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai) = wf%orbital_energies(a + wf%n_o) - wf%orbital_energies(i)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_gs_orbital_differences_ccs
!
!
   subroutine construct_Jacobian_transform_ccs(wf, r_or_l, X, w)
!!
!!    Construct Jacobian transform
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Modified by Rolf H. Myhre, Oct 2019
!!
!!    Constructs R = AX or R = A^T X
!!
!!    Removed calculation of residual, this is now done in the solver
!!
!!    Wrapper for Jacobian transformations
!!
!!    r_or_l: string that should be 'left' or 'right', 
!!            determines if Jacobian or Jacobian transpose is called
!!
!!    X: On input contains the vector to transform, 
!!       on output contains the transformed vector
!!
!!    w: Excitation energy. Only used for debug prints for CCS, CCSD etc.
!!       but is passed to the effective_jacobian_transform for lowmem_CC2 and CC3
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout)   :: X
!
      real(dp), intent(in), optional :: w
!
      if (present(w)) then
         call output%printf('Calling Jacobian (a0) transform with energy: (f19.12)', &
                            pl='debug', chars=[r_or_l], reals=[w])
      endif
!
!     Compute the transformed matrix
      if (r_or_l .eq. "right") then
!
         call wf%jacobian_transformation(X) ! X <- AX
!
      elseif (r_or_l .eq. 'left') then
!
         call wf%jacobian_transpose_transformation(X) ! X <- XA
!
      else
!
         call output%error_msg('Neither left nor right in construct_Jacobian_transform')
!
      endif
!
   end subroutine construct_Jacobian_transform_ccs
!
!
   subroutine get_cvs_projector_ccs(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai
!
      call zero_array(projector, wf%n_es_amplitudes)
!
      do core = 1, n_cores
!
        i = core_MOs(core)
!
        do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai) = one
!
        enddo
     enddo
!
   end subroutine get_cvs_projector_ccs
!
!
   subroutine print_dominant_amplitudes_ccs(wf)
!!
!!    Print dominant amplitudes
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      call wf%print_dominant_x1(wf%t1,'t')
!
   end subroutine print_dominant_amplitudes_ccs
!
!
   subroutine print_dominant_x_amplitudes_ccs(wf, x, tag)
!!
!!    Print dominant amplitudes
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: x
!
      character(len=1) :: tag
!
      call wf%print_dominant_x1(x(1:wf%n_t1),tag)
!
   end subroutine print_dominant_x_amplitudes_ccs
!
!
   subroutine print_dominant_x1_ccs(wf, x1, tag)
!!
!!    Print dominant x1
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Prints the 10 most dominant single amplitudes,
!!    or sorts them if there are fewer than twenty of them.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: x1
      character(len=1), intent(in)                :: tag
!
      real(dp), dimension(:), allocatable :: abs_x1
!
      integer, dimension(:), allocatable :: dominant_indices
      real(dp), dimension(:), allocatable     :: dominant_values
!
      integer :: n_elements, elm, i, a
!
!     Sort according to largest contributions
!
      call mem%alloc(abs_x1, wf%n_t1)
      abs_x1 = abs(x1)
!
      n_elements = 10
      if (n_elements .gt. wf%n_t1) n_elements = wf%n_t1
!
      call mem%alloc(dominant_indices, n_elements)
      call mem%alloc(dominant_values, n_elements)
!
      dominant_indices = 0
      dominant_values  = zero
      call get_n_highest(n_elements, wf%n_t1, abs_x1, dominant_values, dominant_indices)
!
!     Print largest contributions
!
      call output%printf('Largest single amplitudes:', pl='m', fs='(/t6,a)')
      call output%print_separator('m', 35, '-', fs='(t6,a)')
      call output%printf('a       i         ' // tag // '(a,i)', pl='m', fs='(t9,a)')
      call output%print_separator('m', 35, '-', fs='(t6,a)')
!
      do elm = 1, n_elements
!
         call invert_compound_index(dominant_indices(elm), a, i, wf%n_v, wf%n_o)
!
         call output%printf('(i4)    (i3)   (f19.12)', fs='(t6,a)', pl='m', &
                             ints=[a, i], reals=[x1(dominant_indices(elm))])
!
      enddo
!
      call output%print_separator('m', 36, '-', fs='(t6,a)')
!
      call mem%dealloc(dominant_indices, n_elements)
      call mem%dealloc(dominant_values, n_elements)
      call mem%dealloc(abs_x1, wf%n_t1)
!
   end subroutine print_dominant_x1_ccs
!
!
   real(dp) function get_t1_diagnostic_ccs(wf)
!!
!!    Get t1 diagnostic
!!    Written by Eirik F. Kjønstad
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      get_t1_diagnostic_ccs = get_l2_norm(wf%t1, wf%n_t1)
      get_t1_diagnostic_ccs = get_t1_diagnostic_ccs/sqrt(real(wf%system%n_electrons,kind=dp))
!
   end function get_t1_diagnostic_ccs
!
!
   subroutine set_cvs_start_indices_ccs(wf, start_indices)
!!
!!    Set CVS start indices
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      integer, dimension(wf%n_singlet_states) :: start_indices
!
!     Local variables
!
      integer :: a, i, current_root, core
!
      logical :: all_selected
!
!     Calculate start indices using Koopman's
!
      all_selected = .false.
      a =  0
      current_root = 0
!
      do while (.not. all_selected)
!
         a = a + 1
!
         do core = 1, wf%n_core_MOs
!
            i = wf%core_MOs(core)
!
            current_root = current_root + 1
            start_indices(current_root) = wf%n_v*( i - 1) + a
!
            if (current_root .eq. wf%n_singlet_states) then
!
               all_selected = .true.
               exit
!
            endif
!
         enddo
!
      enddo
!
   end subroutine set_cvs_start_indices_ccs
!
!
   real(dp) function L_R_overlap_ccs(wf, L, left_state, R, right_state)
!!
!!    Calculates the overlap between a given left and right state
!!    Written by Josefine Andersen, Apr 2019
!!    Modified by Alexander C. Paul to be overwritable for CC3
!!
      class(ccs), intent(in) :: wf

      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L
      integer, intent(in) :: left_state
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R
      integer, intent(in) :: right_state
!
      real(dp) :: ddot
!
      L_R_overlap_ccs = ddot(wf%n_es_amplitudes, L, 1, R, 1)
!
      call output%printf('Overlap of (i0). left and (i0). right state: (f15.10)',   &
                          fs='(/t6,a)', ints=[left_state, right_state], reals=[L_R_overlap_ccs], pl='debug')
!
   end function L_R_overlap_ccs
!
!
   subroutine form_newton_raphson_t_estimate_ccs(wf, t, dt)
!!
!!    Form Newton-Raphson t estimate 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
!!    Here, t is the full amplitude vector and dt is the correction to the amplitude vector.
!!
!!    The correction is assumed to be obtained from either 
!!    solving the Newton-Raphson equation
!!
!!       A dt = -omega, 
!!
!!    where A and omega are given in the biorthonormal basis,
!!    or from the quasi-Newton equation (A ~ diagonal with diagonal = epsilon) 
!!
!!        dt = -omega/epsilon
!!
!!    Epsilon is the vector of orbital differences. 
!!
!!    On exit, t = t + dt, where the appropriate basis change has been accounted 
!!    for (in particular for the double amplitudes in CCSD wavefunctions). Also,
!!    dt is expressed in the basis compatible with t.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: dt 
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: t 
!
      call daxpy(wf%n_gs_amplitudes, one, dt, 1, t, 1)
!
   end subroutine form_newton_raphson_t_estimate_ccs
!
!
   subroutine make_bath_orbital_ccs(wf)
!!
!!    Make bath orbital
!!    Written by Sarai D. Folkestad
!!
!!    Makes bath orbital with all corresponding integrals 
!!    zero
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_copy
      real(dp), dimension(:), allocatable :: orbital_energies_copy
!
      call mem%alloc(orbital_coefficients_copy, wf%n_ao, wf%n_mo)
!
      call copy_and_scale(one, wf%orbital_coefficients, orbital_coefficients_copy, wf%n_mo*wf%n_ao)
!
      call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
      call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo + wf%n_bath_orbitals)
      call zero_array(wf%orbital_coefficients, wf%n_ao*(wf%n_mo + wf%n_bath_orbitals))
!
      wf%orbital_coefficients(1:wf%n_ao, 1:wf%n_mo) = orbital_coefficients_copy(:,:)
!
      call mem%dealloc(orbital_coefficients_copy, wf%n_ao, wf%n_mo)
!
      call mem%alloc(orbital_energies_copy, wf%n_mo)
!
      call copy_and_scale(one, wf%orbital_energies, orbital_energies_copy, wf%n_mo)
!
      call mem%dealloc(wf%orbital_energies, wf%n_mo)
      call mem%alloc(wf%orbital_energies, wf%n_mo + wf%n_bath_orbitals)
      call zero_array(wf%orbital_energies, wf%n_mo + wf%n_bath_orbitals)
!
      wf%orbital_energies(1:wf%n_mo) = orbital_energies_copy(:)
!
      call mem%dealloc(orbital_energies_copy, wf%n_mo)
!
      wf%n_mo = wf%n_mo + wf%n_bath_orbitals
      wf%n_v = wf%n_v + wf%n_bath_orbitals
!
   end subroutine make_bath_orbital_ccs
!
!
   subroutine set_ip_start_indices_ccs(wf, start_indices)
!!
!!    Set IP start indices
!!    Written by Sarai D. Folkestad, Aug 2019
!!
!!    Sets IP start indices for trial vectors 
!!    from the orbital energies of the 
!!    occupied. (Koopmans)
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      integer, dimension(wf%n_singlet_states), intent(out) :: start_indices
!
      integer :: I
!
      do I = 1, wf%n_singlet_states
!
         start_indices(I) = (wf%n_o - I)*wf%n_v + wf%n_v
!
      enddo
!
   end subroutine set_ip_start_indices_ccs
!
!
   subroutine get_ip_projector_ccs(wf, projector)
!!
!!    Get IP projector 
!!    Written by Sarai D. Folkestad, Aug 2019
!!
!!    Constructs and returns the projector
!!    for an IP calculation (valence).
!!
!!    Only excitations into the last virtual orbital
!!    (the bath orbital) are allowed.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes),intent(out) :: projector
!
      integer :: A, I, AI
!
      call zero_array(projector, wf%n_es_amplitudes)
!
      do I = 1, wf%n_o
         do A = wf%n_v - wf%n_bath_orbitals + 1, wf%n_v
!
            AI = wf%n_v*(I-1) + A
            projector(AI) = one
!
         enddo
      enddo
!
   end subroutine get_ip_projector_ccs
!
!
   subroutine approximate_double_excitation_vectors_ccs(wf, R_ai, R_aibj, omega)
!!
!!    Construct approximate double excitation vectors
!!    Sarai D. Folkestad, May 2019
!!
!!    Approximate double excitation vectors:
!!
!!       R_aibj = 1/Δ_aibj * (P_ai,bj(sum_c R_ci g_bjac - sum_k R_bk g_kjai))/(-ε_aibj + ω^CCS) 
!!              = 1/Δ_aibj * X_aibj/(-ε_aibj + ω^CCS)
!!
!!    Used for CNTOs from CCS. 
!!
!!    For further information see Baudin, P. and Kristensen, K., J. Chem. Phys. 2017, 146, 214114
!!
!!    OBS! Renormalizes R !
!!
!!    OBS! Integrals are in MO basis not T1 basis
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out) :: R_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)               :: R_ai
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(:,:,:), allocatable :: L_Jai, L_Jkj, L_Jac
      real(dp), dimension(:,:,:), allocatable :: X_Jai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kjai
      real(dp), dimension(:), allocatable :: R_aibj_packed
!
      integer :: current_c_batch, req0, req1
!
      type(batching_index) :: batch_c
!
      real(dp) :: R_d_norm_sq, R_s_norm_sq, R_norm, ddot
!
      integer :: a, b, i, j
!
!     Read Cholesky vectors in MO basis
!
      call mem%alloc(L_Jai, wf%integrals%n_J, wf%n_v, wf%n_o)
      call mem%alloc(L_Jkj, wf%integrals%n_J, wf%n_o, wf%n_o)
!
      call wf%integrals%read_cholesky(L_Jai, wf%n_o + 1, wf%n_mo, 1, wf%n_o)
      call wf%integrals%read_cholesky(L_Jkj, 1, wf%n_o, 1, wf%n_o)
!
!     Construct g_kjai integrals in MO basis
!
      call mem%alloc(g_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',            &
                 (wf%n_o**2),         &
                 (wf%n_v)*(wf%n_o),   &
                 wf%integrals%n_J,    &
                 one,                 &
                 L_Jkj,               &
                 wf%integrals%n_J,    &
                 L_Jai,               &
                 wf%integrals%n_J,    &
                 zero,                &
                 g_kjai,              &
                 (wf%n_o**2))
!
      call mem%dealloc(L_Jkj, wf%integrals%n_J, wf%n_o, wf%n_o)
!
!     - sum_k R_bk g_kjai
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  wf%n_o,                 &
                  -one,                   &
                  R_ai,                   & ! R_bk
                  wf%n_v,                 &
                  g_kjai,                 &
                  wf%n_o,                 &
                  zero,                   &
                  R_aibj,                 & ! R_bjai but we will symmetrize anyhow
                  wf%n_v)
!
      call mem%dealloc(g_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_Jai, wf%integrals%n_J, wf%n_v, wf%n_o)
      call zero_array(X_Jai, (wf%integrals%n_J)*(wf%n_v)*(wf%n_o))
!
      req0 = 0
      req1 = (wf%n_v)*(wf%integrals%n_J)
!
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1)
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
         call mem%alloc(L_Jac, wf%integrals%n_J, wf%n_v, batch_c%length)
!
         call wf%integrals%read_cholesky(L_Jac, wf%n_o + 1, wf%n_mo, wf%n_o + batch_c%first, wf%n_o + batch_c%last)
!
!        X_ai_J = sum_c R_ci L_ac_J
!
         call dgemm('N', 'N',                   &
                     wf%n_v*(wf%integrals%n_J), &
                     wf%n_o,                    &
                     batch_c%length,            &
                     one,                       &
                     L_Jac,                     &
                     wf%n_v*(wf%integrals%n_J), &
                     R_ai(batch_c%first, 1),    & ! R_ci
                     wf%n_v,                    &
                     one,                       &
                     X_Jai,                     &
                     wf%n_v*(wf%integrals%n_J))
!
         call mem%dealloc(L_Jac, wf%integrals%n_J, wf%n_v, batch_c%length)
!
      enddo
!
!     sum_J X_ai_J L_bj_J
!
      call dgemm('T', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  wf%integrals%n_J,    &
                  one,                 &
                  X_Jai,               &
                  wf%integrals%n_J,    &
                  L_Jai,               & ! L_Jbj
                  wf%integrals%n_J,    &
                  one,                 &
                  R_aibj,              &
                  (wf%n_v)*(wf%n_o))

!
      call mem%dealloc(L_Jai, wf%integrals%n_J, wf%n_v, wf%n_o)
      call mem%dealloc(X_Jai, wf%integrals%n_J, wf%n_v, wf%n_o)
!
!     Symmetrize
!
      call symmetric_sum(R_aibj, wf%n_o*wf%n_v)
!
!     Binormalization factor
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            R_aibj(a, i, a, i) = half*R_aibj(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do
!
!     Divide by orbital differences and CCS excitation energy
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  R_aibj(a, i, b, j) = R_aibj(a, i, b, j)/(- wf%orbital_energies(a + wf%n_o) - &
                                                            wf%orbital_energies(b + wf%n_o) + &
                                                            wf%orbital_energies(i) + wf%orbital_energies(j) + omega)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Normalize full excitation vector with singles and doubles part.
!
      call mem%alloc(R_aibj_packed, wf%n_o*wf%n_v*(wf%n_o*wf%n_v+1)/2)
      call packin(R_aibj_packed, R_aibj, wf%n_o*wf%n_v)
!
      R_d_norm_sq = ddot(wf%n_o*wf%n_v*(wf%n_o*wf%n_v+1)/2, R_aibj_packed, 1, R_aibj_packed, 1)
      R_s_norm_sq = ddot(wf%n_o*wf%n_v, R_ai, 1, R_ai, 1)
!
      R_norm = sqrt(R_s_norm_sq + R_d_norm_sq)
!
      call dscal(wf%n_o*wf%n_v, one/R_norm, R_ai, 1)
      call dscal((wf%n_o*wf%n_v)**2, one/R_norm, R_aibj, 1)
!
      call mem%dealloc(R_aibj_packed, wf%n_o*wf%n_v*(wf%n_o*wf%n_v+1)/2)
!
   end subroutine approximate_double_excitation_vectors_ccs
!
!
   subroutine read_cvs_settings_ccs(wf)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      if (input%requested_keyword_in_section('core excitation', 'solver cc es')) then 
!  
!        Determine the number of core MOs 
!
         wf%n_core_MOs = input%get_n_elements_for_keyword_in_section('core excitation', 'solver cc es')
!
!        Then read the vector of core MOs
!
         call wf%initialize_core_MOs()
!
         call input%get_array_for_keyword_in_section('core excitation', 'solver cc es', wf%n_core_MOs, wf%core_MOs)
!
      else
!
         call output%error_msg('found no specified core MOs in input for CVS calculation')
!
      endif 
!
      wf%cvs = .true.
!
   end subroutine read_cvs_settings_ccs
!
!
   subroutine mo_preparations_ccs(wf)
!!
!!    MO preparations
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Routine which initializes the MO integral tool,
!!    MO transforms the Cholesky vectors, and does other
!!    preparations related to modifications of the MOs,
!!    such as frozen core, change of basis from canonical
!!    orbitals and shifting of bath orbitals (not implemented).
!!
!!    This routine is not overwritten for 
!!    descendants of standard CC-type (e.g., CCSD, CC2, CC3)
!!    but will be so for MLCC methods.
!!
      implicit none
!
      class(ccs) :: wf
!
      wf%integrals = mo_integral_tool(wf%n_o, wf%n_v, wf%system%n_J)
      call wf%construct_and_write_mo_cholesky(wf%n_mo, wf%orbital_coefficients, wf%integrals%cholesky_mo)
!
   end subroutine mo_preparations_ccs
!
!
   subroutine read_frozen_orbital_terms_ccs(wf)
!!
!!    Read frozen orbital contributions
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    Reads frozen orbital contributions to the
!!    mo Fock matrix.
!!
      implicit none
!
      class(ccs) :: wf
!
      if (wf%frozen_core) then
!      
         call wf%initialize_mo_fock_fc_term()
!
         wf%mo_fock_fc_file = sequential_file('MO_Fock_FC')
!
         call wf%mo_fock_fc_file%open_('read', 'rewind')
         call wf%mo_fock_fc_file%read_(wf%mo_fock_fc_term, wf%n_mo**2)
         call wf%mo_fock_fc_file%close_('keep')
!
      endif
!
      if (wf%frozen_hf_mos) then
!
         call wf%initialize_mo_fock_frozen_hf_term()
!
         wf%mo_fock_frozen_hf_file = sequential_file('MO_frozen_hf_Fock')
!
         call wf%mo_fock_frozen_hf_file%open_('read', 'rewind')
         call wf%mo_fock_frozen_hf_file%read_(wf%mo_fock_frozen_hf_term, wf%n_mo**2)
         call wf%mo_fock_frozen_hf_file%close_('keep')
!
      endif
!
   end subroutine read_frozen_orbital_terms_ccs
!
!
   subroutine check_for_degeneracies_ccs(wf, transformation, threshold)
!!
!!    Check for degeneracies in the excited states
!!    Written by Alexander C. Paul and Rolf H. Myhre, Oct 2019
!!
!!    The solver might not ensure orthogonality of the states (e.g. in DIIS).
!!    Therefore, degenerate roots might be found where the eigenvectors
!!    are equal up to a sign (and within the convergence threshold)
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(:,:), allocatable  :: R
!
      integer :: current_state, state
      integer :: n_degeneracy, p, q
      integer :: reduced_degeneracy
!
      current_state = 1
!
      do while (current_state .le. wf%n_singlet_states)
!
         n_degeneracy = 1
         state = current_state + 1
!
!        Check for degeneracies in the excitation energies
!        NB: energies are assumed to be sorted (cf. DIIS - quicksort)
!
         do while (state .le. wf%n_singlet_states)
!
            if (transformation .eq. 'right') then
!
               if (abs(wf%right_excitation_energies(state) - &
                       wf%right_excitation_energies(current_state)) .gt. threshold) exit
!
            else if (transformation .eq. 'left') then
!
               if (abs(wf%left_excitation_energies(state) - &
                       wf%left_excitation_energies(current_state)) .gt. threshold) exit
!
            else
!
               call output%error_msg('Tried to check for parallel states but argument ' &
                                     // trim(transformation) // ' not recognized.')
!
            end if
!
            n_degeneracy = n_degeneracy + 1
!
            state = state + 1
!
         end do
!
!        If degeneracies are found:
!           - read all degenerate states
!           - loop through all state pairs (p != q): 
!                 -If 2 states are parallel print warning 
!                  and reduce the degree of degeneracy
!
         if (n_degeneracy .gt. 1) then
!
            call mem%alloc(R, wf%n_es_amplitudes, n_degeneracy)
!
            do p = 1, n_degeneracy
               call wf%read_excited_state(R(:,p), current_state + p - 1, trim(transformation))
            end do
!
            reduced_degeneracy = n_degeneracy
!
            do p = 2, n_degeneracy ! if p == 1 then also q == 1 and we cycle anyway
               do q = 1, p
!
                  if(p .eq. q) cycle
!
                  if(are_vectors_parallel(R(:,p), R(:,q), wf%n_es_amplitudes, threshold)) then
!
                     call output%printf('Warning: The (a0) states (i0) and (i0) are parallel.',&
                                         ints=[current_state + p - 1, current_state + q - 1],  &
                                         chars=[trim(transformation)], fs='(//t3,a)', pl='m')
!
                     reduced_degeneracy = reduced_degeneracy - 1
!
                  end if
!
               end do
            end do
!
            if(reduced_degeneracy .gt. 1) then
!
               call output%printf('Found a (i0)-fold degeneracy involving the &
                                  &(i0). state.', fs='(//t6,a)', pl='v',      &
                                  ints=[reduced_degeneracy, current_state])
!
            end if
!
            call mem%dealloc(R, wf%n_es_amplitudes, n_degeneracy)
!
         end if
!
         current_state = current_state + n_degeneracy
!
      end do
!
   end subroutine check_for_degeneracies_ccs
!
!
   subroutine biorthonormalize_L_and_R_ccs(wf, energy_threshold, residual_threshold, skip_states)
!!
!!    Biorthonormalize the left and right eigenvectors
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    For that: 
!!       - check if left and right excitation energies are consistent
!!       - check for degenerate states
!!       - check for and discard parallel states
!!       - If degeneracies are present: biorthonormalize using Gram-Schmidt
!!          following Kohaupt, L., Rocky Mountain J. Math., 44, 1265, (2014)
!!       - Normalize and write to file
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), intent(in) :: energy_threshold, residual_threshold
!
      logical, dimension(wf%n_singlet_states), intent(inout) :: skip_states
!
      real(dp), dimension(:,:), allocatable :: R, R_normalized
      real(dp), dimension(:,:), allocatable :: L, L_normalized
!
      real(dp), dimension(:), allocatable :: overlap_LR
!
      integer, dimension(:), allocatable :: order
!
      logical, dimension(:), allocatable :: parallel
!
      real(dp) :: LT_R, overlap_L_Rnorm, overlap_Lnorm_R, norm_R
      real(dp) :: ddot
!
      integer :: state, current_state, prev_state, p, q
      integer :: n_degeneracy, reduced_degeneracy_r, reduced_degeneracy_l
      integer :: n_overlap_zero, state_nonzero_overlap
!
      logical :: biorthonormalize
!
!     Loop through states look for degenerate states (skip degeneracy)
!     discard parallel states
!     biorthonormalize degenerate states
!
      call output%printf('Biorthonormalization of left and right excited state vectors', &
                          pl='v', fs='(/t3,a)')
!
      current_state = 1
!
      do while (current_state .le. wf%n_singlet_states)
!
         skip_states(current_state) = .false.
!
         if (abs(wf%left_excitation_energies(current_state) &
               - wf%right_excitation_energies(current_state)) .lt. energy_threshold) then
!
            call output%printf('The left and right states corresponding to root (i0) &
                               &are consistent', fs='(/t6,a)', ints=[current_state], pl='v')
!
            n_degeneracy = 1
            state = current_state + 1
!
!           :: Check for degeneracies in both left and right excitation energies ::
!
            do while (state .le. wf%n_singlet_states)
!
               if(abs(wf%left_excitation_energies(state)          &
                    - wf%left_excitation_energies(current_state)) &
                  .gt. energy_threshold) exit
!
               if (abs(wf%right_excitation_energies(state)           &
                     - wf%right_excitation_energies(current_state))  &
                     .lt. energy_threshold) then
!
                  n_degeneracy = n_degeneracy + 1
!
               else 
!
                  call output%error_msg('Different degree of degeneracy in the   &
                  &    left excited states compared to the right excited states')
!
               end if
!
               state = state + 1
!
            end do
!
            call output%printf('Degree of degeneracy: (i0)', fs='(t6,a)', &
                                ints=[n_degeneracy], pl='debug')
!
            if(n_degeneracy .gt. 1) then
!
!              :: Check for parallel "right" states ::
!
               call mem%alloc(R, wf%n_es_amplitudes, n_degeneracy)
!
               do p = 1, n_degeneracy
                  call wf%read_excited_state(R(:,p), current_state + p - 1, 'right')
               end do
!
               reduced_degeneracy_r = n_degeneracy
!
               call mem%alloc(parallel, n_degeneracy)
!
               do p = 1, n_degeneracy
!
                  parallel(p) = .false.
!
                  do q = 1, p
!
                     if(q .eq. p) cycle
!
                     if(are_vectors_parallel(R(:,p), R(:,q), &
                        wf%n_es_amplitudes, residual_threshold)) then
!
                        parallel(p) = .true.
!
                        call output%printf('Warning: The right states (i0) and (i0) &
                                           &are parallel', pl='m', fs='(/t3,a)',   &
                                           ints=[current_state+p-1, current_state+q-1])
!
                        reduced_degeneracy_r = reduced_degeneracy_r - 1
!
                     end if
!
                  end do
               end do
!
!              Remove parallel state from array R
               if(reduced_degeneracy_r .ne. n_degeneracy) then
!
                  state = 0
!
                  do p = 1, n_degeneracy
                     if(parallel(p)) cycle
!
                     state = state + 1
                     call dcopy(wf%n_es_amplitudes, &
                                R(:,p), 1,          &
                                R(:,state), 1)
!
                  end do
!
               end if
!
!              :: Check for parallel "left" states ::
!
               call mem%alloc(L, wf%n_es_amplitudes, n_degeneracy)
!
               do p = 1, n_degeneracy
                  call wf%read_excited_state(L(:,p), current_state + p - 1, 'left')
               end do
!
               reduced_degeneracy_l = n_degeneracy
!
               do p = 1, n_degeneracy
!
                  parallel(p) = .false.
                  skip_states(current_state + p - 1) = .false.
!
                  do q = 1, p
!
                     if(q .eq. p) cycle
!
                     if(are_vectors_parallel(L(:,p), L(:,q), &
                        wf%n_es_amplitudes, residual_threshold)) then
!
                        parallel(p) = .true.
!
                        call output%printf('Warning: The left states (i0) and (i0) &
                                           &are parallel', pl='m', fs='(/t3,a)',   &
                                           ints=[current_state + p - 1, current_state+q-1])
!
                        reduced_degeneracy_l = reduced_degeneracy_l - 1
!
                     end if
!
                  end do
               end do
!
               if(reduced_degeneracy_l .ne. reduced_degeneracy_r) then
                  call output%error_msg('Different degree of degeneracy in the &
                                        &left excited states compared to the &
                                        &right excited states')
               end if
!
!              Remove parallel state from array L
               if(reduced_degeneracy_l .ne. n_degeneracy) then
!
                  state = 0
!
                  do p = 1, n_degeneracy
!
                     if(p .gt. reduced_degeneracy_l) then
                        skip_states(current_state + p - 1) = .true.
                     end if
!
                     if(parallel(p)) cycle
!
                     state = state + 1
!
                     call dcopy(wf%n_es_amplitudes, &
                                L(:,p), 1,          &
                                L(:,state), 1)
!
                  end do
!
               end if
!
               if (reduced_degeneracy_r .gt. 1) then
!
                  call output%printf('Found a degeneracy between:', fs='(/t6,a)', pl='n')
                  call output%print_separator('n', 29,'-', fs='(t6,a)')
                  call output%printf('State     Excitation Energy', fs='(t6,a)', pl='n')
!
                  do p = 1, n_degeneracy
!
                     if(parallel(p)) cycle
!
                     call output%printf(' (i2)     (f19.12)', &
                          fs='(t6,a)', pl='n',                &
                          ints=[current_state + p - 1],       & 
                          reals=[wf%right_excitation_energies(current_state+p-1)])
!
                  end do
!
               end if
!
               call mem%dealloc(parallel, n_degeneracy)
!
!              :: Biorthonormalize states ::
!              -----------------------------
!
!              (following Kohaupt, L., Rocky Mountain J. Math., 44, 1265, (2014))
!
!              non degenerate L/R states are biorthogonal, thus only normalization needed
!
!              degenerate L/R states should be biorthogonal as well
!              if not the k-th state is determined in terms of the 
!              previously biorthogonalized states i
!
!              Intermediate:  L'(k) = L(k) - sum_i < L(k)|R"(i)> * L"(i)
!              Biorthonormal: L"(k) = < L'(k)|R(k)>^(-1) * L'(k)
!
!              Biorthonormal: R"(k) = R(k) - sum_i < R(k)|L"(i)> * R"(i)
!
!              Renormalize R"(k) afterwards and then binormalize L"(k) to R"(k)
!
               call mem%alloc(L_normalized, wf%n_es_amplitudes, reduced_degeneracy_l)
               call mem%alloc(R_normalized, wf%n_es_amplitudes, reduced_degeneracy_r)
!
               do p = 1, reduced_degeneracy_l
!
                  call dcopy(wf%n_es_amplitudes, &
                             L(:,p), 1,          &
                             L_normalized(:,p), 1)
!
               end do
!
               call mem%alloc(overlap_LR, reduced_degeneracy_r)
               call mem%alloc(order, reduced_degeneracy_r)
!
!              We first check if the degenerate states are biorthogonal
!              then we only need to binormalize (biorthonormalize = .false.).
!              If we need to biorthonormalize the logical is set to .true.
!              so that we don't run into the wrong branch of the if-statement
!
               biorthonormalize = .false.
!
               do state = 1, reduced_degeneracy_l
!
                  n_overlap_zero = 0
                  state_nonzero_overlap = 0
!
                  do q = 1, reduced_degeneracy_r
!
!                    Overlap of Singles and Doubles should be enough also for CC3
                     overlap_lr(q) = abs(ddot(wf%n_es_amplitudes, &
                                              L(:,state), 1,      &
                                              R(:,q), 1))
!
                     if(overlap_lr(q) .lt. residual_threshold) then
!
                        n_overlap_zero = n_overlap_zero + 1
!
                     else
!
                        state_nonzero_overlap = q
!
                     end if
!
                  end do
!
!                 Simple binormalization if only 1 right state 
!                 has significant overlap with the left state
!
                  if((.not. biorthonormalize) .and. &
                     (n_overlap_zero .eq. reduced_degeneracy_r-1)) then
!
                     call output%printf('Degenerate states except one are biorthogonal &
                                       &- Thus, only binormalize', pl='v', fs='(/t6,a)')
!
                     LT_R = wf%L_R_overlap(L(:, state),                             &
                                           current_state + state - 1,               &
                                           R(:, state_nonzero_overlap),             &
                                           current_state + state_nonzero_overlap - 1)
!
!                    Sanity check that the left and corresponding right state are not orthogonal
!
                     if(abs(LT_R) .lt. 1.0d-2) then
!
                        call output%printf('Warning: Overlap of (i0). left and &
                                          &right state close to zero.',        &
                                           pl='m', ints=[current_state])
!
                     end if
!
!                    Normalize the new left state to the right state
                     call dscal(wf%n_es_amplitudes,     &
                                one/LT_R,               &
                                L_normalized(:, state), & 
                                1)
!
!                    Copy to have correct ordering in R_normalized
                     call dcopy(wf%n_es_amplitudes,            &
                                R(:,state_nonzero_overlap), 1, &
                                R_normalized(:, state), 1)
!
                  else ! Actual biorthonormalization needed
!
                     call output%printf('Biorthonormalization of degenerate states', &
                                         pl='n',fs='(/t6,a)')
!
!                    Make sure to not run into the other branch of the if statement
                     biorthonormalize = .true.
!
!                    sort overlap_lr and select corresponding R state for L(p)
!                    R(:,order(1)) has the maximal overlap with L(p)
!
                     call quicksort_with_index_descending(overlap_lr, &
                                                          order,      &
                                                          reduced_degeneracy_r)
!
                     call dcopy(wf%n_es_amplitudes,       &
                                R(:,order(1)), 1,         &
                                R_normalized(:,state), 1)
!
                     do prev_state = 1, state - 1
!
!                       :: Construct biorthogonal left state ::
!
!                       L'(k) = L(k) - sum_i < L(k)|R"(i)> * L"(i)
!
                        overlap_L_Rnorm = ddot(wf%n_es_amplitudes,         &
                                               L(:,state),                 &
                                               1,                          &
                                               R_normalized(:,prev_state), &
                                               1)
!
                        call daxpy(wf%n_es_amplitudes,         &
                                   -overlap_L_Rnorm,           &
                                   L_normalized(:,prev_state), &
                                   1,                          &
                                   L_normalized(:,state),      &
                                   1)
!
!                       :: Construct biorthonormal right state ::
!
!                       R"(k) = R(k) - sum_i < R(k)|L"(i)> * R"(i)
!
                        overlap_Lnorm_R = ddot(wf%n_es_amplitudes,          &
                                                R(:,order(1)),              &
                                                1,                          &
                                                L_normalized(:,prev_state), &
                                                1)
!
                        call daxpy(wf%n_es_amplitudes,         &
                                   -overlap_Lnorm_R,           &
                                   R_normalized(:,prev_state), &
                                   1,                          &
                                   R_normalized(:,state),      &
                                   1)
!
                     end do
!
!                    :: Binormalize L'(state) to the corresponding R-state ::
!
                     overlap_Lnorm_R = ddot(wf%n_es_amplitudes,    &
                                            L_normalized(:,state), &
                                            1,                     &
                                            R(:,order(1)),         &
                                            1)
!
                     call dscal(wf%n_es_amplitudes, one/overlap_Lnorm_R, L_normalized(:,state), 1)
!
!                    Zero out R state that has already been used
!                    Thus, overlap_lr(q) will be zero and this state will not be selected again
!
                     call zero_array(R(:,order(1)), wf%n_es_amplitudes)
!
!                    :: Renormalize singles and doubles part of the right vectors ::
!
                     norm_R = ddot(wf%n_es_amplitudes,    &
                                   R_normalized(:,state), &
                                   1,                     &
                                   R_normalized(:,state), &
                                   1)
!
                     call dscal(wf%n_es_amplitudes, one/norm_R, R_normalized(:,state), 1)
!
!                    :: Binormalize the left to the right vectors ::
!                    ::        including triples if present       ::
!
                     LT_R = wf%L_R_overlap(L_normalized(:,state),     &
                                           current_state + state - 1, &
                                           R_normalized(:,state),     &
                                           current_state + state - 1)
!
!                    Sanity check that the left and corresponding right state are not orthogonal
!
                     if(abs(LT_R) .lt. 1.0d-2) then
!
                        call output%printf('Warning: Overlap of (i0). left and right state close to zero.', &
                                          pl='m', ints=[current_state])
!
                     else if(abs(LT_R) .lt. residual_threshold) then
!
                        call output%printf('Overlap of (i0). left and right state less than &
                                          &threshold: (e8.3).', reals=[residual_threshold], &
                                          pl='m', ints=[current_state])
!
                        call output%error_msg('Trying to binormalize nonoverlapping states.')
!
                     end if
!
                     call dscal(wf%n_es_amplitudes, one/LT_R, L_normalized(:, state), 1)
!
                  end if
!
                  call wf%save_excited_state(L_normalized(:, state), current_state+state-1, 'left')
                  call wf%save_excited_state(R_normalized(:, state), current_state+state-1, 'right')
!
               end do
!
               call mem%dealloc(R_normalized, wf%n_es_amplitudes, reduced_degeneracy_r)
               call mem%dealloc(L_normalized, wf%n_es_amplitudes, reduced_degeneracy_l)
!
               call mem%dealloc(R, wf%n_es_amplitudes, n_degeneracy)
               call mem%dealloc(L, wf%n_es_amplitudes, n_degeneracy)
!
               call mem%dealloc(overlap_LR, reduced_degeneracy_r)
               call mem%dealloc(order, reduced_degeneracy_r)
!
            else ! States are not degenerate, thus only binormalize
!
               call mem%alloc(R, wf%n_es_amplitudes, 1)
               call wf%read_excited_state(R, current_state, 'right')
!
               call mem%alloc(L, wf%n_es_amplitudes, 1)
               call wf%read_excited_state(L, current_state, 'left')
!
               LT_R = wf%L_R_overlap(L,             &
                                     current_state, &
                                     R,             &
                                     current_state)
!
!              Sanity check that the left and corresponding right state are not orthogonal
!
               if(abs(LT_R) .lt. 1.0d-2) then
!
                  call output%printf('Warning: Overlap of (i0). left and right state close to zero.', &
                                     pl='m', ints=[current_state])
!
               else if(abs(LT_R) .lt. residual_threshold) then
!
                  call output%printf('Overlap of (i0). left and right state less than &
                                    &threshold: (e8.3).', reals=[residual_threshold], &
                                     pl='m', ints=[current_state])
!
                  call output%error_msg('Trying to binormalize biorthogonal states.')
!
               end if
!
               call mem%dealloc(R, wf%n_es_amplitudes, 1)
!
!              Normalize the new left state to the right state
               call dscal(wf%n_es_amplitudes, 1/LT_R, L, 1)
!
               call wf%save_excited_state(L, current_state, 'left')
!
               call mem%dealloc(L, wf%n_es_amplitudes, 1)
!
            end if
!
         else ! Sanity check failed - roots ordered incorrectly
!
            call output%printf('Eigenvector (i0) is not left-right consistent &
                              & to threshold (e8.3).', fs='(/t6,a)', pl='m',  &
                                ints=[current_state], reals=[energy_threshold])
!
            call output%printf('Energies (left, right): (f19.12) (f19.12)', &
                 reals=[wf%left_excitation_energies(current_state),         &
                 wf%right_excitation_energies(current_state)], fs='(/t6,a)',&
                 pl='m')
!
            call output%error_msg('while biorthonormalizing.')
!
         end if
!
         current_state = current_state + n_degeneracy
!
      end do
!
   end subroutine biorthonormalize_L_and_R_ccs
!
!
   subroutine read_mm_fock_contributions_ccs(wf)
!!
!!    Read MM Fock contributions
!!    Written by Tommaso Giovannini, Oct 2019
!!
!!    Reads HF MM Fock contributions at convergence
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%initialize_mm_matrices()
!      
      if(wf%system%mm%forcefield.eq.'non-polarizable') &
         call mem%alloc(wf%nopol_h_wx,wf%n_ao,wf%n_ao)
!
      call wf%read_mm_matrices()
!
!
   end subroutine read_mm_fock_contributions_ccs
!
!
   subroutine read_pcm_fock_contributions_ccs(wf)
!!
!!    Read MM Fock contributions
!!    Written by Tommaso Giovannini, Oct 2019
!!
!!    Reads HF MM Fock contributions at convergence
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%initialize_pcm_matrices()
!
      call wf%read_pcm_matrices()
!
!
   end subroutine read_pcm_fock_contributions_ccs
!
!
end module ccs_class
