!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
   use parameters
!
   use wavefunction_class, only: wavefunction
!
   use global_in, only: input
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
!
   use stream_file_class, only: stream_file
!
   use range_class, only: range_
!
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
   use abstract_eri_cholesky_c_class, only: abstract_eri_cholesky_c
   use eri_adapter_class, only: eri_adapter
   use eri_adapter_c_class, only: eri_adapter_c
!
!  Need these here to avoid intel segfault
   use eri_memory_tool_c_class, only: eri_memory_tool_c
   use eri_tool_c_class, only: eri_tool_c
   use eri_cholesky_disk_c_class, only: eri_cholesky_disk_c
   use eri_cholesky_memory_c_class, only: eri_cholesky_memory_c
!
   use array_utilities, only: zdot
!
   use batching_index_class, only: batching_index
!
   implicit none
!
   type, extends(wavefunction) :: ccs
!
      real(dp)    :: correlation_energy
      complex(dp) :: correlation_energy_complex
!
      integer :: n_gs_amplitudes
      integer :: n_es_amplitudes
      integer :: n_t1
!
      integer :: n_singlet_states
      integer :: n_bath_orbitals
!
      integer :: n_core_MOs
      logical :: cvs, rm_core
!
      real(dp), dimension(:), allocatable :: left_excitation_energies
      real(dp), dimension(:), allocatable :: right_excitation_energies
!
      logical :: bath_orbital
!
      logical :: need_g_abcd
!
      type(stream_file) :: t_file, tbar_file
!
      type(stream_file), dimension(:), allocatable :: l_files, r_files
!
      type(eri_adapter), allocatable :: eri_t1
      type(eri_adapter_c), allocatable :: eri_t1_c
!
      class(abstract_eri_cholesky), allocatable :: L_mo
      class(abstract_eri_cholesky), allocatable :: L_t1
!
      class(abstract_eri_cholesky_c), allocatable :: L_mo_c
      class(abstract_eri_cholesky_c), allocatable :: L_t1_c
!
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

      real(dp),    dimension(:),   allocatable :: r0
      complex(dp), dimension(:),   allocatable :: r0_complex
!
      integer, dimension(:), allocatable :: core_MOs
!
   contains
!
      procedure :: print_banner                                  => print_banner_ccs
      procedure :: print_amplitude_info                          => print_amplitude_info_ccs
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
      procedure :: destruct_fock                                 => destruct_fock_ccs
      procedure :: destruct_fock_complex                         => destruct_fock_ccs_complex
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
      procedure :: initialize_density_intermediates &
                => initialize_density_intermediates_ccs
      procedure :: destruct_density_intermediates   &
                => destruct_density_intermediates_ccs
!
      procedure :: initialize_right_excitation_energies          => initialize_right_excitation_energies_ccs
      procedure :: destruct_right_excitation_energies            => destruct_right_excitation_energies_ccs
!
      procedure :: initialize_left_excitation_energies           => initialize_left_excitation_energies_ccs
      procedure :: destruct_left_excitation_energies             => destruct_left_excitation_energies_ccs
!
      procedure :: initialize_core_MOs                           => initialize_core_MOs_ccs
      procedure :: destruct_core_MOs                             => destruct_core_MOs_ccs
!
!     File handling procedures
!
      procedure :: initialize_files                              => initialize_files_ccs
      procedure :: initialize_ground_state_files                 => initialize_ground_state_files_ccs
      procedure :: initialize_excited_state_files                => initialize_excited_state_files_ccs
!
      procedure :: read_settings                                 => read_settings_ccs
      procedure :: read_cvs_settings                             => read_cvs_settings_ccs
      procedure :: read_rm_core_settings                         => read_rm_core_settings_ccs
!
      procedure :: read_singles_vector                           => read_singles_vector_ccs
      procedure :: save_singles_vector                           => save_singles_vector_ccs
!
      procedure :: save_amplitudes                               => save_amplitudes_ccs
      procedure :: read_amplitudes                               => read_amplitudes_ccs
      procedure :: save_multipliers                              => save_multipliers_ccs
      procedure :: read_multipliers                              => read_multipliers_ccs
      procedure :: save_excited_state                            => save_excited_state_ccs
      procedure :: read_excited_state                            => read_excited_state_ccs
      procedure :: read_excitation_vector_file                   => read_excitation_vector_file_ccs
      procedure :: save_excitation_vector_on_file                => save_excitation_vector_on_file_ccs
      procedure :: check_and_get_restart_vector                  => check_and_get_restart_vector_ccs
      procedure :: get_restart_vector                            => get_restart_vector_ccs
!
      procedure :: save_tbar_intermediates                       => save_tbar_intermediates_ccs
!
      procedure :: get_density_for_plotting &
                => get_density_for_plotting_ccs
!
!     Print summaries
!
      procedure :: print_gs_summary                              => print_gs_summary_ccs
      procedure :: print_X1_diagnostics                          => print_X1_diagnostics_ccs
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
      procedure :: set_excitation_energies                       => set_excitation_energies_ccs
!
!     Procedures related to the Fock matrix
!
      procedure :: construct_fock                                => construct_fock_ccs
      procedure :: construct_fock_complex                        => construct_fock_ccs_complex
!
      procedure :: add_frozen_fock_terms                         => add_frozen_fock_terms_ccs
      procedure :: add_frozen_fock_terms_complex                 => add_frozen_fock_terms_ccs_complex
!
      procedure :: add_t1_fock_length_dipole_term                => add_t1_fock_length_dipole_term_ccs
      procedure :: add_t1_fock_length_dipole_term_complex        => add_t1_fock_length_dipole_term_ccs_complex
!
      procedure :: construct_fock_ai_t1                          => construct_fock_ai_t1_ccs
      procedure :: construct_fock_ia_t1                          => construct_fock_ia_t1_ccs
      procedure :: construct_fock_ij_t1                          => construct_fock_ij_t1_ccs
      procedure :: construct_fock_ab_t1                          => construct_fock_ab_t1_ccs
!
      procedure :: construct_fock_ai_t1_complex                  => construct_fock_ai_t1_ccs_complex
      procedure :: construct_fock_ia_t1_complex                  => construct_fock_ia_t1_ccs_complex
      procedure :: construct_fock_ij_t1_complex                  => construct_fock_ij_t1_ccs_complex
      procedure :: construct_fock_ab_t1_complex                  => construct_fock_ab_t1_ccs_complex
!
!     Procedures related to the omega vector
!
      procedure :: construct_omega                               => construct_omega_ccs
      procedure :: construct_omega_complex                       => construct_omega_ccs_complex
!
      procedure :: omega_ccs_a1                                  => omega_ccs_a1_ccs
      procedure :: omega_ccs_a1_complex                          => omega_ccs_a1_ccs_complex
!
!     Routines related to the F transformation
!
      procedure :: F_transformation                              => F_transformation_ccs
      procedure :: F_x_transformation                            => F_x_transformation_ccs
      procedure :: F_x_mu_transformation                         => F_x_mu_transformation_ccs
!
      procedure :: F_ccs_a1_0                                    => F_ccs_a1_0_ccs
      procedure :: F_ccs_a1_1                                    => F_ccs_a1_1_ccs
      procedure :: F_ccs_b1_1                                    => F_ccs_b1_1_ccs
      procedure :: F_ccs_c1_1                                    => F_ccs_c1_1_ccs
!
!     Procedures related to the multiplier equation vector
!
      procedure :: construct_multiplier_equation                 => construct_multiplier_equation_ccs
      procedure :: construct_multiplier_equation_complex         => construct_multiplier_equation_ccs_complex
!
      procedure :: construct_eta                                 => construct_eta_ccs
      procedure :: construct_eta_complex                         => construct_eta_ccs_complex
!
!     Procedures related to the multiplier solver
!
      procedure :: print_banner_davidson_cc_multipliers          => print_banner_davidson_cc_multipliers_ccs
      procedure :: get_initial_cc_multipliers                    => get_initial_cc_multipliers_ccs
      procedure :: get_cc_multipliers_preconditioner             => get_cc_multipliers_preconditioner_ccs
      procedure :: cc_multipliers_summary                        => cc_multipliers_summary_ccs
!
!     Procedures related to the Jacobian transformation
!
      procedure :: construct_Jacobian_transform                  => construct_Jacobian_transform_ccs
!
      procedure :: prepare_for_Jacobians                 &
                => prepare_for_Jacobians_ccs
!
      procedure :: prepare_for_approximate_Jacobians     &
                => prepare_for_approximate_Jacobians_ccs
!
      procedure :: approximate_Jacobian_transform        &
                => approximate_Jacobian_transform_ccs
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
      procedure :: compute_eom_transition_moments                => compute_eom_transition_moments_ccs
!
      procedure :: construct_gs_density                          => construct_gs_density_ccs
      procedure :: construct_gs_density_complex                  => construct_gs_density_ccs_complex
      procedure :: mu_ref_density_terms                          => mu_ref_density_terms_ccs
      procedure :: mu_ref_density_terms_complex                  => mu_ref_density_terms_ccs_complex
!
      procedure :: construct_left_transition_density             => construct_left_transition_density_ccs
!
      procedure :: construct_right_transition_density            => construct_right_transition_density_ccs
      procedure :: construct_es_density                          => construct_es_density_ccs
      procedure :: mu_nu_density_terms                           => mu_nu_density_terms_ccs
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
      procedure :: construct_xiX                                 => construct_xiX_ccs
      procedure :: xiX_ccs_a1                                    => xiX_ccs_a1_ccs
      procedure :: etaX_eom_a                                    => etaX_eom_a_ccs
      procedure :: calculate_lr_transition_strength              => calculate_lr_transition_strength_ccs
!
!     Routines related to post-processing excited states
!
      procedure :: biorthonormalize_L_and_R                      => biorthonormalize_L_and_R_ccs
      procedure :: L_R_overlap                                   => L_R_overlap_ccs
      procedure :: get_degenerate_states                         => get_degenerate_states_ccs
      procedure :: check_for_parallel_states                     => check_for_parallel_states_ccs
      procedure :: remove_parallel_states                        => remove_parallel_states_ccs
      procedure :: remove_parallel_states_from_file              => remove_parallel_states_from_file_ccs
!
!     One-electron integrals
!
      procedure :: t1_transform                                  => t1_transform_ccs
      procedure :: t1_transform_complex                          => t1_transform_ccs_complex
!
      procedure :: add_t1_terms                                  => add_t1_terms_ccs
      procedure :: add_t1_terms_complex                          => add_t1_terms_ccs_complex
      procedure :: add_t1_terms_and_transform                    => add_t1_terms_and_transform_ccs
      procedure :: add_t1_terms_and_transform_complex            => add_t1_terms_and_transform_ccs_complex
!
!     One-electron integrals
!
      procedure :: get_t1_oei                                    => get_t1_oei_ccs
      procedure :: get_t1_oei_complex                            => get_t1_oei_ccs_complex
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
      procedure :: prepare_for_properties                           => prepare_for_properties_ccs
      procedure :: prepare_for_properties_complex                   => prepare_for_properties_ccs_complex
!
      procedure :: approximate_double_excitation_vectors         => approximate_double_excitation_vectors_ccs
!
!     Frozen core
!
      procedure :: construct_t1_frozen_fock_terms                => construct_t1_frozen_fock_terms_ccs
      procedure :: construct_t1_frozen_fock_terms_complex        => construct_t1_frozen_fock_terms_ccs_complex
!
!     MO preparations
!
      procedure :: mo_preparations                               => mo_preparations_ccs
      procedure :: construct_MO_screening_for_cd                 => construct_MO_screening_for_cd_ccs
!
      procedure :: integral_preparations                         => integral_preparations_ccs
      procedure :: integral_preparations_complex                 => integral_preparations_ccs_complex
!
!     Core-valence separation procedures
!
      procedure :: get_cvs_projector                             => get_cvs_projector_ccs
      procedure :: set_cvs_start_indices                         => set_cvs_start_indices_ccs
!
      procedure :: get_rm_core_projector                         => get_rm_core_projector_ccs
!
!     Other procedures
!
      procedure :: cleanup                                       => cleanup_ccs
      procedure :: general_cc_preparations                       => general_cc_preparations_ccs
      procedure :: set_variables_from_template_wf                => set_variables_from_template_wf_ccs
!
      procedure :: set_initial_amplitudes_guess                  => set_initial_amplitudes_guess_ccs
      procedure :: set_initial_multipliers_guess                 => set_initial_multipliers_guess_ccs
      procedure :: set_ip_start_indices                          => set_ip_start_indices_ccs
      procedure :: get_ip_projector                              => get_ip_projector_ccs
!
      procedure :: scale_amplitudes                              => scale_amplitudes_ccs
!
      procedure :: print_dominant_x_amplitudes                   => print_dominant_x_amplitudes_ccs
      procedure :: print_dominant_amplitudes                     => print_dominant_amplitudes_ccs
!
      procedure, private :: print_dominant_x1
!
      procedure :: make_bath_orbital                             => make_bath_orbital_ccs
!
!     Procedures related to time dependency
!
      procedure :: make_complex                                  => make_complex_ccs
      procedure :: make_ccs_complex                              => make_ccs_complex_ccs
!
      procedure :: cleanup_complex                               => cleanup_complex_ccs
      procedure :: cleanup_ccs_complex                           => cleanup_ccs_complex_ccs
!
      procedure, private :: complexify_cholesky_vectors
!
      procedure :: construct_complex_time_derivative             => construct_complex_time_derivative_ccs
!
      procedure :: construct_complex_time_derivative_amplitudes  => construct_complex_time_derivative_amplitudes_ccs
      procedure :: construct_complex_time_derivative_multipliers => construct_complex_time_derivative_multipliers_ccs
!
!     Initialize wavefunction
!
      procedure :: initialize => initialize_ccs
!
      procedure :: construct_ntos &
                => construct_ntos_ccs
      procedure :: construct_ntos_or_cntos &
                => construct_ntos_or_cntos_ccs
!
      procedure :: construct_t1_cholesky => construct_t1_cholesky_ccs
      procedure :: construct_t1_cholesky_complex => construct_t1_cholesky_ccs_complex
!
      procedure :: construct_cholesky_t1_oo => construct_cholesky_t1_oo_ccs
      procedure :: construct_cholesky_t1_oo_complex => construct_cholesky_t1_oo_ccs_complex
!
      procedure :: construct_cholesky_t1_vo  => construct_cholesky_t1_vo_ccs
      procedure :: construct_cholesky_t1_vo_complex  => construct_cholesky_t1_vo_ccs_complex
!
      procedure :: construct_cholesky_t1_vv  => construct_cholesky_t1_vv_ccs
      procedure :: construct_cholesky_t1_vv_complex  => construct_cholesky_t1_vv_ccs_complex
!
      procedure, private :: integral_factory
      procedure, private :: integral_factory_complex
      procedure, private :: cholesky_factory
      procedure, private :: cholesky_factory_complex
!
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
      include "F_ccs_interface.F90"
      include "jacobian_ccs_interface.F90"
      include "jacobian_transpose_ccs_interface.F90"
      include "mean_value_ccs_interface.F90"
      include "response_ccs_interface.F90"
      include "oei_ccs_interface.F90"
      include "t1_ccs_interface.F90"
      include "fock_ccs_interface.F90"
!
      include "complex_ccs_interface.F90"
!
      include "generated_complex_files/fock_ccs_complex_interface.F90"
      include "generated_complex_files/initialize_destruct_ccs_complex_interface.F90"
      include "generated_complex_files/jacobian_transpose_ccs_complex_interface.F90"
      include "generated_complex_files/multiplier_equation_ccs_complex_interface.F90"
      include "generated_complex_files/omega_ccs_complex_interface.F90"
      include "generated_complex_files/set_get_ccs_complex_interface.F90"
      include "generated_complex_files/t1_ccs_complex_interface.F90"
      include "generated_complex_files/mean_value_ccs_complex_interface.F90"
!
   end interface
!
contains
!
!
   subroutine initialize_ccs(wf, template_wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'ccs'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%n_t1            = wf%n_o*wf%n_v
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
      wf%need_g_abcd     = .false.
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
   end subroutine initialize_ccs
!
!
   subroutine general_cc_preparations_ccs(wf)
!!
!!    General CC preparations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
      implicit none
!
      class(ccs) :: wf
!
      call wf%initialize_files()
!
      wf%bath_orbital = .false.
      wf%need_g_abcd  = .false.
      wf%cvs          = .false.
      wf%rm_core      = .false.
!
      call wf%read_settings()
!
   end subroutine general_cc_preparations_ccs
!
!
   subroutine print_banner_ccs(wf)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad, Dec 2019
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      character(len=200) :: name_
!
      name_ = trim(convert_to_uppercase(wf%name_)) // ' wavefunction'
!
      call output%printf('m', ':: (a0)', chars=[name_], fs='(//t3,a)')
      call output%print_separator('m', len_trim(name_) + 6, '=')
!
!     Print settings
!
      call output%printf('m', 'Bath orbital(s):         (l0)', &
            logs=[wf%bath_orbital], fs='(/t6, a)')
!
!     Print orbital space info for cc
!
      call output%printf('m', ' - Number of orbitals:', &
                         fs='(/t3,a)')
      call output%printf('m', 'Occupied orbitals:    (i0)', &
                         ints=[wf%n_o], fs='(/t6,a)')
      call output%printf('m', 'Virtual orbitals:     (i0)', &
                         ints=[wf%n_v], fs='(t6,a)')
      call output%printf('m', 'Molecular orbitals:   (i0)', &
                         ints=[wf%n_mo], fs='(t6,a)')
      call output%printf('m', 'Atomic orbitals:      (i0)', &
                         ints=[wf%ao%n], fs='(t6,a)')
!
   end subroutine print_banner_ccs
!
!
   subroutine print_amplitude_info_ccs(wf)
!!
!!    Print amplitude info
!!    Written by Sarai D. Folkestad, Dec 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
!     Print orbital space info for cc
!
      call output%printf('m', ' - Number of ground state amplitudes:', fs='(/t3,a)')
!
      call output%printf('m', 'Single excitation amplitudes:  (i0)', &
            ints=[wf%n_t1], fs='(/t6,a)')

!
   end subroutine print_amplitude_info_ccs
!
!
   subroutine set_variables_from_template_wf_ccs(wf, template_wf)
!!
!!    Set variables from template wavefunction
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Transfers the wavefunction variables needed to start a
!!    CC calculation from either an HF of a CC wavefunction
!!
      implicit none
!
      class(ccs) :: wf

      class(wavefunction), intent(in) :: template_wf
!
      call wf%prepare_ao_tool(template = template_wf%ao)
      call wf%prepare_embedding()
!
      wf%n_o = template_wf%n_o
      wf%n_v = template_wf%n_v
!
      if (wf%n_v == 0) then
         call output%error_msg('Number of virtual orbitals is zero. &
                               &Cannot run post-HF methods')
      end if
!
      wf%ao%n = template_wf%ao%n
      wf%n_mo = template_wf%n_mo
!
#ifdef HAS_32BIT_INTEGERS
      if (wf%n_mo .gt. int32_mo_limit) call output%error_msg('too many MOs for 32bit integer')
#endif
!
      wf%hf_energy = template_wf%hf_energy
!
!     Set orbital coefficients and energies
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call dcopy(wf%n_mo, template_wf%orbital_energies, 1, wf%orbital_energies, 1)
!
      call dcopy(wf%ao%n*wf%n_mo, template_wf%orbital_coefficients, 1, wf%orbital_coefficients, 1)
!
!     Handle changes in the number of MOs as a result of
!     special methods
!
      if (wf%bath_orbital) call wf%make_bath_orbital()
!
      wf%exists_frozen_fock_terms = template_wf%exists_frozen_fock_terms
!
      if (wf%exists_frozen_fock_terms) then
!
         call wf%initialize_mo_fock_frozen()
!
         call dcopy(wf%n_mo**2, template_wf%mo_fock_frozen, 1, wf%mo_fock_frozen, 1)
!
         call wf%initialize_frozen_CCT()
!
         call dcopy(wf%ao%n**2, template_wf%frozen_CCT, 1, wf%frozen_CCT, 1)
!
         wf%frozen_dipole = template_wf%frozen_dipole
         wf%frozen_quadrupole = template_wf%frozen_quadrupole
!
      endif
!
   end subroutine set_variables_from_template_wf_ccs
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
      call wf%destruct_mo_fock_frozen()
      call wf%destruct_fock()
      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
      call wf%destruct_gs_density()
      call wf%destruct_frozen_CCT()
!
      call wf%destruct_core_MOs()
!
      deallocate(wf%ao)
      if (wf%embedded) deallocate(wf%embedding)
!
      if (allocated(wf%L_t1)) call wf%L_t1%remove_observer('eri t1')
!
!     To avoid memory leaks with intel, explicit deallocations
      deallocate(wf%eri_t1)
      deallocate(wf%L_mo)
      deallocate(wf%L_t1)
!
      call output%printf('v', '- Cleaning up ' // trim(wf%name_) // ' wavefunction', &
                         fs='(/t3,a)')
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
      if (input%is_section_present('cc')) then
!
         if (input%is_keyword_present('bath orbital','cc')) then
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
   subroutine set_initial_amplitudes_guess_ccs(wf, restart)
!!
!!    Set initial amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!    Adapted by Alexander C. Paul to use the restart logical, Oct 2020
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ccs), intent(inout) :: wf
      logical, intent(in)       :: restart
!
      if (.not. restart) then
!
         call zero_array(wf%t1, wf%n_t1)
!
      else
!
         if (wf%t_file%exists()) then
!
            call output%printf('m', 'Requested restart. Reading in solution from file.', &
                               fs='(/t3,a)')
!
            call wf%read_amplitudes()
            call wf%construct_t1_cholesky(wf%t1, wf%L_mo, wf%L_t1)
!
         else
!
            call zero_array(wf%t1, wf%n_t1)
!
         end if
!
      end if
!
   end subroutine set_initial_amplitudes_guess_ccs
!
!
   subroutine set_initial_multipliers_guess_ccs(wf, restart)
!!
!!    Set initial multipliers guess
!!    Written by Alexander C. Paul , Oct 2020
!!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      logical, intent(in) :: restart
!
      if (.not. restart) then
!
         call copy_and_scale(one, wf%t1, wf%t1bar, wf%n_t1)
!
      else
!
         if (wf%tbar_file%exists()) then
!
            call output%printf('m', 'Requested restart. Reading multipliers from file.', &
                              fs='(/t3,a)')
!
            call wf%read_multipliers()
!
         else
!
            call copy_and_scale(one, wf%t1, wf%t1bar, wf%n_t1)
!
         end if
!
      end if
!
   end subroutine set_initial_multipliers_guess_ccs
!
!
   subroutine construct_Jacobian_transform_ccs(wf, r_or_l, X, R, w)
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
!!
!!    R: On output contains the transformed vector
!!
!!    w: Excitation energy. Only used for debug prints for CCS, CCSD etc.
!!       but is passed to the effective_jacobian_transform for lowmem_CC2 and CC3
!!
      use warning_suppressor
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: R
!
      real(dp), intent(in), optional :: w
!
!     Suppress unused variable compiler warning for 'w'
      call do_nothing(w)
!
!     Compute the transformed matrix
      if (r_or_l .eq. "right") then
!
         call wf%jacobian_transformation(X, R) ! X <- AX
!
      elseif (r_or_l .eq. 'left') then
!
         call wf%jacobian_transpose_transformation(X, R) ! X <- XA
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
   subroutine prepare_for_Jacobians_ccs(wf, r_or_l)
!!
!!    Approximate Jacobian transform
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for preparations to left and right Jacobian transformations.
!!
!!    r_or_l: 'left', 'right', or 'both'
!!            (prepares for A^T, A, or both A^T and A)
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      if (r_or_l == 'right') then
!
         call wf%prepare_for_jacobian()
!
      else if (r_or_l == 'left') then
!
         call wf%prepare_for_jacobian_transpose()
!
      else if (r_or_l == 'both') then
!
         call wf%prepare_for_jacobian()
         call wf%prepare_for_jacobian_transpose()
!
      else
!
         call output%error_msg('Could not recognize "r_or_l" in prepare_for_Jacobians_ccs.')
!
      end if
!
   end subroutine prepare_for_Jacobians_ccs
!
!
   subroutine approximate_Jacobian_transform_ccs(wf, r_or_l, X, R, w)
!!
!!    Approximate Jacobian transform
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for a lower-level Jacobian transformation that is the best approximation
!!    with a lower computational scaling.
!!
      use array_utilities, only: zero_array
      use warning_suppressor
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: R
!
      real(dp), intent(in), optional :: w
!
      call output%error_msg('Approximate Jacobian transformation for this CC method &
                           & has not yet been implemented.')
!
!     Suppress unused variables
      call do_nothing(w)
      call do_nothing(r_or_l)
      call do_nothing(X)
      call zero_array(R, wf%n_t1)
!
   end subroutine approximate_Jacobian_transform_ccs
!
!
   subroutine prepare_for_approximate_Jacobians_ccs(wf, r_or_l)
!!
!!    Prepare for approximate Jacobians
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for preparations to a lower-level Jacobian transformation that is
!!    the best approximation with a lower computational scaling.
!!
!!    r_or_l: 'left', 'right', or 'both'
!!            (prepares for A^T, A, or both A^T and A)
!!
      use warning_suppressor
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
!     Suppress unused variables
      call do_nothing(wf)
      call do_nothing(r_or_l)
!
      call output%error_msg('Approximate Jacobian transformation for this CC method &
                           & has not yet been implemented. (Error due to call to prepare &
                           & for this transformation.)')
!
   end subroutine prepare_for_approximate_Jacobians_ccs
!
!
   subroutine get_cvs_projector_ccs(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      use array_utilities, only: zero_array
!
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
   subroutine get_rm_core_projector_ccs(wf, projector, n_cores, core_MOs)
!!
!!    Get remove core projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
!
      use array_utilities, only: constant_array
!
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
      call constant_array(projector, wf%n_es_amplitudes, one)
!
      do core = 1, n_cores
!
         i = core_MOs(core)
!
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
            projector(ai) = zero
!
         enddo
      enddo
!
   end subroutine get_rm_core_projector_ccs
!
!
   subroutine print_dominant_amplitudes_ccs(wf)
!!
!!    Print dominant amplitudes
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Prints the dominant amplitudes in the cluster amplitude vector.
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
!!    Prints the dominant amplitudes in the amplitude vector x.
!!
!!    tag specified the printed label for the vector, e.g. tag = "t" for
!!    the cluster amplitudes.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: x
!
      character(len=*) :: tag
!
      call wf%print_dominant_x1(x(1:wf%n_t1),tag)
!
   end subroutine print_dominant_x_amplitudes_ccs
!
!
   subroutine print_dominant_x1(wf, x1, tag)
!!
!!    Print dominant x1
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Prints the 10 most dominant single amplitudes,
!!    or sorts them if there are fewer than twenty of them.
!!
!!    tag specified the printed label for the vector, e.g. tag = "t" for
!!    the cluster amplitudes.
!!
      use array_utilities, only: get_n_highest
      use index_invert, only: invert_compound_index
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: x1
      character(len=*), intent(in)             :: tag
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
      call output%printf('m', 'Largest single amplitudes:', fs='(/t6,a)')
      call output%print_separator('m', 35, '-', fs='(t6,a)')
      call output%printf('m', 'a       i         (a0)(a,i)', fs='(t9,a)',chars=[tag])
      call output%print_separator('m', 35, '-', fs='(t6,a)')
!
      do elm = 1, n_elements
!
         call invert_compound_index(dominant_indices(elm), a, i, wf%n_v, wf%n_o)
!
         call output%printf('m', '(i4)    (i3)   (f19.12)', ints=[a, i], &
                            reals=[x1(dominant_indices(elm))], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('m', 36, '-', fs='(t6,a)')
!
      call mem%dealloc(dominant_indices, n_elements)
      call mem%dealloc(dominant_values, n_elements)
      call mem%dealloc(abs_x1, wf%n_t1)
!
   end subroutine print_dominant_x1
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
   function L_R_overlap_ccs(wf, L, left_state, R, right_state) result(overlap)
!!
!!    Calculates the overlap between a given left and right state
!!    Written by Josefine Andersen, Apr 2019
!!
!!    Modified by Alexander C. Paul to be overwritable for CC3 and low-mem CC2
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp) :: overlap
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L
      integer, intent(in) :: left_state
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R
      integer, intent(in) :: right_state
!
      real(dp) :: ddot
!
      overlap = ddot(wf%n_es_amplitudes, L, 1, R, 1)
!
      call output%printf('debug', 'Overlap of (i0). left and (i0). right state: (f15.10)', &
                         ints=[left_state, right_state], reals=[overlap], fs='(t6,a)')
!
   end function L_R_overlap_ccs
!
!
   subroutine scale_amplitudes_ccs(wf, t)
!!
!!    Scale amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Scales t to conform with the convention used in the wavefunction:
!!
!!       t1 <- t1
!!       t2_aiai <- two * t2_aiai
!!       ...
!!
      use warning_suppressor
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: t
!
      call do_nothing(t)
!
   end subroutine scale_amplitudes_ccs
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
      use array_utilities, only: zero_array, copy_and_scale
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_copy
      real(dp), dimension(:), allocatable :: orbital_energies_copy
!
      call mem%alloc(orbital_coefficients_copy, wf%ao%n, wf%n_mo)
!
      call copy_and_scale(one, wf%orbital_coefficients, orbital_coefficients_copy, wf%n_mo*wf%ao%n)
!
      call wf%destruct_orbital_coefficients()
!
      call mem%alloc(wf%orbital_coefficients, wf%ao%n, wf%n_mo + wf%n_bath_orbitals)
      call zero_array(wf%orbital_coefficients, wf%ao%n*(wf%n_mo + wf%n_bath_orbitals))
!
      wf%orbital_coefficients(1:wf%ao%n, 1:wf%n_mo) = orbital_coefficients_copy(:,:)
!
      call mem%dealloc(orbital_coefficients_copy, wf%ao%n, wf%n_mo)
!
      call mem%alloc(orbital_energies_copy, wf%n_mo)
!
      call copy_and_scale(one, wf%orbital_energies, orbital_energies_copy, wf%n_mo)
!
      call wf%destruct_orbital_energies()
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
      use array_utilities, only: zero_array
!
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
   subroutine approximate_double_excitation_vectors_ccs(wf, R_ai, omega, file_)
!!
!!    Construct approximate double excitation vectors
!!    Sarai D. Folkestad, May 2019
!!
!!    Approximate double excitation vectors:
!!
!!       R_aibj = 1/(1 + delta_aibj) * (P_ai,bj(sum_c R_ci g_bjac
!!                                              - sum_k R_bk g_kjai))/(-eps_aibj + omega^CCS)
!!              = 1/(1 + delta_aibj) * X_aibj/(-eps_aibj + omega^CCS)
!!
!!    Used for CNTOs from CCS.
!!
!!    For further information see Baudin, P. and Kristensen, K., J. Chem. Phys. 2017, 146, 214114
!!
!!    Note: Renormalizes R !
!!    Note: Integrals are in MO basis not T1 basis
!!
!!    Routine is organized in the following way:
!!
!!    1. Construct (term 1)
!!
!!       P_aibj (- sum_k R_bk g_kjai)/(-eps_aibj + omega^CCS)/(1 + delta_ai,bj)
!!
!!    2. Add (term 2)
!!
!!       P_ai,bj(sum_c R_ci g_bjac)/(-eps_aibj + omega^CCS)/(1+delta_ai,bj)
!!
!!       and calculate norm of doubles part
!!
!!    3. Normalize full excitation vector (R_ai, R_aibj)
!!
!!    Steps 1-3 are done in batching loops over a and b,
!!    but note that the same batching structure (sizes and numbers of batches)
!!    are used for all three loops using the maximal memory requirement.
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_12_to_21, sort_123_to_231, sort_1234_to_3412
      use reordering, only: symmetric_sum, packin
      use direct_stream_file_class, only: direct_stream_file
      use array_utilities, only: scale_diagonal, zero_array
      use sequential_file_class, only: sequential_file
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: R_ai
!
      real(dp), intent(in) :: omega
!
      type(batching_index) :: batch_a, batch_b
!
      type(direct_stream_file), intent(inout) :: file_ ! record length o^2v
!
      real(dp), dimension(:,:,:), allocatable :: L_Jai, L_Jbj, L_Jkj, L_kjJ, L_Jac, L_Jbc
      real(dp), dimension(:,:,:), allocatable :: X_Jai, X_Jbj, X_aiJ, X_bjJ
!
      real(dp), dimension(:,:,:,:), allocatable :: R_aibj, R_bjai
      real(dp), dimension(:,:,:,:), allocatable :: R_ibj_a, R_aibj_old, R_aibj_batch
!
      real(dp), dimension(:), allocatable :: R_aibj_packed
!
      integer :: current_a_batch, current_b_batch
      integer :: req0, req1_a, req1_b, req2
      integer :: a, i, b, j
!
      real(dp) :: ddot, R_norm, R_s_norm_sq, R_d_norm_sq
!
      type(sequential_file) :: file_temp_1, file_temp_2
!
!     Initialize temporary files
!
      file_temp_1 = sequential_file('approximate_doubles_temp_1')
      file_temp_2 = sequential_file('approximate_doubles_temp_2')
!
      call file_temp_1%open_('readwrite', 'rewind')
!
      call mem%alloc(L_Jkj, wf%eri_t1%n_J, wf%n_o, wf%n_o)
      call wf%L_mo%get(L_Jkj, 1, wf%n_o, 1, wf%n_o)
!
      call mem%alloc(L_kjJ, wf%n_o, wf%n_o, wf%eri_t1%n_J)
      call sort_123_to_231(L_Jkj, L_kjJ, wf%eri_t1%n_J, wf%n_o, wf%n_o)
      call mem%dealloc(L_Jkj, wf%eri_t1%n_J, wf%n_o, wf%n_o)
!
!     Prepare for batching
!
      req0 = 0
!
      req1_a = max(2*wf%n_o*wf%eri_t1%n_J, & ! X_ai_J + L_ai_J  (term 1)
                  wf%n_v*wf%eri_t1%n_J + wf%n_o*wf%eri_t1%n_J, & ! X_ai_J + L_ac_J term 2)
                  wf%n_o**2*wf%n_v) ! R_aibj no batching on b
!
      req1_b = req1_a ! Not exactly true but we want equal batches of a and b
!
      req2 = 2*(wf%n_o**2) ! R_aibj
!
      batch_a = batching_index(wf%n_v)
      batch_b = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, batch_b, req0, req1_a, req1_b, req2, &
                           tag='approximate_double_excitation_vectors_ccs')
!
!     :: Construct term 1 and store in file_temp_1
!
!        P_aibj (- sum_k R_bk g_kjai)/(-eps_aibj + omega^CCS)/(1 + delta_ai,bj)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(L_Jai, wf%eri_t1%n_J, batch_a%length, wf%n_o)
         call wf%L_mo%get(L_Jai, wf%n_o + batch_a%first, &
                                             wf%n_o + batch_a%get_last(), 1, wf%n_o)
!
         call mem%alloc(X_aiJ, batch_a%length, wf%n_o, wf%eri_t1%n_J)
!
         call dgemm('N', 'N',               &
                     batch_a%length,        &
                     wf%n_o*wf%eri_t1%n_J,  &
                     wf%n_o,                &
                     one,                   &
                     R_ai(batch_a%first,1), &
                     wf%n_v,                &
                     L_kjJ,                 &
                     wf%n_o,                &
                     zero,                  &
                     X_aiJ,                 &
                     batch_a%length)
!
         do current_b_batch = 1, batch_b%num_batches
!
            if (current_b_batch == current_a_batch) cycle
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(X_bjJ, batch_b%length, wf%n_o, wf%eri_t1%n_J)
!
            call dgemm('N', 'N',               &
                        batch_b%length,        &
                        wf%n_o*wf%eri_t1%n_J,  &
                        wf%n_o,                &
                        one,                   &
                        R_ai(batch_b%first,1), &
                        wf%n_v,                &
                        L_kjJ,                 &
                        wf%n_o,                &
                        zero,                  &
                        X_bjJ,                 &
                        batch_b%length)
!
            call mem%alloc(R_bjai, batch_b%length, wf%n_o, batch_a%length, wf%n_o)
!
            call dgemm('N', 'N',               &
                        batch_b%length*wf%n_o, &
                        batch_a%length*wf%n_o, &
                        wf%eri_t1%n_J,         &
                        -one,                  &
                        X_bjJ,                 &
                        batch_b%length*wf%n_o, &
                        L_Jai,                 &
                        wf%eri_t1%n_J,         &
                        zero,                  &
                        R_bjai,                &
                        batch_b%length*wf%n_o)
!
            call mem%dealloc(X_bjJ, batch_b%length, wf%n_o, wf%eri_t1%n_J)
!
            call mem%alloc(R_aibj, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
            call sort_1234_to_3412(R_bjai, R_aibj, batch_b%length, wf%n_o, batch_a%length, wf%n_o)
!
            call mem%dealloc(R_bjai, batch_b%length, wf%n_o, batch_a%length, wf%n_o)
!
            call mem%alloc(L_Jbj, wf%eri_t1%n_J, batch_b%length, wf%n_o)
            call wf%L_mo%get(L_Jbj, wf%n_o + batch_b%first, &
                                                     wf%n_o + batch_b%get_last(),  &
                                                     1, wf%n_o)
!
            call dgemm('N', 'N', &
                        batch_a%length*wf%n_o,  &
                        batch_b%length*wf%n_o,  &
                        wf%eri_t1%n_J,          &
                        -one,                   &
                        X_aiJ,                  &
                        batch_a%length*wf%n_o,  &
                        L_Jbj,                  &
                        wf%eri_t1%n_J,          &
                        one,                    &
                        R_aibj,                 &
                        batch_a%length*wf%n_o)
!
            call mem%dealloc(L_Jbj, wf%eri_t1%n_J, batch_b%length, wf%n_o)
!
!           Divide by orbital differences and CCS excitation energy
!
!$omp parallel do private(a, i, b, j)
            do j = 1, wf%n_o
               do b = 1, batch_b%length
                  do i = 1, wf%n_o
                     do a = 1, batch_a%length
!
                        R_aibj(a, i, b, j) = R_aibj(a, i, b, j)/&
                                    (- wf%orbital_energies(batch_a%first - 1 + a + wf%n_o) &
                                     - wf%orbital_energies(batch_b%first - 1 + b + wf%n_o) &
                                     + wf%orbital_energies(i) + wf%orbital_energies(j) + omega)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
!           Write contribution to file
!
            call file_temp_1%write_(R_aibj, batch_b%length*(wf%n_o**2)*batch_a%length)
!
            call mem%dealloc(R_aibj, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
         enddo
!
!        Calculate a = b block
!
         call mem%alloc(R_aibj, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
         call dgemm('N', 'N',                &
                     batch_a%length*wf%n_o,  &
                     batch_a%length*wf%n_o,  &
                     wf%eri_t1%n_J,          &
                     -one,                   &
                     X_aiJ,                  &
                     batch_a%length*wf%n_o,  &
                     L_Jai,                  &
                     wf%eri_t1%n_J,          &
                     zero,                   &
                     R_aibj,                 &
                     batch_a%length*wf%n_o)
!
         call mem%dealloc(L_Jai, wf%eri_t1%n_J, batch_a%length, wf%n_o)
         call mem%dealloc(X_aiJ, batch_a%length, wf%n_o, wf%eri_t1%n_J)
!
!        Symmetrize
!
         call symmetric_sum(R_aibj, wf%n_o*batch_a%length)
!
!        Binormalization factor
!
         call scale_diagonal(half, R_aibj, batch_a%length*wf%n_o)
!
!        Divide by orbital differences and CCS excitation energy
!
!$omp parallel do private(a, i, b, j)
         do j = 1, wf%n_o
            do b = 1, batch_a%length
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     R_aibj(a, i, b, j) = R_aibj(a, i, b, j)/&
                                 (- wf%orbital_energies(batch_a%first - 1 + a + wf%n_o) &
                                  - wf%orbital_energies(batch_a%first - 1 + b + wf%n_o) &
                                  + wf%orbital_energies(i) + wf%orbital_energies(j) + omega)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
!        Write contribution to file
!
         call file_temp_1%write_(R_aibj, (batch_a%length**2)*(wf%n_o**2))
!
         call mem%dealloc(R_aibj, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(L_kjJ, wf%n_o, wf%n_o, wf%eri_t1%n_J)
!
!     Prepare files
!
      call file_temp_1%rewind_()
      call file_temp_2%open_('readwrite', 'rewind')
!
!     ::  Add term 2, store R_aibj in file_temp_2, and calculate norm of doubles part
!
!        += P_ai,bj(sum_c R_ci g_bjac)/(-eps_aibj + omega^CCS)/(1+delta_ai,bj)
!
!     Accumulate dot product of double vectors for subsequent normalization
!
      R_d_norm_sq = 0
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(L_Jac, wf%eri_t1%n_J, batch_a%length, wf%n_v)
         call wf%L_mo%get(L_Jac,                                         &
                          batch_a%first + wf%n_o, batch_a%get_last() + wf%n_o,  &
                          1 + wf%n_o, wf%n_mo)
!
         call mem%alloc(X_Jai, wf%eri_t1%n_J, batch_a%length, wf%n_o)
!
!        X_ai_J = sum_c R_ci L_ac_J
!
         call dgemm('N', 'N',                         &
                     batch_a%length*(wf%eri_t1%n_J),  &
                     wf%n_o,                          &
                     wf%n_v,                          &
                     one,                             &
                     L_Jac,                           & ! L_Ja_c
                     batch_a%length*(wf%eri_t1%n_J),  &
                     R_ai,                            & ! R_c_i
                     wf%n_v,                          &
                     zero,                            &
                     X_Jai,                           &
                     batch_a%length*(wf%eri_t1%n_J))
!
         call mem%dealloc(L_Jac, wf%eri_t1%n_J, batch_a%length, wf%n_v)
!
         call mem%alloc(L_Jai, wf%eri_t1%n_J, batch_a%length, wf%n_o)
         call wf%L_mo%get(L_Jai, batch_a%first + wf%n_o, &
                       batch_a%get_last() + wf%n_o, 1, wf%n_o)
!
         do current_b_batch = 1, batch_b%num_batches
!
            if (current_b_batch == current_a_batch) cycle
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(R_aibj, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
            call mem%alloc(L_Jbj, wf%eri_t1%n_J, batch_b%length, wf%n_o)
            call wf%L_mo%get(L_Jbj, batch_b%first + wf%n_o, &
                          batch_b%get_last() + wf%n_o, 1, wf%n_o)
!
            call dgemm('T', 'N',                &
                        batch_a%length*wf%n_o,  &
                        batch_b%length*wf%n_o,  &
                        wf%eri_t1%n_J,          &
                        one,                    &
                        X_Jai,                  & ! X_J_ai
                        wf%eri_t1%n_J,          &
                        L_Jbj,                  &
                        wf%eri_t1%n_J,          &
                        zero,                   &
                        R_aibj,                 &
                        batch_a%length*wf%n_o)
!
            call mem%dealloc(L_Jbj, wf%eri_t1%n_J, batch_b%length, wf%n_o)
!
            call mem%alloc(L_Jbc, wf%eri_t1%n_J, batch_b%length, wf%n_v)
            call wf%L_mo%get(L_Jbc,                                         &
                          batch_b%first + wf%n_o, batch_b%get_last() + wf%n_o, &
                          1 + wf%n_o, wf%n_mo)
!
            call mem%alloc(X_Jbj, wf%eri_t1%n_J, batch_b%length, wf%n_o)
!
!           X_bj_J = sum_c R_cj L_bc_J
!
            call dgemm('N', 'N',                         &
                        batch_b%length*(wf%eri_t1%n_J),  &
                        wf%n_o,                          &
                        wf%n_v,                          &
                        one,                             &
                        L_Jbc,                           & ! L_Jb_c
                        batch_b%length*(wf%eri_t1%n_J),  &
                        R_ai,                            & ! R_c_j
                        wf%n_v,                          &
                        zero,                            &
                        X_Jbj,                           &
                        batch_b%length*(wf%eri_t1%n_J))
!
            call mem%dealloc(L_Jbc, wf%eri_t1%n_J, batch_b%length, wf%n_v)
!
            call dgemm('T', 'N',                &
                        batch_a%length*wf%n_o,  &
                        batch_b%length*wf%n_o,  &
                        wf%eri_t1%n_J,          &
                        one,                    &
                        L_Jai,                  & ! L_J_ai
                        wf%eri_t1%n_J,          &
                        X_Jbj,                  & ! X_J_bj
                        wf%eri_t1%n_J,          &
                        one,                    &
                        R_aibj,                 &
                        batch_a%length*wf%n_o)
!
!
            call mem%dealloc(X_Jbj, wf%eri_t1%n_J, batch_b%length, wf%n_o)
!
!           Divide by orbital differences and CCS excitation energy
!
!$omp parallel do private(a, i, b, j)
            do j = 1, wf%n_o
               do b = 1, batch_b%length
                  do i = 1, wf%n_o
                     do a = 1, batch_a%length
!
                        R_aibj(a, i, b, j) = R_aibj(a, i, b, j)/&
                                    (- wf%orbital_energies(batch_a%first - 1 + a + wf%n_o) &
                                     - wf%orbital_energies(batch_b%first - 1 + b + wf%n_o) &
                                     + wf%orbital_energies(i) + wf%orbital_energies(j) + omega)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%alloc(R_aibj_old, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
            call file_temp_1%read_(R_aibj_old, (batch_a%length)*(batch_b%length)*(wf%n_o**2))
!
            call daxpy((batch_a%length)*(batch_b%length)*(wf%n_o**2), one, R_aibj_old, 1, R_aibj, 1)
            call mem%dealloc(R_aibj_old, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
!           Add norm
!
            R_d_norm_sq = R_d_norm_sq + half*ddot(batch_a%length*(wf%n_o**2)*batch_b%length, &
                                             R_aibj, 1, R_aibj, 1)
!
!           Write file
!
            call file_temp_2%write_(R_aibj, batch_a%length*(wf%n_o**2)*batch_b%length)
!
            call mem%dealloc(R_aibj, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
         enddo
!
!        Calculate a = b block
!
         call mem%alloc(R_aibj, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
         call dgemm('T', 'N',                &
                     batch_a%length*wf%n_o,  &
                     batch_a%length*wf%n_o,  &
                     wf%eri_t1%n_J,          &
                     one,                    &
                     L_Jai,                  &
                     wf%eri_t1%n_J,          &
                     X_Jai,                  & ! X_J_bj
                     wf%eri_t1%n_J,          &
                     zero,                   &
                     R_aibj,                 &
                     batch_a%length*wf%n_o)
!
         call mem%dealloc(X_Jai, wf%eri_t1%n_J, batch_a%length, wf%n_o)
         call mem%dealloc(L_Jai, wf%eri_t1%n_J, batch_a%length, wf%n_o)
!
!        Symmetrize
!
         call symmetric_sum(R_aibj, wf%n_o*batch_a%length)
!
!        Binormalization factor
!
         call scale_diagonal(half, R_aibj, batch_a%length*wf%n_o)
!
!        Divide by orbital differences and CCS excitation energy
!
!$omp parallel do private(a, i, b, j)
         do j = 1, wf%n_o
            do b = 1, batch_a%length
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     R_aibj(a, i, b, j) = R_aibj(a, i, b, j)/&
                                 (- wf%orbital_energies(batch_a%first - 1 + a + wf%n_o) &
                                  - wf%orbital_energies(batch_a%first - 1 + b + wf%n_o) &
                                  + wf%orbital_energies(i) + wf%orbital_energies(j) + omega)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%alloc(R_aibj_old, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
         call file_temp_1%read_(R_aibj_old, (batch_a%length**2)*(wf%n_o**2))
!
         call daxpy((batch_a%length**2)*(wf%n_o**2), one, R_aibj_old, 1, R_aibj, 1)
         call mem%dealloc(R_aibj_old, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
!        Pack in
!
         call mem%alloc(R_aibj_packed, batch_a%length*wf%n_o*(batch_a%length*wf%n_o + 1)/2)
!
         call packin(R_aibj_packed, R_aibj, batch_a%length*wf%n_o)
!
!        Add norm
!
         R_d_norm_sq = R_d_norm_sq + ddot(batch_a%length*wf%n_o*(batch_a%length*wf%n_o + 1)/2, &
                                          R_aibj_packed, 1, R_aibj_packed, 1)
!
         call mem%dealloc(R_aibj_packed, batch_a%length*wf%n_o*(batch_a%length*wf%n_o + 1)/2)
!
!        Write file
!
         call file_temp_2%write_(R_aibj, (batch_a%length**2)*(wf%n_o**2))
!
         call mem%dealloc(R_aibj, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
      enddo
!
!     Prepare files
!
      call file_temp_1%close_('delete')
      call file_temp_2%rewind_()
      call file_%open_('write')
!
!     :: Normalization
!
!        Determine norm and rescale
!
      R_s_norm_sq = ddot(wf%n_o*wf%n_v, R_ai, 1, R_ai, 1)
      R_norm = sqrt(R_s_norm_sq + R_d_norm_sq)
!
      call dscal(wf%n_o*wf%n_v, one/R_norm, R_ai, 1)
!
         call mem%alloc(R_ibj_a, wf%n_o, wf%n_v, wf%n_o, batch_a%max_length)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         do current_b_batch = 1, batch_b%num_batches
!
            if (current_b_batch == current_a_batch) cycle
!
            call batch_b%determine_limits(current_b_batch)
!
            call mem%alloc(R_aibj_batch, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
            call file_temp_2%read_(R_aibj_batch, batch_a%length*(wf%n_o**2)*batch_b%length)
            call dscal(batch_a%length*(wf%n_o**2)*batch_b%length, one/R_norm, R_aibj_batch, 1)
!
!$omp parallel do private(a, i, b, j)
            do a = 1, batch_a%length
               do j = 1, wf%n_o
                  do b = 1, batch_b%length
                     do i = 1, wf%n_o
!
                        R_ibj_a(i, batch_b%first - 1 + b, j, a) = R_aibj_batch(a, i, b, j)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(R_aibj_batch, batch_a%length, wf%n_o, batch_b%length, wf%n_o)
!
         enddo
!
         call mem%alloc(R_aibj_batch, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
         call file_temp_2%read_(R_aibj_batch, (batch_a%length**2)*(wf%n_o**2))
!
         call dscal((batch_a%length**2)*(wf%n_o**2), one/R_norm, R_aibj_batch, 1)
!
!$omp parallel do private(a, i, b, j)
            do a = 1, batch_a%length
               do j = 1, wf%n_o
                  do b = 1, batch_a%length
                     do i = 1, wf%n_o
!
                        R_ibj_a(i, batch_a%first - 1 + b, j, a) = R_aibj_batch(a, i, b, j)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
         call mem%dealloc(R_aibj_batch, batch_a%length, wf%n_o, batch_a%length, wf%n_o)
!
         call file_%write_(R_ibj_a, batch_a%first, batch_a%get_last())
!
      enddo

      call mem%dealloc(R_ibj_a, wf%n_o, wf%n_v, wf%n_o, batch_a%max_length)
!
      call file_temp_2%close_('delete')
      call file_%close_('keep')
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
      if (input%is_keyword_present('core excitation', 'solver cc es')) then
!
!        Determine the number of core MOs
!
         wf%n_core_MOs = input%get_n_elements_for_keyword('core excitation', 'solver cc es')
!
!        Then read the vector of core MOs
!
         call wf%initialize_core_MOs()
!
         call input%get_array_for_keyword('core excitation', 'solver cc es', wf%n_core_MOs, wf%core_MOs)
!
      else
!
         call output%error_msg('found no specified core MOs in input for CVS calculation')
!
      endif
!
      wf%cvs = .true.
!
!     Consistency check (given currently implemented CVS capabilities)
!
      if (input%is_keyword_present('core excitation', 'solver cc es') .and. &
          input%is_keyword_present('core', 'frozen orbitals')) then
!
         call output%error_msg('No support for CVS with frozen core yet. Turn off frozen core.')
!
      endif
!
   end subroutine read_cvs_settings_ccs
!
!
   subroutine read_rm_core_settings_ccs(wf)
!!
!!    Read remove core settings
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(ccs) :: wf
!
      if (input%is_keyword_present('remove core', 'solver cc es')) then
!
!        Determine the number of core MOs
!
         wf%n_core_MOs = input%get_n_elements_for_keyword('remove core', 'solver cc es')
!
!        Then read the vector of core MOs
!
         call wf%initialize_core_MOs()
!
         call input%get_array_for_keyword('remove core', 'solver cc es', wf%n_core_MOs, wf%core_MOs)
!
         wf%rm_core = .true.
!
      else
!
         call output%error_msg('found no specified core MOs to remove!')
!
      endif
!
!     Consistency check
!
      if (input%is_keyword_present('remove core', 'solver cc es') .and. &
          input%is_keyword_present('core', 'frozen orbitals')) then
!
         call output%error_msg('No support for removal of core with frozen core yet. Turn off frozen core.')
!
      endif
!
   end subroutine read_rm_core_settings_ccs
!
!
   subroutine integral_preparations_ccs(wf, n_J)
!!
!!    Integral preparations
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_J
!
      call output%printf('m', '- Settings for integral handling:', fs='(/t3,a)')
!
      call wf%cholesky_factory(n_J)
      call wf%integral_factory()
!
   end subroutine integral_preparations_ccs
!
!
   subroutine cholesky_factory(wf, n_J)
!!
!!    Cholesky factory
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      use eri_cholesky_disk_class,     only: eri_cholesky_disk
      use eri_cholesky_memory_class,   only: eri_cholesky_memory
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_J
!
      character(len=200) :: cholesky_storage
!
      logical :: has_room_for_cholesky
      logical :: cholesky_in_memory
!
      has_room_for_cholesky = 2*n_J*wf%n_mo**2*dp < mem%get_available()/5
!
      cholesky_storage = 'memory'
      call input%get_keyword('cholesky storage', 'integrals', cholesky_storage)
!
      if (trim(cholesky_storage) /= 'memory' .and. trim(cholesky_storage) /= 'disk') &
         call output%error_msg('illegal value for cholesky storage: (a0)', &
                               chars=[trim(cholesky_storage)])
!
      cholesky_in_memory = trim(cholesky_storage) == 'memory' .and. has_room_for_cholesky
!
      if (cholesky_in_memory) then
         wf%L_t1 = eri_cholesky_memory()
         wf%L_mo = eri_cholesky_memory()
      else
         wf%L_t1 = eri_cholesky_disk('T1')
         wf%L_mo = eri_cholesky_disk('MO')
      endif
!
      call wf%L_mo%initialize(n_J, 2, [wf%n_o, wf%n_v])
      call wf%L_t1%initialize(n_J, 2, [wf%n_o, wf%n_v])
!
      call output%printf('m', &
                         'Cholesky vectors in memory: (l0)', &
                         logs=[cholesky_in_memory],          &
                         fs='(/t6,a)')
!
   end subroutine cholesky_factory
!
!
   subroutine integral_factory(wf)
!!
!!    Integral factory
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!
      use eri_memory_tool_class,       only: eri_memory_tool
      use eri_tool_class,              only: eri_tool
      use eri_adapter_class,           only: eri_adapter
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=200) :: integral_storage
!
      logical :: has_room_for_eri
!
      logical :: try_to_place_eri_in_memory, eri_in_memory
!
      integral_storage = 'memory'
      call input%get_keyword('eri storage', 'integrals', integral_storage)
!
      if (trim(integral_storage) /= 'memory' .and. trim(integral_storage) /= 'none') &
         call output%error_msg('illegal value for eri storage: (a0)', &
                               chars=[trim(integral_storage)])
!
      try_to_place_eri_in_memory = (integral_storage == 'memory' .and. wf%need_g_abcd)
!
      has_room_for_eri = (wf%n_mo**4)*dp < mem%get_available()/5
!
      eri_in_memory = try_to_place_eri_in_memory .and. has_room_for_eri
!
      if (eri_in_memory) then
         wf%eri_t1 = eri_adapter(eri_memory_tool(wf%L_t1), wf%n_o, wf%n_v)
      else
         wf%eri_t1 = eri_adapter(eri_tool(wf%L_t1), wf%n_o, wf%n_v)
      endif
!
      call wf%L_t1%add_observer('eri t1', wf%eri_t1%eri)
!
      call output%printf('m', &
                         'ERI matrix in memory:       (l0)', &
                         logs=[eri_in_memory],               &
                         fs='(t6,a)')
!
   end subroutine integral_factory
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
      call output%printf('v', 'No MO preparations for (a0)', chars=[wf%name_])
!
   end subroutine mo_preparations_ccs
!
!
   subroutine get_degenerate_states_ccs(wf, state, transformation, threshold, &
                                        n_degeneracy, degenerate)
!!
!!    Get degenerate states
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    Checks if there are degeneracies between the given excited "state"
!!    and returns the number of degenerate states
!!
!!    state:          number of the given state
!!    transformation: left, right or both
!!    threshold:      threshold to compare the energies to
!!    n_degeneracy:   number of states with the same energy as "state" (at least 1)
!!    degenerate:     optional array of logicals
!!                    specifying which state is degenerate to "state"
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in) :: state
!
      character(len=*), intent(in) :: transformation
!
      real(dp), intent(in) :: threshold
!
      integer, intent(out) :: n_degeneracy
!
      logical, dimension(wf%n_singlet_states), intent(out) :: degenerate
!
      real(dp) :: delta_l, delta_r
      integer  :: p, q
      logical  :: state_found
!
      n_degeneracy = 1
      state_found = .true.
      degenerate = .false.
      degenerate(state) = .true.
!
      do while (state_found)
!
         state_found = .false.
!
!        Loop over degenerate states
         do p = 1, wf%n_singlet_states
            if (.not. degenerate(p)) cycle
!
!           Loop over all non-degenerate states
            do q = 1, wf%n_singlet_states
               if (degenerate(q)) cycle
!
               if (trim(transformation) .eq. 'right') then
!
                  delta_r = wf%right_excitation_energies(p) - wf%right_excitation_energies(q)
!
                  if (abs(delta_r) .gt. threshold) cycle
!
               else if (trim(transformation) .eq. 'left') then
!
                  delta_l = wf%left_excitation_energies(p) - wf%left_excitation_energies(q)
!
                  if (abs(delta_l) .gt. threshold) cycle
!
               else if (trim(transformation) .eq. 'both') then
!
                  delta_l = wf%left_excitation_energies(p) - wf%left_excitation_energies(q)
                  delta_r = wf%right_excitation_energies(p) - wf%right_excitation_energies(q)
!
                  if (abs(delta_l) .gt. threshold .and. abs(delta_r) .gt. threshold) then
!
                     cycle
!
                  else if  (abs(delta_l) .gt. threshold .neqv. abs(delta_r) .gt. threshold) then
!
                     call output%error_msg('Different degree of degeneracy in the &
                        &left excited states compared to the right excited states')
!
                  end if
!
               else
!
                  call output%error_msg('Tried to check for degenerate states but argument &
                                       &(a0) not recognized.', chars=[trim(transformation)])
!
               end if
!
               n_degeneracy = n_degeneracy + 1
               degenerate(q) = .true.
               state_found = .true.
!
            end do
         end do
      end do
!
   end subroutine get_degenerate_states_ccs
!
!
   subroutine check_for_parallel_states_ccs(wf, side, threshold, parallel_states)
!!
!!    Check for parallel states
!!    Written by Alexander C. Paul and Rolf H. Myhre, Oct 2019
!!
!!    The solver might not ensure orthogonality of the states (e.g. in DIIS).
!!    Therefore, degenerate roots might be found where the eigenvectors
!!    are equal up to a sign (and within the convergence threshold)
!!
      use array_utilities, only: check_for_parallel_vectors
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: side
!
      real(dp), intent(in) :: threshold
!
      logical, dimension(wf%n_singlet_states), intent(out), optional :: parallel_states
!
      real(dp), dimension(:,:), allocatable  :: vector
      logical, dimension(:), allocatable :: parallel, degenerate, checked
!
      integer, dimension(:), allocatable :: degenerate_states
      integer :: current_state, p, counter
      integer :: n_degeneracy, reduced_degeneracy
!
      character(len=100) :: print_states
!
      parallel_states = .false.
!
      call mem%alloc(degenerate, wf%n_singlet_states)
!
      call mem%alloc(checked, wf%n_singlet_states)
      checked = .false.
!
      do current_state = 1, wf%n_singlet_states
!
!        keep track of which states have been checked,
!        cycle if they were checked already
!
         if (checked(current_state)) cycle
!
!        Reset array indicating degeneracies between current_state and the other states
         degenerate = .false.
!
!        Get the number of states degenerate with current_state (at least 1)
!
         call wf%get_degenerate_states(current_state, side,threshold, &
                                          n_degeneracy, degenerate)
!
!        If degeneracies are found:
!           - read all degenerate states
!           - loop through all state pairs (p != q):
!                 -If 2 states are parallel print warning
!                  and reduce the degree of degeneracy
!
         if (n_degeneracy .gt. 1) then
!
            call mem%alloc(vector, wf%n_es_amplitudes, n_degeneracy)
!
            call mem%alloc(degenerate_states, n_degeneracy)
!
            counter = 0
!
            do p = current_state, wf%n_singlet_states
!
!              Only degenerate states have to be checked
               if (.not. degenerate(p)) cycle
!
               counter = counter + 1
!
!              save state number of the degenerate states
               degenerate_states(counter) = p
!
!              This state won't be checked again
               checked(p) = .true.
!
               call wf%read_excited_state(vector(:, counter), p, p, side)
!
            end do
!
            call mem%alloc(parallel, n_degeneracy)
!
            call check_for_parallel_vectors(vector, wf%n_es_amplitudes, n_degeneracy, &
                                            reduced_degeneracy, threshold, parallel)
!
            if(reduced_degeneracy .gt. 1) then
!
               write(print_states,'(a,i0)') ': ', current_state
!
!              Write string of degenerate but not parallel states.
!              state (p=1) has already been added to the string
               do p = 2, n_degeneracy
!
                  if (parallel(p)) cycle
!
                  write(print_states,'(a,a,i0)') trim(print_states), &
                                                 ', ', degenerate_states(p)
!
               end do
!
               call output%printf('n', 'Found a (i0)-fold degeneracy of the states(a0).', &
                                   ints=[reduced_degeneracy], chars=[print_states], &
                                   fs='(/t6,a)')
!
            end if
!
            do p = 1, n_degeneracy
               if(parallel(p)) then
!
                  call output%warning_msg('(a0) state (i0) is parallel to (a0) state (i0).', &
                                           ints=[current_state,degenerate_states(p)],&
                                           chars=[trim(side),trim(side)])
!
                  if (present(parallel_states)) parallel_states(degenerate_states(p)) = .true.
!
               end if
            end do
!
            call mem%dealloc(degenerate_states, n_degeneracy)
!
            call mem%dealloc(parallel, n_degeneracy)
            call mem%dealloc(vector, wf%n_es_amplitudes, n_degeneracy)
!
         end if
!
      end do
!
      call mem%dealloc(degenerate, wf%n_singlet_states)
      call mem%dealloc(checked, wf%n_singlet_states)
!
   end subroutine check_for_parallel_states_ccs
!
!
   subroutine biorthonormalize_L_and_R_ccs(wf, energy_threshold, residual_threshold)
!!
!!    Biorthonormalize the left and right eigenvectors
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    For that:
!!       - check if left and right excitation energies are consistent
!!       - check for degenerate states
!!       - check for and discard parallel states
!!       - If degeneracies are present: biorthonormalize using Gram-Schmidt
!!         following Kohaupt, L., Rocky Mountain J. Math., 44, 1265, (2014)
!!       - Normalize and write to file
!!
      use array_utilities, only: gram_schmidt_biorthonormalization
      use array_utilities, only: check_for_parallel_vectors
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), intent(in) :: energy_threshold, residual_threshold
!
      real(dp), dimension(:,:), allocatable :: R, L
!
      logical, dimension(:), allocatable :: degenerate, biorthonormalized
!
      real(dp) :: LT_R
!
      integer, dimension(:), allocatable :: degenerate_states
!
      integer :: state, current_state, counter
      integer :: n_degeneracy
!
      type(timings) :: timer
!
      timer = timings('Biorthonormalization of left and right states', pl='n')
      call timer%turn_on()
!
      call output%printf('v', 'Biorthonormalization of left and right excited &
                         &state vectors', fs='(/t3,a)')
!
!     Prepare logicals to keep track of which states have already been
!     biorthonormalized and which states are degenerate to the current_state
!
      call mem%alloc(biorthonormalized, wf%n_singlet_states)
      biorthonormalized = .false.
!
      call mem%alloc(degenerate, wf%n_singlet_states)
!
!     Loop through states look for degeneracies and
!     biorthonormalize degenerate states
!
      do current_state = 1, wf%n_singlet_states
!
!        Keep track of which states have already been biorthonormalized
         if (biorthonormalized(current_state)) cycle
!
!        Reset array indicating degenerate states
         degenerate = .false.
!
         if (abs(wf%left_excitation_energies(current_state) &
               - wf%right_excitation_energies(current_state)) .gt. 2*energy_threshold) then
!
!           roots might be ordered incorrectly
!           if so biorthonormalization will give an error
!
            call output%warning_msg('Eigenvector (i0) is not left-right consistent &
                                    &to threshold (e9.3).', ints=[current_state],  &
                                    reals=[2*energy_threshold], fs='(/t6,a)')
!
            call output%printf('m', 'Energies (left, right): (f19.12) (f19.12)',  &
                               reals=[wf%left_excitation_energies(current_state), &
                               wf%right_excitation_energies(current_state)], fs='(t12,a)')
!
            call output%printf('m', 'Biorthonormalization might lead to an error',  &
                                fs='(t12,a)')
!
         end if
!
!        :: Check for degeneracies in both left and right excitation energies ::
!
         call wf%get_degenerate_states(current_state, 'both', 2*energy_threshold, &
                                       n_degeneracy, degenerate)
!
         call output%printf('v', 'Degree of degeneracy: (i0)', &
                            ints=[n_degeneracy], fs='(/t6,a)')
!
!        Read in degenerate states L/R
!        States 1 to current_state have already been biorthonormalized
!
         call mem%alloc(R, wf%n_es_amplitudes, n_degeneracy)
         call mem%alloc(L, wf%n_es_amplitudes, n_degeneracy)
!
         call mem%alloc(degenerate_states, n_degeneracy)
!
         counter = 0
!
         do state = 1, wf%n_singlet_states
!
            if (.not. degenerate(state)) cycle
!
            counter = counter + 1
!
!           save state number of the degenerate states
            degenerate_states(counter) = state
!
!           This state won't be biorthonormalized again
            biorthonormalized(state) = .true.
!
            call wf%read_excited_state(R(:, counter), state, state,'right')
            call wf%read_excited_state(L(:, counter), state, state,'left')
!
         end do
!
         if (n_degeneracy .gt. 1) then
!
            call output%printf('n', 'Found states that are close in energy:', fs='(/t6,a)')
            call output%print_separator('n', 38,'-', fs='(t6,a)')
            call output%printf('n', 'State     Excitation Energy', fs='(t6,a)')
!
            do state = 1, n_degeneracy
!
               call output%printf('n', ' (i2)     (f19.12)', &
                                  ints=[degenerate_states(state)], &
                                  reals=[wf%right_excitation_energies(current_state+state-1)], &
                                  fs='(t6,a)')
!
            end do
!
         end if
!
!        :: Biorthonormalize states ::
!        -----------------------------
!
         call gram_schmidt_biorthonormalization(L, R, wf%n_es_amplitudes, &
                                                n_degeneracy, residual_threshold)
!
!        We binormalize again to include the doubles/triples contributions
!        in low-memory CC2 and CC3.
!
         do state = 1, n_degeneracy
!
            LT_R = wf%L_R_overlap(L(:,state),               &
                                  degenerate_states(state), &
                                  R(:,state),               &
                                  degenerate_states(state))
!
            call dscal(wf%n_es_amplitudes, one/LT_R, L(:, state), 1)
!
         end do
!
         counter = 0
!
         do state = 1, wf%n_singlet_states
!
!           Only save the states that we currently consider
            if (.not. degenerate(state)) cycle
!
            counter = counter + 1
!
            call wf%save_excited_state(R(:,counter), state, state,'right', &
                                       wf%right_excitation_energies(state))
            call wf%save_excited_state(L(:,counter), state, state,'left', &
                                       wf%left_excitation_energies(state))
!
         end do
!
         call mem%dealloc(degenerate_states, n_degeneracy)
         call mem%dealloc(R, wf%n_es_amplitudes, n_degeneracy)
         call mem%dealloc(L, wf%n_es_amplitudes, n_degeneracy)
!
      end do
!
      call mem%dealloc(biorthonormalized, wf%n_singlet_states)
!
      call mem%dealloc(degenerate, wf%n_singlet_states)
!
      call timer%turn_off()
!
   end subroutine biorthonormalize_L_and_R_ccs
!
!
   subroutine print_gs_summary_ccs(wf)
!!
!!    Print ground state summary
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call output%printf('m', '- Ground state summary:', fs='(/t3,a)')
!
      call output%printf('m', 'Final ground state energy (a.u.): (f18.12)', &
                         reals=[wf%energy], fs='(/t6,a)')
!
      call output%printf('m', 'Correlation energy (a.u.):        (f18.12)', &
                         reals=[wf%correlation_energy], fs='(/t6,a)')
!
      call wf%print_dominant_amplitudes()
!
   end subroutine print_gs_summary_ccs
!
!
   subroutine print_X1_diagnostics_ccs(wf, X, label)
!!
!!    Get X1 diagnostics
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      character(len=1), intent(in) :: label
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: X
!
      real(dp) :: get_X1_diagnostics
!
      get_X1_diagnostics = get_l2_norm(X(1:wf%n_t1),wf%n_t1)&
                           /get_l2_norm(X(1:wf%n_es_amplitudes),wf%n_es_amplitudes)
!
      call output%printf('n', 'Fraction singles (|(a0)1|/|(a0)|):  (f19.12)', &
            reals=[get_X1_diagnostics], chars=[label, label], fs='(t6,a)')
!
   end subroutine print_X1_diagnostics_ccs
!
!
   subroutine construct_MO_screening_for_cd_ccs(wf, screening_vector)
!!
!!    Construct MO screening for CD
!!    Written by Sarai D. Folkestad, Feb 2020
!!
!!    Constructs the screening vector
!!
!!       v_x = max_p(C_xp^2)
!!
!!       screening_vector(x,y) = v_x * v_y
!!
!!    which is used to target accuracy in MO integrals
!!    rather than the AO integrals in the decomposition
!!    of the ERIs
!!
!!    See J. Chem. Phys. 150, 194112 (2019) for further
!!    details
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: screening_vector
!
      integer :: x, p, y
!
      real(dp), dimension(:), allocatable :: v
!
      call mem%alloc(v, wf%ao%n)
      call zero_array(v, wf%ao%n)
!
!$omp parallel do private(x, p)
      do x = 1, wf%ao%n
         do p = 1, wf%n_mo
!
            if (wf%orbital_coefficients(x,p)**2 .gt. v(x)) &
               v(x) = wf%orbital_coefficients(x,p)**2
!
         enddo
!
      enddo
!$omp end parallel do
!
!$omp parallel do private(x, y)
      do x = 1, wf%ao%n
         do y = 1, wf%ao%n
!
            screening_vector(x, y) = v(x)*v(y)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(v, wf%ao%n)
!
   end subroutine construct_MO_screening_for_cd_ccs
!
!
   subroutine construct_ntos_ccs(wf, ntos, k, n_sig_v, n_sig_o, l_or_r,threshold)
!!
!!    Construct NTOs
!!    Written by Sarai D. Folkestad, May 2020
!!
!!    'ntos' : array to store NTO orbital coefficients
!!
!!    'k' : which excited state to construct NTOs for
!!
!!    'n_sig_v' and 'n_sig_o' : The number of significant occupied and virtuals
!!
!!    'l_or_r' : excited state vector type
!!
!!     Constructs the M and N matrices
!!
!!       M_ij += sum_a R_ai R_aj
!!       N_ab += sum_i R_ai R_bi
!!
!!    and diagonalizes them.
!!    The NTOs are then constructed and stored in 'ntos'
!!
!!    This constitutes a singular value decomposition
!!    of the excitation vectors.
!!
!!       N = RR^T = V D1 V^T
!!       M = R^TR = U D2 U^T
!!
!!    If a threshold is given, the significant occupied and virtual
!!    orbitals (n_o and n_v)
!!    are determined by considering the
!!    eigenvalues of M and N, using the criterion
!!
!!       1 - sum_i e_i < threshold,
!!
!!    where the eigenvalues are ordered from largest to smallest.
!    Otherwise n_sig_o = n_sig_v = 1.
!
      use nto_tool_class, only: nto_tool
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(out) :: ntos
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in) :: k
      integer, intent(out) :: n_sig_v, n_sig_o
!
      character(len=*), intent(in) :: l_or_r
!
      type(nto_tool), allocatable :: orbital_tool
!
      real(dp), dimension(:), allocatable :: X
!
      call output%printf('n', '- Constructing NTOs for state (i0)', &
                         ints=[k], ffs='(/t3,a)')
!
      call mem%alloc(X, wf%n_es_amplitudes)
!
      call wf%read_excited_state(X, k, k, trim(l_or_r))
!
      orbital_tool = nto_tool(wf%n_o, wf%n_v, wf%ao%n, wf%n_t1)
!
      call orbital_tool%initialize()
!
      call orbital_tool%add_excited_state(X(1:wf%n_t1))
!
      call orbital_tool%transform_orbitals(wf%orbital_coefficients, ntos)
!
      call orbital_tool%get_n_active_orbitals(n_sig_o, n_sig_v, threshold)
!
      call orbital_tool%cleanup
!
      call mem%dealloc(X, wf%n_es_amplitudes)
!
   end subroutine construct_ntos_ccs
!
!
   subroutine construct_ntos_or_cntos_ccs(wf, orbitals, k, n_sig_v, n_sig_o, &
                                          l_or_r, type_, threshold)
!!
!!    Construct NTOs or CNTOs
!!    Written by Sarai D. Folkestad, May 2020
!!
!!    Wrapper to construct either NTOs or CNTOs
!!
!!    orbitals : NTOs/CNTOs
!!
!!    n_sig_o, n_sig_v : number of occupied and virtual NTOs/CNTOs
!!
!!    l_or_r : use left or right excitation vectors
!!
!!    type_ : 'ntos' or 'cntos'
!!
!!    threshold : threshold to determine n_sig_o and n_sig_v
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(out) :: orbitals
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in) :: k
      integer, intent(out) :: n_sig_v, n_sig_o
!
      character(len=*), intent(in) :: l_or_r
      character(len=*), intent(in) :: type_
!
      if (trim(type_) == 'nto') then
!
         call wf%construct_ntos(orbitals, k, n_sig_v, n_sig_o, l_or_r,threshold)
!
      else if (trim(type_) == 'cnto') then
!
         call output%error_msg('Cannot plot CNTOs for ' // trim(wf%name_))
!
      endif
!
   end subroutine construct_ntos_or_cntos_ccs
!
!
   subroutine cholesky_factory_complex(wf, n_J)
!!
!!    Cholesky factory complex
!!    Written by Sarai D. Folkestad, Sep 2021
!!
!     NOTE:
!     This routine is placed here for gcc-8 compatibility,
!     should be moved to complex_ccs submodule when gcc-8 is no longer supported.
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      integer, intent(in) :: n_J
!
      character(len=200) :: cholesky_storage
!
      logical :: has_room_for_cholesky
      logical :: cholesky_in_memory
!
      has_room_for_cholesky = 2*n_J*wf%n_mo**2*(2*dp) < mem%get_available()/5
!
      cholesky_storage = 'memory'
      call input%get_keyword('cholesky storage', 'integrals', cholesky_storage)
!
      if (trim(cholesky_storage) /= 'memory' .and. trim(cholesky_storage) /= 'disk') &
         call output%error_msg('illegal value for cholesky storage')
!
      cholesky_in_memory = trim(cholesky_storage) == 'memory' .and. has_room_for_cholesky
!
      if (cholesky_in_memory) then
         wf%L_t1_c = eri_cholesky_memory_c()
         wf%L_mo_c = eri_cholesky_memory_c()
      else
         wf%L_t1_c = eri_cholesky_disk_c('T1_complex')
         wf%L_mo_c = eri_cholesky_disk_c('MO_complex')
      endif
!
      call wf%L_mo_c%initialize(n_J, 2, [wf%n_o, wf%n_v])
      call wf%L_t1_c%initialize(n_J, 2, [wf%n_o, wf%n_v])
!
      call output%printf('m', &
                         'Cholesky vectors in memory: (l0)', &
                         logs=[cholesky_in_memory],          &
                         fs='(/t6,a)')
!
   end subroutine cholesky_factory_complex
!
end module ccs_class
