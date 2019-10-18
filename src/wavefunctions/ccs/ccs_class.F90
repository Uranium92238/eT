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
   use global_in, only: input
!
   use mo_integral_tool_class, only : mo_integral_tool
!
   use reordering
!
   use sequential_file_class, only : sequential_file
   use string_utilities, only : convert_to_uppercase
   use array_utilities, only : get_l2_norm, copy_and_scale
   use array_utilities, only : get_abs_max_w_index
   use array_utilities, only : get_n_highest, get_n_lowest
   use index_invert, only : invert_compound_index, invert_packed_index
   use batching_index_class, only : batching_index
   use timings_class, only : timings
!
   implicit none
!
   type, extends(wavefunction) :: ccs
!
      integer :: n_gs_amplitudes
      integer :: n_es_amplitudes
      integer :: n_t1
      integer :: n_excited_states
      integer :: n_bath_orbitals
!
      integer :: n_core_MOs
      logical :: cvs 
!
      real(dp), dimension(:), allocatable :: left_excitation_energies, right_excitation_energies
!
      logical :: bath_orbital
      logical :: frozen_core
!
      type(sequential_file) :: t1_file, t1bar_file
      type(sequential_file) :: r1_file, l1_file
      type(sequential_file) :: excitation_energies_file
      type(sequential_file) :: restart_file
!
      type(mo_integral_tool) :: integrals
!
      real(dp) :: hf_energy
!
      real(dp), dimension(:,:), allocatable  :: t1
      real(dp), dimension(:,:), allocatable  :: t1bar
!
      real(dp), dimension(:,:), allocatable  :: fock_ij
      real(dp), dimension(:,:), allocatable  :: fock_ia
      real(dp), dimension(:,:), allocatable  :: fock_ai
      real(dp), dimension(:,:), allocatable  :: fock_ab
!
      real(dp), dimension(:,:), allocatable :: density
      real(dp), dimension(:,:), allocatable :: left_transition_density
      real(dp), dimension(:,:), allocatable :: right_transition_density
!
      integer, dimension(:), allocatable :: core_MOs
!
!     Frozen core variables (FC)
!
      real(dp), dimension(:,:), allocatable :: mo_fock_fc_contribution 
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_fc
!
      integer :: n_frozen_orbitals
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                                     => cleanup_ccs
      procedure :: general_cc_preparations                     => general_cc_preparations_ccs
!
      procedure :: read_hf                                     => read_hf_ccs
      procedure :: initialize_files                            => initialize_files_ccs
      procedure :: initialize_cc_files                         => initialize_cc_files_ccs
      procedure :: initialize_singles_files                    => initialize_singles_files_ccs
!
!     Routines related to the amplitudes & multipliers
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_ccs
      procedure :: destruct_amplitudes                         => destruct_amplitudes_ccs
      procedure :: set_initial_amplitudes_guess                => set_initial_amplitudes_guess_ccs
      procedure :: t1_transform                                => t1_transform_ccs
      procedure :: ao_to_t1_transformation                     => ao_to_t1_transformation_ccs
      procedure :: set_amplitudes                              => set_amplitudes_ccs
      procedure :: get_amplitudes                              => get_amplitudes_ccs
      procedure :: save_amplitudes                             => save_amplitudes_ccs
      procedure :: read_amplitudes                             => read_amplitudes_ccs
!
      procedure :: print_dominant_x_amplitudes                 => print_dominant_x_amplitudes_ccs
      procedure :: print_dominant_amplitudes                   => print_dominant_amplitudes_ccs
      procedure :: print_dominant_x1                           => print_dominant_x1_ccs
      procedure :: get_t1_diagnostic                           => get_t1_diagnostic_ccs
!
      procedure :: save_singles_vector                         => save_singles_vector_ccs
      procedure :: read_singles_vector                         => read_singles_vector_ccs
!
      procedure :: save_excited_state                          => save_excited_state_ccs
      procedure :: read_excited_state                          => read_excited_state_ccs
      procedure :: get_n_excited_states_on_file                => get_n_excited_states_on_file_ccs
!
      procedure :: initialize_right_excitation_energies        => initialize_right_excitation_energies_ccs
      procedure :: initialize_left_excitation_energies         => initialize_left_excitation_energies_ccs
      procedure :: destruct_right_excitation_energies          => destruct_right_excitation_energies_ccs
      procedure :: destruct_left_excitation_energies           => destruct_left_excitation_energies_ccs
      procedure :: save_excitation_energies                    => save_excitation_energies_ccs
      procedure :: read_excitation_energies                    => read_excitation_energies_ccs
      procedure :: get_n_excitation_energies_on_file           => get_n_excitation_energies_on_file_ccs
!
      procedure :: initialize_multipliers                      => initialize_multipliers_ccs
      procedure :: destruct_multipliers                        => destruct_multipliers_ccs
      procedure :: set_multipliers                             => set_multipliers_ccs
      procedure :: get_multipliers                             => get_multipliers_ccs
      procedure :: save_multipliers                            => save_multipliers_ccs
      procedure :: read_multipliers                            => read_multipliers_ccs
!
      procedure :: save_tbar_intermediates                     => save_tbar_intermediates_ccs
!
      procedure :: is_restart_safe                             => is_restart_safe_ccs
      procedure :: write_cc_restart                            => write_cc_restart_ccs
!
!     Routines related to the Fock matrix
!
      procedure :: set_fock                                    => set_fock_ccs
      procedure :: construct_fock                              => construct_fock_ccs
      procedure :: add_frozen_core_fock_contribution           => add_frozen_core_fock_contribution_ccs
      procedure :: add_molecular_mechanics_fock_contribution   => add_molecular_mechanics_fock_contribution_ccs
!
      procedure :: get_gs_orbital_differences                  => get_gs_orbital_differences_ccs
      procedure :: get_es_orbital_differences                  => get_gs_orbital_differences_ccs
!
!     Routines related to the omega vector
!
      procedure :: construct_omega                             => construct_omega_ccs
      procedure :: omega_ccs_a1                                => omega_ccs_a1_ccs
!
      procedure :: form_newton_raphson_t_estimate              => form_newton_raphson_t_estimate_ccs
!
!     Routines related to the Jacobian transformation
!
      procedure :: prepare_for_jacobian                        => prepare_for_jacobian_ccs
      procedure :: prepare_for_jacobian_transpose              => prepare_for_jacobian_transpose_ccs
      procedure :: prepare_for_multiplier_equation             => prepare_for_multiplier_equation_ccs
!
      procedure :: jacobian_transformation                     => jacobian_transformation_ccs 
      procedure :: jacobian_ccs_a1                             => jacobian_ccs_a1_ccs
      procedure :: jacobian_ccs_b1                             => jacobian_ccs_b1_ccs
!
      procedure :: jacobian_transpose_transformation           => jacobian_transpose_transformation_ccs
      procedure :: jacobian_transpose_ccs_a1                   => jacobian_transpose_ccs_a1_ccs
      procedure :: jacobian_transpose_ccs_b1                   => jacobian_transpose_ccs_b1_ccs
!
      procedure :: construct_excited_state_equation            => construct_excited_state_equation_ccs
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_ccs
      procedure :: construct_eta                               => construct_eta_ccs
!
      procedure :: get_cvs_projector                           => get_cvs_projector_ccs
      procedure :: read_cvs_settings                           => read_cvs_settings_ccs
!
      procedure :: set_cvs_start_indices                       => set_cvs_start_indices_ccs
!
!     Routines to get electron repulsion integrals (ERIs)
!
      procedure :: get_ovov                                    => get_ovov_ccs
      procedure :: get_vovo                                    => get_vovo_ccs
      procedure :: get_vvoo                                    => get_vvoo_ccs
      procedure :: get_voov                                    => get_voov_ccs
      procedure :: get_ovvo                                    => get_ovvo_ccs
      procedure :: get_oovv                                    => get_oovv_ccs
      procedure :: get_oooo                                    => get_oooo_ccs
      procedure :: get_vvvv                                    => get_vvvv_ccs
      procedure :: get_ooov                                    => get_ooov_ccs
      procedure :: get_oovo                                    => get_oovo_ccs
      procedure :: get_ovoo                                    => get_ovoo_ccs
      procedure :: get_vooo                                    => get_vooo_ccs
      procedure :: get_vvvo                                    => get_vvvo_ccs
      procedure :: get_vvov                                    => get_vvov_ccs
      procedure :: get_vovv                                    => get_vovv_ccs
      procedure :: get_ovvv                                    => get_ovvv_ccs
!
      procedure :: get_g_pqrs_required                         => get_g_pqrs_required_ccs
!
      procedure, nopass :: need_g_abcd                         => need_g_abcd_ccs
!
!     Routines to initialize and destruct arrays
!
      procedure :: initialize_fock                             => initialize_fock_ccs
      procedure :: initialize_fock_ij                          => initialize_fock_ij_ccs
      procedure :: initialize_fock_ia                          => initialize_fock_ia_ccs
      procedure :: initialize_fock_ai                          => initialize_fock_ai_ccs
      procedure :: initialize_fock_ab                          => initialize_fock_ab_ccs
      procedure :: initialize_t1                               => initialize_t1_ccs
      procedure :: initialize_t1bar                            => initialize_t1bar_ccs
      procedure :: initialize_gs_density                       => initialize_gs_density_ccs
      procedure :: initialize_transition_densities             => initialize_transition_densities_ccs
      procedure :: initialize_core_MOs                         => initialize_core_MOs_ccs 
!
      procedure :: destruct_fock_ij                            => destruct_fock_ij_ccs
      procedure :: destruct_fock_ia                            => destruct_fock_ia_ccs
      procedure :: destruct_fock_ai                            => destruct_fock_ai_ccs
      procedure :: destruct_fock_ab                            => destruct_fock_ab_ccs
      procedure :: destruct_t1                                 => destruct_t1_ccs
      procedure :: destruct_t1bar                              => destruct_t1bar_ccs
      procedure :: destruct_gs_density                         => destruct_gs_density_ccs
      procedure :: destruct_transition_densities               => destruct_transition_densities_ccs
      procedure :: destruct_core_MOs                           => destruct_core_MOs_ccs 
!
!     Routines related to EOM first order property calculations
!
      procedure :: construct_etaX                              => construct_etaX_ccs
      procedure :: construct_eom_etaX                          => construct_eom_etaX_ccs
      procedure :: etaX_ccs_a1                                 => etaX_ccs_a1_ccs
      procedure :: etaX_ccs_b1                                 => etaX_ccs_b1_ccs
!
      procedure :: construct_csiX                              => construct_csiX_ccs
      procedure :: csiX_ccs_a1                                 => csiX_ccs_a1_ccs
!
      procedure :: etaX_eom_a                                  => etaX_eom_a_ccs
!
      procedure :: calculate_transition_strength               => calculate_transition_strength_ccs
!
!     Routines related to one-electron densities
!
      procedure :: prepare_for_density                         => prepare_for_density_ccs
!
      procedure :: construct_gs_density                        => construct_gs_density_ccs
      procedure :: construct_right_transition_density          => construct_right_transition_density_ccs
      procedure :: construct_left_transition_density           => construct_left_transition_density_ccs
!
      procedure :: gs_one_el_density_ccs_oo                    => gs_one_el_density_ccs_oo_ccs
      procedure :: gs_one_el_density_ccs_vo                    => gs_one_el_density_ccs_vo_ccs
!
      procedure :: right_transition_density_ccs_oo             => right_transition_density_ccs_oo_ccs
      procedure :: right_transition_density_ccs_ov             => right_transition_density_ccs_ov_ccs
      procedure :: right_transition_density_ccs_vv             => right_transition_density_ccs_vv_ccs
      procedure :: right_transition_density_ccs_gs_contr       => right_transition_density_ccs_gs_contr_ccs
!
      procedure :: binormalize_L_wrt_R                         => binormalize_L_wrt_R_ccs
!
!     One-electron operators and mean value
!
      procedure :: construct_h                                 => construct_h_ccs 
      procedure :: construct_mu                                => construct_mu_ccs 
      procedure :: construct_q                                 => construct_q_ccs 
!
      procedure :: calculate_expectation_value                 => calculate_expectation_value_ccs
      procedure :: calculate_energy                            => calculate_energy_ccs
!
      procedure :: construct_molecular_gradient                => construct_molecular_gradient_ccs
!
      procedure :: read_settings                               => read_settings_ccs
      procedure :: make_bath_orbital                           => make_bath_orbital_ccs
!
      procedure :: set_ip_start_indices                        => set_ip_start_indices_ccs
      procedure :: get_ip_projector                            => get_ip_projector_ccs
!
      procedure :: approximate_double_excitation_vectors       => approximate_double_excitation_vectors_ccs
!
!     Frozen core
!
      procedure :: initialize_mo_fock_fc_contribution          => initialize_mo_fock_fc_contribution_ccs                
      procedure :: destruct_mo_fock_fc_contribution            => destruct_mo_fock_fc_contribution_ccs                
      procedure :: initialize_orbital_coefficients_fc          => initialize_orbital_coefficients_fc_ccs
      procedure :: destruct_orbital_coefficients_fc            => destruct_orbital_coefficients_fc_ccs
!
      procedure :: coulomb_contribution_fock_fc                => coulomb_contribution_fock_fc_ccs
      procedure :: exchange_contribution_fock_fc               => exchange_contribution_fock_fc_ccs
      procedure :: remove_core_orbitals                        => remove_core_orbitals_ccs
      procedure :: construct_mo_fock_fc_contribution           => construct_mo_fock_fc_contribution_ccs
      procedure :: construct_t1_fock_fc_contribution           => construct_t1_fock_fc_contribution_ccs
!
!     MO preparations
!
      procedure :: mo_preparations                             => mo_preparations_ccs
!   
!     Debug 
!
      procedure :: omega_for_jacobian_debug                    => omega_for_jacobian_debug_ccs
      procedure :: amplitudes_for_jacobian_debug               => amplitudes_for_jacobian_debug_ccs
      procedure :: normalization_for_jacobian_debug            => normalization_for_jacobian_debug_ccs
      procedure :: numerical_test_jacobian                     => numerical_test_jacobian_ccs
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
      include "fock_ccs_interface.F90"
      include "debug_jacobian_ccs_interface.F90"
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
!
      call wf%write_cc_restart()
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
      wf%system => system
!
!     Initialize CC files
!
      call wf%initialize_files()
!
!     Read necessary information from HF
!
      call wf%read_hf()
!
!     Set orbital coefficients and energies
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%read_orbital_coefficients()
      call wf%read_orbital_energies()
!
!     Logicals for special methods
!
      wf%bath_orbital = .false.
      wf%frozen_core = .false.
      wf%cvs = .false.
!
!     Read CC settings from eT.inp (from cc section)
!
      call wf%read_settings()
!
!     Handle changes in the number of MOs as a result of 
!     special methods
!
      if (wf%bath_orbital) call wf%make_bath_orbital()
      if (wf%frozen_core) call wf%remove_core_orbitals()
!
   end subroutine general_cc_preparations_ccs
!
!
   subroutine cleanup_ccs(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad and
!!    Alexander Paul , 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%destruct_amplitudes()
      call wf%destruct_multipliers()
      call wf%destruct_right_excitation_energies()
      call wf%destruct_left_excitation_energies()
      call wf%destruct_mo_fock_fc_contribution()
!
      write(output%unit, '(/t3,a,a,a)') '- Cleaning up ', trim(convert_to_uppercase(wf%name_)), ' wavefunction'
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
      if (input%requested_section('cc')) then
!
         if (input%requested_keyword_in_section('bath orbital','cc')) then 
!
            wf%bath_orbital = .true.
            wf%n_bath_orbitals = 1 ! May want to expand the number of bath orbitals later on
!
         endif
!
         if (input%requested_keyword_in_section('frozen core', 'cc')) wf%frozen_core = .true.
!
      endif
!
   end subroutine read_settings_ccs
!
!
   subroutine read_hf_ccs(wf)
!!
!!    Read HF file
!!    Written by Rolf Heilemann Myhre, May 2019
!!    Short routine to read the HF information from disk
!!
      implicit none
!
      class(ccs) :: wf
!
      type(sequential_file) :: hf_restart_file 
!
      hf_restart_file = sequential_file('hf_restart_file')
      call hf_restart_file%open_('read', 'rewind')
!
      call hf_restart_file%read_(wf%n_ao)     
      call hf_restart_file%read_(wf%n_mo)     
      call hf_restart_file%skip()     
      call hf_restart_file%read_(wf%n_o)     
      call hf_restart_file%read_(wf%n_v)     
      call hf_restart_file%read_(wf%hf_energy)    
!
      call hf_restart_file%close_()
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
   subroutine construct_excited_state_equation_ccs(wf, X, R, w, r_or_l)
!!
!!    Construct excited state equation
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Constructs R = AX - wX, where w = X^T A X and norm(X) = sqrt(X^T X) = 1
!!
!!    Note I: we assume that X is normalized. If it is not,
!!    please normalize before calling the routine.
!!
!!    Note II: this routine constructs the excited state equation
!!    for standard CC models and the effective (!) excited state
!!    equation in perturbative models. In the CC2 routine, for
!!    instance, X and R will be n_o*n_v vectors and A(w) will
!!    depend on the excitation energy w. See, e.g., Weigend and
!!    Hättig's RI-CC2 paper for more on this topic. This means
!!    that w should be the previous w-value when entering the
!!    routine (so that A(w)X may be constructed approximately)
!!    in perturbative models.
!!
!!    Note III: the routine is used by the DIIS excited state solver.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: R
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), intent(inout) :: w
!
      real(dp), dimension(:), allocatable :: X_copy
!
      real(dp) :: ddot
!
      call mem%alloc(X_copy, wf%n_es_amplitudes)
      call dcopy(wf%n_es_amplitudes, X, 1, X_copy, 1)
!
      if (r_or_l .eq. "right") then
!
         call wf%jacobian_transformation(X_copy) ! X_copy <- AX
!
      elseif (r_or_l .eq. 'left') then
!
         call wf%jacobian_transpose_transformation(X_copy) ! X_copy <- XA
!
      else
!
         call output%error_msg('Neither left nor right in construct_excited_state')
!
      endif
!
      w = ddot(wf%n_es_amplitudes, X, 1, X_copy, 1)
!
      call dcopy(wf%n_es_amplitudes, X_copy, 1, R, 1)
      call daxpy(wf%n_es_amplitudes, -w, X, 1, R, 1)
!
      call mem%dealloc(X_copy, wf%n_es_amplitudes)
!
   end subroutine construct_excited_state_equation_ccs
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
      write(output%unit, '(/t6,a)') 'Largest single amplitudes:'
      write(output%unit, '(t6,a)')  '-----------------------------------'
      write(output%unit, '(t6,a)')  '   a       i         ' // tag // '(a,i)             '
      write(output%unit, '(t6,a)')  '-----------------------------------'
!
      do elm = 1, n_elements
!
         call invert_compound_index(dominant_indices(elm), a, i, wf%n_v, wf%n_o)
!
         write(output%unit, '(t6,i4,4x,i3,3x,f19.12)') a, i, x1(dominant_indices(elm))
!
      enddo
!
      write(output%unit, '(t6,a)')  '------------------------------------'
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
      integer, dimension(wf%n_excited_states) :: start_indices
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
            if (current_root .eq. wf%n_excited_states) then
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
   subroutine binormalize_L_wrt_R_ccs(wf, L, R, state)
!!
!!    Calculates the overlap of the left and right states 
!!    and scales the left amplitudes by it
!!    Written by Josefine Andersen, Apr 2019
!!
!!    Consistency/sanity check: Eirik F. Kjønstad, Aug 2019
!!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: R
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: L
!
      integer, intent(in) :: state
!
      real(dp) :: ddot, LT_R
      real(dp) :: energy_threshold
!
!     Sanity check in case roots are ordered incorrectly
!
      if (input%requested_keyword_in_section('energy threshold', 'solver cc es')) then 
!
        call input%get_keyword_in_section('energy threshold', 'solver cc es', energy_threshold)
!
      elseif (input%requested_keyword_in_section('residual threshold', 'solver cc es')) then 
!
        call input%get_keyword_in_section('residual threshold', 'solver cc es', energy_threshold)
!
      else
!
        call output%printf('Note: assuming default energy threshold (1.0d-6) when testing root consistency.', fs='(t6,a)')
!
        energy_threshold = 1.0d-6
!
      endif 
!
      if (abs(wf%left_excitation_energies(state) - wf%right_excitation_energies(state)) > energy_threshold) then 
!
          call output%printf('Eigenvector (i0) is not left-right consistent to threshold (e8.2).', &
                              ints=[state], reals=[energy_threshold], fs='(/t6,a)')
!
          call output%printf('Energies (left, right): (f19.12) (f19.12)', &
                reals=[wf%left_excitation_energies(state), wf%right_excitation_energies(state)], fs='(/t6,a)')
!
          call output%error_msg('while calculating transition strength.')
!
      !  else
!
         ! Note: consider to add in verbose mode
         ! call output%printf('The left and right states corresponding to root (i0) are consistent', ints=[state])
!
      endif 
!
      LT_R = ddot(wf%n_es_amplitudes, L, 1, R, 1)
      call dscal(wf%n_es_amplitudes, one/LT_R, L, 1)
!
   end subroutine binormalize_L_wrt_R_ccs
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
      integer, dimension(wf%n_excited_states), intent(out) :: start_indices
!
      integer :: I
!
      do I = 1, wf%n_excited_states
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                  :: R_ai
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
   subroutine remove_core_orbitals_ccs(wf)
!!
!!    Remove core orbitals
!!    Written by Sarai D. Folkestad, Sep 2018 
!!
!!    Removes 1s orbitals for C and heavier atoms
!!    from orbital coefficient matrix
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_copy
      real(dp), dimension(:), allocatable :: orbital_energies_copy
!
      integer :: I
!
      logical, dimension(:), allocatable :: freeze_atom
!
      integer :: index_max
      real(dp) :: max_
!
      call mem%alloc(freeze_atom, wf%system%n_atoms)
      freeze_atom = .false.
!
!     Figure out how many core orbitals we have
!
!     Number of atoms heavier than Be (atomic number larger than 4)
!
      wf%n_frozen_orbitals = 0
!
      do I = 1, wf%system%n_atoms
!
         if (wf%system%atoms(I)%number_ .ge. 5 .and. wf%system%atoms(I)%number_ .le. 12) then
!
            wf%n_frozen_orbitals = wf%n_frozen_orbitals + 1
            freeze_atom(I) = .true.
!
         elseif (wf%system%atoms(I)%number_ .ge. 13 .and. wf%system%atoms(I)%number_ .le. 30) then
!
            wf%n_frozen_orbitals = wf%n_frozen_orbitals + 5
            freeze_atom(I) = .true.
!
         elseif (wf%system%atoms(I)%number_ .gt. 30) then
!
            call output%error_msg('No frozen core for Z > 30.')
!
         endif
!
      enddo
!
      call mem%alloc(orbital_coefficients_copy, wf%n_ao, wf%n_mo)
!
      call copy_and_scale(one, wf%orbital_coefficients, orbital_coefficients_copy, wf%n_mo*wf%n_ao)
!
      call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
      call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo - wf%n_frozen_orbitals)
!
      wf%orbital_coefficients(1:wf%n_ao, 1:wf%n_mo - wf%n_frozen_orbitals) = &
            orbital_coefficients_copy(1:wf%n_ao, wf%n_frozen_orbitals + 1:wf%n_mo)
!
!     Keep frozen core orbital orbital coefficients
!
      call wf%initialize_orbital_coefficients_fc()
!
      wf%orbital_coefficients_fc(1:wf%n_ao, 1:wf%n_frozen_orbitals) = &
            orbital_coefficients_copy(1:wf%n_ao, 1:wf%n_frozen_orbitals) 
!
!     Check for crossover:  
!
!     If the largest AO weight on a frozen MO does not belong to an 
!     atom we are supposed to freeze, then we stop.
!
      do i = 1, wf%n_frozen_orbitals
!
         call get_abs_max_w_index(wf%orbital_coefficients_fc(:,i), wf%n_ao, max_, index_max)
!
         if (.not. freeze_atom(wf%system%basis2atom(index_max))) &
            call output%error_msg('Detected crossover in frozen core.')
!
      enddo
!
     call mem%dealloc(orbital_coefficients_copy, wf%n_ao, wf%n_mo)
!
     call mem%alloc(orbital_energies_copy, wf%n_mo)
!
     call copy_and_scale(one, wf%orbital_energies, orbital_energies_copy, wf%n_mo)
!
     call mem%dealloc(wf%orbital_energies, wf%n_mo)
     call mem%alloc(wf%orbital_energies, wf%n_mo - wf%n_frozen_orbitals)
!
     wf%orbital_energies(1:wf%n_mo - wf%n_frozen_orbitals) = &
            orbital_energies_copy(wf%n_frozen_orbitals + 1: wf%n_mo)
!
     call mem%dealloc(orbital_energies_copy, wf%n_mo)
!
     wf%n_mo = wf%n_mo  - wf%n_frozen_orbitals
     wf%n_o = wf%n_o  - wf%n_frozen_orbitals     
!
   end subroutine remove_core_orbitals_ccs
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
      if (wf%frozen_core) call wf%construct_mo_fock_fc_contribution()
!
   end subroutine mo_preparations_ccs
!
!
end module ccs_class
