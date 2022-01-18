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
module cc3_class
!
!!
!!    Coupled cluster CC3 class module
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2018-2019
!!
!
   use triples_class, only: triples
!
   use parameters
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use stream_file_class, only: stream_file
   use direct_stream_file_class, only: direct_stream_file
   use batching_index_class, only: batching_index
   use range_class, only: range_
!
   use eri_adapter_class, only: eri_adapter
   use eri_cholesky_disk_class, only: eri_cholesky_disk
!
   implicit none
!
   type, extends(triples) :: cc3
!
!     Jacobian intermediates files
!
!     X_abid = - sum_jck u^abc_ijk g_kcjd
!     X_ajil = - sum_bck u^abc_ijk g_lbkc
      type(direct_stream_file) :: X_dbai
      type(direct_stream_file) :: X_abid
      type(direct_stream_file) :: X_ajil
!
!     Y_ebck = sum_aij L^abc_ijk t^ae_ij
      type(direct_stream_file) :: Y_ebck
!
!     Density intermediates files
!     Y_cmjk = sum_aij tbar^abc_ijk t^ab_im
      type(direct_stream_file) :: Y_cmjk_tbar, Y_ebck_tbar
!
      real(dp), dimension(:,:), allocatable :: GS_cc3_density_oo
      real(dp), dimension(:,:), allocatable :: GS_cc3_density_vv
!
      real(dp), dimension(:,:), allocatable :: L_cc3_density_oo
      real(dp), dimension(:,:), allocatable :: L_cc3_density_vv
!
      type(eri_cholesky_disk), allocatable :: L_c1
      type(eri_adapter), allocatable :: eri_c1
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                              => cleanup_cc3
      procedure :: delete_intermediate_files            => delete_intermediate_files_cc3
!
      procedure :: one_core_index                       => one_core_index_cc3
      procedure :: two_core_indices                     => two_core_indices_cc3
      procedure :: three_core_indices                   => three_core_indices_cc3
      procedure :: ijk_amplitudes_are_zero              => ijk_amplitudes_are_zero_cc3
!
!     Routines related to omega
!
      procedure :: construct_omega                      => construct_omega_cc3
      procedure :: omega_cc3_a                          => omega_cc3_a_cc3
      procedure :: omega_cc3_contractions               => omega_cc3_contractions_cc3
      procedure :: omega1_cc3_permutation               => omega1_cc3_permutation_cc3
      procedure :: omega2_cc3_permutation               => omega2_cc3_permutation_cc3
      procedure :: omega2_fock_cc3_permutation          => omega2_fock_cc3_permutation_cc3
!
      procedure :: estimate_mem_integral_setup          => estimate_mem_integral_setup_cc3
      procedure :: estimate_mem_c1_integral_setup       => estimate_mem_c1_integral_setup_cc3
!
!     Routines used for prepare, both left and right
!
      procedure :: prepare_for_jacobian                 => prepare_for_jacobian_cc3
      procedure :: prepare_for_jacobian_transpose       => prepare_for_jacobian_transpose_cc3
      procedure :: prepare_cc3_jacobian_intermediates   => prepare_cc3_jacobian_intermediates_cc3
      procedure :: construct_x_intermediates            => construct_x_intermediates_cc3
      procedure :: sort_x_to_abid_and_write             => sort_x_to_abid_and_write_cc3
!
!     Routines related to the jacobian
!
      procedure :: construct_Jacobian_transform         => construct_Jacobian_transform_cc3
!
      procedure :: prepare_for_approximate_Jacobians    => prepare_for_approximate_Jacobians_cc3
      procedure :: approximate_Jacobian_transform       => approximate_Jacobian_transform_cc3
!
!     Right hand side transformation
!
      procedure :: effective_jacobian_transformation    => effective_jacobian_transformation_cc3
!
      procedure :: jacobian_cc3_t3_a2                   => jacobian_cc3_t3_a2_cc3
      procedure :: jacobian_cc3_t3_b2                   => jacobian_cc3_t3_b2_cc3
      procedure :: construct_c1_fock                    => construct_c1_fock_cc3
      procedure :: jacobian_cc3_b2_fock                 => jacobian_cc3_b2_fock_cc3
      procedure :: jacobian_cc3_c3_a                    => jacobian_cc3_c3_a_cc3
      procedure :: rho2_fock_cc3_permutation            => rho2_fock_cc3_permutation_cc3
!
      procedure :: initialize_eri_c1 &
                => initialize_eri_c1_cc3
      procedure :: destruct_eri_c1 &
                => destruct_eri_c1_cc3
!
      procedure :: construct_c1_cholesky     &
                => construct_c1_cholesky_cc3
      procedure :: construct_cholesky_c1_oo  &
                => construct_cholesky_c1_oo_cc3
      procedure :: construct_cholesky_c1_vo  &
                => construct_cholesky_c1_vo_cc3
      procedure :: construct_cholesky_c1_vv  &
                => construct_cholesky_c1_vv_cc3
!
!     Routines related to the transpose of the jacobian
!
      procedure :: effective_jacobian_transpose_transformation  &
                => effective_jacobian_transpose_transformation_cc3
!
      procedure :: jacobian_transpose_cc3_t3_a1         => jacobian_transpose_cc3_t3_a1_cc3
      procedure :: jacobian_transpose_cc3_t3_b1         => jacobian_transpose_cc3_t3_b1_cc3
      procedure :: jacobian_transpose_cc3_b1_x_ai       => jacobian_transpose_cc3_b1_x_ai_cc3
      procedure :: construct_x_ai_intermediate          => construct_x_ai_intermediate_cc3
      procedure :: jacobian_transpose_cc3_c3_a          => jacobian_transpose_cc3_c3_a_cc3
      procedure :: outer_product_terms_l3               => outer_product_terms_l3_cc3
      procedure :: jacobian_transpose_cc3_contractions  => jacobian_transpose_cc3_contractions_cc3
      procedure :: jacobian_transpose_cc3_c3_a1_y_o     => jacobian_transpose_cc3_c3_a1_y_o_cc3
      procedure :: jacobian_transpose_cc3_c3_a1_y_v     => jacobian_transpose_cc3_c3_a1_y_v_cc3
      procedure :: construct_Y_vvvo_permutation &
                => construct_Y_vvvo_permutation_cc3
      procedure :: construct_Y_vooo_permutation &
                => construct_Y_vooo_permutation_cc3
      procedure :: construct_X_vo_permutation &
                => construct_X_vo_permutation_cc3
!
!     Routines related to the multipliers
!
      procedure :: prepare_for_multiplier_equation      => prepare_for_multiplier_equation_cc3
      procedure :: construct_multiplier_equation        => construct_multiplier_equation_cc3
      procedure :: save_tbar_intermediates              => save_tbar_intermediates_cc3
!
!     Routines to construct triples amplitudes in batches of a,b,c
!
      procedure :: setup_oovo_abc                       => setup_oovo_abc_cc3
      procedure :: point_oovo_abc                       => point_oovo_abc_cc3
      procedure :: setup_ooov_abc                       => setup_ooov_abc_cc3
      procedure :: point_ooov_abc                       => point_ooov_abc_cc3
      procedure :: setup_vvvo_abc                       => setup_vvvo_abc_cc3
      procedure :: point_vvvo_abc                       => point_vvvo_abc_cc3
      procedure :: setup_vvov_abc                       => setup_vvov_abc_cc3
      procedure :: point_vvov_abc                       => point_vvov_abc_cc3
      procedure :: setup_ovov_abc                       => setup_ovov_abc_cc3
      procedure :: point_ovov_abc                       => point_ovov_abc_cc3
!
      procedure :: estimate_mem_integral_setup_abc      => estimate_mem_integral_setup_abc_cc3
      procedure :: estimate_mem_c1_integral_setup_abc   => estimate_mem_c1_integral_setup_abc_cc3
!
      procedure :: construct_W_abc                      => construct_W_abc_cc3
      procedure :: outer_product_terms_l3_abc           => outer_product_terms_l3_abc_cc3
      procedure :: divide_by_orbital_differences_abc    => divide_by_orbital_differences_abc_cc3
!
!     Routines related to density matrices
!
      procedure :: initialize_gs_density                => initialize_gs_density_cc3
      procedure :: destruct_gs_density                  => destruct_gs_density_cc3
!
      procedure :: initialize_density_intermediates &
                => initialize_density_intermediates_cc3
      procedure :: destruct_density_intermediates &
                => destruct_density_intermediates_cc3
!
      procedure :: mu_ref_density_terms &
                => mu_ref_density_terms_cc3
      procedure :: density_cc3_mu_ref_blocks &
                => density_cc3_mu_ref_blocks_cc3
!
      procedure :: mu_nu_density_terms   &
                => mu_nu_density_terms_cc3
      procedure :: density_cc3_mu_nu_blocks &
                => density_cc3_mu_nu_blocks_cc3
!
      procedure :: density_cc3_mu_ref_abc               => density_cc3_mu_ref_abc_cc3
      procedure :: density_cc3_mu_ref_oo                => density_cc3_mu_ref_oo_cc3
      procedure :: density_cc3_mu_ref_ijk               => density_cc3_mu_ref_ijk_cc3
      procedure :: density_cc3_mu_ref_vv                => density_cc3_mu_ref_vv_cc3
      procedure :: construct_y_vooo_intermediate &
                => construct_y_vooo_intermediate_cc3
!
      procedure :: density_cc3_mu_nu_ov                 => density_cc3_mu_nu_ov_cc3
      procedure :: density_cc3_mu_nu_ijk                => density_cc3_mu_nu_ijk_cc3
      procedure :: density_cc3_mu_nu_abc                => density_cc3_mu_nu_abc_cc3
      procedure :: density_cc3_mu_nu_ov_Z_term          => density_cc3_mu_nu_ov_Z_term_cc3
      procedure :: density_cc3_Z1_ov                    => density_cc3_Z1_ov_cc3
      procedure :: density_cc3_Z2_oo_vv                 => density_cc3_Z2_oo_vv_cc3
      procedure :: density_cc3_Y_vooo_ov                => density_cc3_Y_vooo_ov_cc3
      procedure :: density_cc3_Y_vvvo_ov                => density_cc3_Y_vvvo_ov_cc3
!
      procedure, nopass :: density_cc3_oo_vv_N7_TN_permutation &
                        => density_cc3_oo_vv_N7_TN_permutation_cc3
      procedure, nopass :: density_cc3_oo_vv_N7_NT_permutation &
                        => density_cc3_oo_vv_N7_NT_permutation_cc3
!
!     Routines related to the biorthonormalization
!
      procedure :: L_R_overlap                          => L_R_overlap_cc3
      procedure :: L_R_overlap_triples                  => L_R_overlap_triples_cc3
!
!     Initialize wavefunction
!
      procedure :: initialize                           => initialize_cc3
!
   end type cc3
!
!
   interface
!
      include "initialize_destruct_cc3_interface.F90"
      include "omega_cc3_interface.F90"
      include "multiplier_equation_cc3_interface.F90"
      include "prepare_jacobian_transform_cc3_interface.F90"
      include "jacobian_cc3_interface.F90"
      include "jacobian_transpose_cc3_interface.F90"
      include "abc_batching_cc3_interface.F90"
      include "mean_value_cc3_interface.F90"
      include "response_cc3_interface.F90"
!
   end interface
!
!
contains
!
!
   subroutine initialize_cc3(wf, template_wf)
!!
!!    Initialize
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2018
!!
      use wavefunction_class,       only: wavefunction
      use citation_class,           only: citation
      use citation_printer_class,   only: eT_citations
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      type(citation), allocatable :: reference
!
      wf%name_ = 'cc3'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%n_t1 = wf%n_o*wf%n_v
      wf%n_t2 = wf%n_o*wf%n_v*(wf%n_o*wf%n_v + 1)/2
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
      wf%need_g_abcd     = .true.
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
      reference = citation(implementation = 'CC3',                                          &
                           journal        = 'J. Chem. Theory Comput.',                      &
                           title_         = 'New and Efficient Implementation of CC3',      &
                           volume         = '17',                                           &
                           issue          = '1',                                            &
                           pages          = '117–126',                                      &
                           year           = '2021',                                         &
                           doi            = '10.1021/acs.jctc.0c00686',                     &
                           authors        = [character(len=25) :: 'Alexander C. Paul',      &
                                                                  'Rolf H. Myhre',          &
                                                                  'Henrik Koch'])
!
      call eT_citations%add(reference)
!
   end subroutine initialize_cc3
!
!
   subroutine cleanup_cc3(wf)
!!
!!    Cleanup
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      call wf%ccsd%cleanup()
!
      call wf%delete_intermediate_files
!
   end subroutine cleanup_cc3
!
!
   subroutine delete_intermediate_files_cc3(wf)
!!
!!    Delete intermediate files
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Check which files exist and delete
!!    (depends on type of calculation)
!!
      implicit none
!
      class(cc3) :: wf
!
!     Delete Intermediate files for the jacobian transformations
      if (wf%X_abid%exists()) then
!
         call wf%X_abid%delete_
         call wf%X_ajil%delete_
!
      end if
!
!     Delete additional files for jacobian transpose transformation
      if (wf%Y_ebck%exists()) then
!
         call wf%Y_ebck%delete_
!
      end if
!
!     Intermediates created for/in mean_value
      if (wf%Y_cmjk_tbar%exists()) then
!
         call wf%Y_cmjk_tbar%delete_
!
      end if
!
!     Intermediates created for response
      if (wf%Y_ebck_tbar%exists()) then
!
         call wf%Y_ebck_tbar%delete_
!
      end if
!
   end subroutine delete_intermediate_files_cc3
!
!
   subroutine construct_Jacobian_transform_cc3(wf, r_or_l, X, R, w)
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
!!    X: On input contains the vector to transform.
!!    R: On output contains the transformed vector
!!
!!    w: Excitation energy. Only used for debug prints for CCS, CCSD etc.
!!       but is passed to the effective_jacobian_transform for lowmem_CC2 and CC3
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: R
!
      real(dp), intent(in), optional :: w
!
      if (present(w)) then
         call output%printf('debug', 'Calling Jacobian (a0) transform with &
                            &energy: (f19.12)', chars=[r_or_l], reals=[w])
      else
!
         call output%error_msg('w is missing in construct_Jacobian_transform for cc3')
!
      endif
!
!     Compute the transformed matrix
      if (r_or_l .eq. "right") then
!
         call wf%effective_jacobian_transformation(w, X, R)
!
      else if (r_or_l .eq. "left") then
!
         call wf%effective_jacobian_transpose_transformation(w, X, R, wf%cvs, wf%rm_core)
!
      else
!
         call output%error_msg('Neither left nor right in construct_Jacobian_transform')
!
      end if
!
   end subroutine construct_Jacobian_transform_cc3
!
!
   subroutine approximate_Jacobian_transform_cc3(wf, r_or_l, X, R, w)
!!
!!    Approximate Jacobian transform
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for a lower-level Jacobian transformation that is the best approximation
!!    with a lower computational scaling.
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: R
!
      real(dp), intent(in), optional :: w
!
      call wf%ccsd%construct_Jacobian_transform(r_or_l, X, R, w)
!
   end subroutine approximate_Jacobian_transform_cc3
!
!
   subroutine prepare_for_approximate_Jacobians_cc3(wf, r_or_l)
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
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      call wf%ccsd%prepare_for_Jacobians(r_or_l)
!
   end subroutine prepare_for_approximate_Jacobians_cc3
!
!
   subroutine save_tbar_intermediates_cc3(wf)
!!
!!    Save tbar intermediates
!!    Written by Alexander C. Paul, August 2019
!!
!!    Modified by Rolf H. Myhre Feb. 2020
!!    copy() should not be used
!!
!!    Copy the intermediate Y_ebck constructed as follows:
!!    Y_ebck = sum_aij tbar^abc_ijk * t^ae_ij
!!    Later used in the right transition density matrix
!!
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(cc3) :: wf
!
      type(batching_index) :: batch_k
      integer :: k_batch, req_0, req_1
!
      real(dp), dimension(:,:), allocatable :: Y_ebck
!
      wf%Y_ebck_tbar = direct_stream_file('Y_ebck_tbar', wf%n_v**3)
!
!     Delete if the file already exists
!     e.g. when restarting a crashed CC3 calculation
!
      if(wf%Y_ebck_tbar%exists()) then
         call wf%Y_ebck_tbar%delete_
      end if
!
!     Want to read and write as much as possible at once:
!
      req_0 = 0
      req_1 = wf%n_v**3
!
      batch_k = batching_index(wf%n_o)
      call mem%batch_setup(batch_k, req_0, req_1, 'save_tbar_intermediate_cc3')
      call mem%alloc(Y_ebck, wf%n_v**3, batch_k%max_length)
!
      call wf%Y_ebck_tbar%open_('write')
      call wf%Y_ebck%open_('read')
!
      do k_batch = 1, batch_k%num_batches
         call batch_k%determine_limits(k_batch)
!
         call wf%Y_ebck%read_range(Y_ebck, batch_k)
         call wf%Y_ebck_tbar%write_range(Y_ebck, batch_k)
!
      enddo
!
      call wf%Y_ebck%close_()
      call wf%Y_ebck_tbar%close_()
!
      call mem%dealloc(Y_ebck, wf%n_v**3, batch_k%max_length)
!
      call mem%batch_finalize()
!
   end subroutine save_tbar_intermediates_cc3
!
!
   function L_R_overlap_cc3(wf, L, left_state, R, right_state) result(L_R_overlap)
!!
!!    Left right overlap
!!    Written by Alexander C. Paul, Aug 2019
!!
!!    Computes L^T * R for full space L and R (singles, doubles, triples)
!!
      use array_utilities, only: copy_and_scale, scale_diagonal
      use reordering, only: squareup_and_sort_1234_to_1324, add_1243_to_1234
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L
      integer, intent(in) :: left_state
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R
      integer, intent(in) :: right_state
!
      real(dp), dimension(:,:,:,:), allocatable :: L2, R2
!
      real(dp) :: ddot, L_R_overlap
!
      L_R_overlap = ddot(wf%n_es_amplitudes, L, 1, R, 1)
!
      call wf%construct_c1_cholesky(R(1:wf%n_t1), wf%L_t1, wf%L_c1)
!
      call mem%alloc(L2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(R2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Construct covariant L2 - use R2 as help array
!
      call squareup_and_sort_1234_to_1324(L(wf%n_t1 + 1 : wf%n_es_amplitudes), R2, &
                                          wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call copy_and_scale(two*third, R2, L2, wf%n_t1**2)
      call add_1243_to_1234(third, R2, L2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call squareup_and_sort_1234_to_1324(R(wf%n_t1 + 1 : wf%n_es_amplitudes), R2, &
                                          wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Scale the right doubles vector by 1 + delta_ai,bj
!
      call scale_diagonal(two, R2, wf%n_v, wf%n_o)
!
      call wf%L_R_overlap_triples(wf%left_excitation_energies(left_state),  &
                                  wf%right_excitation_energies(right_state),&
                                  L(1:wf%n_t1), L2, R2, L_R_overlap, wf%cvs, wf%rm_core)
!
      call mem%dealloc(L2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(R2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call output%printf('debug', 'Overlap of (i0). left and (i0). right state: (f15.10)', &
                         ints=[left_state, right_state], &
                         reals=[L_R_overlap], fs='(/t6,a)')
!
   end function L_R_overlap_cc3
!
!
   subroutine L_R_overlap_triples_cc3(wf, omega_L, omega_R, L1, L2, R2, &
                                      LT_R, cvs, rm_core)
!!
!!    Left right overlap triples
!!    Written by Alexander C. Paul, August 2019
!!
!!    Calculates the overlap of the triples part
!!    of the left and right excitation vectors
!!
!!    Right excitation vector:
!!       R_mu3 = (omega - eps_mu3)^-1 (< mu3| [H,R_2] |HF >
!!                                   + < mu3| [[H,R_1],T_2] |HF >)
!!
!!    Left excitation vector:
!!       L_mu3 = (omega - eps^abc_ijk)^-1 (L_mu1 < mu1| [H,tau_nu3] |R >
!!                                       + L_mu2 < mu2| [H,tau_nu3] |R >)
!!
!
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3) :: wf
!
      logical, intent(in) :: cvs, rm_core
!
      real(dp), intent(in) :: omega_L, omega_R
!
      real(dp), intent(inout) :: LT_R
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: L1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: L2, R2
!
      real(dp), dimension(:,:,:,:), allocatable :: t2
!
      real(dp), dimension(:,:,:), allocatable :: R_abc
      real(dp), dimension(:,:,:), allocatable :: L_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
!     Integrals for the construction of the R3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_c1_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick_c1
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck_c1
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_c1_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_c1_p => null()
!
!     Integrals for the construction of the L3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_klic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_iljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ilkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
      type(batching_index) :: batch_i, batch_j, batch_k
      integer  :: i, j, k, i_rel, j_rel, k_rel
      integer  :: i_batch, j_batch, k_batch
      integer  :: req_0, req_1, req_2, req_3, req_1_eri, req_i
      real(dp) :: ddot
      integer :: req_single_batch
!
      call mem%alloc(t2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_c1_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + 2*wf%n_v**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + 3*wf%n_v**3*wf%n_o &
                       + 3*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 3*wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 6*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           'L_R_overlap_triples',      &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(R_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(L_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_ljci_c1, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      endif
!
!     Loop over the batches in i,j,k
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(wf%eri_t1, g_bdci, g_bdci_p, sorting, batch_i)
         call wf%setup_vvvo(wf%eri_c1, g_bdci_c1, g_bdci_c1_p, sorting, batch_i)
!
         call wf%setup_vvov(g_dbic, g_dbic_p, sorting, batch_i, left=.true.)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(wf%eri_t1, g_ljci, g_ljci_p, sorting, batch_j, batch_i)
            call wf%setup_oovo(wf%eri_c1, g_ljci_c1, g_ljci_c1_p, sorting, batch_j, batch_i)
!
            call wf%setup_ooov(g_jlic, g_jlic_p, sorting, batch_j, batch_i)
!
            call wf%setup_ovov(wf%eri_t1, g_ibjc, g_ibjc_p, sorting, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(wf%eri_t1, g_bdcj, g_bdcj_p, sorting, batch_j)
               call wf%setup_vvvo(wf%eri_c1, g_bdcj_c1, g_bdcj_c1_p, sorting, batch_j)
!
               call wf%setup_vvov(g_dbjc, g_dbjc_p, sorting, batch_j, left=.true.)
!
               call wf%setup_oovo(wf%eri_t1, g_licj, g_licj_p, sorting, batch_i, batch_j)
               call wf%setup_oovo(wf%eri_c1, g_licj_c1, g_licj_c1_p, sorting, batch_i, batch_j)
!
               call wf%setup_ooov(g_iljc, g_iljc_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
               call wf%point_vvvo(g_bdcj_c1_p, g_bdci_c1, batch_j%length)
!
               call wf%point_vvvo(g_dbjc_p, g_dbic, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
               call wf%point_vooo(g_licj_c1_p, g_ljci_c1, batch_i%length, batch_j%length)
!
               call wf%point_vooo(g_iljc_p, g_jlic, batch_i%length, batch_j%length)
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%setup_vvvo(wf%eri_t1, g_bdck, g_bdck_p, sorting, batch_k)
                  call wf%setup_vvvo(wf%eri_c1, g_bdck_c1, g_bdck_c1_p, sorting, batch_k)
!
                  call wf%setup_vvov(g_dbkc, g_dbkc_p, sorting, batch_k, left=.true.)
!
                  call wf%setup_oovo(wf%eri_t1, g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(wf%eri_t1, g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
                  call wf%setup_oovo(wf%eri_c1, g_lick_c1, g_lick_c1_p, sorting, batch_i, batch_k)
!
                  call wf%setup_oovo(wf%eri_c1, g_ljck_c1, g_ljck_c1_p, sorting, batch_j, batch_k)
!
                  call wf%setup_oovo(wf%eri_c1, g_lkci_c1, g_lkci_c1_p, sorting, batch_k, batch_i)
!
                  call wf%setup_oovo(wf%eri_c1, g_lkcj_c1, g_lkcj_c1_p, sorting, batch_k, batch_j)
!
!
                  call wf%setup_ooov(g_ilkc, g_ilkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ooov(g_jlkc, g_jlkc_p, sorting, batch_j, batch_k)
                  call wf%setup_ooov(g_klic, g_klic_p, sorting, batch_k, batch_i)
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ovov(wf%eri_t1, g_ibkc, g_ibkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ovov(wf%eri_t1, g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
                  call wf%point_vvvo(g_bdck_c1_p, g_bdci_c1, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbic, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_lick_c1_p, g_ljci_c1, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_c1_p, g_ljci_c1, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_c1_p, g_ljci_c1, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_c1_p, g_ljci_c1, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_ilkc_p, g_jlic, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_jlic, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_kljc_p, g_jlic, batch_k%length, batch_j%length)
!
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
                  call wf%point_vvoo(g_jbkc_p, g_ibjc, batch_j%length, batch_k%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
                  call wf%point_vvvo(g_bdck_c1_p, g_bdcj_c1, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbjc, batch_k%length)
!
                  call wf%setup_oovo(wf%eri_t1, g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
                  call wf%setup_oovo(wf%eri_c1, g_lkcj_c1, g_lkcj_c1_p, sorting, batch_k, batch_j)
!
                  call wf%point_vooo(g_lick_c1_p, g_licj_c1, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_c1_p, g_lkcj_c1, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_c1_p, g_ljci_c1, batch_k%length, batch_i%length)
!
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_ilkc_p, g_iljc, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_kljc, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
!
                  call wf%setup_ovov(wf%eri_t1, g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
!
               endif
!
               do i = batch_i%first, batch_i%get_last()
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%get_last(), i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%get_last(), j)
!
                        k_rel = k - batch_k%first + 1
!
!                       Check for core orbitals (used for excited states):
!                       cvs: i,j,k cannot all correspond to valence orbitals
!                       rm_core: i,j,k may not contain any core orbital
                        if (wf%ijk_amplitudes_are_zero(i, j, k, cvs, rm_core)) cycle
!
!                       Construct R^{abc}_{ijk} for given i, j, k
!                       Using c1-transformed integrals the terms have the same form
!                       as the omega terms (where t_abc = R_abc)
!
                        call wf%construct_V(i, j, k, sorting,             &
                                            R_abc, t2, R2,                &
                                            g_bdci_p(:,:,:,i_rel),        &
                                            g_bdcj_p(:,:,:,j_rel),        &
                                            g_bdck_p(:,:,:,k_rel),        &
                                            g_bdci_c1_p(:,:,:,i_rel),     &
                                            g_bdcj_c1_p(:,:,:,j_rel),     &
                                            g_bdck_c1_p(:,:,:,k_rel),     &
                                            g_ljci_p(:,:,j_rel,i_rel),    &
                                            g_lkci_p(:,:,k_rel,i_rel),    &
                                            g_lkcj_p(:,:,k_rel,j_rel),    &
                                            g_licj_p(:,:,i_rel,j_rel),    &
                                            g_lick_p(:,:,i_rel,k_rel),    &
                                            g_ljck_p(:,:,j_rel,k_rel),    &
                                            g_ljci_c1_p(:,:,j_rel,i_rel), &
                                            g_lkci_c1_p(:,:,k_rel,i_rel), &
                                            g_lkcj_c1_p(:,:,k_rel,j_rel), &
                                            g_licj_c1_p(:,:,i_rel,j_rel), &
                                            g_lick_c1_p(:,:,i_rel,k_rel), &
                                            g_ljck_c1_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, R_abc, omega_R)
!
!                       Construct covariant W_abc for given i,j,k
!                       L_abc is obtained by dividing the linear combination
!                       4W_abc - 2W_bac - 2W_cba - 2W_acb + W_bca + W_cab
!                       by omega - eps^abc_ijk
!
                        call wf%construct_W(i, j, k, sorting, L_abc, L2, &
                                            g_dbic_p(:,:,:,i_rel),      &
                                            g_dbjc_p(:,:,:,j_rel),      &
                                            g_dbkc_p(:,:,:,k_rel),      &
                                            g_jlic_p(:,:,j_rel,i_rel),  &
                                            g_klic_p(:,:,k_rel,i_rel),  &
                                            g_kljc_p(:,:,k_rel,j_rel),  &
                                            g_iljc_p(:,:,i_rel,j_rel),  &
                                            g_ilkc_p(:,:,i_rel,k_rel),  &
                                            g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%outer_product_terms_l3(i, j, k, L1, L2, L_abc,    &
                                                       wf%fock_ia,                &
                                                       g_ibjc_p(:,:,i_rel,j_rel), &
                                                       g_ibkc_p(:,:,i_rel,k_rel), &
                                                       g_jbkc_p(:,:,j_rel,k_rel))
!
                        call construct_contravariant_t3(L_abc, sorting, wf%n_v)
!
                        call wf%divide_by_orbital_differences(i, j, k, L_abc, omega_L)
!
!                       Overlap:
!                          LT_R = sum_{ai >= bj >= ck} L^abc_ijk R^abc_ijk
!                               = 1/6 sum_abcijk L^abc_ijk R^abc_ijk
!
!                       We don't have the factor of 1/6 due to the restrictions
!                       on the loops. Need to account for doubles counting
!                       if 2 indices are the same e.g.:
!                             L^abc_ijk*R^abc_ijk = L^abc_jik*R^abc_jik
!
                        if (i .ne. j .and. j .ne. k) then
                           LT_R = LT_R + ddot(wf%n_v**3, L_abc, 1 , R_abc, 1)
                        else
                           LT_R = LT_R + half*ddot(wf%n_v**3, L_abc, 1 , R_abc, 1)
                        end if
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      call mem%dealloc(R_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(L_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(t2, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci_c1, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck_c1, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      call mem%batch_finalize()
!
   end subroutine L_R_overlap_triples_cc3
!
!
   subroutine estimate_mem_integral_setup_cc3(wf, req0, req1)
!!
!!    Estimate memory integrals setup
!!    Written by Alexander C. Paul, Dec 2020
!!
!!    Estimate maximum memory needed for cc3 integral setup
!!
!!    get_eri_t1_mem returns the memory needed to construct the requested integral
!!    The dimensions sent in specify if an index is batched (1) or of
!!    full dimension (n_o/n_v)
!!    The memory estimate for the first and second pair of indices
!!    is added to the integers req*.
!!
!!    The memory needed to get vvov and vvvo is identical
!!    The memory needed to get oovo and ooov is identical
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(out) :: req0, req1
!
      integer, dimension(2) :: req_vvvo, req_ovov, req_oovo
!
      req_vvvo = wf%eri_t1%get_memory_estimate('vvvo', wf%n_v, wf%n_v, wf%n_v, 1)
      req_ovov = wf%eri_t1%get_memory_estimate('ovov', 1, wf%n_v, 1, wf%n_v)
      req_oovo = wf%eri_t1%get_memory_estimate('oovo', wf%n_o, 1, wf%n_v, 1)
!
      req0 = req_vvvo(1)
      req1 = max(req_vvvo(2), req_ovov(1) + req_ovov(2), req_oovo(1) + req_ovov(2))
!
   end subroutine estimate_mem_integral_setup_cc3
!
!
   subroutine estimate_mem_c1_integral_setup_cc3(wf, req0, req1)
!!
!!    Estimate memory C1 transformed integrals setup
!!    Written by Alexander C. Paul, Dec 2020
!!
!!    Estimate maximum memory needed for cc3 integral setup
!!    for the C1-transformed integrals
!!
!!    get_eri_c1_mem returns the memory needed to construct the requested
!!    c1-transformed integral
!!    The dimensions sent in specify if an index is batched (1) or of
!!    full dimension (n_o/n_v)
!!
!!    6 memory estimates are returned:
!!    1 for each index added to the first 4 req* variables sent in
!!    1 for the first and 1 for the second pair of indices (last 2 req* variables)
!!
!!    NB: The memory requirement is overestimated by the routines.
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(out) :: req0, req1
!
      integer, dimension(2) :: req_vvvo, req_oovo
!
      req_vvvo = wf%eri_c1%get_memory_estimate('vvvo', wf%n_v, wf%n_v, wf%n_v, 1)
      req_oovo = wf%eri_c1%get_memory_estimate('oovo', wf%n_o, 1, wf%n_v, 1)
!
      req0 = req_vvvo(1)
      req1 = max(req_vvvo(2), req_oovo(1) + req_oovo(2))
!
   end subroutine estimate_mem_c1_integral_setup_cc3
!
!
   pure function one_core_index_cc3(wf, i, j, k, check) result(found)
!!
!!    One core index
!!    Written by Alexander C. Paul, Apr 2021
!!
!!    Checks if at least one of i,j,k is a core index.
!!    The logical check is used to skip the check completely.
!!
!!    This is used to determine if a triple loop over i,j,k
!!    can be cycled, eg. due to cvs or removing core orbitals
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: i, j, k
      logical, intent(in) :: check
      logical :: found
!
      found = .false.
!
      if (check) then
         found =    any(wf%core_MOs .eq. i) &
               .or. any(wf%core_MOs .eq. j) &
               .or. any(wf%core_MOs .eq. k)
      end if
!
   end function one_core_index_cc3
!
!
   pure function two_core_indices_cc3(wf, i, j, k, check) result(found)
!!
!!    Two core indices
!!    Written by Alexander C. Paul, Apr 2021
!!
!!    Checks if at least two of i,j,k are core indices.
!!    The logical check is used to skip the check completely.
!!
!!    This is used to determine if a triple loop over i,j,k
!!    can be cycled, eg. to remove core orbitals
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: i, j, k
      logical, intent(in) :: check
      logical :: found
!
      logical :: i_core, j_core, k_core
!
      found = .false.
!
      if (check) then
!
         i_core = any(wf%core_MOs .eq. i)
         j_core = any(wf%core_MOs .eq. j)
         k_core = any(wf%core_MOs .eq. k)
!
         if (i_core .and. j_core) then
            found = .true.
         else if (j_core .and. k_core) then
            found = .true.
         else if (i_core .and. k_core) then
            found = .true.
         end if
!
      end if
!
   end function two_core_indices_cc3
!
!
   pure function three_core_indices_cc3(wf, i, j, k, check) result(found)
!!
!!    Three core indices
!!    Written by Alexander C. Paul, Apr 2021
!!
!!    Checks if i,j and k are core indices.
!!    The logical check is used to skip the check completely.
!!
!!    This is used to determine if a triple loop over i,j,k
!!    can be cycled, eg. to remove core orbitals
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: i, j, k
      logical, intent(in) :: check
      logical :: found
!
      found = .false.
!
      if (check) then
!
         found =     any(wf%core_MOs .eq. i) &
               .and. any(wf%core_MOs .eq. j) &
               .and. any(wf%core_MOs .eq. k)
!
      end if
!
   end function three_core_indices_cc3
!
!
   pure function ijk_amplitudes_are_zero_cc3(wf, i, j, k, cvs, rm_core) result(are_zero)
!!
!!    ijk amplitudes are zero
!!    Written by Alexander C. Paul, Apr 2021
!!
!!    cvs:     returns true if none of i,j,k is a core orbital
!!    rm_core: returns true if at least one of i,j,k is a core orbital
!!
!!    Used to determine if a triple loop over i,j,k shall be cycled in CC3
!!    if core states or states with removed core contribution shall be calculated.
!!    In most of the i,j,k loops all contributions can be skipped
!!    if the conditions below are true.
!!
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: i, j, k
      logical, intent(in) :: cvs, rm_core
      logical :: are_zero
!
      are_zero = .false.
!
      if (i .eq. j .and. j .eq. k) then ! k == j == i
!
         are_zero = .true.
!
      else if (cvs) then
!
         are_zero = .not. (any(wf%core_MOs .eq. i) &
                    .or.   any(wf%core_MOs .eq. j) &
                    .or.   any(wf%core_MOs .eq. k))
!
      else if (rm_core) then
!
         are_zero =   (any(wf%core_MOs .eq. i) &
                  .or. any(wf%core_MOs .eq. j) &
                  .or. any(wf%core_MOs .eq. k))
!
      end if
!
   end function ijk_amplitudes_are_zero_cc3
!
!
   subroutine construct_c1_cholesky_cc3(wf, c1, L, L_c1)
!!
!!    Construct C1 Cholesky
!!
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Constructs the "C1 Cholesky vector"
!!
      use array_utilities, only: zero_array
      use abstract_eri_cholesky_class, only: abstract_eri_cholesky
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:,:,:), allocatable :: L_Jia
!

      call wf%construct_cholesky_c1_oo(L, L_c1, c1)
      call wf%construct_cholesky_c1_vo(L, L_c1, c1)
      call wf%construct_cholesky_c1_vv(L, L_c1, c1)
!
      call mem%alloc(L_Jia, L%n_J, wf%n_o, wf%n_v)
      call zero_array(L_Jia, L%n_J*wf%n_o*wf%n_v)
      call L_c1%set(L_Jia, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
      call mem%dealloc(L_Jia, L%n_J, wf%n_o, wf%n_v)
!
      call wf%L_c1%notify_observers()
!
   end subroutine construct_c1_cholesky_cc3
!
!
   subroutine construct_cholesky_c1_oo_cc3(wf, L, L_c1, c1)
!!
!!    Construct Cholesky oo C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_J_ij_C1 = sum_b L_J_ib_T1 c_bj
!!
!!    and returns the result in L_J_ij_C1
!!
      use abstract_eri_cholesky_class, only: abstract_eri_cholesky
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:), allocatable :: L_J_oo
      real(dp), dimension(:), allocatable :: L_J_ov
!
      type(batching_index) :: batch_o
!
      integer :: o_batch
!
      batch_o = batching_index(wf%n_o)
      call mem%batch_setup(batch_o, 0, L%n_J*(wf%n_v + wf%n_o), tag='Cholesky c1 oo')
!
      call mem%alloc(L_J_oo, L%n_J*batch_o%max_length*wf%n_o)
      call mem%alloc(L_J_ov, L%n_J*batch_o%max_length*wf%n_v)
!
      do o_batch = 1,batch_o%num_batches
!
         call batch_o%determine_limits(o_batch)
!
         call L%get(L_J_ov, batch_o%first, batch_o%get_last(), wf%n_o + 1, wf%n_mo)
!
         call dgemm('N', 'N',                               &
                    L%n_J*batch_o%length, wf%n_o, wf%n_v,   &
                    one,                                    &
                    L_J_ov, L%n_J*batch_o%length,           &
                    c1, wf%n_v,                             &
                    zero,                                   &
                    L_J_oo, L%n_J*batch_o%length)
!
         call L_c1%set(L_J_oo, batch_o%first, batch_o%get_last(), 1, wf%n_o)
!
      enddo
!
      call mem%dealloc(L_J_oo, L%n_J*batch_o%max_length*wf%n_o)
      call mem%dealloc(L_J_ov, L%n_J*batch_o%max_length*wf%n_v)
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_c1_oo_cc3
!
!
   subroutine construct_cholesky_c1_vo_cc3(wf, L, L_c1, c1)
!!
!!    Construct Cholesky vo C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_J_ai_C1 = sum_b L_J_ab_T1 c_bi - sum_j c_aj L_J_ji_T1
!!
!!    and returns the result in L_J_ai_c1
!!
      use abstract_eri_cholesky_class, only: abstract_eri_cholesky
      use array_utilities, only: zero_array
      use reordering, only: add_132_to_123, sort_123_to_132
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ij
      real(dp), dimension(:,:,:), allocatable :: L_J_ji
      real(dp), dimension(:,:,:), allocatable :: L_J_ai
      real(dp), dimension(:,:,:), allocatable :: L_J_ia
      real(dp), dimension(:,:,:), allocatable :: L_J_ab
!
      type(batching_index) :: batch_a, batch_j
!
      integer :: a_batch, j_batch, req
!
      batch_j = batching_index(wf%n_o)
      batch_a = batching_index(wf%n_v)
      req = L%n_J*(wf%n_o + max(wf%n_v,wf%n_o))
!
      call mem%batch_setup(batch_j, batch_a, 0, req, req, 0, tag='Cholesky c1 vo')
!
      do a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(a_batch)
!
         call mem%alloc(L_J_ab, L%n_J, batch_a%length, wf%n_v)
!
         call L%get(L_J_ab,                        &
                    wf%n_o + batch_a%first,        &
                    wf%n_o + batch_a%get_last(),   &
                    wf%n_o + 1,                    &
                    wf%n_mo)
!
         call mem%alloc(L_J_ai, L%n_J, batch_a%length, wf%n_o)
!
         call dgemm('N', 'N',             &
                    L%n_J*batch_a%length, &
                    wf%n_o,               &
                    wf%n_v,               &
                    one,                  &
                    L_J_ab,               &
                    L%n_J*batch_a%length, &
                    c1,                   &
                    wf%n_v,               &
                    zero,                 &
                    L_J_ai,               &
                    L%n_J*batch_a%length)
!
         call mem%dealloc(L_J_ab, L%n_J, batch_a%length, wf%n_v)
         call mem%alloc(L_J_ia, L%n_J, wf%n_o, batch_a%length)
!
         do j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(j_batch)
!
            call mem%alloc(L_J_ji, L%n_J, batch_j%length, wf%n_o)
            call L%get(L_J_ji, batch_j%first, batch_j%get_last(), 1, wf%n_o)

            call mem%alloc(L_J_ij, L%n_J, wf%n_o, batch_j%length)
            call sort_123_to_132(L_J_ji, L_J_ij, L%n_J, batch_j%length, wf%n_o)
            call mem%dealloc(L_J_ji, L%n_J, batch_j%length, wf%n_o)
!
            call dgemm('N', 'T',       &
                       L%n_J*wf%n_o,   &
                       batch_a%length, &
                       batch_j%length, &
                       -one,           &
                       L_J_ij,         &
                       L%n_J*wf%n_o,   &
                       c1(batch_a%first, batch_j%first), &
                       wf%n_v,         &
                       zero,           &
                       L_J_ia,         &
                       L%n_J*wf%n_o)
!
            call mem%dealloc(L_J_ij, L%n_J, wf%n_o, batch_j%length)
            call add_132_to_123(one, L_J_ia, L_J_ai, L%n_J, batch_a%length, wf%n_o)
!
         enddo
!
         call L_c1%set(L_J_ai, wf%n_o+batch_a%first, wf%n_o+batch_a%get_last(), 1, wf%n_o)
         call mem%dealloc(L_J_ai, L%n_J, batch_a%length, wf%n_o)
         call mem%dealloc(L_J_ia, L%n_J, wf%n_o, batch_a%length)
!
      enddo
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_c1_vo_cc3
!
!
   subroutine construct_cholesky_c1_vv_cc3(wf, L, L_c1, c1)
!!
!!    Construct Cholesky vv C1
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_ab_J_c1= - sum_i c_ai L_ib_J_T1 ,
!!
!!    and returns the result in L_J_ab_c1
!!
      use abstract_eri_cholesky_class, only: abstract_eri_cholesky
      use reordering, only: sort_123_to_132
!
      implicit none

      class(cc3), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ab, L_J_ba, L_J_jb, L_J_bj
!
      type(batching_index) :: batch_b
!
      integer :: b_batch
!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, 0, 2*L%n_J*max(wf%n_v, wf%n_o), tag='Cholesky c1 vv')
!
      do b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(b_batch)
!
         call mem%alloc(L_J_jb, L%n_J, wf%n_o, batch_b%length)
         call L%get(L_J_jb, 1, wf%n_o, wf%n_o + batch_b%first, wf%n_o + batch_b%get_last())
!
         call mem%alloc(L_J_bj, L%n_J, batch_b%length, wf%n_o)
         call sort_123_to_132(L_J_jb, L_J_bj, L%n_J, wf%n_o, batch_b%length)
         call mem%dealloc(L_J_jb, L%n_J, wf%n_o, batch_b%length)
!
         call mem%alloc(L_J_ba, L%n_J, batch_b%length, wf%n_v)
         call dgemm('N', 'T',             &
                    L%n_J*batch_b%length, &
                    wf%n_v,               &
                    wf%n_o,               &
                    -one,                 &
                    L_J_bj,               &
                    L%n_J*batch_b%length, &
                    c1,                   &
                    wf%n_v,               &
                    zero,                 &
                    L_J_ba,               &
                    L%n_J*batch_b%length)
!
         call mem%dealloc(L_J_bj, L%n_J, batch_b%length, wf%n_o)
!
         call mem%alloc(L_J_ab, L%n_J, wf%n_v, batch_b%length)
         call sort_123_to_132(L_J_ba, L_J_ab, L%n_J, batch_b%length, wf%n_v)
         call mem%dealloc(L_J_ba, L%n_J, batch_b%length, wf%n_v)
!
         call L_c1%set(L_J_ab, wf%n_o+1, wf%n_mo, wf%n_o + batch_b%first, wf%n_o + batch_b%get_last())
         call mem%dealloc(L_J_ab, L%n_J, wf%n_v, batch_b%length)
!
      enddo
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_c1_vv_cc3
!
!
end module cc3_class
