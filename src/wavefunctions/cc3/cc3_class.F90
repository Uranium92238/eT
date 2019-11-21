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
module cc3_class
!
!!
!!    Coupled cluster CC3 class module
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2018-2019
!!
!
   use kinds
   use ccsd_class
   use direct_file_class, only : direct_file
   use sequential_file_class, only : sequential_file
   use io_utilities, only : single_record_reader, compound_record_reader 
   use io_utilities, only : single_record_writer, compound_record_writer
   use array_utilities, only: entrywise_product
!
   implicit none
!
   type, extends(ccsd) :: cc3
!
!     Ground state integral files
!
      type(direct_file) :: g_bdck_t
      type(direct_file) :: g_ljck_t
      type(direct_file) :: g_dbkc_t
      type(direct_file) :: g_jlkc_t
      type(direct_file) :: L_jbkc_t
!
!     Right Jacobian integral files
!
      type(direct_file) :: g_bdck_c1
      type(direct_file) :: g_ljck_c1
!
!     Left Jacobian integral files
!
      type(direct_file) :: g_becd_t
      type(direct_file) :: g_mjlk_t
      type(direct_file) :: g_ckld_t
      type(direct_file) :: g_cdlk_t
!
!     Jacobian intermediates files
!
      type(direct_file) :: g_lbkc_t
      type(direct_file) :: X_abdi
      type(direct_file) :: X_abid
      type(direct_file) :: Y_bcek
      type(direct_file) :: X_ajil
!
!     Files for batching of the virtual indices
!
      type(direct_file) :: g_bdck_t_v
      type(direct_file) :: g_ljck_t_v
      type(direct_file) :: g_dbkc_t_v
      type(direct_file) :: g_jlkc_t_v
      type(direct_file) :: L_jbkc_t_v
!
      type(direct_file) :: g_bdck_c1_v
      type(direct_file) :: g_ljck_c1_v
!
!     Density intermediates files
!
      type(direct_file) :: Y_clik_tbar
      type(direct_file) :: Z_bcjk
!
      real(dp), dimension(:,:), allocatable :: GS_cc3_density_oo
      real(dp), dimension(:,:), allocatable :: GS_cc3_density_vv
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                   => cleanup_cc3
      procedure :: delete_intermediate_files => delete_intermediate_files_cc3
!
!     Routines related to omega
!
      procedure :: construct_omega     => construct_omega_cc3
!
      procedure :: omega_cc3_a         => omega_cc3_a_cc3
      procedure :: omega_cc3_integrals => omega_cc3_integrals_cc3
      procedure :: omega_cc3_W_calc    => omega_cc3_W_calc_cc3
      procedure :: omega_cc3_eps       => omega_cc3_eps_cc3
      procedure :: omega_cc3_a_n6      => omega_cc3_a_n6_cc3
      procedure :: omega_cc3_a_n7      => omega_cc3_a_n7_cc3
!
!     Routines used for prepare, both left and right
!
      procedure :: prepare_for_jacobian                 => prepare_for_jacobian_cc3
      procedure :: prepare_for_jacobian_transpose       => prepare_for_jacobian_transpose_cc3
      procedure :: prepare_cc3_g_lbkc_t_file            => prepare_cc3_g_lbkc_t_file_cc3
      procedure :: prepare_cc3_jacobian_intermediates   => prepare_cc3_jacobian_intermediates_cc3
      procedure :: construct_x_intermediates            => construct_x_intermediates_cc3
      procedure :: sort_x_to_abid_and_write             => sort_x_to_abid_and_write_cc3
!
      procedure :: prepare_cc3_jacobian_trans_integrals => prepare_cc3_jacobian_trans_integrals_cc3
!
!     Routines for CVS
      procedure :: get_cvs_projector                   => get_cvs_projector_cc3
      procedure :: get_triples_cvs_projector_abc_batch => get_triples_cvs_projector_abc_batch_cc3
!
!     Routines related to the jacobian
!
      procedure :: construct_Jacobian_transform      => construct_Jacobian_transform_cc3
!
!     Right hand side transformation
!
      procedure :: effective_jacobian_transformation => effective_jacobian_transformation_cc3
!
      procedure :: jacobian_cc3_t3_a2                => jacobian_cc3_t3_a2_cc3
      procedure :: jacobian_cc3_t3_b2                => jacobian_cc3_t3_b2_cc3
      procedure :: construct_c1_fock                 => construct_c1_fock_cc3
      procedure :: jacobian_cc3_b2_fock              => jacobian_cc3_b2_fock_cc3
      procedure :: jacobian_cc3_c3_a                 => jacobian_cc3_c3_a_cc3
      procedure :: construct_c1_integrals            => construct_c1_integrals_cc3
!
!     Routines related to the transpose of the jacobian
!
      procedure :: effective_jacobian_transpose_transformation  &
                                    => effective_jacobian_transpose_transformation_cc3
!
      procedure :: jacobian_transpose_cc3_t3_a1       => jacobian_transpose_cc3_t3_a1_cc3
      procedure :: jacobian_transpose_cc3_t3_b1       => jacobian_transpose_cc3_t3_b1_cc3
      procedure :: construct_x_ai_intermediate        => construct_x_ai_intermediate_cc3
      procedure :: jacobian_transpose_cc3_c3_a        => jacobian_transpose_cc3_c3_a_cc3
      procedure :: jacobian_transpose_cc3_c3_calc     => jacobian_transpose_cc3_c3_calc_cc3
      procedure :: jacobian_transpose_cc3_a_n7        => jacobian_transpose_cc3_a_n7_cc3
      procedure :: construct_y_intermediates          => construct_y_intermediates_cc3
      procedure :: jacobian_transpose_cc3_c3_a1_y_o   => jacobian_transpose_cc3_c3_a1_y_o_cc3
      procedure :: jacobian_transpose_cc3_c3_b1_y_v   => jacobian_transpose_cc3_c3_b1_y_v_cc3
!
!     Routines related to the multipliers
!
      procedure :: prepare_for_multiplier_equation    => prepare_for_multiplier_equation_cc3
      procedure :: construct_multiplier_equation      => construct_multiplier_equation_cc3
      procedure :: save_tbar_intermediates            => save_tbar_intermediates_cc3
!
!     Routines to construct triples amplitudes in batches of a,b,c
!
      procedure :: prepare_cc3_integrals_t3_abc_batch => prepare_cc3_integrals_t3_abc_batch_cc3
      procedure :: prepare_cc3_integrals_R3_abc_batch => prepare_cc3_integrals_R3_abc_batch_cc3
      procedure :: prepare_cc3_integrals_L3_abc_batch => prepare_cc3_integrals_L3_abc_batch_cc3
      procedure :: omega_cc3_W_calc_abc_batch         => omega_cc3_W_calc_abc_batch_cc3
      procedure :: omega_cc3_eps_abc_batch            => omega_cc3_eps_abc_batch_cc3
      procedure :: jacobian_transpose_cc3_c3_calc_abc_batch &
                                                      => jacobian_transpose_cc3_c3_calc_abc_batch_cc3
!
!     Routines related to density matrices
!
      procedure :: initialize_gs_density              => initialize_gs_density_cc3
      procedure :: destruct_gs_density                => destruct_gs_density_cc3
!
      procedure :: prepare_for_density                => prepare_for_density_cc3
      procedure :: construct_gs_density               => construct_gs_density_cc3
      procedure :: construct_left_transition_density  => construct_left_transition_density_cc3
      procedure :: construct_right_transition_density => construct_right_transition_density_cc3
!
      procedure :: density_cc3_mu_ref_abc             => density_cc3_mu_ref_abc_cc3
      procedure :: density_cc3_mu_ref_oo              => density_cc3_mu_ref_oo_cc3
      procedure :: density_cc3_mu_ref_ijk             => density_cc3_mu_ref_ijk_cc3
      procedure :: density_cc3_mu_ref_vv              => density_cc3_mu_ref_vv_cc3
      procedure :: construct_y_intermediate_vo3       => construct_y_intermediate_vo3_cc3
!
      procedure :: density_cc3_mu_nu_ov               => density_cc3_mu_nu_ov_cc3
      procedure :: density_cc3_mu_nu_oo_ov_vv         => density_cc3_mu_nu_oo_ov_vv_cc3
      procedure :: density_cc3_mu3_nu2_ov             => density_cc3_mu3_nu2_ov_cc3
      procedure :: density_cc3_mu_nu_ijk              => density_cc3_mu_nu_ijk_cc3
      procedure :: density_cc3_mu_nu_abc              => density_cc3_mu_nu_abc_cc3
      procedure :: density_cc3_mu_nu_ov_Z_term        => density_cc3_mu_nu_ov_Z_term_cc3
      procedure :: density_cc3_mu_nu_vo               => density_cc3_mu_nu_vo_cc3
      procedure :: construct_Z_intermediate           => construct_Z_intermediate_cc3
!
!     Routines related to the biorthonormalization
!
      procedure :: L_R_overlap                            => L_R_overlap_cc3
      procedure :: L_R_overlap_triples                    => L_R_overlap_triples_cc3
!
      procedure :: scale_triples_biorthonormal_factor     => scale_triples_biorthonormal_factor_cc3
      procedure :: scale_triples_biorthonormal_factor_abc => scale_triples_biorthonormal_factor_abc_cc3
!
      procedure :: biorthonormalize_L_and_R               => biorthonormalize_L_and_R_cc3
!
!     Routines for constructing the R-tdm in a noddy way
!     only for debugging purposes
!
      procedure :: right_tdm_debug       => right_tdm_debug_cc3
      procedure :: left_tdm_debug        => left_tdm_debug_cc3
!
      procedure :: debug_left_oo         => debug_left_oo_cc3
      procedure :: debug_left_ov_N7      => debug_left_ov_N7_cc3
      procedure :: debug_left_ov_N6      => debug_left_ov_N6_cc3
      procedure :: debug_left_vv         => debug_left_vv_cc3
      procedure :: debug_right_ov_t3     => debug_right_ov_t3_cc3
      procedure :: debug_right_ov_Y_term => debug_right_ov_Y_term_cc3
      procedure :: debug_right_oo        => debug_right_oo_cc3
      procedure :: debug_right_vv        => debug_right_vv_cc3
      procedure :: debug_right_ov_R3     => debug_right_ov_R3_cc3
      procedure :: debug_right_vo        => debug_right_vo_cc3
!
!     Construction of the full amplitudes
      procedure :: construct_full_R3     => construct_full_R3_cc3
      procedure :: construct_full_t3     => construct_full_t3_cc3
      procedure :: construct_full_tbar3  => construct_full_tbar3_cc3
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
      include "cc3_batching_abc_interface.F90"
      include "zop_cc3_interface.F90"
      include "fop_cc3_interface.F90"
!
      include "debug_transition_density_cc3_interface.F90"
!
   end interface
!
!
   interface cc3
!
      procedure :: new_cc3 
!
   end interface cc3
!
!
contains
!
!
   function new_cc3(system) result(wf)
!!
!!    New CC3
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none
!
      type(cc3) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%name_ = 'cc3'
!
      call wf%general_cc_preparations(system)
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
      wf%n_t2 = (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
      wf%need_g_abcd     = .true.
!
      call wf%initialize_fock()
!
   end function new_cc3
!
!
   subroutine cleanup_cc3(wf)
!!
!!    Cleanup
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      call output%printf('- Cleaning up (a0) wavefunction', chars=[trim(wf%name_)], fs='(/t3,a)')
!
      call wf%destruct_amplitudes
      call wf%destruct_multipliers
      call wf%destruct_right_excitation_energies()
      call wf%destruct_left_excitation_energies()
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
      type(direct_file) :: Y_bcek_tbar
!
!     Delete files for GS
      if (wf%g_bdck_t%exists()) then
!
         call wf%g_bdck_t%delete_
         call wf%g_ljck_t%delete_
         call wf%g_dbkc_t%delete_
         call wf%g_jlkc_t%delete_
         call wf%L_jbkc_t%delete_
!
      end if
!
!     Delete Intermediate files for the jacobian transformations
      if (wf%X_abid%exists()) then
!
         call wf%X_abid%delete_
         call wf%X_ajil%delete_
!
      end if
!
!     Delete additional files for jacobian transformation
      if (wf%g_bdck_c1%exists()) then
!
         call wf%g_bdck_c1%delete_
         call wf%g_ljck_c1%delete_
!
      end if
!
!     Delete additional files for jacobian transpose transformation
      if (wf%g_becd_t%exists()) then
!
         call wf%g_becd_t%delete_
         call wf%g_mjlk_t%delete_
         call wf%g_ckld_t%delete_
         call wf%g_cdlk_t%delete_
         call wf%Y_bcek%delete_
!
      end if
!
!     Delete integral files needed for t3 and tbar3 in batches of a,b,c
      if (wf%g_bdck_t_v%exists()) then
!
         call wf%g_bdck_t_v%delete_
         call wf%g_ljck_t_v%delete_
         call wf%g_dbkc_t_v%delete_
         call wf%g_jlkc_t_v%delete_
         call wf%L_jbkc_t_v%delete_
!
      end if
!
!     Delete additional integral files needed R3 in batches of a,b,c
      if (wf%g_bdck_c1_v%exists()) then
!
         call wf%g_bdck_c1_v%delete_
         call wf%g_ljck_c1_v%delete_
!
      end if
!
!     Delete intermediates needed for the transition densities
      if(wf%Y_clik_tbar%exists()) then
!
         call wf%Y_clik_tbar%delete_
         call wf%Z_bcjk%delete_
!
         Y_bcek_tbar = direct_file('Y_bcek_tbar', wf%n_v**3)
         call Y_bcek_tbar%open_('write')
         call Y_bcek_tbar%delete_
!
      end if
!
   end subroutine delete_intermediate_files_cc3
!
!
   subroutine construct_Jacobian_transform_cc3(wf, r_or_l, X, w)
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
      class(cc3), intent(in) :: wf
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
      else
!
         call output%error_msg('w is missing in construct_Jacobian_transform for lowmem_cc2')
!
      endif
!
!     Compute the transformed matrix
      if (r_or_l .eq. "right") then
!
         call wf%effective_jacobian_transformation(w, X) ! X <- AX
!
      else if (r_or_l .eq. "left") then
!
         call wf%effective_jacobian_transpose_transformation(w, X, wf%cvs) ! X <- A^TX
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
   subroutine save_tbar_intermediates_cc3(wf)
!!
!!    Save tbar intermediates
!!    Written by Alexander C. Paul, August 2019
!!
!!    jacobian_transpose_transformation is used 
!!    for multipliers and left excited states
!!    Copy the intermediate Y_bcek after solving
!!    for the multipliers to reuse it 
!!    in the right transition densities
!!
      implicit none
!
      class(cc3) :: wf
!
      type(direct_file) :: Y_bcek_tbar
!
!     Copy the intermediate Y_bcek constructed as follows:
!     Y_bcek = sum_aij tbar^abc_ijk * t^ae_ij
!
!     Later used in the right transition density matrix
!
      Y_bcek_tbar = direct_file('Y_bcek_tbar', wf%n_v**3)
!
!     Delete if the file already exists 
!     e.g. when restarting a crashed CC3 calculation
!
      if(Y_bcek_tbar%exists()) then
         call Y_bcek_tbar%delete_
      end if
!
      call wf%Y_bcek%copy('Y_bcek_tbar')
!
   end subroutine save_tbar_intermediates_cc3
!
!
   real(dp) function L_R_overlap_cc3(wf, L, left_state, R, right_state)
!!
!!    Left right overlap
!!    Written by Alexander C. Paul, Aug 2019
!!
!!    Computes L^T * R for full space L and R (singles, doubles, triples)
!!
      class(cc3), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: L
      integer, intent(in) :: left_state
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in) :: R
      integer, intent(in) :: right_state
!
      real(dp), dimension(:,:), allocatable :: L_ai, R_ai
      real(dp), dimension(:,:,:,:), allocatable :: L_abij, R_abij
!
      real(dp) :: ddot
!
      integer :: a, i
!
      L_R_overlap_cc3 = ddot(wf%n_es_amplitudes, L, 1, R, 1)
!
!     :: Triples contribution to the overlap ::
!
!     Allocate unpacked doubles Vectors of the left and right vector
!
      call mem%alloc(L_ai, wf%n_v, wf%n_o)
      call mem%alloc(R_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, L, 1, L_ai, 1)
      call dcopy(wf%n_t1, R, 1, R_ai, 1)
!
!     need to reconstruct the c1-transformed integrals for the respective right state
      if (wf%n_singlet_states .gt. 1) then
         call wf%construct_c1_integrals(R_ai)
      end if
      call mem%dealloc(R_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(R_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call squareup_and_sort_1234_to_1324(L(wf%n_t1 + 1 : wf%n_es_amplitudes), L_abij, &
                                          wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup_and_sort_1234_to_1324(R(wf%n_t1 + 1 : wf%n_es_amplitudes), R_abij, &
                                          wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Scale the right doubles vector by 1 + delta_ai,bj
!
!$omp parallel do schedule(static) private(a,i) collapse(2)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            R_abij(a,a,i,i) = two*R_abij(a,a,i,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%L_R_overlap_triples(wf%left_excitation_energies(left_state),  &
                                                   wf%right_excitation_energies(right_state),&
                                                   L_ai, L_abij, R_abij, L_R_overlap_cc3)
!
      call mem%dealloc(L_ai, wf%n_v, wf%n_o)
      call mem%dealloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(R_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call output%printf('Overlap of (i0). left and (i0). right state: (f15.10)',   &
                          fs='(/t6,a)', ints=[left_state, right_state], reals=[L_R_overlap_cc3], pl='debug')
!
   end function L_R_overlap_cc3
!
!
   subroutine L_R_overlap_triples_cc3(wf, omega_L, omega_R, L_ai, &
                                      L_abij, R_abij, LT_R)
!!
!!    Left right overlap triples contribution
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
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega_L
      real(dp), intent(in) :: omega_R
!
      real(dp), intent(inout) :: LT_R
!
!     Unpacked Left vector
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: L_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: L_abij
!
!     Unpacked R2-amplitudes
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
!     Unpacked t2-amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Arrays for triples amplitudes and multipliers
      real(dp), dimension(:,:,:), allocatable :: R_abc
      real(dp), dimension(:,:,:), allocatable :: L_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
      real(dp), dimension(:,:,:), allocatable :: v_abc
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck ! Transpose of the fock matrix ov block
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
!     C1 transformed integrals
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3
!
!     for CVS
      integer :: i_cvs
      logical :: ijk_core
!
      real(dp) :: ddot
!
!     Set up arrays for amplitudes
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Setup and Batching loops
!
      req_0 = 4*(wf%n_v)**3 + (wf%n_v)*(wf%n_o)
      req_1 = 3*(wf%n_v)**3
      req_2 = 3*(wf%n_o)*(wf%n_v) + (wf%n_v)**2
      req_3 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, buffer_size = zero)
!
!     Allocate integral arrays
!
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      end if
!
!     Arrays for the triples amplitudes and intermediates
      call mem%alloc(R_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(L_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Fock matrix subblock: Resorting for easier contractions later
      call mem%alloc(F_ov_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_ck, wf%n_o, wf%n_v)
!
!     Remaining integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
!        Ordered such that batching indices are at the end
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_dbkc_t%open_('read')
      call wf%g_jlkc_t%open_('read')
      call wf%L_jbkc_t%open_('read')
!
      call wf%g_bdck_c1%open_('read')
      call wf%g_ljck_c1%open_('read')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci, wf%g_dbkc_t, g_dbic, &
                                    wf%g_bdck_c1, g_bdci_c1)
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
         g_bdci_c1_p => g_bdci_c1
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci, wf%g_jlkc_t,  &
                                       g_jlic, wf%g_ljck_c1, g_ljci_c1)
            call compound_record_reader(batch_i, batch_j, wf%L_jbkc_t, L_ibjc)
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_ibjc_p => L_ibjc
!
            g_ljci_c1_p => g_ljci_c1
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj, wf%g_dbkc_t, g_dbjc, &
                                          wf%g_bdck_c1, g_bdcj_c1)
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
               g_bdcj_c1_p => g_bdcj_c1
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj, wf%g_jlkc_t,  &
                                          g_iljc, wf%g_ljck_c1, g_licj_c1)
               g_licj_p => g_licj
               g_iljc_p => g_iljc
!
               g_licj_c1_p => g_licj_c1
!
            else
!
               g_bdcj_p => g_bdci
               g_dbjc_p => g_dbic
!
               g_licj_p => g_ljci
               g_iljc_p => g_jlic
!
               g_bdcj_c1_p => g_bdci_c1
!
               g_licj_c1_p => g_ljci_c1
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck, wf%g_dbkc_t, g_dbkc, &
                                             wf%g_bdck_c1, g_bdck_c1)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
!
                  g_bdck_c1_p => g_bdck_c1
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, wf%g_jlkc_t,  &
                                             g_klic, wf%g_ljck_c1, g_lkci_c1)
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
!
                  g_lkci_c1_p => g_lkci_c1
!
                  call compound_record_reader(batch_i, batch_k, wf%g_ljck_t, g_lick, wf%g_jlkc_t,  &
                                             g_ilkc, wf%L_jbkc_t, L_ibkc, wf%g_ljck_c1, g_lick_c1)
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  g_lick_c1_p => g_lick_c1
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_jlkc_t,  &
                                             g_kljc, wf%g_ljck_c1, g_lkcj_c1)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
!
                  g_lkcj_c1_p => g_lkcj_c1
!
                  call compound_record_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck, wf%g_jlkc_t,  &
                                             g_jlkc, wf%L_jbkc_t, L_jbkc, wf%g_ljck_c1, g_ljck_c1)
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
!
                  g_ljck_c1_p => g_ljck_c1
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
!
                  g_bdck_c1_p => g_bdci_c1
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
!
                  g_lkci_c1_p => g_ljci_c1
!
                  g_lick_p => g_ljci
                  g_ilkc_p => g_jlic
                  L_ibkc_p => L_ibjc
!
                  g_lick_c1_p => g_ljci_c1
!
                  g_lkcj_p => g_ljci
                  g_kljc_p => g_jlic
!
                  g_lkcj_c1_p => g_ljci_c1
!
                  g_ljck_p => g_ljci
                  g_jlkc_p => g_jlic
                  L_jbkc_p => L_ibjc
!
                  g_ljck_c1_p => g_ljci_c1
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
!
                  g_bdck_c1_p => g_bdcj_c1
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
!
                  g_lkci_c1_p => g_ljci_c1
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
!
                  g_lick_c1_p => g_licj_c1
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, wf%g_jlkc_t,  &
                                             g_kljc, wf%L_jbkc_t, L_jbkc, wf%g_ljck_c1, g_lkcj_c1)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
!
                  g_lkcj_c1_p => g_lkcj_c1
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_jbkc
!
                  g_ljck_c1_p => g_lkcj_c1
!
!
               endif
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (k .eq. i) then ! k == j == i
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Check if at least one index i,j,k is a core orbital
                        if(wf%cvs) then
!
                           ijk_core = .false.
!
                           do i_cvs = 1, wf%n_core_MOs
!
                              if(     i .eq. wf%core_MOs(i_cvs)   &
                                 .or. j .eq. wf%core_MOs(i_cvs)   &
                                 .or. k .eq. wf%core_MOs(i_cvs))  then
!
                                 ijk_core = .true.
!
                              end if
!
                           end do
!
!                          Cycle if i,j,k are not core orbitals
                           if (.not. ijk_core) cycle
!
                        end if
!
!                       Construct R^{abc}_{ijk} for given i, j, k
!                       Using c1-transformed integrals the terms have the same form 
!                       as the omega terms (where t_abc = R_abc)
!
!                       Therefore the contributions to the R3-amplitudes can be computed 
!                       using the same routine once for t1-transformed and once for 
!                       c1-transformed integrals
!
                        call wf%omega_cc3_W_calc(i, j, k, R_abc, u_abc, R_abij, &
                                                g_bdci_p(:,:,:,i_rel),          &
                                                g_bdcj_p(:,:,:,j_rel),          &
                                                g_bdck_p(:,:,:,k_rel),          &
                                                g_ljci_p(:,:,j_rel,i_rel),      &
                                                g_lkci_p(:,:,k_rel,i_rel),      &
                                                g_lkcj_p(:,:,k_rel,j_rel),      &
                                                g_licj_p(:,:,i_rel,j_rel),      &
                                                g_lick_p(:,:,i_rel,k_rel),      &
                                                g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_W_calc(i, j, k, R_abc, u_abc, t_abij, &
                                                g_bdci_c1_p(:,:,:,i_rel),       &
                                                g_bdcj_c1_p(:,:,:,j_rel),       &
                                                g_bdck_c1_p(:,:,:,k_rel),       &
                                                g_ljci_c1_p(:,:,j_rel,i_rel),   &
                                                g_lkci_c1_p(:,:,k_rel,i_rel),   &
                                                g_lkcj_c1_p(:,:,k_rel,j_rel),   &
                                                g_licj_c1_p(:,:,i_rel,j_rel),   &
                                                g_lick_c1_p(:,:,i_rel,k_rel),   &
                                                g_ljck_c1_p(:,:,j_rel,k_rel),   &
                                                .true.) ! Do not overwrite R_abc
!
                        call wf%omega_cc3_eps(i, j, k, R_abc, omega_R)
!
                        call wf%scale_triples_biorthonormal_factor(i, j, k, R_abc)
!
!                       c3_calc does not zero out the array
                        call zero_array(L_abc, wf%n_v**3)
!
!                       construct L3 for fixed i,j,k
                        call wf%jacobian_transpose_cc3_c3_calc(i, j ,k, L_ai, L_abij,     &
                                                               L_abc, u_abc, v_abc,       &
                                                               F_ov_ck,                   &
                                                               L_ibjc_p(:,:,i_rel,j_rel), &
                                                               L_ibkc_p(:,:,i_rel,k_rel), &
                                                               L_jbkc_p(:,:,j_rel,k_rel), &
                                                               g_dbic_p(:,:,:,i_rel),     &
                                                               g_dbjc_p(:,:,:,j_rel),     &
                                                               g_dbkc_p(:,:,:,k_rel),     &
                                                               g_jlic_p(:,:,j_rel,i_rel), &
                                                               g_klic_p(:,:,k_rel,i_rel), &
                                                               g_kljc_p(:,:,k_rel,j_rel), &
                                                               g_iljc_p(:,:,i_rel,j_rel), &
                                                               g_ilkc_p(:,:,i_rel,k_rel), &
                                                               g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_eps(i, j, k, L_abc, omega_L)
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
      call wf%g_bdck_t%close_()
      call wf%g_ljck_t%close_()
      call wf%g_dbkc_t%close_()
      call wf%g_jlkc_t%close_()
      call wf%L_jbkc_t%close_()
!
      call wf%g_bdck_c1%close_()
      call wf%g_ljck_c1%close_()
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck_c1, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck_c1, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Deallocate amplitudes arrays and Fock matrix
!
      call mem%dealloc(R_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(L_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(F_ov_ck, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine L_R_overlap_triples_cc3
!
!
   subroutine scale_triples_biorthonormal_factor_cc3(wf, i, j, k, R_abc)
!!
!!    Scale triples amplitudes by biorthonormal factor
!!    Written by Alexander C. Paul, August 2019
!!
!!    Removing the restrictions on the sum over the triples excitations
!!    introduces a factor similar to (1+delta_ai,bj) in the doubles.
!!    Usually this factor cancles but it is e.g. needed in EOM-CC3.
!!
!!       sum_{ai >= bj >= ck} = 1/6 sum_aibjck (1 + delta_ai,bj 
!!                                                + delta_ai,ck
!!                                                + delta_bj,ck
!!                                                + 2 delta_ai,bj delta_ai,ck)
!!
!!    This routine scales an array of amplitudes (for single i,j,k)
!!    with the expression in paranthesis
!!    NB:
!!    k can only be equal to i if i == j == k, for which R_abc = 0
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: R_abc
!
      integer :: a, b
!
      if (i .eq. j .and. i .ne. k) then ! i == j, i != k
!
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               R_abc(a,a,b) = two*R_abc(a,a,b)
!
            enddo
         enddo
!
      else if (j .eq. k .and. i .ne. k) then ! i != j, j == k
!
         do b = 1, wf%n_v
            do a = 1, wf%n_v
!
               R_abc(a,b,b) = two*R_abc(a,b,b)
!
            enddo 
         enddo
!
      end if
!
   end subroutine scale_triples_biorthonormal_factor_cc3
!
!
   subroutine scale_triples_biorthonormal_factor_abc_cc3(wf, a, b, c, R_ijk)
!!
!!    Scale triples by biorthonormal factor (batching in the virtuals)
!!    Written by Alexander C. Paul, August 2019
!!
!!    Removing the restrictions on the sum over the triples excitations
!!    introduces a factor similar to (1+delta_ai,bj) in the doubles.
!!    Usually this factor cancles but it is e.g. needed in EOM-CC3.
!!
!!       sum_{ai >= bj >= ck} = 1/6 sum_aibjck (1 + delta_ai,bj 
!!                                                + delta_ai,ck
!!                                                + delta_bj,ck
!!                                                + 2 delta_ai,bj delta_ai,ck)
!!
!!    This routine scales an array of amplitudes (for single a,b,c)
!!    with the expression in paranthesis
!!    NB:
!!    c can only be equal to a if a == b == c for which R_ijk = 0
!!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(inout) :: R_ijk
!
      integer :: i, j
!
      if (a .eq. b .and. a .ne. c) then
!
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               R_ijk(i,i,j) = two*R_ijk(i,i,j)
!
            enddo
         enddo
!
      else if (b .eq. c .and. a .ne. c) then
!
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               R_ijk(i,j,j) = two*R_ijk(i,j,j)
!
            enddo 
         enddo
!
      end if
!
   end subroutine scale_triples_biorthonormal_factor_abc_cc3
!
!
   subroutine get_cvs_projector_cc3(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai, j, b, bj, aibj
!
      call zero_array(projector, wf%n_es_amplitudes)
!
      do core = 1, n_cores
!
        i = core_MOs(core)
!
!$omp parallel do private (a, ai, j, b, bj, aibj)
        do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai) = one
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!
                  projector(aibj + (wf%n_o)*(wf%n_v)) = one
!
               enddo
            enddo
        enddo
!$omp end parallel do
!
     enddo
!
   end subroutine get_cvs_projector_cc3
!
!
   subroutine biorthonormalize_L_and_R_cc3(wf, energy_threshold, residual_threshold, skip_states)
!!
!!    Biorthonormalize the left and right eigenvectors
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    NB: Folding in the Triples part of the excited states leads to 
!!        non biorthogonal eigenvectors of the jacobian
!!        Therefore, a biorthonormalization would have to include the singles, 
!!        doubles and full triples part which is not feasible due to memory.
!!        
!!        It should be possible though to restart from biorthogonal states 
!!        from a CCSD-calculation. Biorthorgonal CC3 states should be obtained 
!!        from that because DIIS does not tend to mix the degenerate states 
!!
!!    For that: 
!!       - check if left and right excitation energies are consistent
!!       - check for degenerate states
!!       - check for and discard parallel states
!!       - If degeneracies are present: Print Warning and biorthonormalize
!!       - Normalize and write to file
!!
      implicit none
!
      class(cc3) :: wf
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
      real(dp) :: LT_R
      real(dp) :: ddot
!
      integer :: state, current_state, p, q
      integer :: n_degeneracy, reduced_degeneracy_r, reduced_degeneracy_l
      integer :: n_overlap_zero, state_nonzero_overlap
!
      logical :: biorthonormalize
!
!     Loop through states look for degenerate states (skip degeneracy)
!     discard parallel states
!     binormalize degenerate states
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
                  call output%warning_msg('Cannot guarantee that degenerate &
                                          &roots are biorthogonal in cc3.')
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
!                 The biorthonormalization will not work in CC3 
!                 unless we include the full triples
!
!                 Print warning but still binormalize the left 
!                 and right vector with the largest overlap
!
                  else
!
                     call output%printf('Binormalize degenerate left state to the right &
                                        &state with the largest overlap', &
                                         pl='n',fs='(/t6,a)')
!
                     call output%warning_msg('Cannot guarantee that degenerate roots &
                                             &are biorthogonal in cc3.')
                     call output%warning_msg('Oscillator strength cannot be trusted &
                                             &for the degenerate states.')
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
!                    Zero out R state that has already been used
!                    Thus, overlap_lr(q) will be zero and this state will not be selected again
!
                     call zero_array(R(:,order(1)), wf%n_es_amplitudes)
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
   end subroutine biorthonormalize_L_and_R_cc3
!
!
end module cc3_class
