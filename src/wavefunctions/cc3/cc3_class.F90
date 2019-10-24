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
!!    Written by Rolf H. Myhre and Alexander Paul, 2018-2019
!!
!
   use kinds
   use ccsd_class
   use direct_file_class, only : direct_file
   use sequential_file_class, only : sequential_file
   use io_utilities, only : single_record_reader, compound_record_reader 
   use io_utilities, only : single_record_writer, compound_record_writer 
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
!
      real(dp), dimension(:,:), allocatable :: GS_cc3_density_oo
      real(dp), dimension(:,:), allocatable :: GS_cc3_density_vv
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                => cleanup_cc3
!
!     Routines related to omega
!
      procedure :: construct_omega        => construct_omega_cc3
!
      procedure :: omega_cc3_a            => omega_cc3_a_cc3
      procedure :: omega_cc3_integrals    => omega_cc3_integrals_cc3
      procedure :: omega_cc3_W_calc       => omega_cc3_W_calc_cc3
      procedure :: omega_cc3_eps          => omega_cc3_eps_cc3
      procedure :: omega_cc3_a_n6         => omega_cc3_a_n6_cc3
      procedure :: omega_cc3_a_n7         => omega_cc3_a_n7_cc3
!
!     Routines used for prepare, both left and right
!
      procedure :: prepare_for_jacobian            => prepare_for_jacobian_cc3
      procedure :: prepare_for_jacobian_transpose  => prepare_for_jacobian_transpose_cc3
!
!     Routines for CVS
      procedure :: get_cvs_projector   => get_cvs_projector_cc3 
!
!     Both
      procedure :: prep_cc3_g_lbkc_t_file          => prep_cc3_g_lbkc_t_file_cc3
      procedure :: prep_cc3_jacobian_intermediates => prep_cc3_jacobian_intermediates_cc3
      procedure :: construct_x_intermediates       => construct_x_intermediates_cc3
      procedure :: sort_x_to_abid_and_write        => sort_x_to_abid_and_write_cc3
!
!     Only left
      procedure :: prep_cc3_jacobian_trans_integrals  => prep_cc3_jacobian_trans_integrals_cc3
!
!     Routines related to the jacobian
!
      procedure :: construct_Jacobian_transform       => construct_Jacobian_transform_cc3
!
!     Right hand side transformation
!
      procedure :: effective_jacobian_transformation  => effective_jacobian_transformation_cc3
!
      procedure :: jacobian_cc3_t3_a2  => jacobian_cc3_t3_a2_cc3
!
      procedure :: jacobian_cc3_t3_b2     => jacobian_cc3_t3_b2_cc3
      procedure :: construct_c1_fock      => construct_c1_fock_cc3
      procedure :: jacobian_cc3_b2_fock   => jacobian_cc3_b2_fock_cc3
!
      procedure :: jacobian_cc3_c3_a      => jacobian_cc3_c3_a_cc3
      procedure :: construct_c1_integrals => construct_c1_integrals_cc3
!
!     Routines related to the transpose of the jacobian
!
      procedure :: effective_jacobian_transpose_transformation  &
                                                => effective_jacobian_transpose_transformation_cc3
!
      procedure :: jacobian_transpose_cc3_t3_a1 => jacobian_transpose_cc3_t3_a1_cc3
!
      procedure :: jacobian_transpose_cc3_t3_b1 => jacobian_transpose_cc3_t3_b1_cc3
      procedure :: construct_x_ai_intermediate  => construct_x_ai_intermediate_cc3
!
      procedure :: jacobian_transpose_cc3_c3_a        => jacobian_transpose_cc3_c3_a_cc3
      procedure :: jacobian_transpose_cc3_c3_calc     => jacobian_transpose_cc3_c3_calc_cc3
      procedure :: jacobian_transpose_cc3_a_n7        => jacobian_transpose_cc3_a_n7_cc3
      procedure :: construct_y_intermediates          => construct_y_intermediates_cc3
      procedure :: jacobian_transpose_cc3_c3_a1_y_o   => jacobian_transpose_cc3_c3_a1_y_o_cc3
      procedure :: jacobian_transpose_cc3_c3_b1_y_v   => jacobian_transpose_cc3_c3_b1_y_v_cc3
!
!     Routines related to the multipliers
!
      procedure :: prepare_for_multiplier_equation => prepare_for_multiplier_equation_cc3
      procedure :: construct_multiplier_equation   => construct_multiplier_equation_cc3
!
!     Routines to construct triples amplitudes in batches of a,b,c
!
      procedure :: prep_cc3_integrals_t3_abc_batch => prep_cc3_integrals_t3_abc_batch_cc3
      procedure :: prep_cc3_integrals_R3_abc_batch => prep_cc3_integrals_R3_abc_batch_cc3
      procedure :: prep_cc3_integrals_L3_abc_batch => prep_cc3_integrals_L3_abc_batch_cc3
      procedure :: omega_cc3_W_calc_abc_batch      => omega_cc3_W_calc_abc_batch_cc3
      procedure :: omega_cc3_eps_abc_batch         => omega_cc3_eps_abc_batch_cc3
      procedure :: jacobian_transpose_cc3_c3_calc_abc_batch &
                                                   => jacobian_transpose_cc3_c3_calc_abc_batch_cc3
!
!     Routines related to the ground state density matrix
!
      procedure :: initialize_gs_density  => initialize_gs_density_cc3
      procedure :: destruct_gs_density    => destruct_gs_density_cc3
!
      procedure :: prepare_for_density    => prepare_for_density_cc3
      procedure :: construct_gs_density   => construct_gs_density_cc3
!
      procedure :: gs_one_el_density_cc3_abc => gs_one_el_density_cc3_abc_cc3
      procedure :: one_el_density_cc3_oo_N7  => one_el_density_cc3_oo_N7_cc3
!
      procedure :: gs_one_el_density_cc3_ijk     => gs_one_el_density_cc3_ijk_cc3
      procedure :: one_el_density_cc3_vv_N7     => one_el_density_cc3_vv_N7_cc3
!
      procedure :: construct_y_intermediate_vo3 => construct_y_intermediate_vo3_cc3
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
!
      call wf%write_cc_restart()
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
   end subroutine cleanup_cc3
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
!
!     Compute the transformed matrix
      if (r_or_l .eq. "right") then
!
         call wf%effective_jacobian_transformation(w, X) ! X <- AX
!
      elseif (r_or_l .eq. 'left') then
!
         call wf%effective_jacobian_transpose_transformation(w, X) ! X <- AX
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
      projector = zero
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
end module cc3_class
