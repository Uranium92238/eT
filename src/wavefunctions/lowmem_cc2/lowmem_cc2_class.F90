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
module lowmem_cc2_class
!
!!
!!    Low-memory coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto and Alexander C. Paul, 2018
!!
!!    Version of CC2 that has an O(M^2) memory requirement - making it suitable
!!    to treat larger systems than the standard CC2 wavefunction.
!!
!
   use ccs_class, only: ccs
!
   use parameters
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use stream_file_class, only: stream_file
   use direct_stream_file_class, only: direct_stream_file
   use batching_index_class, only: batching_index
!
   implicit none
!
   type, extends(ccs) :: lowmem_cc2
!
      type(stream_file) :: jacobian_a1_intermediate_vv
      type(stream_file) :: jacobian_a1_intermediate_oo
!
   contains
!
      procedure :: construct_fock   => construct_fock_lowmem_cc2
!
!     Omega
!
      procedure :: construct_omega &
                => construct_omega_lowmem_cc2
      procedure :: omega_cc2 &
                => omega_cc2_lowmem_cc2
!
      procedure :: setup_L_Jvo &
                => setup_L_Jvo_lowmem_cc2
      procedure :: setup_L_Jov &
                => setup_L_Jov_lowmem_cc2
      procedure :: setup_L_Joo &
                => setup_L_Joo_lowmem_cc2
!
      procedure :: construct_contravariant_t2_single_ij &
                => construct_contravariant_t2_single_ij_lowmem_cc2
      procedure :: make_contravariant_doubles_single_ij &
                => make_contravariant_doubles_single_ij_lowmem_cc2
!
      procedure :: omega_cc2_fock &
                => omega_cc2_fock_lowmem_cc2
      procedure :: omega_cc2_v2o2J &
                => omega_cc2_v2o2J_lowmem_cc2
      procedure :: omega_cc2_Jv2o &
                => omega_cc2_Jv2o_lowmem_cc2
!
      procedure :: calculate_energy &
                => calculate_energy_lowmem_cc2
      procedure :: get_electronic_dipole &
                => get_electronic_dipole_lowmem_cc2
      procedure :: get_electronic_quadrupole &
                => get_electronic_quadrupole_lowmem_cc2
!
!     Jacobian transformation
!
      procedure :: construct_Jacobian_transform       => construct_Jacobian_transform_lowmem_cc2
!
      procedure :: prepare_for_approximate_Jacobians  &
                => prepare_for_approximate_Jacobians_lowmem_cc2
!
      procedure :: approximate_Jacobian_transform     &
                => approximate_Jacobian_transform_lowmem_cc2
!
      procedure :: effective_jacobian_transformation  => effective_jacobian_transformation_lowmem_cc2
!
      procedure :: prepare_for_jacobian               => prepare_for_jacobian_lowmem_cc2
      procedure :: save_jacobian_a1_2_intermediate    => save_jacobian_a1_2_intermediate_lowmem_cc2
      procedure :: save_jacobian_a1_3_intermediate    => save_jacobian_a1_3_intermediate_lowmem_cc2
      procedure :: jacobian_cc2_a1                    => jacobian_cc2_a1_lowmem_cc2
!
      procedure :: effective_jacobian_cc2_a1 => effective_jacobian_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_cc2_b1 => effective_jacobian_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_cc2_c1 => effective_jacobian_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_cc2_d1 => effective_jacobian_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_cc2_e1 => effective_jacobian_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_cc2_f1 => effective_jacobian_cc2_f1_lowmem_cc2
!
!     Transpose Jacobian transformation
!
      procedure :: effective_jacobian_transpose_transformation &
               => effective_jacobian_transpose_transformation_lowmem_cc2
!
      procedure :: jacobian_transpose_cc2_a1 => jacobian_transpose_cc2_a1_lowmem_cc2
      procedure :: jacobian_transpose_cc2_b1 => jacobian_transpose_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_a1   => effective_jacobian_transpose_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_b1   => effective_jacobian_transpose_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_c1   => effective_jacobian_transpose_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_d1   => effective_jacobian_transpose_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_e1   => effective_jacobian_transpose_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_f1   => effective_jacobian_transpose_cc2_f1_lowmem_cc2
!
!     Initialize wavefunction
!
      procedure :: initialize => initialize_lowmem_cc2
!
   end type lowmem_cc2
!
!
   interface
!
      include "omega_lowmem_cc2_interface.F90"
      include "jacobian_lowmem_cc2_interface.F90"
      include "jacobian_transpose_lowmem_cc2_interface.F90"
      include "mean_value_lowmem_cc2_interface.F90"
      include "fock_lowmem_cc2_interface.F90"
!
   end interface
!
!
contains
!
!
   subroutine initialize_lowmem_cc2(wf, template_wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use wavefunction_class, only: wavefunction
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'low memory cc2'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
      wf%need_g_abcd     = .false.
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
   end subroutine initialize_lowmem_cc2
!
!
   subroutine construct_Jacobian_transform_lowmem_cc2(wf, r_or_l, X, R, w)
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
      class(lowmem_cc2), intent(inout) :: wf
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
         call output%error_msg('w is missing in construct_Jacobian_transform for lowmem_cc2')
!
      endif
!
!     Compute the transformed matrix
!
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
   end subroutine construct_Jacobian_transform_lowmem_cc2
!
!
   subroutine approximate_Jacobian_transform_lowmem_cc2(wf, r_or_l, X, R, w)
!!
!!    Approximate Jacobian transform
!!    Written by Eirik F. Kjønstad, Mar 2021
!!
!!    Wrapper for a lower-level Jacobian transformation that is the best approximation
!!    with a lower computational scaling.
!!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: X
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: R
!
      real(dp), intent(in), optional :: w
!
      call wf%ccs%construct_Jacobian_transform(r_or_l, X, R, w)
!
   end subroutine approximate_Jacobian_transform_lowmem_cc2
!
!
   subroutine prepare_for_approximate_Jacobians_lowmem_cc2(wf, r_or_l)
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
      class(lowmem_cc2), intent(inout) :: wf
!
      character(len=*), intent(in) :: r_or_l
!
      call wf%ccs%prepare_for_Jacobians(r_or_l)
!
   end subroutine prepare_for_approximate_Jacobians_lowmem_cc2
!
!
end module lowmem_cc2_class
