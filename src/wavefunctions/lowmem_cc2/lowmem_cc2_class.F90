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
module lowmem_cc2_class
!
!!
!!    Low-memory coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto and Alexander C. Paul, 2018
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: lowmem_cc2
!
      type(sequential_file) :: jacobian_b1_intermediate_vv
      type(sequential_file) :: jacobian_b1_intermediate_oo
!
   contains
!
      procedure :: construct_omega  => construct_omega_lowmem_cc2
!
      procedure :: omega_cc2_a1     => omega_cc2_a1_lowmem_cc2
      procedure :: omega_cc2_b1     => omega_cc2_b1_lowmem_cc2
      procedure :: omega_cc2_c1     => omega_cc2_c1_lowmem_cc2
!
      procedure :: calculate_energy => calculate_energy_lowmem_cc2
!
      procedure :: construct_Jacobian_transform       => construct_Jacobian_transform_lowmem_cc2
!
      procedure :: effective_jacobian_transformation  => effective_jacobian_transformation_lowmem_cc2
!
      procedure :: prepare_for_jacobian               => prepare_for_jacobian_lowmem_cc2
      procedure :: save_jacobian_b1_2_intermediate    => save_jacobian_b1_2_intermediate_lowmem_cc2
      procedure :: save_jacobian_b1_3_intermediate    => save_jacobian_b1_3_intermediate_lowmem_cc2      
      procedure :: jacobian_cc2_a1                    => jacobian_cc2_a1_lowmem_cc2
      procedure :: jacobian_cc2_b1                    => jacobian_cc2_b1_lowmem_cc2
!
      procedure :: effective_jacobian_cc2_a1 => effective_jacobian_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_cc2_b1 => effective_jacobian_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_cc2_c1 => effective_jacobian_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_cc2_d1 => effective_jacobian_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_cc2_e1 => effective_jacobian_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_cc2_f1 => effective_jacobian_cc2_f1_lowmem_cc2
!
      procedure :: effective_jacobian_transpose_transformation => effective_jacobian_transpose_transformation_lowmem_cc2
      procedure :: jacobian_transpose_ccs_b1 => jacobian_transpose_ccs_b1_lowmem_cc2
      procedure :: jacobian_transpose_cc2_a1 => jacobian_transpose_cc2_a1_lowmem_cc2
      procedure :: jacobian_transpose_cc2_b1 => jacobian_transpose_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_a1   => effective_jacobian_transpose_cc2_a1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_b1   => effective_jacobian_transpose_cc2_b1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_c1   => effective_jacobian_transpose_cc2_c1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_d1   => effective_jacobian_transpose_cc2_d1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_e1   => effective_jacobian_transpose_cc2_e1_lowmem_cc2
      procedure :: effective_jacobian_transpose_cc2_f1   => effective_jacobian_transpose_cc2_f1_lowmem_cc2
!
   end type lowmem_cc2
!
!
   interface
!
      include "omega_lowmem_cc2_interface.F90"
      include "jacobian_lowmem_cc2_interface.F90"
      include "jacobian_transpose_lowmem_cc2_interface.F90"
      include "zop_lowmem_cc2_interface.F90"
!
   end interface
!
!
   interface lowmem_cc2
!
      procedure :: new_lowmem_cc2
!
   end interface lowmem_cc2
!
!
contains
!
!
   function new_lowmem_cc2(system) result(wf)
!!
!!    New lowmem CC2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(lowmem_cc2) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%name_ = 'low memory cc2'
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
   end function new_lowmem_cc2
!
!
   subroutine construct_Jacobian_transform_lowmem_cc2(wf, r_or_l, X, w)
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
      class(lowmem_cc2), intent(in) :: wf
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
!
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
   end subroutine construct_Jacobian_transform_lowmem_cc2
!
!
end module lowmem_cc2_class
