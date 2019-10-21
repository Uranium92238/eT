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
module abstract_doubles_class
!
!!
!!    Abstract doubles class module
!!
!!    Abstract doubles class which both cc2 and ccsd inherit from.
!!
!
   use ccs_class
!
   implicit none
!
   type, abstract, extends(ccs) :: abstract_doubles
!
      real(dp), dimension(:), allocatable :: t2
      real(dp), dimension(:), allocatable :: t2bar
!
      integer :: n_t2
!
      type(sequential_file) :: r2_file
      type(sequential_file) :: l2_file
!
      type(sequential_file) :: jacobian_a1_intermediate_vv
      type(sequential_file) :: jacobian_a1_intermediate_oo
!
   contains
!
!     Omega routines in common for abstract_doubles and abstract_doubles
!
      procedure :: omega_doubles_a1 => omega_doubles_a1_abstract_doubles
      procedure :: omega_doubles_b1 => omega_doubles_b1_abstract_doubles
      procedure :: omega_doubles_c1 => omega_doubles_c1_abstract_doubles
!
!     Jacobian transformation routines in common for abstract_doubles and abstract_doubles
!
      procedure :: jacobian_doubles_a1 => jacobian_doubles_a1_abstract_doubles
      procedure :: jacobian_doubles_b1 => jacobian_doubles_b1_abstract_doubles
      procedure :: jacobian_doubles_c1 => jacobian_doubles_c1_abstract_doubles
      procedure :: jacobian_doubles_d1 => jacobian_doubles_d1_abstract_doubles
      procedure :: jacobian_doubles_a2 => jacobian_doubles_a2_abstract_doubles
!
      procedure :: save_jacobian_a1_intermediates  => save_jacobian_a1_intermediates_abstract_doubles
!
!     Jacobian transpose transformation routines in common for abstract_doubles and abstract_doubles
!
      procedure :: jacobian_transpose_doubles_a1   => jacobian_transpose_doubles_a1_abstract_doubles
      procedure :: jacobian_transpose_doubles_b1   => jacobian_transpose_doubles_b1_abstract_doubles
      procedure :: jacobian_transpose_doubles_a2   => jacobian_transpose_doubles_a2_abstract_doubles
!
!     Initializations and destructions   
!
      procedure :: initialize_t2 => initialize_t2_abstract_doubles
      procedure :: destruct_t2   => destruct_t2_abstract_doubles
!
      procedure :: initialize_t2bar => initialize_t2bar_abstract_doubles
      procedure :: destruct_t2bar   => destruct_t2bar_abstract_doubles
!
!     Projectors for excited and ionized states
!
      procedure :: get_cvs_projector   => get_cvs_projector_abstract_doubles
      procedure :: get_ip_projector    => get_ip_projector_abstract_doubles
!
!     Ground state density
!
      procedure :: construct_gs_density          => construct_gs_density_abstract_doubles
      procedure :: gs_one_el_density_doubles_oo  => gs_one_el_density_doubles_oo_abstract_doubles
      procedure :: gs_one_el_density_doubles_vv  => gs_one_el_density_doubles_vv_abstract_doubles
      procedure :: gs_one_el_density_doubles_ov  => gs_one_el_density_doubles_ov_abstract_doubles
!
!     Transition densities
!
      procedure :: construct_left_transition_density     => construct_left_transition_density_abstract_doubles
      procedure :: construct_right_transition_density    => construct_right_transition_density_abstract_doubles
      procedure :: right_transition_density_doubles_ov   => right_transition_density_doubles_ov_abstract_doubles
      procedure :: right_transition_density_doubles_vo   => right_transition_density_doubles_vo_abstract_doubles
!
!     Eta and cxi
!
      procedure :: construct_eom_etaX  => construct_eom_etaX_abstract_doubles
      procedure :: construct_etaX      => construct_etaX_abstract_doubles 
!
      procedure :: etaX_eom_a       => etaX_eom_a_abstract_doubles    
      procedure :: etaX_doubles_a1  => etaX_doubles_a1_abstract_doubles    
      procedure :: etaX_doubles_a2  => etaX_doubles_a2_abstract_doubles    
      procedure :: etaX_doubles_b2  => etaX_doubles_b2_abstract_doubles    
!     
      procedure :: construct_csiX   => construct_csiX_abstract_doubles  
      procedure :: csiX_doubles_a1  => csiX_doubles_a1_abstract_doubles    
      procedure :: csiX_doubles_a2  => csiX_doubles_a2_abstract_doubles    
!
      procedure :: etaX_eom_doubles_a1 => etaX_eom_doubles_a1_abstract_doubles
!
   end type abstract_doubles
!
   interface
!
      include "omega_abstract_doubles_interface.F90"
      include "jacobian_abstract_doubles_interface.F90"
      include "jacobian_transpose_abstract_doubles_interface.F90"
      include "file_handling_abstract_doubles_interface.F90"
      include "initialize_destruct_abstract_doubles_interface.F90"
      include "zop_abstract_doubles_interface.F90"
      include "fop_abstract_doubles_interface.F90"
!
   end interface
!
contains
!
!
   subroutine get_cvs_projector_abstract_doubles(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(abstract_doubles), intent(inout) :: wf
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
   end subroutine get_cvs_projector_abstract_doubles
!
!
   subroutine get_ip_projector_abstract_doubles(wf, projector)
!!
!!    Get IP projector 
!!    Written by Sarai D. Folkestad, Aug 2019
!!
!!    Constructs and returns the projector
!!    for an IP calculation (valence).
!!
!!
      implicit none
!
      class(abstract_doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes),intent(out) :: projector
!
      integer :: A, I, AI, B, J, BJ, AIBJ
!
      call zero_array(projector, wf%n_es_amplitudes)
!
      do I = 1, wf%n_o
         do A = wf%n_v - wf%n_bath_orbitals + 1, wf%n_v
!
            AI = wf%n_v*(I-1) + A
            projector(AI) = one
!
            do J = 1, wf%n_o
               do B = 1, wf%n_v
!
                  BJ = wf%n_v*(J-1) + B
                  AIBJ = max(AI, BJ)*(max(AI, BJ)-3)/2 + AI + BJ
!
                  projector(wf%n_t1 + AIBJ) = one
!
               enddo
            enddo 
!
         enddo
      enddo
!
   end subroutine get_ip_projector_abstract_doubles
!
!
end module abstract_doubles_class
