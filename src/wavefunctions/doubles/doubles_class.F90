!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module doubles_class
!
!!
!!    Abstract doubles class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019
!!
!!    Doubles class which contains routines that are in common 
!!    for both CCSD and CC2.
!!
!!    Both CCSD and CC2 inherit from this doubles class 
!!    but need specialized routines.
!!
!
   use ccs_class
   use array_utilities, only : scale_diagonal
!
   implicit none
!
   type, abstract, extends(ccs) :: doubles
!
      real(dp),    dimension(:), allocatable :: t2
      complex(dp), dimension(:), allocatable :: t2_complex
!
      real(dp),    dimension(:), allocatable :: t2bar
      complex(dp), dimension(:), allocatable :: t2bar_complex
!
      real(dp),    dimension(:,:,:,:), allocatable :: u_aibj ! 2t_aibj - t_ajbi
      complex(dp), dimension(:,:,:,:), allocatable :: u_aibj_complex
!
      integer :: n_t2
!
      type(sequential_file) :: jacobian_a1_intermediate_vv
      type(sequential_file) :: jacobian_a1_intermediate_oo
!
      type(sequential_file) :: jacobian_transpose_a1_intermediate_vv
      type(sequential_file) :: jacobian_transpose_a1_intermediate_oo
!
   contains
!
!     Omega procedures in common for doubles and doubles
!
      procedure :: omega_doubles_a1                      => omega_doubles_a1_doubles
      procedure :: omega_doubles_a1_complex              => omega_doubles_a1_doubles_complex
!
      procedure :: omega_doubles_b1                      => omega_doubles_b1_doubles
      procedure :: omega_doubles_b1_complex              => omega_doubles_b1_doubles_complex
!
      procedure :: omega_doubles_c1                      => omega_doubles_c1_doubles
      procedure :: omega_doubles_c1_complex              => omega_doubles_c1_doubles_complex
!
!     Jacobian transformation procedures in common for doubles and doubles
!
      procedure :: get_orbital_differences               => get_orbital_differences_doubles
!
      procedure :: jacobian_doubles_a1                   => jacobian_doubles_a1_doubles
      procedure :: jacobian_doubles_b1                   => jacobian_doubles_b1_doubles
      procedure :: jacobian_doubles_c1                   => jacobian_doubles_c1_doubles
      procedure :: jacobian_doubles_d1                   => jacobian_doubles_d1_doubles
      procedure :: jacobian_doubles_a2                   => jacobian_doubles_a2_doubles
!
      procedure :: save_jacobian_a1_intermediates        => save_jacobian_a1_intermediates_doubles
!
!     Jacobian transpose transformation procedures in common for CC2 and CCSD
!
      procedure :: jacobian_transpose_doubles_a1         => jacobian_transpose_doubles_a1_doubles
      procedure :: jacobian_transpose_doubles_a1_complex => jacobian_transpose_doubles_a1_doubles_complex
!
      procedure :: jacobian_transpose_doubles_b1         => jacobian_transpose_doubles_b1_doubles
      procedure :: jacobian_transpose_doubles_b1_complex => jacobian_transpose_doubles_b1_doubles_complex
!
      procedure :: jacobian_transpose_doubles_a2         => jacobian_transpose_doubles_a2_doubles
      procedure :: jacobian_transpose_doubles_a2_complex => jacobian_transpose_doubles_a2_doubles_complex
!
      procedure :: save_jacobian_transpose_a1_intermediates         => save_jacobian_transpose_a1_intermediates_doubles
      procedure :: save_jacobian_transpose_a1_intermediates_complex => save_jacobian_transpose_a1_intermediates_doubles_complex
!
!     Initialization/destruction procedures
!
      procedure :: initialize_t2                         => initialize_t2_doubles
      procedure :: initialize_t2_complex                 => initialize_t2_doubles_complex
!
      procedure :: destruct_t2                           => destruct_t2_doubles
      procedure :: destruct_t2_complex                   => destruct_t2_doubles_complex
!
      procedure :: initialize_t2bar                      => initialize_t2bar_doubles
      procedure :: initialize_t2bar_complex              => initialize_t2bar_doubles_complex
!
      procedure :: destruct_t2bar                        => destruct_t2bar_doubles
      procedure :: destruct_t2bar_complex                => destruct_t2bar_doubles_complex
!
      procedure :: initialize_u_aibj                     => initialize_u_aibj_doubles
      procedure :: initialize_u_aibj_complex             => initialize_u_aibj_doubles_complex
!
      procedure :: destruct_u_aibj                       => destruct_u_aibj_doubles
      procedure :: destruct_u_aibj_complex               => destruct_u_aibj_doubles_complex
!
!     File handling procedures
!
      procedure :: read_doubles_vector                   => read_doubles_vector_doubles
      procedure :: save_doubles_vector                   => save_doubles_vector_doubles
      procedure :: read_excitation_vector_file           => read_excitation_vector_file_doubles
      procedure :: save_excitation_vector_on_file        => save_excitation_vector_on_file_doubles
      procedure :: get_restart_vector                    => get_restart_vector_doubles
!
!     Projectors for excited and ionized states
!
      procedure :: get_cvs_projector                     => get_cvs_projector_doubles
      procedure :: get_rm_core_projector                 => get_rm_core_projector_doubles
      procedure :: get_ip_projector                      => get_ip_projector_doubles
!
!     Properties and densities
!
      procedure :: construct_gs_density                  => construct_gs_density_doubles
      procedure :: construct_gs_density_complex          => construct_gs_density_doubles_complex
!
      procedure :: construct_left_transition_density     => construct_left_transition_density_doubles
      procedure :: construct_right_transition_density    => construct_right_transition_density_doubles
!
      procedure :: density_doubles_mu_ref_oo             => density_doubles_mu_ref_oo_doubles
      procedure :: density_doubles_mu_ref_oo_complex     => density_doubles_mu_ref_oo_doubles_complex

      procedure :: density_doubles_mu_ref_vv             => density_doubles_mu_ref_vv_doubles
      procedure :: density_doubles_mu_ref_vv_complex     => density_doubles_mu_ref_vv_doubles_complex
!
      procedure :: density_doubles_mu_ref_ov             => density_doubles_mu_ref_ov_doubles
      procedure :: density_doubles_mu_ref_ov_complex     => density_doubles_mu_ref_ov_doubles_complex
!
      procedure :: density_doubles_mu_nu_ov              => density_doubles_mu_nu_ov_doubles
      procedure :: density_doubles_mu_nu_vo              => density_doubles_mu_nu_vo_doubles
!
!     Eta and cxi
!
      procedure :: construct_eom_etaX                    => construct_eom_etaX_doubles
!
      procedure :: construct_etaX                        => construct_etaX_doubles
!
      procedure :: etaX_eom_a                            => etaX_eom_a_doubles
!
      procedure :: etaX_doubles_a1                       => etaX_doubles_a1_doubles
      procedure :: etaX_doubles_a2                       => etaX_doubles_a2_doubles
      procedure :: etaX_doubles_b2                       => etaX_doubles_b2_doubles
!
      procedure :: construct_xiX                         => construct_xiX_doubles
!
      procedure :: xiX_doubles_a1                        => xiX_doubles_a1_doubles
      procedure :: xiX_doubles_a2                        => xiX_doubles_a2_doubles
      procedure :: etaX_eom_doubles_a1                   => etaX_eom_doubles_a1_doubles
!
!     Procedures related to time dependency
!
      procedure :: make_doubles_complex                  => make_doubles_complex_doubles
      procedure :: cleanup_doubles_complex               => cleanup_doubles_complex_doubles
!
      procedure :: print_amplitude_info                  => print_amplitude_info_doubles
!
   end type doubles
!
   interface
!
      include "omega_doubles_interface.F90"
      include "jacobian_doubles_interface.F90"
      include "jacobian_transpose_doubles_interface.F90"
      include "file_handling_doubles_interface.F90"
      include "initialize_destruct_doubles_interface.F90"
      include "zop_doubles_interface.F90"
      include "response_doubles_interface.F90"
      include "complex_doubles_interface.F90"
!
      include "autogenerated_complex_files/initialize_destruct_doubles_complex_interface.F90"
      include "autogenerated_complex_files/jacobian_transpose_doubles_complex_interface.F90"
      include "autogenerated_complex_files/omega_doubles_complex_interface.F90"
      include "autogenerated_complex_files/zop_doubles_complex_interface.F90"
!
   end interface
!
contains
!
!
   subroutine print_amplitude_info_doubles(wf)
!!
!!    Print amplitude info
!!    Written by Sarai D. Folkestad, Dec 2019
!!
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      call wf%ccs%print_amplitude_info()
!
      call output%printf('m', 'Double excitation amplitudes:  (i0)', &
            ints=[wf%n_t2], fs='(t6,a)')
!
   end subroutine print_amplitude_info_doubles
!
!
   subroutine get_cvs_projector_doubles(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
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
   end subroutine get_cvs_projector_doubles
!
!
   subroutine get_rm_core_projector_doubles(wf, projector, n_cores, core_MOs)
!!
!!    Get remove core projector
!!    Written by Sarai D. Folkestad, 2021
!!
!
      use array_utilities, only: constant_array
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai, j, b, bj, aibj
!
      call constant_array(projector, wf%n_es_amplitudes, one)
!
      do core = 1, n_cores
!
        i = core_MOs(core)
!
!$omp parallel do private (a, ai, j, b, bj, aibj)
        do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai) = zero
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!
                  projector(aibj + (wf%n_o)*(wf%n_v)) = zero
!
               enddo
            enddo
        enddo
!$omp end parallel do
!
     enddo
!
   end subroutine get_rm_core_projector_doubles
!
!
   subroutine get_ip_projector_doubles(wf, projector)
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
      class(doubles), intent(in) :: wf
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
   end subroutine get_ip_projector_doubles
!
!
   subroutine get_orbital_differences_doubles(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      integer, intent(in) :: N 
!
      real(dp), dimension(N), intent(out) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj
!
      call wf%ccs%get_orbital_differences(orbital_differences, wf%n_t1)
!
      if (N .eq. wf%n_t1) return ! Requested only singles orbital differences
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            do j = 1, wf%n_o 
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j-1) + b 
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     orbital_differences(aibj + wf%n_t1) = wf%orbital_energies(a + wf%n_o)      &
                                                         - wf%orbital_energies(i)               &
                                                         + wf%orbital_energies(b + wf%n_o)      &
                                                         - wf%orbital_energies(j)
!
                  endif
!
               enddo
            enddo  
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_orbital_differences_doubles
!
!
end module doubles_class
