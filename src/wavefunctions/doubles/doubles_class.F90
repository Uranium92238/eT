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
   use ccs_class, only: ccs
!
   use parameters
   use global_out, only: output
   use timings_class, only: timings
   use memory_manager_class, only: mem
   use stream_file_class, only: stream_file
   use amplitude_file_storer_class, only: amplitude_file_storer
!
   implicit none
!
   type, extends(ccs) :: doubles
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
      type(stream_file) :: jacobian_a1_intermediate_vv
      type(stream_file) :: jacobian_a1_intermediate_oo
!
      type(stream_file) :: jacobian_transpose_a1_intermediate_vv
      type(stream_file) :: jacobian_transpose_a1_intermediate_oo
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
      procedure, public :: jacobian_transformation &
                        => jacobian_transformation_doubles
!
      procedure :: prepare_for_jacobian &
                => prepare_for_jacobian_doubles

!
      procedure, private :: jacobian_doubles_a1
      procedure, private :: jacobian_doubles_b1
      procedure, private :: jacobian_doubles_c1
      procedure, private :: jacobian_doubles_d1
      procedure, private :: jacobian_doubles_a2
      procedure, public  :: jacobian_doubles_b2 ! public because of intel
!
      procedure, private :: save_jacobian_a1_intermediates
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
      procedure :: get_restart_vector                    => get_restart_vector_doubles
!
      procedure :: get_full_multipliers &
                => get_multipliers_doubles
!
      procedure :: get_full_multipliers_complex &
                => get_multipliers_doubles_complex
!
!     Projectors for excited and ionized states
!
      procedure :: cvs_projection                        => cvs_projection_doubles
      procedure :: rm_core_projection                    => rm_core_projection_doubles
      procedure :: ip_projection                         => ip_projection_doubles
!
!     Properties and densities
!
      procedure :: construct_gs_density &
                => construct_gs_density_doubles
      procedure :: construct_gs_density_complex &
                => construct_gs_density_doubles_complex
!
      procedure :: mu_ref_density_terms &
                => mu_ref_density_terms_doubles
      procedure :: mu_ref_density_terms_complex &
                => mu_ref_density_terms_doubles_complex
!
      procedure :: construct_right_transition_density &
                => construct_right_transition_density_doubles
!
      procedure :: mu_nu_density_terms &
                => mu_nu_density_terms_doubles
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
!
      procedure :: construct_cntos_for_plotting &
                => construct_cntos_for_plotting_doubles
      procedure :: construct_ntos_or_cntos &
                => construct_ntos_or_cntos_doubles
!
!     F transformation
!
      procedure :: F_doubles_a1_1 => F_doubles_a1_1_doubles
      procedure :: F_doubles_a2_1 => F_doubles_a2_1_doubles
      procedure :: F_doubles_a1_2 => F_doubles_a1_2_doubles
      procedure :: F_doubles_b1_2 => F_doubles_b1_2_doubles
      procedure :: F_doubles_c1_2 => F_doubles_c1_2_doubles
!
      procedure :: F_x_mu_transformation => F_x_mu_transformation_doubles
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
      include "mean_value_doubles_interface.F90"
      include "response_doubles_interface.F90"
      include "complex_doubles_interface.F90"
      include "F_doubles_interface.F90"
!
      include "generated_complex_files/initialize_destruct_doubles_complex_interface.F90"
      include "generated_complex_files/jacobian_transpose_doubles_complex_interface.F90"
      include "generated_complex_files/omega_doubles_complex_interface.F90"
      include "generated_complex_files/mean_value_doubles_complex_interface.F90"
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
   subroutine cvs_projection_doubles(wf, vector, n_cores, core_MOs)
!!
!!    CVS projection
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: vector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: i, a, ai, j, b, bj, aibj
!
      call wf%ccs%cvs_projection(vector, n_cores, core_MOs)
!
!$omp parallel do private(i, a, j, b, ai, bj, aibj) collapse(2)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  if(all(core_MOs /= i) .and. all(core_MOs /= j)) then
!
                     ai = wf%n_v*(i - 1) + a
                     bj = wf%n_v*(j - 1) + b
                     aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!
                     vector(aibj + (wf%n_o)*(wf%n_v)) = zero
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine cvs_projection_doubles
!
!
   subroutine rm_core_projection_doubles(wf, vector, n_cores, core_MOs)
!!
!!    Remove core projection
!!    Written by Sarai D. Folkestad, 2021
!!
!
      use array_initialization, only: constant_array
!
      implicit none
!
      class(doubles), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: vector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: i, a, ai, j, b, bj, aibj, core
!
      call wf%ccs%rm_core_projection(vector(1:wf%n_t1), n_cores, core_MOs)
!
!$omp parallel do private(i, a, j, b, ai, bj, aibj, core) collapse(2)
      do core = 1, n_cores
         do a = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  i = core_MOs(core)
!
                  ai = wf%n_v*(i - 1) + a
                  bj = wf%n_v*(j - 1) + b
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!
                  vector(aibj + (wf%n_o)*(wf%n_v)) = zero
!
               enddo
            enddo
         enddo
      enddo

!$omp end parallel do
!
   end subroutine rm_core_projection_doubles
!
!
   subroutine ip_projection_doubles(wf, vector)
!!
!!    IP projection
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
      real(dp), dimension(wf%n_es_amplitudes),intent(inout) :: vector
!
      integer :: A, I, AI, B, J, BJ, AIBJ
!
      call wf%ccs%ip_projection(vector)
!
!$omp parallel do private(i, a, j, b, ai, bj, aibj) collapse(2)
      do I = 1, wf%n_o
         do A = 1, wf%n_v
            do J = 1, wf%n_o
               do B = 1, wf%n_v
!
                  if ((A <= wf%n_v - wf%n_bath_orbitals) .and. &
                      (B <= wf%n_v - wf%n_bath_orbitals)) then
!
                     AI = wf%n_v*(I-1) + A
                     BJ = wf%n_v*(J-1) + B
                     AIBJ = max(AI, BJ)*(max(AI, BJ)-3)/2 + AI + BJ
!
                     vector(wf%n_t1 + AIBJ) = zero
                  endif
!
               enddo
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine ip_projection_doubles
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
   subroutine construct_cntos_for_plotting_doubles(wf, cntos, k, n_sig_v, &
                                                   n_sig_o, l_or_r, threshold)
!!
!!    Construct CNTOs for plotting
!!    Written by Sarai D. Folkestad, May 2020
!!
!!    'cntos' : array to store CNTO orbital coefficients
!!
!!    'k' : which excited state to construct CNTOs for
!!
!!    'n_sig_v' and 'n_sig_o' : The number of significant occupied and virtuals
!!
!!    'l_or_r' : excited state vector type
!!
!!     Constructs the M and N matrices
!!
!!       M_ij += sum_a R_ai R_aj + 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R^k_aibjR^k_blaj)
!!       N_ab += sum_i R_ai R_bi + 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R^k_aicjR^k_bicj)
!!
!!    and diagonalizes them.
!!    The CNTOs are then constructed and stored in 'cntos'
!!
!!    Significant occupied and virtual orbitals (n_sig_o and n_sig_v)
!!    are determined by considering the
!!    eigenvalues of M and N, using the criterion
!!
!!       1 - sum_i e_i < threshold,
!!
!!    where the eigenvalues are ordered from largest to smallest.
!!
!
      use cnto_tool_class, only: cnto_tool
!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(out) :: cntos
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in)  :: k
      integer, intent(out) :: n_sig_v, n_sig_o
!
      character(len=*), intent(in) :: l_or_r
!
      real(dp), dimension(:), allocatable :: X
!
      type(cnto_tool), allocatable :: orbital_tool
!
      integer :: n_singles_and_doubles
!
      call output%printf('n', '- Constructing CNTOs for state (i0)', &
                         ints=[k], ffs='(/t3,a)')
!
      call mem%alloc(X, wf%n_es_amplitudes)
!
      call wf%read_excited_state(X, k, k, trim(l_or_r))
!
      n_singles_and_doubles = wf%n_t1 + wf%n_t1*(wf%n_t1 + 1)/2
!
      orbital_tool = cnto_tool(wf%n_o, wf%n_v, wf%ao%n, n_singles_and_doubles)
!
      call orbital_tool%initialize()
!
      call orbital_tool%add_excited_state(X(1:n_singles_and_doubles))
!
      call orbital_tool%transform_orbitals(wf%orbital_coefficients, cntos)
!
      call orbital_tool%get_n_active_orbitals(n_sig_o, n_sig_v, threshold)
!
      call orbital_tool%cleanup
!
      call mem%dealloc(X, wf%n_es_amplitudes)
!
   end subroutine construct_cntos_for_plotting_doubles
!
!
   subroutine construct_ntos_or_cntos_doubles(wf, orbitals, k, &
                                              n_sig_v, n_sig_o, l_or_r, type_,threshold)
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
      class(doubles) :: wf
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
         call wf%construct_cntos_for_plotting(orbitals, k, n_sig_v, n_sig_o,l_or_r, threshold)
!
      endif
!
   end subroutine construct_ntos_or_cntos_doubles
!
!
   subroutine get_multipliers_doubles(wf, multipliers)
!!
!!    Get multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes) :: multipliers
!
      call dcopy(wf%n_t1, wf%t1bar, 1, multipliers, 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, multipliers(wf%n_t1 + 1:), 1)
!
   end subroutine get_multipliers_doubles
!
!
   subroutine get_multipliers_doubles_complex(wf, multipliers)
!!
!!    Get multipliers complex
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(doubles), intent(in) :: wf
!
      complex(dp), dimension(wf%n_es_amplitudes) :: multipliers
!
      call zcopy(wf%n_t1, wf%t1bar_complex, 1, multipliers, 1)
      call zcopy(wf%n_t2, wf%t2bar_complex, 1, multipliers(wf%n_t1 + 1:), 1)
!
   end subroutine get_multipliers_doubles_complex
!
!
end module doubles_class
