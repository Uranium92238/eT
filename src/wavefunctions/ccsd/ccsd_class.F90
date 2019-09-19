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
module ccsd_class
!
!!
!!    Coupled cluster singles and doubles (CCSD) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!
   use abstract_doubles_class
!
   implicit none
!
   type, extends(abstract_doubles) :: ccsd
!
      type(file) :: t2_file, t2bar_file
!
!     Intermediate files 
!
      type(sequential_file) :: jacobian_c2_intermediate_oovo_1
      type(sequential_file) :: jacobian_c2_intermediate_oovo_2
      type(sequential_file) :: jacobian_c2_intermediate_oovo_3
      type(sequential_file) :: jacobian_d2_intermediate
      type(sequential_file) :: jacobian_e2_intermediate
      type(sequential_file) :: jacobian_g2_intermediate_vovo
      type(sequential_file) :: jacobian_h2_intermediate_voov_1
      type(sequential_file) :: jacobian_h2_intermediate_voov_2
      type(sequential_file) :: jacobian_j2_intermediate_oooo
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: initialize_files                            => initialize_files_ccsd
      procedure :: initialize_doubles_files                    => initialize_doubles_files_ccsd
!
!     Routines related to the amplitudes
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_ccsd
      procedure :: set_initial_amplitudes_guess                => set_initial_amplitudes_guess_ccsd
      procedure :: set_t2_to_mp2_guess                         => set_t2_to_mp2_guess_ccsd
      procedure :: set_amplitudes                              => set_amplitudes_ccsd
      procedure :: get_amplitudes                              => get_amplitudes_ccsd
      procedure :: read_amplitudes                             => read_amplitudes_ccsd
      procedure :: save_amplitudes                             => save_amplitudes_ccsd
      procedure :: print_dominant_x2                           => print_dominant_x2_ccsd
      procedure :: print_dominant_amplitudes                   => print_dominant_amplitudes_ccsd
      procedure :: print_dominant_x_amplitudes                 => print_dominant_x_amplitudes_ccsd
!
      procedure :: from_biorthogonal_to_biorthonormal          => from_biorthogonal_to_biorthonormal_ccsd
!
!     Routines related to omega
!
      procedure :: construct_omega                             => construct_omega_ccsd
!
      procedure :: omega_ccsd_a2                               => omega_ccsd_a2_ccsd
      procedure :: omega_ccsd_b2                               => omega_ccsd_b2_ccsd
      procedure :: omega_ccsd_c2                               => omega_ccsd_c2_ccsd
      procedure :: omega_ccsd_d2                               => omega_ccsd_d2_ccsd
      procedure :: omega_ccsd_e2                               => omega_ccsd_e2_ccsd
!
      procedure :: form_newton_raphson_t_estimate              => form_newton_raphson_t_estimate_ccsd
!
!     Routines related to Jacobian transformation
!
      procedure :: jacobian_transform_trial_vector             => jacobian_transform_trial_vector_ccsd
      procedure :: jacobian_ccsd_transformation                => jacobian_ccsd_transformation_ccsd
!
      procedure :: jacobian_ccsd_b2                            => jacobian_ccsd_b2_ccsd
      procedure :: jacobian_ccsd_c2                            => jacobian_ccsd_c2_ccsd
      procedure :: jacobian_ccsd_d2                            => jacobian_ccsd_d2_ccsd
      procedure :: jacobian_ccsd_e2                            => jacobian_ccsd_e2_ccsd
      procedure :: jacobian_ccsd_f2                            => jacobian_ccsd_f2_ccsd
      procedure :: jacobian_ccsd_g2                            => jacobian_ccsd_g2_ccsd
      procedure :: jacobian_ccsd_h2                            => jacobian_ccsd_h2_ccsd
      procedure :: jacobian_ccsd_i2                            => jacobian_ccsd_i2_ccsd
      procedure :: jacobian_ccsd_j2                            => jacobian_ccsd_j2_ccsd
      procedure :: jacobian_ccsd_k2                            => jacobian_ccsd_k2_ccsd
!
      procedure :: prepare_for_jacobian                        => prepare_for_jacobian_ccsd
!
      procedure :: save_jacobian_c2_intermediates              => save_jacobian_c2_intermediates_ccsd
      procedure :: save_jacobian_d2_intermediate               => save_jacobian_d2_intermediate_ccsd
      procedure :: save_jacobian_e2_intermediate               => save_jacobian_e2_intermediate_ccsd
      procedure :: save_jacobian_g2_intermediates              => save_jacobian_g2_intermediates_ccsd
      procedure :: save_jacobian_h2_intermediates              => save_jacobian_h2_intermediates_ccsd
      procedure :: save_jacobian_j2_intermediate               => save_jacobian_j2_intermediate_ccsd
!
!     Routines related to Jacobian transpose transformation
!
      procedure :: jacobian_transpose_transform_trial_vector   => jacobian_transpose_transform_trial_vector_ccsd
      procedure :: jacobian_transpose_ccsd_transformation      => jacobian_transpose_ccsd_transformation_ccsd
!
      procedure :: jacobian_transpose_ccsd_d1                  => jacobian_transpose_ccsd_d1_ccsd
      procedure :: jacobian_transpose_ccsd_e1                  => jacobian_transpose_ccsd_e1_ccsd
      procedure :: jacobian_transpose_ccsd_f1                  => jacobian_transpose_ccsd_f1_ccsd
      procedure :: jacobian_transpose_ccsd_g1                  => jacobian_transpose_ccsd_g1_ccsd
!
      procedure :: jacobian_transpose_ccsd_b2                  => jacobian_transpose_ccsd_b2_ccsd
      procedure :: jacobian_transpose_ccsd_c2                  => jacobian_transpose_ccsd_c2_ccsd
      procedure :: jacobian_transpose_ccsd_d2                  => jacobian_transpose_ccsd_d2_ccsd
      procedure :: jacobian_transpose_ccsd_e2                  => jacobian_transpose_ccsd_e2_ccsd
      procedure :: jacobian_transpose_ccsd_f2                  => jacobian_transpose_ccsd_f2_ccsd
      procedure :: jacobian_transpose_ccsd_g2                  => jacobian_transpose_ccsd_g2_ccsd
      procedure :: jacobian_transpose_ccsd_h2                  => jacobian_transpose_ccsd_h2_ccsd
      procedure :: jacobian_transpose_ccsd_i2                  => jacobian_transpose_ccsd_i2_ccsd
!
      procedure :: get_gs_orbital_differences                  => get_gs_orbital_differences_ccsd
      procedure :: get_es_orbital_differences                  => get_gs_orbital_differences_ccsd
      procedure :: calculate_energy                            => calculate_energy_ccsd
!
      procedure :: construct_eta                               => construct_eta_ccsd
!
      procedure :: get_multipliers                             => get_multipliers_ccsd
      procedure :: set_multipliers                             => set_multipliers_ccsd
      procedure :: initialize_multipliers                      => initialize_multipliers_ccsd
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_ccsd
      procedure :: save_multipliers                            => save_multipliers_ccsd
      procedure :: read_multipliers                            => read_multipliers_ccsd
      procedure :: destruct_multipliers                        => destruct_multipliers_ccsd
!
      procedure, nopass :: need_g_abcd                         => need_g_abcd_ccsd
!
      procedure :: normalization_for_jacobian_debug            => normalization_for_jacobian_debug_ccsd
!
   end type ccsd
!
!
   interface
!
      include "file_handling_ccsd_interface.F90"
      include "initialize_destruct_ccsd_interface.F90"
      include "set_get_ccsd_interface.F90"
      include "omega_ccsd_interface.F90"
      include "multiplier_equation_ccsd_interface.F90"
      include "jacobian_ccsd_interface.F90"
      include "jacobian_transpose_ccsd_interface.F90"
      include "zop_ccsd_interface.F90"
      include "debug_jacobian_ccsd_interface.F90"
!
   end interface
!
!
   interface ccsd 
!
      procedure :: new_ccsd 
!
   end interface ccsd 
!
!
contains
!
!
   function new_ccsd(system) result(wf)
!!
!!    New CCSD
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(ccsd) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      wf%name_ = 'ccsd'
!
      wf%system => system
      wf%bath_orbital = .false.
      wf%cvs = .false.
!
      call wf%read_hf()
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%initialize_files()
!
      call wf%read_orbital_coefficients()
      call wf%read_orbital_energies()
!
      wf%bath_orbital = .false.
      wf%frozen_core = .false.
!
      call wf%read_settings()
!
      if (wf%bath_orbital) call wf%make_bath_orbital()
      if (wf%frozen_core) call wf%remove_core_orbitals()
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
      wf%n_t2 = (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
!
      call wf%write_cc_restart()
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
   end function new_ccsd
!
!
   logical function need_g_abcd_ccsd()
!!
!!    Need g_abcd 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Returns whether the vvvv-part of the ERI matrix 
!!    is used to calculate the ground and/or excited state 
!!    equations. If not, there is no need to compute the
!!    entire ERI matrix and store it in memory.
!!
      implicit none 
!
      need_g_abcd_ccsd = .true.
!
   end function need_g_abcd_ccsd
!
!
   subroutine set_initial_amplitudes_guess_ccsd(wf)
!!
!!    Set initial amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      call zero_array(wf%t1, wf%n_o*wf%n_v)
!
      call wf%set_t2_to_mp2_guess()
!
   end subroutine set_initial_amplitudes_guess_ccsd
!
!
   subroutine set_t2_to_mp2_guess_ccsd(wf)
!!
!!    Set t2 amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    t_aibj = - g_aibj/ε_aibj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer :: a, b, i, j, ai, bj, aibj
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%get_vovo(g_aibj)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i-1) + a
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
                     wf%t2(aibj) = g_aibj(a,i,b,j)/(wf%orbital_energies(i) + &
                                                    wf%orbital_energies(j) - &
                                                    wf%orbital_energies(a + wf%n_o) - &
                                                    wf%orbital_energies(b + wf%n_o))
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine set_t2_to_mp2_guess_ccsd
!
!
   subroutine get_gs_orbital_differences_ccsd(wf, orbital_differences, N)
!!
!!    Get orbital differences
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      integer, intent(in) :: N
      real(dp), dimension(N), intent(inout) :: orbital_differences
!
      integer :: a, i, ai, b, j, bj, aibj
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai) = wf%orbital_energies(a + wf%n_o) - wf%orbital_energies(i)
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
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v)) = wf%orbital_energies(a + wf%n_o) &
                                                                   - wf%orbital_energies(i) &
                                                                   +  wf%orbital_energies(b + wf%n_o) &
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
   end subroutine get_gs_orbital_differences_ccsd
!
!
   subroutine print_dominant_amplitudes_ccsd(wf)
!!
!!    Print dominant amplitudes
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      call wf%print_dominant_x1(wf%t1,'t')
      call wf%print_dominant_x2(wf%t2,'t')
!
   end subroutine print_dominant_amplitudes_ccsd
!
!
   subroutine print_dominant_x_amplitudes_ccsd(wf, x, tag)
!!
!!    Print dominant amplitudes  (TODO)
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: x
!
      character(len=1) :: tag
!
      call wf%print_dominant_x1(x(1:wf%n_t1),tag)
      call wf%print_dominant_x2(x(wf%n_t1 + 1:wf%n_gs_amplitudes),tag)
!
   end subroutine print_dominant_x_amplitudes_ccsd
!
!
   subroutine print_dominant_x2_ccsd(wf, x2, tag)
!!
!!    Print dominant x2
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Prints the 20 most dominant double amplitudes,
!!    or sorts them if there are fewer than twenty of them.
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_t2), intent(in) :: x2
      character(len=1), intent(in)    :: tag
!
      real(dp), dimension(:), allocatable :: abs_x2
!
      integer, dimension(:), allocatable :: dominant_indices
      real(dp), dimension(:), allocatable     :: dominant_values
!
      integer :: n_elements, elm, i, a, j, b, ai, bj
!
!     Sort according to largest contributions
!
      call mem%alloc(abs_x2, wf%n_t2)
      abs_x2 = abs(x2)
!
      n_elements = 10
      if (n_elements .gt. wf%n_t2) n_elements = wf%n_t2
!
      call mem%alloc(dominant_indices, n_elements)
      call mem%alloc(dominant_values, n_elements)
!
      dominant_indices = 0
      call zero_array(dominant_values, n_elements)
      call get_n_highest(n_elements, wf%n_t2, abs_x2, dominant_values, dominant_indices)
!
!     Print largest contributions
!
      write(output%unit, '(/t6,a)') 'Largest double amplitudes:'
      write(output%unit, '(t6,a)')  '--------------------------------------------------'
      write(output%unit, '(t6,a)')  '   a      i       b      j         ' // tag // '(ai,bj)             '
      write(output%unit, '(t6,a)')  '--------------------------------------------------'
!
      do elm = 1, n_elements
!
         call invert_packed_index(dominant_indices(elm), ai, bj, (wf%n_o)*(wf%n_v))
         call invert_compound_index(ai, a, i, wf%n_v, wf%n_o)
         call invert_compound_index(bj, b, j, wf%n_v, wf%n_o)
!
         write(output%unit, '(t6,i4,4x,i3,4x,i4,4x,i3,3x,f19.12)') a, i, b, j, x2(dominant_indices(elm))
!
      enddo
!           '
      write(output%unit, '(t6,a)')  '--------------------------------------------------'
!
      call mem%dealloc(dominant_indices, n_elements)
      call mem%dealloc(dominant_values, n_elements)
      call mem%dealloc(abs_x2, wf%n_t2)
!
   end subroutine print_dominant_x2_ccsd
!
!
   subroutine from_biorthogonal_to_biorthonormal_ccsd(wf, X)
!!
!!    From biorthogonal to biorthonormal 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_t2), intent(inout) :: X 
!
      real(dp), dimension(:,:), allocatable :: X_unpacked
!
      integer :: I
!
      call mem%alloc(X_unpacked, wf%n_t1, wf%n_t1)
      call squareup(X, X_unpacked, wf%n_t1)
!
!$omp parallel do private(I)
      do I = 1, wf%n_t1 
!
         X_unpacked(I,I) = X_unpacked(I,I)/two
!
      enddo
!$omp end parallel do 
!
      call packin(X, X_unpacked, wf%n_t1)
      call mem%dealloc(X_unpacked, wf%n_t1, wf%n_t1)
!
   end subroutine from_biorthogonal_to_biorthonormal_ccsd
!
!
   subroutine form_newton_raphson_t_estimate_ccsd(wf, t, dt)
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
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: dt 
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: t 
!
      integer :: ai, aiai 
!
!     Change dt doubles diagonal to match the definition of the 
!     double amplitudes 
!
!$omp parallel do private(ai, aiai)
      do ai = 1, wf%n_t1
!
         aiai = ai*(ai - 3)/2 + 2*ai
         dt(wf%n_t1 + aiai) = two*dt(wf%n_t1 + aiai)
!
      enddo 
!$omp end parallel do 
!
!     Add the dt vector to the t vector 
!
      call daxpy(wf%n_gs_amplitudes, one, dt, 1, t, 1)    
!
   end subroutine form_newton_raphson_t_estimate_ccsd
!
!
end module ccsd_class
