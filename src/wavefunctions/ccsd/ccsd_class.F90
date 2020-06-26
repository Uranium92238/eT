!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   use doubles_class
!
   implicit none
!
   type, extends(doubles) :: ccsd
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
      type(sequential_file) :: jacobian_transpose_d1_intermediate
      type(sequential_file) :: jacobian_transpose_e1_intermediate
      type(sequential_file) :: jacobian_transpose_f1_intermediate
      type(sequential_file) :: jacobian_transpose_g1_intermediate
      type(sequential_file) :: jacobian_transpose_d2_intermediate
      type(sequential_file) :: jacobian_transpose_e2_oo_intermediate
      type(sequential_file) :: jacobian_transpose_e2_vv_intermediate
      type(sequential_file) :: jacobian_transpose_f2_intermediate
      type(sequential_file) :: jacobian_transpose_g2_intermediate
      type(sequential_file) :: jacobian_transpose_g2_intermediate_2
      type(sequential_file) :: jacobian_transpose_i2_intermediate
!
   contains
!
!     Initialization/destruction procedures
!
      procedure :: initialize_amplitudes                      => initialize_amplitudes_ccsd
      procedure :: initialize_amplitudes_complex              => initialize_amplitudes_ccsd_complex
!
      procedure :: destruct_amplitudes                        => destruct_amplitudes_ccsd
      procedure :: destruct_amplitudes_complex                => destruct_amplitudes_ccsd_complex
!
      procedure :: initialize_multipliers                     => initialize_multipliers_ccsd
      procedure :: initialize_multipliers_complex             => initialize_multipliers_ccsd_complex
!
      procedure :: destruct_multipliers                       => destruct_multipliers_ccsd
      procedure :: destruct_multipliers_complex               => destruct_multipliers_ccsd_complex
!
!     Set/get procedures
!
      procedure :: set_amplitudes                             => set_amplitudes_ccsd
      procedure :: set_amplitudes_complex                     => set_amplitudes_ccsd_complex
!
      procedure :: get_amplitudes                             => get_amplitudes_ccsd
      procedure :: get_amplitudes_complex                     => get_amplitudes_ccsd_complex
!
      procedure :: set_multipliers                            => set_multipliers_ccsd
      procedure :: set_multipliers_complex                    => set_multipliers_ccsd_complex
!
      procedure :: get_multipliers                            => get_multipliers_ccsd
      procedure :: get_multipliers_complex                    => get_multipliers_ccsd_complex
!
!     Procedures related to omega
!
      procedure :: construct_omega                            => construct_omega_ccsd
      procedure :: construct_omega_complex                    => construct_omega_ccsd_complex
!
      procedure :: omega_ccsd_a2                              => omega_ccsd_a2_ccsd
      procedure :: omega_ccsd_a2_complex                      => omega_ccsd_a2_ccsd_complex
!
      procedure :: omega_ccsd_b2                              => omega_ccsd_b2_ccsd
      procedure :: omega_ccsd_b2_complex                      => omega_ccsd_b2_ccsd_complex
!
      procedure :: omega_ccsd_c2                              => omega_ccsd_c2_ccsd
      procedure :: omega_ccsd_c2_complex                      => omega_ccsd_c2_ccsd_complex
!
      procedure :: omega_ccsd_d2                              => omega_ccsd_d2_ccsd
      procedure :: omega_ccsd_d2_complex                      => omega_ccsd_d2_ccsd_complex
!
      procedure :: omega_ccsd_e2                              => omega_ccsd_e2_ccsd
      procedure :: omega_ccsd_e2_complex                      => omega_ccsd_e2_ccsd_complex
!
      procedure :: construct_u_aibj                           => construct_u_aibj_ccsd
      procedure :: construct_u_aibj_complex                   => construct_u_aibj_ccsd_complex
!
!     Procedures related to Jacobian transformation
!
      procedure :: jacobian_transformation                    => jacobian_transformation_ccsd
!
      procedure :: jacobian_ccsd_b2                           => jacobian_ccsd_b2_ccsd
      procedure :: jacobian_ccsd_c2                           => jacobian_ccsd_c2_ccsd
      procedure :: jacobian_ccsd_d2                           => jacobian_ccsd_d2_ccsd
      procedure :: jacobian_ccsd_e2                           => jacobian_ccsd_e2_ccsd
      procedure :: jacobian_ccsd_f2                           => jacobian_ccsd_f2_ccsd
      procedure :: jacobian_ccsd_g2                           => jacobian_ccsd_g2_ccsd
      procedure :: jacobian_ccsd_h2                           => jacobian_ccsd_h2_ccsd
      procedure :: jacobian_ccsd_i2                           => jacobian_ccsd_i2_ccsd
      procedure :: jacobian_ccsd_j2                           => jacobian_ccsd_j2_ccsd
      procedure :: jacobian_ccsd_k2                           => jacobian_ccsd_k2_ccsd
!
      procedure :: prepare_for_jacobian                       => prepare_for_jacobian_ccsd
!
      procedure :: save_jacobian_c2_intermediates             => save_jacobian_c2_intermediates_ccsd
      procedure :: save_jacobian_d2_intermediate              => save_jacobian_d2_intermediate_ccsd
      procedure :: save_jacobian_e2_intermediate              => save_jacobian_e2_intermediate_ccsd
      procedure :: save_jacobian_g2_intermediates             => save_jacobian_g2_intermediates_ccsd
      procedure :: save_jacobian_h2_intermediates             => save_jacobian_h2_intermediates_ccsd
      procedure :: save_jacobian_j2_intermediate              => save_jacobian_j2_intermediate_ccsd
!
!     Procedures related to Jacobian transpose transformation
!
      procedure :: prepare_for_jacobian_transpose             => prepare_for_jacobian_transpose_ccsd
      procedure :: prepare_for_jacobian_transpose_complex     => prepare_for_jacobian_transpose_ccsd_complex
!
      procedure :: jacobian_transpose_transformation          => jacobian_transpose_transformation_ccsd
      procedure :: jacobian_transpose_transformation_complex  => jacobian_transpose_transformation_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_d1                 => jacobian_transpose_ccsd_d1_ccsd
      procedure :: jacobian_transpose_ccsd_d1_complex         => jacobian_transpose_ccsd_d1_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_e1                 => jacobian_transpose_ccsd_e1_ccsd
      procedure :: jacobian_transpose_ccsd_e1_complex         => jacobian_transpose_ccsd_e1_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_f1                 => jacobian_transpose_ccsd_f1_ccsd
      procedure :: jacobian_transpose_ccsd_f1_complex         => jacobian_transpose_ccsd_f1_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_g1                 => jacobian_transpose_ccsd_g1_ccsd
      procedure :: jacobian_transpose_ccsd_g1_complex         => jacobian_transpose_ccsd_g1_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_b2                 => jacobian_transpose_ccsd_b2_ccsd
      procedure :: jacobian_transpose_ccsd_b2_complex         => jacobian_transpose_ccsd_b2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_c2                 => jacobian_transpose_ccsd_c2_ccsd
      procedure :: jacobian_transpose_ccsd_c2_complex         => jacobian_transpose_ccsd_c2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_d2                 => jacobian_transpose_ccsd_d2_ccsd
      procedure :: jacobian_transpose_ccsd_d2_complex         => jacobian_transpose_ccsd_d2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_e2                 => jacobian_transpose_ccsd_e2_ccsd
      procedure :: jacobian_transpose_ccsd_e2_complex         => jacobian_transpose_ccsd_e2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_f2                 => jacobian_transpose_ccsd_f2_ccsd
      procedure :: jacobian_transpose_ccsd_f2_complex         => jacobian_transpose_ccsd_f2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_g2                 => jacobian_transpose_ccsd_g2_ccsd
      procedure :: jacobian_transpose_ccsd_g2_complex         => jacobian_transpose_ccsd_g2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_h2                 => jacobian_transpose_ccsd_h2_ccsd
      procedure :: jacobian_transpose_ccsd_h2_complex         => jacobian_transpose_ccsd_h2_ccsd_complex
!
      procedure :: jacobian_transpose_ccsd_i2                 => jacobian_transpose_ccsd_i2_ccsd
      procedure :: jacobian_transpose_ccsd_i2_complex         => jacobian_transpose_ccsd_i2_ccsd_complex
!
      procedure :: save_jacobian_transpose_d1_intermediates           => save_jacobian_transpose_d1_intermediates_ccsd
      procedure :: save_jacobian_transpose_d1_intermediates_complex   => save_jacobian_transpose_d1_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_e1_intermediates           => save_jacobian_transpose_e1_intermediates_ccsd
      procedure :: save_jacobian_transpose_e1_intermediates_complex   => save_jacobian_transpose_e1_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_f1_intermediates           => save_jacobian_transpose_f1_intermediates_ccsd
      procedure :: save_jacobian_transpose_f1_intermediates_complex   => save_jacobian_transpose_f1_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_g1_intermediates           => save_jacobian_transpose_g1_intermediates_ccsd
      procedure :: save_jacobian_transpose_g1_intermediates_complex   => save_jacobian_transpose_g1_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_d2_intermediates           => save_jacobian_transpose_d2_intermediates_ccsd
      procedure :: save_jacobian_transpose_d2_intermediates_complex   => save_jacobian_transpose_d2_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_e2_oo_intermediate         => save_jacobian_transpose_e2_oo_intermediate_ccsd
      procedure :: save_jacobian_transpose_e2_oo_intermediate_complex => save_jacobian_transpose_e2_oo_intermediate_ccsd_complex
!
      procedure :: save_jacobian_transpose_e2_vv_intermediate         => save_jacobian_transpose_e2_vv_intermediate_ccsd
      procedure :: save_jacobian_transpose_e2_vv_intermediate_complex => save_jacobian_transpose_e2_vv_intermediate_ccsd_complex
!
      procedure :: save_jacobian_transpose_f2_intermediates           => save_jacobian_transpose_f2_intermediates_ccsd
      procedure :: save_jacobian_transpose_f2_intermediates_complex   => save_jacobian_transpose_f2_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_g2_intermediates           => save_jacobian_transpose_g2_intermediates_ccsd
      procedure :: save_jacobian_transpose_g2_intermediates_complex   => save_jacobian_transpose_g2_intermediates_ccsd_complex
!
      procedure :: save_jacobian_transpose_i2_intermediates           => save_jacobian_transpose_i2_intermediates_ccsd
      procedure :: save_jacobian_transpose_i2_intermediates_complex   => save_jacobian_transpose_i2_intermediates_ccsd_complex
!
!     Procedures related to multiplier equation
!
      procedure :: prepare_for_multiplier_equation            => prepare_for_multiplier_equation_ccsd
      procedure :: prepare_for_multiplier_equation_complex    => prepare_for_multiplier_equation_ccsd_complex
!
      procedure :: construct_eta                              => construct_eta_ccsd
      procedure :: construct_eta_complex                      => construct_eta_ccsd_complex
!
!     Other procedures
!
      procedure :: set_initial_amplitudes_guess               => set_initial_amplitudes_guess_ccsd
      procedure :: set_t2_to_mp2_guess                        => set_t2_to_mp2_guess_ccsd
!
      procedure :: read_amplitudes                            => read_amplitudes_ccsd
      procedure :: save_amplitudes                            => save_amplitudes_ccsd
!
      procedure :: save_multipliers                           => save_multipliers_ccsd
      procedure :: read_multipliers                           => read_multipliers_ccsd
!
      procedure :: print_dominant_x2                          => print_dominant_x2_ccsd
      procedure :: print_dominant_amplitudes                  => print_dominant_amplitudes_ccsd
      procedure :: print_dominant_x_amplitudes                => print_dominant_x_amplitudes_ccsd
!
      procedure :: form_newton_raphson_t_estimate             => form_newton_raphson_t_estimate_ccsd
!
      procedure :: normalization_for_jacobian_debug           => normalization_for_jacobian_debug_ccsd
!
      procedure :: get_gs_orbital_differences                 => get_gs_orbital_differences_ccsd
!
      procedure :: get_es_orbital_differences                 => get_gs_orbital_differences_ccsd
!
      procedure :: calculate_energy                           => calculate_energy_ccsd
      procedure :: calculate_energy_complex                   => calculate_energy_ccsd_complex
!
!     Procedures related to time dependency
!
      procedure :: make_complex                               => make_complex_ccsd
      procedure :: cleanup_complex                            => cleanup_complex_ccsd
!
      procedure :: construct_complex_time_derivative_amplitudes       => construct_complex_time_derivative_amplitudes_ccsd
      procedure :: construct_complex_time_derivative_multipliers      => construct_complex_time_derivative_multipliers_ccsd
!
!     Initialize wavefunction
!
      procedure :: initialize  => initialize_ccsd
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
      include "complex_ccsd_interface.F90"
!
      include "autogenerated_complex_files/initialize_destruct_ccsd_interface_complex.F90"
      include "autogenerated_complex_files/jacobian_transpose_ccsd_interface_complex.F90"
      include "autogenerated_complex_files/multiplier_equation_ccsd_interface_complex.F90"
      include "autogenerated_complex_files/omega_ccsd_interface_complex.F90"
      include "autogenerated_complex_files/set_get_ccsd_interface_complex.F90"
      include "autogenerated_complex_files/zop_ccsd_interface_complex.F90"
!
   end interface
!
!
contains
!
!
   subroutine initialize_ccsd(wf, template_wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'ccsd'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%need_g_abcd = .true.
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
      wf%n_t2 = (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
      wf%need_g_abcd     = .true.
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
   end subroutine initialize_ccsd
!
!
   subroutine set_initial_amplitudes_guess_ccsd(wf)
!!
!!    Set initial amplitudes guess
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Sets the cluster amplitudes to the MP2 guess, i.e. 
!!
!!       t_ai   = zero 
!!       t_aibj = - g_aibj / epsilon_aibj
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
!!    Sets the doubles amplitudes to the MP2 guess:
!!
!!       t_aibj = - g_aibj/ε_aibj
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp) :: eps_ai
!
      integer :: a, b, i, j, ai, bj, aibj, b_end
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%integrals%construct_g_pqrs_mo(g_aibj,                  &
                                            wf%n_o+1, wf%n_o+wf%n_v, &
                                            1, wf%n_o,               &
                                            wf%n_o+1, wf%n_o+wf%n_v, &
                                            1, wf%n_o)
!
!$omp parallel do schedule(guided) collapse(2) &
!$omp private(a, i, b, j, ai, bj, aibj, b_end, eps_ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i-1) + a
            eps_ai = wf%orbital_energies(i) - wf%orbital_energies(a + wf%n_o)
!
            do j = 1, i
!
               if (j .ne. i) then
                  b_end = wf%n_v
               else
                  b_end = a
               end if
!
               do b = 1, b_end
!
                  bj = wf%n_v*(j-1) + b
!
                  aibj = (ai*(ai-3)/2) + ai + bj
!
                  wf%t2(aibj) = g_aibj(b,j,a,i)/(eps_ai +                     &
                                                 wf%orbital_energies(j) -     &
                                                 wf%orbital_energies(b + wf%n_o))
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
!!    Sets the (ground state) orbital differences vector:
!!
!!       epsilon_ai     = epsilon_a - epsilon_i 
!!       epsilon_aibj   = epsilon_a + epsilon_b - epsilon_i - epsilon_j 
!!
!!    Here, the epsilon vector has dimensionality N = n_t1 + n_t2, i.e. 
!!    the doubles part of the vector is packed. 
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
!!    Prints the dominant amplitudes in the cluster amplitude vector.
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
!!    Print dominant x amplitudes
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
!!    Prints the dominant amplitudes in the amplitude vector x.
!!
!!    tag specified the printed label for the vector, e.g. tag = "t" for 
!!    the cluster amplitudes. 
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: x
!
      character(len=*) :: tag
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
!!    tag specified the printed label for the vector, e.g. tag = "t" for 
!!    the cluster amplitudes. 
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_t2), intent(in) :: x2
      character(len=*), intent(in)    :: tag
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
      call output%printf('m', 'Largest double amplitudes:', fs='(/t6,a)') 
      call output%print_separator('m', 50, '-', fs='(t6,a)')
      call output%printf('m', 'a      i       b      j         (a0)(ai,bj)', &
                         fs='(t9,a)', chars=[tag])
      call output%print_separator('m', 50, '-', fs='(t6,a)')
!
      do elm = 1, n_elements
!
         call invert_packed_index(dominant_indices(elm), ai, bj, (wf%n_o)*(wf%n_v))
         call invert_compound_index(ai, a, i, wf%n_v, wf%n_o)
         call invert_compound_index(bj, b, j, wf%n_v, wf%n_o)
!
         call output%printf('m', '(i4)    (i3)    (i4)    (i3)   (f19.12)', &
                            ints=[a, i, b, j], &
                            reals=[x2(dominant_indices(elm))], fs='(t6,a)')
!
      enddo
!           '
      call output%print_separator('m', 50, '-', fs='(t6,a)')
!
      call mem%dealloc(dominant_indices, n_elements)
      call mem%dealloc(dominant_values, n_elements)
      call mem%dealloc(abs_x2, wf%n_t2)
!
   end subroutine print_dominant_x2_ccsd
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
