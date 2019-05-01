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
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: ccsd
!
      real(dp), dimension(:), allocatable :: t2
      real(dp), dimension(:), allocatable :: t2bar
!
      type(file) :: t2_file, t2bar_file
      type(file) :: r2_file, l2_file
!
      integer :: n_t2  
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: prepare                                     => prepare_ccsd
!
      procedure :: initialize_files                            => initialize_files_ccsd
      procedure :: initialize_doubles_files                    => initialize_doubles_files_ccsd
!
!     Routines related to the amplitudes
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_ccsd
      procedure :: initialize_t2                               => initialize_t2_ccsd
      procedure :: destruct_t2                                 => destruct_t2_ccsd
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
      procedure :: save_doubles_vector                         => save_doubles_vector_ccsd
      procedure :: read_doubles_vector                         => read_doubles_vector_ccsd
!
      procedure :: read_excited_state                          => read_excited_state_ccsd
!
      procedure :: save_excited_state                          => save_excited_state_ccsd 
!
!     Routines related to omega
!
      procedure :: construct_omega                             => construct_omega_ccsd
!
      procedure :: omega_ccsd_a1                               => omega_ccsd_a1_ccsd
      procedure :: omega_ccsd_b1                               => omega_ccsd_b1_ccsd
      procedure :: omega_ccsd_c1                               => omega_ccsd_c1_ccsd
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
      procedure :: jacobian_ccsd_a1                            => jacobian_ccsd_a1_ccsd
      procedure :: jacobian_ccsd_b1                            => jacobian_ccsd_b1_ccsd
      procedure :: jacobian_ccsd_c1                            => jacobian_ccsd_c1_ccsd
      procedure :: jacobian_ccsd_d1                            => jacobian_ccsd_d1_ccsd
!
      procedure :: jacobian_ccsd_a2                            => jacobian_ccsd_a2_ccsd
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
!     Routines related to Jacobian transpose transformation
!
      procedure :: jacobian_transpose_transform_trial_vector   => jacobian_transpose_transform_trial_vector_ccsd
      procedure :: jacobian_transpose_ccsd_transformation      => jacobian_transpose_ccsd_transformation_ccsd
!
      procedure :: jacobian_transpose_ccsd_a1                  => jacobian_transpose_ccsd_a1_ccsd
      procedure :: jacobian_transpose_ccsd_b1                  => jacobian_transpose_ccsd_b1_ccsd
      procedure :: jacobian_transpose_ccsd_c1                  => jacobian_transpose_ccsd_c1_ccsd
      procedure :: jacobian_transpose_ccsd_d1                  => jacobian_transpose_ccsd_d1_ccsd
      procedure :: jacobian_transpose_ccsd_e1                  => jacobian_transpose_ccsd_e1_ccsd
      procedure :: jacobian_transpose_ccsd_f1                  => jacobian_transpose_ccsd_f1_ccsd
      procedure :: jacobian_transpose_ccsd_g1                  => jacobian_transpose_ccsd_g1_ccsd
!
      procedure :: jacobian_transpose_ccsd_a2                  => jacobian_transpose_ccsd_a2_ccsd
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
      procedure :: initialize_t2bar                            => initialize_t2bar_ccsd
      procedure :: get_multipliers                             => get_multipliers_ccsd
      procedure :: set_multipliers                             => set_multipliers_ccsd
      procedure :: initialize_multipliers                      => initialize_multipliers_ccsd
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_ccsd
      procedure :: save_multipliers                            => save_multipliers_ccsd
      procedure :: read_multipliers                            => read_multipliers_ccsd
      procedure :: destruct_multipliers                        => destruct_multipliers_ccsd
      procedure :: destruct_t2bar                              => destruct_t2bar_ccsd
!
      procedure :: get_cvs_projector                           => get_cvs_projector_ccsd
!
!     Routines related to property calculations
!
      procedure :: construct_etaX                              => construct_etaX_ccsd  
      procedure :: etaX_ccsd_a1                                => etaX_ccsd_a1_ccsd    
      procedure :: etaX_ccsd_a2                                => etaX_ccsd_a2_ccsd    
      procedure :: etaX_ccsd_b2                                => etaX_ccsd_b2_ccsd    
!     
      procedure :: construct_csiX                              => construct_csiX_ccsd  
      procedure :: csiX_ccsd_a1                                => csiX_ccsd_a1_ccsd    
      procedure :: csiX_ccsd_a2                                => csiX_ccsd_a2_ccsd    
!
      procedure :: add_etaX_eom_correction                     => add_etaX_eom_correction_ccsd
      procedure :: etaX_eom_ccsd_a1                            => etaX_eom_ccsd_a1_ccsd
!
      procedure, nopass :: need_g_abcd                         => need_g_abcd_ccsd
!
!     One-electron density 
!
      procedure :: construct_density                           => construct_density_ccsd
!
      procedure :: one_el_density_ccsd_oo                      => one_el_density_ccsd_oo_ccsd
      procedure :: one_el_density_ccsd_vv                      => one_el_density_ccsd_vv_ccsd
      procedure :: one_el_density_ccsd_ov                      => one_el_density_ccsd_ov_ccsd
!
   end type ccsd
!
!
   interface
!
      include "fop_ccsd_interface.F90"
      include "files_ccsd_interface.F90"
      include "omega_ccsd_interface.F90"
      include "get_set_ccsd_interface.F90"
      include "jacobian_ccsd_interface.F90"
      include "jacobian_transpose_ccsd_interface.F90"
!
   end interface
!
!
contains
!
!
   subroutine prepare_ccsd(wf, system)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      type(file) :: hf_restart_file 
!
      wf%name_ = 'ccsd'
!
      wf%system => system
!
      call hf_restart_file%init('hf_restart_file', 'sequential', 'unformatted')
!
      call disk%open_file(hf_restart_file, 'read', 'rewind')
!
      read(hf_restart_file%unit) wf%n_ao 
      read(hf_restart_file%unit) wf%n_mo 
      read(hf_restart_file%unit) 
      read(hf_restart_file%unit) wf%n_o  
      read(hf_restart_file%unit) wf%n_v  
      read(hf_restart_file%unit) wf%hf_energy  
!
      call disk%close_file(hf_restart_file)
!
      wf%n_t1 = (wf%n_o)*(wf%n_v)
      wf%n_t2 = (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
!
      wf%n_gs_amplitudes = wf%n_t1 + wf%n_t2
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
!
      call wf%initialize_files()
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%read_orbital_coefficients()
      call wf%read_orbital_energies()
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
   end subroutine prepare_ccsd
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
   subroutine initialize_amplitudes_ccsd(wf)
!!
!!    Initialize amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Allocates the amplitudes. This routine must be overwritten in
!!    descendants which have more amplitudes.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%initialize_t1()
      call wf%initialize_t2()
!
   end subroutine initialize_amplitudes_ccsd
!
!
   subroutine initialize_t2_ccsd(wf)
!!
!!    Initialize t2 amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (.not. allocated(wf%t2)) call mem%alloc(wf%t2, wf%n_t2)
!
   end subroutine initialize_t2_ccsd
!
!
   subroutine destruct_t2_ccsd(wf)
!!
!!    Destruct t2 amplitudes
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2)) call mem%dealloc(wf%t2, wf%n_t2)
!
   end subroutine destruct_t2_ccsd
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
      wf%t1 = zero
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
   subroutine calculate_energy_ccsd(wf)
!!
!!     Calculate energy (CCSD)
!!     Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!     Andreas Skeidsvoll, 2018
!!
!!     Calculates the CCSD energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb ! g_iajb
!
      integer :: a = 0, i = 0, b = 0, j = 0, ai = 0
      integer :: bj = 0, aibj = 0
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iajb)
!
!     Set the initial value of the energy
!
      wf%energy = wf%hf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = (i-1)*wf%n_v + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  aibj = (max(ai,bj)*(max(ai,bj)-3)/2) + ai + bj
!
!                 Add the correlation energy
!
                  wf%energy = wf%energy +                                   &
                                 (wf%t2(aibj) + (wf%t1(a,i))*(wf%t1(b,j)))* &
                                 (two*g_iajb(i,a,j,b) - g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
   end subroutine calculate_energy_ccsd
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
   subroutine construct_eta_ccsd(wf, eta)
!!
!!    Construct eta (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Note: the routine assumes that eta is initialized and that the Fock matrix
!!    has been constructed.
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: eta
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
      real(dp), dimension(:,:,:,:), allocatable :: eta_aibj
!
      integer :: i = 0, a = 0, j = 0, b = 0, aibj = 0
      integer :: bj = 0, ai = 0
!
      eta = zero
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
            eta(ai) = two*(wf%fock_ia(i, a)) ! eta_ai = 2 F_ia
!
         enddo
      enddo
!
!     eta_ai_bj = 2* L_iajb = 4 * g_iajb(i,a,j,b) - 2 * g_iajb(i,b,j,a)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iajb)
!
      call mem%alloc(eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      eta_aibj = zero
!
      call add_2143_to_1234(four, g_iajb, eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-two, g_iajb, eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Pack vector into doubles eta
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  aibj = max(ai, bj)*(max(ai,bj)-3)/2 + ai + bj
!
                  eta(wf%n_t1 + aibj) = eta_aibj(a, i, b, j)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_eta_ccsd
!
!
   subroutine initialize_multipliers_ccsd(wf)
!!
!!    Initialize multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Allocates the multipliers. This routine must be overwritten in
!!    descendants which have more multipliers.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%initialize_t1bar()
      call wf%initialize_t2bar()
!
   end subroutine initialize_multipliers_ccsd
!
!
   subroutine initialize_t2bar_ccsd(wf)
!!
!!    Initialize T2-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (.not. allocated(wf%t2bar)) call mem%alloc(wf%t2bar, wf%n_t2)
!
   end subroutine initialize_t2bar_ccsd
!
!
   subroutine construct_multiplier_equation_ccsd(wf, equation)
!!
!!    Construct multiplier equation
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Constructs
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation
!
      real(dp), dimension(:), allocatable :: eta
!
!     Copy the multipliers, eq. = t-bar
!
      call dcopy(wf%n_t1, wf%t1bar, 1, equation, 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, equation(wf%n_t1 + 1), 1)
!
!     Transform the multipliers by A^T, eq. = t-bar^T A
!
      call wf%jacobian_transpose_ccsd_transformation(equation)
!
!     Add eta, eq. = t-bar^T A + eta
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_gs_amplitudes, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
   end subroutine construct_multiplier_equation_ccsd
!
!
   subroutine destruct_multipliers_ccsd(wf)
!!
!!    Destruct multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
!!    Deallocates the multipliers. This routine must be overwritten in
!!    descendants which have more multipliers.
!!
      implicit none
!
      class(ccsd) :: wf
!
      call wf%destruct_t1bar()
      call wf%destruct_t2bar()
!
   end subroutine destruct_multipliers_ccsd
!
!
   subroutine destruct_t2bar_ccsd(wf)
!!
!!    Destruct T2-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2bar)) call mem%dealloc(wf%t2bar, wf%n_gs_amplitudes)
!
   end subroutine destruct_t2bar_ccsd
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
      dominant_values  = zero
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
   subroutine get_cvs_projector_ccsd(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
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
   end subroutine get_cvs_projector_ccsd
!
!
   subroutine one_el_density_ccsd_oo_ccsd(wf)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad
!!
!!    D_ij -= sum_abk t_akb,i tbar_akb,j 
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_akbj, t_akbi
!
      call mem%alloc(tbar_akbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(t_akbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2bar, tbar_akbj, (wf%n_v)*(wf%n_o))
      call squareup(wf%t2, t_akbi, (wf%n_v)*(wf%n_o))
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_v**2)*(wf%n_o),   &
                  -one,                   &
                  t_akbi,                 & ! t_akb_i
                  (wf%n_v**2)*(wf%n_o),   &
                  tbar_akbj,              & ! tbar_akb_j
                  (wf%n_v**2)*(wf%n_o),   &
                  one,                    &
                  wf%density,             &
                  wf%n_mo)
!
      call mem%dealloc(tbar_akbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_akbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine one_el_density_ccsd_oo_ccsd
!
!
   subroutine one_el_density_ccsd_vv_ccsd(wf)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad
!!
!!    D_ab += sum_jci tbar_a,jci t_b,jci
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ajci, t_bjci
!
      call mem%alloc(tbar_ajci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(t_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(wf%t2bar, tbar_ajci, (wf%n_v)*(wf%n_o))
      call squareup(wf%t2, t_bjci, (wf%n_v)*(wf%n_o))
!
      call dgemm('N', 'T',                            &
                  wf%n_v,                             &
                  wf%n_v,                             &
                  (wf%n_o**2)*(wf%n_v),               &
                  one,                                &
                  tbar_ajci,                          & ! tbar_a_jci
                  wf%n_v,                             &
                  t_bjci,                             & ! t_b_jci
                  wf%n_v,                             &
                  one,                                &
                  wf%density(wf%n_o + 1, wf%n_o + 1), &
                  wf%n_mo)
!
      call mem%dealloc(tbar_ajci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine one_el_density_ccsd_vv_ccsd
!
!
   subroutine one_el_density_ccsd_ov_ccsd(wf)
!!
!!    One electron density vo
!!    Written by Sarai D. Folkestad
!!
!!    D_ia += sum_bj u^{ab}_ij tbar_bj = sum_bj u_ia,bj tbar_bj 
!!
!!    u^{ab}_ij = 2t_aibj - t_ajbi = u_iabj 
!!
      implicit none
!
      class(ccsd) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: u_iabj, t_aibj
      real(dp), dimension(:,:), allocatable :: D_ia
!
      integer :: i, a
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(u_iabj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      u_iabj = zero
      call add_2134_to_1234(two, t_aibj, u_iabj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call add_2431_to_1234(-one, t_aibj, u_iabj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(D_ia, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  u_iabj,              & ! u_ia_bj
                  (wf%n_o)*(wf%n_v),   &
                  wf%t1bar,            & ! tbar_bj
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  D_ia,                &
                  (wf%n_o)*(wf%n_v))

      call mem%dealloc(u_iabj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            wf%density(i, wf%n_o + a) = wf%density(i, wf%n_o + a) + D_ia(i, a)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_ia, wf%n_o, wf%n_v)
!
   end subroutine one_el_density_ccsd_ov_ccsd
!
!
   subroutine construct_density_ccsd(wf)
!!
!!    Construct one-electron density
!!    Written by Sarai Dery Folkestad
!!
!!    Constructs the one-electron density 
!!    matrix in the T1 basis
!!
      implicit none
!
      class(ccsd) :: wf
!
      wf%density = zero
!
      call wf%one_el_density_ccs_oo()
      call wf%one_el_density_ccs_vo()
!      
      call wf%one_el_density_ccsd_oo()
      call wf%one_el_density_ccsd_vv()
      call wf%one_el_density_ccsd_ov()
!
   end subroutine construct_density_ccsd
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
      do I = 1, wf%n_t1 
!
         X_unpacked(I,I) = X_unpacked(I,I)/two
!
      enddo 
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
      do ai = 1, wf%n_t1
!
         aiai = ai*(ai - 3)/2 + 2*ai
         dt(wf%n_t1 + aiai) = two*dt(wf%n_t1 + aiai)
!
      enddo  
!
!     Add the dt vector to the t vector 
!
      call daxpy(wf%n_gs_amplitudes, one, dt, 1, t, 1)    
!
   end subroutine form_newton_raphson_t_estimate_ccsd
!
!
end module ccsd_class
