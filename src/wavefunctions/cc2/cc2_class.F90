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
module cc2_class
!
!!
!!    Coupled cluster singles and perturbative doubles (CC2) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!
   use doubles_class
!
   implicit none
!
   type, extends(doubles) :: cc2
!
   contains
!
!     Ground state
!
      procedure :: construct_omega                             => construct_omega_cc2
      procedure :: calculate_energy                            => calculate_energy_cc2
!
!     Amplitudes and multipliers
!
      procedure :: construct_t2                                => construct_t2_cc2
      procedure :: construct_u_aibj                            => construct_u_aibj_cc2
!
      procedure :: construct_t2bar                             => construct_t2bar_cc2
!
!     Excited state 
!
      procedure :: get_es_orbital_differences                  => get_es_orbital_differences_cc2
!
!     Jacobian
!
      procedure :: prepare_for_jacobian                        => prepare_for_jacobian_cc2
!
      procedure :: jacobian_transformation                     => jacobian_transformation_cc2
      procedure :: jacobian_cc2_b2                             => jacobian_cc2_b2_cc2
!
!     Jacobian transpose
!
      procedure :: prepare_for_jacobian_transpose              => prepare_for_jacobian_transpose_cc2
!
      procedure :: jacobian_transpose_transformation           => jacobian_transpose_transformation_cc2
      procedure :: jacobian_transpose_cc2_b2                   => jacobian_transpose_cc2_b2_cc2
!
!     Multiplier equation
!
      procedure :: prepare_for_multiplier_equation             => prepare_for_multiplier_equation_cc2
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_cc2
!
!     Initialize and destruct
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_cc2 
      procedure :: destruct_amplitudes                         => destruct_amplitudes_cc2 
      procedure :: destruct_multipliers                        => destruct_multipliers_cc2 
!
!     Ground state density 
!
      procedure :: prepare_for_density                         => prepare_for_density_cc2
!
!     Debug 
!
      procedure :: omega_for_jacobian_debug                    => omega_for_jacobian_debug_cc2
      procedure :: amplitudes_for_jacobian_debug               => amplitudes_for_jacobian_debug_cc2
      procedure :: normalization_for_jacobian_debug            => normalization_for_jacobian_debug_cc2
      procedure :: construct_omega2                            => construct_omega2_cc2
!
!     Initialize wavefunction
!
      procedure :: initialize                                  => initialize_cc2
!
   end type cc2
!
   interface
!
      include "initialize_destruct_cc2_interface.F90"
      include "omega_cc2_interface.F90"
      include "multiplier_equation_cc2_interface.F90"
      include "jacobian_cc2_interface.F90"
      include "jacobian_transpose_cc2_interface.F90"
      include "zop_cc2_interface.F90"
      include "debug_jacobian_cc2_interface.F90"
!
   end interface
!
!
contains
!
!
   subroutine initialize_cc2(wf, template_wf)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'cc2'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_t2            = wf%n_t1*(wf%n_t1+1)/2
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1 + wf%n_t2
      wf%need_g_abcd     = .false.
!
      call wf%initialize_fock()
!
      call wf%print_amplitude_info()
!
   end subroutine initialize_cc2
!
!
   subroutine construct_t2_cc2(wf)
!!
!!    Construct t2 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    Construct
!!
!!       t_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, t_aibj
!
      integer :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)  
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)  
!
      call wf%eri%get_eri_t1('vovo', g_aibj)
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o 
         do b = 1, wf%n_v 
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  t_aibj(a, i, b, j) = (g_aibj(a, i, b, j))/ &
                                   (wf%orbital_energies(i) &
                                  + wf%orbital_energies(j) &
                                  - wf%orbital_energies(wf%n_o + a) &
                                  - wf%orbital_energies(wf%n_o + b))

!
               enddo
            enddo
         enddo 
      enddo
!$omp end parallel do
!
      call packin(wf%t2, t_aibj, wf%n_t1)    
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
   end subroutine construct_t2_cc2
!
!
   subroutine construct_t2bar_cc2(wf)
!!
!!    Construct t2bar
!!    Written by Sarai D. Folkestad, May, 2019
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: t2bar, g_iajb
!
      integer :: a, i, b, j
!
      call mem%alloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(t2bar, (wf%n_o*wf%n_v)**2)
!
!     t2bar = sum_ai tbar_ai A_ai,aibj
!
      call wf%jacobian_transpose_doubles_a2(t2bar, wf%t1bar)
      call symmetric_sum(t2bar, wf%n_t1)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov', g_iajb)
!
!     t2bar += η_aibj
!
      call add_2143_to_1234(four, g_iajb, t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-two, g_iajb, t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     t2bar = t2bar/(-ε_aibj)
!
!$omp parallel do private(a, b, i, j)
      do b = 1, wf%n_v
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  t2bar(a, i, b, j) = t2bar(a, i, b, j)/(- wf%orbital_energies(a + wf%n_o) &
                                                         -  wf%orbital_energies(b + wf%n_o) &
                                                         +  wf%orbital_energies(i) &
                                                         +  wf%orbital_energies(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call packin(wf%t2bar, t2bar, wf%n_t1)
!
      call mem%dealloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_t2bar_cc2
!
!
   subroutine construct_u_aibj_cc2(wf)
!!
!!    Construct u_aibj
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    Construct
!!
!!       u_aibj = 2t_aibj - t_ajbi
!!
!!    with
!!
!!       t_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      integer :: a, i, b, j
      real(dp) :: eps
!
      type(timings) :: timer
!
      timer = timings('Construct u_aibj CC2')
      call timer%turn_on()
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)  
!
      call wf%eri%get_eri_t1('vovo', g_aibj)
!
!$omp parallel do private(a, i, b, j)
      do j = 1, wf%n_o 
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  eps =   wf%orbital_energies(i)          &
                        + wf%orbital_energies(j)          &
                        - wf%orbital_energies(wf%n_o + a) &
                        - wf%orbital_energies(wf%n_o + b) 
!
                  wf%u_aibj(a, i, b, j) = (two*g_aibj(a, i, b, j) - g_aibj(a, j, b, i)) / eps
!
               enddo
            enddo
         enddo 
      enddo
!$omp end parallel do
!    
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
      call timer%turn_off()
!
   end subroutine construct_u_aibj_cc2
!
!
   subroutine get_es_orbital_differences_cc2(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(cc2), intent(in) :: wf
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
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v)) = wf%orbital_energies(a + wf%n_o)     &
                                                                     - wf%orbital_energies(i) &
                                                                     +  wf%orbital_energies(b + wf%n_o)  &
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
   end subroutine get_es_orbital_differences_cc2
!
!
end module cc2_class
