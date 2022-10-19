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
submodule (cc2_class) multiplier_equation_cc2
!
!!
!!    Multiplier equation
!!
!!    Routines for calculation of the multiplier equation,
!!
!!       t-bar^T A + eta = 0,
!!
!!    where t-bar is the multiplier vector, and
!!
!!       A_mu,nu = < mu | exp(-T) [H, τ_nu] exp(T) | R >
!!       eta_mu  = < R | exp(-T) [H, τ_mu] exp(T) | R >.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_multiplier_equation_cc2(wf)
!!
!!    Prepare for jacobian transpose transformation
!!    Written by Tor S. Haugland, Oct 2019
!!
!!    Based on prepare_for_multiplier_equation_cc3 by Alex Paul
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      type(timings) :: prep_timer
!
      prep_timer = timings("Prepare for multiplier equation", pl='normal')
      call prep_timer%turn_on()
!
      call output%printf('v', 'Preparing for (a0) multiplier equations', &
                         chars=[trim(wf%name_)], fs='(/t3,a)')
!
      call wf%prepare_for_jacobian_transpose()
!
      call prep_timer%turn_off()
!
   end subroutine prepare_for_multiplier_equation_cc2
!
!
   module subroutine construct_multiplier_equation_cc2(wf, equation)
!!
!!    Construct multiplier equation
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Constructs
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
!!    Solves analytically for tbar_aibj
!!
!!       tbar_aibj = - (η_aibj + sum_ai tbar_ai A_ai,aibj)/ε_aibj
!!
!!    where
!!
!!       η_aibj = 2 L_iajb
!!
!!    and uses this to set up 'equation'
!!
!!       η_ai + sum_bj tbar_bj A_bj,ai + sum_bjck tbar_bjck A_{bjck,ai}
!!
      use array_initialization, only: zero_array
      use reordering, only: symmetric_sum, add_2143_to_1234, add_2341_to_1234
!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation
!
      real(dp), dimension(:), allocatable :: eta
      real(dp), dimension(:,:,:,:), allocatable :: t2bar
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      integer :: a, b, i, j
!
!     Construct t2bar
!
      call mem%alloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o, set_zero=.true.)
!
!     t2bar = sum_ai tbar_ai A_ai,aibj
!
      call wf%jacobian_transpose_doubles_a2(t2bar, wf%t1bar)
      call symmetric_sum(t2bar, wf%n_t1)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri_t1%get('ovov', g_iajb)
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
!     Set up the multipliers equation
!
!     equation = sum_bj tbar_bj A_bj,ai
!
      call zero_array(equation, wf%n_gs_amplitudes)
!
      call wf%jacobian_transpose_ccs_a1(equation, wf%t1bar)
      call wf%jacobian_transpose_ccs_b1(equation, wf%t1bar)
      call wf%jacobian_transpose_doubles_a1(equation, wf%t1bar, wf%u_aibj)
!
!     equation += sum_bjck tbar_bjck A_{bjck,ai}
!
      call wf%jacobian_transpose_doubles_b1(equation, t2bar)
!
      call mem%dealloc(t2bar, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Add eta, equation = t-bar^T A + eta
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_gs_amplitudes, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
   end subroutine construct_multiplier_equation_cc2
!
!
end submodule multiplier_equation_cc2
