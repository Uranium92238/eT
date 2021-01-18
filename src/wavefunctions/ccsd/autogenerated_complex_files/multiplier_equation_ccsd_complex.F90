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
submodule (ccsd_class) multiplier_equation_ccsd_complex
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
   module subroutine prepare_for_multiplier_equation_ccsd_complex(wf)
!!
!!    Prepare for multiplier equation
!!    Written by Tor S. Haugland, 2019
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      call output%printf('v', '- Prepare for multiplier equation')
      call wf%prepare_for_jacobian_transpose_complex()
!
   end subroutine prepare_for_multiplier_equation_ccsd_complex
!
!
   module subroutine construct_eta_ccsd_complex(wf, eta)
!!
!!    Construct eta (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Note: the routine assumes that eta is initialized and that the Fock matrix
!!    has been constructed.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: eta
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_iajb
      complex(dp), dimension(:,:,:,:), allocatable :: eta_aibj
!
      integer :: i, a, ai
!
      call zero_array_complex(eta, wf%n_gs_amplitudes)
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
            eta(ai) = two_complex*(wf%fock_ia_complex(i, a)) ! eta_ai = 2 F_ia
!
         enddo
      enddo
!$omp end parallel do
!
!     eta_ai_bj = 2* L_iajb = 4 * g_iajb(i,a,j,b) - 2 * g_iajb(i,b,j,a)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri_complex%get_eri_t1('ovov', g_iajb)
!
      call mem%alloc(eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array_complex(eta_aibj, (wf%n_o*wf%n_v)**2)
!
      call add_2143_to_1234(four_complex, g_iajb, eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-two_complex, g_iajb, eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Pack vector into doubles eta
!
      call packin(eta(wf%n_t1 + 1 : wf%n_gs_amplitudes), eta_aibj, wf%n_t1)
!
      call mem%dealloc(eta_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_eta_ccsd_complex
!
!
end submodule multiplier_equation_ccsd_complex
