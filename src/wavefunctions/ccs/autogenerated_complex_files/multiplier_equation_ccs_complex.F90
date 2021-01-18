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
submodule (ccs_class) multiplier_equation_ccs_complex
!
!!
!!    Multiplier equation submodule
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
   module subroutine prepare_for_multiplier_equation_ccs_complex(wf)
!!
!!    Prepare for the construction of the multipliers
!!    Written by Alexander C. Paul, July 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
!     For now, do nothing.
!
      call output%printf('v', '- No preparations for the ' // trim(wf%name_) // &
                         ' multiplier equation.', fs='(/t3,a)')
!
   end subroutine prepare_for_multiplier_equation_ccs_complex
!
!
   module subroutine construct_multiplier_equation_ccs_complex(wf, equation)
!!
!!    Construct multiplier equation
!!    Written by Eirik F. Kjønstad, Oct 2018
!!
!!    Constructs
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation
!
      complex(dp), dimension(:), allocatable :: eta
!
!     Copy the multipliers, eq. = t-bar
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%get_multipliers_complex(eta)
!
!     Transform the multipliers by A^T, eq. = t-bar^T A
!
      call wf%jacobian_transpose_transformation_complex(eta, equation)
!
!     Add eta, eq. = t-bar^T A + eta
!
      call wf%construct_eta_complex(eta)
!
      call zaxpy(wf%n_gs_amplitudes, one_complex, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
   end subroutine construct_multiplier_equation_ccs_complex
!
!
   module subroutine construct_eta_ccs_complex(wf, eta)
!!
!!    Construct eta
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: eta
!
      integer :: i, a, ai
!
!$omp parallel do private(a, i, ai)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = (wf%n_v)*(i - 1) + a
            eta(ai) = two_complex*(wf%fock_ia_complex(i, a))
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_eta_ccs_complex
!
!
end submodule multiplier_equation_ccs_complex
