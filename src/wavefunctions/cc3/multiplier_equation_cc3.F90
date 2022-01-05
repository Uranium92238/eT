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
submodule (cc3_class) multiplier_equation_cc3
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
   module subroutine prepare_for_multiplier_equation_cc3(wf)
!!
!!    Prepare for Mutliplier equation
!!    Written by Alexander C. Paul, July 2019
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
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
   end subroutine prepare_for_multiplier_equation_cc3
!
!
   module subroutine construct_multiplier_equation_cc3(wf, equation)
!!
!!    Construct multiplier equation
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Adapted by Alexander C. Paul, June 2019
!!
!!    Constructs
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation
!
      real(dp), dimension(:), allocatable :: eta
!
!     Copy the multipliers, eq. = t-bar
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%get_multipliers(eta)
!
!     Transform the multipliers by A^T, eq. = t-bar^T A
!
!     Same as A^T transformation but with zero frequency and without cvs and "remove core"
      call wf%effective_jacobian_transpose_transformation(zero, eta, equation, &
                                                          cvs=.false., rm_core=.false.)
!
!     No triples contributions to η
!     Construct eta(CCSD) and add, eq. = t-bar^T A + eta
!
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_gs_amplitudes, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
   end subroutine construct_multiplier_equation_cc3
!
!
end submodule multiplier_equation_cc3
