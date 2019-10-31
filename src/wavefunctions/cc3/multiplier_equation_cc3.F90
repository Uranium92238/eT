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
submodule (cc3_class) multiplier_equation_cc3
!
!!
!!    Multiplier equation (CC3)
!!    Set up by Andreas Skeidsvoll, Aug 2019
!!
!!    Equation used for the construction of CC3 multipliers.
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
      prep_timer = timings("Time preparing for multiplier equation")
      call prep_timer%turn_on()
!
      call output%printf('Preparing for (a0) multiplier equations',   &
                          chars=[trim(wf%name_)], pl='v', fs='(/t3,a)')
!
      if (.not. wf%X_ajil%exists()) call wf%prep_cc3_jacobian_intermediates()
      if (.not. wf%g_cdlk_t%exists()) call wf%prep_cc3_jacobian_trans_integrals()
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
      class(cc3), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation
!
      real(dp), dimension(:), allocatable :: eta
!
!     Copy the multipliers, eq. = t-bar
!
      call wf%get_multipliers(equation)
!
!     Transform the multipliers by A^T, eq. = t-bar^T A
!
!     Same as A^T transformation but with zero frequency and cvs is switched off
      call wf%effective_jacobian_transpose_transformation(omega = zero, c = equation, cvs = .false.)
!
!     No triples contributions to η
!     Construct eta(CCSD) and add, eq. = t-bar^T A + eta
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
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
