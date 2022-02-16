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
module cc_eta_xi_calculator_class
!
!!
!! CC eta and xi calculator class (abstract)
!! Written by Eirik F. Kjønstad, Jan 2022
!!
!
   use parameters
   use ccs_class, only: ccs
!
   implicit none
!
   type, abstract :: cc_eta_xi_calculator
!
      class(ccs), pointer :: wf
!
   contains
!
      procedure(calculate), public, deferred :: calculate
!
   end type cc_eta_xi_calculator
!
!
   abstract interface
!
      subroutine calculate(this, xiX, etaX)
!
         import :: cc_eta_xi_calculator, dp
!
         implicit none
!
         class(cc_eta_xi_calculator) :: this
!
         real(dp), dimension(this%wf%n_es_amplitudes, 3), intent(out) :: xiX
         real(dp), dimension(this%wf%n_es_amplitudes, 3), intent(out) :: etaX
!
      end subroutine calculate
!
   end interface
!
!
end module cc_eta_xi_calculator_class
