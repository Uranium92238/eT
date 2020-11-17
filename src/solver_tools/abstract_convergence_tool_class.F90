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
module abstract_convergence_tool_class
!
!!
!!    Abstract convergence tool class module
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!    Rolf H. Myhre and Alexander C. Paul, 2020
!! 
!!    Defines the interface to the convergence tool
!
   use kinds
!
   implicit none
!
   type, abstract :: abstract_convergence_tool
!
      real(dp) :: energy_threshold
      real(dp) :: residual_threshold
      logical  :: energy_convergence
!
   contains
!
      procedure (has_converged),  deferred :: has_converged
      procedure (print_settings), deferred :: print_settings
!
      procedure, non_overridable :: set_energy_threshold
      procedure, non_overridable :: set_residual_threshold
!
   end type abstract_convergence_tool
!
   abstract interface
!
      pure function has_converged(this, residual_norm, dE, iteration) result(converged)
!
         use kinds
         import abstract_convergence_tool
!
         implicit none
!
         class(abstract_convergence_tool), intent(in)           :: this
         real(dp),                         intent(in)           :: residual_norm
         real(dp),                         intent(in)           :: dE
         integer,                          intent(in), optional :: iteration
!
         logical :: converged
!
      end function has_converged
!
      subroutine print_settings(this)
!
         import abstract_convergence_tool
!
         implicit none
!
         class(abstract_convergence_tool), intent(in) :: this
         logical :: converged
!
      end subroutine print_settings
!
   end interface
!
contains
!
   subroutine set_energy_threshold(this, energy_threshold)
!!
!!    Set energy thresholds
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(abstract_convergence_tool), intent(inout) :: this
      real(dp),                         intent(in)    :: energy_threshold
!
      this%energy_threshold = energy_threshold
      this%energy_convergence = .true.
!
   end subroutine set_energy_threshold
!
   subroutine set_residual_threshold(this, residual_threshold)
!!
!!    Set residual thresholds
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(abstract_convergence_tool), intent(inout) :: this
      real(dp),                         intent(in)    :: residual_threshold
!
      this%residual_threshold = residual_threshold
!
   end subroutine set_residual_threshold
!
end module abstract_convergence_tool_class
