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
module convergence_tool_class
!
!!
!!    Convergence tool class module
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre and Alexander C. Paul, 2020
!
   use kinds
   use global_out,                      only: output
!
   implicit none
!
   type :: convergence_tool
!
      real(dp) :: energy_threshold
      real(dp) :: residual_threshold
      logical  :: energy_convergence
!
   contains
!
      procedure, public :: has_converged
      procedure, public :: print_settings
      procedure, public :: set_energy_threshold
      procedure, public :: set_residual_threshold
      procedure, public :: get_energy_threshold
      procedure, public :: get_residual_threshold

!
   end type convergence_tool
!
   interface convergence_tool
!
      procedure :: new_convergence_tool
!
   end interface convergence_tool
!
contains
!
   function new_convergence_tool(energy_threshold, residual_threshold, energy_convergence) result(this)
!!
!!    New convergence tool
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre and Alexander C. Paul, 2020
!!
      implicit none
!
      type(convergence_tool) :: this
!
      real(dp), intent(in) :: residual_threshold
      real(dp), intent(in) :: energy_threshold
      logical,  intent(in) :: energy_convergence
!
      this%residual_threshold = residual_threshold
      this%energy_threshold   = energy_threshold
      this%energy_convergence = energy_convergence
!
   end function new_convergence_tool
!
!
   pure function has_converged(this, residual_norm, dE, iteration) result(converged)
!!
!!    Has converged
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Rolf H. Myhre and Alexander C. Paul, 2020
!!
!!    Returns true if equations are converged according to the residual
!!    norm and, if requested, the energy change (dE)
!!
!!    If the iteration number is passed, convergence can be
!!    reached in the first iteration if the residual norm is
!!    below the threshold, even though dE > energy threshold.
!!
!
      implicit none
!
      class(convergence_tool),   intent(in)           :: this
      real(dp),                  intent(in)           :: residual_norm
      real(dp),                  intent(in)           :: dE
      integer,                   intent(in), optional :: iteration
!
      logical :: converged
!
      logical :: converged_energy
      logical :: converged_residual
!
      converged_energy   = (abs(dE) .lt. this%energy_threshold)
      converged_residual = residual_norm .lt. this%residual_threshold
!
      if (this%energy_convergence) then
!
         converged  = converged_energy .and. converged_residual
!
      else
!
         converged  = converged_residual
!
      endif
!
      if (present(iteration)) then
         if (converged_residual .and. iteration .eq. 1) converged = .true.
      endif
!
   end function has_converged
!
!
   subroutine print_settings(this)
!!
!!    Print
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2020
!!
      implicit none
!
      class(convergence_tool), intent(in) :: this
!
      call output%printf('m', '- Convergence thresholds', fs='(/t3,a)')
!
      call output%printf('m', 'Residual threshold:           (e11.4)', &
                         reals=[this%residual_threshold], fs='(/t6,a)')
!
      if (this%energy_convergence) &
         call output%printf('m', 'Energy threshold:             (e11.4)', &
                         reals=[this%energy_threshold], fs='(t6,a)')
!
   end subroutine print_settings
!
!
   subroutine get_energy_threshold(this, energy_threshold)
!!
!!    Get energy thresholds
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(convergence_tool), intent(inout) :: this
      real(dp),                intent(out)    :: energy_threshold
!
      energy_threshold = this%energy_threshold
!
   end subroutine get_energy_threshold
!
!
   subroutine get_residual_threshold(this, residual_threshold)
!!
!!    Get residual thresholds
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(convergence_tool), intent(inout) :: this
      real(dp),                intent(out)    :: residual_threshold
!
      residual_threshold = this%residual_threshold
!
   end subroutine get_residual_threshold
!
!
   subroutine set_energy_threshold(this, energy_threshold)
!!
!!    Set energy thresholds
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(convergence_tool), intent(inout) :: this
      real(dp),                intent(in)    :: energy_threshold
!
      this%energy_threshold = energy_threshold
      this%energy_convergence = .true.
!
   end subroutine set_energy_threshold
!
!
   subroutine set_residual_threshold(this, residual_threshold)
!!
!!    Set residual thresholds
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none
!
      class(convergence_tool), intent(inout) :: this
      real(dp),                intent(in)    :: residual_threshold
!
      this%residual_threshold = residual_threshold
!
   end subroutine set_residual_threshold
!
!
end module convergence_tool_class
