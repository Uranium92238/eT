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
module function_class
!
!!
!!    Function
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Class that represents a function, where one can evaluate function values 
!!    (and values of derivatives) given set of specifiable parameters
!!
!
   use kinds
!
   implicit none
!
   type, abstract :: function_
!
      integer :: n_parameters
!
   contains 
!
      procedure(get_parameters), deferred    :: get_parameters
      procedure(set_parameters), deferred    :: set_parameters
      procedure(evaluate), deferred          :: evaluate 
      procedure(evaluate_gradient), deferred :: evaluate_gradient 
      procedure(initialize), deferred        :: initialize 
!
   end type function_
!
!
   abstract interface
!
      function get_parameters(this) result(x)
!
         import :: function_, dp 
!
         implicit none 
!
         class(function_) :: this 
!
         real(dp), dimension(this%n_parameters) :: x 
!
      end function get_parameters
!
      subroutine set_parameters(this, x)
!
         import :: function_, dp 
!
         implicit none 
!
         class(function_) :: this 
!
         real(dp), dimension(this%n_parameters), intent(in) :: x 
!
      end subroutine set_parameters
!
      function evaluate(this) result(energy)
!
         import :: function_, dp 
!
         implicit none 
!
         class(function_) :: this 
!
         real(dp) :: energy
!
      end function evaluate
!
      function evaluate_gradient(this) result(g)
!
         import :: function_, dp 
!
         implicit none 
!
         class(function_) :: this 
!
         real(dp), dimension(this%n_parameters) :: g
!
      end function evaluate_gradient
!
      subroutine initialize(this)
!
         import :: function_ 
!
         implicit none 
!
         class(function_) :: this 
!
      end subroutine initialize
!
   end interface 
!
!
end module function_class
