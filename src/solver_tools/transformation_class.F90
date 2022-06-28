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
module transformation_class
!
!!
!!    Abstract transformation tool class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!!    Defines the interface to the transformation tool classes
!!
!!    These classes perform the transformation of a vector
!!    by some matrix, that is not necessarily stored in memory.
!!
!
   use kinds
!
   implicit none
!
   type, abstract :: transformation
!
      integer :: n_parameters
!
   contains
!
      procedure (transform_transformation),  deferred :: transform
      procedure (initialize_transformation), deferred :: initialize
!
   end type  transformation
!
   abstract interface
!
      subroutine transform_transformation(this, trial, transform)
!
         use parameters
!
         import transformation
!
         implicit none
!
         class(transformation), intent(in) :: this
         real(dp), dimension(this%n_parameters) :: trial, transform
!
      end subroutine
!
!
      subroutine initialize_transformation(this)
!
         use parameters
!
         import transformation
!
         implicit none
!
         class(transformation), intent(in) :: this
!
      end subroutine
   end interface
!
!
end module transformation_class
