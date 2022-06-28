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
module frequency_dependent_transformation_class
!
!!
!!    Frequency-dependent transformation class
!!    Written by Sarai D. Folkestad, May 2022
!!
!
   use parameters
   use transformation_class,   only: transformation
!
   implicit none
!
   type, abstract, extends(transformation) :: frequency_dependent_transformation
!
      real(dp), private :: frequency
!
   contains
!
      procedure  :: set_frequency
      procedure  :: get_frequency
!
   end type  frequency_dependent_transformation
!
contains
!
   subroutine set_frequency(this, frequency)
!
      implicit none
!
      class(frequency_dependent_transformation), intent(inout) :: this
      real(dp), intent(in) :: frequency
!
      this%frequency = frequency
!
   end subroutine set_frequency
!
!
   function get_frequency(this) result(frequency)
!
      implicit none
!
      class(frequency_dependent_transformation), intent(in) :: this
      real(dp) :: frequency
!
      frequency = this%frequency
!
   end function get_frequency
!
!
end module frequency_dependent_transformation_class
