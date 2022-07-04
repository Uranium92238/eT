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
module null_preconditioner_getter_class
!
!!
!!    Null preconditioner getter class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!
   use parameters
   use preconditioner_getter_class, only: preconditioner_getter
!
   implicit none
!
   type, extends(preconditioner_getter) :: null_preconditioner_getter
!
   contains
!
      procedure, public :: get => get_null_preconditioner_getter
!
   end type null_preconditioner_getter
!
   interface  null_preconditioner_getter
!
      procedure :: new_null_preconditioner_getter
!
   end interface  null_preconditioner_getter
!
contains
!
   function new_null_preconditioner_getter(n_parameters) result(this)
!!
!!    New null preconditioner getter
!!    Written by Sarai D. Folkestad, May 2021
!!
      implicit none
!
      type(null_preconditioner_getter) :: this
      integer, intent(in) :: n_parameters
!
      this%n_parameters = n_parameters
!
   end function new_null_preconditioner_getter
!
   subroutine get_null_preconditioner_getter(this, preconditioner)
!!
!!    Get
!!    Written by Sarai D. Folkestad, May 2021
!!
!
      use array_utilities, only: constant_array
!
      implicit none
!
      class(null_preconditioner_getter), intent(in) :: this
      real(dp), dimension(this%n_parameters), intent(out) :: preconditioner
!
      call constant_array(preconditioner, this%n_parameters, one)
!
   end subroutine get_null_preconditioner_getter
!
end module null_preconditioner_getter_class
