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
module preconditioner_getter_class
!
!!
!!    Preconditioner getter class module
!!    Written by Sarai D. Folkestad, May 2021
!!
!
   use parameters
!
   implicit none
!
   type, abstract :: preconditioner_getter
!
      integer :: n_parameters
!
   contains
!
      procedure(get_preconditioner_getter), deferred, public :: get
!
   end type preconditioner_getter
!
!
   abstract interface
!
      subroutine get_preconditioner_getter(this, preconditioner)
!!
!!       Get
!!       Written by Sarai D. Folkestad, May 2021
!!
         use parameters
!
         import preconditioner_getter
!
         implicit none
!
         class(preconditioner_getter), intent(in) :: this
         real(dp), dimension(this%n_parameters), intent(out) :: preconditioner
!
      end subroutine
!
   end interface
!
end module preconditioner_getter_class
