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
module abstract_solver_class
!
!!
!!    Abstract solver class module
!!    Written by Eirik F. Kjønstad, 2021
!!
!
   implicit none
!
   type, abstract :: abstract_solver
!
   contains
!
      procedure(run), deferred, public     :: run
!
      procedure, public :: cleanup &
                        => cleanup_abstract_solver
!
   end type abstract_solver
!
!
   abstract interface
!
      subroutine run(this)
!
         import :: abstract_solver
!
         implicit none
!
         class(abstract_solver), intent(inout) :: this
!
      end subroutine run
!
   end interface
!
!
contains
!
!
   subroutine cleanup_abstract_solver(this)
!!
!!    Cleanup
!!    Written by Eirik F. Kjønstad, 2022
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(abstract_solver), intent(inout) :: this
!
      call do_nothing(this)
!
   end subroutine cleanup_abstract_solver
!
!
end module abstract_solver_class
