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
module linear_storage_tool_class
!
!!
!!    Abstract linear storage tool class module
!!    Written by Regina Matveeva, Sep 2021
!!
!!    Based on eigen_storage_tool class by Sarai D. Folkestad
!!
!!    Handles the storage of the solutions to an eigenvalue problem.
!!    Stores the eigenvectors and eigenvalues.
!!
!
   use parameters
!
   implicit none
!
   type, abstract :: linear_storage_tool
!
      integer :: n_parameters
!
   contains
!
      procedure (store_linear_storage_tool), deferred, public :: store
!
end type  linear_storage_tool
!
   abstract interface
!
      subroutine store_linear_storage_tool(this, solution_vector, n)
!
         use parameters
!
         import linear_storage_tool
!
         implicit none
!
         class(linear_storage_tool),               intent(inout) :: this
         integer,                                  intent(in) :: n
         real(dp), dimension(this%n_parameters),   intent(in) :: solution_vector
!
      end subroutine
!
   end interface
!
end module linear_storage_tool_class
