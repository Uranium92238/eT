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
module file_linear_equation_storage_tool_class
!
!!
!!    File linear equation storage tool
!!    Written by Eirik Kjønstad, Oct 2022
!!
!
   use parameters
   use stream_file_class, only: stream_file
   use linear_equation_storage_tool_class, only: linear_equation_storage_tool
!
   implicit none
!
   type, extends(linear_equation_storage_tool) :: file_linear_equation_storage_tool
!
      integer, private :: n_solutions 
      class(stream_file), dimension(:), pointer, private :: files
!
   contains
!
      procedure, public :: store => store_file_linear_equation_storage_tool
!
   end type  file_linear_equation_storage_tool
!
!
   interface  file_linear_equation_storage_tool
!
      procedure :: new_file_linear_equation_storage_tool
!
   end interface  file_linear_equation_storage_tool
!
!
contains
!
   function new_file_linear_equation_storage_tool(n_parameters, n_solutions, files) result(this)
!!
!!    File linear equation storage tool
!!    Written by Eirik Kjønstad, Oct 2022
!!
      implicit none
!
      type(file_linear_equation_storage_tool) :: this
!
      integer, intent(in) :: n_parameters, n_solutions
!
      class(stream_file), dimension(n_solutions), intent(in), target :: files 
!
      this%n_parameters = n_parameters
      this%n_solutions = n_solutions
      this%files => files 
!
   end function new_file_linear_equation_storage_tool
!
!
   subroutine store_file_linear_equation_storage_tool(this, solution_vector, n)
!!
!!    Store
!!    Written by Eirik Kjønstad, Oct 2022
!!
      implicit none
!
      class(file_linear_equation_storage_tool),   intent(inout) :: this
      integer,                                    intent(in) :: n
      real(dp), dimension(this%n_parameters),     intent(in) :: solution_vector
!
      call this%files(n)%open_('rewind')
      call this%files(n)%write_(solution_vector, this%n_parameters)
      call this%files(n)%close_()
!
   end subroutine store_file_linear_equation_storage_tool
!
!
end module file_linear_equation_storage_tool_class
