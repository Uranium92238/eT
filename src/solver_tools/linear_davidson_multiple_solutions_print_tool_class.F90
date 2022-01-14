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
module linear_davidson_multiple_solutions_print_tool_class
!
!!
!!    Linear Davidson multiple solutions print tool class module
!!    Written by Regina Matveeva, Dec 2021
!!
!
   use linear_davidson_print_tool_class, only: linear_davidson_print_tool
!
   implicit none
!
   type, extends(linear_davidson_print_tool) :: linear_davidson_multiple_solutions_print_tool
!
   contains
!
      procedure, public, nopass :: print_iteration_header &
                                => print_iter_header_linear_davidson_multiple_solutions_print_tool
      procedure, public, nopass :: print_iteration &
                                => print_iteration_linear_davidson_multiple_solutions_print_tool
!
   end type linear_davidson_multiple_solutions_print_tool
!
   contains
!
   subroutine print_iter_header_linear_davidson_multiple_solutions_print_tool()
!!
!!    Print iteration header
!!    Written by Regina Matveeva, Dec 2021
!!
      implicit none
!
   end subroutine print_iter_header_linear_davidson_multiple_solutions_print_tool
!
!
   subroutine print_iteration_linear_davidson_multiple_solutions_print_tool(n_solutions,   &
                                                                            residual_norm, &
                                                                            iteration,     &
                                                                            reduced_dimension)
!!
!!    Print iteration
!!    Written by Regina Matveeva, Sept 2021
!!
      use parameters
      use global_out, only: output
!
      implicit none
!
      integer, intent(in) :: n_solutions, iteration, reduced_dimension
      real(dp), dimension(n_solutions), intent(in) :: residual_norm
!
      integer :: n
!
      call output%printf('n', 'Iteration:               (i4)', &
                        ints=[iteration], fs='(/t3,a)')
      call output%printf('n', 'Reduced space dimension: (i4)', &
                        ints=[reduced_dimension])
!
      call output%printf('n', ' Solution       Residual norm', fs='(/t3,a)')
      call output%print_separator('n', 29,'-')
!
      do n = 1, n_solutions
!
         call output%printf('n', ' (i3)              (e11.4)', &
                            ints=[n], reals=[residual_norm(n)])
!
      enddo
!
      call output%print_separator('n', 29,'-')
!
   end subroutine print_iteration_linear_davidson_multiple_solutions_print_tool
!
end module linear_davidson_multiple_solutions_print_tool_class
