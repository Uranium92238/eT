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
module linear_davidson_single_solution_print_tool_class
!
!!
!!    Linear Davidson single solution print tool class module
!!    Written by Regina Matveeva, Dec 2021
!!
!
   use linear_davidson_print_tool_class, only: linear_davidson_print_tool
   use global_out, only: output
!
   implicit none
!
   type, extends(linear_davidson_print_tool) :: linear_davidson_single_solution_print_tool
!
   contains
!
      procedure, public, nopass :: print_iteration_header &
                                => print_iter_header_linear_davidson_single_solution_print_tool
!
      procedure, public, nopass :: print_iteration &
                                => print_iteration_linear_davidson_single_solution_print_tool
!
      procedure, public, nopass :: print_summary &
                                => print_summary_linear_davidson_single_solution_print_tool
!
   end type linear_davidson_single_solution_print_tool
!
!
   contains
!
!
   subroutine print_iter_header_linear_davidson_single_solution_print_tool()
!!
!!    Print iteration header
!!    Written by Regina Matveeva, Dec 2021
!!
      implicit none
!
      call output%printf('n', ' Iteration       Residual norm', fs='(/t3,a)')
      call output%print_separator('n', 31, '-')
!
   end subroutine print_iter_header_linear_davidson_single_solution_print_tool
!
!
   subroutine print_iteration_linear_davidson_single_solution_print_tool(n_solutions,   &
                                                                         residual_norm, &
                                                                         iteration,     &
                                                                         reduced_dimension)
!!
!!    Print iteration
!!    Written by Regina Matveeva, Sept 2021
!!
      use parameters
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      integer, intent(in) :: n_solutions, iteration, reduced_dimension
      real(dp), dimension(n_solutions), intent(in) :: residual_norm
!
      call do_nothing(reduced_dimension)
!
      call output%printf('n', ' (i3)              (e11.4)', &
                            ints=[iteration], reals=[residual_norm(1)])
!
   end subroutine print_iteration_linear_davidson_single_solution_print_tool
!
!
   subroutine print_summary_linear_davidson_single_solution_print_tool(n_solutions, &
                                                                       converged,   &
                                                                       iteration)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      integer, intent(in) :: n_solutions, iteration
      logical, dimension(n_solutions), intent(in) :: converged
!
      call output%print_separator('n', 31, '-')
!
            if (.not. all(converged)) then
!
         call output%error_msg("Did not converge in the max number of &
                               &iterations in Davidson solver.")
!
      else
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
      endif
!
   end subroutine print_summary_linear_davidson_single_solution_print_tool
!
!
end module linear_davidson_single_solution_print_tool_class
