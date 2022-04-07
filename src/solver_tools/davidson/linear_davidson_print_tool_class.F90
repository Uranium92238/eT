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
module linear_davidson_print_tool_class
!
!!
!!    Linear Davidson print tool class module
!!    Written by Regina Matveeva, Dec 2021
!!
!
   use parameters
!
   implicit none
!
   type, abstract :: linear_davidson_print_tool
!
   contains
!
      procedure, public, nopass :: print_banner &
                                => print_banner_linear_davidson_print_tool
      procedure, public, nopass :: print_settings &
                                => print_settings_linear_davidson_print_tool
      procedure, public, nopass :: print_summary &
                                => print_summary_linear_davidson_print_tool
!
      procedure(print_iteration_header), deferred, public, nopass  :: print_iteration_header
      procedure(print_iteration), deferred, public, nopass  :: print_iteration
!
   end type linear_davidson_print_tool
!
   abstract interface
!
   subroutine print_iteration_header()
!!
!!    Print iteration header
!!    Written by Regina Matveeva, Sept 2021
!!
      implicit none
!
   end subroutine print_iteration_header
!
!
   subroutine print_iteration(n_solutions, residual_norm, iteration, reduced_dimension)
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
   end subroutine print_iteration
!
!
   end interface
!
   contains
!
   subroutine print_banner_linear_davidson_print_tool()
!!
!!    Print banner
!!    Written by Regina Matveeva, Dec 2021
!!
      use global_out, only: output
!
      implicit none
!
      character(len=100) :: name_
      character(len=500) :: description
!
      name_ = 'Davidson linear equation solver'
!
      description = 'A Davidson solver that solves a linear equation: (A - freq_k) X_k = F. &
                    &This equation is solved in a reduced space. &
                    &A description of the algorithm can be &
                    &found in E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call output%printf('m', ' - ' // trim(name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(name_)) + 6, '-')
!
      call output%printf('n', description, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_linear_davidson_print_tool
!
!
   subroutine print_settings_linear_davidson_print_tool(residual_threshold, &
                                                        max_iterations)
!!
!!    Print settings
!!    Written by Regina Matveeva, Sept 2021
!!
      use global_out, only: output
!
      implicit none
!
      real(dp), intent(in) :: residual_threshold
      integer, intent(in)  :: max_iterations
!
      call output%printf('m', '- Davidson solver settings', fs='(/t3,a)')
!
      call output%printf('m', 'Residual threshold:             (e9.2)', &
                   reals=[residual_threshold], &
                   fs='(/t6,a)')
!
      call output%printf('m', 'Max number of iterations:     (i11)', &
                         ints=[max_iterations], fs='(t6,a)')
!
   end subroutine print_settings_linear_davidson_print_tool
!
!
   subroutine print_summary_linear_davidson_print_tool(n_solutions, &
                                                       converged,   &
                                                       iteration)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad
!!
      use global_out, only: output
!
      implicit none
!
      integer, intent(in) :: n_solutions, iteration
      logical, dimension(n_solutions), intent(in) :: converged
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
   end subroutine print_summary_linear_davidson_print_tool
!
end module linear_davidson_print_tool_class
