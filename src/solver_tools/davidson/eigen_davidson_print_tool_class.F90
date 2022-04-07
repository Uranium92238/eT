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
module eigen_davidson_print_tool_class
!
!!
!!    Eigen Davidson print tool class module
!!    Written by Regina Matveeva, Dec 2021
!!
!
   use eigen_davidson_tool_class, only: eigen_davidson_tool
!
   implicit none
!
   type :: eigen_davidson_print_tool
!
   contains
!
      procedure, public, nopass :: print_banner &
                                => print_banner_eigen_davidson_print_tool
!
      procedure, public, nopass :: print_iteration &
                                => print_iteration_eigen_davidson_print_tool
!
      procedure, public, nopass :: print_settings &
                                => print_settings_eigen_davidson_print_tool
!
      procedure, public :: print_summary &
                           => print_summary_eigen_davidson_print_tool
!
   end type eigen_davidson_print_tool
!
   contains
!
!
   subroutine print_banner_eigen_davidson_print_tool()
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
      name_ = 'Davidson eigenvalue equation solver'
!
      description = 'A Davidson solver that solves an eigenvalue equation: M x = omega x. &
                    &This equation is solved in a reduced space. &
                    &A description of the algorithm can be &
                    &found in E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call output%printf('m', trim(name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(name_)), '-')
!
      call output%printf('n', description, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_eigen_davidson_print_tool
!
!
   subroutine print_settings_eigen_davidson_print_tool(n_solutions, max_iterations)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad, May 2021
!!
      use global_out, only: output
!
      implicit none
!
      integer, intent(in) :: n_solutions, max_iterations
!
      call output%printf('m', '- Davidson solver settings', fs='(/t3,a)')
!
      call output%printf('m', 'Number of singlet states:     (i11)', &
                         ints=[n_solutions], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations:     (i11)', &
                         ints=[max_iterations], fs='(t6,a)')
!
   end subroutine print_settings_eigen_davidson_print_tool
!
   subroutine print_iteration_eigen_davidson_print_tool(n_solutions,   &
                                                        omega,         &
                                                        new_omega,     &
                                                        residual_norm, &
                                                        iteration,     &
                                                        reduced_dimension)
!!
!!    Print iteration
!!    Written by Sarai D. Folkestad, May 2021
!!
      use parameters
      use global_out, only: output
!
      implicit none
!
      integer,  intent(in) :: n_solutions, iteration, reduced_dimension
      real(dp), dimension(n_solutions), intent(in) :: residual_norm, omega
      complex(dp), dimension(n_solutions), intent(in) :: new_omega
!
      integer :: n
!
      call output%printf('n', 'Iteration:               (i4)', &
                         ints=[iteration], fs='(/t3,a)')
      call output%printf('n', 'Reduced space dimension: (i4)', ints=[reduced_dimension])
!
      call output%printf('n', &
                         'Root  Eigenvalue (Re)   Eigenvalue (Im)    |residual|   Delta E (Re)', &
                         fs='(/t3,a)', ll=120)
!
      call output%print_separator('n', 73,'-')
!
      do n = 1, n_solutions
!
         call output%printf('n', '(i4) (f16.12)  (f16.12)    (e11.4)  (e11.4) ', &
                            ints=[n],                            &
                            reals=[real(new_omega(n)),           &
                                   aimag(new_omega(n)),          &
                                   residual_norm(n),             &
                                   abs(real(new_omega(n)) - omega(n))], &
                                   ll=120)
!
      enddo
!
      call output%print_separator('n', 73,'-')
!
   end subroutine print_iteration_eigen_davidson_print_tool
!
!
   subroutine print_summary_eigen_davidson_print_tool(this, n_solutions, converged, iteration)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad, May 2021
!!
      use global_out, only: output
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(eigen_davidson_print_tool), intent(in) :: this
      integer, intent(in) :: n_solutions, iteration
      logical, dimension(n_solutions), intent(in) :: converged
!
      call do_nothing(this)
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
   end subroutine print_summary_eigen_davidson_print_tool
!
!
end module eigen_davidson_print_tool_class
