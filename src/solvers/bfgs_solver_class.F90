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
module bfgs_solver_class
!
!!
!!    BFGS solver
!!    Written by Eirik F. Kjønstad, 2019 and 2021
!!
!
   use parameters
   use memory_manager_class,              only: mem
   use global_out,                        only: output
   use function_class,                    only: function_
   use redundant_internal_coords_class,   only: redundant_internal_coords
   use bfgs_tool_class,                   only: bfgs_tool
!
   use convergence_tool_class,            only: convergence_tool 
!
   implicit none
!
   type :: bfgs_solver
!
      class(function_), allocatable, private :: energy_function 
!
      class(redundant_internal_coords), pointer, private :: internals 
!
      type(bfgs_tool), private :: bfgs 
!
      class(convergence_tool), allocatable, private :: convergence_checker
!
      integer, private :: n_cartesian, n_internal, iteration, max_iterations
      real(dp), private :: max_step 
!
      real(dp), dimension(:,:), allocatable, private :: geometries 
      real(dp), dimension(:,:), allocatable, private :: gradients 
      real(dp), dimension(:), allocatable, private   :: energies 
! 
      real(dp), dimension(:), allocatable, private :: energy_differences  
      real(dp), dimension(:), allocatable, private :: gradient_maxs  
!
   contains
!
      procedure, public :: run &
                        => run_bfgs_solver
!
      procedure, public :: initialize &
                        => initialize_bfgs_solver
!
      procedure, public :: finalize &
                        => finalize_bfgs_solver
!
      procedure, private :: calculate_gradient_max
      procedure, private :: calculate_energy_difference
      procedure, private :: has_converged
      procedure, private :: update_geometry
      procedure, private :: update_bfgs_hessian_and_calculate_step
      procedure, private :: print_iteration
      procedure, private :: print_summary
      procedure, private :: project_bfgs_hessian
!
   end type bfgs_solver
!
!
   interface bfgs_solver 
!
      procedure :: new_bfgs_solver
!
   end interface bfgs_solver
!
!
contains
!
!
   function new_bfgs_solver(energy_function,    &
                            internals,          &
                            max_iterations,     &
                            max_step,           &
                            energy_threshold,   &
                            gradient_threshold) result(this)
!!
!!    New BFGS solver
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(function_) :: energy_function
      class(redundant_internal_coords), target :: internals 
!
      type(bfgs_solver) :: this
!
      integer, intent(in)  :: max_iterations
      real(dp), intent(in) :: max_step 
      real(dp), intent(in) :: energy_threshold 
      real(dp), intent(in) :: gradient_threshold 
!
      this%energy_function = energy_function 
!
      this%internals => internals
!
      this%n_cartesian = this%energy_function%n_parameters
      this%max_iterations = max_iterations
      this%max_step = max_step
!
      this%convergence_checker = convergence_tool(energy_threshold, &
                                                  gradient_threshold, &
                                                  energy_convergence=.true.)
!
   end function new_bfgs_solver
!
!
   subroutine initialize_bfgs_solver(this)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      real(dp), dimension(:), allocatable :: Hessian_diagonal
!
      call mem%alloc(this%geometries, this%n_cartesian, this%max_iterations)
      call mem%alloc(this%gradients, this%n_cartesian, this%max_iterations)
      call mem%alloc(this%energies, this%max_iterations)
      call mem%alloc(this%energy_differences, this%max_iterations)
      call mem%alloc(this%gradient_maxs, this%max_iterations)
!
      this%n_internal = this%internals%get_n_internal()
!
      this%bfgs = bfgs_tool(this%n_internal, this%max_step)
      call this%bfgs%initialize_arrays()
!
      call mem%alloc(Hessian_diagonal, this%n_internal)
!
      call this%internals%get_approximate_diagonal_Hessian(Hessian_diagonal)
      call this%bfgs%set_initial_Hessian_diagonal(Hessian_diagonal)
!
      call this%project_bfgs_hessian()
!
      call mem%dealloc(Hessian_diagonal, this%n_internal)
!
   end subroutine initialize_bfgs_solver
!
!
   subroutine project_bfgs_hessian(this)
!! 
!!    Project BFGS Hessian
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver), intent(inout) :: this 
!
      real(dp), dimension(:,:), allocatable :: Hessian
!
      call mem%alloc(Hessian, this%n_internal, this%n_internal)
!
      call this%bfgs%get_hessian(Hessian)
!
      call this%internals%Wilson_project_matrix(Hessian)
!
      call this%bfgs%set_hessian(Hessian)
!
      call mem%dealloc(Hessian, this%n_internal, this%n_internal)
!
   end subroutine project_bfgs_hessian
!
!
   subroutine run_bfgs_solver(this)
!!
!!    Run 
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(bfgs_solver) :: this
!
      logical :: converged 
!
      converged = .false.
      this%iteration = 0
!
      call this%energy_function%initialize()
!
      do while (.not. converged .and. this%iteration < this%max_iterations)
!
         this%iteration = this%iteration + 1
!
         this%geometries(:,this%iteration) = this%energy_function%get_parameters()
         this%energies(this%iteration)     = this%energy_function%evaluate()
         this%gradients(:, this%iteration) = this%energy_function%evaluate_gradient() 
!
         this%gradient_maxs(this%iteration)      = this%calculate_gradient_max()
         this%energy_differences(this%iteration) = this%calculate_energy_difference()
!
         call this%print_iteration()
!
         converged = this%has_converged()
!
         if (.not. converged) then 
!
            call this%update_geometry()
!
         else 
!
            call output%printf('m', 'Geometry converged in (i0) iterations!', &
                               ints=[this%iteration], fs='(t6,a)')
!
         end if 
!
      enddo
!
      if (.not. converged) &
         call output%error_msg('Was not able to converge the geometry &
                               &in (i0) iterations!', ints=[this%iteration])
!
      call this%print_summary()
!
   end subroutine run_bfgs_solver
!
!
   function has_converged(this) result(converged)
!!
!!    Has converged?
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      logical :: converged 
!
      converged = this%convergence_checker%has_converged(this%gradient_maxs(this%iteration), &
                                                         this%energy_differences(this%iteration))
!
   end function has_converged 
!
!
   function calculate_gradient_max(this) result(gradient_max)
!!
!!    Calculate gradient_max 
!!    Written by Eirik F. Kjønstad, 2021
!!
      use array_utilities, only: get_abs_max
!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      real(dp) :: gradient_max 
!
      gradient_max = get_abs_max(this%gradients(:, this%iteration), &
                                 this%n_cartesian)
!
   end function calculate_gradient_max
!
!
   function calculate_energy_difference(this) result(energy_difference)
!!
!!    Calculate energy difference
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      real(dp) :: energy_difference
!
      real(dp) :: current_energy, previous_energy 
!
      current_energy = this%energies(this%iteration)
!
      previous_energy = zero
! 
      if (this%iteration .gt. 1) then 
!
         previous_energy = this%energies(this%iteration - 1)
!
      endif       
!
      energy_difference = abs(current_energy - previous_energy)
!
   end function calculate_energy_difference
!
!
   subroutine update_geometry(this)
!!
!!    Update geometry
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      real(dp), dimension(:), allocatable :: s_q, x 
!
      call output%printf('m', 'Geometry not yet converged. Updating geometry via BFGS-RFO step.', fs='(/t6,a)')
!
      call mem%alloc(x, this%n_cartesian)
      call mem%alloc(s_q, this%n_internal)
!
      call this%update_bfgs_hessian_and_calculate_step(s_q)
!
      x = this%geometries(:, this%iteration)
      call this%internals%compute_next_geometry(x, s_q)      
!
      call output%printf('m', 'New and updated geometry identified!', fs='(/t6,a)')
      call this%internals%print_geometry()
!
      call this%energy_function%set_parameters(x)
!
      call mem%dealloc(x, this%n_cartesian)
      call mem%dealloc(s_q, this%n_internal)
!
   end subroutine update_geometry
!
!
   subroutine update_bfgs_hessian_and_calculate_step(this, s_q)
!!
!!    Update BFGS Hessian and calculate step
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      real(dp), dimension(this%n_internal), intent(out) :: s_q 
!
      real(dp), dimension(:), allocatable :: g_x, g_q, dq 
!
      call mem%alloc(g_x, this%n_cartesian)
      call mem%alloc(g_q, this%n_internal)
      call mem%alloc(dq, this%n_internal)
!
      g_x = this%gradients(:, this%iteration)
!
      call this%internals%calculate_internal_gradient(g_x, g_q)
!
      call this%internals%Wilson_project_vector(g_q)
!
      call this%internals%get_internals_difference(dq)
!
      call this%bfgs%update_hessian(dq, g_q)
!
      call this%project_bfgs_hessian()
!
      call this%bfgs%get_step(g_q, s_q)
!
      call this%internals%Wilson_project_vector(s_q)
      call this%internals%Wilson_project_vector(s_q)
!
      call mem%dealloc(g_x, this%n_cartesian)
      call mem%dealloc(g_q, this%n_internal)
      call mem%dealloc(dq, this%n_internal)
!
   end subroutine update_bfgs_hessian_and_calculate_step
!
!
   subroutine print_iteration(this)
!!
!!    Print iteration 
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      call output%printf('m', '- Geometry optimization iteration (i0):', ints=[this%iteration], fs='(/t3,a)')
!
      call output%printf('m', '                             max(gradient)     dE', fs='(/t6,a)')
      call output%print_separator('m', 60,'-', fs='(t6,a)')
!
      call output%printf('m', 'Convergence criterion:      (e11.4)       (e11.4)', &
            reals=[this%convergence_checker%residual_threshold, this%convergence_checker%energy_threshold], fs='(t6,a)')
      call output%printf('m', 'Current value:              (e11.4)       (e11.4)', &
            reals=[this%gradient_maxs(this%iteration), this%energy_differences(this%iteration)], fs='(t6,a)')
!
      call output%print_separator('m', 60,'-', fs='(t6,a)')
!
   end subroutine print_iteration
!
!
   subroutine print_summary(this)
!!
!!    Print summary
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      integer :: i
!
      call output%printf('m', '- Summary of geometry optimization:', fs='(/t3,a/)')
!
      call output%printf('m', 'Iteration     Energy               max(gradient)     dE', fs='(t6,a)')
      call output%print_separator('m', 65,'-', fs='(t6,a)')
!
      do i = 1, this%iteration
!
         call output%printf('m', '(i2)         (f19.12)    (e11.4)       (e11.4)', &
                        ints=[i], reals=[this%energies(i), this%gradient_maxs(i), this%energy_differences(i)], &
                        fs='(t6,a)')
!
      end do 
!
      call output%print_separator('m', 65,'-', fs='(t6,a)')
!
      call this%internals%print_geometry()
!
   end subroutine print_summary
!
!
   subroutine finalize_bfgs_solver(this)
!!
!!    Finalize
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(bfgs_solver) :: this 
!
      call mem%dealloc(this%geometries, this%n_cartesian, this%max_iterations)
      call mem%dealloc(this%gradients, this%n_cartesian, this%max_iterations)
      call mem%dealloc(this%energies, this%max_iterations)
      call mem%dealloc(this%energy_differences, this%max_iterations)
      call mem%dealloc(this%gradient_maxs, this%max_iterations)
      call this%bfgs%cleanup()
!
   end subroutine finalize_bfgs_solver
!
!
end module bfgs_solver_class
