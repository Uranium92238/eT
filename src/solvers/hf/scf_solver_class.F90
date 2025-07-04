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
module scf_solver_class
!
!!
!!		Self-consistent field solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad  2018-2020
!!
!!    Solves the self-consistent field equations
!!
!!       F(C) C = C e  (*)
!!
!!    The equations are solved by
!!
!!    1. Setting an initial guess for C
!!
!!    --------- ITERATIVE LOOP ----------
!!
!!    2. Constructing F and gradient
!!
!!    3. Convergence test on gradient
!!
!!    4. Extrapolation (convergence acceleration)
!!
!!    5. Solving (*)
!!
!!    -----------------------------------
!!
!
   use parameters
!
   use hf_class,                  only: hf
   use memory_manager_class,      only: mem
   use timings_class,             only: timings
   use global_out,                only: output
!
   use accelerator_factory_class, only: accelerator_factory
   use accelerator_tool_class,    only: accelerator_tool
   use convergence_tool_class,    only: convergence_tool
!
   implicit none
!
   type :: scf_solver
!
      logical, private :: skip
!
!     Iterative loop
      integer,  private :: max_iterations
!
!     Equation specifications
      integer, private :: dim_
      integer, private :: n_equations
      integer, private :: gradient_dimension
      integer, private :: packed_F_dimension
!
!     Convergence acceleration (none, DIIS, or CROP)
      character(len=200),                     private :: acceleration_type
      class(accelerator_tool), allocatable,   private :: accelerator
      type(accelerator_factory), allocatable, private :: accelerator_creator
!
      class(convergence_tool), allocatable :: convergence_checker
!
   contains
!
      procedure :: run &
                => run_scf_solver
!
!
      procedure, private :: do_scf_step &
                         => do_scf_step_scf_solver
!
      procedure, private :: prepare_accelerator &
                         => prepare_accelerator_scf_solver
!
      procedure, private :: print_settings &
                         => print_settings_scf_solver
!
      procedure, nopass, private :: print_iteration &
                                 => print_iteration_scf_solver
!
      procedure, nopass, private :: print_summary &
                                 => print_summary_scf_solver
!
      procedure, nopass, private :: print_iteration_banner &
                                 => print_iteration_banner_scf_solver
!
      procedure, private :: skip_scf &
                         => skip_scf_scf_solver
!
   end type scf_solver
!
   interface scf_solver
!
      procedure :: new_scf_solver
      procedure :: new_scf_solver_from_parameters
!
   end interface scf_solver
!
!
contains
!
   function new_scf_solver(acceleration_type,   &
                           max_iterations,      &
                           skip,                &
                           dim_,                &
                           n_equations,         &
                           gradient_dimension,  &
                           convergence_checker) result(solver)
!!
!!    New SCF solver
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(scf_solver)             :: solver
      character(len=*),          intent(in) :: acceleration_type
      logical,                   intent(in) :: skip
      type(convergence_tool),    intent(in) :: convergence_checker
      integer, intent(in) :: dim_
      integer, intent(in) :: n_equations
      integer, intent(in) :: gradient_dimension, max_iterations
!
!     Set defaults
!
      solver%max_iterations      = max_iterations
      solver%acceleration_type   = acceleration_type
      solver%skip                = skip
      solver%dim_                = dim_
      solver%n_equations         = n_equations
      solver%gradient_dimension  = gradient_dimension
!
!     Initialize convergence checker with default threshols
!
      solver%convergence_checker = convergence_checker
!
      solver%accelerator_creator = accelerator_factory('solver scf')
!
      call solver%print_settings()
!
   end function new_scf_solver
!
!
   function new_scf_solver_from_parameters(max_iterations,        &
                                           gradient_threshold,    &
                                           acceleration_type,     &
                                           skip,                  &
                                           energy_threshold, gradient_dimension) result(solver)
!!
!!    New SCF solver from parameters
!!    Written by Tor S. Haugland and Sarai D. Folkestad, 2019-2020
!!
      implicit none
!
      type(scf_solver) :: solver
!
      integer,            intent(in)           :: max_iterations
      real(dp),           intent(in)           :: gradient_threshold
      character(len=*),   intent(in)           :: acceleration_type
      logical,            intent(in)           :: skip
      real(dp),           intent(in), optional :: energy_threshold
      integer, intent(in) :: gradient_dimension
!
      solver%max_iterations      = max_iterations
      solver%acceleration_type   = acceleration_type
      solver%skip                = skip
      solver%gradient_dimension  = gradient_dimension
!
      solver%accelerator_creator = accelerator_factory('solver scf')
!
      if (present(energy_threshold)) then
!
         solver%convergence_checker = convergence_tool(energy_threshold  = energy_threshold,   &
                                                       residual_threshold = gradient_threshold, &
                                                       energy_convergence = .true.)
!
      else
!
         solver%convergence_checker = convergence_tool(energy_threshold  = gradient_threshold, &
                                                       residual_threshold = gradient_threshold, &
                                                       energy_convergence = .false.)
      endif
!
      call solver%print_settings()
!
   end function new_scf_solver_from_parameters
!
!
   subroutine run_scf_solver(solver, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018-2020
!!
      use array_utilities, only: get_abs_max
!
      implicit none
!
      class(scf_solver) :: solver
      class(hf)         :: wf
!
      logical  :: converged
      real(dp) :: energy, previous_energy, max_gradient
      integer  :: iteration
!
      real(dp), dimension(:), allocatable       :: F, gradient
      real(dp), dimension(:,:), allocatable     :: e
      real(dp), dimension(:,:,:), allocatable   :: C
!
      type(timings) :: timer, iteration_timer
!
      timer = timings('SCF solver', pl='normal')
      call timer%turn_on()
!
      iteration_timer = timings('SCF iteration', pl='normal')
!
!     If there is only one AO, equations are solved
!
      if (solver%dim_ == 1) return
!
      if (solver%skip) then
!
         call solver%skip_scf(wf)
         call timer%turn_off()
         return
!
      endif
!
      solver%packed_F_dimension = ((solver%dim_)*(solver%dim_ + 1)/2)*solver%n_equations
!
      call solver%prepare_accelerator()
!
      converged = .false.
      iteration = 0
      energy    = zero
!
      call mem%alloc(F, solver%packed_F_dimension)
      call mem%alloc(C, solver%dim_, solver%dim_, solver%n_equations)
      call mem%alloc(gradient, solver%gradient_dimension)
      call mem%alloc(e, solver%dim_, solver%n_equations)
!
      call solver%print_iteration_banner()
!
      do while ((.not. converged) .and. (iteration .lt. solver%max_iterations))
!
         call iteration_timer%turn_on()
!
         iteration = iteration + 1
!
         call wf%get_F(F)
         call wf%get_gradient(gradient)
!
         previous_energy = energy
         energy          = wf%get_energy()
         max_gradient    = get_abs_max(gradient, solver%gradient_dimension)
!
         converged =  solver%convergence_checker%has_converged(max_gradient,              &
                                                               energy - previous_energy,  &
                                                               iteration)

         call solver%print_iteration(iteration, energy, previous_energy, max_gradient)
!
         if (.not. converged) then
!
            call solver%accelerator%do_(F, gradient)
            call solver%do_scf_step(F, C, e)
            call wf%set_C_and_e(C, e)
!
         endif
!
         call iteration_timer%turn_off()
         call iteration_timer%reset()
!
      enddo
!
      call solver%print_summary(iteration, converged)
!
      call mem%dealloc(F, solver%packed_F_dimension)
      call mem%dealloc(C, solver%dim_, solver%dim_, solver%n_equations)
      call mem%dealloc(gradient, solver%gradient_dimension)
      call mem%dealloc(e, solver%dim_, solver%n_equations)
!
      call solver%accelerator%finalize()
!
      call timer%turn_off()
!
   end subroutine run_scf_solver
!
!
   subroutine do_scf_step_scf_solver(solver, F_packed, C, e)
!!
!!    Do SCF step
!!    Written by Sarai D. Folkestad
!!
!!    Solves the equation
!!
!!       F(C_i) C_i+1 = e C_i+1
!
      use reordering, only: squareup
      use array_utilities, only: diagonalize_symmetric
!
      implicit none
!
      class(scf_solver),                                                 intent(in)  :: solver
      real(dp), dimension(solver%packed_F_dimension),                    intent(in)  :: F_packed
      real(dp), dimension(solver%dim_, solver%dim_, solver%n_equations), intent(out) :: C
      real(dp), dimension(solver%dim_, solver%n_equations),              intent(out) :: e
!
      integer                             :: i, offset
!
      type(timings) :: timer
!
      timer = timings('SCF diagonalization', pl='normal')
      call timer%turn_on()
!
      do i = 1, solver%n_equations
!
         offset = ((solver%dim_)*(solver%dim_ + 1)/2)*(i-1)
!
         call squareup(F_packed(offset + 1 :), C(:,:,i), solver%dim_)
         call diagonalize_symmetric(C(:,:,i), solver%dim_, e(:,i))
!
      enddo
!
      call timer%turn_off()
!
   end subroutine do_scf_step_scf_solver
!
!
   subroutine prepare_accelerator_scf_solver(solver)
!!
!!    Prepare accelerator
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(scf_solver), intent(inout) :: solver
!
      solver%accelerator = solver%accelerator_creator%create(solver%acceleration_type,   &
                                                             solver%packed_F_dimension,  &
                                                             solver%gradient_dimension)
      call solver%accelerator%initialize()
!
   end subroutine prepare_accelerator_scf_solver
!
!
   subroutine print_iteration_scf_solver(iteration, energy, previous_energy, max_gradient)
!!
!!    Print iteration
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      integer,  intent(in) :: iteration
      real(dp), intent(in) :: energy
      real(dp), intent(in) :: previous_energy
      real(dp), intent(in) :: max_gradient
!
      call output%printf('n', '(i4)  (f25.12)    (e11.4)    (e11.4)', &
                           ints=[iteration], reals=[energy, max_gradient, &
                           abs(energy - previous_energy)], fs='(t3,a)')
!
   end subroutine print_iteration_scf_solver
!
!
   subroutine print_summary_scf_solver(iteration, converged)
!!
!!    Print summary
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      integer, intent(in) :: iteration
      logical, intent(in) :: converged
!
      call output%print_separator('n', 63,'-')
!
      if (converged) then
!
         call output%printf('n', 'Convergence criterion met in (i0) iterations!', ints=[iteration])
!
      else
!
         call output%error_msg('Was not able to converge the equations in &
                                 &the given number of maximum iterations.')
!
      endif
!
   end subroutine print_summary_scf_solver
!
!
   subroutine print_iteration_banner_scf_solver()
!!
!!    Print iteration banner
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      call output%printf('n', 'Iteration       Energy (a.u.)      Max(grad.)    Delta E (a.u.)', &
                        fs='(/t3,a)')
      call output%print_separator('n', 63,'-')
!
   end subroutine print_iteration_banner_scf_solver
!
!
   subroutine print_settings_scf_solver(solver)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(scf_solver), intent(in) :: solver
!
      call output%printf('m', '- SCF solver settings:',fs='(/t3,a)')
!
      call output%printf('m', 'Maximum iterations:           (i11)', &
                         ints=[solver%max_iterations], fs='(/t6,a)')
!
      call output%printf('m', 'Acceleration type:            (a11)', &
                         chars=[trim(solver%acceleration_type)], fs='(t6,a)')
!
      call solver%convergence_checker%print_settings()
!
   end subroutine print_settings_scf_solver
!
!
   subroutine skip_scf_scf_solver(solver, wf)
!!
!!    Skip SCF
!!    Written by Sarai D. Folkestad.
!!
!!    Prints the energy and maximal gradient element when SCF is skipped.
!!
      use array_utilities, only: get_abs_max
!
      implicit none
!
      class(scf_solver), intent(in) :: solver
      class(hf), intent(in)         :: wf
!
      real(dp)                            :: energy, max_gradient
      real(dp), dimension(:), allocatable :: gradient
!
      call output%warning_msg('skipping SCF solver!')
!
      call mem%alloc(gradient, solver%gradient_dimension)
      call wf%get_gradient(gradient)
!
      energy          = wf%get_energy()
      max_gradient    = get_abs_max(gradient, solver%gradient_dimension)
!
      call mem%dealloc(gradient, solver%gradient_dimension)

      call solver%print_iteration_banner()
      call solver%print_iteration(iteration = 1,            &
                                  energy = energy,          &
                                  previous_energy = zero,   &
                                  max_gradient = max_gradient)
      call solver%print_summary(iteration = 1, converged = .true.)
!
   end subroutine skip_scf_scf_solver
!
!
end module scf_solver_class
