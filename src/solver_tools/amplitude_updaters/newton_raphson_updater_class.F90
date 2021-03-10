!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
module newton_raphson_updater_class
!!
!!    Newton-Raphson updater class module
!!    Written by Eirik F. Kjønstad, 2020
!! 
!!    Updates the amplitudes according to a Newton-Raphson estimate,
!!
!!       amplitudes_mu = amplitudes_mu - (A^-1)_mu,nu residual_nu,
!!
!!    where A is the Jacobian matrix.
!!
   use parameters
   use amplitude_updater_class,    only: amplitude_updater
   use ccs_class,                  only: ccs
   use memory_manager_class,       only: mem
   use global_out,                 only: output
   use linear_davidson_tool_class, only: linear_davidson_tool
!
   implicit none
!
   type, extends(amplitude_updater) :: newton_raphson_updater
!
      real(dp), private :: relative_threshold 
      logical, private  :: records_in_memory
      integer, private  :: max_dim_red
      integer, private  :: max_iterations
!
      type(linear_davidson_tool), allocatable :: davidson 
!
   contains
!
      procedure :: calculate_update => calculate_update_newton_raphson_updater
!
   end type newton_raphson_updater
!
!
   interface newton_raphson_updater
!
      procedure :: new_newton_raphson_updater
!
   end interface newton_raphson_updater
!
!
contains 
!
!
   function new_newton_raphson_updater(n_amplitudes,             &
                                            scale_amplitudes,    &
                                            scale_residual,      &
                                            relative_threshold,  &
                                            records_in_memory,   &
                                            max_iterations) result(this)
!!
!!    New Newton-Raphson updator
!!    Written by Eirik F. Kjønstad, 2020
!!

!!    scale_amplitudes / scale_residual: if true, scales the diagonal of the double
!!                                       amplitudes by a factor of two

!!    relative_threshold:                convergence threshold relative to residual norm
!!    records_in_memory:                 in Davidson, place trials/transforms in memory; 
!!                                       otherwise, they are on disk
!!    max_iterations:                    error if reached without convergence
!!
      implicit none 
!
      integer, intent(in) :: n_amplitudes, max_iterations
!
      logical, intent(in) :: records_in_memory, scale_amplitudes, scale_residual
!
      real(dp), intent(in) :: relative_threshold
!
      type(newton_raphson_updater) :: this 
!
      this%n_amplitudes       = n_amplitudes
      this%scale_amplitudes   = scale_amplitudes
      this%scale_residual     = scale_residual
      this%relative_threshold = relative_threshold
      this%records_in_memory  = records_in_memory
      this%max_dim_red        = 50
      this%max_iterations     = max_iterations
!
      this%davidson = linear_davidson_tool(name_        = 'newton_raphson_amplitude_updator', &
                                      n_parameters      = this%n_amplitudes,                  &
                                      lindep_threshold  = 1.0d-11,                            &
                                      max_dim_red       = this%max_dim_red,                   &
                                      n_equations       = 1,                                  &
                                      records_in_memory = this%records_in_memory)
!
   end function new_newton_raphson_updater
!
!
   subroutine calculate_update_newton_raphson_updater(this, wf, residual, update)
!!
!!    Calculate update
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018-2020
!!
!!    Solves the equation A (update) = - (residual) for the update vector.
!!
!
      use timings_class, only: timings
      use array_utilities, only: get_l2_norm
!
      implicit none 
!
      class(newton_raphson_updater), intent(inout) :: this  
!
      class(ccs), intent(inout) :: wf 
!
      real(dp), dimension(this%n_amplitudes), intent(inout) :: residual
      real(dp), dimension(this%n_amplitudes), intent(inout) :: update
!
      real(dp), dimension(:), allocatable :: preconditioner, c, x
!
      real(dp) :: rho_norm, rho_threshold
!
      integer :: iteration, trial
      logical :: converged
!
      type(timings), allocatable :: iteration_timer, timer 
!
      timer = timings('Newton-Raphson total update time', pl='n')
      call timer%turn_on()
!
!     Preparations for iterative loop
!
      rho_threshold = get_l2_norm(residual, this%n_amplitudes) * this%relative_threshold
!
      call this%davidson%set_lindep_threshold(min(1.0d-11, rho_threshold))
!
      call this%davidson%initialize()
      call this%davidson%set_rhs(residual)
!
      call mem%alloc(preconditioner, this%n_amplitudes)
      call wf%get_orbital_differences(preconditioner, this%n_amplitudes)
!
      call this%davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, this%n_amplitudes)
!
      call this%davidson%set_trials_to_preconditioner_guess()
!
      call wf%prepare_for_jacobian() 
!
      call output%printf('n', 'Iteration     Residual norm', fs='(/t6,a)')
      call output%print_separator('n', 29, '-', fs='(t6,a)')
!
      iteration = 0
      converged = .false.
!
      call mem%alloc(c, this%davidson%n_parameters)
      call mem%alloc(x, this%davidson%n_parameters)
!
      iteration_timer = timings('Newton-Raphson iteration time', pl='v')
!
!     Iterative loop
!
      do while (.not. converged .and. (iteration .le. this%max_iterations))
!
         iteration = iteration + 1
         call iteration_timer%turn_on()
!
!        Reduced space preparations 
!
         if (this%davidson%red_dim_exceeds_max()) call this%davidson%set_trials_to_solutions()
!
         call this%davidson%update_reduced_dim()
!
         call this%davidson%orthonormalize_trial_vecs() 
!
!        Transform new trial vector and set transform
!
         trial = this%davidson%first_new_trial() ! only one new trial, so first = last
!
         call this%davidson%get_trial(c, trial)
!
         call wf%construct_Jacobian_transform('right', c, x, w = zero)
!
         call this%davidson%set_transform(x, trial)
!
!        Solve equation in the trial space
!
         call this%davidson%solve_reduced_problem()
!
!        Test for convergence - and add residual to trial space if not converged
!
         call this%davidson%construct_residual(x, 1)
         rho_norm = get_l2_norm(x, this%n_amplitudes)
!
         converged = .true.
         if (rho_norm >= rho_threshold) then 
!
            converged = .false.
            call this%davidson%add_new_trial(x, 1)
!
         endif 
!
         call output%printf('n', '(i3)          (e11.4)', &
                            ints=[iteration], reals=[rho_norm], fs='(t6,a)')
!
         call iteration_timer%turn_off()
         call iteration_timer%reset()
!
      enddo
!
      call output%print_separator('n', 29, '-', fs='(t6,a/)')
!
      if (.not. converged) then
!
         call output%error_msg('was not able to converge the equations in the ' // &
                               'specified max number of iterations!')
!
      endif
!
      call mem%dealloc(x, this%davidson%n_parameters)
      call mem%dealloc(c, this%davidson%n_parameters)
!
      call this%davidson%construct_solution(update, 1)
      call dscal(this%n_amplitudes, -one, update, 1)
!
      call this%davidson%cleanup()
!
      call timer%turn_off()
!
   end subroutine calculate_update_newton_raphson_updater
!
!
end module newton_raphson_updater_class
