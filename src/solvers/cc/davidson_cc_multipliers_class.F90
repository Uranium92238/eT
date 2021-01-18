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
module davidson_cc_multipliers_class
!
!!
!! Davidson coupled cluster multipliers solver class module
!! Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!  
!! Solves the multiplier or left CC ground state equation 
!!
!!    tbar^T A = - eta^T
!!
!! for t-bar. Here, A is the coupled cluster Jacobian 
!!
!!    A_mu,nu = < mu | [H-bar, tau_nu] | HF >,    H-bar = e-T H eT,
!!
!! and 
!!
!!    eta_mu = < HF | [H-bar, tau_nu] | HF >.
!!
!! The multipliers are tbar and give the left CC ground state as 
!!
!!    < Lambda | = < HF | e-T + sum_mu tbar_mu < mu | e-T 
!!
!! The solutions are determined using the Davidson
!! reduced space algorithm, where the eigenvalue problem 
!! is solved in a subspace generated from the residuals* obtained
!! in the preceding iterations. See E. R. Davidson, J. Comput. Phys. 
!! 17, 87 (1975) for more details.
!!
!! * Preconditioned using orbital differences approximation of A 
!!   (See davidson_cc_es solver for more details.)
!!
!
   use kinds
   use ccs_class
   use linear_davidson_tool_class
   use array_utilities, only: get_l2_norm
!
   implicit none
!
   type :: davidson_cc_multipliers
!
      character(len=100) :: name_ = 'Davidson coupled cluster multipliers solver'
!
      character(len=500) :: description = 'A Davidson solver that solves the multiplier equation: t-bar^T A = -η. This linear &
                                            &equation is solved in a reduced space. A description of the algorithm can be &  
                                            &found in E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      integer :: max_iterations, max_dim_red
!
      real(dp) :: residual_threshold
!
      character(len=200) :: storage 
      logical :: restart, records_in_memory
!
      type(timings) :: timer
!
   contains
!
      procedure :: cleanup                         => cleanup_davidson_cc_multipliers
!
      procedure :: print_banner                    => print_banner_davidson_cc_multipliers
      procedure :: print_settings                  => print_settings_davidson_cc_multipliers
!
      procedure :: read_settings                   => read_settings_davidson_cc_multipliers
!
      procedure :: run                             => run_davidson_cc_multipliers
!
      procedure, nopass :: set_precondition_vector => set_precondition_vector_davidson_cc_multipliers     
!
   end type davidson_cc_multipliers
!
!
   interface davidson_cc_multipliers
!
      procedure :: new_davidson_cc_multipliers
!
   end interface davidson_cc_multipliers
!
!
contains
!
!
   function new_davidson_cc_multipliers(wf, restart) result(solver)
!!
!!    New Davidson CC multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(davidson_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
      logical, intent(in) :: restart
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' multipliers')
      call solver%timer%turn_on()
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set default settings, read, & print 
!
      solver%max_iterations      = 100
      solver%residual_threshold  = 1.0d-5
      solver%restart             = restart
      solver%max_dim_red         = 50
      solver%records_in_memory   = .false.
      solver%storage             = 'disk'
!
      call solver%read_settings()
      call solver%print_settings()
!
!     Determine whether to store records in memory or on file
!
      if (trim(solver%storage) == 'memory') then 
!
         solver%records_in_memory = .true.
!
      elseif (trim(solver%storage) == 'disk') then 
!
         solver%records_in_memory = .false.
!
      else 
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // &
                                 trim(solver%storage))
!
      endif 
!
   end function new_davidson_cc_multipliers
!
!
   subroutine print_settings_davidson_cc_multipliers(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(davidson_cc_multipliers) :: solver 
!
      call output%printf('m', '- Davidson CC multipliers solver settings:', fs='(/t3,a)')
!
      call output%printf('m', 'Residual threshold:       (e9.2)', &
                         reals=[solver%residual_threshold], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations: (i9)', &
                         ints=[solver%max_iterations], fs='(t6,a)')
!
   end subroutine print_settings_davidson_cc_multipliers
!
!
   subroutine run_davidson_cc_multipliers(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
      type(linear_davidson_tool), allocatable :: davidson
!
      logical :: converged_residual
!
      real(dp), dimension(:), allocatable :: eta, c, X, multipliers
!
      integer :: iteration 
!
      real(dp) :: residual_norm
!
      call wf%prepare_for_multiplier_equation()
      call wf%initialize_multipliers()
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%construct_eta(eta)
!
      if (get_l2_norm(eta, wf%n_gs_amplitudes) < solver%residual_threshold) then 
!
         call output%printf('m', 'Right hand side is zero to within threshold (e6.1).', &
                                  reals=[solver%residual_threshold], fs='(/t3,a)')
         call output%printf('m', 'Equations solved with multipliers set to zero!')
!
         call mem%alloc(multipliers, wf%n_gs_amplitudes)
         call zero_array(multipliers, wf%n_gs_amplitudes)
!
         call wf%set_multipliers(multipliers)
!
         call wf%print_dominant_x_amplitudes(multipliers, 'tbar')
!
         call mem%dealloc(multipliers, wf%n_gs_amplitudes)
         call mem%dealloc(eta, wf%n_gs_amplitudes)
!
         return
!
      endif
!
!     :: Initialize solver tool and set preconditioner ::
!
      call dscal(wf%n_gs_amplitudes, -one, eta, 1)
!
      davidson = linear_davidson_tool('multipliers',                                  &
                                       wf%n_gs_amplitudes,                            &
                                       min(1.0d-11, tenth*solver%residual_threshold), &
                                       solver%max_dim_red,                            &
                                       eta, 1)
!
      call davidson%initialize(solver%records_in_memory)
!
      call solver%set_precondition_vector(wf, davidson) 
!
!     :: Set start vector / initial guess ::
!
      call wf%set_initial_multipliers_guess(solver%restart)
      call mem%alloc(X, wf%n_gs_amplitudes)
      call wf%get_multipliers(X)
      call davidson%set_trial(X, 1)
!
!     :: Iterative loop ::
!
      call output%printf('n', 'Iteration     Residual norm', fs='(/t3,a)')
      call output%print_separator('n', 27,'-')
!
      iteration = 0
      converged_residual = .false.
!
      call mem%alloc(c, davidson%n_parameters)
!
      do while (.not. converged_residual .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1       
!
!        Reduced space preparations 
!
         if (davidson%red_dim_exceeds_max()) call davidson%set_trials_to_solutions()
!
         call davidson%update_reduced_dim()
!
         call davidson%orthonormalize_trial_vecs() 
!
!        Transform new trial vector and write to file
!
         call davidson%get_trial(c, davidson%dim_red)
!
         call wf%construct_Jacobian_transform('left', c, X, zero)
!
         call davidson%set_transform(X, davidson%dim_red)
!
!        Solve problem in reduced space
!
         call davidson%solve_reduced_problem()
!
!        Check if convergence criterion on residual is satisfied,
!        and, if not, construct a new trial vector using the residual vector
!
         call davidson%construct_residual(X, 1)
         residual_norm = get_l2_norm(X, wf%n_gs_amplitudes)
!
         converged_residual = .true.
         if (residual_norm >= solver%residual_threshold) then 
!
            converged_residual = .false.
            call davidson%add_new_trial(X, 1)
!
         endif 
!
!        Print iteration and residual norm
!
         call output%printf('n', '(i3)            (e11.4)', ints=[iteration], reals=[residual_norm])
!
         call davidson%construct_solution(X, 1)
         call wf%set_multipliers(X)
         call wf%save_multipliers()
!
      enddo ! End of iterative loop 
!
      call mem%dealloc(c, wf%n_gs_amplitudes)
      call mem%dealloc(X, wf%n_gs_amplitudes)
!
      call output%print_separator('n', 27,'-')
!
!     :: Summary of calculation :: 
!
      if (converged_residual) then
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(/t3,a)')
!
         call output%printf('n', '- Davidson CC multipliers solver summary:', fs='(/t3,a)')
!
         call mem%alloc(multipliers, wf%n_gs_amplitudes)
!
         call wf%get_multipliers(multipliers)
         call wf%print_dominant_x_amplitudes(multipliers, 'tbar')
!
         call mem%dealloc(multipliers, wf%n_gs_amplitudes)
!
      else
!
         call output%error_msg('was not able to converge the equations in the given ' //&
                                 'number of iterations in run_davidson_cc_multipliers.')
!
      endif
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
      call davidson%cleanup()
!
   end subroutine run_davidson_cc_multipliers
!
!
   subroutine cleanup_davidson_cc_multipliers(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
      class(davidson_cc_multipliers) :: solver
!
      call wf%save_multipliers()  
!
      call solver%timer%turn_off()
!
      call output%printf('m', '- Finished solving the ' //  &
                         trim(convert_to_uppercase(wf%name_)) // ' multipliers equations', &
                         fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec):  (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6,a)')
!
      call output%printf('m', 'Total cpu time (sec):   (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_davidson_cc_multipliers
!
!
   subroutine print_banner_davidson_cc_multipliers(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(davidson_cc_multipliers) :: solver 
!
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')

!
      call output%printf('n', solver%description, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_davidson_cc_multipliers
!
!
   subroutine set_precondition_vector_davidson_cc_multipliers(wf, davidson)
!!
!!    Set precondition vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets precondition vector to orbital differences 
!!
      implicit none
!
      class(ccs) :: wf
!
      type(linear_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(preconditioner, wf%n_gs_amplitudes)
      call davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_gs_amplitudes)
!
   end subroutine set_precondition_vector_davidson_cc_multipliers
!
!
   subroutine read_settings_davidson_cc_multipliers(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_multipliers) :: solver 
!
      call input%get_keyword_in_section('threshold', 'solver cc multipliers', solver%residual_threshold)
      call input%get_keyword_in_section('max iterations', 'solver cc multipliers', solver%max_iterations)
      call input%get_keyword_in_section('max reduced dimension', 'solver cc multipliers', solver%max_dim_red)
!
      call input%get_keyword_in_section('storage', 'solver cc multipliers', solver%storage)
!
   end subroutine read_settings_davidson_cc_multipliers
!
!
end module davidson_cc_multipliers_class
