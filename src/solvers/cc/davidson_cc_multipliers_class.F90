!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
!!    Davidson coupled cluster multipliers solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!  
!!    Solves the multiplier or left CC ground state equation 
!!
!!       t-bar^T A = - eta^T
!!
!!    for t-bar. The solutions are determined using the 
!!    Davidson reduced space algorithm, where the eigenvalue problem 
!!    is solved in a subspace generated from the residuals obtained
!!    in the preceding iterations. See E. R. Davidson, J. Comput. Phys. 
!!    17, 87 (1975) for more details.
!!
!
   use kinds
   use ccs_class
   use linear_davidson_tool_class
!
   implicit none
!
   type :: davidson_cc_multipliers
!
      character(len=100) :: name_ = 'Davidson coupled cluster multipliers solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
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
   function new_davidson_cc_multipliers(wf) result(solver)
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
      solver%residual_threshold  = 1.0d-6
      solver%restart             = .false.
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
      call output%printf('- Davidson CC multipliers solver settings:', pl='m', fs='(/t3,a)')
!
      call output%printf('Residual threshold:       (e9.2)', pl='m', fs='(/t6,a)', reals=[solver%residual_threshold])
      call output%printf('Max number of iterations: (i9)', pl='m', fs='(t6,a)', ints=[solver%max_iterations])
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
      type(linear_davidson_tool) :: davidson
!
      logical :: converged_residual
!
      real(dp), dimension(:), allocatable :: eta, c, multipliers, residual 
!
      integer :: iteration
!
      real(dp) :: residual_norm, ddot, norm_trial
!
!     :: Initialize solver tool and set preconditioner ::
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%construct_eta(eta)
!
      davidson = linear_davidson_tool('multipliers', wf%n_gs_amplitudes, &
                  solver%residual_threshold, solver%max_dim_red, -eta, solver%records_in_memory)
!
      call solver%set_precondition_vector(wf, davidson)
!
!     :: Set start vector / initial guess ::
!
      if (solver%restart) then ! Read multiplier vector from file
!
         call output%printf('Requested restart. Reading multipliers from file.', pl='m', fs='(/t3,a)')
!
         call wf%is_restart_safe('ground state')
         call wf%initialize_multipliers()
         call wf%read_multipliers()
!
         call mem%alloc(multipliers, wf%n_gs_amplitudes)
         call wf%get_multipliers(multipliers)
         call wf%destruct_multipliers()
!
         norm_trial = sqrt(ddot(wf%n_gs_amplitudes, multipliers, 1, multipliers, 1))
         call dscal(wf%n_gs_amplitudes, one/norm_trial, multipliers, 1)
!
         call davidson%set_trial(multipliers, 1)
         call mem%dealloc(multipliers, wf%n_gs_amplitudes)
!
      else ! Use - eta_mu / eps_mu as first trial 
!
         call dscal(wf%n_gs_amplitudes, -one, eta, 1)
         call davidson%precondition(eta)
!
         norm_trial = sqrt(ddot(wf%n_gs_amplitudes, eta, 1, eta, 1))
         call dscal(wf%n_gs_amplitudes, one/norm_trial, eta, 1)
!
         call davidson%set_trial(eta, 1)
!
      endif 
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
!     :: Iterative loop ::
!
      call output%printf('Iteration     Residual norm', pl='n', fs='(/t3,a)')
      call output%printf('---------------------------', pl='n')
!
      iteration = 0
      converged_residual = .false.
!
      do while (.not. converged_residual .and. (iteration .le. solver%max_iterations))
!
         iteration = iteration + 1       
         call davidson%iterate()
!
!        Transform new trial vector and write to file
!
         call mem%alloc(c, davidson%n_parameters)
!
         call davidson%get_trial(c, davidson%dim_red)
!
         call wf%construct_Jacobian_transform('left', c, zero)
!
         call davidson%set_transform(c, davidson%dim_red)
!
         call mem%dealloc(c, davidson%n_parameters)
!
!        Solve problem in reduced space
!
         call davidson%solve_reduced_problem()
!
!        Check if convergence criterion on residual is satisfied,
!        and, if not, construct a new trial vector using the residual vector
!
         call mem%alloc(residual, wf%n_gs_amplitudes)
!
         call davidson%construct_residual(residual, 1)
         residual_norm = get_l2_norm(residual, wf%n_gs_amplitudes)
!
         converged_residual = .true.
         if (residual_norm >= solver%residual_threshold) then 
!
            converged_residual = .false.
            call davidson%construct_next_trial(residual)
!
         endif 
!
         call mem%dealloc(residual, wf%n_gs_amplitudes)
!
!        Print iteration and residual norm
!
         call output%printf('(i3)            (e11.4)', pl='n', &
                              ints=[iteration], reals=[residual_norm])
!
      enddo ! End of iterative loop 
!
      call output%printf('---------------------------', pl='n')
!
!     :: Summary of calculation :: 
!
      if (converged_residual) then
!
         call output%printf('Convergence criterion met in (i0) iterations!', pl='m', &
                              ints=[iteration], fs='(/t3,a)')
!
         call output%printf('- Davidson CC multipliers solver summary:', pl='n', fs='(/t3,a)')
!
         call mem%alloc(multipliers, wf%n_gs_amplitudes)
         call davidson%construct_solution(multipliers, 1)
!
         call wf%print_dominant_x_amplitudes(multipliers, 'tbar')
!
         call wf%initialize_multipliers()
         call wf%set_multipliers(multipliers)
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
      call output%printf('- Finished solving the ' // trim(convert_to_uppercase(wf%name_)) // &
                           ' multipliers equations', pl='m', fs='(/t3,a)')
!
      call output%printf('Total wall time (sec):  (f20.5)', pl='m', &
                           reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6,a)')
!
      call output%printf('Total cpu time (sec):   (f20.5)', pl='m', &
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
      call output%printf(':: ' // solver%name_, pl='m', fs='(//t3,a)')
      call output%printf(':: ' // solver%author, pl='m', fs='(t3,a)')
      call output%printf(solver%description, pl='n', ffs='(/t3,a)', fs='(t3,a)', lfs='(t3,a)')
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
!
      if (input%requested_keyword_in_section('restart', 'solver cc multipliers')) solver%restart = .true.    
!
      call input%get_keyword_in_section('storage', 'solver cc multipliers', solver%storage)
!
   end subroutine read_settings_davidson_cc_multipliers
!
!
end module davidson_cc_multipliers_class
