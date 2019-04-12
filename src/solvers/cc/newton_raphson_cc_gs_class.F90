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
module newton_raphson_cc_gs_class
!
!!
!!		Newton-Raphson coupled cluster ground state solver class module
!!		Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!  
!
   use kinds
   use file_class
   use ccs_class
   use diis_tool_class
!
   implicit none
!
   type :: newton_raphson_cc_gs
!
      character(len=100) :: tag = 'Newton-Raphson coupled cluster ground state solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2019'
!
      character(len=500) :: description1 = 'A Newton-Raphson CC ground state equations solver. Solves the &
                                             &ground state equation using updates Δt based on the Newton equation, &
                                             &A Δt = -Ω, where the A is the Jacobian and Ω the omega vector.'
!
      integer :: max_iterations
      integer :: diis_dimension
!
      real(dp) :: omega_threshold
      real(dp) :: energy_threshold
!
      integer  :: max_micro_iterations
      real(dp) :: micro_residual_threshold
      real(dp) :: relative_micro_residual_threshold
!
      logical :: restart
!
   contains
!     
      procedure, nopass :: cleanup          => cleanup_newton_raphson_cc_gs
      procedure :: prepare                  => prepare_newton_raphson_cc_gs
      procedure :: run                      => run_newton_raphson_cc_gs
!
      procedure :: do_micro_iterations      => do_micro_iterations_newton_raphson_cc_gs
!
      procedure :: read_settings            => read_settings_newton_raphson_cc_gs
      procedure :: print_banner             => print_banner_newton_raphson_cc_gs
      procedure :: print_settings           => print_settings_newton_raphson_cc_gs
      procedure, nopass :: print_summary            => print_summary_newton_raphson_cc_gs
!
   end type newton_raphson_cc_gs
!
!
contains
!
!
   subroutine prepare_newton_raphson_cc_gs(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(newton_raphson_cc_gs) :: solver
!
      class(ccs) :: wf
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%max_iterations            = 100
      solver%diis_dimension            = 100
      solver%max_micro_iterations      = 100
      solver%relative_micro_residual_threshold = 1.0d-2
      solver%energy_threshold          = 1.0d-6
      solver%omega_threshold           = 1.0d-6
      solver%restart                   = .false.
!
!     Read & print settings (thresholds, etc.)
!
      call solver%read_settings()
      call solver%print_settings()
!
!     Set the amplitudes to the initial guess or read if restart
!
      call wf%initialize_amplitudes()
!
!     Prepare restart information file 
!
      if (solver%restart) then
!
         call wf%read_amplitudes()
         call wf%integrals%write_t1_cholesky(wf%t1) 
! 
      else
!
         call wf%integrals%write_t1_cholesky(wf%t1) 
         call wf%set_initial_amplitudes_guess()
!
      endif
!
   end subroutine prepare_newton_raphson_cc_gs
!
!
   subroutine print_settings_newton_raphson_cc_gs(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(newton_raphson_cc_gs) :: solver 
!
      write(output%unit, '(/t3,a)')      '- DIIS accelerated Newton-Raphson CC ground state solver settings:'
!
      write(output%unit, '(/t6,a32,e9.2)') 'Omega threshold:                ', solver%omega_threshold
      write(output%unit, '(t6,a32,e9.2)')  'Energy threshold:               ', solver%energy_threshold
      write(output%unit, '(t6,a32,e9.2)')  'Relative micro threshold:       ', solver%relative_micro_residual_threshold

      write(output%unit, '(/t6,a32,i9)')   'DIIS dimension:                 ', solver%diis_dimension
      write(output%unit, '(t6,a32,i9)')    'Max number of iterations:       ', solver%max_iterations
      write(output%unit, '(t6,a32,i9)')    'Max number of micro-iterations: ', solver%max_micro_iterations
!
      flush(output%unit)
!
   end subroutine print_settings_newton_raphson_cc_gs
!
!
   subroutine run_newton_raphson_cc_gs(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(newton_raphson_cc_gs) :: solver
!
      class(ccs) :: wf
!
      type(diis_tool) :: diis_manager
!
      real(dp), dimension(:), allocatable :: omega, dt, t  
!
      integer :: iteration
!
      logical :: converged, converged_omega, converged_energy
!
      real(dp) :: energy, prev_energy, omega_norm
!
      call diis_manager%init('cc_gs_diis', wf%n_gs_amplitudes, wf%n_gs_amplitudes, 8)
!
      iteration = 0
      prev_energy = zero
!
      converged_energy = .false.
      converged_omega = .false.
      converged = .false.
!
      call mem%alloc(omega, wf%n_gs_amplitudes)
      call mem%alloc(dt, wf%n_gs_amplitudes)
      call mem%alloc(t, wf%n_gs_amplitudes)
!
      do while (.not. converged .and. iteration .le. solver%max_iterations) 
!
         iteration = iteration + 1
!
         write(output%unit, '(/t3,a)') 'Iteration    Energy (a.u.)        |omega|       Delta E (a.u.) '
         write(output%unit, '(t3,a)')  '---------------------------------------------------------------'
         flush(output%unit)
!
!        Construct Fock, calculate energy, and construct omega 
!
         call wf%construct_fock()
!
         call wf%calculate_energy()
         energy = wf%energy 
!
         call wf%construct_omega(omega)
!
         omega_norm = get_l2_norm(omega, wf%n_gs_amplitudes)
!
!        Print energy, energy difference and residual, then test convergence 
!
         write(output%unit, '(t3,i3,10x,f17.12,4x,e11.4,4x,e11.4)') iteration, wf%energy, &
                                          omega_norm, abs(energy-prev_energy)
         flush(output%unit)
!
         converged_energy   = abs(energy-prev_energy) .lt. solver%energy_threshold
         converged_omega    = omega_norm              .lt. solver%omega_threshold
!
         converged = converged_omega .and. converged_energy
!
         if (iteration .eq. 1 .and. converged_omega) converged = .true. ! Exception to the rule
!
!        If not converged, perform micro-iterations to get an estimate for the next amplitudes 
!
         if (.not. converged) then 
!
            solver%micro_residual_threshold = omega_norm*solver%relative_micro_residual_threshold
            call solver%do_micro_iterations(wf, omega, dt)
!
            call wf%get_amplitudes(t)
!
            call daxpy(wf%n_gs_amplitudes, one, dt, 1, t, 1)
!
            call diis_manager%update(omega, t)
!
            call wf%set_amplitudes(t)
!
!           Compute the new T1 transformed Cholesky vectors,
!           and store in memory the entire ERI-T1 matrix if possible and necessary 
!
            call wf%integrals%write_t1_cholesky(wf%t1)
            if (wf%need_g_abcd()) call wf%integrals%can_we_keep_g_pqrs_t1()
!
         endif 
!
         prev_energy = energy 
!
      enddo
!
      write(output%unit, '(t3,a)')  '---------------------------------------------------------------'
!
      if (.not. converged) then 
!   
         write(output%unit, '(/t3,a)')  'Warning: was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
!
      else
!
         call solver%print_summary(wf)
!
      endif 
!
      call mem%dealloc(omega, wf%n_gs_amplitudes)
      call mem%dealloc(dt, wf%n_gs_amplitudes)
      call mem%dealloc(t, wf%n_gs_amplitudes)
!
!
   end subroutine run_newton_raphson_cc_gs
!
!
   subroutine do_micro_iterations_newton_raphson_cc_gs(solver, wf, omega, dt)
!!
!!    Do micro iterations 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      use linear_davidson_tool_class
!
      implicit none 
!
      class(newton_raphson_cc_gs), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in)  :: omega 
      real(dp), dimension(wf%n_gs_amplitudes), intent(out) :: dt  
!
      real(dp), dimension(:), allocatable :: epsilon, first_trial, c_i
!
      real(dp) :: norm_trial, residual_norm, ddot
!
      type(linear_davidson_tool) :: davidson 
!
      integer :: micro_iteration 
!
      logical :: converged_residual 
!
!     Initialize solver tool and set preconditioner 
!
      call mem%alloc(epsilon, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(epsilon, wf%n_gs_amplitudes)
!
      call davidson%prepare('cc_gs_newton_raphson', wf%n_gs_amplitudes, solver%micro_residual_threshold, -omega)
!
      call davidson%set_preconditioner(epsilon)
      call mem%dealloc(epsilon, wf%n_gs_amplitudes)
!
!     Set start vector / initial guess 
!
!     Use - omega_mu / eps_mu as first guess 
!
      call mem%alloc(first_trial, wf%n_gs_amplitudes)
      call dcopy(wf%n_gs_amplitudes, omega, 1, first_trial, 1)
      call dscal(wf%n_gs_amplitudes, -one, first_trial, 1) 
!
      call davidson%precondition(first_trial)
!
      norm_trial = sqrt(ddot(wf%n_gs_amplitudes, first_trial, 1, first_trial, 1))
      call dscal(wf%n_gs_amplitudes, one/norm_trial, first_trial, 1)
!
      call davidson%write_trial(first_trial, 'rewind')
      call mem%dealloc(first_trial, wf%n_gs_amplitudes)
!
!     Enter iterative loop
!
      write(output%unit,'(/t3,a)') 'Micro-iteration     Residual norm'
      write(output%unit,'(t3,a)')  '---------------------------------'
      flush(output%unit)
!
      micro_iteration = 1
      converged_residual = .false.
!
      do while (.not. converged_residual .and. (micro_iteration .le. solver%max_micro_iterations))
!
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, davidson%n_parameters)
!
         call davidson%read_trial(c_i, davidson%dim_red)
         call wf%jacobian_transform_trial_vector(c_i) 
!
         if (micro_iteration == 1) then
!
            call davidson%write_transform(c_i, 'rewind')
!
         else
!
            call davidson%write_transform(c_i, 'append')
!
         endif
!
         call mem%dealloc(c_i, davidson%n_parameters)
!
!        Solve problem in reduced space
!
         call davidson%construct_reduced_matrix()
         call davidson%construct_reduced_gradient()
         call davidson%solve_reduced_problem()
!
!        Construct new trials and check if convergence criterion on residual is satisfied
!
         davidson%n_new_trials = 0
!
         call davidson%construct_next_trial_vec(residual_norm)
!
         write(output%unit,'(t3,i3,16x,e11.4)') micro_iteration, residual_norm
         flush(output%unit)
!
         converged_residual = .true.
!
         if (residual_norm .gt. solver%micro_residual_threshold) converged_residual = .false.
!   
         davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
         micro_iteration = micro_iteration + 1       
!
      enddo
!
      write(output%unit,'(t3,a/)')  '---------------------------------'
      flush(output%unit)
!
      if (.not. converged_residual) then
!
         write(output%unit, '(/t6,a)')  'Warning: was not able to converge the equations in the given'
         write(output%unit, '(t6,a/)')  'number of maximum micro-iterations.'
         flush(output%unit)
!
      endif
!
      call davidson%construct_X(dt, 1)
      call davidson%cleanup()
!
   end subroutine do_micro_iterations_newton_raphson_cc_gs
!
!
   subroutine cleanup_newton_raphson_cc_gs(wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%save_amplitudes()  
!
   end subroutine cleanup_newton_raphson_cc_gs
!
!
   subroutine print_banner_newton_raphson_cc_gs(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(newton_raphson_cc_gs) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description1,'(t3,a)',.false.)
!
   end subroutine print_banner_newton_raphson_cc_gs
!
!
   subroutine read_settings_newton_raphson_cc_gs(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(newton_raphson_cc_gs) :: solver    
!
      call input%get_keyword_in_section('omega threshold', 'solver cc gs', solver%omega_threshold)
      call input%get_keyword_in_section('energy threshold', 'solver cc gs', solver%energy_threshold)
      call input%get_keyword_in_section('diis dimension', 'solver cc gs', solver%diis_dimension)
      call input%get_keyword_in_section('max iterations', 'solver cc gs', solver%max_iterations)
      call input%get_keyword_in_section('rel micro threshold', 'solver cc gs', solver%relative_micro_residual_threshold)
      call input%get_keyword_in_section('max micro iterations', 'solver cc gs', solver%max_micro_iterations)
!
      if (input%requested_keyword_in_section('restart', 'solver cc gs')) solver%restart = .true.
!
   end subroutine read_settings_newton_raphson_cc_gs
!
!
   subroutine print_summary_newton_raphson_cc_gs(wf)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp) :: t1_diagnostic 
!
      write(output%unit, '(/t3,a)') '- DIIS accelerated Newton-Raphson CC ground state solver summary:'
!
      write(output%unit, '(/t6,a33,f18.12)') 'Final ground state energy (a.u.):', wf%energy 
      call wf%print_dominant_amplitudes()
!
      t1_diagnostic = wf%get_t1_diagnostic() 
      write(output%unit, '(/t6,a32,f14.12)') 'T1 diagnostic (|T1|/sqrt(N_e)): ', t1_diagnostic
!
   end subroutine print_summary_newton_raphson_cc_gs
!
!
end module newton_raphson_cc_gs_class
