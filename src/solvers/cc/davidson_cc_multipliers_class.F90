module davidson_cc_multipliers_class
!
!!
!!    Davidson coupled cluster multipliers solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!  
!
   use kinds
   use file_class
   use ccs_class
   use linear_davidson_tool_class
!
   implicit none
!
   type :: davidson_cc_multipliers
!
      character(len=100) :: tag = 'Davidson coupled cluster multipliers solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
      character(len=500) :: description = 'A Davidson CC multiplier equations solver.'
!
      integer :: max_iterations
!
      real(dp) :: residual_threshold
!
      logical :: restart
!
   contains
!
      procedure :: prepare                         => prepare_davidson_cc_multipliers
      procedure, nopass :: cleanup                 => cleanup_davidson_cc_multipliers
!
      procedure :: print_banner                    => print_banner_davidson_cc_multipliers
      procedure :: print_settings                  => print_settings_davidson_cc_multipliers
!
      procedure, nopass :: print_summary           => print_summary_davidson_cc_multipliers
!
      procedure :: read_settings                   => read_settings_davidson_cc_multipliers
!
      procedure :: run                             => run_davidson_cc_multipliers
!
      procedure, nopass :: set_precondition_vector => set_precondition_vector_davidson_cc_multipliers     
!
      procedure, nopass :: transform_trial_vector  => transform_trial_vector_davidson_cc_multipliers 
!
   end type davidson_cc_multipliers
!
!
contains
!
!
   subroutine prepare_davidson_cc_multipliers(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set default settings
!
      solver%max_iterations = 20
!
      solver%residual_threshold  = 1.0d-6
!
      solver%restart = .false.
!
      call solver%read_settings()
!
      call solver%print_settings()
!
      call wf%construct_fock()
!
      call wf%initialize_multipliers()
!
   end subroutine prepare_davidson_cc_multipliers
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
      write(output%unit, '(/t3,a)')      '- Davidson CC multipliers solver settings:'
!
      write(output%unit, '(/t6,a26,e9.2)') 'Residual threshold:       ', solver%residual_threshold
!
      write(output%unit, '(t6,a26,i9)')    'Max number of iterations: ', solver%max_iterations
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
      real(dp), dimension(:), allocatable :: eta, c_i, multipliers
!
      integer :: iteration
!
      real(dp) :: residual_norm, ddot, norm_trial
!
!     Set eta
!
      call mem%alloc(eta, wf%n_gs_amplitudes)
      call wf%construct_eta(eta)
      call dscal(wf%n_gs_amplitudes, -one, eta, 1)
!
!     Initialize solver tool
!
      call davidson%prepare('multipliers', wf%n_gs_amplitudes, solver%residual_threshold, eta)
!
!     Set start vector
!
!     Precondition eta and use as trial
!
      call solver%set_precondition_vector(wf, davidson)
      call davidson%precondition(eta)
!
      norm_trial = sqrt(ddot(wf%n_gs_amplitudes, eta, 1, eta, 1))
      call dscal(wf%n_gs_amplitudes, one/norm_trial, eta, 1)
!
      call davidson%write_trial(eta, 'rewind')
!
      call mem%dealloc(eta, wf%n_gs_amplitudes)
!
!     Enter iterative loop
!
      iteration = 1
!
      write(output%unit,'(/t3,a)') 'Iteration     Residual norm'
      write(output%unit,'(t3,a)')  '---------------------------'
      flush(output%unit)
!
      converged_residual = .false.
!
      do while (.not. converged_residual .and. (iteration .le. solver%max_iterations))
!
!        Transform new trial vectors and write to file
!
         call mem%alloc(c_i, davidson%n_parameters)
!
         call davidson%read_trial(c_i, davidson%dim_red)
         call solver%transform_trial_vector(wf, c_i)
!
         if (iteration == 1) then
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
         write(output%unit,'(t3,i3,12x,e11.4)') iteration, residual_norm
         flush(output%unit)
!
         converged_residual = .true.
!
         if (residual_norm .gt. solver%residual_threshold) converged_residual = .false.
!
         if (davidson%dim_red .ge. davidson%max_dim_red) then
!
           call davidson%set_trials_to_solutions()
!
         else
!
            davidson%dim_red = davidson%dim_red + davidson%n_new_trials
!
         endif
!
         iteration = iteration + 1       
!
      enddo
!
      write(output%unit,'(t3,a)') '---------------------------'
      flush(output%unit)
!
      if (converged_residual) then
!
         call solver%print_summary(davidson, wf)
!
      else
!
         write(output%unit, '(/t3,a)')  'Warning: was not able to converge the equations in the given'
         write(output%unit, '(t3,a/)')  'number of maximum iterations.'
         flush(output%unit)
!
      endif
!
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
!
      call davidson%construct_X(multipliers, 1)
!
      call wf%set_multipliers(multipliers)
!
      call mem%dealloc(multipliers, wf%n_gs_amplitudes)
!
   end subroutine run_davidson_cc_multipliers
!
!
   subroutine cleanup_davidson_cc_multipliers(wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      call wf%save_multipliers()  
      call wf%destruct_multipliers()
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
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
!
   end subroutine print_banner_davidson_cc_multipliers
!
!
   subroutine transform_trial_vector_davidson_cc_multipliers(wf, c_i)
!!
!!    Transform trial vector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Transforms the trial vector according to specified transformation routine.
!!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: c_i
!
      call wf%jacobian_transpose_transform_trial_vector(c_i) 
!
   end subroutine transform_trial_vector_davidson_cc_multipliers
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
   subroutine print_summary_davidson_cc_multipliers(davidson, wf)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(linear_davidson_tool) :: davidson
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(:), allocatable :: X
!
      write(output%unit, '(/t3,a)') '- Multipliers vector amplitudes:'
!
      call mem%alloc(X, wf%n_gs_amplitudes)
!
      call davidson%construct_X(X, 1)         
!
      call wf%print_dominant_x_amplitudes(X, 'r')
!
   end subroutine print_summary_davidson_cc_multipliers
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
      integer :: n_specs, i
      character(len=100) :: line
!
      if (.not. requested_section('multipliers')) return
!
      call move_to_section('multipliers', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:10) == 'threshold:' ) then
!
            read(line(11:100), *) solver%residual_threshold
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         endif
!
      enddo
!
   end subroutine read_settings_davidson_cc_multipliers
!
!
end module davidson_cc_multipliers_class
