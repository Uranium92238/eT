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
module newton_raphson_cc_gs_class
!
!!
!!	Newton-Raphson coupled cluster ground state solver class module
!!	Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!  
!! Solves the ground state coupled cluster equation (or, amplitude
!! equations)
!!
!!    Omega_mu = < mu | H-bar | HF > = 0,     H-bar = e-T H eT,
!!
!! for the cluster amplitudes t_mu in T = sum_mu t_mu tau_mu. The 
!! cluster amplitude give the (right) coupled cluster ground state as 
!!
!!    | CC > = e^T | HF >.
!!
!! The equation is solved using the direct inversion of the iterative 
!! subspace (DIIS) algorithm. See Pulay, P. Convergence acceleration 
!! of iterative sequences. The case of SCF iteration. Chem. Phys. Lett. 
!! 1980, 73, 393−398. This algorithm combines estimates for the parameters 
!! and the errors and finds a least-squares solution to the error being 
!! zero (see diis_tool solver tool for more details). 
!!
!! The amplitude estimates that are used to DIIS-extrapolate are the 
!! exact Newton-Raphson update estimates: 
!!
!!    t_mu <- t_mu + Delta t_mu = t_mu - (A^-1 Omega)_mu.
!!
!! The equation for the update estimate,
!!
!!    A (Delta t) = - Omega,
!!
!! is solved approximately in each macro-iteration (giving 
!! a number of micro-iterations, each involving a transformation
!! by the Jacobian matrix A). Micro-iterations are performed 
!! using the Davidson algorithm for linear CC equations - see 
!! davidson_cc_multipliers solver for more details.
!!
!
   use kinds
   use ccs_class
   use diis_tool_class

   use abstract_convergence_tool_class,   only : abstract_convergence_tool
   use convergence_tool_class,            only : convergence_tool
!
   implicit none
!
   type :: newton_raphson_cc_gs
!
      character(len=100) :: name_ = 'DIIS Newton-Raphson coupled cluster ground state solver'
!
      character(len=500) :: description1 = 'A Newton-Raphson CC ground state equations solver. Solves the &
                                             &ground state equation using updates Δt based on the Newton equation, &
                                             &A Δt = -Ω, where the A is the Jacobian and Ω the omega vector.'
!
      integer :: max_iterations
      integer :: diis_dimension
!
      integer  :: max_micro_iterations
      integer  :: max_micro_dim_red
      real(dp) :: micro_residual_threshold
      real(dp) :: relative_micro_residual_threshold
!
      logical :: crop ! Standard DIIS if false; CROP variant of DIIS if true
!
      character(len=200) :: storage 
      logical :: restart, records_in_memory 
!
      type(timings), allocatable :: timer 
!
      class(abstract_convergence_tool), allocatable :: convergence_checker
!
   contains
!     
      procedure :: cleanup                  => cleanup_newton_raphson_cc_gs
      procedure :: run                      => run_newton_raphson_cc_gs
!
      procedure :: do_micro_iterations      => do_micro_iterations_newton_raphson_cc_gs
!
      procedure :: read_settings            => read_settings_newton_raphson_cc_gs
      procedure :: print_banner             => print_banner_newton_raphson_cc_gs
      procedure :: print_settings           => print_settings_newton_raphson_cc_gs
!
   end type newton_raphson_cc_gs
!
!
   interface newton_raphson_cc_gs
!
      procedure :: new_newton_raphson_cc_gs
!
   end interface newton_raphson_cc_gs
!
!
contains
!
!
   function new_newton_raphson_cc_gs(wf, restart) result(solver)
!!
!!    New Newton-Rahpson GS 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(newton_raphson_cc_gs) :: solver
!
      class(ccs) :: wf
!
      logical, intent(in) :: restart
!
      solver%timer = timings('Newton-Raphson CC GS solver time', pl='minimal')
      call solver%timer%turn_on()
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set standard settings 
!
      solver%max_iterations                     = 100
      solver%diis_dimension                     = 8
      solver%max_micro_iterations               = 100
      solver%relative_micro_residual_threshold  = 1.0d-2
      solver%restart                            = restart
      solver%max_micro_dim_red                  = 50
      solver%storage                            = 'disk'
      solver%crop                               = .false.
!
!     Initialize convergence checker with default threshols
!
      solver%convergence_checker = convergence_tool(energy_threshold   = 1.0d-5,   &
                                                    residual_threshold = 1.0d-5,   &
                                                    energy_convergence = .false.)
!
!     Read & print settings (thresholds, etc.)
!
      call solver%read_settings()
!
      call solver%print_settings()
!
!     Set the amplitudes to the initial guess or read if restart
!
      call wf%initialize_amplitudes()
!
      call wf%eri%set_t1_to_mo()
      call wf%eri%place_g_mo_in_memory()
!
      call wf%set_initial_amplitudes_guess(solver%restart)
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
   end function new_newton_raphson_cc_gs
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
      call output%printf('m', '- DIIS accelerated Newton-Raphson CC ground &
                         &state solver settings:', fs='(/t3,a)')
!
      call solver%convergence_checker%print_settings()
!
      call output%printf('m', 'Relative micro threshold: (e14.3)', &
                         reals=[solver%relative_micro_residual_threshold], fs='(t6,a)')
!
      call output%printf('m', 'DIIS dimension:           (i14)', &
                         ints=[solver%diis_dimension], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations: (i14)', &
                         ints=[solver%max_iterations], fs='(t6,a)')
      call output%printf('m', 'Max number of micro-iterations: (i8)', &
                         ints=[solver%max_micro_iterations], fs='(t6,a)')
!
      if (solver%crop) then 
!
         call output%printf('m', 'Enabled CROP in the DIIS algorithm.', fs='(/t6,a)')
!
      endif
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
      type(diis_tool) :: diis
!
      real(dp), dimension(:), allocatable :: omega, dt, t  
!
      integer :: iteration, micro_iterations
!
      logical :: converged
!
      real(dp) :: energy, prev_energy, omega_norm
!
      type(timings), allocatable :: macro_iteration_timer
!
      diis = diis_tool('cc_gs_diis',                        &
                        wf%n_gs_amplitudes,                 &
                        wf%n_gs_amplitudes,                 &
                        dimension_=solver%diis_dimension,   &
                        crop=solver%crop)
!
      call diis%initialize_storers(solver%records_in_memory)
!
      iteration = 0
      micro_iterations = 0
      prev_energy = zero
!
      converged = .false.
!
      macro_iteration_timer = timings('Newton-Raphson macro-iteration time ' // &
                                       '(includes micro-iteration time)', pl='normal')
!
      call mem%alloc(omega, wf%n_gs_amplitudes)
      call mem%alloc(dt, wf%n_gs_amplitudes)
      call mem%alloc(t, wf%n_gs_amplitudes)
!
      do while (.not. converged .and. iteration .le. solver%max_iterations) 
!
         iteration = iteration + 1
         call macro_iteration_timer%turn_on()
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
         call output%printf('n', 'Macro-iter.    Energy (a.u.)        |omega|   &
                            &    Delta E (a.u.)', fs='(/t3,a)')
!
         call output%print_separator('n', 64, '-', fs='(t3,a)')
!
         call output%printf('n', '(i5)         (f17.12)    (e11.4)    (e11.4)', &
                            ints=[iteration], reals=[wf%energy, omega_norm, abs(energy-prev_energy)])
!
         call output%print_separator('n', 64, '-', fs='(t3,a)')
!
         converged = solver%convergence_checker%has_converged(omega_norm, energy-prev_energy, iteration)
!
!        If not converged, perform micro-iterations to get an estimate for the next amplitudes 
!
         if (.not. converged) then 
!
            solver%micro_residual_threshold = omega_norm*solver%relative_micro_residual_threshold
            call solver%do_micro_iterations(wf, omega, dt, micro_iterations)
!
            call wf%get_amplitudes(t)
!
            call wf%form_newton_raphson_t_estimate(t, dt)
!
            call diis%update(omega, t)
!
            call wf%set_amplitudes(t)
!
!           Update the Cholesky (and electron repulsion integrals, if in memory) 
!           to new T1 amplitudes 
!
            call wf%eri%update_t1_integrals(wf%t1)
!
         endif 
!
         prev_energy = energy 
!
         call wf%save_amplitudes()
!
         call macro_iteration_timer%turn_off()
         call macro_iteration_timer%reset()
!
      enddo
!
      if (.not. converged) then 
!
         call output%warning_msg('Was not able to converge the equations      &
                                 &in the given number of maximum iterations.')
!
      else
!
         call output%printf('m', 'Convergence criterion met in (i0) iterations!', &
                            ints=[iteration], fs='(t3,a)')
!
         call wf%print_gs_summary()
!
      endif 
!
      call mem%dealloc(omega, wf%n_gs_amplitudes)
      call mem%dealloc(dt, wf%n_gs_amplitudes)
      call mem%dealloc(t, wf%n_gs_amplitudes)
!
      call diis%finalize_storers()
!
   end subroutine run_newton_raphson_cc_gs
!
!
   subroutine do_micro_iterations_newton_raphson_cc_gs(solver, wf, omega, dt, final_iteration)
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
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in)  :: omega 
      real(dp), dimension(wf%n_gs_amplitudes), intent(out) :: dt  
!
      integer, intent(out) :: final_iteration
!
      real(dp), dimension(:), allocatable :: preconditioner, c, residual
!
      real(dp), dimension(:), allocatable :: minus_omega  
!
      real(dp) :: residual_norm
!
      type(linear_davidson_tool) :: davidson 
!
      integer :: micro_iteration 
!
      logical :: converged_residual 
!
      type(timings), allocatable :: micro_iteration_timer, timer 
!
      timer = timings('Newton-Raphson total micro-iteration time', pl='normal')
      call timer%turn_on()
!
!     Initialize solver tool and set preconditioner 
!
      call mem%alloc(preconditioner, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(preconditioner, wf%n_gs_amplitudes)
!
      call mem%alloc(minus_omega, wf%n_gs_amplitudes)
      call copy_and_scale(-one, omega, minus_omega, wf%n_gs_amplitudes)
!
      davidson = linear_davidson_tool('cc_gs_newton_raphson',                          &
                                       wf%n_gs_amplitudes,                             &
                                       min(1.0d-11, solver%micro_residual_threshold),  &
                                       solver%max_micro_dim_red,                       &
                                       minus_omega, 1)
!
      call davidson%initialize(solver%records_in_memory)
!
      call davidson%set_preconditioner(preconditioner)
      call mem%dealloc(preconditioner, wf%n_gs_amplitudes)
!
!     Set start vector / initial guess 
!
!     Use - omega_mu / eps_mu as first guess 
!
      call davidson%set_trials_to_preconditioner_guess()
!
!     Prepare intermediates for Jacobian transformation
!
      call wf%prepare_for_jacobian()
!
!     Enter iterative loop
!
      call output%printf('n', 'Micro-iter.  Residual norm', fs='(/t6,a)')
      call output%print_separator('n', 26, '-', fs='(t6,a)')
!
      micro_iteration = 0
      converged_residual = .false.
!
      call mem%alloc(c, davidson%n_parameters)
      call mem%alloc(residual, davidson%n_parameters)
!
      micro_iteration_timer = timings('Newton-Raphson micro-iteration time', pl='verbose')
!
      do while (.not. converged_residual .and. (micro_iteration .le. solver%max_micro_iterations))
!
         micro_iteration = micro_iteration + 1
         call micro_iteration_timer%turn_on()
!
!        Reduced space preparations 
!
         if (davidson%red_dim_exceeds_max()) call davidson%set_trials_to_solutions()
!
         call davidson%update_reduced_dim()
!
         call davidson%orthonormalize_trial_vecs() 
!
!        Transform new trial vectors and write to file
!
         call davidson%get_trial(c, davidson%dim_red)
         call wf%construct_Jacobian_transform('right', c, w = zero)
         call davidson%set_transform(c, davidson%dim_red)
!
!        Solve problem in reduced space
!
         call davidson%solve_reduced_problem()
!
!        Construct new trials and check if convergence criterion on residual is satisfied
!
         call davidson%construct_residual(residual, 1)
         residual_norm = get_l2_norm(residual, wf%n_gs_amplitudes)
!
         converged_residual = .true.
         if (residual_norm >= solver%micro_residual_threshold) then 
!
            converged_residual = .false.
            call davidson%add_new_trial(residual, 1)
!
         endif 
!
         call output%printf('n', '(i3)          (e11.4)', &
                            ints=[micro_iteration], reals=[residual_norm], fs='(t6,a)')
!
         call micro_iteration_timer%turn_off()
         call micro_iteration_timer%reset()
!
      enddo
!
      call mem%dealloc(residual, davidson%n_parameters)
      call mem%dealloc(c, davidson%n_parameters)
!
      call output%print_separator('n', 26, '-', fs='(t6,a)')
!
      if (.not. converged_residual) then
!
         call output%error_msg('was not able to converge the equations in the given ' // &
                                'max number of iterations in do_micro_iterations_newton_raphson_cc_gs.')
!
      endif
!
      call davidson%construct_solution(dt, 1)
      call davidson%cleanup()
!
      final_iteration = micro_iteration
!
      call mem%dealloc(minus_omega, wf%n_gs_amplitudes)
!
      call timer%turn_off()
!
   end subroutine do_micro_iterations_newton_raphson_cc_gs
!
!
   subroutine cleanup_newton_raphson_cc_gs(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(newton_raphson_cc_gs) :: solver 
!
      class(ccs) :: wf
!
      call output%printf('m', '- Finished solving the (a0) ground state equations', &
                         chars=[convert_to_uppercase(wf%name_)], fs='(/t3, a)')
!
      call wf%save_amplitudes()
!
      call solver%timer%turn_off() 
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
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('m', solver%description1, ffs='(/t3,a)')
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
      real(dp) :: energy_threshold, omega_threshold
!
      if (input%requested_keyword_in_section('energy threshold', 'solver cc gs')) then
!
         call input%get_keyword_in_section('energy threshold', 'solver cc gs', energy_threshold)
         call solver%convergence_checker%set_energy_threshold(energy_threshold)
!
      endif
!
      if (input%requested_keyword_in_section('omega threshold', 'solver cc gs')) then
!
         call input%get_keyword_in_section('omega threshold', 'solver cc gs', omega_threshold)
         call solver%convergence_checker%set_residual_threshold(omega_threshold)
!
      endif
      call input%get_keyword_in_section('diis dimension', 'solver cc gs', solver%diis_dimension)
      call input%get_keyword_in_section('max iterations', 'solver cc gs', solver%max_iterations)
      call input%get_keyword_in_section('rel micro threshold', 'solver cc gs', solver%relative_micro_residual_threshold)
      call input%get_keyword_in_section('max micro iterations', 'solver cc gs', solver%max_micro_iterations)
!
      call input%get_keyword_in_section('storage', 'solver cc gs', solver%storage)
!
      if (input%requested_keyword_in_section('crop', 'solver cc gs')) then 
!
         solver%crop = .true.
!
      endif
!
   end subroutine read_settings_newton_raphson_cc_gs
!
!
end module newton_raphson_cc_gs_class
