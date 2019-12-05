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
      real(dp) :: omega_threshold
      real(dp) :: energy_threshold
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
   contains
!     
      procedure, nopass :: cleanup          => cleanup_newton_raphson_cc_gs
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
      solver%energy_threshold                   = 1.0d-6
      solver%omega_threshold                    = 1.0d-6
      solver%restart                            = restart
      solver%max_micro_dim_red                  = 50
      solver%storage                            = 'disk'
      solver%crop                               = .false.
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
!
         call wf%integrals%write_t1_cholesky(wf%t1) 
!
         if(wf%need_g_abcd .and. wf%integrals%room_for_g_pqrs_t1()) &
            call wf%integrals%place_g_pqrs_t1_in_memory()
! 
      else
!
         call wf%integrals%write_t1_cholesky(wf%t1) 
!
         if(wf%need_g_abcd .and. wf%integrals%room_for_g_pqrs_t1()) &
            call wf%integrals%place_g_pqrs_t1_in_memory()
!
         call wf%set_initial_amplitudes_guess()
!
      endif
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
      call output%printf('- DIIS accelerated Newton-Raphson CC ground state solver settings:', &
                          pl='m', fs='(/t3,a)')
!
      call output%printf('Omega threshold:          (e14.3)', reals=[solver%omega_threshold], pl='m', fs='(/t6,a)')
      call output%printf('Energy threshold:         (e14.3)', reals=[solver%energy_threshold], pl='m', fs='(t6,a)')
      call output%printf('Relative micro threshold: (e14.3)',               &
                          reals=[solver%relative_micro_residual_threshold], &
                          pl='m', fs='(t6,a)')
!
      call output%printf('DIIS dimension:           (i14)', ints=[solver%diis_dimension], &
                          pl='m', fs='(/t6,a)')
      call output%printf('Max number of iterations: (i14)', ints=[solver%max_iterations], &
                          pl='m', fs='(t6,a)')
      call output%printf('Max number of micro-iterations: (i8)', &
                          ints=[solver%max_micro_iterations],    &
                          pl='m', fs='(t6,a)')
!
      if (solver%crop) then 
!
         call output%printf('Enabled CROP in the DIIS algorithm.', pl='minimal', fs='(/t6,a)')
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
      logical :: converged, converged_omega, converged_energy
!
      real(dp) :: energy, prev_energy, omega_norm
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
         call output%printf('Macro-iter.    Energy (a.u.)        |omega|       Delta E (a.u.)', &
                             fs='(/t3,a)', pl='n')
!
         call output%print_separator('n', 64, '-', fs='(t3,a)')
!
         call output%printf('(i5)         (f17.12)    (e11.4)    (e11.4)',           &
                             ints=[iteration],                                       &
                             reals=[wf%energy, omega_norm, abs(energy-prev_energy)], &
                             pl='n')
!
         call output%print_separator('n', 64, '-', fs='(t3,a)')
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
!           Compute the new T1 transformed Cholesky vectors,
!           and store in memory the entire ERI-T1 matrix if possible and necessary 
!
            call wf%integrals%write_t1_cholesky(wf%t1)
            if (wf%integrals%get_eri_t1_mem()) &
               call wf%integrals%update_g_pqrs_t1_in_memory()
!
         endif 
!
         prev_energy = energy 
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
         call output%printf('Convergence criterion met in (i0) iterations!', &
                             ints=[iteration], pl='m', fs='(t3,a)')
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
!     Initialize solver tool and set preconditioner 
!
      call mem%alloc(preconditioner, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(preconditioner, wf%n_gs_amplitudes)
!
      call mem%alloc(minus_omega, wf%n_gs_amplitudes)
      call copy_and_scale(-one, omega, minus_omega, wf%n_gs_amplitudes)
!
      davidson = linear_davidson_tool('cc_gs_newton_raphson', wf%n_gs_amplitudes, &
         solver%micro_residual_threshold, solver%max_micro_dim_red, minus_omega, 1)
!
      call davidson%initialize_trials_and_transforms(solver%records_in_memory)
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
      call output%printf('Micro-iter.  Residual norm', pl='n', fs='(/t6,a)')
      call output%print_separator('n', 26, '-', fs='(t6,a)')
!
      micro_iteration = 0
      converged_residual = .false.
!
      call mem%alloc(c, davidson%n_parameters)
      call mem%alloc(residual, davidson%n_parameters)
!
      do while (.not. converged_residual .and. (micro_iteration .le. solver%max_micro_iterations))
!
         micro_iteration = micro_iteration + 1
         call davidson%iterate()
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
            call davidson%construct_next_trial(residual, 1)
!
         endif 
!
         call output%printf('(i3)          (e11.4)', pl='n', ints=[micro_iteration], &
                            reals=[residual_norm], fs='(t6,a)')
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
      call davidson%finalize_trials_and_transforms()
!
      final_iteration = micro_iteration
!
      call mem%dealloc(minus_omega, wf%n_gs_amplitudes)
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
      call output%printf('- Finished solving the (a0) ground state equations', & 
                        chars=[convert_to_uppercase(wf%name_)], &
                        fs='(/t3, a)', pl='minimal')
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
      call output%printf(' - ' // trim(solver%name_), pl='m', fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf(solver%description1, pl='m', ffs='(/t3,a)')
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
