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
module diis_cc_multipliers_class
!
!!
!! DIIS coupled cluster multipliers solver class module
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
!! The equation is solved using the direct inversion of the iterative 
!! subspace (DIIS) algorithm. See Pulay, P. Convergence acceleration 
!! of iterative sequences. The case of SCF iteration. Chem. Phys. Lett. 
!! 1980, 73, 393−398. This algorithm combines estimates for the parameters 
!! and the errors and finds a least-squares solution to the error being 
!! zero (see diis_tool solver tool for more details). 
!!
!! The update estimates used in DIIS extrapolation are the ones 
!! resulting from the orbital differences approximation of A 
!! (See davidson_cc_es solver for more details.),
!!
!!    tbar_mu <- tbar_mu - R_mu/epsilon_mu,
!!
!! where R_mu = (tbar^T A + eta^T)_mu is the residual vector.
!!
!
   use kinds
   use parameters 
   use global_out, only: output
   use timings_class, only: timings
   use global_in, only: input
   use array_utilities, only: get_l2_norm
   use string_utilities, only: convert_to_uppercase
   use ccs_class, only: ccs 
   use diis_tool_class, only: diis_tool 
   use precondition_tool_class, only: precondition_tool 
   use memory_manager_class, only: mem
!
   implicit none
!
   type :: diis_cc_multipliers
!
      character(len=100) :: tag = 'DIIS multipliers solver'
      character(len=100) :: author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
!
      character(len=500) :: description1 = 'A DIIS CC multiplier equations solver. It combines a quasi-Newton &
                                           &perturbation theory estimate of the next multipliers, using &
                                           &least square fitting to find an an optimal combination of &  
                                           &previous estimates such that the update is minimized.'
!
      character(len=500) :: description2 = 'See Helgaker et al., Molecular Electronic Structure Theory, &
                                           &Chapter 13, for the more details on this algorithm.'
!
      integer :: diis_dimension
!
      integer :: max_iterations
!
      real(dp) :: residual_threshold
!
      character(len=200) :: storage 
      logical :: restart, records_in_memory 
!
      type(timings) :: timer
!
      class(precondition_tool), allocatable :: preconditioner 
!
   contains
!     
      procedure :: run                      => run_diis_cc_multipliers
      procedure :: cleanup                  => cleanup_diis_cc_multipliers
!
      procedure :: read_settings            => read_settings_diis_cc_multipliers
!
      procedure :: print_banner             => print_banner_diis_cc_multipliers
      procedure, nopass :: print_summary    => print_summary_diis_cc_multipliers
!
      procedure :: print_settings           => print_settings_diis_cc_multipliers
!
   end type diis_cc_multipliers
!
!
   interface diis_cc_multipliers
!
      procedure :: new_diis_cc_multipliers
!
   end interface diis_cc_multipliers
!
!
contains
!
!
   function new_diis_cc_multipliers(wf) result(solver)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(diis_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable :: eps
!
      solver%timer = timings(trim(convert_to_uppercase(wf%name_)) // ' multipliers')
      call solver%timer%turn_on()
!
!     Print solver banner
!
      call solver%print_banner()
!
!     Set default settings
!
      solver%diis_dimension = 8
      solver%max_iterations = 100
      solver%residual_threshold  = 1.0d-6
      solver%restart = .false.
      solver%storage = 'disk'
!
      call solver%read_settings()
!
      call solver%print_settings()
!
      call wf%construct_fock()
!
      call wf%initialize_multipliers()
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
!     Initialize preconditioner 
!
      call mem%alloc(eps, wf%n_gs_amplitudes)
      call wf%get_gs_orbital_differences(eps, wf%n_gs_amplitudes)
!
      solver%preconditioner = precondition_tool(eps, wf%n_gs_amplitudes)
!
      call mem%dealloc(eps, wf%n_gs_amplitudes)
!
   end function new_diis_cc_multipliers
!
!
   subroutine print_settings_diis_cc_multipliers(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(diis_cc_multipliers) :: solver 
!
      call output%printf('- DIIS CC multipliers solver settings:', pl='m', fs='(/t3,a)')
!
      call output%printf('Residual threshold: (e15.2)', pl='m', fs='(/t6,a)', &
                          reals=[solver%residual_threshold])

      call output%printf('DIIS dimension: (i19)', pl='m', fs='(/t6,a)', &
                          ints=[solver%diis_dimension])
      call output%printf('Max number of iterations: (i9)', pl='m', fs='(t6,a)', &
                          ints=[solver%max_iterations])
!
   end subroutine print_settings_diis_cc_multipliers
!
!
   subroutine run_diis_cc_multipliers(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
      type(diis_tool) :: diis
!
      logical :: converged_residual 
!
      real(dp) :: residual_norm
!
      real(dp), dimension(:), allocatable :: residual  
      real(dp), dimension(:), allocatable :: multipliers
      real(dp), dimension(:), allocatable :: epsilon  
!
      integer :: iteration
!
      call wf%prepare_for_multiplier_equation
!
      diis = diis_tool('cc_multipliers_diis', wf%n_gs_amplitudes, wf%n_gs_amplitudes, &
                  solver%records_in_memory, dimension_=solver%diis_dimension)
!
      call mem%alloc(residual, wf%n_gs_amplitudes)
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
      call mem%alloc(epsilon, wf%n_gs_amplitudes)
!
      call wf%get_gs_orbital_differences(epsilon, wf%n_gs_amplitudes)
!
      if (solver%restart) then 
!
         call output%printf('Requested restart. Reading multipliers from file.', &
                             pl='m', fs='(/t3,a)')
!
         call wf%is_restart_safe('ground state')
         call wf%read_multipliers()
         call wf%get_multipliers(multipliers) 
!
      else
!
         multipliers = zero
         call wf%set_multipliers(multipliers)
!
      endif
!
      converged_residual = .false.
!
      call output%printf('Iteration    Norm residual  ', pl='n', fs='(/t3,a)')
      call output%print_separator('n', 28,'-', fs='(t3,a)')
!
      iteration   = 1
!
      do while (.not. converged_residual .and. iteration .le. solver%max_iterations)         
!
!        Construct the multiplier equations 
!
         call wf%construct_multiplier_equation(residual)
         residual_norm = get_l2_norm(residual, wf%n_gs_amplitudes)
!
         call output%printf('(i3)         (e11.4)', ints=[iteration], reals=[residual_norm], pl='n')
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         converged_residual = residual_norm .lt. solver%residual_threshold
!
         if (converged_residual) then
!
            call output%print_separator('n', 28,'-', fs='(t3,a)')
            call output%printf('Convergence criterion met in (i0) iterations!', ints=[iteration], pl='n', fs='(/t3,a)')
!
         else
!
!           Precondition residual, shift multipliers by preconditioned residual, 
!           then ask for the DIIS update of the multipliers 
!
            call solver%preconditioner%do_(residual)
!
            call wf%get_multipliers(multipliers)
            call daxpy(wf%n_gs_amplitudes, one, residual, 1, multipliers, 1)
!
            call diis%update(residual, multipliers)
            call wf%set_multipliers(multipliers)
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      if (.not. converged_residual) then 
!   
         call output%print_separator('m', 63,'-', fs='(t3,a)')
         call output%warning_msg('was not able to converge the equations      &
                                 &in the given number of maximum iterations.')
!
      else
!
         call solver%print_summary(wf, multipliers)
!
      endif 
!
      call mem%dealloc(residual, wf%n_gs_amplitudes)
      call mem%dealloc(multipliers, wf%n_gs_amplitudes)
      call mem%dealloc(epsilon, wf%n_gs_amplitudes)
!
   end subroutine run_diis_cc_multipliers
!
!
   subroutine cleanup_diis_cc_multipliers(solver, wf)
!!
!! 	Cleanup 
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(diis_cc_multipliers) :: solver
      class(ccs) :: wf
!
      call wf%save_multipliers()
!
      call solver%timer%turn_off()
!
      call output%printf('- Finished solving the ' // trim(wf%name_) // &
                         ' multipliers equations', pl='n', fs='(/t3, a)')
!
      call output%printf('Total wall time (sec): (f20.5)',              &
                          reals=[solver%timer%get_elapsed_time('wall')], &
                          pl='n', fs='(/t6,a)')
      call output%printf('Total cpu time (sec):  (f20.5)',             &
                          reals=[solver%timer%get_elapsed_time('cpu')], &
                          pl='n', fs='(t6,a)')
!
   end subroutine cleanup_diis_cc_multipliers
!
!
   subroutine print_banner_diis_cc_multipliers(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(diis_cc_multipliers) :: solver 
!
      call output%long_string_print(solver%tag,'(//t3,a)',.true.)
      call output%long_string_print(solver%author,'(t3,a/)',.true.)
      call output%long_string_print(solver%description1,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
      call output%long_string_print(solver%description2)
!
   end subroutine print_banner_diis_cc_multipliers
!
!
   subroutine print_summary_diis_cc_multipliers(wf, X)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: X
!
      call output%printf('- DIIS CC multipliers solver summary:', pl='m', fs='(/t3,a)')
!
      call wf%print_dominant_x_amplitudes(X, 'r')
!
   end subroutine print_summary_diis_cc_multipliers
!
!
   subroutine read_settings_diis_cc_multipliers(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cc_multipliers) :: solver 
!
      call input%get_keyword_in_section('threshold', 'solver cc multipliers', solver%residual_threshold)
      call input%get_keyword_in_section('max iterations', 'solver cc multipliers', solver%max_iterations)
!
      if (input%requested_keyword_in_section('restart', 'solver cc multipliers')) solver%restart = .true.    
!
      call input%get_keyword_in_section('storage', 'solver cc multipliers', solver%storage)
!
   end subroutine read_settings_diis_cc_multipliers
!
!
end module diis_cc_multipliers_class
!
