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
   use parameters
   use ccs_class, only: ccs
   use global_out, only: output
   use timings_class, only: timings
   use diis_tool_class, only: diis_tool
   use precondition_tool_class, only: precondition_tool
   use memory_manager_class, only: mem
   use amplitude_updater_class, only: amplitude_updater
!
   implicit none
!
   type :: diis_cc_multipliers
!
      character(len=100) :: name_ = 'DIIS coupled cluster multipliers solver'
!
      character(len=500) :: description1 = 'A DIIS CC multiplier equations solver. It combines a quasi-Newton &
                                           &perturbation theory estimate of the next multipliers, using &
                                           &least square fitting to find an an optimal combination of &
                                           &previous estimates such that the update is minimized.'
!
      character(len=500) :: description2 = 'See Helgaker et al., Molecular Electronic Structure Theory, &
                                           &Chapter 13, for the more details on this algorithm.'
!
      integer :: max_iterations
!
      real(dp) :: residual_threshold
!
      logical :: restart
!
      type(timings) :: timer
!
      class(amplitude_updater), allocatable :: tbar_updater
!
      class(precondition_tool), allocatable :: preconditioner
!
      type(diis_tool), allocatable :: diis
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
   function new_diis_cc_multipliers(wf, restart, tbar_updater) result(solver)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      type(diis_cc_multipliers) :: solver
!
      class(ccs) :: wf
!
      logical, intent(in) :: restart
!
      class(amplitude_updater), intent(in) :: tbar_updater
!
      logical :: records_in_memory, crop
      integer :: diis_dimension
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
      solver%max_iterations      = 100
      solver%residual_threshold  = 1.0d-5
      solver%restart             = restart
!
      records_in_memory   = .false.
      diis_dimension      = 8
      crop                = .false.
!
      call solver%read_settings(records_in_memory, crop, diis_dimension)
      call solver%print_settings()
!
      solver%tbar_updater = tbar_updater
!
      solver%diis = diis_tool('cc_multipliers_diis',      &
                               wf%n_gs_amplitudes,        &
                               wf%n_gs_amplitudes,        &
                               dimension_=diis_dimension, &
                               crop=crop,                 &
                               records_in_memory=records_in_memory)
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
      call output%printf('m', '- DIIS CC multipliers solver settings:', fs='(/t3,a)')
!
      call output%printf('m', 'Residual threshold:       (e9.2)', &
                         reals=[solver%residual_threshold], fs='(/t6,a)')
!
      call output%printf('m', 'Max number of iterations: (i9)', &
                         ints=[solver%max_iterations], fs='(t6,a)')
!
   end subroutine print_settings_diis_cc_multipliers
!
!
   subroutine run_diis_cc_multipliers(solver, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(diis_cc_multipliers) :: solver
!
      class(ccs) :: wf
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
      call wf%initialize_multipliers()
      call wf%prepare_for_multiplier_equation()
!
      call solver%diis%initialize()
!
      call mem%alloc(residual, wf%n_gs_amplitudes)
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
      call mem%alloc(epsilon, wf%n_gs_amplitudes)
!
      call wf%get_orbital_differences(epsilon, wf%n_gs_amplitudes)
!
!     Initialize preconditioner
!
      solver%preconditioner = precondition_tool(wf%n_gs_amplitudes)
      call solver%preconditioner%initialize_and_set_precondition_vector(epsilon)
!
      call wf%set_initial_multipliers_guess(solver%restart)
      call wf%get_multipliers(multipliers)
!
      converged_residual = .false.
!
      call output%printf('n', 'Iteration    Norm residual  ', fs='(/t3,a)')
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
         call output%printf('n', '(i3)         (e11.4)', ints=[iteration], reals=[residual_norm])
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         converged_residual = residual_norm .lt. solver%residual_threshold
!
         if (converged_residual) then
!
            call output%print_separator('n', 28,'-', fs='(t3,a)')
            call output%printf('n', 'Convergence criterion met in (i0) iterations!', &
                               ints=[iteration], fs='(t3,a)')
!
         else
!
!           Get next guess for multipliers, and perform DIIS extrapolation on it
!
            call wf%get_multipliers(multipliers)
!
            call solver%tbar_updater%update(wf, multipliers, residual)
!
            call solver%diis%update(residual, multipliers)
!
            call wf%set_multipliers(multipliers)
            call wf%save_multipliers()
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
         call output%error_msg('was not able to converge the equations      &
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
      call solver%diis%finalize()
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
      call solver%preconditioner%destruct_precondition_vector()
      call wf%save_multipliers()
!
      call solver%timer%turn_off()
!
      call output%printf('n', '- Finished solving the ' // trim(wf%name_) // &
                         ' multipliers equations', fs='(/t3, a)')
!
      call output%printf('n', 'Total wall time (sec): (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('n', 'Total cpu time (sec):  (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('cpu')], fs='(t6,a)')
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
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('m', solver%description1, ffs='(/t3,a)')
      call output%printf('m', solver%description2, ffs='(/t3,a)')
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
      call output%printf('m', '- DIIS CC multipliers solver summary:', fs='(/t3,a)')
!
      call wf%print_dominant_x_amplitudes(X, 'tbar')
!
   end subroutine print_summary_diis_cc_multipliers
!
!
   subroutine read_settings_diis_cc_multipliers(solver, records_in_memory, crop, diis_dimension)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      use global_in, only: input
!
      implicit none
!
      class(diis_cc_multipliers) :: solver
      logical, intent(inout)     :: records_in_memory, crop
      integer, intent(inout)     :: diis_dimension
!
      call input%get_keyword('threshold',                         &
                                        'solver cc multipliers',  &
                                        solver%residual_threshold)
!
      call input%get_keyword('max iterations',                    &
                                        'solver cc multipliers',  &
                                        solver%max_iterations)
!
      call input%get_keyword('diis dimension',                    &
                                        'solver cc multipliers',  &
                                        diis_dimension)
!
      crop = input%is_keyword_present('crop', 'solver cc multipliers')
!
!     Determine whether to store records in memory or on file
!
      call input%place_records_in_memory('solver cc multipliers', records_in_memory)
!
   end subroutine read_settings_diis_cc_multipliers
!
!
end module diis_cc_multipliers_class
!
