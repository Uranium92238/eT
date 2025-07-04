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
module diis_cc_gs_class
!
!!
!!	DIIS coupled cluster ground state solver class module
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
!!    |CC > = e^T | HF >.
!!
!! The equation is solved using the direct inversion of the iterative
!! subspace (DIIS) algorithm. See Pulay, P. Convergence acceleration
!! of iterative sequences. The case of SCF iteration. Chem. Phys. Lett.
!! 1980, 73, 393−398. This algorithm combines estimates for the parameters
!! and the errors and finds a least-squares solution to the error being
!! zero (see diis_tool solver tool for more details).
!!
!
   use parameters
!
   use global_out,              only: output
   use memory_manager_class,    only: mem
   use ccs_class,               only: ccs
   use diis_tool_class,         only: diis_tool
   use timings_class,           only: timings
   use convergence_tool_class,  only: convergence_tool
!
   use amplitude_updater_class, only: amplitude_updater
!
   use abstract_solver_class, only: abstract_solver
!
   implicit none
!
   type, extends(abstract_solver) :: diis_cc_gs
!
      character(len=100) :: name_ = 'DIIS coupled cluster ground state solver'
!
      character(len=500) :: description1 = 'A DIIS CC ground state amplitude equations solver. It uses &
                                           &an extrapolation of previous quasi-Newton perturbation theory &
                                           &estimates of the next amplitudes. See Helgaker et al., Molecular &
                                           &Electronic Structure Theory, Chapter 13.'
!
      integer :: max_iterations
!
      logical :: restart
!
      class(amplitude_updater), allocatable :: t_updater
!
      class(convergence_tool), allocatable :: convergence_checker
!
      type(diis_tool), allocatable :: diis
!
      class(ccs), pointer :: wf
!
   contains
!
      procedure :: run                      => run_diis_cc_gs
      procedure :: cleanup                  => cleanup_diis_cc_gs
!
      procedure :: print_banner             => print_banner_diis_cc_gs
      procedure :: read_settings            => read_settings_diis_cc_gs
!
      procedure :: print_settings           => print_settings_diis_cc_gs
!
   end type diis_cc_gs
!
!
   interface diis_cc_gs
!
      procedure :: new_diis_cc_gs
!
   end interface diis_cc_gs
!
!
contains
!
!
   function new_diis_cc_gs(wf, restart, t_updater, storage) result(this)
!!
!!    New DIIS CC GS
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    restart:   read amplitudes instead of generating an initial guess
!!    t_updater: object for calculating the t-update (i.e. dt_mu), see 'amplitude_updater' class.
!!
      implicit none
!
      type(diis_cc_gs) :: this
!
      class(ccs), target :: wf
!
      logical, intent(in) :: restart
!
      class(amplitude_updater), intent(in) :: t_updater
!
      character(len=*), intent(in) :: storage
!
      integer :: diis_dimension
      logical :: crop ! Standard DIIS if false; CROP variant of DIIS if true
      logical :: records_in_memory
!
      this%wf => wf
!
      this%iteration_timer = timings('DIIS CC GS iteration', pl='normal')
      this%total_timer = timings('CC GS solver total', pl='minimal')
      call this%total_timer%turn_on()
!
!     Print solver banner
!
      call this%print_banner()
!
!     Set standard settings
!
      this%max_iterations   = 100
      this%restart          = restart
      this%t_updater        = t_updater
!
      diis_dimension    = 8
      crop              = .false.
      records_in_memory = .false.
!
!     Initialize convergence checker with default threshols
!
      this%convergence_checker = convergence_tool(energy_threshold   = 1.0d-5,   &
                                                    residual_threshold = 1.0d-5,   &
                                                    energy_convergence = .false.)
!
!     Read & print settings (thresholds, etc.)
!
      call this%read_settings(crop, diis_dimension)
      call this%print_settings()
!
!     Set the amplitudes to the initial guess or read if restart
!
      call this%wf%initialize_amplitudes()
!
      call this%wf%set_initial_amplitudes_guess(this%restart)
!
!     Determine whether to store records in memory or on file
!
      if (trim(storage) == 'memory') then
!
         records_in_memory = .true.
!
      elseif (trim(storage) == 'disk') then
!
         records_in_memory = .false.
!
      else
!
         call output%error_msg('Could not recognize keyword storage in solver: ' // trim(storage))
!
      endif
!
      this%diis = diis_tool('cc_gs_diis',            &
                        this%wf%n_gs_amplitudes,     &
                        this%wf%n_gs_amplitudes,     &
                        dimension_=diis_dimension,   &
                        crop=crop,                   &
                        records_in_memory=records_in_memory)
!
   end function new_diis_cc_gs
!
!
   subroutine print_settings_diis_cc_gs(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(diis_cc_gs) :: this
!
      call output%printf('m', '- CC ground state solver settings:', fs='(/t3,a)')
!
      call this%convergence_checker%print_settings()
!
      call output%printf('m', 'Max number of iterations: (i8)', &
                         ints=[this%max_iterations], fs='(t6, a)')
!
   end subroutine print_settings_diis_cc_gs
!
!
   subroutine run_diis_cc_gs(this)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use array_utilities, only: get_l2_norm
!
      implicit none
!
      class(diis_cc_gs), intent(inout) :: this
!
      logical :: converged
!
      real(dp) :: energy, previous_energy
      real(dp) :: omega_norm
!
      real(dp), dimension(:), allocatable :: omega
      real(dp), dimension(:), allocatable :: amplitudes
!
      integer :: iteration
!
      call this%diis%initialize()
!
      call mem%alloc(omega, this%wf%n_gs_amplitudes)
      call mem%alloc(amplitudes, this%wf%n_gs_amplitudes)
!
      converged = .false.
!
      call output%printf('n', 'Iteration    Energy (a.u.)        |omega|       &
                         &Delta E (a.u.) ', fs='(/t3,a)')
      call output%print_separator('n', 63,'-')
!
      iteration = 0
      previous_energy = zero
!
      call this%wf%save_amplitudes()
!
      do while (.not. converged .and. iteration .le. this%max_iterations)
!
         iteration = iteration + 1
         call this%iteration_timer%turn_on()
!
!        Calculate the energy and error vector omega
!
         call this%wf%construct_fock(task = 'gs')
!
         call this%wf%calculate_energy()
         energy = this%wf%energy
!
         call this%wf%construct_omega(omega)
!
         omega_norm = get_l2_norm(omega, this%wf%n_gs_amplitudes)
!
         call output%printf('n', '(i3)  (f25.12)    (e11.4)    (e11.4)', &
                            ints=[iteration], reals=[this%wf%energy, omega_norm, &
                            abs(this%wf%energy-previous_energy)], fs='(t3, a)')
!
!        Test for convergence & prepare for next iteration if not yet converged
!
         converged = this%convergence_checker%has_converged(omega_norm,                    &
                                                              this%wf%energy - previous_energy,   &
                                                              iteration)
!
         if (converged) then
!
            call output%print_separator('n', 63,'-')
!
            call output%printf('n', 'Convergence criterion met in (i0) iterations!', &
                               ints=[iteration], fs='(t3,a)')
!
         else
!
!           Calculate the next guess for the amplitudes,
!           and improve it with DIIS extrapolation
!
            call this%wf%get_amplitudes(amplitudes)
            call this%t_updater%update(this%wf, amplitudes, omega)
!
            call this%diis%update(omega, amplitudes)
!
            call this%wf%set_amplitudes(amplitudes)
!
!           Update t1-transformed electron repulsion integrals
!
            call this%wf%construct_t1_cholesky(this%wf%t1, this%wf%L_mo, this%wf%L_t1)
!
         endif
!
         previous_energy = energy
         call this%wf%save_amplitudes()
!
         call this%iteration_timer%turn_off()
         call this%iteration_timer%reset()
!
      enddo
!
      call mem%dealloc(omega, this%wf%n_gs_amplitudes)
      call mem%dealloc(amplitudes, this%wf%n_gs_amplitudes)
!
      if (.not. converged) &
         call output%error_msg('Did not converge in the max number of iterations.')
!
      call this%wf%print_gs_summary()
!
      call this%diis%finalize()
!
   end subroutine run_diis_cc_gs
!
!
   subroutine cleanup_diis_cc_gs(this)
!!
!! 	Cleanup
!! 	Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(diis_cc_gs), intent(inout) :: this
!
      call this%wf%save_amplitudes()
!
      call this%total_timer%turn_off()
!
      call output%printf('m', '- Finished solving the (a0) ground state equations', &
                         chars=[convert_to_uppercase(this%wf%name_)], fs='(/t3, a)')
!
      call output%printf('m', 'Total wall time (sec): (f20.5)', &
                         reals=[this%total_timer%get_elapsed_time('wall')], fs='(/t6, a)')
      call output%printf('m', 'Total cpu time (sec):  (f20.5)', &
                         reals=[this%total_timer%get_elapsed_time('cpu')], fs='(t6, a)')
!
   end subroutine cleanup_diis_cc_gs
!
!
   subroutine print_banner_diis_cc_gs(this)
!!
!!    Print banner
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none
!
      class(diis_cc_gs) :: this
!
      call output%printf('m', ' - ' // trim(this%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(this%name_)) + 6, '-')
!
      call output%printf('m', this%description1, ffs='(/t3,a)')
!
   end subroutine print_banner_diis_cc_gs
!
!
   subroutine read_settings_diis_cc_gs(this, crop, diis_dimension)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      use global_in, only: input
!
      implicit none
!
      class(diis_cc_gs)      :: this
!
      logical, intent(inout) :: crop
      integer, intent(inout) :: diis_dimension
!
      real(dp) :: energy_threshold, omega_threshold
!
      if (input%is_keyword_present('energy threshold', 'solver cc gs')) then
!
         call input%get_keyword('energy threshold', 'solver cc gs', energy_threshold)
         call this%convergence_checker%set_energy_threshold(energy_threshold)
!
      endif
!
      if (input%is_keyword_present('omega threshold', 'solver cc gs')) then
!
         call input%get_keyword('omega threshold', 'solver cc gs', omega_threshold)
         call this%convergence_checker%set_residual_threshold(omega_threshold)
!
      endif
!
      call input%get_keyword('max iterations', 'solver cc gs', this%max_iterations)
!
      crop = input%is_keyword_present('crop', 'solver cc gs')
      call input%get_keyword('diis dimension', 'solver cc gs', diis_dimension)
!
   end subroutine read_settings_diis_cc_gs
!
!
end module diis_cc_gs_class
