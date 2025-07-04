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
module cc_lr_transition_moments_task_class
!
!!
!! CC linear response transition moments task class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use parameters
   use ccs_class, only: ccs
   use cc_transition_moments_task_class, only: cc_transition_moments_task
   use memory_manager_class, only: mem
   use global_out, only: output
   use stream_file_class, only: stream_file
   use cc_eta_xi_calculator_class, only: cc_eta_xi_calculator
   use cc_lr_eta_xi_calculator_class, only: cc_lr_eta_xi_calculator
!
   implicit none
!
   type, extends(cc_transition_moments_task) :: cc_lr_transition_moments_task
!
      logical, private :: dipole_length
!
      real(dp), dimension(:,:), allocatable, private :: xiX
      real(dp), dimension(:,:), allocatable, private :: etaX
!
      type(stream_file), dimension(:), allocatable :: M_vectors
!
      class(cc_eta_xi_calculator), allocatable, private :: eta_xi_calculator
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_lr_transition_moments_task
!
      procedure, private         :: determine_M_vectors
      procedure, private, nopass :: solve_M_vector_equations
!
   end type cc_lr_transition_moments_task
!
   interface cc_lr_transition_moments_task
!
      procedure :: new_cc_lr_transition_moments_task
!
   end interface cc_lr_transition_moments_task
!
!
contains
!
!
   function new_cc_lr_transition_moments_task() result(this)
!!
!!    New CC LR transition moments task
!!    Written by Alexander C. Paul, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      type(cc_lr_transition_moments_task) :: this
!
      this%name_ = 'Determining CC LR transition moments'
!
      this%dipole_length = input%is_keyword_present('dipole length','cc response')
!
      if (.not. this%dipole_length) &
         call output%error_msg('No operator selected in response calculation')
!
   end function new_cc_lr_transition_moments_task
!
!
   subroutine execute_cc_lr_transition_moments_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, Nov 2019
!!
!!    Adapted from routines written by Josefine H. Andersen, spring 2019.
!!
!!    Routine that performs three tasks:
!!
!!       - Calls routine to determine the M vectors, which are needed to compute the
!!         transition moments in response theory
!!
!!       - Calls wavefunction routine that computes the moments and strengths from
!!         the etaX, xiX and M vectors.
!!
!!       - Calls routine to print the transition moments / strengths to output
!!
!!    This routine was based on a deprecated routine for EOM moments
!!    (Josefine H. Andersen, Sarai D. Folkestad, June 2019). Reintroduced
!!    with M vector contributions and making use of skip_states array
!!    (Eirik F. Kjønstad, Nov 2019).
!!    Removed skip_states array as the parallel states are removed by the engine
!!    (Alexander C. Paul Sep 2020)
!!
      use timings_class, only: timings
!
      implicit none
!
      class(cc_lr_transition_moments_task), intent(inout) :: this
      class(ccs), intent(inout), target :: wf
!
      real(dp), dimension(:,:,:), allocatable :: transition_moments
!
      integer :: state, component
!
      real(dp), dimension(:), allocatable :: M ! M vector associated with state
!
      call this%print_header()
      call this%start_timer()
!
      call wf%prepare_for_properties
!
!     Dipole matrix: only first row and column will be computed for TM from the GS
      call mem%alloc(transition_moments, wf%n_singlet_states+1, wf%n_singlet_states+1, 3)
!
      call mem%alloc(this%xiX, wf%n_es_amplitudes, 3)
      call mem%alloc(this%etaX, wf%n_es_amplitudes, 3)
!
      this%eta_xi_calculator = cc_lr_eta_xi_calculator(wf)
      call this%eta_xi_calculator%calculate(this%xiX, this%etaX)
!
      call this%determine_M_vectors(wf)
!
      call mem%alloc(M, wf%n_es_amplitudes)
!
      do state = 1, wf%n_singlet_states
!
         call this%M_vectors(state)%open_('rewind')
         call this%M_vectors(state)%read_(M, wf%n_es_amplitudes)
         call this%M_vectors(state)%close_()
!
         do component = 1, 3
!
            call wf%calculate_lr_transition_strength(this%etaX(:,component), &
                                                     this%xiX(:,component), &
                                                     state, &
                                                     transition_moments(state+1, 1, component), &
                                                     transition_moments(1, state+1, component), &
                                                     M)
!
         enddo
!
      enddo
!
      call mem%dealloc(M, wf%n_es_amplitudes)
!
      call mem%dealloc(this%xiX, wf%n_es_amplitudes, 3)
      call mem%dealloc(this%etaX, wf%n_es_amplitudes, 3)
!
      call this%print_transition_moment_summary(wf, transition_moments, &
                                                initial_states=[0], &
                                                calculation_type='LR')
!
      call mem%dealloc(transition_moments, wf%n_singlet_states+1, wf%n_singlet_states+1, 3)
!
      call this%end_timer()
!
   end subroutine execute_cc_lr_transition_moments_task
!
!
   subroutine determine_M_vectors(this, wf)
!!
!!    Determine M vectors
!!    Written by Eirik F. Kjønstad and Josefine H. Andersen, 2019
!!
!!    Determines the solutions to the equations
!!
!!       (A^T + omega_k I) M_k = -F R_k,
!!
!!    where F is the so-called F transformation,
!!    (F_mu,nu = < Lambda | [[H-bar,tau_mu],tau_nu] | HF > ), and
!!    the omega_k are the excitation energies.
!!
!!    The M vectors are required to compute the transition moments
!!    using CC response theory.
!!
!!    On exit, the converged M vectors are stored in the file array
!!    this%M_vectors.
!!
!!    This solver wrapper routine, and the file storage for M vectors, was made by
!!    Eirik F. Kjønstad, Nov 2019. It is adapted/based on the general structure set up
!!    to determine M vectors originally written by Josefine H. Andersen, spring 2019.
!!
      implicit none
!
      class(cc_lr_transition_moments_task), intent(inout) :: this
!
      class(ccs) :: wf
!
      integer :: k
!
      real(dp), dimension(:), allocatable :: R           ! Stores R_k temporarily
      real(dp), dimension(:,:), allocatable :: minus_FR  ! [-F R_k], k = 1, 2, ..., n_singlet_states
!
      character(len=200) :: file_name
!
!     Build the matrix FR of right-hand-sides, -FR = [-F R_k], k = 1, 2, ..., n_singlet_states
!
      call mem%alloc(R, wf%n_es_amplitudes)
      call mem%alloc(minus_FR, wf%n_es_amplitudes, wf%n_singlet_states)
!
      do k = 1, wf%n_singlet_states
!
         call wf%read_excited_state(R, k, k, 'right')
!
         call wf%F_transformation(R, minus_FR(:, k))
!
      enddo
!
      call mem%dealloc(R, wf%n_es_amplitudes)
!
      call dscal(wf%n_es_amplitudes*wf%n_singlet_states, -one, minus_FR, 1)
!
!     Initialize file array to store the converged M_vectors
!
      allocate(this%M_vectors(wf%n_singlet_states))
!
      do k = 1, wf%n_singlet_states
!
         write(file_name, '(a, i3.3)') 'M_vector_state_', k
!
         this%M_vectors(k) = stream_file(trim(file_name))
!
      enddo
!
!     Converge and store M vectors
!
      call this%solve_M_vector_equations(wf, minus_FR, this%M_vectors)
!
      call mem%dealloc(minus_FR, wf%n_es_amplitudes, wf%n_singlet_states)
!
   end subroutine determine_M_vectors
!
!
   subroutine solve_M_vector_equations(wf, minus_FR, M_files)
!!
!!    Solve M vector equations
!!    Written by Eirik Kjønstad, Oct 2022
!!
      use linear_davidson_tool_class,                          only: linear_davidson_tool
      use file_linear_equation_storage_tool_class,             only: file_linear_equation_storage_tool
      use vector_getter_rhs_tool_class,                        only: vector_getter_rhs_tool
      use linear_davidson_multiple_solutions_print_tool_class, only: linear_davidson_multiple_solutions_print_tool
!
      use cc_jacobian_transformation_class,        only: cc_jacobian_transformation
      use linear_davidson_solver_class,            only: linear_davidson_solver
      use cc_jacobian_preconditioner_getter_class, only: cc_jacobian_preconditioner_getter
      use cc_response_solver_settings_class,       only: cc_response_solver_settings
!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes, wf%n_singlet_states), intent(in), target :: minus_FR
!
      class(stream_file), dimension(wf%n_singlet_states), intent(inout), target :: M_files
!
      class(linear_davidson_solver),                        allocatable :: solver
      class(linear_davidson_tool),                          allocatable :: davidson
      class(cc_jacobian_transformation),                    allocatable :: transformer
      class(file_linear_equation_storage_tool),             allocatable :: storer
      class(vector_getter_rhs_tool),                        allocatable :: rhs_getter
      class(cc_jacobian_preconditioner_getter),             allocatable :: preconditioner
      class(linear_davidson_multiple_solutions_print_tool), allocatable :: printer
      class(cc_response_solver_settings),                   allocatable :: settings
!
      real(dp), dimension(:), allocatable :: frequencies
!
      settings = cc_response_solver_settings()
!
      printer = linear_davidson_multiple_solutions_print_tool()
!
      davidson = linear_davidson_tool(name_             = 'M_vectors_davidson',              &
                                      n_parameters      = wf%n_es_amplitudes,                &
                                      lindep_threshold  = min(1.0d-11, settings%threshold),  &
                                      max_dim_red       = 100,                               &
                                      n_rhs             = wf%n_singlet_states,               &
                                      n_frequencies     = wf%n_singlet_states,               &
                                      records_in_memory = settings%records_in_memory)
!
      transformer = cc_jacobian_transformation(wf,                         &
                                               side         = 'left',      &
                                               n_parameters = wf%n_es_amplitudes)
!
      storer = file_linear_equation_storage_tool(n_parameters = wf%n_es_amplitudes,    &
                                                 n_solutions  = wf%n_singlet_states,   &
                                                 files        = M_files)
!
      preconditioner = cc_jacobian_preconditioner_getter(wf, wf%n_es_amplitudes)
!
      rhs_getter = vector_getter_rhs_tool(n_parameters = wf%n_es_amplitudes,  &
                                          n_rhs        = wf%n_singlet_states, &
                                          rhs          = minus_FR)
!
      frequencies = -wf%right_excitation_energies
      solver = linear_davidson_solver(transformer        = transformer,             &
                                      davidson           = davidson,                &
                                      storer             = storer,                  &
                                      preconditioner     = preconditioner,          &
                                      rhs_getter         = rhs_getter,              &
                                      printer            = printer,                 &
                                      n_rhs              = wf%n_singlet_states,     &
                                      n_solutions        = wf%n_singlet_states,     &
                                      max_iterations     = settings%max_iterations, &
                                      residual_threshold = settings%threshold,      &
                                      frequencies        = frequencies)
!
      call solver%run()
      call solver%cleanup()
!
   end subroutine solve_M_vector_equations
!
!
end module cc_lr_transition_moments_task_class
