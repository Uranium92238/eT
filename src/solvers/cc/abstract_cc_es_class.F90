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
module abstract_cc_es_class
!
!!
!! Abstract coupled cluster excited state solver class module
!! Written by Eirik F. Kjønstad, Sarai D. Folkestad, June 2019
!!
!! Class that gathers functionality common the CC excited state
!! solvers (in particular, DIIS and Davidson). These solvers determine
!! states L and R and excitation energies omega that satisfy the
!! eigenvalue equation
!!
!!    A R = omega R        L^T A = omega L^T,
!!
!! where A is the coupled cluster Jacobian matrix,
!!
!!    A_mu,nu = < mu | [H-bar, tau_nu] | HF >,    H-bar = e-T H eT.
!!
!
   use parameters
!
   use global_out, only : output
!
   use memory_manager_class, only : mem
   use timings_class, only : timings
!
   use ccs_class, only : ccs
!
   use es_start_vector_tool_class, only : es_start_vector_tool
!
   use precondition_tool_class, only: precondition_tool
!
   use es_projection_tool_class, only : es_projection_tool
   use es_valence_projection_tool_class, only: es_valence_projection_tool
   use es_cvs_projection_tool_class, only: es_cvs_projection_tool
   use es_ip_projection_tool_class, only: es_ip_projection_tool
   use es_rm_core_projection_tool_class, only: es_rm_core_projection_tool
!
   use convergence_tool_class, only: convergence_tool
!
   use abstract_solver_class, only: abstract_solver
!
   implicit none
!
   type, extends(abstract_solver), abstract :: abstract_cc_es
!
      character(len=100) :: name_
      character(len=100) :: tag
!
      character(len=500) :: description1
      character(len=500) :: description2
!
      integer :: max_iterations
!
      class(convergence_tool), allocatable :: convergence_checker
!
      logical :: restart
!
      integer :: n_singlet_states
!
      character(len=40) :: transformation
      character(len=40) :: restart_transformation
      character(len=40) :: es_type
!
      real(dp), dimension(:), allocatable :: energies
!
      type(timings) :: timer
!
      class(es_start_vector_tool), allocatable  :: start_vectors
      class(es_projection_tool), allocatable    :: projector
      class(precondition_tool), allocatable     :: preconditioner
!
      class(ccs), pointer :: wf
!
   contains
!
      procedure :: print_banner                     => print_banner_abstract_cc_es
!
      procedure :: read_es_settings                 => read_es_settings_abstract_cc_es
      procedure :: print_es_settings                => print_es_settings_abstract_cc_es
!
      procedure :: cleanup                          => cleanup_abstract_cc_es
      procedure :: print_summary                    => print_summary_abstract_cc_es
!
      procedure :: prepare_wf_for_excited_state     => prepare_wf_for_excited_state_abstract_cc_es
!
      procedure :: initialize_start_vector_tool     => initialize_start_vector_tool_abstract_cc_es
      procedure :: initialize_projection_tool       => initialize_projection_tool_abstract_cc_es
!
      procedure :: initialize_energies              => initialize_energies_abstract_cc_es
      procedure :: destruct_energies                => destruct_energies_abstract_cc_es
!
      procedure :: set_initial_guesses              => set_initial_guesses_abstract_cc_es
!
   end type abstract_cc_es
!
!
contains
!
!
   subroutine initialize_energies_abstract_cc_es(this)
!!
!!    Initialize energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initialize excitation energies
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      if (.not. allocated(this%energies)) &
            call mem%alloc(this%energies, this%n_singlet_states)
!
   end subroutine initialize_energies_abstract_cc_es
!
!
   subroutine destruct_energies_abstract_cc_es(this)
!!
!!    Destruct energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Destruct excitation energies
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      if (allocated(this%energies)) &
            call mem%dealloc(this%energies, this%n_singlet_states)
!
   end subroutine destruct_energies_abstract_cc_es
!
!
   subroutine print_banner_abstract_cc_es(this)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      call output%printf('m', ' - ' // trim(this%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(this%name_)) + 6, '-')
!
      call output%printf('n', this%description1, ffs='(/t3,a)', fs='(t3,a)')
      call output%printf('n', this%description2, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_abstract_cc_es
!
!
   subroutine read_es_settings_abstract_cc_es(this, records_in_memory)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
      use global_in, only: input
!
      implicit none
!
      class(abstract_cc_es) :: this
!
      logical, intent(inout) :: records_in_memory
!
      real(dp) :: eigenvalue_threshold, residual_threshold
!
      if (input%is_keyword_present('energy threshold', 'solver cc es')) then
!
         call input%get_keyword('energy threshold', 'solver cc es', eigenvalue_threshold)
         call this%convergence_checker%set_energy_threshold(eigenvalue_threshold)
!
      endif
!
      if (input%is_keyword_present('residual threshold', 'solver cc es')) then
!
         call input%get_keyword('residual threshold', 'solver cc es', residual_threshold)
         call this%convergence_checker%set_residual_threshold(residual_threshold)
!
      endif
!
      call input%get_keyword('max iterations', 'solver cc es', this%max_iterations)
!
      call input%get_required_keyword('singlet states', 'solver cc es', this%n_singlet_states)
!
      if (input%is_keyword_present('core excitation', 'solver cc es') .and. .not. &
          input%is_keyword_present('ionization', 'solver cc es')) this%es_type = 'core'
!
      if (input%is_keyword_present('ionization', 'solver cc es') .and. .not. &
          input%is_keyword_present('core excitation', 'solver cc es')) this%es_type = 'ionize'
!
      if (input%is_keyword_present('remove core', 'solver cc es')) this%es_type = 'remove core'
!
      call input%place_records_in_memory('solver cc es', records_in_memory)
!
   end subroutine read_es_settings_abstract_cc_es
!
!
   subroutine print_es_settings_abstract_cc_es(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      call output%printf('m', '- Settings for coupled cluster excited state &
                         &solver (' //trim(this%tag) // '):', fs='(/t3,a)')
!
      call output%printf('m', 'Calculation type:    (a0)', &
                         chars=[trim(this%es_type)], fs='(/t6,a)')
      call output%printf('m', 'Excitation vectors:  (a0)', &
                         chars=[trim(this%transformation)], fs='(t6,a)')
!
      call this%convergence_checker%print_settings()
!
      call output%printf('m', 'Number of singlet states:     (i11)', &
                         ints=[this%n_singlet_states], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations:     (i11)', &
                         ints=[this%max_iterations], fs='(t6,a)')
!
   end subroutine print_es_settings_abstract_cc_es
!
!
   subroutine cleanup_abstract_cc_es(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(abstract_cc_es), intent(inout) :: this
!
      call this%destruct_energies()
      call this%projector%destruct_projection_vector()
!
      call this%timer%turn_off()
!
     call output%printf('m', '- Finished solving the ' //  &
                        trim(convert_to_uppercase(this%wf%name_)) // ' excited state &
                        &equations ('// trim(this%transformation) //')', fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec): (f20.5)', &
                         reals=[this%timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'Total cpu time (sec):  (f20.5)', &
                         reals=[this%timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_abstract_cc_es
!
!
   subroutine print_summary_abstract_cc_es(this, X, X_order)
!!
!!    Print summary
!!    Written by Eirik F. Kjønstad, Dec 2018
!!    Modified by Eirik F. Kjønstad, Mar 2020
!!
!!    Prints summary of excited states. Lists the dominant amplitudes and
!!    the energies, along with fraction of singles.
!!
!!    X:         array with excited states stored in the columns
!!
!!    X_order:   optional index list giving the ordering of states from low
!!               to high energy. Default is to assume that X is already
!!               ordered from low to high energies.
!!
!!    Warning: it is assumed on entry that the energies are already ordered according
!!             to energy. It is just the states that are not necessarily ordered.
!!
!!    Eirik F. Kjønstad, Mar 2020: added X_order to avoid duplicate copy of states
!!                                 in DIIS solver.
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(abstract_cc_es), intent(in) :: this
!
      real(dp), dimension(this%wf%n_es_amplitudes, this%n_singlet_states), intent(in) :: X
!
      integer, dimension(this%n_singlet_states), optional, intent(in) :: X_order
!
      integer, dimension(:), allocatable :: X_order_local
!
      integer :: state, state_index
!
      character(len=1) :: label ! R or L, depending on whether left or right transformation
!
!     Set up list that gives ordering of energies from low to high
!
      call mem%alloc(X_order_local, this%n_singlet_states)
!
      if (present(X_order)) then
!
         X_order_local = X_order
!
      else
!
         do state = 1, this%n_singlet_states
!
            X_order_local(state) = state
!
         enddo
!
      endif
!
!     Print excited state vectors
!
      label = trim(adjustl(convert_to_uppercase(this%transformation(1:1))))
!
      call output%printf('n', '- Excitation vector amplitudes:', fs='(/t3,a)')
!
      do state = 1, this%n_singlet_states
!
         state_index = X_order_local(state)
!
         call output%printf('n', 'Electronic state nr. (i0)', ints=[state], fs='(/t6,a)')
!
         call output%printf('n', 'Energy (Hartree):             (f19.12)', &
                            reals=[this%energies(state_index)], fs='(/t6,a)')
!
         call this%wf%print_X1_diagnostics(X(:,state_index), label)
         call this%wf%print_dominant_x_amplitudes(X(1, state_index), label)
!
      enddo
!
      call mem%dealloc(X_order_local, this%n_singlet_states)
!
!     Print excited state energies
!
      call output%printf('m', '- Electronic excitation energies:', fs='(/t6,a)')
!
      call output%printf('m', 'Excitation energy', fs='(/t39,a)')
      call output%print_separator('m', 42, '-', fs='(t27,a)')
      call output%printf('m', ' State                (Hartree)             (eV)', fs='(t6,a)')
      call output%print_separator('m', 63, '-', fs='(t6,a)')
!
      do state = 1, this%n_singlet_states
!
         call output%printf('m', '(i4)             (f19.12)   (f19.12)',         &
                            ints=[state], reals=[this%energies(state),   &
                            this%energies(state)*Hartree_to_eV], fs='(t6,a)')
!
      enddo
!
      call output%print_separator('m', 63, '-', fs='(t6,a)')
      call output%printf('m', 'eV/Hartree (CODATA 2014): (f11.8)', &
                         reals=[Hartree_to_eV], fs='(t6,a)')
!
   end subroutine print_summary_abstract_cc_es
!
!
   subroutine prepare_wf_for_excited_state_abstract_cc_es(this)
!!
!!    Prepare wf for excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2019
!!
      implicit none
!
      class(abstract_cc_es), intent(in) :: this
!
      call this%wf%initialize_excited_state_files()
      call this%wf%prepare_for_Jacobians(this%transformation)
!
      if (this%transformation == 'right') then
!
         call this%wf%initialize_right_excitation_energies()
!
      else if (this%transformation == 'left') then
!
         call this%wf%initialize_left_excitation_energies()
!
      else if (this%transformation == 'both') then
!
         call this%wf%initialize_right_excitation_energies()
         call this%wf%initialize_left_excitation_energies()
!
      end if
!
   end subroutine prepare_wf_for_excited_state_abstract_cc_es
!
!
   subroutine initialize_start_vector_tool_abstract_cc_es(this)
!!
!!    Initialize start vector tool
!!    Written by Eirik F. Kjønstad, Sep 2019
!!
      use global_in, only: input
!
      use es_manual_start_vector_tool_class,    only: es_manual_start_vector_tool
      use es_valence_start_vector_tool_class,   only: es_valence_start_vector_tool
      use es_cvs_start_vector_tool_class,       only: es_cvs_start_vector_tool
      use es_ip_start_vector_tool_class,        only: es_ip_start_vector_tool
!
      implicit none
!
      class(abstract_cc_es) :: this
!
      if (trim(this%es_type) == 'core') then
!
         call this%wf%read_cvs_settings()

      elseif (trim(this%es_type) == 'remove core') then
!
         call this%wf%read_rm_core_settings()
!
      endif
!
      if (input%is_keyword_present('state guesses', 'solver cc es')) then
!
         this%start_vectors = es_manual_start_vector_tool(this%wf)
!
      else
!
         if (trim(this%es_type) == 'valence') then
!
            this%start_vectors = es_valence_start_vector_tool(this%wf)
!
         elseif (trim(this%es_type) == 'core') then
!
            this%start_vectors = es_cvs_start_vector_tool(this%wf)
!
         elseif (trim(this%es_type) == 'ionize') then
!
            this%start_vectors = es_ip_start_vector_tool(this%wf)
!
         elseif (trim(this%es_type) == 'remove core') then
!
            this%start_vectors = es_valence_start_vector_tool(this%wf)
!
         else
!
            call output%error_msg('could not recognize excited state type in abstract_cc_es')
!
         endif
!
      endif
!
   end subroutine initialize_start_vector_tool_abstract_cc_es
!
!
   subroutine initialize_projection_tool_abstract_cc_es(this)
!!
!!    Initialize projection tool
!!    Written by Eirik F. Kjønstad, Sep 2019
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      if (trim(this%es_type) == 'valence') then
!
         this%projector = es_valence_projection_tool()
!
      elseif (trim(this%es_type) == 'core') then
!
         this%projector = es_cvs_projection_tool(this%wf)
!
      elseif (trim(this%es_type) == 'ionize') then
!
         this%projector = es_ip_projection_tool(this%wf)
!
      elseif (trim(this%es_type) == 'remove core') then
!
         this%projector = es_rm_core_projection_tool(this%wf)
!
      else
!
         call output%error_msg('could not recognize excited state type in abstract_cc_es')
!
      endif
!
   end subroutine initialize_projection_tool_abstract_cc_es
!
!
   subroutine set_initial_guesses_abstract_cc_es(this, X, first, last)
!!
!!    Set initial guesses
!!    Written by Alexander C. Paul, Oct 2020
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      integer, intent(in) :: first, last
!
      real(dp), dimension(this%wf%n_es_amplitudes,first:last), intent(out) :: X
!
      integer :: state
!
      do state = first, last
!
         call this%start_vectors%get(this%wf, state, X(:,state),  &
                                       this%energies(state), &
                                       this%transformation,  &
                                       this%restart)
!
         if (this%projector%active) call this%projector%do_(X(:,state))
!
      enddo
!
   end subroutine set_initial_guesses_abstract_cc_es
!
!
end module abstract_cc_es_class
