!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
   use global_in,  only :input
!
   use memory_manager_class, only : mem
   use timings_class, only : timings
   use array_utilities, only : get_l2_norm
   use string_utilities, only : convert_to_uppercase
!
   use ccs_class, only : ccs
!
   use es_start_vector_tool_class, only : es_start_vector_tool
   use es_valence_start_vector_tool_class, only: es_valence_start_vector_tool
   use es_cvs_start_vector_tool_class, only: es_cvs_start_vector_tool
   use es_ip_start_vector_tool_class, only: es_ip_start_vector_tool
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
   implicit none
!
   type, abstract :: abstract_cc_es
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
   contains
!
      procedure(run_abstract_cc_es), deferred :: run
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
   abstract interface
      subroutine run_abstract_cc_es(solver, wf)
!
         import abstract_cc_es, ccs
!
         implicit none
!
         class(abstract_cc_es) :: solver
!
         class(ccs) :: wf
!
      end subroutine
   end interface
!
!
contains
!
!
   subroutine initialize_energies_abstract_cc_es(solver)
!!
!!    Initialize energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initialize excitation energies
!!
      implicit none
!
      class(abstract_cc_es) :: solver
!
      if (.not. allocated(solver%energies)) &
            call mem%alloc(solver%energies, solver%n_singlet_states)
!
   end subroutine initialize_energies_abstract_cc_es
!
!
   subroutine destruct_energies_abstract_cc_es(solver)
!!
!!    Destruct energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Destruct excitation energies
!!
      implicit none
!
      class(abstract_cc_es) :: solver
!
      if (allocated(solver%energies)) &
            call mem%dealloc(solver%energies, solver%n_singlet_states)
!
   end subroutine destruct_energies_abstract_cc_es
!
!
   subroutine print_banner_abstract_cc_es(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(abstract_cc_es) :: solver 
!
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('n', solver%description1, ffs='(/t3,a)', fs='(t3,a)')
      call output%printf('n', solver%description2, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_abstract_cc_es
!
!
   subroutine read_es_settings_abstract_cc_es(solver, records_in_memory)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(abstract_cc_es) :: solver 
!
      logical, intent(inout) :: records_in_memory
!
      real(dp) :: eigenvalue_threshold, residual_threshold
!
      if (input%is_keyword_present('energy threshold', 'solver cc es')) then
!
         call input%get_keyword('energy threshold', 'solver cc es', eigenvalue_threshold)
         call solver%convergence_checker%set_energy_threshold(eigenvalue_threshold)
!
      endif
!
      if (input%is_keyword_present('residual threshold', 'solver cc es')) then
!
         call input%get_keyword('residual threshold', 'solver cc es', residual_threshold)
         call solver%convergence_checker%set_residual_threshold(residual_threshold)
!
      endif
!
      call input%get_keyword('max iterations', 'solver cc es', solver%max_iterations)
!               
      call input%get_required_keyword('singlet states', 'solver cc es', solver%n_singlet_states)
!  
      if (input%is_keyword_present('core excitation', 'solver cc es') .and. .not. &
          input%is_keyword_present('ionization', 'solver cc es')) solver%es_type = 'core'
!
      if (input%is_keyword_present('ionization', 'solver cc es') .and. .not. &
          input%is_keyword_present('core excitation', 'solver cc es')) solver%es_type = 'ionize'
!
      if (input%is_keyword_present('remove core', 'solver cc es')) solver%es_type = 'remove core'
!
      call input%place_records_in_memory('solver cc es', records_in_memory)
!
   end subroutine read_es_settings_abstract_cc_es
!
!
   subroutine print_es_settings_abstract_cc_es(solver)
!!
!!    Print settings    
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(abstract_cc_es) :: solver 
!
      call output%printf('m', '- Settings for coupled cluster excited state &
                         &solver (' //trim(solver%tag) // '):', fs='(/t3,a)')
!
      call output%printf('m', 'Calculation type:    (a0)', &
                         chars=[trim(solver%es_type)], fs='(/t6,a)')
      call output%printf('m', 'Excitation vectors:  (a0)', &
                         chars=[trim(solver%transformation)], fs='(t6,a)')
!
      call solver%convergence_checker%print_settings()
!
      call output%printf('m', 'Number of singlet states:     (i11)', &
                         ints=[solver%n_singlet_states], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations:     (i11)', &
                         ints=[solver%max_iterations], fs='(t6,a)')
!
   end subroutine print_es_settings_abstract_cc_es
!
!
   subroutine cleanup_abstract_cc_es(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(abstract_cc_es) :: solver
      class(ccs), intent(in) :: wf
!
      call solver%destruct_energies()
      call solver%projector%destruct_projection_vector()
!
      call solver%timer%turn_off()
!
     call output%printf('m', '- Finished solving the ' //  &
                        trim(convert_to_uppercase(wf%name_)) // ' excited state &
                        &equations ('// trim(solver%transformation) //')', fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec): (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'Total cpu time (sec):  (f20.5)', &
                         reals=[solver%timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_abstract_cc_es
!
!
   subroutine print_summary_abstract_cc_es(solver, wf, X, X_order)
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
      implicit none 
!
      class(abstract_cc_es), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes, solver%n_singlet_states), intent(in) :: X
!
      integer, dimension(solver%n_singlet_states), optional, intent(in) :: X_order
!
      integer, dimension(:), allocatable :: X_order_local
!
      integer :: state, state_index
!
      character(len=1) :: label ! R or L, depending on whether left or right transformation 
!
!     Set up list that gives ordering of energies from low to high
!
      call mem%alloc(X_order_local, solver%n_singlet_states)
!
      if (present(X_order)) then 
!
         X_order_local = X_order
!
      else 
!
         do state = 1, solver%n_singlet_states
!
            X_order_local(state) = state
!
         enddo
!
      endif 
!
!     Print excited state vectors 
!
      label = trim(adjustl(convert_to_uppercase(solver%transformation(1:1))))
!
      call output%printf('n', '- Excitation vector amplitudes:', fs='(/t3,a)')
!
      do state = 1, solver%n_singlet_states 
!
         state_index = X_order_local(state) 
!
         call output%printf('n', 'Electronic state nr. (i0)', ints=[state], fs='(/t6,a)')
!
         call output%printf('n', 'Energy (Hartree):             (f19.12)', &
                            reals=[solver%energies(state_index)], fs='(/t6,a)')
!
         call wf%print_X1_diagnostics(X(:,state_index), label)
         call wf%print_dominant_x_amplitudes(X(1, state_index), label)
!
      enddo 
!
      call mem%dealloc(X_order_local, solver%n_singlet_states)
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
      do state = 1, solver%n_singlet_states
!
         call output%printf('m', '(i4)             (f19.12)   (f19.12)',         &
                            ints=[state], reals=[solver%energies(state),   &
                            solver%energies(state)*Hartree_to_eV], fs='(t6,a)')
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
   subroutine prepare_wf_for_excited_state_abstract_cc_es(solver, wf)
!!
!!    Prepare wf for excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2019
!!
      implicit none
!
      class(abstract_cc_es), intent(in) :: solver 
      class(ccs), intent(inout)         :: wf
!
      call wf%initialize_excited_state_files()
      call wf%prepare_for_Jacobians(solver%transformation)
!
      if (solver%transformation == 'right') then 
!
         call wf%initialize_right_excitation_energies()
!
      else if (solver%transformation == 'left') then 
!
         call wf%initialize_left_excitation_energies()
!
      else if (solver%transformation == 'both') then 
!
         call wf%initialize_right_excitation_energies()
         call wf%initialize_left_excitation_energies()
!
      end if
!
   end subroutine prepare_wf_for_excited_state_abstract_cc_es
!
!
   subroutine initialize_start_vector_tool_abstract_cc_es(solver, wf) 
!!
!!    Initialize start vector tool 
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none
!
      class(abstract_cc_es) :: solver 
!
      class(ccs) :: wf 
!
      if (trim(solver%es_type) == 'valence') then 
!
         solver%start_vectors = es_valence_start_vector_tool(wf)
!
      elseif (trim(solver%es_type) == 'core') then 
!
         call wf%read_cvs_settings()
         solver%start_vectors = es_cvs_start_vector_tool(wf)
!
      elseif (trim(solver%es_type) == 'ionize') then 
!
         solver%start_vectors = es_ip_start_vector_tool(wf)
!
      elseif (trim(solver%es_type) == 'remove core') then
!
         call wf%read_rm_core_settings()
         solver%start_vectors = es_valence_start_vector_tool(wf)
!
      else 
!
         call output%error_msg('could not recognize excited state type in abstract_cc_es')
!
      endif
!
   end subroutine initialize_start_vector_tool_abstract_cc_es
!
!
   subroutine initialize_projection_tool_abstract_cc_es(solver, wf)
!!
!!    Initialize projection tool 
!!    Written by Eirik F. Kjønstad, Sep 2019 
!!
      implicit none
!
      class(abstract_cc_es) :: solver 
!
      class(ccs) :: wf 
!
      if (trim(solver%es_type) == 'valence') then 
!
         solver%projector = es_valence_projection_tool()
!
      elseif (trim(solver%es_type) == 'core') then 
!
         solver%projector = es_cvs_projection_tool(wf)
!
      elseif (trim(solver%es_type) == 'ionize') then 
!
         solver%projector = es_ip_projection_tool(wf)
!
      elseif (trim(solver%es_type) == 'remove core') then 
!
         solver%projector = es_rm_core_projection_tool(wf)
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
   subroutine set_initial_guesses_abstract_cc_es(solver, wf, X, first, last)
!!
!!    Set initial guesses
!!    Written by Alexander C. Paul, Oct 2020
!!
      implicit none
!
      class(abstract_cc_es) :: solver
!
      class(ccs) :: wf
!
      integer, intent(in) :: first, last
!
      real(dp), dimension(wf%n_es_amplitudes,first:last), intent(out) :: X
!
      integer :: state
!
      do state = first, last
!
         call solver%start_vectors%get(wf, state, X(:,state),  &
                                       solver%energies(state), &
                                       solver%transformation,  &
                                       solver%restart)
!
         if (solver%projector%active) call solver%projector%do_(X(:,state))
!
      enddo 
!
   end subroutine set_initial_guesses_abstract_cc_es
!
!
end module abstract_cc_es_class
