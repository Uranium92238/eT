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
module abstract_cc_es_class
!
!!
!!    Abstract coupled cluster excited state solver class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, June 2019
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
   use es_projection_tool_class, only : es_projection_tool
   use es_valence_projection_tool_class, only: es_valence_projection_tool
   use es_cvs_projection_tool_class, only: es_cvs_projection_tool
   use es_ip_projection_tool_class, only: es_ip_projection_tool
!
   implicit none
!
   type, abstract :: abstract_cc_es
!
      character(len=100) :: name_ 
      character(len=100) :: tag 
      character(len=100) :: author 
!
      character(len=500) :: description1
      character(len=500) :: description2 
!
      integer :: max_iterations
!
      real(dp) :: eigenvalue_threshold  
      real(dp) :: residual_threshold  
!
      character(len=200) :: storage 
      logical :: records_in_memory 
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
      class(es_start_vector_tool), allocatable  :: start_vector_tool
      class(es_projection_tool), allocatable    :: projection_tool
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
      procedure :: determine_restart_transformation => determine_restart_transformation_abstract_cc_es 
!
      procedure :: initialize_start_vector_tool     => initialize_start_vector_tool_abstract_cc_es
      procedure :: initialize_projection_tool       => initialize_projection_tool_abstract_cc_es
!
      procedure :: initialize_energies              => initialize_energies_abstract_cc_es
      procedure :: destruct_energies                => destruct_energies_abstract_cc_es
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
      call output%printf(':: ' // solver%name_, pl='m', fs='(//t3,a)')
      call output%printf(':: ' // solver%author, pl='m', fs='(t3,a)')
      call output%printf(solver%description1, pl='n', ffs='(/t3,a)', fs='(t3,a)')
      call output%printf(solver%description2, pl='n', ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_abstract_cc_es
!
!
   subroutine read_es_settings_abstract_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(abstract_cc_es) :: solver 
!
      call input%get_keyword_in_section('residual threshold', 'solver cc es', solver%residual_threshold)
      call input%get_keyword_in_section('energy threshold', 'solver cc es', solver%eigenvalue_threshold)
      call input%get_keyword_in_section('max iterations', 'solver cc es', solver%max_iterations)
!               
      call input%get_required_keyword_in_section('singlet states', 'solver cc es', solver%n_singlet_states)
!
      if (input%requested_keyword_in_section('restart', 'solver cc es')) solver%restart = .true.    
      if (input%requested_keyword_in_section('left eigenvectors', 'solver cc es')) solver%transformation = 'left'    
      if (input%requested_keyword_in_section('right eigenvectors', 'solver cc es')) solver%transformation = 'right'    
!
      if (input%requested_keyword_in_section('core excitation', 'solver cc es') .and. .not. &
          input%requested_keyword_in_section('ionization', 'solver cc es')) solver%es_type = 'core'
!
      if (input%requested_keyword_in_section('ionization', 'solver cc es') .and. .not. &
          input%requested_keyword_in_section('core excitation', 'solver cc es')) solver%es_type = 'ionize'
!
      call input%get_keyword_in_section('storage', 'solver cc es', solver%storage)
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
      call output%printf('- Settings for coupled cluster excited state solver (' //trim(solver%tag) // '):', &
                              pl='m', fs='(/t3,a)')
!
      call output%printf('Energy threshold:     (e9.2)', reals=[solver%eigenvalue_threshold], pl='m', fs='(/t6,a)')
      call output%printf('Residual threshold:   (e9.2)', reals=[solver%residual_threshold], pl='m', fs='(/t6,a)')
!
      call output%printf('Number of singlet states:  (i6)', ints=[solver%n_singlet_states], pl='m', fs='(/t6,a)')
      call output%printf('Max number of iterations:  (i6)', ints=[solver%max_iterations], pl='m', fs='(t6,a)')
      call output%printf('Calculation type:          (a0)', chars=[trim(solver%es_type)], pl='m', fs='(t6,a)')
!
      call output%printf('Solving for the (a0) eigenvectors.', chars=[trim(solver%transformation)], pl='m', fs='(/t6,a)')
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
!
      call solver%timer%turn_off()
!
     call output%printf('- Finished solving the ' // trim(convert_to_uppercase(wf%name_)) &
                 // ' excited state equations ('// trim(solver%transformation) //')', pl='m', fs='(/t3,a)')
!
      call output%printf('Total wall time (sec): (f20.5)', reals=[solver%timer%get_elapsed_time('wall')], pl='m', fs='(/t6,a)')
      call output%printf('Total cpu time (sec):  (f20.5)', reals=[solver%timer%get_elapsed_time('cpu')], pl='m', fs='(t6,a)')
!
   end subroutine cleanup_abstract_cc_es
!
!
   subroutine print_summary_abstract_cc_es(solver, wf, X)
!!
!!    Print summary 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(abstract_cc_es), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes, solver%n_singlet_states), intent(in) :: X
!
      integer :: state 
!
      character(len=1) :: label ! R or L, depending on whether left or right transformation 
!
      label = trim(adjustl(convert_to_uppercase(solver%transformation(1:1))))
!
      call output%printf('- Excitation vector amplitudes:', pl='n', fs='(/t3,a)')
!
      do state = 1, solver%n_singlet_states
!
         call output%printf('Electronic state nr. (i0)', pl='n', ints=[state], fs='(/t6,a)')
!
         call output%printf('Energy (Hartree):             (f19.12)', pl='n', &
                                 reals=[solver%energies(state)], fs='(/t6,a)')
!
         call output%printf('Fraction singles (|(a0)1|/|(a0)|):  (f19.12)', pl='n', &
            reals=[get_l2_norm(X(1:wf%n_t1,state),wf%n_t1)/get_l2_norm(X(:,state),wf%n_es_amplitudes)], &
            chars=[label, label], fs='(t6,a)')
!
         call wf%print_dominant_x_amplitudes(X(1, state), label)
!
      enddo 
!
      call output%printf('- Electronic excitation energies:', pl='m', fs='(/t6,a)')
!
      call output%printf('                                 Excitation energy             ', pl='m', fs='(/t6,a)')
      call output%printf('                     ------------------------------------------', pl='m', fs='(t6,a)')
      call output%printf('State                (Hartree)             (eV)                ', pl='m', fs='(t6,a)')
      call output%printf('---------------------------------------------------------------', pl='m', fs='(t6,a)')
!
      do state = 1, solver%n_singlet_states
!
         call output%printf('(i2)              (f19.12)    (f19.12)', pl='m', fs='(t6,a)', &
                     ints=[state], reals=[solver%energies(state), solver%energies(state)*Hartree_to_eV])
!
      enddo 
!
      call output%printf('---------------------------------------------------------------', pl='m', fs='(t6,a)')
      call output%printf('eV/Hartree (CODATA 2014): (f11.8)', pl='m', fs='(t6,a)', reals=[Hartree_to_eV])
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
      class(abstract_cc_es), intent(in)   :: solver 
      class(ccs), intent(inout)           :: wf
!
      if (solver%transformation == 'right') then 
!
         call wf%initialize_excited_state_files('right')
         call wf%prepare_for_jacobian()
!
      else if (solver%transformation == 'left') then 
!
         call wf%initialize_excited_state_files('left')
         call wf%prepare_for_jacobian_transpose()
!
      else if (solver%transformation == 'both') then 
!
         call wf%initialize_excited_state_files('both')
         call wf%prepare_for_jacobian()
         call wf%prepare_for_jacobian_transpose()
!
      end if
!
   end subroutine prepare_wf_for_excited_state_abstract_cc_es
!
!
   subroutine determine_restart_transformation_abstract_cc_es(solver, wf)
!!
!!    Determine number of states on file 
!!    Written by Eirik F. Kjønstad, June 2019
!!
      implicit none 
!
      class(abstract_cc_es), intent(inout) :: solver
!
      class(ccs), intent(in) :: wf 
!
      integer :: n_left_vectors_on_file, n_right_vectors_on_file
!
      n_left_vectors_on_file = wf%get_n_excited_states_on_file('left')
      n_right_vectors_on_file = wf%get_n_excited_states_on_file('right')
!
      if (solver%transformation == 'right') then 
!
         if (n_right_vectors_on_file > 0) then 
!
            solver%restart_transformation = 'right'
!
         elseif (n_right_vectors_on_file == 0 .and. n_left_vectors_on_file > 0) then 
!
            solver%restart_transformation = 'left'
!
         else
!
            call output%error_msg('Could not restart excited state calculation.')
!
         endif 
!
      elseif (solver%transformation == 'left') then
!
         if (n_left_vectors_on_file > 0) then 
!
            solver%restart_transformation = 'left'
!
         elseif (n_left_vectors_on_file == 0 .and. n_right_vectors_on_file > 0) then 
!
            solver%restart_transformation = 'right'
!
         else
!
            call output%error_msg('Could not restart excited state calculation.')
!
         endif          
!
      else 
!
         call output%error_msg('Could not restart excited state calculation.')
!
      endif  
!
   end subroutine determine_restart_transformation_abstract_cc_es
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
         solver%start_vector_tool = es_valence_start_vector_tool(wf)
!
      elseif (trim(solver%es_type) == 'core') then 
!
         call wf%read_cvs_settings()
         solver%start_vector_tool = es_cvs_start_vector_tool(wf)
!
      elseif (trim(solver%es_type) == 'ionize') then 
!
         solver%start_vector_tool = es_ip_start_vector_tool(wf)
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
         solver%projection_tool = es_valence_projection_tool()
!
      elseif (trim(solver%es_type) == 'core') then 
!
         solver%projection_tool = es_cvs_projection_tool(wf)
!
      elseif (trim(solver%es_type) == 'ionize') then 
!
         solver%projection_tool = es_ip_projection_tool(wf)
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
end module abstract_cc_es_class
