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
module cc_eom_transition_moments_task_class
!
!!
!! CC EOM transition moments task class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use parameters
   use ccs_class, only: ccs
   use cc_transition_moments_task_class, only: cc_transition_moments_task
   use memory_manager_class, only: mem
   use global_out, only: output
!
   implicit none
!
   type, extends(cc_transition_moments_task) :: cc_eom_transition_moments_task
!
      integer, private :: n_initial_states
      integer, dimension(:), allocatable, private :: initial_states
      logical, private :: transition_moments, permanent_moments, dipole_length
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_eom_transition_moments_task
!
      procedure, public :: cleanup &
                        => cleanup_cc_eom_transition_moments_task
!
      procedure, private :: print_permanent_moments_summary
      procedure, private :: print_operator_for_initial_states
      procedure, private :: set_initial_states
      procedure, private :: read_settings
!
   end type cc_eom_transition_moments_task
!
!
   interface cc_eom_transition_moments_task
!
      procedure :: new_cc_eom_transition_moments_task
!
   end interface cc_eom_transition_moments_task
!
!
contains
!
!
   function new_cc_eom_transition_moments_task() result(this)
!!
!!    New CC EOM transition moments task
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(cc_eom_transition_moments_task) :: this
!
      this%name_ = 'Determining CC EOM transition moments'
!
   end function new_cc_eom_transition_moments_task
!
!
   subroutine execute_cc_eom_transition_moments_task(this, wf)
!!
!!    Execute
!!    Written by Josefine H. Andersen, Sarai D. Folkestad
!!    and Alexander C. Paul, June 2019
!!
!!    Computes the EOM dipole transition moments using transition densities
!!    and dipole moment integrals.
!!
!!    Restructured by Alexander C. Paul and Sarai D. Folkestad, Apr 2020
!!    Restructured for general (GS and ES) transition momnets and state
!!    properties and to allow for the calculation of transition moments
!!    of states generated in another calculation.
!!
      use timings_class, only: timings
!
      implicit none
!
      class(cc_eom_transition_moments_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      real(dp), dimension(:,:,:), allocatable :: operator_
      real(dp), dimension(:,:,:), allocatable :: dipole_matrix
!
      call this%print_header()
      call this%start_timer()
!
      call this%read_settings()
!
      call wf%prepare_for_properties()
      call wf%initialize_gs_density()
      call wf%construct_gs_density()
      call wf%initialize_density_intermediates()
!
      call output%printf('m', ':: EOM properties calculation', fs='(/t3,a)')
!
      if (this%dipole_length) then
!
         call mem%alloc(operator_, wf%n_mo, wf%n_mo, 3)
         call wf%get_t1_oei('dipole', operator_)
!
      endif
!
      call mem%alloc(dipole_matrix, wf%n_singlet_states+1, wf%n_singlet_states+1, 3)
!
      call wf%compute_eom_transition_moments(operator_,               &
                                             dipole_matrix,           &
                                             this%n_initial_states,   &
                                             this%initial_states,     &
                                             this%transition_moments, &
                                             this%permanent_moments)
!
      if (this%transition_moments) then
!
         call this%print_transition_moment_summary(wf, dipole_matrix, &
                                                   this%initial_states, 'EOM')
!
      else
!
!        Print transition moments from GS to ES also
!        if only permanent moments are requested
         call this%print_transition_moment_summary(wf, dipole_matrix, &
                                                   [0], 'EOM')
      end if
!
      if (this%permanent_moments) then
         call this%print_permanent_moments_summary(wf, 'dipole', dipole_matrix, 3)
      end if
!
      call mem%dealloc(dipole_matrix, wf%n_singlet_states+1, wf%n_singlet_states+1, 3)
      call mem%dealloc(operator_, wf%n_mo, wf%n_mo, 3)
!
      call wf%destruct_density_intermediates
!
      call this%end_timer()
!
   end subroutine execute_cc_eom_transition_moments_task
!
!
   subroutine set_initial_states(this, n_states)
!!
!!    Set initial states
!!    Written Alexander C. Paul and Sarai D. Folkestad, Apr 2021
!!
!!    Determine the states for which properties shall be calculated
!!    if we are doing transition moments or permanent moments:
!!
!!       Example of usage:
!!       initial states: {0, 1, 2}
!!       or: [0,2]
!!
!!       gives all transitions from the ground state
!!          0 -> 1, 2, 3, ..
!!       and excitations between excited state 1 and 2 (ordered according to energy)
!!       to all other excited states
!!          1 -> 2, 3, ...
!!          2 -> 3, 4, ...
!!       or properties of the excited states 1 and 2
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_eom_transition_moments_task), intent(inout) :: this
      integer, intent(in) :: n_states
!
      integer, dimension(:), allocatable :: initial_states
!
      if (input%is_keyword_present('initial states','cc response')) then
!
         this%n_initial_states = input%get_n_elements_for_keyword('initial states', &
                                                                    'cc response')
!
         call mem%alloc(initial_states, this%n_initial_states)
!
         call input%get_array_for_keyword('initial states', 'cc response', &
                                          this%n_initial_states, initial_states)
!
!        Transitions from the ground state shall always be considered
         if (any(initial_states == 0)) then
!
            call mem%alloc(this%initial_states, this%n_initial_states)
!
            this%initial_states = initial_states
!
            call mem%dealloc(initial_states, this%n_initial_states)
!
         else
!
!           Copy and add 0 as initial state
!
            call mem%alloc(this%initial_states, this%n_initial_states + 1)
!
            this%initial_states(1)  = 0
            this%initial_states(2:) = initial_states
!
            call mem%dealloc(initial_states, this%n_initial_states)
!
            this%n_initial_states = this%n_initial_states + 1
!
         endif
!
         if (this%n_initial_states > n_states + 1) then
            call output%error_msg('Requested properties for more states &
                                  &than requested in the calculation.')
         end if
!
         if (any(this%initial_states > n_states)) then
            call output%error_msg('Requested properties for state with a number &
                                  &larger than the number of states')
         end if
!
         if (any(this%initial_states < 0)) then
            call output%error_msg('Requested properties of state with a negative number.')
         end if
!
      else
!
!        Only calculate transitions from ground state (state 0)
!
         this%n_initial_states = 1
!
         call mem%alloc(this%initial_states, this%n_initial_states)
!
         this%initial_states(1) = 0
!
      endif
!
   end subroutine set_initial_states
!
!
   subroutine print_permanent_moments_summary(this, wf, operator_type, &
                                              electronic, n_components)
!!
!!    Print permanent moments summary
!!    Written by Alexander C. Paul, May 2020
!!
      use global_out, only: output
!
      implicit none
!
      class(cc_eom_transition_moments_task), intent(in) :: this
!
      class(ccs), intent(in) :: wf
!
      character(len=*), intent(in) :: operator_type
      integer, intent(in) :: n_components
!
      real(dp), dimension(wf%n_singlet_states+1, wf%n_singlet_states+1, n_components), &
                                                                           intent(in) :: electronic
!
      real(dp), dimension(:), allocatable :: nuclear
!
      character(len=5), dimension(n_components) :: components
!
      call output%printf('m', '- Summary of EOM permanent moments calculation:', fs='(/t3,a)')
!
      call mem%alloc(nuclear, n_components)
!
      if (operator_type == 'dipole') then
!
         call output%printf('m', 'Total permanent dipole moments in [a.u.]:', fs='(/t6,a)')
         call output%print_separator('m', 41, '=', fs='(t6,a)')
!
         call output%printf('m', 'Conversion factor from au to Debye: (f11.9)', &
                             reals=[au_to_debye], fs='(/t6,a)')
!
         nuclear = wf%get_nuclear_dipole()
!
         components = ['X', 'Y', 'Z']
!
      end if
!
      call this%print_operator_for_initial_states(wf, nuclear, electronic, &
                                                  components, n_components)
!
      call mem%dealloc(nuclear, n_components)
!
   end subroutine print_permanent_moments_summary
!
!
   subroutine print_operator_for_initial_states(this, wf, nuclear, electronic, &
                                                components, n_components)
!!
!!    Print operator for initial states
!!    Written by Alexander C. Paul, May 2022
!!
      implicit none
!
      class(cc_eom_transition_moments_task), intent(in) :: this
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(wf%n_singlet_states+1, wf%n_singlet_states+1, n_components), &
                                                                           intent(in) :: electronic
!
      real(dp), dimension(n_components), intent(in) :: nuclear
!
      character(len=5), dimension(n_components), intent(in) :: components
!
      character(len=18), dimension(this%n_initial_states + 1) :: column_labels
!
      real(dp), dimension(:,:), allocatable :: data_
!
      integer :: state_i, i, c
!
      call mem%alloc(data_, n_components, this%n_initial_states+1)
!
      column_labels(1) = "Nuclear"
      data_(:,1) = nuclear
!
!     Select only the states that shall be printed
!
      do i = 1, this%n_initial_states
!
         state_i = this%initial_states(i)
!
         write(column_labels(i + 1), '(a,1x,i0)')  "State", state_i
!
         do c = 1, n_components
!
            data_(c, i+1) = nuclear(c) + electronic(state_i+1, state_i+1, c)
!
         end do
!
      end do
!
      call output%print_table("Comp.", components, column_labels, data_, &
                              n_components, this%n_initial_states+1)
!
      call mem%dealloc(data_, n_components, this%n_initial_states+1)
!
   end subroutine print_operator_for_initial_states
!
!
   subroutine read_settings(this)
!!
!!    Read settings
!!    Written by Alexander C. Paul, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_eom_transition_moments_task) :: this
!
      integer :: n_states
!
      this%dipole_length = input%is_keyword_present('dipole length','cc response')
      this%permanent_moments = input%is_keyword_present('permanent moments','cc response')
      this%transition_moments = input%is_keyword_present('transition moments','cc response')
!
      call input%get_keyword('singlet states','solver cc es', n_states)
      call this%set_initial_states(n_states)
!
      if (.not. this%dipole_length) &
         call output%error_msg('no operator selected in response calculation')
!
   end subroutine read_settings
!
!
   subroutine cleanup_cc_eom_transition_moments_task(this)
!!
!!    Cleanup
!!    Written by Alexander C. Paul, Jan 2022
!!
      implicit none
!
      class(cc_eom_transition_moments_task), intent(inout) :: this
!
      call mem%dealloc(this%initial_states, this%n_initial_states)
!
   end subroutine cleanup_cc_eom_transition_moments_task
!
!
end module cc_eom_transition_moments_task_class
