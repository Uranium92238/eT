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
module cc_amplitude_response_task_class
!
!!
!! CC amplitude response task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use parameters
   use global_in,             only: input
   use global_out,            only: output
   use ccs_class,             only: ccs
   use cc_task_class,         only: cc_task
   use sequential_file_class, only: sequential_file
   use memory_manager_class,  only: mem
!
   implicit none
!
   type, extends(cc_task) :: cc_amplitude_response_task
!
      integer, private :: n_frequencies
      real(dp), dimension(:), allocatable, private :: frequencies
!
      type(sequential_file), dimension(:,:,:), allocatable :: t_responses
      real(dp), dimension(:,:), allocatable :: xiX
!
      logical, dimension(3), private :: calculate_component
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_amplitude_response_task
!
      procedure, private :: read_and_initialize_frequencies
      procedure, private :: solve_response_equations
!
   end type cc_amplitude_response_task
!
!
   interface cc_amplitude_response_task
!
      procedure :: new_cc_amplitude_response_task
!
   end interface cc_amplitude_response_task
!
!
contains
!
!
   function new_cc_amplitude_response_task(calculate_component) result(this)
!!
!!    New CC amplitude response
!!    Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
      implicit none
!
      logical, dimension(3), intent(in) :: calculate_component
!
      type(cc_amplitude_response_task) :: this
!
      this%name_ = 'Determining CC amplitude response'
!
      this%calculate_component = calculate_component
!
   end function new_cc_amplitude_response_task
!
!
   subroutine read_and_initialize_frequencies(this)
!!
!!    Read and intialize frequencies
!!    Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
      implicit none
!
      class(cc_amplitude_response_task) :: this
!
      if (input%is_keyword_present('frequencies', 'cc response')) then
!
         this%n_frequencies = input%get_n_elements_for_keyword('frequencies', 'cc response')
!
         call mem%alloc(this%frequencies, this%n_frequencies)
!
         call input%get_array_for_keyword('frequencies', 'cc response', &
                                          this%n_frequencies, this%frequencies)
!
      else
!
         call output%error_msg('Requested polarizabilities, but no &
                               &frequencies were provided!')
!
      endif
!
   end subroutine read_and_initialize_frequencies
!
!
   subroutine execute_cc_amplitude_response_task(this, wf)
!!
!!    Calculate
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_amplitude_response_task), intent(inout) :: this
!
      class(ccs), intent(inout), target :: wf
!
      real(dp), dimension(:,:,:), allocatable :: X
!
      integer :: k
!
      call this%print_header()
      call this%start_timer()
!
      call this%read_and_initialize_frequencies()
!
      call mem%alloc(X, wf%n_mo, wf%n_mo, 3)
!
      call wf%get_t1_oei('dipole', X)
!
      call mem%alloc(this%xiX, wf%n_es_amplitudes, 3)
!
      do k = 1, 3
         call wf%construct_xiX(X(:,:,k), this%xiX(:,k))
      end do
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo, 3)
!
      call wf%prepare_for_jacobian()
!
      call this%solve_response_equations(wf)
!
      call mem%dealloc(this%xiX, wf%n_es_amplitudes, 3)
!
      call this%end_timer()
!
   end subroutine execute_cc_amplitude_response_task
!
!
   subroutine solve_response_equations(this, wf)
!!
!!    Solve response equations
!!    Written by Eirik F. Kjønstad and Josefine H. Andersen, 2019
!!
!!    Makes preparations for determining the amplitude response t^X(omega_k) and t^X(-omega_k),
!!    where X is the dipole components (x,y,z) and omega_k the set of frequencies for which
!!    the polarizability is to be evaluated.
!!
!!    These responses are solutions of the equations
!!
!!       (A - omega_k I) t^X(omega_k) = -xi^X
!!
!!    The amplitude response t^X(+- omega_k) are stored to file. In particular,
!!    these are stored in a sequential file array:
!!
!!       this%t_responses(freq, component, sign_), where:
!!
!!          freq:       1,2,3,...; denotes the frequency number as specified in input
!!                      e.g., if freq = 2, the file contains an amplitude response vector
!!                      for the second frequency specified on input
!!
!!          component:  1,2,3 (x,y,z); the component of the dipole moment vector
!!
!!          sign_:      1,2 (+, -); whether it is t^X(+omega_k) or t^X(-omega_k)
!!
!!    This wrapper routine, and the file storage for amplitude response vectors, was made by
!!    Eirik F. Kjønstad, Nov 2019. It is adapted/based on the general structure set up to solve
!!    amplitude response originally written by Josefine H. Andersen, spring 2019.
!!
      use davidson_cc_linear_equations_class, only: davidson_cc_linear_equations
!
      implicit none
!
      class(cc_amplitude_response_task) :: this
!
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable     :: rhs ! right-hand-side
!
      integer :: k, freq, sign_
!
      real(dp) :: prefactor
!
      character(len=200) :: file_name
      character(len=1)   :: sign_character
!
      class(davidson_cc_linear_equations), allocatable :: t_response_solver
!
!     Initialize amplitude response files for storing solutions
!
      allocate(this%t_responses(this%n_frequencies, 3, 2))
!
      do k = 1, 3
         do freq = 1, this%n_frequencies
            do sign_ = 1, 2
!
               write(file_name, '(a, i3.3, a, i3.3, a, i3.3)') 'dipole_t_response_component_', &
                                                            k, '_frequency_', freq, '_sign_', sign_
!
               this%t_responses(freq,k,sign_) = sequential_file(file_name)
!
            enddo
         enddo
      enddo
!
!     Solve response equations for one component and all frequencies
!
      call mem%alloc(rhs, wf%n_es_amplitudes)
!
      do k = 1, 3
!
!        kth component not needed for requested polarizability calculations
         if (.not. this%calculate_component(k)) cycle
!
         call copy_and_scale(-one, this%xiX(:,k), rhs, wf%n_es_amplitudes)
!
         do sign_ = 1, 2
!
            sign_character = '-'
            if (sign_ == 2) sign_character = '+'
!
            call output%printf('v', 'Determining amplitude response...', fs='(/t3,a)')
!
            call output%printf('v', 'Asking solver to get response for ' //  &
                               'component (i0) and ((a1))-frequencies.', &
                               ints=[k], chars=[sign_character], fs='(t3,a)')
!
            t_response_solver = davidson_cc_linear_equations(wf,                                   &
                                                             section='cc response',                &
                                                             eq_description='Solving for the       &
                                                             &amplitude response vectors in CC     &
                                                             &response theory.',                   &
                                                             n_frequencies=this%n_frequencies,   &
                                                             n_rhs=1)
!
            prefactor = real((-1)**sign_, kind=dp) ! for frequencies
!
            call t_response_solver%run(wf, rhs, prefactor*this%frequencies, &
                              this%t_responses(:,k,sign_), 'right')
!
            call t_response_solver%cleanup(wf)
!
         enddo
      enddo
!
      call mem%dealloc(rhs, wf%n_es_amplitudes)
!
   end subroutine solve_response_equations
!
!
end module cc_amplitude_response_task_class
