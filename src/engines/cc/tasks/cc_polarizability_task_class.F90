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
module cc_polarizability_task_class
!
!!
!! CC polarizability task class
!! Written by Alexander C. Paul, Sarai D. Folkestad, and Eirik F. Kjønstad, 2022
!!
!
   use parameters
   use global_out,                     only: output
   use memory_manager_class,           only: mem
   use ccs_class,                      only: ccs
   use cc_task_class,                  only: cc_task
   use stream_file_class,              only: stream_file
   use cc_eta_xi_calculator_class,     only: cc_eta_xi_calculator
   use cc_eom_eta_xi_calculator_class, only: cc_eom_eta_xi_calculator
   use cc_lr_eta_xi_calculator_class,  only: cc_lr_eta_xi_calculator
   use cc_F_transformation_class,      only: cc_F_transformation
   use cc_null_F_transformation_class, only: cc_null_F_transformation
   use cc_lr_F_transformation_class,   only: cc_lr_F_transformation
!
   implicit none
!
   type, extends(cc_task) :: cc_polarizability_task
!
      logical, private :: eom
!
      real(dp), dimension(:,:), allocatable, private :: xiX
      real(dp), dimension(:,:), allocatable, private :: etaX
!
      type(stream_file), dimension(:,:,:), allocatable, private :: t_response_files
!
      integer, private :: n_frequencies = 0
      real(dp), dimension(:), allocatable, private :: frequencies
!
      logical, dimension(:), allocatable, private   :: compute_t_response
      logical, dimension(:,:), allocatable, private :: compute_polarizability
!
      class(cc_eta_xi_calculator), allocatable, private :: eta_xi_calculator
      class(cc_F_transformation), allocatable, private :: F_transformation
!
   contains
!
      procedure, public :: execute &
                        => execute_cc_polarizability_task
!
      procedure, private :: read_and_initialize_frequencies
!
      procedure, private :: initialize_t_response_files
      procedure, private :: determine_t_responses
!
      procedure, private :: set_polarizability_components_to_compute
      procedure, private :: set_t_response_components_to_compute
!
      procedure, private :: calculate_polarizability
      procedure, private :: read_amplitude_response_vector
!
      final :: destructor
!
   end type cc_polarizability_task
!
!
   interface cc_polarizability_task
!
      procedure :: new_cc_polarizability_task
!
   end interface cc_polarizability_task
!
!
contains
!
!
   function new_cc_polarizability_task() result(this)
!!
!!    New CC polarizability task
!!    Written by Alexander C. Paul, Jan 2022
!!
      use global_in, only: input
!
      implicit none
!
      type(cc_polarizability_task) :: this
!
      this%name_ = 'Determining CC polarizabilities'
!
      this%eom = input%is_keyword_present('eom', 'cc response')
!
   end function new_cc_polarizability_task
!
!
   subroutine execute_cc_polarizability_task(this, wf)
!!
!!    Execute
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(cc_polarizability_task), intent(inout) :: this
      class(ccs), target,            intent(inout) :: wf
!
      call this%print_header()
      call this%start_timer()
!
      call this%read_and_initialize_frequencies()
!
      call this%initialize_t_response_files()
!
      allocate(this%compute_polarizability(3,3))
      allocate(this%compute_t_response(3))
!
      call this%set_polarizability_components_to_compute()
      call this%set_t_response_components_to_compute()
!
      call mem%alloc(this%xiX, wf%n_es_amplitudes, 3)
      call mem%alloc(this%etaX, wf%n_es_amplitudes, 3)
!
      call wf%prepare_for_properties()
!
      if (this%eom) then
!
         this%eta_xi_calculator = cc_eom_eta_xi_calculator(wf)
         this%F_transformation  = cc_null_F_transformation(wf)
!
      else
!
         this%eta_xi_calculator = cc_lr_eta_xi_calculator(wf)
         this%F_transformation  = cc_lr_F_transformation(wf)
!
      endif
!
      call this%eta_xi_calculator%calculate(this%xiX, this%etaX)
!
      call this%determine_t_responses(wf)
      call this%calculate_polarizability(wf)
!
      call mem%dealloc(this%xiX, wf%n_es_amplitudes, 3)
      call mem%dealloc(this%etaX, wf%n_es_amplitudes, 3)
!
      call this%end_timer()
!
   end subroutine execute_cc_polarizability_task
!
!
   subroutine read_and_initialize_frequencies(this)
!!
!!    Read and initialize frequencies
!!    Written Eirik F. Kjønstad, 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_polarizability_task) :: this
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
   subroutine initialize_t_response_files(this)
!!
!!    Initialize t response files
!!    Written Eirik F. Kjønstad, 2022
!!
      implicit none
!
      class(cc_polarizability_task) :: this
!
      integer :: q, freq, sign_
      character(len=200) :: file_name
!
!     We have t responses for a given frequency, component (x,y,z), and sign (+,-)
!
      allocate(this%t_response_files(this%n_frequencies, 3, 2))
!
      do q = 1, 3
         do freq = 1, this%n_frequencies
            do sign_ = 1, 2
!
               write(file_name, '(3(a, i3.3))') 'dipole_t_response_component_', &
                                                q, '_frequency_', freq, '_sign_', sign_
!
               this%t_response_files(freq, q, sign_) = stream_file(file_name)
!
            enddo
         enddo
      enddo
!
   end subroutine initialize_t_response_files
!
!
   subroutine determine_t_responses(this, wf)
!!
!!    Determine t responses
!!    Written Eirik F. Kjønstad, 2022
!!
      use array_initialization, only: copy_and_scale
      use davidson_cc_linear_equations_class, only: davidson_cc_linear_equations
!
      implicit none
!
      class(cc_polarizability_task) :: this
!
      class(ccs) :: wf
!
      real(dp), dimension(:), allocatable :: rhs
!
      integer :: q, sign_
!
      real(dp) :: prefactor
!
      character(len=1) :: sign_character
!
      type(davidson_cc_linear_equations), allocatable :: t_response_solver
!
      call wf%prepare_for_jacobian()
!
      call mem%alloc(rhs, wf%n_es_amplitudes)
!
      do q = 1, 3
!
         if (.not. this%compute_t_response(q)) cycle
!
         call copy_and_scale(-one, this%xiX(:,q), rhs, wf%n_es_amplitudes)
!
         do sign_ = 1, 2
!
            sign_character = '-'
            if (sign_ == 2) sign_character = '+'
!
            call output%printf('v', 'Determining amplitude response', fs='(/t3,a)')
!
            call output%printf('v', 'Asking solver to get response for ' //  &
                               'component (i0) and ((a1))-frequencies.', &
                               ints=[q], chars=[sign_character], fs='(t3,a)')
!
            t_response_solver = davidson_cc_linear_equations(wf,                                   &
                                                             section='cc response',                &
                                                             eq_description='Solving for the       &
                                                             &amplitude response vectors in CC     &
                                                             &response theory.',                   &
                                                             n_frequencies=this%n_frequencies,     &
                                                             n_rhs=1)
!
            prefactor = real((-1)**sign_, kind=dp) ! for frequencies
!
            call t_response_solver%run(wf,                                     &
                                       rhs,                                    &
                                       prefactor*this%frequencies,             &
                                       this%t_response_files(:, q, sign_),     &
                                       'right')
!
            call t_response_solver%cleanup(wf)
!
         enddo
      enddo
!
      call mem%dealloc(rhs, wf%n_es_amplitudes)
!
   end subroutine determine_t_responses
!
!
   subroutine set_polarizability_components_to_compute(this)
!!
!!    Set polarizability components to compute
!!    Written by Alexander C. Paul, May 2021
!!
      use global_in, only: input
!
      implicit none
!
      class(cc_polarizability_task), intent(inout) :: this
!
      integer, dimension(:), allocatable :: polarizabilities
      integer :: k, n_polarizabilities
!
      this%compute_polarizability = .false.
!
      n_polarizabilities = input%get_n_elements_for_keyword('polarizabilities', 'cc response')
!
      if (n_polarizabilities == 0) then
!
!        Not specified means you get all components (xx, xy, xz, ...)
!        of the polarizability
!
         this%compute_polarizability(:,:) = .true.
!
      else
!
!        Figure out which polarizabilities-components have been requested
!        to set true/false in the compute_polarizability array
!
         call mem%alloc(polarizabilities, n_polarizabilities)
!
         call input%get_array_for_keyword('polarizabilities', 'cc response', &
                                          n_polarizabilities, polarizabilities)
!
         do k = 1, n_polarizabilities
!
            if (polarizabilities(k) == 11) then
!
               this%compute_polarizability(1,1) = .true.
!
            elseif (polarizabilities(k) == 12 .or. polarizabilities(k) == 21) then
!
               this%compute_polarizability(1,2) = .true.
               this%compute_polarizability(2,1) = .true.
!
            elseif (polarizabilities(k) == 13 .or. polarizabilities(k) == 31) then
!
               this%compute_polarizability(1,3) = .true.
               this%compute_polarizability(3,1) = .true.
!
            elseif (polarizabilities(k) == 22) then
!
               this%compute_polarizability(2,2) = .true.
!
            elseif (polarizabilities(k) == 23 .or. polarizabilities(k) == 32) then
!
               this%compute_polarizability(2,3) = .true.
               this%compute_polarizability(3,2) = .true.
!
            elseif (polarizabilities(k) == 33) then
!
               this%compute_polarizability(3,3) = .true.
!
            endif
!
         enddo
!
         call mem%dealloc(polarizabilities, n_polarizabilities)
!
      endif
!
   end subroutine set_polarizability_components_to_compute
!
!
   subroutine set_t_response_components_to_compute(this)
!!
!!    Set t response components
!!    Written by Alexander C. Paul, May 2021
!!
!!    Determine which components of the response amplitudes will be
!!    necessary to compute the requested polarizabilities,
!!
      implicit none
!
      class(cc_polarizability_task), intent(inout) :: this
!
      integer :: k, l
!
      this%compute_t_response = .false.
!
      do k = 1, 3
         do l = 1, k
!
            if (this%compute_polarizability(k,l)) then
!
               this%compute_t_response(k) = .true.
               this%compute_t_response(l) = .true.
!
            endif
!
         enddo
      enddo
!
   end subroutine set_t_response_components_to_compute
!
!
   subroutine calculate_polarizability(this, wf)
!!
!!    Calculate polarizability
!!    Written by Josefine H. Andersen, spring 2019
!!
!!       << X, Y >>(omega) = 1/2*[ eta^X (t^Y(omega) + t^Y(-omega))
!!                               + eta^Y (t^X(omega) + t^X(-omega))
!!                               + t^X(-omega) F t^Y(omega)
!!                               + t^X(omega) F t^Y(-omega) ], (*)
!!
!!    for the set of requested frequencies omega. Note that only components
!!    requested by the user are computed; see keyword "polarizabilities".
!!
!!    This routine is called after the amplitude response
!!    vectors (t^X, t^Y)  needed for the polarizability components
!!    have been converged. These are solutions to the equations
!!
!!       (A - omega_k I) t^X(omega_k) = -xi^X
!!
!!    Currently, X and Y are components of the dipole operator.
!!
!!    The terms involving the F matrix (terms 3 and 4 in (*)) are only
!!    added in LR theory. These are thus ignored if EOM polarizabilities
!!    are requested.
!!
!!    New storage of t response and adapted to late-2019
!!    program structure by Eirik F. Kjønstad, Nov 2019.
!!
!!    Revision of cross terms June 2021 by Anna Kristina Schnack-Petersen
!!
      implicit none
!
      class(cc_polarizability_task) :: this
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: tk, tl ! holds amplitude response, components k,l
                                                      ! columns: [+, -] (positive/negative frequency)
!
      real(dp), dimension(:,:), allocatable :: Ftk ! holds F tk
                                                   ! columns: [+ -]
!
      integer :: k, l, freq, sign_
!
      real(dp) :: ddot, polarizability
!
      character(len=4), dimension(3) :: operator_ = ['mu_x', 'mu_y', 'mu_z']
!
      call output%printf('m', 'The convention applied here defines the polarizabilities as &
                              &the response functions, without negative sign.', fs='(t6,a)')
!
!     Allocate arrays to hold amplitude response vectors as well as
!     F-transformed response vectors (if LR)
!
      call mem%alloc(tk, wf%n_es_amplitudes, 2)
      call mem%alloc(tl, wf%n_es_amplitudes, 2)
      call mem%alloc(Ftk, wf%n_es_amplitudes, 2)
!
!     Compute polarizabilities
!
      do freq = 1, this%n_frequencies
!
         do k = 1, 3
!
            if (.not. any(this%compute_polarizability(k,:))) cycle
!
            call this%read_amplitude_response_vector(wf, tk, freq, k)
!
            do sign_ = 1, 2 ! +, -
!
               call this%F_transformation%transform(tk(:,sign_), Ftk(:,sign_))
!
            enddo
!
            do l = 1, k
!
               if (this%compute_polarizability(k,l)) then
!
                  polarizability = zero
!
                  call this%read_amplitude_response_vector(wf, tl, freq, l)
!
                  do sign_ = 1, 2 ! +, -
!
                     polarizability = polarizability +                            &
                         half*ddot(wf%n_es_amplitudes, this%etaX(:,k), 1, tl(:,sign_), 1) + &
                         half*ddot(wf%n_es_amplitudes, this%etaX(:,l), 1, tk(:,sign_), 1)
!
                  enddo
!
                  polarizability = polarizability +                  &
                        half*ddot(wf%n_es_amplitudes, Ftk(:,1), 1, tl(:,2), 1) + &
                        half*ddot(wf%n_es_amplitudes, Ftk(:,2), 1, tl(:,1), 1)
!
                  call output%printf('m', '<< ' // operator_(k) // ', ' //  &
                                     operator_(l) // ' >>' // '((e8.2)): (f19.12)', &
                                     reals=[this%frequencies(freq), &
                                     polarizability], fs='(t6,a)')
!
               endif
!
            enddo
!
         enddo
!
      enddo
!
      call mem%dealloc(tk, wf%n_es_amplitudes, 2)
      call mem%dealloc(tl, wf%n_es_amplitudes, 2)
      call mem%dealloc(Ftk, wf%n_es_amplitudes, 2)
!
   end subroutine calculate_polarizability
!
!
   subroutine read_amplitude_response_vector(this, wf, t, freq, k)
!
      implicit none
!
      class(cc_polarizability_task) :: this
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes, 2) :: t
      integer :: k, sign_,freq
!
      do sign_ = 1, 2
!
         call this%t_response_files(freq, k, sign_)%open_('rewind')
         call this%t_response_files(freq, k, sign_)%read_(t(:,sign_), wf%n_es_amplitudes)
         call this%t_response_files(freq, k, sign_)%close_()
!
      enddo
!
   end subroutine read_amplitude_response_vector
!
!
   subroutine destructor(this)
!!
!!    Destructor
!!    Written by Eirik F. Kjønstad, 2022
!!
      implicit none
!
      type(cc_polarizability_task) :: this
!
      if (allocated(this%frequencies)) &
         call mem%dealloc(this%frequencies, this%n_frequencies)
!
   end subroutine destructor
!
!
end module cc_polarizability_task_class
