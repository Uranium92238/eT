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
module davidson_cc_response_class
!
!!
!!    Davidson coupled cluster response solver class module
!!    Written by Josefine H. Andersen, 2019
!!  
!
   use kinds
   use file_class
   use ccs_class
   use linear_davidson_tool_class
!
   implicit none
!
   type :: davidson_cc_response
!
      character(len=100) :: tag = 'Davidson coupled cluster response solver'
      character(len=100) :: author = 'Josefine H. Andersen, 2019'
      character(len=500) :: description = 'A Davidson CC response equations solver.'
!
      integer :: max_iterations
!
      real(dp) :: residual_threshold
!
      logical :: restart
!
   contains
!
      procedure :: prepare                         => prepare_davidson_cc_response
      procedure, nopass :: cleanup                 => cleanup_davidson_cc_response
!
      !procedure :: print_banner                    => print_banner_davidson_cc_response
      !procedure :: print_settings                  => print_settings_davidson_cc_response
!
      !procedure, nopass :: print_summary           => print_summary_davidson_cc_response
!
      procedure :: read_settings                   => read_settings_davidson_cc_response
!
      procedure :: run                             => run_davidson_cc_response
!
      procedure, nopass :: set_precondition_vector => set_precondition_vector_davidson_cc_response
!
      !procedure, nopass :: transform_trial_vector  => transform_trial_vector_davidson_cc_response
!
   end type davidson_cc_response
!
!
contains
!
!
   subroutine prepare_davidson_cc_response(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      class(ccs) :: wf
!
!     Print solver banner
!
      !call solver%print_banner()
!
!     Set default settings
!
      solver%max_iterations      = 100
      solver%residual_threshold  = 1.0d-6
      solver%restart             = .false.
!
      !call solver%read_settings()
!
      !call solver%print_settings()
!
      !call wf%initialize_response()
!
   end subroutine prepare_davidson_cc_response
!
!
   subroutine print_settings_davidson_cc_response(solver)
!!
!!    Print settings    
!!    Written by Josefine H. Andersen, March 2019
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      write(output%unit, '(/t3,a)')      '- Davidson CC response solver settings:'
!
      write(output%unit, '(/t6,a26,e9.2)') 'Residual threshold:       ', solver%residual_threshold
!
      write(output%unit, '(t6,a26,i9)')    'Max number of iterations: ', solver%max_iterations
!
   end subroutine print_settings_davidson_cc_response
!
!
   subroutine run_davidson_cc_response(solver, wf)
!!
!!    Run solver
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      class(ccs) :: wf
!
      type(linear_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: RHS, c_i, Xn, frequencies
!

!
   end subroutine run_davidson_cc_response
!
!
   subroutine cleanup_davidson_cc_response(solver)
!!
!!    Cleanup
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      ! save the produced vectors: t^x or M^f
      ! destruct the produced vectors: t^x or M^f
!
   end subroutine cleanup_davidson_cc_response
!
!
   subroutine print_banner_davidson_cc_response(solver)
!!
!!    Print banner
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description,'(t3,a)',.false.,'(t3,a)','(t3,a/)')
!
   end subroutine print_banner_davidson_cc_response
!
!
   subroutine transform_trial_vector_davidson_cc_response(wf, c_i)
!!
!!    Transform trial vector 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Transforms the trial vector according to specified transformation routine.
!!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes, 1), intent(inout) :: c_i
!
      call wf%jacobian_transpose_transform_trial_vector(c_i)
!
   end subroutine transform_trial_vector_davidson_cc_response
!
!
   subroutine set_precondition_vector_davidson_cc_response(wf, davidson)
!!
!!    Set precondition vector
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      type(linear_davidson_tool) :: davidson
!
      real(dp), dimension(:,:), allocatable :: preconditioner
!
      call mem%alloc(preconditioner, wf%n_gs_amplitudes, 1)
      call wf%get_gs_orbital_differences(preconditioner, wf%n_gs_amplitudes)
!
      ! call to a function that subtracts omegas
!
   end subroutine set_precondition_vector_davidson_cc_response
!
!
   subroutine read_settings_davidson_cc_response(solver)
!!
!!    Read settings 
!!    Written by Josefine H. Andersen, 2019
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      integer :: n_specs, i
      character(len=100) :: line
!
      if (.not. requested_section('cc response')) return
!
      call move_to_section('cc response', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:10) == 'threshold:' ) then
!
            read(line(11:100), *) solver%residual_threshold
!
         elseif (line(1:7) == 'restart' ) then
!
            solver%restart = .true.
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         endif
!
      enddo
!
   end subroutine read_settings_davidson_cc_response
!
!
end module davidson_cc_response_class
