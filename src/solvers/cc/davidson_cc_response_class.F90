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
      logical :: moments, polarizability
      logical :: restart
!
      integer :: dim_rhs = 1, n_freq = 0
!
   contains
!
      procedure :: prepare                         => prepare_davidson_cc_response
      procedure :: cleanup                         => cleanup_davidson_cc_response
!
      procedure :: print_banner                    => print_banner_davidson_cc_response
      procedure :: print_settings                  => print_settings_davidson_cc_response
!
      !procedure, nopass :: print_summary           => print_summary_davidson_cc_response
!
      procedure :: read_settings                   => read_settings_davidson_cc_response
!
      procedure :: run                             => run_davidson_cc_response
!
      procedure, nopass :: set_precondition_matrix => set_precondition_matrix_davidson_cc_response
!
      !procedure, nopass :: transform_trial_vector  => transform_trial_vector_davidson_cc_response
!
      procedure :: construct_rhs                   => construct_rhs_davidson_cc_response
      procedure :: build_fr_matrix                 => build_fr_matrix_davidson_cc_response
!
      procedure :: get_frequencies                 => get_frequencies_davidson_cc_response
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
      call solver%print_banner()
!
!     Set default settings
!
      solver%max_iterations      = 100
      solver%residual_threshold  = 1.0d-6
      solver%restart             = .false.
      solver%moments             = .false.
      solver%polarizability      = .false.
!
      call solver%read_settings()
!
      call solver%print_settings()
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
      write(output%unit, '(t6,a26,i9)')    'Number of frequencies:    ', solver%n_freq
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
      logical :: converged_residual
!
      real(dp), dimension(:,:), allocatable :: rhs, c_i, frequencies
!
      integer :: iteration
!
      real(dp) :: residual_norm, ddot, norm_trial
!
!     Get frequencies from file or input
!
      call mem%alloc(frequencies, solver%n_freq, 1)
      call solver%get_frequencies(wf, frequencies)
!
!     Get right-hand-side vector
!
      call mem%alloc(rhs, wf%n_es_amplitudes, solver%dim_rhs)
      call solver%construct_rhs(wf, rhs)
!
      call davidson%prepare_response('response', wf%n_es_amplitudes, solver%residual_threshold, rhs, solver%dim_rhs)
!
      call solver%set_precondition_matrix(wf, davidson, frequencies, solver%n_freq)
!
!     Set initial guess: RHS_mu / (eps_mu - omega)
!
      call davidson%precondition_response(rhs, solver%dim_rhs, solver%n_freq)
!
      !norm_trial = sqrt(ddot(wf%n_es_amplitudes, eta, 1, eta, 1)) 
!
      call mem%dealloc(frequencies, solver%n_freq, 1)
      call mem%dealloc(rhs, wf%n_es_amplitudes, solver%dim_rhs)
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
      write(output%unit, '(/t3,a,a,a)') 'Cleaning up ', trim(solver%tag), '.'
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
   subroutine construct_rhs_davidson_cc_response(solver, wf, rhs)
!!
!!    Construct right-hand-side vector.
!!    Written by Josefine H. Andersen, March 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      class(davidson_cc_response) :: solver
!
      real(dp), dimension(wf%n_es_amplitudes, solver%dim_rhs), intent(inout) :: rhs
!
      if (solver%moments) then
!
         !call solver%build_fr_matrix(wf, rhs)
         call wf%build_fr_matrix(rhs, solver%dim_rhs)
!
      elseif (solver%polarizability) then
!
         call wf%construct_csiX('dipole_length', rhs)
!
      endif
!
   end subroutine construct_rhs_davidson_cc_response
!
!
   subroutine get_frequencies_davidson_cc_response(solver, wf, freq)
!!
!!    Get frequencies
!!    Written by Josefine H. Andersen, Mar 2019
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      class(ccs) :: wf
!
      real(dp), dimension(solver%n_freq, 1), intent(out) :: freq
!
      if (solver%moments) then
!
         call wf%read_excitation_energies(solver%n_freq, freq)
!
         call dscal(wf%n_es_amplitudes, -one, freq, 1)
!
      elseif (solver%polarizability) then
!
         !call solver%read_freq_from_input(freq)
!
      endif
!
   end subroutine get_frequencies_davidson_cc_response
!
!
   subroutine build_fr_matrix_davidson_cc_response(solver, wf, fr)
!!
!!    Build matrix with F-transformed right vectors i as columns
!!    Written by Josefine H. Andersen
!!
      implicit none
!
      class(davidson_cc_response) :: solver
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes, solver%dim_rhs), intent(out) :: fr
!
      !real(dp), dimension(:,:), allocatable :: r_n
      real(dp), dimension(:,:), allocatable :: R_n
!
      real(dp) :: ddot
!
      integer :: i = 1, state = 1
!
      call mem%alloc(R_n, wf%n_es_amplitudes, 1)
      R_n = zero
      !call mem%alloc(r_n, wf%n_es_amplitudes, 1)

      call wf%read_excited_state(R_n, state, 'right')
!
      write(output%unit,'(t6,a,f19.10)') 'R vector norm = ', &
      ddot(wf%n_es_amplitudes, R_n, 1, R_n, 1)
      flush(output%unit)
!
      call wf%F_transform_vector(R_n)

!
      !do i = 1, solver%dim_rhs
!
         !call wf%read_excited_state(r_n, 1, 'right')
!
         !call wf%F_transform_vector(r_n)
!
         !call daxpy(wf%n_es_amplitudes, one, r_n, 1, fr(:,i), 1)
!
      !enddo
!
      call mem%dealloc(r_n, wf%n_es_amplitudes, 1)
!
   end subroutine build_fr_matrix_davidson_cc_response
!
!
   subroutine transform_trial_vector_davidson_cc_response(wf, c_i)
!!
!!    Transform trial vector 
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
   subroutine set_precondition_matrix_davidson_cc_response(wf, davidson, freq_vec, n_freq)
!!
!!    Set precondition matrix  
!!    Written by Josefine H. Andersen, 2019
!!    
!!       preconditioner = (epsilon_a - epsilon_i) - omega_j
!!
!!    One column per frequency
!!
      implicit none
!
      class(ccs) :: wf
!
      type(linear_davidson_tool) :: davidson
!
      integer, intent(in) :: n_freq
!
      real(dp), dimension(n_freq, 1), intent(in) :: freq_vec
!
      real(dp), dimension(:,:), allocatable :: preconditioner
!
      real(dp), dimension(:,:), allocatable :: orb_diff
!
      integer :: i
!
      call mem%alloc(orb_diff, wf%n_gs_amplitudes, 1)
      call wf%get_gs_orbital_differences(preconditioner, wf%n_gs_amplitudes)
!
!     Loop through frequencies to generate n_freq precondition vectors
!
      call mem%alloc(preconditioner, wf%n_gs_amplitudes, n_freq)
!
      do i = 1, n_freq
!
         call wf%get_gs_orbital_differences(preconditioner(:,i), wf%n_gs_amplitudes)
!
         preconditioner(:, i) = orb_diff(:, 1) - freq_vec(i, 1)
!
      enddo
!
      call davidson%set_response_preconditioner(preconditioner, n_freq)
      call mem%dealloc(preconditioner, wf%n_gs_amplitudes, n_freq)
!
   end subroutine set_precondition_matrix_davidson_cc_response
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
!     Read response section
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
         if (line(1:19) == 'residual threshold:' ) then
!
            read(line(20:100), *) solver%residual_threshold
!
         elseif (line(1:15) == 'max iterations:' ) then
!
            read(line(16:100), *) solver%max_iterations
!
         endif
!
      enddo
!
!     Read property section to get requested property
!
      if (.not. requested_section('cc properties')) return
!
      call move_to_section('cc properties', n_specs)
!
      do i = 1, n_specs
!
         read(input%unit, '(a100)') line
         line = remove_preceding_blanks(line)
!
         if (line(1:18) == 'transition moments' ) then
!
            solver%moments = .true.
!
         elseif (line(1:14) == 'polarizability' ) then
!
            solver%polarizability = .true.
!
         endif
!
      enddo
!
      if ( .not. solver%moments .and. .not. solver%polarizability) then
!
            call output%error_msg('No RHS specified for linear response.')
!
      endif
! 
!     Read excited state section to get n excited states 
!     OBS: only necessary when 'moments' are true
!
      if (solver%moments) then
!
         if (.not. requested_section('cc excited state')) then
!
            call output%error_msg('number of excitations must be specified.')
!
         endif
!
         call move_to_section('cc excited state', n_specs)
!
         do i = 1, n_specs
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:15) == 'singlet states:' ) then
!
               read(line(16:100), *) solver%n_freq
!
               solver%dim_rhs = solver%n_freq
!
            endif
!
         enddo
!
      endif
!
   end subroutine read_settings_davidson_cc_response
!
!
end module davidson_cc_response_class
