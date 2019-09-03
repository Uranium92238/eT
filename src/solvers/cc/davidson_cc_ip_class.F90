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
module davidson_cc_ip_class
!
!!
!!    Davidson coupled cluster ionized state solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    
!
   use davidson_cc_es_class
!
   implicit none
!
   type, extends(davidson_cc_es) :: davidson_cc_ip
!
   contains
!
      procedure :: set_start_vectors      => set_start_vectors_davidson_cc_ip
      procedure :: set_projection_vector  => set_projection_vector_davidson_cc_ip
!
   end type davidson_cc_ip
!
!
   interface davidson_cc_ip
!
      procedure :: new_davidson_cc_ip
!
   end interface davidson_cc_ip
!
!
contains
!
   function new_davidson_cc_ip(transformation, wf) result(solver)
!!
!!    New Davidson CC IP
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(davidson_cc_ip) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' ionized state')
      call solver%timer%turn_on()
!
      solver%tag = 'Davidson IP'
!
      solver%name_ = 'Davidson coupled cluster ionized state solver'
      solver%author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
      solver%description1 = 'A Davidson solver that calculates ionization energies and the &
                            &corresponding right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                            &problem is solved in a reduced space, the dimension of which is expanded &
                            &until the convergence criteria are met.'
      solver%description2 = 'A complete description of the Davidson algorithm can be found in &
                            &E. R. Davidson, J. Comput. Phys. 17, 87 (1975).'
!
      call solver%print_banner()
!
!     Set defaults, then read possible non-defaults
!
      solver%n_singlet_states     = 0
      solver%max_iterations       = 100
      solver%eigenvalue_threshold = 1.0d-6
      solver%residual_threshold   = 1.0d-6
      solver%transformation       = 'right'
      solver%restart              = .false.
      solver%max_dim_red          = 100 
      solver%transformation       = trim(transformation)
!
      call solver%read_settings()
      call solver%print_settings()
!
      call solver%initialize_energies()
      solver%energies = zero
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      write(output%unit, '(/t3,a,a,a)') 'Solving for the ', trim(solver%transformation), ' eigenvectors.'
      flush(output%unit)
!
      wf%n_excited_states = solver%n_singlet_states
!
   end function new_davidson_cc_ip
!
!
   subroutine set_start_vectors_davidson_cc_ip(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cc_ip) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: c_i
!
      integer, dimension(:), allocatable :: start_indices
!
      integer :: trial, n_solutions_on_file
!
      if (.not. wf%bath_orbital) call output%error_msg('Calculation of IPs requires bath orbitals.')
!
      if (solver%restart) then 
!
!        Read the solutions from file & set as initial trial vectors 
!
         call wf%is_restart_safe('excited state')
!
         call solver%determine_restart_transformation(wf) ! Read right or left?
         n_solutions_on_file = wf%get_n_excited_states_on_file(solver%restart_transformation)
!
         call output%printf('Requested restart - there are (i0) (a0) eigenvectors on file.', &
                              ints=[n_solutions_on_file], chars=[solver%restart_transformation])
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         do trial = 1, n_solutions_on_file
!
            call wf%read_excited_state(c_i, trial, solver%restart_transformation)
            call davidson%write_trial(c_i)
!
         enddo 
!
         call mem%dealloc(c_i, wf%n_es_amplitudes)
!
      else
!
         n_solutions_on_file = 0
!
      endif 
!
      if (n_solutions_on_file .lt. solver%n_singlet_states) then ! Koopman for the rest

         call mem%alloc(start_indices, solver%n_singlet_states)
!
         call wf%set_ip_start_indices(start_indices, solver%n_singlet_states)
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         do trial = n_solutions_on_file + 1, solver%n_singlet_states
!
            call zero_array(c_i, wf%n_es_amplitudes)
            c_i(start_indices(trial)) = one
            call davidson%write_trial(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes)
         call mem%dealloc(start_indices, solver%n_singlet_states)
!
      endif
!
      call davidson%orthonormalize_trial_vecs()
!
   end subroutine set_start_vectors_davidson_cc_ip
!
!
   subroutine set_projection_vector_davidson_cc_ip(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cc_ip) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: projector
!
      if (.false.) then
!
         call output%printf(solver%name_) ! Stupid hack to disable warning
!
      endif
!
      call mem%alloc(projector, wf%n_es_amplitudes)
      davidson%do_projection = .true.
!
      call wf%get_ip_projector(projector)
!
      call davidson%set_projector(projector)
      call mem%dealloc(projector, wf%n_es_amplitudes)
!
   end subroutine set_projection_vector_davidson_cc_ip
!
!
end module davidson_cc_ip_class
