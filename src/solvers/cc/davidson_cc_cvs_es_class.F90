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
module davidson_cvs_cc_es_class
!
!!
!!    Davidson CVS coupled cluster excited state solver class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    
!
   use davidson_cc_es_class
!
   implicit none
!
   type, extends(davidson_cc_es) :: davidson_cvs_cc_es
!
      integer :: n_core_Mos
!
      integer, dimension(:), allocatable :: core_MOs
!
   contains
!
      procedure :: read_settings          => read_settings_davidson_cvs_cc_es
      procedure :: read_cvs_settings      => read_cvs_settings_davidson_cvs_cc_es
!
      procedure :: set_start_vectors      => set_start_vectors_davidson_cvs_cc_es
      procedure :: set_projection_vector  => set_projection_vector_davidson_cvs_cc_es
!
      procedure :: initialize_core_MOs    => initialize_core_MOs_davidson_cvs_cc_es
!
      procedure :: destruct_core_MOs      => destruct_core_MOs_davidson_cvs_cc_es
!
   end type davidson_cvs_cc_es
!
!
   interface davidson_cvs_cc_es
!
      procedure :: new_davidson_cvs_cc_es
!
   end interface davidson_cvs_cc_es
!
!
contains
!
   function new_davidson_cvs_cc_es(transformation, wf) result(solver)
!!
!!    New Davidson CVS CC ES 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(davidson_cvs_cc_es) :: solver
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: transformation
!
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' core excited state')
      call solver%timer%turn_on()
!
      solver%tag = 'Davidson CVS'
!
      solver%name_ = 'Davidson coupled cluster core excited state solver'
      solver%author = 'E. F. Kjønstad, S. D. Folkestad, 2018'
      solver%description1 = 'A Davidson CVS solver that calculates core excitation energies and the &
                            &corresponding right eigenvectors of the Jacobian matrix, A. The eigenvalue &
                            &problem is solved in a reduced space, the dimension of which is expanded &
                            &until the convergence criteria are met. In addition the CVS aproximation is &
                            &used to obtain the core excitations'
      solver%description2 = 'A complete description of the Davidson algorithm can be found in &
                            &E. R. Davidson, J. Comput. Phys. 17, 87 (1975). &
                            &A description of the CVS approximation can be found in &
                            &S. Coriani & H. Koch, J. Chem. Phys. 143, 181103 (2015).'
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
   end function new_davidson_cvs_cc_es
!
!
   subroutine read_settings_davidson_cvs_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cvs_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_davidson_settings()
      call solver%read_cvs_settings()
!
   end subroutine read_settings_davidson_cvs_cc_es
!
!
   subroutine read_cvs_settings_davidson_cvs_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cvs_cc_es) :: solver 
!
      if (input%requested_keyword_in_section('core excitation', 'solver cc es')) then 
!  
!        Determine the number of core MOs 
!
         solver%n_core_MOs = input%get_n_elements_for_keyword_in_section('core excitation', 'solver cc es')
!
!        Then read the vector of core MOs for CVS
!
         call solver%initialize_core_MOs()
!
         call input%get_array_for_keyword_in_section('core excitation', 'solver cc es', solver%n_core_MOs, solver%core_MOs)
!
      else
!
         call output%error_msg('found no specified core MOs in input for cvs calculation')
!
      endif 
!
   end subroutine read_cvs_settings_davidson_cvs_cc_es
!
!
   subroutine set_start_vectors_davidson_cvs_cc_es(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cvs_cc_es) :: solver
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
         call wf%set_cvs_start_indices(solver%n_core_MOs, solver%core_MOs, solver%n_singlet_states, start_indices)
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         do trial = n_solutions_on_file + 1, solver%n_singlet_states
!
            c_i = zero
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
   end subroutine set_start_vectors_davidson_cvs_cc_es
!
!
   subroutine set_projection_vector_davidson_cvs_cc_es(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cvs_cc_es) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: projector
!
      call mem%alloc(projector, wf%n_es_amplitudes)
      davidson%do_projection = .true.
!
      call wf%get_cvs_projector(projector, solver%n_core_MOs, solver%core_MOs)
!
      call davidson%set_projector(projector)
      call mem%dealloc(projector, wf%n_es_amplitudes)
!
   end subroutine set_projection_vector_davidson_cvs_cc_es
!
!
   subroutine initialize_core_MOs_davidson_cvs_cc_es(solver)
!!
!!    Initialize core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(davidson_cvs_cc_es) :: solver
!
      if (.not. allocated(solver%core_MOs)) call mem%alloc(solver%core_MOs, solver%n_core_MOs)
!
   end subroutine initialize_core_MOs_davidson_cvs_cc_es
!
!
   subroutine destruct_core_MOs_davidson_cvs_cc_es(solver)
!!
!!    Destruct core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(davidson_cvs_cc_es) :: solver
!
      if (allocated(solver%core_MOs)) call mem%dealloc(solver%core_MOs, solver%n_core_MOs)
!
   end subroutine destruct_core_MOs_davidson_cvs_cc_es
!
!
end module davidson_cvs_cc_es_class
