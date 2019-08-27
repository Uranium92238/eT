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
module diis_cvs_cc_es_class
!
!!
!!    DIIS CVS coupled cluster excited state solver class module
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Based on the diis solver
!! 
!
   use diis_cc_es_class, only: diis_cc_es
   use kinds
   use file_class
   use ccs_class
!
   implicit none
!
   type, extends(diis_cc_es) :: diis_cvs_cc_es
!
      integer :: n_core_Mos
!
      integer, dimension(:), allocatable :: core_MOs
!
   contains
!
      procedure :: read_settings          => read_settings_diis_cvs_cc_es
      procedure :: read_cvs_settings      => read_cvs_settings_diis_cvs_cc_es
!
      procedure :: set_start_vectors      => set_start_vectors_diis_cvs_cc_es
      procedure :: set_projection_vector  => set_projection_vector_diis_cvs_cc_es
!
      procedure :: initialize_core_MOs    => initialize_core_MOs_diis_cvs_cc_es
!
      procedure :: destruct_core_MOs      => destruct_core_MOs_diis_cvs_cc_es
!
   end type diis_cvs_cc_es
!
!
   interface diis_cvs_cc_es
!
      procedure :: new_diis_cvs_cc_es
!
   end interface diis_cvs_cc_es
!
!
contains
!
!
   function new_diis_cvs_cc_es(transformation, wf) result(solver)
!!
!!    New DIIS CVS CC ES 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(diis_cvs_cc_es) :: solver
      class(ccs), intent(in) :: wf
!
      character(len=*), intent(in) :: transformation
!
      solver%timer = new_timer(trim(convert_to_uppercase(wf%name_)) // ' core excited state')
      call solver%timer%turn_on()
!
      solver%tag = 'DIIS CVS'
!
      solver%name_         = 'DIIS coupled cluster core excited state solver'
      solver%author        = 'E. F. Kjønstad, S. D. Folkestad, 2018'
      solver%description1  = 'A DIIS CVS solver that calculates core excitation energies and the &
                            &corresponding right eigenvectors of the Jacobian matrix, A.  The eigenvalue &
                            &problem is solved by DIIS extrapolation of residuals for each &
                            &eigenvector until the convergence criteria are met.'
      solver%description2  = 'More on the DIIS algorithm can be found in &
                           &P. Pulay, Chemical Physics Letters, 73(2), 393-398 (1980).'
!
      call solver%print_banner()
!
!     Set defaults
!
      solver%n_singlet_states     = 0
      solver%max_iterations       = 100
      solver%eigenvalue_threshold = 1.0d-6
      solver%residual_threshold   = 1.0d-6
      solver%transformation       = 'right'
      solver%diis_dimension       = 20
      solver%restart              = .false.
      solver%transformation       = trim(transformation)
!
      solver%do_projection = .false.
!
      call solver%read_settings()
      call solver%print_settings()
!
      if (solver%n_singlet_states == 0) call output%error_msg('number of excitations must be specified.')
!
      call mem%alloc(solver%energies, solver%n_singlet_states)
      solver%energies = zero
!
   end function new_diis_cvs_cc_es
!
!
   subroutine read_settings_diis_cvs_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cvs_cc_es) :: solver 
!
      call solver%read_es_settings()
      call solver%read_diis_settings()
      call solver%read_cvs_settings()
!
   end subroutine read_settings_diis_cvs_cc_es
!
!
   subroutine read_cvs_settings_diis_cvs_cc_es(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(diis_cvs_cc_es) :: solver 
!
      if (input%requested_keyword_in_section('core excitation', 'solver cc es')) then 
!  
!        Determine the number of core MOs 
!
         solver%n_core_MOs = input%get_n_elements_for_keyword_in_section('core excitation', 'solver cc es')
!
!        Then read the vector of core MOs
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
   end subroutine read_cvs_settings_diis_cvs_cc_es
!
!
   subroutine initialize_core_MOs_diis_cvs_cc_es(solver)
!!
!!    Initialize core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(diis_cvs_cc_es) :: solver
!
      if (.not. allocated(solver%core_MOs)) call mem%alloc(solver%core_MOs, solver%n_core_MOs)
!
   end subroutine initialize_core_MOs_diis_cvs_cc_es
!
!
   subroutine destruct_core_MOs_diis_cvs_cc_es(solver)
!!
!!    Destruct core MOs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
      implicit none
!
      class(diis_cvs_cc_es) :: solver
!
      if (allocated(solver%core_MOs)) call mem%dealloc(solver%core_MOs, solver%n_core_MOs)
!
   end subroutine destruct_core_MOs_diis_cvs_cc_es
!
!
   subroutine set_projection_vector_diis_cvs_cc_es(solver, wf)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(diis_cvs_cc_es) :: solver
!
      class(ccs) :: wf
!
      call mem%alloc(solver%projector, wf%n_es_amplitudes)
      solver%do_projection = .true.
!
      call wf%get_cvs_projector(solver%projector, solver%n_core_MOs, solver%core_MOs)
!
   end subroutine set_projection_vector_diis_cvs_cc_es
!
!
   subroutine set_start_vectors_diis_cvs_cc_es(solver, wf, R, orbital_differences)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(diis_cvs_cc_es), intent(in) :: solver 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes, solver%n_singlet_states), intent(inout) :: R 
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: orbital_differences 
!
      integer, dimension(:), allocatable :: cvs_start_indices
!
      integer :: state
!
      if (.false.) write (output%unit) orbital_differences(1) ! Hack to avoid warning
!
      call mem%alloc(cvs_start_indices, solver%n_singlet_states)
      call wf%set_cvs_start_indices(solver%n_core_MOs, solver%core_MOs, solver%n_singlet_states, cvs_start_indices)
!
      call zero_array(R,(wf%n_es_amplitudes)*(solver%n_singlet_states))
!
      do state = 1, solver%n_singlet_states
!
         R(cvs_start_indices(state), state) = one
!
      enddo 
!
      call mem%dealloc(cvs_start_indices, solver%n_singlet_states)
!
   end subroutine set_start_vectors_diis_cvs_cc_es
!
end module diis_cvs_cc_es_class
