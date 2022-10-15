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
!!    A_mu,nu = < mu |[H-bar, tau_nu] | HF >,    H-bar = e-T H eT.
!!
!
   use parameters
!
   use global_out, only : output
!
   use memory_manager_class, only : mem
   use timings_class, only : timings
!
   use ccs_class, only : ccs
!
   use start_vector_tool_class, only : start_vector_tool
!
   use precondition_tool_class, only: precondition_tool
!
   use abstract_projection_tool_class, only: abstract_projection_tool
!
   use convergence_tool_class, only: convergence_tool
!
   use abstract_solver_class, only: abstract_solver
!
   implicit none
!
   type, extends(abstract_solver), abstract :: abstract_cc_es
!
      character(len=100) :: name_
      character(len=100) :: tag
!
      character(len=500) :: description1
      character(len=500) :: description2
!
      integer :: max_iterations, n_singlet_states
!
      logical :: restart
!
      character(len=40) :: transformation
!
      real(dp), dimension(:), allocatable :: energies
!
      class(convergence_tool), allocatable         :: convergence_checker
      class(start_vector_tool), allocatable        :: start_vectors
      class(abstract_projection_tool), allocatable :: projector
      class(precondition_tool), allocatable        :: preconditioner
!
      class(ccs), pointer :: wf
!
   contains
!
      procedure :: print_banner                     => print_banner_abstract_cc_es
      procedure :: print_es_settings                => print_es_settings_abstract_cc_es
!
      procedure :: cleanup                          => cleanup_abstract_cc_es
!
      procedure :: prepare_wf_for_excited_state     => prepare_wf_for_excited_state_abstract_cc_es
!
      procedure :: initialize_energies              => initialize_energies_abstract_cc_es
      procedure :: destruct_energies                => destruct_energies_abstract_cc_es
!
   end type abstract_cc_es
!
!
contains
!
!
   subroutine initialize_energies_abstract_cc_es(this)
!!
!!    Initialize energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Initialize excitation energies
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      if (.not. allocated(this%energies)) &
         call mem%alloc(this%energies, this%n_singlet_states)
!
   end subroutine initialize_energies_abstract_cc_es
!
!
   subroutine destruct_energies_abstract_cc_es(this)
!!
!!    Destruct energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Destruct excitation energies
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      if (allocated(this%energies)) &
         call mem%dealloc(this%energies, this%n_singlet_states)
!
   end subroutine destruct_energies_abstract_cc_es
!
!
   subroutine print_banner_abstract_cc_es(this)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      call output%printf('m', ' - ' // trim(this%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(this%name_)) + 6, '-')
!
      call output%printf('n', this%description1, ffs='(/t3,a)', fs='(t3,a)')
      call output%printf('n', this%description2, ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_abstract_cc_es
!
!
   subroutine print_es_settings_abstract_cc_es(this)
!!
!!    Print settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(abstract_cc_es) :: this
!
      call output%printf('m', '- Settings for coupled cluster excited state &
                         &solver (' //trim(this%tag) // '):', fs='(/t3,a)')
      call output%printf('m', 'Excitation vectors:  (a0)', &
                         chars=[trim(this%transformation)], fs='(t6,a)')
!
      call this%convergence_checker%print_settings()
!
      call output%printf('m', 'Number of singlet states:     (i11)', &
                         ints=[this%n_singlet_states], fs='(/t6,a)')
      call output%printf('m', 'Max number of iterations:     (i11)', &
                         ints=[this%max_iterations], fs='(t6,a)')
!
   end subroutine print_es_settings_abstract_cc_es
!
!
   subroutine cleanup_abstract_cc_es(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(abstract_cc_es), intent(inout) :: this
!
      call this%destruct_energies()
!
      call this%total_timer%turn_off()
!
      call output%printf('m', '- Finished solving the ' //  &
                        trim(convert_to_uppercase(this%wf%name_)) // ' excited state &
                        &equations ('// trim(this%transformation) //')', fs='(/t3,a)')
!
      call output%printf('m', 'Total wall time (sec): (f20.5)', &
                         reals=[this%total_timer%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'Total cpu time (sec):  (f20.5)', &
                         reals=[this%total_timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine cleanup_abstract_cc_es
!
!
   subroutine prepare_wf_for_excited_state_abstract_cc_es(this)
!!
!!    Prepare wf for excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2019
!!
      implicit none
!
      class(abstract_cc_es), intent(in) :: this
!
      call this%wf%prepare_for_Jacobians(this%transformation)
!
   end subroutine prepare_wf_for_excited_state_abstract_cc_es
!
!
end module abstract_cc_es_class
