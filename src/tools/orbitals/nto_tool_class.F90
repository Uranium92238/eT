!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module nto_tool_class
!
!!
!!    Natural transition orbital tool class
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Handles construction and diagonalization of M and N
!!
!!    M_ij +=  sum_a X_ai X_aj
!!    N_ab +=  sum_i X_ai X_bi
!!
!
   use parameters
   use global_out,           only: output
   use memory_manager_class, only: mem
   use array_utilities,      only: zero_array, diagonalize_symmetric
   use stream_file_class,    only: stream_file
!
   implicit none
!
   type :: nto_tool
!
      integer :: n_ao, n_o, n_v, X_length, n_states
!
      real(dp), dimension(:),   allocatable :: eigenvalues_o, eigenvalues_v
      real(dp), dimension(:,:), allocatable :: M, N
!
      type(stream_file), allocatable :: M_file, N_file
!
   contains
!
         procedure, public :: initialize &
                           => initialize_nto_tool
!
         procedure, public :: read_M_and_N &
                           => read_M_and_N_nto_tool
!
         procedure, public :: write_M_and_N &
                           => write_M_and_N_nto_tool
!
         procedure, public :: add_excited_state &
                           => add_excited_state_nto_tool
!
         procedure, public :: add_singles_to_M_and_N &
                           => add_singles_to_M_and_N_nto_tool
!
         procedure, public :: add_contributions_to_M_and_N &
                           => add_contributions_to_M_and_N_nto_tool
!
         procedure, public :: diagonalize_M_and_N &
                           => diagonalize_M_and_N_nto_tool
!
         procedure, public :: transform_orbitals &
                           => transform_orbitals_nto_tool
!
         procedure, public :: get_n_active_orbitals &
                           => get_n_active_orbitals_nto_tool
!
         procedure, public :: cleanup &
                           => cleanup_nto_tool
!
         procedure, public :: prepare &
                           => prepare_nto_tool
!
         procedure, private, nopass :: print_n_active_orbitals
!
   end type  nto_tool
!
   interface  nto_tool
!
      procedure :: new_nto_tool
!
   end interface  nto_tool
!
contains
!
!
   function new_nto_tool(n_o, n_v, n_ao, X_length) result(this)
!!
!!    New nto tool
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      type(nto_tool) :: this
!
      integer, intent(in) :: n_o, n_v, n_ao, X_length
!
      this%n_ao = n_ao
      this%n_o = n_o
      this%n_v = n_v
      this%X_length = X_length
!
      call this%prepare('nto')
!
   end function new_nto_tool
!
!
   subroutine prepare_nto_tool(this, tag)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      class(nto_tool) :: this
!
      character(len=*) :: tag
!
      this%n_states = 0
!
      this%M_file = stream_file(tag // '_M_transformation')
      this%N_file = stream_file(tag // '_N_transformation')
!
      if (this%n_ao .le. 0) then
         call output%error_msg('n_ao less or equal 0 in (a0)-tool', chars=[tag])
      else if (this%n_o .le. 0) then
         call output%error_msg('n_o less or equal 0 in (a0)-tool', chars=[tag])
      else if (this%n_v .le. 0) then
         call output%error_msg('n_v less or equal 0 in (a0)-tool', chars=[tag])
      else if (this%X_length .le. 0) then
         call output%error_msg('Length of vector less or equal 0 in (a0)-tool', chars=[tag])
      end if
!
   end subroutine prepare_nto_tool
!
!
   subroutine initialize_nto_tool(this)
!!
!!    Initialize
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      call mem%alloc(this%eigenvalues_o, this%n_o)
      call mem%alloc(this%eigenvalues_v, this%n_v)
!
      call mem%alloc(this%M, this%n_o, this%n_o)
      call mem%alloc(this%N, this%n_v, this%n_v)
!
      call zero_array(this%M, this%n_o**2)
      call zero_array(this%N, this%n_v**2)
!
      this%n_states = 0
!
   end subroutine initialize_nto_tool
!
!
   subroutine read_M_and_N_nto_tool(this, files_found)
!!
!!    Read M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Read the M and N matrices from file for restart
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      logical, intent(out) :: files_found
!
      integer :: n_states
!
      files_found = this%M_file%exists() .and. this%N_file%exists()
!
      if (files_found) then
!
         call this%M_file%open_('read', 'rewind')
!
         call this%M_file%read_(this%n_states)
         call this%M_file%read_(this%M, this%n_o**2)
!
         call this%M_file%close_()
!
         call this%N_file%open_('read', 'rewind')
!
         call this%N_file%read_(n_states)
!
         if (n_states .ne. this%n_states) then
            call output%error_msg('M and N do not correspond to the same ntos. M &
                                  &consists of (i0) while N consists of (i0) states.', &
                                  ints=[this%n_states, n_states])
         end if
!
         call this%N_file%read_(this%N, this%n_v**2)
!
         call this%N_file%close_()
!
      end if
!
   end subroutine read_M_and_N_nto_tool
!
!
   subroutine write_M_and_N_nto_tool(this)
!!
!!    Write M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    For restart
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      call this%M_file%open_('write', 'rewind')
!
      call this%M_file%write_(this%n_states)
      call this%M_file%write_(this%M, this%n_o**2)
!
      call this%M_file%close_()
!
      call this%N_file%open_('write', 'rewind')
!
      call this%N_file%write_(this%n_states)
      call this%N_file%write_(this%N, this%n_v**2)
!
      call this%N_file%close_()
!
   end subroutine write_M_and_N_nto_tool
!
!
   subroutine add_excited_state_nto_tool(this, X)
!!
!!    Add excited state
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Add the contribution of an excited state X to M and N
!!
!
      use array_utilities, only: copy_and_scale
!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      real(dp), dimension(this%X_length), intent(in) :: X
!
      real(dp) :: norm_, ddot
!
      real(dp), dimension(:), allocatable :: X_copy
!
      this%n_states = this%n_states + 1
!
      call mem%alloc(X_copy, this%X_length)
!
      norm_ = sqrt(ddot(this%X_length, X, 1, X, 1))
      call copy_and_scale(one/norm_, X, X_copy, this%X_length)
!
      call this%add_contributions_to_M_and_N(X_copy)
!
      call mem%dealloc(X_copy, this%X_length)
!
   end subroutine add_excited_state_nto_tool
!
!
   subroutine add_contributions_to_M_and_N_nto_tool(this, X)
!!
!!    Add contributions to M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      real(dp), dimension(this%X_length), intent(in) :: X
!
      call this%add_singles_to_M_and_N(X)
!
   end subroutine add_contributions_to_M_and_N_nto_tool
!
!
   subroutine add_singles_to_M_and_N_nto_tool(this, X)
!!
!!    Add singles to M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Add the contribution of an singly excited state X to M and N
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_v, this%n_o), intent(in) :: X
!
      call dgemm('T', 'N',    &
                  this%n_o,   &
                  this%n_o,   &
                  this%n_v,   &
                  one,        &
                  X,          & ! X_ai
                  this%n_v,   &
                  X,          & ! X_aj
                  this%n_v,   &
                  one,        &
                  this%M,     & ! M_ij
                  this%n_o)
!
      call dgemm('N', 'T',    &
                  this%n_v,   &
                  this%n_v,   &
                  this%n_o,   &
                  one,        &
                  X,          & ! X_ai
                  this%n_v,   &
                  X,          & ! X_bi
                  this%n_v,   &
                  one,        &
                  this%N,     & ! N_ab
                  this%n_v)
!
   end subroutine add_singles_to_M_and_N_nto_tool
!
!
   subroutine diagonalize_M_and_N_nto_tool(this)
!!
!!    Add singles to M and N
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      real(dp) :: normalization_factor
!
!     Scale by -1 to get eigenvalues in the correct order.
!     Scale by 1/n_states to get eigenvalues that add up to one.
!
      normalization_factor = -one/real(this%n_states, kind=dp)
!
      call dscal(this%n_o**2, normalization_factor, this%M, 1)
      call dscal(this%n_v**2, normalization_factor, this%N, 1)
!
      call diagonalize_symmetric(this%M, this%n_o, this%eigenvalues_o)
      call diagonalize_symmetric(this%N, this%n_v, this%eigenvalues_v)
!
      call dscal(this%n_o, -one, this%eigenvalues_o, 1)
      call dscal(this%n_v, -one, this%eigenvalues_v, 1)
!
   end subroutine diagonalize_M_and_N_nto_tool
!
!
   subroutine get_n_active_orbitals_nto_tool(this, n_nto_o, n_nto_v, threshold)
!!
!!    Get number of active orbitals
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
!!    Number of NTOs/CNTOs selected according to
!!    1 - sum_i lambda_i < threshold, lambda is eigenvalue of M or N
!!
      use math_utilities, only: get_n_values_greater_sum
!
      implicit none
!
      class(nto_tool), intent(in) :: this
!
      integer, intent(out) :: n_nto_o, n_nto_v
      real(dp), intent(in) :: threshold
!
      n_nto_o = get_n_values_greater_sum(this%n_o, this%eigenvalues_o, one-threshold)
      call this%print_n_active_orbitals(n_nto_o, this%eigenvalues_o, 'occupied')
!
      n_nto_v = get_n_values_greater_sum(this%n_v, this%eigenvalues_v, one-threshold)
      call this%print_n_active_orbitals(n_nto_v, this%eigenvalues_v, 'virtual')
!
   end subroutine get_n_active_orbitals_nto_tool
!
!
   subroutine print_n_active_orbitals(n_nto, eigenvalues, tag)
!!
!!    Print n active orbitals
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      integer, intent(in) :: n_nto
!
      real(dp), dimension(n_nto), intent(in) :: eigenvalues
!
      character(len=*), intent(in) :: tag
!
      integer :: p
!
      call output%printf('n', 'Number of significant (a0) orbitals: (i0)', &
                         ints = [n_nto], chars=[trim(tag)], ffs='(/t6,a)')
!
      call output%printf('n', 'Significant eigenvalues:', ffs='(/t8,a)')

      do p = 1, n_nto
!
         call output%printf('n', '(e16.8)', reals = [eigenvalues(p)], ffs='(t8,a)')
!
      enddo
!
   end subroutine print_n_active_orbitals
!
!
   subroutine transform_orbitals_nto_tool(this, C, C_ntos)
!!
!!    Transform orbitals
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      real(dp), dimension(this%n_ao, this%n_o+this%n_v), intent(in)  :: C
      real(dp), dimension(this%n_ao, this%n_o+this%n_v), intent(out) :: C_ntos
!
      call this%diagonalize_M_and_N()
!
      call zero_array(C_ntos, this%n_ao * (this%n_o+this%n_v))
!
      call dgemm('N', 'N',    &
                  this%n_ao,  &
                  this%n_o,   &
                  this%n_o,   &
                  one,        &
                  C,          &
                  this%n_ao,  &
                  this%M,     &
                  this%n_o,   &
                  zero,       &
                  C_ntos,     &
                  this%n_ao)
!
      call dgemm('N', 'N',                &
                  this%n_ao,              &
                  this%n_v,               &
                  this%n_v,               &
                  one,                    &
                  C(:, this%n_o+1),       &
                  this%n_ao,              &
                  this%N,                 &
                  this%n_v,               &
                  one,                    &
                  C_ntos(:, this%n_o+1),  &
                  this%n_ao)
!
   end subroutine transform_orbitals_nto_tool
!
!
   subroutine cleanup_nto_tool(this)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Alexander C. Paul, May 2021
!!
      implicit none
!
      class(nto_tool), intent(inout) :: this
!
      call mem%dealloc(this%eigenvalues_o, this%n_o)
      call mem%dealloc(this%eigenvalues_v, this%n_v)
!
      call mem%dealloc(this%M, this%n_o, this%n_o)
      call mem%dealloc(this%N, this%n_v, this%n_v)
!
   end subroutine cleanup_nto_tool
!
!
end module nto_tool_class

