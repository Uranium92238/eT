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
module ri_basis_class
!
!!
!!    RI basis class
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!
   use kinds
   use iso_c_binding
!
   use shell_class, only : shell
!
   implicit none
!
!
   type :: ri_basis
!
      type(shell), dimension(:), allocatable, private :: shells
!
      integer :: n_shells
      integer :: n_ao
      integer :: max_dim
!
      character(len=200) :: basis_set
!
   contains
!
      procedure, public :: initialize_shells
      procedure, public :: cleanup_shells
      procedure, public :: get_shell_length
      procedure, public :: get_shell_first
      procedure, public :: get_shell_last
      procedure, public :: get_shell
!
   end type ri_basis
!
   include "../libint/ri_basis_cdef.F90"
   include "../libint/libint_initialization_cdef.F90"
!
   interface ri_basis
!
      procedure :: new_ri_basis
!
   end interface ri_basis
!
!
contains
!
!
   function new_ri_basis(basis_set, eri_precision) result(this)
!!
!!    New ri basis
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!
      use iso_c_binding, only: c_char, c_null_char, c_double, c_int
!
      implicit none
!
      character(len=200), intent(in) :: basis_set
      real(dp), intent(in) :: eri_precision
!
      type(ri_basis) :: this
!
      character(len=200, kind=c_char) :: basis_set_c
      real(c_double) :: eri_precision_c
      integer(c_int), dimension(1) :: n_shells_c, n_ao_c
!
      this%basis_set = basis_set
!
      basis_set_c = trim(basis_set)//c_null_char
      eri_precision_c = real(eri_precision, kind=c_double)
!
      call prepare_for_ri_c(eri_precision_c, basis_set_c)
!
      call get_n_ri_ao_c(n_ao_c(1))
!
      call get_n_ri_shells_c(n_shells_c(1))
!
      this%n_ao = int(n_ao_c(1))
      this%n_shells = int(n_shells_c(1))
!
   end function new_ri_basis
!
!
   subroutine initialize_shells(this)
!!
!!    Initialize ri shells
!!    Written by Sarai D. Folkestad, Aug 2021
!!
!
      use iso_c_binding, only: c_int
!
      implicit none
!
      class(ri_basis) :: this

      integer(c_int), dimension(1) :: n_ao_c, shell_
!
      integer :: offset, i
!
      allocate(this%shells(this%n_shells))
!
      offset = 0
      this%max_dim = 0
!
      do i = 1, this%n_shells
!
      shell_(1) = int(i, kind=c_int)
      call get_ri_shell_size_c(shell_, n_ao_c)
      this%shells(i) = shell(first=offset + 1,       &
                             length=int(n_ao_c(1)),  &
                             number_=i,              &
                             cartesian=.false.) ! Solid harmonics are default for l>1 in Libint
!
         offset = offset + int(n_ao_c(1))
!
         if (int(n_ao_c(1)) > this%max_dim) this%max_dim = int(n_ao_c(1))
!
      enddo
!
   end subroutine initialize_shells
!
!
   subroutine cleanup_shells(this)
!!
!!    Cleanup shells
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      class(ri_basis) :: this
!
      deallocate(this%shells)
!
   end subroutine cleanup_shells
!
   pure function get_shell_length(this, I) result(length)
!!
!!    Get shell length
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ri_basis), intent(in) :: this
      integer, intent(in) :: I
      integer :: length
!
      length = this%shells(I)%length
!
   end function get_shell_length
!
   pure function get_shell_first(this, I) result(first)
!!
!!    Get shell first
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ri_basis), intent(in) :: this
      integer, intent(in) :: I
      integer :: first
!
      first = this%shells(I)%first
!
   end function get_shell_first
!
   pure function get_shell_last(this, I) result(last)
!!
!!    Get shell last
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ri_basis), intent(in) :: this
      integer, intent(in) :: I
      integer :: last
!
      last = this%shells(I)%get_last()
!
   end function get_shell_last
!
   pure function get_shell(this, I) result(shell_I)
!!
!!    Get shell
!!    Written by Sarai D. Folkestad
!!
      implicit none
!
      class(ri_basis), intent(in) :: this
      integer, intent(in) :: I
      type(shell) :: shell_I
!
      shell_I = this%shells(I)
!
   end function get_shell
!
end module ri_basis_class
