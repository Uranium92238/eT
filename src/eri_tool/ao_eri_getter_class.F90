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
module ao_eri_getter_class
!!
!! AO ERI getter class
!! Written by Sarai D. Folkestad, 2022
!!
!! Tool to extract ERIs from the AO-tool
!!
!!
   use parameters
!
   use global_out,                     only : output
   use memory_manager_class,           only : mem
   use ao_tool_class,                  only : ao_tool
!
   implicit none
!
   type :: ao_eri_getter
!
      type(ao_tool), pointer :: ao
!
   contains
!
      procedure, public :: get_eri => get_eri_ao_eri_getter
!
   end type ao_eri_getter
!
   interface ao_eri_getter
!
      procedure :: new_ao_eri_getter
!
   end interface ao_eri_getter
!
contains
!
!
   function new_ao_eri_getter(ao) result(this)
!!
!!    New
!!    Written by Sarai D. Folkestad, Aug 2021
!!
      implicit none
!
      type(ao_eri_getter) :: this
      class(ao_tool), intent(in), target :: ao
!
      this%ao => ao
!
   end function new_ao_eri_getter
!
!
   subroutine get_eri_ao_eri_getter(this, g, A, B, C, D, precision_, skip)
!!
!!    Get ERI
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2020
!!
!!    Wrapper for ao tool routine
!!
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(ao_eri_getter) :: this
!
      integer :: A, B, C, D
!
      real(dp), dimension(*), intent(out) :: g
!
      real(dp), optional, intent(in) :: precision_
      integer, optional, intent(out) :: skip
!
      call this%ao%get_eri(g, A, B, C, D, precision_, skip)
!
   end subroutine get_eri_ao_eri_getter
!
end module ao_eri_getter_class
