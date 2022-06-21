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
module abstract_G_adder_class
!!
!! Abstract G adder class
!! Written by Sarai D. Folkestad, 2021
!!
!! Adds contributions to the two-electron contribution to
!! the Fock matrix (G)
!!
   use kinds
   use ao_tool_class, only: ao_tool
   use memory_manager_class, only: mem
!
   implicit none
!
   type, abstract :: abstract_G_adder
!
   contains
!
      procedure (add_abstract), deferred, nopass ::  add
!
   end type abstract_G_adder
!
!
   abstract interface
!
!
      subroutine add_abstract(ao, D, G, eri, s1, s2, s3, s4, thread, n_threads, degeneracy)
!!
      use parameters
!
      import abstract_G_adder
      import ao_tool
!
      implicit none
!
      class(ao_tool), intent(in)  :: ao
      integer,                                    intent(in)    :: n_threads, thread
      real(dp), dimension(ao%n, ao%n, n_threads), intent(inout) :: G
      real(dp), dimension(ao%n, ao%n),            intent(in)    :: D
      real(dp), dimension(ao%max_sh_size**4),     intent(inout) :: eri
      real(dp), intent(in) :: degeneracy
!
      integer, intent(in) :: s1, s2, s3, s4
!
      end subroutine add_abstract
!
!
   end interface
!
end module abstract_G_adder_class
