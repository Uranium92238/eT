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
module J_adder_class
!!
!! J adder class
!! Written by Sarai D. Folkestad, 2021
!!
   use parameters
   use ao_tool_class, only: ao_tool
   use abstract_G_adder_class, only: abstract_G_adder
!
   implicit none
!
   type, extends(abstract_G_adder) :: J_adder
!
   contains
!
      procedure, nopass :: add => add_J_adder
!
   end type J_adder
!
contains
!
!
   subroutine add_J_adder(ao, D, G, eri, s1, s2, s3, s4, thread, n_threads, degeneracy)
!!
!!    Add
!!    Written by Sarai D. Folkestad, 2021
!!
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
      integer :: w, x, y, z, w_red, x_red, y_red, z_red, wxyz
      real(dp) :: D_yz, D_wx
!
      do z = ao%shells(s4)%first, ao%shells(s4)%get_last()
         z_red = z - ao%shells(s4)%first + 1
!
         do y = ao%shells(s3)%first, ao%shells(s3)%get_last()
            y_red = y - ao%shells(s3)%first + 1
!
            D_yz = D(y, z)
!
            do x = ao%shells(s2)%first, ao%shells(s2)%get_last()
               x_red = x - ao%shells(s2)%first + 1
!
               do w = ao%shells(s1)%first, ao%shells(s1)%get_last()
                  w_red = w - ao%shells(s1)%first + 1
!
                  D_wx = D(w, x)
!
                  wxyz = ao%shells(s1)%length &
                       *(ao%shells(s2)%length &
                       *(ao%shells(s3)%length &
                       *(z_red-1)+y_red-1)+x_red-1)+w_red
!
                  eri(wxyz) = degeneracy * eri(wxyz)
!
                  G(w, x, thread + 1) = G(w, x, thread + 1) + eri(wxyz) * D_yz
                  G(y, z, thread + 1) = G(y, z, thread + 1) + eri(wxyz) * D_wx
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine add_J_adder
!
!
end module J_adder_class
