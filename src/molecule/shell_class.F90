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
module shell_class
!
!!
!!    Shell class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use interval_class
   use global_out, only : output
   use shell_details_class, only : shell_details
!
   implicit none
!
   type, extends(interval) :: shell ! interval: AO index range of the shell 
!
      type(shell_details) :: basis_details
!
      integer    :: size_cart = -1 ! The number of basis functions in cartesian
      integer    :: l         = -1 ! The angular momentum
      integer    :: number_   = -1 ! The shell number (according to the ordering given by Libint)
!
   contains
!
      procedure :: determine_angular_momentum => determine_angular_momentum_shell
      procedure :: determine_last_ao_index    => determine_last_ao_index_shell
!
   end type shell
!
!
   interface shell 
!
      procedure :: new_shell 
!
   end interface shell 
!
!
contains
!
!
   function new_shell(first, length, number_) result(sh)
!!
!!    New shell 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    first:   the first AO index of the shell 
!!    length:  the number of AOs in the shell 
!!    number_: the shell number in the full list of shells (according to Libint)
!!
      implicit none 
!
      integer, intent(in) :: first, length 
!
      integer, intent(in) :: number_ 
!
      type(shell) :: sh 
!
      sh%first = first 
      sh%length = length 
!
      call sh%determine_last_ao_index()
      call sh%determine_angular_momentum()
!
      sh%number_ = number_
!
   end function new_shell
!
!
   subroutine determine_angular_momentum_shell(sh)
!!
!!    Determine angular momentum
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the angular momentum by counting the number of basis
!!    functions in the shell (l < 10)
!!
      implicit none
!
      class(shell) :: sh
!
      integer :: i
!
      i = 0
      sh%l = -1
!
      do while (i .lt. 10)
!
         if ((2*i + 1) .eq. sh%length .or. (((i+1)*(i+2))/2) .eq. sh%length) then
!
            sh%l = i
!
         endif
!
         i = i + 1
!
      enddo
!
      if (sh%l .eq. -1) then
!
         call output%error_msg('could not determine angular momentum of shell.')
!
      endif
!
   end subroutine determine_angular_momentum_shell
!
!
   subroutine determine_last_ao_index_shell(sh)
!!
!!    Determine last AO index
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(shell) :: sh
!
      sh%last = sh%first + sh%length - 1
!
   end subroutine determine_last_ao_index_shell
!
!
end module shell_class
