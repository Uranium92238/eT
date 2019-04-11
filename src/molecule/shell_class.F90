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
   use output_file_class
   use disk_manager_class
!
   implicit none
!
   type :: shell
!
      integer    :: size  = -1  ! The number of basis functions
      integer    :: first = -1  ! The first AO index
      integer    :: last  = -1  ! The last AO index
      integer    :: l     = -1  ! The angular momentum
      integer    :: number_ = -1 ! The shell number (according to the ordering given by Libint)
!
   contains
!
      procedure :: determine_angular_momentum => determine_angular_momentum_shell
      procedure :: determine_last_ao_index    => determine_last_ao_index_shell
!
   end type shell
!
contains
!
   subroutine determine_angular_momentum_shell(sh)
!!
!!    Determine angular momentum
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Determines the angular momentum by counting the number of basis
!!    functions in the shell (2l + 1). If l >= 10 or the shell size is
!!    not 2l + 1 for some l = 1, 2, 3, ..., 9, the routine will print an
!!    error and stop.
!!
      implicit none
!
      class(shell) :: sh
!
      integer :: i
!
      i = 0
!
      do while (i .lt. 10)
!
         if ((2*i + 1) .eq. sh%size) then
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
      sh%last = sh%first + sh%size - 1
!
   end subroutine determine_last_ao_index_shell
!
!
end module
