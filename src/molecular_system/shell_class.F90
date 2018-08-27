module shell_class
!
!!
!!    Shell class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
!
   implicit none
!
   type :: shell
!
      integer(i15)    :: size  = -1  ! The number of basis functions
      integer(i15)    :: first = -1  ! The first AO index
      integer(i15)    :: last  = -1  ! The last AO index
      integer(i15)    :: l     = -1  ! The angular momentum
      integer(kind=4) :: number = -1 ! The shell number (according to the ordering given by Libint)
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
      integer(i15) :: i
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
         write(output%unit, '(/t3,a)') 'Error: could not determine angular momentum of shell.'
         stop
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
