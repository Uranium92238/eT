module atomic_class
!
!!
!!    Atom class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use shell_class
   use disk_manager_class
!
   implicit none
!
   type :: atomic
!
      character(len=2) :: symbol
      integer(i15)     :: number
!
      character(len=100) :: basis
!
      integer(i15) :: n_shells
      type(shell), dimension(:,:), allocatable :: shells
!
      real(dp) :: x
      real(dp) :: y
      real(dp) :: z
!
   contains
!
      procedure          :: set_number      => set_number_atom
      procedure, private :: symbol_2_number => symbol_2_number_atom
!
   end type atomic
!
contains
!
!
   subroutine set_number_atom(atom)
!!
!!    Set number
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Sets the atomic number from atomic symbol
!!
      implicit none
!
      class(atomic) :: atom
!
      call atom%symbol_2_number()
!
   end subroutine set_number_atom
!
!
   subroutine symbol_2_number_atom(atom)
!!
!!    Symbol to number
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Uses the periodic table to determine atomic number
!!    from atomic symbol
!!
      use periodic_table
!
      implicit none
!
      class(atomic) :: atom
!
      integer(i15) :: i = 0
!
      atom%number = 0
!
      do i = 1, size_periodic_table
         if (atomic_symbol(i) == atom%symbol) then
!
            atom%number = i
!
            return
!
         endif
      enddo
!
      if (atom%number == 0) then
!
         write(output%unit,'(/t3,a)') 'Error: illegal atomic symbol, check the eT.inp file '
         stop
!
      endif
!
   end subroutine symbol_2_number_atom
!
!
end module atomic_class
