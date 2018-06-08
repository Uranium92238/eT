module atom_class
!
!!
!!    Atom class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
!
   implicit none
!
   type :: atom
!
      character(len=2)  :: symbol
      integer(i15)      :: number 
!
      character(len=100) :: basis
!
      real(dp) :: x
      real(dp) :: y
      real(dp) :: z
!
   contains
!
      procedure            :: set_number      => set_number_atom
      procedure, private   :: symbol_2_number => symbol_2_number_atom
!
   end type atom
!
contains
!
!
   subroutine set_number_atom(a)
!!
!!    Set number
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Sets the atomic number from atomic symbol
!!
      implicit none
!
      class(atom) :: a 
!
      call a%symbol_2_number()
      write(output%unit,*) a%number
!
   end subroutine set_number_atom
!
!
   subroutine symbol_2_number_atom(a)
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
      class(atom) :: a
!
      integer(i15) :: i = 0
!
      do i = 1, size_periodic_table
         if (atomic_symbol(i) == a%symbol) then
!
            a%number = i
!
            return
!
         endif
      enddo
!
      if (a%number == 0) then
!
         write(output%unit,'(/t3,a)') 'Error: illegal atomic symbol, check the eT.inp file '
         stop
!
      endif
!
   end subroutine symbol_2_number_atom
!
!
end module atom_class
