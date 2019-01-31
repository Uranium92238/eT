module atom_init
!
   use kinds
   use iso_c_binding
!
   include "atom_init_cdef.F90"
!
contains
!
   subroutine get_shell_numbers(atom, sn)
!
      implicit none
!
      integer :: atom
      integer(i6), dimension(1,1) :: sn
      integer(i6) :: atom4 
!
      atom4 = int(atom,i6)
!
      call get_shell_numbers_c(atom4, sn)
!
   end subroutine get_shell_numbers
!
   subroutine get_first_ao_in_shells(atom, faois)
!
      implicit none
!
      integer :: atom
      integer(i6), dimension(1,1) :: faois
      integer(i6) :: atom4 
!
      atom4 = int(atom,i6)
!
      call get_first_ao_in_shells_c(atom4, faois)
!
   end subroutine get_first_ao_in_shells
!
   subroutine get_n_basis_in_shells(atom, nbis)
!
      implicit none
!
      integer :: atom
      integer(i6), dimension(1,1) :: nbis
      integer(i6) :: atom4 
!
      atom4 = int(atom,i6)
!
      call get_n_basis_in_shells_c(atom4, nbis)
!
   end subroutine get_n_basis_in_shells
!
   subroutine get_n_shells_on_atoms(nsoa)
!
      implicit none
!
      integer(i6), dimension(1,1) :: nsoa
!
      call get_n_shells_on_atoms_c(nsoa)
!
   end subroutine get_n_shells_on_atoms
!
end module atom_init
