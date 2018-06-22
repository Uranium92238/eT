module atom_init
!
   use kinds
   use iso_c_binding
!
   include "atom_init_cdef.F90"
!
contains
!
   subroutine get_first_ao_in_shells(atom, faois)
!
      implicit none
!
      integer(kind=4) :: atom
      integer(kind=4), dimension(1,1) :: faois
!
      call get_first_ao_in_shells_c(atom, faois)
!
   end subroutine get_first_ao_in_shells
!
   subroutine get_n_basis_in_shells(atom, nbis)
!
      implicit none
!
      integer(kind=4) :: atom
      integer(kind=4), dimension(1,1) :: nbis
!
      call get_n_basis_in_shells_c(atom, nbis)
!
   end subroutine get_n_basis_in_shells
!
   subroutine get_n_shells_on_atoms(nsoa)
!
      implicit none
!
      integer(kind=4), dimension(1,1) :: nsoa
!
      call get_n_shells_on_atoms_c(nsoa)
!
   end subroutine get_n_shells_on_atoms
!
end module atom_init
