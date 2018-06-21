interface
!
   subroutine get_n_basis_in_shells_c(atom, nbis) bind(C, name='get_n_basis_in_shells')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int) :: atom
      integer(c_int), dimension(1,1) :: nbis
!
   end subroutine get_n_basis_in_shells_c
!
   subroutine get_n_shells_on_atoms_c(nsoa) bind(C, name='get_n_shells_on_atoms')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int), dimension(1,1) :: nsoa
!
   end subroutine get_n_shells_on_atoms_c
!
end interface
