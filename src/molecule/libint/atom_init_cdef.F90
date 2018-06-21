interface
!
   subroutine get_n_shells_on_atom_c(nsoa) bind(C, name='get_n_shells_on_atom')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int), dimension(1,1) :: nsoa
!
   end subroutine get_n_shells_on_atom_c
!
end interface
