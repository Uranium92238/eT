interface
!
   integer function get_n_shells_on_atom_c(i) bind(C, name='get_n_shells_on_atom')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int64_t) :: i
!
   end function get_n_shells_on_atom_c
!
end interface
