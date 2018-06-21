module atom_init
!
   use kinds
   use iso_c_binding
!
   include "atom_init_cdef.F90"
!
contains
!
   integer(i15) function get_n_shells_on_atom(i)
!
      implicit none
!
      integer(i15) :: i
!
      get_n_shells_on_atom = get_n_shells_on_atom_c(i)
!
   end function get_n_shells_on_atom
!
end module atom_init
