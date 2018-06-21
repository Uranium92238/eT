module atom_init
!
   use kinds
   use iso_c_binding
!
   include "atom_init_cdef.F90"
!
contains
!
   subroutine get_n_shells_on_atom(nsoa)
!
      implicit none
!
      integer(kind=4), dimension(1,1) :: nsoa
!
      call get_n_shells_on_atom_c(nsoa)
!
   end subroutine get_n_shells_on_atom
!
end module atom_init
