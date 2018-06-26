module g_wxyz
!
   use kinds
   use iso_c_binding
!
   include "g_wxyz_cdef.F90"
!
contains
!
   subroutine get_ao_g_wxyz(g)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: g
!
      call get_ao_g_wxyz_c(g)
!
   end subroutine get_ao_g_wxyz
!
end module g_wxyz
