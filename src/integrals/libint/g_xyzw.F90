module g_xyzw
!
   use kinds
   use iso_c_binding
!
   include "g_xyzw_cdef.F90"
!
contains
!
   subroutine get_ao_g_xyzw(g)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: g
!
      call get_ao_g_xyzw_c(g)
!
   end subroutine get_ao_g_xyzw
!
end module g_xyzw
