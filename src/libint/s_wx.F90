module s_wx
!
   use kinds
   use iso_c_binding
!
   include "s_wx_cdef.F90"
!
contains
!
   subroutine construct_ao_s_wx(s, s1, s2)
!
      implicit none
!
      real(dp), dimension(1,1) :: s
      integer(i6) :: s1, s2
!
      call construct_ao_s_wx_c(s, s1, s2)
!
   end subroutine construct_ao_s_wx
!
end module s_wx
