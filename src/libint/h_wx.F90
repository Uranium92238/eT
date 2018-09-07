module h_wx
!
   use kinds
   use iso_c_binding
!
   include "h_wx_cdef.F90"
!
contains
!
   subroutine construct_ao_h_wx(h, s1, s2)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: h
      integer(kind=8) :: s1, s2
!
      call construct_ao_h_wx_c(h, s1, s2)
!
   end subroutine construct_ao_h_wx
!
end module h_wx
