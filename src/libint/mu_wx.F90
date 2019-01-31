module mu_wx
!
   use kinds
   use iso_c_binding
!
   include "mu_wx_cdef.F90"
!
contains
!
   subroutine construct_ao_mu_wx(mu_X, mu_Y, mu_Z, s1, s2)
!
      implicit none
!
      real(dp), dimension(1,1) :: mu_X, mu_Y, mu_Z
      integer(i6) :: s1, s2
!
      call construct_ao_mu_wx_c(mu_X, mu_Y, mu_Z, s1, s2)
!
   end subroutine construct_ao_mu_wx
!
end module mu_wx
