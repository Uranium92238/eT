interface
!
   subroutine construct_ao_mu_wx_c(mu_X, mu_Y, mu_Z, s1, s2) bind(C, name='construct_ao_mu_wx')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(1,1) :: mu_X, mu_Y, mu_Z
      integer(c_long) :: s1, s2
!
   end subroutine construct_ao_mu_wx_c
!
end interface
