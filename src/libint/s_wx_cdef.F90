interface
!
   subroutine construct_ao_s_wx_c(s, s1, s2) bind(C, name='construct_ao_s_wx')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(1,1) :: s
      integer(c_long) :: s1, s2
!
   end subroutine construct_ao_s_wx_c
!
end interface
