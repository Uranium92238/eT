interface
!
   subroutine construct_ao_h_wx_c(h, s1, s2) bind(C, name='construct_ao_h_wx')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(1,1) :: h
      integer(c_int) :: s1, s2
!
   end subroutine construct_ao_h_wx_c
!
end interface
