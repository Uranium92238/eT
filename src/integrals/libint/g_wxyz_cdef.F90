interface
!
   subroutine get_ao_g_wxyz_c(g, s1, s2, s3, s4) bind(C, name='get_ao_g_wxyz')
!
      use iso_c_binding
      implicit none
!
      real(c_double), dimension(1,1) :: g
      integer(c_long) :: s1, s2, s3, s4
!
   end subroutine get_ao_g_wxyz_c
!
end interface
