interface
!
   subroutine get_ao_L_wxyz_c(L, s1, s3) bind(C, name='get_ao_L_wxyz')
!
      use iso_c_binding
      implicit none
!
      integer(c_long) :: s1
      integer(c_long) :: s3
!
      real(c_double), dimension(1,1) :: L
!
   end subroutine get_ao_L_wxyz_c
!
end interface
