interface
!
   subroutine get_ao_s_xy_c(s) bind(C, name='get_ao_s_xy')
!
      use iso_c_binding
      implicit none
!
      real(c_double), dimension(1,1) :: s
!
   end subroutine get_ao_s_xy_c
!
end interface
