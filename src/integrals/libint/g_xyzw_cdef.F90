interface
!
   subroutine get_ao_g_xyzw_c(g) bind(C, name='get_ao_g_xyzw')
!
      use iso_c_binding
      implicit none
!
      real(c_double), dimension(1,1) :: g
!
   end subroutine get_ao_g_xyzw_c
!
end interface
