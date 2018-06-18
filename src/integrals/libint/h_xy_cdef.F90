interface
!
   subroutine get_ao_h_xy_c(h) bind(C, name='get_ao_h_xy')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(1,1) :: h
!
   end subroutine get_ao_h_xy_c
!
   subroutine get_n_aos_c(n_ao) bind(C, name='get_n_aos')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int64_t) :: n_ao
!
   end subroutine get_n_aos_c
!
end interface
