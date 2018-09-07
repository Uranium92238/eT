interface
!
   subroutine get_n_shells_c(ns) bind(C, name='get_n_shells')
!
      use iso_c_binding
!
      implicit none
!
      integer(c_int64_t) :: ns
!
   end subroutine get_n_shells_c
!
   subroutine get_ao_s_xy_c(s) bind(C, name='get_ao_s_xy')
!
      use iso_c_binding
!
      implicit none
!
      real(c_double), dimension(1,1) :: s
!
   end subroutine get_ao_s_xy_c
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
