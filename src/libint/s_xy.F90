module s_xy
!
   use kinds
   use iso_c_binding
!
   include "s_xy_cdef.F90"
!
contains
!
   subroutine get_n_shells(ns)
!
      implicit none
!
      integer(i15) :: ns
!
      call get_n_shells_c(ns)
!
   end subroutine get_n_shells
!
   subroutine get_ao_s_xy(s)
!
      implicit none
!
      real(dp), dimension(1,1) :: s
!
      call get_ao_s_xy_c(s)
!
   end subroutine get_ao_s_xy
!
   subroutine construct_ao_s_wx(s, s1, s2)
!
      implicit none
!
      real(dp), dimension(1,1) :: s
      integer(i15) :: s1, s2
!
      call construct_ao_s_wx_c(s, s1, s2)
!
   end subroutine construct_ao_s_wx
!
end module s_xy
