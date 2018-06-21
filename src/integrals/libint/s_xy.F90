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
      integer(kind=8) :: ns
!
      call get_n_shells_c(ns)
!
   end subroutine get_n_shells
!
   subroutine get_ao_s_xy(s)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: s
!
      call get_ao_s_xy_c(s)
!
   end subroutine get_ao_s_xy
!
end module s_xy
