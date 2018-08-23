module h_xy
!
   use kinds
   use iso_c_binding
!
   include "h_xy_cdef.F90"
!
contains
!
   subroutine get_ao_h_xy(h)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: h
!
      call get_ao_h_xy_c(h)
!
   end subroutine get_ao_h_xy
!
   subroutine get_ao_h_xy_sp(h, s1, s2)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: h
      integer(kind=8) :: s1, s2
!
      call get_ao_h_xy_sp_c(h, s1, s2)
!
   end subroutine get_ao_h_xy_sp
!
   subroutine get_n_aos(n_ao)
!
      use iso_c_binding
!
      implicit none
!
      integer(i15) :: n_ao
!
      call get_n_aos_c(n_ao)
!
   end subroutine get_n_aos
!
end module h_xy
