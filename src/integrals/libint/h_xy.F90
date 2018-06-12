module h_xy
!
   use kinds
   use iso_c_binding
!
   include "h_xy_cdef.F90"
!
contains
!
   subroutine get_ao_xy(h)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: h
!
      call get_ao_xy_c(h)
!
   end subroutine get_ao_xy
!
   subroutine get_n_aos(n_ao)
!
      use iso_c_binding
      implicit none
!
      integer(i15) :: n_ao
!
      call get_n_aos_c(n_ao)
!
   end subroutine get_n_aos
!
end module h_xy
