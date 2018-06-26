module L_wxyz
!
   use kinds
   use iso_c_binding
!
   include "L_wxyz_cdef.F90"
!
contains
!
   subroutine get_ao_L_wxyz(L, s1, s3)
!
      implicit none
!
      integer(kind=4) :: s1
      integer(kind=4) :: s3
!
      real(kind=8), dimension(1,1) :: L
!
      call get_ao_L_wxyz_c(L, s1, s3)
!
   end subroutine get_ao_L_wxyz
!
end module L_wxyz
