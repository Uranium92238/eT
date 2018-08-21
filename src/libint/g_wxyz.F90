module g_wxyz
!
   use kinds
   use iso_c_binding
!
   include "g_wxyz_cdef.F90"
!
contains
!
   subroutine get_ao_g_wxyz(g, s1, s2, s3, s4)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: g
!
      integer(kind=8) :: s1, s2, s3, s4
!
      call get_ao_g_wxyz_c(g, s1, s2, s3, s4)
!
   end subroutine get_ao_g_wxyz
!
   subroutine get_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, epsilon, thread)
!
      implicit none
!
      real(kind=8), dimension(1,1) :: g
      real(kind=8) :: epsilon
!
      integer(kind=8) :: s1, s2, s3, s4, thread 
!
      call get_ao_g_wxyz_epsilon_c(g, s1, s2, s3, s4, epsilon, thread)
!
   end subroutine get_ao_g_wxyz_epsilon
!
end module g_wxyz
