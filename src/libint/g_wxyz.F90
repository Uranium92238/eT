module g_wxyz
!
   use kinds
   use iso_c_binding
!
   include "g_wxyz_cdef.F90"
!
contains
!
   subroutine construct_ao_g_wxyz(g, s1, s2, s3, s4)
!
      implicit none
!
      real(dp), dimension(1,1) :: g
!
      integer(i6) :: s1, s2, s3, s4
!
      call construct_ao_g_wxyz_c(g, s1, s2, s3, s4)
!
   end subroutine construct_ao_g_wxyz
!
   subroutine construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, epsilon, thread, skip, n1, n2, n3, n4)
!
      implicit none
!
      real(dp), dimension(1,1) :: g
      real(dp) :: epsilon
!
      integer(i6) :: s1, s2, s3, s4, thread, skip, n1, n2, n3, n4 
!
      call construct_ao_g_wxyz_epsilon_c(g, s1, s2, s3, s4, epsilon, thread, skip, n1, n2, n3, n4)
!
   end subroutine construct_ao_g_wxyz_epsilon
!
end module g_wxyz
