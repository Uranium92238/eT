module L_xyzw
!
   use kinds
   use iso_c_binding
!
   include "L_xyzw_cdef.F90"
!
contains
!
   subroutine initialize_coulomb()
!
      implicit none
!
      call initialize_coulomb_c()
!
   end subroutine initialize_coulomb
!
   subroutine initialize_basis()
!
      implicit none
!
      call initialize_basis_c()
!
   end subroutine initialize_basis
!
   subroutine get_ao_L_xyzw(L, s1, s3)
!
      implicit none
!
      integer(kind=4) :: s1
      integer(kind=4) :: s3
!
      real(kind=8), dimension(1,1) :: L
!
      call get_ao_L_xyzw_c(L, s1, s3)
!
   end subroutine get_ao_L_xyzw
!
end module L_xyzw
