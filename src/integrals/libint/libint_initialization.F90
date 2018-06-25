module libint_initialization
!
   use kinds
   use iso_c_binding
!
   include "libint_initialization_cdef.F90"
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
end module libint_initialization
