module libint_initialization
!
   use kinds
   use iso_c_binding
!
   include "libint_initialization_cdef.F90"
!
contains
!
   subroutine initialize_overlap()
!
      implicit none
!
      call initialize_overlap_c()
!
   end subroutine initialize_overlap
!
   subroutine initialize_coulomb()
!
      implicit none
!
      call initialize_coulomb_c()
!
   end subroutine initialize_coulomb
!
   subroutine initialize_kinetic()
!
      implicit none
!
      call initialize_kinetic_c()
!
   end subroutine initialize_kinetic
!
   subroutine initialize_nuclear()
!
      implicit none
!
      call initialize_nuclear_c()
!
   end subroutine initialize_nuclear
!
   subroutine initialize_basis(basis_set)
!
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      implicit none
!
      character(len=*) :: basis_set
      character(len=40) :: basis_set_temp
      character(kind=c_char),dimension(40) :: cpp_temp
!
      integer(kind=4) :: j
!
      basis_set_temp = trim(basis_set)//c_null_char
      basis_set_temp = trim(basis_set_temp)
!
      write(*,*)basis_set_temp, len_trim(basis_set_temp)
!
      do j=1,len_trim(basis_set_temp)
!
         cpp_temp(j) = basis_set_temp(j:j)  
!
      enddo
!
      call initialize_basis_c(cpp_temp)
!
   end subroutine initialize_basis
!
end module libint_initialization
