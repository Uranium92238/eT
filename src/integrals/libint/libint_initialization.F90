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
   subroutine initialize_basis(basis_set, mol_name)
!
      use iso_c_binding, only: C_CHAR, C_NULL_CHAR
      implicit none
!
      character(len=*) :: basis_set
      character(len=40) :: basis_set_temp
      character(kind=c_char),dimension(40) :: cpp_temp_basis
!
      character(len=*) :: mol_name
      character(len=40) :: mol_name_temp
      character(kind=c_char),dimension(40) :: cpp_temp_mol_name
!
      integer(kind=4) :: j
!
      basis_set_temp = trim(basis_set)//c_null_char
      mol_name_temp = trim(mol_name)//c_null_char
!
      write(*,*)basis_set_temp, len_trim(basis_set_temp)
      write(*,*)mol_name_temp, len_trim(mol_name_temp)
!
      do j=1,len_trim(basis_set_temp)
!
         cpp_temp_basis(j) = basis_set_temp(j:j)  
!
      enddo
!
      do j=1,len_trim(mol_name_temp)
!
         cpp_temp_mol_name(j) = mol_name_temp(j:j)  
!
      enddo
!
      call initialize_basis_c(cpp_temp_basis, cpp_temp_mol_name)
!
   end subroutine initialize_basis
!
end module libint_initialization
