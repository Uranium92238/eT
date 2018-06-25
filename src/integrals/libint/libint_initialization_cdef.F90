interface
!
   subroutine initialize_basis_c() bind(C, name='initialize_basis')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_basis_c
!
   subroutine initialize_coulomb_c() bind(C, name='initialize_coulomb')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_coulomb_c
!
end interface
