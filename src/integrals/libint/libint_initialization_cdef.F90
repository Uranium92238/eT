interface
!
   subroutine initialize_libint_c() bind(C, name='initialize_libint')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_libint_c
!
   subroutine finalize_libint_c() bind(C, name='finalize_libint')
!
      use iso_c_binding
      implicit none
!
   end subroutine finalize_libint_c
!
   subroutine initialize_basis_c() bind(C, name='initialize_basis')
!
      use iso_c_binding, only: C_CHAR
      implicit none
!
      character(kind = c_char) :: basisset(*)
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
   subroutine initialize_kinetic_c() bind(C, name='initialize_kinetic')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_kinetic_c
!
   subroutine initialize_overlap_c() bind(C, name='initialize_overlap')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_overlap_c
!
   subroutine initialize_nuclear_c() bind(C, name='initialize_nuclear')
!
      use iso_c_binding
      implicit none
!
   end subroutine initialize_nuclear_c
!
end interface
