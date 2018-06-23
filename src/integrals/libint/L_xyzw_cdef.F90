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
   subroutine get_ao_L_xyzw_c(L, s1, s3) bind(C, name='get_ao_L_xyzw')
!
      use iso_c_binding
      implicit none
!
      integer(c_int) :: s1
      integer(c_int) :: s3
!
      real(c_double), dimension(1,1) :: L
!
   end subroutine get_ao_L_xyzw_c
!
end interface
