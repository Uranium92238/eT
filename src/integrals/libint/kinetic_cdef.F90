interface
!
   subroutine get_ao_xy_kinetic_c(h) bind(C, name='get_ao_xy_kinetic')
!
      use iso_c_binding
      implicit none
!
      real(c_double), value :: h
!
   end subroutine get_ao_xy_kinetic_c
!
end interface
