module kinetic
!
   use kinds
   use iso_c_binding
!
   include "kinetic_cdef.F90"
!
contains
!
   subroutine get_ao_xy_kinetics(h)
!
      implicit none
!
      real(kind=8), dimension(2,1), intent(inout) :: h
!
      call get_ao_xy_kinetic_c(h)
!
   end subroutine get_ao_xy_kinetics
!
end module kinetic
