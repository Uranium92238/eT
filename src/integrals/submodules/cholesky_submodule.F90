submodule (integral_manager_class) cholesky
!
!!
!!
!
   implicit none
!
!
contains
!
!
   module subroutine cholesky_decompose_integral_manager(integral, molecule)
!!
!!
!!
!!
      implicit none
!
      class(integral_manager) :: integral
      class(molecular_system) :: molecule
!
      write(output%unit,*) ' Hei hei! '
!
   end subroutine cholesky_decompose_integral_manager

!
!
end submodule cholesky
