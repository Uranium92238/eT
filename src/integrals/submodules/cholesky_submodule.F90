submodule (integral_manager_class) cholesky
!
!!
!!    Cholesky submodule
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
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
!!    Cholesky decompose
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(integral_manager) :: integral
      class(molecular_system) :: molecule
!
      real(dp), dimension(:,:), allocatable :: D_xy
!
      call integral%construct_two_electron_diagonal(D_xy)
!
   end subroutine cholesky_decompose_integral_manager
!
!
   module subroutine construct_two_electron_diagonal_integral_manager(integral, D_xy)
!!
!!    Construct two electron diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(integral_manager) :: integral
!
      real(dp), dimension(:,:) :: D_xy
!
   end subroutine construct_two_electron_diagonal_integral_manager
!
!
end submodule cholesky
