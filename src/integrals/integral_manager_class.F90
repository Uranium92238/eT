module integral_manager_class
!
!!
!!    Integral_manager class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
!  Fortran interfaces to C++ routines
!
   use h_xy
   use s_xy
   use g_wxyz
   use L_wxyz
   use libint_initialization
!
!  Disk & memory class modules
!
   use file_class
   use memory_manager_class
!
   use molecular_system_class
!
   implicit none
!
!  Class definition
!
   type :: integral_manager
!
!     No attributes yet
!
   contains
!
      procedure :: cholesky_decompose => cholesky_decompose_integral_manager
!
      procedure :: get_ao_h_xy   => get_ao_h_xy_integral_manager   ! h_αβ
      procedure :: get_ao_s_xy   => get_ao_s_xy_integral_manager   ! s_αβ
      procedure :: get_ao_g_wxyz => get_ao_g_wxyz_integral_manager ! g_αβγδ
!
   end type integral_manager
!
!
!   interface
!
!
!      module subroutine cholesky_decompose_integral_manager(integrals, molecule)
!!
!!       Cholesky decompose
!!       Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!         implicit none
!
!         class(integral_manager) :: integrals
!         class(molecular_system) :: molecule
!
!      end subroutine cholesky_decompose_integral_manager
!
!
!   end interface
!
!
contains
!
!
   subroutine get_ao_h_xy_integral_manager(int, h)
!!
!!    Get h_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the h_αβ integral_manager in the array h.
!!
      implicit none
!
      class(integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: h
!
      call get_ao_h_xy(h)
!
   end subroutine get_ao_h_xy_integral_manager
!
!
   subroutine get_ao_s_xy_integral_manager(int, s)
!!
!!    Get s_αβ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the s_αβ integral_manager in the array s.
!!
      implicit none
!
      class(integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: s
!
      call get_ao_s_xy(s)
!
   end subroutine get_ao_s_xy_integral_manager
!
!
   subroutine get_ao_g_wxyz_integral_manager(int, g, s1, s2, s3, s4)
!!
!!    Get g_αβγδ integrals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Fortran wrapper for the C++ routine that calculates and
!!    saves the g_αβγδ integral_manager in the array g.
!!
!!    s1 is first shell index and so on.
!!
      implicit none
!
      class(integral_manager) :: int
!
      real(kind=8), dimension(:,:), intent(inout) :: g
!
      integer(kind=8), intent(in) :: s1, s2, s3, s4
!
      call get_ao_g_wxyz(g, s1, s2, s3, s4)
!
   end subroutine get_ao_g_wxyz_integral_manager
!
!
end module integral_manager_class
