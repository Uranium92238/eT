module mo_integral_tool_class
!
!!
!!    MO integral tool class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use file_class
   use memory_manager_class
!
   implicit none
!
!  Class definition
!
   type :: mo_integral_tool
!
      logical, private :: eri_file           = .false.
      logical, private :: eri_t1_file        = .false. 
      logical, private :: cholesky_file      = .true.
      logical, private :: cholesky_t1_file   = .false.
!
   contains
!
      procedure :: need_t1 => need_t1_mo_integral_tool
!
   end type mo_integral_tool
!
!
contains
!
!
   logical function need_t1_mo_integral_tool(integrals)
!!
!!    Need t1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    Returns true if neither t1-transformed integrals nor t1-transformed cholesky
!!    vectors are on file.
!!
      implicit none
!
      class(mo_integral_tool) :: integrals 
!
      need_t1_mo_integral_tool = ((.not. integrals%eri_t1_file) .and. (.not. integrals%cholesky_t1_file))
!
   end function need_t1_mo_integral_tool
!
!
end module mo_integral_tool_class

