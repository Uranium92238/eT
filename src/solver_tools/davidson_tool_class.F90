module davidson_tool_class
!
!!
!!    Abstract Davidson solver class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2018
!!
!!
!
   use kinds
   use parameters
!
   use file_class
   use disk_manager_class
   use memory_manager_class
!
!
   type, abstract :: davidson_tool
!
   contains
!
!     Read and write routines
!
      procedure, non_overridable :: read    => read_davidson_tool
      procedure, non_overridable :: write   => write_davidson_tool
!
      procedure, non_overridable :: read_trials    => read_trials_davidson_tool
      procedure, non_overridable :: write_trials   => write_trials_davidson_tool
!
      procedure, non_overridable :: read_transforms  => read_transforms_davidson_tool
      procedure, non_overridable :: write_transforms => write_transforms_davidson_tool  
!
!     Other procedures
!
      procedure, non_overridable :: construct_reduced_matrix => construct_reduced_matrix_davidson_tool
!
!     Deferred routines
!
      procedure(essential_davidson), deferred :: solve_reduced_problem
      procedure(essential_davidson), deferred :: construct_residual
      procedure(essential_davidson), deferred :: construct_solution
!
   end type davidson_tool
!
!
   abstract interface
!
      subroutine essential_davidson(davidson)
!
         import :: davidson_tool
!
         implicit none 
!
         class(davidson_tool) :: davidson 
!
      end subroutine essential_davidson
!
   end interface
!
contains
!
!
   subroutine read_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine read_davidson_tool
!
!
   subroutine write_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine write_davidson_tool
!
!
   subroutine read_trials_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine read_trials_davidson_tool
!
!
   subroutine write_trials_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine write_trials_davidson_tool
!
!
   subroutine read_transforms_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine read_transforms_davidson_tool
!
!
   subroutine write_transforms_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine write_transforms_davidson_tool
!
!
   subroutine construct_reduced_matrix_davidson_tool(davidson)
!!
!!
!!
      implicit none
!
      class(davidson_tool) :: davidson
!
   end subroutine construct_reduced_matrix_davidson_tool
!
end module davidson_tool_class
