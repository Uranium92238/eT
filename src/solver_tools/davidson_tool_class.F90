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
      procedure :: read    => read_davidson_tool
      procedure :: write   => write_davidson_tool
!
      procedure :: read_trials    => read_trials_davidson_tool
      procedure :: write_trials   => write_trials_davidson_tool
!
      procedure :: read_transforms 
      procedure :: write_transforms      
!
!     
!
    !  procedure(essential_davidson), deferred :: construct_reduced_matrix
!
      procedure(essential_davidson), deferred :: solve_reduced_problem
!
      procedure(essential_davidson), deferred :: construct_residual
!
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
end module davidson_tool_class
