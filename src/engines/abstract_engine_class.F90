module abstract_engine_class
!
!!
!!    Abstract engine class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
!
   implicit none
!
   type, abstract :: abstract_engine
!
   contains
!
      procedure(essential_engine), deferred :: initialize 
      procedure(essential_engine), deferred :: finalize   
      procedure(essential_engine), deferred :: run        
!
      procedure(essential_engine), deferred :: print_banner  
      procedure(essential_engine), deferred :: print_summary 
!
   end type abstract_engine
!
!
   abstract interface
!
      subroutine essential_engine(engine)
!
         import :: abstract_engine
!
         implicit none 
!
         class(abstract_engine) :: engine
!
      end subroutine essential_engine
!
   end interface
!
contains
!
!
end module abstract_engine_class
!