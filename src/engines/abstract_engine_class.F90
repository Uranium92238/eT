module abstract_engine_class
!
!!
!!    Abstract engine class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use ccs_class
!
   implicit none
!
   type, abstract :: abstract_engine
!
      character(len=100) :: name_
!
   contains
!
      procedure(essential_engine), deferred :: prepare 
      procedure(essential_engine), deferred :: cleanup   
      procedure(essential_engine_w_wf), deferred :: run        
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
!
      subroutine essential_engine_w_wf(engine, wf)
!
         import :: abstract_engine, ccs
!
         implicit none 
!
         class(abstract_engine) :: engine
!
         class(ccs) :: wf
!
      end subroutine essential_engine_w_wf
!
!
   end interface
!
contains
!
!
end module abstract_engine_class
!
