module cc2_class
!
!!
!!    Coupled cluster singles and perturbative doubles (cc2) 
!!    class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: cc2
!
      real(dp), dimension(:,:), allocatable :: t2    
!
      integer(i15) :: n_t2  
!
   contains
!
!
   end type cc2
!
!
   interface
!
!
   end interface 
!
!
contains
!
end module cc2_class
