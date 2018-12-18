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
   contains
!
      procedure :: construct_omega  => construct_omega_cc2
      procedure :: omega_cc2_a1     => omega_cc2_a1_cc2
      procedure :: omega_cc2_b1     => omega_cc2_b1_cc2
      procedure :: omega_cc2_c1     => omega_cc2_c1_cc2
!
   end type cc2
!
   interface
!
      include "../submodules/cc2/omega_cc2_interface.F90"
!
   end interface 
!
!
contains
!
end module cc2_class
