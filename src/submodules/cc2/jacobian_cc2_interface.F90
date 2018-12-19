
   subroutine effective_jacobian_transformation_cc2(wf, b)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none 
!
      class(cc2) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1) :: b
!
   end subroutine effective_jacobian_transformation_cc2