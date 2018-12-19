submodule (cc2_class) jacobian
!
!!
!!    Jacobian submodule (CC2)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Routines for the linear transform of trial
!!    vectors by the transpose of the Jacobian matrix 
!!
!!    σ_i = A^T * b_i,
!!
!!    where
!!   
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | ν >.
!! 
!! 
!
   implicit none 
!
!
contains
!
!
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
!
!
end submodule jacobian
