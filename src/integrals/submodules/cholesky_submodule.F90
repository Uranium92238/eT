submodule (integral_manager_class) cholesky
!
!!
!!    Cholesky submodule
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
   use index
   use reordering
   use array_utilities
   use array_analysis
!
   implicit none
!
!
contains
!
!
!

!
!
end submodule cholesky
!
     !call mem%alloc(approximate_diagonal, size_AB, 1)
!
     !approximate_diagonal = zero
!
     !do I = 1, size_AB
     !   do K = 1, n_cholesky
     !      approximate_diagonal(I,1) = approximate_diagonal(I,1) + L_K_yz(K, I)**2
     !   enddo
     !enddo
!
     !call mem%alloc(difference, size_AB, 1)
     !difference = zero
     !do I = 1, size_AB
     !      difference(I, 1) = approximate_diagonal(I,1) - first_sig_D_xy(I, 1)
     !enddo
!
     !max_diff = zero
     ! do I = 1, size_AB
     !   if (abs(difference(I, 1)) .gt. max_diff) max_diff = abs(difference(I, 1))
     !enddo
     !write(output%unit, '(a60, e12.4)')'Maximal difference between approximate and actual diagonal: ', max_diff
     !call mem%dealloc(L_K_yz, n_cholesky, size_AB)
