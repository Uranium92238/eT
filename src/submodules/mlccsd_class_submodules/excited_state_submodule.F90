submodule (mlccsd_class) excited_state
!
!!
!!    Excited state submodule (MLCCSD) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Aug 2017
!!
!! 
!
   implicit none 
!
   logical :: debug   = .false.
   logical :: timings = .false.
!
!
contains
!
   subroutine initialize_excited_states_mlccsd(wf)
!!
!!    Initialize excited states
!!    Written by Sarai D. Folkestad, Aug 2017
!!
!!    Calculates and sets n_s2am, and updates n_parameters
!!    for excited state calculation
!!
      implicit none 
!    
      class(mlccsd) :: wf
!
      integer(i15) :: n_o
      integer(i15) :: n_v
!     
!     Add packed number of double amplitudes 
!
      n_o = wf%n_CC2_o + wf%n_CCSD_o
      n_v = wf%n_CC2_v + wf%n_CCSD_v
!
      wf%n_x2am = (n_o*n_v)*(n_o*n_v + 1)/2
!
      wf%n_parameters = wf%n_t1am + wf%n_x2am                     
!
   end subroutine initialize_excited_states_mlccsd
!
   subroutine calculate_orbital_differences_mlccsd(wf, orbital_diff)
!!
!!    Calculate Orbital Differences (MLCCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
!!    Calculates orbital differences
!!
!!       1) ε_I^A = ε_A - ε_I
!!       2) ε_ij^ab = ε_a + ε_b - ε_i - ε_j (for active spaces only)
!!
!!    and puts them in orbital_diff, which is a vector of length n_parameters.        
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
!     Active space variables
!
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! last active occupied index 
      integer(i15) :: last_active_v ! last active virtual index
      integer(i15) :: n_active_o
      integer(i15) :: n_active_v         
!
      integer(i15) :: active_space
!
      integer(i15) :: offset = 0
!
      integer(i15) :: A = 0, I = 0, b = 0, j = 0
      integer(i15) :: AI = 0, bj = 0
      integer(i15) :: aibj = 0
!
      do I = 1, wf%n_o
         do A = 1, wf%n_v
!
            AI = index_two(A, I, wf%n_v)
!
            orbital_diff(AI, 1) = wf%fock_diagonal(A + wf%n_o, 1) - wf%fock_diagonal(I, 1)
!
         enddo
      enddo
!
!     Calculate first/last indices
! 
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      do i = 1, n_active_o
!
         do a = 1, n_active_v
!
            ai = index_two(a, i, n_active_v)
!
            do j = 1, n_active_o
!
               do b = 1, n_active_v
!
                  bj = index_two(b, j, n_active_v)
!
                  aibj = index_packed(ai, bj)
!
                  orbital_diff((wf%n_o)*(wf%n_v) + aibj, 1) &
                                                 = wf%fock_diagonal(wf%n_o + a + first_active_v - 1, 1) &
                                                 - wf%fock_diagonal(i + first_active_o - 1, 1) &
                                                 + wf%fock_diagonal(wf%n_o + b + first_active_v - 1, 1) &
                                                 - wf%fock_diagonal(j + first_active_o - 1, 1)
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine calculate_orbital_differences_mlccsd
!
end submodule excited_state