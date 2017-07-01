submodule (mlcc2_class) excited_state
!
!!
!!    Excited state submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!    inititialize_excited_states: Initializes number of s2 amplitudes (n_s2am), and adds it n_parameters
!!    calculate_orbital_differences: Calculates the orbital differences, including the double excitation differences 
!!                                   in the active CC2 spaces
!!    transform_trial_vectors: Transforms the new trial vectors. rho = Ac
!!
!!    Upper case indices are general indices, lower case indices are restricted
!!    to the CC2 orbital space.
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
!
   module subroutine initialize_excited_states_mlcc2(wf)
!!
!!    Initialize excited states
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Calculates and sets n_s2am, and updates n_parameters
!!    for excited state calculation
!!
      implicit none 
!    
      class(mlcc2) :: wf
!
      integer(i15) :: active_space
!     
!     Add packed number of double amplitudes 
!
      wf%n_total_active_o = 0
      wf%n_total_active_v = 0
!
      wf%n_s2am = 0
!
      do active_space = 1, wf%n_active_spaces
!
         wf%n_s2am = wf%n_s2am &
                  + ((wf%n_CC2_v(active_space,1))*(wf%n_CC2_o(active_space,1)))&
                   *((wf%n_CC2_v(active_space,1) )*(wf%n_CC2_o(active_space,1))+1)/2 
!
      enddo
!
      
      wf%n_parameters = wf%n_parameters + wf%n_s2am
                       
!
   end subroutine initialize_excited_states_mlcc2
!
!
   module subroutine calculate_orbital_differences_mlcc2(wf, orbital_diff)
!!
!!    Calculate Orbital Differences (CCSD)
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
!     Adding double terms only for active spaces
!
      offset = 0
      do active_space = 1, wf%n_active_spaces
!
!        Calculate first/last indices
! 
         call wf%get_CC2_active_indices(first_active_o, first_active_v, active_space)
!
         n_active_o = wf%n_CC2_o(active_space,1) 
         n_active_v = wf%n_CC2_v(active_space,1)
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
                     orbital_diff((wf%n_o)*(wf%n_v) + aibj + offset, 1) &
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
         offset = offset &
               + (n_active_o*n_active_v)*(n_active_o*n_active_v+1)/2
!
      enddo
!
   end subroutine calculate_orbital_differences_mlcc2
!
!
   module subroutine transform_trial_vectors_mlcc2(wf, first_trial, last_trial)
!!
!!    Transformation of Trial Vectors (MLCC2)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Each trial vector in first_trial to last_trial is read from file and
!!    transformed before the transformed vector is written to file.
!!
!!    Singles and doubles part of the transformed vectors are written to 
!!    the same record in file transformed_vec, record length is n_parameters long.
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      real(dp), dimension(:,:), allocatable :: c_a_i
      real(dp), dimension(:,:), allocatable :: c_aibj
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
      integer(i15) :: trial = 0 
!
!
!     Allocate c_a_i and c_aibj
!
      call allocator(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero 
!
      call allocator(c_aibj, wf%n_s2am, 1)
      c_aibj = zero 
!
!     Open trial vector- and transformed vector files
!
      call generate_unit_identifier(unit_trial_vecs)
      open(unit=unit_trial_vecs, file='trial_vec', action='read', status='unknown', &
           access='direct', form='unformatted', recl=dp*wf%n_parameters, iostat=ioerror)
!
      call generate_unit_identifier(unit_rho)
      open(unit=unit_rho, file='transformed_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*wf%n_parameters, iostat=ioerror)
!
!     For each trial vector: read, transform and write  
!  
      do trial = first_trial, last_trial
!
         read(unit_trial_vecs, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
         call wf%jacobian_mlcc2_transformation(c_a_i, c_aibj)
!
         write(unit_rho, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
      enddo
!
!     Close files
!
      close(unit_trial_vecs) 
      close(unit_rho)                                
!
!     Deallocate c_a_i and c_aibj
!
      call deallocator(c_a_i, wf%n_v, wf%n_o)
      call deallocator(c_aibj, wf%n_s2am, 1)
!
   end subroutine transform_trial_vectors_mlcc2
!
!
end submodule