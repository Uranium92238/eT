submodule (ccsd_class) excited_state
!
!!
!!    Excited state  submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    calculate_orbital_differences: calculates the orbital energy differences.
!!    transform_trial_vectors:       transform trial vectors by the Jacobian (or Jacobian^T).
!!
!
contains 
!
!
   module subroutine calculate_orbital_differences_ccsd(wf,orbital_diff)
!!
!!       Calculate Orbital Differences (CCSD)
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!       Calculates orbital differences
!!
!!          1) ε_i^a = ε_a - ε_i
!!          2) ε_ij^ab = ε_a + ε_b - ε_i - ε_j
!!
!!       and puts them in orbital_diff, which is a vector of length n_parameters.        
!!
         implicit none
!
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
         integer(i15) :: a = 0, i = 0, b = 0, j = 0
         integer(i15) :: ai = 0, bj = 0
         integer(i15) :: aibj = 0
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
!
               ai = index_two(a, i, wf%n_v)
!
               orbital_diff(ai, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1)
!
               do j = 1, wf%n_o
                  do b = 1, wf%n_v
!
                     bj = index_two(b, j, wf%n_v)
!
                     aibj = index_packed(ai, bj)
!
                     orbital_diff((wf%n_o)*(wf%n_v)+aibj, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1) &
                                                               + wf%fock_diagonal(b + wf%n_o, 1) - wf%fock_diagonal(j, 1)

!
                  enddo
               enddo
            enddo
         enddo
!
   end subroutine calculate_orbital_differences_ccsd
!
!
   module subroutine transform_trial_vectors_ccsd(wf, first_trial, last_trial)
!!
!!    Transformation Trial Vectors (CCSD)
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
      class(ccsd) :: wf
!
      integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
      real(dp), dimension(:,:), allocatable :: c_a_i
      real(dp), dimension(:,:), allocatable :: c_aibj
!
      integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
      integer(i15) :: trial = 0 
!
!     Allocate c_a_i and c_aibj
!
      call allocator(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero 
!
      call allocator(c_aibj, wf%n_t2am, 1)
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
         if (wf%current_task == 'excited_state') then
!
            if (wf%excited_state_task =='right_valence') then
!
               call wf%jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='right_core') then
!
               !call wf%cvs_jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='left_valence') then
!               
               call wf%jacobian_transpose_ccsd_transformation(c_a_i, c_aibj)
!
            else
!
               write(unit_output,*) 'Error: Excited state task not recognized'
               stop
!
            endif
!
         elseif (wf%current_task == 'response') then
!
            if (wf%response_task == 'left_eigenvectors') then
!
               call wf%jacobian_transpose_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%response_task == 'multipliers') then 
!
               call wf%jacobian_transpose_ccsd_transformation(c_a_i, c_aibj)
!
            else
!
               write(unit_output,*) 'Error: Response task not recognized'
               stop
!
            endif
!
         else
!
            write(unit_output,*) 'Error: Current task not recognized'
            stop
!
         endif
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
      call deallocator(c_aibj, wf%n_t2am, 1)
!
   end subroutine transform_trial_vectors_ccsd
!
!
   module subroutine excited_state_preparations_ccsd(wf)
!!
!!    Excited State Preparations (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      class(ccsd) :: wf 
!
!     Store vvvv-electronic repulsion integrals to file if there is space 
!
      call wf%store_t1_vv_vv_electronic_repulsion
!
!     Store voov-electronic repulsion integrals to file if there is space
!
      call wf%store_t1_vo_ov_electronic_repulsion
!
!     Store vvvo-electronic repulsion integrals to file if there is space
!
      call wf%store_t1_vv_vo_electronic_repulsion
!
   end subroutine excited_state_preparations_ccsd
!
!
end submodule excited_state