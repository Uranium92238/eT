submodule (cc2_class) excited_state
!
!!
!!    Excited state submodule (CC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Contains the following family of procedures of the CC2 class:
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
   module subroutine initialize_excited_states_cc2(wf)
!!
!!    Initialize excited states
!!    Written by Sarai D. Folkestad, June 2017
!!
!!    Calculates and sets n_s2am, and updates n_parameters
!!    for excited state calculation
!!
      implicit none 
!    
      class(cc2) :: wf
!
      wf%n_s2am = ((wf%n_v)*(wf%n_o))*((wf%n_v )*(wf%n_o)+1)/2 
!
      
      wf%n_parameters = wf%n_parameters + wf%n_s2am                
!
   end subroutine initialize_excited_states_cc2
!
!
   module subroutine calculate_orbital_differences_cc2(wf, orbital_diff)
!!
!!    Calculate Orbital Differences (CC2)
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
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
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
      do i = 1, wf%n_o
!
         do a = 1, wf%n_v
!
            ai = index_two(a, i, wf%n_v)
!
            do j = 1, wf%n_o
!
               do b = 1, wf%n_v
!
                  bj = index_two(b, j, wf%n_v)
!
                  aibj = index_packed(ai, bj)
!
                  orbital_diff((wf%n_o)*(wf%n_v) + aibj, 1) &
                                                 = wf%fock_diagonal(wf%n_o + a , 1) &
                                                 - wf%fock_diagonal(i, 1) &
                                                 + wf%fock_diagonal(wf%n_o + b, 1) &
                                                 - wf%fock_diagonal(j, 1)
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine calculate_orbital_differences_cc2
!
!
   module subroutine transform_trial_vectors_cc2(wf, first_trial, last_trial)
!!
!!    Transformation of Trial Vectors (CC2)
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
      class(cc2) :: wf
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
         if (wf%excited_state_task=='right_valence') then
!
               call wf%jacobian_cc2_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='right_core') then
!
               call wf%cvs_jacobian_cc2_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='left_valence') then
!
      !         call wf%jacobian_transpose_cc2_transformation(c_a_i, c_aibj)
!
            else
!
               write(unit_output,*) 'Error: Excited state task not recognized'
               stop
!
         endif
! 
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
   end subroutine transform_trial_vectors_cc2
!
!
   module subroutine cvs_residual_projection_cc2(wf, residual)
!!
!!    Residual projection (CC2), 
!!    Written by Sarai D. Folkestad Aug. 2017    
!!
      implicit none
!
      class(cc2) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital

      do i = 1, wf%n_o
!
         core_orbital = .false.
         do core = 1, wf%tasks%n_cores
!
            if (i .eq. wf%tasks%index_core_mo(core, 1)) core_orbital = .true.
!
         enddo
!
         if (.not. core_orbital) then
            do a = 1, wf%n_v
               ai = index_two(a, i, wf%n_v)
               residual(ai, 1) = zero
            enddo
         endif
!
      enddo
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            core_orbital = .false.
            do core = 1, wf%tasks%n_cores
!
               if ((i .eq. wf%tasks%index_core_mo(core, 1)) .or. &
                  (j .eq. wf%tasks%index_core_mo(core, 1))) core_orbital = .true.
!
            enddo
!
            if (.not. core_orbital) then
               do a = 1, wf%n_v
                  do b = 1, wf%n_v
                     ai = index_two(a, i, wf%n_v)
                     bj = index_two(b, j, wf%n_v)
                     aibj = index_packed(ai, bj)

                     residual(wf%n_t1am + aibj, 1) = zero
                  enddo
               enddo
            endif
         enddo
      enddo
!
    end subroutine cvs_residual_projection_cc2
!
!
end submodule