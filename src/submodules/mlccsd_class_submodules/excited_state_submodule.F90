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
   module subroutine initialize_excited_states_mlccsd(wf)
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
!
   module subroutine transform_trial_vectors_mlccsd(wf, first_trial, last_trial)
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
      class(mlccsd) :: wf
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
      call allocator(c_aibj, wf%n_x2am, 1)
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
         if (wf%excited_state_task=='right_valence') then
!
               call wf%jacobian_mlccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='right_core') then
!
               call wf%cvs_jacobian_mlccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='left_valence') then
!
      !         call wf%jacobian_transpose_mlccsd_transformation(c_a_i, c_aibj)
!
            else
!
               write(unit_output,*) 'Error: Excited state task not recognized'
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
      call deallocator(c_aibj, wf%n_x2am, 1)
!
   end subroutine transform_trial_vectors_mlccsd
!
!
end submodule excited_state