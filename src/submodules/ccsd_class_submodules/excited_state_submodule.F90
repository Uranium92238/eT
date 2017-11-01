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
!        Test for excited state or ionized state
!
         if (wf%tasks%excited_state .or. wf%tasks%core_excited_state ) then ! Excited state
!
!           Left, right core or valence
!
            if (wf%excited_state_task=='right_eigenvectors') then
!
               call wf%jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='left_eigenvectors') then
!
               call wf%jacobian_transpose_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%response_task=='multipliers') then 
!
               call wf%jacobian_transpose_ccsd_transformation(c_a_i, c_aibj)
!
             elseif (wf%excited_state_task=='right_core') then
               write(unit_output,*)'core transformation'
!
               call wf%cvs_jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            else
!
               write(unit_output,*) 'Error: Response task not recognized'
               stop
!
            endif
!
            write(unit_rho, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
      elseif (wf%tasks%ionized_state  .or. wf%tasks%core_ionized_state ) then ! Ionized state
!
            if (wf%excited_state_task=='right_valence') then
!
               call wf%ionization_jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='right_core') then
!
               !call wf%ionized_cvs_jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            else
!
               write(unit_output,*) 'Error: Response task not recognized'
               stop
!
            endif
!
            write(unit_rho, rec=trial, iostat=ioerror) c_a_i, c_aibj
!
         endif
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
      module subroutine print_excitation_vector_ccsd(wf, vec, unit_id)
!!
         implicit none
!  
         class(ccsd) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id     
!
         integer(i15) :: a = 0, i = 0, ai = 0, b = 0, j = 0, bj = 0, aibj = 0
!
!        Print singles part
!
         write(unit_id,'(2a6,a12)')'a', 'i', 'coeff'
         write(unit_id,'(t3,a)')'-------------------------'
!
         do a = 1, wf%n_v
            do i = 1, wf%n_o
!  
               ai = index_two(a, i, wf%n_v)
               if (abs(vec(ai, 1)) .gt. 1.0D-03) then
                  write(unit_id,'(2i6,f12.4)') a, i, vec(ai, 1)
               endif
!
         enddo
      enddo
      flush(unit_id)
!
!     Print doubles part
!
      write(unit_id,'(/4a6, a11)')'a', 'i', 'b', 'j', 'coeff'
      write(unit_id,'(t3,a)')'---------------------------------'
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  ai = index_two(a, i, wf%n_v)
                  bj = index_two(b, j, wf%n_v)
                  aibj = index_packed(ai, bj)
!
                  if (abs(vec((wf%n_o)*(wf%n_v) + aibj, 1)) .gt. 1.0D-03 .and. ai .ge. bj) then
                     write(unit_id,'(4i6,f12.4)') a, i, b, j, vec((wf%n_o)*(wf%n_v) + aibj, 1)
                  endif
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine print_excitation_vector_ccsd
!
!
   module subroutine cvs_residual_projection_ccsd(wf, residual)
!!
!!    Residual projection (CCSD), 
!!    Written by Sarai D. Folkestad Aug. 2017    
!!
      implicit none
!
      class(ccsd) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital
!
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
               do a = 1,  wf%n_v
                  do b = 1,  wf%n_v
                     ai = index_two(a, i,  wf%n_v)
                     bj = index_two(b, j,  wf%n_v)
                     aibj = index_packed(ai, bj)

                     residual(wf%n_t1am + aibj, 1) = zero
                  enddo
               enddo
            endif
         enddo
      enddo
!
    end subroutine cvs_residual_projection_ccsd
!
!
end submodule excited_state