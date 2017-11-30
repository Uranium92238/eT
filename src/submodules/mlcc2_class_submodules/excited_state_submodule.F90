submodule (mlcc2_class) excited_state
!
!!
!!    Excited state submodule (MLCC2) 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Contains the following family of procedures of the MLCC2 class:
!!
!!    inititialize_excited_states: Initializes number of s2 amplitudes (n_x2am), and adds it n_parameters
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
   module subroutine excited_state_preparations_mlcc2(wf)
!!
!!    Excited State Preparations (MLCC2)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      class(mlcc2) :: wf 
!
!     Do nothing for mlcc2
!
   end subroutine excited_state_preparations_mlcc2
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
      wf%n_x2am = ((wf%n_CC2_v)*(wf%n_CC2_o))&
                   *((wf%n_CC2_v )*(wf%n_CC2_o)+1)/2 
!
      
      wf%n_parameters = wf%n_parameters + wf%n_x2am
                       
!
   end subroutine initialize_excited_states_mlcc2
!
!
   module subroutine calculate_orbital_differences_mlcc2(wf, orbital_diff)
!!
!!    Calculate Orbital Differences (MLCC2)
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
!     Calculate active space indices
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
         if (wf%excited_state_task=='right_valence' .or. wf%excited_state_task=='right_core') then
!
               call wf%jacobian_mlcc2_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_task=='left_valence') then
!
               write(unit_output,*)'Error: Jacobian transpose not implemented for mlcc2'
               stop
!
            else
!
               write(unit_output,*) 'Error: Excited state task not recognized'
               stop
!
         endif
!
!        -::- Projections -::-
!
!        Test for core calculation 
!
         if (wf%tasks%core_excited_state .or. wf%tasks%core_ionized_state) then
!  
!           Project out contamination from valence contributions
!
            call wf%cvs_rho_a_i_projection(c_a_i)
            call wf%cvs_rho_aibj_projection(c_aibj)
!
         endif
!
!        Test for ionization calculation
!
         if (wf%tasks%ionized_state .or. wf%tasks%core_ionized_state) then
!
!           Project out contamination from regular excitations
!
            write(unit_output,*)'Error: Ionized state not implemented for mlcc2'
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
   end subroutine transform_trial_vectors_mlcc2
!
!
   module subroutine cvs_residual_projection_mlcc2(wf, residual)
!!
!!    Residual projection (MLCC2), 
!!    Written by Sarai D. Folkestad Aug. 2017    
!!
      implicit none
!
      class(mlcc2) :: wf
      real(dp), dimension(wf%n_parameters, 1) :: residual
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, core = 0, ai = 0, bj = 0, aibj = 0
!
      logical :: core_orbital
!
!     Active space variables
!     
      integer(i15) :: n_active_o = 0, n_active_v = 0
      integer(i15) :: first_active_o ! first active occupied index 
      integer(i15) :: first_active_v ! first active virtual index
      integer(i15) :: last_active_o ! first active occupied index 
      integer(i15) :: last_active_v ! first active virtual index
!   
!     Calculate first/last indeces
!     
      call wf%get_CC2_active_indices(first_active_o, first_active_v)
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!   
      last_active_o = first_active_o + n_active_o - 1
      last_active_v = first_active_v + n_active_v - 1   
      do i = 1, wf%n_o
!
         core_orbital = .false.
         do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
            if (i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) core_orbital = .true.
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
      do i = 1, n_active_o
         do j = 1, n_active_o
!
            core_orbital = .false.
            do core = 1, wf%core_excited_state_specifications%n_equivalent_cores
!
               if ((i .eq. wf%core_excited_state_specifications%index_core_mo(core, 1)) .or. &
                  (j .eq. wf%core_excited_state_specifications%index_core_mo(core, 1))) core_orbital = .true.
!
            enddo
!
            if (.not. core_orbital) then
               do a = 1, n_active_v
                  do b = 1, n_active_v
                     ai = index_two(a, i, n_active_v)
                     bj = index_two(b, j, n_active_v)
                     aibj = index_packed(ai, bj)

                     residual(wf%n_t1am + aibj, 1) = zero
                  enddo
               enddo
            endif
         enddo
      enddo
!
    end subroutine cvs_residual_projection_mlcc2
!
!
      module subroutine print_excitation_vector_mlcc2(wf, vec, unit_id)
!!
!!
!!
         implicit none
!  
         class(mlcc2) :: wf
!
         real(dp), dimension(wf%n_parameters, 1) :: vec
!
         integer(i15) :: unit_id     
!
         integer(i15) :: a = 0, i = 0, ai = 0, b = 0, j = 0, bj = 0, aibj = 0
         integer(i15) :: n_active_o, n_active_v
         real(dp)     :: a_active_i_active, a_active_i_inactive, a_inactive_i_active, a_inactive_i_inactive, total
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
      call wf%get_CC2_n_active(n_active_o, n_active_v)
!
      do a = 1, n_active_v
         do i = 1, n_active_o
            do b = 1, n_active_v
               do j = 1, n_active_o
!
                  ai = index_two(a, i, n_active_v)
                  bj = index_two(b, j, n_active_v)
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
!    
!
      a_active_i_active       = 0
      a_inactive_i_active     = 0
      a_active_i_inactive     = 0
      a_inactive_i_inactive   = 0
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = index_two(a, i, wf%n_v)
!
            if (a .le. n_active_v .and. i .le. n_active_o) then
!
               a_active_i_active = a_active_i_active + (vec(ai, 1))**2
!
            elseif (a .le. n_active_v .and. i .gt. n_active_o) then
!
               a_active_i_inactive = a_active_i_inactive + (vec(ai, 1))**2
!
            elseif (a .gt. n_active_v .and. i .le. n_active_o) then
!
               a_inactive_i_active = a_inactive_i_active + (vec(ai, 1))**2
!
            elseif (a .gt. n_active_v .and. i .gt. n_active_o) then
!
               a_inactive_i_inactive = a_inactive_i_inactive + (vec(ai, 1))**2
!
            endif

         enddo
      enddo
!
!     Print active space stats:
!
      total = (a_active_i_active + a_active_i_inactive + a_inactive_i_active + a_inactive_i_inactive)
      write(unit_id,'(/a10, 3a12)')'T->T:', 'S->T:', 'T->S:', 'S->S:'
      write(unit_id,'(a50)')'--------------------------------------------------'
      write(unit_id,'(4f12.5)') a_active_i_active/total,&
         a_active_i_inactive/total, &
         a_inactive_i_active/total,&
         a_inactive_i_inactive/total
      write(unit_id,'(/t3,a40,f12.5/)')'Singles contribution to excitation is', total
   end subroutine print_excitation_vector_mlcc2
!
end submodule