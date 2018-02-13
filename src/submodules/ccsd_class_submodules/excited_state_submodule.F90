submodule (ccsd_class) excited_state
!
!!
!!    Excited state submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    calculate_orbital_differences:    calculates the orbital energy differences.
!!    transform_trial_vectors:          transform trial vectors by the Jacobian (or Jacobian^T).
!!
!!    print_excitation_vector:          description missing
!!    analyze_double_excitation_vector: -/-
!!    summary_excited_state_info:       -/-
!!    
!
contains 
!
!
   module subroutine calculate_orbital_differences_ccsd(wf,orbital_diff)
!!
!!    Calculate orbital differences (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates orbital differences
!!
!!          1) ε_i^a = ε_a - ε_i
!!          2) ε_ij^ab = ε_a + ε_b - ε_i - ε_j
!!
!!    and puts them in orbital_diff, which is a vector of length n_parameters.        
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
!!    Transformation trial vectors (CCSD)
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
      call wf%mem%alloc(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero 
!
      call wf%mem%alloc(c_aibj, wf%n_t2am, 1)
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
!        Test for left or right transformation
!
         if (wf%tasks%current == 'excited_state') then
!
            if (wf%excited_state_specifications%right) then
!
               call wf%jacobian_ccsd_transformation(c_a_i, c_aibj)
!
            elseif (wf%excited_state_specifications%left) then
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
         elseif (wf%tasks%current == 'multipliers') then
!
            call wf%jacobian_transpose_ccsd_transformation(c_a_i, c_aibj)
!
         else
!
            write(unit_output,*) 'Error: Current task not recognized'
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
            call wf%ionization_rho_a_i_projection(c_a_i)
            call wf%ionization_rho_aibj_projection(c_aibj)
!
         endif
!
!        Write transformed vector to file
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
      call wf%mem%dealloc(c_a_i, wf%n_v, wf%n_o)
      call wf%mem%dealloc(c_aibj, wf%n_t2am, 1)
!
   end subroutine transform_trial_vectors_ccsd
!
!
   module subroutine print_excitation_vector_ccsd(wf, vec, unit_id)
!!
!!    Print excitation vector 
!!    Written by Sarai D. Folkestad, Oct 2017
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
!     Print singles part
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
   module subroutine excited_state_preparations_ccsd(wf)
!!
!!    Excited state preparations (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      class(ccsd) :: wf 
!
      if (wf%settings%print_level == 'developer') write(unit_output,'(/t3,a/)') 'Preparing for excited state calculation:'
!
!     Read single amplitudes 
!
      call wf%read_single_amplitudes
!
!     Store vvvv-electronic repulsion integrals to file if there is space 
!
      call wf%store_t1_vv_vv_electronic_repulsion
!
!     Store voov-electronic repulsion integrals to file if there is space
!
      call wf%store_t1_vo_ov_electronic_repulsion
!
!     Store vvvo and vvov-electronic repulsion integrals to file if there is space
!
      call wf%store_t1_vv_vo_electronic_repulsion
      call wf%store_t1_vv_ov_electronic_repulsion
!
!     Set current task to excited state calculation 
! 
      wf%tasks%current = 'excited_state'
!
!     Set filename for solution vectors
!
      if (wf%tasks%core_excited_state .or. wf%tasks%core_ionized_state) then   ! Core excitation
!
         if (wf%excited_state_specifications%right) then ! Right vectors
!
            wf%excited_state_specifications%solution_file = 'right_core'
!
         else ! Left vectors 
!                             
            wf%excited_state_specifications%solution_file = 'left_core'
!
         endif
!
      else ! Valence excitation
!
         if (wf%excited_state_specifications%left) then ! Left vectors
!
            wf%excited_state_specifications%solution_file = 'left_valence'
!
         else ! Right vectors
!
            wf%excited_state_specifications%solution_file = 'right_valence'
!
         endif
!
      endif
!
   end subroutine excited_state_preparations_ccsd
!
!
   module subroutine analyze_double_excitation_vector_ccsd(wf, vec, n, sorted_short_vec, index_list)
!!
!!    Analyze double excitation vector 
!!    Written by Sarai D. Folkestad, Oct 2017
!!
!!    Sorts the double excitation part of the excited state according to size,
!!    in order for the print routine to show the largest contributions 
!!    in the printout information.
!!
      implicit none
!  
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_t2am, 1) :: vec    
!
      integer(i15) :: a = 0, i = 0, ai = 0, b = 0, j = 0, bj = 0, aibj = 0, k = 0
!
      integer(i15) :: n    ! Number of elements wanted
!  
      real(dp), dimension(n, 1)    :: sorted_short_vec
!  
      integer(i15), dimension(n, 4) ::index_list
!  
!     Variables for sorting
!  
      real(dp)     :: min
      integer(i15) :: min_pos
!  
      real(dp)     :: swap     = zero
      integer(i15) :: swap_i = 0, swap_a = 0, swap_j = 0, swap_b = 0
!
!     Placing the n first elements of vec into sorted_short_vec
!
      index_list = 0 
      sorted_short_vec(1,1) = vec(1,1)
      index_list(1,1) = 1
      index_list(1,2) = 1
      index_list(1,3) = 1
      index_list(1,4) = 1
!
      min = abs(sorted_short_vec(1,1))
      min_pos = 1
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  ai = index_two(a,i, wf%n_v)
                  bj = index_two(b,j, wf%n_v)
!
                  if (ai .ge. bj) then
!
                     aibj = index_packed(ai,bj)
!
                     if (aibj .le. n) then
                        sorted_short_vec(aibj,1)   = vec(aibj,1)
                        index_list(min_pos,1)      = a
                        index_list(min_pos,2)      = i
                        index_list(min_pos,3)      = b
                        index_list(min_pos,4)      = j
!
                        if (abs(sorted_short_vec(i,1)) .le. min) then
!
                           min = abs(sorted_short_vec(i,1))
                           min_pos = i
!
                        endif
                     else
!
                     if (abs(vec(aibj,1)) .ge. min) then
!
                        sorted_short_vec(min_pos,1) = vec(aibj,1)
                        index_list(min_pos,1) = a
                        index_list(min_pos,2) = i
                        index_list(min_pos,3) = b
                        index_list(min_pos,4) = j
                        min = abs(vec(aibj,1))
!
                     endif
                  endif
!
                  do k = 1, n
                     if (abs(sorted_short_vec(k, 1)) .lt. min) then
!
                        min = abs(sorted_short_vec(k, 1))
                        min_pos = k
!
                     endif
                  enddo
               endif
!
               enddo
            enddo
         enddo
      enddo  
! 
!      Sorting sorted_short_vec
! 
       do i = 1, n
          do j = 1, n - 1
             if (abs(sorted_short_vec(j,1)) .lt. abs(sorted_short_vec(j+1, 1))) then
! 
                swap = sorted_short_vec(j,1)
                sorted_short_vec(j,1) = sorted_short_vec(j+1, 1)
                sorted_short_vec(j+1, 1) = swap
! 
                swap_a = index_list(j, 1)
                swap_i = index_list(j, 2)
                swap_b = index_list(j, 3)
                swap_j = index_list(j, 4)
!
                index_list(j,1) = index_list(j + 1,1)
                index_list(j,2) = index_list(j + 1,2)
                index_list(j,3) = index_list(j + 1,3)
                index_list(j,4) = index_list(j + 1,4)
                index_list(j + 1,1) = swap_a
                index_list(j + 1,2) = swap_i
                index_list(j + 1,3) = swap_b
                index_list(j + 1,4) = swap_j
! 
             endif
          enddo
       enddo
! 
!
   end subroutine analyze_double_excitation_vector_ccsd
!
!
   module subroutine summary_excited_state_info_ccsd(wf, energies)
!!
!!    Summary of excited state info 
!!    Written by Sarai D. Folkestad, Oct 2017
!!
!!    Prints the analysis of the excitation vectors to the main 
!!    output file.
!!
      implicit none
!  
      class(ccsd) :: wf
!
      real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: energies
!
      integer(i15) :: unit_solution = -1, ioerror = 0
      integer(i15) :: state = 0, i = 0
      real(dp), dimension(:,:), allocatable :: solution_ai, solution_aibj
      real(dp), dimension(:,:), allocatable :: sorted_max_vec_singles, sorted_max_vec_doubles
      integer(i15), dimension(:,:), allocatable :: index_list_singles, index_list_doubles
      real(dp) :: norm, ddot
!  
!     Open solution vector file
!  
      call generate_unit_identifier(unit_solution)
!
      open(unit=unit_solution, file=wf%excited_state_specifications%solution_file, &
            action='read', status='unknown', &
            access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
      if (ioerror .ne. 0) write(unit_output,*) 'Error while opening solution file'
!
!     Allocations
!
      call wf%mem%alloc(solution_ai, wf%n_t1am, 1)
      call wf%mem%alloc(solution_aibj, wf%n_t2am, 1)
!
      call wf%mem%alloc(sorted_max_vec_singles, min(wf%n_t1am, 20), 1)
      call wf%mem%alloc(sorted_max_vec_doubles, min(wf%n_t2am, 20), 1)
!
      call wf%mem%alloc_int(index_list_singles, min(wf%n_t1am, 20), 2)
      call wf%mem%alloc_int(index_list_doubles, min(wf%n_t2am, 20), 4)
!
      do state = 1, wf%excited_state_specifications%n_singlet_states
!
         write(unit_output,'(/t3,a30,i3,a1/)')'Analysis of excitation vector ',state, ':' 
         write(unit_output,'(t6, a, f14.8)')'Excitation energy [a.u.]:   ', energies(state,1)
         write(unit_output,'(t6, a, f14.8)')'Excited state energy [a.u.]:', wf%energy + energies(state,1)
! 
!        Read the solution
!
         solution_ai    = zero
         solution_aibj  = zero
         read(unit_solution, rec=state) solution_ai, solution_aibj
!
!        Calculate the contribution from single excitations
!
         norm = sqrt(ddot(wf%n_t1am, solution_ai,1,solution_ai,1))
         write(unit_output,'(/t6,a,f6.4)')'Single excitation contribution to excitation vector: ', norm
!
!        Analysis of excitation vectors
!
         write(unit_output,'(/t6,a)') 'Largest contributions to excitation vector:' 
!
         write(unit_output,'(t6,a32)')'------------------------------------------------------'
         write(unit_output,'(t6,a3, 8x, a3, 8x, a10)')'a', 'i', 'amplitude'
         write(unit_output,'(t6,a32)')'------------------------------------------------------'
!
!        Get 20 highest amplitudes
!
         call wf%analyze_single_excitation_vector(solution_ai, min(wf%n_t1am, 20), &
                         sorted_max_vec_singles, index_list_singles)
         call wf%analyze_double_excitation_vector(solution_aibj, min(wf%n_t2am, 20), &
                         sorted_max_vec_doubles, index_list_doubles)
!
!        And print them
!
         do i = 1, min(wf%n_t1am, 20)
!
            if    (abs(sorted_max_vec_singles(i, 1)) .lt. 1.0D-03) then
!
               exit 
!
            else
!
               write(unit_output,'(t6,i3, 8x,i3, 10x, f8.5)') &
                                                               index_list_singles(i, 1),&
                                                               index_list_singles(i, 2),&
                                                               sorted_max_vec_singles(i, 1) 
            endif
         enddo
!
         write(unit_output,'(t6,a32)')'------------------------------------------------------'
         write(unit_output,'(t6,a54)')'------------------------------------------------------'
         write(unit_output,'(t6,a3, 8x, a3, 8x, a3, 8x, a3, 8x, a10)')'a','i','b','j', 'amplitude'
         write(unit_output,'(t6,a54)')'------------------------------------------------------'
!
         do i = 1, min(wf%n_t2am, 20)
!
            if    (abs(sorted_max_vec_doubles(i, 1)) .lt. 1.0D-03) then
!
               exit 
!
            else
!
               write(unit_output,'(t6,i3, 8x,i3, 8x,i3, 8x, i3, 10x, f8.5)')&
                                                               index_list_doubles(i, 1),&
                                                               index_list_doubles(i, 2),&
                                                               index_list_doubles(i, 3),&
                                                               index_list_doubles(i, 4),&
                                                               sorted_max_vec_doubles(i, 1) 
            endif
         enddo
         write(unit_output,'(t6,a54)')'------------------------------------------------------'
!  
      enddo
!
!     Deallocations
!
      call wf%mem%dealloc(solution_ai, wf%n_t1am, 1)
      call wf%mem%dealloc(solution_aibj, wf%n_t2am, 1)
!
      call wf%mem%dealloc(sorted_max_vec_singles, min(wf%n_t1am, 20), 1)
      call wf%mem%dealloc(sorted_max_vec_doubles, min(wf%n_t2am, 20), 1)
!
      call wf%mem%dealloc_int(index_list_singles, min(wf%n_t1am, 20), 2)
      call wf%mem%dealloc_int(index_list_doubles, min(wf%n_t2am, 20), 4)
!
      close(unit_solution)
!
   end subroutine summary_excited_state_info_ccsd
!
!
end submodule excited_state