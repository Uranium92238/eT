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
      module subroutine excited_state_preparations_mlccsd(wf)
!!
!!    Excited State Preparations (MLCCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for preparation tasks (if any). Can be overwritten
!!    in descendants if other preparations prove necessary.    
!!
      class(mlccsd) :: wf 
!
      integer(i15) :: n_o
      integer(i15) :: n_v
!
!     Set current task to excited state calculation 
! 
      wf%tasks%current = 'excited_state'
!     
!     Set n_parameters
!
      n_o = wf%n_CC2_o + wf%n_CCSD_o
      n_v = wf%n_CC2_v + wf%n_CCSD_v
!
      wf%n_x2am = (n_o*n_v)*(n_o*n_v + 1)/2
!
      wf%n_parameters = wf%n_t1am + wf%n_x2am                     
!
      call wf%initialize_single_amplitudes
      call wf%read_single_amplitudes
!
!     Set filename for solution vectors
!
      if (wf%tasks%core_excited_state .or. wf%tasks%core_ionized_state) then   ! Core excitation
!
         if (wf%excited_state_specifications%right) then                         ! Right vectors
            wf%excited_state_specifications%solution_file = 'right_core'
         else                                                                    ! Left vectors
            write(unit_output,*)'Error: Jacobian transpose transformation not implemented for core excitations' ! S: should be able to get these with the same projections however so...
            stop
         endif
!
      else                                                                    ! Valence excitation
!
         if (wf%excited_state_specifications%left) then                          ! Right vectors
            wf%excited_state_specifications%solution_file = 'left_valence'
         else                                                                    ! Left vectors
            wf%excited_state_specifications%solution_file = 'right_valence'
         endif
!
      endif
!
   end subroutine excited_state_preparations_mlccsd
!
!
   module subroutine excited_state_cleanup_mlccsd(wf)
!!
!!    Excited State Cleanup (MLCCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Oct 2017
!!
!!    A routine for cleanup tasks (if any). Can be overwritten
!!    in descendants if other cleanups prove necessary.    
!!
      implicit none
!
      class(mlccsd) :: wf 
!
!     Deallocate the amplitudes 
!
      call wf%destruct_single_amplitudes
      call wf%destruct_cc2_double_amplitudes
!
   end subroutine excited_state_cleanup_mlccsd
!
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
!     Allocate c_a_i and c_aibj
!
      call wf%mem%alloc(c_a_i, wf%n_v, wf%n_o)
      c_a_i = zero 
!
      call wf%mem%alloc(c_aibj, wf%n_x2am, 1)
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
         if (wf%excited_state_specifications%right) then
!
            call wf%jacobian_mlccsd_transformation(c_a_i, c_aibj)
!
         elseif (wf%excited_state_specifications%left) then
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
      call wf%mem%dealloc(c_a_i, wf%n_v, wf%n_o)
      call wf%mem%dealloc(c_aibj, wf%n_x2am, 1)
!
   end subroutine transform_trial_vectors_mlccsd
!
!
      module subroutine print_excitation_vector_mlccsd(wf, vec, unit_id)
!!
!!
!!
         implicit none
!  
         class(mlccsd) :: wf
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
      call wf%get_CCSD_n_active(n_active_o, n_active_v)
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
!
!
   end subroutine print_excitation_vector_mlccsd
!
   module subroutine summary_excited_state_info_mlccsd(wf, energies)
!!
!!
!!
      implicit none
!  
      class(mlccsd) :: wf
!
      real(dp), dimension(wf%excited_state_specifications%n_singlet_states,1) :: energies
!
      integer(i15) :: unit_solution = -1, ioerror = 0
      integer(i15) :: state = 0, i = 0, a  = 0, ai = 0
!
      real(dp), dimension(:,:), allocatable :: solution_ai, solution_aibj
      real(dp), dimension(:,:), allocatable :: sorted_max_vec_singles, sorted_max_vec_doubles
 !
      integer(i15), dimension(:,:), allocatable :: index_list_singles, index_list_doubles
!
      real(dp) :: norm, ddot
      real(dp) :: a_active_i_active, a_active_i_cc2_inactive, a_active_i_ccs_inactive 
      real(dp) :: a_cc2_inactive_i_active, a_ccs_inactive_i_active, a_cc2_inactive_i_cc2_inactive
      real(dp) :: a_cc2_inactive_i_ccs_inactive, a_ccs_inactive_i_ccs_inactive, a_ccs_inactive_i_cc2_inactive
      real(dp) :: total
!
      integer(i15) :: n_active_o, n_active_v, n_active_o_cc2, n_active_v_cc2
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
      call wf%mem%alloc(solution_aibj, wf%n_x2am, 1)
!
      call wf%mem%alloc(sorted_max_vec_singles, 20, 1)
      call wf%mem%alloc(sorted_max_vec_doubles, 20, 1)
!
      call wf%mem%alloc_int(index_list_singles, 20, 2)
      call wf%mem%alloc_int(index_list_doubles, 20, 4)
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
         call wf%analyze_single_excitation_vector(solution_ai, 20, sorted_max_vec_singles, index_list_singles)
         call wf%analyze_double_excitation_vector(solution_aibj, 20, sorted_max_vec_doubles, index_list_doubles)
!
!        And print them
!
         do i = 1, 20
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
         do i = 1, 20
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
!        MLCC Specific print
!
         
         a_active_i_active       = 0 ! T->T
!
         a_active_i_cc2_inactive = 0 ! T->S
         a_active_i_ccs_inactive = 0 ! T->R
!
         a_cc2_inactive_i_active = 0 ! S->T
         a_ccs_inactive_i_active = 0 ! R->T
!
         a_cc2_inactive_i_cc2_inactive   = 0 ! S->S
         a_cc2_inactive_i_ccs_inactive   = 0 ! S->R
!
         a_ccs_inactive_i_ccs_inactive   = 0 ! R->R 
         a_ccs_inactive_i_cc2_inactive   = 0 ! R->S
!
         call wf%get_CCSD_n_active(n_active_o, n_active_v)

         call wf%get_CC2_n_active(n_active_o_CC2, n_active_v_cc2)
!
         do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = index_two(a, i, wf%n_v)
!  
            if (a .le. n_active_v) then ! T->
!
               if (i .le. n_active_o) then ! ->T
!
                  a_active_i_active = a_active_i_active + (solution_ai(ai, 1))**2
!
               elseif (i .gt. n_active_o .and. i .le. n_active_o_cc2) then ! ->S
!
                  a_active_i_cc2_inactive = a_active_i_cc2_inactive + (solution_ai(ai, 1))**2
!
               else ! ->R
!
                  a_active_i_ccs_inactive = a_active_i_ccs_inactive + (solution_ai(ai, 1))**2
!
               endif
!
            elseif (a .gt. n_active_v .and. a .le. n_active_v_cc2) then ! S->
!
               if (i .le. n_active_o) then ! ->T
!
                  a_cc2_inactive_i_active = a_cc2_inactive_i_active + (solution_ai(ai, 1))**2
!
               elseif (i .gt. n_active_o .and. i .le. n_active_o_cc2) then ! ->S
!
                  a_cc2_inactive_i_cc2_inactive = a_cc2_inactive_i_cc2_inactive + (solution_ai(ai, 1))**2
!
               else ! ->R
!
                  a_cc2_inactive_i_ccs_inactive = a_cc2_inactive_i_ccs_inactive + (solution_ai(ai, 1))**2
!
               endif
!
            else ! R->
!
               if (i .le. n_active_o) then ! ->T
!
                  a_ccs_inactive_i_active = a_ccs_inactive_i_active + (solution_ai(ai, 1))**2
!
               elseif (i .gt. n_active_o .and. i .le. n_active_o_cc2) then ! ->S
!
                  a_ccs_inactive_i_cc2_inactive = a_ccs_inactive_i_cc2_inactive + (solution_ai(ai, 1))**2
!
               else ! ->R
!
                  a_ccs_inactive_i_ccs_inactive = a_ccs_inactive_i_ccs_inactive + (solution_ai(ai, 1))**2
!
               endif
            endif
!
         enddo
      enddo
!
!     Print active space stats:
!
      total = (a_active_i_active + a_active_i_cc2_inactive + a_active_i_ccs_inactive + a_cc2_inactive_i_active &
               + a_ccs_inactive_i_active + a_cc2_inactive_i_cc2_inactive + a_cc2_inactive_i_ccs_inactive &
               + a_ccs_inactive_i_ccs_inactive + a_ccs_inactive_i_cc2_inactive)
!
      write(unit_output,'(/t6, a7, 8a9)')'T->T:', 'T->S:', 'T->R:', 'S->T:', 'S->S:',&
       'S->R', 'R->T:', 'R->S:', 'R->R:'
      write(unit_output,'(t6, a83)')'-----------------------------------------------------------------------------------'
      write(unit_output,'(t6, 9f9.5)')a_active_i_active/total, a_active_i_cc2_inactive/total, &
                                     a_active_i_ccs_inactive/total, a_cc2_inactive_i_active/total, &
                                     a_ccs_inactive_i_active/total, a_cc2_inactive_i_cc2_inactive/total, &
                                     a_cc2_inactive_i_ccs_inactive/total, &
                                     a_ccs_inactive_i_ccs_inactive/total , a_ccs_inactive_i_cc2_inactive/total
      write(unit_output,'(t6, a83)')'-----------------------------------------------------------------------------------'

      enddo
!
!     Deallocations
!
      call wf%mem%dealloc(solution_ai, wf%n_t1am, 1)
      call wf%mem%dealloc(solution_aibj, wf%n_x2am, 1)
!
      call wf%mem%dealloc(sorted_max_vec_singles, 20, 1)
      call wf%mem%dealloc(sorted_max_vec_doubles, 20, 1)
!
      call wf%mem%dealloc_int(index_list_singles, 20, 2)
      call wf%mem%dealloc_int(index_list_doubles, 20, 4)
!
      close(unit_solution)
!
   end subroutine summary_excited_state_info_mlccsd
!
end submodule excited_state
!            if (a .le. n_active_v .and. i .le. n_active_o) then
!!
!               a_active_i_active = a_active_i_active + (solution_ai(ai, 1))**2
!!
!            elseif (a .le. n_active_v .and. i .gt. n_active_o) then
!!
!               a_active_i_inactive = a_active_i_inactive + (solution_ai(ai, 1))**2
!!
!            elseif (a .gt. n_active_v .and. i .le. n_active_o) then
!!
!               a_inactive_i_active = a_inactive_i_active + (solution_ai(ai, 1))**2
!!
!            elseif (a .gt. n_active_v .and. i .gt. n_active_o) then
!!
!               a_inactive_i_inactive = a_inactive_i_inactive + (solution_ai(ai, 1))**2
!!
!            endif