submodule (ccs_class) jacobian
!
!!
!!    Jacobian submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    initialize_trial_vectors_ccs:
!!    find_lowest_orbital_diff:
!!    jacobian_transformation:
!!    rho_ccs_a1:
!!    rho_ccs_b1:
!
   implicit none 
!
!
contains
!
!
      module subroutine initialize_trial_vectors_ccs(wf)
!!
!!       Initialize trial vectors
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Initializes start trial vectors for the calculation of 
!!       singlet excited states and writes them to file 'trial_vecs'.
!!
!!       n start vectors are constructed by finding the n lowest orbital differences,      
!!       where n = n_singlet_states. Vector i has a 1.0D0 at the element corresponding to the i'th lowest
!!       orbital difference and 0.0d0 everywhere else
!!
         implicit none
!
         class(ccs) :: wf
!       
         integer(i15), dimension(:,:), allocatable :: index_lowest_obital_diff
!
         real(dp), dimension(:,:), allocatable :: c
! 
         integer(i15) :: i = 0, j = 0
!
         integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
!
!        Allocate array for the indices of the lowest orbital differences
!
         call allocator_int( index_lowest_obital_diff, wf%tasks%n_singlet_states, 1)
         index_lowest_obital_diff = zero
!
!        Find indecies of lowest orbital differences
!
         call wf%find_start_trial_indices(index_lowest_obital_diff)
!
!        Generate start trial vectors c and write to file
!
         call allocator(c, wf%n_parameters, 1)
!
!        Prepare for writing trial vectors to file
!
         call generate_unit_identifier(unit_trial_vecs)
         open(unit=unit_trial_vecs, file='trial_vec', action='write', status='unknown', &
           access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
         do i = 1, (wf%tasks%n_singlet_states)
            c = zero
            c(index_lowest_obital_diff(i,1),1) = one
            write(unit_trial_vecs, rec=i, iostat=ioerror) (c(j,1), j = 1, wf%n_parameters)
         enddo
!
!        Close file
!     
         close(unit_trial_vecs)
!
!        Deallocate c
!
         call deallocator(c, wf%n_parameters, 1)
!
!        Deallocate index_lowest_obital_diff
!
         call deallocator_int( index_lowest_obital_diff, wf%tasks%n_singlet_states, 1)
!
      end subroutine initialize_trial_vectors_ccs
!
!
      module subroutine trial_vectors_from_stored_solutions_ccs(wf)
!!
!!
!!
         implicit none
!
         class(ccs) :: wf
!
         logical      :: solution_exists = .false.
         logical      :: more_trials = .true.
!
         integer(i15) :: ioerror = 0, unit_solution = 0, unit_trial_vecs = 0
         integer(i15) :: number_of_solutions = 0
!
         integer(i15) :: i = 0, j = 0
!
         real(dp), dimension(:,:), allocatable :: c_i, c_j
!
         real(dp) :: ddot, dot_prod = zero, norm = zero 
!
!        Open solution vector file - if it does not exist return
!
         inquire(file='solution_vectors', exist=solution_exists)
!
!        If no solution vector file, return and use orbital differences.
!
         if (.not. solution_exists) return
!
!        Open files
!
         call generate_unit_identifier(unit_trial_vecs)
         open(unit=unit_trial_vecs, file='trial_vec', action='readwrite', status='old', &
         access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror)
!
         call generate_unit_identifier(unit_solution)
         open(unit=unit_solution, file='solution_vectors', action='read', status='unknown', &
         access='direct', form='unformatted', recl=dp*(wf%n_parameters), iostat=ioerror) 
!
!        Allocate c_i
!
         call allocator(c_i, wf%n_parameters, 1)
         c_i = zero
!
         i = 1
         do while ((i .le. wf%tasks%n_singlet_states) .and. more_trials)
!
!        Read old solutions and count them
!
            read(unit_solution, rec=i, iostat=ioerror) c_i
            if (ioerror .ne. 0) write(unit_output,*) 'Error reading solution vecs'
!
            if (ioerror .eq. 0) then
!
               write(unit_trial_vecs, rec = i)c_i
!
            else
!
               more_trials = .false.
!
            endif
!
            i = i + 1
!
         enddo
!
!        Deallocate c_i 
!
         call deallocator(c_i, wf%n_parameters, 1)
!
!        Close solution file
!
         close(unit_solution)
!
!        Allocate c_i and c_j
!
         call allocator(c_i, wf%n_parameters, 1)
         call allocator(c_j, wf%n_parameters, 1)
         c_i = zero
         c_j = zero
!
!        Reorthonormalize trial vectors
!
         do i = 1, wf%tasks%n_singlet_states
!
            read(unit_trial_vecs, rec=i, iostat=ioerror) c_i
!
            do j = 1, i-1
!
               read(unit_trial_vecs, rec=j, iostat=ioerror) c_j
               dot_prod = ddot(wf%n_parameters, c_j, 1, c_i, 1)
               call daxpy(wf%n_parameters, -dot_prod, c_j, 1, c_i, 1)
!  
               norm = sqrt(ddot(wf%n_parameters, c_i, 1, c_i, 1))
               call dscal(wf%n_parameters, one/norm, c_i, 1)
!
            enddo
            write(unit_trial_vecs, rec = i)c_i
         enddo
!  
         call deallocator(c_i, wf%n_parameters, 1)
         call deallocator(c_j, wf%n_parameters, 1)  
!
!        Close trial vector file
!
         close(unit_trial_vecs)     
!
      end subroutine trial_vectors_from_stored_solutions_ccs
!
!
      module subroutine find_start_trial_indices_ccs(wf, index_list)
!!
!!       Find indices for lowest orbital differences
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!
         implicit none
!
         class(ccs) :: wf
         integer(i15), dimension(wf%tasks%n_singlet_states,1), intent(inout) :: index_list
!
         real(dp), dimension(:,:), allocatable :: orbital_diff
         real(dp), dimension(:,:), allocatable :: lowest_orbital_diff
!
         integer(i15) :: a = 0, i = 0, j = 0
!
         integer(i15) :: ai = 0
!
         real(dp)     :: max
         integer(i15) :: max_pos 
!
         real(dp)     :: swap     = zero
         integer(i15) :: swap_int = 0
!
!        Allocate orbital_diff
!
         call allocator(orbital_diff,wf%n_parameters,1)
         orbital_diff = zero
!
!        Calculate orbital differences
!
         call wf%calculate_orbital_differences(orbital_diff)
!
!        Finding lowest orbital differences
!
         call allocator(lowest_orbital_diff, wf%tasks%n_singlet_states, 1)
         
         lowest_orbital_diff = zero
!
         call get_n_lowest(wf%tasks%n_singlet_states, wf%n_parameters, orbital_diff, lowest_orbital_diff, index_list)
!
         call deallocator(orbital_diff,wf%n_parameters,1)
!
         call deallocator(lowest_orbital_diff, wf%tasks%n_singlet_states, 1)
!
!
      end subroutine find_start_trial_indices_ccs
!
!
      module subroutine calculate_orbital_differences_ccs(wf,orbital_diff)
!!
!!       Calculate and return orbital differences
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad May 2017
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_parameters, 1) :: orbital_diff
!
         integer(i15) :: a = 0, i = 0
         integer(i15) :: ai = 0
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
               ai = index_two(a, i, wf%n_v)
               orbital_diff(ai, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1)
            enddo
         enddo
!
      end subroutine calculate_orbital_differences_ccs
!
!
      module subroutine transform_trial_vectors_ccs(wf, first_trial, last_trial)
!!
!!       Construct Jacobian Transformation of trial vectors
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Each trial vector in first_trial to last_trial is read from file and
!!       transformed before the transformed vector is written to file.
!!
         implicit none
!
         class(ccs) :: wf
!
         integer(i15), intent(in) :: first_trial, last_trial ! Which trial_vectors we are to transform
!
         real(dp), dimension(:,:), allocatable :: c_a_i
!
         integer(i15) :: unit_trial_vecs = 0, unit_rho = 0, ioerror = 0
         integer(i15) :: trial = 0 
!
!        Allocate c_a_i
!
         call allocator(c_a_i, wf%n_v, wf%n_o)
         c_a_i = zero 
!
!        Open trial vector and transformed vector files
!
         call generate_unit_identifier(unit_trial_vecs)
         open(unit=unit_trial_vecs, file='trial_vec', action='read', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
         call generate_unit_identifier(unit_rho)
         open(unit=unit_rho, file='transformed_vec', action='write', status='old', &
           access='direct', form='unformatted', recl=dp*(wf%n_v)*(wf%n_o), iostat=ioerror)
!
!        For each trial vector: Read, transform and write
!               
         do trial = first_trial, last_trial
!
            read(unit_trial_vecs, rec=trial, iostat=ioerror) c_a_i
!
            call wf%jacobian_transformation(c_a_i)
!
!           Write transformed vector to file
!
            write(unit_rho, rec=trial, iostat=ioerror) c_a_i
          
         enddo
         close(unit_trial_vecs) 
         close(unit_rho)                                
!
!        Deallocate c_a_i
!
         call deallocator(c_a_i, wf%n_v, wf%n_o)
!
      end subroutine transform_trial_vectors_ccs
!
!
      module subroutine jacobian_transformation_ccs(wf, c_a_i)
!!
!!       Jacobian transformation
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!
         implicit none
!
         class(ccs) :: wf 
         real(dp), dimension(wf%n_v, wf%n_o)   :: c_a_i       
!
         real(dp), dimension(:,:), allocatable :: rho_a_i
!
         call allocator(rho_a_i, wf%n_v, wf%n_o)
         rho_a_i = zero
!
!        A1-term
!
         call wf%jacobian_ccs_a1(c_a_i,rho_a_i)
!
!        B1-term
!
         call wf%jacobian_ccs_b1(c_a_i,rho_a_i)
!
         call dcopy((wf%n_o)*(wf%n_v), rho_a_i, 1, c_a_i, 1)
!
         call deallocator(rho_a_i, wf%n_v, wf%n_o)
!
      end subroutine jacobian_transformation_ccs
!
      module subroutine jacobian_ccs_a1_ccs(wf,rho,c1)
!!
!!       A1 contribution to right transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Calculates the A1 term of the right transform of the
!!       Jacobian,
!!
!!       A1: sum_b F_ab*c_bi - sum_j F_ji*c_aj
!!
!!       and adds it to the rho vector.
!!
         implicit none
!
         class(ccs) :: wf
         real(dp), dimension(wf%n_v,wf%n_o) :: c1
         real(dp), dimension(wf%n_v,wf%n_o) :: rho
!   
!
!        sum_b F_a_b * c_b_i
!
         call dgemm('N', 'N',     &
                     wf%n_v,      &
                     wf%n_o,      &
                     wf%n_v,      &
                     one,         &
                     wf%fock_ab,  &
                     wf%n_v,      &
                     c1,          &
                     wf%n_v,      &
                     one,         &
                     rho,         &
                     wf%n_v)
!
!        - sum_j c_a_j * F_j_i
!
         call dgemm('N','N',      &
                     wf%n_v,      &
                     wf%n_o,      &
                     wf%n_o,      &
                     -one,        &
                     c1,          &
                     wf%n_v,      &
                     wf%fock_ij,  &
                     wf%n_o,      &
                     one,         &
                     rho,         &
                     wf%n_v)
!
      end subroutine jacobian_ccs_a1_ccs
!
!
      module subroutine jacobian_ccs_b1_ccs(wf,rho,c1)
!!
!!       B1 contribution to right transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Calculates the B1 term of the right transform of the
!!       Jacobian,
!!
!!       B1: sum_bj L_aijb*c_bj = sum_bj (2*g_aijb-g_abji)c_bj
!!
!!       and adds it to the rho vector.
!!
         implicit none
!
         class(ccs) :: wf   
         real(dp), dimension(wf%n_v,wf%n_o) :: c1
         real(dp), dimension(wf%n_v,wf%n_o) :: rho       

!
!        Integrals
!
         real(dp), dimension(:,:), allocatable :: L_ai_J
         real(dp), dimension(:,:), allocatable :: L_jb_J
         real(dp), dimension(:,:), allocatable :: L_ji_J
         real(dp), dimension(:,:), allocatable :: L_ab_J
!
         real(dp), dimension(:,:), allocatable :: g_ai_jb
         real(dp), dimension(:,:), allocatable :: g_ab_ji
         real(dp), dimension(:,:), allocatable :: L_ai_jb
!
!        Reorderings of c and rho
!
         real(dp), dimension(:,:), allocatable :: c_jb
!
!        Batching variables
!
         integer(i15) :: b_batch = 0, b_first = 0, b_last = 0, b_length = 0
         integer(i15) :: required = 0, available = 0, n_batch = 0, batch_dimension = 0
         integer(i15) :: max_batch_length = 0
!
!        Indices
!
         integer(i15) :: a = 0, b = 0
         integer(i15) :: i = 0, j = 0
!
         integer(i15) :: ab = 0
         integer(i15) :: ai = 0, jb = 0
         integer(i15) :: ji = 0
!
         logical :: reorder
!
!
!        Preparing for batching over b
!
         required = max(((wf%n_o)**2)*((wf%n_v)**2)                              &  !
                        + 2*(wf%n_J)*(wf%n_v)**2 + 2*(wf%n_J)*(wf%n_v)*(wf%n_o), &  ! Needed to get L_ab_J
                        2*((wf%n_o)**2)*((wf%n_v)**2)                            &  !
                        + (wf%n_J)*(wf%n_v)**2 + (wf%n_J)*(wf%n_o)**2,           &  ! Needed to get g_ab_ij
                        3*((wf%n_o)**2)*((wf%n_v)**2))                              ! Needed to get L_ai_jb
!
         required = 4*required ! In words
         available = get_available()
!
         batch_dimension  = wf%n_v ! Batch over the virtual index b
         max_batch_length = 0      ! Initilization of unset variables 
         n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do b_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(b_first, b_last, b_batch, max_batch_length, batch_dimension)
         b_length = b_last - b_first + 1 
!            
!        Allocate and get L_ai_J
! 
         call allocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
!
         call wf%get_cholesky_ai(L_ai_J)
!
!        Allocate and get L_jb_J
!
         call allocator(L_jb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
         call wf%get_cholesky_ia(L_jb_J)
!
!        Allocate and construct g_ai_jb constrained to the batch
!
         call allocator(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
!
         call dgemm('N', 'T',                                 &
                     (wf%n_v)*(wf%n_o),                       &
                     (wf%n_o)*b_length,                       &
                     wf%n_J,                                  &
                     one,                                     &
                     L_ai_J,                                  &
                     (wf%n_v)*(wf%n_o),                       &
                     L_jb_J(index_two(1, b_first, wf%n_o),1), &
                     (wf%n_v)*(wf%n_o),                       &
                     zero,                                    &
                     g_ai_jb,                                 &
                     (wf%n_v)*(wf%n_o))
!
!        Deallocate L_ai_J and L_jb_J
!
         call deallocator(L_ai_J, (wf%n_v)*(wf%n_o), wf%n_J)
         call deallocator(L_jb_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
!        Allocate and get L_ab_J
!
         call allocator(L_ab_J, (wf%n_v)*b_length, wf%n_J)
!
         reorder = .false.
         call wf%get_cholesky_ab(L_ab_J, b_first, b_last, reorder, 1, wf%n_v)
!
!        Allocate and get L_ji_J
!
         call allocator(L_ji_J, (wf%n_o)**2, wf%n_J)
!
         call wf%get_cholesky_ij(L_ji_J)
!
!        Allocate and construct g_ab_ji
!
         call allocator(g_ab_ji, (wf%n_v)*b_length, (wf%n_o)**2)
!
         call dgemm('N', 'T',           &
                     (wf%n_v)*b_length, &
                     (wf%n_o)**2,       &
                     wf%n_J,            &
                     one,               &
                     L_ab_J,            &
                     (wf%n_v)*b_length, &
                     L_ji_J,            &
                     (wf%n_o)**2,       &
                     zero,              &
                     g_ab_ji,           &
                     (wf%n_v)*b_length)
!
!        Deallocate and get L_ab_J and L_ji_J
!
         call deallocator(L_ab_J, (wf%n_v)*b_length, wf%n_J)
         call deallocator(L_ji_J, (wf%n_o)**2, wf%n_J)
!
!        Allocate L_ai_jb
!
         call allocator(L_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
         L_ai_jb = zero
!
!        Construct L_ai_jb = 2*g_ai_jb - g_ab_ij
!
         do i = 1, wf%n_o
            do b = 1, b_length
               do j = 1, wf%n_o
                  ji = index_two(j, i , wf%n_o)
                  jb = index_two(j, b, wf%n_o)
                  do a = 1, wf%n_v
                     ai = index_two(a, i, wf%n_v)
                     ab = index_two(a, b, wf%n_v)
!
                     L_ai_jb(ai, jb) = two*g_ai_jb(ai, jb) - g_ab_ji(ab, ji)
!
                  enddo
               enddo
            enddo
         enddo

!
!        Deallocate g_ai_jb and g_ab_ji
!
         call deallocator(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
         call deallocator(g_ab_ji, (wf%n_v)*b_length, (wf%n_o)**2)
!
!        Allocate reordering of c and rho
!
         call allocator(c_jb, (wf%n_o)*b_length, 1)
!
!        reordered c amplitudes
!
         c_jb = zero
         do j = 1, wf%n_o
            do b = 1, b_length
               jb = index_two(j, b, wf%n_o)
!
               c_jb(jb, 1) = c1(b+b_first-1, j)
!
            enddo
         enddo

!
!        Create rho contribution from the batch
!
         call dgemm('N','N',            &
                     (wf%n_v)*(wf%n_o), &
                     1,                 &
                     b_length*(wf%n_o), &
                     one,               &
                     L_ai_jb,           &
                     (wf%n_v)*(wf%n_o), &
                     c_jb,              &
                     b_length*(wf%n_o), &
                     one,               &
                     rho,               &
                     (wf%n_v)*(wf%n_o))         
!
!        Deallocate L_ai_jb
! 
         call deallocator(L_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
!
!        Deallocate c reordered
!
         call deallocator(c_jb, (wf%n_o)*b_length, 1)
!
      enddo ! Looping over batches
!
!
      end subroutine jacobian_ccs_b1_ccs
!
!
end submodule