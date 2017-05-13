submodule (ccs_class) jacobian
!
!!
!!    Jacobian submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!
!
   implicit none 
!
!
contains
!
!
      subroutine construct_right_transform_ccs(wf)
!!
!!       Construct Right Transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!
         implicit none
!
         class(ccs) :: wf      
!
!        Allocate rho and c vectors
!
         call allocator(wf%c1am, wf%n_v, wf%n_o)
         call allocator(wf%rho1_a_i, wf%n_v, wf%n_o)
!
         wf%c1am = zero
         wf%rho1_a_i = zero
!
!        Construct rho
!
!        A1-term
!
         call wf%rho_ccs_a1()
!
!        B1-term
!
         call wf%rho_ccs_b1()
!  
      end subroutine construct_right_transform_ccs
!
!
      subroutine rho_ccs_a1_ccs(wf)
!!
!!       A1 contribution to right transform of Jacobian
!!       Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!       Calculates the A1 term of the right transform of the
!!       Jacobian,
!!
!!       A1: sum_b F_ab*c_bi + sum_j F_ji*c_aj
!!
!!       and adds it to the rho vector.
!!
         implicit none
!
         class(ccs) :: wf        
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
                     wf%c1am,     &
                     wf%n_v,      &
                     one,         &
                     wf%rho1_a_i, &
                     wf%n_v)
!
!        - sum_j c_a_j * F_j_i
!
         call dgemm('N','N',      &
                     wf%n_v,      &
                     wf%n_o,      &
                     wf%n_o,      &
                     -one,        &
                     wf%c1am,     &
                     wf%n_v,      &
                     wf%fock_ij,  &
                     wf%n_o,      &
                     one,         &
                     wf%rho1_a_i, &
                     wf%n_v)
!
      end subroutine rho_ccs_a1_ccs
!
!
      subroutine rho_ccs_b1_ccs(wf)
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
         real(dp), dimension(:,:), allocatable :: rho_ai
!
!        Batching variables
!
         integer(i15) :: b_batch, b_first, b_last, b_length
         integer(i15) :: required, available, n_batch, batch_dimension, max_batch_length
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
         call wf%get_cholesky_ab(L_ab_J, b_first, b_last, (wf%n_v)*b_length, reorder)
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
         call allocator(rho_ai, (wf%n_o)*(wf%n_v), 1)
!
!        reordered c amplitudes
!
         do j = 1, wf%n_o
            do b = 1, b_length
               jb = index_two(j, b, wf%n_o)
!
               c_jb(jb,1) = wf%c1am(b, j)
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
                     zero,              &
                     rho_ai,            &
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
!        Add contribution to rho      
!
         do i = 1, wf%n_o
            do a = 1, wf%n_v
               ai = index_two(a, i, wf%n_v)
!
               wf%rho1_a_i = wf%rho1_a_i + rho_ai(ai, 1)
            enddo
         enddo
!
!        Deallocate rho reordered
!
         call deallocator(rho_ai, (wf%n_o)*(wf%n_v), 1)
!
      enddo ! Looping over batches
!
!
      end subroutine rho_ccs_b1_ccs
!
!
end submodule