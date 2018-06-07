submodule (ccs_class) jacobian
!
!!
!!    Jacobian submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    jacobian_ccs_transformation: directs the transformation by A.
!!    jacobian_ccs_a1:             adds the A1 term to the transformed vector. 
!!    jacobian_ccs_b1:             adds the B1 term to the transformed vector.
!!
   implicit none 
!
   character(len=40) :: integral_type
!
contains
!
!
   module subroutine jacobian_ccs_transformation_ccs(wf, c_a_i)
!!
!!    Jacobian CCS transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >. 
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!! 
!!    On exit, c is overwritten by rho. 
!!
      implicit none
!
      class(ccs) :: wf 
      real(dp), dimension(wf%n_v, wf%n_o)   :: c_a_i       
!
      real(dp), dimension(:,:), allocatable :: rho_a_i
!
      call wf%mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
!     A1-term
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
!
!     B1-term
!
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
!
!     Place rho_a_i in c_a_i
!
      c_a_i = zero
!
      call dcopy((wf%n_o)*(wf%n_v), rho_a_i, 1, c_a_i, 1)
!
      call wf%mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
   end subroutine jacobian_ccs_transformation_ccs
!
!
   module subroutine jacobian_ccs_a1_ccs(wf, rho, c1)
!!
!!    Jacobian CCS A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Calculates the A1 term,
!!
!!       sum_b F_ab*c_bi - sum_j F_ji*c_aj
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_o) :: c1
      real(dp), dimension(wf%n_v,wf%n_o) :: rho
!
!     sum_b F_a_b * c_b_i
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
!     - sum_j c_a_j * F_j_i
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
!
   end subroutine jacobian_ccs_a1_ccs
!
!
   module subroutine jacobian_ccs_b1_ccs(wf,rho,c1)
!!
!!    Jacobian CCS B1 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the B1 term,
!!
!!       sum_bj L_aijb*c_bj = sum_bj (2*g_aijb-g_abji)c_bj,
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
!   
      real(dp), dimension(wf%n_v,wf%n_o) :: c1
      real(dp), dimension(wf%n_v,wf%n_o) :: rho       

!
!     Integrals
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
!     Reorderings of c and rho
!
      real(dp), dimension(:,:), allocatable :: c_jb
!
!     Batching variables
!
      integer(i15) :: b_batch = 0, b_first = 0, b_last = 0, b_length = 0
      integer(i15) :: required = 0, available = 0, n_batch = 0, batch_dimension = 0
      integer(i15) :: max_batch_length = 0
!
!     Indices
!
      integer(i15) :: a = 0, b = 0
      integer(i15) :: i = 0, j = 0
!
      integer(i15) :: ab = 0
      integer(i15) :: ai = 0, jb = 0, jb_full = 0
      integer(i15) :: ji = 0
!
      logical :: reorder
!
!     Allocate and construct g_ai_jb
!
      call wf%mem%alloc(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_ov(integral_type, g_ai_jb)
!
!     Preparing for batching over b
!
      required = max(((wf%n_o)**2)*((wf%n_v)**2)                                 &  !
                        + 2*(wf%n_J)*(wf%n_v)**2 + 2*(wf%n_J)*(wf%n_v)*(wf%n_o), &  ! Needed to get L_ab_J
                        2*((wf%n_o)**2)*((wf%n_v)**2)                            &  !
                        + (wf%n_J)*(wf%n_v)**2 + (wf%n_J)*(wf%n_o)**2,           &  ! Needed to get g_ab_ij
                        3*((wf%n_o)**2)*((wf%n_v)**2))                              ! Needed to get L_ai_jb
!
      required = 4*required ! In words
!
      batch_dimension  = wf%n_v ! Batch over the virtual index b
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, wf%mem%available, max_batch_length, n_batch, batch_dimension)           
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
!        Allocate and construct g_ab_ji
!
         call wf%mem%alloc(g_ab_ji, (wf%n_v)*b_length, (wf%n_o)**2)
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_oo(integral_type, &
                           g_ab_ji,       &
                           1,             &
                           wf%n_v,        &
                           b_first,       &
                           b_last,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)        
!
!        Allocate L_ai_jb
!
         call wf%mem%alloc(L_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
         L_ai_jb = zero
!
!        Construct L_ai_jb = 2*g_ai_jb - g_ab_ij
!
         do i = 1, wf%n_o
            do b = 1, b_length
               do j = 1, wf%n_o
!
                  ji = index_two(j, i, wf%n_o)
                  jb = index_two(j, b, wf%n_o)
                  jb_full = index_two(j, b + b_first - 1, wf%n_o)
!
                  do a = 1, wf%n_v
!
                     ai = index_two(a, i, wf%n_v)
                     ab = index_two(a, b, wf%n_v)
!
                     L_ai_jb(ai, jb) = two*g_ai_jb(ai, jb_full) - g_ab_ji(ab, ji)
!
                  enddo
               enddo
            enddo
         enddo

!
!        Deallocate g_ai_jb and g_ab_ji
!
         call wf%mem%dealloc(g_ab_ji, (wf%n_v)*b_length, (wf%n_o)**2)
!
!        Allocate reordering of c and rho
!
         call wf%mem%alloc(c_jb, (wf%n_o)*b_length, 1)
!
!        reordered c amplitudes
!
         c_jb = zero
         do j = 1, wf%n_o
            do b = 1, b_length
!
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
         call wf%mem%dealloc(L_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*b_length)
!
!        Deallocate c reordered
!
         call wf%mem%dealloc(c_jb, (wf%n_o)*b_length, 1)
!
      enddo ! Looping over batches
!
      call wf%mem%dealloc(g_ai_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_ccs_b1_ccs
!
!
end submodule jacobian
