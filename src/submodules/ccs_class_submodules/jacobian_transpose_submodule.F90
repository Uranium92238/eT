submodule (ccs_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (CCS)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, June 2017
!!
!!    Contains the following family of procedures of the CCS class:
!!
!!    jacobian_transpose_ccs_transformation: performs the transposed Jacobian transformation (A^T)
!!    jacobian_transpose_ccs_a1:             adds the CCS A1 term for the A^T transformation 
!!    jacobian_transpose_ccs_b1:             adds the CCS B1 term for the A^T transformation 
!! 
!
   implicit none 
!
   character(len=40) :: integral_type
!
contains
!
!
   module subroutine jacobian_transpose_ccs_transformation_ccs(wf, b_a_i)
!!
!!    Jacobian transpose transformation (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation 
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    The transformation is performed as sigma^T = b^T A, where b is the vector
!!    sent to the routine. On exit, the vector b is equal to sigma (the transformed
!!    vector).
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
!
      real(dp), dimension(:,:), allocatable :: sigma_a_i
!
!     Allocate the transformed vector 
!
      call allocator(sigma_a_i, wf%n_v, wf%n_o)
      sigma_a_i = zero 
!
!     Calculate and add the CCS contributions 
!
      call wf%jacobian_transpose_ccs_a1(sigma_a_i, b_a_i)
      call wf%jacobian_transpose_ccs_b1(sigma_a_i, b_a_i)
!
!     Copy the result into the incoming vector 
!
      call dcopy((wf%n_o)*(wf%n_v), sigma_a_i, 1, b_a_i, 1)
!
!     Deallocate the transformed vector 
!
      call deallocator(sigma_a_i, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccs_transformation_ccs
!
!
   module subroutine jacobian_transpose_ccs_a1_ccs(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose A1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_c b_ci F_ca - sum_k b_ak F_ik,
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
!
!     Add sum_c F_ca b_ci = sum_c F_ac^T b_ci     
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ab, &
                  wf%n_v,     &
                  b_a_i,      &
                  wf%n_v,     &
                  one,        &
                  sigma_a_i,  &
                  wf%n_v)
!
!     Add - sum_k b_ak F_ik = - sum_k b_ak F_ki^T 
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  b_a_i,      &
                  wf%n_v,     &
                  wf%fock_ij, &
                  wf%n_o,     &
                  one,        &
                  sigma_a_i,  &
                  wf%n_v)
!
   end subroutine jacobian_transpose_ccs_a1_ccs
!
!
   module subroutine jacobian_transpose_ccs_b1_ccs(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the B1 term,
!!
!!       sum_ck L_ckia b_ck
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
!
      real(dp), dimension(:,:), allocatable :: L_ck_J ! L_ck^J 
      real(dp), dimension(:,:), allocatable :: L_ia_J ! L_ia^J 
!
      real(dp), dimension(:,:), allocatable :: g_ck_ia ! g_ckia 
!
      real(dp), dimension(:,:), allocatable :: L_ik_J ! L_ik^J 
      real(dp), dimension(:,:), allocatable :: L_ca_J ! L_ca^J 
!
      real(dp), dimension(:,:), allocatable :: g_ca_ik ! g_caik 
!
      real(dp), dimension(:,:), allocatable :: L_ai_ck ! L_ckia = 2 * g_ckia - g_caik
!
!     Batching variables 
!
      integer(i15) :: required = 0, available = 0 
      integer(i15) :: batch_dimension = 0, max_batch_length = 0, n_batch = 0
!
      integer(i15) :: a_batch = 0, a_length = 0, a_first = 0, a_last = 0
!
!     Indices 
!
      integer(i15) :: a = 0, i = 0, c = 0, k = 0, ck = 0, ik = 0, ai = 0
      integer(i15) :: ia = 0, ca = 0
!
!     Form the integral g_ck_ia = g_ckia  
!
      call allocator(g_ck_ia, (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
!
      integral_type = 'electronic_repulsion'
      call wf%get_vo_ov(integral_type, g_ck_ia)
!
!     :: Form L_ai_ck = L_ckia in batches over a ::
!
      call allocator(L_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      L_ai_ck = zero
!
!     Prepare for batching over index a
! 
      required = 1 ! Not a correct estimate - needs to be set!
!     
      required  = 4*required ! In words
      available = get_available()
!
      batch_dimension  = wf%n_v ! Batch over the virtual index a
      max_batch_length = 0      ! Initilization of unset variables 
      n_batch          = 0
!
      call num_batch(required, available, max_batch_length, n_batch, batch_dimension)           
!
!     Loop over the number of a batches 
!
      do a_batch = 1, n_batch
!
!        For each batch, get the limits for the a index 
!
         call batch_limits(a_first, a_last, a_batch, max_batch_length, batch_dimension)
         a_length = a_last - a_first + 1 
!
!        Form g_ca_ik
!
         call allocator(g_ca_ik, (wf%n_v)*a_length, (wf%n_o)**2)
!
         integral_type = 'electronic_repulsion'
         call wf%get_vv_oo(integral_type, &
                           g_ca_ik,       &
                           1,             &
                           wf%n_v,        &
                           a_first,       &
                           a_last,        &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)
!
!        Set L_ai_ck = L_ckia = 2 * g_ckia - g_caik 
!        for the current batch over a 
!
         do k = 1, wf%n_o
            do c = 1, wf%n_v
!
               ck = index_two(c, k, wf%n_v)
!
               do i = 1, wf%n_o
!
                  ik = index_two(i, k, wf%n_o)
!
                  do a = 1, a_length
!
                     Ai = index_two(a + a_first - 1, i, wf%n_v) ! Full space a 
                     iA = index_two(i, a + a_first - 1, wf%n_o) ! Full space a 
                     ca = index_two(c, a, wf%n_v)
!
                     L_ai_ck(Ai, ck) = two*g_ck_ia(ck, iA) - g_ca_ik(ca, ik)
!
                  enddo
               enddo
            enddo
         enddo
!
         call deallocator(g_ca_ik, (wf%n_v)*a_length, (wf%n_o)**2)
!
      enddo ! End of batches over a 
!
      call deallocator(L_ik_J, (wf%n_o)**2, wf%n_J)
      call deallocator(g_ck_ia, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Add sum_ck L_ckia b_ck = sum_ck L_ai_ck b_ck 
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  L_ai_ck,           &
                  (wf%n_v)*(wf%n_o), &
                  b_a_i,             & ! "b_ai"
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  sigma_a_i,         & ! "sigma_ai"
                  (wf%n_v)*(wf%n_o))
!
      call deallocator(L_ai_ck, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccs_b1_ccs
!
!
end submodule jacobian_transpose