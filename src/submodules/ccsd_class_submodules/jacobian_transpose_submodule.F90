submodule (ccsd_class) jacobian_transpose
!
!!
!!    Jacobian transpose submodule (CCSD)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad, May 2017
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    jacobian_transpose_ccsd_transformation: performs the transposed Jacobian transformation (A^T)
!! 
!
   implicit none 
!
!
contains
!
!
   module subroutine jacobian_transpose_ccsd_transformation_ccsd(wf, b_a_i, b_aibj)
!!
!!    Jacobian transpose transformation (CCSD)
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
      class(ccsd) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
      real(dp), dimension(wf%n_t2am, 1)   :: b_aibj 
!
      real(dp), dimension(:,:), allocatable :: sigma_a_i
      real(dp), dimension(:,:), allocatable :: sigma_ai_bj
!
!     Allocate the transformed singles vector 
!
      call allocator(sigma_a_i, wf%n_v, wf%n_o)
      sigma_a_i = zero 
!
!     Calculate and add the CCS contributions to the 
!     singles transformed vector 
!
      call wf%jacobian_transpose_ccs_a1(sigma_a_i, b_a_i)
      call wf%jacobian_transpose_ccs_b1(sigma_a_i, b_a_i)
!
!     Calculate and add the CCSD contributions to the
!     singles transformed vector 
!
      call wf%jacobian_transpose_ccsd_a1(sigma_a_i, b_a_i)
!
   end subroutine jacobian_transpose_ccsd_transformation_ccsd
!
!
   module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD A1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_ckdl b_ck L_iald u_kl^cd,
!! 
!!    abd adds it to the transformed vector sigma_a_i.
!!
      implicit none 
!
      class(ccsd) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o) :: b_a_i 
      real(dp), dimension(wf%n_v, wf%n_o) :: sigma_a_i 
!
      real(dp), dimension(:,:), allocatable :: u_ld_ck
!
      real(dp), dimension(:,:), allocatable :: X_ld ! An intermediate, see below
!
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J 
      real(dp), dimension(:,:), allocatable :: g_ia_ld ! g_iald 
      real(dp), dimension(:,:), allocatable :: L_ai_ld ! L_iald 
!
      integer(i15) :: k = 0, c = 0, d = 0, l = 0, ck = 0, dk = 0, dl = 0
      integer(i15) :: ld = 0, cl = 0, ckdl = 0, cldk = 0, i = 0, a = 0, ld = 0
      integer(i15) :: id = 0, la = 0, ia = 0, ai = 0
!
!     Read the amplitudes from disk 
!
      call wf%initialize_amplitudes 
      call wf%read_double_amplitudes 
!
!     Form u_ld_ck = u_kl^cd = 2 * t_kl^cd - t_lk^cd 
!
      call allocator(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            ck = index_two(c, k, wf%n_v)
!
            do d = 1, wf%n_v
!
               dk = index_two(d, k, wf%n_v)
!
               do l = 1, wf%n_o
!
                  dl = index_two(d, l, wf%n_v)
                  ld = index_two(l, d, wf%n_o)
                  cl = index_two(c, l, wf%n_v)
!
                  ckdl = index_packed(ck, dl)
                  cldk = index_packed(cl, dk)
!
                  u_ld_ck(ld, ck) = two*(wf%t2am(ckdl,1)) - wf%t2am(cldk,1)
!
               enddo
            enddo
         enddo
      enddo
!
      call wf%destruct_amplitudes
!
!     Form the intermediate X_ld = sum_ck u_ld_ck b_ck  
!
      call allocator(X_ld, (wf%n_v)*(wf%n_o), 1)
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  u_ld_ck,           &
                  (wf%n_v)*(wf%n_o), &
                  b_a_i,             & ! "b_ck"
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_ld,              &
                  (wf%n_v)*(wf%n_o))
!
      call deallocator(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form L_ai_ld = L_iald = 2 * g_iald - g_idla 
!                           = 2 * g_ia_ld(ia,ld) - g_ia_ld(id,la)
!
      call allocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call wf%get_cholesky_ia(L_ia_J)
!
      call allocator(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call dgemm('N','T',            &
                  (wf%n_o)*(wf%n_v), & 
                  (wf%n_o)*(wf%n_v), &
                  wf%n_J,            &
                  one,               &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  L_ia_J,            &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  g_ia_ld,           &
                  (wf%n_o)*(wf%n_v))
!
      call deallocator(L_ia_J, (wf%n_o)*(wf%n_v), wf%n_J)
!
      call allocator(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ai_ld = zero 
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o
!
            ld = index_two(l, d, wf%n_o)
!
            do i = 1, wf%n_o
!
               id = index_two(i, d, wf%n_o)
!
               do a = 1, wf%n_v
!
                  la = index_two(l, a, wf%n_o)
                  ia = index_two(i, a, wf%n_o)
                  ai = index_two(a, i, wf%n_v)
!
                  L_ai_ld(ai, ld) = two*g_ia_ld(ia, ld) - g_ia_ld(id, la)
!
               enddo
            enddo
         enddo
      enddo
!
      call deallocator(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Calculate and add sum_ckdl b_ck L_iald u_kl^cd 
!                       = sum_ld L_ai_ld X_ld
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_ai_ld,           &
                  (wf%n_o)*(wf%n_v), &
                  X_ld,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  sigma_a_i,         & ! "sigma_ai"
                  (wf%n_o)*(wf%n_v))
!
   end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
end submodule jacobian_transpose