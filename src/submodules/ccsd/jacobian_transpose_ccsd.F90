submodule (ccsd_class) jacobian_transpose_ccsd
!
!!
!!    Jacobian transpose submodule (CCSD)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
!!    Contains the following family of procedures of the CCSD class:
!!
!!    jacobian_transpose_ccsd_transformation: performs the transformation by the CCSD
!!                                            Jacobian transpose matrix A^T, placing the result in the
!!                                            incoming vector. 
!!
!!    jacobian_transpose_ccsd_x1:             adds the X1 term to the transformed singles vector; 
!!                                            x = a, b, c, ..., g 
!!    jacobian_transpose_ccsd_x2:             adds the X2 term to the transformed doubles vector; 
!!                                            x = a, b, ..., i
!! 
!
   implicit none 
!
!
contains
!
!
   module subroutine jacobian_transpose_ccsd_a1_ccsd(wf, sigma_a_i, b_a_i)
!!
!!    Jacobian transpose CCSD A1 
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Andreas Skeidsvoll, 2018
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
      real(dp), dimension(:,:), allocatable :: t_ld_ck
!
      real(dp), dimension(:,:), allocatable :: X_ld ! An intermediate, see below
!
      real(dp), dimension(:,:), allocatable :: g_ia_ld ! g_iald 
      real(dp), dimension(:,:), allocatable :: L_ai_ld ! L_iald 
!
      integer(i15) :: k = 0, c = 0, d = 0, l = 0, ck = 0, dk = 0, dl = 0
      integer(i15) :: ld = 0, cl = 0, ckdl = 0, cldk = 0, i = 0, a = 0
      integer(i15) :: id = 0, la = 0, ia = 0, ai = 0
!
!     Read the amplitudes from disk 
! 
      !call wf%read_double_amplitudes 
!
!     Form u_ld_ck = u_kl^cd = 2 * t_kl^cd - t_lk^cd = 2 * t_ld_ck(ld, ck) - t_ld_ck(kd, cl)
!
      call mem%alloc(t_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call squareup_and_sort_1234_to_4312(wf%t2, t_ld_ck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      u_ld_ck = zero
!
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, t_ld_ck, 1, u_ld_ck, 1)
      call add_4231_to_1234(-one, t_ld_ck, u_ld_ck, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      !call wf%destruct_double_amplitudes
!
!     Form the intermediate X_ld = sum_ck u_ld_ck b_ck  
!
      call mem%alloc(X_ld, (wf%n_v)*(wf%n_o), 1)
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
      call mem%dealloc(u_ld_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     Form L_ai_ld = L_iald = 2 * g_iald - g_idla 
!                           = 2 * g_ia_ld(ia,ld) - g_ia_ld(id,la)
!
      call mem%alloc(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_ia_ld)
!
      call mem%alloc(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_ai_ld = zero
!
      call daxpy((wf%n_v)**2*(wf%n_o)**2, two, g_ia_ld, 1, L_ai_ld, 1)
      call add_1432_to_1234(-one, g_ia_ld, L_ai_ld, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ia_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
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
      call mem%dealloc(L_ai_ld, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(X_ld, (wf%n_o)*(wf%n_v), 1)
!
   end subroutine jacobian_transpose_ccsd_a1_ccsd
!
!
end submodule jacobian_transpose_ccsd