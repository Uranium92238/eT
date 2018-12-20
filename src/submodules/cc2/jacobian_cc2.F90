submodule (cc2_class) jacobian
!

!
!!
!!    Jacobian submodule (CC2)
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | ν >.
!!
!
   implicit none
!
!
contains
!
!
   subroutine effective_jacobian_transformation_cc2(wf, c)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c
!
   end subroutine effective_jacobian_transformation_cc2
!
   module subroutine construct_jacobian_cc2_A1_cc2(wf, rho_a_i, c_b_j)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai^A1 = sum_bj (2 g_aijb - g_abji) * c_bj
!!
!!    Separate calculation of both terms due to batching
!!
      implicit none
!
      class(cc2) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_b_j
      real(dp), dimension(:,:), allocatable :: c_j_b
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ai_jb
      real(dp), dimension(:,:), allocatable :: g_ab_ji
!
!     Indices
!
      integer(i15) :: b, j
!
!
!     Explicit reordering of c_b_j
!
      call mem%alloc(c_j_b, wf%n_o, wf%n_v)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            c_j_b(j, b) = c_b_j(b, j)
!
         enddo
      enddo
!
      call mem%alloc(g_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
      call wf%get_voov(g_ai_jb)
!
!     rho_a_i = rho_a_i + sum_bj 2 g_aijb * c_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  g_ai_jb,           &
                  (wf%n_v)*(wf%n_o), &
                  c_j_b,             &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_a_i,           &
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(g_ab_ji, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call wf%get_vvoo(g_ab_ji)
!
      call sort_1234_to_1432(g_ab_ji, g_ai_jb, (wf%n_v), (wf%n_v), (wf%n_o), (wf%n_o))
!
      call mem%dealloc(g_ab_ji, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
!     rho_a_i = rho_a_i - sum_bj g_abji * c_jb
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  g_ai_jb,           &
                  (wf%n_v)*(wf%n_o), &
                  c_j_b,             &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_a_i,           &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(g_ai_jb, (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
      call mem%dealloc(c_j_b, wf%n_o, wf%n_v)
!
   end subroutine construct_jacobian_cc2_A1_cc2
!
!
   module subroutine construct_jacobian_cc2_B1_cc2(wf, rho_a_i, c_b_j)
!!
!!    Jacobian CC2 B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai^B1 = 2 L_ckjb * c_bj * sum_bj (2 t^ac_ik - t^ac_ki)
!!
!!    Separate calculation of both terms due to batching
!!
      implicit none
!
      class(cc2) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_b_j
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_a_i
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: I_kc
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kc_jb
      real(dp), dimension(:,:), allocatable :: g_ai_ck
      real(dp), dimension(:,:), allocatable :: L_kc_bj
      real(dp), dimension(:,:), allocatable :: u_ai_kc
!
!     Indices
!
      integer(i15) :: b, c, j, k, jb, kc, jc, kb, bj, a, i, ai, ck, ak, ci
!
!     Construct L_kc_jb
!
      call mem%alloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_kc_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_jb)
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  kc = wf%n_o*(c-1) + k
                  jb = wf%n_o*(b-1) + j
                  jc = wf%n_o*(c-1) + j
                  kb = wf%n_o*(b-1) + k
                  bj = wf%n_v*(j-1) + b
!
                  L_kc_bj(kc, bj) = two * g_kc_jb(kc, jb) - g_kc_jb(kb, jc)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     I_kc = sum_jb L_kcjb * c_bj
!
      call mem%alloc(I_kc, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  L_kc_bj,           &
                  (wf%n_o)*(wf%n_v), &
                  c_b_j,             &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  I_kc,              &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_kc_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     t_ai_ck = - g_aick/ε^{ac}_{ik}
!
      call mem%alloc(g_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%get_vovo(g_ai_ck)
!
      call mem%alloc(u_ai_kc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do c = 1, wf%n_v
               do k = 1, wf%n_o
!
                  kc = wf%n_o*(c-1) + k
                  ai = wf%n_v*(i-1) + a
                  ak = wf%n_v*(k-1) + a
                  ck = wf%n_v*(k-1) + c
                  ci = wf%n_v*(i-1) + c
!
                  u_ai_kc(ai,kc) = - (two*g_ai_ck(ai,ck)- g_ai_ck(ak,ci))&
                                          /(wf%fock_diagonal(a + wf%n_o, 1) &
                                          + wf%fock_diagonal(c + wf%n_o, 1)&
                                          - wf%fock_diagonal(i, 1) - wf%fock_diagonal(k, 1))
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_a_i = rho_a_i + sum_ck 2 t_aick I_kc
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  u_ai_kc,           &
                  (wf%n_v)*(wf%n_o), &
                  I_kc,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_a_i,           &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(I_kc, wf%n_o, wf%n_v)
      call mem%alloc(u_ai_kc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!
   end subroutine
!
!
end submodule jacobian
