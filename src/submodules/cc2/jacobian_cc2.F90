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
      real(dp), dimension(:,:), allocatable :: c_a_i
!
      real(dp), dimension(:,:), allocatable :: rho_a_i
!
      integer(i15) :: i, j, a, b, ai ! Index
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_a_i, wf%n_v, wf%n_o)
      rho_a_i = zero
!
      call mem%alloc(c_a_i, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c_a_i(a, i) = c(ai, 1)
!
         enddo
      enddo
!$omp end parallel do
!
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_a_i, c_a_i)
      call wf%jacobian_ccs_b1(rho_a_i, c_a_i)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call mem%dealloc(c_a_i, wf%n_v, wf%n_o)
      call mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
!
   end subroutine effective_jacobian_transformation_cc2
!
!
   module subroutine jacobian_cc2_A1_cc2(wf, rho_a_i, c_b_j)
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
      class(cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_b_j
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_a_i
!
      real(dp), dimension(:,:), allocatable :: c_j_b
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
   end subroutine jacobian_cc2_A1_cc2
!
!
   module subroutine jacobian_cc2_B1_cc2(wf, rho_a_i, c_b_j, eps_o, eps_v)
!!
!!    Jacobian CC2 B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai^B1 = 2 L_kcjb c_bj (2 t^ac_ik - t^ac_ki)
!!                - L_kcjb t^cb_ki c_aj - L_kcjb t^ca_kj c_bi
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_b_j
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_a_i
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: I_kc
      real(dp), dimension(:,:), allocatable :: I_ji
      real(dp), dimension(:,:), allocatable :: I_ab
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kc_jb
      real(dp), dimension(:,:), allocatable :: g_ck_bi
      real(dp), dimension(:,:), allocatable :: g_ck_aj
      real(dp), dimension(:,:), allocatable :: g_ai_ck
      real(dp), dimension(:,:), allocatable :: L_kc_bj
      real(dp), dimension(:,:), allocatable :: L_kc_jb
      real(dp), dimension(:,:), allocatable :: L_jc_kb
      real(dp), dimension(:,:), allocatable :: u_ai_kc
!
!     Indices
!
      integer(i15) :: a, i, b, c, j, k
      integer(i15) :: jb, kc, jc, kb, bj, ai, ck, ak, ci, bi, aj, cj
!
!     :: Term 1: 2 L_kcjb * c_bj * (2 t^ac_ik - t^ac_ki)  ::
!
!     Construct L_kcjb = 2 g_kc_jb - g_kb_jc ordered as
!
!               L_kc_bj = 2 g_kc_jb - g_kb_jc
!                 1234        1243      1342
!
      call mem%alloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_kc_bj, (wf%n_o)*(wf%n_v), (wf%n_v)*(wf%n_o))
      L_kc_bj = zero
!
      call wf%get_ovov(g_kc_jb)
!
      call add_1243_to_1234(two, g_kc_jb, L_kc_bj, (wf%n_o), (wf%n_v), (wf%n_v), (wf%n_o))
      call add_1342_to_1234(-one, g_kc_jb, L_kc_bj, (wf%n_o), (wf%n_v), (wf%n_v), (wf%n_o))
!
      call mem%dealloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     I_kc = sum_jb L_kcjb * c_bj = sum_jb L_kc_bj * c_bj
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
!     u_ai_ck = 2 t^ac_ik - t^ac_ki = - (2 g_ai_ck - g_ak_ci)/ε^{ac}_{ik}
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
                                          /(eps_v(a) + eps_v(c) &
                                          - eps_o(i) - eps_o(k))
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ai_ck, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_a_i = rho_a_i + sum_ck u^ac_ik I_kc = rho_a_i + sum_ck u_ai_kc I_kc
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
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
!     :: Term 2: L_kcjb t^cb_ki c_aj ::
!
!     Construct L_kcjb = 2 g_kc_jb - g_kb_jc ordered as
!
!               L_jc_kb = 2 g_kc_jb - g_kb_jc
!                 1234        3214      3412
!
      call mem%alloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_jc_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_jc_kb = zero
!
      call wf%get_ovov(g_kc_jb)
!
      call add_3214_to_1234(two, g_kc_jb, L_jc_kb, (wf%n_o), (wf%n_v), (wf%n_o), (wf%n_v))
      call add_3412_to_1234(-one, g_kc_jb, L_jc_kb, (wf%n_o), (wf%n_v), (wf%n_o), (wf%n_v))
!
      call mem%dealloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call mem%alloc(g_ck_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%get_vovo(g_ck_bi)
!
!     t_ck_bi = t^cb_ki = - g_ck_bi/ε^{cb}_{ik}
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do i = 1, wf%n_o
!
                  ck = (wf%n_v)*(k-1) + c
                  bi = (wf%n_v)*(i-1) + b
!
                  g_ck_bi(ck,bi) = - g_ck_bi(ck,bi) &
                                          /(wf%fock_diagonal(c + wf%n_o, 1) &
                                          + wf%fock_diagonal(b + wf%n_o, 1) &
                                          - wf%fock_diagonal(i, 1) - wf%fock_diagonal(k, 1))
               enddo
            enddo
         enddo
      enddo
!
!     Intermediat I_ji = L_kcjb t^cb_ki = L_j_ckb g_ckb_i
!
      call mem%alloc(I_ji, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',                    &
                  (wf%n_o),                   &
                  (wf%n_o),                   &
                  (wf%n_v)*(wf%n_o)*(wf%n_v), &
                  one,                        &
                  L_jc_kb,                    &
                  (wf%n_o),                   &
                  g_ck_bi,                    &
                  (wf%n_v)*(wf%n_o)*(wf%n_v), &
                  zero,                       &
                  I_ji,                       &
                  (wf%n_o))
!
      call mem%dealloc(L_jc_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(g_ck_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_a_i = rho_a_i - c_a_j I_ji = rho_a_i - c_b_j I_ji
!
      call dgemm('N', 'N',  &
                  (wf%n_v), &
                  (wf%n_o), &
                  (wf%n_o), &
                  one,      &
                  c_b_j,    &
                  (wf%n_v), &
                  I_ji,     &
                  (wf%n_o), &
                  -one,     &
                  rho_a_i,  &
                  (wf%n_v))
!
      call mem%dealloc(I_ji, wf%n_o, wf%n_o)
!
!
!     :: Term 3: L_kcjb t^ca_kj c_bi ::
!
!     Construct L_kcjb = 2 g_kc_jb - g_kb_jc ordered as
!
!     L_kc_jb = 2 g_kc_jb - g_kb_jc
!       1234        1234      1432
!
      call mem%alloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_kc_jb)
!
      L_kc_jb = two * g_kc_jb
!
      call add_1432_to_1234(-one, g_kc_jb, L_kc_jb, (wf%n_o), (wf%n_v), (wf%n_o), (wf%n_v))
!
      call mem%dealloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     t_ak_cj = t^ca_kj = - g_ck_aj/ε^{ca}_{jk}
!
      call mem%alloc(g_ck_aj,(wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%get_vovo(g_ck_aj)
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do a = 1, wf%n_v
               do j = 1, wf%n_o
!
                  ck = (wf%n_v)*(k-1) + c
                  aj = (wf%n_v)*(j-1) + a
                  ak = (wf%n_v)*(k-1) + a
                  cj = (wf%n_v)*(j-1) + c
!
                  g_ck_aj(ak,cj) = - g_ck_aj(ck,aj) &
                                    /(wf%fock_diagonal(a + wf%n_o, 1) &
                                    + wf%fock_diagonal(c + wf%n_o, 1) &
                                    - wf%fock_diagonal(j, 1) - wf%fock_diagonal(k, 1))
               enddo
            enddo
         enddo
      enddo
!
!     Intermediat I_ab = L_kcjb t^ac_cj = t_ak_cj L_kc_jb
!
      call mem%alloc(I_ab, (wf%n_v), (wf%n_o))
!
      call dgemm('N', 'N',                    &
                  (wf%n_v),                   &
                  (wf%n_v),                   &
                  (wf%n_o)*(wf%n_v)*(wf%n_o), &
                  one,                        &
                  g_ck_aj,                    &
                  (wf%n_v),                   &
                  L_kc_jb,                    &
                  (wf%n_o)*(wf%n_v)*(wf%n_o), &
                  zero,                       &
                  I_ab,                       &
                  (wf%n_v))
!
      call mem%dealloc(L_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(g_ck_aj,(wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_a_i = rho_a_i - I_ab c_b_i = rho_a_i - I_ab c_b_j
!
      call dgemm('N', 'N',  &
                  (wf%n_v), &
                  (wf%n_o), &
                  (wf%n_v), &
                  one,      &
                  I_ab,     &
                  (wf%n_v), &
                  c_b_j,    &
                  (wf%n_v), &
                  -one,     &
                  rho_a_i,  &
                  (wf%n_v))
!
      call mem%dealloc(I_ab, (wf%n_v), (wf%n_o))
!
!
   end subroutine jacobian_cc2_B1_cc2
!
!
   module subroutine effective_jacobian_cc2_a1_cc2(wf, omega, rho_a_i, c_c_j, eps_o, eps_v)
!!
!!    Jacobian CC2 E1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai^E1 = - L_kijb  (g_akbc * c_cj + g_bjac * c_ck) /(omega - ε_akbj) * 1/(1 + delta_ak,bj)
!!
!!    Every term will be done separately due to different batching
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Sent to the routine
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_a_i
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_c_j
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp) :: omega
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_ak_bc
      real(dp), dimension(:,:), allocatable :: g_bj_ac
      real(dp), dimension(:,:), allocatable :: g_jb_ki
      real(dp), dimension(:,:), allocatable :: g_ji_kb
      real(dp), dimension(:,:), allocatable :: L_jb_ki
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: I_ak_bj
      real(dp), dimension(:,:), allocatable :: I_bj_ak
!
!     Indices
!
      integer(i15) :: i, j, k, a, b, c
      integer(i15) :: ak, bj, aj, bk, ji, kb, jb, ki
!
!     :: Term 1: - L_kijb * g_akbc * c_cj /(omega - ε_akbj) * 1/(1 + delta_ak,bj)  ::
!
!     I_ak_bj = sum_c g_akbc * c_cj = sum_c g_akb_c * c_c_j
!
      call mem%alloc(g_ak_bc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_v))
      call mem%alloc(I_ak_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%get_vovv(g_ak_bc)
!
      call dgemm('N', 'N',                   &
                  (wf%n_v)**2*(wf%n_o),      &
                  (wf%n_o),                  &
                  (wf%n_v),                  &
                  one,                       &
                  g_ak_bc,                   &
                  (wf%n_v)**2*(wf%n_o),      &
                  c_c_j,                     &
                  (wf%n_v),                  &
                  zero,                      &
                  I_ak_bj,                   &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_ak_bc, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     Reordering of I_ak_bj as aj,bk and scaling by 1/(omega - ε_akbj) * 1/(1 + delta_ak,bj)
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
!
            ak = wf%n_v*(k-1) + a
            I_ak_bj(ak,ak) = half*I_ak_bj(ak,ak)
!
               do b = 1, wf%n_v
                  do j = 1, wf%n_o
!
                     bj = wf%n_v*(j-1) + b
                     aj = wf%n_v*(j-1) + a
                     bk = wf%n_v*(k-1) + b
!
                     I_ak_bj(aj,bk) = - I_ak_bj(ak,bj)&
                                       /(omega - eps_v(a) - eps_v(b) &
                                       + eps_o(j) + eps_o(k))
!
                  enddo
               enddo
         enddo
      enddo
!
!     Construct L_kijb = 2 g_ki_jb - g_kb_ji ordered as
!
!               L_jb_ki = L_ki_jb = 2 g_jb_ki - g_ji_kb
!
      call mem%alloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%alloc(g_ji_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%alloc(L_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call wf%get_ovoo(g_jb_ki)
      call wf%get_ooov(g_ji_kb)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do i = 1, wf%n_o
!
                  jb = wf%n_o*(b - 1) + j
                  ki = wf%n_o*(i - 1) + k
                  ji = wf%n_o*(i - 1) + j
                  kb = wf%n_o*(b - 1) + k
!
                  L_jb_ki(jb,ki) = two * g_jb_ki(jb,ki) - g_ji_kb(ji,kb)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%dealloc(g_ji_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
!     rho_a_i = rho_a_i - sum_jbk I_a_jbk * L_jbk_i
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  I_ak_bj,              &
                  (wf%n_v),             &
                  L_jb_ki,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_a_i,              &
                  (wf%n_v))
!
      call mem%dealloc(I_ak_bj, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call mem%dealloc(L_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
!     :: Term 2: - L_kijb * g_bjac * c_ck /(omega - ε_akbj) * 1/(1 + delta_ak,bj)  ::
!
!     I_bj_ak = sum_c g_bjac * c_ck = sum_c g_bja_c * c_c_k
!
      call mem%alloc(g_bj_ac, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_v))
      call mem%alloc(I_bj_ak, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call wf%get_vovv(g_bj_ac)
!
      call dgemm('N', 'N',                   &
                  (wf%n_v)**2*(wf%n_o),      &
                  (wf%n_o),                  &
                  (wf%n_v),                  &
                  one,                       &
                  g_bj_ac,                   &
                  (wf%n_v)**2*(wf%n_o),      &
                  c_c_j,                     &
                  (wf%n_v),                  &
                  zero,                      &
                  I_bj_ak,                   &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_bj_ac, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_v))
!
!     Reordering of I_bj_ak as aj,bk and scaling by 1/(omega - ε_akbj) * 1/(1 + delta_ak,bj)
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
!
            ak = wf%n_v*(k-1) + a
            I_bj_ak(ak,ak) = half*I_bj_ak(ak,ak)
!
               do b = 1, wf%n_v
                  do j = 1, wf%n_o
!
                     bj = wf%n_v*(j-1) + b
                     aj = wf%n_v*(j-1) + a
                     bk = wf%n_v*(k-1) + b
!
                     I_bj_ak(aj,bk) = - I_bj_ak(bj,ak)&
                                       /(omega - eps_v(a) - eps_v(b) &
                                       + eps_o(j) + eps_o(k))
!
                  enddo
               enddo
         enddo
      enddo
!
!     Due to different batching L has to be reconstructed
!
!     Construct L_kijb = 2 g_ki_jb - g_kb_ji ordered as
!
!               L_jb_ki = L_ki_jb = 2 g_jb_ki - g_ji_kb
!
      call mem%alloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%alloc(g_ji_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%alloc(L_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
      call wf%get_ovoo(g_jb_ki)
      call wf%get_ooov(g_ji_kb)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do i = 1, wf%n_o
!
                  jb = wf%n_o*(b - 1) + j
                  ki = wf%n_o*(i - 1) + k
                  ji = wf%n_o*(i - 1) + j
                  kb = wf%n_o*(b - 1) + k
!
                  L_jb_ki(jb,ki) = two * g_jb_ki(jb,ki) - g_ji_kb(ji,kb)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
      call mem%dealloc(g_ji_kb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
!     rho_a_i = rho_a_i - sum_jbk I_a_jbk * L_jbk_i
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  I_bj_ak,              &
                  (wf%n_v),             &
                  L_jb_ki,              &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_a_i,              &
                  (wf%n_v))
!
      call mem%dealloc(I_bj_ak, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_jb_ki, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_o))
!
   end subroutine effective_jacobian_cc2_a1_cc2
!
!
   module subroutine effective_jacobian_cc2_b1_cc2(wf, omega, rho_a_i, c_a_i, eps_o, eps_v)
!!
!!    Effective jacobian B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_a_i
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_a_i
!
      real(dp), dimension(wf%n_o), intent(in)  :: eps_o
      real(dp), dimension(wf%n_v), intent(in)  :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_lkai, X_ckai
!
      integer(i15) :: a, c, i, k
!
!     Construct X_c_kai = sum_l c_c_l g_l_kai
!
      call mem%alloc(g_lkai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_lkai)
!
      call mem%alloc(X_ckai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_a_i,               & ! c_c_l
                  wf%n_v,              &
                  g_lkai,              & ! g_l_kai
                  wf%n_o,              &
                  zero,                &
                  X_ckai,              &  ! X_c_kai
                  wf%n_v)
!
      call mem%dealloc(g_lkai, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     Divide X_ckai by 1/(1 + delta_ai,ck)*1/(-epsilon_aick + omega)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            X_ckai(a, i, a, i) = X_ckai(a, i, a, i)/two
!
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  X_ckai(c, k, a, i) = X_ckai(c, k, a, i)/(eps_o(k) + eps_o(i) &
                                                         - eps_v(a) - eps_v(c) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
!     Add all four terms while considering batched indices i and k
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  rho_a_i(a, i) = rho_a_i(a, i) - two*X_ckai(c, k, a, i)*wf%fock_ia(k, c) &
                                                + X_ckai(a, k, c, i)*wf%fock_ia(k, c)
!
                  rho_a_i(a, k) = rho_a_i(a, k)  - two*X_ckai(a, k, c, i)*wf%fock_ia(i, c) &
                                                 + X_ckai(c, k, a, i)*wf%fock_ia(i, c)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_ckai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine effective_jacobian_cc2_b1_cc2
!
!
end submodule jacobian
