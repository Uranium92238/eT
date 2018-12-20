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
!!    rho_ai^B1 = 2 L_kcjb c_bj (2 t^ac_ik - t^ac_ki)
!!                - L_kcjb t^cb_ki c_aj - L_kcjb t^ca_kj c_bi
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
      real(dp), dimension(:,:), allocatable :: I_ji
!
!     Integrals
!
      real(dp), dimension(:,:), allocatable :: g_kc_jb
      real(dp), dimension(:,:), allocatable :: g_ck_bi
      real(dp), dimension(:,:), allocatable :: g_ai_ck
      real(dp), dimension(:,:), allocatable :: L_kc_bj
      real(dp), dimension(:,:), allocatable :: L_jb_kc
      real(dp), dimension(:,:), allocatable :: u_ai_kc
!
!     Indices
!
      integer(i15) :: b, c, j, k, jb, kc, jc, kb, bj, a, i, ai, ck, ak, ci, bi
!
!     :: Term 1: 2 L_kcjb * c_bj * sum_bj (2 t^ac_ik - t^ac_ki)  ::
!
!
!     L_kc_bj = L_kcjb = 2 g_kc_jb - g_kb_jc
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
!
!     L_jb_kc = L_kcjb = 2 g_kc_jb - g_kb_jc
!
      call mem%alloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      L_jb_kc = zero
!
      call wf%get_ovov(g_kc_jb)
!
      call add_3412_to_1234(two, g_kc_jb, L_jb_kc, (wf%n_o), (wf%n_v), (wf%n_o), (wf%n_v))
      call add_1432_to_1234(-one, g_kc_jb, L_jb_kc, (wf%n_o), (wf%n_v), (wf%n_o), (wf%n_v))
!
      call mem%dealloc(g_kc_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call mem%alloc(g_ck_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
      call wf%get_vovo(g_ck_bi)
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
!     Intermediat I_ji = L_kcjb t^cb_ki = L_jb_kc g_ck_bi
!
      call mem%alloc(I_ji, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',                    &
                  (wf%n_o),                   &
                  (wf%n_o),                   &
                  (wf%n_v)*(wf%n_o)*(wf%n_v), &
                  one,                        &
                  L_jb_kc,                    &
                  (wf%n_o),                   &
                  g_ck_bi,                    &
                  (wf%n_v)*(wf%n_o)*(wf%n_v), &
                  zero,                       &
                  I_ji,                       &
                  (wf%n_o))
!
      call mem%dealloc(L_jb_kc, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%dealloc(g_ck_bi, (wf%n_v)*(wf%n_o), (wf%n_v)*(wf%n_o))
!
!     rho_a_i = rho_a_i + c_a_j I_ji = rho_a_i + c_b_j I_ji
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
   end subroutine
!
!
end submodule jacobian
