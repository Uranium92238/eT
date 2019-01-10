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
   subroutine effective_jacobian_transformation_cc2(wf, omega, c)
!!
!!    Effective jacobian transformation
!!    Written by Eirik F. Kjønstad and Sarai Dery Folkestad
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), intent(in) :: omega
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_a_i
!
      real(dp), dimension(:,:), allocatable :: rho_a_i
!
      real(dp), dimension(:), allocatable :: eps_o
      real(dp), dimension(:), allocatable :: eps_v
!
      integer(i15) :: i, a, ai ! Index
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
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call mem%alloc(eps_o, wf%n_o)
      call mem%alloc(eps_v, wf%n_v)
!
      eps_o = wf%fock_diagonal(1:wf%n_o,1)
      eps_v = wf%fock_diagonal(wf%n_o + 1 : wf%n_mo, 1)
!
      call wf%jacobian_cc2_a1(rho_a_i, c_a_i)
      call wf%jacobian_cc2_b1(rho_a_i, c_a_i, eps_o, eps_v)
!
      call wf%effective_jacobian_cc2_a1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_b1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_c1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_d1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_e1(omega, rho_a_i, c_a_i, eps_o, eps_v)
      call wf%effective_jacobian_cc2_f1(omega, rho_a_i, c_a_i, eps_o, eps_v)
!
      call dcopy(wf%n_amplitudes, rho_a_i, 1, c, 1)
!
      call mem%dealloc(c_a_i, wf%n_v, wf%n_o)
      call mem%dealloc(rho_a_i, wf%n_v, wf%n_o)
!
   end subroutine effective_jacobian_transformation_cc2
!
!
   module subroutine jacobian_cc2_a1_cc2(wf, rho_ai, c_bj)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai =+ sum_bj (2 g_aijb - g_abji) * c_bj
!!
!!    Separate calculation of both terms due to batching
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(:,:), allocatable :: c_jb
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aijb
      real(dp), dimension(:,:,:,:), allocatable :: g_abji
!
!     Indices
!
      integer(i15) :: a, b, i, j
!
!     Explicit reordering of c_bj
!
      call mem%alloc(c_jb, wf%n_o, wf%n_v)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            c_jb(j, b) = c_bj(b, j)
!
         enddo
      enddo
!
      call mem%alloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_voov(g_aijb,      &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
!     rho_a_i = rho_a_i + sum_bj 2 g_aijb * c_bj
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  g_aijb,            &
                  (wf%n_v)*(wf%n_o), &
                  c_jb,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(g_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_vvoo(g_abji,      &
                        1, wf%n_v,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_o)
!
!     Sort g_abji(a,b,j,i) as g_abji(a,i,j,b)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  g_aijb(a,i,j,b) = g_abji(a,b,j,i)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_abji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     rho_a_i = rho_a_i - sum_bj g_aijb * c_jb
!
      call dgemm('N', 'N',           &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  g_aijb,            &
                  (wf%n_v)*(wf%n_o), &
                  c_jb,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  rho_ai,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(c_jb, wf%n_o, wf%n_v)
!
   end subroutine jacobian_cc2_a1_cc2
!
!
   module subroutine jacobian_cc2_b1_cc2(wf, rho_ai, c_bj, eps_o, eps_v)
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c_bj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: X_kc, X_ji, X_ab
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcjb, g_ckbi, g_ckaj, g_aick
      real(dp), dimension(:,:,:,:), allocatable :: L_kcbj, L_kcjb, L_jckb
      real(dp), dimension(:,:,:,:), allocatable :: u_aikc, t_akcj
!
!     Indices
!
      integer(i15) :: a, i, b, c, j, k
!
!     :: Term 1: 2 L_kcjb * c_bj * (2 t^ac_ik - t^ac_ki)  ::
!
!     Construct L_kcjb = 2 g_kc_jb - g_kb_jc ordered as
!
!               L_kc_bj = 2 g_kc_jb - g_kb_jc
!                 1234        1243      1342
!
      call mem%alloc(g_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcjb,      &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call mem%alloc(L_kcbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  L_kcbj(k,c,b,j) = two * g_kcjb(k,c,j,b) - g_kcjb(k,b,j,c)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     X_kc = sum_jb L_kcjb * c_bj = sum_jb L_kc_bj * c_bj
!
      call mem%alloc(X_kc, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',           &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  L_kcbj,            &
                  (wf%n_o)*(wf%n_v), &
                  c_bj,              &
                  (wf%n_v)*(wf%n_o), &
                  zero,              &
                  X_kc,              &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_kcbj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     u_aick = 2 t^acik - t^acki = - (2 g_aick - g_akci)/ε^{ac}_{ik}
!
      call mem%alloc(g_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_vovo(g_aick,      &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o)
!
      call mem%alloc(u_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            do c = 1, wf%n_v
               do k = 1, wf%n_o
!
                  u_aikc(a,i,k,c) = - (two*g_aick(a,i,c,k)- g_aick(a,k,c,i))&
                                       /(eps_v(a) + eps_v(c) - eps_o(i) - eps_o(k))
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai = rho_ai + sum_ck u^ac_ik X_kc = rho_ai + sum_ck u_aikc X_kc
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  u_aikc,              &
                  (wf%n_v)*(wf%n_o),   &
                  X_kc,                &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_kc, wf%n_o, wf%n_v)
      call mem%dealloc(u_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!
!     :: Term 2: L_kcjb t^cb_ki c_aj ::
!
!     Construct L_kcjb = 2 g_kcjb - g_kbjc ordered as
!
!               L_jckb = 2 g_kcjb - g_kbjc
!
!     Batching over k and b
!
      call mem%alloc(g_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(L_jckb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcjb,      &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
    ! call mem%alloc(g_kbjc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
    ! call wf%get_ovov(g_kbjc,      &
    !                   1, wf%n_o,  &
    !                   1, wf%n_v,  &
    !                   1, wf%n_o,  &
    !                   1, wf%n_v)
!
      do j = 1, wf%n_o
         do c = 1, wf%n_v
            do k = 1, wf%n_o
               do b = 1, wf%n_v
!
                  L_jckb(j,c,k,b) = two * g_kcjb(k,c,j,b) - g_kcjb(j,c,k,b)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
   ! call mem%dealloc(g_kbjc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Construct t_ckbi = - g_ckbi/ε^{cb}_{ik}
!
      call mem%alloc(g_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_vovo(g_ckbi,      &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o)
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do i = 1, wf%n_o
!
                  g_ckbi(c,k,b,i) = - g_ckbi(c,k,b,i) &
                                    /(eps_v(c) + eps_v(b)- eps_o(i) - eps_o(k))
               enddo
            enddo
         enddo
      enddo
!
!     Intermediat X_ji = sum_ckb L_kcjb t^cb_ki = sum_ckb L_jckb g_ckbi
!
      call mem%alloc(X_ji, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',                    &
                  (wf%n_o),                   &
                  (wf%n_o),                   &
                  (wf%n_v)**2*(wf%n_o),       &
                  one,                        &
                  L_jckb,                     &
                  (wf%n_o),                   &
                  g_ckbi,                     &
                  (wf%n_v)**2*(wf%n_o),       &
                  zero,                       &
                  X_ji,                       &
                  (wf%n_o))
!
      call mem%dealloc(L_jckb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(g_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai = rho_ai - c_aj X_ji = rho_ai - c_bj X_ji
!
      call dgemm('N', 'N',  &
                  (wf%n_v), &
                  (wf%n_o), &
                  (wf%n_o), &
                  -one,     &
                  c_bj,     & ! c_a_j
                  (wf%n_v), &
                  X_ji,     &
                  (wf%n_o), &
                  one,      &
                  rho_ai,   &
                  (wf%n_v))
!
      call mem%dealloc(X_ji, wf%n_o, wf%n_o)
!
!
!     :: Term 3: L_kcjb t^ca_kj c_bi ::
!
!     Construct L_kcjb = 2 g_kcjb - g_kbjc ordered as
!
!     L_kc_jb = 2 g_kc_jb - g_kb_jc
!
!     Batching over k and j
!
      call mem%alloc(g_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(L_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcjb,      &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  L_kcjb(k,c,j,b) = two * g_kcjb(k,c,j,b) - g_kcjb(k,b,j,c)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Construct t_akcj = - g_ckaj/ε^{ca}_{jk}
!
      call mem%alloc(g_ckaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_vovo(g_ckaj,      &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o)
!
      call mem%alloc(t_akcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
            do a = 1, wf%n_v
               do j = 1, wf%n_o
!
                  t_akcj(a,k,c,j) = - g_ckaj(c,k,a,j) &
                                  /(eps_v(a) + eps_v(c)- eps_o(j) - eps_o(k))
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ckaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Intermediat X_ab = L_kcjb t^ac_cj = t_akcj L_kcjb
!
      call mem%alloc(X_ab, (wf%n_v), (wf%n_v))
!
      call dgemm('N', 'N',                    &
                  (wf%n_v),                   &
                  (wf%n_v),                   &
                  (wf%n_o)*(wf%n_v)*(wf%n_o), &
                  one,                        &
                  t_akcj,                     &
                  (wf%n_v),                   &
                  L_kcjb,                     &
                  (wf%n_o)*(wf%n_v)*(wf%n_o), &
                  zero,                       &
                  X_ab,                       &
                  (wf%n_v))
!
      call mem%dealloc(t_akcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_kcjb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     rho_ai = rho_ai - X_ab c_bi
!
      call dgemm('N', 'N',  &
                  (wf%n_v), &
                  (wf%n_o), &
                  (wf%n_v), &
                  -one,     &
                  X_ab,     &
                  (wf%n_v), &
                  c_bj,     & ! c_b_i
                  (wf%n_v), &
                  one,      &
                  rho_ai ,  &
                  (wf%n_v))
!
      call mem%dealloc(X_ab, (wf%n_v), (wf%n_v))
!
   end subroutine jacobian_cc2_b1_cc2
!
!
   module subroutine effective_jacobian_cc2_a1_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 A1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    rho_ai =+ F_kc * (-eps_ai,ck + w)^-1 * (1 + delta_ai,ck)^-1 * (2 g_aicd c_dk + 2 g_ckad c_di - g_akcd c_di - g_ciad c_dk)
!!           =+ F_kc * (-eps_ai,ck + w)^-1 * (1 + delta_ai,ck)^-1 * (2 X_aick - X_akci + 2 X_ckai - X_ciak)
!!           =+ F_kc * (Y_aick + Y_ckai)
!!
!!    The term is calculated in batches over the a and c indices.
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aick
      real(dp), dimension(:,:,:,:), allocatable :: Y_aick
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aicd
!
      real(dp), dimension(:,:), allocatable :: F_ck
!
      integer(i15) :: a, i, c, k
!
!     Construct X_aick = sum_d g_aicd c_dk
!
      call mem%alloc(g_aicd, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      call wf%get_vovv(g_aicd,     &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_v)
!
      call mem%alloc(X_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',          &
                  wf%n_o*wf%n_v**2, &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  g_aicd,           &
                  wf%n_o*wf%n_v**2, &
                  c_ai,             & ! c_dk
                  wf%n_v,           &
                  zero,             &
                  X_aick,           &
                  wf%n_o*wf%n_v**2)
!
      call mem%dealloc(g_aicd, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
!     Scale X: Y_aick = (-eps_ai,ck + w)^-1 * (1 + delta_ai,ck)^-1 * (2 X_aick - X_akci)
!
      call mem%alloc(Y_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  Y_aick(a, i, c, k) = (two*X_aick(a, i, c, k) - X_aick(a, k, c, i))/&
                        (- eps_v(a) - eps_v(c) + eps_o(i) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            Y_aick(a, i, a, i) = Y_aick(a, i, a, i)/two
!
         enddo
      enddo
!
      call mem%dealloc(X_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder occ-vir Fock matrix, then compute & add the term,
!     2 F_ck * (Y_aick + Y_ckai) to the transformed vector rho_ai
!
      call mem%alloc(F_ck, wf%n_v, wf%n_o)
!
      do k = 1, wf%n_o
         do c = 1, wf%n_v
!
            F_ck(c, k) = wf%fock_ia(k, c)
!
         enddo
      enddo
!
      call dgemm('N', 'N',       &
                  wf%n_o*wf%n_v, &
                  1,             &
                  wf%n_o*wf%n_v, &
                  one,           &
                  Y_aick,        &
                  wf%n_o*wf%n_v, &
                  F_ck,          &
                  wf%n_o*wf%n_v, &
                  one,           &
                  rho_ai,        &
                  wf%n_o*wf%n_v)
!
      call dgemm('T', 'N',       &
                  wf%n_o*wf%n_v, &
                  1,             &
                  wf%n_o*wf%n_v, &
                  one,           &
                  Y_aick,        &
                  wf%n_o*wf%n_v, &
                  F_ck,          &
                  wf%n_o*wf%n_v, &
                  one,           &
                  rho_ai,        &
                  wf%n_o*wf%n_v)
!
      call mem%dealloc(F_ck, wf%n_v, wf%n_o)
      call mem%dealloc(Y_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine effective_jacobian_cc2_a1_cc2
!
!
   module subroutine effective_jacobian_cc2_b1_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Effective jacobian B1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
!!    Effective B1 = - 2 sum_{kcl} F_kc (1/Δ_{ai,ck})*(1/(ε_{aick} + ω)) * (g_ailk c_cl + g_ckli c_al)
!!                     + sum_{kcl} F_kc (1/Δ_{ak,ci})*(1/(ε_{akci} + ω)) * (g_akli c_cl + g_cilk c_al)
!!                 =   2 sum_{kcl} F_kc (- 2*X_ckai - 2*X_aick + X_ciak + X_akci)
!!
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in)  :: eps_o
      real(dp), dimension(wf%n_v), intent(in)  :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_lkai, X_ckai
!
      integer(i15) :: a, c, i, k
!
      integer(i15) :: req0, req1_i, req1_k, req2
!
      integer(i15) :: current_i_batch, current_k_batch
!
      type(batching_index) :: batch_i, batch_k
!
      req0 = 0
!
      req1_i = (wf%integrals%n_J)*(wf%n_v)
      req1_k = (wf%integrals%n_J)*(wf%n_o)
!
      req2 = 2*(wf%n_o)*(wf%n_v)
!
      call batch_i%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_k, req0, req1_i, req1_k, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
!           Construct X_c_kai = sum_l c_c_l g_l_kai
!
            call mem%alloc(g_lkai, wf%n_o, batch_k%length, wf%n_v, batch_i%length)
!
            call wf%get_oovo(g_lkai,                        &
                              1, wf%n_o,                    &
                              batch_k%first, batch_k%last,  &
                              1, wf%n_v,                    &
                              batch_i%first, batch_i%last)
!
            call mem%alloc(X_ckai, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
            call dgemm('N', 'N',                                  &
                        wf%n_v,                                   &
                        (batch_i%length)*(batch_k%length)*wf%n_v, &
                        wf%n_o,                                   &
                        one,                                      &
                        c_ai,                                     & ! c_c_l
                        wf%n_v,                                   &
                        g_lkai,                                   & ! g_l_kai
                        wf%n_o,                                   &
                        zero,                                     &
                        X_ckai,                                   &  ! X_c_kai
                        wf%n_v)
!
            call mem%dealloc(g_lkai, wf%n_o, batch_k%length, wf%n_v, batch_i%length)
!
!           Divide X_ckai by 1/(1 + delta_ai,ck)*1/(-epsilon_aick + omega)
!
            do i = 1, batch_i%length
               do a = 1, wf%n_v
!
                  X_ckai(a, i, a, i) = X_ckai(a, i, a, i)/two
!
                  do k = 1, batch_k%length
                     do c = 1, wf%n_v
!
                        X_ckai(c, k, a, i) = X_ckai(c, k, a, i)/(eps_o(k + batch_k%first - 1) &
                                                               + eps_o(i + batch_i%first - 1) &
                                                               - eps_v(a) - eps_v(c) + omega)
!
                     enddo
                  enddo
               enddo
            enddo
!
!           Add all four terms while considering batched indices i and k
!
            do i = 1, batch_i%length
               do a = 1, wf%n_v
                  do k = 1, batch_k%length
                     do c = 1, wf%n_v
!
                        rho_ai(a, i + batch_i%first - 1) = rho_ai(a, i + batch_i%first - 1)                                &
                                                            - two*X_ckai(c, k, a, i)*(wf%fock_ia(k + batch_k%first - 1, c))&
                                                            + X_ckai(a, k, c, i)*(wf%fock_ia(k + batch_k%first - 1, c))
!
                        rho_ai(a, k + batch_k%first - 1) = rho_ai(a, k + batch_k%first - 1)                                & !(k <-> i)
                                                            - two*X_ckai(a, k, c, i)*(wf%fock_ia(i + batch_i%first - 1, c))&
                                                            + X_ckai(c, k, a, i)*(wf%fock_ia(i + batch_i%first - 1, c))
!
                     enddo
                  enddo
               enddo
            enddo
!
            call mem%dealloc(X_ckai, wf%n_v, batch_k%length, wf%n_v, batch_i%length)
!
         enddo
      enddo
!
   end subroutine effective_jacobian_cc2_b1_cc2
!
!
   module subroutine effective_jacobian_cc2_c1_cc2(wf, omega, rho_ai, c_cj, eps_o, eps_v)
!!
!!    Jacobian CC2 C1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai^C1 =+ - L_kijb  (g_akbc * c_cj + g_bjac * c_ck) (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1
!!              =+ - L_kijb  (X_akbj + X_bjak)
!!
!!    Every term will be done separately due to different batching (k,b and abi or jbi)
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Sent to the routine
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_cj
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_akbc, g_bjac, g_jbki, g_jikb
      real(dp), dimension(:,:,:,:), allocatable :: L_jbki
!
!     Intermediates
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akbj
      real(dp), dimension(:,:,:,:), allocatable :: X_bjak, X_ajbk
!
!     Indices
!
      integer(i15) :: i, j, k, a, b
!
!     :: Term 1: - L_kijb * g_akbc * c_cj (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1  ::
!
!     Construct X_akbj = sum_c g_akbc * c_cj
!
      call mem%alloc(g_akbc, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      call wf%get_vovv(g_akbc,     &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_v)
!
      call mem%alloc(X_akbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_o)*(wf%n_v)**2, &
                  (wf%n_o),             &
                  (wf%n_v),             &
                  one,                  &
                  g_akbc,               &
                  (wf%n_v)**2*(wf%n_o), &
                  c_cj,                 &
                  (wf%n_v),             &
                  zero,                 &
                  X_akbj,               &
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(g_akbc, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
!     Reordering of X_akbj as ajbk and scaling by (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1
!
      call mem%alloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
!
            X_akbj(a,k,a,k) = half*X_akbj(a,k,a,k)
!
         enddo
      enddo
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  X_ajbk(a,j,b,k) = X_akbj(a,k,b,j) &
                                 /(- eps_v(a) - eps_v(b) + eps_o(j) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_akbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct L_kijb = 2 g_kijb - g_kbji ordered as
!
!               L_jbki = L_ki_jb = 2 g_jbki - g_jikb
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ovoo(g_jbki,     &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_o)
!
      call wf%get_ooov(g_jikb,     &
                        1, wf%n_o, &
                        1, wf%n_o, &
                        1, wf%n_o, &
                        1, wf%n_v)
!

      call mem%alloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do i = 1, wf%n_o
!
                  L_jbki(j,b,k,i) = two * g_jbki(j,b,k,i) - g_jikb(j,i,k,b)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_ai = rho_ai - sum_jbk X_ajbk * L_jbki
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  X_ajbk,               &
                  (wf%n_v),             &
                  L_jbki,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               &
                  (wf%n_v))
!

      call mem%dealloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2: - L_kijb * g_bjac * c_ck (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1  ::
!
!     X_bjak = sum_c g_bjac * c_ck = sum_c g_bjac * c_ck
!
      call mem%alloc(g_bjac, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      call wf%get_vovv(g_bjac,      &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_v)
!
      call mem%alloc(X_bjak, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_v)**2*(wf%n_o), &
                  (wf%n_o),             &
                  (wf%n_v),             &
                  one,                  &
                  g_bjac,               &
                  (wf%n_v)**2*(wf%n_o), &
                  c_cj,                 &
                  (wf%n_v),             &
                  zero,                 &
                  X_bjak,               &
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(g_bjac, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
!     Reordering of X_bjak as aj,bk and scaling by (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1
!
      call mem%alloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
!
            X_bjak(a,k,a,k) = half*X_bjak(a,k,a,k)
!
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  X_ajbk(a,j,b,k) = X_bjak(b,j,a,k)&
                                    /(- eps_v(a) - eps_v(b) + eps_o(j) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_bjak, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Due to different batching L has to be reconstructed
!
!     Construct L_kijb = 2 g_ki_jb - g_kb_ji ordered as
!
!               L_jb_ki = L_ki_jb = 2 g_jb_ki - g_ji_kb
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ovoo(g_jbki,      &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_o)
!
      call wf%get_ooov(g_jikb,      &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call mem%alloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do i = 1, wf%n_o
!
                  L_jbki(j,b,k,i) = two * g_jbki(j,b,k,i) - g_jikb(j,i,k,b)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_ai = rho_ai - sum_jbk X_ajbk * L_jbki
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_v)*(wf%n_o)**2, &
                  -one,                 &
                  X_ajbk,               &
                  (wf%n_v),             &
                  L_jbki,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               &
                  (wf%n_v))
!
      call mem%dealloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine effective_jacobian_cc2_c1_cc2
!
!
   module subroutine effective_jacobian_cc2_d1_cc2(wf, omega, rho_ai, c_bl, eps_o, eps_v)
!!
!!    Jacobian CC2 D1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai =+ - L_kijb  (- g_aklj * c_bl - g_bjlk * c_al) (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1
!!           =+ - L_kijb  (X_bjak + X_akbj)
!!
!!    Every term will be done separately due to different batching kj, bk
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Sent to the routine
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_bl
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ljak, g_jbki, g_jikb, g_lkbj
      real(dp), dimension(:,:,:,:), allocatable :: L_jbki
!
!     Intermediates
!
      real(dp), dimension(:,:,:,:), allocatable :: X_akbj, X_bjak, X_ajbk
!
!     Indices
!
      integer(i15) :: i, j, k, a, b
!
!     :: Term 1: L_kijb * g_aklj * c_bl (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1  ::
!
!     Construct X_bjak = sum_l g_aklj * c_bl
!
      call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_ljak,     &
                        1, wf%n_o, &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_o)
!
      call mem%alloc(X_bjak, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o)**2*(wf%n_v), &
                  (wf%n_o),             &
                  one,                  &
                  c_bl,                 &
                  (wf%n_v),             &
                  g_ljak,               &
                  (wf%n_o),             &
                  zero,                 &
                  X_bjak,               &
                  (wf%n_v))
!
      call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reordering of X_bjak as ajbk and scaling by (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1
!
      call mem%alloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
!
            X_bjak(a,k,a,k) = half*X_bjak(a,k,a,k)
!
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  X_ajbk(a,j,b,k) = X_bjak(b,j,a,k) &
                                 /(- eps_v(a) - eps_v(b) + eps_o(j) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_bjak, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct L_kijb = 2 g_kijb - g_kbji ordered as
!
!               L_jbki = L_ki_jb = 2 g_jbki - g_jikb
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_ovoo(g_jbki,      &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_o)
!
      call mem%alloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_jikb,      &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call mem%alloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do i = 1, wf%n_o
!
                  L_jbki(j,b,k,i) = two * g_jbki(j,b,k,i) - g_jikb(j,i,k,b)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_ai = rho_ai + sum_jbk X_bjak(a,j,b,k) * L_jbki
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_ajbk,               &
                  (wf%n_v),             &
                  L_jbki,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               &
                  (wf%n_v))
!
      call mem%dealloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2: L_kijb * g_bjlk * c_al (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1  ::
!
!     Construct X_akbj = sum_c g_bjlk * c_al
!
      call mem%alloc(g_lkbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_lkbj,     &
                        1, wf%n_o, &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_o)
!
      call mem%alloc(X_akbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o)**2*(wf%n_v), &
                  (wf%n_o),             &
                  one,                  &
                  c_bl,                 & ! c_al
                  (wf%n_v),             &
                  g_lkbj,               &
                  (wf%n_o),             &
                  zero,                 &
                  X_akbj,               &
                  (wf%n_v))
!
      call mem%dealloc(g_lkbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reordering of X_akbj as ajbk and scaling by (omega - ε_akbj)^-1 * (1 + delta_ak,bj)^-1
!
      call mem%alloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do a = 1, wf%n_v
         do k = 1, wf%n_o
!
            X_akbj(a,k,a,k) = half*X_akbj(a,k,a,k)
!
            do b = 1, wf%n_v
               do j = 1, wf%n_o
!
                  X_ajbk(a,j,b,k) = X_akbj(a,k,b,j) &
                                 /(- eps_v(a) - eps_v(b) + eps_o(j) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_akbj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct L_kijb = 2 g_kijb - g_kbji ordered as
!
!               L_jbki = L_ki_jb = 2 g_jbki - g_jikb
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_ovoo(g_jbki,     &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_o)
!
      call mem%alloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_jikb,     &
                        1, wf%n_o, &
                        1, wf%n_o, &
                        1, wf%n_o, &
                        1, wf%n_v)
!
      call mem%alloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do i = 1, wf%n_o
!
                  L_jbki(j,b,k,i) = two * g_jbki(j,b,k,i) - g_jikb(j,i,k,b)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     rho_ai = rho_ai + sum_jbk X_bjak(a,j,b,k) * L_jbki
!
      call dgemm('N', 'N',                &
                  (wf%n_v),               &
                  (wf%n_o),               &
                  (wf%n_v)*(wf%n_o)**2,   &
                  one,                    &
                  X_ajbk,                 &
                  (wf%n_v),               &
                  L_jbki,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  one,                    &
                  rho_ai,                 &
                  (wf%n_v))
!
      call mem%dealloc(X_ajbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine effective_jacobian_cc2_d1_cc2
!
!
   module subroutine effective_jacobian_cc2_e1_cc2(wf, omega, rho_ai, c_dk, eps_o, eps_v)
!!
!!    Jacobian CC2 F1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
!!    Implicit calculation of the doubles vector
!!    rho_ai^E1 =+ L_abkc  (g_bicd * c_dk + g_ckbd * c_di) (omega - ε_bick)^-1 * (1 + delta_bi,ck)^-1
!!              =+ L_abkc  (X_bick + X_ckbi)
!!
!!    Every term will be done separately due to different batching (k,b and )
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
!     Sent to the routine
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_dk
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_bicd, g_ckbd, g_abkc, g_kbac
      real(dp), dimension(:,:,:,:), allocatable :: L_abkc
!
!     Intermediates
!
      real(dp), dimension(:,:,:,:), allocatable :: X_bick, X_ckbi, X_bkci
!
!     Indices
!
      integer(i15) :: a, b, c, i, k
!
!     :: Term 1:  L_abkc * g_bicd * c_dk (omega - ε_bick)^-1 * (1 + delta_bi,ck)^-1  ::
!
!     Construct X_bick = sum_d g_bicd * c_dk
!
      call mem%alloc(g_bicd, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      call wf%get_vovv(g_bicd,     &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_v)
!
      call mem%alloc(X_bick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_o)*(wf%n_v)**2, &
                  (wf%n_o),             &
                  (wf%n_v),             &
                  one,                  &
                  g_bicd,               &
                  (wf%n_v)**2*(wf%n_o), &
                  c_dk,                 &
                  (wf%n_v),             &
                  zero,                 &
                  X_bick,               &
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(g_bicd, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
!     Reordering of X_bick as bkci and scaling by (omega - ε_bick)^-1 * (1 + delta_bi,ck)^-1
!
      call mem%alloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do b = 1, wf%n_v
         do i = 1, wf%n_o
!
            X_bick(b,i,b,i) = half*X_bick(b,i,b,i)
!
            do c = 1, wf%n_v
               do k = 1, wf%n_o
!
                  X_bkci(b,k,c,i) = X_bick(b,i,c,k) &
                                 /(- eps_v(b) - eps_v(c) + eps_o(i) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_bick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct L_abkc = 2 g_abkc - g_kbac
!
      call mem%alloc(g_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vvov(g_abkc,     &
                        1, wf%n_v, &
                        1, wf%n_v, &
                        1, wf%n_o, &
                        1, wf%n_v)
!

      call mem%alloc(g_kbac, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
!
      call wf%get_ovvv(g_kbac,     &
                        1, wf%n_o, &
                        1, wf%n_v, &
                        1, wf%n_v, &
                        1, wf%n_v)
!

      call mem%alloc(L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do k = 1, wf%n_o
               do c = 1, wf%n_v
!
                  L_abkc(a,b,k,c) = two * g_abkc(a,b,k,c) - g_kbac(k,b,a,c)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(g_kbac, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
!
!     rho_ai = rho_ai + sum_bkc L_abkc * X_bkci
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  L_abkc,               &
                  (wf%n_v),             &
                  X_bkci,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  rho_ai,               &
                  (wf%n_v))
!
      call mem%dealloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2: L_abkc * g_ckbd * c_di (omega - ε_bick)^-1 * (1 + delta_bi,ck)^-1  ::
!
!     Construct X_ckbi = sum_d g_ckbd * c_di
!
      call mem%alloc(g_ckbd, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      call wf%get_vovv(g_ckbd,      &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_v)
!
      call mem%alloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_v)**2*(wf%n_o), &
                  (wf%n_o),             &
                  (wf%n_v),             &
                  one,                  &
                  g_ckbd,               &
                  (wf%n_v)**2*(wf%n_o), &
                  c_dk,                 & ! c_di
                  (wf%n_v),             &
                  zero,                 &
                  X_ckbi,               &
                  (wf%n_o)*(wf%n_v)**2)
!
      call mem%dealloc(g_ckbd, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
!     Reordering of X_ckbi as bkci and scaling by (omega - ε_bick)^-1 * (1 + delta_bi,ck)^-1
!
      call mem%alloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do c = 1, wf%n_v
         do k = 1, wf%n_o
!
            X_ckbi(c,k,c,k) = half*X_ckbi(c,k,c,k)
!
            do b = 1, wf%n_v
               do i = 1, wf%n_o
!
                  X_bkci(b,k,c,i) = X_ckbi(c,k,b,i)&
                                 /(- eps_v(b) - eps_v(c) + eps_o(i) + eps_o(k) + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai = rho_ai - sum_bkc L_abkc * X_bkci
!
      call dgemm('N', 'N',              &
                  (wf%n_v),             &
                  (wf%n_o),             &
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  L_abkc,               &
                  (wf%n_v),             &
                  X_bkci,               & ! X_bkci
                  (wf%n_o)*(wf%n_v)**2, &
                  one,                  &
                  rho_ai,               &
                  (wf%n_v))
!
      call mem%dealloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
   end subroutine effective_jacobian_cc2_e1_cc2
!
!
   module subroutine effective_jacobian_cc2_f1_cc2(wf, omega, rho_ai, c_ai, eps_o, eps_v)
!!
!!    Jacobian CC2 effective F1
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Linda Goletto, and Alexander Paul, Jan 2019
!!
!!    Effective F1 = - L_abkc ((1/Δ_{bi,ck})*(1/(ε_{bick} + ω)) * (g_lkbi c_cl + g_lick c_bl))
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
!
      real(dp), dimension(wf%n_o), intent(in) :: eps_o
      real(dp), dimension(wf%n_v), intent(in) :: eps_v
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ckbi, X_bkci, X_bick
      real(dp), dimension(:,:,:,:), allocatable :: L_abkc, g_ackb, g_lkbi, g_lick
!
      integer(i15) :: b, c, i, k
!
!     Term 1: - L_abkc (1/Δ_{bi,ck})*(1/(ε_{bick} + ω)) * (g_lkbi c_cl)
!
      call mem%alloc(g_lkbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_lkbi,      &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o)
!
      call mem%alloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_o,              &
                  one,                 &
                  c_ai,                & ! c_c_l
                  wf%n_v,              &
                  g_lkbi,              & ! g_l_kbi
                  wf%n_o,              &
                  zero,                &
                  X_ckbi,              & ! X_c_kbi
                  wf%n_v)
!
      call mem%dealloc(g_lkbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder and scale
!
      call mem%alloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do i = 1, wf%n_o
         do c = 1, wf%n_v
            do k = 1, wf%n_o
               do b = 1, wf%n_v
!
                  X_bkci(b, k, c, i) = X_ckbi(c, k, b, i)/( eps_o(i) + eps_o(k)  &
                                                          - eps_v(b) - eps_v(c)  &
                                                          + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      do k = 1, wf%n_o
         do b = 1, wf%n_v
!
            X_bkci(b, k, b, k) = X_bkci(b, k, b, k)/two
!
         enddo
      enddo
!
      call mem%dealloc(X_ckbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vvov(L_abkc,      &
                        1, wf%n_v,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call dscal((wf%n_v**3)*(wf%n_o), -two, L_abkc, 1)
!
      call mem%alloc(g_ackb, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vvov(g_ackb,      &
                        1, wf%n_v,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call add_1432_to_1234(one, g_ackb, L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ackb, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o*(wf%n_v**2),  &
                  one,                 &
                  L_abkc,              &
                  wf%n_v,              &
                  X_bkci,              &
                  wf%n_o*(wf%n_v**2),  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2: - L_abkc ((1/Δ_{bi,ck})*(1/(ε_{bick} + ω)) * (g_lick c_bl))
!
      call mem%alloc(g_lick, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_lick,      &
                        1, wf%n_o,  &
                        1, wf%n_o,  &
                        1, wf%n_v,  &
                        1, wf%n_o)
!
       call mem%alloc(X_bick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_o,              &
                  one,                 &
                  c_ai,                & ! c_b_l
                  wf%n_v,              &
                  g_lick,              & ! g_l_ick
                  wf%n_o,              &
                  zero,                &
                  X_bick,              & ! X_b_ick
                  wf%n_v)
!
      call mem%dealloc(g_lick, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder and scale
!
      call mem%alloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do i = 1, wf%n_o
         do c = 1, wf%n_v
            do k = 1, wf%n_o
               do b = 1, wf%n_v
!
                  X_bkci(b, k, c, i) = X_bick(b, i, c, k)/( eps_o(i) + eps_o(k)  &
                                                          - eps_v(b) - eps_v(c)  &
                                                          + omega)
!
               enddo
            enddo
         enddo
      enddo
!
      do k = 1, wf%n_o
         do b = 1, wf%n_v
!
            X_bkci(b, k, b, k) = X_bkci(b, k, b, k)/two
!
         enddo
      enddo
!
      call mem%dealloc(X_bick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vvov(L_abkc,      &
                        1, wf%n_v,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call dscal((wf%n_v**3)*(wf%n_o), -two, L_abkc, 1)
!
      call mem%alloc(g_ackb, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_vvov(g_ackb,      &
                        1, wf%n_v,  &
                        1, wf%n_v,  &
                        1, wf%n_o,  &
                        1, wf%n_v)
!
      call add_1432_to_1234(one, g_ackb, L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_ackb, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o*(wf%n_v**2),  &
                  one,                 &
                  L_abkc,              &
                  wf%n_v,              &
                  X_bkci,              &
                  wf%n_o*(wf%n_v**2),  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(L_abkc, wf%n_v, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(X_bkci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine effective_jacobian_cc2_f1_cc2
!
!
end submodule jacobian
