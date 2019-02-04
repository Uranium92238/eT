submodule (ccsd_class) F_ccsd
!
!!
!!    F-transformation submodule (ccsd)
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    Routines for the linear transform of
!!    vectors by the F matrix 
!!
!!    ρ = F * c,
!!
!!    where
!!   
!!    F_μ,ν = < Λ' | [[ exp(T) H_0 exp(T), τ_μ ], τ_ν ] | R >,
!!
!!    Where < Λ' | = < R | + sum_μ tbar_μ < μ |
!!  
! 
   implicit none
!
!
contains
!
!
   module subroutine F_transform_vector_ccsd(wf, c)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018     
!!
!!    Directs the transformation by the F matrix.
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj
!
      real(dp), dimension(:,:), allocatable :: rho_ai
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj
!
      integer(i15) :: a, i, b, j, ai, bj, aibj
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      rho_ai = zero
!
!$omp parallel do private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c_ai(a, i) = c(ai, 1)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%F_ccs_a1_0(c_ai, rho_ai)
!  
      call wf%F_ccs_a1_1(c_ai, rho_ai)
      call wf%F_ccs_b1_1(c_ai, rho_ai)
      call wf%F_ccs_c1_1(c_ai, rho_ai)
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do private(a, i, ai, bj, aibj)
      do b = 1, wf%n_v
         do j = 1, wf%n_o
!
            bj = wf%n_v*(j - 1) + b
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  aibj = max(ai, bj)*(max(ai, bj) - 3)/2 + ai + bj
!
                  c_aibj(a, i, b, j) = c((wf%n_v)*(wf%n_o) + aibj, 1)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            c_aibj(a, i, a, i) = two*c_aibj(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%F_ccsd_a1_1(c_aibj, rho_ai)
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      rho_aibj = zero
!
!     Terms to rho_aibj that should afterwards be permuted (ai,bj):
!
      call wf%F_ccsd_a2_1(c_ai, rho_aibj)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_transform_vector_ccsd
!
!
   module subroutine F_ccsd_a1_1_ccsd(wf, c_aibj, rho_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1_ai = 2 L_iajb c_bjck tbar_ck
!!                - (L_jbic tbar_ak + L_iajc tbar_bk + L_jbka tbar_ci) c_bjck 
!!
!!    Eqns (59) and (60)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: rho_ai
!
!     Local variables      
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb 
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj 
      real(dp), dimension(:,:,:,:), allocatable :: c_cjbk 
!
      real(dp), dimension(:,:), allocatable :: X_bj
      real(dp), dimension(:,:), allocatable :: X_ac
      real(dp), dimension(:,:), allocatable :: X_ki
      real(dp), dimension(:,:), allocatable :: X_cj
!
!     L_iajb = 2 g_iajb - g_jaib (ordered as L_aibj)
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iajb)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      L_aibj = zero
      call add_2143_to_1234(two, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(two, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Term 1: 2 L_iajb c_bjck tbar_ck
!
      call mem%alloc(X_bj, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  c_aibj,              & ! c_bj_ck
                  (wf%n_v)*(wf%n_o),   &
                  wf%t1bar,            &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_bj,                &
                  (wf%n_v)*(wf%n_o))
!
!     rho_ai += L_aibj X_bj ( = c_bjck L_ia )
!
      call dgemm('N', 'N', &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  two,                 &
                  L_aibj,              & ! L_ai_bj
                  (wf%n_v)*(wf%n_o),   &
                  X_bj,                &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_bj, wf%n_v, wf%n_o)
!
!     Term 2: - L_jbic tbar_ak c_bjck 
!
!     NOTE: We pretend that L_aibj(a, i, b, j) = L_iabj is
!           L_bjci(b, j, c, i) = L_jbic
!
!     X_ki = c_bjck L_bjci ( = c_bjck L_jbic )
!          
      call mem%alloc(X_ki, wf%n_o, wf%n_o) 
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  (wf%n_v**2)*wf%n_o,  & 
                  one,                 &
                  c_aibj,              & ! c_bjc_k
                  (wf%n_v**2)*wf%n_o,  & 
                  L_aibj,              & ! L_bjc_i
                  (wf%n_v**2)*wf%n_o,  & 
                  zero,                &
                  X_ki,                &
                  wf%n_o)
!
      call dgemm('N', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  wf%t1bar,   & ! tbar_ak
                  wf%n_v,     &
                  X_ki,       &
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)

!
      call mem%dealloc(X_ki, wf%n_o, wf%n_o) 
!
!     Term 3: -  L_iajc tbar_bk c_bjck 
!
!     NOTE: we will now pretend that L_aibj(a, i, b, j) = L_iabj is
!           L_aicj(a, i, c, j) = L_iajc
!
!     Reorder c_bjck to c_cjbk
!
      call mem%alloc(c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(c_aibj, c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_cj, wf%n_v, wf%n_o)
!
!     X_cj = c_cjbk t_bk ( = c_bjck t_bk)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  c_cjbk,              & ! c_cj_bk
                  (wf%n_v)*(wf%n_o),   &
                  wf%t1bar,            & ! t_bk
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_cj,                &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += L_aicj X_cj = (L_iajc X_cj)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  -one,                &
                  L_aibj,              & ! L_ai_ck
                  (wf%n_v)*(wf%n_o),   &
                  X_cj,                &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))

!
      call mem%alloc(X_cj, wf%n_v, wf%n_o)

!
!     Term 4_: - L_jbka tbar_ci c_bjck = - L_kajb tbar_ci c_bjck
!
!     NOTE: we will now pretend that L_aibj(a, i, b, j) = L_iabj is
!           L_akbj(a, k, b, j) = L_kajb
! 
!     X_ac = L_akbj c_ckbj ( = L_kajb c_ckbj)
!
      call mem%alloc(X_ac, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  L_aibj,              &  ! L_a_kbj
                  wf%n_v,              &
                  c_aibj,              &  ! c_c_kbj
                  wf%n_v,              &
                  zero,                &
                  X_ac,                &
                  wf%n_v)                    
!
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
!     rho_ai += X_ac tbar_ci
!
      call dgemm('N', 'N',       &
                  wf%n_v,        &
                  wf%n_o,        &
                  wf%n_v,        &
                  -one,          &
                  X_ac,          &
                  wf%n_v,        &
                  wf%t1bar,      & ! tbar_ci
                  wf%n_v,        &
                  one,           &
                  rho_ai,        &
                  wf%n_v)
!
      call mem%dealloc(X_ac, wf%n_v, wf%n_v)
!
   end subroutine F_ccsd_a1_1_ccsd
!
!
   module subroutine F_ccsd_a2_1_ccsd(wf, c_ai, rho_aibj)
!!
!!    F transformation A2,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A2,1_aibj = 2 L_iakc tbar_bj c_ck
!!                   - (L_jbic tbar_ak + L_jbka tbar_ci + L_kcib tbar_aj) c_ck
!!
!!    Eqns (62) and (63)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: L_iakc, g_iakc 
      real(dp), dimension(:,:,:,:), allocatable :: rho_jbia, rho_iajb
      real(dp), dimension(:,:,:,:), allocatable :: X_jbik
!
      real(dp), dimension(:,:), allocatable :: X_ia, X_ib, X_ik, c_kc 
!
      integer(i15) :: a, i, b, j
!
!     Construct L_iakc = 2 g_iakc - g_icka, which will be used in all terms  
!
      call mem%alloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iakc)
!
      L_iakc = two*g_iakc 
      call add_1432_to_1234(-one, g_iakc, L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 1, 2 L_iakc tbar_bj c_ck
!        
!     X_ia = 2 L_iakc c_ck 
!
      call mem%alloc(c_kc, wf%n_o, wf%n_v)
      call sort_12_to_21(c_ai, c_kc, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ia, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  two,               &
                  L_iakc,            &
                  (wf%n_o)*(wf%n_v), & 
                  c_kc,              &
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ia,              &
                  (wf%n_o)*(wf%n_v))
!
!     rho_aibj =+ X_ia tbar_bj
!
!$omp parallel do private(j,b,i,a)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v 
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + X_ia(i,a)*wf%t1bar(b,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(X_ia, wf%n_o, wf%n_v)
!
!     :: Term 2, - L_jbic tbar_ak c_ck 
!     
!     X_jbik = - L_jbic c_ck 
! 
      call mem%alloc(X_jbik, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  wf%n_v,               &
                  -one,                 &
                  L_iakc,               & ! L_jbi,c 
                  (wf%n_v)*(wf%n_o)**2, &
                  c_ai,                 & ! c_c,k
                  wf%n_v,               &
                  zero,                 &
                  X_jbik,               & ! X_jbi,k
                  (wf%n_v)*(wf%n_o)**2)
!
!     rho_jbia = X_jbik tbar_ak = X_jbik tbar^T_ka
!
      call mem%alloc(rho_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N','T',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  wf%n_o,               &
                  one,                  &
                  X_jbik,               & ! X_jbi,k
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%t1bar,             & ! tbar_a,k 
                  wf%n_v,               &
                  zero,                 &
                  rho_jbia,             & ! rho_jbi,a 
                  (wf%n_v)*(wf%n_o)**2)
!
!     rho_aibj =+ rho_jbia 
!
      call add_4321_to_1234(one, rho_jbia, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 3, - L_kajb tbar_ci c_ck
!
!     Construct X_ik = -tbar_ci c_ck = - tbar^T_ic c_ck 
!
      call mem%alloc(X_ik, wf%n_o, wf%n_o)
!
      call dgemm('T','N',   &
                  wf%n_o,   &
                  wf%n_o,   &
                  wf%n_v,   &
                  -one,     &
                  wf%t1bar, & ! tbar_c,i
                  wf%n_v,   &
                  c_ai,     & ! c_c,k 
                  wf%n_v,   &
                  zero,     &
                  X_ik,     & 
                  wf%n_o)
!
!     rho_iajb = X_ik L_kajb 
! 
      call mem%alloc(rho_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  wf%n_o,               &
                  one,                  &
                  X_ik,                 & ! X_i,k
                  wf%n_o,               &
                  L_iakc,               & ! L_k,ajb 
                  wf%n_o,               &
                  zero,                 &
                  rho_iajb,             & ! rho_i,ajb 
                  wf%n_o)
!
!     rho_aibj =+ rho_iajb 
!  
      call add_2143_to_1234(one, rho_iajb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 4, - L_kcib tbar_aj c_ck
!
!     X_ib = - c_ck L_kcib = - c_kc L_kcib 
!
      call mem%alloc(X_ib, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_kc,              &
                  (wf%n_o)*(wf%n_v), &
                  L_iakc,            & ! L_kc,ib 
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ib,              & ! X_1,ib 
                  1)
!
!     rho_aibj =+ X_ib tbar_aj
!
!$omp parallel do private(j,b,i,a)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v 
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + X_ib(i,b)*wf%t1bar(a,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(X_ib, wf%n_o, wf%n_v)
      call mem%dealloc(c_kc, wf%n_o, wf%n_v)
      call mem%dealloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
   end subroutine F_ccsd_a2_1_ccsd
!
!
   module subroutine F_ccsd_a1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation A1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1 = - (g_ibck tbar_ajck + g_jack tbar_bick) c_bj
!!
!!    First two terms of (65)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ckib
      real(dp), dimension(:,:,:,:), allocatable :: X_ckij, X_jcki, X_jick, X_ickj
!
      real(dp), dimension(:,:), allocatable :: rho_ia 
!
!     :: Term 1, - g_ckib tbar_ajck c_bj 
!
      call mem%alloc(g_ckib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov(g_ckib)
!
!     X_ckij = - g_ckib c_bj 
!
      call mem%alloc(X_ckij, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  wf%n_v,               &
                  -one,                 &
                  g_ckib,               & ! g_cki,b
                  (wf%n_v)*(wf%n_o)**2, &
                  c_ai,                 & ! c_b,j 
                  wf%n_v,               &
                  zero,                 &
                  X_ckij,               & ! X_cki,j
                  (wf%n_v)*(wf%n_o)**2)
!
!     rho_ai = tbar_ajck X_ckij = tbar_ajck X_jcki 
!
      call mem%alloc(X_jcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_4123(X_ckij, X_jcki, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(X_ckij, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,jck
                  wf%n_v,               &
                  X_jcki,               & ! X_jck,i 
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_jcki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2, - g_ckja tbar_bick c_bj
!
!     X_jick = - c_bj tbar_bick = - c^T_j,b tbar_b,ick 
!
      call mem%alloc(X_jick, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  c_ai,                 & ! c_b,j 
                  wf%n_v,               &
                  tbar_aibj,            & ! tbar_b,ick 
                  wf%n_v,               &
                  zero,                 &
                  X_jick,               & ! X_j,ick 
                  wf%n_o)
!
!     rho_ia = X_jick g_ckja = X_ickj g_ckja
!
      call mem%alloc(X_ickj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_2341(X_jick, X_ickj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_jick, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ia, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  X_ickj,               & ! X_i,ckj 
                  wf%n_o,               &
                  g_ckib,               & ! g_ckj,a
                  (wf%n_v)*(wf%n_o)**2, &
                  zero,                 &
                  rho_ia,               & ! rho_i,a
                  wf%n_o)
!
      call add_21_to_12(one, rho_ia, rho_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_ickj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_ckib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(rho_ia, wf%n_o, wf%n_v)
!
   end subroutine F_ccsd_a1_2_ccsd
!
!
   module subroutine F_ccsd_b1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation B1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_B1,1 = - (g_ikcb tbar_akcj + g_jkca tbar_bkci) c_bj
!!
!!    Term 3 & 4 of (65)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
!     Local variables 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ikcb
      real(dp), dimension(:,:,:,:), allocatable :: X_ikcj
      real(dp), dimension(:,:,:,:), allocatable :: X_jkci
!
!     Term 1: - g_ikcb tbar_akcj c_bj
!
      call mem%alloc(g_ikcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call wf%get_oovv(g_ikcb)
!
      call mem%alloc(X_ikcj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',                &
                  (wf%n_o**2)*(wf%n_v),   &
                  wf%n_o,                 &
                  wf%n_v,                 &
                  one,                    &
                  g_ikcb,                 & ! g_ikc_b
                  (wf%n_o**2)*(wf%n_v),   &
                  c_ai,                   & ! c_b_j 
                  wf%n_v,                 &
                  zero,                   &
                  X_ikcj,                 &
                  (wf%n_o**2)*(wf%n_v))
!
      call dgemm('N', 'T',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  -one,                   &
                  tbar_aibj,              & ! tbar_a_kcj
                  wf%n_v,                 &
                  X_ikcj,                 & ! X_i_kcj
                  wf%n_o,                 &
                  one,                    &
                  rho_ai,                 &
                  wf%n_v)
!
      call mem%dealloc(X_ikcj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2:  - g_jkca tbar_bkci c_bj
!
      call mem%alloc(X_jkci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_b_j
                  wf%n_v,              &
                  tbar_aibj,           & ! t_b_kci
                  wf%n_v,              &
                  zero,                &
                  X_jkci,              & ! X_j_kci
                  wf%n_o)
!
!     NOTE: we will now pretend that g_ikcb is g_jkca    
!
      call dgemm('T', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  -one,                &
                  g_ikcb,              & ! g_jkc_a
                  (wf%n_o**2)*wf%n_v,  &
                  X_jkci,              & ! X_jkc_i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(g_ikcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(X_jkci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_b1_2_ccsd
!
!
   module subroutine F_ccsd_c1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj)
!!
!!    F transformation C1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_C1,1 = - (g_ikjl tbar_akbl + g_cadb tbar_cidj) c_bj
!!
!!    Last two terms of (65)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: X_jlak
      real(dp), dimension(:,:,:,:), allocatable :: X_akjl
      real(dp), dimension(:,:,:,:), allocatable :: g_ikjl
      real(dp), dimension(:,:,:,:), allocatable :: X_bcdi
      real(dp), dimension(:,:,:,:), allocatable :: X_dbci
      real(dp), dimension(:,:,:,:), allocatable :: tbar_jcdi
      real(dp), dimension(:,:,:,:), allocatable :: g_dbca
!
      type(batching_index) :: batch_c, batch_d
!
      integer(i15) :: req0, req1_c, req1_d, req2, current_c_batch, current_d_batch
      integer(i15) :: c, d, i, j
!
!     Term 1:  - g_ikjl tbar_blak c_bj
!
      call mem%alloc(X_jlak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_b_j
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_b_lak
                  wf%n_v,              &
                  zero,                &
                  X_jlak,              & ! X_j_lak
                  wf%n_o)
!
!     Reorder X_jlak to X_akjl
!
      call mem%alloc(X_akjl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call sort_1234_to_3412(X_jlak, X_akjl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(X_jlak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_oooo(g_ikjl)
!
      call dgemm('N', 'T',       &
                  wf%n_v,        &
                  wf%n_o,        &
                  wf%n_o**3,     &
                  -one,          &
                  X_akjl,        &
                  wf%n_o,        &
                  g_ikjl,        &
                  wf%n_o,        &
                  one,           &
                  rho_ai,        &
                  wf%n_v)
!
      call mem%dealloc(X_akjl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Term 2: - g_cadb tbar_cidj c_bj
!
      req0 = 0
!
      req1_c =  wf%integrals%n_J*wf%n_v
      req1_d =  wf%integrals%n_J*wf%n_v
!
      req2 = (wf%n_o)*(wf%n_v) + wf%n_v**2
!
      call batch_c%init(wf%n_v)
      call batch_d%init(wf%n_v)
!
      call mem%batch_setup(batch_c, batch_d, req0, req1_c, req1_d, req2)
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
         do current_d_batch = 1, batch_d%num_batches
!
            call batch_d%determine_limits(current_d_batch)
!
            call mem%alloc(tbar_jcdi, wf%n_o, batch_c%length, batch_d%length, wf%n_o)
!
!           Reorder tbar_cidj to tbar_jcdi
!
!$omp parallel do private(i, d, c, j)
            do i = 1, wf%n_o
               do d = 1, batch_d%length
                  do c = 1, batch_c%length
                     do j = 1, wf%n_o
!
                        tbar_jcdi(j, c, d, i) = tbar_aibj(c + batch_c%first - 1, i, d + batch_d%first - 1, j)
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
!           X_bcdi = tbar_jcdi c_bj ( = tbar_cicj c_bj )
!
            call mem%alloc(X_bcdi, wf%n_v, batch_c%length, batch_d%length, wf%n_o)
!
            call dgemm('N', 'N',                                     &
                        wf%n_v,                                      &
                        (batch_d%length)*(batch_c%length)*(wf%n_o),  &
                        wf%n_o,                                      &
                        one,                                         &
                        c_ai,                                        & ! c_bj 
                        wf%n_v,                                      &
                        tbar_jcdi,                                   & ! tbar_j_cdi
                        wf%n_o,                                      &
                        zero,                                        &
                        X_bcdi,                                      & ! X_b_cdi
                        wf%n_v)
!
            call mem%dealloc(tbar_jcdi, wf%n_o, batch_c%length, batch_d%length, wf%n_o)
!
!           Reorder X_bcdi to X_dbci
!
            call mem%alloc(X_dbci, batch_d%length, wf%n_v, batch_c%length, wf%n_o)
            call sort_1234_to_3124(X_bcdi, X_dbci, wf%n_v, batch_c%length, batch_d%length, wf%n_o)
            call mem%dealloc(X_bcdi, wf%n_v, batch_c%length, batch_d%length, wf%n_o)
!
            call mem%alloc(g_dbca, batch_d%length, wf%n_v, batch_c%length, wf%n_v)
            call wf%get_vvvv(g_dbca,                        &
                              batch_d%first, batch_d%last,  &
                              1, wf%n_v,                    &
                              batch_c%first, batch_c%last,  &
                              1, wf%n_v)
!
            call dgemm('T', 'N',    &
                        wf%n_v,     &
                        wf%n_o,     &
                        wf%n_v**3,  &
                        -one,       &
                        g_dbca,     & ! g_dbc_a
                        wf%n_v**3,  &
                        X_dbci,     & ! X_dbc_i
                        wf%n_v**3,  &
                        one,        &
                        rho_ai,     &
                        wf%n_v)
!
            call mem%dealloc(g_dbca, batch_d%length, wf%n_v, batch_c%length, wf%n_v)
            call mem%dealloc(X_dbci, batch_d%length, wf%n_v, batch_c%length, wf%n_o)
!
         enddo ! batch_d
      enddo ! batch_c

   end subroutine F_ccsd_c1_2_ccsd
!
!
   module subroutine F_ccsd_d1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation D1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_D1,1 = (g_iljc tbar_albk + g_jlic tbar_blak + g_klja tbar_clbi) c_bjck
!!
!!    First three terms of (66)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: c_bkjc
      real(dp), dimension(:,:,:,:), allocatable :: c_jcbk
      real(dp), dimension(:,:,:,:), allocatable :: c_kjcb
      real(dp), dimension(:,:,:,:), allocatable :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable :: g_jcli
      real(dp), dimension(:,:,:,:), allocatable :: tbar_cbli
      real(dp), dimension(:,:,:,:), allocatable :: X_bkli
      real(dp), dimension(:,:,:,:), allocatable :: X_ilbk
      real(dp), dimension(:,:,:,:), allocatable :: X_kbli
      real(dp), dimension(:,:,:,:), allocatable :: X_kjli
      real(dp), dimension(:,:,:,:), allocatable :: X_klji
!
!     Term 1: g_iljc tbar_albk c_bjck
!
!     Reorder c_ckbj to c_jcbk 
!
      call mem%alloc(c_jcbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4132(c_aibj, c_jcbk , wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_iljc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_ooov(g_iljc)
!
      call mem%alloc(X_ilbk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  g_iljc,              & ! g_il_jc
                  wf%n_o**2,           &
                  c_jcbk,              & ! c_jc_bk
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_ilbk,              & ! X_il_bk
                  wf%n_o**2)
!
      call dgemm('N', 'T',                &
                  (wf%n_v),               &
                  (wf%n_o),               &
                  (wf%n_o**2)*(wf%n_v),   &
                  one,                    &
                  tbar_aibj,              & ! tbar_a_lbk
                  wf%n_v,                 &
                  X_ilbk,                 & ! X_i_lbk
                  wf%n_o,                 &
                  one,                    &
                  rho_ai,                 &
                  wf%n_v)
!
      call mem%dealloc(X_ilbk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2: g_jlic tbar_blak c_bjck
!
!     NOTE: We now pretend that g_iljc is g_jlic
!
!     Reorder g_jlic to g_jcli
!
      call mem%alloc(g_jcli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1423(g_iljc, g_jcli, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%alloc(X_bkli, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  wf%n_o**2,           &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  c_jcbk,              & ! c_jc_bk
                  (wf%n_o)*(wf%n_v),   &
                  g_jcli,              & ! g_jc_li
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_bkli,              & ! X_bk_li 
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(g_jcli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Reorder X_bkli to X_kbli
!
      call mem%alloc(X_kbli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_2134(X_bkli, X_kbli, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_bkli, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  tbar_aibj,           & ! t_a_kbl
                  wf%n_v,              &
                  X_kbli,              & ! X_kbl_i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)

!
      call mem%alloc(X_kbli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(c_bkjc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     Term 3:  g_klja tbar_clbi c_bjck
!
!     Reorder tbar_clbi to tbar_cbli and c_ckbj to c_kjcb
!
      call mem%alloc(tbar_cbli, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(tbar_aibj, tbar_cbli, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_kjcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2413(c_aibj, c_kjcb, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_kjli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',       &
                  wf%n_o**2,     &
                  wf%n_o**2,     &
                  wf%n_v**2,     &
                  one,           &
                  c_kjcb,        & ! c_kj_cb
                  wf%n_o**2,     &
                  tbar_cbli,     & ! tbar_ck_li
                  wf%n_v**2,     &
                  zero,          &
                  X_kjli,        &
                  wf%n_o**2)
!
      call mem%dealloc(tbar_cbli, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(c_kjcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     Reorder X_kjli to X_klji
!
      call mem%alloc(X_klji, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_1234_to_1324(X_kjli, X_klji, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_kjli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     NOTE: we now pretend that g_iljc is g_klja
!
      call dgemm('T', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o**3,  &
                  one,        &
                  g_iljc,     & ! g_klj_a
                  wf%n_o**3,  &
                  X_klji,     & ! X_klj_i
                  wf%n_o**3,  &
                  one,        &
                  rho_ai,     &
                  wf%n_v) 
!
      call mem%dealloc(g_iljc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_klji, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine F_ccsd_d1_2_ccsd
!
!
   module subroutine F_ccsd_e1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation E1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_E1,1 = -(g_ibdc tbar_ajdk + g_jadc tbar_bidk + g_kbda tbar_cjdi) c_bjck
!!
!!    Last three terms of (66)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_jkid, X_jdki, X_idcj, X_kbdi, X_kbdi_red, X_idcj_red
      real(dp), dimension(:,:,:,:), allocatable :: c_jkbc, c_bkcj, c_kbcj
      real(dp), dimension(:,:,:,:), allocatable :: tbar_idbk
      real(dp), dimension(:,:,:,:), allocatable :: g_ibdc, g_bcid, g_dcja, g_kbda
!
      integer(i15) :: req0, req1, current_d_batch
      type(batching_index) :: batch_d
!
!     :: Term 1, - g_ibdc tbar_ajdk c_bjck
!
!     X_jkid = - sum_bc g_ibdc c_bjck
!            = - sum_bc c_jk,bc g_bc,id  
!
!     Construct intermediate in batches over d 
!
      call mem%alloc(c_jkbc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2413(c_aibj, c_jkbc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_jkid, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      X_jkid = zero 
!
      req0 = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%integrals%n_J)*(wf%n_v) + 2*(wf%n_o)*(wf%n_v)**2
!
      call batch_d%init(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1)
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_ibdc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%get_ovvv(g_ibdc,                       &
                           1, wf%n_o,                   &
                           1, wf%n_v,                   &
                           batch_d%first, batch_d%last, &
                           1, wf%n_v) 
!
         call mem%alloc(g_bcid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
         call sort_1234_to_2413(g_ibdc, g_bcid, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call mem%dealloc(g_ibdc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call dgemm('N','N',                      &
                     (wf%n_o)**2,                 &
                     (wf%n_o)*(batch_d%length),   &
                     (wf%n_v)**2,                 &
                     -one,                        &
                     c_jkbc,                      & ! c_jk,bc 
                     (wf%n_o)**2,                 &
                     g_bcid,                      & ! g_bc,i#d
                     (wf%n_v)**2,                 &
                     one,                         &
                     X_jkid(1,1,1,batch_d%first), & ! X_jk,id
                     (wf%n_o)**2)
!
         call mem%dealloc(g_bcid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      enddo
!
      call mem%dealloc(c_jkbc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     rho_ai =+ tbar_ajdk X_jkid = tbar_ajdk X_jdki
!
      call mem%alloc(X_jdki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1423(X_jkid, X_jdki, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_jkid, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,jdk
                  wf%n_v,               &
                  X_jdki,               & ! X_jdk,i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_jdki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2, - g_jadc tbar_bidk c_bjck
!
!     X_idjc = - tbar_bidk c_bjck = - tbar_id,bk c_bk,cj
! 
      call mem%alloc(tbar_idbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2314(tbar_aibj, tbar_idbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_bkcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(c_aibj, c_bkcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_idcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  tbar_idbk,         & ! tbar_id,bk
                  (wf%n_o)*(wf%n_v), &
                  c_bkcj,            & ! c_bk,cj
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_idcj,            & ! X_id,cj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(tbar_idbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(c_bkcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai =+ X_idcj g_jadc = X_idcj g_dcja 
!
!     We construct constributions in batches over d
!
      req0 = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%integrals%n_J)*(wf%n_v) + (wf%n_o)*(wf%n_v)**2 + (wf%n_v)*(wf%n_o)**2
!
      call batch_d%init(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1)
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_dcja, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_dcja,                       &
                           batch_d%first, batch_d%last, &
                           1, wf%n_v,                   &
                           1, wf%n_o,                   &
                           1, wf%n_v)
!
         call mem%alloc(X_idcj_red, wf%n_o, batch_d%length, wf%n_v, wf%n_o)
!
         X_idcj_red(:,:,:,:) = X_idcj(:,batch_d%first:batch_d%last,:,:)
!
         call dgemm('T','T',                             &
                     wf%n_v,                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     one,                                &
                     g_dcja,                             & ! g_dcj,a 
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     X_idcj_red,                         & ! g_i,dcj
                     wf%n_o,                             &
                     one,                                &
                     rho_ai,                             & ! rho_a,i
                     wf%n_v)
!
         call mem%dealloc(g_dcja, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
         call mem%dealloc(X_idcj_red, wf%n_o, batch_d%length, wf%n_v, wf%n_o)
!
      enddo
!
!     :: Term 3, - g_kbda tbar_cjdi c_bjck
!
!     X_kbdi = - tbar_cjdi c_bjck = - c_kb,cj tbar_cj,di
!
      call mem%alloc(c_kbcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4132(c_aibj, c_kbcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_kbdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  c_kbcj,            & ! c_kb,cj 
                  (wf%n_o)*(wf%n_v), &
                  tbar_aibj,         & ! tbar_cj,di
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_kbdi,            & ! X_kb,di
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_kbcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     rho_ai =+ X_kbdi g_kbda 
!
!     Construct term in batches over d 
!
      req0 = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%integrals%n_J)*(wf%n_v) + (wf%n_o)*(wf%n_v)**2 + (wf%n_v)*(wf%n_o)**2
!
      call batch_d%init(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1)
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_kbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%get_ovvv(g_kbda,                       &
                           1, wf%n_o,                   &
                           1, wf%n_v,                   &
                           batch_d%first, batch_d%last, &
                           1, wf%n_v)
!
         call mem%alloc(X_kbdi_red, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
         X_kbdi_red(:,:,:,:) = X_kbdi(:,:,batch_d%first:batch_d%last,:)
!
         call dgemm('T','N',                             &
                     wf%n_v,                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     one,                                &
                     g_kbda,                             & ! g_kbd,a
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     X_kbdi_red,                         & ! X_kbd,i
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     one,                                &
                     rho_ai,                             & ! rho_a,i
                     wf%n_v)
!
         call mem%dealloc(X_kbdi_red, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
         call mem%dealloc(g_kbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
      enddo 
!
      call mem%dealloc(X_kbdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_e1_2_ccsd
!
!
   module subroutine F_ccsd_f1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation F1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_F1,2 = -(F_ib tbar_ckaj + F_ja tbar_ckbi) c_bjck 
!!
!!    The first two terms of (67)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ijck 
!
      real(dp), dimension(:,:), allocatable :: X_ij
!
!     :: Term 1, - F_ib tbar_ajck c_bjck 
!
!     X_ijck = - F_ib c_bjck 
!
      call mem%alloc(X_ijck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  -one,                 &
                  wf%fock_ia,           & ! F_i,b 
                  wf%n_o,               &
                  c_aibj,               & ! c_b,jck
                  wf%n_v,               &
                  zero,                 &
                  X_ijck,               & ! X_i,jck
                  wf%n_o)
!
!     rho_ai =+ tbar_a,jck X_i,jck = tbar_a,jck X^T_jck,i 
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,jck
                  wf%n_v,               &
                  X_ijck,               & ! X_i,jck
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_ijck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2, - F_ja tbar_ckbi c_ckbj 
!
!     X_ij = - tbar_ckb,i c_ckb,j 
!
      call mem%alloc(X_ij, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  -one,                 &
                  tbar_aibj,            & ! tbar_ckb,i
                  (wf%n_o)*(wf%n_v)**2, &
                  c_aibj,               & ! c_ckb,j)
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_ij,                 & ! X_i,j
                  wf%n_o)
!
!     rho_ai =+ F_j,a X_i,j
!
      call dgemm('T','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  one,        & 
                  wf%fock_ia, & ! F_j,a 
                  wf%n_o,     &
                  X_ij,       & ! X_i,j
                  wf%n_o,     &
                  one,        &
                  rho_ai,     & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_ij, wf%n_o, wf%n_o)
!
   end subroutine F_ccsd_f1_2_ccsd
!
!
   module subroutine F_ccsd_g1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation G1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_G1,2 = -(L_iljb tbar_ckal + L_jlia tbar_ckbl) c_bjck 
!!
!!    The third and fourth terms of (67)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iljb, L_ilbj 
      real(dp), dimension(:,:,:,:), allocatable :: X_ilck
!
      real(dp), dimension(:,:), allocatable :: X_jl 
!
!     :: Term 1, - L_iljb tbar_alck c_bjck 
!
!     Construct L_iljb ordered as L_ilbj = 2 g_iljb - g_jlib
!
      call mem%alloc(L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_iljb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%get_ooov(g_iljb)
!
      L_ilbj = zero 
!
      call add_1243_to_1234(two, g_iljb, L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4213_to_1234(-one, g_iljb, L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iljb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     X_ilck = - L_iljb c_bjck = - L_il,bj c_bj,ck 
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)**2,       &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  L_ilbj,            & ! L_il,bj 
                  (wf%n_o)**2,       &
                  c_aibj,            & ! c_bj,ck 
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_ilck,            &
                  (wf%n_o)**2)
!
!     rho_ai =+ tbar_a,lck X_i,lck 
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,lck
                  wf%n_v,               &
                  X_ilck,               & ! X_i,lck 
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               & ! rho_a,i 
                  wf%n_v)
!
      call mem%dealloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2, - L_jlia tbar_ckbl c_ckbj 
!
!     X_jl = - c_ckb,j tbar_ckb,l 
!
      call mem%alloc(X_jl, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &   
                  wf%n_o,               &
                  wf%n_o,               &
                  (wf%n_o)*(wf%n_v)**2, &
                  -one,                 &
                  c_aibj,               & ! c_ckb,j
                  (wf%n_o)*(wf%n_v)**2, &
                  tbar_aibj,            & ! tbar_ckb,l)
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_jl,                 & ! X_j,l
                  wf%n_o)
!
!     rho_ai = X_jl L_jlia = L_jl,ai X_jl 
!
!     Note: recall that we have L_ljia ordered as L_ljai 
!
      call dgemm('T','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)**2,       &
                  one,               &
                  L_ilbj,            & ! L_jl,ai
                  (wf%n_o)**2,       &
                  X_jl,              & 
                  (wf%n_o)**2,       &    
                  one,               &
                  rho_ai,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_jl, wf%n_o, wf%n_o)
!
   end subroutine F_ccsd_g1_2_ccsd
!
!
   module subroutine F_ccsd_h1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation H1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_H1,2 = (L_dajb tbar_ckdi + L_dbia tbar_ckdj) c_bjck 
!!
!!    The last two terms of (67)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_bjdi, X_jbdi, red_X_jbdi
!
      real(dp), dimension(:,:), allocatable :: X_db, red_X_db, rho_ia
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbda, L_jbda
!
      integer(i15)         :: req0, req1, current_d_batch
      type(batching_index) :: batch_d
!
!     Note: we construct both terms at once 
!
!     Make intermediates for both: 
!
!        X_bjdi = c_bjck tbar_ckdi = c_bj,ck tbar_ck,di 
!        X_db   = tbar_djck c_bjck = tbar_d,jck c_b,jck 
!
      call mem%alloc(X_bjdi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  c_aibj,            & ! c_bj,ck
                  (wf%n_o)*(wf%n_v), &
                  tbar_aibj,         & ! tbar_ck,di
                  (wf%n_o)*(wf%n_v), &
                  zero,              &
                  X_bjdi,            & ! X_bj,di
                  (wf%n_o)*(wf%n_v))
!
      call mem%alloc(X_db, wf%n_v, wf%n_v)
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_d,jck
                  wf%n_v,               &
                  c_aibj,               & ! c_b,jck
                  wf%n_v,               &
                  zero,                 &
                  X_db,                 & ! X_d,b
                  wf%n_v)
!
!     rho_ai =+ L_jbda X_bjdi = L_jbd,a X_jbd,i 
!     rho_ia =+ L_ia,db X_db (accumulate in temporary rho_ia)
!
!     Note: we batch over d 
!
      call mem%alloc(X_jbdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2134(X_bjdi, X_jbdi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(X_bjdi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ia, wf%n_o, wf%n_v)
      rho_ia = zero 
!
      req0 = (wf%integrals%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%integrals%n_J)*(wf%n_v) + 2*(wf%n_o)*(wf%n_v)**2 &
                  + (wf%n_v)*(wf%n_o)**2 + wf%n_v
!
      call batch_d%init(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1)
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%get_ovvv(g_jbda,                       &
                           1, wf%n_o,                   &
                           1, wf%n_v,                   &
                           batch_d%first, batch_d%last, &
                           1, wf%n_v)
!
         call mem%alloc(L_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         L_jbda = two*g_jbda
         call add_1432_to_1234(-one, g_jbda, L_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call mem%dealloc(g_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
!        rho_ai =+ L_jbd,a X_jbd,i 
!
         call mem%alloc(red_X_jbdi, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
         red_X_jbdi(:,:,:,:) = X_jbdi(:,:,batch_d%first:batch_d%last,:)
!
         call dgemm('T','N',                             &
                     wf%n_v,                             &
                     wf%n_o,                             &
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     one,                                &
                     L_jbda,                             & ! L_jbd,a 
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     red_X_jbdi,                         & ! X_jbd,i
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     rho_ai,                             & ! rho_a,i
                     wf%n_v)
!
         call mem%dealloc(red_X_jbdi, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
!        rho_ia =+ L_iadb X_db      
!
         call mem%alloc(red_X_db, batch_d%length, wf%n_v)
!
         red_X_db(:,:) = X_db(batch_d%first:batch_d%last,:)
!
         call dgemm('N','N',                    &
                     (wf%n_o)*(wf%n_v),         &
                     1,                         &
                     (wf%n_v)*(batch_d%length), &
                     one,                       & 
                     L_jbda,                    & ! L_ia,db
                     (wf%n_o)*(wf%n_v),         &
                     red_X_db,                  & ! X_db 
                     (wf%n_v)*(batch_d%length), &
                     one,                       &
                     rho_ia,                    & ! rho_ia 
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(L_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
         call mem%dealloc(red_X_db, batch_d%length, wf%n_v)
!
      enddo 
!
      call add_21_to_12(one, rho_ia, rho_ai, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ia, wf%n_o, wf%n_v)
!
   end subroutine F_ccsd_h1_2_ccsd
!
!
   module subroutine F_ccsd_i1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj, t_aibj)
!!
!!    F transformation I1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_I1,2 = - L_ldic tbar_bjak (t_bjck c_dl + t_bjdl c_ck) 
!!               - L_ldka tbar_bjci (t_bjck c_dl + t_bjdl c_ck)
!!               - L_ialc tbar_bjdk (t_bjck c_dl + t_bjdl c_ck)
!!
!!    Equation (23) in Sarai's document
!!    Equation (69) in Eirik's document
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable :: X_cd, X_ci, X_kl
      real(dp), dimension(:,:), allocatable :: Y_ac, Y_cl, Y_ki
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldic
      real(dp), dimension(:,:,:,:), allocatable :: L_dlci
      real(dp), dimension(:,:,:,:), allocatable :: X_kibj, X_kidl
      real(dp), dimension(:,:,:,:), allocatable :: Y_kbji, Y_kdli, Y_kibj, Y_kidl
!
!     Term 1: - L_ldic tbar_bjak t_bjck c_dl  
!
!     Construct L_ldic = 2 g_ldic - g_lcid (ordered as L_dlci)
!
      call mem%alloc(g_ldic, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_ldic)
!
      call mem%alloc(L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_dlci = zero
!
      call add_2143_to_1234(two, g_ldic, L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_ldic, L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ldic, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!  
!     X_ci = c_dl L_dlci (c_dl L_ldic)
!
      call mem%alloc(X_ci, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',          &
                  (wf%n_o)*(wf%n_v),&
                  1,                &
                  (wf%n_o)*(wf%n_v),&
                  one,              &
                  L_dlci,           & ! L_dl_ci
                  (wf%n_o)*(wf%n_v),&
                  c_ai,             & ! c_dl
                  (wf%n_o)*(wf%n_v),&
                  zero,             &
                  X_ci,             &
                  (wf%n_o)*(wf%n_v))
!
!     Y_ac = tbar_akbj t_ckbj
!
      call mem%alloc(Y_ac, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             & 
                  wf%n_v,              &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  tbar_aibj,           & ! tbar_a_kbj
                  wf%n_v,              &
                  t_aibj,              & ! t_c_kbj
                  wf%n_v,              &
                  zero,                &
                  Y_ac,                &
                  wf%n_v)
!
!     rho_ai -= Y_ac X_ci
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  Y_ac,    & ! Y_a_c
                  wf%n_v,  &
                  X_ci,    & ! X_c_i
                  wf%n_v,  &
                  one,     &
                  rho_ai,  & ! rho_a_i
                  wf%n_v)
!
      call mem%dealloc(Y_ac, wf%n_v, wf%n_v)
!
!     Term 2: - L_ldic tbar_bjak t_bjdl c_ck
!
!     Using L_dlci = L_cidl
!     X_kidl = L_cidl c_ck (L_ldic c_ck)
!
      call mem%alloc(X_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_c_k 
                  wf%n_v,              &
                  L_dlci,              & ! L_c_idl
                  wf%n_v,              &
                  zero,                &
                  X_kidl,              &
                  wf%n_o)
!
!     Y_kibj =  X_kidl t_dlbj
!
      call mem%alloc(Y_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  X_kidl,              & ! X_ki_dl
                  wf%n_o**2,           &
                  t_aibj,              & ! t_dl_bj
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  Y_kibj,              &
                  wf%n_o**2)
!
      call mem%dealloc(X_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder Y_kibj to Y_kbji
!
      call mem%alloc(Y_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1342(Y_kibj, Y_kbji, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(Y_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai -= tbar_akbj Y_kbji 
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  -one,                &
                  tbar_aibj,           & ! t_a_kbj 
                  wf%n_v,              &
                  Y_kbji,              & ! Y_kbj_i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(Y_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Term 3: - L_ldka tbar_bjci t_bjck c_dl 
!
!     Y_ki = t_bjck tbar_bjci 
!
      call mem%alloc(Y_ki, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_v**2)*(wf%n_o),   &
                  one,                    &
                  t_aibj,                 & ! t_bjc_k
                  (wf%n_v**2)*(wf%n_o),   &
                  tbar_aibj,              & ! tbar_bjc_i
                  (wf%n_v**2)*(wf%n_o),   &
                  zero,                   &
                  Y_ki,                   & ! Y_k_i
                  wf%n_o)

!     Pretend that X_ci = c_dl L_ldic is X_ak = c_dl L_ldka
!
!     rho_ai -= X_ak Y_ki
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_o,  &
                  -one,    &
                  X_ci,    &
                  wf%n_v,  &
                  Y_ki,    &
                  wf%n_o,  &
                  one,     &
                  rho_ai,  &
                  wf%n_v)
!
      call mem%dealloc(X_ci, wf%n_v, wf%n_o)
!
!     Term 4:  - L_ldka tbar_bjci t_bjdl c_ck
!
!     X_kibj = c_ck tbar_cibj
!
      call mem%alloc(X_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  wf%n_v,                 &
                  one,                    &
                  c_ai,                   & ! c_c_k
                  wf%n_v,                 &
                  tbar_aibj,              & ! tbar_c_ibj
                  wf%n_v,                 &
                  zero,                   &
                  X_kibj,                 &
                  wf%n_o)
!
!     Y_kidl = X_kibj t_bjdl
!
      call mem%alloc(Y_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  X_kibj,              & ! X_ki_bj
                  wf%n_o**2,           &
                  t_aibj,              & ! t_bj_dl
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  Y_kidl,              & ! Y_ki_dl
                  wf%n_o**2)
!
      call mem%dealloc(X_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Pretend that L_dlci is L_dlak = L_akdl (L_ldka = 2 g_ldka - g_lakd)
!  
!     Reorder Y_kidl to Y_kdli
!
      call mem%alloc(Y_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1342(Y_kidl, Y_kdli, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(Y_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  -one,                &
                  L_dlci,              & ! L_a_kdl
                  wf%n_v,              &
                  Y_kdli,              & ! Y_kdl_i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%alloc(Y_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Term 5: - L_ialc t_ckbj tbar_dkbj c_dl
!
!     X_cd = t_ckbj tbar_dkbj
!
      call mem%alloc(X_cd, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o**2,  &
                  one,        &
                  t_aibj,     & ! t_c_kbj
                  wf%n_v,     &
                  tbar_aibj,  & ! tbar_d_kbj
                  wf%n_v,     &
                  zero,       &
                  X_cd,       &
                  wf%n_v)
!
!     Y_cl =  X_cd c_dl
!
      call mem%alloc(Y_cl, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N', &
                  wf%n_v,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  one,     &
                  X_cd,    & ! X_c_d
                  wf%n_v,  &
                  c_ai,    & ! c_d_l
                  wf%n_v,  &
                  zero,    &
                  Y_cl,    &
                  wf%n_v)
!
      call mem%dealloc(X_cd, wf%n_v, wf%n_v)
!
!     rho_ai -= L_ialc Y_cl
!
!     NOTE: Pretend that L_dlci (= L_ldic) is L_aicl 
!
      call dgemm('N', 'N',       &
                  wf%n_v*wf%n_o, &
                  1,             &
                  wf%n_v*wf%n_o, &
                  -one,          &
                  L_dlci,        &
                  wf%n_v*wf%n_o, &
                  Y_cl,          &
                  wf%n_v*wf%n_o, &
                  one,           &
                  rho_ai,        &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(Y_cl, wf%n_v, wf%n_o)
!
!     Term 6:  - L_ialc tbar_bjdk t_bjdl c_ck 
!
!     X_kl = tbar_bjdk t_bjdl
!
      call mem%alloc(X_kl, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  (wf%n_v**2)*wf%n_o,  &
                  one,                 &
                  tbar_aibj,           & ! tbar_bjd_k
                  (wf%n_v**2)*wf%n_o,  &
                  t_aibj,              & ! t_bjd_l
                  (wf%n_v**2)*wf%n_o,  &
                  zero,                &
                  X_kl,                &
                  wf%n_o)
!
!     Y_cl = c_ck X_kl
!
      call mem%alloc(Y_cl, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  one,        &
                  c_ai,       & ! c_c_k
                  wf%n_v,     &
                  X_kl,       & ! X_k_l
                  wf%n_o,     &
                  zero,       &
                  Y_cl,       & ! Y_c_l
                  wf%n_v)
! 
      call mem%dealloc(X_kl, wf%n_o, wf%n_o)
!
!     rho_ai -= L_ialc Y_cl 
!
!     NOTE: Still pretending L_dlci (= L_ldic) is L_aicl 
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  -one,                &
                  L_dlci,              & ! L_aicl
                  (wf%n_o)*(wf%n_v),   &
                  Y_cl,                &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(Y_cl, wf%n_v, wf%n_o)
      call mem%dealloc(L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_i1_2_ccsd
!
!
   module subroutine F_ccsd_j1_2_ccsd(wf, c_ai, rho_ai, tbar_aibj, t_aibj)
!!
!!    F transformation J1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_J1,2 = (g_kbid tbar_cjal + g_jcid tbar_bkal + g_kajd tbar_cibl 
!!                g_jakd tbar_bicl + g_ibkd tbar_ajcl + g_lakb tbar_dicj) t_ckdl c_bj
!!
!!    Equation (68)
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
   end subroutine F_ccsd_j1_2_ccsd
!
!
end submodule F_ccsd
