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
end submodule F_ccsd
