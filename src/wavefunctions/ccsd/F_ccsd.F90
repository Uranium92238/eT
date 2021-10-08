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
   module subroutine F_x_mu_transformation_ccsd(wf, c, rho, x)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018     
!!
!!    Directs the transformation by the F matrix, defined as
!!
!!       F(X)_mu,nu = < X | [[H-bar,tau_mu],tau_nu] | HF >
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in) :: c
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(out) :: rho
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in) :: x
!
      real(dp), dimension(:,:), allocatable :: c_ai, x_ai, rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj, c_abij
      real(dp), dimension(:,:,:,:), allocatable :: x_aibj, x_abij
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj, rho_abij
!
      type(timings), allocatable :: timer 
!
      timer = timings('F transformation CCSD', 'm')
      call timer%turn_on()
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(x_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, x, 1, x_ai, 1)
!
      call zero_array(rho_ai, wf%n_v*wf%n_o)
!
      call dcopy(wf%n_t1, c, 1, c_ai, 1)
!
      call wf%F_ccs_a1_1(c_ai, rho_ai, x_ai)
      call wf%F_ccs_b1_1(c_ai, rho_ai, x_ai)
      call wf%F_ccs_c1_1(c_ai, rho_ai, x_ai)
!
      call mem%alloc(x_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(x(wf%n_t1 + 1:), x_aibj, wf%n_t1)
!
      call wf%F_ccsd_a1_2(c_ai, rho_ai, x_aibj)
      call wf%F_ccsd_b1_2(c_ai, rho_ai, x_aibj)
      call wf%F_ccsd_c1_2(c_ai, rho_ai, x_aibj)
!
      call mem%alloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(c(wf%n_t1 + 1:), c_aibj, wf%n_t1)
!
      call scale_diagonal(two, c_aibj, wf%n_t1)
!
      call wf%F_ccsd_a1_1(c_aibj, rho_ai, x_ai)
!
      call wf%F_ccsd_d1_2(c_aibj, rho_ai, x_aibj)
      call wf%F_ccsd_e1_2(c_aibj, rho_ai, x_aibj)
      call wf%F_ccsd_f1_2(c_aibj, rho_ai, x_aibj)
      call wf%F_ccsd_g1_2(c_aibj, rho_ai, x_aibj)
      call wf%F_ccsd_h1_2(c_aibj, rho_ai, x_aibj)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_t1)
!  
      call wf%F_ccsd_i1_2(c_ai, rho_ai, x_aibj, t_aibj)
      call wf%F_ccsd_j1_2(c_ai, rho_ai, x_aibj, t_aibj)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, rho_ai, 1, rho, 1)
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(rho_aibj, (wf%n_v**2)*(wf%n_o**2))
!
      call wf%F_ccsd_a2_1(c_ai, rho_aibj, x_ai)
!
      call mem%dealloc(x_ai, wf%n_v, wf%n_o)
!
      call wf%F_ccsd_a2_2(c_ai, rho_aibj, x_aibj)
      call wf%F_ccsd_b2_2(c_ai, rho_aibj, x_aibj)
      call wf%F_ccsd_c2_2(c_ai, rho_aibj, x_aibj)
      call wf%F_ccsd_d2_2(c_ai, rho_aibj, x_aibj)
      call wf%F_ccsd_e2_2(c_ai, rho_aibj, x_aibj)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
      call wf%F_ccsd_f2_2(c_aibj, rho_aibj, x_aibj)
      call wf%F_ccsd_h2_2(c_aibj, rho_aibj, x_aibj)
!
      call symmetric_sum(rho_aibj, wf%n_t1)
!
      call mem%alloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(c_aibj, c_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(x_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(x_aibj, x_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(x_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(rho_aibj, rho_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%F_ccsd_g2_2(c_abij, rho_abij, x_abij)
!
      call mem%dealloc(x_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(c_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call zero_array(rho(wf%n_t1+1:), wf%n_t2)
!
      call packin_and_add_from_1324_order_real(rho(wf%n_t1+1:), rho_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine F_x_mu_transformation_ccsd
!
!
   module subroutine F_ccsd_a1_1_ccsd(wf, c_aibj, rho_ai, tbar_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1_ai = 2 L_iajb * c_bjck * tbar_ck
!!                - (L_jbic * tbar_ak + L_iajc * tbar_bk + L_jbka * tbar_ci) * c_bjck 
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: tbar_ai 
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
      call wf%eri%get_eri_t1('ovov',g_iajb)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(L_aibj, (wf%n_v**2)*(wf%n_o**2))
      call add_2143_to_1234(two, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_iajb, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Term 1: 2 * L_iajb * c_bjck * tbar_ck
!
!     X_bj = c_bjck * tbar_ck
!
      call mem%alloc(X_bj, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  c_aibj,              & ! c_bj,ck
                  (wf%n_v)*(wf%n_o),   &
                  tbar_ai,             &
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_bj,                &
                  (wf%n_v)*(wf%n_o))
!
!     rho_ai += 2 * L_aibj * X_bj 
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  two,                 &
                  L_aibj,              & ! L_ai,bj
                  (wf%n_v)*(wf%n_o),   &
                  X_bj,                &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_bj, wf%n_v, wf%n_o)
!
!     Term 2: - L_jbic * tbar_ak * c_bjck 
!
!     X_ki = c_bjck L_bjci
!          
      call mem%alloc(X_ki, wf%n_o, wf%n_o) 
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  (wf%n_v**2)*wf%n_o,  & 
                  one,                 &
                  c_aibj,              & ! c_bjc,k
                  (wf%n_v**2)*wf%n_o,  & 
                  L_aibj,              & ! L_bjc,i
                  (wf%n_v**2)*wf%n_o,  & 
                  zero,                &
                  X_ki,                &
                  wf%n_o)
!
!     rho_ai -=  tbar_ak * X_ki 
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  tbar_ai,             & ! tbar_a,k
                  wf%n_v,              &
                  X_ki,                &
                  wf%n_o,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)

!
      call mem%dealloc(X_ki, wf%n_o, wf%n_o) 
!
!     Term 3: -  L_iajc * tbar_bk * c_bjck 
!
!     Reorder c_bjck to c_cjbk
!
      call mem%alloc(c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(c_aibj, c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_cj, wf%n_v, wf%n_o)
!
!     X_cj = c_cjbk t_bk
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  c_cjbk,              & ! c_cj,bk
                  (wf%n_v)*(wf%n_o),   &
                  tbar_ai,             & ! t_bk
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_cj,                &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(c_cjbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += L_aicj * X_cj 
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  1,                   &
                  (wf%n_v)*(wf%n_o),   &
                  -one,                &
                  L_aibj,              & ! L_ai,ck
                  (wf%n_v)*(wf%n_o),   &
                  X_cj,                &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ai,              &
                  (wf%n_v)*(wf%n_o))

!
      call mem%dealloc(X_cj, wf%n_v, wf%n_o)

!
!     Term 4: - L_jbka * tbar_ci * c_bjck 
!
!     X_ac = L_jbak * c_bjck
!
      call mem%alloc(X_ac, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  L_aibj,              &  ! L_a,kbj
                  wf%n_v,              &
                  c_aibj,              &  ! c_c,kbj
                  wf%n_v,              &
                  zero,                &
                  X_ac,                &
                  wf%n_v)                    
!
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
!     rho_ai += X_ac * tbar_ci
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v,              &
                  -one,                &
                  X_ac,                &
                  wf%n_v,              &
                  tbar_ai,             & ! tbar,ci
                  wf%n_v,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(X_ac, wf%n_v, wf%n_v)
!
   end subroutine F_ccsd_a1_1_ccsd
!
!
   module subroutine F_ccsd_a2_1_ccsd(wf, c_ai, rho_aibj, tbar_ai)
!!
!!    F transformation A2,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A2,1_aibj = 2 * L_iakc * tbar_bj * c_ck
!!                   - (L_jbic * tbar_ak + L_jbka * tbar_ci + L_kcib * tbar_aj) c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: tbar_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: L_iakc, g_iakc 
      real(dp), dimension(:,:,:,:), allocatable :: rho_jbia, rho_iajb
      real(dp), dimension(:,:,:,:), allocatable :: X_jbik
!
      real(dp), dimension(:,:), allocatable :: X_ia, X_ib, X_ik, c_kc 
!
      integer :: a, i, b, j
!
!     Construct L_iakc = 2 g_iakc - g_icka
!
      call mem%alloc(L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ovov',g_iakc)
!
      call copy_and_scale(two, g_iakc, L_iakc, wf%n_o**2*wf%n_v**2)
      call add_1432_to_1234(-one, g_iakc, L_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%dealloc(g_iakc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 1: 2 * L_iakc * tbar_bj * c_ck
!        
!     X_ia = 2 * L_iakc * c_ck 
!
      call mem%alloc(c_kc, wf%n_o, wf%n_v)
      call sort_12_to_21(c_ai, c_kc, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ia, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  (wf%n_o)*(wf%n_v),    &
                  1,                    &
                  (wf%n_o)*(wf%n_v),    &
                  two,                  &
                  L_iakc,               &
                  (wf%n_o)*(wf%n_v),    & 
                  c_kc,                 &
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  X_ia,                 &
                  (wf%n_o)*(wf%n_v))
!
!     rho_aibj =+ X_ia * tbar_bj
!
!$omp parallel do private(j,b,i,a)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v 
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + X_ia(i,a)*tbar_ai(b,j)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(X_ia, wf%n_o, wf%n_v)
!
!     :: Term 2: - L_jbic * tbar_ak * c_ck 
!     
!     X_jbik = - L_jbic * c_ck 
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
!     rho_jbia = X_jbik * tbar_ak
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
                  tbar_ai,              & ! tbar_a,k 
                  wf%n_v,               &
                  zero,                 &
                  rho_jbia,             & ! rho_jbi,a 
                  (wf%n_v)*(wf%n_o)**2)
!
      call mem%dealloc(X_jbik, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     rho_aibj =+ rho_jbia 
!
      call add_4321_to_1234(one, rho_jbia, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_jbia, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 3: - L_kajb * tbar_ci * c_ck
!
!     X_ik = -tbar_ci * c_ck
!
      call mem%alloc(X_ik, wf%n_o, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  wf%n_o,               &
                  wf%n_v,               &
                  -one,                 &
                  tbar_ai,              & ! tbar_c,i
                  wf%n_v,               &
                  c_ai,                 & ! c_c,k 
                  wf%n_v,               &
                  zero,                 &
                  X_ik,                 & 
                  wf%n_o)
!
!     rho_iajb = X_ik * L_kajb 
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
      call mem%dealloc(X_ik, wf%n_o, wf%n_o)
!
!     rho_aibj =+ rho_iajb 
!  
      call add_2143_to_1234(one, rho_iajb, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     :: Term 4: - L_kcib * tbar_aj * c_ck
!
!     X_ib = - c_ck * L_kcib  
!
      call mem%alloc(X_ib, wf%n_o, wf%n_v)
!
      call dgemm('N','N',               &
                  1,                    &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)*(wf%n_v),    &
                  -one,                 &
                  c_kc,                 & ! c_1,kc
                  1,                    &
                  L_iakc,               & ! L_kc,ib 
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  X_ib,                 & ! X_1,ib 
                  1)
!
!     rho_aibj =+ X_ib * tbar_aj
!
!$omp parallel do private(j,b,i,a)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o
               do a = 1, wf%n_v 
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + X_ib(i,b)*tbar_ai(a,j)
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
!!    rho_A1,1 = - (g_ibck * tbar_ajck + g_jack * tbar_bick) * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
!     :: Term 1: - g_ckib * tbar_ajck * c_bj 
!
      call mem%alloc(g_ckib, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('voov',g_ckib)
!
!     X_ckij = - g_ckib * c_bj 
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
!     rho_ai = tbar_ajck * X_ckij  
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
!     :: Term 2: - g_ckja * tbar_bick * c_bj
!
!     X_jick = - c_bj * tbar_bick 
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
!     rho_ia = X_jick * g_ckja 
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
!!    rho_B1,1 = - (g_ikcb * tbar_akcj + g_jkca * tbar_bkci) * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
!     Term 1: - g_ikcb * tbar_akcj * c_bj
!
      call mem%alloc(g_ikcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call wf%eri%get_eri_t1('oovv',g_ikcb)
!
!     X_ikcj = g_ikcb * c_bj
!
      call mem%alloc(X_ikcj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',              &
                  (wf%n_o**2)*(wf%n_v), &
                  wf%n_o,               &
                  wf%n_v,               &
                  one,                  &
                  g_ikcb,               & ! g_ikc,b
                  (wf%n_o**2)*(wf%n_v), &
                  c_ai,                 & ! c_b,j 
                  wf%n_v,               &
                  zero,                 &
                  X_ikcj,               &
                  (wf%n_o**2)*(wf%n_v))
!
!     rho_ai -= tbar_akcj * X_ikcj^T
!
      call dgemm('N', 'T',              &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o**2)*(wf%n_v), &
                  -one,                 &
                  tbar_aibj,            & ! tbar_a,kcj
                  wf%n_v,               &
                  X_ikcj,               & ! X_i,kcj
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               &
                  wf%n_v)
!
      call mem%dealloc(X_ikcj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2:  - g_jkca * tbar_bkci * c_bj
!
!     X_jkci = c_bj^T * tbar_bkci
!
      call mem%alloc(X_jkci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',              &
                  wf%n_o,               &
                  (wf%n_o**2)*wf%n_v,   &
                  wf%n_v,               &
                  one,                  &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  tbar_aibj,            & ! t_b,kci
                  wf%n_v,               &
                  zero,                 &
                  X_jkci,               & ! X_j,kci
                  wf%n_o)
!
!     rho_ai -= g_jkca^T * X_jkci    
!
      call dgemm('T', 'N',              &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o**2)*wf%n_v,   &
                  -one,                 &
                  g_ikcb,               & ! g_jkc,a
                  (wf%n_o**2)*wf%n_v,   &
                  X_jkci,               & ! X_jkc,i
                  (wf%n_o**2)*wf%n_v,   &
                  one,                  &
                  rho_ai,               &
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
!!    rho_C1,1 = (g_ikjl * tbar_akbl - g_cadb * tbar_cidj) * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
      integer :: req0, req1_c, req1_d, req2, current_c_batch, current_d_batch
      integer :: c, d, i, j
!
!     Term 1:  g_ikjl * tbar_blak * c_bj
!
!     X_jlak = c_bj^T * tbar_blak
!
      call mem%alloc(X_jlak, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_b,j
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_b,lak
                  wf%n_v,              &
                  zero,                &
                  X_jlak,              & ! X_j,lak
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
      call wf%eri%get_eri_t1('oooo',g_ikjl)
!
!     rho_ai += X_akjl * g_ikjl^T
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o**3,           &
                  one,                 & 
                  X_akjl,              &
                  wf%n_v,              &
                  g_ikjl,              &
                  wf%n_o,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(X_akjl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(g_ikjl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Term 2: - g_cadb * tbar_cidj * c_bj
!
      req0 = 0
!
      req1_c =  wf%eri%n_J*wf%n_v
      req1_d =  wf%eri%n_J*wf%n_v
!
      req2 = (wf%n_o)*(wf%n_v) + wf%n_v**2
!
      batch_c = batching_index(wf%n_v)
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, batch_d, req0, req1_c, req1_d, req2, &
                           tag='F_ccsd_c1_2_ccsd')
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
!           X_bcdi = c_bj * tbar_jcdi ( = tbar_cidj * c_bj )
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
                        tbar_jcdi,                                   & ! tbar_j,cdi
                        wf%n_o,                                      &
                        zero,                                        &
                        X_bcdi,                                      & ! X_b,cdi
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
            call wf%eri%get_eri_t1('vvvv',g_dbca,                               &
                                   batch_d%first, batch_d%get_last(),           &
                                   1, wf%n_v,                                   &
                                   batch_c%first, batch_c%get_last(),           &
                                   1, wf%n_v)
!
            call dgemm('T', 'N',                                     &
                        wf%n_v,                                      &
                        wf%n_o,                                      &
                        (batch_d%length)*(batch_c%length)*(wf%n_v),  &
                        one,                                         & 
                        g_dbca,                                      & ! g_dbc,a
                        (batch_d%length)*(batch_c%length)*(wf%n_v),  &
                        X_dbci,                                      & ! X_dbc,i
                        (batch_d%length)*(batch_c%length)*(wf%n_v),  &
                        one,                                         &
                        rho_ai,                                      &
                        wf%n_v)
!
            call mem%dealloc(g_dbca, batch_d%length, wf%n_v, batch_c%length, wf%n_v)
            call mem%dealloc(X_dbci, batch_d%length, wf%n_v, batch_c%length, wf%n_o)
!
         enddo ! batch_d
      enddo ! batch_c
!
      call mem%batch_finalize()
!
   end subroutine F_ccsd_c1_2_ccsd
!
!
   module subroutine F_ccsd_d1_2_ccsd(wf, c_aibj, rho_ai, tbar_aibj)
!!
!!    F transformation D1,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_D1,1 = (g_iljc * tbar_albk + g_jlic * tbar_blak + g_klja * tbar_clbi) * c_bjck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: tbar_aibj
!
!     Local variables
!
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
!     Term 1: g_iljc * tbar_albk * c_bjck
!
!     Reorder c_ckbj to c_jcbk 
!
      call mem%alloc(c_jcbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4132(c_aibj, c_jcbk , wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_iljc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ooov',g_iljc)
!
!     X_ilbk = g_iljc * c_jcbk
!
      call mem%alloc(X_ilbk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',                &
                  wf%n_o**2,              &
                  (wf%n_o)*(wf%n_v),      &
                  (wf%n_o)*(wf%n_v),      &
                  one,                    &
                  g_iljc,                 & ! g_il,jc
                  wf%n_o**2,              &
                  c_jcbk,                 & ! c_jc,bk
                  (wf%n_o)*(wf%n_v),      &
                  zero,                   &
                  X_ilbk,                 & ! X_il,bk
                  wf%n_o**2)
!
!     rho_ai += tbar_albk * X_ilbk^T
!
      call dgemm('N', 'T',                &
                  (wf%n_v),               &
                  (wf%n_o),               &
                  (wf%n_o**2)*(wf%n_v),   &
                  one,                    &
                  tbar_aibj,              & ! tbar_a,lbk
                  wf%n_v,                 &
                  X_ilbk,                 & ! X_i,lbk
                  wf%n_o,                 &
                  one,                    &
                  rho_ai,                 &
                  wf%n_v)
!
      call mem%dealloc(X_ilbk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2: g_jlic * tbar_blak * c_bjck
!
!     Reorder g_jlic to g_jcli
!
      call mem%alloc(g_jcli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1423(g_iljc, g_jcli, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     X_bkli = c_jcbk^T * g_jcli
!
      call mem%alloc(X_bkli, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',                &
                  (wf%n_o)*(wf%n_v),      &
                  wf%n_o**2,              &
                  (wf%n_o)*(wf%n_v),      &
                  one,                    &
                  c_jcbk,                 & ! c_jc,bk
                  (wf%n_o)*(wf%n_v),      &
                  g_jcli,                 & ! g_jc,li
                  (wf%n_o)*(wf%n_v),      &
                  zero,                   &
                  X_bkli,                 & ! X_bk,li 
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
!     rho_ai += tbar_blak * X_kbli 
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_o**2)*wf%n_v,     &
                  one,                    &
                  tbar_aibj,              & ! t_a,kbl
                  wf%n_v,                 &
                  X_kbli,                 & ! X_kbl,i
                  (wf%n_o**2)*wf%n_v,     &
                  one,                    &
                  rho_ai,                 &
                  wf%n_v)

!
      call mem%dealloc(X_kbli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(c_jcbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Term 3:  g_klja * tbar_clbi * c_bjck
!
!     Reorder tbar_clbi to tbar_cbli and c_ckbj to c_kjcb
!
      call mem%alloc(tbar_cbli, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(tbar_aibj, tbar_cbli, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_kjcb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2413(c_aibj, c_kjcb, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_kjli = c_kjcb * tbar_cbli
!
      call mem%alloc(X_kjli, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',                &
                  wf%n_o**2,              &
                  wf%n_o**2,              &
                  wf%n_v**2,              &
                  one,                    &
                  c_kjcb,                 & ! c_kj,cb
                  wf%n_o**2,              &
                  tbar_cbli,              & ! tbar_cb,li
                  wf%n_v**2,              &
                  zero,                   &
                  X_kjli,                 &
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
!     rho_ai += g_klja^T * X_klji
!
      call dgemm('T', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  wf%n_o**3,              &
                  one,                    &
                  g_iljc,                 & ! g_klj,a
                  wf%n_o**3,              &
                  X_klji,                 & ! X_klj,i
                  wf%n_o**3,              &
                  one,                    &
                  rho_ai,                 &
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
!!    rho_E1,1 = -(g_ibdc * tbar_ajdk + g_jadc * tbar_bidk + g_kbda * tbar_cjdi) * c_bjck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
      integer :: req0, req1, current_d_batch
      integer :: i, d, c, j, b, k
      type(batching_index) :: batch_d
!
!     :: Term 1: - g_ibdc * tbar_ajdk * c_bjck
!
!     X_jkid = - g_ibdc * c_bjck  
!
      call mem%alloc(c_jkbc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2413(c_aibj, c_jkbc, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_jkid, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(X_jkid, wf%n_v*(wf%n_o**3))
!
      req0 = (wf%eri%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%eri%n_J)*(wf%n_v) + 2*(wf%n_o)*(wf%n_v)**2
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_e1_2_ccsd 1')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_ibdc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv',g_ibdc,                       &
                                1, wf%n_o,                           &
                                1, wf%n_v,                           &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v) 
!
         call mem%alloc(g_bcid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
         call sort_1234_to_2413(g_ibdc, g_bcid, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call mem%dealloc(g_ibdc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call dgemm('N','N',                             &
                     (wf%n_o)**2,                        &
                     (wf%n_o)*(batch_d%length),          &
                     (wf%n_v)**2,                        &
                     -one,                               &
                     c_jkbc,                             & ! c_jk,bc 
                     (wf%n_o)**2,                        &
                     g_bcid,                             & ! g_bc,id
                     (wf%n_v)**2,                        &
                     one,                                &
                     X_jkid(1,1,1,batch_d%first),        & ! X_jk,id
                     (wf%n_o)**2)
!
         call mem%dealloc(g_bcid, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(c_jkbc, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     rho_ai += tbar_ajdk * X_jkid
!
      call mem%alloc(X_jdki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_1423(X_jkid, X_jdki, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(X_jkid, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N','N',                                &
                  wf%n_v,                                &
                  wf%n_o,                                &
                  (wf%n_v)*(wf%n_o)**2,                  &
                  one,                                   &
                  tbar_aibj,                             & ! tbar_a,jdk
                  wf%n_v,                                &
                  X_jdki,                                & ! X_jdk,i
                  (wf%n_v)*(wf%n_o)**2,                  &
                  one,                                   &
                  rho_ai,                                & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_jdki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2: - g_jadc * tbar_bidk * c_bjck
!
!     X_idjc = - tbar_bidk * c_bjck
! 
      call mem%alloc(tbar_idbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2314(tbar_aibj, tbar_idbk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(c_bkcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(c_aibj, c_bkcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_idcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N','N',                                &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  -one,                                  &
                  tbar_idbk,                             & ! tbar_id,bk
                  (wf%n_o)*(wf%n_v),                     &
                  c_bkcj,                                & ! c_bk,cj
                  (wf%n_o)*(wf%n_v),                     &
                  zero,                                  &
                  X_idcj,                                & ! X_id,cj
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(tbar_idbk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(c_bkcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += X_idcj * g_jadc
!
      req0 = (wf%eri%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%eri%n_J)*(wf%n_v) + (wf%n_o)*(wf%n_v)**2 + (wf%n_v)*(wf%n_o)**2
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_e1_2_ccsd 2')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_dcja, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%eri%get_eri_t1('vvov',g_dcja,                       &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v,                           &
                                1, wf%n_o,                           &
                                1, wf%n_v)
!
         call mem%alloc(X_idcj_red, wf%n_o, batch_d%length, wf%n_v, wf%n_o)
!
!$omp parallel do private(j,c,d,i)
         do j = 1, wf%n_o
            do c = 1, wf%n_v
               do d = 1, batch_d%length
                  do i = 1, wf%n_o
                     X_idcj_red(i, d, c, j) = X_idcj(i, batch_d%first-1+d, c, j)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
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
      call mem%batch_finalize()
!
      call mem%dealloc(X_idcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     :: Term 3: - g_kbda * tbar_cjdi * c_bjck
!
!     X_kbdi = - tbar_cjdi c_bjck

      call mem%alloc(c_kbcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_4132(c_aibj, c_kbcj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_kbdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call dgemm('N','N',                                &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  -one,                                  &
                  c_kbcj,                                & ! c_kb,cj 
                  (wf%n_o)*(wf%n_v),                     &
                  tbar_aibj,                             & ! tbar_cj,di
                  (wf%n_o)*(wf%n_v),                     &
                  zero,                                  &
                  X_kbdi,                                & ! X_kb,di
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(c_kbcj, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     rho_ai += X_kbdi * g_kbda 
!
      req0 = (wf%eri%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%eri%n_J)*(wf%n_v) + (wf%n_o)*(wf%n_v)**2 + (wf%n_v)*(wf%n_o)**2
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_e1_2_ccsd 3')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_kbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv',g_kbda,                       &
                                1, wf%n_o,                           &
                                1, wf%n_v,                           &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v)
!
         call mem%alloc(X_kbdi_red, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
!$omp parallel do private(k,b,d,i)
         do i = 1, wf%n_o
            do d = 1, batch_d%length
               do b = 1, wf%n_v
                  do k = 1, wf%n_o
                     X_kbdi_red(k, b, d, i) = X_kbdi(k, b, batch_d%first-1+d, i)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
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
      call mem%batch_finalize()
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
!!    rho_F1,2 = -(F_ib * tbar_ckaj + F_ja * tbar_ckbi) * c_bjck 
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ijck 
!
      real(dp), dimension(:,:), allocatable :: X_ij
!
!     :: Term 1: - F_ib * tbar_ajck * c_bjck 
!
!     X_ijck = - F_ib * c_bjck 
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
!     rho_ai += tbar_ajck * X_ijck
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
!     :: Term 2: - F_ja * tbar_ckbi * c_ckbj 
!
!     X_ij = - tbar_ckbi * c_ckbj 
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
                  c_aibj,               & ! c_ckb,j
                  (wf%n_o)*(wf%n_v)**2, &
                  zero,                 &
                  X_ij,                 & ! X_i,j
                  wf%n_o)
!
!     rho_ai += F_ja * X_ij
!
      call dgemm('T','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  wf%n_o,               &
                  one,                  & 
                  wf%fock_ia,           & ! F_j,a 
                  wf%n_o,               &
                  X_ij,                 & ! X_i,j
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               & ! rho_a,i
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
!!    rho_G1,2 = -(L_iljb * tbar_ckal + L_jlia * tbar_ckbl) * c_bjck 
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
!     :: Term 1: - L_iljb * tbar_alck * c_bjck 
!
!     Construct L_iljb ordered as L_ilbj = 2 g_iljb - g_jlib
!
      call mem%alloc(L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(g_iljb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ooov',g_iljb)
!
      call zero_array(L_ilbj, wf%n_v*(wf%n_o**3))
!
      call add_1243_to_1234(two, g_iljb, L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4213_to_1234(-one, g_iljb, L_ilbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_iljb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     X_ilck = - L_iljb * c_bjck
!
      call mem%alloc(X_ilck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_o)**2,          &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)*(wf%n_v),    &
                  -one,                 &
                  L_ilbj,               & ! L_il,bj 
                  (wf%n_o)**2,          &
                  c_aibj,               & ! c_bj,ck 
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  X_ilck,               &
                  (wf%n_o)**2)
!
!     rho_ai += tbar_alck * X_i,lck 
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
!     :: Term 2: - L_jlia * tbar_ckbl * c_ckbj 
!
!     X_jl = - c_ckbj * tbar_ckbl 
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
!     rho_ai += X_jl * L_jlia 
!
      call dgemm('T','N',               &
                  (wf%n_o)*(wf%n_v),    &
                  1,                    &
                  (wf%n_o)**2,          &
                  one,                  &
                  L_ilbj,               & ! L_jl,ai
                  (wf%n_o)**2,          &
                  X_jl,                 & 
                  (wf%n_o)**2,          &    
                  one,                  &
                  rho_ai,               &
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
!!    rho_H1,2 = (L_dajb * tbar_ckdi + L_dbia * tbar_ckdj) * c_bjck 
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
      integer         :: req0, req1, current_d_batch
      integer         :: j, b, d, i
      type(batching_index) :: batch_d
!
!     X_bjdi = c_bjck * tbar_ckdi 
!
      call mem%alloc(X_bjdi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',                                &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  one,                                   &
                  c_aibj,                                & ! c_bj,ck
                  (wf%n_o)*(wf%n_v),                     &
                  tbar_aibj,                             & ! tbar_ck,di
                  (wf%n_o)*(wf%n_v),                     &
                  zero,                                  &
                  X_bjdi,                                & ! X_bj,di
                  (wf%n_o)*(wf%n_v))
!
!     X_db   = tbar_djck * c_bjck
!
      call mem%alloc(X_db, wf%n_v, wf%n_v)
!
      call dgemm('N','T',                                &
                  wf%n_v,                                &
                  wf%n_v,                                &
                  (wf%n_v)*(wf%n_o)**2,                  &
                  one,                                   &
                  tbar_aibj,                             & ! tbar_d,jck
                  wf%n_v,                                &
                  c_aibj,                                & ! c_b,jck
                  wf%n_v,                                &
                  zero,                                  &
                  X_db,                                  & ! X_d,b
                  wf%n_v)
!
!     rho_ai += L_jbda * X_bjdi
!     rho_ia += L_iadb * X_db
!
      call mem%alloc(X_jbdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call sort_1234_to_2134(X_bjdi, X_jbdi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_bjdi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ia, wf%n_o, wf%n_v)
      call zero_array(rho_ia, wf%n_o*wf%n_v)
!
      req0 = (wf%eri%n_J)*(wf%n_o)*(wf%n_v)
      req1 = (wf%eri%n_J)*(wf%n_v) + 2*(wf%n_o)*(wf%n_v)**2          &
            + (wf%n_v)*(wf%n_o)**2 + wf%n_v
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_h1_2_ccsd')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv',g_jbda,                       &
                                1, wf%n_o,                           &
                                1, wf%n_v,                           &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v)
!
         call mem%alloc(L_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call copy_and_scale(two, g_jbda, L_jbda, wf%n_o*wf%n_v**2*batch_d%length)
         call add_1432_to_1234(-one, g_jbda, L_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call mem%dealloc(g_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
!        rho_ai from X_jbdi 
!
         call mem%alloc(red_X_jbdi, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
!$omp parallel do private(j,b,d,i)
         do i = 1, wf%n_o
            do d = 1, batch_d%length
               do b = 1, wf%n_v
                  do j = 1, wf%n_o
                     red_X_jbdi(j, b, d, i) = X_jbdi(j, b, batch_d%first-1+d, i)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
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
                     one,                                &
                     rho_ai,                             & ! rho_a,i
                     wf%n_v)
!
         call mem%dealloc(red_X_jbdi, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
!        rho_ia from X_db      
!
         call mem%alloc(red_X_db, batch_d%length, wf%n_v)
!
!$omp parallel do private(b,d)
         do b = 1, wf%n_v
            do d = 1, batch_d%length
               red_X_db(d, b) = X_db(batch_d%first-1+d, b)
            enddo
         enddo
!$omp end parallel do
!
         call dgemm('N','N',                             &
                     (wf%n_o)*(wf%n_v),                  &
                     1,                                  &
                     (wf%n_v)*(batch_d%length),          &
                     one,                                & 
                     L_jbda,                             & ! L_ia,db
                     (wf%n_o)*(wf%n_v),                  &
                     red_X_db,                           & ! X_db 
                     (wf%n_v)*(batch_d%length),          &
                     one,                                &
                     rho_ia,                             & ! rho_ia 
                     (wf%n_o)*(wf%n_v))
!
         call mem%dealloc(L_jbda, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
         call mem%dealloc(red_X_db, batch_d%length, wf%n_v)
!
      enddo 
!
      call mem%batch_finalize()
!
      call mem%dealloc(X_db, wf%n_v, wf%n_v)
      call mem%dealloc(X_jbdi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
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
!!    rho_I1,2 = - L_ldic * tbar_bjak * (t_bjck * c_dl + t_bjdl * c_ck) 
!!               - L_ldka * tbar_bjci * (t_bjck * c_dl + t_bjdl * c_ck)
!!               - L_ialc * tbar_bjdk * (t_bjck * c_dl + t_bjdl * c_ck)
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
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
!     Term 1: - L_ldic * tbar_bjak * t_bjck * c_dl  
!
!     Construct L_ldic = 2 * g_ldic - g_lcid (ordered as L_dlci)
!
      call mem%alloc(g_ldic, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov',g_ldic)
!
      call mem%alloc(L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_dlci, (wf%n_v**2)*(wf%n_o**2))
!
      call add_2143_to_1234(two, g_ldic, L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_ldic, L_dlci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ldic, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!  
!     X_ci = c_dl * L_dlci 
!
      call mem%alloc(X_ci, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  L_dlci,              & ! L_dl,ci
                  (wf%n_o)*(wf%n_v),   &
                  c_ai,                & ! c_dl
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_ci,                &
                  (wf%n_o)*(wf%n_v))
!
!     Y_ac = tbar_akbj * t_ckbj
!
      call mem%alloc(Y_ac, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             & 
                  wf%n_v,              &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  tbar_aibj,           & ! tbar_a,kbj
                  wf%n_v,              &
                  t_aibj,              & ! t_c,kbj
                  wf%n_v,              &
                  zero,                &
                  Y_ac,                &
                  wf%n_v)
!
!     rho_ai -= Y_ac * X_ci
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v,              &
                  -one,                &
                  Y_ac,                & ! Y_a,c
                  wf%n_v,              &
                  X_ci,                & ! X_c,i
                  wf%n_v,              &
                  one,                 &
                  rho_ai,              & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(Y_ac, wf%n_v, wf%n_v)
!
!     Term 2: - L_ldic * tbar_bjak * t_bjdl * c_ck
!
!     X_kidl = L_cidl * c_ck
!
      call mem%alloc(X_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_c,k 
                  wf%n_v,              &
                  L_dlci,              & ! L_c,idl
                  wf%n_v,              &
                  zero,                &
                  X_kidl,              &
                  wf%n_o)
!
!     Y_kibj =  X_kidl * t_dlbj
!
      call mem%alloc(Y_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  X_kidl,              & ! X_ki,dl
                  wf%n_o**2,           &
                  t_aibj,              & ! t_dl,bj
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
!     rho_ai -= tbar_akbj * Y_kbji 
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  -one,                &
                  tbar_aibj,           & ! t_a,kbj 
                  wf%n_v,              &
                  Y_kbji,              & ! Y_kbj,i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(Y_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Term 3: - L_ldka * tbar_bjci * t_bjck * c_dl 
!
!     Y_ki = t_bjck * tbar_bjci 
!
      call mem%alloc(Y_ki, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  (wf%n_v**2)*(wf%n_o),&
                  one,                 &
                  t_aibj,              & ! t_bjc,k
                  (wf%n_v**2)*(wf%n_o),&
                  tbar_aibj,           & ! tbar_bjc,i
                  (wf%n_v**2)*(wf%n_o),&
                  zero,                &
                  Y_ki,                & ! Y_k,i
                  wf%n_o)
!
!     rho_ai -= X_ak * Y_ki
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  X_ci,                & ! X_a,k
                  wf%n_v,              &
                  Y_ki,                &
                  wf%n_o,              &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(Y_ki, wf%n_o, wf%n_o)
      call mem%dealloc(X_ci, wf%n_v, wf%n_o)
!
!     Term 4:  - L_ldka * tbar_bjci * t_bjdl * c_ck
!
!     X_kibj = c_ck * tbar_cibj
!
      call mem%alloc(X_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*(wf%n_v),&
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_c,k
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_c,ibj
                  wf%n_v,              &
                  zero,                &
                  X_kibj,              &
                  wf%n_o)
!
!     Y_kidl = X_kibj * t_bjdl
!
      call mem%alloc(Y_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  X_kibj,              & ! X_ki,bj
                  wf%n_o**2,           &
                  t_aibj,              & ! t_bj,dl
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  Y_kidl,              & ! Y_ki,dl
                  wf%n_o**2)
!
      call mem%dealloc(X_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder Y_kidl to Y_kdli
!
      call mem%alloc(Y_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1342(Y_kidl, Y_kdli, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(Y_kidl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     rho -= L_akdl * Y_kdli
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  -one,                &
                  L_dlci,              & ! L_a,kdl
                  wf%n_v,              &
                  Y_kdli,              & ! Y_kdl,i
                  (wf%n_o**2)*wf%n_v,  &
                  one,                 &
                  rho_ai,              &
                  wf%n_v)
!
      call mem%dealloc(Y_kdli, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     Term 5: - L_ialc * t_ckbj * tbar_dkbj * c_dl
!
!     X_cd = t_ckbj * tbar_dkbj
!
      call mem%alloc(X_cd, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  wf%n_o**2 * wf%n_v,  &
                  one,                 &
                  t_aibj,              & ! t_c,kbj
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_d,kbj
                  wf%n_v,              &
                  zero,                &
                  X_cd,                &
                  wf%n_v)
!
!     Y_cl =  X_cd * c_dl
!
      call mem%alloc(Y_cl, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_v,              &
                  one,                 &
                  X_cd,                & ! X_c,d
                  wf%n_v,              &
                  c_ai,                & ! c_d,l
                  wf%n_v,              &
                  zero,                &
                  Y_cl,                &
                  wf%n_v)
!
      call mem%dealloc(X_cd, wf%n_v, wf%n_v)
!
!     rho_ai -= L_ialc * Y_cl
!
      call dgemm('N', 'N',             &
                  wf%n_v*wf%n_o,       &
                  1,                   &
                  wf%n_v*wf%n_o,       &
                  -one,                &
                  L_dlci,              & ! L_ai,cl
                  wf%n_v*wf%n_o,       &
                  Y_cl,                &
                  wf%n_v*wf%n_o,       &
                  one,                 &
                  rho_ai,              &
                  wf%n_v*wf%n_o)
!
      call mem%dealloc(Y_cl, wf%n_v, wf%n_o)
!
!     Term 6:  - L_ialc * tbar_bjdk * t_bjdl * c_ck 
!
!     X_kl = tbar_bjdk * t_bjdl
!
      call mem%alloc(X_kl, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  (wf%n_v**2)*wf%n_o,  &
                  one,                 &
                  tbar_aibj,           & ! tbar_bjd,k
                  (wf%n_v**2)*wf%n_o,  &
                  t_aibj,              & ! t_bjd,l
                  (wf%n_v**2)*wf%n_o,  &
                  zero,                &
                  X_kl,                &
                  wf%n_o)
!
!     Y_cl = c_ck * X_kl
!
      call mem%alloc(Y_cl, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v,              &
                  wf%n_o,              &
                  wf%n_o,              &
                  one,                 &
                  c_ai,                & ! c_c,k
                  wf%n_v,              &
                  X_kl,                & ! X_k,l
                  wf%n_o,              &
                  zero,                &
                  Y_cl,                & ! Y_c,l
                  wf%n_v)
! 
      call mem%dealloc(X_kl, wf%n_o, wf%n_o)
!
!     rho_ai -= L_ialc * Y_cl 
!
      call dgemm('N', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  -one,                &
                  L_dlci,              & ! L_ai,cl
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
!!    rho_J1,2 = (g_kbid * tbar_cjal + g_jcid * tbar_bkal + g_kajd * tbar_cibl 
!!              + g_jakd * tbar_bicl + g_ibkd * tbar_ajcl + g_lakb * tbar_dicj) * t_ckdl * c_bj
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                 :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: rho_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: t_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_idkb, g_ijcd, tbar_cdji, X_ajkl, X_aklj, X_ijcl
      real(dp), dimension(:,:,:,:), allocatable :: X_idkj, t_cdkl, t_cldk, t_dkcl, t_kdcl, t_klcd
      real(dp), dimension(:,:,:,:), allocatable :: X_ijdk, X_jicl, X_jkal, X_jlci, X_kdij
      real(dp), dimension(:,:,:,:), allocatable :: Y_ijcl, Y_ijdk, Y_ijkl, Y_kdji, Y_lcji, Y_klji
!
      call mem%alloc(g_idkb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov',g_idkb)
!
!     :: Term 1: g_idkb * tbar_alcj * t_ckdl * c_bj
!
!     X_idk,j = g_idkb * c_bj 
!
      call mem%alloc(X_idkj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  wf%n_v,               &
                  one,                  &
                  g_idkb,               & ! g_idk,b
                  (wf%n_v)*(wf%n_o)**2, &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  zero,                 &
                  X_idkj,               & ! X_idk,j
                  (wf%n_v)*(wf%n_o)**2)
!
!     Y_ijcl = X_idkj * t_ckdl
!
      call mem%alloc(X_ijdk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1423(X_idkj, X_ijdk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X_idkj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(t_dkcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3214(t_aibj, t_dkcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_o)**2,          &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)*(wf%n_v),    &
                  one,                  &
                  X_ijdk,               & ! X_ij,dk 
                  (wf%n_o)**2,          &
                  t_dkcl,               & ! t_dk,cl
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  Y_ijcl,               & ! Y_ij,cl 
                  (wf%n_o)**2)
!
      call mem%dealloc(X_ijdk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += tbar_alcj * Y_ijcl
!
      call mem%alloc(Y_lcji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_4321(Y_ijcl, Y_lcji, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(Y_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,lcj
                  wf%n_v,               &
                  Y_lcji,               & ! Y_lcj,i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(Y_lcji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2: g_jcid * tbar_bkal * t_ckdl * c_bj
!
!     X_jkal = c_bj * tbar_bkal 
!
      call mem%alloc(X_jkal, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  c_ai,                 & ! c_b,j 
                  wf%n_v,               &
                  tbar_aibj,            & ! tbar_b,kal 
                  wf%n_v,               &
                  zero,                 &
                  X_jkal,               & ! X_j,kal 
                  wf%n_o)
!
!     Y_ijkl = g_jcid * t_ckdl
!
      call mem%alloc(g_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3124(g_idkb, g_ijcd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(t_cdkl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(t_aibj, t_cdkl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_o)**2,          &
                  (wf%n_o)**2,          &
                  (wf%n_v)**2,          &
                  one,                  &
                  g_ijcd,               & ! g_ij,cd 
                  (wf%n_o)**2,          &
                  t_cdkl,               & ! t_cd,kl
                  (wf%n_v)**2,          &
                  zero,                 &
                  Y_ijkl,               & ! Y_ij,kl 
                  (wf%n_o)**2)
!
      call mem%dealloc(t_cdkl, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_ijcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     rho_ai += X_jkal * Y_ijkl
!
      call mem%alloc(X_ajkl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_1234_to_3124(X_jkal, X_ajkl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_jkal, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o)**3,          &
                  one,                  &
                  X_ajkl,               & ! X_a,jkl
                  wf%n_v,               &
                  Y_ijkl,               & ! Y_i,jkl
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_ajkl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(Y_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     :: Term 3: g_jdka * tbar_blci * t_ckdl * c_bj
!
!     X_jlci = c_bj * tbar_blci 
!
      call mem%alloc(X_jlci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  c_ai,                 & ! c_b,j 
                  wf%n_v,               &
                  tbar_aibj,            & ! tbar_b,lci 
                  wf%n_v,               &
                  zero,                 &
                  X_jlci,               & ! X_j,lci
                  wf%n_o)
!
!     Y_ijdk = t_ckdl * X_jlci 
!
      call mem%alloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(t_aibj, t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(X_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_4132(X_jlci, X_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_jlci, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_ijdk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_o)**2,          &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)*(wf%n_v),    &
                  one,                  &
                  X_ijcl,               & ! X_ij,cl
                  (wf%n_o)**2,          &
                  t_cldk,               & ! t_cl,dk
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  Y_ijdk,               & ! Y_ij,dk
                  (wf%n_o)**2)
!
      call mem%dealloc(t_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += g_jdka * Y_ijdk 
!
      call dgemm('T','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  g_idkb,               & ! g_jdk,a 
                  (wf%n_v)*(wf%n_o)**2, &
                  Y_ijdk,               & ! Y_i,jdk 
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               & ! rho_a,i 
                  wf%n_v)
!
      call mem%dealloc(Y_ijdk, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 4: g_kdja * t_ckdl * c_bj * tbar_bicl
!
!     X_jicl = c_bj * tbar_bicl 
!
      call mem%alloc(X_jicl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_v,               &
                  one,                  &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  tbar_aibj,            & ! tbar_b,icl
                  wf%n_v,               &
                  zero,                 &
                  X_jicl,               & ! X_j,icl
                  wf%n_o)
!
!     Y_kdji = t_ckdl * X_jicl
!
      call mem%alloc(Y_kdji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(t_kdcl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call sort_1234_to_2314(t_aibj, t_kdcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','T',               &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)**2,          &
                  (wf%n_o)*(wf%n_v),    &
                  one,                  &
                  t_kdcl,               & ! t_kd,cl
                  (wf%n_o)*(wf%n_v),    &
                  X_jicl,               & ! X_ji,cl
                  (wf%n_o)**2,          &
                  zero,                 &
                  Y_kdji,               & ! Y_kd,ji
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_jicl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai += g_kdja * Y_kdji
!
      call dgemm('T','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  g_idkb,               & ! g_kdj,a
                  (wf%n_v)*(wf%n_o)**2, &
                  Y_kdji,               & ! Y_kdj,i
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(Y_kdji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 5: g_kdib * tbar_ajcl * t_ckdl * c_bj
!
!     X_kdij = g_kdib * c_bj 
!
      call mem%alloc(X_kdij, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  wf%n_v,               &
                  one,                  &
                  g_idkb,               & ! g_kdi,b
                  (wf%n_v)*(wf%n_o)**2, &
                  c_ai,                 & ! c_b,j
                  wf%n_v,               &
                  zero,                 &
                  X_kdij,               & ! X_kdi,j
                  (wf%n_v)*(wf%n_o)**2)
!
!     Y_ijcl = X_kdij * t_ckdl
!
      call mem%alloc(Y_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',               &
                  (wf%n_o)**2,          &
                  (wf%n_o)*(wf%n_v),    &
                  (wf%n_o)*(wf%n_v),    &
                  one,                  &
                  X_kdij,               & ! X_kd,ij
                  (wf%n_o)*(wf%n_v),    &
                  t_kdcl,               & ! t_kd,cl 
                  (wf%n_o)*(wf%n_v),    &
                  zero,                 &
                  Y_ijcl,               & ! Y_ij,cl 
                  (wf%n_o)**2)
!
      call mem%dealloc(t_kdcl, wf%n_o, wf%n_v, wf%n_v, wf%n_o) 
!
!     rho_ai += tbar_ajcl * Y_ijcl
!
      call dgemm('N','T',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  one,                  &
                  tbar_aibj,            & ! tbar_a,jcl
                  wf%n_v,               &
                  Y_ijcl,               & ! Y_i,jcl
                  wf%n_o,               &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(Y_ijcl, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 6: g_lakb * tbar_dicj * t_ckdl * c_bj
!
!     X_lakj = g_lakb * c_bj  (X kept from previous term, X_kdij)
!
!     Y_klji = t_ckdl * tbar_dicj
!
      call mem%alloc(Y_klji, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%alloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_2413(t_aibj, t_klcd, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(tbar_cdji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_3142(tbar_aibj, tbar_cdji, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',               &
                  (wf%n_o)**2,          &
                  (wf%n_o)**2,          &
                  (wf%n_v)**2,          &
                  one,                  &
                  t_klcd,               & ! t_kl,cd
                  (wf%n_o)**2,          &
                  tbar_cdji,            & ! tbar_cd,ji
                  (wf%n_v)**2,          &
                  zero,                 &
                  Y_klji,               & ! Y_kl,ji
                  (wf%n_o)**2)
!
      call mem%dealloc(tbar_cdji, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_klcd, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!     rho_ai += X_lakj * Y_klji
!
      call mem%alloc(X_aklj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_1234_to_2314(X_kdij, X_aklj, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X_kdij, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',               &
                  wf%n_v,               &
                  wf%n_o,               &
                  (wf%n_o)**3,          &
                  one,                  &
                  X_aklj,               & ! X_a,klj
                  wf%n_v,               &
                  Y_klji,               & ! Y_klj,i
                  (wf%n_o)**3,          &
                  one,                  &
                  rho_ai,               & ! rho_a,i
                  wf%n_v)
!
      call mem%dealloc(X_aklj, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(Y_klji, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(t_dkcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_idkb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
   end subroutine F_ccsd_j1_2_ccsd
!
!
   module subroutine F_ccsd_a2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation a2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A2,2 = (g_iljc * tbar_albk + g_klja * tbar_clbi + g_ilkb * tbar_cjal) * c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iljc, g_bikl
      real(dp), dimension(:,:,:,:), allocatable :: X_iljk, X_lkij, X_klaj, X_klbi
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ablk
      real(dp), dimension(:,:,:,:), allocatable :: rho_abij, rho_jabi, rho_biaj
!
!     Term 1: g_iljc * tbar_albk * c_ck
!
      call mem%alloc(g_iljc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ooov',g_iljc)
!
!     X_iljk = g_iljc * c_ck
!
      call mem%alloc(X_iljk, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call dgemm('N', 'N',             &
                  wf%n_o**3,           &
                  wf%n_o,              &
                  wf%n_v,              &
                  one,                 &
                  g_iljc,              & ! g_ilj,c
                  wf%n_o**3,           &
                  c_ai,                & ! c_c,k
                  wf%n_v,              &
                  zero,                &
                  X_iljk,              &
                  wf%n_o**3)
!
!     Reorder X_iljk to X_lkij
!
      call mem%alloc(X_lkij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call sort_1234_to_2413(X_iljk, X_lkij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(X_iljk, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Reorder tbar_albk to tbar_ablk
!
      call mem%alloc(tbar_ablk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(tbar_aibj, tbar_ablk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += tbar_lkab * X_lkij
!
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_v**2,           &
                  wf%n_o**2,           &
                  wf%n_o**2,           &
                  one,                 &
                  tbar_ablk,           &
                  wf%n_v**2,           &
                  X_lkij,              &
                  wf%n_o**2,           &
                  zero,                &
                  rho_abij,            &
                  wf%n_v**2) 
!
      call mem%dealloc(tbar_ablk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(X_lkij, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call add_1324_to_1234(one, rho_abij, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Term 2: g_klja * tbar_clbi * c_ck
!
!     X_klbi = c_ck * tbar_clbi
!
      call mem%alloc(X_klbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_c,k
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_c,lbi
                  wf%n_v,              &
                  zero,                &
                  X_klbi,              & ! X_k,lbi
                  wf%n_o)
!
      call mem%alloc(rho_jabi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     rho_aibj += g_klja * X_klbi
!
      call dgemm('T', 'N',             &
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  wf%n_o**2,           &
                  one,                 &
                  g_iljc,              & ! g_kl,ja
                  wf%n_o**2,           &
                  X_klbi,              & ! X_kl,bi
                  wf%n_o**2,           &
                  zero,                &
                  rho_jabi,            &
                  (wf%n_o)*(wf%n_v))
!
      call add_4132_to_1234(one, rho_jabi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_jabi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Term 3:  g_ilkb * tbar_alcj * c_ck
!
!     Reorder X_kjal to X_klaj
!
      call mem%alloc(X_klaj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(X_klbi, X_klaj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_klbi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder g_ilkb to g_bikl
!
      call mem%alloc(g_bikl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call sort_1234_to_4132(g_iljc, g_bikl, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(g_iljc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%alloc(rho_biaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += g_bikl * X_klaj
!
      call dgemm('N', 'N',             & 
                  (wf%n_o)*(wf%n_v),   &
                  (wf%n_o)*(wf%n_v),   &
                  wf%n_o**2,           &
                  one,                 &
                  g_bikl,              & ! g_bi,kl
                  (wf%n_o)*(wf%n_v),   &
                  X_klaj,              & ! X_kl,aj
                  (wf%n_o)**2,         &
                  zero,                &
                  rho_biaj,            &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_klaj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_bikl, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call add_3214_to_1234(one, rho_biaj, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_biaj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_a2_2_ccsd
!
!
   module subroutine F_ccsd_b2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation b2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_B2,2 = -(g_ibdc * tbar_ajdk + g_kbda * tbar_cjdi + g_dbic * tbar_akdj) * c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ibdk, g_dbic, g_ibdc, Y_dbik, red_X_ibdk
      real(dp), dimension(:,:,:,:), allocatable :: g_abkd, red_Y_dbik, red_Z_kdij, Z_kdij, rho_abij
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi, tbar_ajdk, X_dkbi, Y_dkbi, Z_kjdi
!
      integer         :: req0, req1, current_d_batch
      integer         :: i, b, d, k, j
      type(batching_index) :: batch_d
!
!     Construct intermediate for term 2: Z_kjdi = - c_ck * tbar_cjdi
!
      call mem%alloc(Z_kjdi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T','N',                                &
                  wf%n_o,                                &
                  (wf%n_v)*(wf%n_o)**2,                  &
                  wf%n_v,                                &
                  -one,                                  &
                  c_ai,                                  & ! c_c,k
                  wf%n_v,                                &
                  tbar_aibj,                             & ! tbar_c,jdi
                  wf%n_v,                                &
                  zero,                                  &
                  Z_kjdi,                                & ! Z_k,jdi 
                  wf%n_o)
!
      call mem%alloc(Z_kdij, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1342(Z_kjdi, Z_kdij, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(Z_kjdi, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct intermediates for terms 1 and 3:
!
      call mem%alloc(X_ibdk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%alloc(Y_dbik, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%alloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call zero_array(X_ibdk, (wf%n_v**2)*(wf%n_o**2))
      call zero_array(Y_dbik, (wf%n_v**2)*(wf%n_o**2))
      call zero_array(rho_abij, (wf%n_v**2)*(wf%n_o**2))
!
      req0 = (wf%n_o)*(wf%n_v)*(wf%eri%n_J)
      req1 = max((wf%n_v)*(wf%eri%n_J) + (wf%n_o)*(wf%n_v)**2,       & ! L_dc^J, g_ibdc
                  2*(wf%n_o)*(wf%n_v)**2,                            & ! 2 x g_ibdc
                  (wf%n_o)*(wf%n_v)**2 + (wf%n_v)*(wf%n_o)**2)         ! g_ibdc, Y_dbik
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_b2_2_ccsd')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
         call mem%alloc(g_ibdc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv',g_ibdc,                       & 
                                1, wf%n_o,                           &
                                1, wf%n_v,                           &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v)
!
!        X_ibdk += - g_ibdc * c_ck
!
         call mem%alloc(red_X_ibdk, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     -one,                               &
                     g_ibdc,                             & ! g_ibd,c
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     c_ai,                               & ! c_c,k
                     wf%n_v,                             &
                     zero,                               &
                     red_X_ibdk,                         & ! X_ibd,k 
                     (wf%n_o)*(wf%n_v)*(batch_d%length))
!
!$omp parallel do private(i,b,d,k)
         do k = 1, wf%n_o
            do d = 1, batch_d%length
               do b = 1, wf%n_v
                  do i = 1, wf%n_o
                     X_ibdk(i, b, batch_d%first-1+d, k) = X_ibdk(i, b, batch_d%first-1+d, k) + red_X_ibdk(i, b, d, k)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(red_X_ibdk, wf%n_o, wf%n_v, batch_d%length, wf%n_o)
!
!        Y_dbik += - g_icdb * c_ck   
!
         call mem%alloc(g_dbic, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
         call sort_1234_to_3412(g_ibdc, g_dbic, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
         call mem%dealloc(g_ibdc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call mem%alloc(red_Y_dbik, batch_d%length, wf%n_v, wf%n_o, wf%n_o)
!
         call dgemm('N','N',                             &
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     -one,                               &
                     g_dbic,                             & ! g_dbi,c 
                     (wf%n_o)*(wf%n_v)*(batch_d%length), &
                     c_ai,                               & ! c_c,k
                     wf%n_v,                             &
                     zero,                               &
                     red_Y_dbik,                         & ! X_dbi,k
                     (wf%n_o)*(wf%n_v)*(batch_d%length))
!
!$omp parallel do private(d,b,i,k)
         do k = 1, wf%n_o
            do i = 1, wf%n_o
               do b = 1, wf%n_v
                  do d = 1, batch_d%length
                     Y_dbik(batch_d%first-1+d, b, i, k) = Y_dbik(batch_d%first-1+d, b, i, k) + red_Y_dbik(d, b, i, k)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(red_Y_dbik, batch_d%length, wf%n_v, wf%n_o, wf%n_o)
!
!        For term 2: contruct rho_abij += g_kbda * Z_kjdi      
!
         call mem%alloc(g_abkd, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
         call sort_1234_to_2431(g_dbic, g_abkd, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
         call mem%dealloc(g_dbic, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%alloc(red_Z_kdij, wf%n_o, batch_d%length, wf%n_o, wf%n_o)
!$omp parallel do private(k, d, i, j)
         do j = 1, wf%n_o
            do i = 1, wf%n_o
               do d = 1, batch_d%length
                  do k = 1, wf%n_o
                     red_Z_kdij(k, d, i, j) = Z_kdij(k, batch_d%first-1+d, i, j)
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call dgemm('N','N',                             &
                     (wf%n_v)**2,                        &
                     (wf%n_o)**2,                        &
                     (wf%n_o)*(batch_d%length),          &
                     one,                                &
                     g_abkd,                             & ! g_ab,kd 
                     (wf%n_v)**2,                        &
                     red_Z_kdij,                         & ! Z_kd,ij 
                     (wf%n_o)*(batch_d%length),          &
                     one,                                &
                     rho_abij,                           & ! rho_ab,ij 
                     (wf%n_v)**2)
!
         call mem%dealloc(red_Z_kdij, wf%n_o, batch_d%length, wf%n_o, wf%n_o)
         call mem%dealloc(g_abkd, wf%n_v, wf%n_v, wf%n_o, batch_d%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(Z_kdij, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 2: rho_aibj += rho_abij 
!
      call add_1324_to_1234(one, rho_abij, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     :: Term 1: rho_aibj += tbar_ajdk * X_ibdk
!
      call mem%alloc(X_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_3421(X_ibdk, X_dkbi, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(X_ibdk, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N','N',                                &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  one,                                   &
                  tbar_aibj,                             & ! tbar_aj,dk
                  (wf%n_o)*(wf%n_v),                     &
                  X_dkbi,                                & ! X_dk,bi 
                  (wf%n_o)*(wf%n_v),                     &
                  zero,                                  &
                  rho_ajbi,                              & ! rho_aj,bi
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(X_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 3: rho_aibj += tbar_akdj * Y_dbik
!
      call mem%alloc(tbar_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(tbar_aibj, tbar_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(Y_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1423(Y_dbik, Y_dkbi, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(Y_dbik, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dgemm('N','N',                                &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  (wf%n_o)*(wf%n_v),                     &
                  one,                                   &
                  tbar_ajdk,                             & ! tbar_aj,dk
                  (wf%n_o)*(wf%n_v),                     &
                  Y_dkbi,                                & ! Y_dk,bi
                  (wf%n_o)*(wf%n_v),                     &
                  one,                                   &
                  rho_ajbi,                              & ! rho_aj,bi
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(Y_dkbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(tbar_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_b2_2_ccsd
!
!
   module subroutine F_ccsd_c2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation c2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_C2,2 = -(F_jc * tbar_aibk + F_ka * tbar_bjci) * c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable     :: X_jk
      real(dp), dimension(:,:,:,:), allocatable :: X_kibj
!
!     Term 1: -F_jc * tbar_aibk * c_ck
!
!     X_jk = F_jc * c_ck
!
      call mem%alloc(X_jk, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  wf%n_v,              &
                  one,                 &
                  wf%fock_ia,          & ! F_j,c
                  wf%n_o,              &
                  c_ai,                & ! c_c,k
                  wf%n_v,              &
                  zero,                &
                  X_jk,                &
                  wf%n_o)
!
!     rho_aibj -= tbar_aibk * X_jk
!
      call dgemm('N', 'T',             &
                  (wf%n_v**2)*wf%n_o,  &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  tbar_aibj,           & ! tbar_aib,k
                  (wf%n_v**2)*wf%n_o,  &
                  X_jk,                &
                  wf%n_o,              &
                  one,                 &
                  rho_aibj,            &
                  (wf%n_v**2)*wf%n_o)
!
      call mem%dealloc(X_jk, wf%n_o, wf%n_o)
!
!     Term 2: - F_ka * tbar_bjci * c_ck
!
!     X_kibj = c_ck * tbar_bjci
!
      call mem%alloc(X_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                 wf%n_o,               &
                 (wf%n_o**2)*wf%n_v,   &
                 wf%n_v,               &
                 one,                  &
                 c_ai,                 & ! c_c,k 
                 wf%n_v,               &
                 tbar_aibj,            & ! tbar_c,ibj
                 wf%n_v,               &
                 zero,                 &
                 X_kibj,               &
                 wf%n_o)
!
!     rho_aibj -= F_ka * X_kibj
!
      call dgemm('T', 'N',             &
                  wf%n_v,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_o,              &
                  -one,                &
                  wf%fock_ia,          & ! F_k,a
                  wf%n_o,              &
                  X_kibj,              &
                  wf%n_o,              &
                  one,                 &
                  rho_aibj,            &
                  wf%n_v)

      call mem%dealloc(X_kibj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_c2_2_ccsd
!
!
   module subroutine F_ccsd_d2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation d2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_D2,2 = -(L_jlkc * tbar_aibl + L_kljb * tbar_aicl) * c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable     :: X_lj
      real(dp), dimension(:,:,:,:), allocatable :: X_klai
      real(dp), dimension(:,:,:,:), allocatable :: g_jlkc
      real(dp), dimension(:,:,:,:), allocatable :: L_ljck, L_klbj
!
!     Term 1: - L_jlkc * tbar_aibl * c_ck
!
!     L_jlkc = 2 g_jlkc - g_kljc (ordered as L_ljck)
!
      call mem%alloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ooov',g_jlkc)
!
      call mem%alloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(L_ljck, wf%n_v*(wf%n_o**3))
!
      call add_2143_to_1234(two, g_jlkc, L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call add_4123_to_1234(-one, g_jlkc, L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
!     X_lj = L_ljck * c_ck 
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemm('N', 'N',             &
                  wf%n_o**2,           &
                  1,                   &
                  (wf%n_o)*(wf%n_v),   &
                  one,                 &
                  L_ljck,              & ! L_lj,ck
                  wf%n_o**2,           &
                  c_ai,                & ! c_ck
                  (wf%n_o)*(wf%n_v),   &
                  zero,                &
                  X_lj,                &
                  wf%n_o**2)
!
!     rho_aibj -= tbar_aibl * X_lj
!
      call dgemm('N', 'N',             &
                  (wf%n_v**2)*wf%n_o,  &
                  wf%n_o,              &
                  wf%n_o,              &
                  -one,                &
                  tbar_aibj,           & ! t_aib,l
                  (wf%n_v**2)*wf%n_o,  &
                  X_lj,                &
                  wf%n_o,              &
                  one,                 &
                  rho_aibj,            & ! rho_aib,j
                  (wf%n_v**2)*wf%n_o) 
!
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
!     Term 2: - L_kljb * tbar_aicl * c_ck
!
!     X_klai = c_ck * tbar_aicl
!
      call mem%alloc(X_klai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  (wf%n_o**2)*wf%n_v,  &
                  wf%n_v,              &
                  one,                 &
                  c_ai,                & ! c_c,k
                  wf%n_v,              &
                  tbar_aibj,           & ! tbar_c,lai
                  wf%n_v,              &
                  zero,                &
                  X_klai,              & ! X_k,lai
                  wf%n_o)
!
!     Reorder L_lkbj to L_klbj
!
      call mem%alloc(L_klbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_2134(L_ljck, L_klbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(L_ljck, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj -= X_klai * L_klbj
!
      call dgemm('T', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  wf%n_o**2,           &
                  -one,                &
                  X_klai,              &
                  wf%n_o**2,           &
                  L_klbj,              &
                  wf%n_o**2,           &
                  one,                 &
                  rho_aibj,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_klai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_klbj, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_d2_2_ccsd
!
!
   module subroutine F_ccsd_e2_2_ccsd(wf, c_ai, rho_aibj, tbar_aibj)
!!
!!    F transformation e2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_E2,2 = (L_dakc * tbar_bjdi + Ldcia * tbar_bjdk) * c_ck
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:), allocatable :: X_ad
      real(dp), dimension(:,:,:,:), allocatable :: X_aidk
      real(dp), dimension(:,:,:,:), allocatable :: L_adck, L_aidc
      real(dp), dimension(:,:,:,:), allocatable :: g_dakc, g_iadc
!
      integer         :: req0, req1, current_d_batch
      type(batching_index) :: batch_d
!
!     Term 1: L_dakc * tbar_bjdi * c_ck
!
!     X_ad = L_dakc * c_ck 
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
      call zero_array(X_ad, wf%n_v**2)
!
      req0 = (wf%eri%n_J)*(wf%n_v)*(wf%n_o)
      req1 = max((wf%eri%n_J)*(wf%n_v) + (wf%n_v**2)*(wf%n_o), 2*(wf%n_v**2)*(wf%n_o))
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_e2_2_ccsd 1')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
!        L_dakc = 2 g_dakc - g_dcka (ordered as L_adck)
!
         call mem%alloc(g_dakc, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%eri%get_eri_t1('vvov',g_dakc,                       &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v,                           &
                                1, wf%n_o,                           &
                                1, wf%n_v)
!
         call mem%alloc(L_adck, wf%n_v, batch_d%length, wf%n_v, wf%n_o)
         call zero_array(L_adck, (wf%n_v**2)*batch_d%length*wf%n_o)
         call add_2143_to_1234(two, g_dakc, L_adck, wf%n_v, batch_d%length, wf%n_v, wf%n_o)
         call add_2341_to_1234(-one, g_dakc, L_adck, wf%n_v, batch_d%length, wf%n_v, wf%n_o)
!
         call mem%dealloc(g_dakc, batch_d%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('N', 'N',                          &
                     wf%n_v*batch_d%length,            &
                     1,                                &
                     (wf%n_o)*(wf%n_v),                &
                     one,                              &
                     L_adck,                           & ! L_ad,ck
                     wf%n_v*batch_d%length,            &
                     c_ai,                             & ! c_ck
                     (wf%n_o)*(wf%n_v),                &
                     one,                              &
                     X_ad(1, batch_d%first),           &
                     wf%n_v*batch_d%length)
!
         call mem%dealloc(L_adck, wf%n_v, batch_d%length, wf%n_v, wf%n_o)
!
      enddo ! batch_d
!
      call mem%batch_finalize()
!
!     rho_aibj += X_ad * tbar_dibj
!
      call dgemm('N', 'N',                             &
                  wf%n_v,                              &
                  (wf%n_o**2)*(wf%n_v),                &
                  wf%n_v,                              &
                  one,                                 &
                  X_ad,                                & ! X_a,d
                  wf%n_v,                              &
                  tbar_aibj,                           & ! tbar_d,ibj
                  wf%n_v,                              &
                  one,                                 &
                  rho_aibj,                            & ! rho_a,ibj
                  wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
!
!     Term 2: L_dcia * tbar_bjdk * c_ck
!
!     X_aidk = L_dcia * c_ck
!
      call mem%alloc(X_aidk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(X_aidk, (wf%n_v**2)*(wf%n_o**2))
!
      req0 = (wf%eri%n_J)*(wf%n_v)*(wf%n_o)
      req1 = max((wf%eri%n_J)*(wf%n_v) + (wf%n_v**2)*(wf%n_o), 2*(wf%n_v**2)*(wf%n_o))
!
      batch_d = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_d, req0, req1, tag='F_ccsd_e2_2_ccsd 2')
!
      do current_d_batch = 1, batch_d%num_batches
!
         call batch_d%determine_limits(current_d_batch)
!
!        L_iadc = 2 g_iadc - g_icda (ordered as L_aidc)
!
         call mem%alloc(g_iadc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call wf%eri%get_eri_t1('ovvv',g_iadc,                       &
                                1, wf%n_o,                           &
                                1, wf%n_v,                           &
                                batch_d%first, batch_d%get_last(),   &
                                1, wf%n_v)
!
         call mem%alloc(L_aidc, wf%n_v, wf%n_o, batch_d%length, wf%n_v)
         call zero_array(L_aidc, (wf%n_v**2)*wf%n_o*batch_d%length)
         call add_2134_to_1234(two, g_iadc, L_aidc, wf%n_v, wf%n_o, batch_d%length, wf%n_v)
         call add_2431_to_1234(-one, g_iadc, L_aidc, wf%n_v, wf%n_o, batch_d%length, wf%n_v)
!
         call mem%dealloc(g_iadc, wf%n_o, wf%n_v, batch_d%length, wf%n_v)
!
         call dgemm('N',  'N',                         &
                     (wf%n_v)*(wf%n_o)*batch_d%length, &
                     wf%n_o,                           &
                     wf%n_v,                           &
                     one,                              &
                     L_aidc,                           & ! L_aid,c
                     (wf%n_v)*(wf%n_o)*batch_d%length, &
                     c_ai,                             & ! c_c,k
                     wf%n_v,                           &
                     one,                              &
                     X_aidk(1,1,batch_d%first,1),      & ! X_aid,k
                     (wf%n_v**2)*(wf%n_o))
!
         call mem%dealloc(L_aidc, wf%n_v, wf%n_o, batch_d%length, wf%n_v)
!
      enddo ! batch_d
!
      call mem%batch_finalize()
!
!     rho_aibj += X_aidk * tbar_bjdk
!
      call dgemm('N', 'N',                             &
                  (wf%n_v)*(wf%n_o),                   &
                  (wf%n_v)*(wf%n_o),                   &
                  (wf%n_v)*(wf%n_o),                   &
                  one,                                 &
                  X_aidk,                              & ! X_ai,dk
                  (wf%n_v)*(wf%n_o),                   &
                  tbar_aibj,                           & ! tbar_dk,bj
                  (wf%n_v)*(wf%n_o),                   &
                  one,                                 &
                  rho_aibj,                            & ! rho_ai,bj
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_aidk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_e2_2_ccsd
!
!
   module subroutine F_ccsd_f2_2_ccsd(wf, c_aibj, rho_aibj, tbar_aibj)
!!
!!    F transformation f2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_F2,2 = (g_ibkd * tbar_ajcl + g_kbid * tbar_cjal) * c_ckdl
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout)   :: rho_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)      :: tbar_aibj
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: X_ajdk, X_bicl
      real(dp), dimension(:,:,:,:), allocatable :: c_cldk
      real(dp), dimension(:,:,:,:), allocatable :: g_ibkd, g_bidk, g_bidk_2
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ajcl
      real(dp), dimension(:,:,:,:), allocatable :: rho_ajbi
!
!     Reorder c_ckdl to c_cldk
!
      call mem%alloc(c_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(c_aibj, c_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 1: g_ibkd * tbar_ajcl * c_ckdl
!
!     X_ajdk = tbar_ajcl * c_cldk 
!
      call mem%alloc(X_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  tbar_aibj,           & ! tbar_aj,cl
                  (wf%n_v)*(wf%n_o),   &
                  c_cldk,              & ! c_cl,dk
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_ajdk,              &
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(g_ibkd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov',g_ibkd)
!
!     Reorder g_ibkd to g_bidk
!
      call mem%alloc(g_bidk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_2143(g_ibkd, g_bidk, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(g_ibkd, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += X_ajdk * g_bidk
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  X_ajdk,              & ! X_aj,dk
                  (wf%n_v)*(wf%n_o),   &
                  g_bidk,              & ! g_bi,dk
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  rho_ajbi,            &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_ajdk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2:  g_kbid * tbar_cjal * c_ckdl
!
!     Reorder g_bkdi to g_bidk
!
      call mem%alloc(g_bidk_2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(g_bidk, g_bidk_2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_bidk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_bicl = g_bidk * c_dkcl
!
      call mem%alloc(X_bicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  g_bidk_2,            & ! g_bi,dk
                  (wf%n_v)*(wf%n_o),   &
                  c_cldk,              & ! c_dk,cl
                  (wf%n_v)*(wf%n_o),   &
                  zero,                &
                  X_bicl,              &
                  (wf%n_v)*(wf%n_o))
!  
      call mem%dealloc(g_bidk_2, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(c_cldk, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Reorder tbar_alcj to t_ajcl
!
      call mem%alloc(tbar_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(tbar_aibj, tbar_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += tbar_ajcl * X_bicl
!
      call dgemm('N', 'T',             &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  tbar_ajcl,           & ! tbar_aj,cl
                  (wf%n_v)*(wf%n_o),   &
                  X_bicl,              & ! X_bi,cl
                  (wf%n_v)*(wf%n_o),   &
                  one,                 &
                  rho_ajbi,            & ! rho_aj,bi
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(tbar_ajcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_bicl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_f2_2_ccsd
!
!
   module subroutine F_ccsd_g2_2_ccsd(wf, c_abij, rho_abij, tbar_abij)
!!
!!    F transformation g2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_G2,2 = (tbar_bkal * g_jcid + tbar_cjdi * g_lakb) * c_ckdl
!!
!!    Observe that this term is already symmetrized ! 
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)      :: c_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: rho_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)      :: tbar_abij
!
!     Local variables
!
      real(dp), dimension(:,:,:,:), allocatable :: g_abij
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
      real(dp), dimension(:,:,:,:), allocatable :: X_ijkl
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%eri%get_eri_t1('ovov',g_iajb)
!
!     Reorder g_aibj to g_abij
!
      call mem%alloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_2413(g_iajb, g_abij, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Term 1: tbar_akbl * g_icjd * c_ckdl
!
!     X_ijkl = g_icjd * c_ckdl
!
      call mem%alloc(X_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',    &
                  wf%n_o**2,  &
                  wf%n_o**2,  &
                  wf%n_v**2,  &
                  one,        &
                  g_abij,     & ! g_cd,ij
                  wf%n_v**2,  &
                  c_abij,     & ! c_cd,kl 
                  wf%n_v**2,  &
                  zero,       & 
                  X_ijkl,     &
                  wf%n_o**2)
!
!     rho_abij += tbar_abkl * X_ijkl
!
      call dgemm('N', 'T',    &
                  wf%n_v**2,  &
                  wf%n_o**2,  &
                  wf%n_o**2,  &
                  one,        &
                  tbar_abij,  & ! tbar_ab,kl
                  wf%n_v**2,  &
                  X_ijkl,     & ! X_ij,kl
                  wf%n_o**2,  &
                  one,        &
                  rho_abij,   & ! rho_ab,ij
                  wf%n_v**2)
!
      call mem%dealloc(X_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Term 2: tbar_djci * g_kalb * c_ckdl
!
!     X_ijkl = tbar_cidj * c_ckdl
!
      call mem%alloc(X_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',    &
                  wf%n_o**2,  &
                  wf%n_o**2,  &
                  wf%n_v**2,  &
                  one,        &
                  tbar_abij,  & ! tbar_cd,ij
                  wf%n_v**2,  &
                  c_abij,     & ! c_cd,kl 
                  wf%n_v**2,  &
                  zero,       & 
                  X_ijkl,     &
                  wf%n_o**2)
!
!     rho_abij += g_abkl * X_ijkl
!
      call dgemm('N', 'T',    &
                  wf%n_v**2,  &
                  wf%n_o**2,  &
                  wf%n_o**2,  &
                  one,        &
                  g_abij,     & ! g_ab,kl
                  wf%n_v**2,  &
                  X_ijkl,     & ! X_ij,kl
                  wf%n_o**2,  &
                  one,        &
                  rho_abij,   & ! rho_ab,ij
                  wf%n_v**2)
!
      call mem%dealloc(X_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine F_ccsd_g2_2_ccsd
!
!
   module subroutine F_ccsd_h2_2_ccsd(wf, c_aibj, rho_aibj, tbar_aibj)
!!
!!    F transformation h2,2 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_H2,2 = 2 * tbar_aick * L_jbld * c_ckdl 
!!               - tbar_aicl * L_jbkd * c_ckdl
!!               - tbar_aibl * L_kcjd * c_ckdl
!!               - tbar_aidj * L_kclb * c_ckdl
!!               - tbar_ckdj * L_ialb * c_ckdl
!!               - tbar_ckal * L_jbid * c_ckdl
!!               - tbar_ckaj * L_ldib * c_ckdl
!!
!!    Edited by A. K. Schnack-Petersen og E. F. Kjønstad, Sep 2021 
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: tbar_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: rho_aibj
!
!     Local variables     
!
      real(dp), dimension(:,:,:,:), allocatable :: L_bjdl, g_jbld, X_bjck
      real(dp), dimension(:,:,:,:), allocatable :: c_dkcl, X_bjcl
      real(dp), dimension(:,:), allocatable :: X_lj
      real(dp), dimension(:,:), allocatable :: X_bd
      real(dp), dimension(:,:), allocatable :: X_ad
      real(dp), dimension(:,:,:,:), allocatable :: X_ajdl, rho_ajbi
!
!     Construct L_jbld = 2 g_jbld - g_jdlb (ordered as L_bjdl)
!
      call mem%alloc(g_jbld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call mem%alloc(L_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%eri%get_eri_t1('ovov',g_jbld)
!
      call zero_array(L_bjdl, (wf%n_v**2)*(wf%n_o**2))
!
      call add_2143_to_1234(two, g_jbld, L_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_jbld, L_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_jbld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     Term 1: 2 * tbar_aick * L_jbld * c_ckdl 
!
!     X_bjck = L_jbld * c_ckdl
!
      call mem%alloc(X_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',                &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  two,                    &
                  L_bjdl,                 &
                  (wf%n_v)*(wf%n_o),      &
                  c_aibj,                 & ! c_dl,ck
                  (wf%n_v)*(wf%n_o),      &
                  zero,                   &
                  X_bjck,                 &
                  (wf%n_v)*(wf%n_o))
!
!     rho_aibj += tbar_aick * X_bjck
!
      call dgemm('N', 'T',                &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  one,                    &
                  tbar_aibj,              & ! tbar_ai,ck
                  (wf%n_v)*(wf%n_o),      &
                  X_bjck,                 &
                  (wf%n_v)*(wf%n_o),      &
                  one,                    &
                  rho_aibj,               &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_bjck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 2: - tbar_aicl * L_jbkd * c_ckdl
!
      call mem%alloc(c_dkcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call sort_1234_to_1432(c_aibj, c_dkcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     X_bjcl = - L_jbkd * c_ckdl
!
      call mem%alloc(X_bjcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',                &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  -one,                   &
                  L_bjdl,                 & ! L_bj,dk
                  (wf%n_v)*(wf%n_o),      &
                  c_dkcl,                 &
                  (wf%n_v)*(wf%n_o),      &
                  zero,                   &
                  X_bjcl,                 &
                  (wf%n_v)*(wf%n_o))

!
      call mem%dealloc(c_dkcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += tbar_aicl * X_bjcl
!
      call dgemm('N', 'T',                &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  one,                    &
                  tbar_aibj,              & ! tbar_ai,cl
                  (wf%n_v)*(wf%n_o),      &
                  X_bjcl,                 &
                  (wf%n_v)*(wf%n_o),      &
                  one,                    &
                  rho_aibj,               &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(X_bjcl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Term 3: - tbar_aibl * L_kcjd * c_ckdl
!
!     X_lj = - c_ckdl * L_kcjd
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_v**2)*(wf%n_o),   &
                  -one,                   &
                  c_aibj,                 & ! c_ckd,l
                  (wf%n_v**2)*(wf%n_o),   &
                  L_bjdl,                 & ! L_ckd,j
                  (wf%n_v**2)*(wf%n_o),   &
                  zero,                   &
                  X_lj,                   &
                  wf%n_o)
!
!     rho_aibj += tbar_aibl * X_lj
!
      call dgemm('N', 'N',                &
                  (wf%n_v**2)*(wf%n_o),   &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  one,                    &
                  tbar_aibj,              & ! tbar_aib,l
                  (wf%n_v**2)*(wf%n_o),   &
                  X_lj,                   &
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               &
                  (wf%n_v**2)*(wf%n_o))
!
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
!     Term 4: - tbar_aidj * L_kclb * c_ckdl
!
!     X_bd = L_kclb * c_ckdl
!
      call mem%alloc(X_bd, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',                &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  -one,                   &
                  L_bjdl,                 & ! L_b,lck
                  wf%n_v,                 &
                  c_aibj,                 & ! c_d,lck
                  wf%n_v,                 &
                  zero,                   &
                  X_bd,                   &
                  wf%n_v)
!
!     rho_aibj += rho_bjai = X_bd * tbar_aidj (we symmetrize all later)
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  wf%n_v,                 &
                  one,                    &
                  X_bd,                   &
                  wf%n_v,                 &
                  tbar_aibj,              & ! tbar_d,jai
                  wf%n_v,                 &
                  one,                    &
                  rho_aibj,               & ! rho_b,jai
                  wf%n_v)
!
      call mem%dealloc(X_bd, wf%n_v, wf%n_v)
!
!     Term 5: - tbar_ckdj * L_ialb * c_ckdl
!
!     X_lj = - c_ckdl * tbar_ckdj
!
      call mem%alloc(X_lj, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  (wf%n_v**2)*wf%n_o,     &
                  -one,                   &
                  c_aibj,                 & ! c_ckd,l
                  (wf%n_v**2)*wf%n_o,     &
                  tbar_aibj,              & ! tbar_ckd,j
                  (wf%n_v**2)*wf%n_o,     &
                  zero,                   &
                  X_lj,                   &
                  wf%n_o)
!
!     rho_aibj += L_aibl * X_lj
!
      call dgemm('N', 'N',                &
                  (wf%n_v**2)*(wf%n_o),   &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  one,                    &
                  L_bjdl,                 & ! L_aib,l
                  (wf%n_v**2)*(wf%n_o),   &
                  X_lj,                   &
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               &
                  (wf%n_v**2)*(wf%n_o))
!
      call mem%dealloc(X_lj, wf%n_o, wf%n_o)
!
!     Term 6: - tbar_ckal * L_jbid * c_ckdl
!
!     X_ad = -tbar_ckal * c_ckdl
!
      call mem%alloc(X_ad, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',                &
                  wf%n_v,                 &
                  wf%n_v,                 &
                  (wf%n_o**2)*wf%n_v,     &
                  -one,                   &
                  tbar_aibj,              & ! tbar_a,lck
                  wf%n_v,                 &
                  c_aibj,                 & ! c_d,lck
                  wf%n_v,                 &
                  zero,                   &
                  X_ad,                   &
                  wf%n_v)  
!
!     rho_aibj += X_ad * L_dibj
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  (wf%n_o**2)*(wf%n_v),   &
                  wf%n_v,                 &
                  one,                    &
                  X_ad,                   &
                  wf%n_v,                 &
                  L_bjdl,                 & ! L_d,ibj
                  wf%n_v,                 &
                  one,                    &
                  rho_aibj,               &
                  wf%n_v)
!
      call mem%dealloc(X_ad, wf%n_v, wf%n_v)
!
!     Term 7: - tbar_ckaj * L_ldib * c_ckdl
!
!     X_ajdl = - tbar_ckaj * c_ckdl
!
      call mem%alloc(X_ajdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call dgemm('N', 'N',                &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  -one,                   &
                  tbar_aibj,              & ! t_aj,ck
                  (wf%n_v)*(wf%n_o),      &
                  c_aibj,                 & ! c_ck,dl
                  (wf%n_v)*(wf%n_o),      &
                  zero,                   &
                  X_ajdl,                 &
                  (wf%n_v)*(wf%n_o))
!
      call mem%alloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_aibj += X_ajdl * L_dlbi
!
      call dgemm('N', 'N',                &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  (wf%n_v)*(wf%n_o),      &
                  one,                    &
                  X_ajdl,                 &
                  (wf%n_v)*(wf%n_o),      &
                  L_bjdl,                 & ! L_dl,bi
                  (wf%n_v)*(wf%n_o),      &
                  zero,                   &
                  rho_ajbi,               &
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_bjdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(X_ajdl, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call add_1432_to_1234(one, rho_ajbi, rho_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(rho_ajbi, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine F_ccsd_h2_2_ccsd
!
!
end submodule F_ccsd
