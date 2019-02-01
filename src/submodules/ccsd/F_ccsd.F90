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
   end subroutine F_ccsd_a2_1_ccsd
!
!
   module subroutine F_ccsd_a1_2_ccsd(wf, c_ai, rho_ai)
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)   :: rho_ai

!
   end subroutine F_ccsd_a1_2_ccsd
!
!
   module subroutine F_ccsd_b1_2_ccsd(wf, c_ai, rho_ai)
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)   :: rho_ai
!
   end subroutine F_ccsd_b1_2_ccsd
!
!
   module subroutine F_ccsd_c1_2_ccsd(wf, c_ai, rho_ai)
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
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)      :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)   :: rho_ai
!
   end subroutine F_ccsd_c1_2_ccsd
!
!
   module subroutine F_ccsd_d1_2_ccsd(wf, c_aibj, rho_ai)
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
!
   end subroutine F_ccsd_d1_2_ccsd
!
!
   module subroutine F_ccsd_e1_2_ccsd(wf, c_aibj, rho_ai)
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
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                :: rho_ai
!
   end subroutine F_ccsd_e1_2_ccsd
!
end submodule F_ccsd
