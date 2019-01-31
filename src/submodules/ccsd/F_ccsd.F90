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
   end subroutine F_ccsd_a2_1_ccsd
!
!
end submodule F_ccsd
