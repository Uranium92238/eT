submodule (ccs_class) F_ccs
!
!!
!!    F-transformation submodule (ccs)
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
   module subroutine F_transform_vector_ccs(wf, c)
!!
!!    F transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018     
!!
!!    Directs the transformation by the F matrix.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: c
!
      real(dp), dimension(:,:), allocatable :: rho
!
      integer(i15) :: a, i, ai
!
      call mem%alloc(rho, wf%n_v, wf%n_o)
      rho = zero
!
      call wf%F_ccs_a1_0(c, rho)
!  
      call wf%F_ccs_a1_1(c, rho)
      call wf%F_ccs_b1_1(c, rho)
      call wf%F_ccs_c1_1(c, rho)
!
!$omp parallel do private(a, i, ai)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai, 1) = rho(a, i)
!
         enddo
      enddo
!$omp end parallel do      
!
      call mem%dealloc(rho, wf%n_v, wf%n_o)
!
   end subroutine F_transform_vector_ccs
!
!
   module subroutine F_ccs_a1_0_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,0 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,0_ai = 2 L_iajb c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
   end subroutine F_ccs_a1_0_ccs
!
!
   module subroutine F_ccs_a1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation A1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_A1,1_ai = - (F_ib tbar_aj + F_ja tbar_bi) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!
   end subroutine F_ccs_a1_1_ccs
!
!
   module subroutine F_ccs_b1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation B1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_B1,1_ai = - (L_ikjb tbar_ak + L_jkia tbar_bk) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!
   end subroutine F_ccs_b1_1_ccs
!
!
   module subroutine F_ccs_c1_1_ccs(wf, c_ai, rho_ai)
!!
!!    F transformation C1,1 term
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Feb 2018
!!
!!    rho_C1,1_ai = - (L_cajb tbar_ci + L_cbia tbar_cj) c_bj
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho_ai
!
!
   end subroutine F_ccs_c1_1_ccs
!
end submodule F_ccs
