submodule (cc2_class) jacobian_cc2
!
!!
!!    Jacobian submodule (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    and Andreas Skeidsvoll, 2018
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
   module subroutine jacobian_transform_trial_vector_cc2(wf, c_i)
!!
!!    Jacobian transform trial vector 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Sep 2018
!!
      class(cc2), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1) :: c_i
!
      call wf%jacobian_cc2_transformation(c_i)
!
   end subroutine jacobian_transform_trial_vector_cc2
!
!
   module subroutine jacobian_cc2_transformation_cc2(wf, c)
!!
!!    Jacobian transformation (cc2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
!!
!!    Directs the transformation by the cc2 Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >,
!!
!!    where the basis employed for the brackets is biorthonormal.
!!    The transformation is rho = A c, i.e.,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck
!!                  + 1/2 sum_ckdl A_mu,ckdl c_ckdl (1 + delta_ck,dl).
!!
!!    On exit, c is overwritten by rho. That is, c_ai = rho_ai,
!!    and c_aibj = rho_aibj.
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1)   :: c
!
      real(dp), dimension(:,:), allocatable     :: c_ai
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj
!
      real(dp), dimension(:,:), allocatable     :: rho_ai   
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj
!
!
      integer(i15) :: i, j, a, b, ai, bj, aibj ! Index
!
!     Allocate and zero the transformed vector (singles part)
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      rho_ai = zero
!
      call mem%alloc(c_ai, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, ai) 
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
!     :: CCS contributions to the singles c vector ::
!
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     :: CC2 contributions to the transformed singles vector ::
!
      call wf%jacobian_cc2_a1(rho_ai, c_ai)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!  
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c_aibj(a, i, b, j) = c(wf%n_o*wf%n_v + aibj, 1)
                     c_aibj(b, j, a, i) = c(wf%n_o*wf%n_v + aibj, 1)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Scale the doubles vector by 1 + delta_ai,bj, i.e.
!     redefine to c_ckdl = c_ckdl (1 + delta_ck,dl)
!
!$omp parallel do schedule(static) private(ai) 
      do a =1, wf%n_v
         do i = 1, wf%n_o
!
            c_aibj(a, i,a, i) = two*c_aibj(a, i, a, i)
!
         enddo
      enddo
!$omp end parallel do
!
      call wf%jacobian_cc2_b1(rho_ai, c_aibj)
      call wf%jacobian_cc2_c1(rho_ai, c_aibj)
!
!     Done with singles vector c; overwrite it with
!     transformed vector for exit
!
!$omp parallel do schedule(static) private(a, i, ai) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            c(ai, 1) = rho_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do

      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
!     :: CC2 contributions to the transformed doubles vector ::
!
!     Allocate unpacked transformed vector
!
      call mem%alloc(rho_aibj, (wf%n_o), (wf%n_v), (wf%n_o), (wf%n_v))
      rho_aibj = zero
!
!     Contributions from singles vector c
!
      call wf%jacobian_cc2_a2(rho_aibj, c_ai)
!
      call mem%dealloc(c_ai, wf%n_v, wf%n_o)
!
!     Contributions from doubles vector c      
!
      call wf%jacobian_cc2_b2(rho_aibj, c_aibj)
!
      call mem%dealloc(c_aibj, (wf%n_v), (wf%n_o), (wf%n_v), (wf%n_o))
!
!     Overwrite the incoming doubles c vector & pack in
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!  
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  if (ai .ge. bj) then
!
                     aibj = ai*(ai-3)/2 + ai + bj
!
                     c((wf%n_o)*(wf%n_v) + aibj, 1) = rho_aibj(a, i, b, j)
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine jacobian_cc2_transformation_cc2
!
!
   module subroutine jacobian_cc2_a1_cc2(wf, rho_ai, c_ai)
!!
!!    Jacobian CC2 A1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_ai^A1 = sum_bjck (L_kcjb u_aick c_bj - g_kcjb (u_ckbi c_aj + u_ckaj c_bi))
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)   :: rho_ai   
!
   end subroutine jacobian_cc2_a1_cc2
!
!
   module subroutine jacobian_cc2_b1_cc2(wf, rho_ai, c_aibj)
!!
!!    Jacobian CC2 B1
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_ai^B1 = 2 sum_bj F_jb c_aibj - F_jb c_ajbi ) - sum_kjb L_kijb c_akbj + sum_bkc L_abkc c_bick
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                  :: rho_ai   
!
   end subroutine jacobian_cc2_b1_cc2
!
!
   module subroutine jacobian_cc2_a2_cc2(wf, rho_aibj, c_ai)
!!
!!    Jacobian CC2 A2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^A2 = /(1/Δ_aibj)P_aibj sum_c g_abkc c_cj - sum_k g_aikj c_bk, 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)  :: rho_aibj   
!
   end subroutine jacobian_cc2_a2_cc2
!
!
   module subroutine jacobian_cc2_b2_cc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai D. Folkestad Eirik F. Kjønstad Jan 2019
!!
!!    rho_aibj^B2 = ε_aibj c_aibj/(1/Δ_aibj) 
!!
      implicit none
!
      class(cc2) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)   :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(out)  :: rho_aibj   
!
   end subroutine jacobian_cc2_b2_cc2
!
!
end submodule jacobian_cc2
