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
!     ! Explicit reordering of c_b_j
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
      call sort_1234_to_1432(g_ab_ji, g_ai_jb, (wf%n_v)*(wf%n_v), (wf%n_o)*(wf%n_o), &
                                              (wf%n_v)*(wf%n_o), (wf%n_o)*(wf%n_v))
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

end submodule jacobian
