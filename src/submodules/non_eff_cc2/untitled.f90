submodule (non_eff_cc2_class) omega_cc2
!
!!
!!    Omega submodule (CC2)
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto and Alexander Paul, 2018
!!
!!    Routines to construct 
!!
!!    Ω =  < mu | exp(-T) H exp(T) | R >
!!
!
   implicit none
!
!
contains
!
   module subroutine construct_omega_cc2(wf, omega)
!!
!!    Construct omega 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Direqts the construction of the projection vector < mu | exp(-T) H exp(T) | R >
!!    for the current wavefunction amplitudes.
!!
      implicit none
!
      class(non_eff_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: u_aibj
!
      omega = zero
!
      call wf%omega_ccs_a1(omega)
!
      call mem%alloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%construct_u(u_aibj)
!
      call wf%omega_cc2_a1(omega, u_aibj, omega)
      call wf%omega_cc2_b1(omega, u_aibj, omega)
      call wf%omega_cc2_c1(omega, u_aibj, omega)
!
      call mem%dealloc(u_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine construct_omega_cc2
!
!
   module subroutine omega_cc2_a1_cc2(wf, u_bjci, omega)
!!
!!    Omega CC2 A1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Calculates the A1 term,
!!
!!       A1: sum_ckd u_bj_ci * g_abjc,
!!
!!    with 
!!       
!!       u_bj_ci = 2*t_bj_ci - t_bi_cj
!!
!!    and adds it to the projection vector omega
!!
      implicit none
!
      class(non_eff_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_bjci
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout)          :: omega
!
   end subroutine omega_cc2_a1_cc2
!
!
   module subroutine omega_cc2_b1_cc2(wf, u_ajbk, omega)
!!
!!    Omega CC2 B1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, 
!!    Linda Goletto, and Alexander Paul, Dec 2018
!!
!!    Calculates the B1 term,
!!
!!       B1: - sum_ckl g_kb,ji * u_aj,bk,
!!
!!    with
!!
!!      u_aj_bk = 2t_aj,bk - t_ak,bj
!!
!!    and adds it to the projection vector (omega) of
!!    the wavefunction object wf.
!!
      implicit none
!
      class(non_eff_cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_ajbk
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout)          :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kbji
      real(dp), dimension(:,:,:,:), allocatable :: g_jbki
!
!     g_kbji ordered as g_jbki
!
      call mem%alloc(g_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%get_ovoo(g_kbji)
!
      call mem%alloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call sort_1234_to_3214(g_kbji, g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_kbji, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
!     omega_ai += - sum_ckl g_kb,ji * u_aj,bk
!
      call dgemm('N', 'N',            &
                  wf%n_v,             &
                  wf%n_o,             &
                  (wf%n_o**2)*wf%n_v, &
                  -one,               &
                  u_ajbk,             & ! u_a_jbk
                  wf%n_v,             &
                  g_jbki,             & ! g_jbk_i
                  (wf%n_o**2)*wf%n_v, &
                  one,                &
                  omega,              & ! omega_a_i
                  wf%n_v)
!
      call mem%dealloc(g_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine omega_cc2_b1_cc2
!
!
   module subroutine omega_cc2_c1_cc2(wf, u_aibj, omega)
!!
!!    Omega CC2 C1 term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,  Jan 2019
!!
!!    Calculates the C1 term,
!!
!!       C1: - sum_bj u_ai,bj * F_jb,
!!
!!    with 
!!       
!!       u_ai_bj = 2*t_ai_bj - t_aj_bi
!!
!!
      implicit none
!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: u_aibj
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout)          :: omega
!
    end subroutine omega_cc2_c1_cc2
!
!
    module subroutine construct_u_cc2(wf, u)
!!
!!    Construct U 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    Construct
!!
!!       u_aibj = 2t_aibj - t_ajbi
!!
!!    with
!!
!!       t_aibj = - g_aibj/ε_aibj
!!
!!    where
!!
!!       ε_aibj = ε_a - ε_i + ε_b - ε_j 
!!
!!    and ε_r is the r'th orbital energy.
!!
      class(cc2), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(inout) :: u
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)  
!
      call wf%get_vovo(g_aibj)
!
      do b = 1, wf%n_v 
         do j = 1, wf%n_o 
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  u_aibj(a, i, b, j) = (two* g_aibj(a, i, b, j) - g_aibj(a, j, i, b))/ &
                                       (wf%fock_diagonal(i, 1) + fock_diagonal(j, 1) &
                                        - wf%fock_diagonal(wf%n_o + a, 1) &
                                        - wf%fock_diagonal(wf%n_o + b, 1) )
!
               enddo
            enddo
         enddo 
      enddo
!    
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
    end subroutine construct_u_cc2
!
!
end submodule omega_cc2
