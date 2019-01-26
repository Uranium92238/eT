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
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)              :: omega
!
      type(batching_index) :: batch_a 
!
      integer(i15) :: req0, req1
!
      call batch_a%init(wf%n_v)
!
      req0 = (wf%n_o)*(wf%n_v)*(wf%integrals%n_J)
      req1 = (wf%n_v)**2*(wf%n_o) + (wf%n_v)*(wf%integrals%n_J)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%integrals%get_vvov(g_abjc,                      &
                                    batch_a%first, batch_a%last, &
                                    1, wf%n_v,                   &
                                    1, wf%n_o,                   &
                                    1, wf%n_v)
!
         call mem%alloc(omega_ai, batch_a%length, wf%n_o)
!
         call dgemm('N','N',               &
                     batch_a%length,       &
                     wf%n_o,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     one,                  &
                     g_abjc,               &
                     batch_a%length,       &
                     u_bjci,               &
                     (wf%n_o)*(wf%n_v)**2, &
                     zero,                 &
                     omega_ai,             &
                     batch_a%length)
!
!$omp parallel do private(i, a)
         do i = 1, wf%n_o
            do a = 1, batch_a%length
!  
               omega(a + batch_a%first - 1, i) = omega(a + batch_a%first - 1, i) 
                                                & + omega_ai(a, i)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(omega_ai, batch_a%length, wf%n_o)
!
      enddo 
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
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad, Jan 2019
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
      real(dp), dimension(:,:), allocatable :: F_bj
!
      call mem%alloc(F_bj, wf%n_o, wf%n_v)
      call sort_12_to_21(wf%fock_ia, F_bj, wf%n_o, wf%n_v)
!
      call dgemm('N','N',            &
                  (wf%n_o)*(wf%n_v), &
                  1,                 &
                  (wf%n_o)*(wf%n_v), &
                  -one,              &
                  u_aibj,            &
                  (wf%n_o)*(wf%n_v), &
                  F_bj,              &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  omega,             &
                  (wf%n_o)*(wf%n_v))
!
      call mem%dealloc(F_bj, wf%n_o, wf%n_v)
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
!$omp parallel do private(b, j, i, a)
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
!$omp end parallel do
!    
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)      
!
    end subroutine construct_u_cc2
!
!
   module subroutine calculate_energy_cc2(wf)
!!
!!    Calculate energy 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    E = E_HF + sum_aibj (t_i^a*t_j^b + t_ij^ab) L_iajb
!!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, g_iajb 
!
      real(dp) :: correlation_energy
!
      integer(i15) :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%integrals%get_vovo(g_aibj)
      call wf%integrals%get_ovov(g_iajb)
!
      correlation_energy = zero 
!
!$omp parallel do private(a,i,b,j) reduction(+:correlation_energy)
      do b = 1, wf%n_v
         do i = 1, wf%n_o 
            do j = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  correlation_energy = correlation_energy +                                &
                                       (wf%t1(a, i)*wf%t1(b, j) -                          &
                                       (g_aibj(a,i,b,j))/(wf%fock_diagonal(wf%n_o + a, 1)  &
                                                         + wf%fock_diagonal(wf%n_o + b, 1) &
                                                         - wf%fock_diagonal(i,1)           &
                                                         - wf%fock_diagonal(j,1)))         &
                                       *(two*g_iajb(i,a,j,b)-g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_cc2
!
!
end submodule omega_cc2
