!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
submodule (mlcc2_class) jacobian_mlcc2
!
!!
!!    Jacobian submodule 
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!   
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | ν >.
!!  
!!    NOTE: All routines adapted from the CC2 routines written by Sarai D. Folkestad 
!!    and Eirik F. Kjønstad
!!
! 
   implicit none
!
!
contains
!
!
   module subroutine prepare_for_jacobian_mlcc2(wf)
!!
!!    Prepare for Jacobian
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
      implicit none 
!
      class(mlcc2), intent(inout) :: wf 
!
      type(timings), allocatable :: timer
!
      timer = timings('Prepare for Jacobian transformation MLCC2', pl='normal')
      call timer%turn_on()
!
      call wf%save_jacobian_a1_intermediates(wf%n_cc2_o, wf%n_cc2_v, 1, 1)
!
      call timer%turn_off()
!
   end subroutine prepare_for_jacobian_mlcc2
!
!
   module subroutine jacobian_transformation_mlcc2(wf, c, rho)
!!
!!    Jacobian transformation (mlcc2)
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
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
      use array_utilities, only: scale_diagonal
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)  :: c
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: rho
!
      real(dp), dimension(:,:,:,:), allocatable :: c_aibj
      real(dp), dimension(:,:,:,:), allocatable :: rho_aibj
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian transformation MLCC2', pl='normal')
      call timer%turn_on()
!
!     Zero the transformed vector
!
      call zero_array(rho, wf%n_es_amplitudes)
!
!     CCS contributions to the singles rho vector
!
      call wf%ccs%jacobian_transformation(c(1 : wf%n_t1), rho(1 : wf%n_t1))
!
!     CC2 contributions to the transformed singles vector  
!
      call wf%jacobian_cc2_a1(rho(1 : wf%n_t1), c(1 : wf%n_t1), wf%n_cc2_o, wf%n_cc2_v, &
                              1, 1)
!
!     Allocate the incoming unpacked doubles vector
!
      call mem%alloc(c_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      call squareup(c(wf%n_t1 + 1 : wf%n_es_amplitudes), &
                    c_aibj, wf%n_cc2_v*wf%n_cc2_o)
!
!     Scale the doubles vector by (1 + delta_ai,bj)
!
      call scale_diagonal(two, c_aibj, wf%n_cc2_o*wf%n_cc2_v)
!     
      call wf%jacobian_cc2_b1(rho(1 : wf%n_t1), c_aibj, wf%n_cc2_o, wf%n_cc2_v, &
                              1, 1, wf%n_cc2_o, wf%n_cc2_v)
!
!     CC2 contributions to the transformed doubles vector
!
      call mem%alloc(rho_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
      call zero_array(rho_aibj, (wf%n_cc2_o**2)*(wf%n_cc2_v**2))
!
!     Contributions from singles vector c
!
      call wf%jacobian_cc2_a2(rho_aibj, c(1 : wf%n_t1), wf%n_cc2_o, wf%n_cc2_v, &
                              1, 1, wf%n_cc2_o, wf%n_cc2_v)
! 
!     Biorthonormalize
!
      call scale_diagonal(half, rho_aibj, wf%n_cc2_o*wf%n_cc2_v)
!
!     Symmetrize (b2 term is already symmetric)
!
      call symmetric_sum(rho_aibj, wf%n_cc2_o*wf%n_cc2_v)
!
      call wf%jacobian_cc2_b2(rho_aibj, c_aibj)
!
      call mem%dealloc(c_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
!     Pack in
!
      call packin(rho(wf%n_t1 + 1:), rho_aibj, wf%n_cc2_v*wf%n_cc2_o)
!
      call mem%dealloc(rho_aibj, wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_transformation_mlcc2
!
!
   module subroutine jacobian_cc2_a1_mlcc2(wf, rho_ai, c_ai, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v)
!!
!!    Jacobian CC2 A1
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Calculates the A1 term
!!
!!       A1:   sum_bjck L_kcjb u_aick c_bj 
!!           - sum_bjck g_jbkc u_ckbi c_aj 
!!           - sum_bjck g_kcjb u_ckaj c_bi
!!
!!    and adds it to rho_ai.
!!
!!    Index restrictions:
!!
!!       Term 1:
!!       
!!          a, i, c, k : CC2 orbitals
!!
!!          b, j : unrestricted
!!
!!       Term 2:
!!       
!!          b, i, c, k : CC2 orbitals
!!
!!          a, j : unrestricted
!!
!!       Term 3:
!!       
!!          a, j, c, k : CC2 orbitals
!!
!!          b, i : unrestricted
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)   :: rho_ai
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
      real(dp), dimension(:,:), allocatable     :: X_ck, rho_ai_copy, X_kc
      real(dp), dimension(:,:), allocatable     :: Y_ab, Y_ji, c_jb
      real(dp), dimension(:,:,:), allocatable   :: X_Jkj, X_Jjk, L_Jov, L_Jkc
      real(dp), dimension(:), allocatable       :: X_J 
!
      integer :: a, i, c, k, J
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCC2 A1 transformation', pl='verbose')
      call timer%turn_on()
!
!     Term 1: sum_bjck L_kcjb u_aick c_bj 
!       
!        a, i, c, k : CC2 orbitals
!
!        b, j : unrestricted
!
!     2 sum_bjck g_kcjb u_aick c_bj = L_Jkc L_Jjb c_bj u_aick      
!
      call mem%alloc(L_Jov, wf%eri%n_J, wf%n_o, wf%n_v)
      call wf%eri%get_cholesky_t1(L_Jov, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
!
      call mem%alloc(c_jb, wf%n_o, wf%n_v)
      call sort_12_to_21(c_ai, c_jb, wf%n_v, wf%n_o)
!
      call mem%alloc(X_J, wf%eri%n_J)
      call dgemv('N',            &
                  wf%eri%n_J,    &
                  wf%n_v*wf%n_o, &
                  one,           &
                  L_Jov,         &
                  wf%eri%n_J,    &
                  c_jb,          &
                  1,             &
                  zero,          &
                  X_J,           &
                  1)
!
      call mem%dealloc(c_jb, wf%n_o, wf%n_v)
!
      call mem%alloc(L_Jkc, wf%eri%n_J, n_cc2_o, n_cc2_v)
!
!$omp parallel do private (c, k, J)
      do c = 1, n_cc2_v
         do k = 1, n_cc2_o
            do J = 1, wf%eri%n_J
!
               L_Jkc(J, k, c) = L_Jov(J, k + first_cc2_o - 1, c + first_cc2_v - 1)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(X_kc, n_cc2_o, n_cc2_v)
!
      call dgemv('T',              & 
                  wf%eri%n_J,      &
                  n_cc2_v*n_cc2_o, &                  
                  two,             &
                  L_Jkc,           &
                  wf%eri%n_J,      &
                  X_J,             &
                  1,               &
                  zero,            &
                  X_kc,            &
                  1)
!
      call mem%dealloc(X_J, wf%eri%n_J)
      call mem%dealloc(L_Jkc, wf%eri%n_J, n_cc2_o, n_cc2_v)
!
      call mem%alloc(X_ck, n_cc2_v, n_cc2_o)
      call sort_12_to_21(X_kc, X_ck, n_cc2_o, n_cc2_v)
      call mem%dealloc(X_kc, n_cc2_o, n_cc2_v)
!
!     - sum_bjck g_kbjc u_aick c_bj = L_Jkb L_Jjc c_bj u_aick      
!
!     Note: L_J_jc_t1 = L_J_jc = L_J_cj
!
      call mem%alloc(X_Jkj, wf%eri%n_J, n_cc2_o, wf%n_o)
!
      call dgemm('N', 'N',            &
                  wf%eri%n_J*n_cc2_o, &
                  wf%n_o,             &
                  wf%n_v,             &
                  one,                &
                  L_Jov,              &
                  wf%eri%n_J*wf%n_o,  &
                  c_ai,               &
                  wf%n_v,             &
                  zero,               &
                  X_Jkj,              &
                  wf%eri%n_J*n_cc2_o)
!
      call mem%alloc(X_Jjk, wf%eri%n_J, wf%n_o, n_cc2_o)
      call sort_123_to_132(X_Jkj, X_Jjk, wf%eri%n_J, n_cc2_o, wf%n_o)
      call mem%dealloc(X_Jkj, wf%eri%n_J, n_cc2_o, wf%n_o)
!
      call dgemm('T', 'N',                  &
                  n_cc2_v,                  &
                  n_cc2_o,                  &
                  wf%eri%n_J*wf%n_o,        &
                  -one,                     &
                  L_Jov(1, 1, first_cc2_v), &
                  wf%eri%n_J*wf%n_o,        &
                  X_Jjk,                    &
                  wf%eri%n_J*wf%n_o,        &
                  one,                      &
                  X_ck,                     &
                  n_cc2_v)
!
      call mem%dealloc(X_Jjk, wf%eri%n_J, wf%n_o, n_cc2_o)
      call mem%dealloc(L_Jov, wf%eri%n_J, wf%n_o, wf%n_v)
      call mem%alloc(rho_ai_copy, n_cc2_v, n_cc2_o)
!
      call dgemm('N', 'N',          &
                  n_cc2_o*n_cc2_v,  &
                  1,                &
                  n_cc2_o*n_cc2_v,  &
                  one,              &
                  wf%u_aibj,        & ! u_ai_ck
                  n_cc2_o*n_cc2_v,  &
                  X_ck,             &
                  n_cc2_o*n_cc2_v,  &
                  zero,             &
                  rho_ai_copy,      &
                  n_cc2_o*n_cc2_v)
!
      call mem%dealloc(X_ck, n_cc2_v, n_cc2_o)
!
!$omp parallel do private(a, i) collapse(2)
      do i = 1, n_cc2_o
         do a = 1, n_cc2_v
!
            rho_ai(a + first_cc2_v - 1, i + first_cc2_o - 1) = &
            rho_ai(a + first_cc2_v - 1, i + first_cc2_o - 1) + rho_ai_copy(a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_ai_copy, n_cc2_v, n_cc2_o)
!
!     Term 2: - sum_bjck g_jbkc u_ckbi c_aj = - sum_j c_aj Y_ji 
!      
!        b, i, c, k : CC2 orbitals
!
!        a, j : unrestricted
!
      call mem%alloc(Y_ji, wf%n_o, n_cc2_o)
!
      call wf%jacobian_a1_intermediate_oo%open_('read', 'rewind')
!
      call wf%jacobian_a1_intermediate_oo%read_(Y_ji, n_cc2_o*wf%n_o)
!
      call wf%jacobian_a1_intermediate_oo%close_('keep')
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  n_cc2_o,                &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai,                   & 
                  wf%n_v,                 &
                  Y_ji,                   & 
                  wf%n_o,                 &
                  one,                    &
                  rho_ai(1, first_cc2_o), &
                  wf%n_v)
!
      call mem%dealloc(Y_ji, wf%n_o, n_cc2_o)
!
!     Term 3:  - sum_bjck g_kcjb u_ajck c_bi = - sum_b Y_ab c_bi 
!
!        a, j, c, k : CC2 orbitals
!
!        b, i : unrestricted
!
      call mem%alloc(Y_ab, n_cc2_v, wf%n_v)
!
      call wf%jacobian_a1_intermediate_vv%open_('read', 'rewind')
!
      call wf%jacobian_a1_intermediate_vv%read_(Y_ab, n_cc2_v*wf%n_v)
!
      call wf%jacobian_a1_intermediate_vv%close_('keep')
!
      call dgemm('N', 'N', &
                  n_cc2_v, &
                  wf%n_o,  &
                  wf%n_v,  &
                  -one,    &
                  Y_ab,    &
                  n_cc2_v, &
                  c_ai,    &
                  wf%n_v,  &
                  one,     &
                  rho_ai(first_cc2_v, 1), &
                  wf%n_v)
!
      call mem%dealloc(Y_ab, n_cc2_v, wf%n_v)

      call timer%turn_off()
!
   end subroutine jacobian_cc2_a1_mlcc2
!
!
   module subroutine jacobian_cc2_a2_mlcc2(wf, rho_aibj, c_ai, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Jacobian CC2 A2
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Calculates the A2 term 
!!
!!       A2: sum_c g_aibc c_cj - sum_k g_aikj c_bk, 
!!
!!    and adds it to rho_aibj.
!!
!!    The first term is calculated in batches over c.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
!!       Term 1: 
!!
!!          c : unrestricted
!!
!!       Term 2: 
!!
!!          k : unrestricted
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                       :: c_ai
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(out)  :: rho_aibj   
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kjai 
      real(dp), dimension(:,:,:,:), allocatable :: g_aibc
!
      type(batching_index) :: batch_c 
      integer              :: req0, req1, current_c_batch 
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCC2 A2 transformation', pl='verbose')
      call timer%turn_on()
!
!     Term 2: - sum_k c_bk g_kjai
!
!        a, i, b, j : CC2 orbitals
!
!        k : unrestricted
!
      call mem%alloc(g_kjai, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call wf%eri%get_eri_t1('oovo', g_kjai, 1, wf%n_o, first_cc2_o, last_cc2_o, &
                                             first_cc2_v, last_cc2_v, first_cc2_o, last_cc2_o)
!
      call dgemm('N','N',                 &
                  n_cc2_v,                &
                  (n_cc2_v)*(n_cc2_o)**2, &
                  wf%n_o,                 &
                  -one,                   &
                  c_ai,                   & ! c_b_k
                  wf%n_v,                 &
                  g_kjai,                 & ! g_k_jai
                  wf%n_o,                 &
                  one,                    &
                  rho_aibj,               & ! rho_b_jai but we will symmetrize later
                  n_cc2_v)
!
      call mem%dealloc(g_kjai, wf%n_o, n_cc2_o, n_cc2_v, n_cc2_o)
!
!     Term 1: sum_c g_aibc c_cj
!
!        a, i, b, j : CC2 orbitals
!
!        c : unrestricted
!
      req0 = (wf%eri%n_J)*(n_cc2_o)*(n_cc2_v)
      req1 = (n_cc2_o)*(n_cc2_v)**2 + (wf%eri%n_J)*(n_cc2_v)
!
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_c, req0, req1)
!
      do current_c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(current_c_batch)
!
         call mem%alloc(g_aibc, n_cc2_v, n_cc2_o, n_cc2_v, batch_c%length)
!
         call wf%eri%get_eri_t1('vovv', g_aibc, first_cc2_v, last_cc2_v, &
                                                first_cc2_o, last_cc2_o, &
                                                first_cc2_v, last_cc2_v, &
                                                batch_c%first, batch_c%last)
!
         call dgemm('N','N',                             &
                     (n_cc2_o)*(n_cc2_v)**2,             &
                     n_cc2_o,                            &
                     batch_c%length,                     &
                     one,                                &
                     g_aibc,                             &! g_aib_c 
                     (n_cc2_o)*(n_cc2_v)**2,             &
                     c_ai(batch_c%first, first_cc2_o),   & ! c_c_j
                     wf%n_v,                             &
                     one,                                &
                     rho_aibj,                           & ! rho_aib_j 
                     (n_cc2_o)*(n_cc2_v)**2)
!
         call mem%dealloc(g_aibc, n_cc2_v, n_cc2_o, n_cc2_v, batch_c%length)
!
      enddo  ! batch over c 
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_a2_mlcc2
!
!
   module subroutine jacobian_cc2_b1_mlcc2(wf, rho_ai, c_aibj, n_cc2_o, n_cc2_v, &
                                          first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v)
!!
!!    Jacobian CC2 B1
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Calculates the B1 term
!!
!!       B1: 2 sum_bj F_jb c_aibj - F_jb c_ajbi  
!!           - sum_kjb L_kijb c_akbj + sum_bkc L_abkc c_bick
!!
!!    And adds it to rho_ai.
!!
!!    The fourth term is calculated in batches of index a.
!!    
!!    Index restrictions:
!!
!!       Terms 1 and 2:
!!
!!          a, i, b, j : CC2 orbitals
!!
!!       Term 3:
!!
!!          a, k, b, j : CC2 orbitals
!!
!!          i : unrestricted
!!
!!       Term 4:
!!
!!          c, k, b, i : CC2 orbitals
!!
!!          a : unrestricted
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v, last_cc2_o, last_cc2_v
!
      real(dp), dimension(n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o), intent(in)  :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(out)                     :: rho_ai   
!
      real(dp), dimension(:,:,:,:), allocatable :: X_aijb 
      real(dp), dimension(:,:,:,:), allocatable :: g_jbki 
      real(dp), dimension(:,:,:,:), allocatable :: g_abkc
      real(dp), dimension(:,:,:,:), allocatable :: L_abkc
      real(dp), dimension(:,:,:,:), allocatable :: L_kbji
      real(dp), dimension(:,:,:,:), allocatable :: c_bkci 
!
      real(dp), dimension(:,:), allocatable :: F_jb, rho_ai_copy
!
      type(batching_index) :: batch_a 
      integer              :: req0, req1, current_a_batch
      integer              :: a, b, i, j
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCC2 B1 transformation', pl='verbose')
      call timer%turn_on()
!
!     Terms 1 and 2: 2 sum_bj F_jb c_aibj - F_jb c_ajbi = sum_bj X_aijb F_jb 
!
!        a, i, b, j : CC2 orbitals
!
!     Make X_aijb = 2 c_aibj - c_ajbi 
!
      call mem%alloc(X_aijb, n_cc2_v, n_cc2_o, n_cc2_o, n_cc2_v)
      call zero_array(X_aijb, (n_cc2_v**2)*(n_cc2_o**2))
!
      call add_1243_to_1234(two, c_aibj, X_aijb, n_cc2_v, n_cc2_o, n_cc2_o, n_cc2_v)
      call add_1342_to_1234(-one, c_aibj, X_aijb, n_cc2_v, n_cc2_o, n_cc2_o, n_cc2_v)
!
      call mem%alloc(F_jb, n_cc2_o, n_cc2_v)
!
!$omp parallel do private(b, j)
      do b = 1, n_cc2_v
         do j = 1, n_cc2_o
!
            F_jb(j, b) = wf%fock_ia(j + first_cc2_o - 1, b + first_cc2_v - 1)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(rho_ai_copy, n_cc2_v, n_cc2_o)
!
      call dgemm('N','N',              &
                  (n_cc2_o)*(n_cc2_v), &
                  1,                   &
                  (n_cc2_o)*(n_cc2_v), &
                  one,                 &
                  X_aijb,              & ! X_ai_jb 
                  (n_cc2_o)*(n_cc2_v), &
                  F_jb,                & 
                  (n_cc2_o)*(n_cc2_v), &
                  zero,                &
                  rho_ai_copy,         &
                  (n_cc2_o)*(n_cc2_v))
!
      call mem%dealloc(X_aijb, n_cc2_v, n_cc2_o, n_cc2_o, n_cc2_v)
      call mem%dealloc(F_jb, n_cc2_o, n_cc2_v)
!
!$omp parallel do private(a, i) collapse(2)
      do i = 1, n_cc2_o
         do a = 1, n_cc2_v
!
            rho_ai(a + first_cc2_v - 1, i + first_cc2_o - 1) = &
            rho_ai(a + first_cc2_v - 1, i + first_cc2_o - 1) + rho_ai_copy(a,i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(rho_ai_copy, n_cc2_v, n_cc2_o)
!
!     Term 3: - sum_kjb c_akbj L_jbki
!
!        a, k, b, j : CC2 orbitals
!
!        i : unrestricted
!
!     Make L_jbki ordered as L_kbji (= 2 g_jbki - g_kbji)
!
      call mem%alloc(g_jbki, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
      call wf%eri%get_eri_t1('ovoo', g_jbki, first_cc2_o, last_cc2_o, first_cc2_v, last_cc2_v, &
                                             first_cc2_o, last_cc2_o, 1, wf%n_o)
!
      call mem%alloc(L_kbji, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call copy_and_scale(-one, g_jbki, L_kbji, (n_cc2_o**2)*n_cc2_v*wf%n_o)
      call add_3214_to_1234(two, g_jbki, L_kbji, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call mem%dealloc(g_jbki, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
      call dgemm('N','N',                 &
                  n_cc2_v,                &
                  wf%n_o,                 &
                  (n_cc2_v)*(n_cc2_o**2), &
                  -one,                   &
                  c_aibj,                 & ! c_a,kbj
                  n_cc2_v,                &
                  L_kbji,                 & ! L_kbj,i 
                  (n_cc2_v)*(n_cc2_o**2), &
                  one,                    &
                  rho_ai(first_cc2_v, 1), &
                  wf%n_v)
!
      call mem%dealloc(L_kbji, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_o)
!
!     Term 4: sum_bkc L_abkc c_ckbi, batch over a
!
!        c, k, b, i : CC2 orbitals
!
!        a : unrestricted
!
!     Order c_ckbi as c_bkci
!
      call mem%alloc(c_bkci, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
      call sort_1234_to_3214(c_aibj, c_bkci, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
      req0 = (n_cc2_o)*(n_cc2_v)*(wf%eri%n_J)
      req1 = max((n_cc2_v)*(wf%eri%n_J) + (n_cc2_o)*(n_cc2_v)**2, 2*(n_cc2_o)*(n_cc2_v)**2)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_abkc, batch_a%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call wf%eri%get_eri_t1('vvov', g_abkc, batch_a%first, batch_a%last, &
                                                first_cc2_v, last_cc2_v,     &   
                                                first_cc2_o, last_cc2_o,     &
                                                first_cc2_v, last_cc2_v)
!
         call mem%alloc(L_abkc, batch_a%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call copy_and_scale(two, g_abkc, L_abkc, batch_a%length*(n_cc2_v**2)*n_cc2_o)
         call add_1432_to_1234(-one, g_abkc, L_abkc, batch_a%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call mem%dealloc(g_abkc, batch_a%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
         call dgemm('N','N',                             &
                     batch_a%length,                     &
                     n_cc2_o,                            &
                     (n_cc2_o)*(n_cc2_v)**2,             &
                     one,                                &
                     L_abkc,                             & ! L_a,bkc 
                     batch_a%length,                     &
                     c_bkci,                             & ! c_bkc,i
                     (n_cc2_o)*(n_cc2_v)**2,             &
                     one,                                &
                     rho_ai(batch_a%first, first_cc2_o), &
                     wf%n_v)
!
         call mem%dealloc(L_abkc, batch_a%length, n_cc2_v, n_cc2_o, n_cc2_v)
!
      enddo ! batch_a
!
      call mem%dealloc(c_bkci, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_b1_mlcc2
!
!
   module subroutine jacobian_cc2_b2_mlcc2(wf, rho_aibj, c_aibj)
!!
!!    Jacobian CC2 B2
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Constructs the B2 term
!!
!!       B2: ε_aibj c_aibj/(1/Δ_aibj) 
!!
!!    and adds it to rho_aibj.
!!
!!    Index restrictions:
!!
!!       a, i, b, j : CC2 orbitals
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)     :: c_aibj
      real(dp), dimension(wf%n_cc2_v, wf%n_cc2_o, wf%n_cc2_v, wf%n_cc2_o), intent(inout)     :: rho_aibj   
!
      integer :: i, j, a, b
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian MLCC2 B2 transformation', pl='verbose')
      call timer%turn_on()
!
!     c_aibj/(1/Δ_aibj) 
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_cc2_v
         do i = 1, wf%n_cc2_o
!
            c_aibj(a, i, a, i) = half*c_aibj(a, i, a, i) 
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(a, i, b, j)
      do b = 1, wf%n_cc2_v
         do j = 1, wf%n_cc2_o
            do i = 1, wf%n_cc2_o
               do a = 1, wf%n_cc2_v
!
                  rho_aibj(a,i,b,j) = rho_aibj(a,i,b,j) + c_aibj(a,i,b,j)*&
                                          (- wf%orbital_energies(i) &
                                           - wf%orbital_energies(j) &
                                           + wf%orbital_energies(wf%n_o + a) &
                                           + wf%orbital_energies(wf%n_o + b) )
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call timer%turn_off()
!
   end subroutine jacobian_cc2_b2_mlcc2
!
!
   module subroutine save_jacobian_a1_intermediates_mlcc2(wf, n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v)
!!
!!    Save jacobian a1 intermediates
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adapted from jacobian_cc2.F90
!!    written by Eirik F. Kjønstad, Sarai D. Folkestad, 2018
!!
!!    Constructs the intermediates 
!!
!!       Y_ji = - sum_bck g_jbkc u_ckbi
!!       Y_ab = - sum_jck u_ajck g_kcjb
!!      
!!    Which are constructed in save_jacobian_a1_intermediates
!!    and stored on files
!!
!!       jacobian_a1_intermediate_oo
!!       jacobian_a1_intermediate_vv
!!
!!    which are wavefunction variables
!!
!!    Index restrictions:
!!
!!       oo intermediate:
!!    
!!          c, k, b, i : CC2 orbitals
!!
!!          j : unrestricted
!!
!!       vv intermediate:
!!
!!          a, j, k, c : CC2 orbitals
!!
!!          b : unrestricted 
!!
      implicit none
!
      class(mlcc2) :: wf
!
      integer, intent(in) :: n_cc2_o, n_cc2_v, first_cc2_o, first_cc2_v
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jbkc, u_bkci, g_jckb
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:,:), allocatable :: Y_Ab
      real(dp), dimension(:,:), allocatable :: Y_jI
!
      integer :: c, k, b, j
!
      timer = timings('Jacobian mlcc2 A1 intermediate construction', pl='verbose')
      call timer%turn_on()
!
!     Construct integrals for both intermediates
!
      call mem%alloc(g_jbkc, wf%n_o, n_cc2_v, n_cc2_o, wf%n_v)  
      call wf%eri%get_eri_t1('ovov', g_jbkc,                         &
                             1, wf%n_o,                              &
                             first_cc2_v, first_cc2_v + n_cc2_v - 1, &
                             first_cc2_o, first_cc2_o + n_cc2_o - 1, &
                             1, wf%n_v)
!
!     Y_ji = - sum_bck g_jbkc u_ckbi
!
!     Reorder u_ckbi as u_bkci
!
      call mem%alloc(u_bkci, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
      call sort_1234_to_3214(wf%u_aibj, u_bkci, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
      call mem%alloc(Y_ji, wf%n_o, n_cc2_o)
!
      call dgemm('N', 'N',                &
                  wf%n_o,                 &
                  n_cc2_o,                &
                  (n_cc2_v**2)*(n_cc2_o), &
                  one,                    &
                  g_jbkc,                 & ! g_j_bkc
                  wf%n_o,                 &
                  u_bkci,                 & ! u_bkc_i 
                  (n_cc2_v**2)*(n_cc2_o), &
                  zero,                   &
                  Y_ji,                   &
                  wf%n_o)
!
      call mem%dealloc(u_bkci, n_cc2_v, n_cc2_o, n_cc2_v, n_cc2_o)
!
      wf%jacobian_a1_intermediate_oo = sequential_file('jacobian_a1_intermediate_oo_doubles')
      call wf%jacobian_a1_intermediate_oo%open_('write', 'rewind')
!
      call wf%jacobian_a1_intermediate_oo%write_(Y_ji, n_cc2_o*wf%n_o)
!
      call mem%dealloc(Y_ji, wf%n_o, n_cc2_o)
!
      call wf%jacobian_a1_intermediate_oo%close_('keep')
!
!     Term 3: Y_ab  - sum_jck u_ajck g_kcjb 
!
      call mem%alloc(g_jckb, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_v)
!
!$omp parallel do private(b, j, c, k)
      do b = 1, wf%n_v
         do k = 1, n_cc2_o
            do c = 1, n_cc2_v
               do j = 1, n_cc2_o
!
                  g_jckb(j, c, k, b) = g_jbkc(k + first_cc2_o - 1, c, j, b)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_jbkc, wf%n_o, n_cc2_v, n_cc2_o, wf%n_v)  
!
      call mem%alloc(Y_Ab, n_cc2_v, wf%n_v)
!
      call dgemm('N', 'N',                &
                  n_cc2_v,                &
                  wf%n_v,                 &
                  n_cc2_v*(n_cc2_o**2),   &
                  one,                    &
                  wf%u_aibj,              & ! u_ajck
                  n_cc2_v,                &
                  g_jckb,                 & 
                  n_cc2_v*(n_cc2_o**2),   &
                  zero,                   &
                  Y_Ab,                   &
                  n_cc2_v)
!
      call mem%dealloc(g_jckb, n_cc2_o, n_cc2_v, n_cc2_o, wf%n_v)
!
      wf%jacobian_a1_intermediate_vv = sequential_file('jacobian_a1_intermediate_vv_doubles')
      call wf%jacobian_a1_intermediate_vv%open_('write', 'rewind')
!
      call wf%jacobian_a1_intermediate_vv%write_(Y_Ab, n_cc2_v*wf%n_v)
!
      call mem%dealloc(Y_Ab, n_cc2_v, wf%n_v)
!
      call wf%jacobian_a1_intermediate_vv%close_('keep')
!
      call timer%turn_off()
!
   end subroutine save_jacobian_a1_intermediates_mlcc2
!
!
end submodule jacobian_mlcc2
