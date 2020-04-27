!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
submodule (doubles_class) jacobian_doubles
!
!!
!!    Jacobian submodule
!!
!!    Routines for the linear transform of trial
!!    vectors by the Jacobian matrix 
!!
!!    ρ_i = A * c_i,
!!
!!    where
!!   
!!    A_μ,ν = < μ | exp(-T) [H, τ_ν] exp(T) | R >.
!!  
! 
   implicit none
!
!
contains
!
!
   module subroutine save_jacobian_a1_intermediates_doubles(wf)
!!
!!    Save jacobian a1 intermediates
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2019
!!
!!    Constructs the intermediates 
!!
!!       Y_jl   = t_ckdj L_kcld 
!!       Y_bd   = t_blck L_kcld 
!!
!!    Which are constructed in save_jacobian_a1_intermediates
!!    and stored on files
!!
!!       jacobian_a1_intermediate_oo
!!       jacobian_a1_intermediate_vv
!!
!!    which are wavefunction variables
!!
      implicit none
!
      class(doubles) :: wf
!
      type(timings), allocatable :: jacobian_a1_intermediate_timer
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ldkc
      real(dp), dimension(:,:,:,:), allocatable :: L_dlck
      real(dp), dimension(:,:,:,:), allocatable :: t_blck
!
      real(dp), dimension(:,:), allocatable :: Y_bd
      real(dp), dimension(:,:), allocatable :: Y_jl
!
      jacobian_a1_intermediate_timer = &
         timings('Jacobian doubles G2 intermediate construction', pl='verbose')
!
      call jacobian_a1_intermediate_timer%turn_on()
!
      call mem%alloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
      call wf%get_ovov(g_ldkc)
!
      call mem%alloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call zero_array(L_dlck, (wf%n_v**2)*(wf%n_o**2))
      call add_2143_to_1234(two, g_ldkc, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-one, g_ldkc, L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ldkc, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call mem%alloc(t_blck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_blck, wf%n_t1)
!
!     Y_bd = t_blck L_ldkc
!
      call mem%alloc(Y_bd, wf%n_v, wf%n_v)
!
      call dgemm('N', 'T',             &
                  wf%n_v,              &
                  wf%n_v,              &
                  wf%n_v*(wf%n_o**2),  &
                  one,                 &
                  t_blck,              & ! t_b_lck
                  wf%n_v,              &
                  L_dlck,              & ! L_d_lck
                  wf%n_v,              &
                  zero,                &
                  Y_bd,                &
                  wf%n_v)
!
      wf%jacobian_a1_intermediate_vv = sequential_file('jacobian_a1_intermediate_vv_doubles')
      call wf%jacobian_a1_intermediate_vv%open_('write', 'rewind')
!
      call wf%jacobian_a1_intermediate_vv%write_(Y_bd, wf%n_v**2)
!
      call mem%dealloc(Y_bd, wf%n_v, wf%n_v)
!
      call wf%jacobian_a1_intermediate_vv%close_('keep')
!
!     Y_jl = t_ckdj L_kcld 
!
!     Note: pretend that t_blck is t_ckdj
!     Using symmetry L_dlck = L_ckdl (L_kcld = L_ldkc)
!
      call mem%alloc(Y_jl, wf%n_o, wf%n_o)
!
      call dgemm('T', 'N',             &
                  wf%n_o,              &
                  wf%n_o,              &
                  wf%n_o*(wf%n_v**2),  &
                  one,                 &
                  t_blck,              & ! t_ckd_j
                  wf%n_o*(wf%n_v**2),  &
                  L_dlck,              & ! L_ckd_l
                  wf%n_o*(wf%n_v**2),  &
                  zero,                &
                  Y_jl,                &
                  wf%n_o)
!
      call mem%dealloc(t_blck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_dlck, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      wf%jacobian_a1_intermediate_oo = sequential_file('jacobian_a1_intermediate_oo_doubles')
      call wf%jacobian_a1_intermediate_oo%open_('write', 'rewind')
!
      call wf%jacobian_a1_intermediate_oo%write_(Y_jl, wf%n_o**2)
!
      call mem%dealloc(Y_jl, wf%n_o, wf%n_o)
!
      call wf%jacobian_a1_intermediate_oo%close_('keep')
!
      call jacobian_a1_intermediate_timer%turn_off()
!
   end subroutine save_jacobian_a1_intermediates_doubles
!
!
   module subroutine jacobian_doubles_a1_doubles(wf, rho_ai, c_ai)
!!
!!    Jacobian doubles A1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^A1 = sum_ckdl L_kcld (u_ki^ca c_dl - t_kl^ad c_ci  - t_ki^cd c_al)
!!              = sum_ckdl L_kcld u_ki^ca c_dl - Y_ac c_ci - Y_il c_al)
!!
      implicit none
!
      class(doubles) :: wf
!
!     Vectors sent to the routine
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o)                :: rho_ai
!
!     Integrals
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kcld
      real(dp), dimension(:,:,:,:), allocatable :: L_kcdl
!
!     Intermediates
!
      real(dp), dimension(:,:), allocatable :: X_kc
      real(dp), dimension(:,:), allocatable :: Y_il
      real(dp), dimension(:,:), allocatable :: Y_ac
!
!     Amplitudes
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ckai
      real(dp), dimension(:,:,:,:), allocatable :: u_aikc
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian doubles A1', pl='verbose')
      call timer%turn_on()
!
!     Term 1: sum_ckdl L_kcld u_ki^ca c_dl ::
!
      call mem%alloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_kcld)
!
      call mem%alloc(L_kcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call zero_array(L_kcdl, wf%n_t1**2)
!
!     L_kcdl(lc,dk) = L_kcld = 2*g_kcld - g_ldkc 
!
      call add_1243_to_1234(two, g_kcld, L_kcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
      call add_1342_to_1234(-one, g_kcld, L_kcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_kcld, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
!     X_kc = sum_kd L_kcld c_dl 
!
      call mem%alloc(X_kc, wf%n_o, wf%n_v)
!
      call dgemv('N',                &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  L_kcdl,            &
                  (wf%n_o)*(wf%n_v), &
                  c_ai,              & ! c_dl
                  1,                 &
                  zero,              &
                  X_kc,              &
                  1)
!
      call mem%dealloc(L_kcdl, wf%n_o, wf%n_v, wf%n_v, wf%n_o)
!
!     Form u_ai_kc = u_ki^ca = 2 * t_ki^ca - t_ik^ca 
!
      call mem%alloc(t_ckai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_ckai, (wf%n_o)*(wf%n_v))
!
!     u_ai_kc(ai, lc) = 2 * t_ki^ca - t_ik^ca 
!
      call mem%alloc(u_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(u_aikc, wf%n_t1**2)
!
      call add_4312_to_1234(two, t_ckai, u_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_4213_to_1234(-one, t_ckai, u_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call mem%dealloc(t_ckai, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     rho_ai =+ sum_lc u_ai_lc X_lc
!
      call dgemv('N',                &
                  (wf%n_v)*(wf%n_o), &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  u_aikc,            & 
                  (wf%n_v)*(wf%n_o), &
                  X_kc,              &
                  1, &
                  one,               &
                  rho_ai,            &
                  1)
!
      call mem%dealloc(u_aikc, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call mem%dealloc(X_kc, wf%n_o, wf%n_v)
!
!     Term 2: - Y_ac c_ci = - t_kl^ad L_kcld c_ci
!
      call mem%alloc(Y_ac, wf%n_v, wf%n_v)
!
      call wf%jacobian_a1_intermediate_vv%open_('read', 'rewind')
!
      call wf%jacobian_a1_intermediate_vv%read_(Y_ac, wf%n_v**2)
!
      call wf%jacobian_a1_intermediate_vv%close_('keep')
!
      call dgemm('N', 'N',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  -one,       &
                  Y_ac,       &
                  wf%n_v,     &
                  c_ai,       & ! c_ci
                  wf%n_v,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)
!
      call mem%dealloc(Y_ac, wf%n_v, wf%n_v)
!
!     Term 3: - Y_il c_al = - t_ki^cd c_al = - t_ki^cd L_kcld c_al
!
      call mem%alloc(Y_il, wf%n_o, wf%n_o)
!
      call wf%jacobian_a1_intermediate_oo%open_('read', 'rewind')
!
      call wf%jacobian_a1_intermediate_oo%read_(Y_il, wf%n_o**2)
!
      call wf%jacobian_a1_intermediate_oo%close_('keep')
!
      call dgemm('N', 'T',    &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  c_ai,       & ! c_ak 
                  wf%n_v,     &
                  Y_il,       &
                  wf%n_o,     &
                  one,        &
                  rho_ai,     &
                  wf%n_v)
!
      call mem%dealloc(Y_il, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_doubles_a1_doubles
!
!
 module subroutine jacobian_doubles_b1_doubles(wf, rho_ai, c_aibj)
!!
!!    Jacobian doubles B1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!!    rho_ai^B1 = sum_bj F_jb (2*c_aibj - c_ajbi)
!!              = sum_bj F_jb v_aijb
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: v_aijb
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian doubles B1', pl='verbose')
      call timer%turn_on()
!
!     Construct v_aibj = 2*c_aibj - c_ajbi ordered as
!
!        v_ai_jb(a,i,j,b) = 2*c_aibj(a,i,b,j) - c_aibj(a,j,b,i)
!
!     and do the matrix multiplication with F_jb
!
      call mem%alloc(v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call zero_array(v_aijb, wf%n_t1**2)
!
      call add_1243_to_1234(two, c_aibj, v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call add_1342_to_1234(-one, c_aibj, v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemv('N',                &
                  (wf%n_o)*(wf%n_v), &
                  (wf%n_o)*(wf%n_v), &
                  one,               &
                  v_aijb,            & ! v_ai,jb
                  (wf%n_o)*(wf%n_v), &
                  wf%fock_ia,        & ! F_jb
                  1,                 &
                  one,               &
                  rho_ai,            & ! rho_ai
                  1)
!
      call mem%dealloc(v_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call timer%turn_off()
!
   end subroutine jacobian_doubles_b1_doubles
!
!
   module subroutine jacobian_doubles_c1_doubles(wf, rho_ai, c_aibj)
!!
!!    Jacobian doubles C1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^C1 = - sum_bjk L_jikb c_ajbk
!!              = - sum_bjk (2*g_jikb - g_kijb) c_ajbk
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_aibj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_jikb
      real(dp), dimension(:,:,:,:), allocatable :: L_jbki
!
      type(timings), allocatable :: timer
!
      timer = timings('Jacobian doubles C1', pl='verbose')
      call timer%turn_on()
!
!     Construct L_jikb = 2*g_jikb - g_kijb as
!
!        L_jb_ki(jb,ki) = 2*g_ji_kb(ji,kb) - g_ji_kb(ki,jb)
!
!     and then contract with c_ajbk = c_aibj(aj,bk).
!
      call mem%alloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_ooov(g_jikb)
!
      call mem%alloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(L_jbki, (wf%n_o**3)*wf%n_v)
!
      call add_1432_to_1234(two, g_jikb, L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call add_3412_to_1234(-one, g_jikb, L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_jikb, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call dgemm('N', 'N',                &
                  wf%n_v,                 &
                  wf%n_o,                 &
                  (wf%n_v)*((wf%n_o)**2), &
                  -one,                   &
                  c_aibj,                 & ! c_a,jbk
                  wf%n_v,                 &
                  L_jbki,                 & ! L_jbk,i
                  (wf%n_v)*((wf%n_o)**2), &
                  one,                    &
                  rho_ai,                 & ! rho_ai
                  wf%n_v)
!
      call mem%dealloc(L_jbki, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_doubles_c1_doubles
!
!
   module subroutine jacobian_doubles_d1_doubles(wf, rho_ai, c_bicj)
!!
!!    Jacobian doubles D1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_ai^D1 =  sum_bcj L_abjc c_bicj
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in) :: c_bicj
      real(dp), dimension(wf%n_v, wf%n_o) :: rho_ai
!
      integer :: current_a_batch
!
      type(batching_index) :: batch_a
!
      real(dp), dimension(:,:,:,:), allocatable :: c_bjci
      real(dp), dimension(:,:,:,:), allocatable :: g_abjc
      real(dp), dimension(:,:,:,:), allocatable :: L_abjc
!
      type(timings), allocatable :: timer
!
      integer :: req0, req1
!
      timer = timings('Jacobian doubles D1', pl='verbose')
      call timer%turn_on()
!
!     Reorder c_bicj to c_bjci
!
      call mem%alloc(c_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call sort_1234_to_1432(c_bicj, c_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Prepare for batching over index a
!
      req0 = wf%n_o*wf%integrals%n_J*wf%n_v
!
      req1 = wf%n_v*wf%integrals%n_J + 2*(wf%n_v**2)*(wf%n_o)
!
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches
!
!        Determine the limits for the current a-batch
!
         call batch_a%determine_limits(current_a_batch)
!
!        Construct L_abjc = 2 g_abjc - g_acjb
!
         call mem%alloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call wf%get_vvov(g_abjc,                        &
                           batch_a%first, batch_a%last,  &
                           1, wf%n_v,                    &
                           1, wf%n_o,                    &
                           1, wf%n_v)
!
         call mem%alloc(L_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call copy_and_scale(two, g_abjc, L_abjc, (wf%n_v**2)*wf%n_o*batch_a%length)
         call add_1432_to_1234(-one, g_abjc, L_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call mem%dealloc(g_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
         call dgemm('N', 'N',                   &
                     batch_a%length,            &
                     wf%n_o,                    &
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     L_abjc,                    & ! L_a,bjc
                     batch_a%length,            &
                     c_bjci,                    & ! c_bjc,i
                     (wf%n_o)*(wf%n_v)**2,      &
                     one,                       &
                     rho_ai(batch_a%first, 1),  &
                     wf%n_v)
!
         call mem%dealloc(L_abjc, batch_a%length, wf%n_v, wf%n_o, wf%n_v)
!
      enddo ! End batching over a
!
      call mem%dealloc(c_bjci, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine jacobian_doubles_d1_doubles
!
!
   module subroutine jacobian_doubles_a2_doubles(wf, rho_aibj, c_ai)
!!
!!    Jacobian doubles A2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2017-2018
!!
!!    rho_aibj^A2 = sum_c g_aibc c_cj - sum_k g_aikj c_bk
!!
      implicit none
!
      class(doubles) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)          :: c_ai
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o)      :: rho_aibj
!
      real(dp), dimension(:,:,:,:), allocatable :: g_kjai
      real(dp), dimension(:,:,:,:), allocatable :: g_aibc
!
      integer :: current_b_batch
!
      type(batching_index) :: batch_b
!
      type(timings), allocatable :: timer
!
      integer :: req0, req1
!
      timer = timings('Jacobian doubles A2', pl='verbose')
      call timer%turn_on()
!
!     :: Term 1. - sum_k g_aikj c_bk ::
!
!     rho_bjai =- c_bk g_kjai 
!
      call mem%alloc(g_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%get_oovo(g_kjai)
!
      call dgemm('N', 'N',              &
                  wf%n_v,               &
                  (wf%n_v)*(wf%n_o)**2, &
                  wf%n_o,               &
                  -one,                 &
                  c_ai,                 & ! c_b,k
                  wf%n_v,               &
                  g_kjai,               & ! g_k,jai
                  wf%n_o,               &
                  one,                  &
                  rho_aibj,             & ! rho_b,jai -> will be (ai,bj)-symmetrized 
                  wf%n_v)
!
      call mem%dealloc(g_kjai, wf%n_o, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Term 2. rho_aibj =+ sum_c g_aibc c_cj ::
!
!     We do the matrix multiplication as g_aib_c c_cj, batching over b.
!
      req0 = wf%n_o*wf%integrals%n_J*wf%n_v
!
      req1 = wf%n_v*wf%integrals%n_J + (wf%n_v**2)*(wf%n_o)
!
      batch_b = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_b, req0, req1)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
!        Calculate the batch contribution to rho_aib_j = sum_c g_aibc c_cj
!
         call mem%alloc(g_aibc, wf%n_v, wf%n_o, wf%n_v, batch_b%length)
!
         call wf%get_vovv(g_aibc,        &
                          1,             &
                          wf%n_v,        &
                          1,             &
                          wf%n_o,        &
                          batch_b%first, &
                          batch_b%last,  &
                          1,             &
                          wf%n_v)
!
         call dgemm('N', 'N',                            &
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     wf%n_o,                             &
                     wf%n_v,                             &
                     one,                                &
                     g_aibc,                             & ! g_aib,c
                     (wf%n_v)*(wf%n_o)*(batch_b%length), &
                     c_ai,                               & ! c_c,j
                     wf%n_v,                             &
                     one,                                &
                     rho_aibj(1,1, batch_b%first,1),     &
                     (wf%n_o)*(wf%n_v)**2)
!
         call mem%dealloc(g_aibc, wf%n_v, wf%n_o, wf%n_v, batch_b%length)
!
      enddo ! End of batches over b
!
      call timer%turn_off()
!
   end subroutine jacobian_doubles_a2_doubles
!
!
end submodule jacobian_doubles
