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
submodule (cc3_class) batching_abc
!
!!
!!    Batching abc submodule
!!
!!    Routines that construct the triples T-amplitudes and
!!    components of the excitation vectors in batches 
!!    of the virtual indices a,b,c
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_cc3_integrals_t3_abc_cc3(wf)
!!
!!    Prepare integral files t3 amplitudes in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    (bd|ck) ordered as dk,bc
!!    (lj|ck) ordered as ljk,c
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      type(batching_index) :: batch_c
!
      integer :: req_0, req_c
      integer :: c_batch
!
      batch_c = batching_index(wf%n_v)
!
!     (bd|ck) stored dk,bc
!
      req_0 = wf%integrals%n_J*wf%n_v**2
      req_c = 2*wf%n_o*wf%n_v**2 + wf%integrals%n_J*wf%n_o
!
      call mem%batch_setup(batch_c, req_0, req_c)
!
      wf%g_bdck_t_v = direct_stream_file('g_bdck_t_v', wf%n_v*wf%n_o)
      call wf%g_bdck_t_v%open_('write')
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, batch_c%length, wf%n_o)
         call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_v, batch_c%length)
!
         call wf%get_vvvo(g_pqrs,                        &
                           1, wf%n_v,                    &
                           1, wf%n_v,                    &
                           batch_c%first, batch_c%last,  &
                           1, wf%n_o)
!
         call sort_1234_to_2413(g_pqrs, h_pqrs, wf%n_v, wf%n_v, batch_c%length, wf%n_o)
!
         call wf%g_bdck_t_v%write_compound_full_batch(h_pqrs, wf%n_v, batch_c)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, batch_c%length, wf%n_o)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_v, batch_c%length)
!
      enddo
!
      call wf%g_bdck_t_v%close_()
!
!     (lj|ck) stored ljk,c
!
      req_0 = wf%integrals%n_J*wf%n_o**2
      req_c = 2*wf%n_o**3 + wf%integrals%n_J*wf%n_o
!
      call mem%batch_setup(batch_c, req_0, req_c)
!
      wf%g_ljck_t_v = direct_stream_file('g_ljck_t_v', wf%n_o**3)
      call wf%g_ljck_t_v%open_('write')
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, batch_c%length, wf%n_o)
         call mem%alloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
!
         call wf%get_oovo(g_pqrs,                        &
                           1, wf%n_o,                    &
                           1, wf%n_o,                    &
                           batch_c%first, batch_c%last,  &
                           1, wf%n_o)
!
         call sort_1234_to_1243(g_pqrs, h_pqrs, wf%n_o, wf%n_o, batch_c%length, wf%n_o)
!
         call wf%g_ljck_t_v%write_interval(h_pqrs, batch_c)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, batch_c%length, wf%n_o)
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
!
      enddo
!
      call wf%g_ljck_t_v%close_()
!
   end subroutine prepare_cc3_integrals_t3_abc_cc3
!
!
   module subroutine prepare_cc3_integrals_R3_abc_cc3(wf, R_ai)
!!
!!    Prepare integral files R3 amplitudes in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck')   ordered as dk,bc
!!    g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k)   ordered as ljk,c
!!
!!    NB: The integrals (bd|ck) and (lj|ck) constructed in 
!!        prepare_cc3_integrals_t3_abc_cc3 are also needed
!!
!!    Based on construct_c1_integrals_cc3 
!!    written by Rolf H. Myhre and Alexander C. Paul
!!
!
      implicit none
!
      class(cc3) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ck_c1      ! c1 transformed Cholesky vector
      real(dp), dimension(:,:,:), allocatable :: L_J_ck, L_J_bd ! Cholesky vectors
!
      type(batching_index) :: batch_c
!
      integer :: req_0, req_c
      integer :: c_batch
!
      batch_c = batching_index(wf%n_v)
!
!     g'_bdck = (b'd|ck) + (bd|c'k) + (bd|ck') ordered dk,bc
!
      req_0 = 2*wf%integrals%n_J*wf%n_v**2
      req_c = max((wf%integrals%n_J)*(wf%n_o) + (wf%n_v)**3, 2*(wf%n_v)**2*wf%n_o)
!
      call mem%alloc(L_J_bd, wf%integrals%n_J, wf%n_v, wf%n_v)
!
      wf%g_bdck_c_v = direct_stream_file('g_bdck_c_v', wf%n_v*wf%n_o)
      call wf%g_bdck_c_v%open_('write')
!
      call mem%batch_setup(batch_c, req_0, req_c)
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
!        :: Term 1: g_b'dck += sum_J L_J_b'd L_J_ck ::
!
!        here L_J_bd is used for the R1-transformed cholesky-vector 
         call wf%integrals%construct_cholesky_ab_c1(L_J_bd, R_ai, 1, wf%n_v, 1,  wf%n_v)
!
         call mem%alloc(L_J_ck, wf%integrals%n_J, batch_c%length, wf%n_o)
         call wf%integrals%get_cholesky_t1(L_J_ck, wf%n_o + batch_c%first, &
                                            wf%n_o + batch_c%last, 1, wf%n_o)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, batch_c%length, wf%n_o)
!
         call dgemm('T', 'N',                   &
                     (wf%n_v)**2,               &
                     (batch_c%length)*(wf%n_o), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_bd,                    & ! L_J_bd b is R1-transformed
                     wf%integrals%n_J,          &
                     L_J_ck,                    & ! L_J_ck
                     wf%integrals%n_J,          & 
                     zero,                      &
                     g_pqrs,                    & ! (b'd|ck)
                     (wf%n_v)**2)
!
         call mem%dealloc(L_J_ck, wf%integrals%n_J, batch_c%length, wf%n_o)
!
!        :: Term 2: g_bdc'k += sum_J L_J_bd L_J_c'k_c1 ::
!
         call mem%alloc(L_J_ck_c1, wf%integrals%n_J, batch_c%length, wf%n_o)
         call wf%integrals%construct_cholesky_ai_a_c1(L_J_ck_c1, R_ai, batch_c%first,  &
                                                      batch_c%last, 1, wf%n_o)
!
         call wf%integrals%get_cholesky_t1(L_J_bd, wf%n_o + 1, wf%n_o + wf%n_v, &
                                            wf%n_o + 1, wf%n_o + wf%n_v)
!
         call dgemm('T', 'N',                   &
                     (wf%n_v)**2,               &
                     (batch_c%length)*(wf%n_o), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_bd,                    & ! L_J_bd
                     wf%integrals%n_J,          &
                     L_J_ck_c1,                 & ! L_J_ck c is R1-transformed
                     wf%integrals%n_J,          & 
                     one,                       &
                     g_pqrs,                    & ! (bd|c'k)
                     (wf%n_v)**2)
!
!        :: Term 3: g_bdck' += sum_J L_J_bd L_J_ck'_c1 ::
!
         call wf%integrals%construct_cholesky_ai_i_c1(L_J_ck_c1, R_ai, batch_c%first,  &
                                                      batch_c%last, 1, wf%n_o)
!
         call dgemm('T', 'N',                   &
                     (wf%n_v)**2,               &
                     (batch_c%length)*(wf%n_o), &
                     wf%integrals%n_J,          &
                     one,                       &
                     L_J_bd,                    & ! L_J_bd
                     wf%integrals%n_J,          &
                     L_J_ck_c1,                 & ! L_J_ck k is R1-transformed
                     wf%integrals%n_J,          & 
                     one,                       &
                     g_pqrs,                    & ! (bd|c'k)
                     (wf%n_v)**2)
!
         call mem%dealloc(L_J_ck_c1, wf%integrals%n_J, batch_c%length, wf%n_o)
!
!        Sort from g_pqrs = (b'd|ck) + (bd|c'k) + (bd|ck') to h_pqrs ordered as dk,bc
!
         call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_v, batch_c%length)
!
         call sort_1234_to_2413(g_pqrs, h_pqrs, wf%n_v, wf%n_v, batch_c%length, wf%n_o)
!
         call wf%g_bdck_c_v%write_compound_full_batch(h_pqrs, wf%n_v, batch_c)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, batch_c%length, wf%n_o)
         call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_v, batch_c%length)
!
      enddo ! c_batch
!
      call mem%dealloc(L_J_bd, wf%integrals%n_J, wf%n_v, wf%n_v)
!
      call wf%g_bdck_c_v%close_()
!
!
!     g'_ljck = (lj'|ck) + (lj|ck') + (lj|c'k)   ordered as ljk,c
!
!     assume it's possible to hold n_v*n_o**3 and g_ljck_c1 is already on file
!
      wf%g_ljck_c_v = direct_stream_file('g_ljck_c_v', wf%n_o**3)
      call wf%g_ljck_c_v%open_('write')
!
      call wf%g_ljck_c%open_('read')
!
      call mem%alloc(g_pqrs, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call wf%g_ljck_c%read_(g_pqrs, 1, wf%n_o*wf%n_o) ! g'_ljck ordered lcjk
!
      call mem%alloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
      call sort_1234_to_1342(g_pqrs, h_pqrs, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_pqrs, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%g_ljck_c_v%write_(h_pqrs, 1, wf%n_v)
!
      call mem%dealloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%g_ljck_c_v%close_()
      call wf%g_ljck_c%close_()
!
   end subroutine prepare_cc3_integrals_R3_abc_cc3
!
!
   module subroutine prepare_cc3_integrals_L3_abc_cc3(wf)
!!
!!    Prepare integral files for L3 amplitudes in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    (db|kc) ordered as dk,bc
!!    (jl|kc) ordered as ljk,c
!!    L_jbkc  ordered as jk,bc
!!
!!    Based on omega_cc3_integrals_cc3 written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
      real(dp), dimension(:,:), allocatable     :: o2_help ! Help array for constructing L_jbkc
!
      integer :: b, c
      type(batching_index) :: batch_c
!
      integer :: req_0, req_c
      integer :: c_batch
!
      batch_c = batching_index(wf%n_v)
!
!     (db|kc) stored as kd,bc
!
      req_0 = wf%integrals%n_J*wf%n_v**2
      req_c = 2*wf%n_o*wf%n_v**2 + wf%integrals%n_J*wf%n_o
!
      call mem%batch_setup(batch_c, req_0, req_c)
!
      wf%g_dbkc_t_v = direct_stream_file('g_dbkc_t_v', wf%n_v*wf%n_o)
      call wf%g_dbkc_t_v%open_('write')
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
         call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_c%length)
         call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_v, batch_c%length)
!
         call wf%get_vvov(g_pqrs,      &
                           1, wf%n_v,  &
                           1, wf%n_v,  &
                           1, wf%n_o,  &
                           batch_c%first, batch_c%last)
!
         call sort_1234_to_3124(g_pqrs, h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_c%length)
!
         call wf%g_dbkc_t_v%write_compound_full_batch(h_pqrs, wf%n_v, batch_c)
!
         call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_c%length)
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_v, batch_c%length)
!
      enddo
!
      call wf%g_dbkc_t_v%close_()
!
!     (jl|kc) stored jkl,c
!
      req_0 = wf%integrals%n_J*wf%n_o**2
      req_c = 2*wf%n_o**3 + wf%integrals%n_J*wf%n_o
!
      call mem%batch_setup(batch_c, req_0, req_c)
!
      wf%g_jlkc_t_v = direct_stream_file('g_jlkc_t_v', wf%n_o**3)
      call wf%g_jlkc_t_v%open_('write')
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
         call mem%alloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
!
         call wf%get_ooov(g_pqrs,      &
                           1, wf%n_o,  &
                           1, wf%n_o,  &
                           1, wf%n_o,  &
                           batch_c%first, batch_c%last)
!
         call sort_1234_to_1324(g_pqrs, h_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
!
         call wf%g_jlkc_t_v%write_interval(h_pqrs, batch_c)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
         call mem%dealloc(h_pqrs, wf%n_o, wf%n_o, wf%n_o, batch_c%length)
!
      enddo
!
      call wf%g_jlkc_t_v%close_()
!
!     L_jbkc ordered jk,bc
!
      call mem%alloc(o2_help, wf%n_o, wf%n_o)
!
      req_0 = wf%integrals%n_J*wf%n_o*wf%n_v
      req_c = 2*wf%n_o**2*wf%n_v + wf%integrals%n_J*wf%n_o
!
      call mem%batch_setup(batch_c, req_0, req_c)
!
      wf%L_jbkc_t_v = direct_stream_file('L_jbkc_t_v', wf%n_o**2)
      call wf%L_jbkc_t_v%open_('write')
!
      call batch_c%determine_limits(1)
      call mem%alloc(h_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_c%length)
!
      do c_batch = 1, batch_c%num_batches
!
         call batch_c%determine_limits(c_batch)
!
         call mem%alloc(g_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_c%length)
!
         call wf%get_ovov(g_pqrs,      &
                           1, wf%n_o,  &
                           1, wf%n_v,  &
                           1, wf%n_o,  &
                           batch_c%first, batch_c%last)
!
         call sort_1234_to_1324(g_pqrs, h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_c%length)
!
         do c = 1, batch_c%length
            do b = 1, wf%n_v
!
               call sort_12_to_21(h_pqrs(:,:,b,c), o2_help, wf%n_o, wf%n_o)
!
               call dscal(wf%n_o**2, two, h_pqrs(:,:,b,c), 1)
!
               call daxpy(wf%n_o**2, -one, o2_help, 1, h_pqrs(:,:,b,c), 1)
!
            enddo
         enddo
!
         call wf%L_jbkc_t_v%write_compound_full_batch(h_pqrs, wf%n_v, batch_c)
!
         call mem%dealloc(g_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_c%length)
!
      enddo
!
      call batch_c%determine_limits(1)
      call mem%dealloc(h_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_c%length)
!
      call mem%dealloc(o2_help, wf%n_o, wf%n_o)
!
      call wf%L_jbkc_t_v%close_()
!
   end subroutine prepare_cc3_integrals_L3_abc_cc3
!
!
   module subroutine omega_cc3_W_calc_abc_cc3(wf, a, b, c,            &
                                              t_ijk, u_ijk, t_ijab,   &
                                              g_ljak, g_ljbk, g_ljck, &
                                              g_bdak, g_cdak, g_cdbk, &
                                              g_adbk, g_adck, g_bdck, &
                                              overwrite)
!!
!!    Omega CC3 intermediate W_ijk for fixed a,b,c
!!    Written by Rolf H. Myhre and Alexander C. Paul July 2019
!!
!!    Contributions to W
!!     W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il(lj|ck))
!!
!!    based on omega_cc3_W_calc written by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_ijk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_ijk
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in)   :: t_ijab
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljak
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljbk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)           :: g_ljck
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_bdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_cdak
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_cdbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_adbk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_adck
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                   :: g_bdck
!
!     If present and false, t_abc is not overwritten by first dgemm
      logical, optional, intent(in) :: overwrite
!
      real(dp) :: alpha
!
      if (.not. present(overwrite)) then
         alpha = zero
      elseif (overwrite) then
         alpha = zero
      else
         alpha = one
      endif
!
!
!     u_ijk terms
!     -----------
!
!     t^db_ij*(ad|ck)
!     ---------------
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,b),  & ! t_ij_d,b
                  wf%n_o**2,        &
                  g_adck,           & ! g_dk,ac
                  wf%n_v,           &
                  alpha,            &
                  t_ijk,            &
                  wf%n_o**2)
!
!     -t^ab_il*(lj|ck)
!     ----------------
!
      call dgemm('N', 'N',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  t_ijab(:,:,a,b),  & ! t_i_l,ab
                  wf%n_o,           &
                  g_ljck,           & ! g_l_jk,c
                  wf%n_o,           &
                  one,              &
                  t_ijk,            &
                  wf%n_o)
!
!     t^dc_jk*(bd|ai)
!     ---------------
!
      call dgemm('T', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_v,           &
                  one,              &
                  g_bdak,           & ! g_d_i,ba
                  wf%n_v,           &
                  t_ijab(:,:,:,c),  & ! t_jk_d,c
                  wf%n_o**2,        &
                  one,              &
                  t_ijk,            &
                  wf%n_o)
!
!     -t^ca_kl*(li|bj)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljbk,           & ! g_l_ij,b
                  wf%n_o,           &
                  t_ijab(:,:,c,a),  & ! t_k_l,ca
                  wf%n_o,           &
                  one,              &
                  t_ijk,            &
                  wf%n_o**2)
!
!
!     u_ikj terms
!     -----------
!
!     t^dc_ik*(ad|bj)
!     ---------------
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,c),  & ! t_ik_d,c
                  wf%n_o**2,        &
                  g_adbk,           & ! g_dj,ab
                  wf%n_v,           &
                  zero,             &
                  u_ijk,            &
                  wf%n_o**2)
!
!     -t^ac_il*(lk|bj)
!     ----------------
!
      call dgemm('N', 'N',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  t_ijab(:,:,a,c),  & ! t_i_l,ac
                  wf%n_o,           &
                  g_ljbk,           & ! g_l_kj,b
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o)
!
!     t^db_kj*(cd|ai)
!     ---------------
!
      call dgemm('T', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_v,           &
                  one,              &
                  g_cdak,           & ! g_d_i,ca
                  wf%n_v,           &
                  t_ijab(:,:,:,b),  & ! t_kj_d,b
                  wf%n_o**2,        &
                  one,              &
                  u_ijk,            &
                  wf%n_o)
!
!     -t^ba_jl*(li|ck)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljck,           & ! g_l_ik,c
                  wf%n_o,           &
                  t_ijab(:,:,b,a),  & ! t_j_l,ba
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
!
      call sort_123_to_132_and_add(u_ijk, t_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!
!     u_jik terms
!     -----------
!
!     t^da_ji*(bd|ck)
!     ---------------
!
!
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,a),  & ! t_ji_d,a
                  wf%n_o**2,        &
                  g_bdck,           & ! g_dk,bc
                  wf%n_v,           &
                  zero,             &
                  u_ijk,            &
                  wf%n_o**2)
!
!     -t^cb_kl*(lj|ai)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljak,           & ! g_l_ji,a
                  wf%n_o,           &
                  t_ijab(:,:,c,b),  & ! t_k_l,cb
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
      call sort_123_to_213_and_add(u_ijk, t_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!
!     u_kij terms
!     -----------
!
!     t^da_ki*(cd|bj)
!     ---------------
!      
      call dgemm('N', 'N',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  t_ijab(:,:,:,a),  & ! t_ki_d,a
                  wf%n_o**2,        &
                  g_cdbk,           & ! g_dj,cb
                  wf%n_v,           &
                  zero,             &
                  u_ijk,            &
                  wf%n_o**2)
!
!     -t^bc_jl*(lk|ai)
!     ----------------
!
      call dgemm('T', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_o,           &
                  -one,             &
                  g_ljak,           & ! g_l_ki,a
                  wf%n_o,           &
                  t_ijab(:,:,b,c),  & ! t_j_l,bc
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
      call sort_123_to_231_and_add(u_ijk, t_ijk, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine omega_cc3_W_calc_abc_cc3
!
!
   module subroutine omega_cc3_eps_abc_cc3(wf, a, b, c, t_ijk, omega)
!!
!!    Omega CC3 epsilon denominator in batches of a,b,c
!!    Written by Alexander C. Paul, July 2019    
!!
!!    Divide W^abc_ijk by -epsilon^abc_ijk to obtain T^abc_ijk
!!    Optional argument omega for jacobian transformations
!!
!!    t^abc_ijk = -W^abc_ijk/epsilon^abc_ijk
!!
!!    based on omega_cc3_eps_cc3 by Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(inout) :: t_ijk
!
      real(dp), optional :: omega
!
      integer i, j, k
!
      real(dp) :: epsilon_abc, epsilon_k, epsilon_kj
!
!
      if (present(omega)) then 
!
         epsilon_abc =  omega - wf%orbital_energies(wf%n_o + a)   &
                        - wf%orbital_energies(wf%n_o + b)         &
                        - wf%orbital_energies(wf%n_o + c)
                     else
!
         epsilon_abc =  - wf%orbital_energies(wf%n_o + a)   &
                        - wf%orbital_energies(wf%n_o + b)   &
                        - wf%orbital_energies(wf%n_o + c)
                     end if 
!
!$omp parallel do schedule(static) private(i)
      do i = 1, wf%n_o
!
         t_ijk(i,i,i) = zero
!
      enddo
!$omp end parallel do
!
!$omp parallel do schedule(static) private(i,j,k,epsilon_k,epsilon_kj)
      do k = 1, wf%n_o
!
         epsilon_k = epsilon_abc + wf%orbital_energies(k)
!
         do j = 1, wf%n_o
!
            epsilon_kj = epsilon_k + wf%orbital_energies(j)
!
            do i = 1, wf%n_o
!
               t_ijk(i,j,k) = t_ijk(i,j,k)*one/(epsilon_kj + wf%orbital_energies(i))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine omega_cc3_eps_abc_cc3
!
!
   module subroutine jacobian_transpose_cc3_c3_calc_abc_cc3(wf, a, b ,c, c_ia, c_ijab, &
                                                            c_ijk, u_ijk, v_ijk, F_kc, &
                                                            L_jakb, L_jakc, L_jbkc,    &
                                                            g_jlka, g_jlkb, g_jlkc,    &
                                                            g_dbka, g_dcka, g_dckb,    &
                                                            g_dakb, g_dakc, g_dbkc)
!!
!!    Calculate the L3 amplitudes for fixed indices a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    C^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + C_abij*F_kc - L_abik*F_jc)
!!    + sum_l (C_ablk g_iljc - C_abil L_jlkc) - sum_d (C_adjk g_ibdc - C_adij L_dbkc)
!!
!!    Contibutions from outer products:
!!    P^abc_ijk (C_ai*L_jbkc - C_ak*L_jbic + C_abij*F_kc - C_abik*F_jc)
!!
!!    Contibutions from matrix multiplication:
!!      sum_l P^abc_ijk (C_ablk g_iljc + C_abil g_jckl - 2 C_abil g_jlkc) 
!!    - sum_d P^abc_ijk (C_adjk g_ibdc + C_adij g_dckb - 2 C_adij g_dbkc)
!!
!!    based on jacobian_transpose_cc3_c3_calc written by A.Paul and Rolf H. Myhre
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: c_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: c_ijab
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: c_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: u_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: v_ijk
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: F_kc
!
!     L_ibjc ordered ij,bc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jakb
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jakc
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: L_jbkc
!
!     g_jlkc ordered jkl,c
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlka
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlkb
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in) :: g_jlkc
!
!     g_dbkc ordered kd,bc
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dbka
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dcka
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dckb
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dakb
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dakc
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: g_dbkc
!
!
!     :: Contribution 1 ::
!
!     c_ijk <- u_ijk,ikj,kji = -sum_l (2* c^ab_il g_jlkc - c^ab_il g_kljc - c^ab_kl g_jlic)
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  c_ijab(:,:,a,b),  & ! c_i_l,ab c_k_l,ab
                  wf%n_o,           &
                  g_jlkc,           & ! g_jk_l,c g_kj_l,c g_ji_l,c
                  wf%n_o**2,        &
                  zero,             &
                  u_ijk,            &
                  wf%n_o)
!
!     c_ijk <- u_ijk,ikj,kji = sum_d (2* c^db_ij g_dakc - c^db_kj g_daic - c^db_ik g_dajc)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  c_ijab(:,:,:,b),  & ! c_ij_d,b c_kj_d,b c_ik_d,b
                  wf%n_o**2,        &
                  g_dakc,           & ! g_k_d,ac g_i_d,ac g_j_d,ac
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
!     c_ijk <- u_ijk,ikj,kji = 2*L_iajb*c^c_k - L_iakb*c^c_j - L_kajb*c^c_i
!
      call dger(wf%n_o**2, &
               wf%n_o,     &
               one,        &
               L_jakb,     & ! L_ij,ab L_ik,ab L_kj,ab
               1,          &
               c_ia(:,c),  & ! c_k,c c_j,c c_i,c
               1,          &
               u_ijk,      &
               wf%n_o**2)
!
!     c_ijk <- u_ijk,ikj,kji = 2*c^ab_ij*F_kc - c^ab_ik*F_jc - c^ab_kj*F_ic
!
      call dger(wf%n_o**2,       &
               wf%n_o,           &
               one,              &
               c_ijab(:,:,a,b),  & ! c_ij,ab c_ik,ab c_kj,ab
               1,                &
               F_kc(:,c),        & ! F_k,c F_j,c F_i,c
               1,                &
               u_ijk,            &
               wf%n_o**2)
!
!     c_ijk <- u_ijk,ikj,kji <- v_jik,jki,kij 
!     = - sum_l (2*c^ba_jl g_ilkc - c^ba_jl g_klic - c^ba_kl g_iljc)
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  c_ijab(:,:,b,a),  & ! c_j_l,ba c_k_l,ba
                  wf%n_o,           &
                  g_jlkc,           & ! g_ik_l,c g_ki_l,c g_ij_l,c
                  wf%n_o**2,        &
                  zero,             &
                  v_ijk,            &
                  wf%n_o)
!
!     c_ijk <- u_ijk,ikj,kji <- v_jik,jki,kij
!     = sum_d (2* c^da_ji g_dbkc - c^da_jk g_dbic - c^da_ki g_dbjc)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  c_ijab(:,:,:,a),  & ! c_ji_d,a c_ki_d,a c_jk_d,a
                  wf%n_o**2,        &
                  g_dbkc,           & ! g_k_d,bc g_j_d,bc g_i_d,bc
                  wf%n_o,           &
                  one,              &
                  v_ijk,            &
                  wf%n_o**2)
!
      call sort_123_to_213_and_add(v_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
      call add_two_123_min_132_min_321(u_ijk, c_ijk, wf%n_o)
!
!
!     :: Contribution 2 ::
!
!     c_ijk <- u_ikj,ijk,jki = -sum_l (2* c^ac_il g_kljb - c^ac_il g_jlkb - c^ac_jl g_klib)
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  c_ijab(:,:,a,c),  & ! c_i_l,ac c_j_l,ac
                  wf%n_o,           &
                  g_jlkb,           & ! g_kj_l,b g_jk_l,b g_ki_l,b
                  wf%n_o**2,        &
                  zero,             &
                  u_ijk,            &
                  wf%n_o)
!
!     c_ijk <- u_ikj,ijk,jki = sum_d (2* c^dc_ik g_dajb - c^dc_jk g_daib - c^dc_ij g_dakb)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  c_ijab(:,:,:,c),  & ! c_ik_d,c c_jk_d,c c_ij_d,c
                  wf%n_o**2,        &
                  g_dakb,           & ! g_j_d,ab g_i_d,ab g_k_d,ab
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
!     c_ijk <- u_ikj,ijk,jki = 2*L_iakc*c^b_j - L_iajc*c^b_k - L_jakc*c^b_i
!
      call dger(wf%n_o**2, &
               wf%n_o,     &
               one,        &
               L_jakc,     & ! L_ik,ac L_ij,ac L_jk,ac
               1,          &
               c_ia(:,b),  & ! c_j,b c_k,b c_i,b
               1,          &
               u_ijk,      &
               wf%n_o**2)
!
!     c_ijk <- u_ikj,ijk,jki = 2*c^ac_ik*F_jb - c^ac_ij*F_kb - c^ac_jk*F_ib
!
      call dger(wf%n_o**2,       &
               wf%n_o,           &
               one,              &
               c_ijab(:,:,a,c),  & ! c_ik,ac c_ij,ac c_jk,ac
               1,                &
               F_kc(:,b),        & ! F_j,b F_k,b F_i,b
               1,                &
               u_ijk,            &
               wf%n_o**2)
!
!     c_ijk <- u_ikj,ijk,jki <- v_kij,kji,jik 
!     = - sum_l (2*c^ca_kl g_iljb - c^ca_kl g_jlib - c^ca_jl g_ilkb)
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  c_ijab(:,:,c,a),  & ! c_k_l,ca c_j_l,ca
                  wf%n_o,           &
                  g_jlkb,           & ! g_ij_l,b g_ji_l,b g_ik_l,b
                  wf%n_o**2,        &
                  zero,             &
                  v_ijk,            &
                  wf%n_o)
!
!     c_ijk <- u_ikj,ijk,jki <- v_kij,kji,jik
!     = sum_d (2* c^da_ki g_dcjb - c^da_ji g_dckb - c^da_kj g_dcib)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  c_ijab(:,:,:,a),  & ! c_ki_d,a c_ji_d,a c_kj_d,a
                  wf%n_o**2,        &
                  g_dckb,           & ! g_j_d,cb g_k_d,cb g_i_d,cb
                  wf%n_o,           &
                  one,              &
                  v_ijk,            &
                  wf%n_o**2)
!
      call sort_123_to_213_and_add(v_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
      call add_two_132_min_231_min_123(u_ijk, c_ijk, wf%n_o)
!
!
!     :: Contribution 3 ::
!
!     c_ijk <- u_jki,jik,ikj = -sum_l (2* c^bc_jl g_klia - c^bc_jl g_ilka - c^bc_il g_klja)
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  c_ijab(:,:,b,c),  & ! c_j_l,bc c_j_l,bc c_i_l,bc
                  wf%n_o,           &
                  g_jlka,           & ! g_ki_l,a g_ik_l,a g_kj_l,a
                  wf%n_o**2,        &
                  zero,             &
                  u_ijk,            &
                  wf%n_o)
!
!     c_ijk <- u_jki,jik,ikj = sum_d (2* c^dc_jk g_dbia - c^dc_ik g_dbja - c^dc_ji g_dbka)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  c_ijab(:,:,:,c),  & ! c_jk_d,c c_ik_d,c c_ji_d,c
                  wf%n_o**2,        &
                  g_dbka,           & ! g_i_d,ba g_j_d,ba g_k_d,ba
                  wf%n_o,           &
                  one,              &
                  u_ijk,            &
                  wf%n_o**2)
!
!     c_ijk <- u_jki,jik,ikj = 2*L_jbkc*c^a_i - L_jbic*c^a_k - L_ibkc*c^a_j
!
      call dger(wf%n_o**2, &
               wf%n_o,     &
               one,        &
               L_jbkc,     & ! L_jk,bc L_ji,bc L_ik,bc
               1,          &
               c_ia(:,a),  & ! c_i,a c_k,a c_j,a
               1,          &
               u_ijk,      &
               wf%n_o**2)
!
!     c_ijk <- u_jki,jik,ikj = 2*c^bc_jk*F_ia - c^bc_ji*F_ka - c^bc_ik*F_ja
!
      call dger(wf%n_o**2,       &
               wf%n_o,           &
               one,              &
               c_ijab(:,:,b,c),  & ! c_jk,bc c_ji,bc c_ik,bc
               1,                &
               F_kc(:,a),        & ! F_i,a F_k,a F_j,a
               1,                &
               u_ijk,            &
               wf%n_o**2)
!
!     c_ijk <- u_jki,jik,ikj <- v_kji,kij,ijk
!     = - sum_l (2*c^cb_kl g_jlia - c^cb_kl g_ilja - c^cb_il g_jlka)
!
      call dgemm('N', 'T',          &
                  wf%n_o,           &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  -one,             &
                  c_ijab(:,:,c,b),  & ! c_k_l,cb c_i_l,cb
                  wf%n_o,           &
                  g_jlka,           & ! g_ji_l,a g_ij_l,a g_jk_l,a
                  wf%n_o**2,        &
                  zero,             &
                  v_ijk,            &
                  wf%n_o)
!
!     c_ijk <- u_jki,jik,ikj <- v_kji,kij,ijk
!     = sum_d (2* c^db_kj g_dcia - c^db_ij g_dcka - c^db_ki g_dcja)
!
      call dgemm('N', 'T',          &
                  wf%n_o**2,        &
                  wf%n_o,           &
                  wf%n_v,           &
                  one,              &
                  c_ijab(:,:,:,b),  & ! c_kj_d,b c_ij_d,b c_ki_d,b
                  wf%n_o**2,        &
                  g_dcka,           & ! g_i_d,ca g_k_d,ca g_j_d,ca
                  wf%n_o,           &
                  one,              &
                  v_ijk,            &
                  wf%n_o**2)
!
      call sort_123_to_213_and_add(v_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
      call add_two_231_min_213_min_132(u_ijk, c_ijk, wf%n_o)
!
   end subroutine jacobian_transpose_cc3_c3_calc_abc_cc3
!
!
   module subroutine get_triples_cvs_projector_abc_cc3(wf, projector_ijk)
!!
!!    Get triples cvs projector for fixed a,b,c
!!    Written by Alexander C. Paul and Rolf H. Myhre, September 2019
!!
!!    Set up projector for cvs for the triples amplitudes 
!!    if we batch over the virtual indices
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: projector_ijk
!
      integer :: i, j, k
!
      call zero_array(projector_ijk, wf%n_o**3)
!
      do k = 1, wf%n_o
         do j = 1, k
            do i = 1, wf%n_core_MOs
!
               projector_ijk(wf%core_MOs(i),j,k) = one
               projector_ijk(j,wf%core_MOs(i),k) = one
               projector_ijk(j,k,wf%core_MOs(i)) = one
!
               if (j .ne. k) then
!
                  projector_ijk(wf%core_MOs(i),k,j) = one
                  projector_ijk(k,wf%core_MOs(i),j) = one
                  projector_ijk(k,j,wf%core_MOs(i)) = one
!
               end if
!
            end do
         end do
      end do
!
   end subroutine get_triples_cvs_projector_abc_cc3
!
!
end submodule batching_abc
