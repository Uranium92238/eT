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
submodule (cc3_class) omega_cc3
!
!!
!!    Omega submodule
!!
!!    Routines to construct
!!
!!    Omega =  < mu| exp(-T) H exp(T) |R >
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_omega_cc3(wf, omega)
!!
!!    Construct omega (CC3)
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Directs the construction of the projection vector < mu| exp(-T) H exp(T) |R >
!!    for the current amplitudes of the object wfn
!!
      use array_utilities, only: scale_diagonal
!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      real(dp), dimension(:,:), allocatable     :: omega1
      real(dp), dimension(:,:,:,:), allocatable :: omega_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj, t_abij, omega_aibj
!
      type(timings), allocatable :: cc3_timer
      type(timings), allocatable :: ccsd_timer
      type(timings), allocatable :: timer 
!
      timer       = timings('Construct omega CC3', pl='normal')
      ccsd_timer  = timings('Omega CC3 (CCSD contribution)', pl='normal')
      cc3_timer   = timings('Omega CC3 (CC3 contribution)', pl='normal')
!
      call timer%turn_on()
!
      call mem%alloc(omega1, wf%n_v, wf%n_o)
      call zero_array(omega1, wf%n_t1)
!
!     Construct CCSD singles contributions
!
      call ccsd_timer%turn_on()
!
      call wf%omega_ccs_a1(omega1)
!  
      call wf%construct_u_aibj()
!
      call wf%omega_doubles_a1(omega1, wf%u_aibj)
      call wf%omega_doubles_b1(omega1, wf%u_aibj)
      call wf%omega_doubles_c1(omega1, wf%u_aibj)
!
!     Construct CCSD doubles contributions
!
      call mem%alloc(omega_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call zero_array(omega_aibj, wf%n_t1**2)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, wf%n_t1)

      call wf%omega_ccsd_c2(omega_aibj, t_aibj)
      call wf%omega_ccsd_d2(omega_aibj, t_aibj)
      call wf%omega_ccsd_e2(omega_aibj, t_aibj)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(t_aibj, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call ccsd_timer%freeze()
!    
      call mem%alloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(omega_aibj, omega_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(omega_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call cc3_timer%turn_on()
      call wf%omega_cc3_a(omega1, omega_abij)
      call cc3_timer%turn_off()
!
      call symmetrize_12_and_34(omega_abij, wf%n_v, wf%n_o)
!
      call ccsd_timer%turn_on()
      call wf%omega_ccsd_a2(omega_abij, t_abij)
      call wf%omega_ccsd_b2(omega_abij, t_abij)
      call scale_diagonal(half, omega_abij, wf%n_v, wf%n_o)
      call ccsd_timer%turn_off()
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call dcopy(wf%n_t1, omega1, 1, omega, 1)
!
      call mem%dealloc(omega1, wf%n_v, wf%n_o)
!
      call packin(omega(wf%n_t1+1 : wf%n_gs_amplitudes), omega_abij, wf%n_v, wf%n_o)
!
      call mem%dealloc(omega_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine construct_omega_cc3
!
!
   module subroutine omega_cc3_a_cc3(wf, omega1, omega2)
!!
!!    CC3 Omega terms
!!    Written by Rolf H. Myhre, January 2019
!!
!!    t_mu3 = -< mu3|{U,T2}|HF > (epsilon_mu3)^-1
!!
!!    omega_mu1 += < mu1|[H,T3]|HF >
!!
!!    omega_mu2 += < mu2|[H,T3]|HF >
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout) :: omega2
!
!     Arrays for triples amplitudes
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
      real(dp), dimension(:,:,:), allocatable :: v_abc
!
!     Unpacked doubles amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck !Transpose the fock matrix ov block
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lkcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_licj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_lick
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lkcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_licj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_lick_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_klic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_kljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_iljc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ilkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_klic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_kljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_iljc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ilkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jlkc_p => null()
!
!     L_kbjc and L_jbkc each other's transpose,
!     but utilising this makes the code more complicated and
!     error prone without any huge advantages
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_kbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_kbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_kbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_kbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3
!
!     Set up required integrals on disk
      call wf%omega_cc3_integrals()
!
      call mem%alloc(t_abij,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2,t_abij,wf%n_v,wf%n_o,wf%n_v,wf%n_o)
!
      req_0 = 3*wf%n_v**3 + wf%n_v*wf%n_o
      req_1 = 2*wf%n_v**3
      req_2 = 2*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3)
!
!     Allocate integral arrays and assign pointers.
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,wf%n_o,wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%alloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%alloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
!
         call mem%alloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
!
      endif
!
!     Resorting of the Fock-Matrix for easier contractions later
      call mem%alloc(F_ov_ck,wf%n_v,wf%n_o)
      call sort_12_to_21(wf%fock_ia,F_ov_ck,wf%n_o,wf%n_v)
!
!     Arrays for the triples amplitudes
      call mem%alloc(t_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%alloc(u_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%alloc(v_abc,wf%n_v,wf%n_v,wf%n_v)
!
!     Remaining integrals
!
      if (batch_i%num_batches .eq. 1) then !no batching
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         call mem%alloc(L_jbic,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%alloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%alloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
!
         call mem%alloc(g_jlic,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%alloc(g_klic,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%alloc(g_kljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%alloc(g_iljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%alloc(g_ilkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%alloc(g_jlkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
!
         call mem%alloc(L_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(L_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(L_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(L_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(L_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%alloc(L_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
!
      endif
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_dbkc_t%open_('read')
      call wf%g_jlkc_t%open_('read')
      call wf%L_jbkc_t%open_('read')
!
      do i_batch = 1,batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%g_bdck_t%read_interval(g_bdci, batch_i)
         call wf%g_dbkc_t%read_interval(g_dbic, batch_i)
!
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
         do j_batch = 1,i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%g_ljck_t%read_compound(g_ljci, batch_j, batch_i)
            call wf%g_jlkc_t%read_compound(g_jlic, batch_j, batch_i)
            call wf%L_jbkc_t%read_compound(L_jbic, batch_j, batch_i)
!
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_jbic_p => L_jbic
!
            if (j_batch .ne. i_batch) then
!
               call wf%g_bdck_t%read_interval(g_bdcj, batch_j)
               call wf%g_dbkc_t%read_interval(g_dbjc, batch_j)
!
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
               call wf%g_ljck_t%read_compound(g_licj, batch_i, batch_j)
               call wf%g_jlkc_t%read_compound(g_iljc, batch_i, batch_j)
               call wf%L_jbkc_t%read_compound(L_ibjc, batch_i, batch_j)
!
               g_licj_p => g_licj
               g_iljc_p => g_iljc
               L_ibjc_p => L_ibjc
!
            else
!
               g_bdcj_p => g_bdci
               g_dbjc_p => g_dbic
!
               g_licj_p => g_ljci
               g_iljc_p => g_jlic
               L_ibjc_p => L_jbic
!
            endif
!
            do k_batch = 1,j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%g_bdck_t%read_interval(g_bdck, batch_k)
                  call wf%g_dbkc_t%read_interval(g_dbkc, batch_k)
!
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
!
                  call wf%g_ljck_t%read_compound(g_lkci, batch_k, batch_i)
                  call wf%g_jlkc_t%read_compound(g_klic, batch_k, batch_i)
                  call wf%L_jbkc_t%read_compound(L_kbic, batch_k, batch_i)
!
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
                  L_kbic_p => L_kbic
!
                  call wf%g_ljck_t%read_compound(g_lick, batch_i, batch_k)
                  call wf%g_jlkc_t%read_compound(g_ilkc, batch_i, batch_k)
                  call wf%L_jbkc_t%read_compound(L_ibkc, batch_i, batch_k)
!
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  call wf%g_ljck_t%read_compound(g_lkcj, batch_k, batch_j)
                  call wf%g_jlkc_t%read_compound(g_kljc, batch_k, batch_j)
                  call wf%L_jbkc_t%read_compound(L_kbjc, batch_k, batch_j)
!
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  call wf%g_ljck_t%read_compound(g_ljck, batch_j, batch_k)
                  call wf%g_jlkc_t%read_compound(g_jlkc, batch_j, batch_k)
                  call wf%L_jbkc_t%read_compound(L_jbkc, batch_j, batch_k)
!
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
!
               else if (k_batch .eq. i_batch) then ! k_batch = j_batch = i_batch
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
                  L_kbic_p => L_jbic
!
                  g_lick_p => g_ljci
                  g_ilkc_p => g_jlic
                  L_ibkc_p => L_jbic
!
                  g_lkcj_p => g_ljci
                  g_kljc_p => g_jlic
                  L_kbjc_p => L_jbic
!
                  g_ljck_p => g_ljci
                  g_jlkc_p => g_jlic
                  L_jbkc_p => L_jbic
!
               else !k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
                  L_kbic_p => L_jbic
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
!
                  call wf%g_ljck_t%read_compound(g_lkcj, batch_k, batch_j)
                  call wf%g_jlkc_t%read_compound(g_kljc, batch_k, batch_j)
                  call wf%L_jbkc_t%read_compound(L_kbjc, batch_k, batch_j)
!
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
                  L_kbjc_p => L_kbjc
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_kbjc
!
               endif
!
               do i = batch_i%first, batch_i%last
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%last, i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%last, j)
!
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
                        k_rel = k - batch_k%first + 1
!
                        call wf%construct_W(i, j, k, t_abc, u_abc, t_abij, &
                                            g_bdci_p(:,:,:,i_rel),         &
                                            g_bdcj_p(:,:,:,j_rel),         &
                                            g_bdck_p(:,:,:,k_rel),         &
                                            g_ljci_p(:,:,j_rel,i_rel),     &
                                            g_lkci_p(:,:,k_rel,i_rel),     &
                                            g_lkcj_p(:,:,k_rel,j_rel),     &
                                            g_licj_p(:,:,i_rel,j_rel),     &
                                            g_lick_p(:,:,i_rel,k_rel),     &
                                            g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, t_abc)
!
                        call wf%omega_cc3_a_n7(i, j, k, t_abc, u_abc, v_abc, omega2, &
                                               g_dbic_p(:,:,:,i_rel),                &
                                               g_dbjc_p(:,:,:,j_rel),                &
                                               g_dbkc_p(:,:,:,k_rel),                &
                                               g_jlic_p(:,:,j_rel,i_rel),            &
                                               g_klic_p(:,:,k_rel,i_rel),            &
                                               g_kljc_p(:,:,k_rel,j_rel),            &
                                               g_iljc_p(:,:,i_rel,j_rel),            &
                                               g_ilkc_p(:,:,i_rel,k_rel),            &
                                               g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_a_n6(i, j, k, t_abc, u_abc,      &
                                               omega1, omega2, F_ov_ck,    &
                                               L_jbic_p(:,:,j_rel,i_rel),  &
                                               L_kbic_p(:,:,k_rel,i_rel),  &
                                               L_kbjc_p(:,:,k_rel,j_rel),  &
                                               L_ibjc_p(:,:,i_rel,j_rel),  &
                                               L_ibkc_p(:,:,i_rel,k_rel),  &
                                               L_jbkc_p(:,:,j_rel,k_rel))
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!     Close files
!
      call wf%g_bdck_t%close_()
      call wf%g_ljck_t%close_()
      call wf%g_dbkc_t%close_()
      call wf%g_jlkc_t%close_()
      call wf%L_jbkc_t%close_()
!
!     Deallocate the integral arrays
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%dealloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,wf%n_o)
         call mem%dealloc(g_ljci,wf%n_o,wf%n_v,wf%n_o,wf%n_o)
         call mem%dealloc(g_jlic,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
         call mem%dealloc(L_jbic,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
      else
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%dealloc(g_bdcj,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%dealloc(g_bdck,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
!
         call mem%dealloc(g_dbic,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%dealloc(g_dbjc,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
         call mem%dealloc(g_dbkc,wf%n_v,wf%n_v,wf%n_v,batch_i%length)
!
         call mem%dealloc(g_ljci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_lkci,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_lkcj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_licj,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_lick,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(g_ljck,wf%n_o,wf%n_v,batch_i%length,batch_i%length)
!
         call mem%dealloc(g_jlic,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%dealloc(g_klic,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%dealloc(g_kljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%dealloc(g_iljc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%dealloc(g_ilkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
         call mem%dealloc(g_jlkc,wf%n_v,wf%n_o,batch_i%length,batch_i%length)
!
         call mem%dealloc(L_jbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(L_kbic,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(L_kbjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(L_ibjc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(L_ibkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
         call mem%dealloc(L_jbkc,wf%n_v,wf%n_v,batch_i%length,batch_i%length)
!
      endif
!
      call mem%dealloc(F_ov_ck,wf%n_v,wf%n_o)
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%dealloc(u_abc,wf%n_v,wf%n_v,wf%n_v)
      call mem%dealloc(v_abc,wf%n_v,wf%n_v,wf%n_v)
!
      call mem%dealloc(t_abij,wf%n_v,wf%n_v,wf%n_o,wf%n_o)
!
   end subroutine omega_cc3_a_cc3
!
!
   module subroutine omega_cc3_integrals_cc3(wf)
!!
!!    Construct integrals need in CC3 Omega and store on disk
!!    (bd|ck) ordered as dbc,k
!!    (db|kc) ordered as bcd,k
!!    (lj|ck) ordered as lc,jk
!!    (jl|kc) ordered as cl,jk
!!    (jb|kc) stored as L_jbkc = 2g_jbkc - g_jckb ordered as bc,jk
!!
!!    Written by Rolf H. Myhre, January 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_pqrs ! Array for constructed integrals
      real(dp), dimension(:,:,:,:), allocatable :: h_pqrs ! Array for sorted integrals
!
      type(batching_index) :: batch_k
!
      integer :: req_0, req_k
      integer :: current_k_batch
!
      batch_k = batching_index(wf%n_o)
!
!     (bd|ck)
!
      req_0 = 0
      req_k = wf%n_v**3
      call wf%eri%get_eri_t1_mem('vvvo', req_0, req_k, wf%n_v, wf%n_v, wf%n_v, 1, qp=.true.)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      call mem%alloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
!
      wf%g_bdck_t = direct_stream_file('g_bdck_t',wf%n_v**3)
      call wf%g_bdck_t%open_('write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call wf%eri%get_eri_t1('vvvo', g_pqrs, &
                                first_s=batch_k%first, last_s=batch_k%last, qp=.true.)
!
         call wf%g_bdck_t%write_interval(g_pqrs, batch_k)
!
      enddo
!
      call wf%g_bdck_t%close_()
      call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
!
!     (db|kc)
!
      req_0 = 0
      req_k = wf%n_v**3
      call wf%eri%get_eri_t1_mem('vvov', req_0, req_k, wf%n_v, wf%n_v, 1, wf%n_v)
!
      call mem%alloc(g_pqrs, wf%n_v, wf%n_v, batch_k%max_length, wf%n_v)
      call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
!
      wf%g_dbkc_t = direct_stream_file('g_dbkc_t',wf%n_v**3)
      call wf%g_dbkc_t%open_('write')
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call wf%eri%get_eri_t1('vvov', g_pqrs, first_r=batch_k%first, last_r=batch_k%last)
!
         call sort_1234_to_2413(g_pqrs,h_pqrs,wf%n_v,wf%n_v,batch_k%length,wf%n_v)
!
         call wf%g_dbkc_t%write_interval(h_pqrs, batch_k)
!
      enddo
!
      call wf%g_dbkc_t%close_()
!
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_v, batch_k%max_length)
      call mem%dealloc(g_pqrs, wf%n_v, wf%n_v, batch_k%max_length, wf%n_v)
!
!     (lj|ck)
!
      req_0 = 0
      req_k = 2*wf%n_o**2*wf%n_v
      call wf%eri%get_eri_t1_mem('oovo', req_0, req_k, wf%n_o, wf%n_o, wf%n_v, 1)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      wf%g_ljck_t = direct_stream_file('g_ljck_t',wf%n_v*wf%n_o)
      call wf%g_ljck_t%open_('write')
!
      call mem%alloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%max_length)
      call mem%alloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o ,batch_k%max_length)
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call wf%eri%get_eri_t1('oovo', g_pqrs, first_s=batch_k%first,last_s=batch_k%last)
!
         call sort_1234_to_1324(g_pqrs, h_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%length)
!
         call wf%g_ljck_t%write_compound_full_batch(h_pqrs, wf%n_o, batch_k)
!
      enddo
!
      call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, wf%n_v, batch_k%max_length)
      call mem%dealloc(h_pqrs, wf%n_o, wf%n_v, wf%n_o, batch_k%max_length)
!
      call wf%g_ljck_t%close_()
!
!     (jl|kc)
!
      req_0 = 0
      req_k = 2*wf%n_o**2*wf%n_v
      call wf%eri%get_eri_t1_mem('ooov', req_0, req_k, wf%n_o, wf%n_o, 1, wf%n_v)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      wf%g_jlkc_t = direct_stream_file('g_jlkc_t',wf%n_v*wf%n_o)
      call wf%g_jlkc_t%open_('write')
!
      call mem%alloc(g_pqrs, wf%n_o, wf%n_o, batch_k%max_length, wf%n_v)
      call mem%alloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_k%max_length)
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call wf%eri%get_eri_t1('ooov', g_pqrs, first_r=batch_k%first, last_r=batch_k%last)
!
         call sort_1234_to_4213(g_pqrs, h_pqrs, wf%n_o, wf%n_o, batch_k%length, wf%n_v)
!
         call wf%g_jlkc_t%write_compound_full_batch(h_pqrs, wf%n_o, batch_k)
!
      enddo
!
      call mem%dealloc(g_pqrs, wf%n_o, wf%n_o, batch_k%max_length, wf%n_v)
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_o, wf%n_o, batch_k%max_length)
!
      call wf%g_jlkc_t%close_()
!
!     (jb|kc)
!
      req_0 = 0
      req_k = 2*wf%n_v**2*wf%n_o
      call wf%eri%get_eri_t1_mem('ovov', req_0, req_k, wf%n_o, wf%n_v, 1, wf%n_v)
!
      call mem%batch_setup(batch_k,req_0,req_k)
!
      wf%L_jbkc_t = direct_stream_file('L_jbkc_t',wf%n_v**2)
      call wf%L_jbkc_t%open_('write')
!
      call mem%alloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%max_length)
      call mem%alloc(g_pqrs, wf%n_o, wf%n_v, batch_k%max_length, wf%n_v)
!
      do current_k_batch = 1,batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call wf%eri%get_eri_t1('ovov', g_pqrs, first_r=batch_k%first, last_r=batch_k%last)
!
         call sort_1234_to_1432(g_pqrs, h_pqrs, wf%n_o, wf%n_v, batch_k%length, wf%n_v)
         call dscal(wf%n_v**2*wf%n_o*batch_k%length, two, g_pqrs, 1)
         call daxpy(wf%n_v**2*wf%n_o*batch_k%length, -one, h_pqrs, 1, g_pqrs, 1)
         call sort_1234_to_2413(g_pqrs, h_pqrs, wf%n_o, wf%n_v, batch_k%length, wf%n_v)
!
         call wf%L_jbkc_t%write_compound_full_batch(h_pqrs, wf%n_o, batch_k)
!
      enddo
!
      call mem%dealloc(g_pqrs, wf%n_o, wf%n_v, batch_k%max_length, wf%n_v)
      call mem%dealloc(h_pqrs, wf%n_v, wf%n_v, wf%n_o, batch_k%max_length)
!
      call wf%L_jbkc_t%close_()
!
   end subroutine omega_cc3_integrals_cc3
!
!
   module subroutine omega_cc3_a_n6_cc3(wf, i, j, k, t_abc, u_abc, &
                                        omega1, omega2, F_ov_ck,   &
                                        L_jbic, L_kbic, L_kbjc,    &
                                        L_ibjc, L_ibkc, L_jbkc)
!!
!!    omega cc3 n6 terms
!!
!!    Calculate the triples contribution to omega1 and
!!    the Fock contribution to omega2 scaling as n^6
!!
!!    Written by Rolf H. Myhre, January 2019
!!
!!    omega^a_i += sum_bcjk (t^abc_ijk - t^cba_ijk) L_jbkc
!!    
!!    omega^ab_ij += P^{ab}_{ij}sum_ck (t^abc_ijk - t^cba_ijk) F_kc
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout)                   :: omega1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: F_ov_ck
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbic
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_kbjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibjc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_ibkc
      real(dp), dimension(wf%n_v, wf%n_v), intent(in)                      :: L_jbkc
!
!     Construct u_abc = t_abc - t_cba
!     Zero if i == k, but this is never true
!
      call construct_123_minus_321(t_abc, u_abc, wf%n_v)
!
!     omega_ai += sum_bc (t^abc - t^cba) L_jbkc
!
      call dgemv('N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 one, &
                 u_abc, &
                 wf%n_v, &
                 L_jbkc, &
                 1, &
                 one, &
                 omega1(:,i), &
                 1)
!
!     omega_abij += sum_c (t^abc - t^cba) F_kc
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 one, &
                 u_abc, &
                 wf%n_v**2, &
                 F_ov_ck(:,k), &
                 1, &
                 one, &
                 omega2(:,:,i,j), &
                 1)
!
!     omega_ak += sum_cb (t^cba - t^abc) L_jbic
!
      call dgemv('N', &
                 wf%n_v, &
                 wf%n_v**2, &
                 -one, &
                 u_abc, &
                 wf%n_v, &
                 L_jbic, &
                 1, &
                 one, &
                 omega1(:,k), &
                 1)
!
!     omega_abkj += sum_c (t^cba - t^abc) F_ic
!
      call dgemv('N', &
                 wf%n_v**2, &
                 wf%n_v, &
                 -one, &
                 u_abc, &
                 wf%n_v**2, &
                 F_ov_ck(:,i), &
                 1, &
                 one, &
                 omega2(:,:,k,j), &
                 1)
!
!
!
      if (i .ne. j .and. j .ne. k) then
!
!        Construct u_abc = t_acb - t_cab
!
         call construct_132_minus_312(t_abc, u_abc, wf%n_v)
!
!        omega_ai += sum_cb (t^acb - t^cab) L_kbjc
!
         call dgemv('N', &
                    wf%n_v, &
                    wf%n_v**2, &
                    one, &
                    u_abc, &
                    wf%n_v, &
                    L_kbjc, &
                    1, &
                    one, &
                    omega1(:,i), &
                    1)
!
!        omega_abik += sum_c (t^acb - t^cab) F_jc
!
         call dgemv('N', &
                    wf%n_v**2, &
                    wf%n_v, &
                    one, &
                    u_abc, &
                    wf%n_v**2, &
                    F_ov_ck(:,j), &
                    1, &
                    one, &
                    omega2(:,:,i,k), &
                    1)
!
!           omega_aj += sum_cb (t^cab - t^acb) L_kbic
!
            call dgemv('N', &
                       wf%n_v, &
                       wf%n_v**2, &
                       -one, &
                       u_abc, &
                       wf%n_v, &
                       L_kbic, &
                       1, &
                       one, &
                       omega1(:,j), &
                       1)
!
!           omega_abjk += sum_c (t^cab - t^acb) F_ic
!
            call dgemv('N', &
                       wf%n_v**2, &
                       wf%n_v, &
                       -one, &
                       u_abc, &
                       wf%n_v**2, &
                       F_ov_ck(:,i), &
                       1, &
                       one, &
                       omega2(:,:,j,k), &
                       1)
!
!        Construct u_abc = t_bac - t_bca
!        This is zero if j == k
!
         call construct_213_minus_231(t_abc, u_abc, wf%n_v)
!
!        omega_aj += sum_cb (t^bac - t^bca) L_ibkc
!
         call dgemv('N', &
                    wf%n_v, &
                    wf%n_v**2, &
                    one, &
                    u_abc, &
                    wf%n_v, &
                    L_ibkc, &
                    1, &
                    one, &
                    omega1(:,j), &
                    1)
!
!        omega_abji += sum_c (t^bac - t^bca) F_kc
!
         call dgemv('N', &
                    wf%n_v**2, &
                    wf%n_v, &
                    one, &
                    u_abc, &
                    wf%n_v**2, &
                    F_ov_ck(:,k), &
                    1, &
                    one, &
                    omega2(:,:,j,i), &
                    1)
!
!        omega_ak += sum_cb (t^bca - t^bac)*L_ibjc
!
         call dgemv('N', &
                    wf%n_v, &
                    wf%n_v**2, &
                    -one, &
                    u_abc, &
                    wf%n_v, &
                    L_ibjc, &
                    1, &
                    one, &
                    omega1(:,k), &
                    1)
!
!        omega_abki += sum_c (t^bca - t^bac) F_jc
!
         call dgemv('N', &
                    wf%n_v**2, &
                    wf%n_v, &
                    -one, &
                    u_abc, &
                    wf%n_v**2, &
                    F_ov_ck(:,j), &
                    1, &
                    one, &
                    omega2(:,:,k,i), &
                    1)
!
!
      end if
!
   end subroutine omega_cc3_a_n6_cc3
!
!
   module subroutine omega_cc3_a_n7_cc3(wf, i, j, k, t_abc, u_abc, v_abc, omega2, &
                                        g_dbic, g_dbjc, g_dbkc, &
                                        g_jlic, g_klic, g_kljc, g_iljc, g_ilkc, g_jlkc)
!!
!!    omega cc3 n7 terms
!!
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Calculate the triples contribution to omega2. Scaling as n^7
!!
!!    omega_ablj -= P^ab_lj sum_cik(2t^abc_ijk - t^cba_ijk - t^acb_ijk) g_ilkc
!!
!!    omega_adij += P^ad_ij sum_bck(2t^abc_ijk - t^cba_ijk - t^acb_ijk) g_dbkc
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)             :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(inout)   :: omega2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbic
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbjc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)              :: g_dbkc
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_jlic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_klic
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_kljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_iljc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_ilkc
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                      :: g_jlkc
!
!     construct u_abc = 2t_abc - t_acb - t_cba
!
      call construct_123_min_132_min_321(t_abc, u_abc, wf%n_v)
!
!     omega_adij += sum_bc (2t_abc - t_acb - t_cba) g_dbkc
!
      call dgemm('N','N', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v**2, &
                 one, &
                 u_abc, &
                 wf%n_v, &
                 g_dbkc, &
                 wf%n_v**2, &
                 one, &
                 omega2(:,:,i,j), &
                 wf%n_v)
!
!     omega_ablj += sum_c (2t_abc - t_acb - t_cba) g_ilkc
!
      call dgemm('N','N', &
                 wf%n_v**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 u_abc, &
                 wf%n_v**2, &
                 g_ilkc, &
                 wf%n_v, &
                 one, &
                 omega2(:,:,:,j), &
                 wf%n_v**2)
!
!
      if (j .ne. i) then
!
!        resort to u_bac = 2t_bac - t_bca - t_cab
!      
         call sort_123_to_213(u_abc,v_abc,wf%n_v,wf%n_v,wf%n_v)
!      
!        omega_adji += sum_bc (2t_bac - t_bca - t_cab) g_dbkc
!      
         call dgemm('N','N', &
                    wf%n_v, &
                    wf%n_v, &
                    wf%n_v**2, &
                    one, &
                    v_abc, &
                    wf%n_v, &
                    g_dbkc, &
                    wf%n_v**2, &
                    one, &
                    omega2(:,:,j,i), &
                    wf%n_v)
!      
!        omega_abli += sum_c (2t_bac - t_bca - t_cab) g_jlkc
!      
         call dgemm('N','N', &
                    wf%n_v**2, &
                    wf%n_o, &
                    wf%n_v, &
                    -one, &
                    v_abc, &
                    wf%n_v**2, &
                    g_jlkc, &
                    wf%n_v, &
                    one, &
                    omega2(:,:,:,i), &
                    wf%n_v**2)
!
      end if ! j .ne. i
!
!     construct u_cba = 2t_cba - t_abc - t_bca
!
      call construct_321_min_231_min_123(t_abc, u_abc, wf%n_v)
!
!     omega_adkj += sum_bc (2t_cba - t_abc - t_bca) g_dbic
!
      call dgemm('N','N', &
                 wf%n_v, &
                 wf%n_v, &
                 wf%n_v**2, &
                 one, &
                 u_abc, &
                 wf%n_v, &
                 g_dbic, &
                 wf%n_v**2, &
                 one, &
                 omega2(:,:,k,j), &
                 wf%n_v)
!
!     omega_ablj += sum_c (2t_cba - t_abc - t_bca) g_klic
!
      call dgemm('N','N', &
                 wf%n_v**2, &
                 wf%n_o, &
                 wf%n_v, &
                 -one, &
                 u_abc, &
                 wf%n_v**2, &
                 g_klic, &
                 wf%n_v, &
                 one, &
                 omega2(:,:,:,j), &
                 wf%n_v**2)
!
      if (j .ne. k) then
!
!        resort to u_cab = 2t_cab - t_acb - t_bac
!
         call sort_123_to_213(u_abc,v_abc,wf%n_v,wf%n_v,wf%n_v)
!
!        omega_adjk += sum_bc (2t_cab - t_acb - t_bac) g_dbic
!
         call dgemm('N','N', &
                    wf%n_v, &
                    wf%n_v, &
                    wf%n_v**2, &
                    one, &
                    v_abc, &
                    wf%n_v, &
                    g_dbic, &
                    wf%n_v**2, &
                    one, &
                    omega2(:,:,j,k), &
                    wf%n_v)
!
!        omega_ablk += sum_c (2t_cab - t_acb - t_bac) g_jlic
!
         call dgemm('N','N', &
                    wf%n_v**2, &
                    wf%n_o, &
                    wf%n_v, &
                    -one, &
                    v_abc, &
                    wf%n_v**2, &
                    g_jlic, &
                    wf%n_v, &
                    one, &
                    omega2(:,:,:,k), &
                    wf%n_v**2)
!
         if (j .ne. i) then
!
!        construct u_acb = 2t_acb - t_abc - t_cab
!
            call construct_132_min_123_min_312(t_abc, u_abc, wf%n_v)
!
!           omega_adik += sum_bc (2t_acb - t_abc - t_cab) g_dbjc
!
            call dgemm('N','N', &
                       wf%n_v, &
                       wf%n_v, &
                       wf%n_v**2, &
                       one, &
                       u_abc, &
                       wf%n_v, &
                       g_dbjc, &
                       wf%n_v**2, &
                       one, &
                       omega2(:,:,i,k), &
                       wf%n_v)
!
!              omega_ablk += sum_c (2t_acb - t_abc - t_cab) g_iljc
!
               call dgemm('N','N', &
                          wf%n_v**2, &
                          wf%n_o, &
                          wf%n_v, &
                          -one, &
                          u_abc, &
                          wf%n_v**2, &
                          g_iljc, &
                          wf%n_v, &
                          one, &
                          omega2(:,:,:,k), &
                          wf%n_v**2)
!
!           resort to u_abc = 2t_bca - t_bac - t_cba
!
            call sort_123_to_213(u_abc,v_abc,wf%n_v,wf%n_v,wf%n_v)
!
!           omega_adki += sum_bc (2t_bca - t_bac - t_cba) g_dbjc
!
            call dgemm('N','N', &
                       wf%n_v, &
                       wf%n_v, &
                       wf%n_v**2, &
                       one, &
                       v_abc, &
                       wf%n_v, &
                       g_dbjc, &
                       wf%n_v**2, &
                       one, &
                       omega2(:,:,k,i), &
                       wf%n_v)
!
!           omega_abli += sum_c (2t_bca - t_bac - t_cba) g_kljc
!
            call dgemm('N','N', &
                       wf%n_v**2, &
                       wf%n_o, &
                       wf%n_v, &
                       -one, &
                       v_abc, &
                       wf%n_v**2, &
                       g_kljc, &
                       wf%n_v, &
                       one, &
                       omega2(:,:,:,i), &
                       wf%n_v**2)
!
         end if ! j .ne. i
      end if ! k .ne. j
!
   end subroutine omega_cc3_a_n7_cc3
!
!
end submodule
