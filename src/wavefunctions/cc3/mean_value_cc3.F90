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
submodule (cc3_class) mean_value_cc3
!
!!
!!    Mean-value submodule
!!
!!    Contains routines related to the mean values, i.e.
!!    the construction of density matrices as well as expectation
!!    value calculation.
!!
!!    The ground state density is constructed as follows:
!!
!!          D_pq = < Lambda| E_pq |CC >
!!    where:
!!          < Lambda| = < HF| + sum_mu tbar_mu < mu| exp(-T)
!!
!!
!!    In general a CC density matrix can be written as:
!!
!!          D_pq = < X| e^(-T) E_pq e^T |Y >
!!
!!    where X and Y are left and right vectors with contributions from
!!    a reference determinant and excited determinants (< mu|, |nu >):
!!
!!          D_pq =             X_ref < HF| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_mu  < mu| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_ref < HF| e^(-T) E_pq e^T |mu >  Y_mu
!!                 + sum_mu,nu X_mu  < mu| e^(-T) E_pq e^T |nu >  Y_nu
!!
!!    Depending on the type of density matrix (Ground state, transition ,
!!    excited state, interstate transition) different states and thus different
!!    amplitudes X_ref, X_mu, Y_ref and Y_mu will contribute.
!!
!!    In EOM theory the states can be written as the following vectors:
!!
!!          |CC >     = R_0 = (1, 0)
!!          |Lambda > = L_0 = (1, tbar_mu)
!!          |R_k >    = R_k = (-sum_mu(tbar_mu*R_mu), R_mu)
!!          |L_k >    = L_k = (0, L_mu)
!!
!!    The routine names derive from the contribution of the vectors:
!!
!!       ref_ref: first component of the vector for the left and right state
!!
!!       mu_ref:  second component of the vector for the left and
!!                first component of the vector for the right state
!!
!!       ref_mu:  first component of the vector for the left and
!!                second component of the vector for the right state
!!
!!       mu_nu:   second component of the vector for the left and right state
!!
!
   implicit none
!
!
contains
!
!
   module subroutine mu_ref_density_terms_cc3(wf, density, state, L)
!!
!!    Density mu ref terms
!!    Written by Alexander C. Paul, May 2021
!!
!!    Constructs terms of the form:
!!       sum_mu L_mu < mu| E_pq |HF >
!!
!!    corresponding to terms of the ground state density
!!    and the left transition density.
!!
!!    Directs the construction of the CC3 part of the GS and
!!    left transition density, as the equations are the same, but the
!!    required information is different for the left GS and the left ES.
!!
      use array_utilities, only: add_to_subblock, zero_array
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
!
      integer, intent(in) :: state
!
!     L might only be contiguous in the ranges (1:n_t1) and (1+t1: n_t1+n_t2)
!     as L can also be a combined array of t1bar + t2bar
      real(dp), dimension(wf%n_t1+wf%n_t2), intent(in) :: L
!
      real(dp), dimension(:,:), allocatable :: d_ov, d_vo
!
      logical  :: cvs, rm_core, ground_state
      real(dp) :: omega
!
      type(range_) :: o_range, v_range
!
      type(timings)     :: timer
      character(len=40) :: timer_name
!
      write(timer_name, '(a,i0,a)') 'CC3 contribution to <', state, '|E_pq|0>'
      timer = timings(trim(timer_name), pl='v')
!
      call zero_array(density, wf%n_mo**2)
      call wf%ccsd%mu_ref_density_terms(density, state, L)
!
      call timer%turn_on() ! CC3 contribution
!
      call mem%alloc(d_ov, wf%n_o, wf%n_v)
      call mem%alloc(d_vo, wf%n_v, wf%n_o)
!
      o_range = range_(1, wf%n_o)
      v_range = range_(wf%n_o+1, wf%n_v)
!
      if (state == 0) then
!
         omega = zero
         cvs = .false.
         rm_core = .false.
         ground_state = .true.
!
         call wf%density_cc3_mu_ref_blocks(wf%GS_cc3_density_oo, d_ov, d_vo, &
                                           wf%GS_cc3_density_vv, L, omega,   &
                                           cvs, rm_core, ground_state)
!
         call add_to_subblock(one, wf%GS_cc3_density_oo, density, o_range, o_range)
         call add_to_subblock(one, wf%GS_cc3_density_vv, density, v_range, v_range)
!
      else
!
         omega = wf%left_excitation_energies(state)
         ground_state = .false.
!
         call wf%density_cc3_mu_ref_blocks(wf%L_cc3_density_oo, d_ov, d_vo, &
                                           wf%L_cc3_density_vv, L, omega,   &
                                           wf%cvs, wf%rm_core, ground_state)
!
         call add_to_subblock(one, wf%L_cc3_density_oo, density, o_range, o_range)
         call add_to_subblock(one, wf%L_cc3_density_vv, density, v_range, v_range)
!
      end if
!
      call add_to_subblock(one, d_ov, density, o_range, v_range)
      call add_to_subblock(one, d_vo, density, v_range, o_range)
!
      call mem%dealloc(d_ov, wf%n_o, wf%n_v)
      call mem%dealloc(d_vo, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine mu_ref_density_terms_cc3
!
!
   module subroutine density_cc3_mu_ref_blocks_cc3(wf, d_oo, d_ov, d_vo, d_vv, &
                                                   L, omega, cvs, rm_core, gs)
!!
!!    Density CC3 mu ref blocks
!!    Written by Alexander C. Paul, May 2021
!!
!!    Constructs terms of the ground state density
!!    and the left transition density in blocks.
!!
!!       sum_mu L_mu < mu| E_pq |HF >
!!
      use reordering, only: squareup_and_sort_1234_to_1324
      use reordering, only: sort_12_to_21, sort_1234_to_3412
      use reordering, only: construct_covariant_1324
      use array_utilities, only: zero_array
!
      implicit none
!
      class(cc3), intent(in) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: d_oo
      real(dp), dimension(wf%n_o, wf%n_v), intent(out) :: d_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(out) :: d_vo
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: d_vv
!
      real(dp), intent(in) :: omega
      logical, intent(in)  :: cvs, rm_core, gs
!
      real(dp), dimension(wf%n_t1 + wf%n_t2), intent(in) :: L
!
      real(dp), dimension(:,:), allocatable :: L_ia
      real(dp), dimension(:,:,:,:), allocatable :: L_abij, t_abij
      real(dp), dimension(:,:,:,:), allocatable :: L_ijab, t_ijab
!
      call zero_array(d_oo, wf%n_o**2)
      call zero_array(d_ov, wf%n_o*wf%n_v)
      call zero_array(d_vo, wf%n_v*wf%n_o)
      call zero_array(d_vv, wf%n_v**2)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Construct covariant L2 = 1/3 (2L^ab_ij + L^ba_ij)
      call mem%alloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call construct_covariant_1324(L(wf%n_t1+1:), L_abij, wf%n_v, wf%n_o)
!
!     CC3 contribution in batches of i,j,k
!
      call wf%density_cc3_mu_ref_ijk(d_ov, d_vv, omega, L(1:wf%n_t1), L_abij, t_abij, &
                                     cvs, rm_core, keep_Y=gs)
!
!     Need t2, L2 in ijab ordering for the d_oo term
      call mem%alloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(t_abij, t_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(L_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(L_abij, L_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(L_ia, wf%n_o, wf%n_v)
      call sort_12_to_21(L(1:wf%n_t1), L_ia, wf%n_v, wf%n_o)
!
!     CC3 contribution in batches of a,b,c
!
      call wf%density_cc3_mu_ref_abc(d_oo, omega, L_ia, L_ijab, t_ijab, cvs, rm_core)
!
      call mem%dealloc(L_ia, wf%n_o, wf%n_v)
      call mem%dealloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(L_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
   end subroutine density_cc3_mu_ref_blocks_cc3
!
!
   module subroutine density_cc3_mu_ref_abc_cc3(wf, density_oo, omega, tbar_ia, &
                                                tbar_ijab, t_ijab, cvs, rm_core)
!!
!!    One electron density excited-determinant/reference term
!!    in batches of the virtual orbitals a,b,c
!!    Written by Alexander C. Paul, July 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    Construct t^abc_ijk and tbar^abc_ijk in batches of a,b,c.
!!    Calls routines calculating the actual density contributions
!!    to the oo-block
!!
!!    t_mu3 = -< mu3|{U,T2}|HF > (eps_mu3)^-1
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1| [H,tau_nu3] |R >
!!                             + tbar_mu2 < mu2| [H,tau_nu3] |R >
!!
!!    oo-block:
!!          D_kl += -1/2 sum_ij,abc t^abc_ijk tbar^abc_ijl
!!
      use omp_lib
      use array_utilities, only: copy_and_scale, zero_array
      use reordering, only: construct_contravariant_t3
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), intent(in) :: omega
!
      logical, intent(in)  :: cvs, rm_core
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: tbar_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: tbar_ijab
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: t_ijab
!
      real(dp), dimension(:,:,:,:), allocatable :: t_ijk
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ijk
      real(dp), dimension(:,:,:,:), allocatable :: u_ijk
      real(dp), dimension(:,:,:,:), allocatable :: v_ijk
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
!     Integrals for the construction of the T3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljak
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljbk
      real(dp), dimension(:,:,:,:), allocatable, target :: g_ljck
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljak_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljbk_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_ljck_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target :: g_bdak
      real(dp), dimension(:,:,:,:), allocatable, target :: g_cdak
      real(dp), dimension(:,:,:,:), allocatable, target :: g_cdbk
      real(dp), dimension(:,:,:,:), allocatable, target :: g_adbk
      real(dp), dimension(:,:,:,:), allocatable, target :: g_adck
      real(dp), dimension(:,:,:,:), allocatable, target :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_bdak_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_cdak_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_cdbk_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_adbk_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_adck_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_bdck_p => null()
!
!     Integrals for the construction of the multipliers
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jlka
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jlkb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_jlkc
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jlka_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jlkb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_jlkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dbka
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dcka
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dckb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dakb
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dakc
      real(dp), dimension(:,:,:,:), allocatable, target :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dbka_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dcka_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dckb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dakb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dakc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer :: g_dbkc_p => null()
!
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jakb
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jakc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jakb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jakc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
!     Temporary density array for each thread
      real(dp), dimension(:,:,:), allocatable :: density_oo_thread
!
      integer              :: a, b, c, a_rel, b_rel, c_rel
      type(batching_index) :: batch_a, batch_b, batch_c
      integer              :: a_batch, b_batch, c_batch
      integer              :: req_0, req_1, req_2, req_3, req_1_eri, req_a
      integer              :: n_threads, thread_n, i
      integer              :: req_single_batch
!
      type(timings) :: abc_timer
!
      abc_timer = timings('CC3 density mu ref abc', pl='v')
      call abc_timer%turn_on
!
      thread_n  = 1
      n_threads = 1
!
!$    n_threads = omp_get_max_threads()
      call mem%alloc(density_oo_thread, wf%n_o, wf%n_o, n_threads)
      call zero_array(density_oo_thread, wf%n_o**2*n_threads)
!
      batch_a = batching_index(wf%n_v)
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup_abc(req_0, req_1_eri)
      req_1_eri = req_1_eri + max(wf%n_v**2*wf%n_o, wf%n_o**3)
      req_0 = req_0 + 4*wf%n_o**3*n_threads
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_v + 2*wf%n_v**3*wf%n_o &
                       + 2*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 2*wf%n_o**3
      req_a = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 4*wf%n_o*wf%n_v + wf%n_o**2
      req_3 = 0
!
      call mem%batch_setup(batch_a, batch_b, batch_c,  &
                           req_0, req_a, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           'density_cc3_mu_ref_abc',   &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(u_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(v_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%alloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
!
         call mem%alloc(g_bdak, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
         call mem%alloc(g_dbka, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
!
         call mem%alloc(g_jakb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
         if (wf%n_o .lt. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_ljbk, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%alloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_jlkb, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%alloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%alloc(g_bdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_cdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_cdbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_adbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_adck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         call mem%alloc(g_dbka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dcka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dckb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dakb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dakc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_dbkc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
!
         call mem%alloc(g_jakb, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_jakc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%alloc(g_jbkc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_o, wf%n_v, wf%n_v, batch_a%max_length)
         else
            call mem%alloc(sorting, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         end if
!
      endif
!
!     Loop over the batches in a,b,c
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(a_batch)
!
         call wf%setup_oovo_abc(g_ljak, g_ljak_p, sorting, batch_a)
         call wf%setup_ooov_abc(g_jlka, g_jlka_p, sorting, batch_a)
!
         do b_batch = 1, a_batch
!
            call batch_b%determine_limits(b_batch)
!
            call wf%setup_vvvo_abc(g_bdak, g_bdak_p, sorting, batch_b, batch_a)
            call wf%setup_vvov_abc(g_dbka, g_dbka_p, sorting, batch_b, batch_a)
            call wf%setup_ovov_abc(g_jakb, g_jakb_p, sorting, batch_a, batch_b)
!
            if (b_batch .ne. a_batch) then
!
               call wf%setup_oovo_abc(g_ljbk, g_ljbk_p, sorting, batch_b)
               call wf%setup_ooov_abc(g_jlkb, g_jlkb_p, sorting, batch_b)
               call wf%setup_vvvo_abc(g_adbk, g_adbk_p, sorting, batch_a, batch_b)
               call wf%setup_vvov_abc(g_dakb, g_dakb_p, sorting, batch_a, batch_b)
!
            else
!
               call wf%point_oovo_abc(g_ljbk_p, g_ljak, batch_b%length)
               call wf%point_ooov_abc(g_jlkb_p, g_jlka, batch_b%length)
               call wf%point_vvvo_abc(g_adbk_p, g_bdak, batch_a%length, batch_b%length)
               call wf%point_vvov_abc(g_dakb_p, g_dbka, batch_a%length, batch_b%length)
!
            endif
!
            do c_batch = 1, b_batch
!
               call batch_c%determine_limits(c_batch)
!
               if (c_batch .ne. b_batch) then ! c_batch != b_batch and c_batch != a_batch
!
                  call wf%setup_oovo_abc(g_ljck, g_ljck_p, sorting, batch_c)
!
                  call wf%setup_ooov_abc(g_jlkc, g_jlkc_p, sorting, batch_c)
!
                  call wf%setup_vvvo_abc(g_cdak, g_cdak_p, sorting, batch_c, batch_a)
                  call wf%setup_vvvo_abc(g_adck, g_adck_p, sorting, batch_a, batch_c)
                  call wf%setup_vvvo_abc(g_cdbk, g_cdbk_p, sorting, batch_c, batch_b)
                  call wf%setup_vvvo_abc(g_bdck, g_bdck_p, sorting, batch_b, batch_c)
!
                  call wf%setup_vvov_abc(g_dcka, g_dcka_p, sorting, batch_c, batch_a)
                  call wf%setup_vvov_abc(g_dakc, g_dakc_p, sorting, batch_a, batch_c)
                  call wf%setup_vvov_abc(g_dbkc, g_dbkc_p, sorting, batch_b, batch_c)
                  call wf%setup_vvov_abc(g_dckb, g_dckb_p, sorting, batch_c, batch_b)
!
                  call wf%setup_ovov_abc(g_jakc, g_jakc_p, sorting, batch_a, batch_c)
                  call wf%setup_ovov_abc(g_jbkc, g_jbkc_p, sorting, batch_b, batch_c)
!
               else if (c_batch .eq. a_batch) then ! c_batch = b_batch = a_batch
!
                  call wf%point_oovo_abc(g_ljck_p, g_ljak, batch_c%length)
                  call wf%point_ooov_abc(g_jlkc_p, g_jlka, batch_c%length)
!
                  call wf%point_vvvo_abc(g_cdak_p, g_bdak, batch_c%length, batch_a%length)
                  call wf%point_vvvo_abc(g_adck_p, g_bdak, batch_a%length, batch_c%length)
                  call wf%point_vvvo_abc(g_cdbk_p, g_bdak, batch_c%length, batch_b%length)
                  call wf%point_vvvo_abc(g_bdck_p, g_bdak, batch_b%length, batch_c%length)
!
                  call wf%point_vvov_abc(g_dcka_p, g_dbka, batch_c%length, batch_a%length)
                  call wf%point_vvov_abc(g_dakc_p, g_dbka, batch_a%length, batch_c%length)
                  call wf%point_vvov_abc(g_dckb_p, g_dbka, batch_c%length, batch_b%length)
                  call wf%point_vvov_abc(g_dbkc_p, g_dbka, batch_b%length, batch_c%length)
!
                  call wf%point_ovov_abc(g_jakc_p, g_jakb, batch_a%length, batch_c%length)
                  call wf%point_ovov_abc(g_jbkc_p, g_jakb, batch_b%length, batch_c%length)
!
               else ! c_batch == b_batch != a_batch
!
                  call wf%point_oovo_abc(g_ljck_p, g_ljbk, batch_c%length)
                  call wf%point_ooov_abc(g_jlkc_p, g_jlkb, batch_c%length)
!
                  call wf%setup_vvvo_abc(g_cdbk, g_cdbk_p, sorting, batch_c, batch_b)
                  call wf%point_vvvo_abc(g_bdck_p, g_cdbk, batch_b%length, batch_c%length)
                  call wf%point_vvvo_abc(g_cdak_p, g_bdak, batch_c%length, batch_a%length)
                  call wf%point_vvvo_abc(g_adck_p, g_adbk, batch_a%length, batch_c%length)
!
                  call wf%setup_vvov_abc(g_dckb, g_dckb_p, sorting, batch_c, batch_b)
                  call wf%point_vvov_abc(g_dcka_p, g_dbka, batch_c%length, batch_a%length)
                  call wf%point_vvov_abc(g_dakc_p, g_dakb, batch_a%length, batch_c%length)
                  call wf%point_vvov_abc(g_dbkc_p, g_dckb, batch_b%length, batch_c%length)
!
                  call wf%setup_ovov_abc(g_jbkc, g_jbkc_p, sorting, batch_c, batch_b)
                  call wf%point_ovov_abc(g_jakc_p, g_jakb, batch_a%length, batch_c%length)
!
               endif
!
!$omp parallel do private(a, a_rel, b, b_rel, c, c_rel, thread_n) &
!$omp schedule(guided, 1)
               do a = batch_a%first, batch_a%get_last()
!
                  a_rel = a - batch_a%first + 1
!
!$                thread_n = omp_get_thread_num() + 1
!
                  do b = batch_b%first, min(batch_b%get_last(), a)
!
                     b_rel = b - batch_b%first + 1
!
                     do c = batch_c%first, min(batch_c%get_last(), b)
!
                        if (a .eq. b .and. a .eq. c) then
                           cycle
                        end if
!
                        c_rel = c - batch_c%first + 1
!
!                       Construct t^abc_ijk for given a,b,c
                        call wf%construct_W_abc(a, b, c,                   &
                                                t_ijk(:,:,:,thread_n),     &
                                                u_ijk(:,:,:,thread_n),     &
                                                t_ijab,                    &
                                                g_ljak_p(:,:,:,a_rel),     &
                                                g_ljbk_p(:,:,:,b_rel),     &
                                                g_ljck_p(:,:,:,c_rel),     &
                                                g_bdak_p(:,:,b_rel,a_rel), &
                                                g_cdak_p(:,:,c_rel,a_rel), &
                                                g_cdbk_p(:,:,c_rel,b_rel), &
                                                g_adbk_p(:,:,a_rel,b_rel), &
                                                g_adck_p(:,:,a_rel,c_rel), &
                                                g_bdck_p(:,:,b_rel,c_rel))
!
                        call wf%divide_by_orbital_differences_abc(a, b, c, &
                                                      t_ijk(:,:,:,thread_n))
!
                        call wf%construct_W_abc(a, b, c,                   &
                                                tbar_ijk(:,:,:,thread_n),  &
                                                v_ijk(:,:,:,thread_n),     &
                                                tbar_ijab,                 &
                                                g_jlka_p(:,:,:,a_rel),     &
                                                g_jlkb_p(:,:,:,b_rel),     &
                                                g_jlkc_p(:,:,:,c_rel),     &
                                                g_dbka_p(:,:,b_rel,a_rel), &
                                                g_dcka_p(:,:,c_rel,a_rel), &
                                                g_dckb_p(:,:,c_rel,b_rel), &
                                                g_dakb_p(:,:,a_rel,b_rel), &
                                                g_dakc_p(:,:,a_rel,c_rel), &
                                                g_dbkc_p(:,:,b_rel,c_rel))
!
                        call wf%outer_product_terms_l3_abc(a, b, c, tbar_ia,          &
                                                           tbar_ijab,                 &
                                                           tbar_ijk(:,:,:,thread_n),  &
                                                           wf%fock_ia,                &
                                                           g_jakb_p(:,:,a_rel,b_rel), &
                                                           g_jakc_p(:,:,a_rel,c_rel), &
                                                           g_jbkc_p(:,:,b_rel,c_rel))
!
                        call construct_contravariant_t3(tbar_ijk(:,:,:,thread_n), &
                                                        v_ijk(:,:,:,thread_n), wf%n_o)
!
                        call wf%divide_by_orbital_differences_abc(a, b, c, &
                                                                  tbar_ijk(:,:,:,thread_n), &
                                                                  omega, cvs, rm_core)
!
                        call wf%density_cc3_mu_ref_oo(a, b, c, density_oo_thread(:,:,thread_n), &
                                                      t_ijk(:,:,:,thread_n),                    &
                                                      u_ijk(:,:,:,thread_n),                    &
                                                      tbar_ijk(:,:,:,thread_n),                 &
                                                      v_ijk(:,:,:,thread_n))
!
                     enddo ! loop over c
                  enddo ! loop over b
               enddo ! loop over a
!$omp end parallel do
            enddo ! batch_c
         enddo ! batch_b
      enddo ! batch_a
!
!     Collect contributions to the oo-block from all threads
!
      do i = 1, n_threads
         call daxpy(wf%n_o**2, one, density_oo_thread(:,:,i), 1, density_oo, 1)
      enddo
!
      call mem%dealloc(density_oo_thread, wf%n_o, wf%n_o, n_threads)
!
      call mem%dealloc(t_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(u_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(v_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_bdak, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
         call mem%dealloc(g_dbka, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
         call mem%dealloc(g_jakb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
         if (wf%n_o .lt. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_ljbk, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%dealloc(g_bdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_cdak, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_cdbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_adbk, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_adck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         call mem%dealloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_jlkb, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         call mem%dealloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
!
         call mem%dealloc(g_dbka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dcka, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dckb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dakb, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dakc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_dbkc, wf%n_o, wf%n_v, batch_a%max_length, batch_a%max_length)
!
         call mem%dealloc(g_jakb, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_jakc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
         call mem%dealloc(g_jbkc, wf%n_o, wf%n_o, batch_a%max_length, batch_a%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_o, wf%n_v, wf%n_v, batch_a%max_length)
         else
            call mem%dealloc(sorting, wf%n_o, wf%n_o, wf%n_o, batch_a%max_length)
         end if
!
      endif
!
      call mem%batch_finalize()
!
      call abc_timer%turn_off
!
   end subroutine density_cc3_mu_ref_abc_cc3
!
!
   module subroutine density_cc3_mu_ref_oo_cc3(wf, a, b, c, density_oo, t_ijk, &
                                                   u_ijk, tbar_ijk, v_ijk)
!!
!!    One electron density excited-determinant/reference oo-term
!!    Written by Alexander C. Paul, August 2019
!!
!!    Calculates CC3 terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    explicit Term:
!!          D_kl += -1/2 sum_ij,abc t^abc_ijk tbar^abc_ijl
!!
!!    All permutations for a,b,c have to be considered
!!    due to the restrictions in the a,b,c loops
!!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout)  :: density_oo
!
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)  :: t_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: u_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(in)  :: tbar_ijk
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_o), intent(out) :: v_ijk
!
      real(dp) :: factor_ab, factor_bc
!
      if (a .ne. b .and. b .ne. c) then
         factor_ab = -one
         factor_bc = -one
      else if (b .eq. c) then
         factor_ab = -one
         factor_bc = -half
      else
         factor_ab = -half
         factor_bc = -one
      end if
!
!     D_kl -= 1/2 sum_ij t^abc_ijk tbar^abc_ijl
      call wf%density_cc3_oo_vv_N7_TN_permutation(wf%n_o, t_ijk, tbar_ijk, &
                                                  density_oo, factor_ab)
!
!     D_kl -= 1/2 sum_ij t^abc_kij tbar^abc_lij
      call wf%density_cc3_oo_vv_N7_NT_permutation(wf%n_o, t_ijk, tbar_ijk, &
                                                  density_oo, factor_bc)
!
      if (a .ne. b .and. b .ne. c) then
!
!        ikj -> ijk
!        ilj -> ijl
         call sort_123_to_132(t_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
         call sort_123_to_132(tbar_ijk, v_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!        D_kl -= 1/2 sum_ij t^abc_ikj tbar^abc_ilj
         call wf%density_cc3_oo_vv_N7_TN_permutation(wf%n_o, u_ijk, v_ijk, &
                                                     density_oo, -one)
!
      end if
!
   end subroutine density_cc3_mu_ref_oo_cc3
!
!
   module subroutine density_cc3_mu_ref_ijk_cc3(wf, density_ov, density_vv, omega, &
                                                tbar_ai, tbar_abij, t_abij,        &
                                                cvs, rm_core, keep_Y)
!!
!!    One electron density excited-determinant/reference term
!!    in batches of the occupied orbitals i,j,k
!!    Written by Alexander C. Paul, July 2019
!!    Adapted to first construct a covariant intermediate for L3 from
!!    which the full L3 is obtained by a linear combination
!!    by Alexander C. Paul and Rolf H. Myhre, Sep 2020
!!    For that c^ab_ij = 1/3 (2 L^ab_ij + L^ba_ij)
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    Construct t^abc_ijk and tbar^abc_ijk in batches of a,b,c.
!!    Calls routines calculating the actual density contributions
!!    to the ov- and vv-blocks.
!!
!!    t_mu3 = -< mu3| [U,T2] |HF > (eps_mu3)^-1
!!    tbar_mu3 = (- eps_mu3)^-1 (tbar_mu1 < mu1 | [H,tau_nu3] |R >
!!                             + tbar_mu2 < mu2 | [H,tau_nu3] |R >
!!
!!    ov-block:
!!       rho^L_kc += sum_{abij} tbar^ab_ij (t^abc_ijk - t^bac_ijk)
!!       rho^L_ld -= sum_{abijk} tbar^abc_ijk t^ab_lj t^dc_ik
!!
!!    vv-block:
!!       rho^L_cd += 1/2 sum_{abijk} tbar^abc_ijk t^abd_ijk
!!
      use reordering, only: construct_contravariant_t3, add_21_to_12
      use array_utilities, only: zero_array
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      logical, intent(in)  :: cvs, rm_core
!
!     If present and true the intermediate Y_clik will be stored on file
      logical, intent(in), optional :: keep_Y
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t_abij
!
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: tbar_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     Help array used for sorting integrals and amplitudes
      real(dp), dimension(:,:,:,:), allocatable :: sorting
!
!     Integrals for the construction of the T3-amplitudes
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdcj
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_bdck
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdci_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdcj_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_bdck_p => null()
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
!     Integrals for the construction of the triples multipliers
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_dbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbic_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_dbkc_p => null()
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: g_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: g_jbkc_p => null()
!
!     Intermediate ov-density ordered vo
      real(dp), dimension(:,:), allocatable :: density_ai
!
!     Intermediate Y_clik
      real(dp), dimension(:,:,:,:), allocatable :: Y_clik
!
      type(timings) :: ijk_timer
!
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3, req_i, req_1_eri
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: req_single_batch
      logical :: skip
!
      ijk_timer = timings('CC3 density mu ref ijk', pl='v')
      call ijk_timer%turn_on
!
      skip = .false.
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
!     Memory for sorting array and getting the integrals
      call wf%estimate_mem_integral_setup(req_0, req_1_eri)
      req_0 = req_0 + 3*wf%n_v**3 + wf%n_v*wf%n_o + wf%n_v*wf%n_o**3
      req_1_eri = req_1_eri + max(wf%n_v**3, wf%n_o**2*wf%n_v)
!
!     Need less memory if we don't need to batch, so we overwrite the maximum
!     required memory in batch_setup
!
      req_single_batch = req_0 + req_1_eri*wf%n_o + 2*wf%n_v**3*wf%n_o &
                       + 2*wf%n_v*wf%n_o**3 + (wf%n_v*wf%n_o)**2
!
      req_1 = 2*wf%n_v**3
      req_i = req_1 + req_1_eri ! Mem for integral setup only needed for 1 index.
      req_2 = 4*wf%n_o*wf%n_v + wf%n_v**2
      req_3 = 0
!
      call mem%batch_setup(batch_i, batch_j, batch_k,  &
                           req_0, req_i, req_1, req_1, &
                           req_2, req_2, req_2, req_3, &
                           'density_cc3_mu_ref_ijk',   &
                           req_single_batch=req_single_batch)
!
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(tbar_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%alloc(density_ai, wf%n_v, wf%n_o)
      call zero_array(density_ai, wf%n_t1)
!
      call mem%alloc(Y_clik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(Y_clik, wf%n_v*wf%n_o**3)
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%alloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%alloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%alloc(sorting, wf%n_v, wf%n_v, wf%n_v,batch_i%max_length)
         else
            call mem%alloc(sorting, wf%n_v, wf%n_o, wf%n_o,batch_i%max_length)
         end if
!
      end if
!
!     Loop over the batches in i,j,k
!     Read integrals and assign pointers
!     Without pointers we'll have to use three times as much
!     memory for the non-batching case
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call wf%setup_vvvo(g_bdci, g_bdci_p, sorting, batch_i)
!
         call wf%setup_vvov(g_dbic, g_dbic_p, sorting, batch_i, left=.true.)
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call wf%setup_oovo(g_ljci, g_ljci_p, sorting, batch_j, batch_i)
!
            call wf%setup_ooov(g_jlic, g_jlic_p, sorting, batch_j, batch_i)
!
            call wf%setup_ovov(g_ibjc, g_ibjc_p, sorting, batch_i, batch_j)
!
            if (j_batch .ne. i_batch) then
!
               call wf%setup_vvvo(g_bdcj, g_bdcj_p, sorting, batch_j)
!
               call wf%setup_vvov(g_dbjc, g_dbjc_p, sorting, batch_j, left=.true.)
!
               call wf%setup_oovo(g_licj, g_licj_p, sorting, batch_i, batch_j)
!
               call wf%setup_ooov(g_iljc, g_iljc_p, sorting, batch_i, batch_j)
!
            else
!
               call wf%point_vvvo(g_bdcj_p, g_bdci, batch_j%length)
!
               call wf%point_vvvo(g_dbjc_p, g_dbic, batch_j%length)
!
               call wf%point_vooo(g_licj_p, g_ljci, batch_i%length, batch_j%length)
!
               call wf%point_vooo(g_iljc_p, g_jlic, batch_i%length, batch_j%length)
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call wf%setup_vvvo(g_bdck, g_bdck_p, sorting, batch_k)
!
                  call wf%setup_vvov(g_dbkc, g_dbkc_p, sorting, batch_k, left=.true.)
!
                  call wf%setup_oovo(g_lick, g_lick_p, sorting, batch_i, batch_k)
                  call wf%setup_oovo(g_ljck, g_ljck_p, sorting, batch_j, batch_k)
                  call wf%setup_oovo(g_lkci, g_lkci_p, sorting, batch_k, batch_i)
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ooov(g_ilkc, g_ilkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ooov(g_jlkc, g_jlkc_p, sorting, batch_j, batch_k)
                  call wf%setup_ooov(g_klic, g_klic_p, sorting, batch_k, batch_i)
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
!
                  call wf%setup_ovov(g_ibkc, g_ibkc_p, sorting, batch_i, batch_k)
                  call wf%setup_ovov(g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdci, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbic, batch_k%length)
!
                  call wf%point_vooo(g_lick_p, g_ljci, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_ljci, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_lkcj_p, g_ljci, batch_k%length, batch_j%length)
!
                  call wf%point_vooo(g_ilkc_p, g_jlic, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_jlic, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
                  call wf%point_vooo(g_kljc_p, g_jlic, batch_k%length, batch_j%length)
!
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
                  call wf%point_vvoo(g_jbkc_p, g_ibjc, batch_j%length, batch_k%length)
!
               else ! k_batch == j_batch != i_batch
!
                  call wf%point_vvvo(g_bdck_p, g_bdcj, batch_k%length)
!
                  call wf%point_vvvo(g_dbkc_p, g_dbjc, batch_k%length)
!
                  call wf%setup_oovo(g_lkcj, g_lkcj_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_lick_p, g_licj, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_ljck_p, g_lkcj, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_lkci_p, g_ljci, batch_k%length, batch_i%length)
!
                  call wf%setup_ooov(g_kljc, g_kljc_p, sorting, batch_k, batch_j)
                  call wf%point_vooo(g_ilkc_p, g_iljc, batch_i%length, batch_k%length)
                  call wf%point_vooo(g_jlkc_p, g_kljc, batch_j%length, batch_k%length)
                  call wf%point_vooo(g_klic_p, g_jlic, batch_k%length, batch_i%length)
!
                  call wf%setup_ovov(g_jbkc, g_jbkc_p, sorting, batch_j, batch_k)
                  call wf%point_vvoo(g_ibkc_p, g_ibjc, batch_i%length, batch_k%length)
!
               endif
!
               do i = batch_i%first, batch_i%get_last()
!
                  i_rel = i - batch_i%first + 1
!
                  do j = batch_j%first, min(batch_j%get_last(), i)
!
                     j_rel = j - batch_j%first + 1
!
                     do k = batch_k%first, min(batch_k%get_last(), j)
!
                        if (i .eq. j .and. i .eq. k) cycle
!
!                       Check for core orbitals:
!                       cvs: i,j,k cannot all correspond to valence orbitals
!
!                       rm_core: in the contraction X_ck = L^ab_ij t^abc_ijk
!                                one index (k) can correspond to a general orbital
!                                and we can only cycle if 2 or more indices
!                                are core orbitals
!
                        if(cvs .and. .not. wf%one_core_index(i, j, k, cvs)) cycle
!
                        if (rm_core) then
                           if(wf%two_core_indices(i, j, k, rm_core)) then
                              cycle
                           else
                              skip = wf%one_core_index(i, j, k, rm_core)
                           end if
                        end if
!
                        k_rel = k - batch_k%first + 1
!
!                       Construct covariant W_abc for given i,j,k
!                       L_abc is obtained by dividing the linear combination
!                       4V_abc - 2W_bac - 2W_cba - 2W_acb + W_bca + W_cab
!                       by omega - eps^abc_ijk
!
                        if (.not. skip) then
!
                           call wf%construct_W(i, j, k, sorting,          &
                                               tbar_abc, tbar_abij,       &
                                               g_dbic_p(:,:,:,i_rel),     &
                                               g_dbjc_p(:,:,:,j_rel),     &
                                               g_dbkc_p(:,:,:,k_rel),     &
                                               g_jlic_p(:,:,j_rel,i_rel), &
                                               g_klic_p(:,:,k_rel,i_rel), &
                                               g_kljc_p(:,:,k_rel,j_rel), &
                                               g_iljc_p(:,:,i_rel,j_rel), &
                                               g_ilkc_p(:,:,i_rel,k_rel), &
                                               g_jlkc_p(:,:,j_rel,k_rel))
!
                           call wf%outer_product_terms_l3(i, j, k, tbar_ai, tbar_abij, &
                                                          tbar_abc, wf%fock_ia,        &
                                                          g_ibjc_p(:,:,i_rel,j_rel),   &
                                                          g_ibkc_p(:,:,i_rel,k_rel),   &
                                                          g_jbkc_p(:,:,j_rel,k_rel))
!
                           call construct_contravariant_t3(tbar_abc, sorting, wf%n_v)
!
                           call wf%divide_by_orbital_differences(i, j, k, tbar_abc, omega)
!
                        end if
!
!                       construct t3 for fixed i,j,k
                        call wf%construct_W(i, j, k, u_abc, t_abc, t_abij, &
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
                        if (.not. skip) &
                           call wf%density_cc3_mu_ref_vv(i, j, k, density_vv, t_abc, &
                                                         u_abc, tbar_abc, sorting)
!
!                       Need contravariant t3 for X_ai intermediate
                        call construct_contravariant_t3(t_abc, u_abc, wf%n_v)
!
                        call wf%construct_x_ai_intermediate(i, j, k, t_abc, sorting,  &
                                                            tbar_abij, density_ai)
!
                        if (.not. skip) &
                           call wf%construct_y_vooo_intermediate(i, j, k, tbar_abc, sorting, &
                                                                 t_abij, Y_clik)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(tbar_abc, wf%n_v, wf%n_v, wf%n_v)
!
      if (batch_i%num_batches .eq. 1) then
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         end if
!
      else ! batching
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
!
         call mem%dealloc(g_ljci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkci, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lkcj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_licj, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_lick, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ljck, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%max_length, batch_i%max_length)
!
         call mem%dealloc(g_ibjc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_ibkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
         call mem%dealloc(g_jbkc, wf%n_v, wf%n_v, batch_i%max_length, batch_i%max_length)
!
         if (wf%n_o .le. wf%n_v) then
            call mem%dealloc(sorting, wf%n_v, wf%n_v, wf%n_v, batch_i%max_length)
         else
            call mem%dealloc(sorting, wf%n_v, wf%n_o, wf%n_o, batch_i%max_length)
         end if
!
      endif
!
      call add_21_to_12(one, density_ai, density_ov, wf%n_o, wf%n_v)
!
      call mem%dealloc(density_ai, wf%n_v, wf%n_o)
!
!     Write Y_clik to file if keep_Y present and true
!     Reused for the right transition density matrices
!
      if (present(keep_Y)) then
         if (keep_Y) then
!
            wf%Y_cmjk_tbar = direct_stream_file('Y_cmjk_tbar', wf%n_v*wf%n_o**2)
            call wf%Y_cmjk_tbar%open_('write')
!
            call wf%Y_cmjk_tbar%write_(Y_clik, 1, wf%n_o)
!
            call wf%Y_cmjk_tbar%close_()
!
         end if
      end if
!
      call mem%batch_finalize()
!
!     D_ld -= sum_cik Y_clik t^dc_ik
      call wf%density_cc3_Y_vooo_ov(density_ov, t_abij, Y_clik)
      call mem%dealloc(Y_clik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call ijk_timer%turn_off
!
   end subroutine density_cc3_mu_ref_ijk_cc3
!
!
   module subroutine density_cc3_mu_ref_vv_cc3(wf, i, j, k, density_vv, t_abc, &
                                               u_abc, tbar_abc, v_abc)
!!
!!    One electron density excited-determinant/reference vv-term
!!    Written by Alexander C. Paul, August 2019
!!
!!    Calculates CC3 terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    explicit Term:
!!          D_cd += 1/2 sum_ab,ijk tbar^abc_ijk t^abd_ijk
!!
!!    All permutations for i,j,k have to be considered
!!    due to the restrictions in the i,j,k loops
!!
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) ::  i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout)  :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
      real(dp) :: factor_ij, factor_jk
!
      if (i .ne. j .and. j .ne. k) then
         factor_ij = one
         factor_jk = one
      else if (j .eq. k) then
         factor_ij = one
         factor_jk = half
      else
         factor_ij = half
         factor_jk = one
      end if
!
!     D_cd += 1/2 sum_ab tbar^abc_ijk t^abd_ijk
      call wf%density_cc3_oo_vv_N7_TN_permutation(wf%n_v, tbar_abc, t_abc, &
                                                  density_vv, factor_ij)
!
!     D_cd += 1/2 sum_ab tbar^cab_ijk t^dab_ijk
      call wf%density_cc3_oo_vv_N7_NT_permutation(wf%n_v, tbar_abc, t_abc, &
                                                  density_vv, factor_jk)
!
      if (i .ne. j .and. j .ne. k) then
!
!        acb -> abc
!        adb -> abd
!
         call sort_123_to_132(t_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
         call sort_123_to_132(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        D_cd += 1/2 sum_ab tbar^acb_ijk t^adb_ijk
         call wf%density_cc3_oo_vv_N7_TN_permutation(wf%n_v, v_abc, u_abc, &
                                                     density_vv, one)
!
!
      end if
!
   end subroutine density_cc3_mu_ref_vv_cc3
!
!
   module subroutine construct_y_vooo_intermediate_cc3(wf, i, j, k, u_abc, &
                                                      v_abc, t2, Y_clik)
!!
!!    Construct Y_vooo intermediate
!!    Written by Alexander C. Paul and Rolf H. Myhre, April 2019
!!
!!    Y_clik = sum_abj tbar^abc_ijk * t^ab_lj
!!    used to compute the tbar3 contributions to the ov part of the GS density
!!
!!    All permutations for i,j,k have to be considered
!!    due to the restrictions in the i,j,k loops
!!
!!    NB: u_abc changes from in to out,
!!        so this should be the last routine called in a i,j,k loop
!!
      use reordering, only: sort_123_to_213, sort_123_to_312
!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: u_abc, v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t2
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out) :: Y_clik
!
!     Y_clik += sum_ab,j u_abc t^ab_lj
      call wf%construct_Y_vooo_permutation(i, k, u_abc, t2(:,:,:,j), Y_clik, one)
!
!     bca -> abc
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Y_cmkj += sum_ab,i v^abc t^ab_mi
      call wf%construct_Y_vooo_permutation(k, j, v_abc, t2(:,:,:,i), Y_clik, one)
!
!     cab -> bca -> abc
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Y_cmji += sum_ab,k u_abc t^ab_mk
      call wf%construct_Y_vooo_permutation(j, i, u_abc, t2(:,:,:,k), Y_clik, one)
!
      if (i .ne. j .and. j .ne. k) then
!
!        acb -> cab -> bca -> abc
         call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        Y_cmki = sum_ab,j v_abc t^ab_mj
         call wf%construct_Y_vooo_permutation(k, i, v_abc, t2(:,:,:,j), Y_clik, one)
!
!        cba -> acb -> cab -> bca -> abc
         call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        Y_cmij = sum_ab,k u_abc t^ab_mk
         call wf%construct_Y_vooo_permutation(i, j, u_abc, t2(:,:,:,k), Y_clik, one)
!
!        bac -> cba -> acb -> cab -> bca -> abc
         call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        Y_cljk = sum_ab,i v_abc t^ab_mi
         call wf%construct_Y_vooo_permutation(j, k, v_abc, t2(:,:,:,i), Y_clik, one)
!
      end if
!
   end subroutine construct_y_vooo_intermediate_cc3
!
!
   module subroutine density_cc3_oo_vv_N7_TN_permutation_cc3(dim_, u3, v3, &
                                                             density_block, factor)
!!
!!    Density CC3 oo vv N7 T,N permutation
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Constructs one permutation of the contraction
!!    rho_pq = sum_rs u_rsp v_rsq
!!    used for N7 terms of the oo- and vv-block of CC3 densities
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_, dim_), intent(in) :: u3, v3
!
      real(dp), dimension(dim_, dim_), intent(inout) :: density_block
!
      real(dp), intent(in) :: factor
!
      call dgemm('T','N',        &
                  dim_,          &
                  dim_,          &
                  dim_**2,       &
                  factor,        &
                  u3,            & ! u_rs_p
                  dim_**2,       &
                  v3,            & ! v_rs_q
                  dim_**2,       &
                  one,           &
                  density_block, &
                  dim_)
!
   end subroutine density_cc3_oo_vv_N7_TN_permutation_cc3
!
!
   module subroutine density_cc3_oo_vv_N7_NT_permutation_cc3(dim_, u3, v3, density_block, factor)
!!
!!    Density CC3 oo vv N7 N,T permutation
!!    Written by Alexander C. Paul, Jan 2021
!!
!!    Constructs one permutation of the contraction
!!    rho_pq = sum_rs u_prs v_qrs
!!    used for N7 terms of the oo- and vv-block of CC3 densities
!!
      implicit none
!
      integer, intent(in) :: dim_
!
      real(dp), dimension(dim_, dim_, dim_), intent(in) :: u3, v3
!
      real(dp), dimension(dim_, dim_), intent(inout) :: density_block
!
      real(dp), intent(in) :: factor
!
      call dgemm('N','T',        &
                  dim_,          &
                  dim_,          &
                  dim_**2,       &
                  factor,        &
                  u3,            & ! u_p_rs
                  dim_,          &
                  v3,            & ! v_q_rs
                  dim_,          &
                  one,           &
                  density_block, &
                  dim_)
!
   end subroutine density_cc3_oo_vv_N7_NT_permutation_cc3
!
!
end submodule mean_value_cc3
