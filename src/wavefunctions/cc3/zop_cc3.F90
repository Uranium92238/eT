!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
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
submodule (cc3_class) zop_cc3
!
!!
!!    Zeroth order properties submodule 
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
   module subroutine prepare_for_density_cc3(wf)
!!
!!    Prepare for ground state density
!!    Written by Alexander C. Paul and Rolf H. Myhre, March 2019
!!
!!    Sets up the integral files needed to construct t3 and tbar3
!!    batching over the virtual indices.
!!    g_vvvo ordered 2413
!!    g_oovo ordered 1243
!!    g_vvov ordered 1324
!!    g_ooov ordered 2134
!!    L_ovov ordered 1324
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      call wf%prepare_cc3_integrals_t3_abc_batch
      call wf%prepare_cc3_integrals_L3_abc_batch
!
   end subroutine prepare_for_density_cc3
!
!
   module subroutine construct_gs_density_cc3(wf)
!!
!!    Construct one-electron density
!!    Written by Alexander C. Paul
!!    based on construct_gs_density_ccsd by Sarai D. Folkestad
!!
!!    Constructs the one-electron density matrix in the T1 basis
!!
!!    D_pq = < Lambda| E_pq |CC >
!!
!!    Contributions to the density are split up as follows:
!!    D_pq = D_pq(ref-ref) + sum_mu tbar_mu D_pq(mu-ref)
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(:,:), allocatable :: tbar_ia
      real(dp), dimension(:,:,:,:), allocatable :: t_aibj, tbar_aibj
      real(dp), dimension(:,:,:,:), allocatable :: t_abij, tbar_abij
      real(dp), dimension(:,:,:,:), allocatable :: t_ijab, tbar_ijab
!
      real(dp), dimension(:,:), allocatable :: density_ov
!
      integer :: i, j, a, b
!
      type(timings) :: cc3_timer, ccsd_timer
      type(timings) :: cc3_ijk_timer, cc3_abc_timer
!
      cc3_ijk_timer = timings('Time in CC3 GS density ijk batching')
      cc3_abc_timer = timings('Time in CC3 GS density abc batching')
      cc3_timer     = timings('Total CC3 contribution GS density')
      ccsd_timer    = timings('Total CCSD contribution GS density')
!
      call ccsd_timer%turn_on()
!
      call zero_array(wf%density, wf%n_mo**2)
!
!     CCS contrubtions
      call wf%density_ccs_ref_ref_oo(wf%density)
      call wf%density_ccs_mu_ref_vo(wf%density, wf%t1bar)
!
!     CCSD contriubutions
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%density_doubles_mu_ref_ov(wf%density, wf%t1bar, t_aibj)
!
      call mem%alloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tbar_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%density_doubles_mu_ref_oo(wf%density, tbar_aibj, t_aibj)
      call wf%density_doubles_mu_ref_vv(wf%density, tbar_aibj, t_aibj)
!
      call ccsd_timer%turn_off()
!
!     :: CC3 contributions ::
!     -----------------------
!
      call cc3_timer%turn_on()
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(t_aibj, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(tbar_aibj, tbar_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: CC3 contribution to ov- and vv-part ::   
!     ::         in batches of i,j,k         ::
!
      call mem%alloc(density_ov, wf%n_o, wf%n_v)
      call zero_array(density_ov, wf%n_o*wf%n_v)
!
      call zero_array(wf%GS_cc3_density_vv, wf%n_v**2)
!
!     Call with omega = zero and cvs=.false. for GS density
!     Save the Y_vooo intermediate for FOP
!
      call cc3_ijk_timer%turn_on()
      call wf%density_cc3_mu_ref_ijk(density_ov, wf%GS_cc3_density_vv,  &
                                     zero, wf%t1bar, tbar_abij, t_abij, &
                                     cvs=.false., keep_Y=.true.)
      call cc3_ijk_timer%turn_off()
!
!     Copy CC3 ov and vv contributions to the density matrix
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            wf%density(i, wf%n_o+a) = wf%density(i, wf%n_o+a) + density_ov(i, a)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(density_ov, wf%n_o, wf%n_v)
!
!$omp parallel do private(a, b)
      do b = 1, wf%n_v
         do a = 1, wf%n_v
!
            wf%density(wf%n_o+a, wf%n_o+b) = wf%density(wf%n_o+a, wf%n_o+b) &
                                             + wf%GS_cc3_density_vv(a, b)
!
         enddo
      enddo
!$omp end parallel do
!
!     :: CC3 contribution to oo-part ::
!     ::     in batches of a,b,c     ::
!
!     Need t2, tbar2 in ijab ordering for the oo term
      call mem%alloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(t_abij, t_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(tbar_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_1234_to_3412(tbar_abij, tbar_ijab, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(tbar_ia, wf%n_o, wf%n_v)
      call sort_12_to_21(wf%t1bar, tbar_ia, wf%n_v, wf%n_o)
!
      call zero_array(wf%GS_cc3_density_oo, wf%n_o**2)
!
      call cc3_abc_timer%turn_on()
!
!     Call with omega = zero and cvs=.false. for GS density
!
      call wf%density_cc3_mu_ref_abc(wf%GS_cc3_density_oo,     &
                                     zero, tbar_ia, tbar_ijab, &
                                     t_ijab, cvs=.false.)
      call cc3_abc_timer%turn_off()
!
      call mem%dealloc(tbar_ia, wf%n_o, wf%n_v)
      call mem%dealloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call mem%dealloc(tbar_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
!$omp parallel do private(i, j)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
!
            wf%density(i, j) = wf%density(i, j) + wf%GS_cc3_density_oo(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
      call cc3_timer%turn_off()
!
   end subroutine construct_gs_density_cc3
!
!
   module subroutine density_cc3_mu_ref_abc_cc3(wf, density_oo, omega, tbar_ia, &
                                                tbar_ijab, t_ijab, cvs)
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
      use array_utilities, only: entrywise_product
!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), intent(in) :: omega
!
      logical, intent(in)  :: cvs
!
!     Unpacked Multipliers
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: tbar_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: tbar_ijab
!
!     Unpacked t2-amplitudes
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: t_ijab
!
!     Arrays for triples amplitudes and multipliers
      real(dp), dimension(:,:,:,:), allocatable :: t_ijk
      real(dp), dimension(:,:,:,:), allocatable :: tbar_ijk
      real(dp), dimension(:,:,:,:), allocatable :: u_ijk
      real(dp), dimension(:,:,:,:), allocatable :: v_ijk
!
      real(dp), dimension(:,:,:), allocatable :: projector_ijk
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jakb
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jakc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jakb_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jakc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
!     Temporary density array for each thread
!
      real(dp), dimension(:,:,:), allocatable :: density_oo_thread
!
      integer              :: a, b, c, a_rel, b_rel, c_rel
      type(batching_index) :: batch_a, batch_b, batch_c
      integer              :: a_batch, b_batch, c_batch
      integer              :: req_0, req_1, req_2, req_3
      integer              :: n_threads, thread_n, i
!
      thread_n  = 1
      n_threads = 1
!
!     Allocate and construct the projector for triples amplitudes if CVS is requested
      if(cvs) then
!
         call mem%alloc(projector_ijk, wf%n_o, wf%n_o, wf%n_o)
         call wf%get_triples_cvs_projector_abc_batch(projector_ijk)
!
      end if
!
!
!$    n_threads = omp_get_max_threads()
      call mem%alloc(density_oo_thread, wf%n_o, wf%n_o, n_threads)
      call zero_array(density_oo_thread, wf%n_o**2*n_threads)
      
!
!     Setup and Batching loops
      req_0 = 4*(wf%n_o)**3*n_threads
      req_1 = 2*(wf%n_o)**3
      req_2 = 2*(wf%n_o)*(wf%n_v) + wf%n_o**2!
      req_3 = 0
!
      batch_a = batching_index(wf%n_v)
      batch_b = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)
!
      call mem%batch_setup_ident(batch_a, batch_b, batch_c, &
                           req_0, req_1, req_2, req_3, buffer_size = zero)
!
!     Allocate integral arrays for T3-amplitudes
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%alloc(g_bdak, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
      else ! batching
!
         call batch_a%determine_limits(1)
!
!        Ordered such that batching indices are at the end
!
         call mem%alloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%alloc(g_ljbk, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
!
         call mem%alloc(g_bdak, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(g_cdak, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(g_cdbk, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(g_adbk, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(g_adck, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
!
      endif
!
!     Arrays for the triples amplitudes
      call mem%alloc(t_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(u_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(v_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%alloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
!
!
!     Allocate integral arrays for Tbar3 multipliers
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%alloc(g_dbka, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
         call mem%alloc(L_jakb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
!
      else ! batching
!
         call batch_a%determine_limits(1)
!
!        Ordered such that batching indices are at the end
!
         call mem%alloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%alloc(g_jlkb, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%alloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
!
         call mem%alloc(g_dbka, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%alloc(g_dcka, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%alloc(g_dckb, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%alloc(g_dakb, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%alloc(g_dakc, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%alloc(g_dbkc, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
!
         call mem%alloc(L_jakb, wf%n_o, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(L_jakc, wf%n_o, wf%n_o, batch_a%length, batch_a%length)
         call mem%alloc(L_jbkc, wf%n_o, wf%n_o, batch_a%length, batch_a%length)
!
      endif
!
      call wf%g_bdck_t_v%open_('read')
      call wf%g_ljck_t_v%open_('read')
      call wf%g_dbkc_t_v%open_('read')
      call wf%g_jlkc_t_v%open_('read')
      call wf%L_jbkc_t_v%open_('read')
!
!     Loops over the batches
!
      do a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(a_batch)
!
         call single_record_reader(batch_a, wf%g_ljck_t_v, g_ljak, wf%g_jlkc_t_v, g_jlka)
         g_ljak_p => g_ljak
         g_jlka_p => g_jlka
!
         do b_batch = 1, a_batch
!
            call batch_b%determine_limits(b_batch)
!
            call compound_record_reader(batch_b, batch_a, wf%g_bdck_t_v, g_bdak, &
                                       wf%g_dbkc_t_v, g_dbka)
            call compound_record_reader(batch_a, batch_b, wf%L_jbkc_t_v, L_jakb)
            g_bdak_p => g_bdak
            g_dbka_p => g_dbka
            L_jakb_p => L_jakb
!
            if (b_batch .ne. a_batch) then
!
               call single_record_reader(batch_b, wf%g_ljck_t_v, g_ljbk, wf%g_jlkc_t_v, g_jlkb)
               g_ljbk_p => g_ljbk
               g_jlkb_p => g_jlkb
!
               call compound_record_reader(batch_a, batch_b, wf%g_bdck_t_v, g_adbk, &
                                          wf%g_dbkc_t_v, g_dakb)
               g_adbk_p => g_adbk
               g_dakb_p => g_dakb
!
            else
!
               g_ljbk_p => g_ljak
               g_jlkb_p => g_jlka
!
               g_adbk_p => g_bdak
               g_dakb_p => g_dbka
!
            endif
!
            do c_batch = 1, b_batch
!
               call batch_c%determine_limits(c_batch)
!
               if (c_batch .ne. b_batch) then ! c_batch != b_batch and c_batch != a_batch
!
                  call single_record_reader(batch_c, wf%g_ljck_t_v, g_ljck, wf%g_jlkc_t_v, g_jlkc)
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
!
                  call compound_record_reader(batch_c, batch_a, wf%g_bdck_t_v, g_cdak, &
                                             wf%g_dbkc_t_v, g_dcka)
                  g_cdak_p => g_cdak
                  g_dcka_p => g_dcka
!
                  call compound_record_reader(batch_a, batch_c, wf%g_bdck_t_v, g_adck, &
                                             wf%g_dbkc_t_v, g_dakc, wf%L_jbkc_t_v, L_jakc)
                  g_adck_p => g_adck
                  g_dakc_p => g_dakc
                  L_jakc_p => L_jakc
!
                  call compound_record_reader(batch_c, batch_b, wf%g_bdck_t_v, g_cdbk, &
                                             wf%g_dbkc_t_v, g_dckb)
                  g_cdbk_p => g_cdbk
                  g_dckb_p => g_dckb
!
                  call compound_record_reader(batch_b, batch_c, wf%g_bdck_t_v, g_bdck, &
                                             wf%g_dbkc_t_v, g_dbkc, wf%L_jbkc_t_v, L_jbkc)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
                  L_jbkc_p => L_jbkc
!
               else if (c_batch .eq. a_batch) then ! c_batch = b_batch = a_batch
!
                  g_ljck_p => g_ljak
                  g_jlkc_p => g_jlka
!
                  g_cdak_p => g_bdak
                  g_dcka_p => g_dbka
!                  
                  g_adck_p => g_bdak
                  g_dakc_p => g_dbka
                  L_jakc_p => L_jakb
!
                  g_cdbk_p => g_bdak
                  g_dckb_p => g_dbka
!
                  g_bdck_p => g_bdak
                  g_dbkc_p => g_dbka
                  L_jbkc_p => L_jakb
!
               else ! c_batch == b_batch != a_batch
!
                  g_ljck_p => g_ljbk
                  g_jlkc_p => g_jlkb
!
                  g_cdak_p => g_bdak
                  g_dcka_p => g_dbka
!
                  g_adck_p => g_adbk
                  g_dakc_p => g_dakb
                  L_jakc_p => L_jakb
!
!                 b == c, thus L_jbkc == L_jckb
                  call compound_record_reader(batch_c, batch_b, wf%g_bdck_t_v, g_cdbk, &
                                             wf%g_dbkc_t_v, g_dckb, wf%L_jbkc_t_v, L_jbkc)
                  g_cdbk_p => g_cdbk
                  g_dckb_P => g_dckb
!
                  g_bdck_p => g_cdbk
                  g_dbkc_p => g_dckb
                  L_jbkc_P => L_jbkc
!
               endif
!
!$omp parallel do private(a, a_rel, b, b_rel, c, c_rel, thread_n) & 
!$omp schedule(guided, 1)
               do a = batch_a%first, batch_a%last
!
                  a_rel = a - batch_a%first + 1
!
!$                thread_n = omp_get_thread_num() + 1
!
                  do b = batch_b%first, min(batch_b%last, a)
!
                     b_rel = b - batch_b%first + 1
!
                     do c = batch_c%first, min(batch_c%last, b)
!
                        if (a .eq. b .and. a .eq. c) then
                           cycle
                        end if
!
                        c_rel = c - batch_c%first + 1
!
!                       Construct t^abc_ijk for given a,b,c
                        call wf%omega_cc3_W_calc_abc_batch(a, b, c,                       &
                                                            t_ijk(:,:,:,thread_n),        &
                                                            u_ijk(:,:,:,thread_n),        &
                                                            t_ijab,                       &
                                                            g_ljak_p(:,:,:,a_rel),        &
                                                            g_ljbk_p(:,:,:,b_rel),        &
                                                            g_ljck_p(:,:,:,c_rel),        &
                                                            g_bdak_p(:,:,b_rel,a_rel),    &
                                                            g_cdak_p(:,:,c_rel,a_rel),    &
                                                            g_cdbk_p(:,:,c_rel,b_rel),    &
                                                            g_adbk_p(:,:,a_rel,b_rel),    &
                                                            g_adck_p(:,:,a_rel,c_rel),    &
                                                            g_bdck_p(:,:,b_rel,c_rel))
!
                        call wf%omega_cc3_eps_abc_batch(a, b, c, t_ijk(:,:,:,thread_n))
!
!                       c3_calc does not zero out the array
                        call zero_array(tbar_ijk(:,:,:,thread_n), wf%n_o**3)
!
!                       construct tbar3 for fixed a,b,c
                        call wf%jacobian_transpose_cc3_c3_calc_abc_batch(a, b, c, tbar_ia,         &
                                                                        tbar_ijab,                 &
                                                                        tbar_ijk(:,:,:,thread_n),  &
                                                                        u_ijk(:,:,:,thread_n),     &
                                                                        v_ijk(:,:,:,thread_n),     &
                                                                        wf%fock_ia,                &
                                                                        L_jakb_p(:,:,a_rel,b_rel), &
                                                                        L_jakc_p(:,:,a_rel,c_rel), &
                                                                        L_jbkc_p(:,:,b_rel,c_rel), &
                                                                        g_jlka_p(:,:,:,a_rel),     &
                                                                        g_jlkb_p(:,:,:,b_rel),     &
                                                                        g_jlkc_p(:,:,:,c_rel),     &
                                                                        g_dbka_p(:,:,b_rel,a_rel), &
                                                                        g_dcka_p(:,:,c_rel,a_rel), &
                                                                        g_dckb_p(:,:,c_rel,b_rel), &
                                                                        g_dakb_p(:,:,a_rel,b_rel), &
                                                                        g_dakc_p(:,:,a_rel,c_rel), &
                                                                        g_dbkc_p(:,:,b_rel,c_rel))
!                       ! omega = 0
                        call wf%omega_cc3_eps_abc_batch(a, b, c, tbar_ijk(:,:,:,thread_n), omega) 
!
                        if (cvs) then
                           call entrywise_product(wf%n_o**3,tbar_ijk(:,:,:,thread_n), projector_ijk)
                        end if
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
!     Close files
!
      call wf%g_bdck_t_v%close_()
      call wf%g_ljck_t_v%close_()
      call wf%g_dbkc_t_v%close_()
      call wf%g_jlkc_t_v%close_()
      call wf%L_jbkc_t_v%close_()
!
!     Deallocate the integral arrays
!
      if (batch_a%num_batches .eq. 1) then ! no batching
!
         call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_bdak, wf%n_v, wf%n_o, wf%n_v, wf%n_v)
!
         call mem%dealloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, wf%n_v)
         call mem%dealloc(g_dbka, wf%n_o, wf%n_v, wf%n_v, wf%n_v)
         call mem%dealloc(L_jakb, wf%n_o, wf%n_o, wf%n_v, wf%n_v)         
!
      else ! batching
!
         call batch_a%determine_limits(1)
!
!        Ordered such that batching indices are at the end
!
         call mem%dealloc(g_ljak, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%dealloc(g_ljbk, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
!
         call mem%dealloc(g_bdak, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(g_cdak, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(g_cdbk, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(g_adbk, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(g_adck, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_o, batch_a%length, batch_a%length)
!
         call mem%dealloc(g_jlka, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%dealloc(g_jlkb, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
         call mem%dealloc(g_jlkc, wf%n_o, wf%n_o, wf%n_o, batch_a%length)
!
         call mem%dealloc(g_dbka, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%dealloc(g_dcka, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%dealloc(g_dckb, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%dealloc(g_dakb, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%dealloc(g_dakc, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
         call mem%dealloc(g_dbkc, wf%n_o, wf%n_v, batch_a%length, batch_a%length)
!
         call mem%dealloc(L_jakb, wf%n_o, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(L_jakc, wf%n_o, wf%n_o, batch_a%length, batch_a%length)
         call mem%dealloc(L_jbkc, wf%n_o, wf%n_o, batch_a%length, batch_a%length)
!
      endif
!
!     Deallocate amplitude arrays
!
      call mem%dealloc(t_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(u_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(v_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
      call mem%dealloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o, n_threads)
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
!     D_kl -= 1/2 sum_ij t^abc_ijk tbar^abc_ijl
!
      call dgemm('T','N',     &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_o**2,  &
                  -half,      &
                  t_ijk,      & ! t_ij_k
                  wf%n_o**2,  &
                  tbar_ijk,   & ! tbar_ij_l
                  wf%n_o**2,  &
                  one,        & 
                  density_oo, &
                  wf%n_o)
!
!     D_kl -= 1/2 sum_ij t^abc_kij tbar^abc_lij
!
      call dgemm('N','T',     &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_o**2,  &
                  -half,      &
                  t_ijk,      & ! t_k_ij
                  wf%n_o,     &
                  tbar_ijk,   & ! tbar_l_ij
                  wf%n_o,     &
                  one,        & 
                  density_oo, &
                  wf%n_o)
!
!     sort_123_231(jki) -> kij
!
      call sort_123_to_231(t_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
      call sort_123_to_231(tbar_ijk, v_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!     D_kl -= 1/2 sum_ij t^abc_jki tbar^abc_jli
!
      call dgemm('N','T',     &
                  wf%n_o,     &
                  wf%n_o,     &
                  wf%n_o**2,  &
                  -half,      &
                  u_ijk,      & ! t_k_ij
                  wf%n_o,     &
                  v_ijk,      & ! tbar_l_ij
                  wf%n_o,     &
                  one,        & 
                  density_oo, &
                  wf%n_o)
!
      if (a .ne. b .and. b .ne. c) then
!
!        sort_123_to_213(jik) -> ijk
!        sort_123_to_213(ikj) -> kij
!
         call sort_123_to_213(t_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
         call sort_123_to_213(tbar_ijk, v_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!
!        D_kl -= 1/2 sum_ij t^abc_jik tbar^abc_jil
!
         call dgemm('T','N',     &
                     wf%n_o,     &
                     wf%n_o,     &
                     wf%n_o**2,  &
                     -half,      &
                     u_ijk,      & ! t_ij_k
                     wf%n_o**2,  &
                     v_ijk,      & ! tbar_ij_l
                     wf%n_o**2,  &
                     one,        & 
                     density_oo, &
                     wf%n_o)
!
!        D_kl -= 1/2 sum_ij t^abc_ikj tbar^abc_ilj
!
         call dgemm('N','T',     &
                     wf%n_o,     &
                     wf%n_o,     &
                     wf%n_o**2,  &
                     -half,      &
                     u_ijk,      & ! t_k_ij
                     wf%n_o,     &
                     v_ijk,      & ! tbar_l_ij
                     wf%n_o,     &
                     one,        & 
                     density_oo, &
                     wf%n_o)
!
!        sort_123_to_321(kji) -> ijk
!
         call sort_123_to_321(t_ijk, u_ijk, wf%n_o, wf%n_o, wf%n_o)
         call sort_123_to_321(tbar_ijk, v_ijk, wf%n_o, wf%n_o, wf%n_o)
!
!        D_kl -= 1/2 sum_ij t^abc_kji tbar^abc_lji
!
         call dgemm('T','N',     &
                     wf%n_o,     &
                     wf%n_o,     &
                     wf%n_o**2,  &
                     -half,      &
                     u_ijk,      & ! t_ij_k
                     wf%n_o**2,  &
                     v_ijk,      & ! tbar_ij_l
                     wf%n_o**2,  &
                     one,        & 
                     density_oo, &
                     wf%n_o)
!
      end if
!
   end subroutine density_cc3_mu_ref_oo_cc3
!
!
   module subroutine density_cc3_mu_ref_ijk_cc3(wf, density_ov, density_vv, omega, &
                                                tbar_ai, tbar_abij, t_abij,        &
                                                cvs, keep_Y)
!!
!!    One electron density excited-determinant/reference term 
!!    in batches of the occupied orbitals i,j,k
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
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
!
      logical, intent(in)  :: cvs
!
!     If present and true the intermediate Y_clik will be stored on file
      logical, intent(in), optional :: keep_Y
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
!     Unpacked Multipliers
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
!
!     Unpacked t2-amplitudes
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t_abij
!
!     Arrays for triples amplitudes and multipliers
      real(dp), dimension(:,:,:), allocatable :: t_abc
      real(dp), dimension(:,:,:), allocatable :: tbar_abc
      real(dp), dimension(:,:,:), allocatable :: u_abc
      real(dp), dimension(:,:,:), allocatable :: v_abc
!
      real(dp), dimension(:,:), allocatable :: F_ov_ck ! Transpose the fock matrix ov block
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
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibjc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_ibkc
      real(dp), dimension(:,:,:,:), allocatable, target  :: L_jbkc
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibjc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_ibkc_p => null()
      real(dp), dimension(:,:,:,:), contiguous, pointer  :: L_jbkc_p => null()
!
!     Intermediate ov-density ordered vo
      real(dp), dimension(:,:), allocatable :: density_ai
!
!     Intermediate Y_clik
      real(dp), dimension(:,:,:,:), allocatable :: Y_clik, Y_lcik
!
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3
      integer :: i, j, k, i_rel, j_rel, k_rel
      integer :: a
!
!     CVS
      integer :: i_cvs
      logical :: ijk_core
!
!     Setup and Batching loops
!
      req_0 = 4*(wf%n_v)**3 + 2*(wf%n_v)*(wf%n_o) + wf%n_v*wf%n_o**3
      req_1 = 2*(wf%n_v)**3
      req_2 = 2*(wf%n_o)*(wf%n_v) + (wf%n_v)**2
      req_3 = 0
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, buffer_size = zero)
!
!     Allocate integral arrays
!
!     Split up so that the integral and amplitude arrays are closer in mem
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
         call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
      end if
!
!     Arrays for the triples amplitudes and intermediates
      call mem%alloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Array for the ov-density in vo order
      call mem%alloc(density_ai, wf%n_v, wf%n_o)
      call zero_array(density_ai, wf%n_o*wf%n_v)
!
      call mem%alloc(tbar_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%alloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     Fock matrix subblock: Resorting for easier contractions later
      call mem%alloc(F_ov_ck, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_ov_ck, wf%n_o, wf%n_v)
!
!     vooo Intermediate
      call mem%alloc(Y_clik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(Y_clik, wf%n_v*wf%n_o**3)
!
!     Remaining integral arrays
!
      if (batch_i%num_batches .eq. 1) then ! no batching
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else ! batching
!
         call batch_i%determine_limits(1)
!
!        Ordered such that batching indices are at the end
!
         call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%alloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%alloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%alloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%alloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%alloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_dbkc_t%open_('read')
      call wf%g_jlkc_t%open_('read')
      call wf%L_jbkc_t%open_('read')
!
      do i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(i_batch)
!
         call single_record_reader(batch_i, wf%g_bdck_t, g_bdci, wf%g_dbkc_t, g_dbic)
         g_bdci_p => g_bdci
         g_dbic_p => g_dbic
!
         do j_batch = 1, i_batch
!
            call batch_j%determine_limits(j_batch)
!
            call compound_record_reader(batch_j, batch_i, wf%g_ljck_t, g_ljci, &
                                        wf%g_jlkc_t, g_jlic)
            call compound_record_reader(batch_i, batch_j, wf%L_jbkc_t, L_ibjc)
!
            g_ljci_p => g_ljci
            g_jlic_p => g_jlic
            L_ibjc_p => L_ibjc
!
            if (j_batch .ne. i_batch) then
!
               call single_record_reader(batch_j, wf%g_bdck_t, g_bdcj, wf%g_dbkc_t, g_dbjc)
               g_bdcj_p => g_bdcj
               g_dbjc_p => g_dbjc
!
               call compound_record_reader(batch_i, batch_j, wf%g_ljck_t, g_licj, &
                                           wf%g_jlkc_t, g_iljc)
               g_licj_p => g_licj
               g_iljc_p => g_iljc
!
            else
!
               g_bdcj_p => g_bdci
               g_dbjc_p => g_dbic
!
               g_licj_p => g_ljci
               g_iljc_p => g_jlic
!
            endif
!
            do k_batch = 1, j_batch
!
               call batch_k%determine_limits(k_batch)
!
               if (k_batch .ne. j_batch) then ! k_batch != j_batch, k_batch != i_batch
!
                  call single_record_reader(batch_k, wf%g_bdck_t, g_bdck, wf%g_dbkc_t, g_dbkc)
                  g_bdck_p => g_bdck
                  g_dbkc_p => g_dbkc
! 
                  call compound_record_reader(batch_k, batch_i, wf%g_ljck_t, g_lkci, &
                                              wf%g_jlkc_t, g_klic)
!
                  g_lkci_p => g_lkci
                  g_klic_p => g_klic
!
                  call compound_record_reader(batch_i, batch_k, wf%g_ljck_t, g_lick, & 
                                              wf%g_jlkc_t, g_ilkc, wf%L_jbkc_t, L_ibkc)
!
                  g_lick_p => g_lick
                  g_ilkc_p => g_ilkc
                  L_ibkc_p => L_ibkc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, &
                                              wf%g_jlkc_t, g_kljc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
!
                  call compound_record_reader(batch_j, batch_k, wf%g_ljck_t, g_ljck, & 
                                             wf%g_jlkc_t, g_jlkc, wf%L_jbkc_t, L_jbkc)
!
                  g_ljck_p => g_ljck
                  g_jlkc_p => g_jlkc
                  L_jbkc_p => L_jbkc
!
               else if (k_batch .eq. i_batch) then ! k_batch == j_batch == i_batch
!
                  g_bdck_p => g_bdci
                  g_dbkc_p => g_dbic
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
!
                  g_lick_p => g_ljci
                  g_ilkc_p => g_jlic
                  L_ibkc_p => L_ibjc
!
                  g_lkcj_p => g_ljci
                  g_kljc_p => g_jlic
!
                  g_ljck_p => g_ljci
                  g_jlkc_p => g_jlic
                  L_jbkc_p => L_ibjc
!
               else ! k_batch == j_batch != i_batch
!
                  g_bdck_p => g_bdcj
                  g_dbkc_p => g_dbjc
!
                  g_lkci_p => g_ljci
                  g_klic_p => g_jlic
!
                  g_lick_p => g_licj
                  g_ilkc_p => g_iljc
                  L_ibkc_p => L_ibjc
!
                  call compound_record_reader(batch_k, batch_j, wf%g_ljck_t, g_lkcj, &
                                             wf%g_jlkc_t, g_kljc, wf%L_jbkc_t, L_jbkc)
                  g_lkcj_p => g_lkcj
                  g_kljc_p => g_kljc
!
                  g_ljck_p => g_lkcj
                  g_jlkc_p => g_kljc
                  L_jbkc_p => L_jbkc
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
                        if (i .eq. j .and. i .eq. k) cycle
!
!
!                       Check if at least one index i,j,k is a core orbital
                        if(cvs) then
!
                           ijk_core = .false.
!
                           do i_cvs = 1, wf%n_core_MOs
!
                              if(     i .eq. wf%core_MOs(i_cvs)   &
                                 .or. j .eq. wf%core_MOs(i_cvs)   &
                                 .or. k .eq. wf%core_MOs(i_cvs))  then
!
                              ijk_core = .true.
!
                              end if
!
                           end do
!
!                       Cycle if i,j,k are not core orbitals
                        if (.not. ijk_core) cycle
!
                        end if
!
!                       c3_calc does not zero out the array
                        call zero_array(tbar_abc, wf%n_v**3)
!
                        k_rel = k - batch_k%first + 1
!
!                       construct t3 for fixed i,j,k
                        call wf%construct_W(i, j, k, t_abc, v_abc, t_abij, &
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
                        call wf%construct_x_ai_intermediate(i, j, k, t_abc, v_abc,  &
                                                            density_ai, tbar_abij)
!
!                       construct tbar3 for fixed i,j,k
                        call wf%jacobian_transpose_cc3_c3_calc(i, j ,k, tbar_ai, tbar_abij,  &
                                                               tbar_abc, u_abc, v_abc,       &
                                                               F_ov_ck,                      &
                                                               L_ibjc_p(:,:,i_rel,j_rel),    &
                                                               L_ibkc_p(:,:,i_rel,k_rel),    &
                                                               L_jbkc_p(:,:,j_rel,k_rel),    &
                                                               g_dbic_p(:,:,:,i_rel),        &
                                                               g_dbjc_p(:,:,:,j_rel),        &
                                                               g_dbkc_p(:,:,:,k_rel),        &
                                                               g_jlic_p(:,:,j_rel,i_rel),    &
                                                               g_klic_p(:,:,k_rel,i_rel),    &
                                                               g_kljc_p(:,:,k_rel,j_rel),    &
                                                               g_iljc_p(:,:,i_rel,j_rel),    &
                                                               g_ilkc_p(:,:,i_rel,k_rel),    &
                                                               g_jlkc_p(:,:,j_rel,k_rel))
!
                        call wf%divide_by_orbital_differences(i, j, k, tbar_abc, omega)
!
                        call wf%construct_y_intermediate_vo3(i, j, k, tbar_abc, u_abc,   &
                                                            t_abij, Y_clik)
!
                        call wf%density_cc3_mu_ref_vv(i, j, k, density_vv, t_abc,   &
                                                         u_abc, tbar_abc, v_abc)
!
                     enddo ! loop over k
                  enddo ! loop over j
               enddo ! loop over i
            enddo ! batch_k
         enddo ! batch_j
      enddo ! batch_i
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
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      else
         call batch_i%determine_limits(1)
!
         call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdcj, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_bdck, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbjc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
         call mem%dealloc(g_dbkc, wf%n_v, wf%n_v, wf%n_v, batch_i%length)
!
         call mem%dealloc(g_ljci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkci, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lkcj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_licj, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_lick, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(g_ljck, wf%n_o, wf%n_v, batch_i%length, batch_i%length)
!
         call mem%dealloc(g_jlic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_klic, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_kljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_iljc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_ilkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
         call mem%dealloc(g_jlkc, wf%n_v, wf%n_o, batch_i%length, batch_i%length)
!
         call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_ibkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
         call mem%dealloc(L_jbkc, wf%n_v, wf%n_v, batch_i%length, batch_i%length)
!
      endif
!
!     Deallocate amplitudes arrays and Fock matrix
!
      call mem%dealloc(t_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
      call mem%dealloc(v_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(tbar_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(F_ov_ck, wf%n_v, wf%n_o)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            density_ov(i,a) = density_ov(i, a) + density_ai(a, i)
!
         enddo 
      enddo
!
      call mem%dealloc(density_ai, wf%n_v, wf%n_o)
!
!
!     :: D_ld -= sum_cik Y_clik t^dc_ik ::
!
!     Write Y_clik to file if keep_Y present and true
!     Reused for the right transition density matrices
!
      if (present(keep_Y) .and. keep_Y) then
!
         wf%Y_clik_tbar = direct_file('Y_clik_tbar', wf%n_v*wf%n_o**2)
         call wf%Y_clik_tbar%open_('write')
!
         call single_record_writer(wf%n_o, wf%Y_clik_tbar, Y_clik)
!
         call wf%Y_clik_tbar%close_()
!
      end if
!
      call mem%alloc(Y_lcik, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_2134(Y_clik, Y_lcik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(Y_clik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call dgemm('N','T',                 &
                  wf%n_o,                 &
                  wf%n_v,                 &
                  (wf%n_v)*(wf%n_o)**2,   &
                  -one,                   &
                  Y_lcik,                 & ! Y_l_cik
                  wf%n_o,                 & 
                  t_abij,                 & ! t_d_cik
                  wf%n_v,                 &
                  one,                    &
                  density_ov,             &
                  wf%n_o)
!
      call mem%dealloc(Y_lcik, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
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
!     D_cd += 1/2 sum_ab tbar^abc_ijk t^abd_ijk
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_v,     &
                  wf%n_v**2,  &
                  half,       &
                  tbar_abc,   & ! tbar_ab_c
                  wf%n_v**2,  &
                  t_abc,      & ! t_ab_d
                  wf%n_v**2,  &
                  one,        & 
                  density_vv, &
                  wf%n_v)
!
!     D_cd += 1/2 sum_ab tbar^cab_ijk t^dab_ijk
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_v,     &
                  wf%n_v**2,  &
                  half,       &
                  tbar_abc,   & ! tbar_c_ab
                  wf%n_v,     &
                  t_abc,      & ! t_d_ab
                  wf%n_v,     &
                  one,        & 
                  density_vv, &
                  wf%n_v)
!
!     sort_123_231(bca) -> cab
!
      call sort_123_to_231(t_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
      call sort_123_to_231(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     D_cd += 1/2 sum_ab tbar^bca_ijk t^bda_ijk
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_v,     &
                  wf%n_v**2,  &
                  half,       &
                  v_abc,      & ! tbar_c_ab
                  wf%n_v,     &
                  u_abc,      & ! t_d_ab
                  wf%n_v,     &
                  one,        & 
                  density_vv, &
                  wf%n_v)
!
      if (i .ne. j .and. j .ne. k) then
!
!        sort_123_to_213(bac) -> abc
!        sort_123_to_213(acb) -> cab
!
         call sort_123_to_213(t_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
         call sort_123_to_213(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     D_cd += 1/2 sum_ab tbar^bac_ijk t^bad_ijk
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_v,     &
                  wf%n_v**2,  &
                  half,       &
                  v_abc,      & ! tbar_ab_c
                  wf%n_v**2,  &
                  u_abc,      & ! t_ab_d
                  wf%n_v**2,  &
                  one,        & 
                  density_vv, &
                  wf%n_v)
!
!     D_cd += 1/2 sum_ab tbar^acb_ijk t^adb_ijk
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_v,     &
                  wf%n_v**2,  &
                  half,       &
                  v_abc,      & ! tbar_c_ab
                  wf%n_v,     &
                  u_abc,      & ! t_d_ab
                  wf%n_v,     &
                  one,        & 
                  density_vv, &
                  wf%n_v)
!
!        sort_123_to_321(cba) -> abc
!
         call sort_123_to_321(t_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
         call sort_123_to_321(tbar_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     D_cd += 1/2 sum_ab tbar^cba_ijk t^dba_ijk
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_v,     &
                  wf%n_v**2,  &
                  half,       &
                  v_abc,      & ! tbar_ab_c
                  wf%n_v**2,  &
                  u_abc,      & ! t_ab_d
                  wf%n_v**2,  &
                  one,        & 
                  density_vv, &
                  wf%n_v)
!
      end if
!
   end subroutine density_cc3_mu_ref_vv_cc3
!
!
   module subroutine construct_y_intermediate_vo3_cc3(wf, i, j, k, tbar_abc, &
                                                      u_abc, t_abij, Y_clik)
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
!!    equal to some of the terms in jacobian_transpose construct_y_intermediates
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abij
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(out)  :: Y_clik
!
!     Y_clik = sum_ab,j tbar^abc_ijk t^ab_lj
!
      call dgemm('T','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  one,              &
                  tbar_abc,         & ! tbar_ab_c
                  wf%n_v**2,        &
                  t_abij(:,:,:,j),  & ! t_ab_l,j
                  wf%n_v**2,        &
                  one,              &
                  Y_clik(:,:,i,k),  & ! Y_c_l,ik
                  wf%n_v)
!
!     Y_clji = sum_ab,k tbar^cab_ijk t^ab_lk
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  one,              &
                  tbar_abc,         & ! tbar_c_ab
                  wf%n_v,           &
                  t_abij(:,:,:,k),  & ! t_ab_l,k
                  wf%n_v**2,        &
                  one,              &
                  Y_clik(:,:,j,i),  & ! Y_c_l,ji
                  wf%n_v)
!
!
      call sort_123_to_231(tbar_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     123->231(bca) = cab
!
!     Y_clkj = sum_ab,i tbar^bca_ijk t^ab_li
!
      call dgemm('N','N',           &
                  wf%n_v,           &
                  wf%n_o,           &
                  wf%n_v**2,        &
                  one,              &
                  u_abc,            & ! u_c_ab
                  wf%n_v,           &
                  t_abij(:,:,:,i),  & ! t_ab_l,i
                  wf%n_v**2,        &
                  one,              &
                  Y_clik(:,:,k,j),  & ! Y_c_l,kj
                  wf%n_v)
!
      if (k .ne. j .and. j .ne. i) then
!
         call sort_123_to_213(tbar_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        123->213(bac) = abc
!        123->213(acb) = cab
!
!        Y_clik = sum_ab,i tbar^bac_ijk t^ab_li
!
         call dgemm('T','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     one,              &
                     u_abc,            & ! tbar_ab_c
                     wf%n_v**2,        &
                     t_abij(:,:,:,i),  & ! t_ab_l,i
                     wf%n_v**2,        &
                     one,              &
                     Y_clik(:,:,j,k),  & ! Y_c_l,jk
                     wf%n_v)
!
!        Y_clij = sum_ab,k tbar^acb_ijk t^ab_lk
!
         call dgemm('N','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     one,              &
                     u_abc,            & ! tbar_c_ab
                     wf%n_v,           &
                     t_abij(:,:,:,k),  & ! t_ab_l,k
                     wf%n_v**2,        &
                     one,              &
                     Y_clik(:,:,i,j),  & ! Y_c_l,ij
                     wf%n_v)
!
!
         call sort_123_to_132(tbar_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!        123->132(cba) = cab
!
!        Y_clki = sum_ab,j tbar^cba_ijk t^ab_lj
!
         call dgemm('N','N',           &
                     wf%n_v,           &
                     wf%n_o,           &
                     wf%n_v**2,        &
                     one,              &
                     u_abc,            & ! tbar_c_ab
                     wf%n_v,           &
                     t_abij(:,:,:,j),  & ! t_ab_l,j
                     wf%n_v**2,        &
                     one,              &
                     Y_clik(:,:,k,i),  & ! Y_c_l,ki
                     wf%n_v)
!
      end if
!
   end subroutine construct_y_intermediate_vo3_cc3
!
!
end submodule zop_cc3
