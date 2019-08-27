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
submodule (cc3_class) gs_density
!
!!
!!    Constructs the ground state one-electron density matrix (cc3)
!!    Written by Alexander Paul, July 2019
!!
!!    D_pq = < Λ | E_pq | CC >
!!
!!    Routines can be reused for the construction of the transition density
!!
!
   implicit none
!
!
contains
!
!
   module subroutine initialize_gs_density_cc3(wf)
!!
!!    Initialize density and CC3 corrections to the GS-density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      if (.not. allocated(wf%density)) call mem%alloc(wf%density, wf%n_mo, wf%n_mo)
!
!     CC3 corrections to the GS-density are needed for the right transition density
      if (.not. allocated(wf%GS_cc3_density_oo)) call mem%alloc(wf%GS_cc3_density_oo, wf%n_o, wf%n_o)
      if (.not. allocated(wf%GS_cc3_density_vv)) call mem%alloc(wf%GS_cc3_density_vv, wf%n_v, wf%n_v)
!
   end subroutine initialize_gs_density_cc3
!
!
   module subroutine destruct_gs_density_cc3(wf)
!!
!!    Destruct density
!!    Written by Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      if (allocated(wf%density)) call mem%dealloc(wf%density, wf%n_mo, wf%n_mo)
!
!     CC3 corrections to the GS-density are needed for the right transition density
      if (.not. allocated(wf%GS_cc3_density_oo)) call mem%dealloc(wf%GS_cc3_density_oo, wf%n_o, wf%n_o)
      if (.not. allocated(wf%GS_cc3_density_vv)) call mem%dealloc(wf%GS_cc3_density_vv, wf%n_v, wf%n_v)
!
   end subroutine destruct_gs_density_cc3
!
!
   module subroutine prepare_for_density_cc3(wf)
!!
!!    Prepare the construction of the CC3 contribution to the GS density
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
!!
!!    Sets up the file containing an intermediate needed for the GS and
!!    transition density matrices.
!!
!!    Z_alik = sum_bcj tbar^abc_ijk t^bc_jl
!!
!!    And
!!    Sets up the integral files needed to construct t3 and tbar3
!!    batching over the virtual indices.
!!    g_vvvo ordered 2413
!!    g_oovo ordered 1243
!!    g_vvov ordered 3124
!!    g_ooov ordered 1324
!!    L_ovov ordered 1324
!!
      implicit none
!
      class(cc3), intent(inout) :: wf
!
      call wf%prep_cc3_integrals_t3_abc_batch
      call wf%prep_cc3_integrals_L3_abc_batch
!
   end subroutine prepare_for_density_cc3
!
!
   module subroutine construct_gs_density_cc3(wf)
!!
!!    Construct one-electron density
!!    Written by Alexander Paul
!!    based on construct_gs_density_ccsd by Sarai D. Folkestad
!!
!!    Constructs the one-electron density matrix in the T1 basis
!!
!!    D_pq = < Λ | E_pq | CC >
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
      cc3_ijk_timer   = new_timer('Time in CC3 GS density ijk batching')
      cc3_abc_timer   = new_timer('Time in CC3 GS density abc batching')
      cc3_timer      = new_timer('Total CC3 contribution GS density')
      ccsd_timer     = new_timer('Total CCSD contribution GS density')
!
      call ccsd_timer%turn_on()
!
      call zero_array(wf%density, wf%n_mo**2)
!
!     CCS contrubtions
      call wf%gs_one_el_density_ccs_oo(wf%density)
      call wf%gs_one_el_density_ccs_vo(wf%density, wf%t1bar)
!
!     CCSD contriubutions
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%gs_one_el_density_ccsd_ov(wf%density, wf%t1bar, t_aibj)
!
      call mem%alloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tbar_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%gs_one_el_density_ccsd_oo(wf%density, tbar_aibj, t_aibj)
      call wf%gs_one_el_density_ccsd_vv(wf%density, tbar_aibj, t_aibj)
!
      call ccsd_timer%turn_off()
!
!     CC3 contributions
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
      call cc3_ijk_timer%turn_on()
      call wf%gs_one_el_density_cc3_ijk(density_ov, wf%GS_cc3_density_vv, zero, &
                                        wf%t1bar, tbar_abij, t_abij, .true.)
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
!
!     :: CC3 contribution to oo-part ::
!     ::     in batches of a,b,c     ::
!
!     Need t2, tbar2 in ijab ordering for the oo term
      call mem%alloc(t_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_12_to_21(t_abij, t_ijab, wf%n_v**2, wf%n_o**2)
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(tbar_ijab, wf%n_o, wf%n_o, wf%n_v, wf%n_v)
      call sort_12_to_21(tbar_abij, tbar_ijab, wf%n_v**2, wf%n_o**2)
      call mem%dealloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%alloc(tbar_ia, wf%n_o, wf%n_v)
      call sort_12_to_21(wf%t1bar, tbar_ia, wf%n_v, wf%n_o)
!
      call zero_array(wf%GS_cc3_density_oo, wf%n_o**2)
!
      call cc3_abc_timer%turn_on()
      call wf%gs_one_el_density_cc3_abc(wf%GS_cc3_density_oo, zero, &
                                        tbar_ia, tbar_ijab, t_ijab)
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
   module subroutine gs_one_el_density_cc3_abc_cc3(wf, density_oo, omega,     &
                                                   tbar_ia, tbar_ijab, t_ijab)
!!
!!    Construct t^abc_ijk and tbar^abc_ijk in batches of a,b,c and compute
!!    the contributions to the oo-part of the ground state density matrix.
!!
!!    t_μ3 = -< μ3 |{U,T2}| HF >/ε_μ3
!!    tbar_μ3 = (- ε_μ3)^-1 (tbar_μ1 < μ1 | [H,τ_ν3] | R > + tbar_μ2 < μ2 | [H,τ_ν3] | R >
!!
!!    ρ^L_kl -= 1/2 sum_{abc}{ij} tbar^abc_ijl t^abc_ijk
!!
!!    Written by Alexander Paul, July 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), intent(in) :: omega
!
!     Unpacked Multipliers
      real(dp), dimension(wf%n_o, wf%n_v), intent(in) :: tbar_ia
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: tbar_ijab
!
!     Unpacked t2-amplitudes
      real(dp), dimension(wf%n_o, wf%n_o, wf%n_v, wf%n_v), intent(in) :: t_ijab
!
!     Arrays for triples amplitudes and multipliers
      real(dp), dimension(:,:,:), allocatable :: t_ijk
      real(dp), dimension(:,:,:), allocatable :: tbar_ijk
      real(dp), dimension(:,:,:), allocatable :: u_ijk
      real(dp), dimension(:,:,:), allocatable :: v_ijk
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
      integer              :: a, b, c, a_rel, b_rel, c_rel
      type(batching_index) :: batch_a, batch_b, batch_c
      integer              :: a_batch, b_batch, c_batch
      integer              :: req_0, req_1, req_2, req_3
!
!     Setup and Batching loops
      req_0 = 4*(wf%n_o)**3
      req_1 = 2*(wf%n_o)**3
      req_2 = 2*(wf%n_o)*(wf%n_v) + wf%n_o**2!
      req_3 = 0
!
      call batch_a%init(wf%n_v)
      call batch_b%init(wf%n_v)
      call batch_c%init(wf%n_v)
!
      call mem%batch_setup_ident(batch_a, batch_b, batch_c, &
                           req_0, req_1, req_2, req_3, zero)
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
      call mem%alloc(t_ijk, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(u_ijk, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(v_ijk, wf%n_o, wf%n_o, wf%n_o)
      call mem%alloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o)
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
               do a = batch_a%first, batch_a%last
!
                  a_rel = a - batch_a%first + 1
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
!                       c3_calc does not zero out the array
                        call zero_array(tbar_ijk, wf%n_o**3)
!
                        c_rel = c - batch_c%first + 1
!
!                       Construct t^abc_ijk for given a,b,c
                        call wf%omega_cc3_W_calc_abc_batch(a, b, c, t_ijk, u_ijk, t_ijab, &
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
                        call wf%omega_cc3_eps_abc_batch(a, b, c, t_ijk)
!
!                       construct tbar3 for fixed a,b,c
                        call wf%jacobian_transpose_cc3_c3_calc_abc_batch(a, b, c, tbar_ia, tbar_ijab, &
                                                                        tbar_ijk, u_ijk, v_ijk,       &
                                                                        wf%fock_ia,                   &
                                                                        L_jakb_p(:,:,a_rel,b_rel),    &
                                                                        L_jakc_p(:,:,a_rel,c_rel),    &
                                                                        L_jbkc_p(:,:,b_rel,c_rel),    &
                                                                        g_jlka_p(:,:,:,a_rel),        &
                                                                        g_jlkb_p(:,:,:,b_rel),        &
                                                                        g_jlkc_p(:,:,:,c_rel),        &
                                                                        g_dbka_p(:,:,b_rel,a_rel),    &
                                                                        g_dcka_p(:,:,c_rel,a_rel),    &
                                                                        g_dckb_p(:,:,c_rel,b_rel),    &
                                                                        g_dakb_p(:,:,a_rel,b_rel),    &
                                                                        g_dakc_p(:,:,a_rel,c_rel),    &
                                                                        g_dbkc_p(:,:,b_rel,c_rel))
!
                        call wf%omega_cc3_eps_abc_batch(a, b, c, tbar_ijk, omega) ! omega = 0
!
                        call wf%one_el_density_cc3_oo_N7(a, b, c, density_oo, t_ijk,   &
                                                         u_ijk, tbar_ijk, v_ijk)
!
                     enddo ! loop over a
                  enddo ! loop over b
               enddo ! loop over c
            enddo ! batch_c
         enddo ! batch_b
      enddo ! batch_a
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
      call mem%dealloc(t_ijk, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(u_ijk, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(v_ijk, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(tbar_ijk, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine gs_one_el_density_cc3_abc_cc3
!
!
   module subroutine one_el_density_cc3_oo_N7_cc3(wf, a, b, c, density_oo, t_ijk, &
                                                   u_ijk, tbar_ijk, v_ijk)
!!
!!    Calculates triples contribution to the oo-part of the GS-density
!!
!!    D_kl += -1/2 sum_ij,abc t^abc_ijk tbar^abc_ijl
!!
!!    All permutations for a,b,c have to be considered 
!!    due to the restrictions in the a,b,c loops   
!!
!!    Written by Alexander Paul, August 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: a, b, c
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out)  :: density_oo
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
   end subroutine one_el_density_cc3_oo_N7_cc3
!
!
   module subroutine gs_one_el_density_cc3_ijk_cc3(wf, density_ov, density_vv, omega, &
                                                   tbar_ai, tbar_abij, t_abij, keep_Y)
!!
!!    Construct t^abc_ijk and tbar^abc_ijk in batches of i,j,k and compute
!!    the contributions to the ov- and vv-part of the ground state density matrix.
!!
!!    t_μ3 = -< μ3 |{U,T2}| HF >/ε_μ3
!!    tbar_μ3 = (- ε_μ3)^-1 (tbar_μ1 < μ1 | [H,τ_ν3] | R > + tbar_μ2 < μ2 | [H,τ_ν3] | R >
!!
!!    ov-part:
!!       ρ^L_kc += sum_{ab}{ij} tbar^ab_ij (t^abc_ijk - t^bac_ijk)
!!
!!    vv-part:
!!       ρ^L_cd += 1/2 sum_{ab}{ijk} tbar^abc_ijk t^abd_ijk
!!       ρ^L_ld -= sum_{ab}{ijk} tbar^abc_ijk t^ab_lj t^dc_ik
!!
!!    Written by Alexander Paul, August 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), intent(in) :: omega
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
      real(dp), dimension(:,:), allocatable :: F_kc ! Transpose the fock matrix sub block
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
      integer :: i, j, k, i_rel, j_rel, k_rel
      type(batching_index) :: batch_i, batch_j, batch_k
      integer :: i_batch, j_batch, k_batch
      integer :: req_0, req_1, req_2, req_3
      integer :: a
!
!     Setup and Batching loops
!
      req_0 = 4*(wf%n_v)**3 + 2*(wf%n_v)*(wf%n_o) + wf%n_v*wf%n_o**3
      req_1 = 2*(wf%n_v)**3
      req_2 = 2*(wf%n_o)*(wf%n_v) + (wf%n_v)**2
      req_3 = 0
!
      call batch_i%init(wf%n_o)
      call batch_j%init(wf%n_o)
      call batch_k%init(wf%n_o)
!
      call mem%batch_setup_ident(batch_i, batch_j, batch_k, &
                           req_0, req_1, req_2, req_3, zero)
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
      call mem%alloc(F_kc, wf%n_v, wf%n_o)
      call sort_12_to_21(wf%fock_ia, F_kc, wf%n_o, wf%n_v)
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
                        if (i .eq. j .and. i .eq. k) then
                           cycle
                        end if
!
!                       c3_calc does not zero out the array
                        call zero_array(tbar_abc, wf%n_v**3)
!
                        k_rel = k - batch_k%first + 1
!
!                       construct t3 for fixed i,j,k
                        call wf%omega_cc3_W_calc(i, j, k, t_abc, v_abc, t_abij, &
                                                 g_bdci_p(:,:,:,i_rel), &
                                                 g_bdcj_p(:,:,:,j_rel), &
                                                 g_bdck_p(:,:,:,k_rel), &
                                                 g_ljci_p(:,:,j_rel,i_rel), &
                                                 g_lkci_p(:,:,k_rel,i_rel), &
                                                 g_lkcj_p(:,:,k_rel,j_rel), &
                                                 g_licj_p(:,:,i_rel,j_rel), &
                                                 g_lick_p(:,:,i_rel,k_rel), &
                                                 g_ljck_p(:,:,j_rel,k_rel))
!
                        call wf%omega_cc3_eps(i, j, k, t_abc)
!
                        call wf%construct_x_ai_intermediate(i, j, k, t_abc, v_abc,  &
                                                            density_ai, tbar_abij)
!
!                       construct tbar3 for fixed i,j,k
                        call wf%jacobian_transpose_cc3_c3_calc(i, j ,k, tbar_ai, tbar_abij,  &
                                                               tbar_abc, u_abc, v_abc,       &
                                                               F_kc,                         &
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
                        call wf%omega_cc3_eps(i, j, k, tbar_abc, omega)
!
                        call wf%construct_y_intermediate_vo3(i, j, k, tbar_abc, u_abc,   &
                                                            t_abij, Y_clik)
!
                        call wf%one_el_density_cc3_vv_N7(i, j, k, density_vv, t_abc,   &
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
      call mem%dealloc(F_kc, wf%n_v, wf%n_o)
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
   end subroutine gs_one_el_density_cc3_ijk_cc3
!
!
   module subroutine one_el_density_cc3_vv_N7_cc3(wf, i, j, k, density_vv, t_abc, &
                                                   u_abc, tbar_abc, v_abc)
!!
!!    Calculates triples contribution to the vv-part of the GS-density
!!
!!    D_cd += 1/2 sum_ab,ijk tbar^abc_ijk t^abd_ijk
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops   
!!
!!    Written by Alexander Paul, August 2019
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) ::  i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(out)  :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)  :: tbar_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
!     D_cd += 1/2 sum_ab tbar^abc_ijk tbar^abd_ijk
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
!     D_cd += 1/2 sum_ab tbar^cab_ijk tbar^dab_ijk
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
!     D_cd += 1/2 sum_ab tbar^bca_ijk tbar^bda_ijk
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
!     D_cd += 1/2 sum_ab tbar^bac_ijk tbar^bad_ijk
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
!     D_cd += 1/2 sum_ab tbar^acb_ijk tbar^adb_ijk
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
!     D_cd += 1/2 sum_ab tbar^cba_ijk tbar^dba_ijk
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
   end subroutine one_el_density_cc3_vv_N7_cc3
!
!
   module subroutine construct_y_intermediate_vo3_cc3(wf, i, j, k, tbar_abc, u_abc, t_abij, Y_clik)
!!
!!    Constructs the intermediate Y_clik used to compute 
!!    the triples multipliers contributions to the ov part of the GS density
!!
!!    Y_clik = sum_abj tbar^abc_ijk * t^ab_lj
!!
!!    All permutations for i,j,k have to be considered 
!!    due to the restrictions in the i,j,k loops
!!
!!    equal to some of the terms in jacobian_transpose construct_y_intermediates
!!    Written by Alexander Paul and Rolf H. Myhre, April 2019
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
end submodule
