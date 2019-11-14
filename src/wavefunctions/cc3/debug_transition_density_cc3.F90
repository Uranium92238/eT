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
submodule (cc3_class) debug_transition_density_cc3
!
!!
!!    Debug CC3 submodule
!!    Written by Alexander C. Paul September 2019
!!
!!    Constructs full v3o3 arrays for triples amplitudes using unwrapped loops
!!    
!!    Construction of the transition densities is also done with unwrapped loops
!!
!!    NB:
!!    Used to compare results from the "normal" CC3 code
!!
!
   implicit none
!
!
contains
!
!
   module subroutine left_tdm_debug_cc3(wf, state, l_tdm)
!!
!!    Left transition density matrix (debug)
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Constructs CC3 left transition density matrix using unwrapped loops
!!
!!    NB: only for debugging
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: L_k
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: l_tdm
!
      real(dp), dimension(:,:), allocatable :: L_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: L_abij
!
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj, t_aibj
!
      real(dp), dimension(:,:,:,:,:,:), allocatable :: L_abcijk
      real(dp), dimension(:,:,:,:,:,:), allocatable :: t_abcijk
!
      real(dp), dimension(:,:), allocatable :: density_oo, density_ov, density_vv
!
      integer :: i, j, a, b
!
      call zero_array(l_tdm, wf%n_mo**2)
!
      call mem%alloc(L_k, wf%n_es_amplitudes)
      call wf%read_excited_state(L_k, state, 'left')
!
!     Allocate the singles part of the excitation vector
!
      call mem%alloc(L_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, L_k(1:wf%n_t1), 1, L_ai, 1)
!
      call wf%density_ccs_mu_ref_vo(l_tdm, L_ai)
!
      call mem%alloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2, t_aibj, (wf%n_v)*(wf%n_o))
!
      call wf%density_doubles_mu_ref_ov(l_tdm, L_ai, t_aibj)
!
!     Allocate and unpack doubles part of the excitation vector
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(L_k(wf%n_t1 + 1 : wf%n_es_amplitudes), L_aibj, wf%n_t1)
!
      call mem%dealloc(L_k, wf%n_es_amplitudes)
!
      call wf%density_doubles_mu_ref_oo(l_tdm, L_aibj, t_aibj)
      call wf%density_doubles_mu_ref_vv(l_tdm, L_aibj, t_aibj)
!
      call mem%dealloc(t_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: CC3 contributions ::
!     -----------------------
!
      call mem%alloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(L_aibj, L_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Construct L3 ::
!
      call mem%alloc(L_abcijk, wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call wf%construct_full_tbar3(L_ai, L_abij, L_abcijk, wf%left_excitation_energies(state), wf%cvs)
!
      call mem%dealloc(L_ai, wf%n_v, wf%n_o)
!
!     :: ov-part of the density ::
!
      call mem%alloc(density_ov, wf%n_o, wf%n_v)
      call zero_array(density_ov, wf%n_o*wf%n_v)
!
      call wf%debug_left_ov_N7(density_ov, L_abcijk)
!
!     construct full t3
!
      call mem%alloc(t_abcijk, wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call wf%construct_full_t3(t_abcijk)
!
      call wf%debug_left_ov_N6(density_ov, L_abij, t_abcijk)
!
      call mem%dealloc(L_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            l_tdm(i, wf%n_o + a) = l_tdm(i, wf%n_o + a) + density_ov(i, a)
         enddo
      enddo
!
      call mem%dealloc(density_ov, wf%n_o, wf%n_v)
!
!     :: oo-part of the density ::
!
      call mem%alloc(density_oo, wf%n_o, wf%n_o)
      call zero_array(density_oo, wf%n_o**2)
!
      call wf%debug_left_oo(density_oo, L_abcijk, t_abcijk)
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            l_tdm(i,j) = l_tdm(i,j) + density_oo(i,j)
         enddo
      enddo
!
      call mem%dealloc(density_oo, wf%n_o, wf%n_o)
!
!     :: vv-part of the density ::
!
      call mem%alloc(density_vv, wf%n_v, wf%n_v)
      call zero_array(density_vv, wf%n_v**2)
!
      call wf%debug_left_vv(density_vv, L_abcijk, t_abcijk)
!
      do b = 1, wf%n_v
         do a = 1, wf%n_v
!
            l_tdm(wf%n_o + a, wf%n_o + b) = l_tdm(wf%n_o + a, wf%n_o + b) + density_vv(a, b)
!
         enddo
      enddo
!
      call mem%dealloc(density_vv, wf%n_v, wf%n_v)
!
      call mem%dealloc(t_abcijk, wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(L_abcijk, wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
   end subroutine left_tdm_debug_cc3
!
!
   module subroutine right_tdm_debug_cc3(wf, state, r_tdm)
!!
!!    Right transition density matrix (debug)
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    Constructs CC3 right transition density matrix using unwrapped loops
!!
!!    NB: only for debugging
!!
      implicit none
!
      class(cc3) :: wf
!
      integer, intent(in) :: state
!
      real(dp), dimension(:), allocatable :: R_k
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: r_tdm
!
      real(dp), dimension(:,:), allocatable :: R_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_aibj, R_aibj
      real(dp), dimension(:,:,:,:), allocatable :: tbar_abij, R_abij
!
      real(dp), dimension(:,:,:,:,:,:), allocatable :: R_abcijk
      real(dp), dimension(:,:,:,:,:,:), allocatable :: t_abcijk
      real(dp), dimension(:,:,:,:,:,:), allocatable :: tbar_abcijk
!
      real(dp), dimension(:,:), allocatable :: density_oo, density_ov 
      real(dp), dimension(:,:), allocatable :: density_vo, density_vv
!
      real(dp) :: tbar_R_overlap, ddot
!
      integer :: i, j, a, b
!
      call zero_array(r_tdm, (wf%n_mo)**2)
!
      call mem%alloc(R_k, wf%n_es_amplitudes)
      call wf%read_excited_state(R_k, state, 'right')
!
      call mem%alloc(R_ai, wf%n_v, wf%n_o)
!
      call dcopy(wf%n_t1, R_k(1:wf%n_t1), 1, R_ai, 1)
!
!     make sure that the correct integrals are on file
      call wf%construct_c1_integrals(R_ai)
!
      tbar_R_overlap = ddot(wf%n_v*wf%n_o, R_ai, 1, wf%t1bar, 1)
!
      call wf%density_ccs_mu_nu_oo(r_TDM, wf%t1bar, R_ai)
      call wf%density_ccs_ref_mu_ov(r_TDM, R_ai)
      call wf%density_ccs_mu_nu_vv(r_TDM, wf%t1bar, R_ai)
!
      call mem%alloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call squareup(wf%t2bar, tbar_aibj, wf%n_v*wf%n_o)
!
      call wf%density_doubles_mu_nu_ov(r_tdm, tbar_aibj, R_ai)
      call wf%density_doubles_mu_nu_vo(r_tdm, tbar_aibj, R_ai)
!
!     Allocate and unpack doubles part of the excitation vector
!
      call mem%alloc(R_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call squareup(R_k(wf%n_t1 + 1 : wf%n_es_amplitudes), R_aibj, wf%n_t1)
!
      call mem%dealloc(R_k, wf%n_es_amplitudes)
!
!     Scale the doubles vector by 1 + δ_ai,bj
!
!$omp parallel do schedule(static) private(a,i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            R_aibj(a,i,a,i) = two*R_aibj(a,i,a,i)
!
         enddo
      enddo
!$omp end parallel do

!
      tbar_R_overlap = tbar_R_overlap + half*ddot((wf%n_v*wf%n_o)**2, &
                                                   R_aibj,            &
                                                   1,                 &
                                                   tbar_aibj,         &
                                                   1)
!
!
      call wf%density_doubles_mu_ref_ov(r_tdm, wf%t1bar, R_aibj)
!
      call wf%density_doubles_mu_ref_oo(r_tdm, tbar_aibj, R_aibj)
      call wf%density_doubles_mu_ref_vv(r_tdm, tbar_aibj, R_aibj)
!
!     :: CC3 contributions ::
!     -----------------------
!
      call mem%alloc(R_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(R_aibj, R_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(R_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%alloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call sort_1234_to_1324(tbar_aibj, tbar_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(tbar_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     :: Construct tbar3 ::
!
      call mem%alloc(tbar_abcijk,wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
!
      call wf%construct_full_tbar3(wf%t1bar, tbar_abij, tbar_abcijk, omega = zero, cvs=.false.)
!
      call mem%dealloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     Contribution of the singles part
!
!     :: vo-part of the density ::
!
      call mem%alloc(density_vo, wf%n_v, wf%n_o)
      call zero_array(density_vo, wf%n_o*wf%n_v)
!
      call wf%debug_right_vo(density_vo, R_abij, tbar_abcijk)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            r_tdm(wf%n_o + a, i) = r_tdm(wf%n_o + a, i) + density_vo(a, i)
         enddo
      enddo
!
!     :: ov-part of the density ::
!
      call mem%alloc(t_abcijk, wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
      call wf%construct_full_t3(t_abcijk)
!
      call mem%alloc(density_ov, wf%n_o, wf%n_v)
      call zero_array(density_ov, wf%n_o*wf%n_v)
!
      call wf%debug_right_ov_t3(density_ov, density_vo, R_ai, tbar_abcijk, t_abcijk)
!
      call mem%dealloc(density_vo, wf%n_v, wf%n_o)
!
      call mem%dealloc(t_abcijk, wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
!
      call wf%debug_right_ov_Y_term(density_ov, R_abij, tbar_abcijk)
!
!     :: Construct  R3 ::
      call mem%alloc(R_abcijk, wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
      call wf%construct_full_R3(R_abij, R_abcijk, wf%right_excitation_energies(state))
!
      call wf%debug_right_ov_R3(density_ov, R_abcijk)
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
            r_tdm(i, wf%n_o + a) = r_tdm(i, wf%n_o + a) + density_ov(i, a)
         enddo
      enddo
!
      call mem%dealloc(density_ov, wf%n_o, wf%n_v)
!
!     :: oo-part of the density ::
!
!     rho_lk -= sum_abcij tbar^abc_ijk*R^a_i*t^bc_jl 
!     rho_lk -= half sum_abcij tbar^abc_ijk*R^abc_ijl*Δ_ai,bj,cl
!
      call mem%alloc(density_oo, wf%n_o, wf%n_o)
      call zero_array(density_oo, wf%n_o**2)
!
      call wf%debug_right_oo(density_oo, R_ai, tbar_abcijk, R_abcijk)
!
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            r_tdm(i, j) = r_tdm(i, j) + density_oo(i, j)
         enddo
      enddo
!
      call mem%dealloc(density_oo, wf%n_o, wf%n_o)
!
!     :: vv-part of the density
!     rho_cd += sum_abijk tbar^abc_ijk*R^a_i*t^bd_jk 
!     rho_cd += half sum_abcijk tbar^abc_ijk*R^abd_ijk*Δ_ai,bj,dk
!
      call mem%alloc(density_vv, wf%n_v, wf%n_v)
      call zero_array(density_vv, wf%n_v**2)
!
      call wf%debug_right_vv(density_vv, R_ai, tbar_abcijk, R_abcijk)
!
      do b = 1, wf%n_v
         do a = 1, wf%n_v
            r_tdm(wf%n_o + a, wf%n_o + b) = r_tdm(wf%n_o + a, wf%n_o + b) + density_vv(a, b)
         enddo
      enddo
!
      call mem%dealloc(density_vv, wf%n_v, wf%n_v)
!
!     tbar3 R3 overlap
      tbar_R_overlap = tbar_R_overlap + sixth*ddot((wf%n_v*wf%n_o)**3, tbar_abcijk, 1 , R_abcijk, 1)
!
      call mem%dealloc(R_abcijk,wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
      call mem%dealloc(tbar_abcijk,wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o)
!
      call mem%dealloc(R_ai, wf%n_v, wf%n_o)
!
!     Contribution of the ground state density scaled by 
!     the right-hand side reference term (- sum_mu tbar_mu*R_mu)
!
      call wf%density_mu_mu_oo(r_tdm, tbar_R_overlap)
!
      call wf%density_mu_ref(r_tdm,          &
                             wf%density,     &
                             tbar_R_overlap)
!
   end subroutine right_tdm_debug_cc3
!
!
   module subroutine construct_full_R3_cc3(wf, R_abij, R_abcijk, omega)
!!
!!    Construct full R3 debug
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    constructs full v3o3 array containing the triples excitation vector 
!!    using unwrapped loops
!!
!!    R^abc = (ω - ε^abc_ijk)^-1 * P^abc_ijk (sum_d c^ad_ij g_ckbd - sum_l c^ab_il g_cklj
!!             + sum_d t^ad_ij g'_bdck - sum_l t^ab_il g'_cklj
!!
!!    where g'_ are c1-transformed integrals
!!
!!     NB: the amplitudes are scaled by the biorthonormal factor
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_o,wf%n_o), intent(in) :: R_abij
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(inout) :: R_abcijk
!
      real(dp), intent(in) :: omega
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     integrals
      real(dp), dimension(:,:,:,:), allocatable :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable :: g_ljci
      real(dp), dimension(:,:,:,:), allocatable :: g_bdci_c1
      real(dp), dimension(:,:,:,:), allocatable :: g_ljci_c1
!
      integer :: i,j,k,l,a,b,c,d
!
!     CVS
      integer :: i_cvs
      logical :: ijk_core
!
!     Prepare t2 amplitudes and read in integrals
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
      call wf%g_bdck_c1%open_('read')
      call wf%g_ljck_c1%open_('read')
!
      call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call single_record_reader(wf%n_o, wf%g_bdck_t, g_bdci)
!
      call mem%alloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call single_record_reader(wf%n_o, wf%g_bdck_c1, g_bdci_c1)
!
      call zero_array(R_abcijk, (wf%n_v*wf%n_o)**3)
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. i .eq. k) cycle
!
!              Check if at least one index i,j,k is a core orbital
               if(wf%cvs) then
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
!              Cycle if i,j,k are not core orbitals
               if (.not. ijk_core) cycle
!
               end if
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
                        do d = 1, wf%n_v
!
                           R_abcijk(a,b,c,i,j,k) = R_abcijk(a,b,c,i,j,k)   &
                           + R_abij(a,d,i,j)*g_bdci(d,b,c,k)               &
                           + R_abij(b,d,j,i)*g_bdci(d,a,c,k)               &
                           + R_abij(c,d,k,j)*g_bdci(d,b,a,i)               &
                           + R_abij(a,d,i,k)*g_bdci(d,c,b,j)               &
                           + R_abij(b,d,j,k)*g_bdci(d,c,a,i)               &
                           + R_abij(c,d,k,i)*g_bdci(d,a,b,j)               &
!
                           + t_abij(a,d,i,j)*g_bdci_c1(d,b,c,k)            &
                           + t_abij(b,d,j,i)*g_bdci_c1(d,a,c,k)            &
                           + t_abij(c,d,k,j)*g_bdci_c1(d,b,a,i)            &
                           + t_abij(a,d,i,k)*g_bdci_c1(d,c,b,j)            &
                           + t_abij(b,d,j,k)*g_bdci_c1(d,c,a,i)            &
                           + t_abij(c,d,k,i)*g_bdci_c1(d,a,b,j)
!
                        end do
                     end do
                  end do
               end do
            end do 
         end do
      end do
!
      call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(g_bdci_c1, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call wf%g_bdck_t%close_()
      call wf%g_bdck_c1%close_()
!
      call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call compound_record_reader(wf%n_o, wf%n_o, wf%g_ljck_t, g_ljci)
!
      call mem%alloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call compound_record_reader(wf%n_o, wf%n_o, wf%g_ljck_c1, g_ljci_c1)
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. i .eq. k) cycle
!
!              Check if at least one index i,j,k is a core orbital
               if(wf%cvs) then
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
!              Cycle if i,j,k are not core orbitals
               if (.not. ijk_core) cycle
!
               end if
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
                        do l = 1, wf%n_o
!
                           R_abcijk(a,b,c,i,j,k) = R_abcijk(a,b,c,i,j,k)   &
                           - R_abij(a,b,i,l)*g_ljci(l,c,j,k)               &
                           - R_abij(b,a,j,l)*g_ljci(l,c,i,k)               &
                           - R_abij(c,b,k,l)*g_ljci(l,a,j,i)               &
                           - R_abij(a,c,i,l)*g_ljci(l,b,k,j)               &
                           - R_abij(b,c,j,l)*g_ljci(l,a,k,i)               &
                           - R_abij(c,a,k,l)*g_ljci(l,b,i,j)               &
!
                           - t_abij(a,b,i,l)*g_ljci_c1(l,c,j,k)            &
                           - t_abij(b,a,j,l)*g_ljci_c1(l,c,i,k)            &
                           - t_abij(c,b,k,l)*g_ljci_c1(l,a,j,i)            &
                           - t_abij(a,c,i,l)*g_ljci_c1(l,b,k,j)            &
                           - t_abij(b,c,j,l)*g_ljci_c1(l,a,k,i)            &
                           - t_abij(c,a,k,l)*g_ljci_c1(l,b,i,j)
!
                        end do
!
                        R_abcijk(a,b,c,i,j,k) = R_abcijk(a,b,c,i,j,k)   &
                        * one/(omega - wf%orbital_energies(wf%n_o + a)  &
                                     - wf%orbital_energies(wf%n_o + b)  &
                                     - wf%orbital_energies(wf%n_o + c)  &
                                     + wf%orbital_energies(i)           &
                                     + wf%orbital_energies(j)           &
                                     + wf%orbital_energies(k))
!
                        R_abcijk(a,b,c,i,i,i) = zero
!
                     end do
                  end do
!
                  R_abcijk(c,c,c,i,j,k) = zero
!
               end do
            end do
         end do
      end do
!
!     Biorthonormal factor
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            if (i .eq. j) cycle
!
            do b = 1, wf%n_v
               do a = 1, wf%n_v
!
                  if(a .eq. b) cycle
!
                  R_abcijk(a,a,b,i,i,j) = two*R_abcijk(a,a,b,i,i,j)
                  R_abcijk(a,b,a,i,j,i) = two*R_abcijk(a,b,a,i,j,i)
                  R_abcijk(b,a,a,j,i,i) = two*R_abcijk(b,a,a,j,i,i)
!
               end do
            end do
         end do
      end do
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_ljci_c1, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%g_ljck_t%close_()
      call wf%g_ljck_c1%close_()
!
   end subroutine construct_full_R3_cc3
!
!
   module subroutine construct_full_t3_cc3(wf, t_abcijk)
!!
!!    Construct full t3
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    constructs full v3o3 array containing 
!!    the triples ground state amplitudes
!!
!!    Contributions to W
!!    W^abc_ijk = P^abc_ijk(sum_d t^ad_ij g_bdck - sum_l t^ab_il g_ljck)
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout) :: t_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     help array
      real(dp), dimension(:,:,:), allocatable :: u_abc
!
!     integrals
      real(dp), dimension(:,:,:,:), allocatable :: g_bdci
      real(dp), dimension(:,:,:,:), allocatable :: g_ljci
!
      integer :: i,j,k
!
!     Prepare t2 amplitudes and read in integrals
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call wf%g_bdck_t%open_('read')
      call wf%g_ljck_t%open_('read')
!
      call mem%alloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call single_record_reader(wf%n_o, wf%g_bdck_t, g_bdci)
!
      call mem%alloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call compound_record_reader(wf%n_o, wf%n_o, wf%g_ljck_t, g_ljci)
!
      call mem%alloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call zero_array(t_abcijk, (wf%n_v*wf%n_o)**3)
!
!     Construct full t3 in analogy to the loops in omega_cc3
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
            do k = 1, wf%n_o
!
               if (i .eq. j .and. i .eq. k) cycle
!
               call wf%omega_cc3_W_calc(i, j, k,               &
                                       t_abcijk(:,:,:,i,j,k),  &
                                       u_abc, t_abij,          &
                                       g_bdci(:,:,:,i),        &
                                       g_bdci(:,:,:,j),        &
                                       g_bdci(:,:,:,k),        &
                                       g_ljci(:,:,j,i),        &
                                       g_ljci(:,:,k,i),        &
                                       g_ljci(:,:,k,j),        &
                                       g_ljci(:,:,i,j),        &
                                       g_ljci(:,:,i,k),        &
                                       g_ljci(:,:,j,k))
!
                  call wf%omega_cc3_eps(i, j, k, t_abcijk(:,:,:,i,j,k))
!
            end do 
         end do
      end do
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(u_abc, wf%n_v, wf%n_v, wf%n_v)
!
      call mem%dealloc(g_bdci, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(g_ljci, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
!
      call wf%g_bdck_t%close_()
      call wf%g_ljck_t%close_()
!
   end subroutine construct_full_t3_cc3
!
!
   module subroutine construct_full_tbar3_cc3(wf, tbar_ai, tbar_abij, tbar_abcijk, omega, cvs)
!!
!!    Construct full tbar3
!!    Written by Alexander C. Paul, Sep 2019
!!
!!    constructs full v3o3 array containing the triples multipliers
!!    using unwrapped loops
!!
!!    tbar_^abc_ijk 
!!    = (ω - ε^abc_ijk)^-1 P^abc_ijk (tbar__ai*L_jbkc - tbar__ak*L_jbic 
!!                                  + tbar_abij*F_kc - tbar__abik*F_jc
!!                                  + sum_l (tbar__ablk g_iljc - tbar__abil L_jlkc) 
!!                                  - sum_d (tbar__adjk g_ibdc - tbar__adij L_dbkc))
!!
!!    NB: can also be used for left excitation vector
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: tbar_ai
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: tbar_abij
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(inout) :: tbar_abcijk
!
      real(dp), intent(in) :: omega
!
      logical, intent(in) :: cvs
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     integrals
      real(dp), dimension(:,:,:,:), allocatable :: g_dbic
      real(dp), dimension(:,:,:,:), allocatable :: g_jlic
      real(dp), dimension(:,:,:,:), allocatable :: L_ibjc
!
      integer :: i,j,k,l,a,b,c,d
!
!     CVS
      integer :: i_cvs
      logical :: ijk_core
!
!     Prepare t2 amplitudes and read in integrals
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!     Integrals
      call wf%g_dbkc_t%open_('read')
      call wf%g_jlkc_t%open_('read')
      call wf%L_jbkc_t%open_('read')
!
      call mem%alloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call single_record_reader(wf%n_o, wf%g_dbkc_t, g_dbic)
!
      call mem%alloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call compound_record_reader(wf%n_o, wf%n_o, wf%g_jlkc_t, g_jlic)
!
      call mem%alloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call compound_record_reader(wf%n_o, wf%n_o, wf%L_jbkc_t, L_ibjc)
!
!
      call zero_array(tbar_abcijk, (wf%n_v*wf%n_o)**3)
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. i .eq. k) cycle
!
!              Check if at least one index i,j,k is a core orbital
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
!              Cycle if i,j,k are not core orbitals
               if (.not. ijk_core) cycle
!
               end if
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
                        do l = 1, wf%n_o
!
                           tbar_abcijk(a,b,c,i,j,k)                  &
                           = tbar_abcijk(a,b,c,i,j,k)                &
                           + tbar_abij(a,b,k,l)*g_jlic(c,l,j,i)      &
                           + tbar_abij(b,a,k,l)*g_jlic(c,l,i,j)      &
                           + tbar_abij(c,b,i,l)*g_jlic(a,l,j,k)      &
                           + tbar_abij(a,c,j,l)*g_jlic(b,l,k,i)      &
                           + tbar_abij(b,c,i,l)*g_jlic(a,l,k,j)      &
                           + tbar_abij(c,a,j,l)*g_jlic(b,l,i,k)      &
!
                           - two*tbar_abij(a,b,i,l)*g_jlic(c,l,j,k)  &
                           - two*tbar_abij(b,a,j,l)*g_jlic(c,l,i,k)  &
                           - two*tbar_abij(c,b,k,l)*g_jlic(a,l,j,i)  &
                           - two*tbar_abij(a,c,i,l)*g_jlic(b,l,k,j)  &
                           - two*tbar_abij(b,c,j,l)*g_jlic(a,l,k,i)  &
                           - two*tbar_abij(c,a,k,l)*g_jlic(b,l,i,j)  &
!
                           + tbar_abij(a,b,i,l)*g_jlic(c,l,k,j)      &
                           + tbar_abij(b,a,j,l)*g_jlic(c,l,k,i)      &
                           + tbar_abij(c,b,k,l)*g_jlic(a,l,i,j)      &
                           + tbar_abij(a,c,i,l)*g_jlic(b,l,j,k)      &
                           + tbar_abij(b,c,j,l)*g_jlic(a,l,i,k)      &
                           + tbar_abij(c,a,k,l)*g_jlic(b,l,j,i)
!
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. i .eq. k) cycle
!
!              Check if at least one index i,j,k is a core orbital
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
!              Cycle if i,j,k are not core orbitals
               if (.not. ijk_core) cycle
!
               end if
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
                        do d = 1, wf%n_v
!
                           tbar_abcijk(a,b,c,i,j,k)                  &
                           = tbar_abcijk(a,b,c,i,j,k)                &
                           - tbar_abij(a,d,k,j)*g_dbic(b,c,d,i)      &
                           - tbar_abij(b,d,k,i)*g_dbic(a,c,d,j)      &
                           - tbar_abij(c,d,i,j)*g_dbic(b,a,d,k)      &
                           - tbar_abij(a,d,j,k)*g_dbic(c,b,d,i)      &
                           - tbar_abij(b,d,i,k)*g_dbic(c,a,d,j)      &
                           - tbar_abij(c,d,j,i)*g_dbic(a,b,d,k)      &
!
                           + two*tbar_abij(a,d,i,j)*g_dbic(b,c,d,k)  &
                           + two*tbar_abij(b,d,j,i)*g_dbic(a,c,d,k)  &
                           + two*tbar_abij(c,d,k,j)*g_dbic(b,a,d,i)  &
                           + two*tbar_abij(a,d,i,k)*g_dbic(c,b,d,j)  &
                           + two*tbar_abij(b,d,j,k)*g_dbic(c,a,d,i)  &
                           + two*tbar_abij(c,d,k,i)*g_dbic(a,b,d,j)  &
!
                           - tbar_abij(a,d,i,k)*g_dbic(b,c,d,j)      &
                           - tbar_abij(b,d,j,k)*g_dbic(a,c,d,i)      &
                           - tbar_abij(c,d,k,i)*g_dbic(b,a,d,j)      &
                           - tbar_abij(a,d,i,j)*g_dbic(c,b,d,k)      &
                           - tbar_abij(b,d,j,i)*g_dbic(c,a,d,k)      &
                           - tbar_abij(c,d,k,j)*g_dbic(a,b,d,i)
                           
!
                        end do
                     end do
                  end do
               end do
            end do
         end do
      end do
!
!     direct products
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. i .eq. k) cycle
!
!              Check if at least one index i,j,k is a core orbital
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
!              Cycle if i,j,k are not core orbitals
               if (.not. ijk_core) cycle
!
               end if
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
!
                        tbar_abcijk(a,b,c,i,j,k) =   tbar_abcijk(a,b,c,i,j,k)       &
                                                   + tbar_ai(a,i)*L_ibjc(b,c,j,k)   &
                                                   + tbar_ai(b,j)*L_ibjc(a,c,i,k)   &
                                                   + tbar_ai(c,k)*L_ibjc(b,a,j,i)   &
                                                   + tbar_ai(a,i)*L_ibjc(c,b,k,j)   &
                                                   + tbar_ai(b,j)*L_ibjc(c,a,k,i)   &
                                                   + tbar_ai(c,k)*L_ibjc(a,b,i,j)   &
!
                                                   - tbar_ai(a,k)*L_ibjc(b,c,j,i)   &
                                                   - tbar_ai(b,k)*L_ibjc(a,c,i,j)   &
                                                   - tbar_ai(c,i)*L_ibjc(b,a,j,k)   &
                                                   - tbar_ai(a,j)*L_ibjc(c,b,k,i)   &
                                                   - tbar_ai(b,i)*L_ibjc(c,a,k,j)   &
                                                   - tbar_ai(c,j)*L_ibjc(a,b,i,k)   &
!
                                                   + tbar_abij(a,b,i,j)*wf%fock_ia(k,c)   &
                                                   + tbar_abij(b,a,j,i)*wf%fock_ia(k,c)   &
                                                   + tbar_abij(a,c,i,k)*wf%fock_ia(j,b)   &
                                                   + tbar_abij(c,b,k,j)*wf%fock_ia(i,a)   &
                                                   + tbar_abij(b,c,j,k)*wf%fock_ia(i,a)   &
                                                   + tbar_abij(c,a,k,i)*wf%fock_ia(j,b)   &
!
                                                   - tbar_abij(a,b,i,k)*wf%fock_ia(j,c)   &
                                                   - tbar_abij(b,a,j,k)*wf%fock_ia(i,c)   &
                                                   - tbar_abij(a,c,i,j)*wf%fock_ia(k,b)   &
                                                   - tbar_abij(c,b,k,i)*wf%fock_ia(j,a)   &
                                                   - tbar_abij(b,c,j,i)*wf%fock_ia(k,a)   &
                                                   - tbar_abij(c,a,k,j)*wf%fock_ia(i,b)
!
                        tbar_abcijk(a,b,c,i,j,k) = tbar_abcijk(a,b,c,i,j,k)            &
                                                   * one/(omega                        &
                                                   - wf%orbital_energies(wf%n_o + a)   &
                                                   - wf%orbital_energies(wf%n_o + b)   &
                                                   - wf%orbital_energies(wf%n_o + c)   &
                                                   + wf%orbital_energies(i)            &
                                                   + wf%orbital_energies(j)            &
                                                   + wf%orbital_energies(k))
!
                        tbar_abcijk(a,b,c,i,i,i) = zero
!
                     end do
                  end do
!
                  tbar_abcijk(c,c,c,i,j,k) = zero
!
               end do
            end do
         end do
      end do
!
      call wf%g_dbkc_t%close_()
      call wf%g_jlkc_t%close_()
      call wf%L_jbkc_t%close_()
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
      call mem%dealloc(g_dbic, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
      call mem%dealloc(g_jlic, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call mem%dealloc(L_ibjc, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine construct_full_tbar3_cc3
!
!
   module subroutine debug_left_oo_cc3(wf, density_oo, L_abcijk, t_abcijk)
!!
!!    CC3 contribution to the oo-block of the left tdm
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following term using unwrapped loops:
!!
!!       rho_kl -= 1/2 sum_abcijk L^abc_ijk*t^abc_ijk
!!
!!    expects full v3o3 arrays for the amplitudes
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: L_abcijk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: t_abcijk
!
      integer :: a, b, c, i, j, k, l
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           density_oo(k,l) = density_oo(k,l) &
                           - half*L_abcijk(a,b,c,i,j,l)*t_abcijk(a,b,c,i,j,k)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
   end subroutine debug_left_oo_cc3
!
!
   module subroutine debug_left_ov_N7_cc3(wf, density_ov, L_abcijk)
!!
!!    CC3 contribution to the ov-block of the left tdm
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following term using unwrapped loops:
!!
!!       rho_ld -= sum_abcijk L^abc_ijk*t^ad_ik*t^bc_jl
!!       rho_ld -= sum_aik Y_alik*t^ad_ik
!!
!!    expects full v3o3 array for the left amplitudes
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: L_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
      real(dp), dimension(:,:,:,:), allocatable :: Y_alik
!
      integer :: a, b, c, i, j, k, l, d
!
      call mem%alloc(Y_alik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(Y_alik, wf%n_v*wf%n_o**3)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           Y_alik(a,l,i,k) = Y_alik(a,l,i,k) &
                           + L_abcijk(a,b,c,i,j,k)*t_abij(b,c,j,l)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o 
            do a = 1, wf%n_v
               do i = 1, wf%n_o
                  do k = 1, wf%n_o
                     density_ov(l,d) = density_ov(l,d) &
                     - Y_alik(a,l,i,k)*t_abij(a,d,i,k)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(Y_alik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine debug_left_ov_N7_cc3
!
!
   module subroutine debug_left_ov_N6_cc3(wf, density_ov, L_abij, t_abcijk)
!!
!!    CC3 contribution to the ov-block of the left tdm N6 scaling
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following term using unwrapped loops:
!!    expects full v3o3 array for the multipliers
!!
!!       rho_ld += sum_abcijk L^ab_ij*(t^abd_ijl-t^adb_ijl)
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: L_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: t_abcijk
!
      integer :: a, b, d, i, j, l
!
      do l = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. j .eq. l) cycle
!
               do d = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. d) cycle
!
                        density_ov(l,d) = density_ov(l,d) + &
                        L_abij(a,b,i,j)*(t_abcijk(a,b,d,i,j,l)-t_abcijk(a,d,b,i,j,l))
                        
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
   end subroutine debug_left_ov_N6_cc3
!
!
   module subroutine debug_left_vv_cc3(wf, density_vv, L_abcijk, t_abcijk)
!!
!!    CC3 contribution to the vv-block of the left tdm
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following term using unwrapped loops:
!!
!!       rho_cd += 1/2 sum_abcijk L^abc_ijk*t^abd_ijk
!!
!!    expects full v3o3 arrays for the amplitudes
!!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: L_abcijk
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: t_abcijk
!
      integer :: a, b, c, i, j, k, d
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           density_vv(c,d) = density_vv(c,d) &
                           + half*L_abcijk(a,b,c,i,j,k)*t_abcijk(a,b,d,i,j,k)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
   end subroutine debug_left_vv_cc3
!
!
   module subroutine debug_right_vo_cc3(wf, density_vo, R_abij, tbar_abcijk)
!!
!!    CC3 contribution to the vo-block of the right tdm
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following term using unwrapped loops:
!!
!!       rho_ck += sum_abij 1/2 tbar^abc_ijk*R^ab_ij
!!
!!    expects full v3o3 array for the multipliers
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: density_vo
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: tbar_abcijk
!
      integer :: i, j, k, a, b, c
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. j .eq. k) cycle
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
                        density_vo(c,k) = density_vo(c,k) &
                        + half*tbar_abcijk(a,b,c,i,j,k)*R_abij(a,b,i,j)
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
   end subroutine debug_right_vo_cc3
!
!
   module subroutine debug_right_ov_t3_cc3(wf, density_ov, density_vo, R_ai, tbar_abcijk, t_abcijk)
!!
!!    CC3 contribution to the ov-block of the right tdm (t3 contribution)
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following terms using unwrapped loops:
!!
!!       rho_ld += -1/2 tbar^abc_ijk*t^abc_ijl*R^d_k
!!
!!       rho_ld += sum_abcijk 1/2 tbar^abc_ijk R^ab_ij (2t^cd_kl - t^cd_lk)
!!       rho_ld += sum_ck rho_ck (2t^cd_kl - t^cd_lk)
!!
!!       rho_ld += sum_abcijk tbar^abc_ijk*R^a_i*(t^bcd_jkl - t^bcd_jlk)
!!       rho_ld += sum_bcjk Z^bc_jk*(t^bcd_jkl - t^bcd_jlk)
!!
!!    expects full v3o3 arrays for the triples amplitudes
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: density_vo
!
      real(dp), dimension(wf%n_v, wf%n_o),intent(in) :: R_ai
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, wf%n_o, wf%n_o, wf%n_o), intent(in) :: tbar_abcijk, t_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Intermediates
      real(dp), dimension(:,:), allocatable :: gs_oo, gs_vv
      real(dp), dimension(:,:,:,:), allocatable :: Z_bcjk
!
      integer :: a,b,c,i,j,k,l,d
!
!     Construct GS_density_oo as intermediate then contract with R^d_k
!
      call mem%alloc(gs_oo, wf%n_o, wf%n_o)
      call zero_array(gs_oo, wf%n_o**2)
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           gs_oo(l,k) = gs_oo(l,k) &
                           - half*tbar_abcijk(a,b,c,i,j,k)*t_abcijk(a,b,c,i,j,l)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!     Contraction with R^d_k
      do d = 1, wf%n_v
         do l = 1, wf%n_o
            do k = 1, wf%n_o
               density_ov(l,d) = density_ov(l,d) + gs_oo(l,k)*R_ai(d,k)
            enddo
         enddo
      enddo
!
      call mem%dealloc(gs_oo, wf%n_o, wf%n_o)
!
!     rho_ld += -1/2 tbar^abc_ijk*t^abc_ijl*R^d_k
!     Construct GS_density_oo as intermediate then contract with R^d_k
!
      call mem%alloc(gs_vv, wf%n_v, wf%n_v)
      call zero_array(gs_vv, wf%n_v**2)
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           gs_vv(c,d) = gs_vv(c,d) &
                           + half*tbar_abcijk(a,b,c,i,j,k)*t_abcijk(a,b,d,i,j,k)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
!     Contraction with R^d_k
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o
            do c = 1, wf%n_v
               density_ov(l,d) = density_ov(l,d) - gs_vv(c,d)*R_ai(c,l)
            enddo
         enddo
      enddo
!
      call mem%dealloc(gs_vv, wf%n_v, wf%n_v)
!
!     rho_ld += sum_abcijk 1/2 tbar^abc_ijk R^ab_ij (2t^cd_kl - t^cd_lk)
!     rho_ld += sum_ck rho_ck (2t^cd_kl - t^cd_lk)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o 
            do c = 1, wf%n_v
               do k = 1, wf%n_o
                  density_ov(l,d) = density_ov(l,d) &
                  + density_vo(c,k)*(two*t_abij(c,d,k,l) - t_abij(c,d,l,k))
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
!     rho_ld += sum_abcijk tbar^abc_ijk*R^a_i*(t^bcd_jkl - t^bcd_jlk)
!     rho_ld += sum_bcjk Z^bc_jk*(t^bcd_jkl - t^bcd_jlk)
!
      call mem%alloc(Z_bcjk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call zero_array(Z_bcjk, (wf%n_v*wf%n_o)**2)
!
      do k = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (i .eq. j .and. j .eq. k) cycle
!
               do c = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (a .eq. b .and. b .eq. c) cycle
!
                        Z_bcjk(b,c,j,k) = Z_bcjk(b,c,j,k) &
                        + tbar_abcijk(a,b,c,i,j,k)*R_ai(a,i)
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do j = 1, wf%n_o
!
               if (j .eq. k .and. k .eq. l) cycle
!
               do d = 1, wf%n_v
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
!
                        if (b .eq. c .and. c .eq. d) cycle
!
                        density_ov(l,d) = density_ov(l,d) + &
                        Z_bcjk(b,c,j,k)*(t_abcijk(b,c,d,j,k,l)-t_abcijk(b,c,d,j,l,k))
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(Z_bcjk, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine debug_right_ov_t3_cc3
!
!
   module subroutine debug_right_ov_Y_term_cc3(wf, density_ov, R_abij, tbar_abcijk)
!!
!!    CC3 contribution to the ov-block of the right tdm (Y intermediates)
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following terms using unwrapped loops:
!!
!!       rho_ld -= sum_abcijk tbar^abc_ijk*R^ad_ik*t^bc_jl
!!       rho_ld -= sum_aik Y_alik*R^ad_ik
!!
!!    expects full v3o3 array for the multipliers
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R_abij
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
!     Intermediates
      real(dp), dimension(:,:,:,:), allocatable :: Y_alik, Y_acdi
!
      integer :: a,b,c,i,j,k,l,d
!
      call mem%alloc(Y_alik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call zero_array(Y_alik, wf%n_v*wf%n_o**3)
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           Y_alik(a,l,i,k) = Y_alik(a,l,i,k) &
                           + tbar_abcijk(a,b,c,i,j,k)*t_abij(b,c,j,l)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o 
            do a = 1, wf%n_v
               do i = 1, wf%n_o
                  do k = 1, wf%n_o
                     density_ov(l,d) = density_ov(l,d) &
                     - Y_alik(a,l,i,k)*R_abij(a,d,i,k)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(Y_alik, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     rho_ld -= sum_abcijk tbar^abc_ijk*R^ac_il*t^bd_jk
!     rho_ld -= sum_aik Y_alik*R^ac_il
!
      call mem%alloc(Y_acdi, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
      call zero_array(Y_acdi, wf%n_o*wf%n_v**3)
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           Y_acdi(a,c,d,i) = Y_acdi(a,c,d,i) &
                           + tbar_abcijk(a,b,c,i,j,k)*t_abij(b,d,j,k)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      do d = 1, wf%n_v
         do l = 1, wf%n_o 
            do c = 1, wf%n_v
               do a = 1, wf%n_v
                  do i = 1, wf%n_o
                     density_ov(l,d) = density_ov(l,d) &
                     - Y_acdi(a,c,d,i)*R_abij(a,c,i,l)
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(Y_acdi, wf%n_v, wf%n_v, wf%n_v, wf%n_o)
!
   end subroutine debug_right_ov_Y_term_cc3
!
!
   module subroutine debug_right_ov_R3_cc3(wf, density_ov, R_abcijk)
!!
!!    CC3 contribution to the ov-block of the right tdm (R3 contribution)
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following terms using unwrapped loops:
!!
!!       rho_ld += sum_abij tbar^ab_ij*(R^abd_ijl - R^abd_ilj)
!!
!!    expects the full v3o3 arrays for the amplitudes
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(inout) :: density_ov
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: R_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: tbar_abij
!
      integer :: a,b,d,i,j,l
!
      call mem%alloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2bar, tbar_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do l = 1, wf%n_o
         do j = 1, wf%n_o
            do i = 1, wf%n_o
!
               if (j .eq. l .and. i .eq. j) cycle
!
               do d = 1, wf%n_v
                  do b = 1, wf%n_v
                     do a = 1, wf%n_v
!
                        if (d .eq. b .and. a .eq. b) cycle
!
                        density_ov(l,d) = density_ov(l,d) + &
                        tbar_abij(a,b,i,j)*(R_abcijk(a,b,d,i,j,l)-R_abcijk(a,b,d,i,l,j))
                        
!
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(tbar_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine debug_right_ov_R3_cc3
!
!
   module subroutine debug_right_oo_cc3(wf, density_oo, R_ai, tbar_abcijk, R_abcijk)
!!
!!    CC3 contribution to the oo-block of the right tdm
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following terms using unwrapped loops:
!!
!!       rho_lk -= sum_abcij (tbar^abc_ijk*R^a_i*t^bc_jl 
!!                          + half tbar^abc_ijk*R^abc_ijl*delta_ai,bj,cl)
!!
!!    expects the full v3o3 arrays for the amplitudes
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: density_oo
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: R_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      integer :: i,j,k,l,a,b,c
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do l = 1, wf%n_o
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           density_oo(l,k) = density_oo(l,k) &
                           - tbar_abcijk(a,b,c,i,j,k)*R_ai(a,i)*t_abij(b,c,j,l)
!
                           density_oo(l,k) = density_oo(l,k) &
                           - half*tbar_abcijk(a,b,c,i,j,k)*R_abcijk(a,b,c,i,j,l)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine debug_right_oo_cc3
!
!
   module subroutine debug_right_vv_cc3(wf, density_vv, R_ai, tbar_abcijk, R_abcijk)
!!
!!    CC3 contribution to the vv-block of the right tdm
!!    Written by Alexander C. Paul, September 2019
!!
!!    Constructs the following terms using unwrapped loops:
!!
!!       rho_cd += sum_abijk tbar^abc_ijk*R^a_i*t^bd_jk 
!!       rho_cd += half sum_abcijk tbar^abc_ijk*R^abd_ijk*Δ_ai,bj,dk
!!
!!    expects the full v3o3 arrays for the amplitudes
!!
      implicit none
!
      class(cc3) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: density_vv
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: tbar_abcijk
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v,wf%n_o,wf%n_o,wf%n_o), intent(in) :: R_abcijk
!
      real(dp), dimension(:,:,:,:), allocatable :: t_abij
!
      integer :: i,j,k,a,b,c,d
!
      call mem%alloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call squareup_and_sort_1234_to_1324(wf%t2, t_abij, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      do d = 1, wf%n_v
         do k = 1, wf%n_o
            do j = 1, wf%n_o
               do i = 1, wf%n_o
!
                  if (i .eq. j .and. j .eq. k) cycle
!
                  do c = 1, wf%n_v
                     do b = 1, wf%n_v
                        do a = 1, wf%n_v
!
                           if (a .eq. b .and. b .eq. c) cycle
!
                           density_vv(c,d) = density_vv(c,d) &
                           + tbar_abcijk(a,b,c,i,j,k)*R_ai(a,i)*t_abij(b,d,j,k)
!
                           density_vv(c,d) = density_vv(c,d) &
                           + half*tbar_abcijk(a,b,c,i,j,k)*R_abcijk(a,b,d,i,j,k)
!
                        enddo
                     enddo
                  enddo
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(t_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
!
   end subroutine debug_right_vv_cc3
!
!
end submodule debug_transition_density_cc3