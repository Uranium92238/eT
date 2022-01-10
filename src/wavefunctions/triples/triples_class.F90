!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
module triples_class
!
!!
!!    Abstract triples class module
!!    Written by Rolf H. Myhre and Alexander C. Paul, 2018-2019
!!
!!    Triples class containing W_calc - intermediate for triples amplitueds -
!!    used in CCSD(T), CC3 omega and CC3 jacobian.
!!    Also contains the routine dividing by the orbital energy differences
!!
!
   use ccsd_class, only: ccsd
   use parameters
!

   use eri_cholesky_disk_class, only: eri_cholesky_disk
   use abstract_eri_cholesky_class, only: abstract_eri_cholesky
   use eri_1idx_transformed_tool_class, only: eri_1idx_transformed_tool
   use eri_tool_class, only: eri_tool
   use eri_adapter_class, only: eri_adapter
!
   use memory_manager_class, only: mem
!
   implicit none
!
   type, abstract, extends(ccsd) :: triples
!
   contains
!
      procedure :: construct_W                    => construct_W_triples
      procedure :: construct_V                    => construct_V_triples
      procedure :: construct_W_permutation        => construct_W_permutation_triples
!
      procedure :: divide_by_orbital_differences  => divide_by_orbital_differences_triples
!
      procedure :: setup_vvvo                     => setup_vvvo_triples
      procedure :: setup_vvov                     => setup_vvov_triples
      procedure :: setup_ovov                     => setup_ovov_triples
      procedure :: setup_oovo                     => setup_oovo_triples
      procedure :: setup_ooov                     => setup_ooov_triples
!
      procedure :: point_vvvo                     => point_vvvo_triples
      procedure :: point_vooo                     => point_vooo_triples
      procedure :: point_vvoo                     => point_vvoo_triples
!
!
      procedure :: construct_c1_cholesky     => construct_cholesky_c1_triples
      procedure :: construct_cholesky_c1_oo  => construct_cholesky_c1_oo_triples
      procedure :: construct_cholesky_c1_vo  => construct_cholesky_c1_vo_triples
      procedure :: construct_cholesky_c1_vv  => construct_cholesky_c1_vv_triples
!
   end type triples
!
!
contains
!
!
   subroutine construct_W_triples(wf, i, j, k, u_abc, v_abc, t_abij, &
                                  g_bdci, g_bdcj, g_bdck, &
                                  g_ljci, g_lkci, g_lkcj, g_licj, g_lick, g_ljck)
!!
!!    Construct W
!!    Written by Rolf H. Myhre and Alexander C. Paul, Oct 2020
!!
!!    Calculate the intermediate W for occupied indices i,j,k
!!
!!    Contributions to W
!!    W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il (lj|ck))
!!
!!    Each permutation is obtained by a call to construct_W_permutation
!!
!!    T^abc_ijk = W^abc_ijk/eps^abc_ijk
!!
!!    Note: v_abc contains W^abc in the end
!!
      use reordering, only: sort_123_to_312, sort_123_to_213
!
      implicit none
!
      class(triples) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t_abij
!
!     (bd,ci) ordered dbc,i
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdck
!
!     (lj,ck) ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ljci
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lkci
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lkcj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_licj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lick
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ljck
!
!     u_bac terms
!     -----------
!
!     t^bd_ji*(ad|ck)- t^ac_lk*(li|bj)
!     --------------------------------
!
      call wf%construct_W_permutation(j, i, k, u_abc, t_abij, g_bdck, g_licj, zero)
!
!     bac -> cba
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     v_cba terms
!     -----------
!
!     t^cd_kj*(bd|ai) - t^ba_li*(lj|ck)
!     ---------------------------------
!
      call wf%construct_W_permutation(k, j, i, v_abc, t_abij, g_bdci, g_ljck, one)
!
!     cba -> acb
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     u_acb terms
!     -----------
!
!     t^ad_ik*(cd|bj) - t^cb_lj*(lk|ai)
!     ---------------------------------
!
      call wf%construct_W_permutation(i, k, j, u_abc, t_abij, g_bdcj, g_lkci, one)
!
!     acb -> cab
      call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     v_cab terms
!     -----------
!
!     t^cd_ki*(ad|bj) - t^ab_lj*(li|ck)
!     ---------------------------------
!
      call wf%construct_W_permutation(k, i, j, v_abc, t_abij, g_bdcj, g_lick, one)
!
!     cab -> bca
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     u_bca terms
!     -----------
!
!     t^bd_jk*(cd|ai) - t^ca_li*(lk|bj)
!     ---------------------------------
!
      call wf%construct_W_permutation(j, k, i, u_abc, t_abij, g_bdci, g_lkcj, one)
!
!     bca -> abc
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     v_abc terms
!     -----------
!
!     t^ad_ij*(bd|ck) - t^bc_lk*(lj|ai)
!     ---------------------------------
!
      call wf%construct_W_permutation(i, j, k, v_abc, t_abij, g_bdck, g_ljci, one)
!
   end subroutine construct_W_triples
!
!
   subroutine construct_V_triples(wf, i, j, k, u_abc, v_abc, t2, R2, &
                                  g_bdci, g_bdcj, g_bdck, &
                                  g_bdci_c1, g_bdcj_c1, g_bdck_c1, &
                                  g_ljci, g_lkci, g_lkcj, g_licj, g_lick, g_ljck, &
                                  g_ljci_c1, g_lkci_c1, g_lkcj_c1, g_licj_c1, g_lick_c1, g_ljck_c1)
!!
!!    Construct V
!!    Written by Rolf H. Myhre and Alexander C. Paul, Oct 2020
!!
!!    Calculate the intermediate W for occupied indices i,j,k
!!
!!    Contributions to W
!!    V^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck)    - sum_l t^ab_il (lj|ck)
!!                         +sum_d t^ad_ij(bd|ck)_c1 - sum_l t^ab_il (lj|ck)_c1)
!!
!!    R^abc_ijk = V^abc_ijk/(eps^abc_ijk-omega)
!!
!!    Note: v_abc contains V^abc in the end
!!    Note: This routine does the same as calling construct_W twice:
!!          once with t2 and c1-transformed integrals and
!!          once with R2 and t1-transformed integrals
!!
!!    Each permutation is obtained by two calls to construct_W_permutation
!!
      use reordering, only: sort_123_to_312, sort_123_to_213
!
      implicit none
!
      class(triples) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: u_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out) :: v_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t2
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: R2
!
!     (bd,ci) ordered dbc,i
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdck
!
!     c1-transformed (bd,ci) ordered dbc,i
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdci_c1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdcj_c1
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_bdck_c1
!
!     (lj,ck) ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ljci
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lkci
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lkcj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_licj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lick
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ljck
!
!     c1-transformed (lj,ck) ordered cl,jk
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ljci_c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lkci_c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lkcj_c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_licj_c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_lick_c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_ljck_c1
!
!     u_bac terms
!     -----------
!
!     t^bd_ji*(ad|ck)^c1 - t^ac_lk*(li|bj)^c1 + R^bd_ji*(ad|ck) - R^ac_lk*(li|bj)
!     ---------------------------------------------------------------------------
!
      call wf%construct_W_permutation(j, i, k, u_abc, t2, g_bdck_c1, g_licj_c1, zero)
      call wf%construct_W_permutation(j, i, k, u_abc, R2, g_bdck, g_licj, one)
!
!     bac -> cba
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     v_cba terms
!     -----------
!
!     t^cd_kj*(bd|ai)^c1 - t^ba_li*(lj|ck)^c1 + R^cd_kj*(bd|ai) - R^ba_li*(lj|ck)
!     ---------------------------------------------------------------------------
!
      call wf%construct_W_permutation(k, j, i, v_abc, t2, g_bdci_c1, g_ljck_c1, one)
      call wf%construct_W_permutation(k, j, i, v_abc, R2, g_bdci, g_ljck, one)
!
!     cba -> acb
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     u_acb terms
!     -----------
!
!     t^ad_ik*(cd|bj)^c1 - t^cb_lj*(lk|ai)^c1 + R^ad_ik*(cd|bj) - R^cb_lj*(lk|ai)
!     ---------------------------------------------------------------------------
!
      call wf%construct_W_permutation(i, k, j, u_abc, t2, g_bdcj_c1, g_lkci_c1, one)
      call wf%construct_W_permutation(i, k, j, u_abc, R2, g_bdcj, g_lkci, one)
!
!     acb -> cab
      call sort_123_to_213(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     v_cab terms
!     -----------
!
!     t^cd_ki*(ad|bj)^c1 - t^ab_lj*(li|ck)^c1 + R^cd_ki*(ad|bj) - R^ab_lj*(li|ck)
!     ---------------------------------------------------------------------------
!
      call wf%construct_W_permutation(k, i, j, v_abc, t2, g_bdcj_c1, g_lick_c1, one)
      call wf%construct_W_permutation(k, i, j, v_abc, R2, g_bdcj, g_lick, one)
!
!     cab -> bca
      call sort_123_to_312(v_abc, u_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     u_bca terms
!     -----------
!
!     t^bd_jk*(cd|ai)^c1 - t^ca_li*(lk|bj)^c1 + R^bd_jk*(cd|ai) - R^ca_li*(lk|bj)
!     ---------------------------------------------------------------------------
!
      call wf%construct_W_permutation(j, k, i, u_abc, t2, g_bdci_c1, g_lkcj_c1, one)
      call wf%construct_W_permutation(j, k, i, u_abc, R2, g_bdci, g_lkcj, one)
!
!     bca -> abc
      call sort_123_to_312(u_abc, v_abc, wf%n_v, wf%n_v, wf%n_v)
!
!     v_abc terms
!     -----------
!
!     t^ad_ij*(bd|ck)^c1 - t^bc_lk*(lj|ai)^c1 + R^ad_ij*(bd|ck) - R^bc_lk*(lj|ai)
!     ---------------------------------------------------------------------------
!
      call wf%construct_W_permutation(i, j, k, v_abc, t2, g_bdck_c1, g_ljci_c1, one)
      call wf%construct_W_permutation(i, j, k, v_abc, R2, g_bdck, g_ljci, one)
!
   end subroutine construct_V_triples
!
!
   subroutine construct_W_permutation_triples(wf, o1, o2, o3, W, t2, g_vvv_o, g_vo_oo, alpha)
!!
!!    Construct W permutation
!!    Written by Alexander C. Paul and Rolf H. Myhre, Jan 2021
!!
!!    Calculates one permutation of the triples amplitudes:
!!    These contributions are:
!!    1. Contraction of a doubles amplitude and a v^3o integral over a virtual index
!!    W_bac = t^bd_ji*(ad|ck)
!!
!!    2. Contraction of a doubles amplitude and a o^3v integral over an occupied index
!!    W_bac = - t^ac_lk*(li|bj)
!!
      implicit none
!
      class(triples) :: wf
!
      integer, intent(in) :: o1, o2, o3
      real(dp), intent(in) :: alpha
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(inout) :: W
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in) :: t2
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in) :: g_vvv_o
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: g_vo_oo
!
!     W_v1v2v3 += sum_v4 t^v1v4_o1o2 g_v4_v1v2o3
!
      call dgemm('N', 'N',      &
                 wf%n_v,        &
                 wf%n_v**2,     &
                 wf%n_v,        &
                 one,           &
                 t2(:,:,o1,o2), & ! t_b_d,ji
                 wf%n_v,        &
                 g_vvv_o,       & ! g_d_ac,k
                 wf%n_v,        &
                 alpha,         &
                 W, wf%n_v)
!
!     W_v1v2v3 += sum_o4 g_v1o4_o2o2 t^v2v3_o4o3
!
      call dgemm('N', 'T',     &
                 wf%n_v,       &
                 wf%n_v**2,    &
                 wf%n_o,       &
                 -one,         &
                 g_vo_oo,      & ! g_b_l,ij
                 wf%n_v,       &
                 t2(:,:,:,o3), & ! t_ac_l,k
                 wf%n_v**2,    &
                 one,          &
                 W, wf%n_v)
!
   end subroutine construct_W_permutation_triples
!
!
   subroutine divide_by_orbital_differences_triples(wf, i, j, k, t_abc, omega)
!!
!!    Divide by orbital (energy) differences
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Divide an array of triples amplitudes (single i,j,k)
!!    by the respective orbital energy differences eps^abc_ijk
!!
!!    omega: If present divide by eps^abc_ijk - omega instead of eps^abc_ijk
!!
      implicit none
!
      class(triples) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v,wf%n_v,wf%n_v), intent(inout) :: t_abc
!
      real(dp), optional :: omega
!
      integer a, b, c
!
      real(dp) :: epsilon_ijk, epsilon_c, epsilon_cb
!
      epsilon_ijk = wf%orbital_energies(i) + wf%orbital_energies(j)  &
                  + wf%orbital_energies(k)
!
      if (present(omega)) then
!
         epsilon_ijk = epsilon_ijk + omega
!
      end if
!
!$omp parallel do schedule(static) private(c,b,a,epsilon_c,epsilon_cb)
      do c = 1,wf%n_v
!
         epsilon_c = epsilon_ijk - wf%orbital_energies(wf%n_o + c)
!
         do b = 1,wf%n_v
!
            epsilon_cb = epsilon_c - wf%orbital_energies(wf%n_o + b)
!
            do a = 1,wf%n_v
!
               t_abc(a,b,c) = t_abc(a,b,c)/(epsilon_cb - wf%orbital_energies(wf%n_o + a))
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
!
   end subroutine divide_by_orbital_differences_triples
!
!
   subroutine setup_vvvo_triples(wf, g_vvvo, point, unordered_g_vvvo, batch_s, &
                                 eri_c1, mo, left)
!!
!!    Setup vvvo
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    Setups the vvvo integral and the pointer for a triples calculation
!!    The integral is given in batches over the occupied index
!!
!!    The integral is returned in 2134 order:
!!
!!    g_vvvo:           final sorted integral array
!!    point:            pointer to g_vvvo integral
!!    unordered_g_vvvo: help array for reordering
!!    c_ai:             use the c1-transformed integral
!!    mo:               use the "normal" MO integral (for CCSD(T))
!!    left:             return integral in 1324 order (for the transpose Jacobian)
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_1234_to_1324, sort_1234_to_2134
!
      implicit none
!
      class(triples) :: wf
!
      type(batching_index), intent(in) :: batch_s
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_s%length), target, intent(out)  :: g_vvvo
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_vvvo
      type(eri_adapter), optional, intent(inout)  :: eri_c1
!
      logical, optional, intent(in) :: mo, left
      logical :: mo_, left_
!
      type(eri_tool), allocatable :: eri
      type(eri_adapter), allocatable :: adapter
!
      mo_   = .false.
      left_ = .false.
      if (present(mo)) mo_ = mo
      if (present(left)) left_ = left
!
      if (.not. (mo_ .or. present(eri_c1))) then
!
         call wf%eri_t1%get('vvvo', unordered_g_vvvo, &
                                 first_s=batch_s%first, last_s=batch_s%get_last())
!
      else if(present(eri_c1)) then
!
         call eri_c1%get('vvvo', unordered_g_vvvo, &
                        first_s=batch_s%first, last_s=batch_s%get_last())
!
      else if(mo_) then
!
         eri = eri_tool(wf%L_mo)
         adapter = eri_adapter(eri, wf%n_o, wf%n_v)
!
         call adapter%get('vvvo', unordered_g_vvvo, &
                                 first_s=batch_s%first, last_s=batch_s%get_last()) ! UPS NAME !
!
      end if
!
      if (.not. left_) then
!
         call sort_1234_to_2134(unordered_g_vvvo, g_vvvo, wf%n_v, wf%n_v, &
                                                          wf%n_v, batch_s%length)
!
      else
!
         call sort_1234_to_1324(unordered_g_vvvo, g_vvvo, wf%n_v, wf%n_v, &
                                                          wf%n_v, batch_s%length)
!
      end if
!
      point(1:wf%n_v, 1:wf%n_v, 1:wf%n_v, 1:batch_s%length) => g_vvvo
!
   end subroutine setup_vvvo_triples
!
!
   subroutine setup_vvov_triples(wf, g_vvov, point, unordered_g_vvov, batch_r, left)
!!
!!    Setup vvov
!!    Written by Alexander C. Paul, Oct 2020
!!
!!    Setups the vvov integral and the pointer for a triples calculation
!!    The integral is given in batches over the occupied indices
!!
!!    The integral is returned in 2413 order
!!
!!    g_vvov:           final sorted integral array
!!    point:            pointer to g_vvov integral
!!    unordered_g_vvov: help array for reordering
!!    left:             The integral is returned in 1243 order
!!                      (for the transpose Jacobian)
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_1234_to_2413, sort_1234_to_1243
!
      implicit none
!
      class(triples) :: wf
!
      type(batching_index), intent(in) :: batch_r
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, batch_r%length), target, intent(out)  :: g_vvov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), optional, intent(out) :: unordered_g_vvov
!
      logical, optional, intent(in) :: left
      logical :: left_
!
      left_ = .false.
      if(present(left)) left_ = left
!
      if (.not. left_) then ! Omega and Jacobian
!
         call wf%eri_t1%get('vvov', unordered_g_vvov, &
                                 first_r=batch_r%first, last_r=batch_r%get_last())
!
         call sort_1234_to_2413(unordered_g_vvov, g_vvov, &
                                wf%n_v, wf%n_v, batch_r%length, wf%n_v)
!
      else ! transpose Jacobian
!
         call wf%eri_t1%get('vvov', unordered_g_vvov, &
                                 first_r=batch_r%first, last_r=batch_r%get_last())
!
         call sort_1234_to_1243(unordered_g_vvov, g_vvov, &
                                wf%n_v, wf%n_v, batch_r%length, wf%n_v)
!
      end if
!
      point(1:wf%n_v, 1:wf%n_v, 1:wf%n_v, 1:batch_r%length) => g_vvov
!
   end subroutine setup_vvov_triples
!
!
   subroutine setup_ovov_triples(wf, g_ovov, point, unordered_g_ovov, batch_p, batch_r, mo)
!!
!!    Setup ovov
!!    Written by Alexander C. Paul and Rolf H. Myhre, Dec 2020
!!
!!    Setups the ovov integral and the pointer for a triples calculation
!!    The integral is given in batches over the occupied indices
!!
!!    The integral is returned in 2413 order:
!!
!!    g_ovov:           final sorted integral array
!!    point:            pointer to g_ovov integral
!!    unordered_g_ovov: help array for reordering
!!    mo:               if present and true use the "normal" MO integral
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_1234_to_2413
!
      implicit none
!
      class(triples) :: wf
!
      type(batching_index), intent(in) :: batch_p, batch_r
!
      real(dp), dimension(wf%n_v, wf%n_v, batch_p%length, batch_r%length), &
                                                         target, intent(out) :: g_ovov
!
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_ovov
!
      logical, optional, intent(in) :: mo
      logical :: mo_
!
      type(eri_tool), allocatable :: eri
      type(eri_adapter), allocatable :: adapter
!
      mo_ = .false.
      if (present(mo)) mo_ = mo
!
      if(.not. mo_) then
!
         call wf%eri_t1%get('ovov', unordered_g_ovov, &
                                first_p=batch_p%first, last_p=batch_p%get_last(), &
                                first_r=batch_r%first, last_r=batch_r%get_last())
!
      else
!
!
         eri = eri_tool(wf%L_mo)
         adapter = eri_adapter(eri, wf%n_o, wf%n_v)
!
         call adapter%get('ovov', unordered_g_ovov, &
                                first_p=batch_p%first, last_p=batch_p%get_last(), &
                                first_r=batch_r%first, last_r=batch_r%get_last())
!
      end if
!
      call sort_1234_to_2413(unordered_g_ovov, g_ovov, batch_p%length, wf%n_v, &
                                                       batch_r%length, wf%n_v)
!
      point(1:wf%n_v, 1:wf%n_v, 1:batch_p%length, 1:batch_r%length) => g_ovov
!
   end subroutine setup_ovov_triples
!
!
   subroutine setup_oovo_triples(wf, g_oovo, point, unordered_g_oovo, batch_q, batch_s, &
                                 eri_c1, mo)
!!
!!    Setup oovo
!!    Written by Alexander C. Paul and Rolf H. Myhre, Dec 2020
!!
!!    Setups the oovo integral and the pointer for a triples calculation
!!    The integral is given in batches over the occupied indices
!!
!!    The integral is returned in 3124 order:
!!
!!    g_oovo:           final sorted integral array
!!    point:            pointer to g_oovo integral
!!    unordered_g_oovo: help array for reordering
!!    c_ai:             if present and true use the c1-transformed integral
!!    mo:               if present and true use the "normal" MO integral
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_1234_to_3124
!
      implicit none
!
      class(triples) :: wf
!
      type(batching_index), intent(in) :: batch_q, batch_s
!
      real(dp), dimension(wf%n_v, wf%n_o, batch_q%length, batch_s%length), &
                                                         target, intent(out) :: g_oovo
!
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_oovo
!
      type(eri_adapter), optional, intent(inout)  :: eri_c1
!
      logical, optional, intent(in) :: mo
      logical :: mo_
!
      type(eri_tool), allocatable :: eri
      type(eri_adapter), allocatable :: adapter
!
      mo_ = .false.
      if (present(mo)) mo_ = mo
!
      if(.not. (mo_ .or. present(eri_c1))) then
!
         call wf%eri_t1%get('oovo', unordered_g_oovo, &
                                first_q=batch_q%first, last_q=batch_q%get_last(), &
                                first_s=batch_s%first, last_s=batch_s%get_last())
!
      else if(present(eri_c1)) then
!
         call eri_c1%get('oovo', unordered_g_oovo, &
                            first_q=batch_q%first, last_q=batch_q%get_last(), &
                            first_s=batch_s%first, last_s=batch_s%get_last())
!
      else if(mo_) then
!
         eri = eri_tool(wf%L_mo)
         adapter = eri_adapter(eri, wf%n_o, wf%n_v)
!
         call adapter%get('oovo', unordered_g_oovo, &
                                first_q=batch_q%first, last_q=batch_q%get_last(), &
                                first_s=batch_s%first, last_s=batch_s%get_last())
!
      end if
!
      call sort_1234_to_3124(unordered_g_oovo, g_oovo, wf%n_o, batch_q%length, &
                                                       wf%n_v, batch_s%length)
!
      point(1:wf%n_v, 1:wf%n_o, 1:batch_q%length, 1:batch_s%length) => g_oovo
!
   end subroutine setup_oovo_triples
!
!
   subroutine setup_ooov_triples(wf, g_ooov, point, unordered_g_ooov, batch_p, batch_r)
!!
!!    Setup ooov
!!    Written by Alexander C. Paul and Rolf H. Myhre, Dec 2020
!!
!!    Setups the ooov integral and the pointer for a triples calculation
!!    The integral is given in batches over the occupied indices
!!
!!    The integral is returned in 4213 order
!!
!!    g_oovo:           final sorted integral array
!!    point:            pointer to g_ooov integral
!!    unordered_g_ooov: help array for reordering
!!
      use batching_index_class, only: batching_index
      use reordering, only: sort_1234_to_4213
!
      implicit none
!
      class(triples) :: wf
!
      type(batching_index), intent(in) :: batch_p, batch_r
!
      real(dp), dimension(wf%n_v, wf%n_o, batch_p%length, batch_r%length), &
                                                         target, intent(out) :: g_ooov
!
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      real(dp), dimension(:,:,:,:), intent(out) :: unordered_g_ooov
!
      call wf%eri_t1%get('ooov', unordered_g_ooov, &
                             first_p=batch_p%first, last_p=batch_p%get_last(), &
                             first_r=batch_r%first, last_r=batch_r%get_last())
!
      call sort_1234_to_4213(unordered_g_ooov, g_ooov, batch_p%length, wf%n_o, &
                                                       batch_r%length, wf%n_v)
!
      point(1:wf%n_v, 1:wf%n_o, 1:batch_p%length, 1:batch_r%length) => g_ooov
!
   end subroutine setup_ooov_triples
!
!
   subroutine point_vvvo_triples(wf, point, g_vvvo, length)
!!
!!    Point vvvo
!!    Written by Alexander C. Paul and Rolf H. Myhre, Dec 2020
!!
!!    Sets a pointer to a v3o integral
!!    where the occupied index can be batched over
!!    NB: The occupied index needs to be the last index
!!
!!    point:  pointer to g_vvvo integral
!!    g_vvvo: Integral array in vvvo order
!!    length: length of the batch
!!
      implicit none
!
      class(triples), intent(in) :: wf
!
      integer, intent(in) :: length
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v, length), target, intent(in)  :: g_vvvo
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_v, 1:wf%n_v, 1:wf%n_v, 1:length) => g_vvvo
!
   end subroutine point_vvvo_triples
!
!
   subroutine point_vooo_triples(wf, point, g_ooov, length1, length2)
!!
!!    Point vooo
!!    Written by Alexander C. Paul and Rolf H. Myhre, Dec 2020
!!
!!    Sets a pointer to a o^3v integral where two occupied indices are batched
!!    NB: The integral is sorted v, o, batch1, batch2
!!        and the batches are equally long
!!
!!    point:   pointer to g_ooov integral
!!    g_ooov:  Integral array in vooo order
!!    length1: length of the batch for the second to last index
!!    length2: length of the batch for the last index
!!
      implicit none
!
      class(triples), intent(in) :: wf
!
      integer, intent(in) :: length1, length2
!
      real(dp), dimension(wf%n_v, wf%n_o, length1, length2), target, intent(in)  :: g_ooov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_v, 1:wf%n_o, 1:length1, 1:length2) => g_ooov
!
   end subroutine point_vooo_triples
!
!
   subroutine point_vvoo_triples(wf, point, g_ovov, length1, length2)
!!
!!    Point vvoo
!!    Written by Alexander C. Paul and Rolf H. Myhre, Dec 2020
!!
!!    Sets a pointer to a v^2o^2 integral where both occupied indices are batched
!!    NB: The occupied indices need to be the last ones
!!        and the batches are equally long
!!
!!    point:   pointer to g_ovov integral
!!    g_ovov:  Integral array in vvoo order
!!    length1: length of the batch for the second to last index
!!    length2: length of the batch for the last index
!!
      implicit none
!
      class(triples), intent(in) :: wf
!
      integer, intent(in) :: length1, length2
!
      real(dp), dimension(wf%n_v, wf%n_v, length1, length2), target, intent(in)  :: g_ovov
      real(dp), dimension(:,:,:,:), pointer, contiguous, intent(out) :: point
!
      point(1:wf%n_v, 1:wf%n_v, 1:length1, 1:length2) => g_ovov
!
   end subroutine point_vvoo_triples
!
!
   subroutine construct_cholesky_c1_triples(wf, L, L_c1, c1)
!!
!!    Construct Cholesky C1
!!
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Constructs the "C1 Cholesky vector"
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_qp_c1 if .true.
!!
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(triples), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:,:,:), allocatable :: L_Jia
!
      call wf%construct_cholesky_c1_oo(L, L_c1, c1)
      call wf%construct_cholesky_c1_vo(L, L_c1, c1)
      call wf%construct_cholesky_c1_vv(L, L_c1, c1)
!
      call mem%alloc(L_Jia, L%n_J, wf%n_o, wf%n_v)
      call zero_array(L_Jia, L%n_J*wf%n_o*wf%n_v)
      call L_c1%set(L_Jia, 1, wf%n_o, wf%n_o + 1, wf%n_mo)
      call mem%dealloc(L_Jia, L%n_J, wf%n_o, wf%n_v)
!
   end subroutine construct_cholesky_c1_triples
!
!
   subroutine construct_cholesky_c1_oo_triples(wf, L, L_c1, c1)
!!
!!    Construct Cholesky oo C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_J_ij_C1 = sum_b L_J_ib_T1 c_bj
!!
!!    and returns the result in L_J_ij_C1
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_ji_C1 if .true.
!!
!
      use batching_index_class, only: batching_index
      implicit none
!
      class(triples), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:), allocatable :: L_J_oo
      real(dp), dimension(:), allocatable :: L_J_ov
!
      type(batching_index) :: batch_o
!
      integer :: o_batch
!
      batch_o = batching_index(wf%n_o)
      call mem%batch_setup(batch_o, 0, L%n_J*(wf%n_v + wf%n_o), tag='Cholesky c1 oo')
!
      call mem%alloc(L_J_oo, L%n_J*batch_o%max_length*wf%n_o)
      call mem%alloc(L_J_ov, L%n_J*batch_o%max_length*wf%n_v)
!
      do o_batch = 1,batch_o%num_batches
!
         call batch_o%determine_limits(o_batch)
!
         call L%get(L_J_ov, batch_o%first, batch_o%get_last(), wf%n_o + 1, wf%n_mo)
!
         call dgemm('N', 'N',                               &
                    L%n_J*batch_o%length, wf%n_o, wf%n_v,   &
                    one,                                    &
                    L_J_ov, L%n_J*batch_o%length,           &
                    c1, wf%n_v,                             &
                    zero,                                   &
                    L_J_oo, L%n_J*batch_o%length)
!
         call L_c1%set(L_J_oo, batch_o%first, batch_o%get_last(), 1, wf%n_o)
!
      enddo
!
      call mem%dealloc(L_J_oo, L%n_J*batch_o%max_length*wf%n_o)
      call mem%dealloc(L_J_ov, L%n_J*batch_o%max_length*wf%n_v)
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_c1_oo_triples
!
!
   subroutine construct_cholesky_c1_vo_triples(wf, L, L_c1, c1)
!!
!!    Construct Cholesky vo C1
!!    Written by Rolf. H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_J_ai_C1 = sum_b L_J_ab_T1 c_bi - sum_j c_aj L_J_ji_T1
!!
!!    and returns the result in L_J_ai_c1
!!
!!    qp: optional (default = .false) reordering logical, returns L_J_ia_C1 if .true.
!!
      use batching_index_class, only: batching_index
      use array_utilities, only: zero_array
      use reordering, only: add_132_to_123, sort_123_to_132
!
      implicit none
!
      class(triples), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ij
      real(dp), dimension(:,:,:), allocatable :: L_J_ji
      real(dp), dimension(:,:,:), allocatable :: L_J_ai
      real(dp), dimension(:,:,:), allocatable :: L_J_ia
      real(dp), dimension(:,:,:), allocatable :: L_J_ab
!
      type(batching_index) :: batch_a, batch_j
!
      integer :: a_batch, j_batch, req
!
      batch_j = batching_index(wf%n_o)
      batch_a = batching_index(wf%n_v)
      req = L%n_J*(wf%n_o + max(wf%n_v,wf%n_o))
!
      call mem%batch_setup(batch_j, batch_a, 0, req, req, 0, tag='Cholesky c1 vo')
!
      do a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(a_batch)
!
         call mem%alloc(L_J_ab, L%n_J, batch_a%length, wf%n_v)
!
         call L%get(L_J_ab,                        &
                    wf%n_o + batch_a%first,        &
                    wf%n_o + batch_a%get_last(),   &
                    wf%n_o + 1,                    &
                    wf%n_mo)
!
         call mem%alloc(L_J_ai, L%n_J, batch_a%length, wf%n_o)
!
         call dgemm('N', 'N',             &
                    L%n_J*batch_a%length, &
                    wf%n_o,               &
                    wf%n_v,               &
                    one,                  &
                    L_J_ab,               &
                    L%n_J*batch_a%length, &
                    c1,                   &
                    wf%n_v,               &
                    zero,                 &
                    L_J_ai,               &
                    L%n_J*batch_a%length)
!
         call mem%dealloc(L_J_ab, L%n_J, batch_a%length, wf%n_v)
         call mem%alloc(L_J_ia, L%n_J, wf%n_o, batch_a%length)
!
         do j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(j_batch)
!
            call mem%alloc(L_J_ji, L%n_J, batch_j%length, wf%n_o)
            call L%get(L_J_ji, batch_j%first, batch_j%get_last(), 1, wf%n_o)

            call mem%alloc(L_J_ij, L%n_J, wf%n_o, batch_j%length)
            call sort_123_to_132(L_J_ji, L_J_ij, L%n_J, batch_j%length, wf%n_o)
            call mem%dealloc(L_J_ji, L%n_J, batch_j%length, wf%n_o)
!
            call dgemm('N', 'T',       &
                       L%n_J*wf%n_o,   &
                       batch_a%length, &
                       batch_j%length, &
                       -one,           &
                       L_J_ij,         &
                       L%n_J*wf%n_o,   &
                       c1(batch_a%first, batch_j%first), &
                       wf%n_v,         &
                       zero,           &
                       L_J_ia,         &
                       L%n_J*wf%n_o)
!
            call mem%dealloc(L_J_ij, L%n_J, wf%n_o, batch_j%length)
            call add_132_to_123(one, L_J_ia, L_J_ai, L%n_J, batch_a%length, wf%n_o)
!
         enddo
!
         call L_c1%set(L_J_ai, wf%n_o+batch_a%first, wf%n_o+batch_a%get_last(), 1, wf%n_o)
         call mem%dealloc(L_J_ai, L%n_J, batch_a%length, wf%n_o)
         call mem%dealloc(L_J_ia, L%n_J, wf%n_o, batch_a%length)
!
      enddo
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_c1_vo_triples
!
!
   subroutine construct_cholesky_c1_vv_triples(wf, L, L_c1, c1)
!!
!!    Construct Cholesky vv C1
!!    Written by Rolf H. Myhre, Jun 2020
!!
!!    based on routines by Alexander C. Paul
!!
!!    Computes
!!
!!       L_ab_J_c1= - sum_i c_ai L_ib_J_T1 ,
!!
!!    and returns the result in L_J_ab_c1
!!
!
      use batching_index_class, only: batching_index
      use reordering, only: sort_123_to_132
!
      implicit none

      class(triples), intent(inout) :: wf
!
      class(abstract_eri_cholesky), intent(inout) :: L
      class(abstract_eri_cholesky), intent(inout) :: L_c1
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: c1
!
      real(dp), dimension(:,:,:), allocatable :: L_J_ab, L_J_ba, L_J_jb, L_J_bj
!
      type(batching_index) :: batch_b
!
      integer :: b_batch
!
      batch_b = batching_index(wf%n_v)
      call mem%batch_setup(batch_b, 0, 2*L%n_J*max(wf%n_v, wf%n_o), tag='Cholesky c1 vv')
!
      do b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(b_batch)
!
         call mem%alloc(L_J_jb, L%n_J, wf%n_o, batch_b%length)
         call L%get(L_J_jb, 1, wf%n_o, wf%n_o + batch_b%first, wf%n_o + batch_b%get_last())
!
         call mem%alloc(L_J_bj, L%n_J, batch_b%length, wf%n_o)
         call sort_123_to_132(L_J_jb, L_J_bj, L%n_J, wf%n_o, batch_b%length)
         call mem%dealloc(L_J_jb, L%n_J, wf%n_o, batch_b%length)
!
         call mem%alloc(L_J_ba, L%n_J, batch_b%length, wf%n_v)
         call dgemm('N', 'T',             &
                    L%n_J*batch_b%length, &
                    wf%n_v,               &
                    wf%n_o,               &
                    -one,                 &
                    L_J_bj,               &
                    L%n_J*batch_b%length, &
                    c1,                   &
                    wf%n_v,               &
                    zero,                 &
                    L_J_ba,               &
                    L%n_J*batch_b%length)
!
         call mem%dealloc(L_J_bj, L%n_J, batch_b%length, wf%n_o)
!
         call mem%alloc(L_J_ab, L%n_J, wf%n_v, batch_b%length)
         call sort_123_to_132(L_J_ba, L_J_ab, L%n_J, batch_b%length, wf%n_v)
         call mem%dealloc(L_J_ba, L%n_J, batch_b%length, wf%n_v)
!
         call L_c1%set(L_J_ab, wf%n_o+1, wf%n_mo, wf%n_o + batch_b%first, wf%n_o + batch_b%get_last())
         call mem%dealloc(L_J_ab, L%n_J, wf%n_v, batch_b%length)
!
      enddo
!
      call mem%batch_finalize()
!
   end subroutine construct_cholesky_c1_vv_triples
!
!
end module triples_class
