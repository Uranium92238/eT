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
   implicit none
!
   type, abstract, extends(ccsd) :: triples
!
   contains
!
      procedure :: construct_W            => construct_W_triples
      procedure :: divide_by_orbital_differences => divide_by_orbital_differences_triples
!
   end type triples
!
!
contains
!
!
   subroutine construct_W_triples(wf, i, j, k, t_abc, u_abc, t_abij, &
                             g_bdci, g_bdcj, g_bdck, &
                             g_ljci, g_lkci, g_lkcj, g_licj, g_lick, g_ljck, &
                             overwrite)
!!
!!    Construct W intermediate
!!    Written by Rolf H. Myhre, January 2019
!!
!!    Calculate the intermediate W for occupied indices i,j,k
!!
!!    Contributions to W
!!    W^abc_ijk = P^abc_ijk(sum_d t^ad_ij(bd|ck) - sum_l t^ab_il (lj|ck))
!!
!!    T^abc_ijk = W^abc_ijk/eps^abc_ijk
!!    R^abc_ijk = (W^abc_ijk + W^abc_ijk(C1))/(eps^abc_ijk-omega)
!!
      use reordering
!
      implicit none
!
      class(triples) :: wf
!
      integer, intent(in) :: i, j, k
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: t_abc
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(out)          :: u_abc
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_o, wf%n_o), intent(in)   :: t_abij
!
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdci
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdcj
      real(dp), dimension(wf%n_v, wf%n_v, wf%n_v), intent(in)           :: g_bdck
!
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_ljci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lkci
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lkcj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_licj
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_lick
      real(dp), dimension(wf%n_o, wf%n_v), intent(in)                   :: g_ljck
!
!     If present and true, t_abc is not overwritten by first dgemm
!     Needed to calculate triples part of the right excitation vector
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
!     u_abc terms
!     -----------
!
!     t^ad_ij*(bd|ck)
!     ---------------
!
      call dgemm('N', 'N',        &
                 wf%n_v,          &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 one,             &
                 t_abij(:,:,i,j), & !t_a_d,ij
                 wf%n_v,          &
                 g_bdck,          & !g_d_bc,k
                 wf%n_v,          &
                 alpha,           &
                 t_abc,           &
                 wf%n_v)
!
!     t^ab_lj*(li|ck)
!     ---------------
!
      call dgemm('N', 'N',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_o,          &
                 -one,            &
                 t_abij(:,:,:,j), & !t_ab_l,j
                 wf%n_v**2,       &
                 g_lick,          & !g_l_c,ik
                 wf%n_o,          &
                 one,             &
                 t_abc,           &
                 wf%n_v**2)
!
!
!     t^cd_ki*(ad|bj)
!     ---------------
!
      call dgemm('T', 'T',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_v,          &
                 one,             &
                 g_bdcj,          & !g_d_ab,j
                 wf%n_v,          &
                 t_abij(:,:,k,i), & !t_c_d,ki
                 wf%n_v,          &
                 one,             &
                 t_abc,           &
                 wf%n_v**2)
!
!     t^bc_lk*(lj|ai)
!     ---------------
!
      call dgemm('T', 'T',        &
                 wf%n_v,          &
                 wf%n_v**2,       &
                 wf%n_o,          &
                 -one,            &
                 g_ljci,          & !g_l_a,ji
                 wf%n_o,          &
                 t_abij(:,:,:,k), & !t_bc_l,k
                 wf%n_v**2,       &
                 one,             &
                 t_abc,           &
                 wf%n_v)
!
!     u_acb terms
!     -----------
!
!     t^ad_ik*(cd|bj)
!     ---------------
!
      call dgemm('N', 'N',        &
                 wf%n_v,          &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 one,             &
                 t_abij(:,:,i,k), & !t_a_d,ik
                 wf%n_v,          &
                 g_bdcj,          & !g_d_cb,j
                 wf%n_v,          &
                 zero,            &
                 u_abc,           &
                 wf%n_v)
!
!     t^ac_lk*(li|bj)
!     ---------------
!
      call dgemm('N', 'N',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_o,          &
                 -one,            &
                 t_abij(:,:,:,k), & !t_ac_l,k
                 wf%n_v**2,       &
                 g_licj,          & !g_l_b,ij
                 wf%n_o,          &
                 one,             &
                 u_abc,           &
                 wf%n_v**2)
!
!
!     t^bd_ji*(ad|ck)
!     ---------------
!
      call dgemm('T', 'T',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_v,          &
                 one,             &
                 g_bdck,          & !g_d_ac,k
                 wf%n_v,          &
                 t_abij(:,:,j,i), & !t_b_d,ji
                 wf%n_v,          &
                 one,             &
                 u_abc,           &
                 wf%n_v**2)
!
!     t^cb_lj*(lk|ai)
!     ---------------
!
      call dgemm('T', 'T',        &
                 wf%n_v,          &
                 wf%n_v**2,       &
                 wf%n_o,          &
                 -one,            &
                 g_lkci,          & !g_l_a,ki
                 wf%n_o,          &
                 t_abij(:,:,:,j), & !t_cb_l,j
                 wf%n_v**2,       &
                 one,             &
                 u_abc,           &
                 wf%n_v)
!
!
      call sort_123_to_132_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
!     u_bac terms
!     -----------
!
!     t^cd_kj*(bd|ai)
!     ---------------
!
      call dgemm('T', 'T',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_v,          &
                 one,             &
                 g_bdci,          & !g_d_ba,i
                 wf%n_v,          &
                 t_abij(:,:,k,j), & !t_c_d,kj
                 wf%n_v,          &
                 zero,            &
                 u_abc,           &
                 wf%n_v**2)
!
!     t^ba_li*(lj|ck)
!     ---------------
!
      call dgemm('N', 'N',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_o,          &
                 -one,            &
                 t_abij(:,:,:,i), & !t_ba_l,i
                 wf%n_v**2,       &
                 g_ljck,          & !g_l_c,jk
                 wf%n_o,          &
                 one,             &
                 u_abc,           &
                 wf%n_v**2)
!
      call sort_123_to_213_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
!
!     u_cab terms
!     -----------
!
!     t^bd_jk*(cd|ai)
!     ---------------
!
      call dgemm('T', 'T',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_v,          &
                 one,             &
                 g_bdci,          & !g_d_ca,i
                 wf%n_v,          &
                 t_abij(:,:,j,k), & !t_b_d,jk
                 wf%n_v,          &
                 zero,            &
                 u_abc,           &
                 wf%n_v**2)
!
!     t^ca_li*(lk|bj)
!     ---------------
!
      call dgemm('N', 'N',        &
                 wf%n_v**2,       &
                 wf%n_v,          &
                 wf%n_o,          &
                 -one,            &
                 t_abij(:,:,:,i), & !t_ca_l,i
                 wf%n_v**2,       &
                 g_lkcj,          & !g_l_c,kj
                 wf%n_o,          &
                 one,             &
                 u_abc,           &
                 wf%n_v**2)
!
      call sort_123_to_231_and_add(u_abc,t_abc,wf%n_v,wf%n_v,wf%n_v)
!
   end subroutine construct_W_triples
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
!!    If present divide by eps^abc_ijk - omega instead of eps^abc_ijk
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
!$omp parallel do schedule(static) private(a)
      do a = 1,wf%n_v
!
         t_abc(a,a,a) = zero
!
      enddo
!$omp end parallel do
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
end module triples_class
