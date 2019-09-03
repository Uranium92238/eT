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
   module subroutine prepare_for_density_cc3(wf)
!!
!!    Prepare the construction of the CC3 contribution to the GS density
!!    Written by Alexander Paul and Rolf H. Myhre, March 2019
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
!!    tbar_μ3 = (- ε^abc_ijk)^-1 (tbar_μ1 < μ1 | [H,τ_ν3] | R > + tbar_μ2 < μ2 | [H,τ_ν3] | R >
!!
!!    ρ^L_kl += -1/2 sum_{abc}{ij} tbar^abc_ijl t^abc_ijk
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
   end subroutine construct_y_intermediate_vo3_cc3
!