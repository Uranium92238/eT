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
submodule (mp2_class) zop_mp2
!
!!
!!    Zeroth order properties submodule 
!!
!!    Contains routines related to the mean values, i.e. 
!!    the construction of density matrices as well as expectation 
!!    value calculation.
!!
!
   implicit none 
!
!
contains
!
!
   module subroutine calculate_energy_mp2(wf)
!!
!!    Calculate energy
!!    Written by Andreas Skeidsvoll, 2018
!!
!!    Calculates the MP2 energy from HF energy, E_HF, vovo integrals, g_aibj, 
!!    and the orbital energies, eps. The total MP2 energy is calculated as
!!
!!       E = E_HF - sum_aibj g_aibj*L_aibj/(eps(a)+eps(b)-eps(i)-eps(j))
!!
!!    where
!!
!!       L_aibj = 2*g_aibj - g_ajbi.
!!
!!    On entry, it is assumed that the energy is equal to the HF energy 
!!    (i.e. the routine only adds the correction to the energy variable).
!!
      implicit none
!
      class(mp2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
      real(dp), dimension(:), allocatable :: eps
!
      real(dp) :: e2_neg
      integer  :: a, i, b, j
!
      call mem%alloc(eps, wf%n_mo)
      eps = wf%orbital_energies
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%get_vovo(g_aibj)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      L_aibj = two*g_aibj
      call add_1432_to_1234(-one, g_aibj, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      e2_neg = zero 
!
!$omp parallel do schedule(static) private(a, i, b, j) reduction(+:e2_neg)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  e2_neg = e2_neg + g_aibj(a, i, b, j)*L_aibj(a, i, b, j)/(eps(wf%n_o+a)+eps(wf%n_o+b)-eps(i)-eps(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%energy = wf%hf_energy - e2_neg
!
      call mem%dealloc(eps, wf%n_mo)
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine calculate_energy_mp2
!
!
end submodule zop_mp2