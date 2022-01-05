!
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
submodule (mp2_class) mean_value_mp2
!
!!
!!    Mean-value submodule
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
!!    Calculates the MP2 energy as:
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
      use array_utilities, only: copy_and_scale
      use reordering, only: add_1432_to_1234
      use timings_class, only: timings
!
      implicit none
!
      class(mp2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
!
      real(dp) :: e2_neg, eps_bj, eps_ibj
      integer  :: a, i, b, j
!
      type(timings) :: timer
!
      timer = timings('Calculate energy', pl='n')
      call timer%turn_on()
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%eri%get_eri_mo('vovo', g_aibj, 1, wf%n_v, 1, wf%n_o, 1, wf%n_v, 1, wf%n_o)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call copy_and_scale(two, g_aibj, L_aibj, (wf%n_v*wf%n_o)**2)
      call add_1432_to_1234(-one, g_aibj, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      e2_neg = zero
!
!$omp parallel do schedule(static) reduction(+:e2_neg) &
!$omp private(a, i, b, j, eps_bj, eps_ibj)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            eps_bj = wf%orbital_energies(wf%n_o + b) - wf%orbital_energies(j)
!
            do i = 1, wf%n_o
!
               eps_ibj = eps_bj - wf%orbital_energies(i)
!
               do a = 1, wf%n_v
!
                  e2_neg = e2_neg + g_aibj(a, i, b, j)*L_aibj(a, i, b, j) &
                                  /(wf%orbital_energies(wf%n_o + a)+eps_ibj)
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%energy = wf%hf_energy - e2_neg
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call timer%turn_off()
!
   end subroutine calculate_energy_mp2
!
!
end submodule mean_value_mp2
