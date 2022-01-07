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
submodule (lowmem_cc2_class) mean_value_lowmem_cc2
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
   module subroutine calculate_energy_lowmem_cc2(wf)
!!
!!     Calculate energy (lowmem CC2)
!!     Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!     Andreas Skeidsvoll, 2018
!!
!!     Calculates the lowmem CC2 energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(lowmem_cc2), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      real(dp) :: correlation_energy
!
      integer :: i, j, a, b
!
      integer :: req0, req1_i, req1_j, req2
!
      integer :: current_i_batch, current_j_batch
!
      type(batching_index) :: batch_i, batch_j
!
      type(timings) :: timer
!
      timer = timings('Calculate energy', pl='n')
      call timer%turn_on()
!
      req0 = 0
!
      req1_i = (wf%n_v)*(wf%eri_t1%n_J)
      req1_j = (wf%n_v)*(wf%eri_t1%n_J)
!
      req2 = 2*(wf%n_v**2)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2, &
                           tag='calculate_energy_lowmem_cc2')
!
      correlation_energy = zero
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(g_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
            call mem%alloc(g_iajb, batch_i%length, wf%n_v, batch_j%length, wf%n_v)
!
            call wf%eri_t1%get('vovo', g_aibj, &
                                   1, wf%n_v, &
                                   batch_i%first, batch_i%get_last(), &
                                   1, wf%n_v, &
                                   batch_j%first, batch_j%get_last())
!
            call wf%eri_t1%get('ovov', g_iajb, &
                                   batch_i%first, batch_i%get_last(), &
                                   1, wf%n_v, &
                                   batch_j%first, batch_j%get_last(), &
                                   1, wf%n_v)
!
!$omp parallel do private(b,i,j,a) reduction(+:correlation_energy)
            do b = 1, wf%n_v
               do i = 1, batch_i%length
                  do j = 1, batch_j%length
                     do a = 1, wf%n_v
!
                        correlation_energy = correlation_energy + &
                                             (wf%t1(a, i + batch_i%first - 1)*wf%t1(b, j + batch_j%first - 1) &
                                             - (g_aibj(a, i, b, j))/(wf%orbital_energies(wf%n_o + a) &
                                                               + wf%orbital_energies(wf%n_o + b) &
                                                               - wf%orbital_energies(i + batch_i%first - 1) &
                                                               - wf%orbital_energies(j + batch_j%first - 1)))&
                                             *(two*g_iajb(i, a, j, b)-g_iajb(i, b, j, a))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_aibj, wf%n_v, batch_i%length, wf%n_v, batch_j%length)
            call mem%dealloc(g_iajb, batch_i%length, wf%n_v, batch_j%length, wf%n_v)
!
         enddo
      enddo
!
      call mem%batch_finalize()
!
      wf%correlation_energy = correlation_energy
!
      wf%energy = wf%hf_energy + correlation_energy
!
      call timer%turn_off()
!
   end subroutine calculate_energy_lowmem_cc2
!
!
end submodule mean_value_lowmem_cc2
