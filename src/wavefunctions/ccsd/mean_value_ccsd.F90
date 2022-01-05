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
submodule (ccsd_class) mean_value_ccsd
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
   module subroutine calculate_energy_ccsd(wf)
!!
!!    Calculate energy (CCSD)
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Andreas Skeidsvoll, 2018
!!
!!    Calculates the CCSD energy. This is only equal to the actual
!!    energy when the ground state equations are solved, of course.
!!
!!       E = E_hf + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb ! g_iajb
!
      real(dp) :: omp_correlation_energy
!
      integer :: a, i, b, j, ai, bj, aibj
!
      type(timings) :: timer
!
      timer = timings('Calculate energy', pl='n')
      call timer%turn_on()
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('ovov', g_iajb)
!
      omp_correlation_energy = zero
!
!$omp parallel do private(a,i,ai,bj,j,b,aibj) reduction(+:omp_correlation_energy)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = (i-1)*wf%n_v + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
!
                  aibj = (max(ai,bj)*(max(ai,bj)-3)/2) + ai + bj
!
                  omp_correlation_energy = omp_correlation_energy                   &
                                         +(wf%t2(aibj) + (wf%t1(a,i))*(wf%t1(b,j)))* &
                                          (two*g_iajb(i,a,j,b) - g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      wf%correlation_energy = omp_correlation_energy
!
      wf%energy = wf%hf_energy + wf%correlation_energy
!
      call timer%turn_off()
!
   end subroutine calculate_energy_ccsd
!
!
end submodule mean_value_ccsd
