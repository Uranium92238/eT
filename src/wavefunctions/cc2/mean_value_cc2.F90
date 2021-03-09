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
submodule (cc2_class) mean_value_cc2
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
   module subroutine prepare_for_density_cc2(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Sarai D. Folekstad, May 2019
!!
      implicit none
!
      class(cc2), intent(inout) :: wf
!
      call wf%initialize_t2()
      call wf%initialize_t2bar()
!
      call wf%construct_t2()
      call wf%construct_t2bar()
!
   end subroutine prepare_for_density_cc2
!
!
   module subroutine calculate_energy_cc2(wf)
!!
!!    Calculate energy 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Jan 2019
!!
!!    E = E_HF + sum_aibj (t_i^a*t_j^b + t_ij^ab) L_iajb
!!
      class(cc2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj, g_iajb 
!
      real(dp) :: omp_correlation_energy
!
      integer :: a, i, b, j
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%eri%get_eri_t1('vovo', g_aibj)
      call wf%eri%get_eri_t1('ovov', g_iajb)
!
      omp_correlation_energy = zero 
!
!$omp parallel do private(a,i,b,j) reduction(+:omp_correlation_energy)
      do b = 1, wf%n_v
         do i = 1, wf%n_o 
            do j = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  omp_correlation_energy = omp_correlation_energy +                            &
                                           (wf%t1(a, i)*wf%t1(b, j) -                          &
                                           (g_aibj(a,i,b,j))/(wf%orbital_energies(wf%n_o + a)  &
                                                             + wf%orbital_energies(wf%n_o + b) &
                                                             - wf%orbital_energies(i)          &
                                                             - wf%orbital_energies(j)))        &
                                           *(two*g_iajb(i,a,j,b)-g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do 
!
      call mem%dealloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      wf%correlation_energy = omp_correlation_energy
!
      wf%energy = wf%hf_energy + wf%correlation_energy
!
   end subroutine calculate_energy_cc2
!
!
end submodule mean_value_cc2
