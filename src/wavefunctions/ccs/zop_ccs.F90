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
submodule (ccs_class) zop_ccs
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
   module subroutine prepare_for_density_ccs(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
!     For now, do nothing.
!
      call output%printf('- No preparations for the density for ' // trim(wf%name_) // &
                         ' wavefunction.', pl='v', fs='(/t3,a)')
!
   end subroutine prepare_for_density_ccs
!
!
   module subroutine construct_gs_density_ccs(wf)
!!
!!    Construct one-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    Constructs the one-electron density 
!!    matrix in the T1 basis
!!
!!    D_pq = < Λ | E_pq | CC >
!!
      implicit none
!
      class(ccs) :: wf
!
      call zero_array(wf%density, (wf%n_mo)**2)
!
      call wf%gs_one_el_density_ccs_oo(wf%density)
      call wf%gs_one_el_density_ccs_vo(wf%density, wf%t1bar)
!
   end subroutine construct_gs_density_ccs
!
!
   module subroutine gs_one_el_density_ccs_oo_ccs(wf, density)
!!
!!    One electron density oo
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ii = 2  
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, wf%n_o
!
         density(i,i) = density(i,i) + two  
!
      enddo
!$omp end parallel do
!
   end subroutine gs_one_el_density_ccs_oo_ccs
!
!
   module subroutine gs_one_el_density_ccs_vo_ccs(wf, density, tbar_ai)
!!
!!    One electron density vo
!!    Written by Sarai D. Folkestad, 2019
!!
!!    D_ai = tbar_ai 
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      real(dp), dimension(wf%n_v, wf%n_o) :: tbar_ai
!
      integer :: i, a
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!        
            density(wf%n_o + a, i) = density(wf%n_o + a, i) + tbar_ai(a, i)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine gs_one_el_density_ccs_vo_ccs
!
!
   module function calculate_expectation_value_ccs(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Calculate the expectation value of a one-electron
!!    operator Â
!!
!!       < A > = < Λ | Â | CC > = sum_pq A_pq D_pq
!!
!!    where A_pq are the T1-transformed integrals
!!    and D_pq is the a one-electron density matrix
!!    in the T1-basis
!!
      implicit none
!  
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
!
      real(dp) :: expectation_value
!
      real(dp) :: ddot
!
      expectation_value = ddot(wf%n_mo**2, A, 1, density, 1)
!
   end function calculate_expectation_value_ccs
!
!
   module subroutine calculate_energy_ccs(wf)
!!
!!    Calculate energy
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the CCSD energy. This is only equal to the actual
!!    energy when the ground state equations are solved, of course.
!!
!!       E = E_hf + sum_aibj t_i^a t_j^b L_iajb
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      real(dp) :: correlation_energy 
!
      integer :: a, i, b, j, ai, bj, aibj
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov(g_iajb)
!
      correlation_energy = zero 
!
!$omp parallel do private(a,i,ai,bj,j,b,aibj) reduction(+:correlation_energy)
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
                  correlation_energy = correlation_energy + (wf%t1(a,i))*(wf%t1(b,j))* &
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
      wf%energy = wf%hf_energy + correlation_energy
!
   end subroutine calculate_energy_ccs
!
!
   module subroutine calculate_energy_omega_term_ccs(wf)
!!
!!    Calculate energy omega term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds multipliers dot omega to the energy,
!!
!!       energy += Σ_μ tbar_μ Ω_μ,
!!
!!    which appears in the variational energy expression <Λ|H|CC> when Omega ≠ 0.
!!    This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(:), allocatable :: multipliers, omega
!
      real(dp), external :: ddot
!
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
      call mem%alloc(omega, wf%n_gs_amplitudes)
!
      call wf%get_multipliers(multipliers)
      call wf%construct_omega(omega)
!
      wf%energy = wf%energy + ddot(wf%n_gs_amplitudes, multipliers, 1, omega, 1)
!
      call mem%dealloc(multipliers, wf%n_gs_amplitudes)
      call mem%dealloc(omega, wf%n_gs_amplitudes)
!
   end subroutine calculate_energy_omega_term_ccs
!
!
   module subroutine calculate_energy_length_dipole_term_ccs(wf, electric_field)
!!
!!    Calculate energy length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds dipole part of the length gauge electromagnetic potential to the energy,
!!
!!       energy += 2 sum_ii (-μ·E)_ii,
!!
!!    where μ is the vector of electric dipole integral matrices and E is a uniform classical electric
!!    vector field. This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(3), intent(in) :: electric_field
!
      real(dp), dimension(:,:,:), allocatable :: mu
!
      integer :: i
!
!     Construct t1 transformed dipole moment
!
      call mem%alloc(mu, wf%n_mo, wf%n_mo, 3)
      call wf%construct_mu(mu)
!
!     Add one-electron electric field contribution to the diagonal of Fock and one-electron integral terms
!
      do i = 1, wf%n_o
!
         wf%energy = wf%energy - two*(mu(i, i, 1)*electric_field(1)   &
                                      + mu(i, i, 2)*electric_field(2) &
                                      + mu(i, i, 3)*electric_field(3))
!
      enddo
!
   end subroutine calculate_energy_length_dipole_term_ccs
!
!
end submodule zop_ccs