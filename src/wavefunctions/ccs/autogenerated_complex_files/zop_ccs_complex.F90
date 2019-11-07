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
submodule (ccs_class) zop_ccs_complex
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
   module subroutine prepare_for_density_ccs_complex(wf)
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
   end subroutine prepare_for_density_ccs_complex
!
!
   module subroutine construct_gs_density_ccs_complex(wf)
!!
!!    Construct one_complex-electron density
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    Constructs the one_complex-electron density 
!!    matrix in the T1 basis
!!
!!    D_pq = < Λ | E_pq | CC >
!!
      implicit none
!
      class(ccs) :: wf
!
      call zero_array_complex(wf%density_complex, (wf%n_mo)**2)
!
      call wf%gs_one_el_density_ccs_oo_complex(wf%density_complex)
      call wf%gs_one_el_density_ccs_vo_complex(wf%density_complex, wf%t1bar_complex)
!
   end subroutine construct_gs_density_ccs_complex
!
!
   module subroutine gs_one_el_density_ccs_oo_ccs_complex(wf, density)
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
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
      integer :: i
!
!$omp parallel do private(i)
      do i = 1, wf%n_o
!
         density(i,i) = density(i,i) + two_complex  
!
      enddo
!$omp end parallel do
!
   end subroutine gs_one_el_density_ccs_oo_ccs_complex
!
!
   module subroutine gs_one_el_density_ccs_vo_ccs_complex(wf, density, tbar_ai)
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
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      complex(dp), dimension(wf%n_v, wf%n_o) :: tbar_ai
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
   end subroutine gs_one_el_density_ccs_vo_ccs_complex
!
!
   module function calculate_expectation_value_ccs_complex(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Calculate the expectation value of a one_complex-electron
!!    operator Â
!!
!!       < A > = < Λ | Â | CC > = sum_pq A_pq D_pq
!!
!!    where A_pq are the T1-transformed integrals
!!    and D_pq is the a one_complex-electron density matrix
!!    in the T1-basis
!!
      implicit none
!  
      class(ccs), intent(in) :: wf
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
!
      complex(dp) :: expectation_value
!
!
      expectation_value = our_zdotu(wf%n_mo**2, A, 1, density, 1)
!
   end function calculate_expectation_value_ccs_complex
!
!
   module subroutine calculate_energy_ccs_complex(wf)
!!
!!    Calculate energy_complex
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the CCSD energy_complex. This is only equal to the actual
!!    energy_complex when the ground state equations are solved, of course.
!!
!!       E = E_hf + sum_aibj t_i^a t_j^b L_iajb
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_iajb
!
      complex(dp) :: correlation_energy 
!
      integer :: a, i, b, j, ai, bj, aibj
!
      call mem%alloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      call wf%get_ovov_complex(g_iajb)
!
      correlation_energy = zero_complex 
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
                  correlation_energy = correlation_energy + (wf%t1_complex(a,i))*(wf%t1_complex(b,j))* &
                                                      (two_complex*g_iajb(i,a,j,b) - g_iajb(i,b,j,a))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajb, wf%n_o, wf%n_v, wf%n_o, wf%n_v)
!
      wf%energy_complex = wf%hf_energy_complex + correlation_energy
!
   end subroutine calculate_energy_ccs_complex
!
!
   module subroutine calculate_energy_omega_term_ccs_complex(wf)
!!
!!    Calculate energy_complex omega term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds multipliers dot omega to the energy_complex,
!!
!!       energy_complex += Σ_μ tbar_μ Ω_μ,
!!
!!    which appears in the variational energy_complex expression <Λ|H|CC> when Omega ≠ 0.
!!    This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(:), allocatable :: multipliers, omega
!
!
      call mem%alloc(multipliers, wf%n_gs_amplitudes)
      call mem%alloc(omega, wf%n_gs_amplitudes)
!
      call wf%get_multipliers_complex(multipliers)
      call wf%construct_omega_complex(omega)
!
      wf%energy_complex = wf%energy_complex + our_zdotu(wf%n_gs_amplitudes, multipliers, 1, omega, 1)
!
      call mem%dealloc(multipliers, wf%n_gs_amplitudes)
      call mem%dealloc(omega, wf%n_gs_amplitudes)
!
   end subroutine calculate_energy_omega_term_ccs_complex
!
!
   module subroutine calculate_energy_length_dipole_term_ccs_complex(wf, electric_field)
!!
!!    Calculate energy_complex length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds dipole part of the length gauge electromagnetic potential to the energy_complex,
!!
!!       energy_complex += 2 sum_ii (-μ·E)_ii,
!!
!!    where μ is the vector of electric dipole integral matrices and E is a uniform classical electric
!!    vector field. This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      complex(dp), dimension(3), intent(in) :: electric_field
!
      complex(dp), dimension(:,:,:), allocatable :: mu
!
      integer :: i
!
!     Construct t1 transformed dipole moment
!
      call mem%alloc(mu, wf%n_mo, wf%n_mo, 3)
      call wf%construct_mu_complex(mu)
!
!     Add one_complex-electron electric field contribution to the diagonal of Fock and one_complex-electron integral terms
!
      do i = 1, wf%n_o
!
         wf%energy_complex = wf%energy_complex - two_complex*(mu(i, i, 1)*electric_field(1)   &
                                      + mu(i, i, 2)*electric_field(2) &
                                      + mu(i, i, 3)*electric_field(3))
!
      enddo
!
   end subroutine calculate_energy_length_dipole_term_ccs_complex
!
!
end submodule zop_ccs_complex
