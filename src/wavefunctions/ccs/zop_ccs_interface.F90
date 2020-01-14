!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
   module subroutine prepare_for_density_ccs(wf)
!!
!!    Prepare for the construction of density matrices
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
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
!!    D_pq = < Lambda| E_pq |CC >
!!
!!    Contributions to the density are split up as follows:
!!    D_pq = D_pq(ref-ref) + sum_mu tbar_mu D_pq(mu-ref)
!!
      implicit none
!
      class(ccs) :: wf
!
   end subroutine construct_gs_density_ccs
!
!
   module subroutine density_ccs_ref_ref_oo_ccs(wf, density)
!!
!!    One electron density reference-reference oo-term
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Hartree-Fock density contribution:
!!    D_pq += < HF| e^(-T) E_pq e^T |HF >
!!
!!    D_ii = 2  
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
!
   end subroutine density_ccs_ref_ref_oo_ccs
!
!
   module subroutine density_ccs_mu_ref_vo_ccs(wf, density, tbar_ai)
!!
!!    One electron density excited-determinant/reference vo-term
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Computes terms of the form:
!!
!!          D_pq += sum_mu X_mu * < mu| e^(-T) E_pq e^T |HF >
!!
!!    where X_mu is a general amplitude (tbar or L)
!!
!!    explicit term in this routine:
!!          D_ai = tbar_ai 
!!
      implicit none
!
      class(ccs) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: density
      real(dp), dimension(wf%n_v, wf%n_o) :: tbar_ai
!
   end subroutine density_ccs_mu_ref_vo_ccs
!
!
   module function calculate_expectation_value_ccs(wf, A, density) result(expectation_value)
!!
!!    Calculate expectation value
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Calculate the expectation value of a one-electron
!!    operator A
!!
!!       < A > = < Lambda| A | CC > = sum_pq A_pq D_pq
!!
!!    where A_pq are the T1-transformed integrals
!!    and D_pq is the a one-electron density matrix
!!    in the T1-basis
!!
      implicit none
!
      class(ccs), intent(in) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: A
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: density
      real(dp) :: expectation_value
!
   end function calculate_expectation_value_ccs
!
!
   module subroutine calculate_energy_ccs(wf)
!!
!!    Calculate energy
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the CCS energy. This is only equal to the actual
!!    energy when the ground state equations are solved, of course.
!!
!!       E = E_hf + sum_aibj t_i^a t_j^b L_iajb
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
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
!!       energy += sum_mu tbar_mu Omega_mu,
!!
!!    which appears in the energy expression:
!!
!!          < Lambda|H|CC > when Omega != 0.
!!
!!    This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
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
!!       energy += 2 sum_ii (-mu·E)_ii,
!!
!!    where mu is the vector of electric dipole integral matrices 
!!    and E is a uniform classical electric
!!    vector field. This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      real(dp), dimension(3), intent(in) :: electric_field
!
   end subroutine calculate_energy_length_dipole_term_ccs
