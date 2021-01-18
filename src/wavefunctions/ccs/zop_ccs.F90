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
submodule (ccs_class) zop_ccs
!
!!
!!    Zeroth order properties submodule 
!!
!!    Contains routines related to the mean values, i.e. 
!!    the construction of density matrices as well as expectation 
!!    value calculation.
!!
!!    The ground state density is constructed as follows:
!!
!!          D_pq = < Lambda| E_pq |CC >
!!    where: 
!!          < Lambda| = < HF| + sum_mu tbar_mu < mu| exp(-T)
!!
!!
!!    In general a CC density matrix can be written as:
!!
!!          D_pq = < X| e^(-T) E_pq e^T |Y >
!!
!!    where X and Y are left and right vectors with contributions from
!!    a reference determinant and excited determinants (< mu|, |nu >):
!!
!!          D_pq =             X_ref < HF| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_mu  < mu| e^(-T) E_pq e^T |HF >  Y_ref
!!                 + sum_mu    X_ref < HF| e^(-T) E_pq e^T |mu >  Y_mu
!!                 + sum_mu,nu X_mu  < mu| e^(-T) E_pq e^T |nu >  Y_nu
!!
!!    Depending on the type of density matrix (Ground state, transition , 
!!    excited state, interstate transition) different states and thus different
!!    amplitudes X_ref, X_mu, Y_ref and Y_mu will contribute.
!!
!!    In EOM theory the states can be written as the following vectors:
!!
!!          |CC >     = R_0 = (1, 0)
!!          |Lambda > = L_0 = (1, tbar_mu)
!!          |R_k >    = R_k = (-sum_mu(tbar_mu*R_mu), R_mu)
!!          |L_k >    = L_k = (0, L_mu)
!!
!!    The routine names derive from the contribution of the vectors:
!!
!!       ref_ref: first component of the vector for the left and right state
!!
!!       mu_ref:  second component of the vector for the left and 
!!                first component of the vector for the right state
!!
!!       ref_mu:  first component of the vector for the left and 
!!                second component of the vector for the right state
!!
!!       mu_nu:   second component of the vector for the left and right state
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
      call output%printf('v', '- No preparations for the density for ' //  &
                         trim(wf%name_) // ' wavefunction.', fs='(/t3,a)')
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
      call zero_array(wf%density, (wf%n_mo)**2)
!
      call wf%density_ccs_ref_ref_oo(wf%density)
      call wf%density_ccs_mu_ref_vo(wf%density, wf%t1bar)
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
!!    Calculates the CCS energy. This is only equal to the actual
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
      real(dp) :: omp_correlation_energy
!
      integer :: a, i, b, j
!
      integer :: req0, req1_i, req1_j, req2
!
      integer :: current_i_batch, current_j_batch
!
      type(batching_index) :: batch_i, batch_j
!
      req0 = 0
!
      req1_i = (wf%n_v)*(wf%eri%n_J)
      req1_j = (wf%n_v)*(wf%eri%n_J)
!
      req2 = (wf%n_v**2)
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2)
!
      omp_correlation_energy = zero
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(g_iajb, batch_i%length, wf%n_v, batch_j%length, wf%n_v)
!
            call wf%eri%get_eri_t1('ovov', g_iajb, &
                                   batch_i%first, batch_i%last, &
                                   1, wf%n_v, &
                                   batch_j%first, batch_j%last, &
                                   1, wf%n_v)
!
!$omp parallel do private(b,i,j,a) reduction(+:omp_correlation_energy)
            do b = 1, wf%n_v
               do i = 1, batch_i%length
                  do j = 1, batch_j%length
                     do a = 1, wf%n_v
!
                        omp_correlation_energy = omp_correlation_energy +     &
                                              wf%t1(a, i + batch_i%first - 1) &
                                             *wf%t1(b, j + batch_j%first - 1) &
                                             *(two*g_iajb(i, a, j, b)-g_iajb(i, b, j, a))
!
                     enddo
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_iajb, batch_i%length, wf%n_v, batch_j%length, wf%n_v)
!
         enddo
      enddo
!
      wf%correlation_energy = omp_correlation_energy
!
      wf%energy = wf%hf_energy + wf%correlation_energy
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
!!       energy += 2 sum_ii (-mu·E)_ii,
!!
!!    where mu is the vector of electric dipole integral matrices 
!!    and E is a uniform classical electric
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
!     Add one-electron electric field contribution to the diagonal 
!     of Fock and one-electron integral terms
!
      do i = 1, wf%n_o
!
         wf%energy = wf%energy - two*(mu(i, i, 1)*electric_field(1)   &
                                      + mu(i, i, 2)*electric_field(2) &
                                      + mu(i, i, 3)*electric_field(3))
!
      enddo
!
      call mem%dealloc(mu, wf%n_mo, wf%n_mo, 3)
!
   end subroutine calculate_energy_length_dipole_term_ccs
!
!
end submodule zop_ccs
