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
submodule (ccs_class) mean_value_ccs
!
!!
!!    Mean-value submodule
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
   module subroutine prepare_for_properties_ccs(wf)
!!
!!    Prepare for properties
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      call do_nothing(wf)
!
   end subroutine prepare_for_properties_ccs
!
!
   module subroutine construct_gs_density_ccs(wf)
!!
!!    Construct GS density
!!    Written by Sarai Dery Folkestad, 2019
!!
!!    Constructs the one-electron density matrix in the T1 basis
!!
!!    D_pq = < Lambda| E_pq |CC >
!!
!!    Contributions to the density are split up as follows:
!!    D_pq = D_pq(ref-ref) + sum_mu tbar_mu D_pq(mu-ref)
!!
!!    The second term is separated in "mu_ref_density_terms",
!!    as it is used for the left transition density as well.
!!
      implicit none
!
      class(ccs) :: wf
      type(timings) :: timer
!
      timer = timings('Ground state density', pl='m')
      call timer%turn_on
!
      call wf%mu_ref_density_terms(wf%density, 0, wf%t1bar)
      call wf%density_ccs_ref_ref_oo(wf%density)
!
      call timer%turn_off
!
   end subroutine construct_gs_density_ccs
!
!
   module subroutine mu_ref_density_terms_ccs(wf, density, state, L)
!!
!!    Density mu ref terms
!!    Written by Alexander C. Paul, May 2021
!!
!!    Constructs terms of the form:
!!       sum_mu L_mu < mu| E_pq |HF >
!!
!!    corresponding to terms of the ground state density
!!    and the left transition density.
!!
      use array_initialization, only: zero_array
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: density
!
      integer, intent(in) :: state
!
      real(dp), dimension(wf%n_t1), intent(in) :: L
!
      type(timings)     :: timer
      character(len=40) :: timer_name
!
      write(timer_name, '(a,i0,a)') 'CCS contribution to <', state,'|E_pq|0>'
      timer = timings(trim(timer_name), pl='v')
      call timer%turn_on()
!
      call zero_array(density, wf%n_mo**2)
!
      call wf%density_ccs_mu_ref_vo(density, L)
!
      call timer%turn_off()
!
   end subroutine mu_ref_density_terms_ccs
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
!!         = E_hf + sum_aibj 2 t_i^a L^J_ia L^J_jb t_j^b
!!                - sum_aibj t_i^a L^J_ja L^J_ib t_j^b
!!
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(:,:,:), allocatable :: L_Jai, L_Jia, L_Jjb, X_Jij, X_Jji
      real(dp), dimension(:), allocatable :: X_J
!
      integer :: req0, req1_i, req1_j, req2, req_single_batch
!
      integer :: current_i_batch, current_j_batch
!
      type(batching_index) :: batch_i, batch_j
!
      type(timings), allocatable :: timer
!
      real(dp) :: ddot
      real(dp) :: correlation_energy
!
      integer :: i, j, JJ
!
      timer = timings('Calculate energy (ccs)', 'n')
      call timer%turn_on()

      call mem%alloc(X_J, wf%eri_t1%n_J, set_zero=.true.)
!
      req0 = 0
!
      req1_i = wf%eri_t1%n_J*wf%n_v
!
      batch_i = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, req0, req1_i, tag='calculate_energy_ccs 1')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call mem%alloc(L_Jai, wf%eri_t1%n_J, wf%n_v, batch_i%length)
         call wf%L_mo%get(L_Jai, wf%n_o + 1, wf%n_mo, &
                          batch_i%first, batch_i%get_last())
!
         call dgemm('N', 'N',                &
                     wf%eri_t1%n_J,          &
                     1,                      &
                     wf%n_v*batch_i%length,  &
                     one,                    &
                     L_Jai,                  &
                     wf%eri_t1%n_J,          &
                     wf%t1(1,batch_i%first), &
                     wf%n_v*wf%n_o,          &
                     one,                    &
                     X_J,                    &
                     wf%eri_t1%n_J)
!
         call mem%dealloc(L_Jai, wf%eri_t1%n_J, wf%n_v, batch_i%length)
!
      enddo
!
      call mem%batch_finalize()
!
      correlation_energy = two*ddot(wf%eri_t1%n_J, X_J, 1, X_J, 1)
!
      call mem%dealloc(X_J, wf%eri_t1%n_J)
!
      req0 = 0
!
      req1_i = wf%eri_t1%n_J*wf%n_v
      req1_j = wf%eri_t1%n_J*wf%n_v
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
!
      req2 = 3*wf%eri_t1%n_J
!
      req_single_batch = wf%eri_t1%n_J*wf%n_v*wf%n_o + wf%eri_t1%n_J*(wf%n_o**2)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2, &
                           req_single_batch=req_single_batch, tag='calculate_energy_ccs 2')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call mem%alloc(L_Jia, wf%eri_t1%n_J, batch_i%length, wf%n_v)
         call wf%L_mo%get(L_Jia,                                   &
                                     batch_i%first, batch_i%get_last(), &
                                     wf%n_o + 1, wf%n_mo)
 !
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
            call mem%alloc(X_Jij, wf%eri_t1%n_J, batch_i%length, batch_j%length)
!
            call dgemm('N', 'N',                      &
                        wf%eri_t1%n_J*batch_i%length, &
                        batch_j%length,               &
                        wf%n_v,                       &
                        one,                          &
                        L_Jia,                        &
                        wf%eri_t1%n_J*batch_i%length, &
                        wf%t1(1, batch_j%first),      &
                        wf%n_v,                       &
                        zero,                         &
                        X_Jij,                        &
                        wf%eri_t1%n_J*batch_i%length)
!
            if (current_j_batch .ne. current_i_batch) then
!
               call mem%alloc(L_Jjb, wf%eri_t1%n_J, batch_j%length, wf%n_v)
!
               call wf%L_mo%get(L_Jjb,                &
                                batch_j%first,        &
                                batch_j%get_last(),   &
                                wf%n_o + 1, wf%n_mo)
!
               call mem%alloc(X_Jji, wf%eri_t1%n_J, batch_j%length, batch_i%length)
!
               call dgemm('N', 'N',                      &
                           wf%eri_t1%n_J*batch_j%length, &
                           batch_i%length,               &
                           wf%n_v,                       &
                           one,                          &
                           L_Jjb,                        &
                           wf%eri_t1%n_J*batch_j%length, &
                           wf%t1(1, batch_i%first),      &
                           wf%n_v,                       &
                           zero,                         &
                           X_Jji,                        &
                           wf%eri_t1%n_J*batch_j%length)
!
               call mem%dealloc(L_Jjb, wf%eri_t1%n_J, batch_j%length, wf%n_v)
!
!$omp parallel do private(j, i, JJ) reduction(+:correlation_energy)
               do j = 1, batch_j%length
                  do i = 1, batch_i%length
                     do JJ = 1, wf%eri_t1%n_J
!
                        correlation_energy = correlation_energy - X_Jij(JJ,i,j)*X_Jji(JJ,j,i)
!
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
               call mem%dealloc(X_Jji, wf%eri_t1%n_J, batch_j%length, batch_i%length)
!
            else
!
!$omp parallel do private(j, i, JJ) reduction(+:correlation_energy)
               do j = 1, batch_j%length
                  do i = 1, batch_i%length
                     do JJ = 1, wf%eri_t1%n_J
!
                        correlation_energy = correlation_energy - X_Jij(JJ,i,j)*X_Jij(JJ,j,i)
!
                     enddo
                  enddo
               enddo
!$omp end parallel do
!
            endif
!
            call mem%dealloc(X_Jij, wf%eri_t1%n_J, batch_i%length, batch_j%length)
!
         enddo
!
         call mem%dealloc(L_Jia, wf%eri_t1%n_J, batch_i%length, wf%n_v)
!
      enddo
!
      call mem%batch_finalize()
!
      wf%correlation_energy = correlation_energy
!
      wf%energy = wf%hf_energy + wf%correlation_energy
!
      call timer%turn_off()
!
   end subroutine calculate_energy_ccs
!
!
   module function get_electronic_dipole_ccs(wf) result(mu_electronic)
!!
!!    Get electronic dipole
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2019-2021
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(3) :: mu_electronic
!
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: mu_pqk
!
!     Get the integrals mu_pqk for components k = 1, 2, 3 in the T1-transformed basis
!
      call mem%alloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
      call wf%get_t1_oei('dipole', mu_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 3
!
         mu_electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k), wf%density)
!
         if (wf%exists_frozen_fock_terms) then
!
            mu_electronic(k) = mu_electronic(k) + wf%frozen_dipole(k)
!
         endif
!
      enddo
!
      call mem%dealloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
!
   end function get_electronic_dipole_ccs
!
!
   module function get_electronic_quadrupole_ccs(wf) result(q_electronic)
!!
!!    Get electronic quadrupole
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(6) :: q_electronic
!
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: q_pqk
!
!     Get the integrals q_pqk for components k = 1, 2, ..., 6 in the T1-transformed basis
!
      call mem%alloc(q_pqk, wf%n_mo, wf%n_mo, 6)
      call wf%get_t1_oei('quadrupole', q_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 6
!
         q_electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k), wf%density)
!
         if (wf%exists_frozen_fock_terms) then
!
            q_electronic(k) = q_electronic(k) + wf%frozen_quadrupole(k)
!
         endif
!
      enddo
!
      call mem%dealloc(q_pqk, wf%n_mo, wf%n_mo, 6)
!
   end function get_electronic_quadrupole_ccs
!
!
end submodule mean_value_ccs
