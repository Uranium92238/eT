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
submodule (ccs_class) fock_ccs
!
!!
!!    Fock submodule
!!
!!    Submodule containing routines that can be used to construct the t1-transformed Fock matrix.
!!
!
   implicit none
!
!
contains
!
!
   module subroutine construct_fock_ccs(wf, task)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad, Jul 2020
!!
!!    Constructs the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_pq = h_pq + sum_k (2*g_pqkk - g_pkkq) + (effective Fock contributions)
!!
!!    Depending on the 'task' different blocks (ij, ai, ia, ab) will be constructed
!!
      implicit none
!
      class(ccs), intent(inout)              :: wf
      character(len=*), intent(in), optional :: task
      type(timings) :: timer
!
      real(dp), dimension(:,:), allocatable :: h, F_eff
!
      timer = timings('Fock matrix construction (T1 basis)', pl='n')
      call timer%turn_on()
!
      call mem%alloc(h, wf%n_mo, wf%n_mo)
      call mem%alloc(F_eff, wf%n_mo, wf%n_mo, set_zero=.true.)
!
      call wf%get_t1_oei('hamiltonian', h, screening=.true.)
!
      if (wf%exists_frozen_fock_terms) call wf%add_frozen_fock_terms(F_eff)
!
      if (.not. present(task)) then
!
         call wf%construct_fock_ai_t1(h, F_eff)
         call wf%construct_fock_ia_t1(h, F_eff)
         call wf%construct_fock_ab_t1(h, F_eff)
         call wf%construct_fock_ij_t1(h, F_eff)
!
      else
!
         if (trim(task) == 'gs') then
!
            call wf%construct_fock_ai_t1(h, F_eff)
!
         elseif (trim(task) == 'multipliers') then
!
            call wf%construct_fock_ia_t1(h, F_eff)
            call wf%construct_fock_ab_t1(h, F_eff)
            call wf%construct_fock_ij_t1(h, F_eff)
!
         elseif (trim(task) == 'es') then
!
            call wf%construct_fock_ab_t1(h, F_eff)
            call wf%construct_fock_ij_t1(h, F_eff)
!
         else
!
            call output%error_msg('did not recognize task in construct_fock_ccs')
!
         endif
!
      endif
!
      call mem%dealloc(h, wf%n_mo, wf%n_mo)
      call mem%dealloc(F_eff, wf%n_mo, wf%n_mo)
!
      call timer%turn_off()
!
   end subroutine construct_fock_ccs
!
!
   module subroutine add_frozen_fock_terms_ccs(wf, F_pq)
!!
!!    Add frozen Fock terms
!!    Written by Sarai D. Folkestad, 2019
!!
!!    Adds the frozen core contributions to
!!    the effective T1-transformed Fock matrix.
!!
!!    Isolated into subroutine by Eirik F. Kjønstad, 2019
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: F_pq
!
      real(dp), dimension(:,:), allocatable :: F_pq_frozen
!
      call mem%alloc(F_pq_frozen, wf%n_mo, wf%n_mo)
!
      call wf%construct_t1_frozen_fock_terms(F_pq_frozen)
      call daxpy(wf%n_mo**2, one, F_pq_frozen, 1, F_pq, 1)
!
      call mem%dealloc(F_pq_frozen, wf%n_mo, wf%n_mo)
!
   end subroutine add_frozen_fock_terms_ccs
!
!
   module subroutine construct_t1_frozen_fock_terms_ccs(wf, F_pq)
!!
!!    Calculate T1 Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
      call dcopy(wf%n_mo**2, wf%mo_fock_frozen, 1, F_pq, 1)
!
      call wf%t1_transform(F_pq)
!
   end subroutine construct_t1_frozen_fock_terms_ccs
!
!
   module subroutine add_t1_fock_length_dipole_term_ccs(wf, electric_field)
!!
!!    Add t1 Fock length dipole term (CCS)
!!    Written by Andreas Skeidsvoll, Jan 2019
!!
!!    Adds dipole part of the length gauge electromagnetic potential to the Fock matrix,
!!
!!       Fock matrix += -μ·E,
!!
!!    where μ is the vector of electric dipole integral matrices and E is a uniform classical electric
!!    vector field. This routine does not have to be overwritten in descendants.
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
      real(dp), dimension(3), intent(in) :: electric_field
!
      real(dp), dimension(:,:,:), allocatable :: mu, potential
!
      integer :: a, i, b, j
!
      call mem%alloc(mu, wf%n_mo, wf%n_mo, 3)
      call wf%get_t1_oei('dipole', mu)
!
!     Create interaction potential by scaling dipole integrals by minus electric field
!
      call mem%alloc(potential, wf%n_mo, wf%n_mo, 3, set_zero=.true.)
!
      call daxpy((wf%n_mo)**2, -electric_field(1), mu(:,:,1), 1, potential(:,:,1), 1)
      call daxpy((wf%n_mo)**2, -electric_field(2), mu(:,:,2), 1, potential(:,:,2), 1)
      call daxpy((wf%n_mo)**2, -electric_field(3), mu(:,:,3), 1, potential(:,:,3), 1)
!
      call mem%dealloc(mu, wf%n_mo, wf%n_mo, 3)
!
!$omp parallel do private(i,j)
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            wf%fock_ij(i, j) = wf%fock_ij(i, j)     &
                               + potential(i, j, 1) &
                               + potential(i, j, 2) &
                               + potential(i, j, 3)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            wf%fock_ia(i, a) = wf%fock_ia(i, a)              &
                               + potential(i, wf%n_o + a, 1) &
                               + potential(i, wf%n_o + a, 2) &
                               + potential(i, wf%n_o + a, 3)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(i,a)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ai(a, i) = wf%fock_ai(a, i)              &
                               + potential(wf%n_o + a, i, 1) &
                               + potential(wf%n_o + a, i, 2) &
                               + potential(wf%n_o + a, i, 3)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private(a,b)
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            wf%fock_ab(a, b) = wf%fock_ab(a, b)                       &
                               + potential(wf%n_o + a, wf%n_o + b, 1) &
                               + potential(wf%n_o + a, wf%n_o + b, 2) &
                               + potential(wf%n_o + a, wf%n_o + b, 3)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(potential, wf%n_mo, wf%n_mo, 3)
!
   end subroutine add_t1_fock_length_dipole_term_ccs
!
!
   module subroutine construct_fock_ai_t1_ccs(wf, h, F_eff)
!!
!!    Construct Fock ai T1,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the ai block of the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_ai = h_ai + sum_j (2*g_aijj - g_ajji) + (effective Fock contributions)
!!
!!    Effective Fock contributions:
!!
!!       Frozen core by Sarai D. Folkestad, 2019
!!       QM/MM by Tommaso Giovannini, 2019
!!       QM/PCM by Tommaso Giovannini, 2019
!!
!!       Modified by Sarai D. Folkestad, Nov 2019
!!
!!       - Added batching for N^2 memory requirement.
!!
!!       Modified by Sarai D. Folkestad, Dec 2021
!!
!!       - N^4 scaling.
!!
      use batching_index_class, only : batching_index
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: h
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_eff
!
      integer :: i, j, a, jj
!
      real(dp), dimension(:), allocatable :: X_J
      real(dp), dimension(:,:,:), allocatable :: L_Jai
      real(dp), dimension(:,:,:), allocatable :: L_Jaj, L_Jji, L_Jja
!
      integer :: req0, req1_i, current_i_batch
      integer :: req1_j, current_j_batch
!
      type(batching_index) :: batch_i, batch_j
!
!     Set Fock matrix to h + effective Fock contributions
!
!$omp parallel do private(a, i)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ai(a,i) = h(a + wf%n_o, i) + F_eff(a + wf%n_o, i)
!
         enddo
      enddo
!$omp end parallel do
!
!     Add occupied-virtual contributions: F_ai = F_ai + sum_j (2*g_aijj - g_ajji)
!
      call mem%alloc(X_J, wf%L_t1%n_J, set_zero=.true.)
!
!     - Exchange term: L_aj^J L^J_ji
!
      req0 = 0
      req1_j = max(2*wf%L_t1%n_J*wf%n_v, wf%L_t1%n_J*wf%n_v + wf%L_t1%n_J*wf%n_o)
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, tag='construct_fock_ai_t1_ccs_3')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jaj, wf%L_t1%n_J, wf%n_v, batch_j%length)
         call wf%L_t1%get(L_Jaj,                    &
                          wf%n_o + 1, wf%n_mo,      &
                          batch_j%first, batch_j%get_last())
!
         call mem%alloc(L_Jja, wf%L_t1%n_J, batch_j%length, wf%n_v)
         call sort_123_to_132(L_Jaj, L_Jja, wf%L_t1%n_J, wf%n_v, batch_j%length)
         call mem%dealloc(L_Jaj, wf%L_t1%n_J, wf%n_v, batch_j%length)
!
         call mem%alloc(L_Jji, wf%L_t1%n_J, batch_j%length, wf%n_o)
!
         call wf%L_t1%get(L_Jji,                              &
                          batch_j%first, batch_j%get_last(),  &
                          1, wf%n_o)
!
         call dgemm('T', 'N',                      &
                     wf%n_v,                       &
                     wf%n_o,                       &
                     wf%L_t1%n_J*batch_j%length,   &
                     -one,                         &
                     L_Jja,                        &
                     wf%L_t1%n_J*batch_j%length,   &
                     L_Jji,                        &
                     wf%L_t1%n_J*batch_j%length,   &
                     one,                          &
                     wf%fock_ai,                   &
                     wf%n_v)
!
!     Construct X_J = sum_j L^J_jj
!
!$omp parallel do private (J, jj)
         do J = 1, wf%L_t1%n_J
            do jj = 1, batch_j%length
!
               X_J(J) = X_J(J) + L_Jji(J, jj, jj + batch_j%first - 1)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_Jji, wf%L_t1%n_J, batch_j%length, wf%n_o)
         call mem%dealloc(L_Jja, wf%L_t1%n_J, batch_j%length, wf%n_v)
!
      enddo
!
      call mem%batch_finalize()
!
!
!     - Coulomb term: 2 g_aijj = 2 L_ai^J L^J_jj
!
!     F_ai += 2 L_Jai X_J
!
!     batching over i
!
      req0 = 0
      req1_i = wf%L_t1%n_J*wf%n_v
!
      batch_i = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, req0, req1_i, tag='construct_fock_ai_t1_ccs_2')
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         call mem%alloc(L_Jai, wf%L_t1%n_J, wf%n_v, batch_i%length)
!
         call wf%L_t1%get(L_Jai,                 &
                          wf%n_o + 1, wf%n_mo,   &
                          batch_i%first, batch_i%get_last())
!
         call dgemm('T', 'N',                      &
                     wf%n_v*batch_i%length,        &
                     1,                            &
                     wf%L_t1%n_J,                  &
                     two,                          &
                     L_Jai,                        &
                     wf%L_t1%n_J,                  &
                     X_J,                          &
                     wf%L_t1%n_J,                  &
                     one,                          &
                     wf%fock_ai(1,batch_i%first),  &
                     wf%n_v*wf%n_o)
!
         call mem%dealloc(L_Jai, wf%L_t1%n_J, wf%n_v, batch_i%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(X_J, wf%L_t1%n_J)
!
   end subroutine construct_fock_ai_t1_ccs
!
!
   module subroutine construct_fock_ia_t1_ccs(wf, h, F_eff, first_i, last_i, first_a, last_a)
!!
!!    Construct Fock ia T1,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the ia block of the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_ia = h_ia + sum_j (2*g_iajj - g_ijja) + (effective Fock contributions)
!!
!!    Effective Fock contributions:
!!
!!       Frozen core by Sarai D. Folkestad, 2019
!!       QM/MM by Tommaso Giovannini, 2019
!!       QM/PCM by Tommaso Giovannini, 2019
!!
!!       Modified by Sarai D. Folkestad, Nov 2019
!!
!!       Added batching for N^2 memory requirement.
!!
!!       Modified by Sarai D. Folkestad, Jul 2020
!!
!!       Added limits for a and i for
!!       subblock calculation
!!
!
      use batching_index_class, only : batching_index
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: h
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_eff
!
      integer, intent(in), optional :: first_i, last_i, first_a, last_a
!
      type(range_) :: i_range, a_range
!!
      integer :: i, j, a, jj
!
      real(dp), dimension(:), allocatable :: X_J
      real(dp), dimension(:,:), allocatable :: F_ia
      real(dp), dimension(:,:,:), allocatable :: L_Jjj, L_Jia, L_Jij, L_Jji, L_Jja
!
      integer :: req0, req1_a, current_a_batch
      integer :: req1_j, current_j_batch
!
      type(batching_index) :: batch_j, batch_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_a) .and. present(last_a)) then
!
         i_range = range_(first_i, last_i)
         a_range = range_(first_a, last_a)
!
      else
!
         i_range = range_(1, wf%n_o)
         a_range = range_(1, wf%n_v)
!
      endif
!
!     Set Fock matrix to h + effective Fock contributions
!
!$omp parallel do private(a, i)
      do a = a_range%first, a_range%get_last()
         do i = i_range%first, i_range%get_last()
!
            wf%fock_ia(i, a) = h(i, a + wf%n_o) + F_eff(i, a + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
!     Add occupied-virtual contributions: F_ia = F_ia + sum_j (2*g_iajj - g_ijja)
!
!     - Coulomb part: 2 g_iajj = L_Jia L_Jjj
!
!     Construct X_J = sum_j L^J_jj
!
      call mem%alloc(X_J, wf%L_t1%n_J, set_zero=.true.)
!
      req0 = 0
      req1_j = wf%L_t1%n_J*wf%n_o ! NOTE: this is an overestimate but used due to
                                 !       the double j index of L_Jjj
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, tag='construct_fock_ia_t1_ccs_1')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jjj, wf%L_t1%n_J, batch_j%length, batch_j%length)
         call wf%L_t1%get(L_Jjj,                               &
                          batch_j%first, batch_j%get_last(),   &
                          batch_j%first, batch_j%get_last())
!
!$omp parallel do private (J, jj)
         do J = 1, wf%L_t1%n_J
            do jj = 1, batch_j%length
!
               X_J(J) = X_J(J) + L_Jjj(J, jj, jj)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_Jjj, wf%L_t1%n_J, batch_j%length, batch_j%length)
!
      enddo
!
      call mem%batch_finalize()
!
!     F_ia += 2 L_Jia X_J
!
      call mem%alloc(F_ia, i_range%length, a_range%length, set_zero=.true.)
!
      req0 = 0
      req1_a = wf%L_t1%n_J*i_range%length
!
      batch_a = batching_index(a_range%length)
!
      call mem%batch_setup(batch_a, req0, req1_a, tag='construct_fock_ia_t1_ccs_2')
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(L_Jia, wf%L_t1%n_J, i_range%length, batch_a%length)
         call wf%L_t1%get(L_Jia,                                         &
                          i_range%first, i_range%get_last(),             &
                          wf%n_o + batch_a%first + a_range%first - 1,    &
                          wf%n_o + batch_a%get_last() + a_range%first - 1)
!
         call dgemm('T', 'N',                         &
                     i_range%length*batch_a%length,   &
                     1,                               &
                     wf%L_t1%n_J,                     &
                     two,                             &
                     L_Jia,                           &
                     wf%L_t1%n_J,                     &
                     X_J,                             &
                     wf%L_t1%n_J,                     &
                     one,                             &
                     F_ia(1, batch_a%first),          &
                     i_range%length*batch_a%length)
!
         call mem%dealloc(L_Jia, wf%L_t1%n_J, i_range%length, batch_a%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(X_J, wf%L_t1%n_J)
!
!     - Exchange term: - L_ij^J L^J_ja
!
      req0 = 0
      req1_j = max(2*wf%L_t1%n_J*i_range%length, &
                  wf%L_t1%n_J*a_range%length + wf%L_t1%n_J*i_range%length)
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, tag='construct_fock_ia_t1_ccs_3')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jij, wf%L_t1%n_J, i_range%length, batch_j%length)
!
         call wf%L_t1%get(L_Jij,                             &
                          i_range%first, i_range%get_last(), &
                          batch_j%first, batch_j%get_last())
!
         call mem%alloc(L_Jji, wf%L_t1%n_J, batch_j%length, i_range%length)
         call sort_123_to_132(L_Jij, L_Jji, wf%L_t1%n_J, i_range%length, batch_j%length)
         call mem%dealloc(L_Jij, wf%L_t1%n_J, i_range%length, batch_j%length)
!
         call mem%alloc(L_Jja, wf%L_t1%n_J, batch_j%length, a_range%length)
!
         call wf%L_t1%get(L_Jja,                              &
                          batch_j%first, batch_j%get_last(),  &
                          wf%n_o + a_range%first, wf%n_o + a_range%get_last())

!
         call dgemm('T', 'N',                   &
                     i_range%length,            &
                     a_range%length,            &
                     wf%L_t1%n_J*batch_j%length,&
                     -one,                      &
                     L_Jji,                     &
                     wf%L_t1%n_J*batch_j%length,&
                     L_Jja,                     &
                     wf%L_t1%n_J*batch_j%length,&
                     one,                       &
                     F_ia,                      &
                     i_range%length)
!
         call mem%dealloc(L_Jji, wf%L_t1%n_J, batch_j%length, i_range%length)
         call mem%dealloc(L_Jja, wf%L_t1%n_J, batch_j%length, a_range%length)
!
      enddo
!
      call mem%batch_finalize()
!
!$omp parallel do private(a, i)
      do a = 1, a_range%length
         do i = 1, i_range%length
!
            wf%fock_ia(i + i_range%first - 1, a + a_range%first - 1) = &
               wf%fock_ia(i + i_range%first - 1, a + a_range%first - 1) + F_ia(i, a)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(F_ia, i_range%length, a_range%length)
!
   end subroutine construct_fock_ia_t1_ccs
!
!
   module subroutine construct_fock_ab_t1_ccs(wf, h, F_eff, first_a, last_a, first_b, last_b)
!!
!!    Construct Fock ab T1,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the ab block of the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_ab = h_ab + sum_i (2*g_abii - g_aiib) + (effective Fock contributions)
!!
!!    Effective Fock contributions:
!!
!!       Frozen core by Sarai D. Folkestad, 2019
!!       QM/MM by Tommaso Giovannini, 2019
!!       QM/PCM by Tommaso Giovannini, 2019
!!
!!       Modified by Sarai D. Folkestad, Nov 2019
!!
!!       Added batching for N^2 memory requirement.
!!
!!       Modified by Sarai D. Folkestad, Jul 2020
!!
!!       Added limits for a and b for
!!       subblock calculation
!!
!
      use batching_index_class, only : batching_index
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: h
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_eff
!
      integer, intent(in), optional :: first_a, last_a, first_b, last_b
!
      type(range_) :: a_range, b_range
!
      integer :: a, b, jj, J
!
      real(dp), dimension(:,:,:), allocatable :: L_Jab, L_Jjj, L_Jaj, L_Jjb, L_Jja
      real(dp), dimension(:,:), allocatable :: F_ab
      real(dp), dimension(:), allocatable :: X_J
!
      integer :: req0, req1_j, current_j_batch
      integer :: req1_b, current_b_batch
!
      type(batching_index) :: batch_j, batch_b
!
      if (present(first_a) .and. present(last_a) .and. &
          present(first_b) .and. present(last_b)) then
!
         a_range = range_(first_a, last_a)
         b_range = range_(first_b, last_b)
!
      else
!
         a_range = range_(1, wf%n_v)
         b_range = range_(1, wf%n_v)
!
      endif
!
!     Set Fock matrix to h + effective Fock contributions
!
!$omp parallel do private(a, b)
      do b = b_range%first, b_range%get_last()
         do a = a_range%first, a_range%get_last()
!
            wf%fock_ab(a, b) = h(a + wf%n_o, b + wf%n_o) + F_eff(a + wf%n_o, b + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do
!
!     Add virtual-virtual contributions: F_ab = h_ab + sum_i (2*g_abii - g_aiib)
!
!     - Coulomb part: 2 g_abjj = L_Jab L_Jjj
!
!     Construct X_J = sum_j L^J_jj
!
      call mem%alloc(X_J, wf%L_t1%n_J, set_zero=.true.)
!
!     batching over j
!
      req0 = 0
      req1_j = wf%L_t1%n_J*wf%n_o ! NOTE: this estimate is not correct but used due to
                                       !       the double j index of L_Jjj
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, &
                           tag='construct_fock_ab_t1_ccs_1')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jjj, wf%L_t1%n_J, batch_j%length, batch_j%length)
         call wf%L_t1%get(L_Jjj,                              &
                          batch_j%first, batch_j%get_last(),  &
                          batch_j%first, batch_j%get_last())
!
!$omp parallel do private (J, jj)
         do J = 1, wf%L_t1%n_J
            do jj = 1, batch_j%length
!
               X_J(J) = X_J(J) + L_Jjj(J, jj, jj)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_Jjj, wf%L_t1%n_J, batch_j%length, batch_j%length)
!
      enddo
!
      call mem%batch_finalize()
!
!     F_ab += 2 L_Jab X_J
!
!     batching over a
!
      call mem%alloc(F_ab, a_range%length, b_range%length, set_zero=.true.)
!
      req0 = 0
      req1_b = wf%L_t1%n_J*a_range%length
!
      batch_b = batching_index(b_range%length)
!
      call mem%batch_setup(batch_b, req0, req1_b, &
                           tag='construct_fock_ab_t1_ccs_2')
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
         call mem%alloc(L_Jab, wf%L_t1%n_J, a_range%length, batch_b%length)
         call wf%L_t1%get(L_Jab,                                       &
                          wf%n_o + a_range%first,                      &
                          wf%n_o + a_range%get_last(),                 &
                          wf%n_o + batch_b%first + b_range%first - 1,  &
                          wf%n_o + batch_b%get_last() + b_range%first - 1)
!
         call dgemm('T', 'N',                         &
                     a_range%length*batch_b%length,   &
                     1,                               &
                     wf%L_t1%n_J,                     &
                     two,                             &
                     L_Jab,                           &
                     wf%L_t1%n_J,                     &
                     X_J,                             &
                     wf%L_t1%n_J,                     &
                     one,                             &
                     F_ab(1, batch_b%first),          &
                     a_range%length*b_range%length)
!
         call mem%dealloc(L_Jab, wf%L_t1%n_J, a_range%length, batch_b%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(X_J, wf%L_t1%n_J)
!
!     - Exchange term: - L_aj^J L^J_jb
!
!     batching over j
!
      req0 = 0
      req1_j = max(2*wf%L_t1%n_J*a_range%length, &
                  wf%L_t1%n_J*a_range%length + wf%L_t1%n_J*b_range%length)
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, &
                           tag='construct_fock_ab_t1_ccs_3')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jaj, wf%L_t1%n_J, a_range%length, batch_j%length)
!
         call wf%L_t1%get(L_Jaj,                         &
                          wf%n_o + a_range%first,        &
                          wf%n_o + a_range%get_last(),   &
                          batch_j%first, batch_j%get_last())
!
         call mem%alloc(L_Jja, wf%L_t1%n_J, batch_j%length, a_range%length)
         call sort_123_to_132(L_Jaj, L_Jja, wf%L_t1%n_J, a_range%length, batch_j%length)
         call mem%dealloc(L_Jaj, wf%L_t1%n_J, a_range%length, batch_j%length)
!
         call mem%alloc(L_Jjb, wf%L_t1%n_J, batch_j%length, b_range%length)
!
         call wf%L_t1%get(L_Jjb,                              &
                          batch_j%first, batch_j%get_last(),  &
                          wf%n_o + b_range%first,             &
                          wf%n_o + b_range%get_last())

!
         call dgemm('T', 'N',                   &
                     a_range%length,            &
                     b_range%length,            &
                     wf%L_t1%n_J*batch_j%length,&
                     -one,                      &
                     L_Jja,                     &
                     wf%L_t1%n_J*batch_j%length,&
                     L_Jjb,                     &
                     wf%L_t1%n_J*batch_j%length,&
                     one,                       &
                     F_ab,                      &
                     a_range%length)
!
         call mem%dealloc(L_Jjb, wf%L_t1%n_J, batch_j%length, b_range%length)
         call mem%dealloc(L_Jja, wf%L_t1%n_J, batch_j%length, a_range%length)
!
      enddo
!
      call mem%batch_finalize()
!
!$omp parallel do private(a, b)
      do a = 1, a_range%length
         do b = 1, b_range%length
!
            wf%fock_ab(a + a_range%first - 1, b + b_range%first - 1) = &
               wf%fock_ab(a + a_range%first - 1, b + b_range%first - 1) + F_ab(a, b)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(F_ab, a_range%length, b_range%length)
!
   end subroutine construct_fock_ab_t1_ccs
!
!
   module subroutine construct_fock_ij_t1_ccs(wf, h, F_eff, first_i, last_i, first_j, last_j)
!!
!!    Construct Fock ij T1,
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the ij block of the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_ij = h_ij + sum_k (2*g_ijkk - g_ikkj) + (effective Fock contributions)
!!
!!    Effective Fock contributions:
!!
!!       Frozen core by Sarai D. Folkestad, 2019
!!       QM/MM by Tommaso Giovannini, 2019
!!       QM/PCM by Tommaso Giovannini, 2019
!!
!!       Modified by Sarai D. Folkestad, Nov 2019
!!
!!       Added batching for N^2 memory requirement.
!!
!!       Modified by Sarai D. Folkestad, Jul 2020
!!
!!       Added limits for i and j for
!!       subblock calculation
!!
      use batching_index_class, only : batching_index
      use reordering, only: sort_123_to_132
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: h
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_eff
!
      integer, intent(in), optional :: first_i, last_i, first_j, last_j
!
      type(range_) :: i_range, j_range
!
      integer :: i, j, jj
!
      real(dp), dimension(:), allocatable :: X_J
      real(dp), dimension(:,:), allocatable :: F_ik
      real(dp), dimension(:,:,:), allocatable :: L_Jik, L_Jjj, L_Jji, L_Jjk, L_Jij
!
      integer :: req0, req1_k, req1_j, current_k_batch, current_j_batch
!
      type(batching_index) :: batch_k, batch_j
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j)) then
!
         i_range = range_(first_i, last_i)
         j_range = range_(first_j, last_j)
!
      else
!
         i_range = range_(1, wf%n_o)
         j_range = range_(1, wf%n_o)
!
      endif
!
!     Set Fock matrix to h + effective Fock contributions
!
!$omp parallel do
      do j = j_range%first, j_range%get_last()
         do i = i_range%first, i_range%get_last()
!
            wf%fock_ij(i,j) = h(i,j) + F_eff(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
!     Add occupied-occupied contributions: F_ik = F_ik + sum_j (2*g_ikjj - g_ijjk)
!
!     - Coulomb part: 2 g_ikjj = L_Jik L_Jjj
!
!     Construct X_J = sum_j L^J_jj
!
      call mem%alloc(X_J, wf%L_t1%n_J, set_zero=.true.)
!
      req0 = 0
      req1_j = wf%L_t1%n_J*wf%n_o ! NOTE: this estimate is not correct but used due to
                                 !       the double j index of L_Jjj
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, &
                           tag='construct_fock_ij_t1_ccs_1')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jjj, wf%L_t1%n_J, batch_j%length, batch_j%length)
         call wf%L_t1%get(L_Jjj,                              &
                          batch_j%first, batch_j%get_last(),  &
                          batch_j%first, batch_j%get_last())
!
!$omp parallel do private (J, jj)
         do J = 1, wf%L_t1%n_J
            do jj = 1, batch_j%length
!
               X_J(J) = X_J(J) + L_Jjj(J, jj, jj)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(L_Jjj, wf%L_t1%n_J, batch_j%length, batch_j%length)
!
      enddo
!
      call mem%batch_finalize()
!
!     F_ik += 2 L_Jik X_J
!
!     batching over a
!
      call mem%alloc(F_ik, i_range%length, j_range%length, set_zero=.true.)
!
      req0 = 0
      req1_k = wf%L_t1%n_J*i_range%length
!
      batch_k = batching_index(j_range%length)
!
      call mem%batch_setup(batch_k, req0, req1_k, &
                           tag='construct_fock_ij_t1_ccs_2')
!
      do current_k_batch = 1, batch_k%num_batches
!
         call batch_k%determine_limits(current_k_batch)
!
         call mem%alloc(L_Jik, wf%L_t1%n_J, i_range%length, batch_k%length)
         call wf%L_t1%get(L_Jik,                             &
                          i_range%first,                     &
                          i_range%get_last(),                &
                          batch_k%first + j_range%first - 1, &
                          batch_k%get_last() + j_range%first - 1)
!
         call dgemm('T', 'N',                         &
                     i_range%length*batch_k%length,   &
                     1,                               &
                     wf%L_t1%n_J,                     &
                     two,                             &
                     L_Jik,                           &
                     wf%L_t1%n_J,                     &
                     X_J,                             &
                     wf%L_t1%n_J,                     &
                     one,                             &
                     F_ik(1, batch_k%first),          &
                     i_range%length*j_range%length)
!
         call mem%dealloc(L_Jik, wf%L_t1%n_J, i_range%length, batch_k%length)
!
      enddo
!
      call mem%batch_finalize()
!
      call mem%dealloc(X_J, wf%L_t1%n_J)
!
!     - Exchange term: - L_ij^J L^J_jk
!
!     batching over j
!
      req0 = 0
      req1_j = max(2*wf%L_t1%n_J*i_range%length, &
                  wf%L_t1%n_J*i_range%length + wf%L_t1%n_J*j_range%length)
!
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_j, req0, req1_j, &
                           tag='construct_fock_ij_t1_ccs_3')
!
      do current_j_batch = 1, batch_j%num_batches
!
         call batch_j%determine_limits(current_j_batch)
!
         call mem%alloc(L_Jij, wf%L_t1%n_J, i_range%length, batch_j%length)
         call wf%L_t1%get(L_Jij,               &
                          i_range%first,       &
                          i_range%get_last(),  &
                          batch_j%first, batch_j%get_last())
!
         call mem%alloc(L_Jji, wf%L_t1%n_J, batch_j%length, i_range%length)
         call sort_123_to_132(L_Jij, L_Jji, wf%L_t1%n_J, i_range%length, batch_j%length)
         call mem%dealloc(L_Jij, wf%L_t1%n_J, i_range%length, batch_j%length)
!
         call mem%alloc(L_Jjk, wf%L_t1%n_J, batch_j%length, j_range%length)
         call wf%L_t1%get(L_Jjk,                                   &
                          batch_j%first, batch_j%get_last(),  &
                          j_range%first,                      &
                          j_range%get_last())

!
         call dgemm('T', 'N',                   &
                     i_range%length,            &
                     j_range%length,            &
                     wf%L_t1%n_J*batch_j%length,&
                     -one,                      &
                     L_Jji,                     &
                     wf%L_t1%n_J*batch_j%length,&
                     L_Jjk,                     &
                     wf%L_t1%n_J*batch_j%length,&
                     one,                       &
                     F_ik,                      &
                     i_range%length)
!
         call mem%dealloc(L_Jjk, wf%L_t1%n_J, batch_j%length, j_range%length)
         call mem%dealloc(L_Jji, wf%L_t1%n_J, batch_j%length, i_range%length)
!
      enddo
!
      call mem%batch_finalize()
!
!$omp parallel do private(i, j)
      do i = 1, i_range%length
         do j = 1, j_range%length
!
            wf%fock_ij(i + i_range%first - 1, j + j_range%first - 1) = &
               wf%fock_ij(i + i_range%first - 1, j + j_range%first - 1) + F_ik(i, j)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(F_ik, i_range%length, j_range%length)
!
   end subroutine construct_fock_ij_t1_ccs
!
end submodule fock_ccs
