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
   module subroutine construct_fock_ccs(wf)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_pq = h_pq + sum_k (2*g_pqkk - g_pkkq) + (effective Fock contributions)
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
!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_pq
!
      integer :: i, j, k, a, b
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ijkl
      real(dp), dimension(:,:,:,:), allocatable :: g_iajj
      real(dp), dimension(:,:,:,:), allocatable :: g_ijja
      real(dp), dimension(:,:,:,:), allocatable :: g_aijj
      real(dp), dimension(:,:,:,:), allocatable :: g_ajji
      real(dp), dimension(:,:,:,:), allocatable :: g_abii
      real(dp), dimension(:,:,:,:), allocatable :: g_aiib
!
      integer :: req0, req2, req1_i, req1_k, current_i_batch, current_k_batch
      integer :: req1_j, current_j_batch, req1_a, current_a_batch
!
      type(batching_index) :: batch_i, batch_k, batch_j, batch_a
!
!     Set F_pq = h_pq (t1-transformed) 
!
      call mem%alloc(F_pq, wf%n_mo, wf%n_mo)
!
      call wf%construct_h(F_pq)
!
!     Add effective contributions to Fock matrix 
!
      if (wf%exists_frozen_fock_terms) call wf%add_frozen_fock_terms(F_pq)
!
!     Add occupied-occupied contributions: F_ij = F_ij + sum_k (2*g_ijkk - g_ikkj)
!
!     Batching over i and k
!
      req0 = 0
!
      req1_i = (wf%eri%n_J)*(wf%n_o)
      req1_k = (wf%eri%n_J)*(wf%n_o)
!
      req2 =  (wf%n_o**2)
!
      batch_i = batching_index(wf%n_o)
      batch_k = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_k, req0, req1_i, req1_k, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_k_batch = 1, batch_k%num_batches
!
            call batch_k%determine_limits(current_k_batch)
!
            call mem%alloc(g_ijkl, batch_i%length, wf%n_o, batch_k%length, wf%n_o)
!
            call wf%eri%get_eri_t1('oooo', g_ijkl,               &
                                   batch_i%first, batch_i%last,  &
                                   1, wf%n_o,                    &
                                   batch_k%first, batch_k%last,  &
                                   1, wf%n_o)
!
!$omp parallel do private(i,j,k)
            do j = 1, wf%n_o
               do i = 1, batch_i%length
                  do k = 1, batch_k%length
!
                     F_pq(i + batch_i%first - 1, j) = F_pq(i + batch_i%first - 1, j) &
                                 + two*g_ijkl(i, j, k, k + batch_k%first - 1) &
                                 - g_ijkl(i, k + batch_k%first - 1, k, j)
!
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_ijkl, batch_i%length, wf%n_o, batch_k%length, wf%n_o)
!
         enddo
!
      enddo
!
!     Add occupied-virtual contributions: F_ia = F_ia + sum_j (2*g_iajj - g_ijja)
!                                         F_ai = F_ai + sum_j (2*g_aijj - g_ajji)
!
!     batching over i and j
!
      req0 = 0
!
      req1_i = (wf%eri%n_J)*(wf%n_v)
      req1_j = (wf%eri%n_J)*(wf%n_v)
!
      req2 =  2*wf%n_o*wf%n_v
!
      batch_i = batching_index(wf%n_o)
      batch_j = batching_index(wf%n_o)
!
      call mem%batch_setup(batch_i, batch_j, req0, req1_i, req1_j, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_j_batch = 1, batch_j%num_batches
!
            call batch_j%determine_limits(current_j_batch)
!
!           F_ia = F_ia + sum_j (2*g_iajj - g_ijja)
!
            call mem%alloc(g_iajj, batch_i%length, wf%n_v, batch_j%length, batch_j%length)
            call mem%alloc(g_ijja, batch_i%length, batch_j%length, batch_j%length, wf%n_v)
!
            call wf%eri%get_eri_t1('ovoo', g_iajj, batch_i%first, batch_i%last, &
                                                   1, wf%n_v,                   &
                                                   batch_j%first, batch_j%last, &
                                                   batch_j%first, batch_j%last)
!
            call wf%eri%get_eri_t1('ooov', g_ijja, batch_i%first, batch_i%last, &
                                                   batch_j%first, batch_j%last, &
                                                   batch_j%first, batch_j%last, &
                                                   1, wf%n_v)
!
!$omp parallel do private (a, i, j)
            do a = 1, wf%n_v
               do i = 1, batch_i%length
                  do j = 1, batch_j%length
!
                     F_pq(i + batch_i%first - 1, a + wf%n_o) & 
                        = F_pq(i + batch_i%first - 1, a + wf%n_o) &
                        + two*g_iajj(i, a, j, j) - g_ijja(i, j, j, a)
!
                  enddo
               enddo 
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_ijja, batch_i%length, batch_j%length, batch_j%length, wf%n_v)
            call mem%dealloc(g_iajj, batch_i%length, wf%n_v, batch_j%length, batch_j%length)
!
!           F_ai = F_ai + sum_j (2*g_aijj - g_ajji)
!
            call mem%alloc(g_aijj, wf%n_v, batch_i%length, batch_j%length, batch_j%length)
            call mem%alloc(g_ajji, wf%n_v, batch_j%length, batch_j%length, batch_i%length)
!
            call wf%eri%get_eri_t1('vooo', g_aijj, 1, wf%n_v, &
                              batch_i%first, batch_i%last,    &
                              batch_j%first, batch_j%last,    &
                              batch_j%first, batch_j%last)
!
            call wf%eri%get_eri_t1('vooo', g_ajji, 1, wf%n_v, &
                              batch_j%first, batch_j%last,    &
                              batch_j%first, batch_j%last,    &
                              batch_i%first, batch_i%last)
!
!$omp parallel do private(i, a, j)
            do i = 1, batch_i%length
               do a = 1, wf%n_v
                  do j = 1, batch_j%length
!
                     F_pq(a + wf%n_o, i + batch_i%first - 1)   &
                     = F_pq(a + wf%n_o, i + batch_i%first - 1) &
                     + two*g_aijj(a, i, j, j) - g_ajji(a, j, j, i)
!
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_aijj, wf%n_v, batch_i%length, batch_j%length, batch_j%length)
            call mem%dealloc(g_ajji, wf%n_v, batch_j%length, batch_j%length, batch_i%length)
!
         enddo
      enddo
!
!     Add virtual-virtual contributions: F_ab = h_ab + sum_i (2*g_abii - g_aiib) 
!
!     batching over a and i
!
      req0 = 0
!
      req1_i = (wf%eri%n_J)*(wf%n_v)
      req1_a = (wf%eri%n_J)*(wf%n_v)
!
      req2 =  2*wf%n_o*wf%n_v
!
      batch_i = batching_index(wf%n_o)
      batch_a = batching_index(wf%n_v)
!
      call mem%batch_setup(batch_i, batch_a, req0, req1_i, req1_a, req2)
!
      do current_i_batch = 1, batch_i%num_batches
!
         call batch_i%determine_limits(current_i_batch)
!
         do current_a_batch = 1, batch_a%num_batches
!
            call batch_a%determine_limits(current_a_batch)
!
            call mem%alloc(g_abii, batch_a%length, wf%n_v, batch_i%length, batch_i%length)
!
            call wf%eri%get_eri_t1('vvoo', g_abii,              &
                                   batch_a%first, batch_a%last, &
                                   1, wf%n_v,                   &
                                   batch_i%first, batch_i%last, &
                                   batch_i%first, batch_i%last)
!
            call mem%alloc(g_aiib, batch_a%length, batch_i%length, batch_i%length, wf%n_v)
!
            call wf%eri%get_eri_t1('voov', g_aiib,              &
                                   batch_a%first, batch_a%last, &
                                   batch_i%first, batch_i%last, &
                                   batch_i%first, batch_i%last, &
                                   1, wf%n_v)
!
!$omp parallel do private (a, b, i)
            do a = 1, batch_a%length
               do b = 1, wf%n_v
                  do i = 1, batch_i%length
!
                     F_pq(a + batch_a%first - 1 + wf%n_o, b + wf%n_o) &
                     = F_pq(a + batch_a%first - 1 + wf%n_o, b + wf%n_o) &
                     + two*g_abii(a, b, i, i) - g_aiib(a, i, i, b)
!
                  enddo
               enddo
            enddo
!$omp end parallel do
!
            call mem%dealloc(g_aiib, batch_a%length, batch_i%length, batch_i%length, wf%n_v)
            call mem%dealloc(g_abii, batch_a%length, wf%n_v, batch_i%length, batch_i%length)
!
         enddo
      enddo
!
      call wf%set_fock(F_pq) 
      call mem%dealloc(F_pq, wf%n_mo, wf%n_mo)
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq 
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
      call copy_(wf%mo_fock_frozen, F_pq, wf%n_mo, wf%n_mo)
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
      call wf%construct_mu(mu)
!
!     Create interaction potential by scaling dipole integrals by minus electric field
!
      call mem%alloc(potential, wf%n_mo, wf%n_mo, 3)
      call zero_array(potential, wf%n_mo*wf%n_mo*3)
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
end submodule fock_ccs
