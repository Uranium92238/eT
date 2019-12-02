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
submodule (ccs_class) fock_ccs_complex
!
!!
!!    Fock submodule (CCS)
!!    Set up by Andreas Skeidsvoll, Sep 2019
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
   module subroutine construct_fock_ccs_complex(wf)
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
      implicit none
!
      class(ccs) :: wf
!
      complex(dp), dimension(:,:), allocatable :: F_pq
!
      integer :: i, j, k, a, b
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_ijkl
      complex(dp), dimension(:,:,:,:), allocatable :: g_abij
      complex(dp), dimension(:,:,:,:), allocatable :: g_aijb
      complex(dp), dimension(:,:,:,:), allocatable :: g_iajk
      complex(dp), dimension(:,:,:,:), allocatable :: g_aijk
!
!     Set F_pq = h_pq (t1-transformed) 
!
      call mem%alloc(F_pq, wf%n_mo, wf%n_mo)
!
      call wf%construct_h_complex(F_pq)
!
!     Add effective contributions to Fock matrix 
!
      if (wf%frozen_core) call wf%add_frozen_core_fock_term_complex(F_pq)
      if (wf%frozen_hf_mos) call wf%add_frozen_hf_fock_term_complex(F_pq)
      if (wf%mlhf_reference) call wf%add_mlhf_inactive_fock_term_complex(F_pq)
      if (wf%system%mm_calculation) call wf%add_molecular_mechanics_fock_term_complex(F_pq)
      if (wf%system%pcm_calculation) call wf%add_pcm_fock_contribution_complex(F_pq)
!
!     Add occupied-occupied contributions: F_ij = F_ij + sum_k (2*g_ijkk - g_ikkj)
!
      call mem%alloc(g_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_oooo_complex(g_ijkl)
!
!$omp parallel do private(i,j,k)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do k = 1, wf%n_o
!
               F_pq(i, j) = F_pq(i, j) + two_complex*g_ijkl(i,j,k,k) - g_ijkl(i,k,k,j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Add occupied-virtual contributions: F_ia = F_ia + sum_j (2*g_iajj - g_ijja)
!                                         F_ai = F_ai + sum_j (2*g_aijj - g_ajji)
!
      call mem%alloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_ovoo_complex(g_iajk)
!
      call mem%alloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_vooo_complex(g_aijk)
!
!$omp parallel do private(i,a,j)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
               F_pq(i, a + wf%n_o) = F_pq(i, a + wf%n_o) + two_complex*g_iajk(i,a,j,j) - g_iajk(j,a,i,j)
               F_pq(a + wf%n_o, i) = F_pq(a + wf%n_o, i) + two_complex*g_aijk(a,i,j,j) - g_aijk(a,j,j,i)
!
            enddo
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     Add virtual-virtual contributions: F_ab = h_ab + sum_i (2*g_abii - g_aiib) 
!
      call mem%alloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_vvoo_complex(g_abij)
!
      call mem%alloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov_complex(g_aijb)
!
!$omp parallel do private(a,b,i)
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
!
               F_pq(wf%n_o + a, wf%n_o + b) = F_pq(wf%n_o + a, wf%n_o + b) + two_complex*g_abij(a,b,i,i) - g_aijb(a,i,i,b)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%set_fock_complex(F_pq) 
      call mem%dealloc(F_pq, wf%n_mo, wf%n_mo)
!
   end subroutine construct_fock_ccs_complex
!
!
   module subroutine add_frozen_core_fock_term_ccs_complex(wf, F_pq)
!!
!!    Add frozen core Fock contribution 
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
      complex(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq 
!
      complex(dp), dimension(:,:), allocatable :: F_pq_core
!
      call mem%alloc(F_pq_core, wf%n_mo, wf%n_mo)
!
      call wf%construct_t1_fock_fc_term_complex(F_pq_core)
      call zaxpy(wf%n_mo**2, one_complex, F_pq_core, 1, F_pq, 1)
!
      call mem%dealloc(F_pq_core, wf%n_mo, wf%n_mo)      
!
   end subroutine add_frozen_core_fock_term_ccs_complex
!
!
   module subroutine add_frozen_hf_fock_term_ccs_complex(wf, F_pq)
!!
!!    Add frozen HF Fock contribution 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2019 
!!
!!    Adds the contributions from frozen HF orbitals to
!!    the effective T1-transformed Fock matrix.  
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      complex(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq 
!
      complex(dp), dimension(:,:), allocatable :: F_pq_core
!
      call mem%alloc(F_pq_core, wf%n_mo, wf%n_mo)
!
      call wf%construct_t1_fock_frozen_hf_term_complex(F_pq_core)
      call zaxpy(wf%n_mo**2, one_complex, F_pq_core, 1, F_pq, 1)
!
      call mem%dealloc(F_pq_core, wf%n_mo, wf%n_mo)      
!
   end subroutine add_frozen_hf_fock_term_ccs_complex
!
!
   module subroutine add_mlhf_inactive_fock_term_ccs_complex(wf, F_pq)
!!
!!    Add MLHF inactive Fock term
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Nov 2019 
!!
!!    Adds the contribution from MLHF inactive orbitals to
!!    the effective T1-transformed Fock matrix.  
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      complex(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: F_pq 
!
      complex(dp), dimension(:,:), allocatable :: F_pq_mlhf_inactive
!
      call mem%alloc(F_pq_mlhf_inactive, wf%n_mo, wf%n_mo)
!
      call wf%construct_t1_mlhf_inactive_fock_term_complex(F_pq_mlhf_inactive)
      call zaxpy(wf%n_mo**2, one_complex, F_pq_mlhf_inactive, 1, F_pq, 1)
!
      call mem%dealloc(F_pq_mlhf_inactive, wf%n_mo, wf%n_mo)      
!
   end subroutine add_mlhf_inactive_fock_term_ccs_complex
!
!
   module subroutine add_molecular_mechanics_fock_term_ccs_complex(wf, F_pq)
!!
!!    Add molecular mechanics Fock contribution 
!!    Written by Tommaso Giovannini, 2019 
!!
!!    Adds the molecular mechanics contributions to  
!!    the effective T1-transformed Fock matrix. 
!!
!!    Isolated into subroutine by Eirik F. Kjønstad, 2019
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: F_pq 
!
      complex(dp), dimension(:,:), allocatable :: h_mm_t1
!
      call mem%alloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
      if (wf%system%mm%forcefield == 'non-polarizable') then
!
         call wf%ao_to_t1_transformation_complex(wf%nopol_h_wx, h_mm_t1)
         call zaxpy(wf%n_mo**2, half_complex, h_mm_t1, 1, F_pq, 1)
!
      else
!
         call wf%ao_to_t1_transformation_complex(wf%pol_emb_fock, h_mm_t1)
         call zaxpy(wf%n_mo**2, half_complex, h_mm_t1, 1, F_pq, 1)
!
      endif    
!
      call mem%dealloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
   end subroutine add_molecular_mechanics_fock_term_ccs_complex
!
!
   module subroutine add_pcm_fock_contribution_ccs_complex(wf, F_pq)
!!
!!    Add PCM Fock contribution 
!!    Written by Tommaso Giovannini, 2019 
!!
!!    Adds the PCM contributions to  
!!    the effective T1-transformed Fock matrix. 
!!
!!    Isolated into subroutine by Eirik F. Kjønstad, 2019
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: F_pq 
!
      complex(dp), dimension(:,:), allocatable :: h_mm_t1
!
      call mem%alloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
!
      call wf%ao_to_t1_transformation_complex(wf%pcm_fock, h_mm_t1)
      call zaxpy(wf%n_mo**2, half_complex, h_mm_t1, 1, F_pq, 1)
!
!
      call mem%dealloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
   end subroutine add_pcm_fock_contribution_ccs_complex
!
!
   module subroutine construct_t1_fock_fc_term_ccs_complex(wf, F_pq)
!!
!!    Calculate T1 Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
      call copy_(wf%mo_fock_fc_term, F_pq, wf%n_mo, wf%n_mo)
!
      call wf%t1_transform_complex(F_pq)
!
   end subroutine construct_t1_fock_fc_term_ccs_complex
!
!
   module subroutine construct_t1_fock_frozen_hf_term_ccs_complex(wf, F_pq)
!!
!!    Calculate T1 Fock frozen fock contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
      call copy_(wf%mo_fock_frozen_hf_term, F_pq, wf%n_mo, wf%n_mo)
!
      call wf%t1_transform_complex(F_pq)
!
   end subroutine construct_t1_fock_frozen_hf_term_ccs_complex
!
!
   module subroutine construct_t1_mlhf_inactive_fock_term_ccs_complex(wf, F_pq)
!!
!!    Calculate T1 MLHF inactive Fock term
!!    Written by Sarai D. Folkestad and Linda Goletto, Nov 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
      complex(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: F_pq
!
      call copy_(wf%mlhf_inactive_fock_term, F_pq, wf%n_mo, wf%n_mo)
!
      call wf%t1_transform_complex(F_pq)
!
   end subroutine construct_t1_mlhf_inactive_fock_term_ccs_complex
!
!
   module subroutine add_t1_fock_length_dipole_term_ccs_complex(wf, electric_field)
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
      complex(dp), dimension(3), intent(in) :: electric_field
!
      complex(dp), dimension(:,:,:), allocatable :: mu, potential
!
      integer :: a, i, b, j
!
      call mem%alloc(mu, wf%n_mo, wf%n_mo, 3)
      call wf%construct_mu_complex(mu)
!
!     Create interaction potential by scaling dipole integrals by minus electric field
!
      call mem%alloc(potential, wf%n_mo, wf%n_mo, 3)
      call zero_array_complex(potential, wf%n_mo*wf%n_mo*3)
!
      call zaxpy((wf%n_mo)**2, -electric_field(1), mu(:,:,1), 1, potential(:,:,1), 1)
      call zaxpy((wf%n_mo)**2, -electric_field(2), mu(:,:,2), 1, potential(:,:,2), 1)
      call zaxpy((wf%n_mo)**2, -electric_field(3), mu(:,:,3), 1, potential(:,:,3), 1)
!     
      call mem%dealloc(mu, wf%n_mo, wf%n_mo, 3)
!
!$omp parallel do private(i,j)
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            wf%fock_ij_complex(i, j) = wf%fock_ij_complex(i, j)     &
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
            wf%fock_ia_complex(i, a) = wf%fock_ia_complex(i, a)              &
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
            wf%fock_ai_complex(a, i) = wf%fock_ai_complex(a, i)              &
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
            wf%fock_ab_complex(a, b) = wf%fock_ab_complex(a, b)                       &
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
   end subroutine add_t1_fock_length_dipole_term_ccs_complex
!
!
end submodule fock_ccs_complex
