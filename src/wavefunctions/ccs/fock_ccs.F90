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
submodule (ccs_class) fock_ccs
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
   module subroutine construct_fock_ccs(wf)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_pq = h_pq + sum_k (2*g_pqkk - g_pkkq) 
!!                   + (effective Fock contributions)
!!
!!    Since the two-electron ERIs are available already
!!    t1-transformed, our task is to transform the one-
!!    electron term, which we assume is on file in the
!!    MO basis.
!!
!!    Additional effective Fock contributions:
!!
!!    Frozen core by Sarai D. Folkestad, 2019
!!    QM/MM by Tommaso Giovannini, 2019
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_pq
      real(dp), dimension(:,:), allocatable :: F_pq_core
!
      integer :: i, j, k, a, b
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ijkl
      real(dp), dimension(:,:,:,:), allocatable :: g_abij
      real(dp), dimension(:,:,:,:), allocatable :: g_aijb
      real(dp), dimension(:,:,:,:), allocatable :: g_iajk
      real(dp), dimension(:,:,:,:), allocatable :: g_aijk
!
      real(dp), dimension(:,:), allocatable :: h_mm_t1
!
!     Get T1-transformed h integrals, put them in F_pq 
!
      call mem%alloc(F_pq, wf%n_mo, wf%n_mo)
      call wf%construct_h(F_pq)
!
!     Add QM/MM contributions
!
      if(wf%system%mm_calculation .and. wf%system%mm%forcefield.eq.'non-polarizable') then
!
         call mem%alloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
         call wf%mo_transform(half*wf%system%mm%nopol_h_wx, h_mm_t1)
         call wf%t1_transform(h_mm_t1)
!
         F_pq = F_pq + h_mm_t1 
!
         call mem%dealloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
      endif
!
      if(wf%system%mm_calculation .and. wf%system%mm%forcefield.ne.'non-polarizable') then
!
         call mem%alloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
         call wf%mo_transform(half*wf%system%mm%pol_emb_fock, h_mm_t1)
         call wf%t1_transform(h_mm_t1)
!
         F_pq = F_pq + h_mm_t1 
!
         call mem%dealloc(h_mm_t1, wf%n_mo, wf%n_mo) 
!
      endif
!
!     Add frozen core contributions
!
      if (wf%frozen_core) then
!
         call mem%alloc(F_pq_core, wf%n_mo, wf%n_mo)
!
         call wf%construct_t1_fock_fc_contribution(F_pq_core)
         call daxpy(wf%n_mo**2, one, F_pq_core, 1, F_pq, 1)
!
         call mem%dealloc(F_pq_core, wf%n_mo, wf%n_mo)
!
      endif
!
!     Occupied-occupied contributions: F_ij = F_ij + sum_k (2*g_ijkk - g_ikkj)
!
      call mem%alloc(g_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_oooo(g_ijkl)
!
!$omp parallel do private(i,j,k)
      do j = 1, wf%n_o
         do i = 1, wf%n_o
            do k = 1, wf%n_o
!
               F_pq(i, j) = F_pq(i, j) + two*g_ijkl(i,j,k,k) - g_ijkl(i,k,k,j)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Occupied-virtual contributions: F_ia = F_ia + sum_j (2*g_iajj - g_ijja)
!                                     F_ai = F_ai + sum_j (2*g_aijj - g_ajji)
!
      call mem%alloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_ovoo(g_iajk)
!
      call mem%alloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_vooo(g_aijk)
!
!$omp parallel do private(i,a,j)
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
               F_pq(i, a + wf%n_o) = F_pq(i, a + wf%n_o) + two*g_iajk(i,a,j,j) - g_iajk(j,a,i,j)
               F_pq(a + wf%n_o, i) = F_pq(a + wf%n_o, i) + two*g_aijk(a,i,j,j) - g_aijk(a,j,j,i)
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
!     Virtual-virtual contributions: F_ab = h_ab + sum_i (2*g_abii - g_aiib) ::
!
      call mem%alloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_vvoo(g_abij)
!
      call mem%alloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov(g_aijb)
!
!$omp parallel do private(a,b,i)
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
!
               F_pq(wf%n_o + a, wf%n_o + b) = F_pq(wf%n_o + a, wf%n_o + b) + two*g_abij(a,b,i,i) - g_aijb(a,i,i,b)
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%set_fock(F_pq)
      call mem%dealloc(F_pq, wf%n_mo, wf%n_mo)
!
   end subroutine construct_fock_ccs
!
!
   module subroutine coulomb_contribution_fock_fc_ccs(wf)
!!
!!    Coulomb contribution to frozen core fock
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the Coulomb contribution to Fock matrix
!!    from the frozen cores. The contribution is
!!
!!       sum_I 2 g_pqII = sum_I sum_J sum_αβ L_pq^J L_αβ^J C_αI C_βI
!!
!!    where I is core orbital which has been frozen
!!
      implicit none
!
      class(ccs) :: wf 
!
      real(dp), dimension(:,:,:), allocatable :: L_J_xy, L_J_pq
      real(dp), dimension(:,:), allocatable :: D_core_xy
      real(dp), dimension(:), allocatable :: X_J
!
      integer :: x, y, xy
      integer :: req0, req1, current_y_batch, current_q_batch
!
      type(batching_index) :: batch_y, batch_q
!
      call mem%alloc(D_core_xy, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'T',                &
                  wf%n_ao,                &
                  wf%n_ao,                &
                  wf%n_frozen_orbitals,   &
                  one,                    &
                  wf%orbital_coefficients_fc, &
                  wf%n_ao,                &
                  wf%orbital_coefficients_fc, &
                  wf%n_ao,                &
                  zero,                   &
                  D_core_xy,              &
                  wf%n_ao)
!
      call mem%alloc(X_J, wf%system%n_J)
      call zero_array(X_J, wf%system%n_J)
!
      req0 = 0
      req1 = (wf%system%n_J)*(wf%n_ao)
!
      batch_y = batching_index(wf%n_ao)
!
      call mem%batch_setup(batch_y, req0, req1)
!
!     Loop over the number of y batches
!
      do current_y_batch = 1, batch_y%num_batches
!
         call batch_y%determine_limits(current_y_batch) 
!
         call wf%system%ao_cholesky_file%open_('read')
!
         call mem%alloc(L_J_xy, wf%system%n_J, wf%n_ao, batch_y%length) 
!
         do x = 1, wf%n_ao
            do y = 1, batch_y%length
!
               xy = max(x,y+batch_y%first-1)*(max(x,y+batch_y%first-1)-3)/2 + x + y+batch_y%first-1
!
               call wf%system%ao_cholesky_file%read_(L_J_xy(:, x, y),xy)
!
            enddo
         enddo
!
         call wf%system%ao_cholesky_file%close_('keep')
!
!        X_J = L_J_xy * D_core_xy
!
         call dgemv('N', &
                     wf%system%n_J,                &
                     wf%n_ao*batch_y%length,       &
                     one,                          &
                     L_J_xy,                       &
                     wf%system%n_J,                &
                     D_core_xy(1,batch_y%first),   &
                     1,                            &
                     one,                          &
                     X_J,                          &
                     1)
!
         call mem%dealloc(L_J_xy, wf%system%n_J, wf%n_ao,  batch_y%length) 
!
      enddo ! Batch over y
!
      call mem%dealloc(D_core_xy, wf%n_ao, wf%n_ao)
!
      req0 = 0
      req1 = (wf%system%n_J)*(wf%n_mo)
!
      batch_q = batching_index(wf%n_mo)
!
      call mem%batch_setup(batch_q, req0, req1)
!
!     Loop over the number of q batches

      do current_q_batch = 1, batch_q%num_batches
!
         call batch_q%determine_limits(current_q_batch) 
!
         call mem%alloc(L_J_pq, wf%system%n_J, wf%n_mo, batch_q%length)
!
         call wf%integrals%read_cholesky(L_J_pq, 1, wf%n_mo, batch_q%first, batch_q%last)
!
         call dgemv('T',                                             &
                     wf%system%n_J,                                  &
                     wf%n_mo*batch_q%length,                         &
                     two,                                            &
                     L_J_pq,                                         &
                     wf%system%n_J,                                  &
                     X_J,                                            &
                     1,                                              &
                     one,                                            &
                     wf%mo_fock_fc_contribution(1,batch_q%first),    &
                     1)
!
         call mem%dealloc(L_J_pq, wf%system%n_J, wf%n_mo, batch_q%length)
!
      enddo ! Batch over q
!
      call mem%dealloc(X_J, wf%system%n_J)
!
   end subroutine coulomb_contribution_fock_fc_ccs
!
!
   module subroutine exchange_contribution_fock_fc_ccs(wf)
!!
!!    Exchange contribution to frozen core Fock
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Calculates the exchange contribution to Fock matrix
!!    from the frozen cores. The contribution is
!! 
!!       - sum_I g_pIIq = sum_I sum_J sum_αβγδ L_γδ^J L_αβ^J C_γI C_βI C_αp C_δq 
!!
!!       X_Iδ^J = L_γδ^J C_γI
!!
!!       Y_Iq^J = X_Iδ^J C_δq
!!
!!       F_pq -= sum_I sum_J Y_Iq^J Y_Ip^J
!!
!!    where I is core orbital which has been frozen
!!
      implicit none
!
      class(ccs) :: wf 
!
      real(dp), dimension(:,:,:), allocatable :: L_J_xy, X_J_xI, X_J_Ix, Y_J_Ip
!
      integer :: x, y, xy      
      integer :: req0, req1, current_y_batch
!
      type(batching_index) :: batch_y
!
      call mem%alloc(X_J_xI, wf%system%n_J, wf%n_ao, wf%n_frozen_orbitals)
      call zero_array(X_J_xI, (wf%system%n_J)*(wf%n_ao)*(wf%n_frozen_orbitals))
!
      req0 = 0
      req1 = (wf%system%n_J)*(wf%n_ao)
!
      batch_y = batching_index(wf%n_ao)
!
      call mem%batch_setup(batch_y, req0, req1)
!
!     Loop over the number of y batches
!
      do current_y_batch = 1, batch_y%num_batches
!
         call batch_y%determine_limits(current_y_batch) 
!
         call wf%system%ao_cholesky_file%open_('read')
!
         call mem%alloc(L_J_xy, wf%system%n_J, wf%n_ao, batch_y%length) 
!
         do x = 1, wf%n_ao
            do y = 1, batch_y%length
!
               xy = max(x,y + batch_y%first - 1)*(max(x,y + batch_y%first - 1)-3)/2 + x + y + batch_y%first - 1
!
               call wf%system%ao_cholesky_file%read_(L_J_xy(:, x, y),xy)
!
            enddo
         enddo
!
         call wf%system%ao_cholesky_file%close_('keep')
!
!        X_Iδ^J = L_γδ^J C_γI
!
         call dgemm('N', 'N',                               &
                     (wf%system%n_J)*(wf%n_ao),             &
                     wf%n_frozen_orbitals,                  &
                     batch_y%length,                        &
                     one,                                   &
                     L_J_xy,                                &
                     (wf%system%n_J)*(wf%n_ao),             &
                     wf%orbital_coefficients_fc(batch_y%first,1),   &
                     wf%n_ao,                               &
                     one,                                   &
                     X_J_xI,                                &
                     (wf%system%n_J)*(wf%n_ao))
!
         call mem%dealloc(L_J_xy, wf%system%n_J, wf%n_ao, batch_y%length) 
!
      enddo ! Batches of y
!
      call mem%alloc(X_J_Ix, wf%system%n_J, wf%n_frozen_orbitals, wf%n_ao)
!
      call sort_123_to_132(X_J_xI, X_J_Ix, wf%system%n_J, wf%n_ao, wf%n_frozen_orbitals)     
!
      call mem%dealloc(X_J_xI, wf%system%n_J, wf%n_ao, wf%n_frozen_orbitals)
!
!     Y_Iq^J = X_Iδ^J C_δq
!
      call mem%alloc(Y_J_Ip, wf%system%n_J, wf%n_frozen_orbitals, wf%n_mo)
!
      call dgemm('N', 'N',                                  &
                  wf%system%n_J*(wf%n_frozen_orbitals),     &
                  wf%n_mo,                                  &
                  wf%n_ao,                                  &
                  one,                                      &
                  X_J_Ix,                                   &
                  wf%system%n_J*(wf%n_frozen_orbitals),     &
                  wf%orbital_coefficients,                  &
                  wf%n_ao,                                  &
                  zero,                                     &
                  Y_J_Ip,                                   &
                  wf%system%n_J*(wf%n_frozen_orbitals))
!
      call mem%dealloc(X_J_Ix, wf%system%n_J, wf%n_frozen_orbitals, wf%n_ao)
!
!     F_pq -=  Y_J_Ip * Y_J_Iq
!
      call dgemm('T', 'N',                                  &
                  wf%n_mo,                                  &
                  wf%n_mo,                                  &
                  wf%system%n_J*(wf%n_frozen_orbitals),     &
                  -one,                                     &
                  Y_J_Ip,                                   &
                  wf%system%n_J*(wf%n_frozen_orbitals),     &
                  Y_J_Ip,                                   &
                  wf%system%n_J*(wf%n_frozen_orbitals),     &
                  one,                                      &
                  wf%mo_fock_fc_contribution,               &
                  wf%n_mo)
!
      call mem%dealloc(Y_J_Ip, wf%system%n_J, wf%n_frozen_orbitals, wf%n_mo)
!
   end subroutine exchange_contribution_fock_fc_ccs
!
!
   module subroutine construct_mo_fock_fc_contribution_ccs(wf)
!!
!!    Calculate MO Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
      implicit none
!
      class(ccs) :: wf 
!
      
      call output%printf(':: Frozen core approximation is used', &
                        ints=[wf%n_frozen_orbitals], fs='(/t3,a)',pl='minimal')
!
      call output%printf('There are (i0) frozen orbitals.', &
                        ints=[wf%n_frozen_orbitals],pl='minimal',fs='(t6,a)')
!
      call wf%initialize_mo_fock_fc_contribution()
!
      call zero_array(wf%mo_fock_fc_contribution, wf%n_mo**2)
!
      call wf%coulomb_contribution_fock_fc()
      call wf%exchange_contribution_fock_fc()
!
      call wf%destruct_orbital_coefficients_fc()
!
   end subroutine construct_mo_fock_fc_contribution_ccs
!
!
   module subroutine construct_t1_fock_fc_contribution_ccs(wf, F_pq)
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
      call dcopy(wf%n_mo**2, wf%mo_fock_fc_contribution, 1, F_pq, 1)
!
      call wf%t1_transform(F_pq)
!
   end subroutine construct_t1_fock_fc_contribution_ccs
!
!
end submodule fock_ccs
