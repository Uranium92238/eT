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
submodule (mlcc2_class) orbitals_mlcc2
!
!!
!!    MLCC2 orbitals
!!    Written by Sarai D. Folkestad, Apr 2019
!!
!!    This submodule contains routines that handle orbital
!!    transformation and orbital partitioning for MLCC2.
!!
!!    Orbital construction in MLCC calculations consists of 3
!!    steps:
!!
!!    1. Construct the orbital basis specified on input
!!       (e.g., Cholesky, Cholesky-PAO, CNTOs, NTO/Canonical)
!!
!!    2. If the number of active/inactive orbitals not determined
!!       on input, it is determined during orbital construction.
!!       Orbital coefficient matrix is ordered according to 
!!       active and inactive:
!!       
!!          C = [active occupied, inactive occupied, active virtual, inactive virtual]
!!
!!    3. The Fock matrix is constructed, and the F_oo and F_vv blocks are block 
!!       diagonalized, such that the active-active and inactive-inactive
!!       blocks are diagonal. The corresponding basis is the basis used in the 
!!       MLCC calculation.
!!
!!    Note that step 1 and 2 happens at the same time, in the routines 
!!    which constructs the orbitals
!!
!!
!!
!
   implicit none
!
contains
!
!
   module subroutine orbital_partitioning_mlcc2(wf)
!!
!!    Orbital partitioning 
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    This routine drives the orbital 
!!    partitioning for MLCC2 
!!
!!    Based on which orbitals are requested on 
!!    input, we construct this new orbital basis 
!!    and sets, n_cc2_o, n_cc2_v, n_ccs_o and n_ccs_v,
!!    first_cc2_o, last_cc2_o, first_cc2_v, last_cc2_v
!!    first_ccs_o, last_ccs_o, first_ccs_v, last_ccs_v
!!
!!    NOTE: This routine shold always be 
!!    followed by a routine which block diagonalizes
!!    the Fock matrix!
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: T_o, T_v
!
      type(timings) :: timer
!
      logical :: cholesky_o_only
!
      timer = timings('MLCC orbital construction')
      call timer%turn_on()
!
      if (trim(wf%cc2_orbital_type) == 'cholesky') then
!
         call wf%construct_cholesky_orbitals()
!
      elseif (trim(wf%cc2_orbital_type) == 'cnto-approx') then
!
         wf%n_ccs_o = wf%n_o - wf%n_cc2_o
         wf%n_ccs_v = wf%n_v - wf%n_cc2_v
!
         call mem%alloc(T_o, wf%n_o, wf%n_o)
         call mem%alloc(T_v, wf%n_v, wf%n_v)
!
         if (wf%cnto_restart) then
!
            call wf%read_cnto_transformation_matrices(T_o, T_v)
!
         else
!
            call wf%construct_ccs_cnto_transformation_matrices(T_o, T_v)
!
         endif
!
         call wf%construct_cntos(T_o, T_v)
!
         call mem%dealloc(T_o, wf%n_o, wf%n_o)
         call mem%dealloc(T_v, wf%n_v, wf%n_v)
!
      elseif (trim(wf%cc2_orbital_type) == 'nto-canonical') then
!
         call mem%alloc(T_o, wf%n_o, wf%n_o)
!
         call wf%construct_ccs_nto_transformation_matrix(T_o)
!
         call wf%construct_mixed_nto_canonical_orbitals(T_o)
!
         call mem%dealloc(T_o, wf%n_o, wf%n_o)
!
      elseif (trim(wf%cc2_orbital_type) == 'cholesky-pao') then
!
         cholesky_o_only = .true.
         call wf%construct_cholesky_orbitals(cholesky_o_only)
!
         wf%n_cc2_v = 0
         wf%n_ccs_v = 0
!
         call wf%construct_paos()
!
      else
!
         call output%error_msg('could not recognize the selected cc2 orbitals')
!
      endif
!
!     Set orbital partitioning specifications
!
      wf%first_cc2_o = 1
      wf%first_cc2_v = 1
!
      wf%last_cc2_o = wf%n_cc2_o
      wf%last_cc2_v = wf%n_cc2_v
!
      wf%first_ccs_o = wf%last_cc2_o + 1
      wf%first_ccs_v = wf%last_cc2_v + 1
!
      wf%last_ccs_o = wf%n_o
      wf%last_ccs_v = wf%n_v
!
      call timer%turn_off()
!
   end subroutine orbital_partitioning_mlcc2
!
!
   module subroutine construct_cholesky_orbitals_mlcc2(wf, occupied_only)
!!
!!    Construct Cholesky orbitals
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Constructs Cholesky orbitals by
!!    decomposing the HF AO density and the
!!    HF virtual AO density
!!
!!    See A. M. J. Sánchez de Merás, H. Koch, 
!!    I. G. Cuesta, and L. Boman (J. Chem. Phys. 132, 204105 (2010))
!!    for more information on active space generation
!!    using Cholesky decomposition
!!
!!    'occupied_only' : Optional argument that is used to 
!!                      determine if we construct occuped 
!!                      Cholesky orbitals only, or if we also 
!!                      construct the virtual Cholesky orbitals
!!                      DEFAULT: .false.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      logical, intent(in), optional :: occupied_only
!
      logical  :: occupied_only_local
!
      real(dp), dimension(:,:), allocatable :: D
!
      integer, dimension(:), allocatable :: active_aos
!
      integer :: first_ao, last_ao, i, n_active_aos
!
      real(dp), parameter :: full_cd_threshold = 1.0d-4
!
      integer :: mo_offset
!
      occupied_only_local = .false. 
!
      if (present(occupied_only)) occupied_only_local = occupied_only
!
!     Construct active occupied orbitals     
!
!     0. Determine active ao list
!
      call wf%system%first_and_last_ao_active_space('cc2', first_ao, last_ao)
!
      n_active_aos = last_ao - first_ao + 1
!
      call mem%alloc(active_aos, n_active_aos)
!
      do i = 1, n_active_aos
!
         active_aos(i) = first_ao + i - 1
!
      enddo
!
!     1. Set up active occupied density 
!
      call mem%alloc(D, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'T',                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  D,                         &
                  wf%n_ao)
!
!     Construct active occupied orbitals
!
      mo_offset = 0
!
      call wf%construct_orbital_block_by_density_cd(D, wf%n_cc2_o, wf%cholesky_orbital_threshold, mo_offset, active_aos)
!
!     Construct inactive occupied orbitals
!
      mo_offset = wf%n_cc2_o
!
      call wf%construct_orbital_block_by_density_cd(D, wf%n_ccs_o, full_cd_threshold, mo_offset)
!
!     Construct active virtual orbitals     
!
      if (.not. occupied_only_local) then
!
!        1. Set up virtual density     
!
         call dgemm('N', 'T',                                  &
                     wf%n_ao,                                  &
                     wf%n_ao,                                  &
                     wf%n_v,                                   &
                     one,                                      &
                     wf%orbital_coefficients(1, wf%n_o + 1),   &
                     wf%n_ao,                                  &
                     wf%orbital_coefficients(1, wf%n_o + 1),   &
                     wf%n_ao,                                  &
                     zero,                                     &
                     D,                                        &
                     wf%n_ao)
!
!        Construct active virtual orbitals
!
         mo_offset = wf%n_o
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_cc2_v, wf%cholesky_orbital_threshold, mo_offset, active_aos)
!
!        Construct inactive virtual orbitals
!
         mo_offset = wf%n_o + wf%n_cc2_v
!
         call wf%construct_orbital_block_by_density_cd(D, wf%n_ccs_v, full_cd_threshold, mo_offset)
!
      endif
!
      call mem%dealloc(active_aos, n_active_aos)
      call mem%dealloc(D, wf%n_ao, wf%n_ao)
!
   end subroutine construct_cholesky_orbitals_mlcc2
!
!
   module subroutine construct_block_diagonal_fock_orbitals_mlcc2(wf)
!!
!!    Construct block diagonal Fock MOs 
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    This routine block-diagonalizes the occupied-occupied
!!    and virutal-virtual blocks of the Fock matrix s.t.
!!    the active-active, and inactive-inactive blocks are 
!!    diagonal.
!!   
!!    Note that after this routine, the Fock matrix in wf 
!!    corresponds to the old basis.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_oo, F_vv
!
      call mem%alloc(F_oo, wf%n_o, wf%n_o)
      call mem%alloc(F_vv, wf%n_v, wf%n_v)
!
      call dcopy(wf%n_o**2, wf%fock_ij, 1, F_oo, 1)
      call dcopy(wf%n_v**2, wf%fock_ab, 1, F_vv, 1)
!
!     Block diagonal occupied-occupied Fock
!
      call wf%construct_block_diagonal_fock_mos_2_level( wf%n_o, wf%n_cc2_o, F_oo, &
                           wf%orbital_coefficients(:,1:wf%n_o), wf%orbital_energies(1:wf%n_o))

!
!     Block diagonal virtual-virtual Fock
!
      call wf%construct_block_diagonal_fock_mos_2_level( wf%n_v, wf%n_cc2_v, F_vv, &
                           wf%orbital_coefficients(:,wf%n_o + 1 : wf%n_mo), wf%orbital_energies(wf%n_o + 1 : wf%n_mo))
!
      call mem%dealloc(F_oo, wf%n_o, wf%n_o)
      call mem%dealloc(F_vv, wf%n_v, wf%n_v)
!
   end subroutine construct_block_diagonal_fock_orbitals_mlcc2
!
!
   module subroutine construct_block_diagonal_fock_mos_2_level_mlcc2(wf, n_total, n_active, fock, MO_coeff, diagonal)
!!
!!    Construct Fock block diagonal 2 levels
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Construct orbitals that block diagonalizes 
!!    Fock matrix block for two levels (inactive/active)
!!
!!    'n_total' : Total dimmension of Fock matrix block
!!
!!    'n_active' : Dimension of active block of Fock matrix block
!!
!!    'fock' : Fock matrix block to block diagonalize
!!
!!    'MO_coef' : MO coefficients which are updated to 
!!                the new basis which block diagonalizes Fock
!!                matrix block
!!
!!    'diagonal' : Diagonal elements of Fock matrix after block
!!                 diagonalization
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      integer, intent(in) :: n_total, n_active ! Total matrix dimension of block, and number of active orbitals
!
      real(dp), dimension(n_total, n_total), intent(inout) :: fock
      real(dp), dimension(wf%n_ao, n_total), intent(inout) :: mo_coeff
      real(dp), dimension(n_total), intent(out) :: diagonal
!
      real(dp), dimension(:), allocatable :: work
      real(dp), dimension(:), allocatable :: orbital_energies
      real(dp), dimension(:,:), allocatable :: C_active, C_inactive
!
      integer :: info, i, j
!
!     1. Active block     
!
      if (n_active .gt. 0) then
!
!        Diagonalize active block
!
         call mem%alloc(work, 4*n_active)
         call mem%alloc(orbital_energies, n_active)
!
         call dsyev('V','U',           &
                     n_active,         &
                     fock,             &
                     n_total,          &
                     orbital_energies, &
                     work,             &
                     4*n_active,       &
                     info)
!
         call mem%dealloc(work, 4*n_active)
!
         if (info .ne. 0) call output%error_msg('Diagonalization of active fock matrix block failed')
!
!        Setting diagonal (orbital energies)
!
         do i = 1, n_active
!
            diagonal(i) = orbital_energies(i)
!
         enddo
!
         call mem%dealloc(orbital_energies, n_active)
!
      endif
!
!     2. Inactive block
!
      if ((n_total - n_active) .gt. 0) then
!
!        Diagonalize inactive block
!
         call mem%alloc(work, 4*(n_total - n_active))
         call mem%alloc(orbital_energies, (n_total - n_active))
!
         call dsyev('V','U',                                &
                     (n_total - n_active),                  &
                     fock(n_active + 1, n_active + 1),      &
                     n_total,                               &
                     orbital_energies,                      &
                     work,                                  &
                     4*(n_total - n_active),                &
                     info)
!
         call mem%dealloc(work, 4*(n_total - n_active))
!
         if (info .ne. 0) call output%error_msg('Diagonalization of inactive fock matrix block failed')
!
!        Setting diagonal (orbital energies)
!
         do i = 1, (n_total - n_active)
!
            diagonal(n_active + i) = orbital_energies(i)
!
         enddo
!
         call mem%dealloc(orbital_energies, (n_total - n_active))
!
      endif
!
!     Transform orbital coefficients
!
!     1. Active 
!
      if (n_active .gt. 0) then
!
         call mem%alloc(C_active, wf%n_ao, n_active)
!
         call dgemm('N', 'N',    &
                     wf%n_ao,    &
                     n_active,   &
                     n_active,   &
                     one,        &
                     mo_coeff,   &
                     wf%n_ao,    &
                     fock,       &
                     n_total,    &
                     zero,       &
                     C_active,   &
                     wf%n_ao)
!
!$omp parallel do private(i, j)
         do j = 1, n_active
            do i = 1, wf%n_ao
!
              mo_coeff(i, j) = C_active(i, j)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(C_active, wf%n_ao, n_active)
!
      endif
!
!     1. Inactive 
!
      if ((n_total - n_active) .gt. 0) then
!
         call mem%alloc(C_inactive, wf%n_ao, (n_total - n_active))
!
         call dgemm('N', 'N',                            &
                     wf%n_ao,                            &
                     (n_total - n_active),               &
                     (n_total - n_active),               &
                     one,                                &
                     mo_coeff(1, n_active + 1),          &
                     wf%n_ao,                            &
                     fock(n_active + 1, n_active + 1),   &
                     n_total,                            &
                     zero,                               &
                     C_inactive,                         &
                     wf%n_ao)
!
!$omp parallel do private(i, j)
         do j = 1, (n_total - n_active)
            do i = 1, wf%n_ao
!
               mo_coeff(i, j + n_active) = C_inactive(i, j)
!
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(C_inactive, wf%n_ao, (n_total - n_active))
!
      endif    
!
   end subroutine construct_block_diagonal_fock_mos_2_level_mlcc2
!
!
   module subroutine construct_M_and_N_cnto_mlcc2(wf, R_ai, R_aibj, M, N, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij += ( sum_a R_ai R_aj + 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibj R_ajbl )
!!       N_ab += ( sum_i R_ai R_bi + 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    Used to construct CNTOs.
!!
!!    set_to_zero determines if M and N are set to zero initially. This makes it possible for 
!!    the routine ton be used to add to M and N if we are using more than one excitation vector 
!!    to generate the CNTOs.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o, wf%n_v, wf%n_o), intent(in)  :: R_aibj
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)                  :: R_ai
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
!
      logical, intent(in) :: set_to_zero
!
      integer :: a, i
!
      if (set_to_zero) then
!
         call zero_array(M, wf%n_o**2)
         call zero_array(N, wf%n_v**2)
!
      endif
!
!     1. Singles contribution
!
!     M_ij += sum_a R_ai R_aj
!
      call dgemm('T', 'N', &
                  wf%n_o,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  one,     &
                  R_ai,    & ! R_ai
                  wf%n_v,  &
                  R_ai,    & ! R_aj
                  wf%n_v,  &
                  one,     &
                  M,       & ! M_ij
                  wf%n_o)
!
!     N_ab += sum_i R_ai R_aj
!
      call dgemm('N', 'T', &
                  wf%n_v,  &
                  wf%n_v,  &
                  wf%n_o,  &
                  one,     &
                  R_ai,    & ! R_ai
                  wf%n_v,  &
                  R_ai,    & ! R_bi
                  wf%n_v,  &
                  one,     &
                  N,       & ! N_ab
                  wf%n_v)
!
!     2. Doubles contribution
!
!     M_ij += 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R^k_aibj R^k_blaj )
!
      call dgemm('T', 'N',             &
                 wf%n_o,               &
                 wf%n_o,               &
                 (wf%n_v**2)*(wf%n_o), &
                 half,                 &
                 R_aibj,               & ! R_bla_i
                 (wf%n_v**2)*(wf%n_o), &
                 R_aibj,               & ! R_bla_j
                 (wf%n_v**2)*(wf%n_o), &
                 one,                  &
                 M,                    &
                 wf%n_o)
!
!$omp parallel do private(a, i)
   do i = 1, wf%n_o
      do a = 1, wf%n_v
!
            M(i, i) = M(i, i) + half*R_aibj(a,i,a,i)**2
!
      enddo
   enddo
!$omp end parallel do
!
!     N_ab += 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R^k_aicj R^k_bicj)
!
      call dgemm('N', 'T',             &
                 wf%n_v,               &
                 wf%n_v,               &
                 (wf%n_v)*(wf%n_o**2), &
                 half,                 &
                 R_aibj,               &  ! R_a_icj
                 (wf%n_v),             &
                 R_aibj,               &  ! R_b_icj
                 (wf%n_v),             &
                 one,                  &
                 N,                    &
                 wf%n_v)
!
!$omp parallel do private(a, i)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            N(a, a) = N(a, a) + half*R_aibj(a,i,a,i)**2
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_M_and_N_cnto_mlcc2
!
!
   module subroutine construct_cntos_mlcc2(wf, T_o, T_v)
!!
!!    Construct correlated natural transition orbitals
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Transform orbital coefficients to CNTO
!!    basis.
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: T_v
!
      real(dp), dimension(:,:), allocatable     :: MO_coeff
!
      call mem%alloc(MO_coeff, wf%n_ao, wf%n_mo)
      call dcopy((wf%n_ao)*(wf%n_mo), wf%orbital_coefficients, 1, MO_coeff, 1)
!
      call zero_array(wf%orbital_coefficients, (wf%n_ao)*(wf%n_mo))
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_o,                    &
                  wf%n_o,                    &
                  one,                       &
                  MO_coeff,                  &
                  wf%n_ao,                   &
                  T_o,                       &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao)
!
      call dgemm('N', 'N',                                  &
                  wf%n_ao,                                  &
                  wf%n_v,                                   &
                  wf%n_v,                                   &
                  one,                                      &
                  MO_coeff(1, wf%n_o + 1),                  &
                  wf%n_ao,                                  &
                  T_v,                                      &
                  wf%n_v,                                   &
                  one,                                      &
                  wf%orbital_coefficients(1, wf%n_o + 1),   &
                  wf%n_ao)
!
      call mem%dealloc(MO_coeff, wf%n_ao, wf%n_mo)
!
   end subroutine construct_cntos_mlcc2
!
!
   module subroutine read_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Read CNTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Read CNTO transformation matrices.
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
      type(sequential_file) :: transformation_o, transformation_v
!
      transformation_o = sequential_file('cnto_M_transformation', 'unformatted')
      transformation_v = sequential_file('cnto_N_transformation', 'unformatted')
!
      call transformation_o%open_('read', 'rewind')
      call transformation_v%open_('read', 'rewind')
!
      call transformation_o%read_(T_o, wf%n_o**2)
      call transformation_v%read_(T_v, wf%n_v**2)
!
      call transformation_o%close_('keep')
      call transformation_v%close_('keep')
!
   end subroutine read_cnto_transformation_matrices_mlcc2
!
!
   module subroutine construct_ccs_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Construct CCS CNTO transformation matrices
!!    Written by Sarai D. Folkestad, May 2019
!!
!!       - Run CCS calculation
!!
!!       - Construct M and N for CNTOs 
!!
!!       - Diagonalize M and N
!!
!!       - Write transformation matrices to file
!!
!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(out) :: T_v
!
      integer :: n_cnto_states, k
!
      real(dp), dimension(:), allocatable       :: omega_ccs ! CCS excitation energies
      real(dp), dimension(:,:), allocatable     :: R_ai_k
      real(dp), dimension(:,:,:), allocatable   :: R_ai
      real(dp), dimension(:,:,:,:), allocatable :: R_aibj_k
!
      logical :: set_to_zero
!
      type(sequential_file) :: transformation_o, transformation_v
!
      character(len=200) :: r_or_l
!
      n_cnto_states = size(wf%cnto_states)
!
      r_or_l = 'right'
!
      if (input%requested_keyword_in_section('left eigenvectors', 'solver cc es')) r_or_l = 'left'

      call mem%alloc(R_ai, wf%n_v, wf%n_o, n_cnto_states)
      call mem%alloc(omega_ccs, n_cnto_states)
!
!     Run CCS calculation 
!
      call wf%ccs_calculation_for_cntos(r_or_l, n_cnto_states, R_ai, wf%cnto_states, omega_ccs)
!
      set_to_zero = .true.
!
      call mem%alloc(R_ai_k, wf%n_v, wf%n_o)
      call mem%alloc(R_aibj_k, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!   
      do k = 1, n_cnto_states
!
         call dcopy(wf%n_t1, R_ai(1,1,k), 1, R_ai_k, 1)
!
!        Construct approximate double excitation vector
! 
         call wf%approximate_double_excitation_vectors(R_ai_k, R_aibj_k, omega_ccs(k))
!
!        Add contribution to M and N
!
         call wf%construct_M_and_N_cnto(R_ai_k, R_aibj_k, T_o, T_v, set_to_zero)
!
         set_to_zero = .false.
!
      enddo
!
      call mem%dealloc(R_ai_k, wf%n_v, wf%n_o)
      call mem%dealloc(R_aibj_k, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call mem%dealloc(R_ai, wf%n_v, wf%n_o, n_cnto_states)
!
      call wf%diagonalize_M_and_N(T_o, T_v)
!
!     Write eigenvectors of M and N 
!
      transformation_o = sequential_file('cnto_M_transformation', 'unformatted')
      transformation_v = sequential_file('cnto_N_transformation', 'unformatted')
!
      call transformation_o%open_('write', 'rewind')
      call transformation_v%open_('write', 'rewind')
!
      call transformation_o%write_(T_o, wf%n_o**2)
      call transformation_v%write_(T_v, wf%n_v**2)
!
      call transformation_o%close_('keep')
      call transformation_v%close_('keep')
!
   end subroutine construct_ccs_cnto_transformation_matrices_mlcc2
!
!
   module subroutine construct_M_nto_mlcc2(wf, R_ai, M, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij +=  sum_a R_ai R_aj
!!
!!    Used to construct occupied NTOs.
!!
!!    set_to_zero determines if M and N are set to zero initially. This is used such that 
!!    the routine can be used to add to M and N if we are using more than one excitation vector 
!!    to generate the CNTOs.
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
!
      logical, intent(in) :: set_to_zero
!
      if (set_to_zero) then
!
         call zero_array(M, wf%n_o**2)
!
      endif
!
!     1. Singles contribution
!
!     M_ij += sum_a R_ai R_aj
!
      call dgemm('T', 'N', &
                  wf%n_o,  &
                  wf%n_o,  &
                  wf%n_v,  &
                  one,     &
                  R_ai,    & ! R_ai
                  wf%n_v,  &
                  R_ai,    & ! R_aj
                  wf%n_v,  &
                  one,     &
                  M,       & ! M_ij
                  wf%n_o)
!
   end subroutine construct_M_nto_mlcc2
!
!
   module subroutine construct_ccs_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Construct CCS CNTO transformation matrices
!!    Written by Sarai D. Folkestad, May 2019
!!
!!       - Run CCS calculation
!!
!!       - Construct M NTOs 
!!
!!       - Diagonalize M
!!
!!       - Write transformation matrix to file
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
!
      integer :: k, info, n_nto_states
!
      real(dp), dimension(:), allocatable       :: work, eigenvalues
      real(dp), dimension(:,:), allocatable     :: R_ai_k
      real(dp), dimension(:,:,:), allocatable   :: R_ai
!
      logical :: set_to_zero
!
      type(sequential_file) :: transformation_o

      character(len=200) :: r_or_l
!
      n_nto_states = size(wf%nto_states)
!
      r_or_l = 'right'
!
      if (input%requested_keyword_in_section('left eigenvectors', 'solver cc es')) r_or_l = 'left'
!
      call mem%alloc(R_ai, wf%n_v, wf%n_o, n_nto_states)
!
!     Run CCS calculation
!
      call wf%ccs_calculation_for_cntos(r_or_l, n_nto_states, R_ai, wf%nto_states)
!
      set_to_zero = .true.
!
      call mem%alloc(R_ai_k, wf%n_v, wf%n_o)
!   
      do k = 1, n_nto_states
!
         call dcopy(wf%n_t1, R_ai(1,1,k), 1, R_ai_k, 1)
!
!        Add contribution to M 
!
         call wf%construct_M_nto(R_ai_k, T_o, set_to_zero)
!
         set_to_zero = .false.
!
      enddo
!
      call mem%dealloc(R_ai_k, wf%n_v, wf%n_o)
      call mem%dealloc(R_ai, wf%n_v, wf%n_o, n_nto_states)
!
!     Diagonalize -M  to get eigenvectors in correct order
!
      call dscal(wf%n_o**2, -one, T_o, 1)
!
      call mem%alloc(work, 4*wf%n_o)
      call mem%alloc(eigenvalues, wf%n_o)
!
      call dsyev('V','U',        &
                  wf%n_o,        &
                  T_o,           &
                  wf%n_o,        &
                  eigenvalues,   &
                  work,          &
                  4*wf%n_o,      &
                  info)
!
      call mem%dealloc(work, 4*wf%n_o)
      call mem%dealloc(eigenvalues, wf%n_o)
!
      if (info .ne. 0) call output%error_msg('Diagonalization of M failed')
!
!     Write eigenvectors of M 
!
      transformation_o = sequential_file('cnto_M_transformation', 'unformatted')
!
      call transformation_o%open_('write', 'rewind')
!
      call transformation_o%write_(T_o, wf%n_o**2)
!
      call transformation_o%close_('keep')
!
   end subroutine construct_ccs_nto_transformation_matrix_mlcc2
!
!
   module subroutine construct_mixed_nto_canonical_orbitals_mlcc2(wf, T_o)
!!
!!    Construct mixed NTO and canonical orbitals
!!    Written by Sarai D. Folekstad, Jun 2019
!!
!!    Construct occupiued NTOs, leave canonical virtuals
!!
      implicit none 
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
!
      real(dp), dimension(:,:), allocatable     :: MO_coeff
!
      call mem%alloc(MO_coeff, wf%n_ao, wf%n_mo)
      call dcopy((wf%n_ao)*(wf%n_mo), wf%orbital_coefficients, 1, MO_coeff, 1)
!
      call zero_array(wf%orbital_coefficients, (wf%n_ao)*(wf%n_o))
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_o,                    &
                  wf%n_o,                    &
                  one,                       &
                  MO_coeff,                  &
                  wf%n_ao,                   &
                  T_o,                       &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao)
!
      wf%n_ccs_o = wf%n_o - wf%n_cc2_o
      wf%n_ccs_v = wf%n_v - wf%n_cc2_v
!
      call mem%dealloc(MO_coeff, wf%n_ao, wf%n_mo)
!
   end subroutine construct_mixed_nto_canonical_orbitals_mlcc2
!
!
   module subroutine construct_paos_mlcc2(wf)
!!
!!    Construct PAOs
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Construct projected atomic orbitals for
!!    virtual orbitals.
!!
!!    1. Construct PAOs on active atoms
!!
!!    2. Orthonormalize these active virtual orbitals 
!!
!!    3. PAOs for the remaining system by projecting out
!!       both occupied and active virtuals out of all  AOs
!!
!!    4. Orthonormalize these inactive virtual orbitals 
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: D, S, PAO_coeff
!
      integer :: first_ao, last_ao, ao, mo, n_active_aos, rank
!
!     0. Determine active ao list
!
      call wf%system%first_and_last_ao_active_space('cc2', first_ao, last_ao)
!
      n_active_aos = last_ao - first_ao + 1
!
!     1. Set up occupied density 
!
      call mem%alloc(D, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'T',                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  D,                         &
                  wf%n_ao)
!
!     2. Construct PAOs for active atoms
!
      call mem%alloc(PAO_coeff, wf%n_ao, n_active_aos)
!
      call wf%projected_atomic_orbitals(D, PAO_coeff, n_active_aos, first_ao)
!
!     3. Orthonormalize PAOs to get active virtual orbitals
!
      call mem%alloc(S, n_active_aos, n_active_aos)
!
      call wf%get_orbital_overlap(PAO_coeff, n_active_aos, S)
!
      call wf%lovdin_orthonormalization(PAO_coeff, S, n_active_aos, rank)
!
      call mem%dealloc(S, n_active_aos, n_active_aos)
!
      wf%n_cc2_v = rank 
!
!     Set the active virtual orbital coefficients
!
      do mo = wf%n_o + 1, wf%n_o + wf%n_cc2_v
         do ao = 1, wf%n_ao
!
            wf%orbital_coefficients(ao, mo) = PAO_coeff(ao, mo - wf%n_o)
!
         enddo
      enddo
!
      call mem%dealloc(PAO_coeff, wf%n_ao, n_active_aos)
!
      if (rank .lt. wf%n_v) then 
!
!        Set virtual inactive
!
!        4. Construct M = sum_p C_αp C_βp  for p = 1, n_o + n_cc2_v
!
         call dgemm('N', 'T',                   &
                     wf%n_ao,                   &
                     wf%n_ao,                   &
                     wf%n_o + wf%n_cc2_v,       &
                     one,                       &
                     wf%orbital_coefficients,   &
                     wf%n_ao,                   &
                     wf%orbital_coefficients,   &
                     wf%n_ao,                   &
                     zero,                      &
                     D,                         &
                     wf%n_ao)
!
!        Construct PAOs for the remaining virtual orbitals
!
         call mem%alloc(PAO_coeff, wf%n_ao, wf%n_ao)
!
         call wf%projected_atomic_orbitals(D, PAO_coeff, wf%n_ao)
!
!        5. Orthonormalize PAOs to get inactive virtual orbitals
!
         call mem%alloc(S, wf%n_ao, wf%n_ao)
!
         call wf%get_orbital_overlap(PAO_coeff, wf%n_ao, S)
!
         call wf%lovdin_orthonormalization(PAO_coeff, S, wf%n_ao, rank)
!
         call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
         wf%n_ccs_v = rank 
!
!        6. Set inactive virtuals
!
         do mo = wf%n_o + wf%n_cc2_v + 1, wf%n_mo
            do ao = 1, wf%n_ao
!
               wf%orbital_coefficients(ao, mo) = PAO_coeff(ao, mo - wf%n_o - wf%n_cc2_v)
!
            enddo
         enddo
!
         call mem%dealloc(PAO_coeff, wf%n_ao, wf%n_ao)
!
      endif
!
      call mem%dealloc(D, wf%n_ao, wf%n_ao)
!
   end subroutine construct_paos_mlcc2
!
!
   module subroutine diagonalize_M_and_N_mlcc2(wf, T_o, T_v)
!!
!!    Diagonalize M and N
!!    Written by Sarai D. Folkestad, Aug. 2019
!!
!!    Diagonalizes the M and N CNTO matrices.
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: T_v
!
      real(dp), dimension(:), allocatable       :: work, eigenvalues
!
      integer :: info, i, work_size
!
      integer :: count_zero_eigenvalues
!
!     Diagonalize -M and -N to get eigenvectors in correct order
!
      call dscal(wf%n_o**2, -one, T_o, 1)
      call dscal(wf%n_v**2, -one, T_v, 1)
!
      call mem%alloc(work, 1)
      call mem%alloc(eigenvalues, wf%n_o)
!
      call dsyev('V','U',        &
                  wf%n_o,        &
                  T_o,           &
                  wf%n_o,        &
                  eigenvalues,   &
                  work,          &
                  -1,            &
                  info)
!
      work_size = int(work(1))
!
      call mem%dealloc(work, 1)
      call mem%alloc(work, work_size)
!
      call dsyev('V','U',        &
                  wf%n_o,        &
                  T_o,           &
                  wf%n_o,        &
                  eigenvalues,   &
                  work,          &
                  work_size,      &
                  info)
!
      call mem%dealloc(work, work_size)
!
      count_zero_eigenvalues = 0
!
      do i = 1, wf%n_o
!
         if (abs(eigenvalues(i)) .lt. 1.0d-8 ) count_zero_eigenvalues = count_zero_eigenvalues + 1
!
      enddo
!
      if (count_zero_eigenvalues .gt. 0)  &
         call output%printf('Warning: T_o has (i0) zero eigenvalues', ints=[count_zero_eigenvalues],pl='minimal')
!
      call mem%dealloc(eigenvalues, wf%n_o)
!
      if (info .ne. 0) call output%error_msg('Diagonalization of M failed')
!
      call mem%alloc(work, 1)
      call mem%alloc(eigenvalues, wf%n_v)
!
      call dsyev('V','U',        &
                  wf%n_v,        &
                  T_v,           &
                  wf%n_v,        &
                  eigenvalues,   &
                  work,          &
                  -1,            &
                  info)
!
      work_size = int(work(1))
!
      call mem%dealloc(work, 1)
      call mem%alloc(work, work_size)
!
      call dsyev('V','U',        &
                  wf%n_v,        &
                  T_v,           &
                  wf%n_v,        &
                  eigenvalues,   &
                  work,          &
                  work_size,     &
                  info)
!
      call mem%dealloc(work, work_size)
!
      count_zero_eigenvalues = 0
!
      do i = 1, wf%n_v
!
         if (abs(eigenvalues(i)) .lt. 1.0d-8 ) count_zero_eigenvalues = count_zero_eigenvalues + 1
!
      enddo
!
      if (count_zero_eigenvalues .gt. 0)  &
         call output%printf('Warning: T_v has (i0) zero eigenvalues', ints=[count_zero_eigenvalues],pl='minimal')
!
      call mem%dealloc(eigenvalues, wf%n_v)
!
      if (info .ne. 0) call output%error_msg('Diagonalization of N failed')
!
   end subroutine diagonalize_M_and_N_mlcc2
!
!
   module subroutine ccs_calculation_for_cntos_mlcc2(wf, transformation, n_cnto_states, &
                                                         R_ai, cnto_states, omega_ccs)
!!
!!    CCS calculation for CNTOs
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Performs CCS calculation for CNTO (and NTO)
!!    construction.
!!
!
      use diis_cc_gs_class, only: diis_cc_gs
      use davidson_cc_es_class, only: davidson_cc_es
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      character(len=200), intent(in) :: transformation
!
      integer, intent(in) :: n_cnto_states
!
      real(dp), dimension(wf%n_v, wf%n_o, n_cnto_states), intent(out) :: R_ai
!
      integer, dimension(n_cnto_states), intent(in) :: cnto_states
!
      real(dp), dimension(n_cnto_states), intent(out), optional :: omega_ccs
!
!     Local variables
!
      type(ccs), allocatable :: ccs_wf
!
      type(diis_cc_gs) :: cc_gs_solver_diis
      type(davidson_cc_es) :: cc_es_solver_davidson
!
      type(timings) :: timer, timer_gs, timer_es
!
      integer :: n, n_es
!
      real(dp), dimension(:), allocatable :: all_omega_ccs
!
      call output%printf('Running CCS calculation for NTOs/CNTOs.', fs='(/t3,a)',pl='minimal')
!
      if (.not. input%requested_keyword_in_section('print ccs calculation','mlcc')) call output%mute()
!
!     Run CCS calculation
!
      timer = timings('CCS calculation for NTOs/CNTOs')
      call timer%turn_on()
!
!     1. Preparations (Note that cholesky decomposition is already done)
!
      ccs_wf = ccs(wf%system)
!
      ccs_wf%integrals = mo_integral_tool(ccs_wf%n_o, ccs_wf%n_v, ccs_wf%system%n_J)
!
      call ccs_wf%mo_preparations()
!
!     1. Ground state
!
      timer_gs = timings('Ground state CCS calculation for NTOs/CNTOs')
      call timer_gs%turn_on()
!
      cc_gs_solver_diis = diis_cc_gs(ccs_wf)
      call cc_gs_solver_diis%run(ccs_wf)
      call cc_gs_solver_diis%cleanup(ccs_wf)
!
      call timer_gs%turn_off()
!
      call ccs_wf%integrals%write_t1_cholesky(ccs_wf%t1)
!
!     Excited states
!
      timer_es = timings('Excited state CCS calculation for NTOs/CNTOs')
      call timer_es%turn_on()
!
      cc_es_solver_davidson = davidson_cc_es(transformation, ccs_wf)
      call cc_es_solver_davidson%run(ccs_wf)
      call cc_es_solver_davidson%cleanup(ccs_wf)
!
      call timer_es%turn_off()
!
!     Transfer information to mlcc wavefunction
!
!     1. Integrals
!
      wf%integrals = mo_integral_tool(ccs_wf%integrals)
!
!     2. Excitation vectors
!
      n_es = ccs_wf%get_n_excitation_energies_on_file()
!
      if(n_es .lt. n_cnto_states) call output%error_msg('Requested too many CNTO/NTO states')
!
!     Return only the excitation vectors in cnto_states
!
      do n = 1, n_cnto_states
!
         if(n_es .lt. cnto_states(n)) call output%error_msg('Requested non-existent CNTO/NTO state')
!
         call ccs_wf%read_excited_state(R_ai(:,:,n), cnto_states(n), transformation)   
!
      enddo  
!
!     3. Excitation energies (if requested)
!
      if (present(omega_ccs)) then  
!
         call mem%alloc(all_omega_ccs, n_es)
!
!        Read CCS excitaion energies
!
         call wf%read_excitation_energies(n_es, all_omega_ccs)
!
!        Return only the excitation energies in cnto_states
!
         do n = 1, n_cnto_states
!
            omega_ccs(n) = all_omega_ccs(cnto_states(n))   
!
         enddo  
!
      endif
!
!     Cleanup and print
!
      call ccs_wf%cleanup() 
!
      if (.not. input%requested_keyword_in_section('print ccs calculation','mlcc')) call output%unmute()
!
      call timer%turn_off()
!
      call output%printf('Summary of CCS calculation for NTOs/CNTOs:',fs='(/t3,a)',pl='minimal')
!
      call output%printf('Wall time for CCS ground calculation (sec):   (f20.2)', &
             reals=[timer_gs%get_elapsed_time('wall')], fs='(/t6,a)',pl='minimal')
      call output%printf('CPU time for CCS ground calculation (sec):    (f20.2)', &
             reals=[timer_gs%get_elapsed_time('cpu')], fs='(t6,a)',pl='minimal')
!
      call output%printf('Wall time for CCS excited calculation (sec):  (f20.2)', &
             reals=[timer_es%get_elapsed_time('wall')], fs='(/t6,a)',pl='minimal')
      call output%printf('CPU time for CCS excited calculation (sec):   (f20.2)', &
             reals=[timer_es%get_elapsed_time('cpu')], fs='(t6,a)',pl='minimal')
!
   end subroutine ccs_calculation_for_cntos_mlcc2
!
!
   module subroutine check_orthonormality_of_MOs_mlcc2(wf)
!!
!!    Check orthonormality of MOs
!!    Write Sarai D. Folkestad, Sep 2019
!!
!!    Checks that 
!!
!!       C^T S C = I
!!
!!    to ensure that we have orthonormal MOs.
!!
      implicit none

      class(mlcc2) :: wf
!
      real(dp), dimension(:,:), allocatable :: S, I1, I2
!
      integer :: i, j
!
!     Sanity check on the orbitals
!
      call mem%alloc(S, wf%n_ao, wf%n_ao)
      call wf%get_ao_s_wx(S)
!
      call mem%alloc(I1, wf%n_mo, wf%n_ao)
      call mem%alloc(I2, wf%n_mo, wf%n_mo)
!
      call dgemm('T', 'N',                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  S,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  I1,                        &
                  wf%n_mo)
!
      call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'N',                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  I1,                        &
                  wf%n_mo,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  I2,                        &
                  wf%n_mo)
!
      do i = 1, wf%n_mo
         if (abs(I2(i,i) - 1.0d0) .gt. 1.0d-8) call output%error_msg(trim(wf%name_)//' orbitals are not normal')
      enddo
!
      do i = 1, wf%n_mo
         do j = 1, i-1
            if (abs(I2(i,j)) .gt. 1.0d-8) then
            call output%error_msg(trim(wf%name_)//' orbitals are not orthogonal')
            endif
         enddo
      enddo
!
      call mem%dealloc(I1, wf%n_mo, wf%n_ao)
      call mem%dealloc(I2, wf%n_mo, wf%n_mo)

   end subroutine check_orthonormality_of_MOs_mlcc2
!
end submodule orbitals_mlcc2

!
