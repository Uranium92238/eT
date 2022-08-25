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
submodule (mlcc2_class) orbitals_mlcc2
!
!!
!!    MLCC2 orbitals submodule
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
            call output%printf('m', 'Requested restart for CNTOs, &
                                    &reading orbital transformation matrices')
!
            call wf%read_cnto_transformation_matrices(T_o, T_v)
!
         else
!
            call wf%construct_ccs_cnto_transformation_matrices(T_o, T_v)
            call wf%write_cnto_transformation_matrices(T_o, T_v)
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
         if (wf%nto_restart) then
!
            call output%printf('m', 'Requested restart for NTOs, &
                                    &reading orbital transformation matrix')
!
            call wf%read_nto_transformation_matrix(T_o)
!
         else
!
            call wf%construct_ccs_nto_transformation_matrix(T_o)
            call wf%write_nto_transformation_matrix(T_o)
!
         endif
!
         call wf%construct_mixed_nto_canonical_orbitals(T_o)
!
         wf%n_ccs_o = wf%n_o - wf%n_cc2_o
         wf%n_ccs_v = wf%n_v - wf%n_cc2_v
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
!!                      default: .false.
!!
      use cholesky_orbital_tool_class, only: cholesky_orbital_tool
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      logical, intent(in), optional :: occupied_only
!
      logical  :: occupied_only_local
!
      real(dp), dimension(:,:), allocatable :: C
!
      integer :: first_ao, last_ao, n_active_aos
!
      type(cholesky_orbital_tool), allocatable :: cd_tool_o, cd_tool_v
!
      occupied_only_local = .false.
!
      if (present(occupied_only)) occupied_only_local = occupied_only
!
!     Construct occupied orbitals
!
      call wf%ao%get_aos_in_subset('cc2', first_ao, last_ao)
!
      n_active_aos = last_ao - first_ao + 1
!
      call mem%alloc(C, wf%ao%n, wf%ao%n)
!
      cd_tool_o = cholesky_orbital_tool(wf%ao%n, wf%cholesky_orbital_threshold)
      call cd_tool_o%initialize_density()
      call cd_tool_o%set_density_from_orbitals(wf%orbital_coefficients(:,1:wf%n_o), wf%n_o)
!
!     Active
!
      call cd_tool_o%restricted_decomposition(C, wf%n_cc2_o, n_active_aos, first_ao)
      call wf%set_orbital_coefficients(C(:,1:wf%n_cc2_o), wf%n_cc2_o, 1)
!
!     Inactive
!
      call cd_tool_o%full_decomposition(C, wf%n_ccs_o)
      call wf%set_orbital_coefficients(C(:,1:wf%n_ccs_o), wf%n_ccs_o, wf%n_cc2_o + 1)
!
      call cd_tool_o%cleanup()
!
!     Construct virtual orbitals
!
      if (.not. occupied_only_local) then
!
         cd_tool_v = cholesky_orbital_tool(wf%ao%n, wf%cholesky_orbital_threshold)
         call cd_tool_v%initialize_density()
!
         call cd_tool_v%set_density_from_orbitals(wf%orbital_coefficients(:,wf%n_o + 1 : wf%n_mo), wf%n_v)
!
!        Active
!
         call cd_tool_v%restricted_decomposition(C, wf%n_cc2_v, n_active_aos, first_ao)
         call wf%set_orbital_coefficients(C(:,1:wf%n_cc2_v), wf%n_cc2_v, wf%n_o + 1)
!
!        Inactive
!
         call cd_tool_v%full_decomposition(C, wf%n_ccs_v)
         call wf%set_orbital_coefficients(C(:,1:wf%n_ccs_v), wf%n_ccs_v, wf%n_o + wf%n_cc2_v + 1)
!
         call cd_tool_v%cleanup()
!
      endif
!
      call mem%dealloc(C, wf%ao%n, wf%ao%n)
!
   end subroutine construct_cholesky_orbitals_mlcc2
!
!
   module subroutine construct_block_diagonal_fock_orbitals_mlcc2(wf, n_levels, n_occupied_list, &
                                          n_virtual_list, orbital_coefficients, orbital_energies)
!!
!!    Construct block diagonal Fock MOs
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    This routine constructs the MOs which
!!    block-diagonalize the occupied-occupied
!!    and virutal-virtual blocks of the Fock matrix s.t.
!!    the active-active, and inactive-inactive blocks are
!!    diagonal.
!!
!!    Note that after this routine, the Fock matrix in wf
!!    corresponds to the old basis but the MOs are updated.
!!
      use array_utilities, only : block_diagonalize_symmetric, zero_array
!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      integer, intent(in) :: n_levels
!
      integer, dimension(n_levels), intent(in) :: n_occupied_list, n_virtual_list
!
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(inout) :: orbital_coefficients
      real(dp), dimension(wf%n_mo), intent(inout) :: orbital_energies
!
      real(dp), dimension(:,:), allocatable :: F_oo, F_vv, C_copy
!
      integer :: i, offset
!
      call mem%alloc(F_oo, wf%n_o, wf%n_o)
      call mem%alloc(F_vv, wf%n_v, wf%n_v)
!
      call dcopy(wf%n_o**2, wf%fock_ij, 1, F_oo, 1)
      call dcopy(wf%n_v**2, wf%fock_ab, 1, F_vv, 1)
!
!     Block diagonal occupied-occupied Fock
!
      call block_diagonalize_symmetric(F_oo, wf%n_o, n_levels, n_occupied_list, &
                                       orbital_energies(1:wf%n_o))
!
!     Block diagonal virtual-virtual Fock
!
      call block_diagonalize_symmetric(F_vv, wf%n_v, n_levels, n_virtual_list, &
                                    orbital_energies(wf%n_o+1:wf%n_mo))
!
!     Transform blocks
!
      call mem%alloc(C_copy, wf%ao%n, wf%n_mo)
      call dcopy(wf%n_mo*wf%ao%n, orbital_coefficients, 1, C_copy, 1)
!
      call zero_array(orbital_coefficients, wf%n_mo*wf%ao%n)
!
      offset = 1
!
      do i = 1, n_levels
!
         if (n_occupied_list(i) == 0 ) cycle
!
         call dgemm('N', 'N',                      &
                  wf%ao%n,                         &
                  n_occupied_list(i),              &
                  n_occupied_list(i),              &
                  one,                             &
                  C_copy(1, offset),               &
                  wf%ao%n,                         &
                  F_oo(offset, offset),            &
                  wf%n_o,                          &
                  one,                             &
                  orbital_coefficients(1, offset), &
                  wf%ao%n)
!
         offset = offset + n_occupied_list(i)
!
      enddo
!
      offset = 1
!
      do i = 1, n_levels
!
         if (n_virtual_list(i) == 0 ) cycle
!
         call dgemm('N', 'N',                               &
                  wf%ao%n,                                  &
                  n_virtual_list(i),                        &
                  n_virtual_list(i),                        &
                  one,                                      &
                  C_copy(1, offset + wf%n_o),               &
                  wf%ao%n,                                  &
                  F_vv(offset, offset),                     &
                  wf%n_v,                                   &
                  one,                                      &
                  orbital_coefficients(1, offset + wf%n_o), &
                  wf%ao%n)
!
         offset = offset + n_virtual_list(i)
!
      enddo
!
      call mem%dealloc(C_copy, wf%ao%n, wf%n_mo)
      call mem%dealloc(F_oo, wf%n_o, wf%n_o)
      call mem%dealloc(F_vv, wf%n_v, wf%n_v)
!
   end subroutine construct_block_diagonal_fock_orbitals_mlcc2
!
!
   module subroutine construct_M_and_N_cnto_mlcc2(wf, R_ai, R_aibj, M, N, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij += ( sum_a R_ai R_aj + 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibl R_ajbl )
!!       N_ab += ( sum_i R_ai R_bi + 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    Used to construct CNTOs.
!!
!!    set_to_zero determines if M and N are set to zero initially. This makes it possible for
!!    the routine to be used to add to M and N if we are using more than one excitation vector
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
      real(dp) :: zero_or_one
!
      zero_or_one = one
!
      if (set_to_zero) zero_or_one = zero
!
!     1. Singles contribution
!
!     M_ij += sum_a R_ai R_aj
!
      call dgemm('T', 'N',       &
                  wf%n_o,        &
                  wf%n_o,        &
                  wf%n_v,        &
                  one,           &
                  R_ai,          & ! R_ai
                  wf%n_v,        &
                  R_ai,          & ! R_aj
                  wf%n_v,        &
                  zero_or_one,   &
                  M,             & ! M_ij
                  wf%n_o)
!
!     N_ab += sum_i R_ai R_aj
!
      call dgemm('N', 'T',       &
                  wf%n_v,        &
                  wf%n_v,        &
                  wf%n_o,        &
                  one,           &
                  R_ai,          & ! R_ai
                  wf%n_v,        &
                  R_ai,          & ! R_bi
                  wf%n_v,        &
                  zero_or_one,   &
                  N,             & ! N_ab
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
      use array_utilities, only: zero_array
!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: T_v
!
      real(dp), dimension(:,:), allocatable     :: MO_coeff
!
      call mem%alloc(MO_coeff, wf%ao%n, wf%n_mo)
      call dcopy((wf%ao%n)*(wf%n_mo), wf%orbital_coefficients, 1, MO_coeff, 1)
!
      call zero_array(wf%orbital_coefficients, (wf%ao%n)*(wf%n_mo))
!
      call dgemm('N', 'N',                   &
                  wf%ao%n,                   &
                  wf%n_o,                    &
                  wf%n_o,                    &
                  one,                       &
                  MO_coeff,                  &
                  wf%ao%n,                   &
                  T_o,                       &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n)
!
      call dgemm('N', 'N',                                  &
                  wf%ao%n,                                  &
                  wf%n_v,                                   &
                  wf%n_v,                                   &
                  one,                                      &
                  MO_coeff(1, wf%n_o + 1),                  &
                  wf%ao%n,                                  &
                  T_v,                                      &
                  wf%n_v,                                   &
                  one,                                      &
                  wf%orbital_coefficients(1, wf%n_o + 1),   &
                  wf%ao%n)
!
      call mem%dealloc(MO_coeff, wf%ao%n, wf%n_mo)
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
      call wf%T_cnto_o_file%open_('read', 'rewind')
      call wf%T_cnto_v_file%open_('read', 'rewind')
!
      call wf%T_cnto_o_file%read_(T_o, wf%n_o**2)
      call wf%T_cnto_v_file%read_(T_v, wf%n_v**2)
!
      call wf%T_cnto_o_file%close_('keep')
      call wf%T_cnto_v_file%close_('keep')
!
   end subroutine read_cnto_transformation_matrices_mlcc2
!
!
   module subroutine write_cnto_transformation_matrices_mlcc2(wf, T_o, T_v)
!!
!!    Write CNTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Write CNTO transformation matrices.
!!    Used to ensure restart
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(in) :: T_v
!
      call wf%T_cnto_o_file%open_('write', 'rewind')
      call wf%T_cnto_v_file%open_('write', 'rewind')
!
      call wf%T_cnto_o_file%write_(T_o, wf%n_o**2)
      call wf%T_cnto_v_file%write_(T_v, wf%n_v**2)
!
      call wf%T_cnto_o_file%close_('keep')
      call wf%T_cnto_v_file%close_('keep')
!
   end subroutine write_cnto_transformation_matrices_mlcc2
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
      use global_in, only: input
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
!
      logical :: set_to_zero
!
      character(len=200) :: r_or_l
!
      type(direct_stream_file) :: doubles_file
!
      type(timings) :: timer
!
      timer = timings('Construct CNTO matrices')
      call timer%turn_on()
!
      n_cnto_states = size(wf%cnto_states)
!
      r_or_l = 'right'
!
      if (input%is_keyword_present('left eigenvectors', 'solver cc es')) r_or_l = 'left'

      call mem%alloc(R_ai, wf%n_v, wf%n_o, n_cnto_states)
      call mem%alloc(omega_ccs, n_cnto_states)
!
!     Run CCS calculation
!
      call wf%ccs_calculation_for_cntos(r_or_l, n_cnto_states, R_ai, wf%cnto_states, omega_ccs)
!
      set_to_zero = .true.
!
!     Prepare file
!
      doubles_file = direct_stream_file('approximate_doubles', wf%n_v*wf%n_o**2)
!
      call mem%alloc(R_ai_k, wf%n_v, wf%n_o)
!
      do k = 1, n_cnto_states
!
         call dcopy(wf%n_t1, R_ai(1,1,k), 1, R_ai_k, 1)
!
!        Construct approximate double excitation vector
!
         call wf%approximate_double_excitation_vectors(R_ai_k, omega_ccs(k), doubles_file)
!
!        Add contribution to M and N
!
         call wf%construct_M_and_N_singles_cnto(R_ai_k, T_o, T_v, set_to_zero)
         call wf%add_doubles_M_and_N_cnto(T_o, T_v, doubles_file)
!
         set_to_zero = .false.
!
      enddo
!
      call mem%dealloc(omega_ccs, n_cnto_states)
      call mem%dealloc(R_ai_k, wf%n_v, wf%n_o)
      call mem%dealloc(R_ai, wf%n_v, wf%n_o, n_cnto_states)
!
      call wf%diagonalize_M_and_N(T_o, T_v)
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
!!    set_to_zero determines if contributions are added to, or overwrites M and N.
!!
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
      real(dp) :: zero_or_one
!
      zero_or_one = one
!
      if (set_to_zero) zero_or_one = zero
!
!     1. Singles contribution
!
!     M_ij += sum_a R_ai R_aj
!
      call dgemm('T', 'N',       &
                  wf%n_o,        &
                  wf%n_o,        &
                  wf%n_v,        &
                  one,           &
                  R_ai,          & ! R_ai
                  wf%n_v,        &
                  R_ai,          & ! R_aj
                  wf%n_v,        &
                  zero_or_one,   &
                  M,             & ! M_ij
                  wf%n_o)
!
   end subroutine construct_M_nto_mlcc2
!
!
   module subroutine read_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Read NTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Read NTO transformation matrices.
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
!
      call wf%T_nto_o_file%open_('read', 'rewind')
!
      call wf%T_nto_o_file%read_(T_o, wf%n_o**2)
!
      call wf%T_nto_o_file%close_('keep')
!
   end subroutine read_nto_transformation_matrix_mlcc2
!
!
   module subroutine write_nto_transformation_matrix_mlcc2(wf, T_o)
!!
!!    Write NTO transformation matrices
!!    Written by Sarai D. Folkestad, Jun 2019
!!
!!    Write NTO transformation matrices.
!!    Used to ensure restart
!!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
!
      call wf%T_nto_o_file%open_('write', 'rewind')
!
      call wf%T_nto_o_file%write_(T_o, wf%n_o**2)
!
      call wf%T_nto_o_file%close_('keep')
!
   end subroutine write_nto_transformation_matrix_mlcc2
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
!
      use array_utilities, only: diagonalize_symmetric
      use global_in, only: input
!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(out) :: T_o
!
      integer :: k, n_nto_states
!
      real(dp), dimension(:,:), allocatable     :: R_ai_k
      real(dp), dimension(:,:,:), allocatable   :: R_ai
!
      logical :: set_to_zero

      character(len=200) :: r_or_l
!
      n_nto_states = size(wf%nto_states)
!
      r_or_l = 'right'
!
      if (input%is_keyword_present('left eigenvectors', 'solver cc es')) r_or_l = 'left'
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
      call diagonalize_symmetric(T_o, wf%n_o)
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
      use array_utilities, only: zero_array
!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(in) :: T_o
!
      real(dp), dimension(:,:), allocatable     :: MO_coeff
!
      call mem%alloc(MO_coeff, wf%ao%n, wf%n_mo)
      call dcopy((wf%ao%n)*(wf%n_mo), wf%orbital_coefficients, 1, MO_coeff, 1)
!
      call zero_array(wf%orbital_coefficients, (wf%ao%n)*(wf%n_o))
!
      call dgemm('N', 'N',                   &
                  wf%ao%n,                   &
                  wf%n_o,                    &
                  wf%n_o,                    &
                  one,                       &
                  MO_coeff,                  &
                  wf%ao%n,                   &
                  T_o,                       &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n)
!
      call mem%dealloc(MO_coeff, wf%ao%n, wf%n_mo)
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
      call wf%ao%get_aos_in_subset('cc2', first_ao, last_ao)
!
      n_active_aos = last_ao - first_ao + 1
!
!     1. Set up occupied density
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'T',                   &
                  wf%ao%n,                   &
                  wf%ao%n,                   &
                  wf%n_o,                    &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  zero,                      &
                  D,                         &
                  wf%ao%n)
!
      if (wf%exists_frozen_fock_terms) then
!
         call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!
      endif
!
!     2. Construct PAOs for active atoms
!
      call mem%alloc(PAO_coeff, wf%ao%n, n_active_aos)
!
      call wf%project_atomic_orbitals(D, PAO_coeff, n_active_aos, first_ao)
!
!     3. Orthonormalize PAOs to get active virtual orbitals
!
      call mem%alloc(S, n_active_aos, n_active_aos)
!
      call wf%get_orbital_overlap(PAO_coeff, n_active_aos, S)
!
      call wf%lowdin_orthonormalization(PAO_coeff, S, n_active_aos, rank)
!
      call mem%dealloc(S, n_active_aos, n_active_aos)
!
      wf%n_cc2_v = rank
!
!     Set the active virtual orbital coefficients
!
      do mo = wf%n_o + 1, wf%n_o + wf%n_cc2_v
         do ao = 1, wf%ao%n
!
            wf%orbital_coefficients(ao, mo) = PAO_coeff(ao, mo - wf%n_o)
!
         enddo
      enddo
!
      call mem%dealloc(PAO_coeff, wf%ao%n, n_active_aos)
!
      if (rank .lt. wf%n_v) then
!
!        Set virtual inactive
!
!        4. Construct M = sum_p C_αp C_βp  for p = 1, n_o + n_cc2_v
!
         call dgemm('N', 'T',                   &
                     wf%ao%n,                   &
                     wf%ao%n,                   &
                     wf%n_o + wf%n_cc2_v,       &
                     one,                       &
                     wf%orbital_coefficients,   &
                     wf%ao%n,                   &
                     wf%orbital_coefficients,   &
                     wf%ao%n,                   &
                     zero,                      &
                     D,                         &
                     wf%ao%n)
!
         if (wf%exists_frozen_fock_terms) then
!
            call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!
         endif
!
!        Construct PAOs for the remaining virtual orbitals
!
         call mem%alloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
         call wf%project_atomic_orbitals(D, PAO_coeff, wf%ao%n)
!
!        5. Orthonormalize PAOs to get inactive virtual orbitals
!
         call mem%alloc(S, wf%ao%n, wf%ao%n)
!
         call wf%get_orbital_overlap(PAO_coeff, wf%ao%n, S)
!
         call wf%lowdin_orthonormalization(PAO_coeff, S, wf%ao%n, rank)
!
         call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
         wf%n_ccs_v = rank
!
!        6. Set inactive virtuals
!
         do mo = wf%n_o + wf%n_cc2_v + 1, wf%n_mo
            do ao = 1, wf%ao%n
!
               wf%orbital_coefficients(ao, mo) = PAO_coeff(ao, mo - wf%n_o - wf%n_cc2_v)
!
            enddo
         enddo
!
         call mem%dealloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
      endif
!
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
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
!
      use array_utilities, only: diagonalize_symmetric
!
      implicit none
!
      class(mlcc2) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: T_o
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: T_v
!
      real(dp), dimension(:), allocatable :: eigenvalues
!
      integer :: i
!
      integer :: count_zero_eigenvalues
!
!     Diagonalize -M and -N to get eigenvectors in correct order
!
      call dscal(wf%n_o**2, -one, T_o, 1)
      call dscal(wf%n_v**2, -one, T_v, 1)
!
      call mem%alloc(eigenvalues, wf%n_o)
      call diagonalize_symmetric(T_o, wf%n_o, eigenvalues)
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
         call output%warning_msg('T_o has (i0) zero eigenvalues', ints=[count_zero_eigenvalues])
!
      call mem%dealloc(eigenvalues, wf%n_o)
!
      call mem%alloc(eigenvalues, wf%n_v)
      call diagonalize_symmetric(T_v, wf%n_v, eigenvalues)
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
         call output%warning_msg('T_v has (i0) zero eigenvalues', ints=[count_zero_eigenvalues])
!
      call mem%dealloc(eigenvalues, wf%n_v)
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
      use amplitude_updater_class,        only: amplitude_updater
      use quasi_newton_updater_class,     only: quasi_newton_updater
      use diis_cc_gs_class,               only: diis_cc_gs
      use abstract_solver_class,          only: abstract_solver
      use eri_tool_class,                 only: eri_tool
      use eri_cholesky_disk_class,        only: eri_cholesky_disk
      use cc_es_amplitudes_solver_factory_class, only: cc_es_amplitudes_solver_factory
!
      use global_in, only: input
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
      character(len=200) :: storage
      type(cc_es_amplitudes_solver_factory), allocatable :: solver_factory
!
!     Local variables
!
      type(ccs), allocatable :: ccs_wf
!
      type(diis_cc_gs), allocatable :: cc_gs_solver_diis
      class(amplitude_updater), allocatable :: t_updater
!
      class(abstract_solver), allocatable :: cc_es_solver
!
      type(timings) :: timer, timer_gs, timer_es
!
      integer :: n, n_es
!
      real(dp), dimension(:), allocatable :: all_omega_ccs
!
      call output%printf('m', 'Running CCS calculation for NTOs/CNTOs.', fs='(/t3,a)')
!
      if (.not. input%is_keyword_present('print ccs calculation','mlcc')) call output%mute()
!
!     Run CCS calculation
!
      timer = timings('CCS calculation for NTOs/CNTOs')
      call timer%turn_on()
!
!     1. Preparations (Note that cholesky decomposition is already done)
!
      allocate(ccs::ccs_wf)
      call ccs_wf%initialize(wf)
!
      call ccs_wf%integral_preparations(wf%L_mo%n_J)
!
      call ccs_wf%L_mo%set_equal_to(wf%L_mo)
      call ccs_wf%L_t1%set_equal_to(wf%L_mo)
!
      call ccs_wf%mo_preparations()
!
!     1. Ground state
!
      timer_gs = timings('Ground state CCS calculation for NTOs/CNTOs')
      call timer_gs%turn_on()
!
      call ccs_wf%initialize_ground_state_files()
      t_updater = quasi_newton_updater(n_amplitudes     = ccs_wf%n_gs_amplitudes, &
                                       scale_amplitudes = .true., &
                                       scale_residual   = .true.)
!
      storage = 'disk'
      call input%get_keyword('storage', 'solver cc gs', storage)
!
      cc_gs_solver_diis = diis_cc_gs(wf        = ccs_wf,    &
                                     restart   = .false.,   &
                                     t_updater = t_updater, &
                                     storage   = storage)
!
      call cc_gs_solver_diis%run()
      call cc_gs_solver_diis%cleanup()
!
      call timer_gs%turn_off()
!
!     Excited states
!
      timer_es = timings('Excited state CCS calculation for NTOs/CNTOs')
      call timer_es%turn_on()
!
      call ccs_wf%construct_fock('es')
!
      solver_factory = cc_es_amplitudes_solver_factory(ccs_wf%name_, transformation, restart=.false.)
      call solver_factory%create(ccs_wf, cc_es_solver)
!
      call ccs_wf%initialize_excited_state_files()
      call ccs_wf%initialize_excitation_energies()
!
      call cc_es_solver%run()
      call cc_es_solver%cleanup()
!
      call timer_es%turn_off()
!
!     Transfer information to mlcc wavefunction
!
!     1. Excitation vectors
!
      n_es = ccs_wf%n_singlet_states
!
      call mem%alloc(all_omega_ccs, n_es)
!
      if(n_es .lt. n_cnto_states) call output%error_msg('Requested too many CNTO/NTO states')
!
!     Return only the excitation vectors in cnto_states
!
      do n = 1, n_cnto_states
!
         if(n_es .lt. cnto_states(n)) call output%error_msg('Requested non-existent CNTO/NTO state')
!
         call ccs_wf%read_excited_state(R_ai(:,:,n),        &
                                        cnto_states(n),     &
                                        cnto_states(n),     &
                                        transformation,     &
                                        all_omega_ccs(n))
!
      enddo
!
!     2. Excitation energies (if requested)
!
      if (present(omega_ccs)) then
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
      call mem%dealloc(all_omega_ccs, n_es)
!
!     Cleanup and print
!
      call ccs_wf%cleanup()
!
      if (.not. input%is_keyword_present('print ccs calculation','mlcc')) &
         call output%unmute()
!
      call timer%turn_off()
!
      call output%printf('m', '- Summary of CCS calculation for NTOs/CNTOs:',fs='(/t3,a)')
!
      call output%printf('m', 'Wall time for CCS ground calculation (sec):   (f20.2)', &
                         reals=[timer_gs%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'CPU time for CCS ground calculation (sec):    (f20.2)', &
                         reals=[timer_gs%get_elapsed_time('cpu')], fs='(t6,a)')
!
      call output%printf('m', 'Wall time for CCS excited calculation (sec):  (f20.2)', &
                         reals=[timer_es%get_elapsed_time('wall')], fs='(/t6,a)')
      call output%printf('m', 'CPU time for CCS excited calculation (sec):   (f20.2)', &
                         reals=[timer_es%get_elapsed_time('cpu')], fs='(t6,a)')
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
      call mem%alloc(S, wf%ao%n, wf%ao%n)
      call wf%ao%get_oei('overlap', S)
!
      call mem%alloc(I1, wf%n_mo, wf%ao%n)
      call mem%alloc(I2, wf%n_mo, wf%n_mo)
!
      call dgemm('T', 'N',                   &
                  wf%n_mo,                   &
                  wf%ao%n,                   &
                  wf%ao%n,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  S,                         &
                  wf%ao%n,                   &
                  zero,                      &
                  I1,                        &
                  wf%n_mo)
!
      call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'N',                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  wf%ao%n,                   &
                  one,                       &
                  I1,                        &
                  wf%n_mo,                   &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  zero,                      &
                  I2,                        &
                  wf%n_mo)
!
      do i = 1, wf%n_mo
         if (abs(I2(i,i) - one) .gt. 1.0d-6) then
            call output%error_msg(trim(wf%name_)//' orbitals are not normal')
         endif
      enddo
!
      do i = 1, wf%n_mo
         do j = 1, i-1
            if (abs(I2(i,j)) .gt. 1.0d-6) then
               call output%error_msg(trim(wf%name_)//' orbitals are not orthogonal')
            endif
         enddo
      enddo
!
      call mem%dealloc(I1, wf%n_mo, wf%ao%n)
      call mem%dealloc(I2, wf%n_mo, wf%n_mo)

   end subroutine check_orthonormality_of_MOs_mlcc2
!
!
   module subroutine add_doubles_M_and_N_cnto_mlcc2(wf, M, N, doubles_file)
!!
!!    Add doubles M and N CNTO
!!    Written by Sarai D. Folkestad, Feb 2020
!!
!!    The CNTO matrices are defined as
!!
!!       M_ij = ( sum_a R_ai R_aj + 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibl R_ajbl )
!!       N_ab = ( sum_i R_ai R_bi + 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    In this routine the doubles part
!!
!!       M_ij += ( 1/2 sum_abl(1 + δ_ai,bl δ_i,j) R_aibl R_ajbl )
!!       N_ab += ( 1/2 sum_cij(1 + δ_ai,cj δ_a,b) R_aicj R_bicj )
!!
!!    is added for the R_aibj on doubles_file
!!
      use batching_index_class, only: batching_index
!
      implicit none
!
      class(mlcc2) :: wf
!
      type(direct_stream_file) :: doubles_file
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
!
      type(batching_index) :: batch_a, batch_c
!
      integer :: current_a_batch, current_c_batch
      integer :: req0, req1_a, req1_c, req2
!
      integer :: a, i
!
      real(dp), dimension(:,:,:,:), allocatable :: R_ibja, R_kdlc
!
      call doubles_file%open_('read')
!
      req0 = 0
!
      req1_a = wf%n_o**2*wf%n_v
!
      req1_c = req1_a
!
      req2 = 0
!
!     Initialize batching variables
!
      batch_a = batching_index(wf%n_v)
      batch_c = batching_index(wf%n_v)

      call mem%batch_setup(batch_a, batch_c, req0, req1_a, req1_c, req2, &
                           tag='add_doubles_M_and_N_cnto_mlcc2')
!
      call mem%alloc(R_ibja, wf%n_o, wf%n_v, wf%n_o, batch_a%max_length)
      call mem%alloc(R_kdlc, wf%n_o, wf%n_v, wf%n_o, batch_c%max_length)
!
      do current_a_batch = 1, batch_a%num_batches
!
         call batch_a%determine_limits(current_a_batch)
!
         call doubles_file%read_(R_ibja(:,:,:,1:batch_a%length), &
                                 batch_a%first, batch_a%get_last())
!
!        M_ij = 1/2 (1 + delta_ai,bk delta_i,j) R_jbka R_ibka
!
         call dgemm('N', 'T',                      &
                     wf%n_o,                       &
                     wf%n_o,                       &
                     wf%n_v*wf%n_o*batch_a%length, &
                     half,                         &
                     R_ibja,                       &
                     wf%n_o,                       &
                     R_ibja,                       &
                     wf%n_o,                       &
                     one,                          &
                     M,                            &
                     wf%n_o)
!
!$omp parallel do private (i, a)
         do i = 1, wf%n_o
            do a = 1, batch_a%length
!
               M(i,i) = M(i,i) + half*R_ibja(i, a + batch_a%first - 1, i, a)**2
!
            enddo
         enddo
!$omp end parallel do
!
         do current_c_batch = 1, batch_c%num_batches
!
            call batch_c%determine_limits(current_c_batch)
!
            call doubles_file%read_(R_kdlc(:,:,:,1:batch_c%length), &
                                    batch_c%first, batch_c%get_last())
!
!           N_ac += 1/2 sum_dkl(1 + delta_ak,dl delta_a,c) R_akdl R_ckdl
!
            call dgemm('T', 'N',                         &
                        batch_a%length,                  &
                        batch_c%length,                  &
                        wf%n_o**2*wf%n_v,                &
                        half,                            &
                        R_ibja,                          &
                        wf%n_o**2*wf%n_v,                &
                        R_kdlc,                          &
                        wf%n_o**2*wf%n_v,                &
                        one,                             &
                        N(batch_a%first, batch_c%first), &
                        wf%n_v)
!
         enddo
!
!$omp parallel do private (a, i)
         do a = 1, batch_a%length
            do i = 1, wf%n_o
!
               N(a + batch_a%first - 1, a + batch_a%first - 1) &
                  = N(a + batch_a%first - 1, a + batch_a%first - 1) &
                     + half*R_ibja(i, a + batch_a%first - 1, i, a)**2
!
            enddo
         enddo
!$omp end parallel do
!
      enddo
!
!
      call mem%dealloc(R_kdlc, wf%n_o, wf%n_v, wf%n_o, batch_c%max_length)
      call mem%dealloc(R_ibja, wf%n_o, wf%n_v, wf%n_o, batch_a%max_length)
!
      call mem%batch_finalize()
!
      call doubles_file%close_('keep')
!
   end subroutine add_doubles_M_and_N_cnto_mlcc2
!
!
   module subroutine construct_M_and_N_singles_cnto_mlcc2(wf, R_ai, M, N, set_to_zero)
!!
!!    Construct M and N
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Constructs the M and N matrices,
!!
!!       M_ij = ( sum_a R_ai R_aj )
!!       N_ab = ( sum_i R_ai R_bi )
!!
!!    Used to construct CNTOs.
!!
!!    set_to_zero determines if contributions are added to, or overwrites M and N.
!!
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(in) :: R_ai
!
      real(dp), dimension(wf%n_o, wf%n_o), intent(inout) :: M
      real(dp), dimension(wf%n_v, wf%n_v), intent(inout) :: N
!
      logical, intent(in) :: set_to_zero
!
      real(dp) :: zero_or_one
!
      zero_or_one = one
!
      if (set_to_zero) zero_or_one = zero
!
!     1. Singles contribution
!
!     M_ij += sum_a R_ai R_aj
!
      call dgemm('T', 'N',       &
                  wf%n_o,        &
                  wf%n_o,        &
                  wf%n_v,        &
                  one,           &
                  R_ai,          & ! R_ai
                  wf%n_v,        &
                  R_ai,          & ! R_aj
                  wf%n_v,        &
                  zero_or_one,   &
                  M,             & ! M_ij
                  wf%n_o)
!
!     N_ab += sum_i R_ai R_aj
!
      call dgemm('N', 'T',       &
                  wf%n_v,        &
                  wf%n_v,        &
                  wf%n_o,        &
                  one,           &
                  R_ai,          & ! R_ai
                  wf%n_v,        &
                  R_ai,          & ! R_bi
                  wf%n_v,        &
                  zero_or_one,   &
                  N,             & ! N_ab
                  wf%n_v)
!
   end subroutine construct_M_and_N_singles_cnto_mlcc2
!
!
   module subroutine construct_semicanonical_mlcc_orbitals_mlcc2(wf)
!!
!!    Construct semicanonical orbitals
!!    Written by Sarai D. Folkestad
!!
!!    Wrapper to construct orbitals that block diagonalizes the fock matrix
!!    for the different orbitals
!!
      implicit none
!
      class(mlcc2), intent(inout) :: wf
!
      integer :: n_levels
!
      n_levels = 2
!
      call wf%construct_block_diagonal_fock_orbitals(n_levels, [wf%n_cc2_o, wf%n_ccs_o],   &
                                          [wf%n_cc2_v, wf%n_ccs_v], wf%orbital_coefficients, &
                                          wf%orbital_energies)
!
   end subroutine construct_semicanonical_mlcc_orbitals_mlcc2
!
!
end submodule orbitals_mlcc2
!
