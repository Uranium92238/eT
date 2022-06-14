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
submodule (hf_class) frozen_orbital_hf
!
!!
!!    Frozen orbital (HF)
!!
!!    Submodule containing routines for frozen orbitals
!!
!
   implicit none
!
!
contains
!
!
   module subroutine prepare_mos_hf(wf)
!!
!!    Prepare MOs
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
!!    This routine prepares the MOs for coupled cluster
!!    in the cases where there is a reduction in the
!!    number of MOs in CC compared to HF
!!
!!    Examples of this is the frozen core
!!    approximation and if CC is only done
!!    for a localized region of a large molecule
!!    which has been treated at HF level of theory.
!!
      implicit none
!
      class(hf) :: wf
!
!     Destruct MO quantities in the old MO dimension, if they are allocated,
!     before n_mo changes
!
      call wf%destruct_mo_fock()
!
!     Eliminate the core orbitals if frozen core requested
!
!     MO coefficients for core orbitals are placed in
!     wf%orbital_coefficients_fc and removed from wf%orbital_coefficients
!     the number of frozen core orbitals is wf%n_frozen_core_orbitals
!
      if (wf%frozen_core) call wf%remove_core_orbitals()
!
!     Cholesky decomposition of density for reduced space CC calculation
!
!     MO coefficients for frozen hf orbitals now placed in
!     wf%orbital_coefficients_frozen_hf and removed from wf%orbital_coefficients
!     the number of frozen hf orbitals is wf%n_frozen_hf_orbitals
!
      if (wf%frozen_hf_mos) call wf%remove_frozen_hf_orbitals()
!
   end subroutine prepare_mos_hf
!
!
   module subroutine remove_core_orbitals_hf(wf)
!!
!!    Remove core orbitals
!!    Written by Sarai D. Folkestad, Sep 2018
!!
!!    - Removes core orbitals from wf%orbital_coefficients
!!
!!       Removes 1s for C - Mg
!!       Removes 1s, 2s, 2p for Al - Zn
!!
!!    - The core orbitals are stored in wf%orbital_coefficients_fc
!!    - The number of frozen core orbitals is wf%n_frozen_core_orbitals on exit
!!    - On exit wf%n_mo and wf%n_o are updated to not include the core orbitals
!!
!
      use array_utilities, only : get_abs_max_w_index
!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_copy
      real(dp), dimension(:), allocatable :: orbital_energies_copy
!
      integer :: I, ao, mo, K
!
      logical, dimension(:), allocatable :: freeze_atom
!
      integer :: index_max
      real(dp) :: max_
!
!     Figure out how many core orbitals we have, and which atoms to freeze
!
      call mem%alloc(freeze_atom, wf%get_n_active_hf_atoms())
!
      call wf%ao%get_frozen_cores_and_n_frozen_orbitals(frozen      = freeze_atom,                 &
                                                        n_frozen_ao = wf%n_frozen_core_orbitals,   &
                                                        first       = 1,                           &
                                                        last        = wf%get_n_active_hf_atoms())
!
      call mem%alloc(orbital_coefficients_copy, wf%ao%n, wf%n_mo)
!
      call dcopy(wf%n_mo*wf%ao%n, wf%orbital_coefficients, 1, orbital_coefficients_copy, 1)
!
      call wf%destruct_orbital_coefficients()
!
      call mem%alloc(wf%orbital_coefficients, wf%ao%n, wf%n_mo - wf%n_frozen_core_orbitals)
!
!$omp parallel do private (ao, mo)
      do mo = 1, wf%n_mo - wf%n_frozen_core_orbitals
         do ao = 1, wf%ao%n
!
            wf%orbital_coefficients(ao, mo) = &
                  orbital_coefficients_copy(ao, wf%n_frozen_core_orbitals + mo)
!
         enddo
      enddo
!$omp end parallel do
!
!     Keep frozen core orbital orbital coefficients
!
      call wf%initialize_orbital_coefficients_fc()
!
!$omp parallel do private(K, i) collapse(2)
      do K = 1, wf%n_frozen_core_orbitals
         do i = 1, wf%ao%n
!
            wf%orbital_coefficients_fc(i, K) = orbital_coefficients_copy(i, K)
!
         enddo
      enddo
!$omp end parallel do
!
!     Check for crossover:
!
!     If the largest AO weight on a frozen MO does not belong to an
!     atom we are supposed to freeze, then we stop.
!
      do i = 1, wf%n_frozen_core_orbitals
!
         call get_abs_max_w_index(wf%orbital_coefficients_fc(:,i), wf%ao%n, max_, index_max)
!
         if (.not. freeze_atom(wf%ao%ao_to_center(index_max))) &
            call output%error_msg('Frozen MO has large contribution from a non-frozen atom.')
!
      enddo
!
      call mem%dealloc(freeze_atom, wf%get_n_active_hf_atoms())
      call mem%dealloc(orbital_coefficients_copy, wf%ao%n, wf%n_mo)
!
      call mem%alloc(orbital_energies_copy, wf%n_mo)
!
      call dcopy(wf%n_mo, wf%orbital_energies, 1, orbital_energies_copy, 1)
!
      call wf%destruct_orbital_energies()
      call mem%alloc(wf%orbital_energies, wf%n_mo - wf%n_frozen_core_orbitals)
!
!$omp parallel do private (mo)
      do mo = 1, wf%n_mo - wf%n_frozen_core_orbitals
!
         wf%orbital_energies(mo) = orbital_energies_copy(wf%n_frozen_core_orbitals + mo)
!
      enddo
!$omp end parallel do
!
      call mem%dealloc(orbital_energies_copy, wf%n_mo)
!
      wf%n_mo = wf%n_mo  - wf%n_frozen_core_orbitals
      wf%n_o  = wf%n_o  - wf%n_frozen_core_orbitals
!
      call output%printf('m', '- Preparation for frozen core approximation', &
                         fs='(/t3,a)')
!

      call output%printf('m', 'There are (i0) frozen core orbitals.', &
                        ints=[wf%n_frozen_core_orbitals], fs='(/t6,a)')
!
   end subroutine remove_core_orbitals_hf
!
!
   module subroutine remove_frozen_hf_orbitals_hf(wf)
!!
!!    Remove frozen hf orbitals
!!    Written by Sarai D. Folkestad, Feb 2019
!!    Added and modified for HF by Ida-Marie Hoyvik, Oct 2019
!!
!!    - Removes frozen HF molecular orbitals from wf%orbital_coefficients
!!      and places them in wf%orbital_coefficients_frozen_hf
!!    - The number of frozen occupied HF orbitals is wf%n_frozen_hf_o
!!      on exit
!!    - On exit wf%n_mo, wf%n_o and wf%n_v are updated not to include
!!      frozen hf orbitals
!!
      use cholesky_orbital_tool_class, only: cholesky_orbital_tool
!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: D, orbitals_copy, S, PAO_coeff, C_o
!
      integer :: first_ao, last_ao, i, n_active_aos, a, first_inactive_ao
!
      integer :: ao, rank, n_inactive_centers, n_hf_centers
!
      type(cholesky_orbital_tool), allocatable :: cd_tool_o
!
!     Construct active occupied orbitals
!
!     0. Determine active ao list -> all active atoms that are not hf active
!
!     First some sanity checks on the active atoms.
!
!     Are there active atoms spaces?
!
      n_inactive_centers = wf%ao%get_n_centers_in_subset('unclassified')
!
      if (n_inactive_centers .eq. wf%n_atomic_centers) &
         call output%error_msg('frozen HF orbitals requested, but no active atoms on input!')
!
!     If there is only one active atoms space, is it HF active?
!
      n_hf_centers = wf%ao%get_n_centers_in_subset('hf')
!
      if (n_hf_centers + n_inactive_centers .eq. wf%n_atomic_centers) &
         call output%error_msg('frozen HF orbitals requested, but no active CC atoms on input!')
!
!     Get the first and last active AO indices for CC
!
!     Centers in eT are ordered according to levels:
!
!        CC3, CCSD, CC2, CCS, HF, unclassified
!
      first_ao = 1
!
      if (wf%ao%is_center_subset('hf')) then
!
!        First inactive is in the 'hf' subset
!
         call wf%ao%get_aos_in_subset('hf', first=first_inactive_ao)
!
         last_ao = first_inactive_ao - 1
!
      elseif (wf%ao%is_center_subset('unclassified')) then
!
!        First inactive is in the 'unclassified' subset
!
         call wf%ao%get_aos_in_subset('unclassified', first=first_inactive_ao)
!
         last_ao = first_inactive_ao - 1
!
      else
!
!        No inactive; last active AO is equal to the full number of AOs
!
         last_ao = wf%ao%n
!
      endif
!
      n_active_aos = last_ao
!
!     Occupied orbitals
!
!     1. Set up active occupied density
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
!     2. Construct occupied active orbitals
!
      call mem%alloc(C_o, wf%ao%n, wf%ao%n)
!
      cd_tool_o = cholesky_orbital_tool(wf%ao%n, wf%cholesky_orbital_threshold)
      call cd_tool_o%initialize_density()
      call cd_tool_o%set_density_from_density(D, wf%n_o)
!
      call cd_tool_o%restricted_decomposition(C_o, wf%n_o, n_active_aos, 1)
      call wf%set_orbital_coefficients(C_o(:,1:wf%n_o), wf%n_o, 1)
!
      call cd_tool_o%full_decomposition(C_o, wf%n_frozen_hf_o)
      call wf%set_orbital_coefficients(C_o(:,1:wf%n_frozen_hf_o), wf%n_frozen_hf_o, wf%n_o + 1)
!
      call mem%dealloc(C_o, wf%ao%n, wf%ao%n)
!
      call cd_tool_o%cleanup()
!
!     Virtual orbitals
!
      call wf%get_full_idempotent_density(D)
!
!     1. Construct PAOs for active atoms
!
      call mem%alloc(PAO_coeff, wf%ao%n, n_active_aos)
!
      call wf%project_atomic_orbitals(D, PAO_coeff, n_active_aos, first_ao)
!
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
!     2. Orthonormalize PAOs to get active virtual orbitals
!
      call mem%alloc(S, n_active_aos, n_active_aos)
!
      call wf%get_orbital_overlap(PAO_coeff, n_active_aos, S)
!
      call wf%lowdin_orthonormalization(PAO_coeff, S, n_active_aos, rank)
!
      call mem%dealloc(S, n_active_aos, n_active_aos)
!
      wf%n_v = rank
!
!     Update orbital_coefficients and place orbital_coefficients_frozen_hf
!
      call wf%initialize_orbital_coefficients_frozen_hf()
!
      call mem%alloc(orbitals_copy, wf%ao%n, wf%n_mo)
!
      call dcopy(wf%ao%n*wf%n_mo, wf%orbital_coefficients, 1, orbitals_copy, 1)
!
      call wf%destruct_orbital_coefficients()
!
      call mem%alloc(wf%orbital_coefficients, wf%ao%n, wf%n_o + wf%n_v)
!
!$omp parallel do private (ao, i)
      do ao = 1, wf%ao%n
         do i = 1, wf%n_o
!
            wf%orbital_coefficients(ao, i) = orbitals_copy(ao, i)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private (ao, a)
      do ao = 1, wf%ao%n
         do a = 1, wf%n_v
!
            wf%orbital_coefficients(ao, wf%n_o + a) = PAO_coeff(ao, a)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private (ao, i)
      do ao = 1, wf%ao%n
         do i = 1, wf%n_frozen_hf_o
!
            wf%orbital_coefficients_frozen_hf(ao, i) = orbitals_copy(ao, wf%n_o + i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(orbitals_copy, wf%ao%n, wf%n_mo)
      call mem%dealloc(PAO_coeff, wf%ao%n, n_active_aos)
!
      call wf%destruct_orbital_energies()
      call mem%alloc(wf%orbital_energies, wf%n_o + wf%n_v)
!
      wf%n_frozen_hf_orbitals = wf%n_mo - wf%n_o - wf%n_v
!
      wf%n_mo = wf%n_o + wf%n_v
!
!     Diagonalize MO Fock matrix for the new
!     HF basis.
!
      call wf%diagonalize_fock_frozen_hf_orbitals()
!
      call output%printf('m', '- Preparation for frozen Hartree-Fock orbitals', &
                         fs='(/t3,a)')
!
      call output%printf('m', 'There are (i0) frozen occupied orbitals.', &
                        ints=[wf%n_frozen_hf_o], fs='(/t6,a)')
!

      call output%printf('m', 'There are (i0) frozen virtual orbitals.', &
                        ints=[wf%n_frozen_hf_orbitals-wf%n_frozen_hf_o], fs='(t6,a)')
!
!
!     Construct inactive virtuals and add to frozen CC^T. Needed for later PAO construction
!     e.g. to perform CC in reduced space or MLCC with PAOs
!
      call mem%alloc(PAO_coeff, wf%ao%n, wf%ao%n)
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
      call wf%get_full_idempotent_density(D)
!
!     Add active virtual density
!
      call dgemm('N', 'T',                               &
                  wf%ao%n,                               &
                  wf%ao%n,                               &
                  wf%n_v,                                &
                  one,                                   &
                  wf%orbital_coefficients(1,wf%n_o + 1), &
                  wf%ao%n,                               &
                  wf%orbital_coefficients(1,wf%n_o + 1), &
                  wf%ao%n,                               &
                  one,                                   &
                  D,                                     &
                  wf%ao%n)
!
!     Project occupied and active virtual out of AOs
!
      call wf%project_atomic_orbitals(D, PAO_coeff, wf%ao%n)
!
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      call mem%alloc(S, wf%ao%n, wf%ao%n)
!
      call wf%get_orbital_overlap(PAO_coeff, wf%ao%n, S)
!
!     Orthonormalize
!
      call wf%lowdin_orthonormalization(PAO_coeff, S, wf%ao%n, rank)
!
      call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
!     Add to CC^T inactive virtual density to frozen CC^T
!
      call dgemm('N', 'T',       &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  rank,          &
                  one,           &
                  PAO_coeff,     &
                  wf%ao%n,       &
                  PAO_coeff,     &
                  wf%ao%n,       &
                  one,           &
                  wf%frozen_CCT, &
                  wf%ao%n)
!
      call mem%dealloc(PAO_coeff, wf%ao%n, wf%ao%n)
!
   end subroutine remove_frozen_hf_orbitals_hf
!
!
   module subroutine diagonalize_fock_frozen_hf_orbitals_hf(wf)
!!
!!    Diagonalize Fock frozen HF orbitals
!!    Written by Sarai D. Folkestad and Linda Goletto, Nov 2019
!!
!!    Does a diagonalization of the Fock matrix in the
!!    MO basis where the frozen HF orbitals have been removed
!!
!!    Fock matrix is no longer diagonal, because determining
!!    the frozen HF orbitals entails mixing of occupied orbitals
!!    and mixing of virtual orbitals, respectively.
!!
!
      use array_utilities, only: block_diagonalize_symmetric, zero_array
!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:), allocatable  :: C_copy
      integer, dimension(2)                  :: block_dim
      integer                                :: n_blocks
!
!     We do one Roothan-Hall step to get a diagonal Fock matrix
!     (this should only entail occupied-occupied and virtual-virtual orbital mixing.)
!
      call wf%initialize_mo_fock()
      call wf%mo_transform(wf%ao_fock, wf%mo_fock)
!
!     Find active C that diagonalizes Fock in mo basis
!     also updates the orbital energies
!
      n_blocks    = 2
      block_dim   = [wf%n_o, wf%n_v]
!
      call block_diagonalize_symmetric(wf%mo_fock,  wf%n_mo, n_blocks, block_dim, wf%orbital_energies)
!
!     Transform orbitals
!
      call mem%alloc(C_copy, wf%ao%n, wf%n_mo)
      call dcopy(wf%n_mo*wf%ao%n, wf%orbital_coefficients, 1, C_copy, 1)
      call zero_array(wf%orbital_coefficients, wf%n_mo*wf%ao%n)
!
      call dgemm('N', 'N',                &
                  wf%ao%n,                &
                  wf%n_o,                 &
                  wf%n_o,                 &
                  one,                    &
                  C_copy,                 &
                  wf%ao%n,                &
                  wf%mo_fock,             &
                  wf%n_mo,                &
                  one,                    &
                  wf%orbital_coefficients,&
                  wf%ao%n)
!
      call dgemm('N', 'N',                                &
                  wf%ao%n,                                &
                  wf%n_v,                                 &
                  wf%n_v,                                 &
                  one,                                    &
                  C_copy(1, wf%n_o + 1),                  &
                  wf%ao%n,                                &
                  wf%mo_fock(wf%n_o + 1, wf%n_o + 1),     &
                  wf%n_mo,                                &
                  one,                                    &
                  wf%orbital_coefficients(1, wf%n_o + 1), &
                  wf%ao%n)
!
      call mem%dealloc(C_copy, wf%ao%n, wf%n_mo)
!
      call wf%destruct_mo_fock()
!
   end subroutine diagonalize_fock_frozen_hf_orbitals_hf
!
!
   module subroutine construct_mo_fock_fc_term_hf(wf, mo_fc_fock)
!!
!!    Calculate MO Fock frozen core contribution
!!    Written by Sarai D. Folkestad, Sep 2019
!!
!!    Constructs the frozen core contribution to
!!    the fock matrix
!!
!!       F_pq = (2 g_wxyz D^FC_yz - g_wyzx D^FC_xy) C_pw C_qx
!!
!!    in preparation of FC-CC
!!
      implicit none
!
      class(hf) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: mo_fc_fock
!
      real(dp), dimension(:,:), allocatable :: D
      real(dp), dimension(:,:), allocatable :: ao_F_fc
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
!     Add frozen core contribution
!
      call dgemm('N', 'T',                      &
                  wf%ao%n,                      &
                  wf%ao%n,                      &
                  wf%n_frozen_core_orbitals,    &
                  two,                          &
                  wf%orbital_coefficients_fc,   &
                  wf%ao%n,                      &
                  wf%orbital_coefficients_fc,   &
                  wf%ao%n,                      &
                  zero,                         &
                  D,                            &
                  wf%ao%n)
!
      call daxpy(wf%ao%n**2, half, D, 1, wf%frozen_CCT, 1)
!
!     Construct the frozen core contribution to the active Fock matrix
!
      call mem%alloc(ao_F_fc, wf%ao%n, wf%ao%n)
!
      call wf%construct_ao_G(D, ao_F_fc)
      call dscal(wf%ao%n**2, half, ao_F_fc, 1)
!
      call wf%mo_transform(ao_F_fc, mo_fc_fock)
!
      call mem%dealloc(ao_F_fc, wf%ao%n, wf%ao%n)
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      call wf%destruct_orbital_coefficients_fc()
!
   end subroutine construct_mo_fock_fc_term_hf
!
!
   module subroutine construct_mo_fock_frozen_hf_term_hf(wf, mo_frozen_hf_fock)
!!
!!    Construct MO fock frozen hf  contribution
!!    Written by Ida-Marie Høyvik, Oct 2019
!!
!!
!!    Constructs the frozen HF contribution to
!!    the fock matrix
!!
!!       F_pq = (2 g_wxyz D^F_yz - g_wyzx D^F_xy) C_pw C_qx
!!
!!    in preparation of CC in subspace.
!!
      implicit none
!
      class(hf) :: wf
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: mo_frozen_hf_fock
!
      real(dp), dimension(:,:), allocatable :: D
      real(dp), dimension(:,:), allocatable :: ao_F_frozen_hf
!
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
!     add frozen orbitals contribution
!
      call dgemm('N', 'T',                                &
                  wf%ao%n,                                &
                  wf%ao%n,                                &
                  wf%n_frozen_hf_o,                       &
                  two,                                    &
                  wf%orbital_coefficients_frozen_hf,      &
                  wf%ao%n,                                &
                  wf%orbital_coefficients_frozen_hf,      &
                  wf%ao%n,                                &
                  zero,                                   &
                  D,                                      &
                  wf%ao%n)
!
      call daxpy(wf%ao%n**2, half, D, 1, wf%frozen_CCT, 1)
!
!     Construct the frozen core contribution to the active Fock matrix
!
      call mem%alloc(ao_F_frozen_hf, wf%ao%n, wf%ao%n)
!
      call wf%construct_ao_G(D, ao_F_frozen_hf)
      call dscal(wf%ao%n**2, half, ao_F_frozen_hf, 1)
!
      call wf%mo_transform(ao_F_frozen_hf, mo_frozen_hf_fock)
!
      call mem%dealloc(ao_F_frozen_hf, wf%ao%n, wf%ao%n)
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
      call wf%destruct_orbital_coefficients_frozen_hf()
!
   end subroutine construct_mo_fock_frozen_hf_term_hf
!
!
   module function get_n_active_hf_atoms_hf(wf) result(n_active_hf_atoms)
!!
!!    Get number of active hf atoms
!!    Written by Sarai D. Folkestad and Linda Goletto, Dec 2019
!!
!!    Sets the number of active hf atoms in the system
!!
!!    NOTE: modified in mlhf
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer :: n_active_hf_atoms
!
      n_active_hf_atoms = wf%ao%get_n_centers()
!
   end function get_n_active_hf_atoms_hf
!
!
   module subroutine get_full_idempotent_density_hf(wf, D)
!!
!!    Get full idempotent density
!!    Written by Sarai D. Folkestad, Jan 2020
!!
!!    Constructs the full idempotent, occupied density
!!    for determining the frozen HF virtuals
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n,wf%ao%n), intent(out) :: D
!
      call dcopy(wf%ao%n**2, wf%ao_density, 1, D, 1)
      call dscal(wf%ao%n**2, half, D, 1)
!
   end subroutine get_full_idempotent_density_hf
!
!
   module subroutine calculate_frozen_dipole_moment_hf(wf)
!!
!!    Calculate frozen dipole moment
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Calculates the constribution to the dipole
!!    moment from the frozen orbitals
!!
      use array_utilities, only: symmetric_sandwich
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:,:), allocatable :: mu_ao
      real(dp), dimension(:,:), allocatable :: mu_mo
!
      integer :: i, j
!
      wf%frozen_dipole = zero
!
      call mem%alloc(mu_ao, wf%ao%n, wf%ao%n, 3)
!
      call wf%ao%get_oei('dipole', mu_ao)
!
      if (wf%frozen_core) then
!
         call mem%alloc(mu_mo, wf%n_frozen_core_orbitals, wf%n_frozen_core_orbitals)
!
         do i = 1, 3
            call symmetric_sandwich(mu_mo, mu_ao(:,:,i), &
                                    wf%orbital_coefficients_fc, &
                                    wf%ao%n, wf%n_frozen_core_orbitals)
!
            do j = 1, wf%n_frozen_core_orbitals
!
               wf%frozen_dipole(i) = wf%frozen_dipole(i) + two*mu_mo(j,j)
!
            enddo
!
         enddo
!
         call mem%dealloc(mu_mo, wf%n_frozen_core_orbitals, wf%n_frozen_core_orbitals)
!
      endif
!
      if (wf%frozen_hf_mos) then
!
         call mem%alloc(mu_mo, wf%n_frozen_hf_o, wf%n_frozen_hf_o)
!
         do i = 1, 3
!
            call symmetric_sandwich(mu_mo, mu_ao(:,:,i), &
                                    wf%orbital_coefficients_frozen_hf, &
                                    wf%ao%n, wf%n_frozen_hf_o)
!
            do j = 1, wf%n_frozen_hf_o
!
               wf%frozen_dipole(i) = wf%frozen_dipole(i) + two*mu_mo(j,j)
!
            enddo
!
         enddo
!
         call mem%dealloc(mu_mo, wf%n_frozen_hf_o, wf%n_frozen_hf_o)
!
      endif
!
      call mem%dealloc(mu_ao, wf%ao%n, wf%ao%n, 3)
!
   end subroutine calculate_frozen_dipole_moment_hf
!
!
   module subroutine calculate_frozen_quadrupole_moment_hf(wf)
!!
!!    Calculate frozen quadrupole moment
!!    Written by Sarai D. Folkestad, 2020
!!
!!    Calculates the constribution to the quadrupole
!!    moment from the frozen orbitals
!!
      use array_utilities, only: symmetric_sandwich
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:,:), allocatable :: q_ao
      real(dp), dimension(:,:), allocatable :: q_mo
!
      integer :: i, j
!
      wf%frozen_quadrupole = zero
!
      call mem%alloc(q_ao, wf%ao%n, wf%ao%n, 6)
!
      call wf%ao%get_oei('quadrupole', q_ao)
!
      if (wf%frozen_core) then
!
         call mem%alloc(q_mo, wf%n_frozen_core_orbitals, wf%n_frozen_core_orbitals)
!
         do i = 1, 6
!
            call symmetric_sandwich(q_mo, q_ao(:,:,i), &
                                    wf%orbital_coefficients_fc, &
                                    wf%ao%n, wf%n_frozen_core_orbitals)
!
            do j = 1, wf%n_frozen_core_orbitals
!
               wf%frozen_quadrupole(i) = wf%frozen_quadrupole(i) + two*q_mo(j,j)
!
            enddo
!
         enddo
!
         call mem%dealloc(q_mo, wf%n_frozen_core_orbitals, wf%n_frozen_core_orbitals)
!
      endif
!
      if (wf%frozen_hf_mos) then
!
         call mem%alloc(q_mo, wf%n_frozen_hf_o, wf%n_frozen_hf_o)
!
         do i = 1, 6
!
            call symmetric_sandwich(q_mo, q_ao(:,:,i), &
                                    wf%orbital_coefficients_frozen_hf, &
                                    wf%ao%n, wf%n_frozen_hf_o)
!
            do j = 1, wf%n_frozen_hf_o
!
               wf%frozen_quadrupole(i) = wf%frozen_quadrupole(i) + two*q_mo(j,j)
!
            enddo
!
         enddo
!
         call mem%dealloc(q_mo, wf%n_frozen_hf_o, wf%n_frozen_hf_o)
!
      endif
!
      call mem%dealloc(q_ao, wf%ao%n, wf%ao%n, 6)
!
   end subroutine calculate_frozen_quadrupole_moment_hf
!
!
end submodule frozen_orbital_hf
