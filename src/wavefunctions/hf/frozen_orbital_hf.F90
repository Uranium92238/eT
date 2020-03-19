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
!
      use visualization_class, only : visualization
!
      implicit none
!
      class(hf) :: wf
!
      type(visualization), allocatable :: plotter
!
      character(len=200) :: label
!
      real(dp), dimension(:,:), allocatable :: D
!
!     Destruct MO quantities in the old MO dimension, if they are allocated,
!     before n_mo changes
!
      call wf%destruct_mo_fock()
      call wf%destruct_W_mo_update()
!
!     We are done with these and want to delete them before n_mo changes
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
!
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
      if (wf%plot_active_density) then
!
         plotter = visualization(wf%system, wf%n_ao)
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
         label = 'HF_density_for_CC'
!
         call plotter%plot_density(wf%system, D, label)
!
         call mem%dealloc(D, wf%n_ao, wf%n_ao)
!
      endif
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
!!    - The number of frozen core orbitals is wf%n_frozen_core_orbitals 
!!       on exit
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
      call mem%alloc(freeze_atom, wf%get_n_active_hf_atoms())
      freeze_atom = .false.
!
!     Figure out how many core orbitals we have
!
!     Number of atoms heavier than Be (atomic number larger than 4)
!
      wf%n_frozen_core_orbitals = 0
!
      do I = 1, wf%get_n_active_hf_atoms()
!
         if (wf%system%atoms(I)%number_ .ge. 5 .and. wf%system%atoms(I)%number_ .le. 12) then
!
            wf%n_frozen_core_orbitals = wf%n_frozen_core_orbitals + 1
            freeze_atom(I) = .true.
!
         elseif (wf%system%atoms(I)%number_ .ge. 13 .and. wf%system%atoms(I)%number_ .le. 30) then
!
            wf%n_frozen_core_orbitals = wf%n_frozen_core_orbitals + 5
            freeze_atom(I) = .true.
!
         elseif (wf%system%atoms(I)%number_ .gt. 30) then
!
            call output%error_msg('No frozen core for Z > 30.')
!
         endif
!
      enddo
!
      call mem%alloc(orbital_coefficients_copy, wf%n_ao, wf%n_mo)
!
      call dcopy(wf%n_mo*wf%n_ao, wf%orbital_coefficients, 1, orbital_coefficients_copy, 1)
!
      call wf%destruct_orbital_coefficients()
!
      call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo - wf%n_frozen_core_orbitals)
!
!$omp parallel do private (ao, mo)
      do mo = 1, wf%n_mo - wf%n_frozen_core_orbitals
         do ao = 1, wf%n_ao
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
         do i = 1, wf%n_ao
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
         call get_abs_max_w_index(wf%orbital_coefficients_fc(:,i), wf%n_ao, max_, index_max)
!
         if (.not. freeze_atom(wf%system%basis2atom(index_max))) &
            call output%error_msg('Detected crossover in frozen core.')
!
      enddo
!
      call mem%dealloc(freeze_atom, wf%get_n_active_hf_atoms())
      call mem%dealloc(orbital_coefficients_copy, wf%n_ao, wf%n_mo)
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
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: D, orbitals_copy, S, PAO_coeff
!
      integer, dimension(:), allocatable :: active_aos
!
      integer :: first_ao, last_ao, i, n_active_aos, a
!
      real(dp), parameter :: full_cd_threshold = 1.0d-4
!
      integer :: mo_offset, ao, rank
!
      character(len=100) :: last_active_level
!
!     Construct active occupied orbitals
!
!     0. Determine active ao list -> all active atoms that are not hf active
!
!     First some sanity checks on the active atoms.
!     
!     Are there active atoms spaces?
!
      if (wf%system%n_active_atom_spaces .lt. 1) &
         call output%error_msg('frozen hf orbitals requested, but no active atoms on input')
!
!     If there is only one active atoms space, is it HF active?
!
      if ((wf%system%n_active_atom_spaces .eq. 1) .and. &
          (trim(wf%system%active_atom_spaces(1)%level) == 'hf')) &
            call output%error_msg('frozen hf orbitals requested, but no active CC atoms on input')
!
!     Get the first and last ao indices of the last active space.
!
      last_active_level = wf%system%active_atom_spaces(wf%system%n_active_atom_spaces)%level
!
      call wf%system%first_and_last_ao_active_space(trim(last_active_level), first_ao, last_ao)
!
!     Atoms in eT are ordered according to active space level:
!
!        CC3, CCSD, CC2, CCS, HF, inactive.
!
!     For frozen HF orbitals in CC, only the CC active spaces (CC3, CCSD, CC2, CCS) are 
!     kept. We therefore check if the last active space is HF, and in this case the
!     next to last active space is the last CC active space.
!
!     Determine if the last active space is CC active, and if not adjust the last_ao
!     first_ao is always 1
!
      if (trim(last_active_level) == 'hf') then
!
         last_ao = first_ao - 1
!
      endif
!
      first_ao = 1 ! Note: first CC active ao is always 1
!
      n_active_aos = last_ao
!
      call mem%alloc(active_aos, n_active_aos)
!
!$omp parallel do private(i)
      do i = 1, n_active_aos
!
         active_aos(i) =  i
!
      enddo
!$omp end parallel do
!
!     Occupied orbitals
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
!     2. Construct occupied active orbitals
!
      mo_offset = 0
!
      call wf%construct_orbital_block_by_density_cd(D, wf%n_o, wf%cholesky_orbital_threshold, &
                  mo_offset, active_aos)
!
      call mem%dealloc(active_aos, n_active_aos)
!
!     3. Construct occupied orbitals to be frozen 
!
      mo_offset = wf%n_o
!
      call wf%construct_orbital_block_by_density_cd(D, wf%n_frozen_hf_o, full_cd_threshold, &
                  mo_offset)
!
!     Virtual orbitals
!
      call wf%get_full_idempotent_density(D)  
!
!     1. Construct PAOs for active atoms
!
      call mem%alloc(PAO_coeff, wf%n_ao, n_active_aos)
!
      call wf%project_atomic_orbitals(D, PAO_coeff, n_active_aos, first_ao)
!
      call mem%dealloc(D, wf%n_ao, wf%n_ao)
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
      call mem%alloc(orbitals_copy, wf%n_ao, wf%n_mo)
!
      call dcopy(wf%n_ao*wf%n_mo, wf%orbital_coefficients, 1, orbitals_copy, 1)
!
      call wf%destruct_orbital_coefficients()
!
      call mem%alloc(wf%orbital_coefficients, wf%n_ao, wf%n_o + wf%n_v)
!
!$omp parallel do private (ao, i)
      do ao = 1, wf%n_ao
         do i = 1, wf%n_o
!
            wf%orbital_coefficients(ao, i) = orbitals_copy(ao, i)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private (ao, a)
      do ao = 1, wf%n_ao
         do a = 1, wf%n_v
!
            wf%orbital_coefficients(ao, wf%n_o + a) = PAO_coeff(ao, a)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do private (ao, i)
      do ao = 1, wf%n_ao
         do i = 1, wf%n_frozen_hf_o
!
            wf%orbital_coefficients_frozen_hf(ao, i) = orbitals_copy(ao, wf%n_o + i)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(orbitals_copy, wf%n_ao, wf%n_mo)
      call mem%dealloc(PAO_coeff, wf%n_ao, n_active_aos)
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
      call mem%alloc(PAO_coeff, wf%n_ao, wf%n_ao)
      call mem%alloc(D, wf%n_ao, wf%n_ao)
!
      call wf%get_full_idempotent_density(D)
!
!     Add active virtual density
!
      call dgemm('N', 'T',                               &
                  wf%n_ao,                               &
                  wf%n_ao,                               &
                  wf%n_v,                                &
                  one,                                   &
                  wf%orbital_coefficients(1,wf%n_o + 1), &
                  wf%n_ao,                               &
                  wf%orbital_coefficients(1,wf%n_o + 1), &
                  wf%n_ao,                               &
                  one,                                   &
                  D,                                     &
                  wf%n_ao)
!
!     Project occupied and active virtual out of AOs
!
      call wf%project_atomic_orbitals(D, PAO_coeff, wf%n_ao)
!
      call mem%dealloc(D, wf%n_ao, wf%n_ao)
!
      call mem%alloc(S, wf%n_ao, wf%n_ao)
!
      call wf%get_orbital_overlap(PAO_coeff, wf%n_ao, S)
!
!     Orthonormalize
!
      call wf%lowdin_orthonormalization(PAO_coeff, S, wf%n_ao, rank)
!
      call mem%dealloc(S, wf%n_ao, wf%n_ao)
!
!     Add to CC^T inactive virtual density to frozen CC^T
!
     call dgemm('N', 'T',       &
                 wf%n_ao,       &
                 wf%n_ao,       &
                 rank,          &
                 one,           &
                 PAO_coeff,     &
                 wf%n_ao,       &
                 PAO_coeff,     &
                 wf%n_ao,       &
                 one,           &
                 wf%frozen_CCT, &
                 wf%n_ao)
!
     call mem%dealloc(PAO_coeff, wf%n_ao, wf%n_ao)      
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
      implicit none
!
      class(hf) :: wf
!
!     We do one Roothan-Hall step to get a diagonal Fock matrix 
!     (this should only entail occupied-occupied and virtual-virtual orbital mixing.)
!
      call wf%initialize_mo_fock()
      call wf%mo_transform(wf%ao_fock, wf%mo_fock)
!
      call wf%initialize_W_mo_update()
!
!     Find active C that diagonalizes Fock in mo basis
!     also updates the orbital energies
!
      call wf%roothan_hall_update_orbitals_mo()
!
      call wf%destruct_mo_fock()
      call wf%destruct_W_mo_update()
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
      call mem%alloc(D, wf%n_ao, wf%n_ao)
!
!     Add frozen core contribution
!
      call dgemm('N', 'T',                      &
                  wf%n_ao,                      &
                  wf%n_ao,                      &
                  wf%n_frozen_core_orbitals,    &
                  two,                          &
                  wf%orbital_coefficients_fc,   &
                  wf%n_ao,                      &
                  wf%orbital_coefficients_fc,   &
                  wf%n_ao,                      &
                  zero,                         &
                  D,                            &
                  wf%n_ao)
!
      call daxpy(wf%n_ao**2, half, D, 1, wf%frozen_CCT, 1)
!
!     Construct the frozen core contribution to the active Fock matrix
!
      call mem%alloc(ao_F_fc, wf%n_ao, wf%n_ao)
!
      call wf%construct_ao_G(D, ao_F_fc)
!
      call wf%mo_transform(ao_F_fc, mo_fc_fock)
!
      call mem%dealloc(ao_F_fc, wf%n_ao, wf%n_ao)
      call mem%dealloc(D, wf%n_ao, wf%n_ao)
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
      call mem%alloc(D, wf%n_ao, wf%n_ao)
!
!     add frozen orbitals contribution
!
      call dgemm('N', 'T',                                &
                  wf%n_ao,                                &
                  wf%n_ao,                                &
                  wf%n_frozen_hf_o,                       &
                  two,                                    &
                  wf%orbital_coefficients_frozen_hf,      &
                  wf%n_ao,                                &
                  wf%orbital_coefficients_frozen_hf,      &
                  wf%n_ao,                                &
                  zero,                                   &
                  D,                                      &
                  wf%n_ao)
!
      call daxpy(wf%n_ao**2, half, D, 1, wf%frozen_CCT, 1)
!
!     Construct the frozen core contribution to the active Fock matrix
!
      call mem%alloc(ao_F_frozen_hf, wf%n_ao, wf%n_ao)
!
      call wf%construct_ao_G(D, ao_F_frozen_hf)
!
      call wf%mo_transform(ao_F_frozen_hf, mo_frozen_hf_fock)
!
      call mem%dealloc(ao_F_frozen_hf, wf%n_ao, wf%n_ao)
      call mem%dealloc(D, wf%n_ao, wf%n_ao)
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
      n_active_hf_atoms = wf%system%n_atoms
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
      real(dp), dimension(wf%n_ao,wf%n_ao), intent(out) :: D
!
      call dcopy(wf%n_ao**2, wf%ao_density, 1, D, 1)
      call dscal(wf%n_ao**2, half, D, 1)
!
   end subroutine get_full_idempotent_density_hf
!
!
end submodule frozen_orbital_hf
