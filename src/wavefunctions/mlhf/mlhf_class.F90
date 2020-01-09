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
module mlhf_class
!
!!
!!    Multilevel Hartree-Fock (MLHF) class module
!!    Written by Linda Goletto Ida-Marie Høyvik,
!!    and Sarai D. Folkestad, 2019
!!
!!    An initial idempotent density D, is
!!    partitioned into an active and an external part 
!!
!!       D = D^a + D^e
!!
!!    The Hartree-Fock energy expression becomes
!!
!!       E = e(D^a) + e(D^e) + 2Tr[D^aG(D^e)],
!!    
!!    where 
!!
!!       e(D^x) = 2Tr(hD^x) + Tr(D^xG(D^x)).
!!
!!    In the MLHF procedure, we only optimize the 
!!    active density (i.e., only rotate among the active orbitals)
!!
!!    The Fock matrix of the active space is
!!
!!       F = F(D^a) + G(D^e),
!!
!!    and the SCF procedure is performed considering this Fock matrix.
!!
!!    The SCF procedure is always performed in the MO basis, either
!!    with the standard Roothan-Hall procedure or with DIIS acceleration.
!!
!!    For further information, 
!!    see S. Sæther, T. Kjærgaard, H. Koch, and I-M. Høyvik, JCTC 13, no. 11 (2017),
!!    and I-M. Høyvik, Mol. Phys. (2019).
!!
!
   use hf_class
!
   use string_utilities, only : convert_to_uppercase
!
!
   implicit none
!
!  Multilevel Hartree-Fock wavefunction
!
   type, extends(hf) :: mlhf
!
      type(sequential_file) :: inactive_energy_file
      type(sequential_file) :: mlhf_inactive_fock_term_file
      type(sequential_file) :: G_De_ao_file
!
      real(dp), dimension(:,:), allocatable :: G_De
      real(dp), dimension(:,:), allocatable :: G_De_ao
!
      real(dp) :: inactive_energy
      real(dp) :: cholesky_threshold
      real(dp) :: full_space_hf_threshold
!
      logical :: cholesky_virtuals
      logical :: minimal_basis_diagonalization
      logical :: full_space_optimization
      logical :: print_initial_hf
!
   contains
!
      procedure :: prepare                                  => prepare_mlhf
      procedure :: cleanup                                  => cleanup_mlhf
      procedure :: prepare_for_roothan_hall                 => prepare_for_roothan_hall_mlhf
      procedure :: prepare_for_roothan_hall_mo              => prepare_for_roothan_hall_mlhf
      procedure :: print_energy                             => print_energy_mlhf
      procedure :: print_banner                             => print_banner_mlhf
!
      procedure :: append_orbital_info_to_restart           => append_orbital_info_to_restart_mlhf
      procedure :: is_restart_safe                          => is_restart_safe_mlhf
      procedure :: read_for_scf_restart                     => read_for_scf_restart_mlhf
      procedure :: read_for_scf_restart_mo                  => read_for_scf_restart_mlhf
!
      procedure :: roothan_hall_update_orbitals_mo          => roothan_hall_update_orbitals_mo_mlhf
      procedure :: roothan_hall_update_orbitals             => roothan_hall_update_orbitals_mo_mlhf
!
      procedure :: construct_active_mos                     => construct_active_mos_mlhf
      procedure :: set_active_mo_coefficients               => set_active_mo_coefficients_mlhf
      procedure :: construct_virtual_density_from_mo        => construct_virtual_density_from_mo_mlhf
      procedure :: construct_active_paos                    => construct_active_paos_mlhf
      procedure :: print_orbital_space_info                 => print_orbital_space_info_mlhf
!
      procedure :: determine_minimal_basis                  => determine_minimal_basis_mlhf
      procedure :: construct_idempotent_density_from_fock   => construct_idempotent_density_from_fock_mlhf
      procedure :: minimal_basis_fock_diagonalization       => minimal_basis_fock_diagonalization_mlhf
      procedure :: do_initial_full_space_optimization       => do_initial_full_space_optimization_mlhf
!
      procedure :: read_mlhf_settings                       => read_mlhf_settings_mlhf
      procedure :: read_settings                            => read_settings_mlhf
!
      procedure :: initialize_G_De                          => initialize_G_De_mlhf
      procedure :: destruct_G_De                            => destruct_G_De_mlhf
      procedure :: initialize_G_De_ao                       => initialize_G_De_ao_mlhf
      procedure :: destruct_G_De_ao                         => destruct_G_De_ao_mlhf
!
      procedure :: get_active_energy_G_De_term              => get_active_energy_G_De_term_mlhf
!
      procedure :: construct_G_De                           => construct_G_De_mlhf
!
      procedure :: get_max_roothan_hall_gradient            => get_max_roothan_hall_gradient_mlhf
!
      procedure :: update_fock_and_energy                   => update_fock_and_energy_mlhf
      procedure :: update_fock_and_energy_mo                => update_fock_and_energy_mlhf
!
      procedure :: prepare_frozen_fock_terms                => prepare_frozen_fock_terms_mlhf
      procedure :: diagonalize_fock_frozen_hf_orbitals      => diagonalize_fock_frozen_hf_orbitals_mlhf
      procedure :: get_n_active_hf_atoms                    => get_n_active_hf_atoms_mlhf
      procedure :: prepare_mos                              => prepare_mos_mlhf
!
   end type mlhf
!
!
   interface mlhf
!
      procedure :: new_mlhf
!
   end interface mlhf
!
contains
!
!
   function new_mlhf(system) result(wf)
!!
!!    New MLHF
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(mlhf) :: wf
!
      class(molecular_system), target, intent(in) :: system
!
      wf%name_ = 'mlhf'
!
      wf%system => system
!
      wf%cholesky_threshold            = 1.0d-2
      wf%full_space_hf_threshold       = 1.0d-1
      wf%full_space_optimization       = .false.
      wf%minimal_basis_diagonalization = .false.
      wf%cholesky_virtuals             = .false.
      wf%print_initial_hf              = .false.
!
      call wf%read_settings()
!
      call wf%print_banner()
!
      call wf%prepare()
!
   end function new_mlhf
!
!
   subroutine prepare_mlhf(wf)
!!
!!    Prepare
!!    Written by Linda Goletto, Sarai D. Folkestad
!!    and Eirik F. Kjønstad, 2019
!!
!!    Initializes files, writes the restart file used for consistency checks
!!    and constructs screening vectors
!!
      implicit none
!
      class(mlhf) :: wf
!
      wf%n_ao        = wf%system%get_n_aos()
      wf%n_densities = 1
!
      call wf%set_n_mo()
!
      wf%orbital_coefficients_file = sequential_file('orbital_coefficients')
      wf%orbital_energies_file = sequential_file('orbital_energies')
!
!     Initialize restart file, but for MLHF we will write
!     restart information after active space is constructed
!
      wf%restart_file = sequential_file('mlhf_restart_file')
!
!     Initialize other MLHF files
!
      wf%mlhf_inactive_fock_term_file = sequential_file('mlhf_inactive_fock_term')
      wf%G_De_ao_file = sequential_file('G_De_ao')
!
!     Construct screening vectors
!
      call wf%initialize_sp_eri_schwarz()
      call wf%initialize_sp_eri_schwarz_list()
!
      call wf%construct_sp_eri_schwarz()
!
      call wf%initialize_ao_h()
      call wf%get_ao_h_wx(wf%ao_h)
!
   end subroutine prepare_mlhf
!
!
   subroutine prepare_for_roothan_hall_mlhf(wf)
!!
!!    Prepare for Roothan-Hall
!!    Written by Linda Goletto, Ida-Marie Hoyvik
!!    and Sarai D. Folkestad, 2019
!!
!!    Performs the necessary preparations needed to solve
!!    of the Roothan-Hall equation in the iterative cycle
!!    construct Fock - calculate energy - Roothan-Hall-update orbitals -
!!    update the density:
!!
!!    - constructs the ao Fock matrix and performs a Roothan-Hall step
!!      to get the initial idempotent density;
!!    - prints the number of electrons and the energy of the initial guess;
!!    - partitions virtual (with Cholesky or pao) and occupied (with Cholesky)
!!      orbitals into active and inactive;
!!    - allocates active mo specific arrays and constructs them;
!!    - constructs the active density.
!!
!
      use array_utilities, only : identity_array
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp) :: n_electrons
!
      type(timings) :: timer
!
!     Construct AO Fock from SAD density
!
      timer = timings('AO Fock construction', pl='normal')
      call timer%turn_on()
!
!     AO fock construction and energy calculation
!
!     Construct the two electron part of the Fock matrix (G),
!     and add the contribution to the Fock matrix
!
      call wf%construct_ao_G(wf%ao_density, wf%ao_fock)
!
!     Add the one-electron part
!
      call daxpy(wf%n_ao**2, one, wf%ao_h, 1, wf%ao_fock, 1)
!
      call timer%turn_off()
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, wf%ao_h)
!
      call wf%get_n_electrons_in_density(n_electrons)
!
      call output%printf('m', 'Energy of initial guess:      (f25.12)', &
                         reals=[wf%energy], fs='(/t6, a)')
!
      call output%printf('m', 'Number of electrons in guess: (f25.12)', &
                         reals=[n_electrons], fs='(t6, a)')
!
!     Diagonalize to generate idempotent density
!
      call wf%construct_idempotent_density_from_fock()
!
!     Construct active mos
!
      call wf%construct_active_mos()
!
!     Finished with partitioning of density, print info to output
!
      call wf%print_orbital_space_info()
!
!     Allocate active mo specific arrays
!     and construct them
!
      call wf%initialize_W_mo_update()
      call identity_array(wf%w_mo_update, wf%n_mo)
!
      call wf%initialize_G_De_ao()
      call wf%initialize_G_De()
      call wf%initialize_mo_fock()
!
      call wf%construct_G_De()
!
!     Inactive energy contribution
!
      wf%inactive_energy = wf%calculate_hf_energy_from_G(wf%G_De_ao, wf%ao_h)
!
!     Print multilevel orbital information to restart file
!
      call wf%append_orbital_info_to_restart()
!
!     Construct active ao density
!
      call wf%construct_ao_density()
!
      call wf%update_fock_and_energy_mo()
!
      call wf%roothan_hall_update_orbitals_mo()  ! DIIS F => C
      call wf%update_ao_density()
!
   end subroutine prepare_for_roothan_hall_mlhf
!
!
   subroutine append_orbital_info_to_restart_mlhf(wf)
!!
!!    Append orbital info to restart
!!    Written by Linda Goletto and Anders Hutcheson, 2019
!!
      implicit none
!
      class(mlhf) :: wf
!
      integer :: i, n_active_atoms, n_active_spaces
      character(len=100) :: last_active_space_level
!
      n_active_spaces = wf%system%n_active_atom_spaces
      last_active_space_level = wf%system%active_atom_spaces(n_active_spaces)%level
!
      if (trim(last_active_space_level) .ne. 'hf') then
         call output%error_msg('No hf active space found')
      endif
!
      n_active_atoms = wf%system%active_atom_spaces(n_active_spaces)%last_atom
!
      call wf%restart_file%open_('write', 'append')
!
      call wf%restart_file%write_(n_active_atoms)
!
      do i=1, n_active_atoms
         call wf%restart_file%write_(wf%system%atoms(i)%input_number)
      enddo
!
      call wf%restart_file%write_(wf%n_o)
      call wf%restart_file%write_(wf%n_v)
      call wf%restart_file%write_(wf%inactive_energy)
!
      call wf%restart_file%close_
!
   end subroutine append_orbital_info_to_restart_mlhf
!
!
   subroutine roothan_hall_update_orbitals_mo_mlhf(wf)
!!
!!    Roothan-Hall update of orbitals
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Update orbitals after a Roothan-Hall step in the MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      call wf%hf%roothan_hall_update_orbitals_mo()
!
   end subroutine roothan_hall_update_orbitals_mo_mlhf
!
!
   subroutine update_fock_and_energy_mlhf(wf, prev_ao_density)
!!
!!    Update Fock and energy
!!    Written by Linda Goletto and Sarai D. Folkestad, 2019
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock, energy and basis for G_De are to be computed.
!!
!!    Transforms G_De to the new MO basis
!!
!!       G_De_new = W^T G_De_old W
!!
!!    where W is MO_basis_update.
!!
!!    Constructs Fock matrix (active) in the AO basis
!!    and calculates the energy.
!!
!!    Transforms the active Fock matrix to the new MO basis.
!!    and adds the inactive part of the Fock matrix (G_De_new).
!!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in), optional :: prev_ao_density
!
      real(dp), dimension(:,:), allocatable :: Z_pq ! = sum_x G_De_old_wx * w_xq
      real(dp), dimension(:,:), allocatable :: G_De_old
!
!     Update of the MO basis for the G_De
!     G_De =  w^T * G_De_old * w
!
      call mem%alloc(Z_pq, wf%n_mo, wf%n_mo)
      call mem%alloc(G_De_old, wf%n_mo, wf%n_mo)
!
      call dcopy(wf%n_mo**2, wf%G_De, 1, G_De_old, 1)
!
      call dgemm('N', 'N',          &
                  wf%n_mo,          &
                  wf%n_mo,          &
                  wf%n_mo,          &
                  one,              &
                  G_De_old,         &
                  wf%n_mo,          &
                  wf%w_mo_update,   &
                  wf%n_mo,          &
                  zero,             &
                  Z_pq,             &
                  wf%n_mo)
!
      call dgemm('T', 'N',          &
                  wf%n_mo,          &
                  wf%n_mo,          &
                  wf%n_mo,          &
                  one,              &
                  wf%w_mo_update,   &
                  wf%n_mo,          &
                  Z_pq,             &
                  wf%n_mo,          &
                  zero,             &
                  wf%G_De,          &
                  wf%n_mo)
!
      call mem%dealloc(G_De_old, wf%n_mo, wf%n_mo)
      call mem%dealloc(Z_pq, wf%n_mo, wf%n_mo)
!
!     Write the new G(De) to file
!
      call wf%mlhf_inactive_fock_term_file%open_('write', 'rewind')
      call wf%mlhf_inactive_fock_term_file%write_(wf%G_De, wf%n_mo**2)
      call wf%mlhf_inactive_fock_term_file%close_
!
!     AO fock construction and energy calculation
!
      call wf%hf%update_fock_and_energy_mo(prev_ao_density)
!
!     Add the Tr[Da * G(De)] and inactive energy contributions to the energy
!
      wf%energy = wf%energy + wf%inactive_energy + wf%get_active_energy_G_De_term()
!
!     Add G_De to MO fock
!
      call daxpy(wf%n_mo**2, one, wf%G_De, 1, wf%mo_fock, 1)
!
   end subroutine update_fock_and_energy_mlhf
!
!
   subroutine set_active_mo_coefficients_mlhf(wf,C_o,C_v)
!!
!!    Set active MO coefficients
!!    Written by Ida-Marie Høyvik and Linda Goletto, 2019
!!
!!    Places the coefficients in the input
!!    into the coefficient matrix:
!!
!!    - takes the first n_o vectors from C_o (occupied coefficients)
!!      and places into coefficient matrix;
!!    - takes the first n_v vectors from C_v (virtual coefficients)
!!      and places into coefficient matrix.
!!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_ao,wf%n_o), intent(in) :: C_o
      real(dp), dimension(wf%n_ao,wf%n_v), intent(in) :: C_v
!
      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
!
!     new dimension of mo space
      wf%n_mo = wf%n_o + wf%n_v
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
      call wf%initialize_pivot_matrix_ao_overlap()
      call wf%initialize_cholesky_ao_overlap()
!
      call dcopy(wf%n_ao*wf%n_o, C_o, 1, wf%orbital_coefficients, 1)
!
      call dcopy(wf%n_ao*wf%n_v, C_v, 1, wf%orbital_coefficients(1,wf%n_o+1), 1)
!
   end subroutine set_active_mo_coefficients_mlhf
!
!
   subroutine construct_virtual_density_from_mo_mlhf(wf, D_v)
!!
!!    Construct virtual density from MO
!!    Written by Sarai D. Folkestad and Ida-Marie Hoyvik 2018
!!
!!    Construct the virtual density either as
!!
!!       D^V = sum_a C_alpha,a * C_beta,a
!!
!!    or as
!!
!!       D^V = S^inv - D,
!!
!!    if no virtual mo coefficients have been generated.
!!
!
      use array_utilities, only : invert
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out) :: D_v
!
      if (wf%minimal_basis_diagonalization) then
!
         call invert(D_v, wf%ao_overlap, wf%n_ao)
         call daxpy(wf%n_ao**2, -one, wf%ao_density, 1, D_v, 1)
!
      else
!
         call dgemm('N', 'T',                              &
                   wf%n_ao,                                &
                   wf%n_ao,                                &
                   wf%n_v,                                 &
                   one,                                    &
                   wf%orbital_coefficients(1, wf%n_o + 1), &
                   wf%n_ao,                                &
                   wf%orbital_coefficients(1, wf%n_o + 1), &
                   wf%n_ao,                                &
                   zero,                                   &
                   D_v,                                    &
                   wf%n_ao)
!
      endif
!
   end subroutine construct_virtual_density_from_mo_mlhf
!
!
   subroutine construct_idempotent_density_from_fock_mlhf(wf)
!!
!!    Construct idempotent density from Fock
!!    Written by Ida-Marie Høyvik and Sarai D. Folkestad,  Oct 2019.
!!
!!    Constructs the idempotent density.
!!
!!    This is done by either using the standard diagonalization
!!    of the hf class (do_roothan_hall) or by projecting the problem
!!    onto a minimal basis.
!!
!!    Modified by Anders Hutcheson and Linda Goletto, Oct 2019
!!
!!    Added option of an initial optimization in the full space.
!!
!
      implicit none
!
      class(mlhf) :: wf
!
      if (wf%minimal_basis_diagonalization) then
!
!        The problem is projected on a minimal basis
!
         call wf%minimal_basis_fock_diagonalization()
!
      elseif (wf%full_space_optimization) then
!
!        An initial optimization of the density in the full space is performed
!
         call wf%do_initial_full_space_optimization
!
      else
!
!        Standard Roothan-Hall step to generate idempotent density
!
         call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies) ! F => C
         call wf%update_ao_density() ! C => D
!
      endif
!
   end subroutine construct_idempotent_density_from_fock_mlhf
!
!
   subroutine minimal_basis_fock_diagonalization_mlhf(wf)
!!
!!    Minimal basis fock diagonalization
!!    Written by Ida-Marie Høyvik and Sarai D. Folkestad,  Oct 2019.
!!
!!    Projects the problem on a minimal basis and
!!    performs the diagonalization to construct the density
!!
!!    The minimal basis is determined by considering the
!!    magnitude of the diagonal elements of the AO density.
!!
!
      use array_utilities, only : invert, generalized_diagonalization, &
                                  count_n_true, extract_columns_of_matrix
!
      implicit none
!
      class(mlhf) :: wf
!
      logical, dimension(:), allocatable :: significant_AOs
!
      real(dp), dimension(:,:), allocatable :: Q_tr, Q
      real(dp), dimension(:,:), allocatable :: S_minimal, S_minimal_inv
      real(dp), dimension(:,:), allocatable :: T
      real(dp), dimension(:,:), allocatable :: X
      real(dp), dimension(:,:), allocatable :: F_minimal, C_minimal, D_minimal
!
      integer :: I, J
      integer :: n_significant_AOs
!
      call mem%alloc(significant_AOs, wf%n_ao)
!
      call wf%determine_minimal_basis(significant_AOs)
!
      n_significant_AOs = count_n_true(wf%n_ao, significant_AOs)
!
      call mem%alloc(Q_tr, wf%n_ao, n_significant_AOs)
      call extract_columns_of_matrix(n_significant_AOs, wf%n_ao, wf%n_ao, &
                                     wf%ao_overlap, Q_tr, significant_AOs)
!
      call mem%alloc(Q, n_significant_AOs, wf%n_ao)
!
!$omp parallel do private(I, J)
      do I = 1, n_significant_AOs
         do J = 1, wf%n_ao
!
            Q (I,J) = Q_tr(J,I)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(Q_tr, wf%n_ao, n_significant_AOs)
!
      call mem%alloc(S_minimal,  n_significant_AOs,  n_significant_AOs)
      call extract_columns_of_matrix(n_significant_AOs, n_significant_AOs, wf%n_ao, &
                                       Q, S_minimal, significant_AOs)
!
      call mem%dealloc(significant_AOs, wf%n_ao)
!
      call mem%alloc(S_minimal_inv,  n_significant_AOs,  n_significant_AOs)
!
      call invert(S_minimal_inv, S_minimal, n_significant_AOs)
!
      call mem%dealloc(S_minimal,  n_significant_AOs,  n_significant_AOs)
!
      call mem%alloc(T, n_significant_AOs, wf%n_ao)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%n_ao,             &
                  n_significant_AOs,   &
                  one,                 &
                  S_minimal_inv,       &
                  n_significant_AOs,   &
                  Q,                   &
                  n_significant_AOs,   &
                  zero,                &
                  T,                   &
                  n_significant_AOs)
!
      call mem%dealloc(Q, n_significant_AOs, wf%n_ao)
      call mem%dealloc(S_minimal_inv,  n_significant_AOs,  n_significant_AOs)
!
!     Project F down on the minimal basis
!
      call mem%alloc(X, n_significant_AOs, wf%n_ao)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%n_ao,             &
                  wf%n_ao,             &
                  one,                 &
                  T,                   &
                  n_significant_AOs,   &
                  wf%ao_fock,          &
                  wf%n_ao,             &
                  zero,                &
                  X,                   &
                  n_significant_AOs)
!
      call mem%alloc(F_minimal, n_significant_AOs, n_significant_AOs)
!
      call dgemm('N', 'T',             &
                  n_significant_AOs,   &
                  n_significant_AOs,   &
                  wf%n_ao,             &
                  one,                 &
                  X,                   &
                  n_significant_AOs,   &
                  T,                   &
                  n_significant_AOs,   &
                  zero,                &
                  F_minimal,           &
                  n_significant_AOs)
!
      call mem%dealloc(X, n_significant_AOs, wf%n_ao)
!
!     Project S on the minimal basis
!
      call mem%alloc(X, n_significant_AOs, wf%n_ao)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%n_ao,             &
                  wf%n_ao,             &
                  one,                 &
                  T,                   &
                  n_significant_AOs,   &
                  wf%ao_overlap,       &
                  wf%n_ao,             &
                  zero,                &
                  X,                   &
                  n_significant_AOs)
!
      call mem%alloc(S_minimal, n_significant_AOs, n_significant_AOs)
!
      call dgemm('N', 'T',             &
                  n_significant_AOs,   &
                  n_significant_AOs,   &
                  wf%n_ao,             &
                  one,                 &
                  X,                   &
                  n_significant_AOs,   &
                  T,                   &
                  n_significant_AOs,   &
                  zero,                &
                  S_minimal,           &
                  n_significant_AOs)
!
      call mem%dealloc(X, n_significant_AOs, wf%n_ao)
!
      call mem%alloc(C_minimal, n_significant_AOs, n_significant_AOs)
!
!     Diagonalize Fock in the minimal basis
!
      call generalized_diagonalization(n_significant_AOs, F_minimal, S_minimal, C_minimal) ! F => C
!
      call mem%dealloc(S_minimal, n_significant_AOs, n_significant_AOs)
      call mem%dealloc(F_minimal, n_significant_AOs, n_significant_AOs)
!
      call mem%alloc(D_minimal, n_significant_AOs, n_significant_AOs)
!
!     NOTE: Factor two to get the non-idempotent density
!
      call dgemm('N', 'T',             &
                  n_significant_AOs,   &
                  n_significant_AOs,   &
                  wf%n_o,              &
                  two,                 &
                  C_minimal,           &
                  n_significant_AOs,   &
                  C_minimal,           &
                  n_significant_AOs,   &
                  zero,                &
                  D_minimal,           &
                  n_significant_AOs)
!
      call mem%dealloc(C_minimal, n_significant_AOs, n_significant_AOs)
!
      call mem%alloc(X, n_significant_AOs, wf%n_ao)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%n_ao,             &
                  n_significant_AOs,   &
                  one,                 &
                  D_minimal,           &
                  n_significant_AOs,   &
                  T,                   &
                  n_significant_AOs,   &
                  zero,                &
                  X,                   &
                  n_significant_AOs)
!
      call mem%dealloc(D_minimal, n_significant_AOs, n_significant_AOs)
!
      call dgemm('T', 'N',             &
                  wf%n_ao,             &
                  wf%n_ao,             &
                  n_significant_AOs,   &
                  one,                 &
                  T,                   &
                  n_significant_AOs,   &
                  X,                   &
                  n_significant_AOs,   &
                  zero,                &
                  wf%ao_density,       &
                  wf%n_ao)
!
      call mem%dealloc(T, n_significant_AOs, wf%n_ao)
      call mem%dealloc(X, n_significant_AOs, wf%n_ao)
!
   end subroutine minimal_basis_fock_diagonalization_mlhf
!
!
   subroutine determine_minimal_basis_mlhf(wf, significant_AOs)
!!
!!    Determine minimal basis
!!    Written by Ida-Marie Høyvik and Sarai D. Folkestad, Oct 2019
!!
!!    Determines the minimal basis (significant AOs).
!!    An AO alpha is significant if the corresponding diagonal
!!    element of the AO density is above a given threshold:
!!
!!       D_alpha,alpha > occupation_treshold
!!
!!    The output significant_AOs is a logical array where
!!    significant_AOs(I)==.true. if the Ith AO is in the minimal basis.
!!
      implicit none
!
      class(mlhf), intent(in) :: wf
!
      logical, dimension(wf%n_ao), intent(out) :: significant_AOs

      real(dp), parameter :: occupation_treshold = 0.1d0
!
      integer :: ao
!
      significant_AOs = .false.
!
!$omp parallel do private(ao)
      do ao = 1, wf%n_ao
!
         if (wf%ao_density(ao, ao) .gt. occupation_treshold) significant_AOs(ao) = .true.
!
      enddo
!$omp end parallel do
!
   end subroutine determine_minimal_basis_mlhf
!
!
   subroutine read_settings_mlhf(wf)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(mlhf) :: wf
!
      call wf%read_hf_settings()
      call wf%read_frozen_orbitals_settings()
      call wf%read_mlhf_settings()
!
   end subroutine read_settings_mlhf
!
!
   subroutine read_mlhf_settings_mlhf(wf)
!!
!!    Read mlhf settings
!!    Written by Sarai D. Folkestad, Oct 2019
!!
      implicit none
!
      class(mlhf) :: wf               
!
      call input%get_keyword_in_section('cholesky threshold', 'multilevel hf', wf%cholesky_threshold)
      call input%get_keyword_in_section('initial hf threshold', 'multilevel hf', wf%full_space_hf_threshold)
!
      if (input%requested_keyword_in_section('project on minimal basis', 'multilevel hf')) &
            wf%minimal_basis_diagonalization = .true.
!
      if (input%requested_keyword_in_section('cholesky virtuals', 'multilevel hf')) &
            wf%cholesky_virtuals = .true.
!
      if (input%requested_keyword_in_section('initial hf optimization', 'multilevel hf')) &
            wf%full_space_optimization = .true.
!
      if (input%requested_keyword_in_section('print initial hf', 'multilevel hf')) &
            wf%print_initial_hf = .true.
!
!     Sanity checks
!
      if (input%requested_keyword_in_section('initial hf threshold', 'multilevel hf') &
            .and. .not. wf%full_space_optimization) &
            call output%error_msg('Initial hf threshold specified, without specifying initial hf optimization')
!
      if (wf%full_space_optimization .and. wf%minimal_basis_diagonalization) &
            call output%error_msg('Projection on minimal basis and initial hf optimization are not compatible')
!
   end subroutine read_mlhf_settings_mlhf
!
!
   subroutine construct_active_paos_mlhf(wf, C_pao)
!!
!!    Construct active PAOs
!!    Written by Ida-Marie Hoyvik and Sarai D. Folkestad 2019
!!
!!    Constructs the active PAOs and orthogonalizes them,
!!    on exit the first n_v vectors of C_paos are active PAOs
!!    and wf%n_v is updated.
!!
!
      use array_utilities
!
      implicit none
!
      class(mlhf)   :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(out)  :: C_pao
!
      real(dp), dimension(:,:), allocatable :: C_pao_copy
      real(dp), dimension(:,:), allocatable :: S_pao
!
      integer   :: n_active_aos
      integer   :: last_hf_ao
!
      call wf%system%last_ao_active_space('hf', last_hf_ao)
!
      n_active_aos = last_hf_ao ! Because if there are CC active atoms, these are also HF active
!
      call mem%alloc(C_pao_copy, wf%n_ao, n_active_aos)
      call wf%projected_atomic_orbitals(wf%ao_density, C_pao_copy, n_active_aos)
!
      call mem%alloc(S_pao, n_active_aos, n_active_aos)
!
      call wf%get_orbital_overlap(C_pao_copy, n_active_aos, S_pao)
!
      call wf%lowdin_orthonormalization(C_pao_copy, S_pao, n_active_aos, wf%n_v)
!
      call mem%dealloc(S_pao, n_active_aos, n_active_aos)
!
      call dcopy(wf%n_ao*wf%n_v, C_pao_copy, 1, C_pao, 1)
!
      call mem%dealloc(C_pao_copy, wf%n_ao, n_active_aos)
!
   end subroutine construct_active_paos_mlhf
!
!
   real(dp) function get_active_energy_G_De_term_mlhf(wf)
!!
!!    Get active energy contribution from G(De)
!!    Written by Sarai D. Folkestad and Linda Goletto, 2019
!!
!!    Returns the contribution to the active energy from G(D_e)
!!
!!       Tr(D_a * G(D_e))
!!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp) :: E_G_De
!
      integer :: i
!
      E_G_De = zero
!
      do i = 1, wf%n_o
!
         E_G_De = E_G_De + wf%G_De(i,i)
!
      end do
!
      E_G_De = two * E_G_De
!
      get_active_energy_G_De_term_mlhf = E_G_De
!
   end function get_active_energy_G_De_term_mlhf
!
!
   subroutine print_orbital_space_info_mlhf(wf)
!!
!!    Print orbital space information
!!    Written by Linda Goletto and Ida-Marie Hoyvik, 2019
!!
!!    Prints the information on the active orbitals
!!
      implicit none
!
      class(mlhf) :: wf
!
     call output%printf('m', '- Active orbital space:', &
            fs='(/t3, a)')
     call output%printf('m', ' Number of active occupied orbitals: (i8)', &
            ints=[wf%n_o], ffs='(/t6, a)')
     call output%printf('m', ' Number of active virtual orbitals:  (i8)', &
            ints=[wf%n_v], ffs='(t6, a)')
     call output%printf('m', ' Number of active orbitals:          (i8)', &
            ints=[wf%n_v+wf%n_o], ffs='(t6, a)')
!
   end subroutine print_orbital_space_info_mlhf
!
!
   subroutine initialize_G_De_mlhf(wf)
!!
!!    Initialize G(De)
!!    Written by Ida-Marie Hoyvik, Oct 2019
!!
!!    Initializes G_De in the MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (.not. allocated(wf%G_De)) call mem%alloc(wf%G_De, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_G_De_mlhf
!
!
   subroutine initialize_G_De_ao_mlhf(wf)
!!
!!    Initialize G(De)_wx
!!    Written by Ida-Marie Hoyvik and Linda Goletto, 2019
!!
!!    Initializes G_De in the AO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      call mem%alloc(wf%G_De_ao, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_G_De_ao_mlhf
!
!
   subroutine destruct_G_De_mlhf(wf)
!!
!!    Destruct G(De)
!!    Written by Ida-Marie Hoyvik, Oct 2019
!!
!!    Destructs G_De in the MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (allocated(wf%G_De)) call mem%dealloc(wf%G_De, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_G_De_mlhf
!
!
   subroutine destruct_G_De_ao_mlhf(wf)
!!
!!    Destruct G(De)_wx
!!    Written by Ida-Marie Hoyvik and Linda Goletto, 2019
!!
!!    Destructs G_De in the AO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      call mem%dealloc(wf%G_De_ao, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_G_De_ao_mlhf
!
!
   subroutine cleanup_mlhf(wf)
!!
!!    Cleanup
!!    Written by Linda Goletto, Sarai D. Folkestad
!!    and Ida-Marie Hoyvik, Oct 2019
!!
      implicit none
!
      class(mlhf) :: wf
! 
      call wf%save_ao_density()
!
      call wf%destruct_orbital_energies()
      call wf%destruct_orbital_coefficients()
      call wf%destruct_ao_overlap()
      call wf%destruct_ao_fock()
      call wf%destruct_mo_fock()
      call wf%destruct_ao_density()
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
!
      call wf%destruct_W_mo_update()
      call wf%destruct_G_De()
      call wf%destruct_G_De_ao()
!
      call wf%destruct_sp_eri_schwarz()
      call wf%destruct_sp_eri_schwarz_list()
      call wf%destruct_mo_fock_frozen()
!
      call wf%G_De_ao_file%delete_
!
   end subroutine cleanup_mlhf
!
!
   subroutine construct_G_De_mlhf(wf)
!!
!!    Construct G(De)
!!    Written by Linda Goletto, Ida-Marie Høyvik
!!    and Sarai D. Folkestad, 2019
!!
!!    Constructs G(De) in the MO basis
!!
      implicit none
!
      class(mlhf), intent(inout) :: wf
!
      type(timings) :: G_De_construction_timer
!
      G_De_construction_timer = timings('G(De) construction', pl='normal')
      call G_De_construction_timer%turn_on()
!
!     Scale by two to get non-idempotent inactive density De
!
      call dscal(wf%n_ao**2, two, wf%ao_density, 1)
!
      call wf%construct_ao_G(wf%ao_density, wf%G_De_ao)
!
      call wf%G_De_ao_file%open_('write', 'rewind')
      call wf%G_De_ao_file%write_(wf%G_De_ao, wf%n_ao**2)
      call wf%G_De_ao_file%close_
!
      call wf%mo_transform(wf%G_De_ao, wf%G_De)
!
      call wf%mlhf_inactive_fock_term_file%open_('write', 'rewind')
      call wf%mlhf_inactive_fock_term_file%write_(wf%G_De, wf%n_mo**2)
      call wf%mlhf_inactive_fock_term_file%close_
!
      call G_De_construction_timer%turn_off()
!
   end subroutine construct_G_De_mlhf
!
!
   function get_max_roothan_hall_gradient_mlhf(wf) result(max_gradient)
!!
!!    Get max of Roothan-Hall gradient
!!    Written by Sarai D. Folkestad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       E^(1) = - 4 F_ai
!!
!!    and returns the maximum absolute value of E^(1).
!!
      implicit none
!
      class(mlhf), intent(in) :: wf
!
      real(dp) :: max_gradient
!
      real(dp), dimension(:,:), allocatable :: gradient
!
      call mem%alloc(gradient, wf%n_v, wf%n_o)
!
      call wf%hf%get_roothan_hall_mo_gradient(gradient)
!
      max_gradient = get_abs_max(gradient, wf%n_o*wf%n_v)
!
      call mem%dealloc(gradient, wf%n_v, wf%n_o)
!
   end function get_max_roothan_hall_gradient_mlhf
!
!
   subroutine print_energy_mlhf(wf)
!!
!!    Print wavefunction summary
!!    Written by Eirik F. Kjønstad, Sep 2018
!!    Modified by Tommaso Giovannini, March 2019
!!
!!    Prints information related to the wavefunction, most of which is meaningful
!!    only for a properly converged wavefunction.
!!
      implicit none
!
      class(mlhf), intent(inout) :: wf
!
      real(dp) :: nuclear_repulsion, E_G_De
!
      call wf%hf%print_energy()
!
      nuclear_repulsion = wf%system%get_total_nuclear_repulsion()
      E_G_De = wf%get_active_energy_G_De_term()
!
      call output%printf('m', '- Summary of '// &
                         trim(convert_to_uppercase(wf%name_))// &
                         ' active/inactive contributions to electronic energy (a.u.):', &
                         ll=80, fs='(/t3,a)')
!
      call output%printf('m', 'Active energy:             (f19.12)', &
                         reals=[wf%energy - (E_G_De + wf%inactive_energy + nuclear_repulsion)], fs='(/t6,a)')
      call output%printf('m', 'Active-inactive energy:    (f19.12)', &
                         reals=[E_G_De], fs='(t6,a)')
      call output%printf('m', 'Inactive energy:           (f19.12)', &
                         reals=[wf%inactive_energy], fs='(t6,a)')
!
   end subroutine print_energy_mlhf
!
!
   subroutine construct_active_mos_mlhf(wf)
!!
!!    Construct active MOs
!!    Written by Linda Goletto, Ida-Marie Høyvik
!!    and Sarai D. Folkestad, Oct 2019
!!
!!    Partitions virtual (with cholesky or pao) and occupied (with cholesky)
!!    orbitals into active and inactive.
!!
!!    On exit we have updated the wf%n_o, wf%n_v, wf%n_mo, and the
!!    orbital coefficients
!!
!
      use array_utilities, only : cholesky_decomposition_limited_diagonal
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: D_v, C_v, C_o
!
      integer     :: n_vectors
!
      integer     :: n_active_aos
      integer     :: i
      integer     :: last_hf_ao
!
      integer,  dimension(:),   allocatable :: active_aos
!
      type(timings) :: active_mos_construction_timer
!
      active_mos_construction_timer = timings('Active MOs construction', pl='normal')
      call active_mos_construction_timer%turn_on()
!
      n_active_aos = 0
      n_vectors    = 0
!
      call wf%system%last_ao_active_space('hf', last_hf_ao)
!
      n_active_aos =  last_hf_ao ! Because if there are CC active atoms, these are also HF active
!
      call mem%alloc(active_aos, n_active_aos)
!
      do i = 1, n_active_aos
!
         active_aos(i) = i
!
      enddo
!
!     Partition virtual space
!
      call mem%alloc(C_v, wf%n_ao, wf%n_ao)
!
      call dscal(wf%n_ao**2, half, wf%ao_density, 1) ! Idempotent ao density
!
      if (wf%cholesky_virtuals) then
!
!        Cholesky orbitals
!
         call mem%alloc(D_v, wf%n_ao, wf%n_ao)
!
!        Make virtual density before occupied density is changed
!
         call wf%construct_virtual_density_from_mo(D_v)

         call cholesky_decomposition_limited_diagonal(D_v, C_v, wf%n_ao, wf%n_v, &
                                       wf%cholesky_threshold, n_active_aos, active_aos)
!
         call mem%dealloc(D_v, wf%n_ao, wf%n_ao)
!
      else
!
!        PAOs
!
         call wf%construct_active_paos(C_v)
!
      endif
!
!     Partition occupied space (Cholesky orbitals)
!
      call mem%alloc(C_o, wf%n_ao, wf%n_ao)
!
      call cholesky_decomposition_limited_diagonal(wf%ao_density, C_o, wf%n_ao, wf%n_o, &
                                  wf%cholesky_threshold, n_active_aos, active_aos)
!
!     Set molecular orbitals coefficients
!
      call wf%set_active_mo_coefficients(C_o(:,1:wf%n_o), C_v(:,1:wf%n_v))
!
      call mem%dealloc(C_v, wf%n_ao, wf%n_ao)
      call mem%dealloc(C_o, wf%n_ao, wf%n_ao)
!
      call mem%dealloc(active_aos, n_active_aos)
!
      call active_mos_construction_timer%turn_off()
!
   end subroutine construct_active_mos_mlhf
!
!
   subroutine is_restart_safe_mlhf(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Linda Goletto
!!    and Anders Hutcheson, Mar 2019
!!
      implicit none
!
      class(mlhf) :: wf
!
      character(len=*), intent(in) :: task
!
      integer :: n_ao, n_densities, i, n_active_atoms, n_active_spaces, n_electrons
      integer, dimension(:), allocatable :: active_atoms
      character(len=100) :: last_active_space_level
!
      n_active_spaces = wf%system%n_active_atom_spaces
      last_active_space_level = wf%system%active_atom_spaces(n_active_spaces)%level
!
      if (last_active_space_level .ne. 'hf') then
         call output%error_msg('The last active space is not hf, but:' // last_active_space_level)
      endif
!
      n_active_atoms = wf%system%active_atom_spaces(n_active_spaces)%last_atom
!
      call wf%restart_file%open_('read', 'rewind')
!
      call wf%restart_file%read_(n_ao)
      call wf%restart_file%read_(n_densities)
      call wf%restart_file%read_(n_electrons)
      call wf%restart_file%read_(n_active_atoms)
!
      if (n_ao .ne. wf%n_ao) then
         call output%error_msg('attempted to restart MLHF ' // &
                           'with an inconsistent number of atomic orbitals for task ' // trim(task))
      endif
!
      if (n_densities .ne. wf%n_densities) then
         call output%error_msg('attempted to restart MLHF with an inconsistent number ' // &
            'of atomic densities (likely a HF/UHF inconsistency) for task ' // trim(task))
      endif
!
      if (n_electrons .ne. wf%system%get_n_electrons()) then
         call output%error_msg('attempted to restart MLHF with an inconsistent number ' // &
                               'of electrons for task ' // trim(task))
      endif
!
      if (n_active_atoms .ne. wf%system%active_atom_spaces(n_active_spaces)%last_atom) then
         call output%error_msg('attempted to restart MLHF with an inconsistent number ' // &
            'of active atoms for task ' // trim(task))
      endif
!
      call mem%alloc(active_atoms, n_active_atoms)
!
      do i = 1, n_active_atoms
         call wf%restart_file%read_(active_atoms(i))
         if (active_atoms(i) .ne. wf%system%atoms(i)%input_number) then
                     call output%error_msg('attempted to restart MLHF ' // &
                     'with inconsistent active atoms for task ' // trim(task))
         endif
      enddo
!
      call mem%dealloc(active_atoms, n_active_atoms)
!
      call wf%restart_file%close_
!
   end subroutine is_restart_safe_mlhf
!
!
   subroutine read_for_scf_restart_mlhf(wf)
!!
!!    Read for SCF restart
!!    Written by Sarai D. Folkestad, Anders Hutcheson
!!    and Linda Goletto, Oct 2019
!!
      implicit none
!
      class(mlhf) :: wf
!
      integer :: n_active_atoms
!
!     Destruct orbital coeffiecients and orbital energies allocated
!     with the old n_mo
!
      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
!
!     Read n_o, n_v and inactive_energy
!
      call wf%restart_file%open_('read', 'rewind')
!
      call wf%restart_file%skip(3)
      call wf%restart_file%read_(n_active_atoms)
      call wf%restart_file%skip(n_active_atoms)
      call wf%restart_file%read_(wf%n_o)
      call wf%restart_file%read_(wf%n_v)
      call wf%restart_file%read_(wf%inactive_energy)
!
      call wf%restart_file%close_
!
!     Calculate the new n_mo and initialize orbital coefficients and
!     orbital energies with it
!
      wf%n_mo = wf%n_o + wf%n_v
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
      call wf%initialize_pivot_matrix_ao_overlap()
      call wf%initialize_cholesky_ao_overlap()
!
      call wf%read_orbital_coefficients()
      call wf%update_ao_density()
      call wf%read_orbital_energies()
!
!     Allocate active mo specific arrays
!     and read G(De) in the AO and MO basis
!
      call wf%initialize_W_mo_update()
      call wf%initialize_mo_fock()
      call wf%initialize_G_De_ao()
      call wf%initialize_G_De()
!
      call wf%G_De_ao_file%open_('read', 'rewind')
      call wf%G_De_ao_file%read_(wf%G_De_ao, wf%n_ao**2)
      call wf%G_De_ao_file%close_
!
      call wf%mlhf_inactive_fock_term_file%open_('read', 'rewind')
      call wf%mlhf_inactive_fock_term_file%read_(wf%G_De, wf%n_mo**2)
      call wf%mlhf_inactive_fock_term_file%close_
!
      call identity_array(wf%W_mo_update, wf%n_mo)
!
   end subroutine read_for_scf_restart_mlhf
!
!
   subroutine do_initial_full_space_optimization_mlhf(wf)
!!
!!    Do initial full space optimization
!!    Written by Anders Hutcheson and Linda Goletto, 2019
!!
!!    Performs an initial optimization of the full space wavefunction
!!    to a low threshold that can be read from input (otherwise set 
!!    to 1.0d-1), in order to start the multilevelcalculation with
!!    a better starting density.
!!
!
      use scf_diis_hf_class, only: scf_diis_hf
!
      implicit none
!
      class(mlhf) :: wf
!
      type(hf), allocatable            :: full_space_wf
      type(scf_diis_hf), allocatable   :: full_space_solver
!
      integer :: max_iterations, diis_dimension
      character(len=200) :: ao_density_guess, storage
!
!     Set the default settings
!
      max_iterations = 20
      diis_dimension = 8
      ao_density_guess = 'sad'
      storage = 'memory'
!
      call output%printf('m', '- Initial full hf optimization to a gradient threshold of (e9.2)', &
                  reals=[wf%full_space_hf_threshold], &
                  ffs='(/t3,a)', fs='(t6,a)')
!
!     Perform the HF calculation in the full space
!
      if(.not. wf%print_initial_hf) call output%mute()
!
      full_space_wf = hf(wf%system)
!
      full_space_solver = scf_diis_hf(wf=full_space_wf,                 &
                        restart=.false.,                                &
                        diis_dimension=diis_dimension,                  &
                        ao_density_guess=ao_density_guess,              &
                        energy_threshold=wf%full_space_hf_threshold,    &
                        max_iterations=max_iterations,                  &
                        gradient_threshold=wf%full_space_hf_threshold,  &
                        storage=storage,                                &
                        cumulative_threshold=1.0d-2,                    &
                        crop=.false.)
!
      call full_space_solver%run(full_space_wf)
!
!     Update the MLHF orbital coefficients, orbital energies and density
!
      wf%orbital_coefficients = full_space_wf%orbital_coefficients
      wf%orbital_energies = full_space_wf%orbital_energies
!
      call full_space_wf%cleanup()
!
      if(.not. wf%print_initial_hf) call output%unmute()
!
      call wf%update_ao_density()
!
   end subroutine do_initial_full_space_optimization_mlhf
!
!
   subroutine prepare_frozen_fock_terms_mlhf(wf)
!!
!!    Prepare frozen Fock contributions
!!    Written by Sarai D. Folkestad, Oct 2019
!!
!!    This routine prepares the frozen Fock contributions
!!    to coupled cluster. 
!!
!!    Always included
!!
!!       - G(De)
!!
!!    Possible contributions to frozen fock:
!!
!!       - Frozen core
!!
!!       - Frozen HF orbitals
!!
!!       - MM fock (QM/MM)
!!
!!       - PCM fock
!!
!!    Modified by Linda Goletto, Nov 2019
!!
!!    In case of a reduction, the MLHF inactive fock term 
!!    has to be updated to the new MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: mo_fc_fock
      real(dp), dimension(:,:), allocatable :: mo_frozen_hf_fock
      real(dp), dimension(:,:), allocatable :: mo_frozen_mlhf_fock
!
      wf%exists_frozen_fock_terms = .true.
!
      call wf%initialize_mo_fock_frozen()
      call zero_array(wf%mo_fock_frozen, wf%n_mo**2)
!
!     Contribution from frozen core orbitals
!
      if (wf%frozen_core) then
!
         call mem%alloc(mo_fc_fock, wf%n_mo, wf%n_mo)
!
         call wf%construct_mo_fock_fc_term(mo_fc_fock)
         call daxpy(wf%n_mo**2, one, mo_fc_fock, 1, wf%mo_fock_frozen, 1)
!
         call mem%dealloc(mo_fc_fock, wf%n_mo, wf%n_mo)
!
      endif
!
      if (wf%frozen_hf_mos) then
!
         call mem%alloc(mo_frozen_hf_fock, wf%n_mo, wf%n_mo)
!
         call wf%construct_mo_fock_frozen_hf_term(mo_frozen_hf_fock)
         call daxpy(wf%n_mo**2, one, mo_frozen_hf_fock, 1, wf%mo_fock_frozen, 1)
!
         call mem%dealloc(mo_frozen_hf_fock, wf%n_mo, wf%n_mo)
!
      endif
!
!     Add G(De)
!
      call mem%alloc(mo_frozen_mlhf_fock, wf%n_mo, wf%n_mo)
!
      call wf%mo_transform(wf%G_De_ao, mo_frozen_mlhf_fock)
!
      call daxpy(wf%n_mo**2, one, mo_frozen_mlhf_fock, 1, wf%mo_fock_frozen, 1)
!
      call mem%dealloc(mo_frozen_mlhf_fock, wf%n_mo, wf%n_mo)
!
   end subroutine prepare_frozen_fock_terms_mlhf
!
!
   subroutine diagonalize_fock_frozen_hf_orbitals_mlhf(wf)
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
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_effective
!
!     We do one Roothan-Hall step to get a diagonal Fock matrix 
!     (this should only entail occupied-occupied and virtual-virtual orbital mixing.)
!
      call mem%alloc(F_effective, wf%n_ao, wf%n_ao)
      call wf%initialize_mo_fock()
!
!     Add the MLHF inactive Fock term in the AO basis, G_De_ao, to the AO Fock
!
      call dcopy(wf%n_ao**2, wf%ao_fock, 1, F_effective, 1)
      call daxpy(wf%n_ao**2, one, wf%G_De_ao, 1, F_effective, 1)      
!
      call wf%mo_transform(F_effective, wf%mo_fock)
!
      call mem%dealloc(F_effective, wf%n_ao, wf%n_ao)
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
   end subroutine diagonalize_fock_frozen_hf_orbitals_mlhf
!
!
   function get_n_active_hf_atoms_mlhf(wf) result(n_active_hf_atoms)
!!
!!    Get number of active hf atoms
!!    Written by Sarai D. Folkestad and Linda Goletto, Dec 2019
!!
!!    Sets the number of active hf atoms in the system 
!!
      implicit none 
!
      class(mlhf), intent(in)  :: wf
!
      integer :: n_active_hf_atoms, n_active_spaces
!
      n_active_spaces = wf%system%n_active_atom_spaces
!
      n_active_hf_atoms = wf%system%active_atom_spaces(n_active_spaces)%last_atom
!
   end function get_n_active_hf_atoms_mlhf
!
!
   subroutine print_banner_mlhf(wf)
!!
!!    Print banner
!!    Sarai D. Folkestad, Dec 2019
!!
!
      use string_utilities, only : convert_to_uppercase
!
      implicit none
!
      class(mlhf) :: wf
!
      call wf%hf%print_banner()
!
      call output%printf('m', '- MLHF settings:', fs='(/t3,a)')
!
      call output%printf('m', 'Occupied orbitals:    Cholesky', fs='(/t6,a)')
!
      if (wf%cholesky_virtuals) then
!     
         call output%printf('m', 'Virtual orbitals:     Cholesky', fs='(t6,a)')
!
      else
!     
         call output%printf('m', 'Virtual orbitals:     PAOs', fs='(t6,a)')
!
      endif
!
      call output%printf('m', 'Cholesky decomposition threshold: (e9.2)', fs='(/t6,a)', &
            reals=[wf%cholesky_threshold])
!
      if (wf%full_space_optimization) then
!
         call output%printf('m', 'Initial optimization of full AO density enabled', fs='(/t6,a)')
!
      endif
!
   end subroutine print_banner_mlhf
!
!
   subroutine prepare_mos_mlhf(wf)
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
!!
      implicit none
!
      class(mlhf) :: wf
!
!     Destruct MO quantities in the old MO dimension, if they are allocated,
!     before n_mo changes
!
      call wf%destruct_mo_fock()
      call wf%destruct_W_mo_update()
      call wf%destruct_G_De()
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
   end subroutine prepare_mos_mlhf
!
!
end module mlhf_class
