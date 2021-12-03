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
   use stream_file_class, only : stream_file
   use string_utilities, only : convert_to_uppercase
!
!
   implicit none
!
!  Multilevel Hartree-Fock wavefunction
!
   type, extends(hf) :: mlhf
!
      type(stream_file) :: mlhf_file
!
      real(dp), dimension(:,:), allocatable :: G_De_imo
      real(dp), dimension(:,:), allocatable :: G_De_mo
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
      real(dp), dimension(:,:), allocatable :: imo_to_mo ! Unitary transformation between initial MO basis (IMO)
                                                           ! and current MO basis
!
      real(dp), dimension(:,:), allocatable :: imo_fock ! Fock matrix in initial mo basis
!
   contains
!
      procedure :: prepare                                  => prepare_mlhf
      procedure :: cleanup                                  => cleanup_mlhf
      procedure :: prepare_for_roothan_hall                 => prepare_for_roothan_hall_mlhf
      procedure :: print_energy                             => print_energy_mlhf
      procedure :: print_banner                             => print_banner_mlhf
!
      procedure :: read_for_scf_restart                     => read_for_scf_restart_mlhf
!
      procedure :: construct_active_mos                     => construct_active_mos_mlhf
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
      procedure :: initialize_G_De_mo                       => initialize_G_De_mo_mlhf
      procedure :: initialize_G_De_imo                      => initialize_G_De_imo_mlhf
      procedure :: destruct_G_De_mo                         => destruct_G_De_mo_mlhf
      procedure :: destruct_G_De_imo                        => destruct_G_De_imo_mlhf
      procedure :: initialize_G_De_ao                       => initialize_G_De_ao_mlhf
      procedure :: destruct_G_De_ao                         => destruct_G_De_ao_mlhf
!
      procedure :: get_active_energy_G_De_term              => get_active_energy_G_De_term_mlhf
!
      procedure :: construct_G_De                           => construct_G_De_mlhf
!
      procedure :: update_fock_and_energy                   => update_fock_and_energy_mlhf
!
      procedure :: prepare_for_cc                           => prepare_for_cc_mlhf
      procedure :: prepare_frozen_fock_terms                => prepare_frozen_fock_terms_mlhf
      procedure :: diagonalize_fock_frozen_hf_orbitals      => diagonalize_fock_frozen_hf_orbitals_mlhf
      procedure :: get_n_active_hf_atoms                    => get_n_active_hf_atoms_mlhf
      procedure :: prepare_mos                              => prepare_mos_mlhf
!
      procedure :: get_full_idempotent_density              => get_full_idempotent_density_mlhf
!
      procedure :: get_F &
                => get_F_mlhf
!
      procedure :: get_gradient &
                => get_gradient_mlhf
!
      procedure :: set_C_and_e &
                => set_C_and_e_mlhf
!
      procedure :: initialize_imo_to_mo                      => initialize_imo_to_mo_mlhf
      procedure :: destruct_imo_to_mo                        => destruct_imo_to_mo_mlhf
!
      procedure :: initialize_imo_fock                      => initialize_imo_fock_mlhf
      procedure :: destruct_imo_fock                        => destruct_imo_fock_mlhf
!
      procedure :: get_nuclear_dipole                       => get_nuclear_dipole_mlhf
      procedure :: get_nuclear_quadrupole                   => get_nuclear_quadrupole_mlhf
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
   function new_mlhf() result(wf)
!!
!!    New MLHF
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(mlhf) :: wf
!
      wf%name_ = 'mlhf'
!
      wf%cholesky_threshold            = 1.0d-2
      wf%full_space_hf_threshold       = 1.0d-1
      wf%full_space_optimization       = .false.
      wf%minimal_basis_diagonalization = .false.
      wf%cholesky_virtuals             = .false.
      wf%print_initial_hf              = .false.
      wf%cumulative_fock               = .false.
!
      wf%cumulative_fock_threshold = 1.0d0
!
      call wf%read_settings()
!
      call wf%print_banner()
!
   end function new_mlhf
!
!
   subroutine prepare_mlhf(wf, centers, embedding, charge)
!!
!!    Prepare
!!    Written by Linda Goletto, Sarai D. Folkestad
!!    and Eirik F. Kjønstad, 2019
!!
!!    Initializes files, writes the restart file used for consistency checks
!!    and constructs screening vectors
!!
      use atomic_center_class, only: atomic_center
!
      implicit none
!
      class(mlhf) :: wf
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers       
      integer, intent(in), optional :: charge 
!
      logical, intent(in), optional :: embedding 
!
      wf%orbital_file = stream_file('hf_orbitals')
!
      call wf%prepare_ao_tool(centers=centers, charge=charge)
      call wf%prepare_embedding(embedding)
      if (wf%embedded) call wf%embedding%print_description
!
      if (wf%ao%has_ghost_atoms()) &
         call output%warning_msg("Ghosts are experimental in multilevel.")
!
      wf%n_densities = 1
!
      call wf%set_n_mo()
!
!     Initialize other MLHF files
!
      wf%mlhf_file = stream_file('mlhf_restart_file')
!
!     Construct frozen CC^T
!
      call wf%initialize_frozen_CCT()
      call zero_array(wf%frozen_CCT, wf%ao%n**2)
!
      call wf%set_screening_and_precision_thresholds(wf%gradient_threshold)
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
      real(dp), dimension(:,:), allocatable :: h
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
      call dscal(wf%ao%n**2, half, wf%ao_fock, 1)
!
!     Add the one-electron part
! 
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
      call daxpy(wf%ao%n**2, one, h, 1, wf%ao_fock, 1)
!
      call timer%turn_off()
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, h)
!
      n_electrons = wf%get_n_electrons_in_density()
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
      wf%packed_gradient_dimension = (wf%n_mo**2)*wf%n_densities
!
      call wf%initialize_imo_to_mo()
      call identity_array(wf%imo_to_mo, wf%n_mo)
!
      call wf%initialize_G_De_ao()
      call wf%initialize_G_De_mo()
      call wf%initialize_G_De_imo()
      call wf%initialize_mo_fock()
      call wf%initialize_imo_fock()
!
      call wf%construct_G_De()
!
!     Inactive energy contribution
!
      wf%inactive_energy = wf%calculate_hf_energy_from_G(wf%G_De_ao, h)
!
!     Construct active ao density
!
      call wf%construct_ao_density()
!
      call wf%mlhf_file%open_('write', 'rewind')
      call wf%mlhf_file%write_(wf%n_o)
      call wf%mlhf_file%write_(wf%n_v)
      call wf%mlhf_file%write_(wf%inactive_energy)
      call wf%mlhf_file%write_(wf%G_De_ao, wf%ao%n**2)
      call wf%mlhf_file%close_
!
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!
   end subroutine prepare_for_roothan_hall_mlhf
!
!
   subroutine update_fock_and_energy_mlhf(wf, cumulative)
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
!
      use array_utilities, only: symmetric_sandwich_right_transposition, symmetric_sandwich
!
      implicit none
!
      class(mlhf), intent(inout) :: wf
      logical, intent(in) :: cumulative 
!
      real(dp), dimension(:,:), allocatable :: G
      real(dp), dimension(:,:), allocatable :: h
!
      type(timings) :: timer
!
      timer = timings('AO Fock construction', pl='normal')
      call timer%turn_on()
!
      call mem%alloc(h, wf%ao%n, wf%ao%n)
      call wf%get_ao_h(h)
!
      if (cumulative) then
!
         call output%printf('v', 'Fock matrix construction using density differences')
!
!        Construct the two electron part of the Fock matrix (G),
!        and add the contribution to the Fock matrix
!
         call daxpy(wf%ao%n**2, -one, wf%previous_ao_density, 1, wf%ao_density, 1)
!
         call mem%alloc(G, wf%ao%n, wf%ao%n)
!
         call wf%construct_ao_G(wf%ao_density,                    &
                                G,                                &
                                C_screening=.true. )
!
         call daxpy(wf%ao%n**2, half, G, 1, wf%ao_fock, 1)
!
         call mem%dealloc(G, wf%ao%n, wf%ao%n)
!
         call daxpy(wf%ao%n**2, one, wf%previous_ao_density, 1, wf%ao_density, 1)         
!
      else
!
!        AO fock construction and energy calculation
!
!        Construct the two electron part of the Fock matrix (G),
!        and add the contribution to the Fock matrix
!
         call wf%construct_ao_G(wf%ao_density,                    &
                                wf%ao_fock,                       &
                                C_screening=.true.)
         call dscal(wf%ao%n**2, half, wf%ao_fock, 1)
!
!        Add the one-electron part
!
         call daxpy(wf%ao%n**2, one, h, 1, wf%ao_fock, 1)
!
      endif
!
      call timer%turn_off()
!
      wf%energy = wf%calculate_hf_energy_from_fock(wf%ao_fock, h)
!
      call mem%dealloc(h, wf%ao%n, wf%ao%n)
!
!     Transformation of the AO fock in the MO basis
!
      call wf%mo_transform(wf%ao_fock, wf%mo_fock)
!
!     Transform to IMO basis
!
      call symmetric_sandwich_right_transposition(wf%imo_fock, wf%mo_fock, wf%imo_to_mo, wf%n_mo, wf%n_mo)
!
!     Get G_De_mo from G_De_imo
!
      call symmetric_sandwich(wf%G_De_mo, wf%G_De_imo, wf%imo_to_mo, wf%n_mo, wf%n_mo)
!
!     Add G_De to Fock 
!
      call daxpy(wf%n_mo**2, half, wf%G_De_mo, 1, wf%mo_fock, 1)
      call daxpy(wf%n_mo**2, half, wf%G_De_imo, 1, wf%imo_fock, 1)
!
!     Add the Tr[Da * G(De)] and inactive energy contributions to the energy
!
      wf%energy = wf%energy + wf%inactive_energy + wf%get_active_energy_G_De_term()
!
   end subroutine update_fock_and_energy_mlhf
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
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: D_v
!
      if (wf%minimal_basis_diagonalization) then
!
         call invert(D_v, wf%ao%s, wf%ao%n)
         call daxpy(wf%ao%n**2, -one, wf%ao_density, 1, D_v, 1)
!
      else
!
         call dgemm('N', 'T',                              &
                   wf%ao%n,                                &
                   wf%ao%n,                                &
                   wf%n_v,                                 &
                   one,                                    &
                   wf%orbital_coefficients(1, wf%n_o + 1), &
                   wf%ao%n,                                &
                   wf%orbital_coefficients(1, wf%n_o + 1), &
                   wf%ao%n,                                &
                   zero,                                   &
                   D_v,                                    &
                   wf%ao%n)
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
         call wf%diagonalize_fock()
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
      call mem%alloc(significant_AOs, wf%ao%n)
!
      call wf%determine_minimal_basis(significant_AOs)
!
      n_significant_AOs = count_n_true(wf%ao%n, significant_AOs)
!
      call mem%alloc(Q_tr, wf%ao%n, n_significant_AOs)
      call extract_columns_of_matrix(n_significant_AOs, wf%ao%n, wf%ao%n, &
                                     wf%ao%s, Q_tr, significant_AOs)
!
      call mem%alloc(Q, n_significant_AOs, wf%ao%n)
!
!$omp parallel do private(I, J)
      do I = 1, n_significant_AOs
         do J = 1, wf%ao%n
!
            Q (I,J) = Q_tr(J,I)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(Q_tr, wf%ao%n, n_significant_AOs)
!
      call mem%alloc(S_minimal,  n_significant_AOs,  n_significant_AOs)
      call extract_columns_of_matrix(n_significant_AOs, n_significant_AOs, wf%ao%n, &
                                       Q, S_minimal, significant_AOs)
!
      call mem%dealloc(significant_AOs, wf%ao%n)
!
      call mem%alloc(S_minimal_inv,  n_significant_AOs,  n_significant_AOs)
!
      call invert(S_minimal_inv, S_minimal, n_significant_AOs)
!
      call mem%dealloc(S_minimal,  n_significant_AOs,  n_significant_AOs)
!
      call mem%alloc(T, n_significant_AOs, wf%ao%n)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%ao%n,             &
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
      call mem%dealloc(Q, n_significant_AOs, wf%ao%n)
      call mem%dealloc(S_minimal_inv,  n_significant_AOs,  n_significant_AOs)
!
!     Project F down on the minimal basis
!
      call mem%alloc(X, n_significant_AOs, wf%ao%n)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%ao%n,             &
                  wf%ao%n,             &
                  one,                 &
                  T,                   &
                  n_significant_AOs,   &
                  wf%ao_fock,          &
                  wf%ao%n,             &
                  zero,                &
                  X,                   &
                  n_significant_AOs)
!
      call mem%alloc(F_minimal, n_significant_AOs, n_significant_AOs)
!
      call dgemm('N', 'T',             &
                  n_significant_AOs,   &
                  n_significant_AOs,   &
                  wf%ao%n,             &
                  one,                 &
                  X,                   &
                  n_significant_AOs,   &
                  T,                   &
                  n_significant_AOs,   &
                  zero,                &
                  F_minimal,           &
                  n_significant_AOs)
!
      call mem%dealloc(X, n_significant_AOs, wf%ao%n)
!
!     Project S on the minimal basis
!
      call mem%alloc(X, n_significant_AOs, wf%ao%n)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%ao%n,             &
                  wf%ao%n,             &
                  one,                 &
                  T,                   &
                  n_significant_AOs,   &
                  wf%ao%s,       &
                  wf%ao%n,             &
                  zero,                &
                  X,                   &
                  n_significant_AOs)
!
      call mem%alloc(S_minimal, n_significant_AOs, n_significant_AOs)
!
      call dgemm('N', 'T',             &
                  n_significant_AOs,   &
                  n_significant_AOs,   &
                  wf%ao%n,             &
                  one,                 &
                  X,                   &
                  n_significant_AOs,   &
                  T,                   &
                  n_significant_AOs,   &
                  zero,                &
                  S_minimal,           &
                  n_significant_AOs)
!
      call mem%dealloc(X, n_significant_AOs, wf%ao%n)
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
      call mem%alloc(X, n_significant_AOs, wf%ao%n)
!
      call dgemm('N', 'N',             &
                  n_significant_AOs,   &
                  wf%ao%n,             &
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
                  wf%ao%n,             &
                  wf%ao%n,             &
                  n_significant_AOs,   &
                  one,                 &
                  T,                   &
                  n_significant_AOs,   &
                  X,                   &
                  n_significant_AOs,   &
                  zero,                &
                  wf%ao_density,       &
                  wf%ao%n)
!
      call mem%dealloc(T, n_significant_AOs, wf%ao%n)
      call mem%dealloc(X, n_significant_AOs, wf%ao%n)
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
      logical, dimension(wf%ao%n), intent(out) :: significant_AOs

      real(dp), parameter :: occupation_treshold = 0.1d0
!
      integer :: ao
!
      significant_AOs = .false.
!
!$omp parallel do private(ao)
      do ao = 1, wf%ao%n
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
      call input%get_keyword('cholesky threshold', 'multilevel hf', wf%cholesky_threshold)
      call input%get_keyword('initial hf threshold', 'multilevel hf', wf%full_space_hf_threshold)
!
      if (input%is_keyword_present('project on minimal basis', 'multilevel hf')) &
            wf%minimal_basis_diagonalization = .true.
!
      if (input%is_keyword_present('cholesky virtuals', 'multilevel hf')) &
            wf%cholesky_virtuals = .true.
!
      if (input%is_keyword_present('initial hf optimization', 'multilevel hf')) &
            wf%full_space_optimization = .true.
!
      if (input%is_keyword_present('print initial hf', 'multilevel hf')) &
            wf%print_initial_hf = .true.
!
!     Sanity checks
!
      if (input%is_keyword_present('initial hf threshold', 'multilevel hf') &
            .and. .not. wf%full_space_optimization) &
            call output%error_msg('Initial hf threshold specified, without specifying initial hf optimization')
!
      if (wf%full_space_optimization .and. wf%minimal_basis_diagonalization) &
            call output%error_msg('Projection on minimal basis and initial hf optimization are not compatible')
!
   end subroutine read_mlhf_settings_mlhf
!
!
   subroutine construct_active_paos_mlhf(wf, C_pao, n_active_aos)
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

      integer, intent(in)   :: n_active_aos
      real(dp), dimension(wf%ao%n, n_active_aos), intent(out)  :: C_pao
!
      real(dp), dimension(:,:), allocatable :: C_pao_copy
      real(dp), dimension(:,:), allocatable :: S_pao
      real(dp), dimension(:,:), allocatable :: D
!
      integer   :: rank
!
      call dscal(wf%ao%n**2, half, wf%ao_density, 1)
!
      call wf%project_atomic_orbitals(wf%ao_density, C_pao, n_active_aos)
!
      call mem%alloc(S_pao, n_active_aos, n_active_aos)
!
      call wf%get_orbital_overlap(C_pao, n_active_aos, S_pao)
!
      call wf%lowdin_orthonormalization(C_pao, S_pao, n_active_aos, wf%n_v)
!
      call mem%dealloc(S_pao, n_active_aos, n_active_aos)
!
!     Construct inactive virtuals and add to frozen CC^T. Needed for later PAO construction
!     e.g. to perform CC in reduced space or MLCC with PAOs
!
      call mem%alloc(C_pao_copy, wf%ao%n, wf%ao%n)
      call mem%alloc(D, wf%ao%n, wf%ao%n)
!
!     Active virtual density
!
      call dgemm('N', 'T', &
                  wf%ao%n, &
                  wf%ao%n, &
                  wf%n_v,  &
                  one,     &
                  C_pao,   &
                  wf%ao%n, &
                  C_pao,   &
                  wf%ao%n, &
                  zero,    &
                  D,       &
                  wf%ao%n)
!
!     Add AO density (occupied orbitals)
!
      call daxpy(wf%ao%n**2, one, wf%ao_density, 1, D, 1)
!
!     Project occupied and active virtual out of AOs
!
      call wf%project_atomic_orbitals(D, C_pao_copy, wf%ao%n)
!
      call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      call mem%alloc(S_pao, wf%ao%n, wf%ao%n)
!
      call wf%get_orbital_overlap(C_pao_copy, wf%ao%n, S_pao)
!
!     Orthonormalize
!
      call wf%lowdin_orthonormalization(C_pao_copy, S_pao, wf%ao%n, rank)
!
      call mem%dealloc(S_pao, wf%ao%n, wf%ao%n)
!
!     Add to CC^T
!
      call dgemm('N', 'T',       &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  rank,          &
                  one,           &
                  C_pao_copy,    &
                  wf%ao%n,       &
                  C_pao_copy,    &
                  wf%ao%n,       &
                  one,           &
                  wf%frozen_CCT, &
                  wf%ao%n)
!
      call mem%dealloc(C_pao_copy, wf%ao%n, wf%ao%n)      
!
      call dscal(wf%ao%n**2, two, wf%ao_density, 1)
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
!
      use array_utilities, only: symmetric_sandwich
!
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
         E_G_De = E_G_De + wf%G_De_mo(i,i)
!
      end do
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
   subroutine initialize_G_De_mo_mlhf(wf)
!!
!!    Initialize G(De) MO
!!    Written by Sarai D. Folkestad
!!
!!    Initializes G_De in the MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (.not. allocated(wf%G_De_mo)) call mem%alloc(wf%G_De_mo, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_G_De_mo_mlhf
!
!
   subroutine initialize_G_De_imo_mlhf(wf)
!!
!!    Initialize G(De) IMO
!!    Written by Ida-Marie Hoyvik, Oct 2019
!!
!!    Initializes G_De in the IMO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (.not. allocated(wf%G_De_imo)) call mem%alloc(wf%G_De_imo, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_G_De_imo_mlhf
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
      call mem%alloc(wf%G_De_ao, wf%ao%n, wf%ao%n)
!
   end subroutine initialize_G_De_ao_mlhf
!
!
   subroutine destruct_G_De_mo_mlhf(wf)
!!
!!    Destruct G(De) MO
!!    Written by Ida-Marie Hoyvik, Oct 2019
!!
!!    Destructs G_De in the MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (allocated(wf%G_De_mo)) call mem%dealloc(wf%G_De_mo, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_G_De_mo_mlhf
!
!
   subroutine destruct_G_De_imo_mlhf(wf)
!!
!!    Destruct G(De) IMO
!!    Written by Sarai D. Folkestad
!!
!!    Destructs G_De in the IMO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (allocated(wf%G_De_imo)) call mem%dealloc(wf%G_De_imo, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_G_De_imo_mlhf
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
      call mem%dealloc(wf%G_De_ao, wf%ao%n, wf%ao%n)
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
      call wf%destruct_ao_fock()
      call wf%destruct_imo_fock()
      call wf%destruct_mo_fock()
      call wf%destruct_ao_density()
!
      call wf%destruct_imo_to_mo()
      call wf%destruct_G_De_mo()
      call wf%destruct_G_De_imo()
      call wf%destruct_G_De_ao()
!
      call wf%destruct_mo_fock_frozen()
      call wf%destruct_orbital_coefficients_fc()
      call wf%destruct_orbital_coefficients_frozen_hf()
!
      call wf%destruct_frozen_CCT()
!
      deallocate(wf%ao)
!
      if (wf%embedded) deallocate(wf%embedding)
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
      call wf%construct_ao_G(wf%ao_density, wf%G_De_ao, C_screening=.true.)
!
      call wf%mo_transform(wf%G_De_ao, wf%G_De_imo)
      call dcopy(wf%n_mo**2, wf%G_De_imo, 1, wf%G_De_mo, 1)
!
      call G_De_construction_timer%turn_off()
!
   end subroutine construct_G_De_mlhf
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
      nuclear_repulsion = wf%get_nuclear_repulsion()
!
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
      use array_utilities, only : copy_and_scale
      use cholesky_orbital_tool_class, only: cholesky_orbital_tool
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: C_v, C_o
!
      integer     :: n_vectors
!
      integer     :: n_active_aos
      integer     :: last_hf_ao
!
      type(timings) :: active_mos_construction_timer
!
      type(cholesky_orbital_tool), allocatable :: cd_tool_o, cd_tool_v
!
      active_mos_construction_timer = timings('Active MOs construction', pl='normal')
      call active_mos_construction_timer%turn_on()
!
      n_active_aos = 0
      n_vectors    = 0
!
      call wf%ao%get_aos_in_subset('hf', last=last_hf_ao)
!
      n_active_aos =  last_hf_ao ! Because if there are CC active atoms, these are also HF active
!
!     Partition virtual space
!
      call mem%alloc(C_v, wf%ao%n, n_active_aos)
!
      if (wf%cholesky_virtuals) then
!
!        Cholesky orbitals
!
         cd_tool_v = cholesky_orbital_tool(wf%ao%n, wf%cholesky_threshold)
         call cd_tool_v%initialize_density()
         call cd_tool_v%set_density_from_orbitals(wf%orbital_coefficients(:,wf%n_o + 1:wf%n_mo), wf%n_v)
         call cd_tool_v%restricted_decomposition(C_v, wf%n_v, n_active_aos, 1)
!
!        Add inactive virtual density to frozen_CCT to be used if new PAOs are constructed later
!
         call daxpy(wf%ao%n**2, one, cd_tool_v%D, 1, wf%frozen_CCT, 1)
         call cd_tool_v%cleanup()
!
      else
!
!        PAOs
!
         call wf%construct_active_paos(C_v, n_active_aos)
!
      endif
!
!     Partition occupied space (Cholesky orbitals)
!
      call mem%alloc(C_o, wf%ao%n, n_active_aos)
!
      cd_tool_o = cholesky_orbital_tool(wf%ao%n, wf%cholesky_threshold)
      call cd_tool_o%initialize_density()
      call cd_tool_o%set_density_from_density(wf%ao_density, wf%n_o, factor=half)
      call cd_tool_o%restricted_decomposition(C_o, wf%n_o, n_active_aos, 1)
!
!     Add inactive density to frozen_CCT to be used if new PAOs are constructed later
!
      call daxpy(wf%ao%n**2, one, cd_tool_o%D, 1, wf%frozen_CCT, 1)
!
      call copy_and_scale(two, cd_tool_o%D, wf%ao_density, wf%ao%n**2)
      call cd_tool_o%cleanup()
!
!     Set molecular orbitals coefficients
!
      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
!
!     new dimension of mo space
      wf%n_mo = wf%n_o + wf%n_v
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%set_orbital_coefficients(C_o(:,1:wf%n_o), wf%n_o, 1)
      call wf%set_orbital_coefficients(C_v(:,1:wf%n_v), wf%n_v, wf%n_o + 1)
!
      call mem%dealloc(C_v, wf%ao%n, n_active_aos)
      call mem%dealloc(C_o, wf%ao%n, n_active_aos)
!
      call active_mos_construction_timer%turn_off()
!
   end subroutine construct_active_mos_mlhf
!
!
   subroutine read_for_scf_restart_mlhf(wf)
!!
!!    Read for SCF restart
!!    Written by Sarai D. Folkestad, Anders Hutcheson
!!    and Linda Goletto, Oct 2019
!!
!
      use array_utilities, only : symmetric_sandwich
!
      implicit none
!
      class(mlhf) :: wf
      integer(i64) :: n
!
!     Destruct orbital coeffiecients and orbital energies allocated
!     with the old n_mo
!
      call wf%destruct_orbital_coefficients()
      call wf%destruct_orbital_energies()
!
      call wf%orbital_file%open_('read', 'rewind')
      call wf%orbital_file%read_(n, (i64) + 1)
      wf%n_mo = int(n)
      call wf%orbital_file%close_('keep')
!
!     Calculate the new n_mo and initialize orbital coefficients and
!     orbital energies with it
!
      wf%packed_gradient_dimension = (wf%n_mo**2)*wf%n_densities
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%read_orbitals()
!
!     Allocate active mo specific arrays
!     and read G(De) in the AO and MO basis
!
      call wf%initialize_imo_to_mo()
!
      call wf%initialize_mo_fock()
      call wf%initialize_imo_fock()
      call wf%initialize_G_De_ao()
      call wf%initialize_G_De_mo()
      call wf%initialize_G_De_imo()
!
      call wf%mlhf_file%open_('read', 'rewind')
      call wf%mlhf_file%read_(wf%n_o)
      call wf%mlhf_file%read_(wf%n_v)
      call wf%mlhf_file%read_(wf%inactive_energy)
      call wf%mlhf_file%read_(wf%G_De_ao, wf%ao%n**2)
      call wf%mlhf_file%close_
!
      call symmetric_sandwich(wf%G_De_imo, wf%G_De_ao, wf%orbital_coefficients, wf%ao%n, wf%n_mo)
      call dcopy(wf%n_mo**2, wf%G_De_imo, 1, wf%G_De_mo, 1)
!      
      call identity_array(wf%imo_to_mo, wf%n_mo)
      call wf%update_ao_density()
!
      call wf%print_orbital_space_info()
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
      use scf_solver_class, only: scf_solver
      use scf_solver_factory_class, only: scf_solver_factory
!
      implicit none
!
      class(mlhf) :: wf
!
      type(hf), allocatable             :: full_space_wf
      class(scf_solver), allocatable    :: full_space_solver
!
      type(scf_solver_factory) :: factory
!
      call output%printf('m', '- Initial full hf optimization to a gradient threshold of (e9.2)', &
                  reals=[wf%full_space_hf_threshold], &
                  ffs='(/t3,a)', fs='(t6,a)')
!
!     Perform the HF calculation in the full space
!
      if(.not. wf%print_initial_hf) call output%mute()
!
      full_space_wf = hf()
!
      call full_space_wf%prepare()
!
      call full_space_wf%set_gradient_threshold(wf%full_space_hf_threshold)
      call full_space_wf%prepare_for_scf(restart=.false., skip=.false.)
!
      factory = scf_solver_factory(acceleration_type = 'diis', max_iterations = 20, &
                                   energy_threshold = wf%full_space_hf_threshold)
!
      call factory%create(full_space_wf, full_space_solver, restart=.false., skip=.false.)
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
      call daxpy(wf%n_mo**2, half, mo_frozen_mlhf_fock, 1, wf%mo_fock_frozen, 1)
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
!
      use array_utilities, only: block_diagonalize_symmetric
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(:,:), allocatable :: F_effective
!
      real(dp), dimension(:,:), allocatable  :: C_copy
      integer, dimension(2)                  :: block_dim
      integer                                :: n_blocks
!
!     We do one Roothan-Hall step to get a diagonal Fock matrix 
!     (this should only entail occupied-occupied and virtual-virtual orbital mixing.)
!
      call mem%alloc(F_effective, wf%ao%n, wf%ao%n)
      call wf%initialize_mo_fock()
!
!     Add the MLHF inactive Fock term in the AO basis, G_De_ao, to the AO Fock
!
      call dcopy(wf%ao%n**2, wf%ao_fock, 1, F_effective, 1)
      call daxpy(wf%ao%n**2, half, wf%G_De_ao, 1, F_effective, 1)
!
      call wf%mo_transform(F_effective, wf%mo_fock)
!
      call mem%dealloc(F_effective, wf%ao%n, wf%ao%n)
!
      call wf%initialize_imo_to_mo()
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
      call wf%destruct_imo_to_mo()
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
      integer :: n_active_hf_atoms
!
      n_active_hf_atoms = wf%ao%get_n_centers_in_subset('hf')
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
!
      use visualization_class, only : visualization
!
      implicit none
!
      class(mlhf) :: wf
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
      call wf%destruct_imo_fock()
      call wf%destruct_imo_to_mo()
      call wf%destruct_G_De_mo()
      call wf%destruct_G_De_imo()
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
         plotter = visualization(wf%ao)
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
         label = 'MLHF_density_for_CC'
!
         call plotter%plot_density(wf%ao, D, label)
!
         call mem%dealloc(D, wf%ao%n, wf%ao%n)
!
      endif

   end subroutine prepare_mos_mlhf
!
!
   subroutine get_full_idempotent_density_mlhf(wf, D)
!!
!!    Get full idempotent density
!!    Written by Sarai D. Folkestad, Jan 2020
!!
!!    Constructs the full occupied density
!!    for determining the frozen HF virtuals
!!
      implicit none
!
      class(mlhf), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(out) :: D
!
      call wf%hf%get_full_idempotent_density(D)
!
      call daxpy(wf%ao%n**2, one, wf%frozen_CCT, 1, D, 1)
!
   end subroutine get_full_idempotent_density_mlhf
!
!
   subroutine prepare_for_cc_mlhf(wf)
!!
!!    Prepare for CC
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Prepares frozen fock terms, 
!!    and places energy in hf_energy
!!
      implicit none
!
      class(mlhf) :: wf
!
      wf%hf_energy = wf%energy
!
      wf%exists_frozen_fock_terms = .true. ! Always true for MLHF
!
!     Change the MOs if frozen core or frozen hf 
!     is requested
!
      call wf%prepare_mos()
!
!     Prepare frozen Fock terms from frozen core 
!     and frozen HF
!
      call wf%prepare_frozen_fock_terms()
!
   end subroutine prepare_for_cc_mlhf
!
!
   subroutine get_F_mlhf(wf, F_packed)
!!
!!    Get F
!!    Written by Sarai D. Folkestad
!!
!!    Constructs the Fock matrix in the initial MO basis
!!    and returns it packed
!!
      use reordering, only: packin
!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_mo*(wf%n_mo + 1)/2*wf%n_densities), intent(out) :: F_packed
!
      call wf%update_fock_and_energy(cumulative = wf%cumulative_fock)
!
      call packin(F_packed, wf%imo_fock, wf%n_mo)
!
   end subroutine get_F_mlhf
!
!
   subroutine set_C_and_e_mlhf(wf, C, e)
!!
!!    Set C and e
!!    Written by Sarai D. Folkestad, 2020
!! 
!!    Sets the orbital coefficients from the orthonormal MO
!!    update matrix C
!!
!!    Sets the orbital energies (e)
!!
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, wf%n_densities), intent(in)  :: C
      real(dp), dimension(wf%n_mo, wf%n_densities), intent(in)           :: e
!
      real(dp), dimension(:,:), allocatable :: C_old
!
      call mem%alloc(C_old, wf%ao%n, wf%n_mo)
!
!     Back to initial MO basis 
!
      call dgemm('N', 'T',                   &
                  wf%ao%n,                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%ao%n,                   &
                  wf%imo_to_mo,              &
                  wf%n_mo,                   &
                  zero,                      &
                  C_old,                     &
                  wf%ao%n)
!
!     To current MO basis
!
      call dgemm('N', 'N',                   &
                  wf%ao%n,                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  one,                       &
                  C_old,                     &
                  wf%ao%n,                   &
                  C,                         &
                  wf%n_mo,                   &
                  zero,                      &
                  wf%orbital_coefficients,   &
                  wf%ao%n)
!
      call dcopy(wf%n_mo**2, C, 1, wf%imo_to_mo, 1)
      call dcopy(wf%n_mo, e, 1, wf%orbital_energies, 1)
!
      call mem%dealloc(C_old, wf%ao%n, wf%n_mo)
!
      call wf%update_ao_density()
!
      call wf%save_orbitals()
   end subroutine set_C_and_e_mlhf
!
!
   subroutine get_gradient_mlhf(wf, G)
!!
!!    Get gradient
!!    Written by Sarai D. Folkestad
!!
!!    Returns the gradient F_ov in the initial MO basis
!!
!!    If the gradient norm is sufficiently small, 
!!    'cumulative_fock' is enabled
!!
      use reordering, only: packin 
      use array_utilities, only: symmetric_sandwich
      implicit none
!
      class(mlhf) :: wf
!
      real(dp), dimension(wf%packed_gradient_dimension), intent(out) :: G
      real(dp), dimension(:,:), allocatable :: F,  X
!
      integer :: a, i 
!
      call mem%alloc(F, wf%n_o, wf%n_v)
!
!$omp parallel do private(i, a)
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            F(i, a) = -four*wf%mo_fock(i, a + wf%n_o)
!
         enddo
      enddo
!$omp end parallel do 
!
      call mem%alloc(X, wf%n_o, wf%n_mo)
!
      call dgemm('N', 'T',                       &
                  wf%n_o,                        &
                  wf%n_mo,                       &
                  wf%n_v,                        &
                  one,                           &
                  F,                             &
                  wf%n_o,                        &
                  wf%imo_to_mo(1, wf%n_o + 1),   &
                  wf%n_mo,                       &
                  zero,                          &
                  X,                             &
                  wf%n_o)
!
      call dgemm('N', 'N',        &
                  wf%n_mo,        &
                  wf%n_mo,        &
                  wf%n_o,         &
                  one,            &
                  wf%imo_to_mo,   &
                  wf%n_mo,        &
                  X,              &
                  wf%n_o,         &
                  zero,           &
                  G,              &
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_o,wf%n_mo)
      call mem%dealloc(F, wf%n_o, wf%n_v)
!
      wf%cumulative_fock = wf%can_do_cumulative_fock(G)
!
   end subroutine get_gradient_mlhf
!
!
   subroutine initialize_imo_to_mo_mlhf(wf)
!!
!!    Initialize IMO to MO 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Jan 2019
!!
!!    Modified by Ida-Marie Hoyvik, Oct 2019
!!
!!    Initializes the transformation matrix which transforms between
!!    initial and current MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (.not. allocated(wf%imo_to_mo)) call mem%alloc(wf%imo_to_mo, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_imo_to_mo_mlhf
!
!
   subroutine destruct_imo_to_mo_mlhf(wf)
!!
!!    Destruct IMO to MO 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Jan 2019
!!
!!    Destructs the transformation matrix which transforms between
!!    initial and current MO basis
!!
!!    Modified by Ida-Marie Hoyvik, Oct 2019
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (allocated(wf%imo_to_mo)) call mem%dealloc(wf%imo_to_mo, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_imo_to_mo_mlhf
!
!
   subroutine initialize_imo_fock_mlhf(wf)
!!
!!    Initialize IMO to MO 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Jan 2019
!!
!!    Modified by Ida-Marie Hoyvik, Oct 2019
!!
!!    Initializes the transformation matrix which transforms between
!!    initial and current MO basis
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (.not. allocated(wf%imo_fock)) call mem%alloc(wf%imo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_imo_fock_mlhf
!
!
   subroutine destruct_imo_fock_mlhf(wf)
!!
!!    Destruct IMO fock
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad 
!!    and Linda Goletto, Jan 2019
!!
!!    Destructs the transformation matrix which transforms between
!!    initial and current MO basis
!!
!!    Modified by Ida-Marie Hoyvik, Oct 2019
!!
      implicit none
!
      class(mlhf) :: wf
!
      if (allocated(wf%imo_fock)) call mem%dealloc(wf%imo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_imo_fock_mlhf
!
!
   function get_nuclear_dipole_mlhf(wf) result(d)
!!
!!    Get nuclear dipole 
!!    Written by Sarai D. Folkestad, 2021
!!
      use point_charges_class, only: point_charges
!
      implicit none 
!
      class(mlhf), intent(in) :: wf 
!
      real(dp), dimension(3) :: d 
!
      type(point_charges) :: pc
!
      call wf%ao%get_subset_point_charges(pc, 'hf', include_higher_priority = .true.)
      d  = pc%get_dipole()
!
   end function get_nuclear_dipole_mlhf
!
!
   function get_nuclear_quadrupole_mlhf(wf) result(q)
!!
!!    Get nuclear quadrupole 
!!    Written by Sarai D. Folkestad, 2021
!!
      use point_charges_class, only: point_charges
!
      implicit none 
!
      class(mlhf), intent(in) :: wf 
!
      real(dp), dimension(6) :: q 
!
      type(point_charges) :: pc
!
      call wf%ao%get_subset_point_charges(pc, 'hf', include_higher_priority = .true.)
      q  = pc%get_quadrupole()
!
   end function get_nuclear_quadrupole_mlhf
!
!
end module mlhf_class
