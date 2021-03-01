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
module wavefunction_class
!
!!
!!    Wavefunction class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use parameters
   use global_out,            only: output
   use memory_manager_class,  only: mem
   use global_in,             only: input
   use ao_tool_class,         only: ao_tool
   use environment_class,     only: environment
!
   implicit none
!
   type, abstract :: wavefunction
!
      character(len=40) :: name_
!
      real(dp)    :: energy
      real(dp)    :: hf_energy
      complex(dp) :: energy_complex
!
      complex(dp), dimension(3) :: dipole_moment_complex
!
      integer :: n_atomic_centers 
!
      integer :: n_mo ! Number of molecular orbitals 
      integer :: n_o  ! Number of occupied orbitals 
      integer :: n_v  ! Number of virtual orbbitals 
!
      logical :: embedded
!
      class(environment), allocatable :: embedding

      class(ao_tool), allocatable :: ao
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients
      real(dp), dimension(:), allocatable   :: orbital_energies
!
!     Frozen orbital variables. Frozen orbitals are typically frozen core or frozen HF orbitals.
!
      real(dp), dimension(:,:), allocatable :: mo_fock_frozen
!
      real(dp) :: cholesky_orbital_threshold = 1.0D-2
!
      logical :: exists_frozen_fock_terms ! Are there frozen Fock terms?
!
      real(dp), dimension(:,:), allocatable :: frozen_CCT   ! Matrix CC^T where the contraction is
                                                            ! over frozen orbitals. Needed if
                                                            ! PAO construction follows (e.g in mlcc)
                                                            ! other reducrtion of orbitals
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: ao_G
      real(dp), dimension(:,:), allocatable :: mo_fock
!
      real(dp), dimension(3) :: frozen_dipole
      real(dp), dimension(6) :: frozen_quadrupole
!
   contains
!
      procedure :: initialize_orbital_coefficients &
                => initialize_orbital_coefficients_wavefunction
!
      procedure :: initialize_orbital_energies &    
                => initialize_orbital_energies_wavefunction
!
      procedure :: destruct_orbital_coefficients &
                => destruct_orbital_coefficients_wavefunction
!
      procedure :: destruct_orbital_energies &
                => destruct_orbital_energies_wavefunction
!
      procedure :: initialize_frozen_CCT &
                => initialize_frozen_CCT_wavefunction
!
      procedure :: destruct_frozen_CCT &
                => destruct_frozen_CCT_wavefunction
!
      procedure :: initialize_ao_fock &
                => initialize_ao_fock_wavefunction
!
      procedure :: destruct_ao_fock &
                => destruct_ao_fock_wavefunction
!
      procedure :: initialize_mo_fock &
                => initialize_mo_fock_wavefunction
!
      procedure :: destruct_mo_fock &
                => destruct_mo_fock_wavefunction
!
      procedure :: mo_transform &
                => mo_transform_wavefunction
!
      procedure :: project_atomic_orbitals  &
                => project_atomic_orbitals_wavefunction
!
      procedure :: get_orbital_overlap &
                => get_orbital_overlap_wavefunction
!
      procedure :: lowdin_orthonormalization &
                => lowdin_orthonormalization_wavefunction
!
      procedure :: construct_orbital_block_by_density_cd &
                => construct_orbital_block_by_density_cd_wavefunction
!
      procedure :: initialize_mo_fock_frozen &
                => initialize_mo_fock_frozen_wavefunction
!
      procedure :: destruct_mo_fock_frozen &
                => destruct_mo_fock_frozen_wavefunction
!
      procedure :: get_mo_fock &
                => get_mo_fock_wavefunction
!
      procedure :: set_mo_fock &
                => set_mo_fock_wavefunction
!
      procedure :: prepare_ao_tool &
                => prepare_ao_tool_wavefunction
!
      procedure :: set_geometry &
                => set_geometry_wavefunction
!
      procedure :: prepare_embedding &
                => prepare_embedding_wavefunction 
!
      procedure :: contruct_mo_basis_transformation &
                => contruct_mo_basis_transformation_wavefunction
!
      procedure :: get_nuclear_dipole &
                => get_nuclear_dipole_wavefunction
!  
      procedure :: get_nuclear_quadrupole &
                => get_nuclear_quadrupole_wavefunction
!  
      procedure :: get_nuclear_repulsion &
                => get_nuclear_repulsion_wavefunction
!  
      procedure :: get_nuclear_repulsion_1der &
                => get_nuclear_repulsion_1der_wavefunction
!
      procedure :: get_molecular_geometry &
                => get_molecular_geometry_wavefunction
!
   end type wavefunction 
!
!
contains
!
!
   subroutine prepare_ao_tool_wavefunction(wf, centers, template)
!!
!!    Prepare AO tool
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Initializes the AO tool, which handles the atomic orbital integrals in eT.
!!
!!    centers:   Optional set of 'atomic_center' objects that represents the atomic orbital 
!!               basis used by the AO tool and the integral program (Libint). The default 
!!               centers are those specified by the QM geometry in the eT input file.
!!
!!    template: Optional template_ao_tool from which the ao_tool is initialized.
!!               
      use atomic_center_class, only: atomic_center
!
      implicit none 
!
      class(wavefunction) :: wf 
!
      class(atomic_center), dimension(:), optional, intent(in) :: centers
!
      type(ao_tool), optional, intent(in) :: template
!
      wf%ao = ao_tool()
!
      if (present(template)) then
!
         call wf%ao%initialize_ao_tool_from_template(template)
!
      else 
!
         call wf%ao%initialize(centers)
!
      end if
!
      wf%n_atomic_centers = wf%ao%get_n_centers()
!
   end subroutine prepare_ao_tool_wavefunction
!
!
   subroutine set_geometry_wavefunction(wf, R, units)
!!
!!    Set geometry 
!!    Written by Eirik F. Kjønstad, 2020 
!!
!!    Sets the position of the atomic nuclei. 
!!
!!       R:     (3 x n_atoms) array containing the Cartesian coordinates of the atoms 
!!       units: specifies units of R ('angstrom' or 'bohr')
!!
!!    This includes updating the atomic centers of the AO tool 
!!    and the atomic positions of the molecule.
!!
      implicit none
!
      class(wavefunction), intent(inout) :: wf 
!
      character(len=*), intent(in) :: units
!
      real(dp), dimension(3, wf%n_atomic_centers) :: R 
!
      call wf%ao%set_atomic_centers(R, units)
!
   end subroutine set_geometry_wavefunction
!
!
   subroutine initialize_orbital_coefficients_wavefunction(wf)
!!
!!    Initialize MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_coefficients)) call mem%alloc(wf%orbital_coefficients, wf%ao%n, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_wavefunction
!
!
   subroutine destruct_orbital_coefficients_wavefunction(wf)
!!
!!    Destruct MO coefficients
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_coefficients)) call mem%dealloc(wf%orbital_coefficients, wf%ao%n, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_wavefunction
!
!
   subroutine initialize_orbital_energies_wavefunction(wf)
!!
!!    Initialize orbital energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%orbital_energies)) call mem%alloc(wf%orbital_energies, wf%n_mo)
!
   end subroutine initialize_orbital_energies_wavefunction
!
!
   subroutine destruct_orbital_energies_wavefunction(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%orbital_energies)) call mem%dealloc(wf%orbital_energies, wf%n_mo)
!
   end subroutine destruct_orbital_energies_wavefunction
!
!
   subroutine mo_transform_wavefunction(wf, X_wx, Y_pq)
!!
!!    MO transform
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Performs MO transformation of X and saves the result in Y:
!!
!!       Y_pq = sum_wx C_wp X_wx C_xq
!!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in)    :: X_wx 
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Y_pq  
!
      real(dp), dimension(:,:), allocatable :: Z_wq ! = sum_x X_wx C_xq
!
      call mem%alloc(Z_wq, wf%ao%n, wf%n_mo)
!
      call dgemm('N', 'N',                 &
                  wf%ao%n,                 &
                  wf%n_mo,                 &
                  wf%ao%n,                 &
                  one,                     &
                  X_wx,                    &
                  wf%ao%n,                 &
                  wf%orbital_coefficients, & ! C_xq
                  wf%ao%n,                 &
                  zero,                    &
                  Z_wq,                    &
                  wf%ao%n)
!
      call dgemm('T', 'N',                 &
                  wf%n_mo,                 &
                  wf%n_mo,                 &
                  wf%ao%n,                 &
                  one,                     &
                  wf%orbital_coefficients, & ! C_wp 
                  wf%ao%n,                 &
                  Z_wq,                    &
                  wf%ao%n,                 &
                  zero,                    &
                  Y_pq,                    &
                  wf%n_mo)
!
      call mem%dealloc(Z_wq, wf%ao%n, wf%n_mo)
!
   end subroutine mo_transform_wavefunction
!
!
   subroutine project_atomic_orbitals_wavefunction(wf, D, PAO_coeff, n_orbitals, first_ao)
!!
!!    Project atomic orbitals
!!    Written by Linda Goletto and Sarai D. Folkestad, Jun 2019
!!
!!    Projects the orbitals given by D out of the atomic orbitals (PAOs)
!!
!!       C^PAO = I - DS
!!
!!    NOTE: n_orbitals must be equal to n_ao unless
!!    a restricted orbital space is requested.
!!
!!    For restricted space optional argument first_ao
!!    may be used.
!!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      integer, intent(in) :: n_orbitals
!  
      real(dp), dimension(wf%ao%n, wf%ao%n), intent(in)      :: D
      real(dp), dimension(wf%ao%n, n_orbitals), intent(out)  :: PAO_coeff
!
      integer, intent(in), optional :: first_ao
!
      integer :: local_first, local_last
!
      integer :: i
!
      real(dp), dimension(:,:), allocatable :: S
!
!     Set offsets
!
      local_first = 1
      if (present(first_ao)) local_first = first_ao
!
      local_last = local_first + n_orbitals - 1
!
!     Sanity checks and safety guards
!
      if (n_orbitals .gt. wf%ao%n) call output%error_msg('number of orbitals for PAOs exceeds number of AOs')
      if ((local_first .lt. 1) .or. (local_first .gt. wf%ao%n)) call output%error_msg('First PAO exceeds number of AOs')
      if ((local_last .lt. 1) .or. (local_last .gt. wf%ao%n)) call output%error_msg('Last PAO exceeds number of AOs')
!
!     Construct PAO coefficients
!
      call zero_array(PAO_coeff, (wf%ao%n)*(n_orbitals))
!
!$omp parallel do private(i)
      do i = 1, n_orbitals
!
         PAO_coeff(i + local_first - 1, i) = one
!
      enddo
!$omp end parallel do
!  
      call mem%alloc(S, wf%ao%n, wf%ao%n)
!
      call wf%ao%get_oei('overlap', S)
!
      call dgemm('N', 'N',             &
                  wf%ao%n,             &
                  n_orbitals,          &
                  wf%ao%n,             &
                  -one,                &
                  D,                   &
                  wf%ao%n,             &
                  S(1,local_first),    &
                  wf%ao%n,             &
                  one,                 &
                  PAO_coeff,           &
                  wf%ao%n)
!
      call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
   end subroutine project_atomic_orbitals_wavefunction
!
!
   subroutine get_orbital_overlap_wavefunction(wf, orbital_coeff, n_orbitals, S)
!!
!!    Get orbital overlap
!!    Written by Sarai D. Folkestad and Linda Goletto, Jun 2019
!!
!!    Construct the orbital overlap
!!
!!       S = C^T S_AO C
!!
!!    for some set of orbital coefficients.
!!
      implicit none
!
      class(wavefunction) :: wf
!
      integer, intent(in) :: n_orbitals
!
      real(dp), dimension(wf%ao%n, n_orbitals), intent(in) :: orbital_coeff
!
      real(dp), dimension(n_orbitals, n_orbitals), intent(out) :: S
!
      real(dp), dimension(:,:), allocatable :: S_ao, X
!
      call mem%alloc(S_ao, wf%ao%n, wf%ao%n)
      call wf%ao%get_oei('overlap', S_ao)
!
      call mem%alloc(X, n_orbitals, wf%ao%n)
!
      call dgemm('T', 'N',       &
                  n_orbitals,    &
                  wf%ao%n,       &
                  wf%ao%n,       &
                  one,           &
                  orbital_coeff, &
                  wf%ao%n,       &
                  S_ao,          &
                  wf%ao%n,       &
                  zero,          &
                  X,             &
                  n_orbitals)
!
      call mem%dealloc(S_ao, wf%ao%n, wf%ao%n)
!
      call dgemm('N', 'N',       &
                  n_orbitals,    &
                  n_orbitals,    &
                  wf%ao%n,       &
                  one,           &
                  X,             &
                  n_orbitals,    &
                  orbital_coeff, &
                  wf%ao%n,       &
                  zero,          &
                  S,             &
                  n_orbitals)
!
      call mem%dealloc(X, n_orbitals, wf%ao%n)
!
   end subroutine get_orbital_overlap_wavefunction
!
!
   subroutine lowdin_orthonormalization_wavefunction(wf, orbital_coeff, S, n_orbitals, rank)
!!
!!    Lövdin orthonormalization
!!    Written by Linda Goletto and Sarai D. Folkestad, Jun 2019
!!
!!    Orthonormalizes the orbital_coeff using Lövdin
!!    orthonormalization
!!
!!    The orbital overlap
!!
!!       S = C^T S_AO C
!!
!!    is diagonalized
!!
!!       S = U λ U^T
!!
!!    and the orbitals are updated according to
!!
!!       C = C U λ^(-1/2).
!!
!!    Linear dependence is removed by screening on the eigenvalues
!!    with a threshold of 1.0 * 10^-6
!!
!!    Modified by SDF May 2020,
!!
!!    Added flip of orbitals.
!!    Diagonalization using wrapper.
!!
      use array_utilities, only: copy_and_scale, diagonalize_symmetric
!
      implicit none
!
      class(wavefunction) :: wf
!
      integer, intent(in)  :: n_orbitals
      integer, intent(out) :: rank
!
      real(dp), dimension(wf%ao%n, n_orbitals), intent(in)     :: orbital_coeff
      real(dp), dimension(n_orbitals, n_orbitals), intent(inout)  :: S
!
      real(dp), dimension(:), allocatable :: eigenvalues, inv_sqrt_eig
      real(dp), dimension(:,:), allocatable :: orbital_coeff_copy
!
      integer :: i
!
      real(dp), parameter :: threshold = 1.0d-6 ! Threshold for linear dependency.
!
!     Diagonalize S
!
      call mem%alloc(eigenvalues, n_orbitals)
!
      call dscal(n_orbitals**2, -one, S, 1)
      call diagonalize_symmetric(S, n_orbitals, eigenvalues)
!
!     Find the rank of S
!
      rank = 0
!
      do i = 1, n_orbitals
!
         if (abs(eigenvalues(i)) .lt. threshold) exit
!
         rank = rank + 1
!
      enddo
!
!     Invert and square-root the eigenvalues (screening done here)
!
      call mem%alloc(inv_sqrt_eig, rank)
!
      do i = 1, rank
!
         inv_sqrt_eig(i) = 1/(sqrt(abs(eigenvalues(i))))
!
      enddo
!
      call mem%dealloc(eigenvalues, n_orbitals)
!
!     S = U λ^(-1/2)
!
      do i = 1, rank
!
         S(:, i) = S(:, i)*inv_sqrt_eig(i)
!
      enddo
!
      call mem%dealloc(inv_sqrt_eig, rank)
!
      do i = rank + 1, n_orbitals
!
         S(:,i) = zero
!
      enddo
!
!     Transform orbital coefficients
!
      call mem%alloc(orbital_coeff_copy, wf%ao%n, n_orbitals)
!
      call copy_and_scale(one, orbital_coeff, orbital_coeff_copy, (n_orbitals)*(wf%ao%n))
!
      call dgemm('N', 'N',                   &
                  wf%ao%n,                   &
                  rank,                      &
                  n_orbitals,                &
                  one,                       &
                  orbital_coeff_copy,        &
                  wf%ao%n,                   &
                  S,                         &
                  n_orbitals,                &
                  zero,                      &
                  orbital_coeff,             &
                  wf%ao%n)
!
      call mem%dealloc(orbital_coeff_copy, wf%ao%n, n_orbitals)
!
   end subroutine lowdin_orthonormalization_wavefunction
!
!
   subroutine construct_orbital_block_by_density_cd_wavefunction(wf, D, n_vectors, threshold, mo_offset, active_aos)
!!
!!    Construct orbital block by Cholesky decomposition for density
!!    Written by Sarai D. Folkestad, Feb 2019
!!
!!    Cholesky decomposition of density D plus
!!    update of corresponding wavefunction MOs
!!
!!    See A. M. J. Sánchez de Merás, H. Koch,
!!    I. G. Cuesta, and L. Boman (J. Chem. Phys. 132, 204105 (2010))
!!    for more information on active space generation
!!    using Cholesky decomposition
!!
!!    'D' : Density to be decomposed
!!          Given in the AO basis
!!
!!    'n_vectors' : The number of Cholesky
!!                  vectors (orbitals) constructed
!!
!!    'threshold' : The Cholesky decomposition threshold
!!
!!    'mo_offset' : Offset used to set the new MOs
!!
!!    'active_aos' : List of AOs on active atoms (optional)
!!                   This is used in case we want to partially
!!                   decompose density to get active orbitals
!!
!!
!
      use array_utilities, only: cholesky_decomposition_limited_diagonal, full_cholesky_decomposition_effective
!
      implicit none
!
      class(wavefunction), intent(inout) :: wf
!
      real(dp), dimension(wf%ao%n,wf%ao%n), intent(inout) :: D
!
      integer, intent(out) :: n_vectors
!
      real(dp), intent(in) :: threshold
!
      integer, intent(in) :: mo_offset
!
      integer, dimension(:), optional :: active_aos
!
      integer :: n_active_aos, ao, mo
!
      real(dp), dimension(:,:), allocatable :: cholesky_vec
      integer, dimension(:), allocatable  :: keep_vectors
!
      call mem%alloc(cholesky_vec, wf%ao%n, wf%ao%n)
!
      if (present(active_aos)) then
!
!        Active space generation by CD choosing pivots
!        only on active atoms
!
         n_active_aos = size(active_aos)
!
         if (n_active_aos .gt. wf%ao%n) call output%error_msg('More active AOs than total AOs')
!      
         call cholesky_decomposition_limited_diagonal(D, cholesky_vec, wf%ao%n, &
                                                      n_vectors, threshold, &
                                                      n_active_aos, active_aos)
!
!        Set the  MOs to be the ones to be frozen for CC
!
!$omp parallel do private(ao, mo)
         do mo = 1, n_vectors
            do ao = 1, wf%ao%n
!
               wf%orbital_coefficients(ao, mo_offset + mo) = cholesky_vec(ao, mo)
!
            enddo
         enddo
!$omp end parallel do
!
      else
!
         call mem%alloc(keep_vectors, wf%ao%n)
!
!        Full CD of density (used for inactive densities)
!
         call full_cholesky_decomposition_effective(D, cholesky_vec, &
                                             wf%ao%n, n_vectors, &
                                             threshold, keep_vectors)
!
!
         call mem%dealloc(keep_vectors, wf%ao%n)
!
!
!$omp parallel do private(ao, mo)
         do mo = 1, n_vectors
            do ao = 1, wf%ao%n
!
               wf%orbital_coefficients(ao, mo_offset + mo) = cholesky_vec(ao, mo)
!
            enddo
         enddo
!$omp end parallel do
!
      endif
!
      call mem%dealloc(cholesky_vec, wf%ao%n, wf%ao%n)
!
   end subroutine construct_orbital_block_by_density_cd_wavefunction
!
!
   subroutine initialize_mo_fock_frozen_wavefunction(wf)
!!
!!    Initialize MO Fock frozen
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(wavefunction) :: wf
!
      call mem%alloc(wf%mo_fock_frozen, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_mo_fock_frozen_wavefunction
!
!
   subroutine destruct_mo_fock_frozen_wavefunction(wf)
!!
!!    Destruct MO Fock frozen
!!    Written by Sarai D. Folkestad, Sep. 2019
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%mo_fock_frozen)) call mem%dealloc(wf%mo_fock_frozen, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_mo_fock_frozen_wavefunction
!
!
   subroutine initialize_frozen_CCT_wavefunction(wf)
!!
!!    Initialize frozen CC^T
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(wavefunction) :: wf
!
      call mem%alloc(wf%frozen_CCT, wf%ao%n, wf%ao%n)
!
   end subroutine initialize_frozen_CCT_wavefunction
!
!
   subroutine destruct_frozen_CCT_wavefunction(wf)
!!
!!    Destruct frozen CC^T
!!    Written by Sarai D. Folkestad, Jan 2020
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%frozen_CCT)) call mem%dealloc(wf%frozen_CCT, wf%ao%n, wf%ao%n)
!
   end subroutine destruct_frozen_CCT_wavefunction
!
!
   subroutine initialize_ao_fock_wavefunction(wf)
!!
!!    Initialize AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%ao_fock)) call mem%alloc(wf%ao_fock, wf%ao%n, wf%ao%n)
      if (.not. allocated(wf%ao_G)) call mem%alloc(wf%ao_G, wf%ao%n, wf%ao%n)
!
   end subroutine initialize_ao_fock_wavefunction
!
!
   subroutine destruct_ao_fock_wavefunction(wf)
!!
!!    Destruct AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%ao_fock)) call mem%dealloc(wf%ao_fock, wf%ao%n, wf%ao%n)
      if (allocated(wf%ao_G)) call mem%dealloc(wf%ao_G, wf%ao%n, wf%ao%n)
!
   end subroutine destruct_ao_fock_wavefunction
!
!
   subroutine initialize_mo_fock_wavefunction(wf)
!!
!!    Initialize MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (.not. allocated(wf%mo_fock)) call mem%alloc(wf%mo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_mo_fock_wavefunction
!
!
   subroutine destruct_mo_fock_wavefunction(wf)
!!
!!    Destruct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(wavefunction) :: wf
!
      if (allocated(wf%mo_fock)) call mem%dealloc(wf%mo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_mo_fock_wavefunction
!
!
   subroutine get_mo_fock_wavefunction(wf, F)
!!
!!    Get MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Gets the MO Fock
!!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(:,:), intent(out) :: F
!
      call dcopy(wf%n_mo**2, wf%mo_fock, 1, F, 1)
!
   end subroutine get_mo_fock_wavefunction
!
!
   subroutine set_mo_fock_wavefunction(wf, F)
!!
!!    Set MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the MO Fock from input
!!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(:, :), intent(in) :: F
!
      call dcopy(wf%n_mo**2, F, 1, wf%mo_fock, 1)
!
   end subroutine set_mo_fock_wavefunction
!
!
   subroutine prepare_embedding_wavefunction(wf, embedding)
!!
!!    Prepare embedding
!!    Written by Sarai D. Folkestad
!!
!!    Initializes the embedding object
!!    Embedding is either non-polarizable
!!    or polarizable (QM/FQ or PCM)
!!
      use environment_factory_class, only: environment_factory
!      
      implicit none
!
      class(wavefunction), intent(inout)        :: wf
      type(environment_factory), allocatable    :: factory
!
      logical, intent(in), optional :: embedding
!
      if (present(embedding)) then
!
         wf%embedded = embedding
!
      else
!
         wf%embedded = input%is_embedding_on()
!
      endif
!
      if (wf%embedded) then
!
         factory = environment_factory()
!
         wf%embedding = factory%create()
         call wf%embedding%initialize(wf%ao)
!
      end if
!
   end subroutine prepare_embedding_wavefunction
!
!
   subroutine contruct_mo_basis_transformation_wavefunction(wf, C1, C2, T)
!!
!!    Construct MO basis transformation
!!    Written by Sarai D. Folekstad, Nov 2019
!!
!!    Constructs a transformation matrix 'T' which
!!    takes a matrix from one molecular orbital basis
!!    to another.
!!
!!    'C1' : coefficients of the MO basis we end up in
!!
!!    'C2' : coefficients of the MO basis we start out with
!!
!!    The transformation matrix is defined as
!!
!!       T = C1^T S C2
!!
!!    where S is the AO overlap matrix.
!!
      implicit none
!
      class(wavefunction), intent(in) :: wf
!
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(in) :: C1, C2
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(out) :: T
!
      real(dp), dimension(:,:), allocatable :: X, S
!
      call mem%alloc(S, wf%ao%n, wf%ao%n)
      call wf%ao%get_oei('overlap', S)
!
!     X = C1^T S
!
      call mem%alloc(X, wf%n_mo, wf%ao%n)
!
      call dgemm('T', 'N', &
                  wf%n_mo, &
                  wf%ao%n, &
                  wf%ao%n, &
                  one,     &
                  C1,      &
                  wf%ao%n, &
                  S,       &
                  wf%ao%n, &
                  zero,    &
                  X,       &
                  wf%n_mo)
!
      call mem%dealloc(S, wf%ao%n, wf%ao%n)
!
!     T = X C2
!
      call dgemm('N', 'N', &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%ao%n, &
                  one,     &
                  X,       &
                  wf%n_mo, &
                  C2,      &
                  wf%ao%n, &
                  zero,    &
                  T,       &
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_mo, wf%ao%n)
!
   end subroutine contruct_mo_basis_transformation_wavefunction
!
!
   function get_nuclear_dipole_wavefunction(wf) result(d)
!!
!!    Get nuclear dipole 
!!    Written by Eirik F. Kjønstad, Oct 2020
!!
      use point_charges_class,         only: point_charges
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(3) :: d 
!
      type(point_charges) :: pc
!
      call wf%ao%get_point_charges(pc)
      d  = pc%get_dipole()
!
   end function get_nuclear_dipole_wavefunction
!
!
   function get_nuclear_quadrupole_wavefunction(wf) result(q)
!!
!!    Get nuclear quadrupole 
!!    Written by Eirik F. Kjønstad, Oct 2020
!!
      use point_charges_class,         only: point_charges
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(6) :: q 
!
      type(point_charges) :: pc
!
      call wf%ao%get_point_charges(pc)
      q  = pc%get_quadrupole()
!
   end function get_nuclear_quadrupole_wavefunction
!
!
   function get_nuclear_repulsion_wavefunction(wf) result(E)
!!
!!    Get nuclear repulsion 
!!    Written by Eirik F. Kjønstad, Oct 2020
!!
      use point_charges_class, only: point_charges
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp) :: E 
!
      type(point_charges) :: pc
!
      call wf%ao%get_point_charges(pc)
      E = pc%get_coulomb_interaction()
!
   end function get_nuclear_repulsion_wavefunction
!
!
   function get_nuclear_repulsion_1der_wavefunction(wf) result(E)
!!
!!    Get nuclear repulsion 1der
!!    Written by Eirik F. Kjønstad, Oct 2020
!!
!!    Gets first derivative of nuclear repulsion
!!
      use point_charges_class,         only: point_charges
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(3, wf%ao%get_n_centers()) :: E 
!
      type(point_charges) :: pc
!
      call wf%ao%get_point_charges(pc)
      E  = pc%get_coulomb_interaction_1der()
!
   end function get_nuclear_repulsion_1der_wavefunction
!
!
   function get_molecular_geometry_wavefunction(wf) result(R)
!!
!!    Get molecular geometry 
!!    Written by Eirik F. Kjønstad, June 2019 
!!
      use parameters, only: angstrom_to_bohr
!
      implicit none 
!
      class(wavefunction), intent(in) :: wf 
!
      real(dp), dimension(3, wf%n_atomic_centers) :: R
!
      R = wf%ao%get_center_coordinates() 
      R = R * angstrom_to_bohr     
!
   end function get_molecular_geometry_wavefunction
!
!
end module wavefunction_class
