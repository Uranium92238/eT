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
module hf_class
!
!!
!!    Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
!
   use reordering
   use timings_class
   use array_utilities
   use array_analysis
   use interval_class
   use libint_initialization
!
   use omp_lib
!
   implicit none
!
!  Hartree-Fock wavefunction
!
   type, extends(wavefunction) :: hf
!
      real(dp), dimension(:,:), allocatable :: ao_density
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: mo_fock
!
      real(dp), dimension(:,:), allocatable :: ao_overlap
      real(dp), dimension(:,:), allocatable :: cholesky_ao_overlap
      real(dp), dimension(:,:), allocatable :: pivot_matrix_ao_overlap
!
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz
      integer, dimension(:,:), allocatable  :: sp_eri_schwarz_list
!
      real(dp) :: linear_dep_threshold = 1.0D-6
!
      real(dp) :: coulomb_threshold    = 1.0D-12   ! screening threshold
      real(dp) :: exchange_threshold   = 1.0D-10   ! screening threshold
      real(dp) :: libint_epsilon       = 1.0D-20   ! ε for libint, integral precision given
                                                   ! approximately by sqrt(ε)
!
      type(file) :: restart_file
!
      integer :: n_densities
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: cleanup                                  => cleanup_hf
!
      procedure :: is_restart_safe                          => is_restart_safe_hf
!
      procedure :: read_settings                            => read_settings_hf
      procedure :: read_hf_settings                         => read_hf_settings_hf
      procedure :: construct_ao_overlap                     => construct_ao_overlap_hf
      procedure :: decompose_ao_overlap                     => decompose_ao_overlap_hf
      procedure :: print_energy                             => print_energy_hf
!
!     AO Fock and energy related routines
!
      procedure :: construct_ao_fock                        => construct_ao_fock_hf
      procedure :: construct_ao_G                           => construct_ao_G_hf
      procedure :: construct_ao_G_slow                      => construct_ao_G_slow_hf
      procedure :: construct_ao_G_1der                      => construct_ao_G_1der_hf
      procedure :: construct_ao_G_1der_numerical            => construct_ao_G_1der_numerical_hf
!
      procedure :: construct_coulomb_ao_G                   => construct_coulomb_ao_G_hf
      procedure :: construct_exchange_ao_G                  => construct_exchange_ao_G_hf
!
      procedure :: construct_ao_fock_SAD                    => construct_ao_fock_SAD_hf
      procedure :: construct_mo_fock                        => construct_mo_fock_hf
      procedure :: set_ao_fock                              => set_ao_fock_hf
      procedure :: get_ao_fock                              => get_ao_fock_hf
      procedure :: get_fock_ov                              => get_fock_ov_hf
      procedure :: calculate_hf_energy_from_fock            => calculate_hf_energy_from_fock_hf
      procedure :: calculate_hf_energy_from_G               => calculate_hf_energy_from_G_hf
      procedure :: initialize_fock                          => initialize_fock_hf
      procedure :: destruct_fock                            => destruct_fock_hf
      procedure :: update_fock_and_energy                   => update_fock_and_energy_hf
      procedure :: update_fock_and_energy_cumulative        => update_fock_and_energy_cumulative_hf
!
!     AO Density related routines
!
      procedure :: construct_ao_density                     => construct_ao_density_hf
      procedure :: rotate_ao_density                        => rotate_ao_density_hf
      procedure :: purify_ao_density                        => purify_ao_density_hf
      procedure :: decompose_ao_density                     => decompose_ao_density_hf
      procedure :: get_ao_density                           => get_ao_density_hf
      procedure :: set_ao_density                           => set_ao_density_hf
      procedure :: initialize_density                       => initialize_density_hf
      procedure :: update_ao_density                        => update_ao_density_hf
      procedure :: save_ao_density                          => save_ao_density_hf
      procedure :: get_ao_density_sq                        => get_ao_density_sq_hf
      procedure :: set_initial_ao_density_guess             => set_initial_ao_density_guess_hf
      procedure :: set_ao_density_to_sad                    => set_ao_density_to_sad_hf
      procedure :: set_ao_density_to_core_guess             => set_ao_density_to_core_guess_hf
      procedure :: get_n_electrons_in_density               => get_n_electrons_in_density_hf
      procedure :: construct_sp_density_schwarz             => construct_sp_density_schwarz_hf
!
!     MO orbital related routines
!
      procedure :: do_roothan_hall                          => do_roothan_hall_hf
      procedure :: initialize_orbitals                      => initialize_orbitals_hf
      procedure :: roothan_hall_update_orbitals             => roothan_hall_update_orbitals_hf
      procedure :: print_orbital_energies                   => print_orbital_energies_hf
      procedure :: print_orbitals                           => print_orbitals_hf 
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density                    => initialize_ao_density_hf
      procedure :: initialize_ao_fock                       => initialize_ao_fock_hf
      procedure :: initialize_mo_fock                       => initialize_mo_fock_hf
      procedure :: initialize_ao_overlap                    => initialize_ao_overlap_hf
      procedure :: initialize_pivot_matrix_ao_overlap       => initialize_pivot_matrix_ao_overlap_hf
      procedure :: initialize_cholesky_ao_overlap           => initialize_cholesky_ao_overlap_hf
!
      procedure :: destruct_ao_density                      => destruct_ao_density_hf
      procedure :: destruct_ao_fock                         => destruct_ao_fock_hf
      procedure :: destruct_mo_fock                         => destruct_mo_fock_hf
      procedure :: destruct_ao_overlap                      => destruct_ao_overlap_hf
      procedure :: destruct_pivot_matrix_ao_overlap         => destruct_pivot_matrix_ao_overlap_hf
      procedure :: destruct_cholesky_ao_overlap             => destruct_cholesky_ao_overlap_hf
!
!     Gradient and Hessian related routines
!
      procedure :: construct_projection_matrices            => construct_projection_matrices_hf
      procedure :: project_redundant_rotations              => project_redundant_rotations_hf
!
      procedure :: construct_roothan_hall_hessian           => construct_roothan_hall_hessian_hf
      procedure :: construct_roothan_hall_gradient          => construct_roothan_hall_gradient_hf
      procedure :: get_packed_roothan_hall_gradient         => get_packed_roothan_hall_gradient_hf
!
      procedure :: construct_molecular_gradient             => construct_molecular_gradient_hf
!
!     Integral related routines
!
      procedure :: construct_sp_eri_schwarz                 => construct_sp_eri_schwarz_hf
      procedure :: get_n_sig_eri_sp                         => get_n_sig_eri_sp_hf
!
      procedure :: set_n_mo                                 => set_n_mo_hf
!
      procedure :: set_screening_and_precision_thresholds   => set_screening_and_precision_thresholds_hf
      procedure :: print_screening_settings                 => print_screening_settings_hf
!
   end type hf
!
!
   interface hf 
!
      procedure :: new_hf 
!
   end interface hf 
!
!
contains
!
!
   function new_hf(system) result(wf)
!!
!!    New HF
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      type(hf) :: wf
!
      class(molecular_system), target, intent(in) :: system
!
      wf%name_ = 'HF'
!
      wf%system => system
!
      call wf%read_settings()
!
      wf%n_ao        = wf%system%get_n_aos()
      wf%n_densities = 1
!
      call wf%set_n_mo()
!
      call wf%initialize_wavefunction_files()
      call wf%restart_file%init('hf_restart_file', 'sequential', 'unformatted')
!
      call disk%open_file(wf%restart_file, 'readwrite', 'rewind')
!
      write(wf%restart_file%unit) wf%n_ao
      write(wf%restart_file%unit) wf%n_mo
      write(wf%restart_file%unit) wf%n_densities
      write(wf%restart_file%unit) wf%n_o
      write(wf%restart_file%unit) wf%n_v
!
      call disk%close_file(wf%restart_file)
!
   end function new_hf
!
!
   subroutine is_restart_safe_hf(wf, task)
!!
!!    Is restart safe?
!!    Written by Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(hf) :: wf
!
      character(len=*), intent(in) :: task
!
      integer :: n_ao, n_mo, n_densities
!
      call disk%open_file(wf%restart_file, 'read', 'rewind')
!
      read(wf%restart_file%unit) n_ao
      read(wf%restart_file%unit) n_mo
      read(wf%restart_file%unit) n_densities
!
      call disk%close_file(wf%restart_file)
!
      if (n_ao .ne. wf%n_ao) call output%error_msg('attempted to restart HF with an inconsistent number ' // &
                                                   'of atomic orbitals for task ' // trim(task))
!
      if (n_mo .ne. wf%n_mo) call output%error_msg('attempted to restart HF with an inconsistent number ' // &
                                                   'of molecular orbitals for task ' // trim(task))
!
      if (n_densities .ne. wf%n_densities) call output%error_msg('attempted to restart HF with an inconsistent number ' // &
                                                   'of atomic densities (likely a HF/UHF inconsistency) for task ' // trim(task))
!
   end subroutine is_restart_safe_hf
!
!
   subroutine print_energy_hf(wf)
!!
!!    Print wavefunction summary
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints information related to the wavefunction, most of which is meaningful
!!    only for a properly  converged wavefunction. Should be overwritten in descendants
!!    if more or less or other information is present.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp) :: homo_lumo_gap
!
      homo_lumo_gap = wf%orbital_energies(wf%n_o + 1) - wf%orbital_energies(wf%n_o)
!
      write(output%unit, '(/t6,a26,f19.12)') 'HOMO-LUMO gap:            ', homo_lumo_gap
      write(output%unit, '(t6,a26,f19.12)')  'Nuclear repulsion energy: ', wf%system%get_nuclear_repulsion()
      write(output%unit, '(t6,a26,f19.12)')  'Electronic energy:        ', wf%energy - wf%system%get_nuclear_repulsion()
      write(output%unit, '(t6,a26,f19.12)')  'Total energy:             ', wf%energy
!
   end subroutine print_energy_hf
!
!
   subroutine read_settings_hf(wf)
!!
!!    Read settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Designed to be overwritten by descendants with
!!    descendant specific settings (see e.g. UHF).
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%read_hf_settings()
!
   end subroutine read_settings_hf
!
!
   subroutine read_hf_settings_hf(wf)
!!
!!    Read HF settings
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call input%get_keyword_in_section('coulomb threshold', 'solver hf', wf%coulomb_threshold)
      call input%get_keyword_in_section('exchange threshold', 'solver hf', wf%exchange_threshold)
      call input%get_keyword_in_section('integral precision', 'solver hf', wf%libint_epsilon)
!
   end subroutine read_hf_settings_hf
!
!
   subroutine print_orbital_energies_hf(wf, indentation)
!!
!!    Print orbital energies
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Prints the current orbital energies to output
!!    in a hopefully readable way.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      character(len=*), optional :: indentation
!
      character(len=40) :: indent
!
      indent = '6'
      if (present(indentation)) indent = trim(indentation)
!
      write(output%unit, '(/t'//trim(indent)//',a)') '- Molecular orbital energies:'
!
      call print_vector(wf%orbital_energies, wf%n_ao, indent)
!
   end subroutine print_orbital_energies_hf
!
!
   subroutine print_orbitals_hf(wf, n_orbitals, orbital_list)
!!
!!    Print orbitals 
!!    Written by Eirik F. Kjønstad, May 2019
!!
!!    Prints the orbitals with atom & orbital information given.
!!
      implicit none 
!  
      class(hf), intent(in) :: wf 
!
      integer, intent(in) :: n_orbitals
!
      integer, dimension(n_orbitals), intent(in), optional :: orbital_list
!
      integer, parameter :: n_entries  = 4
!
      integer :: atom, mo_offset, first_mo, last_mo, shell, ao, mo, l, mos_to_print
!
      character(len=3) :: adv ! advance = 'yes' or 'no', depending on the circumstances
!
      character(len=1), dimension(6), parameter :: angular_momentum = ['s', 'p', 'd', 'f', 'g', 'h']
!
      integer, dimension(:), allocatable :: orbital_list_local
!
      mos_to_print = min(wf%n_mo, n_orbitals)
!
      write(output%unit, '(/t3,a,i0,a,i0,a)') 'Printing ', mos_to_print, ' MOs (out of the total ', wf%n_mo, ' MOs).'
!
      call mem%alloc(orbital_list_local, mos_to_print)
!
      if (.not. present(orbital_list)) then
!
!        Default is to print MOs [1,n_orbitals]
!
         do mo = 1, mos_to_print
!
            orbital_list_local(mo) = mo
!
         enddo
!
      else
!
!        If list is provided, print only MOs in list
!
         orbital_list_local = orbital_list
!
!        Sanity check
!
         do mo = 1, mos_to_print
!
            if (orbital_list(mo) .lt. 1 .or. orbital_list(mo) .gt. wf%n_mo) &
                     call output%error_msg('Tried to print non-existent orbital')
!
         enddo
!
      endif
!
      do mo_offset = 1, mos_to_print, n_entries
!
         first_mo = mo_offset
         last_mo  = min(mo_offset + n_entries - 1, mos_to_print)
!
         write(output%unit, '(/t3,a)', advance='no') '- Printing molecular orbitals: '
!
         do mo = first_mo, last_mo 
!
            if (mo == first_mo) then
!
               write(output%unit, '(i0)', advance='no') orbital_list_local(mo) 
!
            elseif ((mo == last_mo) ) then
!
            write(output%unit, '(a,1x,i0/)', advance='yes') ',', orbital_list_local(mo) 
!
            else
!
               write(output%unit, '(a,1x,i0)', advance='no') ',', orbital_list_local(mo)
!
            endif
!
         enddo 
!   
         write(output%unit, '(t6,a)', advance='no') 'AO   Atom'
!
         adv = 'no'
         do mo = first_mo, last_mo 
!
            if (mo == last_mo) adv = 'yes'
!
            if (mo == first_mo) then
!
               write(output%unit, '(7x,i4)', advance=adv) orbital_list_local(mo) 
!
            else
!
               write(output%unit, '(8x,i4)', advance=adv) orbital_list_local(mo)
!
            endif
!
         enddo
!
         write(output%unit, '(t5,a)') '---------------------------------------------------------------'
!
         do atom = 1, wf%system%n_atoms
!
            do shell = 1, wf%system%atoms(atom)%n_shells 
!
               l = wf%system%atoms(atom)%shells(shell)%l
!
               do ao = wf%system%atoms(atom)%shells(shell)%first, wf%system%atoms(atom)%shells(shell)%last 
!
                  write(output%unit, '(t3,i4,i4,1x,a2,a1,a1,a1)', advance='no') ao, atom, &
                                                                          wf%system%atoms(atom)%symbol, &
                                                                         '(', angular_momentum(l + 1), ')'
!
                  adv = 'no'
                  do mo = first_mo, last_mo 
!
                     if (mo == last_mo) adv = 'yes'
                     write(output%unit, '(2x, f10.6)', advance=adv) wf%orbital_coefficients(ao, orbital_list_local(mo))
!
                  enddo
!
               enddo
!
            enddo
!
         enddo
!
         write(output%unit, '(t5,a)') '---------------------------------------------------------------'
!
      enddo 
!
      call mem%dealloc(orbital_list_local, n_orbitals)

!
   end subroutine print_orbitals_hf
!
!
   subroutine mo_transform_hf(wf, X_wx, Y_pq)
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
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: X_wx
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Y_pq
!
      real(dp), dimension(:,:), allocatable :: Z_wq ! = sum_x X_wx C_xq
!
      call mem%alloc(Z_wq, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                 &
                  wf%n_ao,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  X_wx,                    &
                  wf%n_ao,                 &
                  wf%orbital_coefficients, & ! C_xq
                  wf%n_ao,                 &
                  zero,                    &
                  Z_wq,                    &
                  wf%n_ao)
!
      call dgemm('T', 'N',                 &
                  wf%n_mo,                 &
                  wf%n_mo,                 &
                  wf%n_ao,                 &
                  one,                     &
                  wf%orbital_coefficients, & ! C_wp
                  wf%n_ao,                 &
                  Z_wq,                    &
                  wf%n_ao,                 &
                  zero,                    &
                  Y_pq,                    &
                  wf%n_mo)
!
      call mem%dealloc(Z_wq, wf%n_ao, wf%n_mo)
!
   end subroutine mo_transform_hf
!
!
   subroutine set_initial_ao_density_guess_hf(wf, guess)
!!
!!    Set initial AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Sets initial AO density (or densities) to the
!!    appropriate initial guess requested by the solver.
!!
      implicit none
!
      class(hf) :: wf
!
      character(len=*) :: guess
!
      real(dp), dimension(:,:), allocatable :: h_wx
!
      if (trim(guess) == 'sad' .or. trim(guess) == 'SAD') then
!
         call wf%set_ao_density_to_sad()
!
      elseif (trim(guess) == 'core' .or. trim(guess) == 'CORE') then
!
         call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
         call wf%get_ao_h_wx(h_wx)
!
         call wf%set_ao_density_to_core_guess(h_wx)
!
         call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      else
!
         call output%error_msg('Guess AO density ' // trim(guess) // ' is currently not supported.')
!
      endif
!
   end subroutine set_initial_ao_density_guess_hf
!
!
   subroutine initialize_orbitals_hf(wf)
!!
!!    Initialize orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the arrays associated with the orbital
!!    coefficients. In spin-unrestricted wavefunctions, this
!!    will include alpha and beta coefficients, though these
!!    are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      wf%orbital_coefficients = zero
      wf%orbital_energies     = zero
!
   end subroutine initialize_orbitals_hf
!
!
   subroutine initialize_density_hf(wf)
!!
!!    Initialize density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO density (or densities).
!!    In spin-unrestricted wavefunctions, this alpha and beta densities,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_ao_density()
!
      wf%ao_density = zero
!
   end subroutine initialize_density_hf
!
!
   subroutine initialize_fock_hf(wf)
!!
!!    Initialize Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Initializes the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%initialize_ao_fock()
!
      wf%ao_fock = zero
!
   end subroutine initialize_fock_hf
!
!
   subroutine destruct_fock_hf(wf)
!!
!!    Destruct Fock
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Destructs the AO Fock matrix (or matrices).
!!    In spin-unrestricted wavefunctions, this alpha and beta matrices,
!!    though these are the same and therefore redundant in restricted
!!    wavefunctions.
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%destruct_ao_fock()
!
   end subroutine destruct_fock_hf
!
!
   subroutine update_fock_and_energy_hf(wf, h_wx)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      call wf%construct_ao_fock(wf%ao_density, wf%ao_fock, h_wx)
!
      call wf%calculate_hf_energy_from_fock(wf%ao_fock, h_wx)
!
   end subroutine update_fock_and_energy_hf
!
!
   subroutine update_fock_and_energy_cumulative_hf(wf, prev_ao_density, h_wx)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Called by solver when a new density has been obtained and
!!    the next Fock and energy is to be computed. By cumulatively
!!    we mean using the density change to build the Fock matrix
!!    in the iterative loop.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(wf%n_ao**2, wf%n_densities), intent(in) :: prev_ao_density
!
      logical :: cumulative
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
!
      cumulative = .true.
      call wf%construct_ao_fock(wf%ao_density, wf%ao_fock, h_wx, cumulative)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)
!
      call wf%calculate_hf_energy_from_fock(wf%ao_fock, h_wx)
!
   end subroutine update_fock_and_energy_cumulative_hf
!
!
   subroutine roothan_hall_update_orbitals_hf(wf)
!!
!!    Roothan-Hall update of orbitals
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    This routine guides the construction of new orbital coefficients
!!    from the current AO Fock matrix (or matrices if the wavefunction
!!    is unrestricted).
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies)
!
   end subroutine roothan_hall_update_orbitals_hf
!
!
   subroutine update_ao_density_hf(wf)
!!
!!    Update AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Updates the AO density (or densities, if unrestricted) based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%construct_ao_density()
!
   end subroutine update_ao_density_hf
!
!
   subroutine save_ao_density_hf(wf)
!!
!!    Save AO density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Save the AO density (or densities, if unrestricted) based
!!    on the current orbital coefficient matrix (or matrices).
!!
      implicit none
!
      class(hf) :: wf
!
      type(file) :: ao_density
!
      call ao_density%init('ao_density', 'sequential', 'formatted')
      call disk%open_file(ao_density, 'readwrite', 'rewind')
!
      write(ao_density%unit, *) wf%ao_density
!
      call disk%close_file(ao_density)
!
   end subroutine save_ao_density_hf
!
!
   subroutine get_ao_density_sq_hf(wf, D)
!!
!!    Get AO density squared
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Returns the unpacked AO density matrix D
!!    (or density matrices in descendants, see overwriting routines)
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao**2,wf%n_densities), intent(inout) :: D
!
      call dcopy(wf%n_ao**2, wf%ao_density, 1, D, 1)
!
   end subroutine get_ao_density_sq_hf
!
!
   subroutine cleanup_hf(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call disk%open_file(wf%restart_file, 'readwrite', 'append')
!
      write(wf%restart_file%unit) wf%energy
!
      call disk%close_file(wf%restart_file)
!
      call wf%destruct_orbital_energies()
      call wf%destruct_ao_overlap()
      call wf%destruct_ao_fock()
      call wf%destruct_ao_density()
      call wf%destruct_pivot_matrix_ao_overlap()
      call wf%destruct_cholesky_ao_overlap()
!
   end subroutine cleanup_hf
!
!
   subroutine initialize_ao_density_hf(wf)
!!
!!    Initialize AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_density)) call mem%alloc(wf%ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_hf
!
!
   subroutine initialize_ao_fock_hf(wf)
!!
!!    Initialize AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_fock)) call mem%alloc(wf%ao_fock, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_fock_hf
!
!
   subroutine initialize_mo_fock_hf(wf)
!!
!!    Initialize MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%mo_fock)) call mem%alloc(wf%mo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_mo_fock_hf
!
!
   subroutine initialize_ao_overlap_hf(wf)
!!
!!    Initialize AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%ao_overlap)) call mem%alloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_overlap_hf
!
!
   subroutine initialize_pivot_matrix_ao_overlap_hf(wf)
!!
!!    Initialize pivot matrix AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%pivot_matrix_ao_overlap)) call mem%alloc(wf%pivot_matrix_ao_overlap, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_pivot_matrix_ao_overlap_hf
!
!
   subroutine initialize_cholesky_ao_overlap_hf(wf)
!!
!!    Initialize cholesky vectors AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%cholesky_ao_overlap)) call mem%alloc(wf%cholesky_ao_overlap, wf%n_mo, wf%n_mo)
!
   end subroutine initialize_cholesky_ao_overlap_hf
!
!
   subroutine destruct_ao_overlap_hf(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_overlap)) call mem%dealloc(wf%ao_overlap, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_overlap_hf
!
!
   subroutine destruct_ao_density_hf(wf)
!!
!!    Destruct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_density)) call mem%dealloc(wf%ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_hf
!
!
   subroutine destruct_ao_fock_hf(wf)
!!
!!    Destruct AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%ao_fock)) call mem%dealloc(wf%ao_fock, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_hf
!
!
   subroutine destruct_mo_fock_hf(wf)
!!
!!    Destruct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%mo_fock)) call mem%dealloc(wf%mo_fock, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_mo_fock_hf
!
!
   subroutine destruct_pivot_matrix_ao_overlap_hf(wf)
!!
!!    Destruct pivot matrix AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%pivot_matrix_ao_overlap)) call mem%dealloc(wf%pivot_matrix_ao_overlap, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_pivot_matrix_ao_overlap_hf
!
!
   subroutine destruct_cholesky_ao_overlap_hf(wf)
!!
!!    Initialize cholesky vectors AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%cholesky_ao_overlap)) call mem%dealloc(wf%cholesky_ao_overlap, wf%n_mo, wf%n_mo)
!
   end subroutine destruct_cholesky_ao_overlap_hf
!
!
   subroutine construct_ao_density_hf(wf)
!!
!!    Construct AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Routine which calculates D_αβ = sum_i C_αi C_βi,
!!    where C are the MO coefficients.
!!
!!    Density is packed
!!
      implicit none
!
      class(hf) :: wf
!
     call dgemm('N', 'T',                   &
                 wf%n_ao,                   &
                 wf%n_ao,                   &
                 wf%n_o,                    &
                 two,                       &
                 wf%orbital_coefficients,   &
                 wf%n_ao,                   &
                 wf%orbital_coefficients,   &
                 wf%n_ao,                   &
                 zero,                      &
                 wf%ao_density,             &
                 wf%n_ao)
!
   end subroutine construct_ao_density_hf
!
!
   subroutine construct_sp_eri_schwarz_hf(wf)
!!
!!    Construct shell-pair electronic-repulsion-integral Schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of g_wxwx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively.
!!
      implicit none
!
      class(hf) :: wf
!
      integer, dimension(:),  allocatable :: sp_eri_schwarz_index_list
      real(dp), dimension(:), allocatable :: sorted_sp_eri_schwarz
!
!     Local variables
!
      integer :: s1, s2, s1s2
!
      real(dp) :: maximum
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      type(interval) :: A_interval, B_interval
!
!     Set the maximum element in each shell pair
!
      call set_coulomb_precision_c(1.0d-50)
!
!$omp parallel do private(s1, s2, s1s2, A_interval, B_interval, g, maximum) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
            s1s2 = (max(s1,s2)*(max(s1,s2)-3)/2) + s1 + s2
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call wf%system%construct_ao_g_wxyz(g, s1, s2, s1, s2)
!
            maximum = get_abs_max(g, ((A_interval%size)*(B_interval%size))**2)
!
            wf%sp_eri_schwarz(s1s2, 1) = sqrt(maximum)
!
            wf%sp_eri_schwarz_list(s1s2, 1) = s1
            wf%sp_eri_schwarz_list(s1s2, 2) = s2
!
         enddo
      enddo
!$omp end parallel do
!
!     Sort the sp_eri_schwarz vector and use the resulting index list
!     to resort the sp_eri_schwarz_list matrix
!
      call mem%alloc(sp_eri_schwarz_index_list, wf%system%n_s*(wf%system%n_s + 1)/2)
      sp_eri_schwarz_index_list = 0
!
      call mem%alloc(sorted_sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
      sorted_sp_eri_schwarz = wf%sp_eri_schwarz(:, 1)
!
      call get_n_highest(wf%system%n_s*(wf%system%n_s + 1)/2, wf%system%n_s*(wf%system%n_s + 1)/2, wf%sp_eri_schwarz,&
                      sorted_sp_eri_schwarz, sp_eri_schwarz_index_list)
!
      wf%sp_eri_schwarz(:, 2) = wf%sp_eri_schwarz(:, 1)
      wf%sp_eri_schwarz(:, 1) = sorted_sp_eri_schwarz
      call mem%dealloc(sorted_sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
      wf%sp_eri_schwarz_list(:,3) = sp_eri_schwarz_index_list
      call mem%dealloc(sp_eri_schwarz_index_list, wf%system%n_s*(wf%system%n_s + 1)/2)
!
   end subroutine construct_sp_eri_schwarz_hf
!
!
   subroutine construct_sp_density_schwarz_hf(wf, sp_density_schwarz, D)
!!
!!    Construct shell-pair density schwarz vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Computes a vector that contains the largest value (in absolute terms)
!!    of D_wx^1/2 for each shell pair (A,B), where w and x is in A and B,
!!    respectively.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s) :: sp_density_schwarz
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), dimension(:,:), allocatable, target :: D_red
      real(dp), dimension(:,:), contiguous, pointer :: D_red_p => null()
!
      type(interval) :: A_interval, B_interval
!
      integer :: s1, s2
      integer :: n_threads = 1, thread = 0
!
      real(dp) :: maximum
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(D_red,wf%system%max_shell_size**2,n_threads)
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D_red_p, maximum, thread) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
!$          thread = omp_get_thread_num()
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            D_red_p(1:A_interval%size,1:B_interval%size) => D_red(1:A_interval%size*B_interval%size,thread+1)
!
            D_red_p = D(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            maximum = get_abs_max(D_red_p, (A_interval%size)*(B_interval%size))
!
            nullify(D_red_p)
!
            sp_density_schwarz(s1, s2) = maximum
            sp_density_schwarz(s2, s1) = maximum
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(D_red,wf%system%max_shell_size**2,n_threads)
!
   end subroutine construct_sp_density_schwarz_hf
!
!
   subroutine get_n_sig_eri_sp_hf(wf, n_sig_sp)
!!
!!    Get number of significant ERI shell-pairs
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the number of significant shell pairs. The threshold
!!    determines how small the largest element of g_wxwx in a shell
!!    pair AB (w in A, x in B) to be ignored completely in the Fock
!!    construction loop.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(inout) :: n_sig_sp
!
      integer :: s1s2
!
      n_sig_sp = 0
!
      do s1s2 = 1, wf%system%n_s*(wf%system%n_s + 1)/2
!
         if (wf%sp_eri_schwarz(s1s2, 1)*wf%sp_eri_schwarz(1, 1) .lt. sqrt(wf%libint_epsilon)) then
!
            exit
!
         else
!
            n_sig_sp = n_sig_sp + 1
!
         endif
!
      enddo
!
   end subroutine get_n_sig_eri_sp_hf
!
!
   subroutine construct_ao_fock_SAD_hf(wf, coulomb, exchange)
!!
!!    Construct AO Fock matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates
!!
!!       F_αβ = h_αβ + sum_γ g_αβγδ D_γδ - 1/2 * sum_γ g_αδγβ D_γδ,
!!
!!    where D is the SAD density, which is block diagonal.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), optional :: coulomb, exchange ! Non-standard screening thresholds
!
      real(dp) :: coulomb_thr, exchange_thr
!
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz, sp_density_schwarz
!
      real(dp), dimension(:,:), allocatable :: h_wx
      integer :: w, x, y, z
!
      integer :: A, B, C, D, atom
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
      type(interval) :: D_interval
!
      real(dp) :: maximum, max_eri, max_density
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      real(dp), dimension(wf%system%max_shell_size**4), target :: g_C, g_K
      real(dp), dimension(:,:,:,:), pointer :: g_C_p, g_K_p
!
      real(dp), dimension(wf%system%max_shell_size**2), target :: D_yz
      real(dp), dimension(:,:), contiguous, pointer :: D_yz_p
!
      integer, dimension(:,:), allocatable :: shells_on_atoms
!
!     Set thresholds to ignore Coulomb and exchange terms,
!     as well as the desired Libint integral precision
!
      coulomb_thr = wf%coulomb_threshold
      if (present(coulomb)) coulomb_thr = coulomb
!
      exchange_thr = wf%exchange_threshold
      if (present(exchange)) exchange_thr = exchange
!
      !n_s = wf%system%get_n_shells()
!
   call mem%alloc(sp_eri_schwarz, wf%system%n_s, wf%system%n_s)
!
!$omp parallel do private(A, B, A_interval, B_interval, g, maximum) schedule(dynamic)
      do A = 1, wf%system%n_s
         do B = 1, A
!
            A_interval = wf%system%shell_limits(A)
            B_interval = wf%system%shell_limits(B)
!
            call wf%system%construct_ao_g_wxyz(g, A, B, A, B)
!
            maximum = get_abs_max(g, ((A_interval%size)*(B_interval%size))**2)
!
            sp_eri_schwarz(A, B) = sqrt(maximum)
            sp_eri_schwarz(B, A) = sqrt(maximum)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(shells_on_atoms, wf%system%n_atoms, 2) ! [first, last]
!
      shells_on_atoms(1, 1) = 1
      shells_on_atoms(1, 2) = wf%system%atoms(1)%n_shells
!
      do atom = 2, wf%system%n_atoms
!
         shells_on_atoms(atom, 1) =  shells_on_atoms(atom - 1, 2) + 1
         shells_on_atoms(atom, 2) =  shells_on_atoms(atom, 1) + wf%system%atoms(atom)%n_shells - 1
!
      enddo
!
      call mem%alloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
      sp_density_schwarz = zero
!
!$omp parallel do private(A, B, A_interval, B_interval, D_yz, D_yz_p, maximum) schedule(dynamic)
      do atom = 1, wf%system%n_atoms
         do A = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
            A_interval = wf%system%shell_limits(A)
!
            do B = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
               B_interval = wf%system%shell_limits(B)
!
               D_yz_p(1 : A_interval%size, 1 : B_interval%size) => D_yz(1 : A_interval%size*B_interval%size)

               D_yz_p = wf%ao_density(A_interval%first : A_interval%last, &
                                    B_interval%first : B_interval%last)
!
               maximum = get_abs_max(D_yz, (A_interval%size)*(B_interval%size))
!
               sp_density_schwarz(A, B) = maximum
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%ao_fock = zero
!
      max_density = get_abs_max(sp_density_schwarz, wf%system%n_s**2)
      max_eri     = get_abs_max(sp_eri_schwarz, wf%system%n_s**2)
!
!$omp parallel do &
!$omp private(A, B, C, D, A_interval, B_interval, C_interval, D_interval, w, x, y, z, &
!$omp g_C, g_C_p) schedule(dynamic)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            if (sp_eri_schwarz(A, B)*max_eri*max_density .lt. coulomb_thr) cycle
!
            do atom = 1, wf%system%n_atoms
!
               do C = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
                  C_interval = wf%system%shell_limits(C)
!
                  do D = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
                     D_interval = wf%system%shell_limits(D)
!
                     if (sp_eri_schwarz(A, B)*sp_eri_schwarz(C, D)*sp_density_schwarz(C, D) &
                                       .lt. coulomb_thr) cycle
!
                     call wf%system%construct_ao_g_wxyz(g_C, A, B, C, D)

                     g_C_p(1 : A_interval%size, 1 : B_interval%size, 1 : C_interval%size, 1 : D_interval%size) &
                                 => g_C(1 : A_interval%size*B_interval%size*C_interval%size*D_interval%size)
!
!                    Add Fock matrix contributions
!
                     if (A .ne. B) then
!
                        do w = A_interval%first, A_interval%last
                           do x = B_interval%first, B_interval%last
                              do y = C_interval%first, C_interval%last
                                 do z = D_interval%first, D_interval%last
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + g_C_p(w, x, y, z)*wf%ao_density(y, z)
!
                                 enddo
                              enddo
                           enddo
                        enddo
!
                     else
!
                        do w = A_interval%first, A_interval%last
                           do x = A_interval%first, w
                              do y = C_interval%first, C_interval%last
                                 do z = D_interval%first, D_interval%last
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + g_C_p(w, x, y, z)*wf%ao_density(y, z)
!
                                 enddo
                              enddo
                           enddo
                        enddo
                     endif
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do &
!$omp private(A, B, C, D, A_interval, B_interval, C_interval, D_interval, w, x, y, z, &
!$omp g_K, g_K_p) schedule(dynamic)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            do atom = 1, wf%system%n_atoms
!
               do C = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
                  C_interval = wf%system%shell_limits(C)
!
                  do D = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
                     D_interval = wf%system%shell_limits(D)
!
                     if (sp_eri_schwarz(A, D)*sp_eri_schwarz(C, B)*sp_density_schwarz(C, D) .lt. exchange_thr) cycle
!
                     call wf%system%construct_ao_g_wxyz(g_K, A, D, C, B)

                     g_K_p(1 : A_interval%size, 1 : D_interval%size, 1 : C_interval%size, 1 : B_interval%size) &
                                    => g_K(1 : A_interval%size*B_interval%size*C_interval%size*D_interval%size)
!
!                   Add Fock matrix contributions
!
                     if (A .ne. B) then
!
                        do w = A_interval%first, A_interval%last
                           do x = B_interval%first, B_interval%last
                              do y = C_interval%first, C_interval%last
                                 do z = D_interval%first, D_interval%last
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + &
                                             (- half*g_K_p(w, z, y, x))*wf%ao_density(y, z)
!
                                 enddo
!
                              enddo
                           enddo
                        enddo
!
                     else
!
                        do w = A_interval%first, A_interval%last
                           do x = A_interval%first, w
                              do y = C_interval%first, C_interval%last
                                 do z = D_interval%first, D_interval%last
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + (- half*g_K_p(w, z, y, x))*wf%ao_density(y, z)
!
                                 enddo
                              enddo
                           enddo
                        enddo
                     endif
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!     Symmetrize
!
!$omp parallel do &
!$omp private(x, y) schedule(static)
      do x = 1, wf%n_ao
         do y = x + 1, wf%n_ao
!
            wf%ao_fock(x,y) = wf%ao_fock(x,y) + wf%ao_fock(y,x)
!
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do &
!$omp private(x, y) schedule(static)
      do y = 1, wf%n_ao
         do x = y + 1, wf%n_ao
!
            wf%ao_fock(x,y) = wf%ao_fock(y,x)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
!
      call wf%get_ao_h_wx(h_wx)
!
      call wf%calculate_hf_energy_from_G(wf%ao_fock, h_wx)

!$omp parallel do &
!$omp private(x, y) schedule(static)
      do y = 1, wf%n_ao
         do x = 1, wf%n_ao
!
            wf%ao_fock(x,y) = wf%ao_fock(x,y) + h_wx(x, y)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
   end subroutine construct_ao_fock_SAD_hf
!
!
   subroutine construct_ao_fock_hf(wf, D, ao_fock, h_wx, cumulative)
!!
!!    Construct AO Fock matrix
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates
!!
!!       F_αβ = h_αβ + sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where D is the AO density. This routine is integral direct, and
!!    it calculates the Hartree-Fock energy by default.
!!
      implicit none
!
      class(hf) :: wf
!
    !  integer :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: D
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(inout) :: ao_fock
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      logical, intent(in), optional :: cumulative
!
      real(dp) :: coulomb_thr, exchange_thr, precision_thr
!
      integer :: thread = 0, n_threads = 1
!
      logical :: local_cumulative
!
      real(dp), dimension(:,:), allocatable :: F, sp_density_schwarz
!
      integer :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      type(timings) :: ao_fock_timer
!
     ! n_s = wf%system%get_n_shells()
!
      ao_fock_timer = new_timer('AO Fock construction')
      call ao_fock_timer%turn_on()
!
!     Set thresholds to ignore Coulomb and exchange terms,
!     as well as the desired Libint integral precision
!
      coulomb_thr   = wf%coulomb_threshold
      exchange_thr  = wf%exchange_threshold
      precision_thr = wf%libint_epsilon
!
      local_cumulative = .false.
      if (present(cumulative)) then
!
         if (cumulative) then
!
            local_cumulative = .true.
!
         else
!
            local_cumulative = .false.
!
         endif
!
      endif
!
!     Construct the density screening vector and the maximum element in the density
!
      call mem%alloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
      call wf%construct_sp_density_schwarz(sp_density_schwarz, D)
      max_D_schwarz = get_abs_max(sp_density_schwarz, wf%system%n_s**2)
!
!     Compute number of significant ERI shell pairs (the Fock construction
!     only loops over these shell pairs) and the maximum element
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
      max_eri_schwarz = get_abs_max(wf%sp_eri_schwarz, wf%system%n_s*(wf%system%n_s + 1)/2)
!
!     Construct the two electron part of the Fock matrix, using the screening vectors
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
!$    n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      F = zero
!
      call wf%construct_ao_G(F, D, n_threads, max_D_schwarz, max_eri_schwarz,     &
                                         sp_density_schwarz,  &
                                         n_sig_sp, coulomb_thr, exchange_thr, precision_thr, &
                                         wf%system%shell_limits)
!
      call mem%dealloc(sp_density_schwarz, wf%system%n_s, wf%system%n_s)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix,
!     and symmetrize the result
!
      if (.not. local_cumulative) ao_fock = zero
      do thread = 1, n_threads
!
         call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, ao_fock, 1)
!
      enddo
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads)
!
      call symmetric_sum(ao_fock, wf%n_ao)
      ao_fock = ao_fock*half
!
      if (.not. local_cumulative) ao_fock = ao_fock + h_wx
!
      call ao_fock_timer%turn_off()
!
   end subroutine construct_ao_fock_hf
!
!
   subroutine construct_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,    &
                                          sp_density_schwarz, n_sig_sp, coulomb_thr, &
                                          exchange_thr, precision_thr, shells)
!!
!!    Construct AO G 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ = sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ (= G(D)_αβ),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads, n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                             &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34,     &
!$omp w, x, y, z, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset,         &
!$omp temp4, temp5, temp6, temp7, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,        &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                                &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4), &
                           sp_density_schwarz_s1s2)
!
               temp8 = max(sp_density_schwarz_s3s2, &
                           sp_density_schwarz_s3s1, &
                           sp_density_schwarz(s4,s2), &
                           sp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr .and. temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4,         &
                  precision_thr/max(temp7,temp8), thread, skip, shells(s1)%size, shells(s2)%size, &
                  shells(s3)%size, shells(s4)%size)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
!              Add Fock matrix contributions
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%size*(shells(s2)%size*(shells(s3)%size*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           temp3 = one_over_eight*temp*d3
                           temp4 = one_over_eight*temp*d4
                           temp5 = one_over_eight*temp*d5
                           temp6 = one_over_eight*temp*d6
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
!
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
!
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call set_coulomb_precision_c(wf%libint_epsilon)
!
   end subroutine construct_ao_G_hf
!
!
   subroutine construct_ao_G_slow_hf(wf, F, D, shells)
!!
!!    Construct AO G slow 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ = sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ (= G(D)_αβ),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
      real(dp), dimension(wf%n_ao, wf%n_ao)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp) :: d1, d2, d3, d4, d5, d6
      real(dp) :: deg, deg_12, deg_34, deg_12_34, temp1, temp2, temp3, temp4, temp5, temp6, temp
!
      integer :: w, x, y, z, s1, s2, s3, s4, s4_max, tot_dim
      integer :: w_red, x_red, y_red, z_red, wxyz
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
      do s1 = 1, wf%system%n_s 
         do s2 = 1, s1 
!
            deg_12 = real(2-s2/s1, kind=dp)
!
            do s3 = 1, s1
!
               s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
               do s4 = 1, s4_max
!
                  deg_34    = real(2-s4/s3, kind=dp)
                  deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

                  deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
                  call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4,         &
                  wf%libint_epsilon, thread, skip, shells(s1)%size, shells(s2)%size, &
                  shells(s3)%size, shells(s4)%size)
!
                  if (skip == 1) cycle
!
                  tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
                  g(1:tot_dim) = deg*g(1:tot_dim)
!
                  do z = shells(s4)%first, shells(s4)%last
!
                     z_red = z - shells(s4)%first + 1
!
                     do y = shells(s3)%first, shells(s3)%last
!
                        y_red = y - shells(s3)%first + 1
!
                        d1 = D(y, z)
!
                        do x = shells(s2)%first, shells(s2)%last
!
                           x_red = x - shells(s2)%first + 1
!
                           d3 = D(x, y)
                           d5 = D(x, z)
!
                           do w = shells(s1)%first, shells(s1)%last
!
                              d2 = D(w, x)
                              d4 = D(w, y)
                              d6 = D(w, z)
!
                              w_red = w - shells(s1)%first + 1
!
                              wxyz = shells(s1)%size*(shells(s2)%size*(shells(s3)%size*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                              temp = g(wxyz)
!
                              temp1 = half*temp*d1
                              temp2 = half*temp*d2
!
                              temp3 = one_over_eight*temp*d3
                              temp4 = one_over_eight*temp*d4
                              temp5 = one_over_eight*temp*d5
                              temp6 = one_over_eight*temp*d6
!
                              F(w, x) = F(w, x) + temp1
                              F(y, x) = F(y, x) - temp6
!
                              F(y, z) = F(y, z) + temp2
                              F(w, z) = F(w, z) - temp3
                              F(x, z) = F(x, z) - temp4
!
                              F(w, y) = F(w, y) - temp5
!
                           enddo
                        enddo
                     enddo
                  enddo
!
               enddo
            enddo
         enddo
      enddo
!
   end subroutine construct_ao_G_slow_hf
!
!
   subroutine construct_coulomb_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,     &
                                                   sp_density_schwarz, &
                                                   n_sig_sp, coulomb_thr, precision_thr, shells)
!!
!!    AO Fock Coulomb construction loop
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the Coulomb two-electron part of the Fock matrix,
!!
!!       F_αβ = F_αβ + sum_γδ g_αβγδ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads,  n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: d1, d2, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp7, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp1, temp2, d1, d2, thread, thread_offset,                            &
!$omp temp7, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,                                &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                            &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. coulomb_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. coulomb_thr) cycle ! Screened out shell pair
!
               temp7 = max(sp_density_schwarz(s3,s4), sp_density_schwarz_s1s2)
!
               if (temp7*temp .lt. coulomb_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp7, thread, skip, shells(s1)%size, shells(s2)%size,    &
                  shells(s3)%size, shells(s4)%size)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
!              Add Fock matrix contributions
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     d1 = D(y, z)
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d2 = D(w, x)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%size*(shells(s2)%size*(shells(s3)%size*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp1 = half*temp*d1
                           temp2 = half*temp*d2
!
                           F(w, thread_offset + x) = F(w, thread_offset + x) + temp1
                           F(y, thread_offset + z) = F(y, thread_offset + z) + temp2
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_coulomb_ao_G_hf
!
!
   subroutine construct_exchange_ao_G_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz, &
                                          sp_density_schwarz, &
                                           n_sig_sp, exchange_thr, precision_thr, shells)
!!
!!    AO Fock exchange construction loop
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβ = F_αβ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      integer, intent(in) :: n_threads,  n_sig_sp
!
      type(interval), dimension(wf%system%n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, exchange_thr, precision_thr
!
      real(dp), dimension(wf%system%n_s, wf%system%n_s), intent(in) :: sp_density_schwarz
!
      real(dp) :: d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp3, temp4, temp5, temp6, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer :: w, x, y, z, s1s2, s1, s2, s3, s4, s4_max, tot_dim
      integer :: s3s4, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(wf%system%max_shell_size**4) :: g
!
      integer :: thread = 0, skip
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, temp3, d3, d4, d5, d6, thread, thread_offset,                           &
!$omp temp4, temp5, temp6, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,           &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, deg_12_34,                            &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
!$       thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix
!
         sp_eri_schwarz_s1s2 = wf%sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = wf%sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = wf%sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = wf%sp_eri_schwarz_list(s1s2_packed, 2)
!
         sp_density_schwarz_s1s2 = sp_density_schwarz(s1, s2)
         if (sp_eri_schwarz_s1s2*(max_D_schwarz)*(max_eri_schwarz) .lt. exchange_thr) cycle
!
         deg_12 = real(2-s2/s1, kind=dp)
!
         do s3 = 1, s1
!
            sp_density_schwarz_s3s2 = sp_density_schwarz(s3, s2)
            sp_density_schwarz_s3s1 = sp_density_schwarz(s3, s1)
!
            s4_max = (s3/s1)*s2 + (1-s3/s1)*s3
!
            do s4 = 1, s4_max
!
               s3s4 = (max(s3,s4)*(max(s3,s4)-3)/2) + s3 + s4
!
               temp = sp_eri_schwarz_s1s2*wf%sp_eri_schwarz(s3s4, 2)

               if (temp*(max_D_schwarz) .lt. exchange_thr) cycle ! Screened out shell pair
!
               temp8 = max(sp_density_schwarz_s3s2,   &
                           sp_density_schwarz_s3s1,   &
                           sp_density_schwarz(s4,s2), &
                           sp_density_schwarz(s1,s4))
!
               if (temp8*temp .lt. exchange_thr) cycle
!
               deg_34    = real(2-s4/s3, kind=dp)
               deg_12_34 = min(1-s3/s1+2-min(s4/s2,s2/s4), 2)

               deg = deg_12*deg_34*deg_12_34 ! Shell degeneracy
!
               call wf%system%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp8, thread, skip, shells(s1)%size, shells(s2)%size,    &
                  shells(s3)%size, shells(s4)%size)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
               g(1:tot_dim) = deg*g(1:tot_dim)
!
!              Add Fock matrix contributions
!
               do z = shells(s4)%first, shells(s4)%last
!
                  z_red = z - shells(s4)%first + 1
!
                  do y = shells(s3)%first, shells(s3)%last
!
                     y_red = y - shells(s3)%first + 1
!
                     do x = shells(s2)%first, shells(s2)%last
!
                        x_red = x - shells(s2)%first + 1
!
                        d3 = D(x, y)
                        d5 = D(x, z)
!
                        do w = shells(s1)%first, shells(s1)%last
!
                           d4 = D(w, y)
                           d6 = D(w, z)
!
                           w_red = w - shells(s1)%first + 1
!
                           wxyz = shells(s1)%size*(shells(s2)%size*(shells(s3)%size*(z_red-1)+y_red-1)+x_red-1)+w_red
!
                           temp = g(wxyz)
!
                           temp3 = one_over_eight*temp*d3
                           temp4 = one_over_eight*temp*d4
                           temp5 = one_over_eight*temp*d5
                           temp6 = one_over_eight*temp*d6
!
                           F(y, thread_offset + x) = F(y, thread_offset + x) - temp6
                           F(w, thread_offset + z) = F(w, thread_offset + z) - temp3
                           F(x, thread_offset + z) = F(x, thread_offset + z) - temp4
                           F(w, thread_offset + y) = F(w, thread_offset + y) - temp5
!
                        enddo
                     enddo
                  enddo
               enddo
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_exchange_ao_G_hf
!
!
   subroutine calculate_hf_energy_from_G_hf(wf, half_GD_wx, h_wx)
!!
!!    Calculate HF energy from G(D)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates the Hartree-Fock energy,
!!
!!       E = Tr(h D) + 1/4 * Tr(D G(D)),
!!
!!    where D is the AO density and
!!
!!       G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ
!!
!!    The traces are calculated as dot products (since B is symmetric here):
!!
!!       Tr(AB) = sum_x (AB)_xx = sum_xy A_xy B_yx = sum_xy A_xy B_xy.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: half_GD_wx
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      wf%energy = wf%system%get_nuclear_repulsion()
!
      wf%energy = wf%energy + ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      wf%energy = wf%energy + half*ddot((wf%n_ao)**2, wf%ao_density, 1, half_GD_wx, 1)

   end subroutine calculate_hf_energy_from_G_hf
!
!
   subroutine calculate_hf_energy_from_fock_hf(wf, F_wx, h_wx)
!!
!!    Calculate HF energy from F
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Calculates the Hartree-Fock energy,
!!
!!       E = Tr(h D) + 1/4 * Tr(D G(D)) = 1/2 * Tr (F D) + 1/2 * Tr(h D),
!!
!!    where D is the AO density and
!!
!!       G(D)_αβ = 2 sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ
!!
!!    The traces are calculated as dot products (since B is symmetric here):
!!
!!       Tr(AB) = sum_x (AB)_xx = sum_xy A_xy B_yx = sum_xy A_xy B_xy = A-dot-B.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp) :: ddot
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: F_wx
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      wf%energy = wf%system%get_nuclear_repulsion()
!
      wf%energy = wf%energy + half*ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      wf%energy = wf%energy + half*ddot((wf%n_ao)**2, wf%ao_density, 1, F_wx, 1)

   end subroutine calculate_hf_energy_from_fock_hf
!
!
   subroutine construct_ao_G_1der_numerical_hf(wf, G_wxqk, dx)
!!
!!    Get AO G(D) 1st derivative numerically
!!    Written by Eirik F. Kjønstad, June 2019 
!!
!!    Uses forward differences to calculate the derivative.
!!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: G_wxqk
!
      real(dp), intent(in) :: dx ! Displacement length in bohr
!
      real(dp), dimension(:,:), allocatable :: G_wx, G_wx_displaced, R_qk, R_qk_displaced
!
      integer :: k, q, w, x 
!
!     Get s at the reference geometry 
!
      call mem%alloc(G_wx, wf%n_ao, wf%n_ao)
      G_wx = zero
      call wf%construct_ao_G_slow(G_wx, wf%ao_density, wf%system%shell_limits)
!
      call symmetric_sum(G_wx, wf%n_ao)
      G_wx = half*G_wx 
!
      call mem%alloc(R_qk, 3, wf%system%n_atoms)
      call mem%alloc(R_qk_displaced, 3, wf%system%n_atoms)
!
      R_qk = wf%system%get_geometry()
      call mem%alloc(G_wx_displaced, wf%n_ao, wf%n_ao)
!
      do k = 1, wf%system%n_atoms
         do q = 1, 3
!
!           Get s at the displaced geometry 
!
            R_qk_displaced = R_qk 
            R_qk_displaced(q,k) = R_qk_displaced(q,k) + dx 
            call wf%system%set_geometry(R_qk_displaced)
!
            G_wx_displaced = zero 
            call wf%construct_ao_G_slow(G_wx_displaced, wf%ao_density, wf%system%shell_limits)
!
            call symmetric_sum(G_wx_displaced, wf%n_ao)
            G_wx_displaced = half*G_wx_displaced 
!
!           Use difference from reference to compute derivative
!
            do w = 1, wf%n_ao
               do x = 1, wf%n_ao 
!
                  G_wxqk(w,x,q,k) = (G_wx_displaced(w,x) - G_wx(w,x))/dx
!
               enddo
            enddo
!
         enddo
      enddo
!    
      call wf%system%set_geometry(R_qk)
!
      call mem%dealloc(G_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(R_qk, 3, wf%system%n_atoms)
      call mem%dealloc(R_qk_displaced, 3, wf%system%n_atoms)
      call mem%dealloc(G_wx_displaced, wf%n_ao, wf%n_ao)
!
   end subroutine construct_ao_G_1der_numerical_hf
!
!
   subroutine construct_ao_G_1der_hf(wf, G_ao, D_ao)
!!
!!    Construct AO G 1der
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    This routine constructs the entire two-electron part of the Fock matrix,
!!
!!       F_αβqk = sum_γδ g_αβγδqk D_γδ - 1/2 * sum_γδ g_αδγβqk D_γδ (= G(D)_αβqk),
!!
!!    where contributions from different threads are gathered column blocks
!!    of the incoming F matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao, 3, wf%system%n_atoms) :: G_ao
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_ao
!
      real(dp), dimension((wf%system%max_shell_size**4)*3*4), target :: g_ABCDqk 
      real(dp), dimension(:,:,:,:,:,:), pointer, contiguous :: g_ABCDqk_p 
!
      integer :: A, B, C, D, D_max, w, x, y, z, n_sig_sp, AB, AB_packed
      integer :: w_red, x_red, y_red, z_red, tot_dim, k, q, n_threads, thread
!
      real(dp) :: d1, d2, d3, d4, d5, d6
!
      integer, dimension(4) :: atoms 
!
      real(dp) :: deg, deg_AB, deg_CD, deg_AB_CD 
!
      real(dp), dimension(3,4) :: temp, temp1, temp2, temp3, temp4, temp5, temp6
!
      real(dp), dimension(:,:,:,:,:), allocatable :: G_ao_t
!
!$    n_threads = omp_get_max_threads()
      call mem%alloc(G_ao_t, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms, n_threads)
      G_ao_t = zero 
!
      call wf%get_n_sig_eri_sp(n_sig_sp)
!
!$omp parallel do private(A, B, C, D, D_max, atoms, deg, deg_CD, deg_AB, deg_AB_CD, g_ABCDqk, g_ABCDqk_p, &
!$omp w, x, y, z, w_red, x_red, y_red, z_red, temp, temp1, temp2, temp3, temp4, temp5, temp6, &
!$omp d1, d2, d3, d4, d5, d6, thread, q, k, tot_dim, AB, AB_packed) schedule(dynamic)
      do AB = 1, n_sig_sp 
!
!$       thread = omp_get_thread_num()
!
         if (wf%sp_eri_schwarz(AB, 1)*wf%sp_eri_schwarz(1, 1) < wf%coulomb_threshold) cycle
         AB_packed = wf%sp_eri_schwarz_list(AB, 3) 
!
         A = wf%sp_eri_schwarz_list(AB_packed, 1)
         B = wf%sp_eri_schwarz_list(AB_packed, 2)
!
         atoms(1) = wf%system%shell2atom(A)
         atoms(2) = wf%system%shell2atom(B)
!
         deg_AB = real(2-B/A, kind=dp)
!
            do C = 1, A
!
               D_max = (C/A)*B + (1-C/A)*C
               atoms(3) = wf%system%shell2atom(C)
!
               do D = 1, D_max
!
                  deg_CD    = real(2-D/C, kind=dp)
                  deg_AB_CD = min(1-C/A+2-min(D/B,B/D), 2)

                  deg = deg_AB*deg_CD*deg_AB_CD ! Shell degeneracy
!
                  atoms(4) = wf%system%shell2atom(D)
!
                  call wf%system%construct_ao_g_wxyz_1der(g_ABCDqk, A, B, C, D)                 
!
                  tot_dim = (wf%system%shell_limits(A)%size)*(wf%system%shell_limits(B)%size)&
                              *(wf%system%shell_limits(C)%size)*(wf%system%shell_limits(D)%size)&
                              *3*4
!
                  g_ABCDqk(1:tot_dim) = deg*g_ABCDqk(1:tot_dim)
!
                  g_ABCDqk_p(1 : wf%system%shell_limits(A)%size, 1 : wf%system%shell_limits(B)%size, &
                             1 : wf%system%shell_limits(C)%size, 1 : wf%system%shell_limits(D)%size, &
                             1 : 3, 1 : 4) => g_ABCDqk(1 : tot_dim)
!
                  do z = wf%system%shell_limits(D)%first, wf%system%shell_limits(D)%last
!
                     z_red = z - wf%system%shell_limits(D)%first + 1
!
                     do y = wf%system%shell_limits(C)%first, wf%system%shell_limits(C)%last
!
                        y_red = y - wf%system%shell_limits(C)%first + 1
!
                        d1 = D_ao(y, z)
!
                        do x = wf%system%shell_limits(B)%first, wf%system%shell_limits(B)%last
!
                           x_red = x - wf%system%shell_limits(B)%first + 1
!
                           d3 = D_ao(x, y)
                           d5 = D_ao(x, z)
!
                           do w = wf%system%shell_limits(A)%first, wf%system%shell_limits(A)%last
!
                              d2 = D_ao(w, x)
                              d4 = D_ao(w, y)
                              d6 = D_ao(w, z)
!
                              w_red = w - wf%system%shell_limits(A)%first + 1
!
                              do k = 1, 4
                                 do q = 1, 3
                                    temp(q, k) = g_ABCDqk_p(w_red, x_red, y_red, z_red, q, k)
                                 enddo
                              enddo
!
                              temp1 = half*temp*d1
                              temp2 = half*temp*d2
!
                              temp3 = one_over_eight*temp*d3
                              temp4 = one_over_eight*temp*d4
                              temp5 = one_over_eight*temp*d5
                              temp6 = one_over_eight*temp*d6
!
                              do k = 1, 4
                                 do q = 1, 3
!
                                    G_ao_t(w, x, q, atoms(k), thread+1) = G_ao_t(w, x, q, atoms(k), thread+1) + temp1(q, k)
                                    G_ao_t(y, x, q, atoms(k), thread+1) = G_ao_t(y, x, q, atoms(k), thread+1) - temp6(q, k)
                                    G_ao_t(y, z, q, atoms(k), thread+1) = G_ao_t(y, z, q, atoms(k), thread+1) + temp2(q, k)
                                    G_ao_t(w, z, q, atoms(k), thread+1) = G_ao_t(w, z, q, atoms(k), thread+1) - temp3(q, k)
                                    G_ao_t(x, z, q, atoms(k), thread+1) = G_ao_t(x, z, q, atoms(k), thread+1) - temp4(q, k)
                                    G_ao_t(w, y, q, atoms(k), thread+1) = G_ao_t(w, y, q, atoms(k), thread+1) - temp5(q, k)
!
                                 enddo
                              enddo
!
                           enddo
                        enddo
                     enddo
                  enddo
!
               enddo
            enddo
         enddo
!$omp end parallel do 
!
      G_ao = zero
!
      do thread = 1, n_threads
!
         call daxpy(3*wf%system%n_atoms*wf%n_ao**2, one, G_ao_t(1,1,1,1,thread), 1, G_ao, 1)
!
      enddo
!
      call mem%dealloc(G_ao_t, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms, n_threads)
!
   end subroutine construct_ao_G_1der_hf
!
!
   subroutine set_ao_density_hf(wf, D)
!!
!!    Set AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO density from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(:,:) :: D ! Packed
!
      call squareup(D, wf%ao_density, wf%n_ao)
!
   end subroutine set_ao_density_hf
!
!
   subroutine set_ao_fock_hf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock from input
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao + 1)/2, wf%n_densities), intent(in) :: F ! Packed
!
      call squareup(F, wf%ao_fock, wf%n_ao)
!
   end subroutine set_ao_fock_hf
!
!
   subroutine get_ao_fock_hf(wf, F)
!!
!!    Set AO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the AO Fock from input
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao+1)/2, wf%n_densities), intent(inout) :: F ! Packed
!
      call packin(F, wf%ao_fock, wf%n_ao)
!
   end subroutine get_ao_fock_hf
!
!
   subroutine construct_ao_overlap_hf(wf)
!!
!!    Construct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      call wf%get_ao_s_wx(wf%ao_overlap)
!
   end subroutine construct_ao_overlap_hf
!
!
   subroutine construct_mo_fock_hf(wf)
!!
!!    Construct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the MO Fock matrix F_pq using the current AO
!!    Fock and the orbital coefficients.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_mo)
!
      call dgemm('N', 'N',                   &
                  wf%n_ao,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%ao_fock,                &
                  wf%n_ao,                   &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  zero,                      &
                  X,                         & ! X = F^ao C
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_mo,                   &
                  wf%n_mo,                   &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  wf%mo_fock,                & ! F = C^T F^ao C
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_ao, wf%n_mo)
!
   end subroutine construct_mo_fock_hf
!
!
   subroutine get_fock_ov_hf(wf, F)
!!
!!    Get HF equations
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the occupied-virtual block of the Fock MO matrix,
!!    and returns the result in the array F.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_o, wf%n_v)   :: F ! F_ia
      real(dp), dimension(:,:), allocatable :: X
!
      call mem%alloc(X, wf%n_ao, wf%n_v)
!
      call dgemm('N', 'N',                                  &
                  wf%n_ao,                                  &
                  wf%n_v,                                   &
                  wf%n_ao,                                  &
                  one,                                      &
                  wf%ao_fock,                               &
                  wf%n_ao,                                  &
                  wf%orbital_coefficients(1, wf%n_o + 1),   &
                  wf%n_ao,                                  &
                  zero,                                     &
                  X,                                        &
                  wf%n_ao)
!
      call dgemm('T', 'N',                   &
                  wf%n_o,                    &
                  wf%n_v,                    &
                  wf%n_ao,                   &
                  one,                       &
                  wf%orbital_coefficients,   &
                  wf%n_ao,                   &
                  X,                         &
                  wf%n_ao,                   &
                  zero,                      &
                  F,                         &
                  wf%n_o)
!
      call mem%dealloc(X, wf%n_ao, wf%n_v)
!
   end subroutine get_fock_ov_hf
!
!
   subroutine get_ao_density_hf(wf, D)
!!
!!    Get AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Packs the AO density into D.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(:,:) :: D
!
      call packin(D, wf%ao_density, wf%n_ao)
!
   end subroutine get_ao_density_hf
!
!
   subroutine decompose_ao_density_hf(wf)
!!
!!    Decompose AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Does a Cholesky decomposition of the AO density matrix,
!!
!!       D^AO_wx = sum_J L_w,J L_J,x^T,
!!
!!    and sets the MO coefficients accordingly.
!!
      implicit none
!
      class(hf) :: wf
!
      integer, dimension(:), allocatable :: used_diag
!
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      integer :: rank, j
!
      call mem%alloc(used_diag, wf%n_ao)
!
      wf%ao_density = half*wf%ao_density
      call full_cholesky_decomposition_system(wf%ao_density, wf%orbital_coefficients, wf%n_ao, rank, &
                                                      1.0d-12, used_diag)
      wf%ao_density = two*wf%ao_density
!
!     Make permutation matrix P
!
      call mem%alloc(perm_matrix, wf%n_ao, wf%n_ao)
!
      perm_matrix = zero
!
      do j = 1, wf%n_ao
!
         perm_matrix(used_diag(j), j) = one
!
      enddo
!
      call mem%dealloc(used_diag, wf%n_ao)
!
   end subroutine decompose_ao_density_hf
!
!
   subroutine decompose_ao_overlap_hf(wf)
!!
!!    Decompose AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs a Cholesky decomposition of the AO overlap matrix S,
!!    to within a given threshold:
!!
!!       P^T S P = L L^T
!!
!!    The routine allocates and sets P and L, the 'permutation matrix'
!!    and the 'cholesky ao overlap'. Moreover, it sets the number of
!!    linearly independent AOs, wf%n_mo, which is sometimes less than
!!    wf%n_ao. From P and L, we can transform equations to the linearly
!!    independent basis and back.
!!
      implicit none
!
      class(hf) :: wf
!
      integer, dimension(:), allocatable :: used_diag
!
      real(dp), dimension(:, :), allocatable :: L
!
      integer :: j
!
      call mem%alloc(used_diag, wf%n_ao)
      used_diag = 0
!
      call mem%alloc(L, wf%n_ao, wf%n_ao) ! Full Cholesky vector
      L = zero
!
      call full_cholesky_decomposition_system(wf%ao_overlap, L, wf%n_ao, wf%n_mo, &
                                                wf%linear_dep_threshold, used_diag)
!
      call wf%initialize_cholesky_ao_overlap()
!
      wf%cholesky_ao_overlap(:,:) = L(1:wf%n_mo, 1:wf%n_mo)
!
      call mem%dealloc(L, wf%n_ao, wf%n_ao)
!
!     Make permutation matrix P
!
      call wf%initialize_pivot_matrix_ao_overlap()
!
      wf%pivot_matrix_ao_overlap = zero
!
      do j = 1, wf%n_mo
!
         wf%pivot_matrix_ao_overlap(used_diag(j), j) = one
!
      enddo
!
      call mem%dealloc(used_diag, wf%n_ao)
!
   end subroutine decompose_ao_overlap_hf
!
!
   subroutine rotate_ao_density_hf(wf, X)
!!
!!    Rotate AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Performs an update of the AO density according to a first-order
!!    truncation of the BCH expansion:
!!
!!       D^AO <- exp(-X S) D^AO exp(S X) ~ D^AO + [D^AO, X]_S = D_AO + C.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: X
!
      real(dp), dimension(:,:), allocatable :: M
      real(dp), dimension(:,:), allocatable :: C
!
!     Construct C = [D^AO, X]_S = D^AO S X - X S D^AO
!
      call mem%alloc(M, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  X,             &
                  wf%n_ao,       &
                  zero,          &
                  M,             & ! M = S X
                  wf%n_ao)
!
      call mem%alloc(C, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  M,             &
                  wf%n_ao,       &
                  zero,          &
                  C,             & ! C = D^AO M = D^AO S X
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  wf%ao_density, &
                  wf%n_ao,       &
                  zero,          &
                  M,             & ! M = S D^AO
                  wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  -one,          &
                  X,             &
                  wf%n_ao,       &
                  M,             &
                  wf%n_ao,       &
                  one,           &
                  C,             & ! C = C - X M = C - X S D^AO = D^AO S X - X S D^AO
                  wf%n_ao)
!
      wf%ao_density = wf%ao_density + C
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
      call mem%dealloc(C, wf%n_ao, wf%n_ao)
!
   end subroutine rotate_ao_density_hf
!
!
   subroutine purify_ao_density_hf(wf, threshold)
!!
!!    Purify AO density
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Purifies a non-idempotent AO density matrix - typically arising from
!!    the non-exact rotation of the density - by the following fixed point algorithm:
!!
!!       D^AO <- 3/4 D^AO S D^AO - 1/2 D^AO S D^AO S D^AO
!!
!!    To check whether the final result is consistent, it is possible to verify that
!!    1/2 D^AO S is idempotent (which was done during debug of routine).
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(:,:), allocatable :: M ! Arrays to temporarily hold matrix products
      real(dp), dimension(:,:), allocatable :: N ! Arrays to temporarily hold matrix products
!
      real(dp), dimension(:,:), allocatable :: prev_ao_density ! Holds previous density matrix
!
      real(dp) :: ddot, error
!
      logical :: pure = .false.
!
      integer :: iteration
      integer, parameter :: max_iterations = 50
!
      iteration = 1
!
      call mem%alloc(M, wf%n_ao, wf%n_ao)
      call mem%alloc(N, wf%n_ao, wf%n_ao)
      call mem%alloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
      pure = .false.
      error = zero
!
      do while (.not. pure .and. iteration .le. max_iterations)
!
         prev_ao_density = wf%ao_density
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_overlap, &
                     wf%n_ao,       &
                     wf%ao_density, &
                     wf%n_ao,       &
                     zero,          &
                     M,             & ! M = S D^AO
                     wf%n_ao)
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_density, &
                     wf%n_ao,       &
                     M,             &
                     wf%n_ao,       &
                     zero,          &
                     N,             & ! N = D^AO M = D^AO S D^AO
                     wf%n_ao)
!
         wf%ao_density = (three/two)*N ! D^AO = 3 D^AO S D^AO
!
         call dgemm('N','N',        &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     wf%n_ao,       &
                     one,           &
                     wf%ao_overlap, &
                     wf%n_ao,       &
                     N,             &
                     wf%n_ao,       &
                     zero,          &
                     M,             & ! M = S N = S D^AO S D^AO
                     wf%n_ao)
!
         call dgemm('N','N',          &
                     wf%n_ao,         &
                     wf%n_ao,         &
                     wf%n_ao,         &
                     one,             &
                     prev_ao_density, &
                     wf%n_ao,         &
                     M,               &
                     wf%n_ao,         &
                     zero,            &
                     N,               & ! N = D^AO M = D^AO S D^AO S D^AO
                     wf%n_ao)
!
         wf%ao_density = wf%ao_density - (one/two)*N
!
         M = wf%ao_density - prev_ao_density
!
         error = sqrt(ddot((wf%n_ao)**2, M, 1, M, 1))
!
         if (error .lt. threshold) then
!
            pure = .true.
!
         endif
!
         iteration = iteration + 1
!
      enddo
!
      if (.not. pure) then
!
         write(output%unit, '(t3,a49,f16.12)') 'Error: could not purify AO density. Final error: ', error
         stop
!
      endif
!
      call mem%dealloc(M, wf%n_ao, wf%n_ao)
      call mem%dealloc(N, wf%n_ao, wf%n_ao)
      call mem%dealloc(prev_ao_density, wf%n_ao, wf%n_ao)
!
   end subroutine purify_ao_density_hf
!
!
   subroutine construct_projection_matrices_hf(wf, Po, Pv, D)
!!
!!    Construct projection matrices
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs Po = 1/2 D S
!!               Pv = 1 - Po
!!
!!    D: AO density matrix, S: AO overlap matrix. Both are assumed
!!    to be allocated and properly set.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao) :: Pv
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      integer :: w, x
!
!     Po = 1/2 D S
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  half,          &
                  D,             &
                  wf%n_ao,       &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  zero,          &
                  Po,            &
                  wf%n_ao)
!
!     Pv = I
!
      Pv = zero
!
      do x = 1, wf%n_ao
!
         Pv(x, x) = one
!
      enddo
!
!     Pv = I - Po
!
      do x = 1, wf%n_ao
         do w = 1, wf%n_ao
!
            Pv(w, x) = Pv(w, x) - Po(w, x)
!
         enddo
      enddo
!
   end subroutine construct_projection_matrices_hf
!
!
   subroutine project_redundant_rotations_hf(wf, X, Po, Pv)
!!
!!    Project redundant rotations
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Here, X is an antisymmetric rotations matrix on entering the routine,
!!    where some parameters are redundant for rotating the AO density. On exit,
!!    the redundant parameters have been projected out of X.
!!
!!    To achieve this, we set
!!
!!       X <- Po X Pv^T + Pv X Po^T,
!!
!!    where Po = 1/2 D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: X
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct
!
!        tmp = X Pv^T   =>   tmp^T = Pv X^T = - Pv X
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N', 'T', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_ao, &
                  Pv,      &
                  wf%n_ao, &
                  zero,    &
                  tmp,     &
                  wf%n_ao)
!
!     X = Po X Pv^T = Po tmp
!
      call dgemm('N', 'N', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one,     &
                  Po,      &
                  wf%n_ao, &
                  tmp,     &
                  wf%n_ao, &
                  zero,    &
                  X,       &
                  wf%n_ao)
!
!     X = X + Pv X Po^T = X - (-Pv X) Po^T = X - tmp^T Po^T
!
      call dgemm('T', 'T', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  -one,    &
                  tmp,     &
                  wf%n_ao, &
                  Po,      &
                  wf%n_ao, &
                  one,     &
                  X,       &
                  wf%n_ao)
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine project_redundant_rotations_hf
!
!
   subroutine construct_roothan_hall_hessian_hf(wf, H, Po, Pv)
!!
!!    Construct Roothan-Hall Hessian
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall Hessian,
!!
!!       H = Fvv - Foo = Pv^T F Pv - Po^T F Po,
!!
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: H
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fvv = Pv^T F Pv and set H = tmp = Fvv
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = wf%ao_fock
      call sandwich(tmp, Pv, Pv, wf%n_ao)
!
      H = tmp
!
!     Construct tmp = Foo = Po^T F Po and set H = H - tmp = Fvv - Foo
!
      tmp = wf%ao_fock
      call sandwich(tmp, Po, Po, wf%n_ao)
!
      H = H - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine construct_roothan_hall_hessian_hf
!
!
   subroutine get_packed_roothan_hall_gradient_hf(wf, G)
!!
!!    Get packed Roothan-Hall gradient
!!    Written by Eirik F. Kjønstad, Nov 2018
!!
!!    Constructs and returns the gradient as an
!!    anti-symmetrized packed vector.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao*(wf%n_ao - 1)/2, wf%n_densities), intent(inout) :: G
!
      real(dp), dimension(:,:), allocatable :: G_sq, Po, Pv
!
      call mem%alloc(Po, wf%n_ao, wf%n_ao)
      call mem%alloc(Pv, wf%n_ao, wf%n_ao)
!
      Po = zero
      Pv = zero
!
      call wf%construct_projection_matrices(Po, Pv, wf%ao_density)
!
      call mem%alloc(G_sq, wf%n_ao, wf%n_ao)
      G_sq = zero
!
      call wf%construct_roothan_hall_gradient(G_sq, Po, Pv, wf%ao_fock)
!
      call mem%dealloc(Po, wf%n_ao, wf%n_ao)
      call mem%dealloc(Pv, wf%n_ao, wf%n_ao)
!
      G = zero
      call packin_anti(G, G_sq, wf%n_ao)
!
      call mem%dealloc(G_sq, wf%n_ao, wf%n_ao)
!
   end subroutine get_packed_roothan_hall_gradient_hf
!
!
   subroutine construct_roothan_hall_gradient_hf(wf, G, Po, Pv, F)
!!
!!    Construct Roothan-Hall gradient
!!    Written by Eirik F. Kjønstad, 2018
!!
!!    Constructs the Roothan-Hall gradient,
!!
!!       G = Fov - Fvo = Po^T F Pv - Pv^T F Po,
!!
!!    where Po = D S and Pv = 1 - Po. In Po, D is the AO density and S
!!    is the AO overlap matrix.
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Po
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: Pv
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: F
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: G
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fov = Po^T F Pv and set G = tmp = Fov
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = F
      call sandwich(tmp, Po, Pv, wf%n_ao)
!
      G = tmp
!
!     Construct tmp = Fvo = Pv^T F Po and set H = H - tmp = Fov - Fvo
!
      tmp = F
      call sandwich(tmp, Pv, Po, wf%n_ao)
!
      G = G - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine construct_roothan_hall_gradient_hf
!
!
   subroutine do_roothan_hall_hf(wf, F, C, e)
!!
!!    Do Roothan-Hall
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Solves the equation F C = S C e for the orbital coefficients C.
!!    More precisely, it solves the equation in a linearly independent
!!    subspace,
!!
!!       P^T F P (P^T C) = P^T S P P^T C e = L L^T (P^T C) e,
!!
!!    which differs from the AO basis if linear dependence is present
!!    to within some given threshold. The components of the equation is
!!    given by the Cholesky decomposition
!!
!!       P^T S P = L L^T,
!!
!!    where P is referred to as the 'permutation matrix' and L the 'cholesky
!!    ao overlap' (note that these are member variables of the solver which
!!    must be set by a call to solver%decompose_ao_overlap(wf)). The number
!!    of linearly independent orbitals is wf%n_mo, whereas the full number
!!    is wf%n_ao.
!!
!!    Default is to not transform the Fock matrix to the MO basis. If
!!    do_mo_transformation is passed and is set true, the MO Fock matrix
!!    is initialized and transformed to the MO basis:
!!
!!       F_mo = C'^T P^T F P C'
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)     :: F
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(inout)  :: C
      real(dp), dimension(wf%n_mo), intent(inout)           :: e
!
      real(dp), dimension(:), allocatable   :: work
      real(dp), dimension(:,:), allocatable :: metric
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: FP
!
      real(dp), dimension(:,:), allocatable :: prev_C
!
      real(dp) :: ddot, orbital_dotprod
!
      integer :: info, p
!
      call mem%alloc(metric, wf%n_mo, wf%n_mo)
!
      call dgemm('N','T',                 &
                  wf%n_mo,                &
                  wf%n_mo,                &
                  wf%n_mo,                &
                  one,                    &
                  wf%cholesky_ao_overlap, &
                  wf%n_mo,                &
                  wf%cholesky_ao_overlap, &
                  wf%n_mo,                &
                  zero,                   &
                  metric,                 & ! metric = L L^T
                  wf%n_mo)
!
!     Allocate reduced space matrices
!
      call mem%alloc(ao_fock, wf%n_mo, wf%n_mo)
!
!     Construct reduced space Fock matrix, F' = P^T F P,
!     which is to be diagonalized over the metric L L^T
!
      call mem%alloc(FP, wf%n_ao, wf%n_mo)
!
      call dgemm('N','N',                       &
                  wf%n_ao,                      &
                  wf%n_mo,                      &
                  wf%n_ao,                      &
                  one,                          &
                  F,                            &
                  wf%n_ao,                      &
                  wf%pivot_matrix_ao_overlap,   &
                  wf%n_ao,                      &
                  zero,                         &
                  FP,                           &
                  wf%n_ao)
!
      call dgemm('T','N',                       &
                  wf%n_mo,                      &
                  wf%n_mo,                      &
                  wf%n_ao,                      &
                  one,                          &
                  wf%pivot_matrix_ao_overlap,   &
                  wf%n_ao,                      &
                  FP,                           &
                  wf%n_ao,                      &
                  zero,                         &
                  ao_fock,                      & ! F' = P^T F P
                  wf%n_mo)
!
      call mem%dealloc(FP, wf%n_mo, wf%n_ao)
!
!     Solve F'C' = L L^T C' e
!
      info = 0
!
      call mem%alloc(work, 4*wf%n_mo)
      work = zero
!
      call dsygv(1, 'V', 'L',       &
                  wf%n_mo,          &
                  ao_fock,          & ! ao_fock on entry, orbital coefficients on exit
                  wf%n_mo,          &
                  metric,           &
                  wf%n_mo,          &
                  e,                &
                  work,             &
                  4*(wf%n_mo),      &
                  info)
!
      call mem%dealloc(metric, wf%n_mo, wf%n_mo)
      call mem%dealloc(work, 4*wf%n_mo)
!
      if (info .ne. 0)  call output%error_msg('Error: could not solve Roothan-Hall equations.')
!
!     Transform back the solutions to original basis, C = P (P^T C) = P C'
!
      call mem%alloc(prev_C, wf%n_ao, wf%n_mo)
      prev_C = C
!
      call dgemm('N','N',                       &
                  wf%n_ao,                      &
                  wf%n_mo,                      &
                  wf%n_mo,                      &
                  one,                          &
                  wf%pivot_matrix_ao_overlap,   &
                  wf%n_ao,                      &
                  ao_fock,                      & ! orbital coefficients
                  wf%n_mo,                      &
                  zero,                         &
                  C,                            &
                  wf%n_ao)
!
      call mem%dealloc(ao_fock, wf%n_mo, wf%n_mo)
!
!     Test for orbitals that were approximately sign-flipped by dsygv,
!     resetting them if this is the case (this changes the orbitals
!     minimally, thus preserving the sweetness of whatever CC guess
!     is on file).
!
      do p = 1, wf%n_mo
!
         orbital_dotprod = ddot(wf%n_ao, prev_C(1, p), 1, C(1, p), 1)/ddot(wf%n_ao, C(1, p), 1, C(1, p), 1)
!
         if (orbital_dotprod .lt. zero) then
!
            C(:,p) = -C(:,p)
!
         endif
!
      enddo
!
   end subroutine do_roothan_hall_hf
!
!
   subroutine set_ao_density_to_sad_hf(wf)
!!
!!    Set AO density to SAD
!!    Written by Eirik F. Kjønstad, Aug-Sep 2018
!!
!!    The routine uses the pre-computed library of atomic densities to
!!    assemble the superposition of atomic densities (SAD) start guess.
!!    The library is assumed to be in the location given by the environ-
!!    ment variable "ET_SAD_DIR". Note that the library is located in the
!!    directory
!!
!!       eT/src/molecular_system/sad/ => set ET_SAD_DIR accordingly
!!
!!    The SAD guess is based on ground state UHF calculations, where
!!    the valence electrons are evenly spread out in the degenerate HOMO
!!    orbitals (to ensure a spherically symmetric AO density, which in
!!    turn ensures a rotationally invariant SAD guess for the molecule).
!!
!!    The algorithm is based on the procedure implemented in Psi4, although
!!    we perform spherical averaging of the densities in the UHF calculations
!!    in order to get a unique rotationally invariant SAD density. This is
!!    especially important for us, since we will use the initial idempotent
!!    density that results from a single Roothan-Hall step in multilevel HF
!!    calculations. The atomic calculations are performed with the ground
!!    state multiplicities, as listed in Griffiths, David J.,
!!    "Introduction to Quantum Mechanics." (1995).
!!
      implicit none
!
      class(hf) :: wf
!
      integer :: I, n_ao_on_atom, first_ao_on_atom, last_ao_on_atom, n_s_on_atom
!
      real(dp), dimension(:,:), allocatable :: atomic_density
!
      wf%ao_density = zero
!
      do I = 1, wf%system%n_atoms
!
         n_ao_on_atom     = wf%system%atoms(I)%n_ao
         n_s_on_atom      = wf%system%atoms(I)%n_shells
!
         first_ao_on_atom = wf%system%atoms(I)%shells(1)%first
         last_ao_on_atom  = wf%system%atoms(I)%shells(n_s_on_atom)%last
!
         call mem%alloc(atomic_density, n_ao_on_atom, n_ao_on_atom)
!
         call wf%system%atoms(I)%read_atomic_density(atomic_density)
!
         wf%ao_density(first_ao_on_atom:last_ao_on_atom, first_ao_on_atom:last_ao_on_atom) &
                                             = atomic_density(1:n_ao_on_atom, 1:n_ao_on_atom)
!
         call mem%dealloc(atomic_density, n_ao_on_atom, n_ao_on_atom)
!
      enddo
!
   end subroutine set_ao_density_to_sad_hf
!
!
   subroutine get_n_electrons_in_density_hf(wf, n_electrons)
!!
!!    Get number of electrons in density
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the number of electrons in the current
!!    AO density, a useful check on the sensibility of the
!!    initial AO density guess. It is by default printed
!!    by some solvers.
!!
      implicit none
!
      class(hf), intent(in)   :: wf
      real(dp), intent(inout) :: n_electrons
!
      integer :: ao
!
      real(dp), dimension(:,:), allocatable :: DS
!
      call mem%alloc(DS, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  wf%ao_overlap, &
                  wf%n_ao,       &
                  zero,          &
                  DS,            &
                  wf%n_ao)
!
      n_electrons = zero
!
      do ao = 1, wf%n_ao
!
         n_electrons = n_electrons + DS(ao, ao)
!
      enddo
!
      call mem%dealloc(DS, wf%n_ao, wf%n_ao)
!
   end subroutine get_n_electrons_in_density_hf
!
!
   subroutine set_ao_density_to_core_guess_hf(wf, h_wx)
!!
!!    Set AO density to core guess
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Solves the Roothan-Hall equation ignoring the two-
!!    electron AO density dependent part of the Fock matrix,
!!    giving an initial density from the resulting orbital
!!    coefficients.
!!
      implicit none
!
      class(hf) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      wf%ao_fock = h_wx
      call wf%do_roothan_hall(wf%ao_fock, wf%orbital_coefficients, wf%orbital_energies)
      call wf%construct_ao_density()
!
   end subroutine set_ao_density_to_core_guess_hf
!
!
   subroutine print_screening_settings_hf(wf)
!!
!!    Print screening settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      write(output%unit, '(/t6,a30,e11.4)') 'Coulomb screening threshold:  ', wf%coulomb_threshold
      write(output%unit, '(t6,a30,e11.4)') 'Exchange screening threshold: ', wf%exchange_threshold
      write(output%unit, '(t6,a30,e11.4)') 'ERI integral precision:       ', wf%libint_epsilon
      flush(output%unit)
!
   end subroutine print_screening_settings_hf
!
!
   subroutine set_n_mo_hf(wf)
!!
!!    Set number of molecular orbitals
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      write(output%unit, '(/t3,a)') '- Cholesky decomposition of AO overlap to get linearly independent orbitals:'
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap()
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      write(output%unit, '(/t6,a30,i8)') 'Number of occupied orbitals:  ', wf%n_o
      write(output%unit, '(t6,a30,i8)')  'Number of virtual orbitals:   ', wf%n_v
      write(output%unit, '(t6,a30,i8)')  'Number of molecular orbitals: ', wf%n_mo
      write(output%unit, '(t6,a30,i8)')  'Number of atomic orbitals:    ', wf%n_ao
!
     if (wf%n_mo .lt. wf%n_ao) &
            write(output%unit, '(/t6, a, i4, a)')'Removed ', wf%n_ao - wf%n_mo, ' AOs due to linear dep.'
!
   end subroutine set_n_mo_hf
!
!
   subroutine set_screening_and_precision_thresholds_hf(wf, gradient_threshold)
!!
!!    Set screening and precision thresholds
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Sets the screening thresholds for Coulomb and exchange
!!    integrals given the convergence threshold for the gradient
!!
!!       coulomb_threshold = gradient_threshold * 1.0d-2
!!       exchange_threshold = gradient_threshold * 1.0d-2
!!
!!       libint_epsilon = (gradient_threshold * 1.0d-3)**2
!!
!!    unless stricter thresholds are already set on input or by default.
!!
      implicit none
!
      class(hf), intent(inout) :: wf
!
      real(dp), intent(in) :: gradient_threshold
!
      if (wf%coulomb_threshold .gt. gradient_threshold*1.0d-3)  wf%coulomb_threshold  = gradient_threshold*1.0d-3
      if (wf%exchange_threshold .gt. gradient_threshold*1.0d-3) wf%exchange_threshold = gradient_threshold*1.0d-3
!
      wf%coulomb_threshold  = min(wf%coulomb_threshold,  1.0d-12)
      wf%exchange_threshold = min(wf%exchange_threshold, 1.0d-10)
!
      if (wf%libint_epsilon .gt. (wf%coulomb_threshold)**2) wf%libint_epsilon = (wf%coulomb_threshold)**2
!
   end subroutine set_screening_and_precision_thresholds_hf
!
!
   subroutine construct_molecular_gradient_hf(wf, E_qk)
!!
!!    Contruct molecular gradient
!!    Written by Åsmund H. Tveten and Eirik F. Kjønstad, 2019
!!
!!    Constructs the molecular gradient,
!! 
!!       E^x = Tr[D h^x] + (1/2)Tr[D G^x(D)] - Tr[D F D S^x] + h_nuc^x.
!!
!!    Here, x denotes the energy in the x direction. In the code, 
!!    x = (q,k), where q denotes the component (x,y, or z) and k 
!!    denotes the atom (k = 1,2,3,...,n_atoms).
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(3, wf%system%n_atoms), intent(inout) :: E_qk ! Molecular gradient
!
      real(dp), dimension(:,:,:,:), allocatable :: h_wxqk
      real(dp), dimension(:,:,:,:), allocatable :: G_wxqk
      real(dp), dimension(:,:,:,:), allocatable :: s_wxqk
!
      real(dp), dimension(:,:), allocatable :: DFD, FD
!
      real(dp), dimension(3, wf%system%n_atoms) :: TrDh_qk, TrDG_qk, TrDFDS_qk
!
      real(dp) :: ddot
!
      integer :: k, q
!
      type(timings) :: s_timer, h_timer, G_timer, non_integral_timer
!
!     Initialize timers 
!
      s_timer = timings('HF gradient - 1st derivative-integrals of S')
      h_timer = timings('HF gradient - 1st derivative-integrals of h')
      G_timer = timings('HF gradient - 1st derivative-integrals of G(D)')
      non_integral_timer = timings('HF gradient - non-integral-time')
!
!     Construct h_nuc^x, and the AO integral derivatives, h^x, S^x, and G^x(D)
!
      E_qk = wf%system%get_nuclear_repulsion_1der_numerical(1.0d-8) ! E_qk = h_nuc_qk
!
      call mem%alloc(h_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%alloc(G_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%alloc(s_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
!
      call s_timer%turn_on()
!
      s_wxqk = zero
      call wf%get_ao_s_wx_1der(s_wxqk)
!
      call s_timer%turn_off()
!
      call G_timer%turn_on()
!
      G_wxqk = zero
      call wf%construct_ao_G_1der(G_wxqk, wf%ao_density)
!
      do k = 1, wf%system%n_atoms

         do q = 1, 3

           call symmetric_sum(G_wxqk(:,:,q,k), wf%n_ao)
           G_wxqk(:,:,q,k) = half*G_wxqk(:,:,q,k)

         enddo

      enddo
!
      call G_timer%turn_off()
!
      call h_timer%turn_on()
!
      h_wxqk = zero
      call wf%get_ao_h_wx_1der(h_wxqk)
!
      call h_timer%turn_off()
!
      call non_integral_timer%turn_on()
!  
!     Construct D F D 
!
      call mem%alloc(FD, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_fock,    & 
                  wf%n_ao,       &
                  wf%ao_density, & 
                  wf%n_ao,       &
                  zero,          &
                  FD,            &  
                  wf%n_ao)
!
      call mem%alloc(DFD, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  one,           &
                  wf%ao_density, &
                  wf%n_ao,       &
                  FD,            &
                  wf%n_ao,       &
                  zero,          &
                  DFD,           & 
                  wf%n_ao)
!
      call mem%dealloc(FD, wf%n_ao, wf%n_ao)
!
!     Perform the traces, adding the contributions to the gradient 
!
      do k = 1, wf%system%n_atoms
         do q = 1, 3
!
            TrDh_qk(q,k) = ddot(wf%n_ao**2, wf%ao_density, 1, h_wxqk(:,:,q,k), 1)
            TrDG_qk(q,k) = ddot(wf%n_ao**2, wf%ao_density, 1, G_wxqk(:,:,q,k), 1)
            TrDFDS_qk(q,k) = ddot(wf%n_ao**2, DFD, 1, s_wxqk(:,:,q,k), 1)
!
            E_qk(q,k) = E_qk(q,k)   &
              + TrDh_qk(q,k)        & 
              + half*TrDG_qk(q,k)   &
              - half*TrDFDS_qk(q,k)
!
         enddo
      enddo 
!
      write(output%unit, '(/t6,a/)') 'Molecular gradient (Hartree/bohr):'
!
      do k = 1, wf%system%n_atoms
!
         write(output%unit, '(t6,a2,f12.6,f12.6,f12.6)') wf%system%atoms(k)%symbol, &
                                          E_qk(1,k), E_qk(2,k), E_qk(3,k) 
!
      enddo
!
      call mem%dealloc(DFD, wf%n_ao, wf%n_ao)
!
      call non_integral_timer%turn_off()
!
      call mem%dealloc(h_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%dealloc(G_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)
      call mem%dealloc(s_wxqk, wf%n_ao, wf%n_ao, 3, wf%system%n_atoms)   
!
   end subroutine construct_molecular_gradient_hf
!
!
end module hf_class
