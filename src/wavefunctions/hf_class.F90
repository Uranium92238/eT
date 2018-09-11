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
   use array_utilities
   use array_analysis
   use interval_class
   use index
!
   implicit none
!
!  Hartree-Fock wavefunction 
!
   type, extends(wavefunction):: hf
!
      real(dp), dimension(:,:), allocatable :: ao_density
      real(dp), dimension(:,:), allocatable :: ao_fock
      real(dp), dimension(:,:), allocatable :: mo_fock
      real(dp), dimension(:,:), allocatable :: orbital_energies
!
      real(dp), dimension(:,:), allocatable :: ao_overlap 
      real(dp), dimension(:,:), allocatable :: cholesky_ao_overlap
      real(dp), dimension(:,:), allocatable :: pivot_matrix_ao_overlap 
!
      real(dp) :: linear_dep_threshold = 1.0D-6
!
      real(dp) :: coulomb_threshold    = 1.0D-11 ! screening threshold 
      real(dp) :: exchange_threshold   = 1.0D-11 ! screening threshold 
      real(dp) :: integral_precision   = 1.0D-14 ! integral engine accuracy
!
	contains
!
!     Preparation and cleanup routines 
!
      procedure :: prepare                            => prepare_hf
      procedure :: cleanup                            => cleanup_hf
      procedure :: read_settings                      => read_settings_hf
      procedure :: read_hf_settings                   => read_hf_settings_hf
      procedure :: construct_ao_overlap               => construct_ao_overlap_hf
      procedure :: decompose_ao_overlap               => decompose_ao_overlap_hf
      procedure :: print_wavefunction_summary         => print_wavefunction_summary_hf
!
!     AO Fock and energy related routines 
!
      procedure :: construct_ao_fock                  => construct_ao_fock_hf            
      procedure :: construct_ao_fock_cumulative       => construct_ao_fock_cumulative_hf 
      procedure :: ao_fock_construction_loop          => ao_fock_construction_loop_hf
      procedure :: ao_fock_coulomb_construction_loop  => ao_fock_coulomb_construction_loop_hf
      procedure :: ao_fock_exchange_construction_loop => ao_fock_exchange_construction_loop_hf
      procedure :: construct_ao_fock_SAD              => construct_ao_fock_SAD_hf
      procedure :: construct_mo_fock                  => construct_mo_fock_hf
      procedure :: set_ao_fock                        => set_ao_fock_hf
      procedure :: get_fock_ov                        => get_fock_ov_hf
      procedure :: calculate_hf_energy_from_fock      => calculate_hf_energy_from_fock_hf
      procedure :: calculate_hf_energy_from_G         => calculate_hf_energy_from_G_hf
      procedure :: initialize_fock                    => initialize_fock_hf
      procedure :: destruct_fock                      => destruct_fock_hf
      procedure :: update_fock_and_energy             => update_fock_and_energy_hf
      procedure :: update_fock_and_energy_cumulative  => update_fock_and_energy_cumulative_hf
!
!     AO Density related routines 
!
      procedure :: construct_ao_density               => construct_ao_density_hf
      procedure :: rotate_ao_density                  => rotate_ao_density_hf
      procedure :: purify_ao_density                  => purify_ao_density_hf
      procedure :: decompose_ao_density               => decompose_ao_density_hf
      procedure :: get_ao_density                     => get_ao_density_hf
      procedure :: set_ao_density                     => set_ao_density_hf
      procedure :: initialize_density                 => initialize_density_hf
      procedure :: update_ao_density                  => update_ao_density_hf
      procedure :: save_ao_density                    => save_ao_density_hf
      procedure :: set_initial_ao_density_guess       => set_initial_ao_density_guess_hf
      procedure :: set_ao_density_to_sad              => set_ao_density_to_sad_hf 
      procedure :: set_ao_density_to_core_guess       => set_ao_density_to_core_guess_hf
      procedure :: get_n_electrons_in_density         => get_n_electrons_in_density_hf
      procedure :: construct_sp_density_schwarz       => construct_sp_density_schwarz_hf
!
!     MO orbital related routines 
!
      procedure :: do_roothan_hall                    => do_roothan_hall_hf
      procedure :: initialize_orbitals                => initialize_orbitals_hf
      procedure :: roothan_hall_update_orbitals       => roothan_hall_update_orbitals_hf
      procedure :: print_orbital_energies             => print_orbital_energies_hf
      procedure :: mo_transform                       => mo_transform_hf
      procedure :: mo_transform_and_save_h            => mo_transform_and_save_h_hf
!
!     Class variable initialize and destruct routines
!
      procedure :: initialize_ao_density              => initialize_ao_density_hf
      procedure :: initialize_ao_fock                 => initialize_ao_fock_hf
      procedure :: initialize_mo_fock                 => initialize_mo_fock_hf
      procedure :: initialize_ao_overlap              => initialize_ao_overlap_hf
      procedure :: initialize_orbital_energies        => initialize_orbital_energies_hf
      procedure :: initialize_pivot_matrix_ao_overlap => initialize_pivot_matrix_ao_overlap_hf
      procedure :: initialize_cholesky_ao_overlap     => initialize_cholesky_ao_overlap_hf
!
      procedure :: destruct_ao_density                => destruct_ao_density_hf
      procedure :: destruct_ao_fock                   => destruct_ao_fock_hf
      procedure :: destruct_mo_fock                   => destruct_mo_fock_hf
      procedure :: destruct_ao_overlap                => destruct_ao_overlap_hf
      procedure :: destruct_orbital_energies          => destruct_orbital_energies_hf
      procedure :: destruct_pivot_matrix_ao_overlap   => destruct_pivot_matrix_ao_overlap_hf
      procedure :: destruct_cholesky_ao_overlap       => destruct_cholesky_ao_overlap_hf
!
!     Gradient and Hessian related routines
!
      procedure :: construct_projection_matrices      => construct_projection_matrices_hf
      procedure :: project_redundant_rotations        => project_redundant_rotations_hf
!
      procedure :: construct_roothan_hall_hessian     => construct_roothan_hall_hessian_hf
      procedure :: construct_roothan_hall_gradient    => construct_roothan_hall_gradient_hf
!
!     Integral related routines 
!
      procedure :: construct_sp_eri_schwarz           => construct_sp_eri_schwarz_hf
      procedure :: get_n_sig_eri_sp                   => get_n_sig_eri_sp_hf
      procedure :: get_ao_h_wx                        => get_ao_h_wx_hf
      procedure :: get_ao_s_wx                        => get_ao_s_wx_hf
      procedure :: get_ao_mu_wx                       => get_ao_mu_wx_hf
!
   end type hf
!
!
contains
!
!
   subroutine prepare_hf(wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      wf%name = 'HF'
      write(output%unit, '(/t3,a)')  ':: Preparing ' // trim(wf%name) // ' wavefunction object'
!
      write(output%unit, '(/t3,a)') 'Reading wavefunction settings.'
!
      call wf%read_settings()
!
      write(output%unit, '(/t6,a30,e10.4)')  'Coulomb screening threshold:  ', wf%coulomb_threshold
      write(output%unit, '(t6,a30,e10.4)')   'Exchange screening threshold: ', wf%exchange_threshold
      write(output%unit, '(t6,a30,e10.4)')   'ERI integral precision:       ', wf%integral_precision
!
      write(output%unit, '(/t3,a)') 'Preparing molecular system object.'
      call wf%system%prepare()
!
      write(output%unit, '(/t6,a14,i1)') 'Charge:       ', wf%system%charge 
      write(output%unit, '(t6,a14,i1)')  'Multiplicity: ', wf%system%multiplicity 
!
      write(output%unit, '(/t6,a35,f18.12)') 'Nuclear repulsion energy (a.u.):   ', wf%system%get_nuclear_repulsion()
      write(output%unit, '(t6,a35,f18.12)')  'Bohr/angstrom value (CODATA 2010): ', bohr_to_angstrom
!
      call wf%system%print_geometry()
!
      wf%n_ao = wf%system%get_n_aos()
!
      call initialize_coulomb()
      call initialize_kinetic()
      call initialize_nuclear()
      call initialize_overlap()
!
      write(output%unit, '(/t3,a)') 'Cholesky decomposing AO overlap matrix to make sure used orbitals are'
      write(output%unit, '(t3,a)')  'linearly independent.'
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap() 
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      write(output%unit, '(/t6,a30,i4)') 'Number of occupied orbitals:  ', wf%n_o 
      write(output%unit, '(t6,a30,i4)')  'Number of virtual orbitals:   ', wf%n_v
      write(output%unit, '(t6,a30,i4)')  'Number of molecular orbitals: ', wf%n_mo 
      write(output%unit, '(t6,a30,i4)')  'Number of atomic orbitals:    ', wf%n_ao 
!
   end subroutine prepare_hf
!
!
   subroutine print_wavefunction_summary_hf(wf)
!!
!!    Print wavefunction summary 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Prints information related to the wavefunction,
!!    most of which is meaningful only for a properly 
!!    converged wavefunction. Should be overwritten in 
!!    descendants if more or less or other information 
!!    is present. 
!!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      real(dp) :: homo_lumo_gap
!
      integer(i15) :: mo 
!
      write(output%unit, '(/t3,a,a,a)') ':: Summary of ', trim(wf%name), ' wavefunction energetics (a.u.)'
!
      homo_lumo_gap = wf%orbital_energies(wf%n_o + 1, 1) - wf%orbital_energies(wf%n_o, 1)
!
      write(output%unit, '(/t3,a26,f19.12)') 'HOMO-LUMO gap:            ', homo_lumo_gap
      write(output%unit, '(t3,a26,f19.12)')  'Nuclear repulsion energy: ', wf%system%get_nuclear_repulsion()
      write(output%unit, '(t3,a26,f19.12)')  'Electronic energy:        ', wf%energy - wf%system%get_nuclear_repulsion()
      write(output%unit, '(t3,a26,f19.12)')  'Total energy:             ', wf%energy
!
      call wf%print_orbital_energies('3')
!
   end subroutine print_wavefunction_summary_hf
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
      integer(i15) :: n_records, i
!
      character(len=100) :: line, value 
!
      if (requested_section('hf')) then ! User has requested something 
!
         call move_to_section('hf', n_records)
!
         do i = 1, n_records
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:18) == 'coulomb threshold:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) wf%coulomb_threshold
               cycle
!
            elseif (line(1:19) == 'exchange threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) wf%exchange_threshold
               cycle
!
            elseif (line(1:19) == 'integral precision:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) wf%integral_precision
               cycle
!
            endif
!
         enddo
!
      endif 
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
      write(output%unit, '(/t'//trim(indent)//',a)') 'Molecular orbital energies:'
!
      call print_vector(wf%orbital_energies, wf%n_ao, indent)
!
   end subroutine print_orbital_energies_hf
!
!
   subroutine mo_transform_and_save_h_hf(wf)
!!
!!    MO transform and save h 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none 
!
      class(hf) :: wf 
!
      real(dp), dimension(:,:), allocatable :: h_wx, h_pq 
!
      type(file) :: h_pq_file
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call mem%alloc(h_pq, wf%n_mo, wf%n_mo)
!
      call wf%get_ao_h_wx(h_wx)
      call wf%mo_transform(h_wx, h_pq)
!
      call h_pq_file%init('h_pq', 'sequential', 'unformatted')
      call disk%open_file(h_pq_file, 'readwrite', 'rewind')
!
      write(h_pq_file%unit) h_pq 
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
      call mem%dealloc(h_pq, wf%n_mo, wf%n_mo)     
!
   end subroutine mo_transform_and_save_h_hf
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
   subroutine update_fock_and_energy_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)
!!
!!    Update Fock and energy
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted 
!!    wavefunctions.
!!
      implicit none 
!
      class(hf) :: wf 
!
      integer(i15), intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      call wf%construct_ao_fock(sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx)      
!
      call wf%calculate_hf_energy_from_fock(wf%ao_fock, h_wx)
!
   end subroutine update_fock_and_energy_hf
!
!
   subroutine update_fock_and_energy_cumulative_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, prev_ao_density, h_wx)
!!
!!    Update Fock and energy cumulatively
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    This routine guides the construction of the Fock matrix (or matrices for
!!    unrestricted wavefunctions) from the current AO density (or densities).
!!    It is called by the solver and is overwritten for unrestricted 
!!    wavefunctions.
!!
      implicit none 
!
      class(hf) :: wf 
!
      integer(i15), intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: prev_ao_density
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      call wf%construct_ao_fock_cumulative(sp_eri_schwarz, sp_eri_schwarz_list, &
                                           n_s, prev_ao_density, h_wx) 
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
      write(ao_density%unit, *) wf%ao_density
      call disk%close_file(ao_density)
!
   end subroutine save_ao_density_hf
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
      call wf%destruct_orbital_energies()
      call wf%destruct_ao_overlap()
      call wf%destruct_orbital_coefficients()
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
   subroutine initialize_orbital_energies_hf(wf)
!!
!!    Initialize orbital energies
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (.not. allocated(wf%orbital_energies)) call mem%alloc(wf%orbital_energies, wf%n_mo, 1)
!
   end subroutine initialize_orbital_energies_hf
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
   subroutine destruct_orbital_energies_hf(wf)
!!
!!    Destruct AO overlap
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf) :: wf
!
      if (allocated(wf%orbital_energies)) call mem%dealloc(wf%orbital_energies, wf%n_mo, 1)
!
   end subroutine destruct_orbital_energies_hf
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
!   subroutine destruct_mo_coefficients_hf(wf)
!!!
!!!    Destruct MO coefficients
!!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!!
!      implicit none
!!
!      class(hf) :: wf
!!
!      if (allocated(wf%orbital_coefficients)) call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!!
!   end subroutine destruct_mo_coefficients_hf
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
   subroutine construct_sp_eri_schwarz_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s)
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
      integer(i15) :: n_s
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2) :: sp_eri_schwarz
!
      integer(i15), dimension(n_s*(n_s + 1)/2, 3) :: sp_eri_schwarz_list
      integer(i15), dimension(:, :), allocatable  :: sp_eri_schwarz_list_copy
!
      integer(i15), dimension(:,:), allocatable :: sp_eri_schwarz_index_list
      real(dp), dimension(:,:),     allocatable :: sorted_sp_eri_schwarz
!
!     Local variables
!
      integer(i15) :: s1, s2, s1s2, counter
!
      real(dp) :: maximum
!
      real(dp), dimension(:,:), allocatable :: g
!
      type(interval) :: A_interval, B_interval
!
!     Set the maximum element in each shell pair 
!
!$omp parallel do private(s1, s2, s1s2, A_interval, B_interval, g, maximum) schedule(dynamic)
      do s1 = 1, n_s
         do s2 = 1, s1
!
            s1s2 = (max(s1,s2)*(max(s1,s2)-3)/2) + s1 + s2
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
            call wf%system%ao_integrals%construct_ao_g_wxyz(g, s1, s2, s1, s2)
!
            maximum = get_abs_max(g, ((A_interval%size)*(B_interval%size))**2)
!
            call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                (A_interval%size)*(B_interval%size))
!
            sp_eri_schwarz(s1s2, 1) = sqrt(maximum)
!
            sp_eri_schwarz_list(s1s2, 1) = s1 
            sp_eri_schwarz_list(s1s2, 2) = s2
!
         enddo
      enddo
!$omp end parallel do
!
!     Sort the sp_eri_schwarz vector and use the resulting index list 
!     to resort the sp_eri_schwarz_list matrix 
!
      call mem%alloc_int(sp_eri_schwarz_index_list, n_s*(n_s + 1)/2, 1)
      sp_eri_schwarz_index_list = 0
!
      call mem%alloc(sorted_sp_eri_schwarz, n_s*(n_s + 1)/2, 1)
      sorted_sp_eri_schwarz = sp_eri_schwarz
!
      call get_n_highest(n_s*(n_s + 1)/2, n_s*(n_s + 1)/2, sp_eri_schwarz, sorted_sp_eri_schwarz, sp_eri_schwarz_index_list)
!
      sp_eri_schwarz(:,2) = sp_eri_schwarz(:,1)
      sp_eri_schwarz(:,1) = sorted_sp_eri_schwarz(:,1)
      call mem%dealloc(sorted_sp_eri_schwarz, n_s*(n_s + 1)/2, 1)
!
      sp_eri_schwarz_list(:,3) = sp_eri_schwarz_index_list(:,1)
      call mem%dealloc_int(sp_eri_schwarz_index_list, n_s*(n_s + 1)/2, 1)
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
      real(dp), dimension(:,:), allocatable :: D_red 
!
      type(interval) :: A_interval, B_interval
!
      integer(i15) :: s1, s2 
!
      real(dp) :: maximum
!
!$omp parallel do private(s1, s2, A_interval, B_interval, D_red, maximum) schedule(dynamic)
      do s1 = 1, wf%system%n_s
         do s2 = 1, s1
!
            A_interval = wf%system%shell_limits(s1)
            B_interval = wf%system%shell_limits(s2)
!
            call mem%alloc(D_red, A_interval%size, B_interval%size)
!
            D_red = D(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
            maximum = get_abs_max(D_red, (A_interval%size)*(B_interval%size))
!
            call mem%dealloc(D_red, A_interval%size, B_interval%size)
!
            sp_density_schwarz(s1, s2) = maximum
            sp_density_schwarz(s2, s1) = maximum
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_sp_density_schwarz_hf
!
!
   subroutine get_n_sig_eri_sp_hf(wf, n_sig_sp, sp_eri_schwarz, threshold)
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
      integer(i15), intent(inout) :: n_sig_sp 
!
      real(dp), intent(in) :: threshold
!
      real(dp), dimension(wf%system%n_s*(wf%system%n_s + 1)/2, 2), intent(in) :: sp_eri_schwarz
!
      integer(i15) :: s1s2
!
      n_sig_sp = 0
      do s1s2 = 1, wf%system%n_s*(wf%system%n_s + 1)/2
!
         if (sp_eri_schwarz(s1s2, 1)**2 .lt. threshold) then
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
   subroutine construct_ao_fock_SAD_hf(wf, coulomb, exchange, precision)
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
      integer(i15) :: n_s
!
      real(dp), optional :: coulomb, exchange, precision ! Non-standard thresholds
!
      real(dp) :: coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(:,:), allocatable :: sp_eri_schwarz, sp_density_schwarz
!
      real(dp), dimension(:,:), allocatable :: h_wx, h_AB
      integer(i15) :: w, x, y, z, wx, yz, wz, yx
!
      integer(i15) :: A, B, C, D, skip, thread, omp_get_thread_num, atom
!
      type(interval) :: A_interval
      type(interval) :: B_interval
      type(interval) :: C_interval
      type(interval) :: D_interval
!
      real(dp) :: maximum, max_eri, max_density
!
      real(dp), dimension(:,:), allocatable :: g, g_C, g_K, D_yz
!
      integer(i15), dimension(:,:), allocatable :: shells_on_atoms
!!
!     Set thresholds to ignore Coulomb and exchange terms,
!     as well as the desired Libint integral precision  
!
      coulomb_thr = 1.0D-11 
      if (present(coulomb)) coulomb_thr = coulomb 
!
      exchange_thr = 1.0D-11
      if (present(exchange)) exchange_thr = exchange 
!
      precision_thr = 1.0D-14
      if (present(precision)) precision_thr = precision 
!
      call set_coulomb_precision(1.0d-25)
!
      n_s = wf%system%get_n_shells()
!
   call mem%alloc(sp_eri_schwarz, n_s, n_s)
!
!$omp parallel do private(A, B, A_interval, B_interval, g, maximum) schedule(dynamic)
      do A = 1, n_s
         do B = 1, A
!
            A_interval = wf%system%shell_limits(A)
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(g, (A_interval%size)*(B_interval%size), &
                              (A_interval%size)*(B_interval%size))
!
            call wf%system%ao_integrals%construct_ao_g_wxyz(g, A, B, A, B)
!
            maximum = get_abs_max(g, ((A_interval%size)*(B_interval%size))**2)
!
            call mem%dealloc(g, (A_interval%size)*(B_interval%size), &
                                (A_interval%size)*(B_interval%size))
!
            sp_eri_schwarz(A, B) = sqrt(maximum)
            sp_eri_schwarz(B, A) = sqrt(maximum)
!
         enddo
      enddo
!$omp end parallel do
!
      call mem%alloc_int(shells_on_atoms, wf%system%n_atoms, 2) ! [first, last]
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
      call mem%alloc(sp_density_schwarz, n_s, n_s)
      sp_density_schwarz = zero
!
!$omp parallel do private(A, B, A_interval, B_interval, D_yz, maximum) schedule(dynamic)
      do atom = 1, wf%system%n_atoms
         do A = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
            A_interval = wf%system%shell_limits(A)
!
            do B = shells_on_atoms(atom, 1), shells_on_atoms(atom, 2)
!
               B_interval = wf%system%shell_limits(B)
!
               call mem%alloc(D_yz, (A_interval%size), (B_interval%size))
!
               D_yz = wf%ao_density(A_interval%first : A_interval%last, B_interval%first : B_interval%last)
!
               maximum = get_abs_max(D_yz, (A_interval%size)*(B_interval%size))
!
               call mem%dealloc(D_yz, (A_interval%size), (B_interval%size))
!
               sp_density_schwarz(A, B) = maximum
!
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call set_coulomb_precision(precision_thr)
!
      wf%ao_fock = zero
!
      max_density = get_abs_max(sp_density_schwarz, n_s**2)
      max_eri     = get_abs_max(sp_eri_schwarz, n_s**2)
!
!      write(output%unit, *)'doing coulomb_'
!      flush(output%unit)
!
!$omp parallel do &
!$omp private(A, B, C, D, A_interval, B_interval, C_interval, D_interval, w, x, y, z, wx, yz, &
!$omp g_C) schedule(dynamic)
      do A = 1, n_s
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
                     if (sp_eri_schwarz(A, B)*sp_eri_schwarz(C, D)*sp_density_schwarz(C, D) .lt. coulomb_thr) cycle               
!
                     call mem%alloc(g_C, (A_interval%size)*(B_interval%size), &
                                    (C_interval%size)*(D_interval%size))
!
                     call wf%system%ao_integrals%construct_ao_g_wxyz(g_C, A, B, C, D)
!     
!                    Add Fock matrix contributions
!     
                     if (A .ne. B) then
!
                        do w = A_interval%first, A_interval%last
                           do x = B_interval%first, B_interval%last
!
                              wx = A_interval%size*(x - B_interval%first) + w - A_interval%first + 1
!
                              do y = C_interval%first, C_interval%last
                                 do z = D_interval%first, D_interval%last
!
                                    yz = C_interval%size*(z - D_interval%first) + y  - C_interval%first + 1
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + g_C(wx, yz)*wf%ao_density(y, z)
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
!
                              wx= A_interval%size*(x - A_interval%first) + w - A_interval%first + 1
!
                              do y = C_interval%first, C_interval%last
                                 do z = D_interval%first, D_interval%last
!
                                    yz = C_interval%size*(z - D_interval%first) + y  - C_interval%first + 1
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + g_C(wx, yz)*wf%ao_density(y, z)
!
                                 enddo
                              enddo
                           enddo
                        enddo
                     endif
!                  
                     call mem%dealloc(g_C, (A_interval%size)*(B_interval%size), &
                                       (C_interval%size)*(D_interval%size))
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!$omp parallel do &
!$omp private(A, B, C, D, A_interval, B_interval, C_interval, D_interval, w, x, y, z, wz, yx, &
!$omp g_K) schedule(dynamic)
      do A = 1, n_s
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
                     call mem%alloc(g_K, (A_interval%size)*(D_interval%size), &
                                    (C_interval%size)*(B_interval%size))
!
                     call wf%system%ao_integrals%construct_ao_g_wxyz(g_K, A, D, C, B)
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
                                    yx = C_interval%size*(x - B_interval%first) + y  - C_interval%first + 1
                                    wz= A_interval%size*(z - D_interval%first) + w - A_interval%first + 1
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + (- half*g_K(wz, yx))*wf%ao_density(y, z)
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
                                    yx = C_interval%size*(x - B_interval%first) + y  - C_interval%first + 1
                                    wz= A_interval%size*(z - D_interval%first) + w - A_interval%first + 1
!
                                    wf%ao_fock(w, x) = wf%ao_fock(w, x) + (- half*g_K(wz, yx))*wf%ao_density(y, z)
!
                                 enddo
                              enddo
                           enddo
                        enddo
                     endif
!                  
                     call mem%dealloc(g_K, (A_interval%size)*(B_interval%size), &
                                       (C_interval%size)*(D_interval%size))
!
                  enddo
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
!
!   write(output%unit, *)'doing symmetrize'
!   flush(output%unit)
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
!      write(output%unit, *)'h and energy'
!      flush(output%unit)
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
!
!$omp parallel do &
!$omp private(A, B, h_AB, A_interval, B_interval, x, y) schedule(static)
      do A = 1, n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, n_s
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(h_AB, A_interval%size, B_interval%size)
            call wf%system%ao_integrals%construct_ao_h_wx(h_AB, A, B)
!!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!!
                   h_wx(A_interval%first - 1 + x, B_interval%first - 1 + y) = h_AB(x, y)
!!
                enddo
             enddo
!
            call mem%dealloc(h_AB, A_interval%size, B_interval%size)
!
         enddo
      enddo
!$omp end parallel do
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
!      write(output%unit, *)'done'
!      flush(output%unit)
!
   end subroutine construct_ao_fock_SAD_hf
!
!
   subroutine construct_ao_fock_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx, coulomb, exchange, precision)
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
      integer(i15), intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      real(dp), optional :: coulomb, exchange, precision   ! Non-standard thresholds, optionals
      real(dp) :: coulomb_thr, exchange_thr, precision_thr ! Actual thresholds 
!
      integer(i15) :: thread, n_threads, omp_get_max_threads
!
      real(dp), dimension(:,:), allocatable :: F, sp_density_schwarz
!
      integer(i15) :: s1s2
!
      integer(i15) :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
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
      precision_thr = wf%integral_precision
      if (present(precision)) precision_thr = precision 
!
!     Construct the density screening vector and the maximum element in the density
!
      call mem%alloc(sp_density_schwarz, n_s, n_s)
      call wf%construct_sp_density_schwarz(sp_density_schwarz, wf%ao_density)
      max_D_schwarz = get_abs_max(sp_density_schwarz, n_s**2)
!
!     Compute number of significant ERI shell pairs (the Fock construction 
!     only loops over these shell pairs) and the maximum element 
!
      call wf%get_n_sig_eri_sp(n_sig_sp, sp_eri_schwarz, 1.0d-20)
      max_eri_schwarz = get_abs_max(sp_eri_schwarz, n_s*(n_s + 1)/2)
!
!     Construct the two electron part of the Fock matrix, using the screening vectors 
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      F = zero 
!
      call wf%ao_fock_construction_loop(F, wf%ao_density, n_threads, max_D_schwarz, max_eri_schwarz,  & 
                                         sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list,     &
                                         n_s, n_sig_sp, coulomb_thr, exchange_thr, precision_thr,     &
                                         wf%system%shell_limits)
!
      call mem%dealloc(sp_density_schwarz, n_s, n_s)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix,
!     and symmetrize the result 
!
      wf%ao_fock = zero
      do thread = 1, n_threads
!
         call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock, 1)
!
      enddo
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads) 
!
      call symmetric_sum(wf%ao_fock, wf%n_ao)
      wf%ao_fock = wf%ao_fock*half
!
      wf%ao_fock = wf%ao_fock + h_wx
!
   end subroutine construct_ao_fock_hf
!
!
  subroutine construct_ao_fock_cumulative_hf(wf, sp_eri_schwarz, sp_eri_schwarz_list, &
                                             n_s, prev_ao_density, h_wx, coulomb, exchange, precision)
!!
!!    Construct AO Fock matrix cumulatively
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Updates the AO Fock matrix incrementally,
!!
!!       F_αβ = F_αβ + sum_γδ g_αβγδ dD_γδ - 1/2 * sum_γδ g_αδγβ dD_γδ,
!!
!!    where dD is the difference in the AO density. This routine is integral direct, and
!!    it calculates the Hartree-Fock energy. It screens on the density difference and is 
!!    therefore appropriate when small changes are made to the density matrix (such as 
!!    in Hartree-Fock iterations). It is assumed that wf%ao_density is the current density
!!    and that the previous density is passed in the array prev_ao_density.
!!
      implicit none
!
      class(hf) :: wf
!
      integer(i15) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: prev_ao_density
!
      real(dp), dimension(n_s*(n_s + 1)/2, 1), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), optional :: coulomb, exchange, precision   ! Non-standard thresholds, optionals
      real(dp) :: coulomb_thr, exchange_thr, precision_thr ! Actual thresholds 
!
      integer(i15) :: n_threads, omp_get_max_threads, thread
!
      real(dp), dimension(:,:), allocatable :: F 
!
      integer(i15) :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: sp_density_schwarz
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
      precision_thr = wf%integral_precision
      if (present(precision)) precision_thr = precision 
!
!     Overwrite the AO density with the density difference,
!     and compute the density screening vector 
!
      call daxpy(wf%n_ao**2, -one, prev_ao_density, 1, wf%ao_density, 1)
!
      call mem%alloc(sp_density_schwarz, n_s, n_s)
      call wf%construct_sp_density_schwarz(sp_density_schwarz, wf%ao_density)     
      max_D_schwarz = get_abs_max(sp_density_schwarz, n_s**2)
!
!     Compute number of significant ERI shell pairs (the Fock construction 
!     only loops over these shell pairs) and the maximum element 
!
      call wf%get_n_sig_eri_sp(n_sig_sp, sp_eri_schwarz, 1.0d-20)
      max_eri_schwarz = get_abs_max(sp_eri_schwarz, n_s*(n_s + 1)/2)
!
!     Construct the two electron part of the Fock matrix, using the screening vectors 
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      call dscal(n_threads*(wf%n_ao)**2, zero, F, 1) 
!
      call wf%ao_fock_construction_loop(F, wf%ao_density, n_threads, max_D_schwarz, max_eri_schwarz, & 
                                         sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list,    &
                                         n_s, n_sig_sp, coulomb_thr, exchange_thr, precision_thr,    &
                                         wf%system%shell_limits)
!
      call daxpy(wf%n_ao**2, one, prev_ao_density, 1, wf%ao_density, 1)
      call mem%dealloc(sp_density_schwarz, n_s, n_s)
!
!     Put the accumulated Fock matrices from each thread into the Fock matrix,
!     and symmetrize the result
!
      do thread = 1, n_threads
!
         call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock, 1)
!
      enddo
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads) 
!
      call symmetric_sum(wf%ao_fock, wf%n_ao)
      call dscal((wf%n_ao)**2, half, wf%ao_fock, 1) 
!
   end subroutine construct_ao_fock_cumulative_hf
!
!
   subroutine ao_fock_construction_loop_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,    & 
                                          sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list, &
                                          n_s, n_sig_sp, coulomb_thr, exchange_thr, precision_thr, shells)
!!
!!    AO Fock construction loop 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    This routine constructs the entire two-electron part of the Fock matrix, 
!!
!!       F_αβ = sum_γδ g_αβγδ D_γδ - 1/2 * sum_γδ g_αδγβ D_γδ,
!!
!!    where contributions from different threads are gathered column blocks 
!!    of the incoming F matrix. 
!!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      integer(i15), intent(in) :: n_threads, n_s, n_sig_sp
!
      type(interval), dimension(n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, exchange_thr, precision_thr
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
      real(dp), dimension(n_s, n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer(i15) :: w, x, y, omp_get_thread_num, z, wx, yz, s1s2, s1, s2, s3, s4, s4_max, tot_dim 
      integer(i15) :: s3s4, s3s4_sorted, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(:,:), allocatable :: g 
!
      integer(i15) :: max_shell_size, thread, skip
!
!     Preallocate the vector that holds the shell quadruple 
!     ERI integrals, then enter the construction loop 
!
      call wf%system%get_max_shell_size(max_shell_size)
      call mem%alloc(g, max_shell_size**4, 1)
!
!$omp parallel do                                                                             &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34,     &
!$omp w, x, y, z, wx, yz, temp1, temp2, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset, &
!$omp temp4, temp5, temp6, temp7, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,        &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, s3s4_sorted, deg_12_34,                   &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
         thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix 
!
         sp_eri_schwarz_s1s2 = sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = sp_eri_schwarz_list(s1s2_packed, 2)
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
               temp = sp_eri_schwarz_s1s2*sp_eri_schwarz(s3s4, 2)

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
               call wf%system%ao_integrals%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4,         &
                  precision_thr/max(temp7,temp8), thread, skip, shells(s1)%size, shells(s2)%size, &
                  shells(s3)%size, shells(s4)%size)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
               g(1:tot_dim,1) = deg*g(1:tot_dim,1)
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
                           temp = g(wxyz, 1)
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
      call mem%dealloc(g, max_shell_size**4, 1)
!
   end subroutine ao_fock_construction_loop_hf
!
!
   subroutine ao_fock_coulomb_construction_loop_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz,     & 
                                                   sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list, &
                                                   n_s, n_sig_sp, coulomb_thr, precision_thr, shells)
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
      integer(i15), intent(in) :: n_threads, n_s, n_sig_sp
!
      type(interval), dimension(n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, coulomb_thr, precision_thr
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
      real(dp), dimension(n_s, n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp7, deg, deg_12, deg_34, deg_12_34
!
      integer(i15) :: w, x, y, omp_get_thread_num, z, wx, yz, s1s2, s1, s2, s3, s4, s4_max, tot_dim 
      integer(i15) :: s3s4, s3s4_sorted, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(:,:), allocatable :: g 
!
      integer(i15) :: max_shell_size, thread, skip
!
!     Preallocate the vector that holds the shell quadruple 
!     ERI integrals, then enter the construction loop 
!
      call wf%system%get_max_shell_size(max_shell_size)
      call mem%alloc(g, max_shell_size**4, 1)
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, wx, yz, temp1, temp2, d1, d2, thread, thread_offset,                    &
!$omp temp7, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,                                &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, s3s4_sorted, deg_12_34,               &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
         thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix 
!
         sp_eri_schwarz_s1s2 = sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = sp_eri_schwarz_list(s1s2_packed, 2)
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
               temp = sp_eri_schwarz_s1s2*sp_eri_schwarz(s3s4, 2)

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
               call wf%system%ao_integrals%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp7, thread, skip, shells(s1)%size, shells(s2)%size,    &
                  shells(s3)%size, shells(s4)%size)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
               g(1:tot_dim,1) = deg*g(1:tot_dim,1)
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
                           temp = g(wxyz, 1)
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
      call mem%dealloc(g, max_shell_size**4, 1)
!
   end subroutine ao_fock_coulomb_construction_loop_hf
!
!
   subroutine ao_fock_exchange_construction_loop_hf(wf, F, D, n_threads, max_D_schwarz, max_eri_schwarz, & 
                                          sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list,       &
                                          n_s, n_sig_sp, exchange_thr, precision_thr, shells)
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
      integer(i15), intent(in) :: n_threads, n_s, n_sig_sp
!
      type(interval), dimension(n_s), intent(in) :: shells
!
      real(dp), dimension(wf%n_ao, wf%n_ao*n_threads)   :: F 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D
!
      real(dp), intent(in) :: max_D_schwarz, max_eri_schwarz, exchange_thr, precision_thr
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
      real(dp), dimension(n_s, n_s), intent(in)               :: sp_density_schwarz
!
      real(dp) :: d1, d2, d3, d4, d5, d6, sp_eri_schwarz_s1s2
      real(dp) :: temp, temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, deg, deg_12, deg_34, deg_12_34
!
      integer(i15) :: w, x, y, omp_get_thread_num, z, wx, yz, s1s2, s1, s2, s3, s4, s4_max, tot_dim 
      integer(i15) :: s3s4, s3s4_sorted, w_red, x_red, y_red, z_red, thread_offset, wxyz, s1s2_packed
!
      real(dp) :: sp_density_schwarz_s1s2, sp_density_schwarz_s3s2, sp_density_schwarz_s3s1
!
      real(dp), dimension(:,:), allocatable :: g 
!
      integer(i15) :: max_shell_size, thread, skip
!
!     Preallocate the vector that holds the shell quadruple 
!     ERI integrals, then enter the construction loop 
!
      call wf%system%get_max_shell_size(max_shell_size)
      call mem%alloc(g, max_shell_size**4, 1)
!
!$omp parallel do                                                                         &
!$omp private(s1, s2, s3, s4, deg, s4_max, temp, s1s2, s1s2_packed, s3s4, deg_12, deg_34, &
!$omp w, x, y, z, wx, yz, temp3, d1, d2, d3, d4, d5, d6, thread, thread_offset,           &
!$omp temp4, temp5, temp6, temp8, w_red, x_red, tot_dim, y_red, z_red, wxyz, g,           &
!$omp sp_eri_schwarz_s1s2, sp_density_schwarz_s1s2, s3s4_sorted, deg_12_34,               &
!$omp sp_density_schwarz_s3s2, sp_density_schwarz_s3s1, skip) schedule(dynamic)
      do s1s2 = 1, n_sig_sp
!
         thread = omp_get_thread_num()
         thread_offset = thread*wf%n_ao ! Start column of thread's Fock matrix 
!
         sp_eri_schwarz_s1s2 = sp_eri_schwarz(s1s2, 1)
!
         s1s2_packed = sp_eri_schwarz_list(s1s2, 3) ! The s1s2-th largest packed index
!
         s1 = sp_eri_schwarz_list(s1s2_packed, 1)
         s2 = sp_eri_schwarz_list(s1s2_packed, 2)
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
               temp = sp_eri_schwarz_s1s2*sp_eri_schwarz(s3s4, 2)

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
               call wf%system%ao_integrals%construct_ao_g_wxyz_epsilon(g, s1, s2, s3, s4, &
                  precision_thr/temp8, thread, skip, shells(s1)%size, shells(s2)%size,    &
                  shells(s3)%size, shells(s4)%size)
!
               if (skip == 1) cycle
!
               tot_dim = (shells(s1)%size)*(shells(s2)%size)*(shells(s3)%size)*(shells(s4)%size)
!
               g(1:tot_dim,1) = deg*g(1:tot_dim,1)
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
                           temp = g(wxyz, 1)
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
      call mem%dealloc(g, max_shell_size**4, 1)
!
   end subroutine ao_fock_exchange_construction_loop_hf
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
      wf%energy = wf%energy + two*(one/four)*ddot((wf%n_ao)**2, wf%ao_density, 1, half_GD_wx, 1)

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
!!       Tr(AB) = sum_x (AB)_xx = sum_xy A_xy B_yx = sum_xy A_xy B_xy.
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
      wf%energy = wf%energy + one/two*ddot((wf%n_ao)**2, h_wx, 1, wf%ao_density, 1)
      wf%energy = wf%energy + one/two*ddot((wf%n_ao)**2, wf%ao_density, 1, F_wx, 1)

   end subroutine calculate_hf_energy_from_fock_hf
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
      real(dp), dimension(:,:) :: F ! Packed
!
      call squareup(F, wf%ao_fock, wf%n_ao)
!
   end subroutine set_ao_fock_hf
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
   subroutine construct_mo_fock_hf(wf, F_pq)
!!
!!    Construct MO Fock
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!!    Constructs the MO Fock matrix F_pq using the current AO 
!!    Fock and the orbital coefficients. 
!!
      implicit none
!
      class(hf), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: F_pq
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
                  F_pq,                      & ! F = C^T F^ao C
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
      integer(kind=4), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(:,:), allocatable :: perm_matrix
!
      real(dp), dimension(:,:), allocatable :: csc
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: rank
      integer(i15) :: i, j
!
      allocate(used_diag(wf%n_ao, 1))
!
      wf%ao_density = half*wf%ao_density
      call full_cholesky_decomposition_system(wf%ao_density, wf%orbital_coefficients, wf%n_ao, rank,&
                                                      1.0D-12, used_diag)
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
         perm_matrix(used_diag(j,1), j) = one
!
      enddo
!
      deallocate(used_diag)
!
!     Sanity check
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      call dgemm('N','N', &
                  wf%n_ao, &
                  wf%n_ao, &
                  wf%n_ao, &
                  one, &
                  perm_matrix, &
                  wf%n_ao, &
                  wf%orbital_coefficients, &
                  wf%n_ao, &
                  zero, &
                  tmp, &
                  wf%n_ao)
!
      wf%orbital_coefficients = tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
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
      integer(kind=4), dimension(:, :), allocatable :: used_diag
!
      real(dp), dimension(:, :), allocatable :: L
!
      integer(i15) :: i, j
!
      allocate(used_diag(wf%n_ao, 1))
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
         wf%pivot_matrix_ao_overlap(used_diag(j, 1), j) = one
!
      enddo
!
      deallocate(used_diag)
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
      integer(i15) :: iteration, p, q
      integer(i15), parameter :: max_iterations = 50
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
   subroutine construct_projection_matrices_hf(wf, Po, Pv)
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
      real(dp), dimension(:,:), allocatable :: tmp
!
      integer(i15) :: w, x
!
!     Po = 1/2 D S
!
      call dgemm('N','N',        &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  wf%n_ao,       &
                  half,          &
                  wf%ao_density, &
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
      integer(i15) :: p, q
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
   subroutine construct_roothan_hall_gradient_hf(wf, G, Po, Pv)
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
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: G
!
      real(dp), dimension(:, :), allocatable :: tmp
!
!     Construct tmp = Fov = Po^T F Pv and set G = tmp = Fov
!
      call mem%alloc(tmp, wf%n_ao, wf%n_ao)
!
      tmp = wf%ao_fock
      call sandwich(tmp, Po, Pv, wf%n_ao)
!
      G = tmp
!
!     Construct tmp = Fvo = Pv^T F Po and set H = H - tmp = Fov - Fvo
!
      tmp = wf%ao_fock
      call sandwich(tmp, Pv, Po, wf%n_ao)
!
      G = G - tmp
!
      call mem%dealloc(tmp, wf%n_ao, wf%n_ao)
!
   end subroutine construct_roothan_hall_gradient_hf
!
!
   subroutine do_roothan_hall_hf(wf, F, C, e, do_mo_transformation)
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
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in)    :: F 
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(inout) :: C 
      real(dp), dimension(wf%n_mo, 1),       intent(inout) :: e  
!
      logical, optional, intent(in) :: do_mo_transformation
!
      real(dp), dimension(:,:), allocatable :: work
      real(dp), dimension(:,:), allocatable :: metric 
      real(dp), dimension(:,:), allocatable :: ao_fock 
      real(dp), dimension(:,:), allocatable :: FP 
!
      real(dp), dimension(:,:), allocatable :: ao_fock_copy 
      real(dp), dimension(:,:), allocatable :: red_orbital_coefficients  
      real(dp), dimension(:,:), allocatable :: tmp  
!
      real(dp) :: ddot, norm
!
      integer(i15) :: info, i
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
      if (present(do_mo_transformation)) then 
!
         if (do_mo_transformation) then 
!
            call mem%alloc(ao_fock_copy, wf%n_mo, wf%n_mo)
            ao_fock_copy = ao_fock 
!
         endif 
!
      endif 
!  
!     Solve F'C' = L L^T C' e
!
      info = 0
!
      call mem%alloc(work, 4*wf%n_mo, 1)
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
      call mem%dealloc(work, 4*wf%n_mo, 1)
!
      if (info .ne. 0) then 
!
         write(output%unit, '(/t3,a/)') 'Error: could not solve Roothan-Hall equations.'
         stop
!
      endif
!
!     If requested MO transformation of Fock matrix, do it 
!
      if (present(do_mo_transformation)) then 
!
         if (do_mo_transformation) then 
!
            call mem%alloc(red_orbital_coefficients, wf%n_mo, wf%n_mo)
            red_orbital_coefficients = ao_fock 
!
            call mem%alloc(tmp, wf%n_mo, wf%n_mo)
!
            call dgemm('N', 'N', &
                        wf%n_mo, &
                        wf%n_mo, &
                        wf%n_mo, &
                        one, &
                        ao_fock_copy, &
                        wf%n_mo, &
                        red_orbital_coefficients, &
                        wf%n_mo, &
                        zero, &
                        tmp, & ! tmp = F' C'  
                        wf%n_mo)
!
            call wf%initialize_mo_fock() ! Allocate if necessary
!
            call dgemm('T', 'N',                  &
                        wf%n_mo,                  &
                        wf%n_mo,                  &
                        wf%n_mo,                  &
                        one,                      &
                        red_orbital_coefficients, &
                        wf%n_mo,                  &
                        tmp,                      &
                        wf%n_mo,                  &
                        zero,                     &
                        wf%mo_fock,               & ! F = C'^T F' C'  
                        wf%n_mo)
!
            call mem%dealloc(ao_fock_copy, wf%n_mo, wf%n_mo)
            call mem%dealloc(red_orbital_coefficients, wf%n_mo, wf%n_mo)
            call mem%dealloc(tmp, wf%n_mo, wf%n_mo)
!
         endif
!
      endif
!
!     Transform back the solutions to original basis, C = P (P^T C) = P C'
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
      integer(i15) :: I, n_ao_on_atom, first_ao_on_atom, last_ao_on_atom, n_s_on_atom
!
      real(dp) :: trace
!
      real(dp), dimension(:,:), allocatable :: atomic_density, temp
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
      integer(i15) :: ao 
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
   subroutine get_ao_h_wx_hf(wf, h)
!!
!!    Get AO h 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the full one-electron h matrix.
!!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: h 
!
      type(interval) :: A_interval, B_interval
!
      integer(i15) :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: h_AB 
!
!$omp parallel do &
!$omp private(A, B, h_AB, A_interval, B_interval, x, y) schedule(static)
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(h_AB, A_interval%size, B_interval%size)
            call wf%system%ao_integrals%construct_ao_h_wx(h_AB, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   h(A_interval%first - 1 + x, B_interval%first - 1 + y) = h_AB(x, y)
                   h(B_interval%first - 1 + y, A_interval%first - 1 + x) = h_AB(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(h_AB, A_interval%size, B_interval%size)
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_ao_h_wx_hf
!
!
   subroutine get_ao_s_wx_hf(wf, s)
!!
!!    Get AO s 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the full one-electron h matrix.
!!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: s 
!
      type(interval) :: A_interval, B_interval
!
      integer(i15) :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: s_AB 
!
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(s_AB, A_interval%size, B_interval%size)
            call wf%system%ao_integrals%construct_ao_s_wx(s_AB, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   s(A_interval%first - 1 + x, B_interval%first - 1 + y) = s_AB(x, y)
                   s(B_interval%first - 1 + y, A_interval%first - 1 + x) = s_AB(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(s_AB, A_interval%size, B_interval%size)
!
         enddo
      enddo
!
   end subroutine get_ao_s_wx_hf
!
!
   subroutine get_ao_mu_wx_hf(wf, mu_X, mu_Y, mu_Z)
!!
!!    Get AO mu
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Uses the integral tool to construct the full dipole integrals
!!    for the X, Y, and Z components.
!!
      implicit none 
!
      class(hf), intent(in) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: mu_X
      real(dp), dimension(wf%n_ao, wf%n_ao) :: mu_Y
      real(dp), dimension(wf%n_ao, wf%n_ao) :: mu_Z
!
      type(interval) :: A_interval, B_interval
!
      integer(i15) :: x, y, A, B
!
      real(dp), dimension(:,:), allocatable :: mu_AB_X 
      real(dp), dimension(:,:), allocatable :: mu_AB_Y 
      real(dp), dimension(:,:), allocatable :: mu_AB_Z 
!
      do A = 1, wf%system%n_s
!
         A_interval = wf%system%shell_limits(A)
!
         do B = 1, A
!
            B_interval = wf%system%shell_limits(B)
!
            call mem%alloc(mu_AB_X, A_interval%size, B_interval%size)
            call mem%alloc(mu_AB_Y, A_interval%size, B_interval%size)
            call mem%alloc(mu_AB_Z, A_interval%size, B_interval%size)
!
            call wf%system%ao_integrals%construct_ao_mu_wx(mu_AB_X, mu_AB_Y, mu_AB_Z, A, B)
!
             do x = 1, A_interval%size
                do y = 1, B_interval%size
!
                   mu_X(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_X(x, y)
                   mu_X(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_X(x, y)
!
                   mu_Y(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_Y(x, y)
                   mu_Y(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_Y(x, y)
!
                   mu_Z(A_interval%first - 1 + x, B_interval%first - 1 + y) = mu_AB_Z(x, y)
                   mu_Z(B_interval%first - 1 + y, A_interval%first - 1 + x) = mu_AB_Z(x, y)
!
                enddo
             enddo
!
            call mem%dealloc(mu_AB_X, A_interval%size, B_interval%size)
            call mem%dealloc(mu_AB_Y, A_interval%size, B_interval%size)
            call mem%dealloc(mu_AB_Z, A_interval%size, B_interval%size)
!
         enddo
      enddo
!
   end subroutine get_ao_mu_wx_hf
!
!
end module hf_class
!!
!      write(output%unit, *)'doing exchange'
!      flush(output%unit)
!!
!!$omp parallel do &
!!$omp private(A, B, C, A_interval, B_interval, C_interval, x, y, z, xy, zz, xz, yz, &
!!$omp g_K) schedule(dynamic)
!      do A = 1, n_s
!!
!         A_interval = wf%system%shell_limits(A)
!!           
!         do C = 1, n_s
!!
!            C_interval = wf%system%shell_limits(C)
!!
!            do B = 1, A
!!
!               B_interval = wf%system%shell_limits(B)
!!
!               if (sp_eri_schwarz(A, C)*sp_eri_schwarz(B, C)*sp_density_schwarz(C, 1) .lt. exchange_thr) cycle
!!
!               call mem%alloc(g_K, (A_interval%size)*(C_interval%size), &
!                                 (B_interval%size)*(C_interval%size))
!!
!               call wf%system%ao_integrals%construct_ao_g_wxyz(g_K, A, C, B, C)
!!
!!              Add Fock matrix contributions
!!
!               if (A .ne. B) then
!!
!                  do x = A_interval%first, A_interval%last
!                     do y = B_interval%first, B_interval%last
!!
!                        xy = A_interval%size*(y - B_interval%first) + x - A_interval%first + 1
!!
!                        do z = C_interval%first, C_interval%last
!!
!                           zz = C_interval%size*(z - C_interval%first) + z  - C_interval%first + 1
!                           xz = A_interval%size*(z - C_interval%first) + x  - A_interval%first + 1
!                           yz = B_interval%size*(z - C_interval%first) + y  - B_interval%first + 1
!!
!                           wf%ao_fock(x, y) = wf%ao_fock(x, y) + (- half*g_K(xz, yz))*wf%ao_density(z, z)
!!
!                        enddo
!                     enddo
!                  enddo
!!
!               else
!!
!                  do x = A_interval%first, A_interval%last
!                     do y = A_interval%first, x
!!
!                        xy = A_interval%size*(y - A_interval%first) + x - A_interval%first + 1
!!
!                        do z = C_interval%first, C_interval%last
!!
!                           zz = C_interval%size*(z - C_interval%first) + z  - C_interval%first + 1
!                           xz = A_interval%size*(z - C_interval%first) + x  - A_interval%first + 1
!                           yz = B_interval%size*(z - C_interval%first) + y  - B_interval%first + 1
!!
!                           wf%ao_fock(x, y) = wf%ao_fock(x, y) + (- half*g_K(xz, yz))*wf%ao_density(z, z)
!!
!                        enddo
!                     enddo
!                  enddo
!               endif
!!                  
!               call mem%dealloc(g_K, (A_interval%size)*(C_interval%size), &
!                                 (B_interval%size)*(C_interval%size))
!!
!            enddo
!         enddo
!      enddo
!!$omp end parallel do
