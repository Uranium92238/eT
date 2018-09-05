module uhf_class
!
!!
!!    Unrestricted Hartree-Fock (HF) class module
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!
   use hf_class
!
   use reordering
   use array_utilities
   use array_analysis
   use interval_class
   use index
!
   implicit none
!
   type, extends(hf):: uhf
!
      integer(i15) :: n_alpha
      integer(i15) :: n_beta 
!
      real(dp), dimension(:,:), allocatable :: ao_density_a 
      real(dp), dimension(:,:), allocatable :: ao_density_b
! 
      real(dp), dimension(:,:), allocatable :: ao_fock_a 
      real(dp), dimension(:,:), allocatable :: ao_fock_b
!
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_a
      real(dp), dimension(:,:), allocatable :: orbital_coefficients_b      
!
      real(dp), dimension(:,:), allocatable :: orbital_energies_a
      real(dp), dimension(:,:), allocatable :: orbital_energies_b
!
	contains
!
      procedure :: prepare => prepare_uhf
!
      procedure :: determine_n_alpha_and_beta_electrons => determine_n_alpha_and_beta_electrons_uhf
!
      procedure :: construct_ao_spin_density => construct_ao_spin_density_uhf
      procedure :: construct_ao_spin_fock    => construct_ao_spin_fock_uhf
!
      procedure :: form_ao_density => form_ao_density_uhf
!
      procedure :: calculate_uhf_energy => calculate_uhf_energy_uhf
      procedure :: set_ao_density_to_core_guess => set_ao_density_to_core_guess_uhf
!
!     Routines that may change in descendants but are required by solvers 
!
      procedure :: initialize_orbitals          => initialize_orbitals_uhf
      procedure :: initialize_density           => initialize_density_uhf
      procedure :: initialize_fock              => initialize_fock_uhf
      procedure :: destruct_fock                => destruct_fock_uhf 
      procedure :: update_fock_and_energy       => update_fock_and_energy_uhf
      procedure :: roothan_hall_update_orbitals => roothan_hall_update_orbitals_uhf
      procedure :: update_ao_density            => update_ao_density_uhf
!
!
!     Initialize and destruct routines
!
      procedure :: initialize_ao_density_a           => initialize_ao_density_a_uhf
      procedure :: initialize_ao_density_b           => initialize_ao_density_b_uhf
!
      procedure :: initialize_ao_fock_a              => initialize_ao_fock_a_uhf
      procedure :: initialize_ao_fock_b              => initialize_ao_fock_b_uhf
!
      procedure :: initialize_orbital_coefficients_a => initialize_orbital_coefficients_a_uhf
      procedure :: initialize_orbital_coefficients_b => initialize_orbital_coefficients_b_uhf
!
      procedure :: initialize_orbital_energies_a     => initialize_orbital_energies_a_uhf
      procedure :: initialize_orbital_energies_b     => initialize_orbital_energies_b_uhf
!
      procedure :: destruct_ao_density_a             => destruct_ao_density_a_uhf
      procedure :: destruct_ao_density_b             => destruct_ao_density_b_uhf
!
      procedure :: destruct_ao_fock_a                => destruct_ao_fock_a_uhf
      procedure :: destruct_ao_fock_b                => destruct_ao_fock_b_uhf
!
      procedure :: destruct_orbital_coefficients_a   => destruct_orbital_coefficients_a_uhf
      procedure :: destruct_orbital_coefficients_b   => destruct_orbital_coefficients_b_uhf
!
      procedure :: destruct_orbital_energies_a       => destruct_orbital_energies_a_uhf
      procedure :: destruct_orbital_energies_b       => destruct_orbital_energies_b_uhf
!
   end type uhf
!
!
contains 
!
!
   subroutine prepare_uhf(wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(uhf) :: wf
!
      wf%name = 'UHF'
!
      call wf%system%prepare()
!
      call get_n_aos(wf%n_ao)
!
      call initialize_coulomb()
      call initialize_kinetic()
      call initialize_nuclear()
      call initialize_overlap()
!
      call wf%initialize_ao_overlap()
      call wf%construct_ao_overlap()
      call wf%decompose_ao_overlap() 
!
      wf%n_o = (wf%system%get_n_electrons())/2
      wf%n_v = wf%n_mo - wf%n_o
!
      call wf%determine_n_alpha_and_beta_electrons()
!
      write(output%unit, *) ' n alpha:', wf%n_alpha
      write(output%unit, *) ' n beta:', wf%n_beta
      write(output%unit, *) 'multiplicity:', wf%system%multiplicity
!
   end subroutine prepare_uhf
!
!
   subroutine initialize_orbitals_uhf(wf)
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
      class(uhf) :: wf 
!
      call wf%initialize_orbital_coefficients()
!
      call wf%initialize_orbital_coefficients_a()
      call wf%initialize_orbital_coefficients_b()
!
   end subroutine initialize_orbitals_uhf
!
!
   subroutine initialize_density_uhf(wf)
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
      class(uhf) :: wf 
!
      call wf%initialize_ao_density()
!
      call wf%initialize_ao_density_a()
      call wf%initialize_ao_density_b()
!
   end subroutine initialize_density_uhf
!
!
   subroutine initialize_fock_uhf(wf)
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
      class(uhf) :: wf 
!
      call wf%initialize_ao_fock()
!
      call wf%initialize_ao_fock_a()
      call wf%initialize_ao_fock_b()
!
   end subroutine initialize_fock_uhf
!
!
   subroutine destruct_fock_uhf(wf)
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
      class(uhf) :: wf 
!
      call wf%destruct_ao_fock()
!
      call wf%destruct_ao_fock_a()
      call wf%destruct_ao_fock_b()
!
   end subroutine destruct_fock_uhf
!
!
   subroutine update_fock_and_energy_uhf(wf, sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx, coulomb, exchange, precision)
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
      class(uhf) :: wf 
!
      integer(i15), intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      real(dp), optional :: coulomb, exchange, precision ! Non-standard thresholds, optionals
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_a, 'alpha',    &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s, &
                                       h_wx, coulomb, exchange, precision)
!
      call wf%construct_ao_spin_fock(wf%ao_density, wf%ao_density_b, 'beta',     &
                                       sp_eri_schwarz, sp_eri_schwarz_list, n_s, &
                                       h_wx, coulomb, exchange, precision)     
!
      call wf%calculate_uhf_energy(h_wx)
!
   end subroutine update_fock_and_energy_uhf
!
!
   subroutine roothan_hall_update_orbitals_uhf(wf)
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
      class(uhf) :: wf 
!
      wf%ao_fock = wf%ao_fock_a 
      call wf%do_roothan_hall()
      wf%orbital_coefficients_a = wf%orbital_coefficients 
!
      wf%ao_fock = wf%ao_fock_b 
      call wf%do_roothan_hall()
      wf%orbital_coefficients_b = wf%orbital_coefficients 
!
   end subroutine roothan_hall_update_orbitals_uhf
!
!
   subroutine update_ao_density_uhf(wf)
!!
!!    Update AO density 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Updates the AO density (or densities, if unrestricted) based 
!!    on the current orbital coefficient matrix (or matrices).
!! 
      implicit none 
!
      class(uhf) :: wf 
!
      call wf%construct_ao_spin_density('alpha') ! Make D_alpha 
      call wf%construct_ao_spin_density('beta')  ! Make D_beta 
!
      call wf%form_ao_density()                  ! Set D = D_alpha + D_beta 
!
   end subroutine update_ao_density_uhf
!
!
   subroutine set_ao_density_to_core_guess_uhf(wf, h_wx)
!!
!!    Set AO density to core guess
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Solves the Roothan-Hall equation ignoring the two-
!!    electron AO density dependent part of the Fock matrix.
!!
!!    Based on the orbital coefficients, the routine constructs 
!!    the associated AO spin densities. 
!!
      implicit none 
!
      class(uhf) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx 
!
      wf%ao_fock = h_wx 
      call wf%do_roothan_hall()
!
      wf%orbital_coefficients_a = wf%orbital_coefficients
      wf%orbital_coefficients_b = wf%orbital_coefficients
!
      call wf%construct_ao_spin_density('alpha')
      call wf%construct_ao_spin_density('beta')
!
   end subroutine set_ao_density_to_core_guess_uhf
!
!
   subroutine determine_n_alpha_and_beta_electrons_uhf(wf)
!!
!!    Determine the number of alpha and beta electrons 
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
!!    Determines the number of alpha and beta elctrons 
!!    from the provided multiplicity and the number of
!!    electrons in the system. We assume that n_alpha 
!!    is greater than or equal to n_beta without loss
!!    of generality.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      integer(i15) :: n_alpha_m_n_beta, n_alpha_p_n_beta 
!
      n_alpha_m_n_beta = wf%system%multiplicity - 1 ! 2 S + 1 - 1 = 2 S
      n_alpha_p_n_beta = wf%system%n_electrons 
!
      wf%n_alpha = (n_alpha_p_n_beta + n_alpha_m_n_beta)/2
      wf%n_beta  = wf%system%n_electrons - wf%n_alpha
!
      if ((wf%n_alpha + wf%n_beta) .ne. wf%system%n_electrons) then 
!
         call output%error_msg('Given multiplicity and number of electrons is inconsistent.')
!
      endif 
!
   end subroutine determine_n_alpha_and_beta_electrons_uhf
!
!
   subroutine construct_ao_spin_density_uhf(wf, sigma)
!!
!!    Construct AO spin density 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Constructs either the alpha or beta density, 
!!    depending on the value of 'spin' (='alpha' or 'beta').
!!    The routine assumes that the associated orbital coefficients 
!!    are allocated and set. Typically, the construction of these
!!    densities will follow a Roothan-Hall step. 
!!
!!    We employ the normalization so that D = D_alpha + D_beta is 
!!    equal to the restricted HF density if the alpha and beta 
!!    densities are equal.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      character(len=*), intent(in) :: sigma
!
      integer(i15) :: x,y,i
!
      if (trim(sigma) == 'alpha') then 
!
         call dgemm('N', 'T',                   &
                     wf%n_ao,                   &
                     wf%n_ao,                   &
                     wf%n_alpha,                &
                     one,                       &
                     wf%orbital_coefficients_a, &
                     wf%n_ao,                   &
                     wf%orbital_coefficients_a, &
                     wf%n_ao,                   &
                     zero,                      &
                     wf%ao_density_a,           &
                     wf%n_ao)
!
      elseif (trim(sigma) == 'beta') then 
!
         call dgemm('N', 'T',                   &
                     wf%n_ao,                   &
                     wf%n_ao,                   &
                     wf%n_beta,                 &
                     one,                       &
                     wf%orbital_coefficients_b, &
                     wf%n_ao,                   &
                     wf%orbital_coefficients_b, &
                     wf%n_ao,                   &
                     zero,                      &
                     wf%ao_density_b,           &
                     wf%n_ao)   
!
      else
!
         call output%error_msg('Did not recognize spin variable in construct_ao_spin_density:' // trim(sigma))
!
      endif  
!
   end subroutine construct_ao_spin_density_uhf
!
!
   subroutine form_ao_density_uhf(wf)
!!
!!    Form AO density 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Given that the spin densities have been constructed,
!!    this routine sets 
!!
!!       D = D_alpha + D_beta     
!!
!!    Note that it is not constructed from the orbital
!!    coefficients and should be preceded by calls to 
!!    "construct_ao_spin_density" to make sure that the 
!!    spin densities are up-to-date.
!!
      implicit none 
!
      class(uhf) :: wf 
!
      wf%ao_density = wf%ao_density_a + wf%ao_density_b
!
   end subroutine form_ao_density_uhf
!
!
   subroutine construct_ao_spin_fock_uhf(wf, D, D_sigma, sigma, &
                     sp_eri_schwarz, sp_eri_schwarz_list, n_s, h_wx, coulomb, exchange, precision)
!!
!!    Construct AO spin Fock 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    The routine computes the alpha or beta Fock matrix, depending 
!!    on the value of the spin 'sigma' (='alpha' or 'beta'):
!!
!!       F_αβ^alpha = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^alpha 
!!       F_αβ^beta  = h_αβ + sum_γδ g_αβγδ D_γδ - sum_γδ g_αδγβ D_γδ^beta 
!!
!!    Here the superscript refers to the spin function, while the subscripts
!!    are AO indices. In contrast to the restricted routine, this one does 
!!    not calculate the energy - a separate call is required to get the 
!!    unrestricted Hartree-Fock energy. 
!!
      implicit none 
!
      class(uhf) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D 
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: D_sigma 
!
      character(len=*), intent(in) :: sigma 
!
      integer(i15), intent(in) :: n_s
!
      real(dp), dimension(wf%n_ao, wf%n_ao), intent(in) :: h_wx
!
      real(dp), dimension(n_s*(n_s + 1)/2, 2), intent(in)     :: sp_eri_schwarz
      integer(i15), dimension(n_s*(n_s + 1)/2, 3), intent(in) :: sp_eri_schwarz_list
!
      integer(i15) :: thread, n_threads, omp_get_max_threads
!
      real(dp), dimension(:,:), allocatable :: F, sp_density_schwarz
!
      real(dp), optional :: coulomb, exchange, precision      ! Non-standard thresholds, optionals
      real(dp) :: coulomb_thr, exchange_thr, precision_thr    ! Actual thresholds 
!
      integer(i15) :: n_sig_sp
!
      real(dp) :: max_D_schwarz, max_eri_schwarz
!
      real(dp), dimension(:,:), allocatable :: scaled_D_sigma ! = 2 * D_sigma
!
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
!     Compute number of significant ERI shell pairs (the Fock construction 
!     only loops over these shell pairs) and the maximum element 
!
      call wf%get_n_sig_eri_sp(n_sig_sp, sp_eri_schwarz, 1.0d-20)
      max_eri_schwarz = get_abs_max(sp_eri_schwarz, n_s*(n_s + 1)/2)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors 
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(sp_density_schwarz, n_s, n_s)
      call wf%construct_sp_density_schwarz(sp_density_schwarz, D)
      max_D_schwarz = get_abs_max(sp_density_schwarz, n_s**2)
!
      n_threads = omp_get_max_threads()
!
      call mem%alloc(F, wf%n_ao, wf%n_ao*n_threads) ! [F(thread 1) F(thread 2) ...]
      F = zero 
!
      call wf%ao_fock_coulomb_construction_loop(F, D, n_threads, max_D_schwarz, max_eri_schwarz,         & 
                                                sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list, &
                                                n_s, n_sig_sp, coulomb_thr, precision_thr,               &
                                                wf%system%shell_limits)
!
!     Construct the Coulomb two electron part of the Fock matrix, using the screening vectors 
!     and parallellizing over available threads (each gets its own copy of the Fock matrix)
!
      call mem%alloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
      scaled_D_sigma = two*D_sigma 
!
      call wf%construct_sp_density_schwarz(sp_density_schwarz, scaled_D_sigma)
      max_D_schwarz = get_abs_max(sp_density_schwarz, n_s**2)      
!
      call wf%ao_fock_exchange_construction_loop(F, scaled_D_sigma, n_threads, max_D_schwarz, max_eri_schwarz, & 
                                                   sp_density_schwarz, sp_eri_schwarz, sp_eri_schwarz_list,    &
                                                   n_s, n_sig_sp, exchange_thr, precision_thr,                 &
                                                   wf%system%shell_limits)      
!
      call mem%dealloc(sp_density_schwarz, n_s, n_s)
      call mem%dealloc(scaled_D_sigma, wf%n_ao, wf%n_ao)
!
!     Add the accumulated Fock matrix F into the correct Fock matrix
!     (i.e., either the alpha or beta Fock matrix )
!
      if (trim(sigma) == 'alpha') then 
!
         wf%ao_fock_a = zero
         do thread = 1, n_threads
!
            call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock_a, 1)
!
         enddo
!
         call symmetric_sum(wf%ao_fock_a, wf%n_ao)
         wf%ao_fock_a = wf%ao_fock_a*half
!
         wf%ao_fock_a = wf%ao_fock_a + h_wx
!
      elseif (trim(sigma) == 'beta') then 
!
         wf%ao_fock_b = zero
         do thread = 1, n_threads
!
            call daxpy(wf%n_ao**2, one, F(1, (thread-1)*wf%n_ao + 1), 1, wf%ao_fock_b, 1)
!
         enddo 
!
         call symmetric_sum(wf%ao_fock_b, wf%n_ao)
         wf%ao_fock_b = wf%ao_fock_b*half
!
         wf%ao_fock_b = wf%ao_fock_b + h_wx
!
      else 
!
         call output%error_msg('Did not recognize spin variable in construct_ao_fock:' // trim(sigma))
!
      endif 
!
      call mem%dealloc(F, wf%n_ao, wf%n_ao*n_threads) 
!
   end subroutine construct_ao_spin_fock_uhf
!
!
   subroutine calculate_uhf_energy_uhf(wf, h_wx)
!!
!!    Calculate UHF energy 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
!!    Calculates the unrestricted Hartree-Fock energy,
!!
!!       E = E_alpha + E_beta = 1/2 Tr (D^alpha h) + 1/2 Tr (D^alpha F) 
!!                            + 1/2 Tr (D^beta h)  + 1/2 Tr (D^beta F)
!!
      implicit none 
!
      class(uhf) :: wf 
!
      real(dp), dimension(wf%n_ao, wf%n_ao) :: h_wx
!
      real(dp) :: ddot 
!
      wf%energy = wf%system%get_nuclear_repulsion()
!
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_a, 1, h_wx, 1)
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_a, 1, wf%ao_fock_a, 1)    
!
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_b, 1, h_wx, 1)
      wf%energy = wf%energy + (one/two)*ddot((wf%n_ao)**2, wf%ao_density_b, 1, wf%ao_fock_b, 1)   
!
   end subroutine calculate_uhf_energy_uhf
!
!
   subroutine initialize_ao_density_a_uhf(wf)
!!
!!    Initialize AO density alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%ao_density_a)) call mem%alloc(wf%ao_density_a, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_a_uhf
!
!
   subroutine destruct_ao_density_a_uhf(wf)
!!
!!    Destruct AO density alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%ao_density_a)) call mem%dealloc(wf%ao_density_a, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_a_uhf
!
!
   subroutine initialize_ao_density_b_uhf(wf)
!!
!!    Initialize AO density beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%ao_density_b)) call mem%alloc(wf%ao_density_b, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_density_b_uhf
!
!
   subroutine destruct_ao_density_b_uhf(wf)
!!
!!    Destruct AO density beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%ao_density_b)) call mem%dealloc(wf%ao_density_b, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_density_b_uhf
!
!
   subroutine initialize_ao_fock_a_uhf(wf)
!!
!!    Initialize AO Fock alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%ao_fock_a)) call mem%alloc(wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_fock_a_uhf
!
!
   subroutine destruct_ao_fock_a_uhf(wf)
!!
!!    Destruct AO Fock alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%ao_fock_a)) call mem%dealloc(wf%ao_fock_a, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_a_uhf
!
!
   subroutine initialize_ao_fock_b_uhf(wf)
!!
!!    Initialize AO Fock beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%ao_fock_b)) call mem%alloc(wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
   end subroutine initialize_ao_fock_b_uhf
!
!
   subroutine destruct_ao_fock_b_uhf(wf)
!!
!!    Destruct AO Fock beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%ao_fock_b)) call mem%dealloc(wf%ao_fock_b, wf%n_ao, wf%n_ao)
!
   end subroutine destruct_ao_fock_b_uhf
!
!
   subroutine initialize_orbital_coefficients_a_uhf(wf)
!!
!!    Initialize orbital coefficients alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%orbital_coefficients_a)) call mem%alloc(wf%orbital_coefficients_a, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_a_uhf
!
!
   subroutine destruct_orbital_coefficients_a_uhf(wf)
!!
!!    Destruct orbital coefficients alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%orbital_coefficients_a)) call mem%dealloc(wf%orbital_coefficients_a, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_a_uhf
!
!
   subroutine initialize_orbital_coefficients_b_uhf(wf)
!!
!!    Initialize orbital coefficients beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%orbital_coefficients_b)) call mem%alloc(wf%orbital_coefficients_b, wf%n_ao, wf%n_mo)
!
   end subroutine initialize_orbital_coefficients_b_uhf
!
!
   subroutine destruct_orbital_coefficients_b_uhf(wf)
!!
!!    Destruct orbital coefficients beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%orbital_coefficients_b)) call mem%dealloc(wf%orbital_coefficients_b, wf%n_ao, wf%n_mo)
!
   end subroutine destruct_orbital_coefficients_b_uhf
!
!
   subroutine initialize_orbital_energies_a_uhf(wf)
!!
!!    Initialize orbital energies alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%orbital_energies_a)) call mem%alloc(wf%orbital_energies_a, wf%n_mo, 1)
!
   end subroutine initialize_orbital_energies_a_uhf
!
!
   subroutine destruct_orbital_energies_a_uhf(wf)
!!
!!    Initialize orbital energies alpha 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%orbital_energies_a)) call mem%dealloc(wf%orbital_energies_a, wf%n_mo, 1)
!
   end subroutine destruct_orbital_energies_a_uhf
!
!
   subroutine initialize_orbital_energies_b_uhf(wf)
!!
!!    Initialize orbital energies beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (.not. allocated(wf%orbital_energies_b)) call mem%alloc(wf%orbital_energies_b, wf%n_mo, 1)
!
   end subroutine initialize_orbital_energies_b_uhf
!
!
   subroutine destruct_orbital_energies_b_uhf(wf)
!!
!!    Destruct orbital energies beta 
!!    Written by Eirik F. Kjønstad, Sep 2018
!!
      implicit none
!
      class(uhf) :: wf 
!
      if (allocated(wf%orbital_energies_b)) call mem%dealloc(wf%orbital_energies_b, wf%n_mo, 1)
!
   end subroutine destruct_orbital_energies_b_uhf
!
!
end module uhf_class
