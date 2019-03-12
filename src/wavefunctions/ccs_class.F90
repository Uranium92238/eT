module ccs_class
!
!!
!!    Coupled cluster singles (ccs) class module
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use wavefunction_class
   use hf_class
!
   use mo_integral_tool_class
!
   use reordering
   use array_utilities
   use array_analysis
   use interval_class
   use index_invert
   use batching_index_class
   use timings_class
!
   implicit none
!
   type, extends(wavefunction) :: ccs
!
      real(dp) :: hf_energy 
!
      integer                                :: n_gs_amplitudes
      integer                                :: n_es_amplitudes
      integer                                :: n_t1
!
      real(dp), dimension(:,:), allocatable  :: t1
      real(dp), dimension(:,:), allocatable  :: t1bar
!
      real(dp), dimension(:,:), allocatable  :: fock_ij
      real(dp), dimension(:,:), allocatable  :: fock_ia
      real(dp), dimension(:,:), allocatable  :: fock_ai
      real(dp), dimension(:,:), allocatable  :: fock_ab
!
      real(dp), dimension(:), allocatable  :: fock_diagonal
!
      type(mo_integral_tool) :: integrals
!
      integer :: n_bath ! Number of bath orbitals (always the last ao/mo indices)
!
   contains
!
!     Preparation and cleanup routines 
!
      procedure :: prepare                                     => prepare_ccs
      procedure :: cleanup                                     => cleanup_ccs
!
!     Routines related to the amplitudes & multipliers
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_ccs 
      procedure :: destruct_amplitudes                         => destruct_amplitudes_ccs 
      procedure :: set_initial_amplitudes_guess                => set_initial_amplitudes_guess_ccs
      procedure :: t1_transform                                => t1_transform_ccs
      procedure :: set_amplitudes                              => set_amplitudes_ccs 
      procedure :: get_amplitudes                              => get_amplitudes_ccs 
      procedure :: save_amplitudes                             => save_amplitudes_ccs
      procedure :: read_amplitudes                             => read_amplitudes_ccs
      procedure :: save_t1                                     => save_t1_ccs 
      procedure :: read_t1                                     => read_t1_ccs 
      procedure :: print_dominant_x_amplitudes                 => print_dominant_x_amplitudes_ccs
      procedure :: print_dominant_amplitudes                   => print_dominant_amplitudes_ccs
      procedure :: print_dominant_x1                           => print_dominant_x1_ccs
      procedure :: get_t1_diagnostic                           => get_t1_diagnostic_ccs
!
      procedure :: initialize_multipliers                      => initialize_multipliers_ccs
      procedure :: destruct_multipliers                        => destruct_multipliers_ccs
      procedure :: set_multipliers                             => set_multipliers_ccs
      procedure :: get_multipliers                             => get_multipliers_ccs
      procedure :: save_multipliers                            => save_multipliers_ccs
      procedure :: read_multipliers                            => read_multipliers_ccs
      procedure :: save_t1bar                                  => save_t1bar_ccs
      procedure :: read_t1bar                                  => read_t1bar_ccs
!
!     Routines related to the Fock matrix 
! 
      procedure :: set_fock                                    => set_fock_ccs
      procedure :: construct_fock                              => construct_fock_ccs
      procedure :: get_gs_orbital_differences                  => get_gs_orbital_differences_ccs
      procedure :: get_es_orbital_differences                  => get_gs_orbital_differences_ccs
      procedure :: calculate_energy                            => calculate_energy_ccs
!
!     Routines related to the omega vector 
!
      procedure :: construct_omega                             => construct_omega_ccs
      procedure :: omega_ccs_a1                                => omega_ccs_a1_ccs
!
!     Routines related to the Jacobian transformation 
!
      procedure :: prepare_for_jacobian                        => prepare_for_jacobian_ccs
!
      procedure :: jacobian_transform_trial_vector             => jacobian_transform_trial_vector_ccs
      procedure :: jacobian_transpose_transform_trial_vector   => jacobian_transpose_transform_trial_vector_ccs
!
      procedure :: jacobian_ccs_transformation                 => jacobian_ccs_transformation_ccs
      procedure :: jacobian_ccs_a1                             => jacobian_ccs_a1_ccs 
      procedure :: jacobian_ccs_b1                             => jacobian_ccs_b1_ccs 
!
      procedure :: jacobian_transpose_ccs_transformation       => jacobian_transpose_ccs_transformation_ccs
      procedure :: jacobian_transpose_ccs_a1                   => jacobian_transpose_ccs_a1_ccs
      procedure :: jacobian_transpose_ccs_b1                   => jacobian_transpose_ccs_b1_ccs
!
      procedure :: construct_excited_state_equation            => construct_excited_state_equation_ccs
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_ccs
      procedure :: construct_eta                               => construct_eta_ccs
!
      procedure :: get_cvs_projector                           => get_cvs_projector_ccs
      procedure :: get_ip_projector                            => get_ip_projector_ccs
!
      procedure :: set_cvs_start_indices                       => set_cvs_start_indices_ccs
!
!     Routines to get electron repulsion integrals (ERIs)
!
      procedure :: get_ovov                                     => get_ovov_ccs
      procedure :: get_vovo                                     => get_vovo_ccs
      procedure :: get_vvoo                                     => get_vvoo_ccs
      procedure :: get_voov                                     => get_voov_ccs
      procedure :: get_ovvo                                     => get_ovvo_ccs
      procedure :: get_oovv                                     => get_oovv_ccs
      procedure :: get_oooo                                     => get_oooo_ccs
      procedure :: get_vvvv                                     => get_vvvv_ccs
      procedure :: get_ooov                                     => get_ooov_ccs
      procedure :: get_oovo                                     => get_oovo_ccs
      procedure :: get_ovoo                                     => get_ovoo_ccs
      procedure :: get_vooo                                     => get_vooo_ccs
      procedure :: get_vvvo                                     => get_vvvo_ccs
      procedure :: get_vvov                                     => get_vvov_ccs
      procedure :: get_vovv                                     => get_vovv_ccs
      procedure :: get_ovvv                                     => get_ovvv_ccs
!
      procedure, nopass :: need_g_abcd                          => need_g_abcd_ccs
!
!     Routines to initialize and destruct arrays 
!
      procedure :: initialize_fock_ij                           => initialize_fock_ij_ccs
      procedure :: initialize_fock_ia                           => initialize_fock_ia_ccs
      procedure :: initialize_fock_ai                           => initialize_fock_ai_ccs
      procedure :: initialize_fock_ab                           => initialize_fock_ab_ccs
      procedure :: initialize_fock_diagonal                     => initialize_fock_diagonal_ccs
      procedure :: initialize_t1                                => initialize_t1_ccs
      procedure :: initialize_t1bar                             => initialize_t1bar_ccs
!
      procedure :: destruct_fock_ij                             => destruct_fock_ij_ccs
      procedure :: destruct_fock_ia                             => destruct_fock_ia_ccs
      procedure :: destruct_fock_ai                             => destruct_fock_ai_ccs
      procedure :: destruct_fock_ab                             => destruct_fock_ab_ccs
      procedure :: destruct_fock_diagonal                       => destruct_fock_diagonal_ccs
      procedure :: destruct_t1                                  => destruct_t1_ccs
      procedure :: destruct_t1bar                               => destruct_t1bar_ccs
!
   end type ccs
!
!
contains
!
!
   subroutine prepare_ccs(wf, ref_wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      class(hf) :: ref_wf
!
      integer :: p
!
      wf%name_ = 'ccs'
!
      wf%system = ref_wf%system
!
      wf%n_ao   = ref_wf%n_ao
      wf%n_mo   = ref_wf%n_mo
      wf%n_o    = ref_wf%n_o
      wf%n_v    = ref_wf%n_v
!
      wf%hf_energy = ref_wf%energy
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
      call wf%initialize_fock_diagonal()
!
      wf%fock_ij(:,:) = ref_wf%mo_fock(1 : wf%n_o, 1 : wf%n_o)
      wf%fock_ia(:,:) = ref_wf%mo_fock(1 : wf%n_o, wf%n_o + 1 : wf%n_mo)
      wf%fock_ai(:,:) = ref_wf%mo_fock(wf%n_o + 1 : wf%n_mo, 1 : wf%n_o)
      wf%fock_ab(:,:) = ref_wf%mo_fock(wf%n_o + 1 : wf%n_mo, wf%n_o + 1 : wf%n_mo)
!
      do p = 1, wf%n_mo
!
         wf%fock_diagonal(p) = ref_wf%mo_fock(p, p)
!
      enddo
!
      call ref_wf%mo_transform_and_save_h()
!
      call wf%initialize_orbital_coefficients()
      wf%orbital_coefficients = ref_wf%orbital_coefficients
!
   end subroutine prepare_ccs
!
!
   subroutine cleanup_ccs(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      write(output%unit, '(/t3,a,a,a)') '- Cleaning up ', trim(wf%name_), ' wavefunction'
!
   end subroutine cleanup_ccs
!
!
   subroutine initialize_amplitudes_ccs(wf)
!!
!!    Initialize amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Allocates the amplitudes. This routine must be overwritten in 
!!    descendants which have more amplitudes. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      call wf%initialize_t1()
!
   end subroutine initialize_amplitudes_ccs
!
!
   logical function need_g_abcd_ccs()
!!
!!    Need g_abcd 
!!    Written by Eirik F. Kjønstad, Mar 2019 
!!
!!    Returns whether the vvvv-part of the ERI matrix 
!!    is used to calculate the ground and/or excited state 
!!    equations. If not, there is no need to compute the
!!    entire ERI matrix and store it in memory.
!!
      implicit none 
!
      need_g_abcd_ccs = .false.
!
   end function need_g_abcd_ccs
!
!
   subroutine destruct_amplitudes_ccs(wf)
!!
!!    Destruct amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Deallocates the amplitudes. This routine must be overwritten in 
!!    descendants which have more amplitudes. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      call wf%destruct_t1()
!
   end subroutine destruct_amplitudes_ccs
!
!
   subroutine set_amplitudes_ccs(wf, amplitudes)
!!
!!    Set amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(ccs) :: wf  
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: amplitudes
!
      call dcopy(wf%n_gs_amplitudes, amplitudes, 1, wf%t1, 1)
!
   end subroutine set_amplitudes_ccs
!
!
   subroutine get_amplitudes_ccs(wf, amplitudes)
!!
!!    Get amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none 
!
      class(ccs), intent(in) :: wf  
!
      real(dp), dimension(wf%n_gs_amplitudes) :: amplitudes
!
      call dcopy(wf%n_gs_amplitudes, wf%t1, 1, amplitudes, 1)
!
   end subroutine get_amplitudes_ccs
!
!
   subroutine save_amplitudes_ccs(wf)
!!
!!    Save amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      call wf%save_t1()
!
   end subroutine save_amplitudes_ccs
!
!
   subroutine read_amplitudes_ccs(wf)
!!
!!    Read amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf
!
      call wf%read_t1()    
!
   end subroutine read_amplitudes_ccs
!
!
   subroutine read_t1_ccs(wf)
!!
!!    Read t1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf
!
      type(file) :: t1_file 
!
      call t1_file%init('t1', 'sequential', 'unformatted')
!
      call disk%open_file(t1_file, 'read', 'rewind')
!
      read(t1_file%unit) wf%t1
!
      call disk%close_file(t1_file)      
!
   end subroutine read_t1_ccs
!
!
   subroutine save_t1_ccs(wf)
!!
!!    Save t1 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(in) :: wf 
!
      type(file) :: t1_file 
!
      call t1_file%init('t1', 'sequential', 'unformatted')
!
      call disk%open_file(t1_file, 'write', 'rewind')
!
      write(t1_file%unit) wf%t1
!
      call disk%close_file(t1_file)
!
   end subroutine save_t1_ccs
!
!
   subroutine save_multipliers_ccs(wf)
!!
!!    Save multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      call wf%save_t1bar()
!
   end subroutine save_multipliers_ccs
!
!
   subroutine read_multipliers_ccs(wf)
!!
!!    Read multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
      call wf%read_t1bar()
!
   end subroutine read_multipliers_ccs
!
!
   subroutine destruct_multipliers_ccs(wf)
!!
!!    Destruct multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Deallocates the multipliers. This routine must be overwritten in 
!!    descendants which have more multipliers. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      call wf%destruct_t1bar()
!
   end subroutine destruct_multipliers_ccs
!
!
   subroutine save_t1bar_ccs(wf)
!!
!!    Save t1bar 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      type(file) :: t1bar_file 
!
      call t1bar_file%init('t1bar', 'sequential', 'unformatted')
!
      call disk%open_file(t1bar_file, 'write', 'rewind')
!
      write(t1bar_file%unit) wf%t1bar
!
      call disk%close_file(t1bar_file)      
!
   end subroutine save_t1bar_ccs
!
!
   subroutine read_t1bar_ccs(wf)
!!
!!    Save t1bar 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
      type(file) :: t1bar_file 
!
      call t1bar_file%init('t1bar', 'sequential', 'unformatted')
!
      call disk%open_file(t1bar_file, 'read', 'rewind')
!
      read(t1bar_file%unit) wf%t1bar
!
      call disk%close_file(t1bar_file)      
!
   end subroutine read_t1bar_ccs
!
!
   subroutine set_initial_amplitudes_guess_ccs(wf)
!!
!!    Set initial amplitudes guess 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      wf%t1 = zero 
!
   end subroutine set_initial_amplitudes_guess_ccs
!
!
   subroutine initialize_multipliers_ccs(wf)
!!
!!    Initialize multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Allocates the multipliers. This routine must be overwritten in 
!!    descendants which have more multipliers. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      call wf%initialize_t1bar()
!
   end subroutine initialize_multipliers_ccs
!
!
   subroutine set_multipliers_ccs(wf, multipliers)
!!
!!    Set multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(ccs) :: wf  
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(in) :: multipliers
!
      call dcopy(wf%n_gs_amplitudes, multipliers, 1, wf%t1bar, 1)
!
   end subroutine set_multipliers_ccs
!
!
   subroutine get_multipliers_ccs(wf, multipliers)
!!
!!    Get multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none 
!
      class(ccs), intent(in) :: wf  
!
      real(dp), dimension(wf%n_gs_amplitudes) :: multipliers
!
      call dcopy(wf%n_gs_amplitudes, wf%t1bar, 1, multipliers, 1)
!
   end subroutine get_multipliers_ccs
!
!
   subroutine calculate_energy_ccs(wf)
!!
!!    Calculate energy 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
      integer :: i = 0
!
      type(file)   :: h_pq_file
!
      real(dp), dimension(:,:), allocatable :: h_pq
!
!     Read MO-transformed h array
!
      call h_pq_file%init('h_pq', 'sequential', 'unformatted')
      call disk%open_file(h_pq_file, 'read')
      rewind(h_pq_file%unit)
!
      call mem%alloc(h_pq, wf%n_mo, wf%n_mo)
      read(h_pq_file%unit) h_pq 
!
      call disk%close_file(h_pq_file)
!
!     Compute energy
!
      wf%energy = wf%system%get_nuclear_repulsion()
!
      do i = 1, wf%n_o
!
         wf%energy = wf%energy + h_pq(i,i) + wf%fock_ij(i,i)
!
      enddo
!
      call mem%dealloc(h_pq, wf%n_mo, wf%n_mo)
!
   end subroutine calculate_energy_ccs
!
!
   subroutine omega_ccs_a1_ccs(wf, omega)
!!
!!    Omega A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, March 2017
!!
!!    Adds the A1 contribution to omega,
!!
!!       Omega_ai^A1 =+ F_ai_T1.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes) :: omega
!
      type(timings) :: omega_ccs_a1_timer
!
      call omega_ccs_a1_timer%init('omega ccs a1')
      call omega_ccs_a1_timer%start()
!
      call daxpy((wf%n_o)*(wf%n_v), one, wf%fock_ai, 1, omega, 1)
!
      call omega_ccs_a1_timer%freeze()
      call omega_ccs_a1_timer%switch_off()
!
   end subroutine omega_ccs_a1_ccs
!
!
   subroutine construct_omega_ccs(wf, omega)
!!
!!    Construct Omega (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: omega
!
      omega = zero
      call wf%omega_ccs_a1(omega)
!
   end subroutine construct_omega_ccs
!
!
   subroutine construct_fock_ccs(wf)
!!
!!    Construct Fock 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 
!!
!!    Constructs the Fock matrix in the t1-transformed MO 
!!    basis using the MO integrals and the current single 
!!    amplitudes:
!!
!!       F_pq = h_pq + sum_k (2*g_pqkk - g_pkkq)
!!
!!    Since the two-electron ERIs are available already 
!!    t1-transformed, our task is to transform the one-
!!    electron term, which we assume is on file in the 
!!    MO basis.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      type(file) :: h_pq_file
!
      real(dp), dimension(:,:), allocatable :: F_pq 
!
      integer :: i, j, k, a, b
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ijkl
      real(dp), dimension(:,:,:,:), allocatable :: g_abij
      real(dp), dimension(:,:,:,:), allocatable :: g_aijb
      real(dp), dimension(:,:,:,:), allocatable :: g_iajk
      real(dp), dimension(:,:,:,:), allocatable :: g_aijk
!
!     Read MO-transformed h integrals into the  
!
      call h_pq_file%init('h_pq', 'sequential', 'unformatted')
      call disk%open_file(h_pq_file, 'read', 'rewind')
!
      call mem%alloc(F_pq, wf%n_mo, wf%n_mo)
      read(h_pq_file%unit) F_pq 
!
      call disk%close_file(h_pq_file)
!
!     Perform t1-transformation of F_pq = h_pq  
!
      call wf%t1_transform(F_pq)
!
!     Occupied-occupied contributions: F_ij = F_ij + sum_k (2*g_ijkk - g_ikkj)
!
      call mem%alloc(g_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_oooo(g_ijkl)
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o 
            do k = 1, wf%n_o
!
               F_pq(i, j) = F_pq(i, j) + two*g_ijkl(i,j,k,k) - g_ijkl(i,k,k,j)
! 
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_ijkl, wf%n_o, wf%n_o, wf%n_o, wf%n_o)
!
!     Occupied-virtual contributions: F_ia = F_ia + sum_j (2*g_iajj - g_ijja)
!                                     F_ai = F_ai + sum_j (2*g_aijj - g_ajji)
!
      call mem%alloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_ovoo(g_iajk)
!
      call mem%alloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
      call wf%get_vooo(g_aijk)
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
            do j = 1, wf%n_o
!
               F_pq(i, a + wf%n_o) = F_pq(i, a + wf%n_o) + two*g_iajk(i,a,j,j) - g_iajk(j,a,i,j)
               F_pq(a + wf%n_o, i) = F_pq(a + wf%n_o, i) + two*g_aijk(a,i,j,j) - g_aijk(a,j,j,i)
!
            enddo
!
         enddo
      enddo
!
      call mem%dealloc(g_iajk, wf%n_o, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_aijk, wf%n_v, wf%n_o, wf%n_o, wf%n_o)
!
!     Virtual-virtual contributions: F_ab = h_ab + sum_i (2*g_abii - g_aiib) ::
!
      call mem%alloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call wf%get_vvoo(g_abij)
!
      call mem%alloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov(g_aijb)
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
            do i = 1, wf%n_o
!
               F_pq(wf%n_o + a, wf%n_o + b) = F_pq(wf%n_o + a, wf%n_o + b) + two*g_abij(a,b,i,i) - g_aijb(a,i,i,b)
!
            enddo
         enddo
      enddo
!
      call mem%dealloc(g_abij, wf%n_v, wf%n_v, wf%n_o, wf%n_o)
      call mem%dealloc(g_aijb, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
      call wf%set_fock(F_pq)
      call mem%dealloc(F_pq, wf%n_mo, wf%n_mo)
!
   end subroutine construct_fock_ccs
!
!
   subroutine set_fock_ccs(wf, F_pq)
!!
!!    Set Fock 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Sets the different blocks of the Fock matrix based on the full 
!!    matrix sent to the routine.
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(in) :: F_pq 
!
      integer :: i, j, a, b 
!
      do i = 1, wf%n_o
         do j = 1, wf%n_o
!
            wf%fock_ij(i,j) = F_pq(i,j)
!
         enddo
      enddo
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            wf%fock_ia(i,a) = F_pq(i, wf%n_o + a)
            wf%fock_ai(a,i) = F_pq(wf%n_o + a, i)
!
         enddo
      enddo
!
      do a = 1, wf%n_v
         do b = 1, wf%n_v
!
            wf%fock_ab(a,b) = F_pq(wf%n_o + a, wf%n_o + b)
!
         enddo
      enddo      
!
   end subroutine set_fock_ccs
!
!
   subroutine t1_transform_ccs(wf, Z_pq)
!!
!!    T1 transform 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
!!    Assumes that Z is in the MO basis and performs the T1 transformation,
!!
!!       Z_pq <- sum_rs X_ps Z_sr Y_qr,    i.e.    Z <- X Z Y^T 
!!
!!    where
!!
!!       X = I - t1 
!!       Y = I + t1^T 
!! 
!!    Here, t1 is a full MO matrix whose only non-zero block is the vir-occ 
!!    part, where it is equal to t_i^a.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo), intent(inout) :: Z_pq 
!
      real(dp), dimension(:,:), allocatable :: X, Y 
!
      real(dp), dimension(:,:), allocatable :: W ! W_sq = sum_r Z_sr Y_rq^T, intermediate 
!
      integer :: p, i, a
!
!     Construct the X and Y arrays 
!
      call mem%alloc(X, wf%n_mo, wf%n_mo)
      call mem%alloc(Y, wf%n_mo, wf%n_mo)
!
      X = zero
      Y = zero
!
      do p = 1, wf%n_mo 
!
         X(p, p) = one 
         Y(p, p) = one 
!
      enddo 
!
      do i = 1, wf%n_o 
         do a = 1, wf%n_v 
!
            X(wf%n_o + a, i) = -wf%t1(a, i)
            Y(i, wf%n_o + a) = wf%t1(a, i) 
!
         enddo
      enddo
!
!     Construct intermediate W = Z Y^T and then use it to do transformation 
!
      call mem%alloc(W, wf%n_mo, wf%n_mo)
!
      call dgemm('N', 'T', &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  Z_pq,    & ! Z_s_r 
                  wf%n_mo, &
                  Y,       & ! Y_q_r 
                  wf%n_mo, &
                  zero,    &
                  W,       & ! W_sq = sum_r Z_sr Y_rq 
                  wf%n_mo) 
!
      call dgemm('N', 'N', &
                  wf%n_mo, &
                  wf%n_mo, &
                  wf%n_mo, &
                  one,     &
                  X,       &
                  wf%n_mo, &
                  W,       &
                  wf%n_mo, &
                  zero,    &
                  Z_pq,    & ! Z_pq = (X W)_pq = sum_s X_ps W_sq = sum_sr X_ps Z_sr Y_rq
                  wf%n_mo)
!
      call mem%dealloc(X, wf%n_mo, wf%n_mo)
      call mem%dealloc(Y, wf%n_mo, wf%n_mo)
      call mem%dealloc(W, wf%n_mo, wf%n_mo)
!
   end subroutine t1_transform_ccs
!
!
   subroutine get_gs_orbital_differences_ccs(wf, orbital_differences, N)
!!
!!    Get orbital differences 
!!    Written by Sarai D. Folkestad, Sep 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      integer, intent(in) :: N 
      real(dp), dimension(N), intent(inout) :: orbital_differences
!
      integer :: a, i, ai
!
      do i = 1, wf%n_o 
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai) = wf%fock_diagonal(a + wf%n_o) - wf%fock_diagonal(i)
!
         enddo
      enddo
!
   end subroutine get_gs_orbital_differences_ccs
!
!
   subroutine initialize_fock_ij_ccs(wf)
!!
!!    Initialize Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ij)) call mem%alloc(wf%fock_ij, wf%n_o, wf%n_o)
!
   end subroutine initialize_fock_ij_ccs
!
!
   subroutine initialize_fock_ia_ccs(wf)
!!
!!    Initialize Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ia)) call mem%alloc(wf%fock_ia, wf%n_o, wf%n_v)
!
   end subroutine initialize_fock_ia_ccs
!
!
   subroutine initialize_fock_ai_ccs(wf)
!!
!!    Initialize Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ai)) call mem%alloc(wf%fock_ai, wf%n_v, wf%n_o)
!
   end subroutine initialize_fock_ai_ccs
!
!
   subroutine initialize_fock_ab_ccs(wf)
!!
!!    Initialize Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_ab)) call mem%alloc(wf%fock_ab, wf%n_v, wf%n_v)
!
   end subroutine initialize_fock_ab_ccs
!
!
   subroutine initialize_fock_diagonal_ccs(wf)
!!
!!    Initialize Fock diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%fock_diagonal)) call mem%alloc(wf%fock_diagonal, wf%n_mo)
!
   end subroutine initialize_fock_diagonal_ccs
!
!
   subroutine initialize_t1_ccs(wf)
!!
!!    Initialize T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%t1)) call mem%alloc(wf%t1, wf%n_v, wf%n_o)
      wf%t1 = zero ! Hack, fix later, for integrals
!
   end subroutine initialize_t1_ccs
!
!
   subroutine initialize_t1bar_ccs(wf)
!!
!!    Initialize T1-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (.not. allocated(wf%t1bar)) call mem%alloc(wf%t1bar, wf%n_v, wf%n_o)
!
   end subroutine initialize_t1bar_ccs
!
!
   subroutine destruct_fock_ij_ccs(wf)
!!
!!    Destruct Fock ij block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ij)) call mem%dealloc(wf%fock_ij, wf%n_o, wf%n_o)
!
   end subroutine destruct_fock_ij_ccs
!
!
   subroutine destruct_fock_ia_ccs(wf)
!!
!!    Destruct Fock ia block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ia)) call mem%dealloc(wf%fock_ia, wf%n_o, wf%n_v)
!
   end subroutine destruct_fock_ia_ccs
!
!
   subroutine destruct_fock_ai_ccs(wf)
!!
!!    Destruct Fock ai block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ai)) call mem%dealloc(wf%fock_ai, wf%n_v, wf%n_o)
!
   end subroutine destruct_fock_ai_ccs
!
!
   subroutine destruct_fock_ab_ccs(wf)
!!
!!    Destruct Fock ab block
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_ab)) call mem%dealloc(wf%fock_ij, wf%n_v, wf%n_v)
!
   end subroutine destruct_fock_ab_ccs
!
!
   subroutine destruct_fock_diagonal_ccs(wf)
!!
!!    Destruct Fock diagonal
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%fock_diagonal)) call mem%dealloc(wf%fock_diagonal, wf%n_mo)
!
   end subroutine destruct_fock_diagonal_ccs
!
!
   subroutine destruct_t1_ccs(wf)
!!
!!    Destruct T1
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1)) call mem%dealloc(wf%t1, wf%n_v, wf%n_o)
!
   end subroutine destruct_t1_ccs
!
!
   subroutine destruct_t1bar_ccs(wf)
!!
!!    Destruct T1-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs) :: wf
!
      if (allocated(wf%t1bar)) call mem%dealloc(wf%t1bar, wf%n_v, wf%n_o)
!
   end subroutine destruct_t1bar_ccs
!
!
   subroutine get_ovov_ccs(wf, g_iajb, first_i, last_i, first_a, last_a, &
                                         first_j, last_j, first_b, last_b)
!!
!!    Get ovov 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iajb 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_a, local_last_a
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
!
      logical :: index_restrictions
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_a) .and. present(last_a) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b)) then
!
         index_restrictions = .true.
!
         local_first_i = first_i 
         local_first_a = first_a 
         local_first_j = first_j 
         local_first_b = first_b 
!
         local_last_i = last_i 
         local_last_a = last_a
         local_last_j = last_j
         local_last_b = last_b
!
      else
!
         index_restrictions = .false.
!
         local_first_i = 1 
         local_first_a = 1 
         local_first_j = 1 
         local_first_b = 1 
!
         local_last_i = wf%n_o 
         local_last_a = wf%n_v
         local_last_j = wf%n_o
         local_last_b = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_iajb, &
                                             local_first_i, local_last_i, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             local_first_j, local_last_j, &
                                             wf%n_o + local_first_b, wf%n_o + local_last_b)
!
   end subroutine get_ovov_ccs
!
!
   subroutine get_oooo_ccs(wf, g_ijkl, first_i, last_i, first_j, last_j, &
                                         first_k, last_k, first_l, last_l)
!!
!!    Get oooo 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijkl 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_l, last_l
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_l, local_last_l
!
      logical :: index_restrictions
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_l) .and. present(last_l)) then
!
         index_restrictions = .true.
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_k = first_k 
         local_first_l = first_l 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_k = last_k
         local_last_l = last_l
!
      else
!
         index_restrictions = .false.
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_k = 1 
         local_first_l = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_l = wf%n_o
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_ijkl, &
                                             local_first_i, local_last_i, &
                                             local_first_j, local_last_j, &
                                             local_first_k, local_last_k, &
                                             local_first_l, local_last_l)
!
   end subroutine get_oooo_ccs
!
!
   subroutine get_ooov_ccs(wf, g_ijka, first_i, last_i, first_j, last_j, &
                                         first_k, last_k, first_a, last_a)
!!
!!    Get ooov 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijka 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_k = first_k 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_k = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_ijka, &
                                             local_first_i, local_last_i, &
                                             local_first_j, local_last_j, &
                                             local_first_k, local_last_k, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a)
!
   end subroutine get_ooov_ccs
!
!
   subroutine get_oovo_ccs(wf, g_ijak, first_i, last_i, first_j, last_j, &
                                         first_a, last_a, first_k, last_k)
!!
!!    Get oovo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijak 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_k = first_k 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_k = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_ijak, &
                                             local_first_i, local_last_i, &
                                             local_first_j, local_last_j, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             local_first_k, local_last_k)
!
   end subroutine get_oovo_ccs
!
!
   subroutine get_ovoo_ccs(wf, g_iajk, first_i, last_i, first_a, last_a, &
                                         first_j, last_j, first_k, last_k)
!!
!!    Get ovoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iajk 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_k = first_k 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_k = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_iajk, &
                                             local_first_i, local_last_i, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             local_first_j, local_last_j, &
                                             local_first_k, local_last_k)
!
   end subroutine get_ovoo_ccs
!
!
   subroutine get_vooo_ccs(wf, g_aijk, first_a, last_a, first_i, last_i, &
                                         first_j, last_j, first_k, last_k)
!!
!!    Get vooo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aijk 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_k, last_k
      integer, optional, intent(in) :: first_a, last_a
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_k, local_last_k
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_k) .and. present(last_k) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_k = first_k 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_k = last_k
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_k = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_k = wf%n_o
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_aijk, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             local_first_i, local_last_i, &
                                             local_first_j, local_last_j, &
                                             local_first_k, local_last_k)
!
   end subroutine get_vooo_ccs
!
!
   subroutine get_vvoo_ccs(wf, g_abij, first_a, last_a, first_b, last_b, &
                                         first_i, last_i, first_j, last_j)
!!
!!    Get vvoo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abij 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_abij, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, wf%n_o + local_last_b, &
                                             local_first_i, local_last_i, &
                                             local_first_j, local_last_j)
!
   end subroutine get_vvoo_ccs
!
!
   subroutine get_vovo_ccs(wf, g_aibj, first_a, last_a, first_i, last_i, &
                                         first_b, last_b, first_j, last_j)
!!
!!    Get vovo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aibj 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      logical :: index_restrictions
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         index_restrictions = .true.
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         index_restrictions = .false.
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_aibj, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             local_first_i, local_last_i, &
                                             wf%n_o + local_first_b, wf%n_o + local_last_b, &
                                             local_first_j, local_last_j)
!
   end subroutine get_vovo_ccs
!
!
   subroutine get_voov_ccs(wf, g_aijb, first_a, last_a, first_i, last_i, &
                                         first_j, last_j, first_b, last_b)
!!
!!    Get voov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aijb 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_aijb, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             local_first_i, local_last_i, &
                                             local_first_j, local_last_j, &
                                             wf%n_o + local_first_b, wf%n_o + local_last_b)
!
   end subroutine get_voov_ccs
!
!
   subroutine get_ovvo_ccs(wf, g_iabj, first_i, last_i, first_a, last_a, &
                                         first_b, last_b, first_j, last_j)
!!
!!    Get ovvo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iabj
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_iabj, &
                                             local_first_i, local_last_i, &
                                             wf%n_o + local_first_a, wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, wf%n_o + local_last_b, &
                                             local_first_j, local_last_j)
!
   end subroutine get_ovvo_ccs
!
!
   subroutine get_oovv_ccs(wf, g_ijab, first_i, last_i, first_j, last_j, &
                                         first_a, last_a, first_b, last_b)
!!
!!    Get oovv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_ijab
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_j, last_j
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_j, local_last_j
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_j) .and. present(last_j) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_j = first_j 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_j = last_j
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_j = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_j = wf%n_o
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_ijab, &
                                             local_first_i, &
                                             local_last_i, &
                                             local_first_j, &
                                             local_last_j, &
                                             wf%n_o + local_first_a, &
                                             wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, &
                                             wf%n_o + local_last_b)
!
   end subroutine get_oovv_ccs
!
!
   subroutine get_vvvo_ccs(wf, g_abci, first_a, last_a, first_b, last_b, &
                                       first_c, last_c, first_i, last_i)
!!
!!    Get vvvo
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abci 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_c = first_c 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_c = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_abci, &
                                             wf%n_o + local_first_a, &
                                             wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, &
                                             wf%n_o + local_last_b, &
                                             wf%n_o + local_first_c, &
                                             wf%n_o + local_last_c,&
                                             local_first_i, &
                                             local_last_i)
!
   end subroutine get_vvvo_ccs
!
!
   subroutine get_vvov_ccs(wf, g_abic, first_a, last_a, first_b, last_b, &
                                         first_i, last_i, first_c, last_c)
!!
!!    Get vvov
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abic 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_c = first_c 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_c = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_abic, &
                                             wf%n_o + local_first_a, &
                                             wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, &
                                             wf%n_o + local_last_b, &
                                             local_first_i, &
                                             local_last_i, &
                                             wf%n_o + local_first_c, &
                                             wf%n_o + local_last_c)
!
   end subroutine get_vvov_ccs
!
!
   subroutine get_vovv_ccs(wf, g_aibc, first_a, last_a, first_i, last_i, &
                                         first_b, last_b, first_c, last_c)
!!
!!    Get vovv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_aibc 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_c = first_c 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_c = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_aibc, &
                                             wf%n_o + local_first_a, &
                                             wf%n_o + local_last_a, &
                                             local_first_i, &
                                             local_last_i, &
                                             wf%n_o + local_first_b, &
                                             wf%n_o + local_last_b, &
                                             wf%n_o + local_first_c, &
                                             wf%n_o + local_last_c)
!
   end subroutine get_vovv_ccs
!
!
   subroutine get_ovvv_ccs(wf, g_iabc, first_i, last_i, first_a, last_a, &
                                         first_b, last_b, first_c, last_c)
!!
!!    Get ovvv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_iabc 
!
      integer, optional, intent(in) :: first_i, last_i 
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_i, local_last_i 
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      if (present(first_i) .and. present(last_i) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         local_first_i = first_i 
         local_first_c = first_c 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_i = last_i 
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         local_first_i = 1 
         local_first_c = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_i = wf%n_o 
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_iabc, &
                                             local_first_i, &
                                             local_last_i, &
                                             wf%n_o + local_first_a, &
                                             wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, &
                                             wf%n_o + local_last_b, &
                                             wf%n_o + local_first_c, &
                                             wf%n_o + local_last_c)
!
   end subroutine get_ovvv_ccs
!
!
   subroutine get_vvvv_ccs(wf, g_abcd, first_a, last_a, first_b, last_b, &
                                         first_c, last_c, first_d, last_d)
!!
!!    Get vvvv
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018
!!
!!    The set of "get pqrs" routines will return the integral as t1-transformed,
!!    with the appropriate index restrictions if passed. If no index restrictions 
!!    are provided, the routines assume that the full integral should be returned.
!!
!!    Note that the MO integral tool controls how the integrals are constructed.
!!    The choice depends on logicals within the tool that knows whether t1-transformed 
!!    Cholesky vectors or the t1-transformed integrals themselves are on file. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(:,:,:,:), contiguous, intent(out) :: g_abcd
!
      integer, optional, intent(in) :: first_d, last_d 
      integer, optional, intent(in) :: first_c, last_c
      integer, optional, intent(in) :: first_a, last_a
      integer, optional, intent(in) :: first_b, last_b
!
      integer :: local_first_d, local_last_d 
      integer :: local_first_c, local_last_c
      integer :: local_first_b, local_last_b
      integer :: local_first_a, local_last_a
!
      logical :: index_restrictions
!
      if (present(first_d) .and. present(last_d) .and. &
          present(first_c) .and. present(last_c) .and. &
          present(first_b) .and. present(last_b) .and. &
          present(first_a) .and. present(last_a)) then
!
         index_restrictions = .true.
!
         local_first_d = first_d 
         local_first_c = first_c 
         local_first_b = first_b 
         local_first_a = first_a 
!
         local_last_d = last_d 
         local_last_c = last_c
         local_last_b = last_b
         local_last_a = last_a
!
      else
!
         index_restrictions = .false.
!
         local_first_d = 1 
         local_first_c = 1 
         local_first_b = 1 
         local_first_a = 1 
!
         local_last_d = wf%n_v 
         local_last_c = wf%n_v
         local_last_b = wf%n_v
         local_last_a = wf%n_v
!
      endif
!
      call wf%integrals%construct_g_pqrs_t1(g_abcd, &
                                             wf%n_o + local_first_a, &
                                             wf%n_o + local_last_a, &
                                             wf%n_o + local_first_b, &
                                             wf%n_o + local_last_b, &
                                             wf%n_o + local_first_c, &
                                             wf%n_o + local_last_c, &
                                             wf%n_o + local_first_d, &
                                             wf%n_o + local_last_d)
!
   end subroutine get_vvvv_ccs
!
!
   subroutine jacobian_transform_trial_vector_ccs(wf, c_i)
!!
!!    Jacobi transform trial vector 
!!    Written by Sarai D. Folkestad, Sep 2018
!!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
      call wf%jacobian_ccs_transformation(c_i)
!
   end subroutine jacobian_transform_trial_vector_ccs
!
!
   subroutine jacobian_transpose_transform_trial_vector_ccs(wf, c_i)
!!
!!    Jacobi transpose transform trial vector 
!!    Written by Sarai D. Folkestad, Sep 2018
!!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes) :: c_i
!
      call wf%jacobian_transpose_ccs_transformation(c_i)
!
   end subroutine jacobian_transpose_transform_trial_vector_ccs
!
!
   subroutine construct_excited_state_equation_ccs(wf, X, R, w)
!!
!!    Construct excited state equation 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Constructs R = AX - wX, where w = X^T A X and norm(X) = sqrt(X^T X) = 1
!!
!!    Note I: we assume that X is normalized. If it is not,
!!    please normalize before calling the routine. 
!!
!!    Note II: this routine constructs the excited state equation
!!    for standard CC models and the effective (!) excited state 
!!    equation in perturbative models. In the CC2 routine, for 
!!    instance, X and R will be n_o*n_v vectors and A(w) will 
!!    depend on the excitation energy w. See, e.g., Weigend and 
!!    Hättig's RI-CC2 paper for more on this topic. This means 
!!    that w should be the previous w-value when entering the 
!!    routine (so that A(w)X may be constructed approximately)
!!    in perturbative models.
!!
!!    Note III: the routine is used by the DIIS excited state solver. 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_es_amplitudes), intent(in)    :: X 
      real(dp), dimension(wf%n_es_amplitudes), intent(inout) :: R
!
      real(dp), intent(inout) :: w 
!
      real(dp), dimension(:), allocatable :: X_copy
!
      real(dp) :: ddot  
!
      call mem%alloc(X_copy, wf%n_es_amplitudes)
      call dcopy(wf%n_es_amplitudes, X, 1, X_copy, 1)
!
      call wf%jacobian_transform_trial_vector(X_copy) ! X_copy <- AX 
!
      w = ddot(wf%n_es_amplitudes, X, 1, X_copy, 1)
!
      call dcopy(wf%n_es_amplitudes, X_copy, 1, R, 1)
      call daxpy(wf%n_es_amplitudes, -w, X, 1, R, 1)
!
      call mem%dealloc(X_copy, wf%n_es_amplitudes)
!
   end subroutine construct_excited_state_equation_ccs
!
!
   subroutine jacobian_ccs_transformation_ccs(wf, c_ai)
!!
!!    Jacobian CCS transformation
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Directs the transformation by the CCSD Jacobi matrix,
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | nu >. 
!!
!!    In particular,
!!
!!       rho_mu = (A c)_mu = sum_ck A_mu,ck c_ck.
!! 
!!    On exit, c is overwritten by rho. 
!!
      implicit none
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: c_ai       
!
      real(dp), dimension(:,:), allocatable :: rho_ai
!
!     Allocate the transformed vector & add the terms to it 
!
      call mem%alloc(rho_ai, wf%n_v, wf%n_o)
      rho_ai = zero
!
      call wf%jacobian_ccs_a1(rho_ai, c_ai)
      call wf%jacobian_ccs_b1(rho_ai, c_ai)
!
!     Then overwrite the c vector with the transformed vector 
!
      call dcopy((wf%n_o)*(wf%n_v), rho_ai, 1, c_ai, 1)
      call mem%dealloc(rho_ai, wf%n_v, wf%n_o)
!
   end subroutine jacobian_ccs_transformation_ccs
!
!
   subroutine jacobian_transpose_ccs_transformation_ccs(wf, b_ai)
!!
!!    Jacobian transpose transformation (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the transpose Jacobian transformation, i.e., the transformation 
!!    by the transpose of the Jacobian matrix
!!
!!       A_mu,nu = < mu | exp(-T) [H, tau_nu] exp(T) | R >.
!!
!!    In particular,
!!
!!       sigma_mu = (b^T A)_mu = sum_ck b_ck A_ck,mu.
!! 
!!    On exit, b is overwritten by sigma. 
!!
      implicit none 
!
      class(ccs) :: wf 
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: b_ai 
!
      real(dp), dimension(:,:), allocatable :: sigma_ai
!
!     Allocate the transformed vector & add the terms to it
!
      call mem%alloc(sigma_ai, wf%n_v, wf%n_o)
      sigma_ai = zero 
!
      call wf%jacobian_transpose_ccs_a1(sigma_ai, b_ai)
      call wf%jacobian_transpose_ccs_b1(sigma_ai, b_ai)
!
!     Then overwrite the b vector with the transformed vector 
!
      call dcopy((wf%n_o)*(wf%n_v), sigma_ai, 1, b_ai, 1)
      call mem%dealloc(sigma_ai, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccs_transformation_ccs
!
!
   subroutine jacobian_ccs_a1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS A1
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_b F_ab c_bi - sum_j F_ji c_aj,
!!
!!    and adds it to the rho vector.
!!
      implicit none
!
      class(ccs) :: wf
!
      real(dp), dimension(wf%n_v*wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1
!
!     sum_b F_a_b c_b_i
!
      call dgemm('N', 'N',     &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_v,      &
                  one,         &
                  wf%fock_ab,  &
                  wf%n_v,      &
                  c1,          &
                  wf%n_v,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
!     - sum_j c_a_j F_j_i
!
      call dgemm('N','N',      &
                  wf%n_v,      &
                  wf%n_o,      &
                  wf%n_o,      &
                  -one,        &
                  c1,          &
                  wf%n_v,      &
                  wf%fock_ij,  &
                  wf%n_o,      &
                  one,         &
                  rho1,        &
                  wf%n_v)
!
   end subroutine jacobian_ccs_a1_ccs
!
!
   subroutine jacobian_ccs_b1_ccs(wf, rho1, c1)
!!
!!    Jacobian CCS B1 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, May 2017
!!
!!    Calculates the B1 term,
!!
!!       sum_bj L_aijb c_bj = sum_bj (2 g_aijb - g_abji) c_bj,
!!
!!    and adds it to the rho1 vector.
!!
      implicit none
!
      class(ccs) :: wf
!   
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: c1
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: rho1      
!
      real(dp), dimension(:,:,:,:), allocatable :: g_abji
      real(dp), dimension(:,:,:,:), allocatable :: L_aijb
!
      real(dp), dimension(:,:), allocatable :: c_jb
!
      type(batching_index) :: batch_b
!
      integer :: req0, req1, j, b, b_red, current_b_batch
!
      call batch_b%init(wf%n_v) 
!
      req0 = wf%n_o*wf%n_v*wf%integrals%n_J ! L_ai^J 
      req1 = 2*wf%n_v*wf%n_o**2 + & ! L_aijb and g_abji 
               wf%n_v*wf%integrals%n_J + & ! L_ab^J 
               wf%n_o ! c_jb 
!
      call mem%batch_setup(batch_b, req0, req1)
!
      do current_b_batch = 1, batch_b%num_batches
!
         call batch_b%determine_limits(current_b_batch)
!
!        Construct L_aijb = 2 g_aijb - g_abji
!
         call mem%alloc(L_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
!
         call wf%get_voov(L_aijb,     &
                           1, wf%n_v,  &
                           1, wf%n_o,  &
                           1, wf%n_o,  &
                           batch_b%first, batch_b%last)   
!
         call dscal(((wf%n_o)**2)*(wf%n_v)*(batch_b%length), two, L_aijb, 1)
!
!        Construct L_aijb = 2 g_aijb - g_abji
!
         call mem%alloc(g_abji, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_abji,                        &
                           1, wf%n_v,                    &
                           batch_b%first, batch_b%last,  &
                           1, wf%n_o,                    &
                           1, wf%n_o)   
!
         call add_1432_to_1234(-one, g_abji, L_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
!
         call mem%dealloc(g_abji, wf%n_v, batch_b%length, wf%n_o, wf%n_o)
!
!        Reorder c1 to do multiply with L_aijb
!
         call mem%alloc(c_jb, (wf%n_o), batch_b%length)
!
         do b = batch_b%first, batch_b%last
            do j = 1, wf%n_o
!
               b_red = b - batch_b%first + 1
!  
               c_jb(j, b_red) = c1(b, j)
!
            enddo
         enddo
    
!
         call dgemm('N', 'N',                   &
                     (wf%n_v)*(wf%n_o),         &
                     1,                         &
                     (wf%n_o)*(batch_b%length), &
                     one,                       &
                     L_aijb,                   &
                     (wf%n_v)*(wf%n_o),         &
                     c_jb,                      &
                     (wf%n_o)*batch_b%length,   &
                     one,                       &
                     rho1,                      &
                     (wf%n_v)*(wf%n_o))      
!
         call mem%dealloc(L_aijb, wf%n_v, wf%n_o, wf%n_o, batch_b%length)
         call mem%dealloc(c_jb, (wf%n_o), (batch_b%length))
!
      enddo
!
   end subroutine jacobian_ccs_b1_ccs
!
!
   subroutine jacobian_transpose_ccs_a1_ccs(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose A1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the A1 term,
!!
!!       sum_c b_ci F_ca - sum_k b_ak F_ik,
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai 
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
!     Add sum_c F_ca b_ci = sum_c F_ac^T b_ci     
!
      call dgemm('T','N',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_v,     &
                  one,        &
                  wf%fock_ab, &
                  wf%n_v,     &
                  b_ai,       &
                  wf%n_v,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
!     Add - sum_k b_ak F_ik = - sum_k b_ak F_ki^T 
!
      call dgemm('N','T',     &
                  wf%n_v,     &
                  wf%n_o,     &
                  wf%n_o,     &
                  -one,       &
                  b_ai,       &
                  wf%n_v,     &
                  wf%fock_ij, &
                  wf%n_o,     &
                  one,        &
                  sigma_ai,   &
                  wf%n_v)
!
   end subroutine jacobian_transpose_ccs_a1_ccs
! 
! 
   subroutine jacobian_transpose_ccs_b1_ccs(wf, sigma_ai, b_ai)
!!
!!    Jacobian transpose B1 (CCS)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Calculates the B1 term,
!!
!!       sum_ck L_ckia b_ck
!!
!!    and adds it to the sigma-vector (b^T -> sigma^T = b^T A).
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_v, wf%n_o), intent(inout) :: sigma_ai 
      real(dp), dimension(wf%n_v, wf%n_o), intent(in)    :: b_ai
!
      real(dp), dimension(:,:,:,:), allocatable :: g_ckia ! g_ckia 
      real(dp), dimension(:,:,:,:), allocatable :: g_caik ! g_caik 
!
      real(dp), dimension(:,:,:,:), allocatable :: L_aick ! L_ckia = 2 * g_ckia - g_caik
!
      integer :: k, c, i, a
!
      integer              :: req0, req1, current_a_batch
      type(batching_index) :: batch_a 
!
!     :: Construct L_aick = L_ckia
!
      call mem%alloc(g_ckia, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
      call wf%get_voov(g_ckia)
!
      call mem%alloc(L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      L_aick = zero
!
      req0 = wf%integrals%n_J*wf%n_o**2 ! L_ik^J
      req1 = wf%n_v*wf%n_o**2 &         ! g_caik
             + wf%integrals%n_J*wf%n_v  ! L_ab^J
!
      call batch_a%init(wf%n_v)
!
      call mem%batch_setup(batch_a, req0, req1)
!
      do current_a_batch = 1, batch_a%num_batches 
!
!        Set part of L_aick = L_ckia = 2 * g_ckia - g_caik for current a batch 
!
         call batch_a%determine_limits(current_a_batch)
!
         call mem%alloc(g_caik, wf%n_v, batch_a%length, wf%n_o, wf%n_o)
!
         call wf%get_vvoo(g_caik,         &
                           1,             &
                           wf%n_v,        &
                           batch_a%first, &
                           batch_a%last,  &
                           1,             &
                           wf%n_o,        &
                           1,             &
                           wf%n_o)
!
!$omp parallel do private(k,c,i,a)
         do k = 1, wf%n_o
            do c = 1, wf%n_v
               do i = 1, wf%n_o
                  do a = 1, batch_a%length
!
                     L_aick(a + batch_a%first - 1,i,c,k) = two*g_ckia(c,k,i,a + batch_a%first - 1) - g_caik(c,a,i,k)
!
                  enddo
               enddo
            enddo
         enddo
!$omp end parallel do
!
         call mem%dealloc(g_caik, wf%n_v, batch_a%length, wf%n_o, wf%n_o)
!
      enddo ! End of batches over a
!
      call mem%dealloc(g_ckia, wf%n_v, wf%n_o, wf%n_o, wf%n_v)
!
!     :: Add sum_ck L_ckia b_ck = sum_ck L_aick b_ck to sigma 
!
      call dgemm('N','N',            &
                  (wf%n_v)*(wf%n_o), &
                  1,                 &
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  L_aick,            &
                  (wf%n_v)*(wf%n_o), &
                  b_ai,              & ! "b_ai"
                  (wf%n_v)*(wf%n_o), &
                  one,               &
                  sigma_ai,          & ! "sigma_ai"
                  (wf%n_v)*(wf%n_o))
!
      call mem%dealloc(L_aick, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
   end subroutine jacobian_transpose_ccs_b1_ccs
!
!
   subroutine construct_eta_ccs(wf, eta)
!!
!!    Construct eta 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
      implicit none
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: eta 
!
      integer :: i, a, ai
!
!$omp parallel do private(a, i, ai)
      do i = 1, wf%n_o 
         do a = 1, wf%n_v 
!
            ai = (wf%n_v)*(i - 1) + a
            eta(ai) = two*(wf%fock_ia(i, a))
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine construct_eta_ccs
!
!
   subroutine construct_multiplier_equation_ccs(wf, equation)
!!
!!    Construct multiplier equation 
!!    Written by Eirik F. Kjønstad, Oct 2018 
!!
!!    Constructs 
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes), intent(inout) :: equation 
!
      real(dp), dimension(:), allocatable :: eta 
!
!     Copy the multipliers, eq. = t-bar 
!
      call dcopy(wf%n_t1, wf%t1bar, 1, equation, 1)
!
!     Transform the multipliers by A^T, eq. = t-bar^T A 
!
      call wf%jacobian_transpose_ccs_transformation(equation)
!
!     Add eta, eq. = t-bar^T A + eta 
!
      call mem%alloc(eta, wf%n_t1)
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_t1, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_t1)
!
   end subroutine construct_multiplier_equation_ccs
!
!
   subroutine add_bath_orbitals_ccs(wf)
!!
!!    Add bath orbitals,
!!    Written by Sarai D. Folkestad, Oct. 2018 
!!
      implicit none
!
      class(ccs) :: wf
!
      integer :: p, q, pq, ao, removed_orbitals
!
      type(file) :: h_pq_file
!
      real(dp), dimension(:,:), allocatable :: orbital_coeff_copy, L_J, h_pq
!
!     Read number and type of bath orbitals (for now only 1 and for ionization)
!
      wf%n_bath = 1
!
!     Add atom X to system (if necessary), with s-type orbitals 
!
      call mem%alloc(orbital_coeff_copy, wf%n_ao, wf%n_mo)
!
      orbital_coeff_copy = wf%orbital_coefficients
!
      call mem%dealloc(wf%orbital_coefficients, wf%n_ao, wf%n_mo)
!
!     Update coefficient matrix
!
      call mem%alloc(wf%orbital_coefficients, wf%n_ao + wf%n_bath, wf%n_mo + wf%n_bath)
      wf%orbital_coefficients = zero
!
      wf%orbital_coefficients(1:wf%n_ao, 1:wf%n_mo) = orbital_coeff_copy(:,:)
!
      call mem%dealloc(orbital_coeff_copy, wf%n_ao, wf%n_mo)
!
      removed_orbitals = wf%n_ao - wf%n_mo ! Due to linear dependancy
!
!     Bath orbitals do not mix with other orbitals
!
      do ao = wf%n_ao + 1, wf%n_ao + wf%n_bath
!
         wf%orbital_coefficients(ao, ao - removed_orbitals) = one
!
      enddo
!
!     Update n_ao, n_mo, n_v
!
      wf%n_ao = wf%n_ao + wf%n_bath
      wf%n_mo = wf%n_mo + wf%n_bath
      wf%n_v  = wf%n_v + wf%n_bath
!
!     Update h_pq matrix
!
      call h_pq_file%init('h_pq', 'sequential', 'unformatted')
      call disk%open_file(h_pq_file, 'readwrite')
      rewind(h_pq_file%unit)
!
      call mem%alloc(h_pq, wf%n_mo, wf%n_mo)
      h_pq = zero
!
      read(h_pq_file%unit) h_pq(1:wf%n_mo - wf%n_bath, 1:wf%n_mo - wf%n_bath)
      rewind(h_pq_file%unit)
      write(h_pq_file%unit) h_pq
!
      call disk%close_file(h_pq_file)
!
!     Update cholesky vectors
!
      call disk%open_file(wf%integrals%cholesky_mo, 'write')
!
      call mem%alloc(L_J, 1, wf%integrals%n_J)
      L_J = zero
!
      do p = wf%n_mo - wf%n_bath, wf%n_mo
         do q = 1, p 
!
            pq = p*(p-3)/2 + p + q
!
            write(wf%integrals%cholesky_mo%unit, rec=pq) L_J
!
         enddo
      enddo
!
      call mem%dealloc(L_J, 1, wf%integrals%n_J)
!
      call disk%close_file(wf%integrals%cholesky_mo)
!
   end subroutine add_bath_orbitals_ccs
!
!
   subroutine get_cvs_projector_ccs(wf, projector, n_cores, core_MOs)
!!
!!    Get CVS projector
!!    Written by Sarai D. Folkestad, Oct 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer, intent(in) :: n_cores
!
      integer, dimension(n_cores), intent(in) :: core_MOs
!
      integer :: core, i, a, ai
!
      projector = zero
!
      do core = 1, n_cores
!
        i = core_MOs(core)
!
        do a = 1, wf%n_v
!
           ai = wf%n_v*(i - 1) + a
           projector(ai) = one
!
        enddo
     enddo
!
   end subroutine get_cvs_projector_ccs
!
!
   subroutine get_ip_projector_ccs(wf, projector)
!!
!!    Get ip projector
!!    Written by Sarai D. Folekstad, Oct 2018
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_es_amplitudes), intent(out) :: projector
!
      integer :: i, a, ai
!
      projector = zero
!
      a = wf%n_v ! Last virtual is bath orbital
!
      do i = 1, wf%n_o
!
         ai = wf%n_v*(i - 1) + a
         projector(ai) = one
!
     enddo
!
   end subroutine get_ip_projector_ccs
!
!
   subroutine print_dominant_amplitudes_ccs(wf)
!!
!!    Print dominant amplitudes 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      call wf%print_dominant_x1(wf%t1,'t')
!
   end subroutine print_dominant_amplitudes_ccs
!
!
   subroutine print_dominant_x_amplitudes_ccs(wf, x, tag)
!!
!!    Print dominant amplitudes
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(wf%n_gs_amplitudes) :: x 
!
      character(len=1) :: tag 
!
      call wf%print_dominant_x1(x(1:wf%n_t1),tag)
!
   end subroutine print_dominant_x_amplitudes_ccs
!
!
   subroutine print_dominant_x1_ccs(wf, x1, tag)
!!
!!    Print dominant x1    
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Prints the 20 most dominant single amplitudes,
!!    or sorts them if there are fewer than twenty of them.
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_t1), intent(in) :: x1 
      character(len=1), intent(in)                :: tag 
!
      real(dp), dimension(:), allocatable :: abs_x1
!
      integer, dimension(:), allocatable :: dominant_indices
      real(dp), dimension(:), allocatable     :: dominant_values
!
      integer :: n_elements, elm, i, a 
!
!     Sort according to largest contributions
!
      call mem%alloc(abs_x1, wf%n_t1)
      abs_x1 = abs(x1)
!
      n_elements = 20
      if (n_elements .gt. wf%n_t1) n_elements = wf%n_t1 
!
      call mem%alloc(dominant_indices, n_elements)
      call mem%alloc(dominant_values, n_elements)
!
      dominant_indices = 0
      dominant_values  = zero 
      call get_n_highest(n_elements, wf%n_t1, abs_x1, dominant_values, dominant_indices)
!
!     Print largest contributions
!
      write(output%unit, '(/t6,a)') 'Largest single amplitudes:'
      write(output%unit, '(t6,a)')  '-----------------------------------------'
      write(output%unit, '(t6,a)')  '  a         i         ' // tag // '(a,i)             '
      write(output%unit, '(t6,a)')  '-----------------------------------------'
!
      do elm = 1, n_elements
!
         call invert_compound_index(dominant_indices(elm), a, i, wf%n_v, wf%n_o)
!
         write(output%unit, '(t6,i3,7x,i3,5x,f19.12)') a, i, x1(dominant_indices(elm))
!
      enddo
!
      write(output%unit, '(t6,a)')  '-----------------------------------------'
!
      call mem%dealloc(dominant_indices, n_elements)
      call mem%dealloc(dominant_values, n_elements)
      call mem%dealloc(abs_x1, wf%n_t1)
!
   end subroutine print_dominant_x1_ccs
!
!
   real(dp) function get_t1_diagnostic_ccs(wf)
!!
!!    Get t1 diagnostic 
!!    Written by Eirik F. Kjønstad 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf 
!
      get_t1_diagnostic_ccs = get_l2_norm(wf%t1, wf%n_t1)
      get_t1_diagnostic_ccs = get_t1_diagnostic_ccs/sqrt(real(wf%system%n_electrons,kind=dp))
!
   end function get_t1_diagnostic_ccs
!
!
   subroutine prepare_for_jacobian_ccs(wf)
!!
!!    Prepare for jacobian
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Jan 2019
!!
      implicit none 
!
      class(ccs), intent(inout) :: wf 
!
!     For now, do nothing.
!
      write(output%unit,'(/t3,a,a,a)') 'No Jacobian preparations for ', &
                                       trim(wf%name_), ' wavefunction.'
!
   end subroutine prepare_for_jacobian_ccs
!
!
   subroutine set_cvs_start_indices_ccs(wf, n_cores, core_MOs, n_start_indices, start_indices)
!!
!!    Set CVS start indices
!!    Written by Sarai D. Folkestad 
!!
      implicit none 
!
      class(ccs), intent(in) :: wf
!
      integer :: n_cores, n_start_indices
!
      integer, dimension(n_cores)          :: core_MOs
      integer, dimension(n_start_indices)  :: start_indices
!
!     Local variables
!
      integer :: a, i, current_root, core
!
      logical :: all_selected
!
!     Calculate start indices using Koopman's
!
      all_selected = .false.
      a =  0 
      current_root = 0
!
      do while (.not. all_selected)
!
         a = a + 1
!
         do core = 1, n_cores
!
            i = core_MOs(core)
!
            current_root = current_root + 1
            start_indices(current_root) = wf%n_v*( i - 1) + a
!
            if (current_root .eq. n_start_indices) then
!
               all_selected = .true.
               exit
!
            endif
!
         enddo
!
      enddo
!
   end subroutine set_cvs_start_indices_ccs
!
!
end module ccs_class
