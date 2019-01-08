module ccsd_class
!
!!
!!    Coupled cluster singles and doubles (ccsd) class module
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad,
!!    Andreas Skeidsvoll and Alice Balbi, 2018
!!
!
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: ccsd
!
      real(dp), dimension(:,:), allocatable :: t2   
      real(dp), dimension(:,:), allocatable :: t2bar   
!
      integer(i15) :: n_t2  
!
   contains
!
!     Preparation and cleanup routines 
!
      procedure :: prepare                                     => prepare_ccsd
      procedure :: cleanup                                     => cleanup_ccsd
!
!     Routines related to the amplitudes 
!
      procedure :: initialize_amplitudes                       => initialize_amplitudes_ccsd 
      procedure :: initialize_t2                               => initialize_t2_ccsd
      procedure :: destruct_t2                                 => destruct_t2_ccsd
      procedure :: set_initial_amplitudes_guess                => set_initial_amplitudes_guess_ccsd
      procedure :: set_t2_to_mp2_guess                         => set_t2_to_mp2_guess_ccsd
      procedure :: set_amplitudes                              => set_amplitudes_ccsd 
      procedure :: get_amplitudes                              => get_amplitudes_ccsd 
      procedure :: read_amplitudes                             => read_amplitudes_ccsd
      procedure :: save_amplitudes                             => save_amplitudes_ccsd
      procedure :: save_t2                                     => save_t2_ccsd
      procedure :: read_t2                                     => read_t2_ccsd
      procedure :: print_dominant_x2                           => print_dominant_x2_ccsd
      procedure :: print_dominant_amplitudes                   => print_dominant_amplitudes_ccsd
      procedure :: print_dominant_x_amplitudes                 => print_dominant_x_amplitudes_ccsd
!
!     Routines related to omega
!
      procedure :: construct_omega                             => construct_omega_ccsd
!
      procedure :: omega_ccsd_a1                               => omega_ccsd_a1_ccsd
      procedure :: omega_ccsd_b1                               => omega_ccsd_b1_ccsd
      procedure :: omega_ccsd_c1                               => omega_ccsd_c1_ccsd
!
      procedure :: omega_ccsd_a2                               => omega_ccsd_a2_ccsd
      procedure :: omega_ccsd_b2                               => omega_ccsd_b2_ccsd
      procedure :: omega_ccsd_c2                               => omega_ccsd_c2_ccsd
      procedure :: omega_ccsd_d2                               => omega_ccsd_d2_ccsd
      procedure :: omega_ccsd_e2                               => omega_ccsd_e2_ccsd
!
!     Routines related to Jacobian transformation
!
      procedure :: jacobian_transform_trial_vector             => jacobian_transform_trial_vector_ccsd
      procedure :: jacobian_ccsd_transformation                => jacobian_ccsd_transformation_ccsd
!
      procedure :: jacobian_ccsd_a1                            => jacobian_ccsd_a1_ccsd
      procedure :: jacobian_ccsd_b1                            => jacobian_ccsd_b1_ccsd
      procedure :: jacobian_ccsd_c1                            => jacobian_ccsd_c1_ccsd
      procedure :: jacobian_ccsd_d1                            => jacobian_ccsd_d1_ccsd
!
      procedure :: jacobian_ccsd_a2                            => jacobian_ccsd_a2_ccsd
      procedure :: jacobian_ccsd_b2                            => jacobian_ccsd_b2_ccsd
      procedure :: jacobian_ccsd_c2                            => jacobian_ccsd_c2_ccsd
      procedure :: jacobian_ccsd_d2                            => jacobian_ccsd_d2_ccsd
      procedure :: jacobian_ccsd_e2                            => jacobian_ccsd_e2_ccsd
      procedure :: jacobian_ccsd_f2                            => jacobian_ccsd_f2_ccsd
      procedure :: jacobian_ccsd_g2                            => jacobian_ccsd_g2_ccsd
      procedure :: jacobian_ccsd_h2                            => jacobian_ccsd_h2_ccsd
      procedure :: jacobian_ccsd_i2                            => jacobian_ccsd_i2_ccsd
      procedure :: jacobian_ccsd_j2                            => jacobian_ccsd_j2_ccsd
      procedure :: jacobian_ccsd_k2                            => jacobian_ccsd_k2_ccsd
!
!     Routines related to Jacobian transpose transformation
!
      procedure :: jacobian_transpose_transform_trial_vector   => jacobian_transpose_transform_trial_vector_ccsd
      procedure :: jacobian_transpose_ccsd_transformation      => jacobian_transpose_ccsd_transformation_ccsd
!
      procedure :: jacobian_transpose_ccsd_a1                  => jacobian_transpose_ccsd_a1_ccsd
      procedure :: jacobian_transpose_ccsd_b1                  => jacobian_transpose_ccsd_b1_ccsd
      procedure :: jacobian_transpose_ccsd_c1                  => jacobian_transpose_ccsd_c1_ccsd
      procedure :: jacobian_transpose_ccsd_d1                  => jacobian_transpose_ccsd_d1_ccsd
      procedure :: jacobian_transpose_ccsd_e1                  => jacobian_transpose_ccsd_e1_ccsd
      procedure :: jacobian_transpose_ccsd_f1                  => jacobian_transpose_ccsd_f1_ccsd
      procedure :: jacobian_transpose_ccsd_g1                  => jacobian_transpose_ccsd_g1_ccsd
!
      procedure :: jacobian_transpose_ccsd_a2                  => jacobian_transpose_ccsd_a2_ccsd
      procedure :: jacobian_transpose_ccsd_b2                  => jacobian_transpose_ccsd_b2_ccsd
      procedure :: jacobian_transpose_ccsd_c2                  => jacobian_transpose_ccsd_c2_ccsd
      procedure :: jacobian_transpose_ccsd_d2                  => jacobian_transpose_ccsd_d2_ccsd
      procedure :: jacobian_transpose_ccsd_e2                  => jacobian_transpose_ccsd_e2_ccsd
      procedure :: jacobian_transpose_ccsd_f2                  => jacobian_transpose_ccsd_f2_ccsd
      procedure :: jacobian_transpose_ccsd_g2                  => jacobian_transpose_ccsd_g2_ccsd
      procedure :: jacobian_transpose_ccsd_h2                  => jacobian_transpose_ccsd_h2_ccsd
      procedure :: jacobian_transpose_ccsd_i2                  => jacobian_transpose_ccsd_i2_ccsd
!
      procedure :: get_orbital_differences                     => get_orbital_differences_ccsd
      procedure :: calculate_energy                            => calculate_energy_ccsd
!
      procedure :: construct_eta                               => construct_eta_ccsd
!
      procedure :: initialize_t2bar                            => initialize_t2bar_ccsd
      procedure :: get_multipliers                             => get_multipliers_ccsd
      procedure :: set_multipliers                             => set_multipliers_ccsd
      procedure :: initialize_multipliers                      => initialize_multipliers_ccsd
      procedure :: construct_multiplier_equation               => construct_multiplier_equation_ccsd
      procedure :: save_multipliers                            => save_multipliers_ccsd
      procedure :: read_multipliers                            => read_multipliers_ccsd
      procedure :: destruct_multipliers                        => destruct_multipliers_ccsd
      procedure :: destruct_t2bar                              => destruct_t2bar_ccsd
      procedure :: read_t2bar                                  => read_t2bar_ccsd
      procedure :: save_t2bar                                  => save_t2bar_ccsd
!
   end type ccsd
!
!
   interface
!
      include "../submodules/ccsd/omega_ccsd_interface.F90"
      include "../submodules/ccsd/jacobian_ccsd_interface.F90"
      include "../submodules/ccsd/jacobian_transpose_ccsd_interface.F90"
!
   end interface 
!
!
contains
!
!
   subroutine prepare_ccsd(wf, ref_wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      class(hf) :: ref_wf
!
      integer(i15) :: p
!
      wf%name = 'ccsd'
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
      wf%n_t1 = (wf%n_o)*(wf%n_v) 
      wf%n_t2 = (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
!
      wf%n_amplitudes = wf%n_t1 + wf%n_t2 
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
         wf%fock_diagonal(p, 1) = ref_wf%mo_fock(p, p)
!
      enddo
!
      call ref_wf%mo_transform_and_save_h()
!
      call wf%initialize_orbital_coefficients()
      wf%orbital_coefficients = ref_wf%orbital_coefficients
!
   end subroutine prepare_ccsd
!
!
   subroutine cleanup_ccsd(wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      write(output%unit, '(/t3,a,a,a)') '- Cleaning up ', trim(wf%name), ' wavefunction'
!
   end subroutine cleanup_ccsd
!
!
   subroutine initialize_amplitudes_ccsd(wf)
!!
!!    Initialize amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Allocates the amplitudes. This routine must be overwritten in 
!!    descendants which have more amplitudes. 
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      call wf%initialize_t1()
      call wf%initialize_t2()
!
   end subroutine initialize_amplitudes_ccsd
!
!
   subroutine initialize_t2_ccsd(wf)
!!
!!    Initialize t2 amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      if (.not. allocated(wf%t2)) call mem%alloc(wf%t2, wf%n_t2, 1)
!
   end subroutine initialize_t2_ccsd
!
!
   subroutine destruct_t2_ccsd(wf)
!!
!!    Destruct t2 amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      if (allocated(wf%t2)) call mem%dealloc(wf%t2, wf%n_t2, 1)
!
   end subroutine destruct_t2_ccsd
!
!
   subroutine set_amplitudes_ccsd(wf, amplitudes)
!!
!!    Set amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(ccsd) :: wf  
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(in) :: amplitudes
!
      call dcopy(wf%n_t1, amplitudes, 1, wf%t1, 1)
      call dcopy(wf%n_t2, amplitudes(wf%n_t1 + 1, 1), 1, wf%t2, 1)
!
   end subroutine set_amplitudes_ccsd
!
!
   subroutine get_amplitudes_ccsd(wf, amplitudes)
!!
!!    Get amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf  
!
      real(dp), dimension(wf%n_amplitudes, 1) :: amplitudes
!
      call dcopy(wf%n_t1, wf%t1, 1, amplitudes, 1)
      call dcopy(wf%n_t2, wf%t2, 1,  amplitudes(wf%n_t1 + 1, 1), 1)
!
   end subroutine get_amplitudes_ccsd
!
!
   subroutine set_initial_amplitudes_guess_ccsd(wf)
!!
!!    Set initial amplitudes guess 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      wf%t1 = zero 
!
      call wf%set_t2_to_mp2_guess()
!
   end subroutine set_initial_amplitudes_guess_ccsd
!
!
   subroutine set_t2_to_mp2_guess_ccsd(wf)
!!
!!    Set t2 amplitudes guess 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    t_aibj = - g_aibj/ε_aibj
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
!
      integer(i15) :: a, b, i, j, ai, bj, aibj
!
      call mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call wf%get_vovo(g_ai_bj)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o 
!
            ai = wf%n_v*(i-1) + a
!
            do j = 1, wf%n_o
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j-1) + b              
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     wf%t2(aibj, 1) = g_ai_bj(ai, bj)/(wf%fock_diagonal(i, 1) + &
                                                       wf%fock_diagonal(j, 1) - &
                                                       wf%fock_diagonal(a + wf%n_o, 1) - &
                                                       wf%fock_diagonal(b + wf%n_o, 1))
!
                  endif
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      call mem%dealloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine set_t2_to_mp2_guess_ccsd
!
!
   subroutine calculate_energy_ccsd(wf)
!!
!!     Calculate energy (CCSD)
!!     Written by Sarai D. Folkestad, Eirik F. Kjønstad, 
!!     Andreas Skeidsvoll, 2018
!!
!!     Calculates the CCSD energy. This is only equal to the actual
!!     energy when the ground state equations are solved, of course.
!!
      implicit none
!
      class(ccsd), intent(inout) :: wf
!
      real(dp), dimension(:,:), allocatable :: g_ia_jb ! g_iajb
!
      integer(i15) :: a = 0, i = 0, b = 0, j = 0, ai = 0
      integer(i15) :: bj = 0, aibj = 0, ia = 0, jb = 0, ib = 0, ja = 0
!
!     Get g_ia_jb = g_iajb
!
      call mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_ia_jb)
!
!     Set the initial value of the energy
!
      wf%energy = wf%hf_energy
!
!     Add the correlation energy E = E + sum_aibj (t_ij^ab + t_i^a t_j^b) L_iajb
!
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = index_two(a, i, wf%n_v)
            ia = index_two(i, a, wf%n_o)
!
            do j = 1, wf%n_o
!
               ja = wf%n_o*(a-1) + j
!
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j - 1) + b
                  jb = wf%n_o*(b - 1) + j
                  ib = wf%n_o*(b - 1) + i
!
                  aibj = (max(ai,bj)*(max(ai,bj)-3)/2) + ai + bj
!
!                 Add the correlation energy
!
                  wf%energy = wf%energy +                                     &
                                 (wf%t2(aibj,1) + (wf%t1(a,i))*(wf%t1(b,j)))* &
                                 (two*g_ia_jb(ia,jb) - g_ia_jb(ib,ja))
!
               enddo
            enddo
         enddo
      enddo
!
!     Deallocate g_ia_jb
!
      call mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine calculate_energy_ccsd
!
!
   subroutine get_orbital_differences_ccsd(wf, orbital_differences)
!!
!!    Get orbital differences 
!!    Written by Eirik F. Kjønstad, Sarai D. Folkestad
!!    and Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: orbital_differences
!
      integer(i15) :: a, i, ai, b, j, bj, aibj
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj, aibj) 
      do a = 1, wf%n_v
         do i = 1, wf%n_o
!
            ai = wf%n_v*(i - 1) + a
!
            orbital_differences(ai, 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1)
!
            do j = 1, wf%n_o 
               do b = 1, wf%n_v
!
                  bj = wf%n_v*(j-1) + b 
!
                  if (ai .ge. bj) then
!
                     aibj = (ai*(ai-3)/2) + ai + bj
!
                     orbital_differences(aibj + (wf%n_o)*(wf%n_v), 1) = wf%fock_diagonal(a + wf%n_o, 1) - wf%fock_diagonal(i, 1) &
                                                                      +  wf%fock_diagonal(b + wf%n_o, 1) - wf%fock_diagonal(j, 1)
!
                  endif
!
               enddo
            enddo  
!
         enddo
      enddo
!$omp end parallel do
!
   end subroutine get_orbital_differences_ccsd
!
!
   subroutine read_amplitudes_ccsd(wf)
!!
!!    Read amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccsd), intent(inout) :: wf
!
      call wf%read_t1()  
      call wf%read_t2()  
!
   end subroutine read_amplitudes_ccsd
!
!
   subroutine save_amplitudes_ccsd(wf)
!!
!!    Read amplitudes 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf
!
      call wf%save_t1()  
      call wf%save_t2()  
!
   end subroutine save_amplitudes_ccsd
!
!
   subroutine read_t2_ccsd(wf)
!!
!!    Read t2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none 
!
      class(ccsd), intent(inout) :: wf
!
      type(file) :: t2_file 
!
      call t2_file%init('t2', 'sequential', 'unformatted')
!
      call disk%open_file(t2_file, 'read', 'rewind')
!
      read(t2_file%unit) wf%t2
!
      call disk%close_file(t2_file)      
!
   end subroutine read_t2_ccsd
!
!
   subroutine save_t2_ccsd(wf)
!!
!!    Save t2
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, May 2017
!!
      implicit none
!
      class(ccsd), intent(in) :: wf 
!
      type(file) :: t2_file 
!
      call t2_file%init('t2', 'sequential', 'unformatted')
!
      call disk%open_file(t2_file, 'write', 'rewind')
!
      write(t2_file%unit) wf%t2
!
      call disk%close_file(t2_file)
!
   end subroutine save_t2_ccsd
!
subroutine construct_eta_ccsd(wf, eta)
!!
!!    Construct eta (CCSD)
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, June 2017
!!
!!    Note: the routine assumes that eta is initialized and that the Fock matrix
!!    has been constructed.
!!
      implicit none
!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: eta 
!
      real(dp), dimension(:,:), allocatable :: g_ia_jb
      real(dp), dimension(:,:), allocatable :: eta_ai_bj
!
      integer(i15) :: i = 0, a = 0, j = 0, b = 0, aibj = 0
      integer(i15) :: bj = 0, ai = 0
!
      eta = zero
!
      do i = 1, wf%n_o
         do a = 1, wf%n_v
!
            ai = wf%n_v*(i - 1) + a
            eta(ai, 1) = two*(wf%fock_ia(i, a)) ! eta_ai = 2 F_ia
!
         enddo
      enddo
!
!     eta_ai_bj = 2* L_iajb = 4 * g_ia_jb(ia,jb) - 2 * g_ia_jb(ib,ja)
!
      call mem%alloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call wf%get_ovov(g_ia_jb)
!
      call mem%alloc(eta_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      eta_ai_bj = zero
!
      call add_2143_to_1234(four, g_ia_jb, eta_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call add_2341_to_1234(-two, g_ia_jb, eta_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      call mem%dealloc(g_ia_jb, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
!     Pack vector into doubles eta
!
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j - 1) + b
!
            do i = 1, wf%n_o
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i - 1) + a
!
                  aibj = max(ai, bj)*(max(ai,bj)-3)/2 + ai + bj
!
                  eta(wf%n_t1 + aibj, 1) = eta_ai_bj(ai, bj)
!
               enddo
            enddo
         enddo
      enddo
!
      call mem%dealloc(eta_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
   end subroutine construct_eta_ccsd
!
!
   subroutine initialize_multipliers_ccsd(wf)
!!
!!    Initialize multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Sep 2018 
!!
!!    Allocates the multipliers. This routine must be overwritten in 
!!    descendants which have more multipliers. 
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      call wf%initialize_t1bar()
      call wf%initialize_t2bar()
!
   end subroutine initialize_multipliers_ccsd
!
!
   subroutine set_multipliers_ccsd(wf, multipliers)
!!
!!    Set multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018 
!!
      implicit none 
!
      class(ccsd) :: wf  
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(in) :: multipliers
!
      call dcopy(wf%n_t1, multipliers, 1, wf%t1bar, 1)
      call dcopy(wf%n_t2, multipliers(wf%n_t1 + 1, 1), 1, wf%t2bar, 1)
!
   end subroutine set_multipliers_ccsd
!
!
   subroutine get_multipliers_ccsd(wf, multipliers)
!!
!!    Get multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf  
!
      real(dp), dimension(wf%n_amplitudes, 1) :: multipliers
!
      call dcopy(wf%n_t1, wf%t1bar, 1, multipliers, 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, multipliers(wf%n_t1 + 1, 1), 1)
!
   end subroutine get_multipliers_ccsd
!
!
   subroutine initialize_t2bar_ccsd(wf)
!!
!!    Initialize T2-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (.not. allocated(wf%t2bar)) call mem%alloc(wf%t2bar, wf%n_t2, 1)
!
   end subroutine initialize_t2bar_ccsd
!
!
   subroutine construct_multiplier_equation_ccsd(wf, equation)
!!
!!    Construct multiplier equation 
!!    Written by Eirik F. Kjønstad, Nov 2018 
!!
!!    Constructs 
!!
!!       t-bar^T A + eta,
!!
!!    and places the result in 'equation'.
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1), intent(inout) :: equation 
!
      real(dp), dimension(:,:), allocatable :: eta 
!
!     Copy the multipliers, eq. = t-bar 
!
      call dcopy(wf%n_t1, wf%t1bar, 1, equation, 1)
      call dcopy(wf%n_t2, wf%t2bar, 1, equation(wf%n_t1 + 1, 1), 1)
!
!     Transform the multipliers by A^T, eq. = t-bar^T A 
!
      call wf%jacobian_transpose_ccsd_transformation(equation)
!
!     Add eta, eq. = t-bar^T A + eta 
!
      call mem%alloc(eta, wf%n_amplitudes, 1)
      call wf%construct_eta(eta)
!
      call daxpy(wf%n_amplitudes, one, eta, 1, equation, 1)
!
      call mem%dealloc(eta, wf%n_amplitudes, 1)
!
   end subroutine construct_multiplier_equation_ccsd
!
!
   subroutine save_multipliers_ccsd(wf)
!!
!!    Save multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf 
!
      call wf%save_t1bar()
      call wf%save_t2bar()
!
   end subroutine save_multipliers_ccsd
!
!
   subroutine read_multipliers_ccsd(wf)
!!
!!    Read multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none 
!
      class(ccsd), intent(inout) :: wf 
!
      call wf%read_t1bar()
      call wf%read_t2bar()
!
   end subroutine read_multipliers_ccsd
!
!
   subroutine destruct_multipliers_ccsd(wf)
!!
!!    Destruct multipliers 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018 
!!
!!    Deallocates the multipliers. This routine must be overwritten in 
!!    descendants which have more multipliers. 
!!
      implicit none 
!
      class(ccsd) :: wf 
!
      call wf%destruct_t1bar()
      call wf%destruct_t2bar()
!
   end subroutine destruct_multipliers_ccsd
!
!
   subroutine save_t2bar_ccsd(wf)
!!
!!    Save t2bar 
!!    Written by Eirik F. Kjønstad, Nov 2018 
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf 
!
      type(file) :: t2bar_file 
!
      call t2bar_file%init('t2bar', 'sequential', 'unformatted')
!
      call disk%open_file(t2bar_file, 'write', 'rewind')
!
      write(t2bar_file%unit) wf%t2bar
!
      call disk%close_file(t2bar_file)      
!
   end subroutine save_t2bar_ccsd
!
!
   subroutine read_t2bar_ccsd(wf)
!!
!!    Save t2bar 
!!    Written by Eirik F. Kjønstad, Nov 2018 
!!
      implicit none 
!
      class(ccsd), intent(inout) :: wf 
!
      type(file) :: t2bar_file 
!
      call t2bar_file%init('t2bar', 'sequential', 'unformatted')
!
      call disk%open_file(t2bar_file, 'read', 'rewind')
!
      read(t2bar_file%unit) wf%t2bar
!
      call disk%close_file(t2bar_file)      
!
   end subroutine read_t2bar_ccsd
!
!
   subroutine destruct_t2bar_ccsd(wf)
!!
!!    Destruct T2-bar
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Nov 2018
!!
      implicit none
!
      class(ccsd) :: wf
!
      if (allocated(wf%t2bar)) call mem%dealloc(wf%t2bar, wf%n_amplitudes, 1)
!
   end subroutine destruct_t2bar_ccsd
!
!
   subroutine print_dominant_amplitudes_ccsd(wf)
!!
!!    Print dominant amplitudes 
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf 
!
      call wf%print_dominant_x1(wf%t1,'t')
      call wf%print_dominant_x2(wf%t2,'t')
!
   end subroutine print_dominant_amplitudes_ccsd
!
!
   subroutine print_dominant_x_amplitudes_ccsd(wf, x, tag)
!!
!!    Print dominant amplitudes  (TODO)
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf 
!
      real(dp), dimension(wf%n_amplitudes, 1) :: x 
!
      character(len=1) :: tag
!
      call wf%print_dominant_x1(x(1:wf%n_t1,1),tag)
      call wf%print_dominant_x2(x(wf%n_t1 + 1:wf%n_amplitudes,1),tag)
!
   end subroutine print_dominant_x_amplitudes_ccsd
!
!
   subroutine print_dominant_x2_ccsd(wf, x2, tag)
!!
!!    Print dominant x2   
!!    Written by Eirik F. Kjønstad, Dec 2018 
!!
!!    Prints the 20 most dominant double amplitudes,
!!    or sorts them if there are fewer than twenty of them.
!!
      implicit none 
!
      class(ccsd), intent(in) :: wf
!
      real(dp), dimension(wf%n_t2, 1) :: x2 
      character(len=1), intent(in)    :: tag 
!
      real(dp), dimension(:,:), allocatable :: abs_x2
!
      integer(i15), dimension(:,:), allocatable :: dominant_indices
      real(dp), dimension(:,:), allocatable     :: dominant_values
!
      integer(i15) :: n_elements, elm, i, a, j, b, ai, bj
!
!     Sort according to largest contributions
!
      call mem%alloc(abs_x2, wf%n_t2, 1)
      abs_x2 = abs(x2)
!
      n_elements = 20
      if (n_elements .gt. wf%n_t2) n_elements = wf%n_t2
!
      call mem%alloc(dominant_indices, n_elements, 1)
      call mem%alloc(dominant_values, n_elements, 1)
!
      dominant_indices = 0
      dominant_values  = zero 
      call get_n_highest(n_elements, wf%n_t2, abs_x2, dominant_values, dominant_indices)
!
!     Print largest contributions
!
      write(output%unit, '(/t6,a)') 'Largest double amplitudes:'
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
      write(output%unit, '(t6,a)')  '  a         i         b         j         ' // tag // '(ai,bj)             '
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
      do elm = 1, n_elements
!
         call invert_packed_index(dominant_indices(elm,1), ai, bj, (wf%n_o)*(wf%n_v))
         call invert_compound_index(ai, a, i, wf%n_v, wf%n_o)
         call invert_compound_index(bj, b, j, wf%n_v, wf%n_o)
!
         write(output%unit, '(t6,i3,7x,i3,7x,i3,7x,i3,5x,f19.12)') a, i, b, j, x2(dominant_indices(elm,1), 1)
!
      enddo
!
      write(output%unit, '(t6,a)')  '---------------------------------------------------------------'
!
      call mem%dealloc(dominant_indices, n_elements, 1)
      call mem%dealloc(dominant_values, n_elements, 1)
      call mem%dealloc(abs_x2, wf%n_t2, 1)
!
   end subroutine print_dominant_x2_ccsd
!
!
end module ccsd_class
