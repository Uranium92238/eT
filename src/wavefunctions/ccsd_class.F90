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
      integer(i15) :: p, i, a
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
!     Nothing here yet
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
      real(dp), dimension(:,:), allocatable :: L_ai_bj
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
      real(dp), dimension(:,:), allocatable :: L_ia_J  ! L_ia^J
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
!
!
end module ccsd_class
