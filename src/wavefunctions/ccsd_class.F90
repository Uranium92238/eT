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
   contains
!
!     Preparation and cleanup routines 
!
      procedure :: prepare                      => prepare_ccsd
      procedure :: cleanup                      => cleanup_ccsd
!
!     Routines related to the amplitudes 
!
      procedure :: initialize_amplitudes        => initialize_amplitudes_ccsd 
      procedure :: initialize_t2                => initialize_t2_ccsd
      procedure :: destruct_t2                  => destruct_t2_ccsd
      procedure :: set_initial_amplitudes_guess => set_initial_amplitudes_guess_ccsd
      procedure :: set_t2_to_mp2_guess          => set_t2_to_mp2_guess_ccsd
      procedure :: set_amplitudes               => set_amplitudes_ccsd 
      procedure :: get_amplitudes               => get_amplitudes_ccsd 
!
!     Routines related to omega
!
      procedure :: construct_omega              => construct_omega_ccsd
!
      procedure :: omega_ccsd_a1                => omega_ccsd_a1_ccsd
      procedure :: omega_ccsd_b1                => omega_ccsd_b1_ccsd
      procedure :: omega_ccsd_c1                => omega_ccsd_c1_ccsd
!
      procedure :: omega_ccsd_a2                => omega_ccsd_a2_ccsd
      procedure :: omega_ccsd_b2                => omega_ccsd_b2_ccsd
      procedure :: omega_ccsd_c2                => omega_ccsd_c2_ccsd
      procedure :: omega_ccsd_d2                => omega_ccsd_d2_ccsd
      procedure :: omega_ccsd_e2                => omega_ccsd_e2_ccsd
!
   end type ccsd
!
!
   interface
!
      include "../submodules/ccsd/omega_ccsd_interface.F90"
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
      wf%n_amplitudes = (wf%n_o)*(wf%n_v) + (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2
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
      if (.not. allocated(wf%t2)) call mem%alloc(wf%t2, (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2, 1)
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
      if (allocated(wf%t2)) call mem%dealloc(wf%t2, (wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v) + 1)/2, 1)
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
      call dcopy((wf%n_o)*(wf%n_v), amplitudes, 1, wf%t1, 1)
      call dcopy((wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v)+1)/2, amplitudes((wf%n_o)*(wf%n_v)+1, 1), 1, wf%t2, 1)
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
      call dcopy((wf%n_o)*(wf%n_v), wf%t1, 1, amplitudes, 1)
      call dcopy((wf%n_o)*(wf%n_v)*((wf%n_o)*(wf%n_v)+1)/2, wf%t2, 1,  amplitudes((wf%n_o)*(wf%n_v)+1, 1), 1)
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
                     wf%t2(aibj, 1) = g_ai_bj(ai, bj)/(wf%orbital_energies(i, 1) + &
                                                       wf%orbital_energies(j, 1) - &
                                                       wf%orbital_energies(a + wf%n_o, 1) - &
                                                       wf%orbital_energies(b + wf%n_o, 1))
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
end module ccsd_class
