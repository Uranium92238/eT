module cc3_class
!
!!
!!    Coupled cluster CC3 class module
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
!
   use ccsd_class
!
   implicit none
!
   type, extends(ccsd) :: cc3
!
!  Integral files
   type(file)  :: g_bdck_t
   type(file)  :: g_ljck_t
   type(file)  :: g_dbkc_t
   type(file)  :: g_jlkc_t
   type(file)  :: L_jbkc_t
   type(file)  :: g_bdck_c
   type(file)  :: g_ljck_c
   type(file)  :: g_dbkc_c
   type(file)  :: g_jlkc_c
!
   contains
!
!     Preparation and cleanup routines
!
      procedure :: prepare                => prepare_cc3
      procedure :: cleanup                => cleanup_cc3
!
!     Routines related to omega
!
      procedure :: construct_omega        => construct_omega_cc3
!
      procedure :: omega_cc3_a            => omega_cc3_a_cc3
      procedure :: omega_cc3_integrals    => omega_cc3_integrals_cc3
      procedure :: omega_cc3_vvv_reader   => omega_cc3_vvv_reader_cc3
      procedure :: omega_cc3_ov_vv_reader => omega_cc3_ov_vv_reader_cc3
      procedure :: omega_cc3_W_calc       => omega_cc3_W_calc_cc3
      procedure :: omega_cc3_eps          => omega_cc3_eps_cc3
      procedure :: omega_cc3_omega1       => omega_cc3_omega1_cc3
      procedure :: omega_cc3_omega2       => omega_cc3_omega2_cc3
!
   end type cc3
!
!
   interface
!
      include "../submodules/cc3/omega_cc3_interface.F90"
!
   end interface
!
!
contains
!
!
   subroutine prepare_cc3(wf, ref_wf)
!!
!!    Prepare
!!    Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      class(hf) :: ref_wf
!
      integer :: p
!
      wf%name_ = 'cc3'
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
   end subroutine prepare_cc3
!
!
   subroutine cleanup_cc3(wf)
!!
!!    Cleanup
!!    Written by Alex C. Paul and Rolf H. Myhre 2018
!!
      implicit none
!
      class(cc3) :: wf
!
      write(output%unit, '(/t3,a,a,a)') '- Cleaning up ', trim(wf%name_), ' wavefunction'
!
   end subroutine cleanup_cc3
!
!
end module cc3_class
