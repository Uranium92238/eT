module ccpt2_class
!
!!
!!  Coupled cluster pertubation theory 2 class module
!!  Written by ?, 2018
!!
!
   use wavefunction_class
   use hf_class
   use ccs_class
!
   implicit none
!
   type, extends(ccs):: ccpt2
!
   contains
!
   procedure :: prepare => prepare_ccpt2
!
   procedure :: calculate_energy => calculate_energy_ccpt2
!
   end type ccpt2
!
contains
!
!
subroutine prepare_ccpt2(wf, ref_wf)
!!
!!    Prepare
!!    Written by ?, 2018
!!
      implicit none
!
      class(ccpt2) :: wf
!
      class(hf) :: ref_wf
!
      integer(i15) :: p
!
      wf%name = 'ccpt2'
!
      write(output%unit, *) 'Se 1'
      flush(output%unit)
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
      wf%n_amplitudes = (wf%n_o)*(wf%n_v)
!
      call wf%initialize_fock_ij()
      call wf%initialize_fock_ia()
      call wf%initialize_fock_ai()
      call wf%initialize_fock_ab()
!
      write(output%unit, *) 'Se 2'
      flush(output%unit)
!
      call wf%initialize_fock_diagonal()
!
      write(output%unit, *) 'Se 3'
      flush(output%unit)
!
      wf%fock_ij(:,:) = ref_wf%mo_fock(1 : wf%n_o , 1 : wf%n_o)
      wf%fock_ia(:,:) = ref_wf%mo_fock(1 : wf%n_o , wf%n_o + 1 : wf%n_v)
      wf%fock_ai(:,:) = ref_wf%mo_fock(1 : wf%n_o , wf%n_o + 1 : wf%n_v)
      wf%fock_ab(:,:) = ref_wf%mo_fock( wf%n_o + 1 : wf%n_v , wf%n_o + 1 : wf%n_v)
!
      write(output%unit, *) 'Se 4'
      flush(output%unit)
!
      do p = 1, wf%n_mo
!
         wf%fock_diagonal(p, 1) = ref_wf%mo_fock(p, p)
!
      enddo
!
      write(output%unit, *) 'Se 4.5'
      flush(output%unit)
!
      call ref_wf%mo_transform_and_save_h()
!
      call wf%initialize_orbital_coefficients()
      wf%orbital_coefficients = ref_wf%orbital_coefficients
!
      write(output%unit, *) 'Se 5'
      flush(output%unit)
!
   end subroutine prepare_ccpt2
!
!
!
   subroutine calculate_energy_ccpt2
!
      implicit none
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: g_aj_bi
!
      call mem%alloc(g_ai_bj, wf%n_o*wf%n_v, wf%n_o*wf%n_v)
      call mem%alloc(g_aj_bi, wf%n_o*wf%n_v, wf%n_o*wf%n_v)

      call wf%get_vovo(g_ai_bj)
   end subroutine calculate_energy_ccpt2