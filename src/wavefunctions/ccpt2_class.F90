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
!   use reordering
!   use array_utilities
!   use array_analysis
!   use interval_class
   use index
!
   implicit none
!
   type, extends(ccs):: ccpt2
!
   real(dp) :: ccpt2_energy
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
   subroutine calculate_energy_ccpt2(wf)
!
      implicit none
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: L_ai_bj
      real(dp), dimension(:,:), allocatable :: eps
      real(dp) :: e2
      integer(i15) :: p, q, r, s, pq, rs, qp, sr
      integer(i15) :: dim_p, dim_q, dim_r, dim_s
!
      call mem%alloc(g_ai_bj, wf%n_amplitudes, wf%n_amplitudes)
      call mem%alloc(L_ai_bj, wf%n_amplitudes, wf%n_amplitudes)
      call mem%alloc(eps, wf%n_o + wf_nv, 1)
!
      eps = wf%orbital_energies
!
      call wf%get_vovo(g_ai_bj)
!
      L_ai_bj = 2.0*g_ai_bj
!
      add_1432_to_1234(-0.5, g_ai_bj, L_ai_bj, wf%n_v, wf%n_o, wf_nv, wf_no)
!
      dim_p = wf%n_v
      dim_q = wf%n_o
      dim_r = wf%n_v
      dim_s = wf%n_o
!
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq) reduce(+:e2)
      do s = 1, dim_s
         do r = 1, dim_r
!
            rs = dim_r*(s-1) + r
            sr = dim_s*(r-1) + s
!
            do q = 1, dim_q
!
               do p = 1, dim_p
!
                  pq = dim_p*(q-1) + p
                  qp = dim_q*(p-1) + q
!
                  e2 = e2 - g_ai_bj(pq, rs)*L_ai_bj(qp, sr)/(eps(p)+eps(r)-eps(q)-eps(s))
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
      wf%ccpt2_energy = wf%hf_energy + e2    
   end subroutine calculate_energy_ccpt2