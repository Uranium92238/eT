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
   use index
!
   implicit none
!
   type, extends(ccs):: ccpt2
!
   contains
!
      procedure :: prepare          => prepare_ccpt2
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
      wf%system = ref_wf%system
!
      wf%n_ao   = ref_wf%n_ao
      wf%n_mo   = ref_wf%n_mo
      wf%n_o    = ref_wf%n_o
      wf%n_v    = ref_wf%n_v
!
      wf%hf_energy = ref_wf%energy
!
      call wf%initialize_orbital_coefficients()
      wf%orbital_coefficients = ref_wf%orbital_coefficients
!
      call wf%initialize_orbital_energies()
      wf%orbital_energies = ref_wf%orbital_energies
!
   end subroutine prepare_ccpt2
!
!
!
   subroutine calculate_energy_ccpt2(wf)
!
!!
!!  Calculates CCPT2 energy from vovo integrals and orbital energies
!!  Written by Andreas Skeidsvoll, 2018
!!
!
      implicit none
!
      class(ccpt2) :: wf 
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: L_ai_bj
      real(dp), dimension(:,:), allocatable :: eps
      real(dp) :: e2_neg
      integer(i15) :: a, i, b, j, ai, bj 
!
      call mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      call mem%alloc(eps, wf%n_o + wf%n_v, 1)
      eps = wf%orbital_energies
!
      call wf%get_vovo(g_ai_bj)
!
      L_ai_bj = two*g_ai_bj
!
      call add_1432_to_1234(-one, g_ai_bj, L_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!!$omp parallel do schedule(static) private(s,r,q,p,pq,rs,ps,rq) reduction(+:e2_neg)
      do j = 1, wf%n_o 
         do b = 1, wf%n_v
!
            bj = wf%n_v*(j-1) + b
!
            do i = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  ai = wf%n_v*(i-1) + a
!
                  e2_neg = e2_neg + g_ai_bj(ai, bj)*L_ai_bj(ai, bj)/(eps(a, 1)+eps(b, 1)-eps(i, 1)-eps(j, 1))
!
               enddo
            enddo
         enddo
      enddo
!!$omp end parallel do
!
      wf%energy = wf%hf_energy - e2_neg
!
   end subroutine calculate_energy_ccpt2
!
end module ccpt2_class
