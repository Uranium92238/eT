module ccpt2_class
!
!!
!!    Coupled cluster pertubation theory 2 (CCPT2) class module
!!    Written by Andreas Skeidsvoll and Eirik F. Kjønstad, 2018
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
      procedure :: prepare                    => prepare_ccpt2
      procedure :: calculate_energy           => calculate_energy_ccpt2
      procedure :: print_wavefunction_summary => print_wavefunction_summary_ccpt2
!
   end type ccpt2
!
contains
!
!
   subroutine prepare_ccpt2(wf, ref_wf)
!!
!!    Prepare
!!    Written by Andreas Skeidsvoll and Eirik F. Kjønstad, 2018, based on prepare_ccs
!!
      implicit none
!
      class(ccpt2) :: wf
!
      class(hf) :: ref_wf
!
      integer(i15) :: p
!
      write(output%unit, *) 'Entered prepare ccpt2'
      flush(output%unit)
!
      wf%name = 'MP2'
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
      call wf%initialize_amplitudes()
      wf%t1 = zero
!
      write(output%unit, *) 'Exiting prepare ccpt2'
      flush(output%unit)
!
   end subroutine prepare_ccpt2
!
!
   subroutine calculate_energy_ccpt2(wf)
!!
!!    Calculate energy
!!    Written by Andreas Skeidsvoll and Eirik Kjønstad, 2018
!!
!!    Calculates MP2 energy from HF energy, E_HF, vovo integrals, g_aibj, and orbital energies, eps.
!!    Total MP2 energy is calculated using
!!
!!       E = E_HF - sum_aibj g_aibj*L_aibj/(eps(a)+eps(b)-eps(i)-eps(j))
!!
!!    where
!!
!!       L_aibj = 2*g_aibj - g_ajbi.
!!
      implicit none
!
      class(ccpt2) :: wf 
!
      real(dp), dimension(:,:), allocatable :: g_ai_bj
      real(dp), dimension(:,:), allocatable :: L_ai_bj
      real(dp), dimension(:,:), allocatable :: eps
      real(dp) :: e2_neg = zero
      integer(i15) :: a, i, b, j, ai, bj 
!
      call mem%alloc(eps, wf%n_mo, 1)
      call mem%alloc(g_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
      call mem%alloc(L_ai_bj, (wf%n_o)*(wf%n_v), (wf%n_o)*(wf%n_v))
!
      eps = wf%orbital_energies
!
      call wf%get_vovo(g_ai_bj)
!
      L_ai_bj = two*g_ai_bj
!
      call add_1432_to_1234(-one, g_ai_bj, L_ai_bj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
!$omp parallel do schedule(static) private(a, i, b, j, ai, bj) reduction(+:e2_neg)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
!
            bj = (wf%n_v)*(j-1) + b
!
            do i = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  ai = (wf%n_v)*(i-1) + a
!
                  e2_neg = e2_neg + g_ai_bj(ai, bj)*L_ai_bj(ai, bj)/(eps(wf%n_o+a,1)+eps(wf%n_o+b,1)-eps(i,1)-eps(j,1))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%energy = wf%hf_energy - e2_neg
!
   end subroutine calculate_energy_ccpt2
!
!
   subroutine print_wavefunction_summary_ccpt2(wf)
!!
!!    Print wavefunction summary 
!!    Written by Andreas Skeidsvoll, 2018, based on print_wavefunction_summary_hf
!!
!!    Prints information related to the wavefunction,
!!    most of which is meaningful only for a properly 
!!    converged wavefunction. Should be overwritten in 
!!    descendants if more or less or other information 
!!    is present. 
!!
      implicit none 
!
      class(ccpt2), intent(in) :: wf
!
      write(output%unit, '(/t3,a,a,a)') ':: Summary of ', trim(wf%name), ' wavefunction energetics (a.u.)'
!
      write(output%unit, '(/t3,a26,f19.12)') 'HF energy:                ', wf%hf_energy
      write(output%unit, '(t3,a26,f19.12)')  'MP2 correction:           ', (wf%energy)-(wf%hf_energy)
      write(output%unit, '(t3,a26,f19.12)')  'MP2 energy:               ', wf%energy
!
      end subroutine print_wavefunction_summary_ccpt2
!
end module ccpt2_class