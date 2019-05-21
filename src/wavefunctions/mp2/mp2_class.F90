!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2019 the authors of eT
!
!  eT is free software: you can redistribute it and/or modify
!  it under the terms of the GNU General Public License as published by
!  the Free Software Foundation, either version 3 of the License, or
!  (at your option) any later version.
!
!  eT is distributed in the hope that it will be useful,
!  but WITHOUT ANY WARRANTY; without even the implied warranty of
!  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
!  GNU General Public License for more details.
!
!  You should have received a copy of the GNU General Public License
!  along with this program. If not, see <https://www.gnu.org/licenses/>.
!
!
module mp2_class
!
!!
!!    Second order Møller-Plesset pertubation theory (MP2/CCPT2) class module
!!    Written by Andreas Skeidsvoll, 2018
!!
!
   use wavefunction_class
   use hf_class
   use ccs_class
!
   implicit none
!
   type, extends(ccs) :: mp2
!
!     No unique variables yet
!
   contains
!
      procedure :: prepare                    => prepare_mp2
      procedure :: calculate_energy           => calculate_energy_mp2
      procedure :: print_wavefunction_summary => print_wavefunction_summary_mp2
!
   end type mp2
!
contains
!
!
   subroutine prepare_mp2(wf, system)
!!
!!    Prepare
!!    Written by Andreas Skeidsvoll, 2018
!!
      implicit none
!
      class(mp2) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      type(file) :: hf_restart_file 
!
      wf%name_ = 'mp2'
!
      wf%system => system
!
      call hf_restart_file%init('hf_restart_file', 'sequential', 'unformatted')
!
      call disk%open_file(hf_restart_file, 'read', 'rewind')
!
      read(hf_restart_file%unit) wf%n_ao 
      read(hf_restart_file%unit) wf%n_mo 
      read(hf_restart_file%unit) 
      read(hf_restart_file%unit) wf%n_o  
      read(hf_restart_file%unit) wf%n_v  
      read(hf_restart_file%unit) wf%hf_energy  
!
      call disk%close_file(hf_restart_file)
!
      call wf%initialize_files()
!
      call wf%initialize_orbital_coefficients()
      call wf%initialize_orbital_energies()
!
      call wf%read_orbital_coefficients()
      call wf%read_orbital_energies()
!
      call wf%initialize_amplitudes()
      wf%t1 = zero
!
   end subroutine prepare_mp2
!
!
   subroutine calculate_energy_mp2(wf)
!!
!!    Calculate energy
!!    Written by Andreas Skeidsvoll, 2018
!!
!!    Calculates the MP2 energy from HF energy, E_HF, vovo integrals, g_aibj, 
!!    and the orbital energies, eps. The total MP2 energy is calculated as
!!
!!       E = E_HF - sum_aibj g_aibj*L_aibj/(eps(a)+eps(b)-eps(i)-eps(j))
!!
!!    where
!!
!!       L_aibj = 2*g_aibj - g_ajbi.
!!
!!    On entry, it is assumed that the energy is equal to the HF energy 
!!    (i.e. the routine only adds the correction to the energy variable).
!!
      implicit none
!
      class(mp2), intent(inout) :: wf 
!
      real(dp), dimension(:,:,:,:), allocatable :: g_aibj
      real(dp), dimension(:,:,:,:), allocatable :: L_aibj
      real(dp), dimension(:), allocatable :: eps
!
      real(dp) :: e2_neg
      integer  :: a, i, b, j
!
      call mem%alloc(eps, wf%n_mo)
      eps = wf%orbital_energies
!
      call mem%alloc(g_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
      call wf%get_vovo(g_aibj)
!
      call mem%alloc(L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      L_aibj = two*g_aibj
      call add_1432_to_1234(-one, g_aibj, L_aibj, wf%n_v, wf%n_o, wf%n_v, wf%n_o)
!
      e2_neg = zero 
!
!$omp parallel do schedule(static) private(a, i, b, j) reduction(+:e2_neg)
      do j = 1, wf%n_o
         do b = 1, wf%n_v
            do i = 1, wf%n_o 
               do a = 1, wf%n_v
!
                  e2_neg = e2_neg + g_aibj(a, i, b, j)*L_aibj(a, i, b, j)/(eps(wf%n_o+a)+eps(wf%n_o+b)-eps(i)-eps(j))
!
               enddo
            enddo
         enddo
      enddo
!$omp end parallel do
!
      wf%energy = wf%hf_energy - e2_neg
!
   end subroutine calculate_energy_mp2
!
!
   subroutine print_wavefunction_summary_mp2(wf)
!!
!!    Print wavefunction summary 
!!    Written by Andreas Skeidsvoll, 2018
!!
!!    Prints information related to the wavefunction,
!!    most of which is meaningful only for a properly 
!!    converged wavefunction. Should be overwritten in 
!!    descendants if more or less or other information 
!!    is present. 
!!
      implicit none 
!
      class(mp2), intent(in) :: wf
!
      write(output%unit, '(/t3,a,a,a)') ':: Summary of ', trim(wf%name_), ' wavefunction energetics (a.u.)'
!
      write(output%unit, '(/t3,a26,f19.12)') 'HF energy:                ', wf%hf_energy
      write(output%unit, '(t3,a26,f19.12)')  'MP2 correction:           ', (wf%energy)-(wf%hf_energy)
      write(output%unit, '(t3,a26,f19.12)')  'MP2 energy:               ', wf%energy
!
   end subroutine print_wavefunction_summary_mp2
!
!
end module mp2_class
