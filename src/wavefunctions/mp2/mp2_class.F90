!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2022 the authors of eT
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
   use ccs_class, only: ccs
!
   use parameters
   use global_out, only: output
   use memory_manager_class, only: mem
!
   implicit none
!
   type, extends(ccs) :: mp2
!
!     No unique variables yet
!
   contains
!
      procedure :: calculate_energy           => calculate_energy_mp2
      procedure :: print_gs_summary           => print_gs_summary_mp2
!
      procedure :: initialize                 => initialize_mp2
!
!
   end type mp2
!
!
   interface
!
      include "mean_value_mp2_interface.F90"
!
   end interface
!
contains
!
!
   subroutine initialize_mp2(wf, template_wf)
!!
!!    Initialize
!!    Written by Andreas Skeidsvoll, 2018
!!
      use wavefunction_class, only: wavefunction
      use array_initialization, only: zero_array
!
      implicit none
!
      class(mp2), intent(inout) :: wf
!
      class(wavefunction), intent(in) :: template_wf
!
      wf%name_ = 'mp2'
!
      call wf%general_cc_preparations()
      call wf%set_variables_from_template_wf(template_wf)
      call wf%print_banner()
!
      wf%n_t1            = (wf%n_o)*(wf%n_v)
      wf%n_gs_amplitudes = wf%n_t1
      wf%n_es_amplitudes = wf%n_t1
      wf%need_g_abcd     = .false.
!
      call wf%initialize_fock()
      call wf%initialize_amplitudes()
      call zero_array(wf%t1, wf%n_o*wf%n_v)
!
      call wf%print_amplitude_info()
!
   end subroutine initialize_mp2
!
!
   subroutine print_gs_summary_mp2(wf)
!!
!!    Print ground state summary
!!    Written by Eirik F. Kjønstad, Dec 2018
!!
      implicit none
!
      class(mp2), intent(inout) :: wf
!
      call output%printf('m', '- Ground state summary:', fs='(/t3,a)')
!
      call output%printf('m', 'HF energy:                (f19.12)', &
                         reals=[wf%hf_energy], fs='(/t6,a)')
!
      call output%printf('m', 'MP2 correction:           (f19.12)', &
                         reals=[ (wf%energy - wf%hf_energy) ], fs='(t6,a)')
!
      call output%printf('m', 'MP2 energy:               (f19.12)', &
                         reals=[wf%energy], fs='(t6,a)')
!
   end subroutine print_gs_summary_mp2
!
!
end module mp2_class
