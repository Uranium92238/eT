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
      procedure :: calculate_energy           => calculate_energy_mp2
!
   end type mp2
!
!
   interface
!
      include "zop_mp2_interface.F90"
!
   end interface
!
!
   interface mp2
!
      procedure :: new_mp2
!
   end interface mp2
!
!
contains
!
!
   function new_mp2(system) result(wf)
!!
!!    New MP2
!!    Written by Andreas Skeidsvoll, 2018
!!
      implicit none
!
      type(mp2) :: wf
!
      class(molecular_system), target, intent(in) :: system 
!
      type(sequential_file) :: hf_restart_file 
!
      wf%name_ = 'mp2'
!
      wf%system => system
!
      hf_restart_file = sequential_file('hf_restart_file')
!
      call hf_restart_file%open_('read', 'rewind')
!
      call hf_restart_file%read_(wf%n_ao)
      call hf_restart_file%read_(wf%n_mo)
      call hf_restart_file%read_()
      call hf_restart_file%read_(wf%n_o)
      call hf_restart_file%read_(wf%n_v)
      call hf_restart_file%read_(wf%hf_energy)
!
      call hf_restart_file%close_()
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
      call zero_array(wf%t1, wf%n_o*wf%n_v)
!
   end function new_mp2
!
!
end module mp2_class
