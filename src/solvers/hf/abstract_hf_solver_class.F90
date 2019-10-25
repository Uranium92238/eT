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
module abstract_hf_solver_class
!!
!!    Abstract HF solver class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use global_in,  only : input
   use global_out, only : output
   use hf_class, only : hf
   use memory_manager_class, only : mem
!
   use parameters   
!
   implicit none 
!
   type, abstract :: abstract_hf_solver 
!
      character(len=100) :: tag
      character(len=100) :: author
      character(len=400) :: description
!
      real(dp) :: energy_threshold  
      real(dp) :: gradient_threshold
!
      integer :: max_iterations
!  
      character(len=200) :: ao_density_guess
!
      integer, dimension(:), allocatable :: orbitals_to_print
!
   contains 
!
      procedure :: print_banner             => print_banner_abstract_hf_solver
!
      procedure :: read_settings            => read_settings_abstract_hf_solver
      procedure :: read_hf_solver_settings  => read_hf_solver_settings_abstract_hf_solver
!
      procedure :: print_hf_solver_settings => print_hf_solver_settings_hf_solver
!
      procedure :: print_summary            => print_summary_abstract_hf_solver
!
      procedure, nopass :: run_single_ao    => run_single_ao_abstract_hf_solver 
!
   end type abstract_hf_solver
!
contains 
!
!
   subroutine run_single_ao_abstract_hf_solver(wf)
!!
!!    Run single AO 
!!    Written by Eirik F. Kjønstad, 2019 
!!
!!    Special case where n_ao == 1 means that the standard 
!!    solvers cannot be used. This is because the packed gradient 
!!    (anti-symmetric) does not exist. 
!!
!!    The equations are already solved, so this routine simply 
!!    constructs the required properties for saving in the cleanup
!!    routine.
!!
      implicit none 
!
      class(hf) :: wf 
!
      real(dp), dimension(:,:), allocatable :: h_wx 
!
      call output%printf('The system contains just one atomic orbital. Just constructing the solutions.', &
                           pl='m', fs='(/t3,a)')
!
      call mem%alloc(h_wx, wf%n_ao, wf%n_ao)
      call wf%get_ao_h_wx(h_wx)
!
      call wf%update_fock_and_energy(h_wx)  
!
      call mem%dealloc(h_wx, wf%n_ao, wf%n_ao)
!
      call wf%roothan_hall_update_orbitals() ! F => C
      call wf%update_ao_density()            ! C => D
!
      call wf%save_orbital_coefficients()
      call wf%save_orbital_energies()
!
   end subroutine run_single_ao_abstract_hf_solver
!
!
   subroutine print_hf_solver_settings_hf_solver(solver)
!!
!!    Print HF solver settings    
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      call output%printf('Energy threshold:             (e11.4)',&
         reals=[solver%energy_threshold], fs='(t6,a)', pl='minimal')
!
      call output%printf('Gradient threshold:           (e11.4)',&
         reals=[solver%gradient_threshold], fs='(t6,a)', pl='minimal')
!
   end subroutine print_hf_solver_settings_hf_solver
!
!
   subroutine read_settings_abstract_hf_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings. This routine is to be overwritten by 
!!    descendants if more settings need to be set. 
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      call solver%read_hf_solver_settings()
!
   end subroutine read_settings_abstract_hf_solver
!
!
   subroutine read_hf_solver_settings_abstract_hf_solver(solver)
!!
!!    Read HF solver settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings specific to this class.
!!
      use memory_manager_class, only: mem
!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      integer :: n_orbitals_to_print
!
      call input%get_keyword_in_section('energy threshold', 'solver scf', solver%energy_threshold)
      call input%get_keyword_in_section('gradient threshold', 'solver scf', solver%gradient_threshold)
      call input%get_keyword_in_section('max iterations', 'solver scf', solver%max_iterations)
      call input%get_keyword_in_section('ao density guess', 'solver scf', solver%ao_density_guess)
!
      if (input%requested_keyword_in_section('print orbitals', 'solver scf')) then
!
         n_orbitals_to_print = input%get_n_elements_for_keyword_in_section('print orbitals', 'solver scf')
!
         call mem%alloc(solver%orbitals_to_print, n_orbitals_to_print)
!
         call input%get_array_for_keyword_in_section('print orbitals', 'solver scf', n_orbitals_to_print, solver%orbitals_to_print)
!
      endif
!
   end subroutine read_hf_solver_settings_abstract_hf_solver
!
!
   subroutine print_banner_abstract_hf_solver(solver)
!!
!!    Print banner
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      call output%long_string_print(solver%tag,'(//t3,a)',.true.)
      call output%long_string_print(solver%author,'(t3,a/)',.true.)
      call output%long_string_print(solver%description)
!
   end subroutine print_banner_abstract_hf_solver
!
!
   subroutine print_summary_abstract_hf_solver(solver, wf)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use string_utilities, only: convert_to_uppercase
      use hf_class, only: hf 
!
      implicit none 
!
      class(abstract_hf_solver) :: solver
      class(hf), intent(inout) :: wf
!
      call output%printf('- Summary of '// trim(convert_to_uppercase(wf%name_))//&
             ' wavefunction energetics (a.u.):', fs='(/t3,a)', pl='minimal')
!
      call wf%print_energy()
      call wf%print_orbital_energies('3')
!
      if (allocated(solver%orbitals_to_print)) then
!      
         call wf%print_orbitals(size(solver%orbitals_to_print), solver%orbitals_to_print)
!
      else
!
         call wf%print_orbitals(wf%n_o + 10)
!
      endif
!
   end subroutine print_summary_abstract_hf_solver
!
end module abstract_hf_solver_class
