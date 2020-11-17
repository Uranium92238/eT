!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
!
   use hf_class,                        only : hf
   use memory_manager_class,            only : mem
   use timings_class,                   only : timings
   use abstract_convergence_tool_class, only : abstract_convergence_tool 
!
   use parameters   
!
   implicit none 
!
   type, abstract :: abstract_hf_solver 
!
      character(len=100) :: name_
      character(len=100) :: tag
      character(len=400) :: description
!
      integer :: max_iterations
!  
      character(len=200) :: ao_density_guess
!
      integer, dimension(:), allocatable :: orbitals_to_print
!
      type(timings), allocatable :: timer 
!
      logical :: energy_convergence
!  
      class(abstract_convergence_tool), allocatable :: convergence_checker
!
   contains 
!
      procedure :: print_banner             => print_banner_abstract_hf_solver
!
      procedure :: read_settings            => read_settings_abstract_hf_solver
      procedure :: read_hf_solver_settings  => read_hf_solver_settings_abstract_hf_solver
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
      call output%printf('m', 'The system contains just one atomic orbital. &
                         &Just constructing the solutions.', fs='(/t3,a)')
!
      call wf%update_fock_and_energy()  
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
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      real(dp) :: energy_threshold, gradient_threshold
!!
      if (input%requested_keyword_in_section('energy threshold', 'solver scf')) then
!
         call input%get_keyword_in_section('energy threshold', 'solver scf', energy_threshold)
         call solver%convergence_checker%set_energy_threshold(energy_threshold)
!
      endif
!
      if (input%requested_keyword_in_section('gradient threshold', 'solver scf')) then
!
         call input%get_keyword_in_section('gradient threshold', 'solver scf', gradient_threshold)
         call solver%convergence_checker%set_residual_threshold(gradient_threshold)
!
      endif
!
      call input%get_keyword_in_section('max iterations', 'solver scf', solver%max_iterations)
      call input%get_keyword_in_section('ao density guess', 'solver scf', solver%ao_density_guess)
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
      call output%printf('m', ' - ' // trim(solver%name_), fs='(/t3,a)')
      call output%print_separator('m', len(trim(solver%name_)) + 6, '-')
!
      call output%printf('n', '(a0)', ffs='(/t3,a)',  chars=[trim(solver%description)])
!
   end subroutine print_banner_abstract_hf_solver
!
!
end module abstract_hf_solver_class
