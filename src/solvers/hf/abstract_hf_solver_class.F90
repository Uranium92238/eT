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
   use kinds 
   use parameters   
   use file_class
   use disk_manager_class
   use io_utilities
!
   implicit none 
!
   type, abstract :: abstract_hf_solver 
!
      character(len=100) :: tag
      character(len=100) :: author
      character(len=400) :: description
!
      real(dp) :: energy_threshold          = 1.0D-6
      real(dp) :: gradient_threshold        = 1.0D-6
!
      integer :: max_iterations        = 100
!  
      character(len=200) :: ao_density_guess = 'SAD'
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
   end type abstract_hf_solver
!
contains 
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
      write(output%unit, '(t6,a30,e11.4)') 'Energy threshold:             ', solver%energy_threshold
      write(output%unit, '(t6,a30,e11.4)') 'Gradient threshold:           ', solver%gradient_threshold
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
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      call input%get_keyword_in_section('energy threshold', 'solver hf', solver%energy_threshold)
      call input%get_keyword_in_section('gradient threshold', 'solver hf', solver%gradient_threshold)
      call input%get_keyword_in_section('max iterations', 'solver hf', solver%max_iterations)
      call input%get_keyword_in_section('ao density guess', 'solver hf', solver%ao_density_guess)
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
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description)
!
   end subroutine print_banner_abstract_hf_solver
!
!
end module abstract_hf_solver_class
