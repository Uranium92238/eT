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
module hf_energy_function_class
!
!!
!!    HF energy function
!!    Written by Eirik F. Kjønstad, 2021
!!
!
   use kinds
   use hf_class,                        only: hf
   use function_class,                  only: function_
   use redundant_internal_coords_class, only: redundant_internal_coords
!
   implicit none
!
   type, extends(function_) :: hf_energy_function
!
      class(hf), pointer, private :: wf 
!
   contains 
!
      procedure, public :: get_parameters &
                        => get_parameters_hf_energy_function
!
      procedure, public :: set_parameters &
                        => set_parameters_hf_energy_function
!
      procedure, public :: evaluate &
                        => evaluate_hf_energy_function
!
      procedure, public :: evaluate_gradient &
                        => evaluate_gradient_hf_energy_function 
!
      procedure, public :: initialize &
                        => initialize_hf_energy_function
!
      procedure, private :: update_wavefunction
!
   end type hf_energy_function
!
!
   interface hf_energy_function
!
      procedure :: new_hf_energy_function
!
   end interface hf_energy_function
!
!
contains 
!
!
   function new_hf_energy_function(wf) result(this)
!!
!!    New HF energy function
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      type(hf_energy_function) :: this 
!
      class(hf), target :: wf 
!
      this%wf => wf 
!
      this%n_parameters = 3 * wf%n_atomic_centers
!
   end function new_hf_energy_function
!  
!
   function get_parameters_hf_energy_function(this) result(x)
!!
!!    Get parameters
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(hf_energy_function) :: this 
!
      real(dp), dimension(this%n_parameters) :: x 
!
      real(dp), dimension(3,this%wf%n_atomic_centers) :: geometry 
!
      geometry = this%wf%get_molecular_geometry()
!
      call dcopy(this%n_parameters, geometry, 1, x, 1)
!
   end function get_parameters_hf_energy_function
!
!
   subroutine set_parameters_hf_energy_function(this, x)
!!
!!    Set geometry
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(hf_energy_function) :: this 
!
      real(dp), dimension(this%n_parameters), intent(in) :: x 
!
      call this%wf%set_geometry(x, 'bohr')
!
      call this%wf%ao%print_centers('angstrom')
      !call this%wf%ao%print_centers('bohr')
!
      call this%update_wavefunction()
!
   end subroutine set_parameters_hf_energy_function
!
!
   function evaluate_hf_energy_function(this) result(energy)
!!
!!    Evaluate
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(hf_energy_function) :: this 
!
      real(dp) :: energy 
!
      energy = this%wf%energy
!
   end function evaluate_hf_energy_function 
!
!
   function evaluate_gradient_hf_energy_function(this) result(g)
!!
!!    Evaluate gradient
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(hf_energy_function) :: this 
!
      real(dp), dimension(this%n_parameters) :: g 
!
      call this%wf%construct_molecular_gradient(g)      
!
   end function evaluate_gradient_hf_energy_function
!
!
   subroutine update_wavefunction(this)
!!
!!    Update wavefunction 
!!    Written by Eirik F. Kjønstad, 2021
!!
      use hf_gs_engine_class, only: hf_gs_engine
!
      implicit none 
!
      class(hf_energy_function) :: this 
!
      type(hf_gs_engine), allocatable :: hf_engine 
!
      call this%wf%set_n_mo() 
!
      hf_engine = hf_gs_engine()
      call hf_engine%ignite(this%wf)
!
   end subroutine update_wavefunction
!
!
   subroutine initialize_hf_energy_function(this)
!!
!!    Initialize
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(hf_energy_function) :: this 
!
      call this%update_wavefunction()
!
   end subroutine initialize_hf_energy_function
!
!
end module hf_energy_function_class
