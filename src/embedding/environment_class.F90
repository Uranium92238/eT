!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2021 the authors of eT
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
module environment_class
!
!!
!!    Abstract environment class module
!!    Written by Sarai D. Folkestad, Sep 2020
!!
!!    Defines the interface to the environment classes
!!
!
   use parameters
!
   use ao_tool_class,      only: ao_tool
   use mm_molecule_class,  only: mm_molecule
!
   implicit none
!
   type, abstract :: environment
!
      character(len=50) :: type_
      integer           :: n_charges
!
   contains
!
!     Public deferred class routines
!
      procedure (initialize_environment),    deferred :: initialize
      procedure (update_environment),        deferred :: update
      procedure (get_energy_environment),    deferred :: get_energy
      procedure (print_energy_environment),  deferred :: print_energy
!
   end type  environment
!
   abstract interface
!
      subroutine initialize_environment(embedding, ao)
!
         import environment
         import ao_tool
!
         implicit none
!
         class(environment),  intent(inout) :: embedding
         type(ao_tool),       intent(inout) :: ao
!
      end subroutine
!
      subroutine update_environment(embedding, ao, density)
!
         use parameters
!
         import environment
         import ao_tool
!
         implicit none
!
         class(environment),              intent(inout)  :: embedding
         type(ao_tool),                   intent(inout)  :: ao
         real(dp), dimension(ao%n, ao%n), intent(in)     :: density
!
      end subroutine
!
      real(dp) function get_energy_environment&
                        (embedding, ao, density) result(embedding_energy)
!
         use parameters
!
         import environment
         import ao_tool
!
         implicit none
!
         class(environment),              intent(in)     :: embedding
         type(ao_tool),                   intent(in)     :: ao
         real(dp), dimension(ao%n, ao%n), intent(in)     :: density
!
      end function
!
      subroutine print_energy_environment(embedding, ao, density)
!
         use parameters
!
         import environment
         import ao_tool
!
         implicit none
!
         class(environment)                           :: embedding
         type(ao_tool),                   intent(in)  :: ao
         real(dp), dimension(ao%n, ao%n), intent(in)  :: density
!
      end subroutine
!
   end interface
!
end module environment_class
