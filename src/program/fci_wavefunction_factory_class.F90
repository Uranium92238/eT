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
module fci_wavefunction_factory_class
!
!!
!!    FCI wavefunction factory class
!!    Written by Eirik F. Kjønstad, May 2022
!!
!
   use hf_class,  only: hf 
   use fci_class, only: fci 
!
   implicit none 
!
   type :: fci_wavefunction_factory
!
      character(len=30), private :: wf_name
!
   contains
!
      procedure, public :: create & 
                        => create_fci_wavefunction_factory
!
      procedure, private :: allocate_wavefunction
!
   end type fci_wavefunction_factory
!
!
   interface fci_wavefunction_factory
!
      procedure :: new_fci_wavefunction_factory
!
   end interface fci_wavefunction_factory
!
!
contains 
!
!
   function new_fci_wavefunction_factory() result(this)
!!
!!    New 
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_in, only: input 
!
      implicit none 
!
      type(fci_wavefunction_factory) :: this 
!
      this%wf_name = input%get_fci_wf()
!
   end function new_fci_wavefunction_factory
!
!
   subroutine create_fci_wavefunction_factory(this, ref_wf, fci_wf)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, May 2022
!!
      implicit none 
!
      class(fci_wavefunction_factory), intent(in) :: this 
!
      class(hf), intent(in) :: ref_wf 
      class(fci), allocatable, intent(out) :: fci_wf
!
      call this%allocate_wavefunction(fci_wf)
!
      call fci_wf%initialize(ref_wf)
!
   end subroutine create_fci_wavefunction_factory
!
!
   subroutine allocate_wavefunction(this, fci_wf)
!!
!!    Allocate wavefunction
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use fci_class,  only: fci
      use global_out, only: output
!
      implicit none
!
      class(fci_wavefunction_factory), intent(in) :: this 
!
      class(fci), allocatable, intent(out) :: fci_wf
!
      select case (trim(this%wf_name))
!
         case ('fci')
!
            allocate(fci::fci_wf)
!
         case default
!
            call output%error_msg('could not recognize FCI method ' // trim(this%wf_name) // '.')
!
      end select 
!
   end subroutine allocate_wavefunction
!
!
end module fci_wavefunction_factory_class
