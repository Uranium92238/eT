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
module cc_wavefunction_factory_class
!
!!
!!    CC wavefunction factory class
!!    Written by Eirik F. Kjønstad, May 2022
!!
!
   use ccs_class, only: ccs 
   use hf_class,  only: hf
!
   implicit none 
!
   type :: cc_wavefunction_factory
!
      character(len=30), private :: wf_name
!
   contains
!
      procedure, public :: create & 
                        => create_cc_wavefunction_factory
!
      procedure, private :: allocate_wavefunction
!
   end type cc_wavefunction_factory
!
!
   interface cc_wavefunction_factory
!
      procedure :: new_cc_wavefunction_factory
!
   end interface cc_wavefunction_factory
!
!
contains 
!
!
   function new_cc_wavefunction_factory() result(this)
!!
!!    New 
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_in, only: input 
!
      implicit none 
!
      type(cc_wavefunction_factory) :: this 
!
      this%wf_name = input%get_cc_wf()
!
   end function new_cc_wavefunction_factory
!
!
   subroutine create_cc_wavefunction_factory(this, ref_wf, cc_wf)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, May 2022
!!
      implicit none 
!
      class(cc_wavefunction_factory), intent(in) :: this 
!
      class(hf), intent(in) :: ref_wf
!
      class(ccs), allocatable, intent(out) :: cc_wf 
!
      call this%allocate_wavefunction(cc_wf)
!
      call cc_wf%initialize(ref_wf)
!
   end subroutine create_cc_wavefunction_factory
!
!
   subroutine allocate_wavefunction(this, cc_wf)
!!
!!    Allocate wavefunction
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use cc2_class,          only: cc2
      use lowmem_cc2_class,   only: lowmem_cc2
      use ccsd_class,         only: ccsd
      use cc3_class,          only: cc3
      use ccsdpt_class,       only: ccsdpt
      use mp2_class,          only: mp2
      use mlcc2_class,        only: mlcc2
      use mlccsd_class,       only: mlccsd
!
      use global_out,         only: output
!
      implicit none 
!
      class(cc_wavefunction_factory), intent(in) :: this 
!
      class(ccs), allocatable, intent(out) :: cc_wf 
!
      select case (trim(this%wf_name))
!
         case ('ccs')
!
            allocate(ccs::cc_wf)
!
         case ('cc2')
!
            allocate(cc2::cc_wf)
!
         case ('lowmem-cc2')
!
            allocate(lowmem_cc2::cc_wf)
!
         case ('ccsd')
!
            allocate(ccsd::cc_wf)
!
         case ('cc3')
!
            allocate(cc3::cc_wf)
!
         case ('ccsd(t)')
!
            allocate(ccsdpt::cc_wf)
!
         case ('mp2')
!
            allocate(mp2::cc_wf)
!
         case ('mlcc2')
!
            allocate(mlcc2::cc_wf)
!
         case ('mlccsd')
!
            allocate(mlccsd::cc_wf)
!
         case default
!
            call output%error_msg('could not recognize CC method ' // trim(this%wf_name) // '.')
!
      end select
!
   end subroutine allocate_wavefunction
!
!
end module cc_wavefunction_factory_class
