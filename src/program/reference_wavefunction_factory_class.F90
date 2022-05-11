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
module reference_wavefunction_factory_class
!
!!
!!    Reference wavefunction factory class
!!    Written by Eirik F. Kjønstad, May 2022
!!
!
   use hf_class,  only: hf
!
   implicit none 
!
   type :: reference_wavefunction_factory
!
      character(len=30), private :: wf_name
!
   contains
!
      procedure, public :: create & 
                        => create_reference_wavefunction_factory
!
      procedure, private :: allocate_wavefunction
!
   end type reference_wavefunction_factory
!
!
   interface reference_wavefunction_factory
!
      procedure :: new_reference_wavefunction_factory
!
   end interface reference_wavefunction_factory
!
!
contains 
!
!
   function new_reference_wavefunction_factory() result(this)
!!
!!    New 
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_in, only: input 
!
      implicit none 
!
      type(reference_wavefunction_factory) :: this 
!
      this%wf_name = input%get_reference_wf()
!
   end function new_reference_wavefunction_factory
!
!
   subroutine create_reference_wavefunction_factory(this, wf)
!!
!!    Create
!!    Written by Eirik F. Kjønstad, May 2022
!!
      implicit none 
!
      class(reference_wavefunction_factory), intent(in) :: this 
!
      class(hf), allocatable, intent(out) :: wf
!
      call this%allocate_wavefunction(wf)
!
      call wf%prepare()
!
   end subroutine create_reference_wavefunction_factory
!
!
   subroutine allocate_wavefunction(this, wf)
!!
!!    Allocate wavefunction
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use uhf_class,    only: uhf
      use mlhf_class,   only: mlhf
      use qed_hf_class, only: qed_hf
      use cuhf_class,   only: cuhf
      use rohf_class,   only: rohf
!
      use global_out,   only: output
!
      implicit none
!
      class(reference_wavefunction_factory), intent(in) :: this 
!
      class(hf), allocatable, intent(out)  :: wf
!
      select case (trim(this%wf_name))
!
         case ('hf')
!
            wf = hf()
!
         case ('uhf')
!
            wf = uhf()
!
         case ('cuhf')
!
            wf = cuhf()
!
         case ('rohf')
!
            wf = rohf()
!
         case ('mlhf')
!
            wf = mlhf()
!
         case ('qed-hf')
!
            wf = qed_hf()
!
         case default
!
            call output%error_msg('did not recognize the reference wavefunction ' &
                                    // trim(this%wf_name) //'.')
!
      end select 
!
   end subroutine allocate_wavefunction
!
!
end module reference_wavefunction_factory_class
