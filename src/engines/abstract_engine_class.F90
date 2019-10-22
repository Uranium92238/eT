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
module abstract_engine_class
!
!!
!!    Abstract engine class module
!!    Written by Tor S. Haugland, 2019
!!
!!    Each descendant must make an ignite-routine with the corresponding wavefunction.
!!
!
   use kinds
!
   use global_out,         only: output
   use timings_class,      only: timings
   use string_utilities,   only: convert_to_uppercase
!
   use wavefunction_class, only: wavefunction
   use hf_class,           only: hf
   use ccs_class,          only: ccs
!
   implicit none
!
   type, abstract :: abstract_engine
!
      character(len=200) :: name_
      character(len=200) :: tag
      character(len=200) :: description  
      character(len=200) :: author
!
      type(timings) :: timer ! Timer for engine. Obs! must be turned on in constructor
!
      character(len=150), dimension(:), allocatable :: tasks   ! The printed tasks of the engine. 
                                                               ! Should be set in constructor
!
   contains
!
      procedure :: print_banner   => print_banner_abstract_engine
      procedure :: print_timings  => print_timings_abstract_engine
!
   end type abstract_engine
!
!
contains
!
!
   subroutine print_banner_abstract_engine(engine, wf)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Prints: 
!!    
!!       - Engine name
!!       - Authors and date
!!       - Wavefunction type
!!       - Engine tasks
!!
!!    Dependancies:
!!
!!       - The printables of the engine must be set for each decendant (constructor)
!!
      implicit none
!
      class(abstract_engine), intent(in) :: engine
      class(wavefunction),    intent(in) :: wf
!
      integer :: task
!
      if (.not. allocated(engine%tasks)) call output%error_msg('Tasks of engine was not set. Do this in prepare.')
!
      call output%printf(":: (a0)", pl='minimal', fs='(//t3,a)', chars=[engine%name_])
      call output%printf(":: (a0)", pl='minimal', chars=[engine%author])
!
      call output%printf("(a0)",    pl='normal', fs='(/t3,a)', chars=[engine%description])
!
      call output%printf('This is a (a0) (a0) calculation. The following tasks will be performed:', &
                         pl='normal', ffs='(/t3,a)', fs='(t3,a)',                                   &
                         chars=[ convert_to_uppercase(wf%name_), engine%tag ] )
!     
      do task = 1, size(engine%tasks)
!
         call output%printf('- (a0)', pl='normal', fs='(t6,a)', chars = [engine%tasks(task)] )
!
      enddo
!
   end subroutine print_banner_abstract_engine
!
!
   subroutine print_timings_abstract_engine(engine, wf)
!!
!!    Print timings
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Prints the timings of the engine.
!!
      implicit none
!
      class(abstract_engine), intent(inout) :: engine
      class(wavefunction),    intent(in)    :: wf
!
      call output%printf('- Timings for the (a0) (a0) calculation', pl='minimal', fs='(/t3, a)', &
                         chars=[ convert_to_uppercase(wf%name_), engine%tag ])
!
      call engine%timer%turn_off()
!
      call output%printf('Total wall time (sec): (f20.5)', pl='normal', fs='(/t6, a)', &
                         reals=[engine%timer%get_elapsed_time('wall')])
      call output%printf('Total cpu time (sec):  (f20.5)', pl='normal', fs='(t6, a)', &
                         reals=[engine%timer%get_elapsed_time('cpu')])
!
   end subroutine print_timings_abstract_engine
!
!
end module abstract_engine_class
