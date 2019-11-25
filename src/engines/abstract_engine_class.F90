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
   use parameters
!
   use global_out,           only: output
   use timings_class,        only: timings
   use string_utilities,     only: convert_to_uppercase
   use memory_manager_class, only: mem
!
   use wavefunction_class, only: wavefunction
!
   implicit none
!
   type, abstract :: abstract_engine
!
      character(len=200) :: name_
      character(len=200) :: tag
      character(len=200) :: description  
      character(len=200) :: author
      logical :: dipole
      logical :: quadrupole
      logical :: plot_density
!
      type(timings) :: timer ! Timer for engine. Obs! must be turned on in constructor
!
      character(len=150), dimension(:), allocatable :: tasks   ! The printed tasks of the engine. 
                                                               ! Should be set in constructor
!
   contains
!
      procedure :: print_banner           => print_banner_abstract_engine
      procedure :: print_timings          => print_timings_abstract_engine
      procedure, nopass :: print_operator => print_operator_abstract_engine
!
      procedure, nopass :: remove_trace   => remove_trace_abstract_engine
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
!!    Dependencies:
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
                         chars=[character(len=500)::convert_to_uppercase(wf%name_), engine%tag] )
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
      call output%printf('- Timings for the (a0) (a0) calculation', pl='normal', fs='(/t3, a)', &
                         chars=[character(len=500)::convert_to_uppercase(wf%name_), engine%tag])
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
   subroutine print_operator_abstract_engine(operator_, electronic, nuclear, total, components, n_components)
!!
!!    Print operator
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!    Modified by Linda Goletto and Tommaso Giovannini, Oct 2019
!!
      use parameters
!!      
      implicit none
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(n_components), intent(in) :: electronic
      real(dp), dimension(n_components), intent(in) :: nuclear
      real(dp), dimension(n_components), intent(in) :: total
      real(dp)  :: dipole_norm
!
      character(len=4), dimension(n_components), intent(in) :: components
!
      character(len=*), intent(in) :: operator_
!
      integer :: k
!
      call output%printf('- Operator: (a0) [a.u.]', pl='minimal', fs='(/t3,a,a,a)', chars=[operator_])
!
      call output%printf(' Comp.   Electronic        Nuclear          Total   ', pl='verbose', fs='(/t6,a)')
      call output%printf('------------------------------------------------------', pl='verbose', fs='(t6,a)')
!
      do k = 1, n_components
!
         call output%printf('(a4)   (E14.7)  (E14.7)  (E14.7)', pl='verbose', fs='(t6,a)', &
                        chars=[components(k)], reals=[electronic(k), nuclear(k), total(k)])
!
      enddo
!
      call output%printf('------------------------------------------------------', pl='verbose', fs='(t6,a)')
!
!     For dipole moments calculate the norm and print Debye units 
! 
      if (index(trim(operator_), 'dipole').gt.0) then
!
         call output%printf('x:     (f14.7)', pl='minimal', fs='(/t6,a)', reals=[total(1)])
         call output%printf('y:     (f14.7)', pl='minimal', fs='(t6,a)', reals=[total(2)])
         call output%printf('z:     (f14.7)', pl='minimal', fs='(t6,a)', reals=[total(3)])
!      
         dipole_norm = dsqrt(total(1)**2 + total(2)**2 + total(3)**2)
!
         call output%printf('|mu|:  (f14.7)', pl='minimal', fs='(/t6,a)', reals=[dipole_norm])
!
         call output%printf('- Operator: (a0) [Debye]', pl='minimal', fs='(/t3,a,a,a)', chars=[operator_])
!         
         call output%printf(' Comp.   Electronic        Nuclear          Total   ', pl='verbose', fs='(/t6,a)')
         call output%printf('------------------------------------------------------', pl='verbose', fs='(t6,a)')
!
         do k = 1, n_components
!
            call output%printf('(a4)   (E14.7)  (E14.7)  (E14.7)', pl='verbose', fs='(t6,a)',   &
                                       chars=[components(k)], reals=[electronic(k)*au_to_debye, &
                                                                     nuclear(k)*au_to_debye,    &
                                                                     total(k)*au_to_debye])
!
         enddo
!
         call output%printf('------------------------------------------------------', pl='verbose', fs='(t6,a)')
!
         call output%printf('x:     (f14.7)', pl='minimal', fs='(/t6,a)', reals=[total(1)*au_to_debye])
         call output%printf('y:     (f14.7)', pl='minimal', fs='(t6,a)', reals=[total(2)*au_to_debye])
         call output%printf('z:     (f14.7)', pl='minimal', fs='(t6,a)', reals=[total(3)*au_to_debye])
!
         call output%printf('|mu|:  (f14.7)', pl='minimal', fs='(/t6,a)', reals=[dipole_norm*au_to_debye])
!
!     For quadrupole moments print Debye*Ang units 
!
      else if (index(trim(operator_), 'quadrupole').gt.0) then
!
         call output%printf('xx:    (f14.7)', pl='minimal', fs='(/t6,a)', reals=[total(1)])
!
         do k = 2, n_components
!
            call output%printf('(a4):    (f14.7)', pl='minimal', fs='(t4,a)', &
                                     chars=[components(k)], reals=[total(k)])
!
         enddo
!
         call output%printf('- Operator: (a0) [Debye*Ang]', pl='minimal', fs='(/t3,a,a,a)', chars=[operator_])
!          
         call output%printf(' Comp.   Electronic        Nuclear          Total   ', pl='verbose', fs='(/t6,a)')
         call output%printf('------------------------------------------------------', pl='verbose', fs='(t6,a)')
!
         do k = 1, n_components
!
            call output%printf('(a4)   (E14.7)  (E14.7)  (E14.7)', pl='verbose', fs='(t6,a)',      &
                          chars=[components(k)], reals=[electronic(k)*au_to_debye*bohr_to_angstrom,&
                                                        nuclear(k)*au_to_debye*bohr_to_angstrom,   &
                                                        total(k)*au_to_debye*bohr_to_angstrom])
!
         enddo
!
         call output%printf('------------------------------------------------------', pl='verbose', fs='(t6,a)')
!
         call output%printf('xx:    (f14.7)', pl='minimal', fs='(/t6,a)', reals=[total(1)*au_to_debye*bohr_to_angstrom])
!
         do k = 2, n_components
!
            call output%printf('(a4):    (f14.7)', pl='minimal', fs='(t4,a)',      &
                          chars=[components(k)], reals=[total(k)*au_to_debye*bohr_to_angstrom])
!
         enddo
!
      endif
!
   end subroutine print_operator_abstract_engine
!
!
   subroutine remove_trace_abstract_engine(M)
!!
!!    Remove trace
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
!!    The assumption here is that M is a 2-tensor ordered as xx, xy, xz, yy, yz, and zz,
!!    where the other elements of the tensor are given by symmetry, such as for the quadrupole
!!    moment. Thus, this routine can be called after a call to "calculate quadrupole moment"
!!    to make the moment trace-free.
!!
      implicit none
!
      real(dp), dimension(6), intent(inout) :: M
!
      real(dp) :: trace_
!
      trace_ = M(1) + M(4) + M(6)
!
      M(1) = (three*M(1) - trace_)/two
      M(4) = (three*M(4) - trace_)/two
      M(6) = (three*M(6) - trace_)/two
!
      M(2) = (three*M(2))/two
      M(3) = (three*M(3))/two
      M(5) = (three*M(5))/two
!
   end subroutine remove_trace_abstract_engine
!
!
end module abstract_engine_class
