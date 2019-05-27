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
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2018
!!
!
   use kinds
   use ccs_class
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
      type(timings) :: timer ! Timer for engine. Obs! must be turned on in prepare, off in cleanup. 
!
      character(len=150), dimension(:), allocatable :: tasks   ! The printed tasks of the engine. 
                                                               ! Should be set in prepare.
!
   contains
!
      procedure :: ignite => ignite_abstract_engine
!
      procedure(essential_engine_w_wf), deferred     :: run 
      procedure(essential_engine), deferred, private :: set_printables 
!
      procedure, non_overridable :: cleanup => cleanup_abstract_engine
!
      procedure, nopass :: do_cholesky => do_cholesky_abstract_engine       
!
      procedure, non_overridable :: print_banner => print_banner_abstract_engine
!
   end type abstract_engine
!
!
   abstract interface
!
      subroutine essential_engine(engine)
!
         import :: abstract_engine
!
         implicit none 
!
         class(abstract_engine) :: engine
!
      end subroutine essential_engine
!
!
      subroutine essential_engine_w_wf(engine, wf)
!
         import :: abstract_engine, ccs
!
         implicit none 
!
         class(abstract_engine) :: engine
!
         class(ccs) :: wf
!
      end subroutine essential_engine_w_wf
!
!
   end interface
!
contains
!
!
   subroutine ignite_abstract_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
!!    Banner, run & cleanup
!!
      implicit none
!
      class(abstract_engine) :: engine
!
      class(ccs) :: wf
!
      call engine%print_banner(wf)
      call engine%run(wf)
      call engine%cleanup(wf)
!
   end subroutine ignite_abstract_engine
!
!
   subroutine do_cholesky_abstract_engine(wf, orbital_coefficients)
!!
!!    Do Cholesky
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
!!    Cholesky decomposition of electronic repiulsion integrals
!!
!
      use eri_cd_class
!
      implicit none
!
      class(ccs), intent(inout) :: wf
!
      real(dp), dimension(wf%n_ao, wf%n_mo), intent(in) :: orbital_coefficients
!
      type(eri_cd) :: eri_chol_solver
!
!     Cholesky decoposition 
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%diagonal_test(wf%system)
      call eri_chol_solver%construct_mo_cholesky_vectors(wf%system, wf%n_mo, orbital_coefficients)
!
      call wf%integrals%prepare(eri_chol_solver%n_cholesky, wf%n_o, wf%n_v)
!
      call eri_chol_solver%cleanup()
!
   end subroutine do_cholesky_abstract_engine
!
!
   subroutine cleanup_abstract_engine(engine, wf)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Prints the timings of the engine. For now
!!    non-overridable procedure
!!
      implicit none
!
      class(abstract_engine), intent(inout)  :: engine
!
      class(ccs), intent(in)                 :: wf
!
      write(output%unit, '(/t3, a)') '- Finalizing the ' // trim(convert_to_uppercase(wf%name_)) // &
                                       ' ' // trim(engine%tag) // ' calculation'
!
      call engine%timer%turn_off()
!
      write(output%unit, '(/t6,a23,f20.5)')  'Total wall time (sec): ', engine%timer%get_elapsed_time('wall')
      write(output%unit, '(t6,a23,f20.5)')   'Total cpu time (sec):  ', engine%timer%get_elapsed_time('cpu')
!
   end subroutine cleanup_abstract_engine
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
!!       - The printables of the engine must be set for each decendant (set_printables and prepare)
!!
      implicit none
!
      class(abstract_engine), intent(in)  :: engine
!
      class(ccs), intent(in)              :: wf
!
      integer :: task
!
      character(len=500) :: calculation_type
!
      call engine%set_printables()
!
      if (.not. allocated(engine%tasks)) call output%error_msg('Tasks of engine was not set. Do this in prepare.')
!
      calculation_type  = 'This is a '// trim(convert_to_uppercase(wf%name_)) // ' ' // trim(engine%tag) // ' calculation.&
                           & The following tasks will be performed:'
!     
      call output%long_string_print(engine%name_,'(//t3,a)',.true.)
      call output%long_string_print(engine%author,'(t3,a/)',.true.)
      call output%long_string_print(engine%description,'(/t3,a)',.false.,'(t3,a)','(t3,a)') 
!
      write(output%unit, '(/t3,a/)') trim(calculation_type)
!
      do task = 1, size(engine%tasks)
!
         write(output%unit, '(t6, a2, a)') '- ', trim(engine%tasks(task))
!
      enddo
!
   end subroutine print_banner_abstract_engine
!
!
end module abstract_engine_class
!
