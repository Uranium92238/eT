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
module zop_engine_class
!!
!!    Zeroth order coupled cluster engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Calculates expectation values < Λ | A | CC > for the operator A,
!!    where A is the dipole or quadrupole operator. 
!!
   use gs_engine_class
!
   type, extends(gs_engine) :: zop_engine
!
      character(len=100) :: tag           = 'Zeroth order coupled cluster properties'
      character(len=100) :: author        = 'E. F. Kjønstad, S. D. Folkestad, 2019'
!
      character(len=500) :: description1  = 'Calculates the time-independent expectation value of&
                                             & one-electron operators A, < A > = < Λ | A | CC >.'
!
      logical :: dipole 
      logical :: quadrupole 
!
   contains
!
      procedure :: prepare                                     => prepare_zop_engine
      procedure :: run                                         => run_zop_engine
!
      procedure :: read_settings                               => read_settings_zop_engine
      procedure :: read_zop_settings                           => read_zop_settings_zop_engine
!
      procedure :: calculate_expectation_values                => calculate_expectation_values_zop_engine
!
      procedure, nopass, private :: print_summary              => print_summary_zop_engine
!
   end type zop_engine
!
contains
!
   subroutine prepare_zop_engine(engine)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(zop_engine) :: engine
!
      engine%name_ = 'Zeroth order properties engine'
!
      engine%dipole                 = .false.
      engine%quadrupole             = .false. 
      engine%multipliers_algorithm  = 'davidson'
!
      call engine%read_settings()
!
   end subroutine prepare_zop_engine
!
!
   subroutine read_settings_zop_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(zop_engine) :: engine 
!
      call engine%read_gs_settings()
      call engine%read_zop_settings()
!
   end subroutine read_settings_zop_engine
!
!
   subroutine read_zop_settings_zop_engine(engine)
!!
!!    Read ZOP settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(zop_engine) :: engine 
!
      if (input%requested_keyword_in_section('dipole','cc zop')) engine%dipole = .true.
      if (input%requested_keyword_in_section('quadrupole','cc zop')) engine%quadrupole = .true.
!
   end subroutine read_zop_settings_zop_engine
!
!
   subroutine run_zop_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs)         :: wf
      class(zop_engine)  :: engine
!
!     Cholesky decoposition of the electron repulsion integrals 
!
      call engine%do_cholesky(wf, wf%orbital_coefficients)
!
!     Determine ground state | CC > 
!
      call engine%do_ground_state(wf)
!
!     Determine multipliers < Λ |
!
      call engine%do_multipliers(wf)
!
!     Compute the one-electron density 
!
      call wf%initialize_density()
      call wf%construct_density()
!
!     Compute the requested expectation values
!
      call long_string_print(engine%tag,'(//t3,a)',.true.)
      call long_string_print(engine%author,'(t3,a/)',.true.)
      call long_string_print(engine%description1,'(t3,a)',.false.,'(t3,a)','(t3,a)')
!
      call engine%calculate_expectation_values(wf)
!
      call wf%destruct_density()
!
   end subroutine run_zop_engine
!
!
   subroutine calculate_expectation_values_zop_engine(engine, wf)
!!
!!    Calculate expectation values 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019 
!!
      implicit none 
!
      class(zop_engine), intent(in) :: engine 
!
      class(ccs), intent(in) :: wf 
!
      real(dp), dimension(3) :: mu_electronic
      real(dp), dimension(3) :: mu_nuclear
      real(dp), dimension(3) :: mu_total
!
      real(dp), dimension(6) :: q_electronic
      real(dp), dimension(6) :: q_nuclear
      real(dp), dimension(6) :: q_total
!
      character(len=4), dimension(:), allocatable :: components
!
      if (engine%dipole) then 
!
         call engine%calculate_dipole_moment(wf, mu_electronic, mu_nuclear, mu_total)
!
         allocate(components(3))
!
         components = (/'x   ',&
                        'y   ',&
                        'z   '/)  
!
         call engine%print_summary('dipole moment', mu_electronic, mu_nuclear, mu_total, &
                                    components, 3)     
!
         deallocate(components)  
!
      endif
!
      if (engine%quadrupole) then 
!
         call engine%calculate_quadrupole_moment(wf, q_electronic, q_nuclear, q_total)
!
         allocate(components(6))
!
         components = (/ 'xx  ',   &
                         'xy  ',   & 
                         'xz  ',   & 
                         'yy  ',   & 
                         'yz  ',   & 
                         'zz  '    /) 
!
         call engine%print_summary('quadrupole moment (with trace)', q_electronic, q_nuclear, q_total, &
                                    components, 6)     
!
         call engine%remove_trace(q_electronic)
         call engine%remove_trace(q_nuclear)
!
         q_total = q_electronic + q_nuclear
!
         call engine%print_summary('traceless quadrupole moment', q_electronic, q_nuclear, q_total, &
                                    components, 6)   
!
         deallocate(components)  
!
      endif 
!
   end subroutine calculate_expectation_values_zop_engine
!
!
   subroutine print_summary_zop_engine(operator_, electronic, nuclear, total, components, n_components)
!!
!!    Print summary
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      integer, intent(in) :: n_components
!
      real(dp), dimension(n_components), intent(in) :: electronic
      real(dp), dimension(n_components), intent(in) :: nuclear
      real(dp), dimension(n_components), intent(in) :: total 
!
      character(len=4), dimension(n_components), intent(in) :: components
!
      character(len=*), intent(in) :: operator_ 
!
      integer :: k
!
      write(output%unit, '(/t3,a,a)') '- Operator: ', trim(operator_)
!
      write(output%unit, '(/t6, a)') 'Cart. comp.  Electronic         Nuclear             Total           '
      write(output%unit, '(t6, a)')  '--------------------------------------------------------------------'
!
      do k = 1, n_components
!
         write(output%unit, '(t9, a4, 3x, f19.12, f19.12, f19.12)') components(k), &
               electronic(k), nuclear(k), total(k)
!
      enddo
!
      write(output%unit, '(t6, a)')  '--------------------------------------------------------------------'
      flush(output%unit)
!
   end subroutine print_summary_zop_engine
!
!
end module zop_engine_class
