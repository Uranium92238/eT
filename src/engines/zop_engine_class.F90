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
   use abstract_engine_class
   use diis_cc_gs_class
   use davidson_cc_multipliers_class
   use diis_cc_multipliers_class
   use ccs_class
   use eri_cd_class
!
   type, extends(abstract_engine) :: zop_engine
!
      character(len=200) :: operator 
!
      integer :: n_components
!
      character(len=4), dimension(:), allocatable :: components
!
   contains
!
      procedure :: prepare                   => prepare_zop_engine
      procedure :: run                       => run_zop_engine
      procedure :: cleanup                   => cleanup_zop_engine
!
      procedure :: read_settings             => read_settings_zop_engine
!
      procedure, private :: set_components   => set_components_zop_engine
      procedure, private :: set_n_components => set_n_components_zop_engine
!
      procedure, private :: construct_operator => construct_operator_zop_engine
      procedure, private :: calculate_nuclear_contribution => calculate_nuclear_contribution_zop_engine
!
      procedure, private :: print_summary => print_summary_zop_engine
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
!     Set standards and then read if nonstandard
!
      call engine%read_settings()
!
      call engine%set_n_components()
!
      call engine%set_components()
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
      call input%get_keyword_in_section('operator', 'cc zop', engine%operator)
!
   end subroutine read_settings_zop_engine
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
!
      class(zop_engine)  :: engine
!
      type(eri_cd)                  :: eri_chol_solver
      type(diis_cc_gs)              :: cc_gs_solver
      type(davidson_cc_multipliers) :: cc_mult_solver
!
      real(dp), dimension(:,:,:), allocatable :: A 
!
      real(dp), dimension(:), allocatable :: expectation_value  
! 
      real(dp), dimension(:), allocatable :: nuclear_contribution
!
      real(dp) :: r2 
!
      integer :: component
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)
!
!     Cholesky decomposition
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test(wf%system)
!
      call eri_chol_solver%construct_mo_cholesky_vecs(wf%system, wf%n_mo, wf%orbital_coefficients)
!
      call wf%integrals%prepare(eri_chol_solver%n_cholesky, wf%n_o, wf%n_v)
!
      call eri_chol_solver%cleanup()
!
!     Ground state solution
!
      call cc_gs_solver%prepare(wf)
      call cc_gs_solver%run(wf)
      call cc_gs_solver%cleanup(wf)
!
!     Multiplier solution 
!
      call wf%initialize_amplitudes()
      call wf%read_amplitudes()
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      call wf%integrals%can_we_keep_g_pqrs_t1()
!
      call cc_mult_solver%prepare(wf)
      call cc_mult_solver%run(wf)
      call cc_mult_solver%cleanup(wf)
!
!     Compute expectation value of A = (A_x A_y A_z) for the operator A 
!
      call wf%initialize_multipliers()
      call wf%read_multipliers()
!
!     Determine the number of components in operator and construct it
!
      call mem%alloc(A, wf%n_mo, wf%n_mo, engine%n_components)
      call mem%alloc(expectation_value, engine%n_components)
      call mem%alloc(nuclear_contribution, engine%n_components)
!
      call engine%construct_operator(wf, A)
      call engine%calculate_nuclear_contribution(wf, nuclear_contribution)
!
      call wf%initialize_density()
      call wf%construct_density()
!
      do component = 1, engine%n_components
!
         expectation_value(component) = wf%calculate_expectation_value(A(:,:,component))
!
      enddo
!
      if (engine%traceless) call engine%remove_trace(expectation_value)
!xx, xy, xz, yy, yz, and zz.
      r2 = expectation_value(1) + expectation_value(4) + expectation_value(6)
!
      write(output%unit, *) 'xx:', (three*expectation_value(1) - r2)/two
      write(output%unit, *) 'yy:', (three*expectation_value(4) - r2)/two
      write(output%unit, *) 'zz:', (three*expectation_value(6) - r2)/two
      write(output%unit, *) 'xy:', (three*expectation_value(2))/two
      write(output%unit, *) 'xz:', (three*expectation_value(3))/two
      write(output%unit, *) 'yz:', (three*expectation_value(5))/two
!
      call engine%print_summary(expectation_value, nuclear_contribution) 
!
      call mem%dealloc(A, wf%n_mo, wf%n_mo, engine%n_components)
      call mem%dealloc(expectation_value, engine%n_components)
!
      call wf%destruct_density()
!
   end subroutine run_zop_engine
!
!
   subroutine remove_trace_zop_engine(engine, M)
!!
!!    Remove trace 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
!!    The assumption here is that M is a 2-tensor ordered as xx, xy, xz, yy, yz, and zz.
!!
      implicit none 
!
      class(zop_engine), intent(in) :: engine 
!
      real(dp), dimension(engine%n_components), intent(inout) :: M
!
      real(dp) :: trace_
!
      real(dp), dimension(engine%n_components) :: temp_M 
!
      if (trim(engine%operator) /= 'quadrupole') then 
!
         call output%error_msg('Cannot remove trace for operator ' // trim(engine%operator))
!
      else 
!
         trace_ = M(1) + M(4) + M(6)
!
         
!
      endif 
!
   end subroutine remove_trace_zop_engine
!
!
   subroutine cleanup_zop_engine(engine)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(zop_engine) :: engine
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(engine%name_)
!
   end subroutine cleanup_zop_engine
!
!
   subroutine set_n_components_zop_engine(engine)
!!
!!    Set n components 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad
!!
!!    Sets the number of components given the operator
!!
      implicit none
!
      class(zop_engine), intent(inout) :: engine
!
      engine%n_components = 0
!
      if (trim(engine%operator) == 'dipole') then
!
         engine%n_components = 3
!
      elseif (trim(engine%operator) == 'quadrupole') then 
!
         engine%n_components = 6
!
      else
!
         call output%error_msg('could not recognize the operator ' &
               // trim(engine%operator) //' in zop engine.')
!
      endif
!
   end subroutine set_n_components_zop_engine
!
!
   subroutine set_components_zop_engine(engine)
!!
!!    Set components 
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
!!    Sets the Cartesian components for given the operator
!!
      implicit none
!
      class(zop_engine), intent(inout) :: engine
!
      if (trim(engine%operator) == 'dipole') then
!
         engine%components = (/'x   ',&
                               'y   ',&
                               'z   '/)
!
      elseif (trim(engine%operator) == 'quadrupole') then 
!
         engine%components = (/ 'xx  ',   &
                                'xy  ',   & 
                                'xz  ',   & 
                                'yy  ',   & 
                                'yz  ',   & 
                                'zz  '    /)
!
      else
!
         call output%error_msg('could not recognize the operator ' &
               // trim(engine%operator) //' in zop engine.')
!
      endif
!
   end subroutine set_components_zop_engine
!
!
   subroutine construct_operator_zop_engine(engine, wf, A)
!!
!!    Construct operator
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(zop_engine), intent(in) :: engine
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(wf%n_mo, wf%n_mo, engine%n_components), intent(out) :: A
!
      if (trim(engine%operator) == 'dipole') then 
!
         call wf%construct_mu(A)
!
      elseif (trim(engine%operator) == 'quadrupole') then 
!
         call wf%construct_q(A)
!
      else      
!
         call output%error_msg('Tried to construct unrecognized one-electron integral matrix '&
                // trim(engine%operator))
!
      endif
!
   end subroutine construct_operator_zop_engine
!
!
   subroutine calculate_nuclear_contribution_zop_engine(engine, wf, nuclear_contribution)
!!
!!    Calculate nuclear contribution
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(zop_engine), intent(in) :: engine
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(engine%n_components), intent(out) :: nuclear_contribution
!
      if (trim(engine%operator) == 'dipole') then 
!
         call wf%system%get_nuclear_dipole(nuclear_contribution)
!
      elseif (trim(engine%operator) == 'quadrupole') then 
!
!        TODO!
!
      else
!
         call output%error_msg('Tried to construct nuclear contribution for unknown operator '&
                // trim(engine%operator))
!
      endif
!
   end subroutine calculate_nuclear_contribution_zop_engine
!
!
   subroutine print_summary_zop_engine(engine, electronic, nuclear)
!!
!!    Print summary
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
      implicit none
!
      class(zop_engine), intent(in) :: engine
!
      real(dp), dimension(engine%n_components), intent(in) :: electronic
      real(dp), dimension(engine%n_components), intent(in) :: nuclear
!
      integer :: component
!
      write(output%unit, '(/t3,a,a)') 'Operator: ', trim(engine%operator)
!
      write(output%unit, '(/t3, a)') 'Cart. comp.  Electronic         Nuclear             Total           '
      write(output%unit, '(t3, a)')  '--------------------------------------------------------------------'
!
      do component = 1, engine%n_components
!
         write(output%unit, '(t6, a4, 3x, f19.12, f19.12, f19.12)') engine%components(component), &
               electronic(component), nuclear(component), electronic(component) + nuclear(component)
!
      enddo
!
      write(output%unit, '(t3, a)')  '--------------------------------------------------------------------'
      flush(output%unit)
!
   end subroutine print_summary_zop_engine
!
!
end module zop_engine_class
