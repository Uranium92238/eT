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
module mean_value_engine_class
!!
!!    Expectation value coupled cluster engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Calculates expectation values < Λ | A | CC > for the operator A,
!!    where A is the dipole or quadrupole operator.
!!
   use kinds
   use global_in,     only: input
   use global_out,    only: output
   use timings_class, only: timings
!
   use gs_engine_class, only: gs_engine
   use ccs_class,       only: ccs
!
   use task_list_class, only: task_list
!
   type, extends(gs_engine) :: mean_value_engine
!
   contains
!
      procedure :: run => run_mean_value_engine
!
      procedure :: read_settings => read_settings_mean_value_engine
!
      procedure :: read_cc_mean_value_settings  => read_cc_mean_value_settings_mean_value_engine
!
      procedure :: calculate_mean_values  => calculate_mean_values_mean_value_engine
!
      procedure :: set_printables  => set_printables_mean_value_engine
!
   end type mean_value_engine
!
!
   interface mean_value_engine
!
      procedure :: new_mean_value_engine
!
   end interface mean_value_engine
!
!
contains
!
   function new_mean_value_engine(wf) result(engine)
!!
!!    New expectation value engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
!     Needed for defaults and sanity checks
      class(ccs), intent(in) :: wf
!
      type(mean_value_engine) :: engine
!
      if (wf%name_ .eq. 'ccsd(t)' .or. &
          wf%name_ .eq. 'low memory cc2' .or. &
          wf%name_ .eq. 'mlccsd' .or. &
          wf%name_ .eq. 'mp2') then
!
         call output%error_msg("Mean value calculations not implemented for (a0)", &
                               chars=[wf%name_])
!
      elseif (wf%name_ .eq. 'mlcc2') then
!
         call output%warning_msg("Mean value calculations not recomended for (a0)!", &
                               chars=[wf%name_])
!
      end if
!
      call engine%set_default_gs_algorithm(wf)
      call engine%set_default_multipliers_algorithm(wf)
!
      engine%dipole               = .false.
      engine%quadrupole           = .false.
      engine%plot_density         = .false.
!
      engine%gs_restart           = .false.
      engine%multipliers_restart  = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_mean_value_engine
!
!
   subroutine read_settings_mean_value_engine(engine)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(mean_value_engine) :: engine
!
      call engine%read_gs_settings()
      call engine%read_cc_mean_value_settings()
!
   end subroutine read_settings_mean_value_engine
!
!
   subroutine read_cc_mean_value_settings_mean_value_engine(engine)
!!
!!    Read expectation value settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(mean_value_engine) :: engine
!
      if (input%is_keyword_present('dipole','cc mean value')) engine%dipole = .true.
      if (input%is_keyword_present('quadrupole','cc mean value')) engine%quadrupole = .true.
!
   end subroutine read_cc_mean_value_settings_mean_value_engine
!
!
   subroutine run_mean_value_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      use visualization_class,  only: visualization
      use memory_manager_class, only: mem
!
      implicit none
!
      class(ccs)               :: wf
      class(mean_value_engine) :: engine
!
      type(visualization), allocatable :: visualizer
!
      real(dp), dimension(:,:), allocatable :: c_D_ct
!
!     Determine ground state | CC >
!
      call engine%gs_engine%run(wf)
!
!     Determine multipliers < Λ |
!
      call engine%do_multipliers(wf)
!
!     Compute the one-electron density
!
      call wf%prepare_for_properties()
      call wf%initialize_gs_density()
      call wf%construct_gs_density()
!
      call engine%calculate_mean_values(wf)
!
      if (engine%plot_density) then
!
!        Transform the density to AO basis and plot
!
         visualizer = visualization(wf%ao)
         call visualizer%initialize(wf%ao)
!
         call mem%alloc(c_D_ct, wf%ao%n, wf%ao%n)
!
         call wf%add_t1_terms_and_transform(wf%density, c_D_ct)
         call visualizer%plot_density(wf%ao, c_D_ct, 'cc_gs_density')
!
         call visualizer%cleanup()
!
         call mem%dealloc(c_D_ct, wf%ao%n, wf%ao%n)
!
      endif
!
      call wf%destruct_gs_density()
!
   end subroutine run_mean_value_engine
!
!
   subroutine set_printables_mean_value_engine(engine)
!!
!!    Set printables
!!    Written by sarai D. Folkestad, May 2019
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(mean_value_engine) :: engine
!
      engine%name_ = 'Mean value coupled cluster engine'
!
      engine%tag   = 'mean value'
!
!     Prepare the list of tasks
!
      engine%tasks = task_list()
!
      call engine%tasks%add(label='cholesky', &
                            description='Cholesky decomposition of the electron &
                                         &repulsion integrals')
!
      call engine%tasks%add(label='mo preparations',                             &
                            description='Preparation of MO basis and integrals')
!
      call engine%tasks%add(label='gs solver',                                &
                            description='Calculation of the ground state ('// &
                           trim((engine%gs_algorithm))//' algorithm)')
!
      call engine%tasks%add(label='multipliers solver',                       &
                            description='Calculation of the multipliers ('    &
                            //trim((engine%multipliers_algorithm))&
                            //' algorithm)')
!
      call engine%tasks%add(label='expectation value',                        &
                            description='Calculation of ground state properties')
!
      if (engine%plot_density) &
         call engine%tasks%add(label='plotting', description='Plot ground state density')
!
      engine%description  = 'Calculates the time-independent expectation value of&
                            & one-electron operators A, < A > = < Λ | A | CC >.'
!
   end subroutine set_printables_mean_value_engine
!
!
   subroutine calculate_mean_values_mean_value_engine(engine, wf)
!!
!!    Calculate expectation values
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(mean_value_engine), intent(in) :: engine
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
      call engine%tasks%print_('expectation value')
!
      if(engine%dipole) then 
!
         call engine%calculate_dipole_moment(wf, mu_electronic, mu_nuclear, mu_total)
!
         allocate(components(3))
!
         components = (/'x   ',&
                        'y   ',&
                        'z   '/)
!
         call engine%print_operator('dipole moment', mu_electronic, mu_nuclear, mu_total, &
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
         call engine%print_operator('quadrupole moment (with trace)', q_electronic, q_nuclear, q_total, &
                                    components, 6)
!
         call engine%remove_trace(q_electronic)
         call engine%remove_trace(q_nuclear)
!
         q_total = q_electronic + q_nuclear
!
         call output%printf('m', 'The traceless quadrupole is calculated as:', fs='(/t6,a)')
         call output%printf('m', 'Q_ij = 1/2[3*q_ij - tr(q)*delta_ij]',fs='(/t9,a)')
         call output%printf('m', 'where q_ij is the non-traceless matrix',fs='(/t6,a)')
!
         call engine%print_operator('traceless quadrupole moment', q_electronic, q_nuclear, q_total, &
                                    components, 6)
!
         deallocate(components)
!
      endif
!
   end subroutine calculate_mean_values_mean_value_engine
!
!
end module mean_value_engine_class
