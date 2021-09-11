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
module es_engine_class
!!
!!    Coupled cluster ground state engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
!
   use kinds
   use global_in,       only: input
   use global_out,      only: output
   use timings_class,   only: timings
   use task_list_class, only: task_list
!
   use gs_engine_class, only: gs_engine
   use ccs_class,       only: ccs
!
   type, extends(gs_engine) :: es_engine
!
      character(len=200) :: es_algorithm
      character(len=200) :: es_type
      character(len=200) :: es_transformation
!
      logical :: es_restart
!
       logical :: plot_ntos, plot_cntos
!
   contains
!
      procedure :: run                       => run_es_engine
!
      procedure :: read_settings             => read_settings_es_engine
      procedure :: read_es_settings          => read_es_settings_es_engine
!
      procedure :: do_excited_state          => do_excited_state_es_engine
!
      procedure :: set_printables            => set_printables_es_engine
!
!
      procedure :: do_nto_visualization      => do_nto_visualization_es_engine
!
      procedure, nopass :: get_thresholds &
                        => get_thresholds_es_engine
!
      procedure, nopass :: plot_orbitals &
                        => plot_orbitals_es_engine
!
      procedure :: set_default_es_algorithm &
                => set_default_es_algorithm_es_engine
!
   end type es_engine
!
!
   interface es_engine
!
      procedure :: new_es_engine
!
   end interface es_engine
!
!
contains
!
!
   function new_es_engine(wf) result(engine)
!!
!!    New ES engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
!     Needed for defaults and sanity checks
      class(ccs), intent(in)       :: wf
!
      type(es_engine) :: engine
!
!     Set standards and then read if nonstandard
!
      call engine%set_default_gs_algorithm(wf)
      call engine%set_default_multipliers_algorithm(wf)
      call engine%set_default_es_algorithm(wf)
!
      engine%es_type                = 'valence'
      engine%es_transformation      = 'right'
!
      engine%gs_restart            = .false.
      engine%multipliers_restart   = .false.
      engine%es_restart            = .false.
      engine%plot_cntos            = .false.
      engine%plot_ntos             = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_es_engine
!
!
   subroutine set_default_es_algorithm_es_engine(engine, wf)
!!
!!    Set default ES algorithm
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none 
!
      class(es_engine), intent(inout) :: engine 
!
      class(ccs), intent(in) :: wf 
!
      if (wf%name_ .eq. 'ccsd(t)' .or. &
          wf%name_ .eq. 'mp2') then
         call output%error_msg("Excited states not implemented for (a0)", &
                               chars=[wf%name_])
      end if
!
      if (wf%name_ .eq. 'cc3' .or. &
          wf%name_ .eq. 'low memory cc2') then
!
         engine%es_algorithm = 'non-linear davidson'
!
      else
!
         engine%es_algorithm = 'davidson'
!
      end if
!
   end subroutine set_default_es_algorithm_es_engine
!
!
   subroutine read_settings_es_engine(engine)
!!
!!    Read es settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(es_engine) :: engine
!
      call engine%read_gs_settings()
      call engine%read_es_settings()
!
   end subroutine read_settings_es_engine
!
!
   subroutine read_es_settings_es_engine(engine)
!!
!!    Read es settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019
!!
      implicit none
!
      class(es_engine) :: engine
!
      call input%get_keyword('algorithm', 'solver cc es', engine%es_algorithm)
!
      if (input%is_keyword_present('core excitation', 'solver cc es') .and. .not. &
          input%is_keyword_present('ionization', 'solver cc es')) engine%es_type = 'core'
!
      if (input%is_keyword_present('ionization', 'solver cc es') .and. .not. &
         input%is_keyword_present('core excitation', 'solver cc es')) engine%es_type = 'ionize'
!
      if (input%is_keyword_present('ionization', 'solver cc es') .and.    &
         input%is_keyword_present('core excitation', 'solver cc es'))     &
            call output%error_msg('XPS still not implemented.')
!
      if (input%is_keyword_present('left eigenvectors', 'solver cc es')) &
                                                engine%es_transformation = 'left'
!
      if (input%is_keyword_present('right eigenvectors', 'solver cc es')) then
         if (engine%es_transformation == 'left') then
!
            engine%es_transformation = 'both'
!
         else
!
            engine%es_transformation = 'right'
!
         end if
      end if
!
      engine%es_restart = input%is_keyword_present('restart', 'solver cc es')
!
!     global restart
      if (input%is_keyword_present('restart', 'do')) engine%es_restart = .true.
!
      engine%plot_ntos  = input%is_keyword_present('plot ntos','visualization')
      engine%plot_cntos = input%is_keyword_present('plot cntos','visualization')
!
   end subroutine read_es_settings_es_engine
!
!
   subroutine run_es_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
!
      real(dp) :: energy_threshold, residual_threshold
!
!     Determine ground state
!
      call engine%gs_engine%run(wf)
!
!     Determine excited states
!
      if (engine%es_transformation .eq. 'both') then
!
         call engine%do_excited_state(wf, 'right', engine%es_restart)
         call engine%do_excited_state(wf, 'left', restart = .true.)
!
      else
!
         call engine%do_excited_state(wf, engine%es_transformation, engine%es_restart)
!
         if (engine%plot_ntos) call engine%do_nto_visualization(wf, 'nto')
!
         if (engine%plot_cntos) call engine%do_nto_visualization(wf, 'cnto')
!
      end if
!
      if (engine%es_algorithm .eq. 'diis') then
!
         call engine%get_thresholds(energy_threshold, residual_threshold)
         call wf%remove_parallel_states(residual_threshold, engine%es_transformation)
!
      end if
!
   end subroutine run_es_engine
!
!
   subroutine do_excited_state_es_engine(engine, wf, transformation, restart)
!!
!!    Do excited state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Solves the excited state (valence or cvs) using
!!    either a DIIS or Davidson solver
!!
!!    or calculates the excitation energies and oscillator strengths
!!    using the asymmetric Lanczos algorithm.
!!
!!    Modified by Torsha Moitra, S. Coriani and Sarai D. Folkestad, Sep-Nov 2019
!!
!!       Added the asymmetric Lanczos solver
!!
      use abstract_cc_es_class, only: abstract_cc_es
      use davidson_cc_es_class, only: davidson_cc_es
      use nonlinear_davidson_cc_es_class, only: nonlinear_davidson_cc_es
      use diis_cc_es_class, only: diis_cc_es
      use asymmetric_lanczos_cc_es_class, only: asymmetric_lanczos_cc_es
!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
!
      character(len=*), intent(in) :: transformation
!
      logical, intent(in) :: restart
!
      class(asymmetric_lanczos_cc_es), allocatable :: cc_es_solver_asymmetric_lanczos

      class(abstract_cc_es), allocatable :: cc_es_solver
!
      call engine%tasks%print_('es solver', &
            append_string='Calculating ' // trim(transformation) //' vectors', append_fs='(t6,a)')
!
!     Prepare for excited state
!
      if (engine%es_algorithm == 'asymmetric lanczos') then
!
         call engine%do_multipliers(wf)
         cc_es_solver_asymmetric_lanczos = asymmetric_lanczos_cc_es(wf)
         call cc_es_solver_asymmetric_lanczos%run(wf)
         call cc_es_solver_asymmetric_lanczos%cleanup(wf)
!
      else
!
         if (engine%es_algorithm == 'diis') then
!
            cc_es_solver = diis_cc_es(transformation, wf, restart)
!
         elseif (engine%es_algorithm == 'davidson') then
!
            if (trim(wf%name_) == 'low memory cc2' .or. trim(wf%name_) == 'cc3') then
!
                call output%error_msg('Davidson not implemented for CC3 and lowmem CC2')
!
            end if
!
            cc_es_solver = davidson_cc_es(transformation, wf, restart)
!
         elseif (engine%es_algorithm == 'non-linear davidson') then
!
            cc_es_solver = nonlinear_davidson_cc_es(transformation, wf, restart)
!
         else
!
            call output%error_msg('Could not start excited state solver. It may be that the &
                                    &algorithm is not implemented for the method specified.')
         endif
!
!        Calculate the terms of the Fock matrix that was not needed for ground state
!
         call wf%construct_fock(task = 'es')
!
         call cc_es_solver%run(wf)
         call cc_es_solver%cleanup(wf)
!
      endif
!
   end subroutine do_excited_state_es_engine
!
!
   subroutine set_printables_es_engine(engine)
!!
!!    Set printables
!!    Written by sarai D. Folkestad, May 2019
!!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(es_engine) :: engine
!
      engine%name_  = 'Excited state coupled cluster engine'
!
      engine%tag = 'excited state'
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
                           description='Calculation of the ground state ('//  &
                           trim(engine%gs_algorithm)//' algorithm)')
!
      if (engine%es_algorithm == 'asymmetric lanczos') then
!
         call engine%tasks%add(label='multipliers solver',                    &
                           description='Calculation of the multipliers ('     &
                           //trim(engine%multipliers_algorithm)&
                           //' algorithm)')
!
         call engine%tasks%add(label='es solver',                             &
                           description='Calculation of the excited state ('// &
                           trim(engine%es_algorithm)//' algorithm)')
!
      else
!
         call engine%tasks%add(label='es solver',                              &
                           description='Calculation of the excited state ('//  &
                           trim((engine%es_algorithm))//' algorithm)')
!
      endif
!
!
      if (engine%plot_ntos .or. engine%plot_cntos) &
         call engine%tasks%add(label='plotting',   &
                            description='Plotting NTOs/CNTOs')
!
      engine%description = 'Calculates the coupled cluster excitation &
                           &vectors and excitation energies'
!
   end subroutine set_printables_es_engine
!
!
   subroutine do_nto_visualization_es_engine(engine, wf, tag)
!!
!!    Do NTO visualization
!!    Written by Sarai D. Folkestad, May 2020
!!
!!    Makes .plt files for NTOs and or CNTOs
!!
      use visualization_class, only : visualization
      use memory_manager_class, only : mem
      use parameters
!
      implicit none
!
      class(es_engine), intent(in) :: engine
      class(ccs), intent(inout) :: wf
!
      character(len=*), intent(in) :: tag
!
      type(visualization), allocatable :: plotter
!
      integer :: state, n_significant_v, n_significant_o
!
      real(dp), dimension(:,:), allocatable :: orbitals
!
      real(dp) :: threshold
!
      threshold = 0.1d0
!
      call input%get_keyword('nto threshold', 'visualization', threshold)
!
      if ((threshold .gt. 1) .or. (threshold .lt. 0)) &
            call output%error_msg('illegal threshold given for NTO/CNTO plotting')
!
      call engine%tasks%print_('plotting')
!
      plotter = visualization(wf%ao)

      call mem%alloc(orbitals, wf%ao%n, wf%n_mo)
!
      do state = 1, wf%n_singlet_states
!
         call wf%construct_ntos_or_cntos(orbitals, state, n_significant_v, n_significant_o, &
                                         engine%es_transformation, tag, threshold)
!
         call engine%plot_orbitals(wf, n_significant_o, n_significant_v, &
                                   orbitals, tag, state, plotter)
!
      enddo
!
      call mem%dealloc(orbitals, wf%ao%n, wf%n_mo)
!
   end subroutine do_nto_visualization_es_engine
!
!
   subroutine get_thresholds_es_engine(energy_threshold, residual_threshold)
!!
!!    Get thresholds
!!    Written by Alexander C. Paul, Oct 2019
!!
!!    Get thresholds from input to perform checks for parallel states and
!!    to check that the left and right states are consistent
!!
      implicit none
!
      real(dp), intent(out) :: energy_threshold, residual_threshold
!
      energy_threshold = 1.0d-3
      residual_threshold = energy_threshold
!
      if (input%is_keyword_present('energy threshold', 'solver cc es') .and. &
          input%is_keyword_present('residual threshold', 'solver cc es')) then
!
        call input%get_keyword('energy threshold', 'solver cc es', energy_threshold)
        call input%get_keyword('residual threshold', 'solver cc es', residual_threshold)
!
      else if (input%is_keyword_present('residual threshold', 'solver cc es')) then
!
        call input%get_keyword('residual threshold', 'solver cc es', residual_threshold)
        energy_threshold = residual_threshold
!
      else if (input%is_keyword_present('energy threshold', 'solver cc es')) then
!
         call input%get_keyword('energy threshold', 'solver cc es', energy_threshold)
         residual_threshold = energy_threshold
!
      endif
!
   end subroutine get_thresholds_es_engine
!
!
   subroutine plot_orbitals_es_engine(wf, n_sig_o, n_sig_v, orbitals, tag, state, plotter)
!!
!!    Plot orbitals
!!    Written by Sarai D. Folkestad, May 2020
!!
!
      use visualization_class, only : visualization
!
      implicit none
!
      class(ccs),                            intent(in)    :: wf
!
      integer,                               intent(in)    :: n_sig_o, n_sig_v, state
      real(dp), dimension(wf%ao%n, wf%n_mo), intent(in)    :: orbitals
      character(len=*),                      intent(in)    :: tag
      type(visualization),                   intent(inout) :: plotter
!
      character(len=200), dimension(:), allocatable :: file_tags
!
      integer :: i
!
      allocate(file_tags(max(n_sig_o, n_sig_v)))
!
      do i = 1, n_sig_o
!
         write(file_tags(i), '(a,i3.3,a,a,a,i3.3)') 'state_', state, '_', tag, '_o_', i
!
      enddo
!
      call plotter%plot_orbitals(wf%ao, orbitals(:, 1:n_sig_o), n_sig_o, file_tags)
!
      do i = 1, n_sig_v
!
         write(file_tags(i), '(a,i3.3,a,a,a,i3.3)') 'state_', state, '_', tag, '_v_', i
!
      enddo
!
      call plotter%plot_orbitals(wf%ao, orbitals(:, wf%n_o+1:n_sig_v), n_sig_v, file_tags)
!
      deallocate(file_tags)
!
   end subroutine plot_orbitals_es_engine
!
!
end module es_engine_class
