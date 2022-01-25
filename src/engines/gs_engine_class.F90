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
module gs_engine_class
!!
!!    Coupled cluster ground state engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use abstract_engine_class, only: abstract_engine
!
   use parameters
!
   use ccs_class,             only: ccs
   use global_in,             only: input
   use global_out,            only: output
   use memory_manager_class,  only: mem
   use timings_class,         only: timings
   use task_list_class,       only: task_list
!
   type, extends(abstract_engine) :: gs_engine
!
      character(len=200) :: multipliers_algorithm
!
      character(len=200) :: gs_algorithm
      character(len=200) :: ri_basis_set
!
      logical :: gs_restart
      logical :: multipliers_restart
      logical :: ri
!
   contains
!
      procedure :: ignite                                => ignite_gs_engine
      procedure :: run                                   => run_gs_engine
!
      procedure :: read_settings                         => read_settings_gs_engine
      procedure :: read_gs_settings                      => read_gs_settings_gs_engine
!
      procedure :: set_printables                        => set_printables_gs_engine
!
      procedure :: do_ground_state                       => do_ground_state_gs_engine
!
      procedure, private :: initialize_amplitude_updater => initialize_amplitude_updater_gs_engine
!
      procedure :: do_multipliers                        => do_multipliers_gs_engine
!
      procedure, nopass :: calculate_dipole_moment       => calculate_dipole_moment_gs_engine
!
      procedure, nopass :: calculate_quadrupole_moment   => calculate_quadrupole_moment_gs_engine
!
      procedure :: do_integral_decomposition             => do_integral_decomposition_gs_engine
      procedure, private :: do_cholesky
      procedure, private :: do_ri
!
      procedure, nopass, private :: enable_multimodel_newton
!
      procedure :: set_default_gs_algorithm &
                => set_default_gs_algorithm_gs_engine
!
      procedure :: set_default_multipliers_algorithm &
                => set_default_multipliers_algorithm_gs_engine
!
   end type gs_engine
!
!
   interface gs_engine
!
      procedure :: new_gs_engine
!
   end interface gs_engine
!
!
contains
!
!
   function new_gs_engine(wf) result(engine)
!!
!!    New GS engine
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
!     Needed for defaults and sanity checks
      class(ccs), intent(in)       :: wf
!
      type(gs_engine) :: engine
!
      call engine%set_default_gs_algorithm(wf)
!
      engine%gs_restart          = .false.
      engine%multipliers_restart = .false.
      engine%ri                  = .false.
!
      call engine%read_settings()
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_gs_engine
!
!
   subroutine set_default_gs_algorithm_gs_engine(engine, wf)
!!
!!    Set default GS algorithm
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(gs_engine), intent(inout) :: engine
!
      class(ccs), intent(in) :: wf
!
      if (wf%name_ == 'cc3') then
!
         engine%gs_algorithm = 'newton-raphson'
!
      else
!
         engine%gs_algorithm = 'diis'
!
      endif
!
   end subroutine set_default_gs_algorithm_gs_engine
!
!
   subroutine set_default_multipliers_algorithm_gs_engine(engine, wf)
!!
!!    Set default multipliers algorithm
!!    Written by Eirik F. Kjønstad, 2021
!!
      implicit none
!
      class(gs_engine), intent(inout) :: engine
!
      class(ccs), intent(in) :: wf
!
      if (wf%name_ == 'cc3') then
!
         engine%multipliers_algorithm = 'newton-raphson'
!
      else if (wf%name_ .eq. 'ccs' .or. &
          wf%name_ .eq. 'cc2' .or. &
          wf%name_ .eq. 'low memory cc2' .or. &
          wf%name_ .eq. 'mlcc2') then
!
         engine%multipliers_algorithm = 'diis'
!
      else
!
         engine%multipliers_algorithm = 'davidson'
!
      end if
!
   end subroutine set_default_multipliers_algorithm_gs_engine
!
!
   subroutine ignite_gs_engine(engine, wf)
!!
!!    Ignite
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019
!!
!!    Banner, run & cleanup
!!
      implicit none
!
      class(gs_engine) :: engine
!
      class(ccs) :: wf
!
      call engine%print_banner(wf)
      call engine%run(wf)
      call engine%print_timings(wf)
!
   end subroutine ignite_gs_engine
!
!
   subroutine run_gs_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(gs_engine) :: engine
      class(ccs)       :: wf
!
      call engine%do_integral_decomposition(wf)
!
      call engine%tasks%print_('mo preparations')
!
      call wf%mo_preparations()
!
      call engine%do_ground_state(wf)
!
   end subroutine run_gs_engine
!
!
   subroutine read_settings_gs_engine(engine)
!!
!!    Read settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(gs_engine) :: engine
!
      call engine%read_gs_settings()
!
   end subroutine read_settings_gs_engine
!
!
   subroutine read_gs_settings_gs_engine(engine)
!!
!!    Read GS engine settings
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      implicit none
!
      class(gs_engine) :: engine
!
      call input%get_keyword('algorithm', 'solver cc multipliers', engine%multipliers_algorithm)
      call input%get_keyword('algorithm', 'solver cc gs',          engine%gs_algorithm)
!
      engine%gs_restart          = input%is_keyword_present('restart', 'solver cc gs')
      engine%multipliers_restart = input%is_keyword_present('restart', 'solver cc multipliers')
!
!     global restart
      if (input%is_keyword_present('restart', 'do')) then
!
         engine%gs_restart = .true.
         engine%multipliers_restart = .true.
!
      end if
!
      if (input%is_keyword_present('ri', 'integrals')) then
!
         engine%ri = .true.
         call input%get_keyword('ri', 'integrals', engine%ri_basis_set)
!
      end if
!
      engine%plot_density = input%is_keyword_present('plot cc density', 'visualization')
!
   end subroutine read_gs_settings_gs_engine
!
!
   subroutine set_printables_gs_engine(engine)
!!
!!    Set Printables
!!    Written by Sarai D. Folkestad, May 2019
!!
!!    Should be overwritten by descendants.
!!
!
      use string_utilities, only: convert_to_uppercase
!
      implicit none
!
      class(gs_engine) :: engine
!
      engine%name_       = 'Ground state coupled cluster engine'
!
      engine%description = 'Calculates the ground state CC wavefunction | CC > = exp(T) | R >'
      engine%tag         = 'ground state'
!
!     Prepare the list of tasks
!
      engine%tasks = task_list()
!
      if (engine%ri) then
!
         call engine%tasks%add(label='ri', &
                               description='RI approximation of the electron &
                                           &repulsion integrals')
!
      else
         call engine%tasks%add(label='cholesky', &
                               description='Cholesky decomposition of the electron &
                                           &repulsion integrals')
      endif
!
      call engine%tasks%add(label='mo preparations',                             &
                            description='Preparation of MO basis and integrals')
!
      call engine%tasks%add(label='gs solver',                                &
                            description='Calculation of the ground state ('// &
                            trim((engine%gs_algorithm))//' algorithm)')
!
   end subroutine set_printables_gs_engine
!
!
   subroutine do_ground_state_gs_engine(engine, wf)
!!
!!    Do ground state
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Performs a ground state coupled cluster calculation. Determines
!!    the amplitudes, i.e. the (right) ground state | CC > = e^T | R >.
!!
!!    After this routine, the wavefunction has the cluster amplitudes
!!    stored in memory.
!!
      use diis_cc_gs_class,            only: diis_cc_gs
!
      use amplitude_updater_class,     only: amplitude_updater
!
      implicit none
!
      class(gs_engine), intent(in) :: engine
!
      class(ccs), intent(inout) :: wf
!
      type(diis_cc_gs), allocatable :: diis_solver
!
      class(amplitude_updater), allocatable :: t_updater
!
      character(len=200) :: storage
!
      call engine%tasks%print_('gs solver')
!
      if (trim(wf%name_) == 'mp2') then
!
         call wf%construct_t1_cholesky(wf%t1, wf%L_mo, wf%L_t1)
!
         call wf%calculate_energy()
!
         call wf%print_gs_summary()

!
      else
!
         storage = 'disk'
         call input%get_keyword('storage', 'solver cc gs', storage)
!
         call engine%initialize_amplitude_updater(wf, t_updater, &
                                                storage, 'solver cc gs')
!
         diis_solver = diis_cc_gs(wf, engine%gs_restart, t_updater, storage)
         call diis_solver%run(wf)
         call diis_solver%cleanup(wf)
!
      endif
!
   end subroutine do_ground_state_gs_engine
!
!
   subroutine initialize_amplitude_updater_gs_engine(engine, wf, &
                                          t_updater, global_storage, section)
!!
!!    Initialize amplitude updater
!!    Written by Eirik F. Kjønstad, 2020
!!
!!    Allocates/initializes the amplitude updater algorithm used by the ground state solver.
!!
      use abstract_jacobian_transformer_class,     only: abstract_jacobian_transformer
      use jacobian_transformer_class,              only: jacobian_transformer
      use approximate_jacobian_transformer_class,  only: approximate_jacobian_transformer
!
      use amplitude_updater_class,      only: amplitude_updater
!
      use quasi_newton_updater_class,   only: quasi_newton_updater
      use newton_raphson_updater_class, only: newton_raphson_updater
!
      use citation_printer_class, only : eT_citations
!
      implicit none
!
      class(gs_engine), intent(in) :: engine
!
      class(ccs), intent(in) :: wf
!
      character(len=*), intent(in) :: global_storage
      character(len=*), intent(in) :: section
!
      logical :: scale_amplitudes, scale_residual
!
      character(len=200) :: side, algorithm
!
      class(amplitude_updater), allocatable :: t_updater
!
      class(abstract_jacobian_transformer), allocatable :: transformer
!
      real(dp) :: relative_threshold
!
      integer :: max_iterations
!
      logical :: records_in_memory, multipliers
      character(len=200) :: storage
!
      if (section == 'solver cc multipliers') then
!
         multipliers = .true.
         algorithm = engine%multipliers_algorithm
!
      else
!
         multipliers = .false.
         algorithm = engine%gs_algorithm
!
      end if
!
      if (algorithm == 'diis') then
!
         if (multipliers) then
!
            scale_amplitudes = .false.
            scale_residual = .false.
!
         else
!
            scale_amplitudes = .true.
            scale_residual = .true.
!
         end if
!
         t_updater = quasi_newton_updater(n_amplitudes     = wf%n_gs_amplitudes, &
                                          scale_amplitudes = scale_amplitudes,   &
                                          scale_residual   = scale_residual)
!
      elseif (algorithm == 'newton-raphson') then
!
         if (trim(wf%name_) == 'cc2') &
            call output%error_msg('Newton-Raphson not implemented for CC2')
!
!        Set defaults for microiterations, and read if non-standard settings
!
         relative_threshold = 1.0d-2
         max_iterations     = 100
         storage            = 'disk'
!
         call input%get_keyword('rel micro threshold',     section, relative_threshold)
         call input%get_keyword('max micro iterations',    section, max_iterations)
         call input%get_keyword('micro iteration storage', section, storage)
!
         records_in_memory = .false.
!
         if (trim(global_storage) == 'memory') records_in_memory = .true.
         if (trim(storage) == 'memory') records_in_memory = .true.
!
         if (multipliers) then
!
            scale_amplitudes = .false.
            scale_residual = .false.
            side = 'left'
!
         else
!
            scale_amplitudes = .true.
            scale_residual = .false.
            side = 'right'
!
         end if
!
!        Determine which transformation to use - exact or approximate Jacobian -
!        and then construct the t_updater
!
         if (engine%enable_multimodel_newton(wf, section)) then
!
            transformer = approximate_jacobian_transformer(side)
!
            call eT_citations%add('Multimodel Newton algorithm')
!
         else
!
            transformer = jacobian_transformer(side)
!
         endif
!
         t_updater = newton_raphson_updater(n_amplitudes       = wf%n_gs_amplitudes, &
                                            scale_amplitudes   = scale_amplitudes,   &
                                            scale_residual     = scale_residual,     &
                                            relative_threshold = relative_threshold, &
                                            records_in_memory  = records_in_memory,  &
                                            max_iterations     = max_iterations,     &
                                            transformer        = transformer)
!
      else
!
         call output%error_msg('(a0) is not a valid ground state algorithm.', &
                                  chars=[engine%gs_algorithm])
!
      endif
!
   end subroutine initialize_amplitude_updater_gs_engine
!
!
   function enable_multimodel_newton(wf, section) result(enable)
!!
!!    Enable multimodel Newton?
!!    Written by Eirik F. Kjønstad, 2021
!!
!!    Default is yes for CC3 and no for other CC methods
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      character(len=*), intent(in) :: section
!
      logical :: enable
!
      character(len=200) :: multimodel_newton
!
      enable = .false.
      if (wf%name_ == 'cc3') enable = .true.
!
      multimodel_newton = 'default'
      call input%get_keyword('multimodel newton', section, multimodel_newton)
!
      if (multimodel_newton == 'on') enable = .true.
      if (multimodel_newton == 'off') enable = .false.
!
   end function enable_multimodel_newton
!
!
   subroutine do_multipliers_gs_engine(engine, wf)
!!
!!    Do multipliers
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Performs a (left) ground state coupled cluster calculation. Determines
!!    the amplitudes, i.e. the (left) ground state < Λ |  = (< R | + sum_mu tbar_mu < mu |) e^-T.
!!
!!    After this routine, the wavefunction has the cluster amplitudes
!!    and the multipliers stored in memory.
!!
!!    Modified by Regina Matveeva, Sep 2021
!!
!!    Replaced the davidson_cc_multipliers_solver with general_linear_davidson_solver
!!    in combination with cc_multipliers_solver_factory
!!
!
      use general_linear_davidson_solver_class, only: general_linear_davidson_solver
      use cc_multipliers_solver_factory_class,  only: cc_multipliers_solver_factory
      use diis_cc_multipliers_class,            only: diis_cc_multipliers
      use string_utilities,                     only: convert_to_uppercase
      use amplitude_updater_class,              only: amplitude_updater
!
      implicit none
!
      class(gs_engine), intent(inout) :: engine
      class(ccs), intent(inout) :: wf
!
      class(general_linear_davidson_solver),     allocatable  :: solver
      type(cc_multipliers_solver_factory)              :: solver_factory
      type(diis_cc_multipliers),          allocatable  :: diis_solver
!
      type(timings), allocatable :: timer
!
      real(dp), dimension(:), allocatable :: frequencies
!
      class(amplitude_updater), allocatable :: tbar_updater
!
      character(len=200) :: storage
!
      call wf%construct_fock(task = 'multipliers')
!
      call engine%tasks%print_('multipliers solver')
!
      timer = timings(trim(convert_to_uppercase(wf%name_)) // ' multipliers')
      call timer%turn_on()
!
      call wf%initialize_multipliers()
!
      if (trim(engine%multipliers_algorithm) == 'davidson') then
!
         if (trim(wf%name_) == 'cc2' .or.            &
             trim(wf%name_) == 'low memory cc2' .or. &
             trim(wf%name_) == 'cc3') then
!
            call output%error_msg('Davidson not implemented for (a0)', &
                                  chars=[wf%name_])
!
         end if
!
         call wf%print_banner_davidson_cc_multipliers()
!
         call solver_factory%create(wf, solver)
!
         call mem%alloc(frequencies, 1)
         frequencies = zero
!
         call solver%run(frequencies)
!
         call mem%dealloc(frequencies, 1)
!
         call wf%cc_multipliers_summary()
!
      elseif (trim(engine%multipliers_algorithm) == 'diis' .or. &
              trim(engine%multipliers_algorithm) == 'newton-raphson') then
!
         storage = 'disk'
         call input%get_keyword('storage', 'solver cc multipliers', storage)
!
         call engine%initialize_amplitude_updater(wf, tbar_updater, &
                                                storage, 'solver cc multipliers')
!
         diis_solver = diis_cc_multipliers(wf, engine%multipliers_restart, tbar_updater)
         call diis_solver%run(wf)
         call diis_solver%cleanup(wf)
!
      else
!
         call output%error_msg('Could not recognize multipliers algorithm (a0)', &
                               chars=[trim(engine%multipliers_algorithm)])
!
      endif
!
      call timer%turn_off()
!
      call output%printf('m', '- Finished solving the (a0) multipliers equations', &
                         fs='(/t3,a)', chars=[trim(convert_to_uppercase(wf%name_))])
!
      call output%printf('m', 'Total wall time (sec):  (f20.5)', &
                         reals=[timer%get_elapsed_time('wall')], fs='(/t6,a)')
!
      call output%printf('m', 'Total cpu time (sec):   (f20.5)', &
                         reals=[timer%get_elapsed_time('cpu')], fs='(t6,a)')
!
   end subroutine do_multipliers_gs_engine
!
!
   subroutine calculate_dipole_moment_gs_engine(wf, electronic, nuclear, total)
!!
!!    Calculate dipole moment
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Determines < mu > = < Λ | mu | CC > for the one-electron operator mu.
!!    The routine assumes that the multipliers and amplitudes have been determined
!!    and that the one-electron density is calculated and in memory.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(3), intent(out) :: electronic
      real(dp), dimension(3), intent(out) :: nuclear
      real(dp), dimension(3), intent(out) :: total
!
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: mu_pqk
!
!     Get the integrals mu_pqk for components k = 1, 2, 3 in the T1-transformed basis
!
      call mem%alloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
      call wf%get_t1_oei('dipole', mu_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 3
!
         electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k), wf%density)
         if (wf%exists_frozen_fock_terms) electronic(k) = electronic(k) + wf%frozen_dipole(k)
!
      enddo
!
      call mem%dealloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
!
!     Get nuclear expectation value contribution, then sum the two
!
      nuclear = wf%get_nuclear_dipole()
!
      total = electronic + nuclear
!
   end subroutine calculate_dipole_moment_gs_engine
!
!
   subroutine calculate_quadrupole_moment_gs_engine(wf, electronic, nuclear, total)
!!
!!    Calculate quadrupole moment
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!!    Determines < q > = < Λ | q | CC > for the one-electron operator q.
!!    The routine assumes that the multipliers and amplitudes have been determined
!!    and that the one-electron density is calculated and in memory.
!!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(6), intent(out) :: electronic
      real(dp), dimension(6), intent(out) :: nuclear
      real(dp), dimension(6), intent(out) :: total
!
      integer :: k
!
      real(dp), dimension(:,:,:), allocatable :: q_pqk
!
!     Get the integrals q_pqk for components k = 1, 2, ..., 6 in the T1-transformed basis
!
      call mem%alloc(q_pqk, wf%n_mo, wf%n_mo, 6)
      call wf%get_t1_oei('quadrupole', q_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 6
!
         electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k), wf%density)
         if (wf%exists_frozen_fock_terms) electronic(k) = electronic(k) + wf%frozen_quadrupole(k)
!
      enddo
!
      call mem%dealloc(q_pqk, wf%n_mo, wf%n_mo, 6)
!
!     Get nuclear expectation value contribution, then sum the two
!
      nuclear = wf%get_nuclear_quadrupole()
!
      total = electronic + nuclear
!
   end subroutine calculate_quadrupole_moment_gs_engine
!
!
   subroutine do_integral_decomposition_gs_engine(engine, wf)
!!
!!    Do integral decomposition
!!    Written by Sarai D. Folkestad, 2021
!!
      implicit none
!
      class(gs_engine), intent(in)     :: engine
      class(ccs),       intent(inout)  :: wf
!
      if (engine%ri) then
!
         call engine%do_ri(wf)
!
      else
!
         call engine%do_cholesky(wf)
!
      endif
!
      call wf%L_t1%set_equal_to(wf%L_mo)
!
   end subroutine do_integral_decomposition_gs_engine
!
!
   subroutine do_ri(engine, wf)
!!
!!    Do RI
!!    Written by Sarai D. Folkestad, Feb 2021
!!
      use eri_ri_class, only: eri_ri
!
      implicit none
!
      class(gs_engine), intent(in)     :: engine
      class(ccs),       intent(inout)  :: wf
!
      type(eri_ri), allocatable :: ri_solver
!
      call engine%tasks%print_('ri')
!
      ri_solver = eri_ri(engine%ri_basis_set, wf%ao%get_libint_epsilon())
      call ri_solver%initialize()
!
      call ri_solver%run(wf%ao)
!
      call wf%integral_preparations(ri_solver%get_n_J())
!
      call ri_solver%construct_cholesky_mo_vectors(wf%ao, wf%ao%n, wf%n_mo, &
                                                   wf%orbital_coefficients, wf%L_mo)
!
      call ri_solver%cleanup()
!
   end subroutine do_ri
!
!
   subroutine do_cholesky(engine, wf)
!!
!!    Do Cholesky
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Feb 2020
!!
!!    Performs Cholesky decomposition of the electron repulsion integral matrix.
!!
!!    For reduced space coupled cluster calculations (frozen HF or MLHF)
!!    the MO screening can be used to reduce the number of Cholesky vectors,
!!    as the accuracy of the integrals in the active MO basis,
!!    rather than the AO basis, is targeted. For details, see
!!    Folkestad, S. D., Kjønstad, E. F., and Koch, H., JCP, 150(19), 194112 (2019).
!
      use eri_cd_class, only: eri_cd
!
      implicit none
!
      class(gs_engine), intent(in) :: engine
!
      class(ccs) :: wf
!
      type(eri_cd), allocatable :: eri_cholesky_solver
!
      logical :: do_MO_screening
!
      real(dp), dimension(:,:), allocatable, target :: screening_vector
!
      call engine%tasks%print_('cholesky')
!
      do_MO_screening = input%is_keyword_present('mo screening', 'solver cholesky')
!
      eri_cholesky_solver = eri_cd(wf%ao)
!
      if (do_MO_screening) then
!
         call output%printf('m', 'Using the MO screening for the Cholesky decomposition', &
                        fs='(/t3,a)')
!
         call mem%alloc(screening_vector, wf%ao%n, wf%ao%n)
         call wf%construct_MO_screening_for_cd(screening_vector)
!
         call eri_cholesky_solver%run(wf%ao, screening_vector)   ! Do the Cholesky decomposition
                                                                 ! using the MO screening
!
         call mem%dealloc(screening_vector, wf%ao%n, wf%ao%n)
!
      else
!
         call eri_cholesky_solver%run(wf%ao) ! Do the Cholesky decomposition
!
      endif
!
      call eri_cholesky_solver%diagonal_test(wf%ao)  ! Determine the largest
                                                     ! deviation in the ERI matrix
!
      call wf%integral_preparations(eri_cholesky_solver%get_n_cholesky())
!
      call eri_cholesky_solver%construct_cholesky_mo_vectors(wf%ao, wf%ao%n, wf%n_mo, &
                                                   wf%orbital_coefficients, wf%L_mo)
!
      call eri_cholesky_solver%cleanup()
!
   end subroutine do_cholesky
!
!
end module gs_engine_class
