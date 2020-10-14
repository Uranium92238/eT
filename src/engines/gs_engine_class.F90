!
!
!  eT - a coupled cluster program
!  Copyright (C) 2016-2020 the authors of eT
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
!
      logical :: restart
      logical :: gs_restart
      logical :: multipliers_restart
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
      procedure :: do_multipliers                        => do_multipliers_gs_engine
!
      procedure, nopass :: calculate_dipole_moment       => calculate_dipole_moment_gs_engine
!
      procedure, nopass :: calculate_quadrupole_moment   => calculate_quadrupole_moment_gs_engine
!
      procedure :: restart_handling                      => restart_handling_gs_engine
!
      procedure, nopass :: do_cholesky                   => do_cholesky_gs_engine
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
      engine%gs_algorithm          = 'diis'
!
      if (wf%name_ .eq. 'cc2' .or. &
          wf%name_ .eq. 'cc3' .or. &
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
      engine%gs_restart            = .false.
      engine%multipliers_restart   = .false.
!
      call engine%read_settings()
!
      engine%restart = engine%gs_restart .or. engine%multipliers_restart
!
      call engine%set_printables()
!
      engine%timer = timings(trim(engine%name_))
      call engine%timer%turn_on()
!
   end function new_gs_engine
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
      call engine%tasks%print_('cholesky')
!
      call engine%do_cholesky(wf)
!
      call engine%tasks%print_('mo preparations')
!
      call wf%mo_preparations() 
!
      call engine%restart_handling(wf)
!
!     Determine ground state
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
      call input%get_keyword_in_section('algorithm', 'solver cc multipliers', &
                                                            engine%multipliers_algorithm)
!
      call input%get_keyword_in_section('algorithm', 'solver cc gs', engine%gs_algorithm )
!
      engine%gs_restart = input%requested_keyword_in_section('restart', 'solver cc gs')
!
      engine%multipliers_restart = &
               input%requested_keyword_in_section('restart', 'solver cc multipliers')
!
      engine%plot_density = &
               input%requested_keyword_in_section('plot cc density', 'visualization')
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
      use newton_raphson_cc_gs_class,  only: newton_raphson_cc_gs
      use string_utilities,            only: convert_to_uppercase
!
      implicit none
!
      class(gs_engine), intent(in) :: engine
!
      class(ccs), intent(inout) :: wf
!
      type(diis_cc_gs), allocatable :: diis_solver
      type(newton_raphson_cc_gs), allocatable :: newton_raphson_solver
!
      call engine%tasks%print_('gs solver')
!
      if (trim(wf%name_) == 'mp2') then
!
         call wf%eri%set_t1_to_mo()
!
         call wf%calculate_energy()
!
         call output%printf('m', ':: Summary of (a0) wavefunction energetics (a.u.)', &
                            chars=[convert_to_uppercase(wf%name_)], fs='(/t3,a)')
!
         call output%printf('m', 'HF energy:                (f19.12)', &
                            reals=[wf%hf_energy], fs='(/t6,a)')
!
         call output%printf('m', 'MP2 correction:           (f19.12)', &
                            reals=[ (wf%energy - wf%hf_energy) ], fs='(t6,a)')
!
         call output%printf('m', 'MP2 energy:               (f19.12)', &
                            reals=[wf%energy], fs='(t6,a)')
!
      else
!
         if (trim(engine%gs_algorithm) == 'diis') then
!
            diis_solver = diis_cc_gs(wf, engine%gs_restart)
            call diis_solver%run(wf)
            call diis_solver%cleanup(wf)
!
         elseif (trim(engine%gs_algorithm) == 'newton-raphson') then 
!
               if (trim(wf%name_) == 'cc2') then
!
                  call output%error_msg('Newton-Raphson not implemented for CC2')
!
               end if
!
            newton_raphson_solver = newton_raphson_cc_gs(wf, engine%gs_restart)
            call newton_raphson_solver%run(wf)
            call newton_raphson_solver%cleanup(wf)
!
         else 
!
            call output%error_msg('(a0) is not a valid ground state algorithm.', &
                                  chars=[engine%gs_algorithm])
!
         endif 
!
      endif
!
   end subroutine do_ground_state_gs_engine
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
      use diis_cc_multipliers_class,      only: diis_cc_multipliers
      use davidson_cc_multipliers_class,  only: davidson_cc_multipliers
!
      implicit none
!
      class(gs_engine), intent(inout) :: engine 
!
      class(ccs), intent(inout) :: wf
!
      type(diis_cc_multipliers), allocatable     :: diis_solver
      type(davidson_cc_multipliers), allocatable :: davidson_solver
!
      call engine%tasks%print_('multipliers solver')
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
         davidson_solver = davidson_cc_multipliers(wf, engine%multipliers_restart)
         call davidson_solver%run(wf)
         call davidson_solver%cleanup(wf)
!
      elseif (trim(engine%multipliers_algorithm) == 'diis') then 
!
         diis_solver = diis_cc_multipliers(wf, engine%multipliers_restart)
         call diis_solver%run(wf)
         call diis_solver%cleanup(wf)
!
      else
!
         call output%error_msg('Could not recognize multipliers algorithm ' // trim(engine%multipliers_algorithm))
!
      endif
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
      call wf%construct_mu(mu_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 3
!
         electronic(k) = wf%calculate_expectation_value(mu_pqk(:,:,k), wf%density)
!
      enddo
!
      call mem%dealloc(mu_pqk, wf%n_mo, wf%n_mo, 3)
!
!     Get nuclear expectation value contribution, then sum the two
!
      call wf%system%get_nuclear_dipole(nuclear)
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
      call wf%construct_q(q_pqk)
!
!     Get electronic expectation value contribution
!
      do k = 1, 6
!
         electronic(k) = wf%calculate_expectation_value(q_pqk(:,:,k), wf%density)
!
      enddo
!
      call mem%dealloc(q_pqk, wf%n_mo, wf%n_mo, 6)
!
!     Get nuclear expectation value contribution, then sum the two
!
      call wf%system%get_nuclear_quadrupole(nuclear)
!
      total = electronic + nuclear
!
   end subroutine calculate_quadrupole_moment_gs_engine
!
!
   subroutine restart_handling_gs_engine(engine, wf)
!!
!!    Restart handling
!!    Written by Sarai D. Folkestad, Nov 2019
!!
!!    Writes the restart information 
!!    if restart is not requested.
!!
!!    If restart is requested performs safety 
!!    checks for restart
!!
      implicit none
!
      class(gs_engine), intent(in) :: engine
      class(ccs), intent(in) :: wf
!
      if (.not. engine%restart) then
!
         call wf%write_cc_restart()
!
      else
!
         if (engine%gs_restart .or. engine%multipliers_restart) &
               call wf%is_restart_safe
!
      endif
!
   end subroutine restart_handling_gs_engine
!
!
   subroutine do_cholesky_gs_engine(wf)
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
!!
!!
      use eri_cd_class, only: eri_cd
      use t1_eri_tool_class, only: t1_eri_tool
!
      implicit none 
!
      class(ccs) :: wf 
!
      type(eri_cd), allocatable :: eri_cholesky_solver 
!
      logical :: do_MO_screening
!
      real(dp), dimension(:,:), allocatable, target :: screening_vector
!
      call output%printf('v', 'Doing Cholesky decomposition of the ERIs.')
!
      do_MO_screening = input%requested_keyword_in_section('mo screening', 'solver cholesky')
!
      eri_cholesky_solver = eri_cd(wf%system)
!
      if (do_MO_screening) then
!
         call output%printf('m', 'Using the MO screening for the Cholesky decomposition', &
                        fs='(/t3,a)')
!
         call mem%alloc(screening_vector, wf%n_ao, wf%n_ao)
         call wf%construct_MO_screening_for_cd(screening_vector)
!
         call eri_cholesky_solver%run(wf%system, screening_vector)   ! Do the Cholesky decomposition
                                                                     ! using the MO screening
!
         call mem%dealloc(screening_vector, wf%n_ao, wf%n_ao) 
!
      else
!
         call eri_cholesky_solver%run(wf%system) ! Do the Cholesky decomposition 
!
      endif
!
      call eri_cholesky_solver%diagonal_test(wf%system)  ! Determine the largest 
                                                         ! deviation in the ERI matrix 
!
      wf%eri = t1_eri_tool(wf%n_o, wf%n_v, eri_cholesky_solver%n_cholesky, wf%need_g_abcd)
!
      call wf%eri%initialize()
!
      call eri_cholesky_solver%construct_cholesky_mo_vectors(wf%system, wf%n_ao, wf%n_mo, &
                                                             wf%orbital_coefficients, wf%eri)
!
      call eri_cholesky_solver%cleanup()
!
   end subroutine do_cholesky_gs_engine
!
!
end module gs_engine_class
