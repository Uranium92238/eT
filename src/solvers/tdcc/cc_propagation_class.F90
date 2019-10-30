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
module cc_propagation_class
!
!!
!! Abstract coupled cluster propagation solver class module.
!! Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!! Module containing the abstract type cc_propagation. This abstract type can, when
!! extended, be used to generate objects that can propagate a time-dependent coupled cluster
!! wavefunction in the prescence of an external time-dependent electric field. During the
!! propagation, the time-dependent amplitudes and multipliers evolve according to the equations
!!
!!    d/dt t_mu(t)    = -i*Omega_mu(t)
!!    d/dt tbar_nu(t) = i*<Lambda(t)|[H(t),tau_nu]|CC(t)>
!!
!! See J. Chem. Phys. 93, 3333 (1990) for more details. The abstract procedure propagation_step
!! needs to be overwritten in descendant types, and should specify a single step of a one-step
!! integration method. The tag and and name should also be set in the descendant type object
!! generator, see e.g. euler_cc_propagation_class.
!!
!! During the propagation, time-dependent properties can be calculated, and written as time series
!! to files. See read_settings for more details.
!!
!
   use parameters
   use memory_manager_class, only: mem
   use global_in, only: input
   use global_out, only: output
   use timings_class, only: timings
   use sequential_file_class, only: sequential_file
   use ccs_class, only: ccs
   use electric_field_class, only: electric_field
!
   implicit none
!
   type, abstract :: cc_propagation
!
      character(len=100) :: tag
      character(len=100) :: name_
      character(len=100) :: author
      character(len=500) :: description
!
      integer :: steps_between_output, vector_length
!
      real(dp) :: ti, tf, h, implicit_threshold
!
      complex(dp), dimension(:,:,:,:), allocatable :: g_pqrs_complex_ti
      complex(dp), dimension(:,:), allocatable     :: t1_complex_ti
!
      logical :: energy_output, dipole_moment_output, electric_field_output, amplitudes_output, multipliers_output
!
      type(sequential_file) :: energy_file
      type(sequential_file) :: dipole_moment_file
      type(sequential_file) :: electric_field_file
      type(sequential_file) :: amplitudes_file
      type(sequential_file) :: multipliers_file
!
   contains
!
      procedure :: new_cc_propagation              => new_cc_propagation_cc_propagation
      procedure :: run                             => run_cc_propagation
      procedure :: cleanup                         => cleanup_cc_propagation
!
      procedure :: read_settings                   => read_settings_cc_propagation
      procedure :: print_settings                  => print_settings_cc_propagation
      procedure :: print_banner                    => print_banner_cc_propagation
      procedure :: print_summary                   => print_summary_cc_propagation
!
      procedure :: update_field_and_wavefunction   => update_field_and_wavefunction_cc_propagation
      procedure :: update_wavefunction             => update_wavefunction_cc_propagation
!
      procedure, nopass :: calculate_energy        => calculate_energy_cc_propagation
      procedure, nopass :: calculate_dipole_moment => calculate_dipole_moment_cc_propagation
!
      procedure :: write_energy_to_file            => write_energy_to_file_cc_propagation
      procedure :: write_dipole_moment_to_file     => write_dipole_moment_to_file_cc_propagation
      procedure :: write_electric_field_to_file    => write_electric_field_to_file_cc_propagation
      procedure :: write_amplitudes_to_file        => write_amplitudes_to_file_cc_propagation
      procedure :: write_multipliers_to_file       => write_multipliers_to_file_cc_propagation
!
      procedure :: open_files                      => open_files_cc_propagation
      procedure :: close_files                     => close_files_cc_propagation
!
      procedure(propagation_step), deferred        :: step ! Single propagation step for the selected method
!
   end type cc_propagation  
!
   abstract interface
!
      subroutine propagation_step(solver, wf, field, ti, dt, ui, uf, n)
!
         use kinds
         use ccs_class, only: ccs
         use electric_field_class, only: electric_field
         import :: cc_propagation
!
         implicit none
!
         class(cc_propagation), intent(inout) :: solver
         class(ccs), intent(inout) :: wf
         class(electric_field), intent(inout) :: field
!
         real(dp), intent(in) :: dt
         real(dp), intent(in) :: ti
!
         integer, intent(in) :: n
!
         complex(dp), dimension(n), intent(in) :: ui
         complex(dp), dimension(n), intent(out) :: uf
!
      end subroutine propagation_step
!
   end interface
!
!
contains
!
!
   subroutine new_cc_propagation_cc_propagation(solver, wf)
!!
!!    New cc propagation
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Should be called by generators in all descendants. Sets printables, reads settings and
!!    allocates arrays.
!!
      implicit none
!
      class(cc_propagation) :: solver
      class(ccs) :: wf
!
!     Set printables
!
      solver%author = 'A. Skeidsvoll, A. Balbi, 2018'
      solver%description = 'A solver that propagates a coupled cluster wavefunction from the ' &
                            // 'ground state using the ' // trim(solver%tag) // ' method.'
!
      call solver%print_banner()
!
!     Set default settings
!
      solver%steps_between_output = 1
!
      solver%implicit_threshold = 1.0e-11_dp
!
      solver%energy_output           = .false.
      solver%dipole_moment_output    = .false.
      solver%electric_field_output   = .false.
      solver%amplitudes_output       = .false.
      solver%multipliers_output      = .false.
!
      call solver%read_settings()
      call solver%print_settings()
!
      solver%vector_length = 2*wf%n_gs_amplitudes
!
!     Initial state t1 and t1 transformed g_pqrs
!
      call mem%alloc(solver%t1_complex_ti, wf%n_v, wf%n_o)!
      call mem%alloc(solver%g_pqrs_complex_ti, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
      call zcopy(wf%n_t1, wf%t1_complex, 1, solver%t1_complex_ti, 1)
      call zcopy(wf%n_mo**4, wf%integrals%g_pqrs_complex, 1, solver%g_pqrs_complex_ti, 1)
!
!     Complex density matrix used to calculate properties
!
      if (solver%energy_output .or. solver%dipole_moment_output) &
         call wf%initialize_gs_density_complex()
!
      call solver%open_files()
!
   end subroutine new_cc_propagation_cc_propagation
!
!
   subroutine run_cc_propagation(solver, wf, field)
!!
!!    Run
!!    Written by Alice Balbi and Andreas Skeidsvoll, Oct 2018
!!
!!    Do the propagation of the coupled cluster wavefunction stepwise until final time is
!!    reached. Properties are calculated and written to file if requested in settings.
!!
      implicit none
!
      class(cc_propagation), intent(inout) :: solver
      class(ccs) :: wf
      class(electric_field) :: field
!
      complex(dp), dimension(:), allocatable :: amplitudes_multipliers, amplitudes_multipliers_next
!
      real(dp) :: t, dt, t_next
!
      integer :: step
!
      logical :: last
!
      type(timings) :: step_timer
!
!     Create vector containing the amplitudes and multipliers
!
      call mem%alloc(amplitudes_multipliers, solver%vector_length)
      call mem%alloc(amplitudes_multipliers_next, solver%vector_length)
!
      call wf%get_amplitudes_complex(amplitudes_multipliers(1:wf%n_gs_amplitudes))
      call wf%get_multipliers_complex( &
         amplitudes_multipliers(wf%n_gs_amplitudes+1:solver%vector_length))
!
!     Propagate until tf
!
      t = solver%ti
      step = 0
      last = .false.
!
      do
!
         call output%printf('Time: (f10.4) au', reals=[t], pl='verbose', fs='(/t3,a)')
!
         call solver%update_field_and_wavefunction(wf, field, t, amplitudes_multipliers)
!
!        Write properties to file
!
         if (mod(step, solver%steps_between_output) == 0) then
!
            if (solver%energy_output .or. solver%dipole_moment_output) &
               call wf%construct_gs_density_complex()
!
            if (solver%energy_output) then
               call solver%calculate_energy(wf, field)
               call solver%write_energy_to_file(wf, t)
            endif
!
            if (solver%dipole_moment_output) then
               call solver%calculate_dipole_moment(wf)
               call solver%write_dipole_moment_to_file(wf, t)
            endif
!
            if (solver%electric_field_output)   call solver%write_electric_field_to_file(field, t)
            if (solver%amplitudes_output)       call solver%write_amplitudes_to_file(wf, t)
            if (solver%multipliers_output)      call solver%write_multipliers_to_file(wf, t)
!
            call output%printf('Properties written at time: (f10.4) au', reals=[t], pl='normal', &
                               fs='(t3,a)')
!
         endif
!
!        Count current step
!
         step = step + 1
!
         if (last) exit
!
!        Make a trial next time value
!
         t_next = solver%ti + real(step, dp)*(solver%h)
!
!        Check if trial next time value goes beyond final time value, adjust size of last time 
!        step if it does
!
         last = (((dt .gt. zero) .and. (t_next .ge. solver%tf)) .or. ((dt .lt. zero) .and. &
                                                                      (t_next .le. solver%tf)))
!
         if (last) then
            dt = solver%tf - t
         else
            dt = solver%h
         endif
!
         if (dt .eq. zero) exit
!
!        Do a propagation step
!
         step_timer = timings('cc propagation step', pl='verbose')
         call step_timer%turn_on()
!
         call solver%step(wf, field, t, dt, amplitudes_multipliers, amplitudes_multipliers_next, &
                          solver%vector_length) 
!
         call step_timer%turn_off()
         call step_timer%reset()
!
!        Update to next amplitudes_multipliers
!
         call zcopy(solver%vector_length, amplitudes_multipliers_next, 1, amplitudes_multipliers, &
                    1)
!
!        Calculate next time value
!
         if (last) then
            t = solver%ti + real(step-1, dp)*(solver%h) + dt
         else
            t = solver%ti + real(step, dp)*(solver%h)
         endif
!
      enddo
!
      call mem%dealloc(amplitudes_multipliers, solver%vector_length)
      call mem%dealloc(amplitudes_multipliers_next, solver%vector_length)
!
      call output%printf('Finished propagating for (f10.4) au', reals=[solver%tf], pl='minimal', &
                         fs='(/t3,a)')
!
   end subroutine run_cc_propagation
!
!
   subroutine cleanup_cc_propagation(solver, wf, field)
!!
!!    Clean up
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Clean up the abstract coupled cluster propagation solver.
!!
      implicit none
!
      class(cc_propagation), intent(inout) :: solver
!
      class(ccs) :: wf
      class(electric_field) :: field
!
      call solver%close_files()
!
!     Initial state t1 and t1 transformed g_pqrs
!
      call mem%dealloc(solver%t1_complex_ti, wf%n_v, wf%n_o)
      call mem%dealloc(solver%g_pqrs_complex_ti, wf%n_mo, wf%n_mo, wf%n_mo, wf%n_mo)
!
!     Complex density matrix used to calculate properties
!
      if (solver%energy_output .or. solver%dipole_moment_output) call wf%destruct_gs_density()
!
      call solver%print_summary(wf, field)
!
   end subroutine cleanup_cc_propagation
!
!
   subroutine update_field_and_wavefunction_cc_propagation(solver, wf, field, t, &
                                                           amplitudes_multipliers)
!!
!!    Update field and wavefunction
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Update the electric field to the field at the given time t, and then the wavefunction
!!    (amplitudes, multipliers, t1 transformed quantities and Fock matrix).
!!
      implicit none 
!
      class(cc_propagation)                            :: solver
      class(ccs), intent(inout)                                 :: wf
      class(electric_field), intent(inout)                      :: field
      real(dp), intent(in)                                      :: t
      complex(dp), dimension(solver%vector_length), intent(in)  :: amplitudes_multipliers
!
!     Update field
!
      call field%update(t)
!
!     Update wavefunction
!
      call solver%update_wavefunction(wf, field, amplitudes_multipliers)
!
   end subroutine update_field_and_wavefunction_cc_propagation
!
!
   subroutine update_wavefunction_cc_propagation(solver, wf, field, amplitudes_multipliers)
!!
!!    Update wavefunction
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Update amplitudes, multipliers, t1 transformed quantities and Fock matrix.
!!
      implicit none 
!
      class(cc_propagation)                            :: solver
      class(ccs), intent(inout)                                 :: wf
      class(electric_field), intent(inout)                      :: field
      complex(dp), dimension(solver%vector_length), intent(in)  :: amplitudes_multipliers
!
      call wf%set_amplitudes_complex(amplitudes_multipliers(1:wf%n_gs_amplitudes))
      call wf%set_multipliers_complex( &
         amplitudes_multipliers(wf%n_gs_amplitudes+1:solver%vector_length))
!
      call wf%t1_transform_4_complex(solver%g_pqrs_complex_ti, wf%integrals%g_pqrs_complex, &
                                     wf%t1_complex - solver%t1_complex_ti)
!
      call wf%construct_fock_complex()
!
!     Contribution to Fock matrix corresponding to the interaction with the electric field.
!
      call wf%add_t1_fock_em_lg_dipole_contribution_complex(cmplx(field%vector, zero, dp))
!
   end subroutine update_wavefunction_cc_propagation
!
!
   subroutine print_banner_cc_propagation(solver)
!!
!!    Print banner
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Print banner of solver.
!!
      implicit none 
!
      class(cc_propagation) :: solver 
!
      call output%printf(':: ' // solver%name_, pl='minimal', fs='(//t3,a)')
      call output%printf(':: ' // solver%author, pl='minimal', fs='(t3,a)')
      call output%printf(solver%description, pl='normal', ffs='(/t3,a)', fs='(t3,a)')
!
   end subroutine print_banner_cc_propagation
!
!
   subroutine calculate_energy_cc_propagation(wf, field)
!!
!!    Calculate energy
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Calculate the energy.
!!
      implicit none 
!
      class(ccs)                     :: wf
      class(electric_field)          :: field
!
!     Total variational energy with interaction, <Λ| e^-T H e^T |R>
!
      call wf%calculate_energy_complex()
      call wf%calculate_energy_omega_contribution_complex()
      call wf%calculate_energy_em_lg_dp_contribution_complex(cmplx(field%vector, zero, dp))
!
   end subroutine calculate_energy_cc_propagation
!
!
   subroutine write_energy_to_file_cc_propagation(solver, wf, t)
!!
!!    Write energy to file
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Write energy to file together with the time t.
!!
      implicit none 
!
      class(cc_propagation) :: solver
      class(ccs)                     :: wf
      real(dp), intent(in)           :: t
!
!     Write energy to file
!
      call solver%energy_file%write_([t, real(wf%energy_complex), aimag(wf%energy_complex)], 3)
!
   end subroutine write_energy_to_file_cc_propagation
!
!
   subroutine calculate_dipole_moment_cc_propagation(wf)
!!
!!    Calculate dipole moment
!!    Written by Andreas Skeidsvoll and Alice Balbi, Oct 2018
!!
!!    Calculate the dipole moment vector.
!!
      implicit none 
!
      class(ccs) :: wf
!
      real(dp), dimension(3) :: nuclear
!
      complex(dp), dimension(:,:,:), allocatable :: mu
!
      complex(dp), external :: zdotu
!
!     Do contraction between electronic dipole integrals and electronic density
!
      call mem%alloc(mu, wf%n_mo, wf%n_mo, 3)
      call wf%construct_mu_complex(mu)
!
      wf%dipole_moment_complex(1) = zdotu(wf%n_mo*wf%n_mo, wf%density_complex, 1, mu(:,:,1), 1)
      wf%dipole_moment_complex(2) = zdotu(wf%n_mo*wf%n_mo, wf%density_complex, 1, mu(:,:,2), 1)
      wf%dipole_moment_complex(3) = zdotu(wf%n_mo*wf%n_mo, wf%density_complex, 1, mu(:,:,3), 1)
!
      call mem%dealloc(mu, wf%n_mo, wf%n_mo, 3)
!
!     Add the nuclear contribution
!
      call wf%system%get_nuclear_dipole(nuclear)
!
      wf%dipole_moment_complex(1) = wf%dipole_moment_complex(1) + cmplx(nuclear(1), zero, dp)
      wf%dipole_moment_complex(2) = wf%dipole_moment_complex(2) + cmplx(nuclear(2), zero, dp)
      wf%dipole_moment_complex(3) = wf%dipole_moment_complex(3) + cmplx(nuclear(3), zero, dp)
!
   end subroutine calculate_dipole_moment_cc_propagation
!
!
   subroutine write_dipole_moment_to_file_cc_propagation(solver, wf, t)
!!
!!    Write dipole moment to file
!!    Written by Andreas Skeidsvoll and Alice Balbi, Oct 2018
!!
!!    Write the dipole moment vector to file together with the time t.
!!
      implicit none 
!
      class(cc_propagation) :: solver
      class(ccs) :: wf
      real(dp), intent(in) :: t
!
!     Write dipole moment to file
!
      call solver%dipole_moment_file%write_([t,                                   &
                                             real(wf%dipole_moment_complex(1)),   &
                                             aimag(wf%dipole_moment_complex(1)),  &
                                             real(wf%dipole_moment_complex(2)),   &
                                             aimag(wf%dipole_moment_complex(2)),  &
                                             real(wf%dipole_moment_complex(3)),   &
                                             aimag(wf%dipole_moment_complex(3))], 7)
!
   end subroutine write_dipole_moment_to_file_cc_propagation
!
!
   subroutine write_electric_field_to_file_cc_propagation(solver, field, t)
!!
!!    Write dipole moment to file
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Write the electric field vector to file together with the time t.
!!
      implicit none 
!
      class(cc_propagation) :: solver
      class(electric_field) :: field
      real(dp), intent(in) :: t
!
      call solver%electric_field_file%write_([t, field%vector(1), zero, &
                                                 field%vector(2), zero, &
                                                 field%vector(3), zero], 7)
!
   end subroutine write_electric_field_to_file_cc_propagation
!
!
   subroutine write_amplitudes_to_file_cc_propagation(solver, wf, t)
!!
!!    Write amplitudes to file
!!    Written by Alice Balbi and Andreas Skeidsvoll, Mar 2019
!!
!!    Write the real part of the vector of amplitudes to file.
!!
      implicit none
!
      class(cc_propagation) :: solver
      class(ccs) :: wf
      real(dp), intent(in) :: t
!
      complex(dp), dimension(:), allocatable :: amplitudes_complex
!
      call mem%alloc(amplitudes_complex, wf%n_gs_amplitudes)
      call wf%get_amplitudes_complex(amplitudes_complex)
!
      call solver%amplitudes_file%write_([t, real(amplitudes_complex)], wf%n_gs_amplitudes + 1)
!
      call mem%dealloc(amplitudes_complex, wf%n_gs_amplitudes)
!
   end subroutine write_amplitudes_to_file_cc_propagation
!
!
   subroutine write_multipliers_to_file_cc_propagation(solver, wf, t)
!!
!!    Write  multipliers to file
!!    Written by Alice Balbi and Andreas Skeidsvoll, Mar 2019
!!
!!    Write the real part of the vector of multipliers to file together with the time t.
!!
      implicit none
!
      class(cc_propagation) :: solver
      class(ccs) :: wf
      real(dp), intent(in) :: t
!
      complex(dp), dimension(:), allocatable :: multipliers_complex
!
      call mem%alloc(multipliers_complex, wf%n_gs_amplitudes)
      call wf%get_multipliers_complex(multipliers_complex)
!
      call solver%multipliers_file%write_([t, real(multipliers_complex)], &
                                          wf%n_gs_amplitudes + 1)
!
      call mem%dealloc(multipliers_complex, wf%n_gs_amplitudes)
!
   end subroutine write_multipliers_to_file_cc_propagation
!
!
   subroutine print_summary_cc_propagation(solver, wf, field)
!!
!!    Print summary 
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Print the summary of the solver.
!!
      implicit none 
!
      class(cc_propagation) :: solver
      class(ccs) :: wf
      class(electric_field) :: field
!
      if (solver%energy_output) then
         call solver%calculate_energy(wf, field)
      endif
!
      if (solver%dipole_moment_output) then
         call solver%calculate_dipole_moment(wf)
      endif
!
      call output%printf('- ' // trim(solver%name_) // ' summary', pl='minimal', fs='(/t3,a)')
      call output%printf('Energy after propagation [au]:', pl='minimal', fs='(/t6,a)')
      call output%print_separator('minimal', 42, '-', fs='(t6,a)')
      call output%printf('    Real part         Imaginary part      ', pl='minimal', fs='(t6,a)')
      call output%print_separator('minimal', 42, '-', fs='(t6,a)')
      call output%printf('    (f16.10)  (f16.10)',                                                &
                         reals=[real(wf%energy_complex), aimag(wf%energy_complex)], pl='minimal', &
                         fs='(t6,a)')
      call output%print_separator('minimal', 42, '-', fs='(t6,a)')
      call output%printf('Dipole moment after propagation [au]:', pl='minimal', fs='(/t6,a)')
      call output%print_separator('minimal', 45, '-', fs='(t6,a)')
      call output%printf('       Real part         Imaginary part      ', pl='minimal', &
                         fs='(t6,a)')
      call output%print_separator('minimal', 45, '-', fs='(t6,a)')
      call output%printf('    x  (f16.10)  (f16.10)', pl='minimal', fs='(t6,a)', &
                         reals=[real(wf%dipole_moment_complex(1)),               &
                                aimag(wf%dipole_moment_complex(1))])
      call output%printf('    y  (f16.10)  (f16.10)', pl='minimal', fs='(t6,a)', &
                         reals=[real(wf%dipole_moment_complex(2)),               &
                                aimag(wf%dipole_moment_complex(2))])
      call output%printf('    z  (f16.10)  (f16.10)', pl='minimal', fs='(t6,a)', &
                         reals=[real(wf%dipole_moment_complex(3)),               &
                                aimag(wf%dipole_moment_complex(3))])
      call output%print_separator('minimal', 42, '-', fs='(t6,a)')
!
   end subroutine print_summary_cc_propagation
!
!
   subroutine read_settings_cc_propagation(solver)
!!
!!    Read settings
!!    Written by Andreas Skeidsvoll, Feb 2019 
!!
!!    Read the settings required to do the propagation.
!!
!!    Settings:
!!
!!    initial time:          the solver starts propagating the coupled cluster wavefunction from
!!                           this time
!!    final time:            the solver stops propagating the coupled cluster wavefunction when
!!                           reaching this time
!!    time step:             the time step used in the integration method
!!    steps between output:  how many time steps should pass between each time the solver writes
!!                           properties to file
!!    implicit treshold:     convergence threshold for the solution of equations in implicit
!!                           integrators
!!    energy output:         writes energy time series to file if this is specified
!!    dipole moment output:  writes dipole moment vector time series to file if this is specified
!!    electric field output: writes electric field vector time series to file if this is specified
!!    density diag output:   writes real part of the density matrix diagonal time series if this is
!!                           specified
!!    amplitudes output:     writes real part of the vector of amplitudes time series to file if
!!                           this is specified
!!    multipliers output:    writes real part of the vector of multipliers time series to file if
!!                           this is specified
!!
      implicit none 
!
      class(cc_propagation) :: solver
!
      call input%get_required_keyword_in_section('initial time', 'solver cc propagation', &
                                                 solver%ti)
      call input%get_required_keyword_in_section('final time', 'solver cc propagation', solver%tf)
      call input%get_required_keyword_in_section('time step', 'solver cc propagation', solver%h)
      call input%get_keyword_in_section('steps between output', 'solver cc propagation', &
                                        solver%steps_between_output)
      call input%get_keyword_in_section('implicit threshold', 'solver cc propagation', &
                                        solver%implicit_threshold)
!
      if (input%requested_keyword_in_section('energy output', 'solver cc propagation')) &
         solver%energy_output = .true.
      if (input%requested_keyword_in_section('dipole moment output', 'solver cc propagation')) &
         solver%dipole_moment_output = .true.
      if (input%requested_keyword_in_section('electric field output', 'solver cc propagation')) &
         solver%electric_field_output = .true.
      if (input%requested_keyword_in_section('amplitudes output', 'solver cc propagation')) &
         solver%amplitudes_output = .true.
      if (input%requested_keyword_in_section('multipliers output', 'solver cc propagation')) &
         solver%multipliers_output = .true.
!
   end subroutine read_settings_cc_propagation
!
!
   subroutine print_settings_cc_propagation(solver)
!!
!!    Print settings    
!!    Written by Andreas Skeidsvoll, Oct 2018
!!
!!    Print the settings of the solver related to the time of propagation
!!
      implicit none 
!
      class(cc_propagation) :: solver 
!
      call output%printf('- Propagation settings:', pl='minimal', fs='(/t3,a)')
!
      call output%printf('Initial time: (f10.4) au', reals=[solver%ti], pl='minimal', fs='(/t6,a)')
      call output%printf('Final time:   (f10.4) au', reals=[solver%tf], pl='minimal', fs='(t6,a)')
      call output%printf('Time step:    (f10.4) au', reals=[solver%h], pl='minimal', fs='(t6,a)')
!
   end subroutine print_settings_cc_propagation
!
!
   subroutine open_files_cc_propagation(solver)
!!
!!    Open files
!!    Written by Andreas Skeidsvoll, Apr 2019
!!
!!    Open the files where the property time series should be written.
!!
      implicit none
!
      class(cc_propagation) :: solver
!
      if (solver%energy_output) then
         solver%energy_file = sequential_file('cc_propagation_energy', 'formatted')
         call solver%energy_file%open_('write', 'rewind')
      endif
!
      if (solver%dipole_moment_output) then
         solver%dipole_moment_file = sequential_file('cc_propagation_dipole_moment', 'formatted')
         call solver%dipole_moment_file%open_('write', 'rewind')
      endif
!
      if (solver%electric_field_output) then
         solver%electric_field_file = sequential_file('cc_propagation_electric_field', 'formatted')
         call solver%electric_field_file%open_('write', 'rewind')
      endif
!
      if (solver%amplitudes_output) then
         solver%amplitudes_file = sequential_file('cc_propagation_amplitudes', 'formatted')
         call solver%amplitudes_file%open_('write', 'rewind')                                           
      endif
!
      if (solver%multipliers_output) then
         solver%multipliers_file = sequential_file('cc_propagation_multipliers', 'formatted')
         call solver%multipliers_file%open_('write','rewind')                                           
      endif
!
   end subroutine open_files_cc_propagation
!
!
   subroutine close_files_cc_propagation(solver)
!!
!!    Close files
!!    Written by Andreas Skeidsvoll, Apr 2019
!!
!!    Close the files where the property time series have been written.
!!
      implicit none
!
      class(cc_propagation) :: solver
!
      if (solver%energy_output) then
         call solver%energy_file%close_
      endif
!
      if (solver%dipole_moment_output) then
         call solver%dipole_moment_file%close_
      endif
!
      if (solver%electric_field_output) then
         call solver%electric_field_file%close_
      endif
!
      if (solver%amplitudes_output) then
         call solver%amplitudes_file%close_
      endif
!
      if (solver%multipliers_output) then
         call solver%multipliers_file%close_
      endif
!
   end subroutine close_files_cc_propagation
!
!
end module cc_propagation_class
