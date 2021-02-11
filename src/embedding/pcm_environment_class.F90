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
module pcm_environment_class
!
!!
!!    PCM environment class module
!!    Written by Tommaso Giovannini, Eirik F. Kjønstad and Sarai D. Folkestad, 2018 and 2020
!!
!
   use parameters
!
   use global_out,           only: output
   use memory_manager_class, only: mem
   use ao_tool_class,        only: ao_tool
   use environment_class,    only: environment
!
#ifdef HAS_PCMSOLVER
   use point_charges_class,  only: point_charges
   use global_in,            only: input 
!
!  External module
!
   use pcmsolver
!
   implicit none
!
   type, extends(environment) :: pcm_environment
!
      type(c_ptr),                      private :: pcm_context
      character(len=200),               private :: input         
      character(len=200),               private :: solvent       
      character(len=200),               private :: solver_type   
      real(dp),                         private :: tesserae_area 
      type(point_charges), allocatable, private :: pcm_points
!
   contains
!
!     Public class routines
!
      procedure :: initialize &
                => initialize_pcm_environment
!
      procedure :: update &
                => update_pcm_environment 
!
      procedure :: get_energy &
                => get_energy_pcm_environment
!
      procedure :: print_energy &
                => print_energy_pcm_environment
!
!     Private class routines
!
      procedure, private :: print_description &
                         => print_description_pcm_environment
!
      procedure, private :: read_parameters &
                         => read_parameters_pcm_environment
!
      procedure, private :: pcmsolver_input &
                         => pcmsolver_input_pcm_environment
!
      procedure, private :: get_nuclei_mm_energy &
                         => get_nuclei_mm_energy_pcm_environment
!
      procedure, private :: get_potential_at_pcm_points &
                         => get_potential_at_pcm_points_pcm_environment
!
   end type  pcm_environment
!
!
   interface  pcm_environment
!
      procedure :: new_pcm_environment
!
   end interface  pcm_environment
!
contains
!
!
   pure function new_pcm_environment() result(embedding)
!!
!!    New pcm environment 
!!    Written by Sarai D. Folkestad, 2020
!!
      implicit none 
!
      type(pcm_environment) :: embedding 
!
      embedding%type_          = 'PCM'
      embedding%n_charges      = 0
      embedding%input          = 'internal'
      embedding%solvent        = '' 
      embedding%tesserae_area  = 0.3d0    ! Angstrom
      embedding%solver_type    = 'iefpcm' ! default IEF
!
   end function new_pcm_environment
!
!
   subroutine initialize_pcm_environment(embedding, ao)
!!
!!    Initialize
!!    Written by Tommaso Giovannini, Oct 2019
!!
!!    Modified by Sarai D. Folkestad, Sep 2020
!!    for environment tool
!!    
      use iso_c_binding,   only: c_char, c_int, c_null_char
      use array_utilities, only: zero_array, copy_and_scale
!
      implicit none
!
      class(pcm_environment), intent(inout)  :: embedding
      type(ao_tool),          intent(inout)  :: ao
!
      real(dp), dimension(:), allocatable    :: r
      type(PCMInput)                         :: host_input
      type(point_charges)                    :: qm_points 
!
      integer(c_int), dimension(4)           :: symmetry_info 
!
      symmetry_info = [int(0,c_int),int(0,c_int),int(0,c_int),int(0,c_int)]
!
      call embedding%print_description()
      call embedding%read_parameters()
!
      if (.not. pcmsolver_is_compatible_library()) &
        call output%error_msg('PCMSolver is not compatible')
!
      host_input = embedding%pcmsolver_input()
! 
!     Extract info on QM system
!
      call ao%get_point_charges(qm_points)
!
      call mem%alloc(r, 3*qm_points%n_charges)
!
      call dcopy(3*qm_points%n_charges, qm_points%r, 1, r, 1)
      call dscal(3*qm_points%n_charges, angstrom_to_bohr, r, 1)
!
!     Initialize pcm solver 
!
      if(trim(embedding%input).eq.'internal') then ! get pcm solver input from eT.inp
!      
         embedding%pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST,               &
                                              int(qm_points%n_charges, kind=c_int), &
                                              qm_points%q, r,                       &
                                              symmetry_info, host_input,            &
                                              c_funloc(host_writer))
!
      else ! get pcm solver input from separate file
!      
         embedding%pcm_context = pcmsolver_new(PCMSOLVER_READER_OWN,                &
                                              int(qm_points%n_charges, kind=c_int), &
                                              qm_points%q, r,                       &
                                              symmetry_info, host_input,            &
                                              c_funloc(host_writer))
!
      endif   
!
      call mem%dealloc(r, 3*qm_points%n_charges) 
!     
!     Collect information on position of pcm point chrges
!                    
      embedding%n_charges = pcmsolver_get_cavity_size(embedding%pcm_context)
!
      embedding%pcm_points = point_charges(embedding%n_charges)
      call embedding%pcm_points%initialize()
!
      call mem%alloc(r, 3*embedding%n_charges)
      call pcmsolver_get_centers(embedding%pcm_context, r)
!
      call dcopy(3*embedding%n_charges, r, 1, embedding%pcm_points%r, 1)
      call dscal(3*embedding%n_charges, bohr_to_angstrom, embedding%pcm_points%r, 1)
!
      call mem%dealloc(r, 3*embedding%n_charges)
!
      call zero_array(embedding%pcm_points%q, embedding%n_charges)
!
      call pcmsolver_print(embedding%pcm_context)
!
      call ao%initialize_oei('electrostatic potential')
!
   end subroutine initialize_pcm_environment
!
!
   subroutine update_pcm_environment(embedding, ao, density)
!!
!!    Update
!!    Written by Tommaso Giovannini
!!
!!    Modified for pcm_environment by Sarai D. Folkestad
!!
!
      use array_utilities, only: zero_array
!
      implicit none
!
      class(pcm_environment),          intent(inout)  :: embedding
      type(ao_tool),                   intent(inout)  :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)     :: density
!
      integer                                         :: I
      real(dp)                                        :: ddot
      real(dp), dimension(:), allocatable             :: potential
      real(dp), dimension(:,:), allocatable           :: v
!
      call mem%alloc(v, ao%n, ao%n)
      call mem%alloc(potential, embedding%n_charges)
!
      call embedding%get_potential_at_pcm_points(potential, ao)
!
      do I = 1, embedding%n_charges
!
         call ao%initialize_external_unit_charges(1, embedding%pcm_points%r(:, I))
         call ao%get_oei('electrostatic potential unit', v)
!
         potential(I) = potential(I) + ddot(ao%n**2, v, 1, density, 1)
!
      enddo
!
      call mem%dealloc(v, ao%n, ao%n)
!
!     Solve PCM equarions to update charges
!
      call pcmsolver_set_surface_function(embedding%pcm_context,                 &
                                          int(embedding%n_charges, kind=c_int),  &
                                          potential,                             &
                                          pcmsolver_fstring_to_carray('NucMEP'))                                    
!
      call pcmsolver_compute_asc(embedding%pcm_context, &
                                 pcmsolver_fstring_to_carray('NucMEP'), &
                                 pcmsolver_fstring_to_carray('NucASC'), &
                                 irrep=0_c_int)                                  
!
      call pcmsolver_get_surface_function(embedding%pcm_context,                 &
                                          int(embedding%n_charges, kind=c_int),  &
                                          embedding%pcm_points%q,                &
                                          pcmsolver_fstring_to_carray('NucASC'))
!
!     Print iteration summary if pl='v'
!
      call output%print_separator('v', 67, fs='(/t3,a)')
      call output%printf('v', 'Atom         PCM ASC            PCM RHS', fs='(t6,a)')
!
      do i = 1, embedding%n_charges
!
         call output%printf('v', '(i4)      (e13.6)      (e13.6)', ints=[i], &
                            reals=[embedding%pcm_points%q(i), &
                            potential(i)], fs='(t6,a)')
!
      enddo
      call output%print_separator('v', 67)
!
      call mem%dealloc(potential, embedding%n_charges)
!
      call ao%initialize_external_charges(embedding%pcm_points)
!
      call ao%construct_stored_oei('electrostatic potential')
!
   end subroutine update_pcm_environment
!
!
   subroutine print_description_pcm_environment(embedding)
!!
!!    Print description
!!    Written by Tommaso Giovannini, April 2019
!!    
      implicit none
!
      class(pcm_environment), intent(in) :: embedding
!
      call output%printf('m', 'This is a PCM calculation', fs='(/t3,a)')
!
      call output%printf('m',  'Polarizable Continuum Solver via PCMSolver', fs='(/t6,a)')
      call output%printf('m', 'For details on PCM, see:', fs='(/t6,a)')
      call output%printf('m', 'Tomasi, Mennucci, Cammi, Chem. Rev. 2005, 105, 2999-3094.', &
                         fs='(t6,a)')
      call output%printf('m', 'For details on PCMSolver, see:', fs='(/t6,a)')
      call output%printf('m', 'Di Remigio et al., IJQC, 2019, 119, e25685', fs='(t6,a)')
!      
!
      if(embedding%input.eq.'external') then 
!
         call output%printf('m', 'PCM Solver was set via external file', fs='(/t6,a)')
!
      else 
!
         call output%printf('m', 'PCM Solver was set via internal parameters', fs='(/t6,a)')
         call output%printf('m', 'Solver Type:   '//trim(embedding%solver_type), fs='(t6,a)')
         call output%printf('m', 'Solvent:       '//trim(embedding%solvent), fs='(t6,a)')
         call output%printf('m', 'Tesserae Area: (f5.3) Å', &
                            reals=[embedding%tesserae_area], fs='(t6,a/)')
!
      endif
!
!
      if(input%requested_cc_calculation()) then 
         call output%printf('m', 'CC calculation: zero-order approximation', fs='(/t6,a)')
         call output%printf('m', 'PCM charges only affect MOs and Fock', fs='(t6,a/)')
      endif
!
   end subroutine print_description_pcm_environment
!
!
   function get_energy_pcm_environment&
                     (embedding, ao, density) result(embedding_energy)
!!
!!    Get energy
!!    Written by Sarai D. Folkestad
!!
!!    Based on PCM routines written by Tommaso Giovannini 
!!
!!    embedding energy = sum_i sum_I q_i Q_I/|r_i - R_I|
!!
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(pcm_environment), intent(in)          :: embedding
      type(ao_tool), intent(in)                   :: ao
      real(dp), dimension(ao%n, ao%n), intent(in) :: density
      real(dp) :: embedding_energy
!
      call do_nothing(density)   
!
      embedding_energy = half*embedding%get_nuclei_mm_energy(ao)
!
   end function get_energy_pcm_environment
!
!
   subroutine print_energy_pcm_environment(embedding, ao, density)
!!
!!    Print energy
!!    Written by Sarai D. Folkestad and Tommaso Giovannini
!!
      implicit none
!
      class(pcm_environment)                       :: embedding
      type(ao_tool),                   intent(in)  :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)  :: density
!
      real(dp) :: scf_energy, electrostatic_energy, ddot
!
!     SCF energy: 1/2 Tr(vD) + sum_i sum_I q_i Q_I/|r_i - R_I|
!
      scf_energy = embedding%get_energy(ao, density)&
                 + half*ddot((ao%n)**2, density, 1, ao%v, 1)
!
!     Electrostatic energy: 
!     Tr(vD) +  sum_i sum_I q_i Q_I/|r_i - R_I|
!
      electrostatic_energy = two*scf_energy
!
      call output%printf('m', '- Summary of QM/PCM energetics:', fs='(/t3,a)')
      call output%printf('m', 'a.u.             eV     kcal/mol', fs='(t42,a)')
!
      call output%printf('m', 'QM/PCM SCF Contribution: (f22.12)(f12.5) (f9.3)', &
                         reals=[scf_energy, &
                                scf_energy*Hartree_to_eV, &
                                scf_energy*Hartree_to_kcalmol], fs='(t6,a)')
!
      call output%printf('m', 'QM/PCM Electrostatic Energy:(f19.12)(f12.5) (f9.3)', &
                         reals=[electrostatic_energy, &
                                electrostatic_energy*Hartree_to_eV, &
                                electrostatic_energy*Hartree_to_kcalmol], fs='(t6,a)')
!
   end subroutine print_energy_pcm_environment
!
!
   subroutine read_parameters_pcm_environment(embedding)
!!
!!    Read parameters PCM
!!    Written by Tommaso Giovannini, October 2019
!!
      implicit none
!
      class(pcm_environment) :: embedding
!
      logical :: file_exists
!      
      call input%get_keyword('input', 'pcm', embedding%input)
!      
      if (trim(embedding%input) .eq. 'internal') then 
!      
         call input%get_keyword('solvent','pcm',embedding%solvent)
         call input%get_keyword('tesserae area','pcm',embedding%tesserae_area)
         call input%get_keyword('solver type','pcm',embedding%solver_type)
!         
         if(len_trim(embedding%solvent) .eq. 0) &
            call output%error_msg('PCM Calculation Internal Input: Solvent Not Specified')
!      
      else if(trim(embedding%input) .eq. 'external') then
!      
         inquire(file='@pcmsolver.inp', exist=file_exists)
!         
         if(.not.file_exists) &
            call output%error_msg('PCM Calculation External Input: Problems with go_pcm.py')
!      
      else
!      
         call output%error_msg('PCM Calculation: neither internal or external input is set')
!      
      endif
!
!
   end subroutine read_parameters_pcm_environment
!
!
   function get_nuclei_mm_energy_pcm_environment(embedding, ao) result(E) 
!!
!!    Get MM energy contribution
!!    Written by Sarai D. Folkestad
!!
!!    Returns the nuclei-MM interaction energy 
!!
!!       E = sum_i sum_I q_i Q_I/|r_i - R_I|
!!
      implicit none
!
      class(pcm_environment), intent(in) :: embedding
      type(ao_tool), intent(in) :: ao
      real(dp) :: E
!
      type(point_charges) :: qm_points
!
      call ao%get_point_charges(qm_points)
!
      E = embedding%pcm_points%get_coulomb_interaction(qm_points)
!
   end function get_nuclei_mm_energy_pcm_environment
!
!
   function pcmsolver_input_pcm_environment(embedding) result(host_input)
!!
!!   PCMSolver input pcm
!!   Adapted by Tommaso Giovannini, Oct 2019
!!
!!   Performs syntactic checks on PCMSolver input
!!   and fills the data structure holding input data
!!
!!   Based on the PCMSolver (LGPL3) routine "pcmsolver_input":
!!
!!      PCMSolver, an API for the Polarizable Continuum Model
!!      Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
!!
!!   Modified to add the possibility of changing values of some properties 
!!   on input (area, solver_type, solvent).
!!
!
      use iso_c_binding, only: c_char, c_int, c_null_char
!
      implicit none
!
      class(pcm_environment) :: embedding
      type(PCMInput) :: host_input
!
!     These parameters would be set by the host input reading subroutine(s)
!     Notice that the strings have a maximum pre-set length.
!     This is to ensure C interoperability
!     Length and area parameters are all assumed to be in Angstrom,
!     the module will convert to Bohr internally
!
      character(kind=c_char, len=7)  :: pcmmod_cavity_type     = 'gepol'//c_null_char
      integer                        :: pcmmod_patch_level     = 2
      real(c_double)                 :: pcmmod_coarsity        = 0.5
      real(c_double)                 :: pcmmod_min_distance    = 0.1
      integer                        :: pcmmod_der_order       = 4
      logical(c_bool)                :: pcmmod_scaling         = .true.
      character(kind=c_char, len=7)  :: pcmmod_radii_set       = 'bondi'//c_null_char
      character(kind=c_char, len=19) :: pcmmod_restart_name    = 'cavity.npz'//c_null_char
      real(c_double)                 :: pcmmod_min_radius      = 100.0
      character(kind=c_char, len=11) :: pcmmod_equation_type   = 'secondkind'//c_null_char
      real(c_double)                 :: pcmmod_correction      = 0.0
      real(c_double)                 :: pcmmod_probe_radius    = 1.0
      character(kind=c_char, len=7)  :: pcmmod_inside_type     = 'vacuum'//c_null_char
      character(kind=c_char, len=21) :: pcmmod_outside_type    = 'uniformdielectric'//c_null_char
      real(c_double)                 :: pcmmod_outside_epsilon = 1.0
!
!     Adjusted dimensions for PCMSolver interface
!
      character(kind=c_char, len=21) :: pcmmod_solvent
      character(kind=c_char, len=7) :: pcmmod_solver_type
!
!
      pcmmod_solver_type = trim(embedding%solver_type)//c_null_char
      pcmmod_solvent     = trim(embedding%solvent)//c_null_char
!
      host_input%cavity_type     = pcmsolver_fstring_to_carray(pcmmod_cavity_type)
!
      host_input%patch_level     = int(pcmmod_patch_level, kind=c_int)
      host_input%coarsity        = pcmmod_coarsity
!      
      host_input%area            = embedding%tesserae_area
      host_input%min_distance    = pcmmod_min_distance
      host_input%der_order       = int(pcmmod_der_order, kind=c_int)
      host_input%scaling         = pcmmod_scaling
      host_input%radii_set       = pcmsolver_fstring_to_carray(pcmmod_radii_set)
      host_input%restart_name    = pcmsolver_fstring_to_carray(pcmmod_restart_name)
!
      host_input%min_radius      = pcmmod_min_radius
!
      host_input%solver_type     = pcmsolver_fstring_to_carray(pcmmod_solver_type)
!
      host_input%solvent         = pcmsolver_fstring_to_carray(pcmmod_solvent)
!
      host_input%equation_type   = pcmsolver_fstring_to_carray(pcmmod_equation_type)
!
      host_input%correction      = pcmmod_correction
      host_input%probe_radius    = pcmmod_probe_radius
!
      host_input%inside_type     = pcmsolver_fstring_to_carray(pcmmod_inside_type)
!
      host_input%outside_epsilon = pcmmod_outside_epsilon
      host_input%outside_type    = pcmsolver_fstring_to_carray(pcmmod_outside_type)
!
   end function pcmsolver_input_pcm_environment
!
!
   subroutine host_writer(message) bind(c, name='host_writer')
!!
!!   Host writer
!!   Adapted by Tommaso Giovannini, Oct 2019
!!
!!   Writer Interface for PCMSolver
!!
!!   Based on the PCMSolver (LGPL3) routine "host_writer":
!!
!!      PCMSolver, an API for the Polarizable Continuum Model
!!      Copyright (C) 2019 Roberto Di Remigio, Luca Frediani and contributors.
!!
!!   Modified to fit with the default indentation of eT messages.
!!
!
     character(kind=c_char), intent(in) :: message(*)
     integer(c_int) :: length, length_2, i, n_messages
     character(len=1000), allocatable :: message_eT(:)
!
     length = 0_c_int
     n_messages = 0_c_int
!
     do
!
       if (message(length + 1_c_int) == c_null_char) exit
!
       if (message(length + 1_c_int) == char(10)) then 
!         
          n_messages = n_messages + 1_c_int
!
       endif
!       
       length = length + 1_c_int
!
     end do
!     
!     
     allocate(message_eT(n_messages))
!
     message_eT = ' '
!     
     length_2 = 1_c_int
     n_messages = 1_c_int
!
     do i = 1_c_int, length
!
        if (message(i) == c_null_char) exit
!
        if (message(i) .eq. char(10)) then
!
            write(message_eT(n_messages),'(1000a)') message(length_2:i-1) 
            length_2 = i + 1_c_int
            n_messages = n_messages + 1_c_int
!
        endif
!
     enddo
!
     do i = 1_c_int, n_messages-1_c_int
!
      call output%printf('n', trim(message_eT(i)), fs='(t6,a)')
!
     enddo
!
     deallocate(message_eT)
!
   end subroutine host_writer
!
!
    subroutine get_potential_at_pcm_points_pcm_environment(embedding, potential, ao)
!!
!!    Get potential at pcm points
!!    Written by Sarai D. Folkestad
!!
!!    Returns the potential at the pcm charges from the 
!!    nuclei
!!
!!      potential(I) = sum_J Q_J/|R_J - r_I|
!!
!!    where Q_J and R_J are the charge and position of the QM nuclei
!!    and r_I is an pcm point 
!!
      implicit none
!
      class(pcm_environment),                   intent(in)  :: embedding
      real(dp), dimension(embedding%n_charges), intent(out) :: potential
      type(ao_tool),                            intent(in)  :: ao
!
      type(point_charges)                                   :: qm_points
!
      call ao%get_point_charges(qm_points)
!
!     Get potential at pcm points
!
      call qm_points%get_potential_at_external_points(             &
                                  embedding%pcm_points%n_charges,  &
                                  embedding%pcm_points%r,          &
                                  potential)
!
    end subroutine get_potential_at_pcm_points_pcm_environment
!
!
#else
!
!  Dummy type to scream errors at you if you compile without PCMSolver
!  but call pcm in input
!
   implicit none
!
   type, extends(environment) :: pcm_environment
!  
!      
   contains
!
      procedure :: initialize &
                => initialize_pcm_environment
!
      procedure :: update &
                => update_pcm_environment 
!
      procedure :: get_energy &
                => get_energy_pcm_environment
!
      procedure :: print_energy &
                => print_energy_pcm_environment
!
   end type pcm_environment
!
   interface pcm_environment
!
      procedure :: new_pcm_environment
!
   end interface pcm_environment
!
contains
!
   function new_pcm_environment() result(embedding)
!
      implicit none
!
      type(pcm_environment) :: embedding
!
      embedding%type_                = ''
      embedding%n_charges            = 0 
!
      call output%error_msg('constructor was called for pcm environment, &
                           &, but eT is not compiled with pcm. &
                           &Use --pcm when running setup to link to PCMSolver.')
!
   end function new_pcm_environment
!
   subroutine initialize_pcm_environment(embedding, ao)
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(pcm_environment), intent(inout) :: embedding
      type(ao_tool),          intent(inout) :: ao
!
      call output%error_msg('initialize was called for pcm environment, &
                           &, but eT is not compiled with pcm. &
                           &Use --pcm when running setup to link to PCMSolver.')
!
      call do_nothing(embedding)
      call do_nothing(ao)
!
   end subroutine
!
!
   subroutine update_pcm_environment(embedding, ao, density)
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(pcm_environment), intent(inout)             :: embedding
      type(ao_tool), intent(inout)                      :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)       :: density
!
      call output%error_msg('update was called for pcm environment, &
                           &, but eT is not compiled with pcm. &
                           &Use --pcm when running setup to link to PCMSolver.')
!
      call do_nothing(embedding)
      call do_nothing(ao)
      call do_nothing(density)
!
   end subroutine
!
!
   function get_energy_pcm_environment&
                     (embedding, ao, density) result(embedding_energy)
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(pcm_environment), intent(in)              :: embedding
      type(ao_tool), intent(in)                       :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)     :: density
!
      real(dp)                                        :: embedding_energy
!
      embedding_energy = zero
!
      call output%error_msg('get_energy was called for pcm environment, &
                           &, but eT is not compiled with pcm. &
                           &Use --pcm when running setup to link to PCMSolver.')
!
      call do_nothing(embedding)
      call do_nothing(ao)
      call do_nothing(density)
!
   end function get_energy_pcm_environment
!
!
   subroutine print_energy_pcm_environment(embedding, ao, density)
!
      use warning_suppressor, only: do_nothing
!
      implicit none
!
      class(pcm_environment)                          :: embedding
      type(ao_tool), intent(in)                       :: ao
      real(dp), dimension(ao%n, ao%n), intent(in)     :: density
!
      call output%error_msg('print_energy was called for pcm environment, &
                           &, but eT is not compiled with pcm. &
                           &Use --pcm when running setup to link to PCMSolver.')
!
      call do_nothing(embedding)
      call do_nothing(ao)
      call do_nothing(density)
!
   end subroutine
!
#endif
!
end module pcm_environment_class
