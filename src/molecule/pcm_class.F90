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
module pcm_class
!
!!
!!    PCM class module
!!    Written by Tommaso Giovannini, Oct 2019
!!
!!    Reference : J. Tomasi, B. Mennucci, R. Cammi.
!!    Chem. Rev. 2005, 105, 2999-3094.
!!
!
   use global_in, only : input
   use global_out, only : output
   use parameters
!
#ifdef HAS_PCMSOLVER
   use atomic_class, only: atomic
   use memory_manager_class, only: mem
   use, intrinsic :: iso_c_binding, only:c_char, c_int, c_null_char
   use pcmsolver
!
   implicit none
!
   type :: pcm
!
      type(c_ptr)    :: pcm_context
!
      integer :: n_tesserae              ! number of tesserae
      character(len=200) :: input        ! input options : internal or external 
      character(len=200) :: solvent      ! solvent type 
      character(len=200) :: solver_type  ! solvent type 
      real(dp) :: tesserae_area          ! area of the tesserae
!      
      real(dp), dimension(:), allocatable   :: grid_coord   ! PCM coordinates
      real(dp), dimension(:), allocatable   :: charges      ! PCM charges 
!      
      real(dp), dimension(:), allocatable   :: pcm_rhs    ! RHS for polarizable PCM
!      
   contains
!
      procedure :: prepare                         => prepare_pcm
      procedure :: cleanup                         => cleanup_pcm
      procedure :: read_parameters                 => read_parameters_pcm
      procedure :: pcmsolver_input                 => pcmsolver_input_pcm
      procedure :: print_description_and_settings  => print_description_and_settings_pcm
!
      procedure :: set_surface_function            => set_surface_function_pcm
      procedure :: get_surface_function            => get_surface_function_pcm
      procedure :: compute_asc                     => compute_asc_pcm
!
   end type pcm
!
contains
!
!
   subroutine prepare_pcm(molecule, n_atoms_qm, qm_coordinates, qm_charges)
!!
!!    Prepare PCM system
!!    Written by Tommaso Giovannini, Oct 2019
!!
      implicit none
!
      class(pcm) :: molecule
!
      integer :: n_atoms_qm
      real(dp), dimension(3*n_atoms_qm) :: qm_coordinates
      real(dp), dimension(n_atoms_qm)   :: qm_charges
!
!     stuff for pcmsolver
!
      type(PCMInput) :: host_input
!      
      integer(c_int) :: symmetry_info(4)
!      
      if (.not. pcmsolver_is_compatible_library()) &
!      
        call output%error_msg('PCMSolver is not compatible')
!
      call molecule%read_parameters()

      symmetry_info = (/0_c_int, 0_c_int, 0_c_int, 0_c_int/)
!      
      host_input = molecule%pcmsolver_input()
!      
      if(trim(molecule%input).eq.'internal') then 
!      
         molecule%pcm_context = pcmsolver_new(PCMSOLVER_READER_HOST, &
                                              int(n_atoms_qm, kind=c_int), qm_charges, qm_coordinates, &
                                              symmetry_info, host_input, &
                                              c_funloc(host_writer))
!      
      else
!      
         molecule%pcm_context = pcmsolver_new(PCMSOLVER_READER_OWN, &
                                              int(n_atoms_qm, kind=c_int), qm_charges, qm_coordinates, &
                                              symmetry_info, host_input, &
                                              c_funloc(host_writer))
!      
      endif
!                                  
      molecule%n_tesserae = pcmsolver_get_cavity_size(molecule%pcm_context)
!               
      call mem%alloc(molecule%grid_coord, 3*molecule%n_tesserae)
!      
      molecule%grid_coord = zero
!      
      call pcmsolver_get_centers(molecule%pcm_context, molecule%grid_coord)
!
      call mem%alloc(molecule%charges, molecule%n_tesserae)
      call mem%alloc(molecule%pcm_rhs, molecule%n_tesserae)
!
      molecule%charges = zero
! 
      call pcmsolver_print(molecule%pcm_context)
!
   end subroutine prepare_pcm
!
!
   subroutine cleanup_pcm(molecule)
!!
!!    Clean Up PCM Molecular System
!!    Written by Tommaso Giovannini, Oct 2019
!!
      implicit none
!
      class(pcm) :: molecule
!
      if (allocated(molecule%pcm_rhs))    call mem%dealloc(molecule%pcm_rhs, molecule%n_tesserae)
!      
      if (allocated(molecule%charges))    call mem%dealloc(molecule%charges, molecule%n_tesserae)
!      
      if (allocated(molecule%grid_coord)) call mem%dealloc(molecule%grid_coord, 3*molecule%n_tesserae)
!
!
   end subroutine cleanup_pcm
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
!
   end subroutine host_writer
!
!
   function pcmsolver_input_pcm(molecule) result(host_input)
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
     class(pcm) :: molecule
     type(PCMInput) :: host_input
!
!    These parameters would be set by the host input reading subroutine(s)
!    Notice that the strings have a maximum pre-set length.
!    This is to ensure C interoperability
!    Length and area parameters are all assumed to be in Angstrom,
!    the module will convert to Bohr internally
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
!    Adjusted dimensions for PCMSolver interface
!
     character(kind=c_char, len=21) :: pcmmod_solvent
     character(kind=c_char, len=7) :: pcmmod_solver_type
!
!    
     pcmmod_solver_type = trim(molecule%solver_type)//c_null_char
     pcmmod_solvent     = trim(molecule%solvent)//c_null_char
!
     host_input%cavity_type     = pcmsolver_fstring_to_carray(pcmmod_cavity_type)
     host_input%patch_level     = int(pcmmod_patch_level, kind=c_int)
     host_input%coarsity        = pcmmod_coarsity
!     
     host_input%area            = molecule%tesserae_area
     host_input%min_distance    = pcmmod_min_distance
     host_input%der_order       = int(pcmmod_der_order, kind=c_int)
     host_input%scaling         = pcmmod_scaling
     host_input%radii_set       = pcmsolver_fstring_to_carray(pcmmod_radii_set)
     host_input%restart_name    = pcmsolver_fstring_to_carray(pcmmod_restart_name)
     host_input%min_radius      = pcmmod_min_radius
!
     host_input%solver_type     = pcmsolver_fstring_to_carray(pcmmod_solver_type)
!
     host_input%solvent         = pcmsolver_fstring_to_carray(pcmmod_solvent)
!
     host_input%equation_type   = pcmsolver_fstring_to_carray(pcmmod_equation_type)
     host_input%correction      = pcmmod_correction
     host_input%probe_radius    = pcmmod_probe_radius
!
     host_input%inside_type     = pcmsolver_fstring_to_carray(pcmmod_inside_type)
     host_input%outside_epsilon = pcmmod_outside_epsilon
     host_input%outside_type    = pcmsolver_fstring_to_carray(pcmmod_outside_type)
!
!
   end function pcmsolver_input_pcm
!
!
   subroutine read_parameters_pcm(molecule)
!!
!!    Read parameters PCM
!!    Written by Tommaso Giovannini, October 2019
!!
      implicit none
!
      class(pcm) :: molecule
      logical :: logical_exists
!
      molecule%input            = 'internal'
      molecule%solvent          = ' ' 
      molecule%tesserae_area    = 0.3d0 ! in Angstrom
      molecule%solver_type      = 'iefpcm' ! default IEF
!      
      call input%get_keyword_in_section('input','pcm',molecule%input)
!      
      if (trim(molecule%input) .eq. 'internal') then 
!      
         call input%get_keyword_in_section('solvent','pcm',molecule%solvent)
!      
         call input%get_keyword_in_section('tesserae area','pcm',molecule%tesserae_area)
!      
         call input%get_keyword_in_section('solver type','pcm',molecule%solver_type)
!         
         if(len_trim(molecule%solvent).eq.0) &
            call output%error_msg('PCM Calculation Internal Input: Solvent Not Specified')
!      
      else if(trim(molecule%input) .eq. 'external') then
!      
         inquire(file='@pcmsolver.inp',exist=logical_exists)
!         
         if(.not.logical_exists) &
            call output%error_msg('PCM Calculation External Input: Problems with go_pcm.py')
!      
      else
!      
         call output%error_msg('PCM Calculation: neither internal or external input is set')
!      
      endif
!
!
   end subroutine read_parameters_pcm
!
!
   subroutine print_description_and_settings_pcm(molecule)
!!
!!    Print system PCM
!!    Written by Tommaso Giovannini, Oct 2019
!!
      implicit none 
!
      class(pcm) :: molecule  
!
      call output%printf('m',  'Polarizable Continuum Solver via PCMSolver', fs='(/t6,a)')
      call output%printf('m', 'For details on PCM, see:', fs='(/t6,a)')
      call output%printf('m', 'Tomasi, Mennucci, Cammi, Chem. Rev. 2005, 105, 2999-3094.', &
                         fs='(t6,a)')
      call output%printf('m', 'For details on PCMSolver, see:', fs='(/t6,a)')
      call output%printf('m', 'Di Remigio et al., IJQC, 2019, 119, e25685', fs='(t6,a)')
!      
!
      if(molecule%input.eq.'external') then 
!
         call output%printf('m', 'PCM Solver was set via external file', fs='(/t6,a)')
!
      else 
!
         call output%printf('m', 'PCM Solver was set via internal parameters', fs='(/t6,a)')
         call output%printf('m', 'Solver Type:   '//trim(molecule%solver_type), fs='(t6,a)')
         call output%printf('m', 'Solvent:       '//trim(molecule%solvent), fs='(t6,a)')
         call output%printf('m', 'Tesserae Area: (f5.3) Å', &
                            reals=[molecule%tesserae_area], fs='(t6,a/)')
!
      endif
!
!
      if(input%requested_cc_calculation()) then 
         call output%printf('m', 'CC calculation: zero-order approximation', fs='(/t6,a)')
         call output%printf('m', 'PCM charges only affect MOs and Fock', fs='(t6,a/)')
      endif
!      
!
   end subroutine print_description_and_settings_pcm
!
!
   subroutine set_surface_function_pcm(molecule, mep_lbl)
!!
!!    Set surface function
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
!!    Wrapper for pcmsolvers set_surface_function
!!
      implicit none 
!
      class(pcm) :: molecule  
      character(len=*), intent(in) :: mep_lbl
!
      call pcmsolver_set_surface_function(molecule%pcm_context, &
                                          int(molecule%n_tesserae, kind=c_int), &
                                          -molecule%pcm_rhs, &
                                          pcmsolver_fstring_to_carray(mep_lbl))
!                                          
   end subroutine set_surface_function_pcm
!
!
   subroutine get_surface_function_pcm(molecule, asc_lbl)
!!
!!    Set surface function
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
!!    Wrapper for pcmsolvers set_surface_function
!!
      implicit none 
!
      class(pcm) :: molecule  
      character(len=*), intent(in) :: asc_lbl
!
      call pcmsolver_get_surface_function(molecule%pcm_context, &
                                          int(molecule%n_tesserae, kind=c_int), &
                                          molecule%charges, &
                                          pcmsolver_fstring_to_carray(asc_lbl))
!
   end subroutine get_surface_function_pcm
!
!
   subroutine compute_asc_pcm(molecule, mep_lbl, asc_lbl)
!!
!!    Compute ASC
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
!!    Wrapper for pcmsolvers compute_asc
!!
      implicit none 
!
      class(pcm) :: molecule  
      character(len=*), intent(in) :: mep_lbl
      character(len=*), intent(in) :: asc_lbl
!
      call pcmsolver_compute_asc(molecule%pcm_context, &
                                 pcmsolver_fstring_to_carray(mep_lbl), &
                                 pcmsolver_fstring_to_carray(asc_lbl), &
                                 irrep=0_c_int)
!                                 
   end subroutine compute_asc_pcm
!
#else
!
!  Dummy type to scream errors at you if you compile without PCMSolver
!  but call pcm in input
!
   implicit none
!
   type :: pcm
!
      real(dp), dimension(:), allocatable   :: charges
      real(dp), dimension(:), allocatable   :: pcm_rhs
      real(dp), dimension(:), allocatable   :: grid_coord
      integer :: n_tesserae             
!      
   contains
!
      procedure :: prepare                         => prepare_pcm
      procedure :: cleanup                         => cleanup_pcm
      procedure :: print_description_and_settings  => print_description_and_settings_pcm
!
      procedure :: set_surface_function            => set_surface_function_pcm
      procedure :: get_surface_function            => get_surface_function_pcm
      procedure :: compute_asc                     => compute_asc_pcm
!
   end type pcm
!
contains
!
!
   subroutine prepare_pcm(molecule, n_atoms_qm, qm_coordinates, qm_charges)
!!
!!    Prepare PCM dummy
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
      implicit none
!
      class(pcm) :: molecule
!
      integer, intent(in)                             :: n_atoms_qm
      real(dp), dimension(3*n_atoms_qm), intent(out)  :: qm_coordinates
      real(dp), dimension(n_atoms_qm), intent(out)    :: qm_charges
!
      call output%error_msg('prepare_pcm called, but eT is not compiled with pcm. &
                            &Use --pcm when running setup to link to PCMSolver.')
!
!     Some assignments to keep the compiler from complaining
      qm_coordinates = zero
      qm_charges = zero
      molecule%n_tesserae = 0
!
   end subroutine prepare_pcm
!
!
   subroutine cleanup_pcm(molecule)
!!
!!    Clean Up dummy
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
      implicit none
!
      class(pcm) :: molecule
!
      call output%error_msg('cleanup_pcm called, but eT is not compiled with pcm. &
                            & Use --pcm when running setup to link to PCMSolver.')
!
!     Some assignment to keep the compiler from complaining
      molecule%n_tesserae = 0
!
   end subroutine cleanup_pcm
!
!
   subroutine print_description_and_settings_pcm(molecule)
!!
!!    Print description and settings dummy
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
      implicit none 
!
      class(pcm) :: molecule  
!
      call output%error_msg('print_description_and_settings_pcm called, &
                            &but eT is not compiled with pcm. &
                            &Use --pcm when running setup to link to PCMSolver.')
!
!     Some assignment to keep the compiler from complaining
      molecule%n_tesserae = 0
!
   end subroutine print_description_and_settings_pcm
!
!
   subroutine set_surface_function_pcm(molecule, mep_lbl)
!!
!!    Set surface function dummy
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
      implicit none 
!
      class(pcm) :: molecule  
      character(len=*), intent(in) :: mep_lbl
!
      call output%error_msg('set_surface_function_pcm called with (a0), &
                            &but eT is not compiled with pcm. &
                            &Use --pcm when running setup to link to PCMSolver.', &
                             chars=[mep_lbl])
!
!     Some assignment to keep the compiler from complaining
      molecule%n_tesserae = 0
!                                          
   end subroutine set_surface_function_pcm
!
!
   subroutine get_surface_function_pcm(molecule, asc_lbl)
!!
!!    Get surface function dummy
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
      implicit none 
!
      class(pcm) :: molecule  
      character(len=*), intent(in) :: asc_lbl
!
      call output%error_msg('get_surface_function_pcm called with (a0), &
                            &but eT is not compiled with pcm. &
                            &Use --pcm when running setup to link to PCMSolver.', &
                             chars=[asc_lbl])
!
!     Some assignment to keep the compiler from complaining
      molecule%n_tesserae = 0
!                                          
   end subroutine get_surface_function_pcm
!
!
   subroutine compute_asc_pcm(molecule, mep_lbl, asc_lbl)
!!
!!    Compute ASC dummy
!!    Written by Rolf Heilemann Myhre, Dec 2019
!!
      implicit none 
!
      class(pcm) :: molecule  
      character(len=*), intent(in) :: mep_lbl
      character(len=*), intent(in) :: asc_lbl
!
      call output%error_msg('get_surface_function_pcm called with (a0) and (a0), &
                            &but eT is not compiled with pcm. &
                            &Use --pcm when running setup to link to PCMSolver.', &
                             chars=[mep_lbl, asc_lbl])
!
!     Some assignment to keep the compiler from complaining
      molecule%n_tesserae = 0
!                                          
   end subroutine compute_asc_pcm
!
#endif
!
end module pcm_class
