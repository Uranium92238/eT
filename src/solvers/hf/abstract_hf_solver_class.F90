module abstract_hf_solver_class
!!
!!    Abstract HF solver class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use kinds 
   use parameters   
   use file_class
   use disk_manager_class
   use io_utilities
!
   implicit none 
!
   type, abstract :: abstract_hf_solver 
!
      character(len=100) :: tag
      character(len=100) :: author
      character(len=400) :: description
!
      real(dp) :: energy_threshold          = 1.0D-6
      real(dp) :: gradient_threshold        = 1.0D-6
!
      integer :: max_iterations        = 100
!  
      character(len=40) :: ao_density_guess = 'SAD'
!
   contains 
!
      procedure :: print_banner             => print_banner_abstract_hf_solver
!
      procedure :: read_settings            => read_settings_abstract_hf_solver
      procedure :: read_hf_solver_settings  => read_hf_solver_settings_abstract_hf_solver
!
      procedure :: print_hf_solver_settings => print_hf_solver_settings_hf_solver
!
   end type abstract_hf_solver
!
contains 
!
!
   subroutine print_hf_solver_settings_hf_solver(solver)
!!
!!    Print HF solver settings    
!!    Written by Eirik F. Kjønstad, Sep 2018 
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      write(output%unit, '(t6,a30,e10.4)') 'Energy threshold:             ', solver%energy_threshold
      write(output%unit, '(t6,a30,e10.4)') 'Gradient threshold:           ', solver%gradient_threshold
!
   end subroutine print_hf_solver_settings_hf_solver
!
!
   subroutine read_settings_abstract_hf_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings. This routine is to be overwritten by 
!!    descendants if more settings need to be set. 
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      call solver%read_hf_solver_settings()
!
   end subroutine read_settings_abstract_hf_solver
!
!
   subroutine read_hf_solver_settings_abstract_hf_solver(solver)
!!
!!    Read HF solver settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings specific to this class.
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      integer :: n_records, i 
!
      character(len=100) :: line, value 
!
      if (requested_section('hf')) then ! User has requested something 
!
         call move_to_section('hf', n_records)
!
         do i = 1, n_records
!
            read(input%unit, '(a100)') line
            line = remove_preceding_blanks(line)
!
            if (line(1:17) == 'energy threshold:') then
!
               value = line(18:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%energy_threshold
               cycle
!
            elseif (line(1:19) == 'gradient threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%gradient_threshold
               cycle
!
            elseif (line(1:15) == 'max iterations:') then 
!
               value = line(16:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%max_iterations
               cycle
!
            elseif (line(1:17) == 'ao density guess:') then 
!
               value = line(18:100)
               value = remove_preceding_blanks(value)
               solver%ao_density_guess = trim(value)
               cycle
!
            endif
!
         enddo
!
      endif 
!
   end subroutine read_hf_solver_settings_abstract_hf_solver
!
!
   subroutine print_banner_abstract_hf_solver(solver)
!!
!!    Print banner
!!    Written by Rolf H. Myhre, 2018
!!
      implicit none 
!
      class(abstract_hf_solver) :: solver 
!
      call long_string_print(solver%tag,'(//t3,a)',.true.)
      call long_string_print(solver%author,'(t3,a/)',.true.)
      call long_string_print(solver%description)
!
   end subroutine print_banner_abstract_hf_solver
!
!
end module abstract_hf_solver_class
