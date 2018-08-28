module hf_solver_class
!!
!!    Abstract HF solver class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use kinds 
   use parameters
   use hf_class
!
   implicit none 
!
   type, abstract :: hf_solver 
!
      real(dp) :: energy_threshold   = 1.0D-6
      real(dp) :: residual_threshold = 1.0D-6
!
      integer(i15) :: max_iterations = 100
!
      real(dp) :: coulomb_thr        = 1.0D-11 ! screening threshold 
      real(dp) :: exchange_thr       = 1.0D-11 ! screening threshold 
      real(dp) :: coulomb_precision  = 1.0D-14 ! integral engine accuracy
!
   contains 
!
      procedure(essential), deferred :: prepare 
      procedure(essential), deferred :: run 
      procedure(essential), deferred :: cleanup
!
      procedure :: read_settings           => read_settings_hf_solver
      procedure :: read_hf_solver_settings => read_hf_solver_settings_hf_solver
!
   end type hf_solver
!
   abstract interface
!
      subroutine essential(solver, wf)
!
         import :: hf, hf_solver
!
         implicit none 
!
         class(hf_solver) :: solver 
!
         class(hf) :: wf 
!
      end subroutine essential
!
   end interface
!
contains 
!
!
   subroutine read_settings_hf_solver(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings. This routine is to be overwritten by 
!!    descendants if more settings need to be set. 
!!
      implicit none 
!
      class(hf_solver) :: solver 
!
      call solver%read_hf_solver_settings()
!
   end subroutine read_settings_hf_solver
!
!
   subroutine read_hf_solver_settings_hf_solver(solver)
!!
!!    Read HF solver settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
!!    Reads the settings specific to this class.
!!
      implicit none 
!
      class(hf_solver) :: solver 
!
      integer(i15) :: n_records, i 
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
            write(output%unit, *) trim(line)
!
            if (line(1:17) == 'energy_threshold:') then
!
               value = line(18:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%energy_threshold
               return
!
            elseif (line(1:19) == 'residual_threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%residual_threshold
               return
!
            elseif (line(1:15) == 'max_iterations:') then 
!
               value = line(16:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%max_iterations
               return
!
            elseif (line(1:18) == 'coulomb_threshold:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%coulomb_thr
               return
!
            elseif (line(1:19) == 'exchange_threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%exchange_thr
               return
!
            elseif (line(1:18) == 'coulomb_precision:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%coulomb_precision
               return
!
            endif
!
         enddo
!
      endif 
!
   end subroutine read_hf_solver_settings_hf_solver
!
!
end module hf_solver_class
