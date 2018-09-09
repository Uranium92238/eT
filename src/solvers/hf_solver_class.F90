module hf_solver_class
!!
!!    Abstract HF solver class module 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use kinds 
   use parameters
   use hf_class
   use uhf_class
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
      procedure :: prepare                 => prepare_hf_solver
      procedure :: run                     => run_hf_solver
      procedure :: cleanup                 => cleanup_hf_solver
!
      procedure :: print_banner            => print_banner_hf_solver
!
      procedure :: read_settings           => read_settings_hf_solver
      procedure :: read_hf_solver_settings => read_hf_solver_settings_hf_solver
!
   end type hf_solver
!
contains 
!
!
   subroutine prepare_hf_solver(solver, wf)
!!
!!    Prepare 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_solver) :: solver
!
      class(hf) :: wf
!
!     Read settings (thresholds, etc.)
!
      call solver%read_settings()
!
   end subroutine prepare_hf_solver
!
!
   subroutine run_hf_solver(solver, wf)
!!
!!    Run 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_solver) :: solver
!
      class(hf) :: wf
!
!     Nothing here 
!
   end subroutine run_hf_solver
!
!
   subroutine cleanup_hf_solver(solver, wf)
!!
!!    Cleanup 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(hf_solver) :: solver
!
      class(hf) :: wf
!
!     Nothing here 
!
   end subroutine cleanup_hf_solver
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
            if (line(1:17) == 'energy threshold:') then
!
               value = line(18:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%energy_threshold
               cycle
!
            elseif (line(1:19) == 'residual threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%residual_threshold
               cycle
!
            elseif (line(1:15) == 'max iterations:') then 
!
               value = line(16:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%max_iterations
               cycle
!
            elseif (line(1:18) == 'coulomb threshold:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%coulomb_thr
               cycle
!
            elseif (line(1:19) == 'exchange threshold:') then 
!
               value = line(20:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%exchange_thr
               cycle
!
            elseif (line(1:18) == 'coulomb precision:') then 
!
               value = line(19:100)
               value = remove_preceding_blanks(value)
               read(value, *) solver%coulomb_precision
               cycle
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
   subroutine print_banner_hf_solver(solver)
!!
!!    Print banner
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none 
!
      class(hf_solver) :: solver 
!
   end subroutine print_banner_hf_solver
!
!
end module hf_solver_class
