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
module eT_class
!
!!
!!    eT class
!!    Written by Eirik F. Kjønstad, May 2022
!!
!
   use timings_class, only: timings
!
   implicit none
!
   type :: eT
!
      type(timings), allocatable, private :: timer
!
   contains
!
      procedure, public :: run
!
      procedure, private :: initializations_for_input_output_and_timing
      procedure, private, nopass :: create_memory_manager
!
      procedure, private, nopass :: get_and_print_n_threads
      procedure, private, nopass :: set_global_print_levels_in_output_and_timing
      procedure, private, nopass :: print_compilation_info
!
      procedure, private :: print_timestamp
      procedure, private, nopass :: get_date_and_time
!
      procedure, private :: run_calculations
!
      procedure, private, nopass :: cholesky_decompose_eris
      procedure, private, nopass :: run_reference_calculation
      procedure, private, nopass :: run_cc_calculation
      procedure, private, nopass :: run_fci_calculation
!
      procedure, private :: print_top_info_to_output_and_timing
      procedure, private :: print_bottom_info_to_output
      procedure, private :: print_total_cpu_and_wall_time
!
   end type eT
!
!
   interface eT
!
      procedure :: new_eT
!
   end interface eT
!
!
contains
!
!
   function new_eT() result(this)
!!
!!    New
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use timings_class, only: timings
!
      implicit none
!
      type(eT) :: this
!
      this%timer = timings("Total time in eT", pl='minimal')
!
   end function new_eT
!
!
   subroutine run(this)
!!
!!    Run
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad,
!!    Alexander C. Paul, and Rolf H. Myhre, 2018-2022
!!
      use parameters
      use global_in
      use global_out
      use timings_class,            only: timings, timing
      use citation_printer_class,   only: eT_citations, citation_printer
      use timings_file_class,       only: timings_file
      use memory_manager_class,     only: mem
!
      implicit none
!
      class(eT), intent(inout) :: this
!
      call this%timer%turn_on()
!
      output = output_file('eT.out')
      timing = timings_file('eT.timing.out')
!
      call output%open_()
      call timing%open_()
!
      input = input_file('eT.inp')
!
      call this%initializations_for_input_output_and_timing()
!
      call this%print_top_info_to_output_and_timing()
!
      call input%check_for_errors()
!
      eT_citations = citation_printer(input)
      mem = this%create_memory_manager(input)
!
      call this%run_calculations()
!
      call mem%check_for_leak()
!
      call input%cleanup()
!
      call this%timer%turn_off()
!
      call output%check_for_warnings()
!
      call this%print_bottom_info_to_output()
!
      call timing%close_()
!
      call output%printf('m', 'eT terminated successfully!', fs='(/t3,a)')
      call output%close_()
!
   end subroutine run
!
!
   subroutine initializations_for_input_output_and_timing(this)
!!
!!    Initializations for input, output, and timing
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_in, only: input
!
      implicit none
!
      class(eT), intent(in) :: this
!
      call input%read_and_process_input()
      call this%set_global_print_levels_in_output_and_timing()
!
   end subroutine initializations_for_input_output_and_timing
!
!
   subroutine print_total_cpu_and_wall_time(this)
!!
!!    Print total CPU and wall time
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_out, only: output
!
      implicit none
!
      class(eT), intent(in) :: this
!
      call output%printf('m', 'Total wall time in eT (sec): (f20.5)', &
                         reals=[this%timer%get_elapsed_time('wall')], fs='(/t3,a)')
!
      call output%printf('m', 'Total cpu time in eT (sec):  (f20.5)', &
                         reals=[this%timer%get_elapsed_time('cpu')], fs='(t3,a)')
!
   end subroutine print_total_cpu_and_wall_time
!
!
   subroutine print_timestamp(this, print_label)
!!
!!    Print timestamp
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_out, only: output
!
      implicit none
!
      class(eT), intent(in) :: this
!
      character(len=*), intent(in) :: print_label
!
      character(len=50) :: timestamp
!
      call this%get_date_and_time(timestamp)
      call output%printf('m', 'Calculation (a0):', chars=[print_label], fs='(/t3,a)', adv=.false.)
      call output%printf('m', "(a0)", chars=[timestamp], fs='(t1,a)')
!
   end subroutine print_timestamp
!
!
   subroutine get_and_print_n_threads()
!!
!!    Get and print number of threads
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use omp_lib
      use global_out, only: output
!
      implicit none
!
      integer :: n_threads
!
      n_threads = 1
!
!$    n_threads = omp_get_max_threads()
!
      if (n_threads .eq. 1) then
!
         call output%printf('m', 'Running on (i0) OMP thread', ints=[n_threads], fs='(/t3,a)')
!
      else
!
         call output%printf('m', 'Running on (i0) OMP threads', ints=[n_threads], fs='(/t3,a)')
!
      endif
!
   end subroutine get_and_print_n_threads
!
!
   function create_memory_manager(input) result(mem_manager)
!!
!!    Create memory manager
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use parameters
      use memory_manager_class, only: memory_manager
      use input_file_class, only: input_file
!
      implicit none
!
      class(input_file), intent(in) :: input
!
      character(len=200) :: mem_unit
      integer(i64)       :: mem_total
!
      type(memory_manager) :: mem_manager
!
      mem_total = 8
      mem_unit  = 'gb'
!
      call input%get_keyword('available', 'memory', mem_total)
      call input%get_keyword('unit',      'memory', mem_unit)
!
      mem_manager = memory_manager(total = mem_total, &
                                   units = mem_unit)
!
   end function create_memory_manager
!
!
   subroutine run_calculations(this)
!!
!!    Run calculations
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use global_out, only: output
      use global_in, only: input
!
      use hf_class, only: hf
      use libint_initialization, only: initialize_libint_c, finalize_libint_c
!
      implicit none
!
      class(eT), intent(in) :: this
!
      class(hf), allocatable :: ref_wf
!
      logical :: requested_cholesky
!
      call initialize_libint_c() ! Safe to use Libint from now on
!
!     Cholesky decomposition of electron repulsion integrals (ERIs)
      requested_cholesky = input%is_keyword_present('cholesky eri', 'do')
      if (requested_cholesky) call cholesky_decompose_eris()
!
      if (input%requested_reference_calculation()) then
!
         call this%run_reference_calculation(ref_wf)
!
         if (input%requested_cc_calculation()) then
!
            call ref_wf%prepare_for_post_HF_method()
            call this%run_cc_calculation(ref_wf)
!
         else if (input%requested_fci_calculation()) then
!
            call ref_wf%prepare_for_post_HF_method()
            call this%run_fci_calculation(ref_wf)
!
         endif
!
         call ref_wf%cleanup()
         deallocate(ref_wf)
!
      else
!
         if (input%requested_cc_calculation()) &
            call output%error_msg('to run CC calculation reference wavefunction must be specified.')
!
         if (.not. requested_cholesky) &
            call output%error_msg('no method nor ERI Cholesky decomposition selected in input.')
!
      endif
!
      call finalize_libint_c() ! No longer safe to use Libint
!
   end subroutine run_calculations
!
!
   subroutine run_reference_calculation(ref_wf)
!!
!!    Run reference calculation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      use hf_class, only: hf
      use hf_engine_class, only: hf_engine
      use reference_wavefunction_factory_class, only: reference_wavefunction_factory
      use reference_engine_factory_class, only: reference_engine_factory
!
      implicit none
!
      class(hf), allocatable, intent(inout)              :: ref_wf
      class(reference_wavefunction_factory), allocatable :: ref_wf_factory
!
      class(hf_engine), allocatable  :: ref_engine
      type(reference_engine_factory) :: ref_engine_factory
!
      ref_wf_factory = reference_wavefunction_factory()
!
      call ref_wf_factory%create(ref_wf)
!
      call ref_engine_factory%create(ref_engine)
!
      call ref_engine%ignite(ref_wf)
!
   end subroutine run_reference_calculation
!
!
   subroutine run_cc_calculation(ref_wf)
!!
!!    Run coupled cluster calculation
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
      use hf_class, only: hf
!
      use ccs_class,                      only: ccs
      use cc_wavefunction_factory_class,  only: cc_wavefunction_factory
!
      use cc_engine_class,          only: cc_engine
      use cc_engine_factory_class,  only: cc_engine_factory
!
      implicit none
!
      class(hf), intent(in) :: ref_wf
!
      class(ccs), allocatable                    :: cc_wf
      type(cc_wavefunction_factory), allocatable :: cc_wf_factory
!
      class(cc_engine), allocatable        :: engine
      type(cc_engine_factory), allocatable :: engine_factory
!
      cc_wf_factory = cc_wavefunction_factory()
      call cc_wf_factory%create(ref_wf, cc_wf)
!
      allocate(cc_engine_factory::engine_factory)
      call engine_factory%create(engine)
!
      call engine%ignite(cc_wf)
!
      call cc_wf%cleanup()
!
   end subroutine run_cc_calculation
!
!
   subroutine run_fci_calculation(ref_wf)
!!
!!    Run FCI calculation
!!    Written by Enrico Ronca, 2020
!!
      use hf_class, only: hf
!
      use fci_class,                      only: fci
      use fci_wavefunction_factory_class, only: fci_wavefunction_factory
!
      use fci_engine_class, only: fci_engine
!
      implicit none
!
      class(hf), intent(in) :: ref_wf
!
      class(fci), allocatable                     :: fci_wf
      type(fci_wavefunction_factory), allocatable :: fci_wf_factory
!
      class(fci_engine), allocatable :: engine
!
      fci_wf_factory = fci_wavefunction_factory()
      call fci_wf_factory%create(ref_wf, fci_wf)
!
      engine = fci_engine()
!
      call engine%ignite(fci_wf)
!
      call fci_wf%cleanup()
!
   end subroutine run_fci_calculation
!
!
   subroutine cholesky_decompose_eris()
!!
!!    Cholesky decompose ERIs
!!    Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019 and Dec 2019
!!
!!    Performs Cholesky decomposition of the electron repulsion integral matrix.
!!
      use eri_cd_class,  only: eri_cd
      use ao_tool_class, only: ao_tool
!
      implicit none
!
      type(eri_cd), allocatable  :: eri_cholesky_solver
      type(ao_tool), allocatable :: ao
!
      ao = ao_tool()
      call ao%initialize()
!
      eri_cholesky_solver = eri_cd(ao)
!
      call eri_cholesky_solver%run(ao)
!
      call eri_cholesky_solver%diagonal_test(ao)  ! Determine the largest
                                                  ! deviation in the ERI matrix
!
      call eri_cholesky_solver%cleanup()
!
   end subroutine cholesky_decompose_eris
!
!
   subroutine set_global_print_levels_in_output_and_timing()
!!
!!    Set global print levels in output and timing
!!    Written by Rolf H. Myhre, Oct. 2019
!!
      use global_out,      only: output
      use global_in,       only: input
      use timings_class,   only: timing
!
      implicit none
!
      character(len=200) :: print_level
!
!     Set default
      print_level = 'normal'
!
!     Overwrite print_level if keyword is present
      call input%get_keyword('output print level', 'print', print_level)
!
!     This is the only place this routine is allowed to be called
      call output%set_global_print_level(print_level)
!
!     Repeat for timing file
!     Set default
      print_level = 'normal'
!
!     Overwrite print_level if keyword is present
      call input%get_keyword('timing print level', 'print', print_level)
!
!     This is the only place this routine is allowed to be called
      call timing%set_global_print_level(print_level)
!
   end subroutine set_global_print_levels_in_output_and_timing
!
!
   subroutine print_top_info_to_output_and_timing(this)
!!
!!    Print to info to output and timing
!!    Written by Eirik F. Kjønstad, 2019
!!
      use parameters
      use global_out, only: output, timing
      use global_in, only: input
!
      implicit none
!
      class(eT), intent(in) :: this
!
      call output%printf('m', 'eT (i0).(i0) - an electronic structure program ', &
                      ints=[major_version, minor_version], fs='(///t22,a)')
!
      call output%print_separator('m',72,'-', fs='(/t3,a)')
!
      call output%printf('m', 'Author list in alphabetical order:', fs='(t4,a)')
!
      call output%print_separator('m',72,'-', fs='(t3,a)')
!
      call output%printf('m', 'J. H. Andersen, '         // &
                           'A. Balbi, '               // &
                           'S. Coriani, '             // &
                           'S. D. Folkestad, '        // &
                           'T. Giovannini, '          // &
                           'L. Goletto, '             // &
                           'T. S. Haugland, '         // &
                           'A. Hutcheson, '           // &
                           'I-M. Høyvik, '            // &
                           'E. F. Kjønstad, '         // &
                           'H. Koch, '                // &
                           'R. Matveeva, '            // &
                           'T. Moitra, '              // &
                           'R. H. Myhre, '            // &
                           'A. C. Paul, '             // &
                           'S. Roet, '                // &
                           'E. Ronca, '               // &
                           'M. Scavino, '             // &
                           'A. K. Schnack-Petersen, ' // &
                           'A. S. Skeidsvoll, '       // &
                           'Å. H. Tveten',               &
                           ffs='(t4,a)', fs='(t4,a)', ll=68)
!
      call output%print_separator('m',72,'-', fs='(t3,a)')
!
      call output%printf('m', 'J. Chem. Phys. 152, 184103 (2020); https://doi.org/10.1063/5.0004713', &
                      fs='(t4,a)')
!
      call output%printf('m', "This is eT (i0).(i0).(i0) (a0)", &
                         ints=[major_version, minor_version, patch_version], &
                         chars = [version_name], fs='(//t4,a)')
!
      call this%print_compilation_info(output)
      call this%print_compilation_info(timing)
!
      call timing%print_banner()
!
      call input%print_input_except_geometry()
!
      call this%print_timestamp(print_label = 'start')
!
      call this%get_and_print_n_threads()
!
   end subroutine print_top_info_to_output_and_timing
!
!
   subroutine print_bottom_info_to_output(this)
!!
!!    Print bottom info to output
!!    Written by Eirik F. Kjønstad, May 2022
!!
      use memory_manager_class, only: mem
      use citation_printer_class, only: eT_citations
      use global_out, only: output
!
      implicit none
!
      class(eT), intent(in) :: this
!
      call mem%print_max_used()
!
      call this%print_total_cpu_and_wall_time()
!
      call this%print_timestamp(print_label = 'end')
!
      call eT_citations%print_(output)
!
   end subroutine print_bottom_info_to_output
!
!
   subroutine get_date_and_time(string)
!!
!!    Get date and time
!!    Written by Rolf H. Myhre, Nov, 2020
!!
!!    Returns a formatted string with date, time and UTC offset
!!
!!    Format: yyyy-mm-dd hh:mm:ss UTC xhh:mm
!!
      implicit none
!
      character(len=*), intent(out) :: string
      character(len=20) :: date, time, zone
!
      call date_and_time(date=date, time=time, zone=zone)
!
      string = ""
      write(string, "(a,a,a,a,a,a)") date(1:4), "-", date(5:6), "-", date(7:8), " "
      write(string(12:), "(a,a,a,a,a)") time(1:2), ":", time(3:4), ":", time(5:6)
      write(string(20:), "(a,a,a,a)") " UTC ", zone(1:3), ":", zone(4:5)
!
   end subroutine get_date_and_time
!
!
   subroutine print_compilation_info(file_)
!!
!!    Print compilation info
!!    Written by Rolf H. Myhre, Nov, 2020
!!
!!    Retrieves compilation information from the
!!    get_compilation_info library generated by CMake
!!    and prints to output
!!
      use output_file_class, only: output_file
!
      implicit none
!
      class(output_file), intent(inout) :: file_
!
      character(len=200) :: string
!
      call file_%print_separator('n',60,'-', fs='(t3,a)')
!
      call get_configuration_time(string)
      call file_%printf("m", "Configuration date: (a0)", chars =[string])
!
      call get_git_branch(string)
      if(len_trim(string) .gt. 0) then
         call file_%printf("m", "Git branch:         (a0)", chars =[string])

         call get_git_hash(string)
         call file_%printf("m", "Git hash:           (a0)", chars =[string])
      endif
!
      call get_fortran_compiler(string)
      call file_%printf("m", "Fortran compiler:   (a0)", chars =[string])
!
      call get_c_compiler(string)
      call file_%printf("m", "C compiler:         (a0)", chars =[string])
!
      call get_cxx_compiler(string)
      call file_%printf("m", "C++ compiler:       (a0)", chars =[string])
!
      call get_lapack_type(string)
      call file_%printf("m", "LAPACK type:        (a0)", chars =[string])
!
      call get_blas_type(string)
      call file_%printf("m", "BLAS type:          (a0)", chars =[string])
!
      call get_int64(string)
      call file_%printf("m", "64-bit integers:    (a0)", chars =[string])
!
      call get_omp(string)
      call file_%printf("m", "OpenMP:             (a0)", chars =[string])
!
      call get_pcm(string)
      call file_%printf("m", "PCM:                (a0)", chars =[string])
!
      call get_forced_batching(string)
      call file_%printf("m", "Forced batching:    (a0)", chars =[string])
!
      call get_runtime_check(string)
      call file_%printf("m", "Runtime checks:     (a0)", chars =[string])
!
      call file_%print_separator("m",60,"-", fs="(t3,a)")
!
   end subroutine print_compilation_info
!
!
end module eT_class
