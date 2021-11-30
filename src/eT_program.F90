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
program eT_program
!
!!
!!  eT - a coupled cluster program
!!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2019
!!
!
   use parameters
   use global_in
   use global_out
   use timings_class,            only: timings, timing
   use memory_manager_class,     only: mem, memory_manager
   use citation_printer_class,   only: eT_citations, citation_printer
   use libint_initialization,    only: initialize_libint_c, finalize_libint_c
   use timings_file_class,       only: timings_file
!
   use hf_class,     only: hf
   use uhf_class,    only: uhf
   use mlhf_class,   only: mlhf
   use cuhf_class,   only: cuhf
   use rohf_class,   only: rohf
!
   use omp_lib
!
   implicit none
!
   integer :: n_threads
!
   type(timings) :: eT_timer
!
   class(hf), allocatable  :: ref_wf
!
!  Memory manager information
!
   character(len=200) :: mem_unit
   integer(i64) :: mem_total
!
!  Interface reference and CC wavefunction calculation,
!  as well as for starting and stopping Libint
!
   character(len=50) :: timestamp
!
   logical :: requested_cholesky
!
!  Interface reference and CC wavefunction calculation
!
   interface
!
      subroutine reference_calculation(ref_wf)
!
         use hf_class, only: hf
!
         implicit none
!
         class(hf), allocatable, intent(inout)  :: ref_wf
!
      end subroutine reference_calculation
!
      subroutine cc_calculation(ref_wf)
!
         use hf_class, only: hf
!
         implicit none
!
         class(hf), intent(in)  :: ref_wf
!
      end subroutine cc_calculation
!
      subroutine print_compilation_info(file_)
!
         use output_file_class, only: output_file
!
         implicit none
!
         class(output_file), intent(inout) :: file_
!
      end subroutine print_compilation_info
!
   end interface
!
!  Prepare input, output and timing file
!
   output = output_file('eT.out')
   call output%open_()
!
   input = input_file('eT.inp')
   call input%open_()
!
   timing = timings_file('eT.timing.out')
   call timing%open_()
!
   eT_timer = timings("Total time in eT", pl='minimal')
   call eT_timer%turn_on()
!
   call print_program_banner()
!
   call output%printf('m', "This is eT (i0).(i0).(i0) (a0)", &
                      ints=[major_version, minor_version, patch_version], &
                      chars = [version_name], fs='(//t4,a)')
!
   call print_compilation_info(output)
   call print_compilation_info(timing)
!
   call timing%print_banner()
!
   call get_date_and_time(timestamp)
   call output%printf('m', "Calculation started: (a0)", chars=[timestamp], fs='(/t3,a)')
!
   call input%read_keywords_and_geometry()
   call input%close_()
!
   eT_citations = citation_printer()
!
   n_threads = 1
!
!$   n_threads = omp_get_max_threads()
!
   if (n_threads .eq. 1) then
      call output%printf('m', 'Running on (i0) OMP thread', ints=[n_threads], fs='(/t3,a)')
   else
      call output%printf('m', 'Running on (i0) OMP threads', ints=[n_threads], fs='(/t3,a)')
   endif
!
!  Set print level in output and timing files
!
   call set_global_print_levels()
!
!  Create memory manager
!
   mem_total = 8
   mem_unit  = 'gb'
!
   call input%get_keyword('available', 'memory', mem_total)
   call input%get_keyword('unit', 'memory', mem_unit)
!
   mem = memory_manager(total=mem_total, &
                        units=mem_unit)
!
   call initialize_libint_c() ! Safe to use Libint from now on
!
!  Cholesky decomposition of electron repulsion integrals (ERIs)
!
   requested_cholesky = input%is_keyword_present('cholesky eri', 'do')
   if (requested_cholesky) call do_eri_cholesky()
!
!  Hartree-Fock calculation
!
   if (input%requested_reference_calculation()) then
!
      call reference_calculation(ref_wf)
!
!     Coupled cluster calculation
!
      if (input%requested_cc_calculation()) then
!
         call ref_wf%prepare_for_cc()
         call cc_calculation(ref_wf)
!
      endif
!
      call ref_wf%cleanup()
!
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
   call timing%printf('m', ":: Total time", fs='(//t3,a)')
   call timing%print_separator('m', 16, '=')
!
   call output%print_separator('m', 60, '-', fs='(/t3,a)')
!
   call mem%check_for_leak()
   call mem%print_max_used()
!
   call output%check_for_warnings()
!
   call eT_timer%turn_off()
   call output%printf('m', 'Total wall time in eT (sec): (f20.5)', &
                      reals=[eT_timer%get_elapsed_time('wall')], fs='(/t3,a)')
   call output%printf('m', 'Total cpu time in eT (sec):  (f20.5)', &
                      reals=[eT_timer%get_elapsed_time('cpu')], fs='(t3,a)')
!
   call get_date_and_time(timestamp)
   call output%printf('m', "Calculation ended: (a0)", chars=[timestamp], fs='(/t3,a)')
!
   call eT_citations%print_(output)
!
   call output%printf('m', 'eT terminated successfully!', fs='(/t3,a)')
!
   call output%close_()
   call timing%close_()
   call input%cleanup_geometry()
   call input%cleanup_keywords()
!
end program eT_program
!
!
subroutine reference_calculation(ref_wf)
!!
!! Reference calculation
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!! Directs the reference state calculation for eT
!!
   use global_in,  only: input
   use global_out, only: output
!
   use reference_engine_class, only: reference_engine
   use hf_geoopt_engine_class, only: hf_geoopt_engine
   use tdhf_engine_class,     only: tdhf_engine
!
   use hf_class,   only: hf
   use uhf_class,  only: uhf
   use mlhf_class, only: mlhf
   use cuhf_class, only: cuhf
   use rohf_class, only: rohf
!
   implicit none
!
   class(hf), allocatable, intent(inout)  :: ref_wf
!
   class(reference_engine), allocatable :: ref_engine
!
   character(len=30) :: ref_wf_name
!
   ref_wf_name = input%get_reference_wf()
!
   if (trim(ref_wf_name) == 'hf') then
!
      ref_wf = hf()
!
   elseif (trim(ref_wf_name) == 'uhf') then
!
      ref_wf = uhf()
!
   elseif (trim(ref_wf_name) == 'cuhf') then
!
      ref_wf = cuhf()
!
   elseif (trim(ref_wf_name) == 'rohf') then
!
      ref_wf = rohf()
!
   elseif (trim(ref_wf_name) == 'mlhf') then
!
      ref_wf = mlhf()
!
   else
!
      call output%error_msg('did not recognize the reference wavefunction ' &
                                    // trim(ref_wf_name) //'.')
!
   endif
!
   call ref_wf%prepare()
!
   if (input%is_keyword_present('ground state geoopt', 'do')) then
!
      ref_engine = hf_geoopt_engine()
!
   elseif (input%is_keyword_present('time dependent hf', 'do')) then
!
      ref_engine = tdhf_engine()
!
   else
!
      ref_engine = reference_engine()
!
   endif
!
   call ref_engine%ignite(ref_wf)
!
end subroutine reference_calculation
!
!
subroutine cc_calculation(ref_wf)
!!
!! Coupled cluster calculation
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!! Directs the coupled cluster calculation for eT
!!
   use global_in,  only: input
   use global_out, only: output
!
   use hf_class, only: hf
!
   use ccs_class, only: ccs
   use cc2_class, only: cc2
   use lowmem_cc2_class, only: lowmem_cc2
   use ccsd_class, only: ccsd
   use cc3_class, only: cc3
   use ccsdpt_class, only: ccsdpt
   use mp2_class, only: mp2
   use mlcc2_class, only: mlcc2
   use mlccsd_class, only: mlccsd
!
   use gs_engine_class, only: gs_engine
   use es_engine_class, only: es_engine
   use response_engine_class, only: response_engine
   use mean_value_engine_class, only: mean_value_engine
   use td_engine_class, only: td_engine
!
   implicit none
!
   class(hf), intent(in)  :: ref_wf
!
   class(ccs), allocatable :: cc_wf
   class(gs_engine), allocatable :: cc_engine
!
   character(len=30) :: cc_wf_name
!
   cc_wf_name = input%get_cc_wf()
!
   select case (trim(cc_wf_name))
!
      case ('ccs')
!
         allocate(ccs::cc_wf)
!
      case ('cc2')
!
         allocate(cc2::cc_wf)
!
      case ('lowmem-cc2')
!
         allocate(lowmem_cc2::cc_wf)
!
      case ('ccsd')
!
         allocate(ccsd::cc_wf)
!
      case ('cc3')
!
         allocate(cc3::cc_wf)
!
      case ('ccsd(t)')
!
         allocate(ccsdpt::cc_wf)
!
      case ('mp2')
!
         allocate(mp2::cc_wf)
!
      case ('mlcc2')
!
         allocate(mlcc2::cc_wf)
!
      case ('mlccsd')
!
         allocate(mlccsd::cc_wf)
!
      case default
!
         call output%error_msg('could not recognize CC method ' // trim(cc_wf_name) // '.')
!
   end select
!
!  initialize wavefunction
!
   call cc_wf%initialize(ref_wf)
!
   if (input%is_keyword_present('response', 'do')) then
!
      cc_engine = response_engine(cc_wf)
!
   elseif (input%is_keyword_present('excited state', 'do')) then
!
      cc_engine = es_engine(cc_wf)
!
   elseif (input%is_keyword_present('mean value', 'do')) then
!
      cc_engine = mean_value_engine(cc_wf)
!
   elseif (input%is_keyword_present('ground state', 'do')) then
!
      cc_engine = gs_engine(cc_wf)
!
   elseif (input%is_keyword_present('time dependent state', 'do')) then
!
      cc_engine = td_engine(cc_wf)
!
   else
!
      call output%error_msg('could not recognize coupled cluster task.')
!
   endif
!
   call cc_engine%ignite(cc_wf)
   call cc_wf%cleanup()
!
end subroutine cc_calculation
!
!
subroutine do_eri_cholesky()
!!
!! Do ERI Cholesky
!! Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019 and Dec 2019
!!
!! Performs Cholesky decomposition of the atomic orbital (AO) electron repulsion
!! integrals.
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
end subroutine do_eri_cholesky
!
!
subroutine set_global_print_levels()
!!
!! Set global print levels
!! Written by Rolf H. Myhre, Oct. 2019
!!
!! Reads and sets the global print levels for the output file
!! and the timing file from input.
!!
   use global_out, only: output
   use global_in, only: input
   use timings_class, only : timing
!
   implicit none
!
   character(len=200) :: print_level
!
!  Set default
   print_level = 'normal'
!
!  Overwrite print_level if keyword is present
   call input%get_keyword('output print level', 'print', print_level)
!
!  This is the only place this routine is allowed to be called
   call output%set_global_print_level(print_level)
!
!  Repeat for timing file
!  Set default
   print_level = 'normal'
!
!  Overwrite print_level if keyword is present
   call input%get_keyword('timing print level', 'print', print_level)
!
!  This is the only place this routine is allowed to be called
   call timing%set_global_print_level(print_level)
!
end subroutine set_global_print_levels
!
!
subroutine print_program_banner()
!!
!! Print program banner
!! Written by Eirik F. Kjønstad, 2019
!!
!! Prints banner, author list, and list of contributors.
!!
   use parameters
   use global_out, only: output
!
   implicit none
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
end subroutine print_program_banner
!
!
subroutine get_date_and_time(string)
!!
!! Get date and time
!! Written by Rolf H. Myhre, Nov, 2020
!!
!! Returns a formatted string with date, time and UTC offset
!!
!! Format: yyyy-mm-dd hh:mm:ss UTC xhh:mm
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
!! Print compilation info
!! Written by Rolf H. Myhre, Nov, 2020
!!
!! Retrieves compilation information from the
!! get_compilation_info library generated by CMake
!! and prints to output
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
