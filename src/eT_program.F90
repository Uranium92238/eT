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
program eT_program
!
!!
!!  eT - a coupled cluster program
!!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2019
!!
!
   use kinds
   use global_in
   use global_out
   use timings_class, only : timings, timing
   use memory_manager_class, only : mem, memory_manager
   use libint_initialization, only : initialize_libint, finalize_libint
   use molecular_system_class, only : molecular_system
!
   use hf_class, only: hf
   use uhf_class, only: uhf
   use mlhf_class, only: mlhf
!
   use omp_lib
!
   implicit none
!
   integer :: n_threads
!
!  Molecular system object
!
   type(molecular_system) :: system
!
!  Timer object
!
   type(timings) :: eT_timer
!
!  Reference wavefunction
!
   class(hf), allocatable  :: ref_wf
!
!  Interface reference and CC wavefunction calculation
!
   interface
!
      subroutine reference_calculation(system, ref_wf)
!
         use molecular_system_class, only: molecular_system
         use hf_class, only: hf
!
         implicit none
!
         type(molecular_system), intent(inout)  :: system
         class(hf), allocatable, intent(inout)  :: ref_wf
!
      end subroutine reference_calculation
!
      subroutine cc_calculation(ref_wf)
!
         use molecular_system_class, only: molecular_system
!
         use hf_class, only: hf
!
         implicit none
!
         class(hf), intent(in)  :: ref_wf
!
      end subroutine cc_calculation
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
   timing = output_file('timing.out')
   call timing%open_()
!
   eT_timer = timings("Total time in eT", pl='minimal')
   call eT_timer%turn_on()
!
!  Print program banner
!
   call print_program_banner()
!
   call input%read_keywords_and_geometry()
   call input%close_()
!
   n_threads = 1
!
!$   n_threads = omp_get_max_threads()
!
   call output%printf('m', 'Running on (i0) OMP thread(s)', ints=[n_threads], fs='(/t3,a)')
!
!  Set print level in output and timing files
!
   call set_global_print_levels()
!
!  Create memory manager
!
   mem = memory_manager()
!
   call initialize_libint() ! Safe to use Libint from now on
!
!  Create molecular system
!
   system = molecular_system()
!
!  Write xyz and deallocate geometry in input
!
   call system%write_xyz_file()
   call input%cleanup_geometry()
!
!  Cholesky decomposition of electron repulsion integrals (ERIs)
!
   if (input%requested_keyword_in_section('cholesky eri', 'do')) call do_eri_cholesky(system) 
!
!  Hartree-Fock calculation
!
   if (input%requested_reference_calculation()) then
!
      call reference_calculation(system, ref_wf)
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
   else
!
      if (input%requested_cc_calculation()) &
         call output%error_msg('to run CC calculation reference wavefunction must be specified.')
!
   endif
!
   call finalize_libint() ! No longer safe to use Libint
!
   call timing%printf('m', ":: Total time", fs='(//t3,a)')
   call timing%print_separator('m', 16, '=')
   call eT_timer%turn_off()
!
   call system%cleanup()
!
   call mem%check_for_leak()
   call mem%print_max_used()
!
   call output%check_for_warnings()
!
   call output%printf('m', 'eT terminated successfully!', fs='(/t3,a)')
!
   call output%close_()
   call timing%close_()
   call input%cleanup_keywords()
!
end program eT_program
!
!
subroutine reference_calculation(system, ref_wf)
!!
!! Reference calculation
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!! Directs the reference state calculation for eT
!!
   use molecular_system_class, only: molecular_system
!
   use global_in,  only: input
   use global_out, only: output
!
   use reference_engine_class, only: reference_engine
   use hf_geoopt_engine_class, only: hf_geoopt_engine
!
   use hf_class, only: hf
   use uhf_class, only: uhf
   use mlhf_class, only: mlhf
!
   implicit none
!
   type(molecular_system), intent(inout) :: system
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
      ref_wf = hf(system)
!
   elseif (trim(ref_wf_name) == 'uhf') then
!
      ref_wf = uhf(system)
!
   elseif (trim(ref_wf_name) == 'mlhf') then
!
      ref_wf = mlhf(system)
!
   else
!
      call output%error_msg('did not recognize the reference wavefunction ' &
                                    // trim(ref_wf_name) //'.')
!
   endif
!
   if (input%requested_keyword_in_section('ground state geoopt', 'do')) then
!
      ref_engine = hf_geoopt_engine()
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
   if (input%requested_keyword_in_section('response', 'do')) then
!
      cc_engine = response_engine(cc_wf)
!
   elseif (input%requested_keyword_in_section('excited state', 'do')) then
!
      cc_engine = es_engine(cc_wf)
!
   elseif (input%requested_keyword_in_section('mean value', 'do')) then
!
      cc_engine = mean_value_engine(cc_wf)
!
   elseif (input%requested_keyword_in_section('ground state', 'do')) then
!
      cc_engine = gs_engine(cc_wf)
!
   elseif (input%requested_keyword_in_section('time dependent state', 'do')) then
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
subroutine do_eri_cholesky(system)
!!
!! Do ERI Cholesky
!! Written by Eirik F. Kjønstad and Sarai D. Folkestad, Apr 2019 and Dec 2019
!!
!! Performs Cholesky decomposition of the atomic orbital (AO) electron repulsion
!! integrals.
!!
   use eri_cd_class,             only: eri_cd
   use molecular_system_class,   only: molecular_system
!
   implicit none
!
   type(molecular_system), intent(inout) :: system
!
   type(eri_cd), allocatable :: eri_cholesky_solver
!
   eri_cholesky_solver = eri_cd(system)
!
   call eri_cholesky_solver%run(system)
!
   call eri_cholesky_solver%diagonal_test(system)  ! Determine the largest 
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
   character(len=200) :: print_level
!
!  Set default
   print_level = 'normal'
!
!  Overwrite print_level if keyword is present
   call input%get_keyword_in_section('output print level', 'print', print_level)
!
!  This is the only place this routine is allowed to be called
   call output%set_global_print_level(print_level)
!
!  Repeat for timing file
!  Set default
   print_level = 'normal'
!
!  Overwrite print_level if keyword is present
   call input%get_keyword_in_section('timing print level', 'print', print_level)
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
   use global_out, only: output
!
   implicit none
!
   call output%printf('m', 'eT - an electronic structure program ', fs='(///t22,a)')
!
   call output%print_separator('m',72,'-', fs='(/t3,a)')
!
   call output%printf('m', 'Author list in alphabetical order:', fs='(t4,a)')
!
   call output%print_separator('m',72,'-', fs='(t3,a)')
!
   call output%printf('m', 'J. H. Andersen, '   // &
                           'A. Balbi, '         // &
                           'S. Coriani, '       // &
                           'S. D. Folkestad, '  // &
                           'T. Giovannini, '    // &
                           'L. Goletto, '       // &
                           'T. S. Haugland, '   // &
                           'A. Hutcheson, '     // &
                           'I-M. Høyvik, '      // &
                           'E. F. Kjønstad, '   // &
                           'H. Koch, '          // &
                           'T. Moitra, '        // &
                           'R. H. Myhre, '      // &
                           'A. C. Paul, '       // &
                           'M. Scavino, '       // &
                           'A. Skeidsvoll, '    // &
                           'Å. H. Tveten',         &
                           ffs='(t4,a)', fs='(t4,a)', ll=70)
!
   call output%print_separator('m',72,'-', fs='(t3,a)')
!
end subroutine print_program_banner
!
