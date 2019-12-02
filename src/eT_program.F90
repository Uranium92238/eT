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
   n_threads = 1
!
!$   n_threads = omp_get_max_threads()
!
   call output%printf('Running on (i0) OMP thread(s)', pl='minimal', fs='(/t3,a)', ints=[n_threads])
!
   call input%check_for_errors()
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
!  Hartree-Fock calculation
!
   if (input%requested_reference_calculation()) call reference_calculation(system)
!
!  Coupled cluster calculation
!
   if (input%requested_cc_calculation()) call cc_calculation(system)
!
   call finalize_libint() ! No longer safe to use Libint
!
   call eT_timer%turn_off()
!
   call system%cleanup()
!
   call mem%check_for_leak()
!
   call output%check_for_warnings()
!
   call output%printf('eT terminated successfully!', pl='minimal', fs='(/t3,a)')
!
   call output%close_()
   call input%close_()
   call timing%close_()
!
end program eT_program
!
!
subroutine reference_calculation(system)
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
   use hf_class, only: hf 
   use uhf_class, only: uhf 
   use mlhf_class, only: mlhf 
!
   use reference_engine_class, only: reference_engine 
   use hf_geoopt_engine_class, only: hf_geoopt_engine 
!
   implicit none
!
   type(molecular_system) :: system
!
   class(hf), allocatable  :: ref_wf
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
   call ref_wf%cleanup()
!
end subroutine reference_calculation
!
!
subroutine cc_calculation(system)
!!
!! Coupled cluster calculation
!! Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
!! Directs the coupled cluster calculation for eT
!!
   use global_in,  only: input
   use global_out, only: output 
!
   use molecular_system_class, only: molecular_system
!
   use ccs_class, only: ccs 
   use cc2_class, only: cc2 
   use lowmem_cc2_class, only: lowmem_cc2
   use ccsd_class, only: ccsd
   use cc3_class, only: cc3
   use ccsdpt_class, only: ccsdpt
   use mp2_class, only: mp2 
   use mlcc2_class, only: mlcc2
!
   use gs_engine_class, only: gs_engine
   use es_engine_class, only: es_engine
   use zop_engine_class, only: zop_engine 
   use fop_engine_class, only: fop_engine 
   use td_engine_class, only: td_engine 
!
   implicit none
!
   type(molecular_system) :: system
!
   class(ccs), allocatable :: cc_wf
   class(gs_engine), allocatable :: cc_engine 
!
   character(len=30) :: cc_wf_name
!
   if (.not. input%requested_reference_calculation()) &
      call output%error_msg('to run CC calculation reference wavefunction must be specified.')
!
   cc_wf_name = input%get_cc_wf()
!
   select case (trim(cc_wf_name))
!
      case ('ccs')
!
         cc_wf = ccs(system)
!
      case ('cc2')
!
         cc_wf = cc2(system)
!
      case ('lowmem-cc2')
!
         cc_wf = lowmem_cc2(system)
!
      case ('ccsd')
!
         cc_wf = ccsd(system)
!
      case ('cc3')
!
         cc_wf = cc3(system)
!
      case ('ccsd(t)')
!
         cc_wf = ccsdpt(system)
!
      case ('mp2')
!
         cc_wf = mp2(system)
!
      case ('mlcc2')
!
         cc_wf = mlcc2(system)
!
      case default
!
         call output%error_msg('could not recognize CC method ' // trim(cc_wf_name) // '.')
!
   end select
!
   if (input%requested_keyword_in_section('fop', 'do')) then
!
      cc_engine = fop_engine()
!
   elseif (input%requested_keyword_in_section('excited state', 'do')) then
!
      cc_engine = es_engine()
!
   elseif (input%requested_keyword_in_section('zop', 'do')) then 
!
      cc_engine = zop_engine()
!
   elseif (input%requested_keyword_in_section('ground state', 'do')) then
!
      cc_engine = gs_engine()
!
   elseif (input%requested_keyword_in_section('time dependent state', 'do')) then
!
      cc_engine = td_engine()
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
   call output%printf('eT - a coupled cluster program ', pl='m', fs='(///t24,a)')
   call output%printf('Original authors: Sarai D. Folkestad, Eirik F. Kjønstad,' &
                     // ' and Henrik Koch', pl='m', fs='(t8,a)', ll=80)
!
   call output%print_separator('m',82,'-', fs='(/t3,a)')
!
   call output%printf('Author:                 Contribution(s):', pl='m', fs='(t4,a)')
!
   call output%print_separator('m',82,'-', fs='(t3,a)')
!
   call output%printf('Josefine H. Andersen   first order properties', pl='m', fs='(t4,a)', ll=82)
   call output%printf('Alice Balbi            CC propagation', pl='m', fs='(t4,a)', ll=82)
!
   call output%printf('Sarai D. Folkestad     &
                      &program design, HF, CCS, CC2, CCSD, Libint-interface, &
                      &Cholesky decomposition, Davidson-tool, CVS, DIIS-tool, &
                      &zeroth order properties, first order properties, IP, &
                      &frozen core, MLCC2, MLHF, visualization', &
                      pl='m', ffs='(t4,a)', fs='(t27,a)', ll=50)
!
   call output%printf('Tommaso Giovannini     QM/MM, PCM, Libint-interface', pl='m', fs='(t4,a)', ll=82)
   call output%printf('Linda Goletto          CC2, MLHF', pl='m', fs='(t4,a)', ll=82)
   call output%printf('Tor S. Haugland        SAD', pl='m', fs='(t4,a)', ll=82)
   call output%printf('Anders Hutcheson       frozen HF orbitals', pl='m', fs='(t4,a)', ll=82)
   call output%printf('Ida-Marie Høyvik       MLHF, frozen HF orbitals', pl='m', fs='(t4,a)', ll=82)
!
   call output%printf('Eirik F. Kjønstad      &
                      &program design, HF, UHF, CCS, CC2, CCSD, DIIS-tool, &
                      &Cholesky decomposition, Libint-interface, Davidson-tool, &
                      &zeroth order properties, first order properties, SAD, &
                      &BFGS-tool', pl='m', ffs='(t4,a)',fs='(t27,a)')
!
   call output%printf('Rolf H. Myhre          CC3, Runtest-interface, launch script, file system', &
                        pl='m', fs='(t4,a)', ll=82)
   call output%printf('Alexander C. Paul      CC2, CC3, first order properties', pl='m', fs='(t4,a)', ll=82)
   call output%printf('Andreas Skeidsvoll     MP2, CC propagation, FFT solver, visualization', &
                      pl='m', fs='(t4,a)', ll=82)
   call output%printf('Åsmund H. Tveten       HF', pl='m', fs='(t4,a)', ll=82)
!
   call output%print_separator('m', 82, '-', fs='(t3,a)')
   call output%printf('Other contributors: M. Scavino', pl='m', fs='(t4,a)', ll=82) 
!
end subroutine print_program_banner
!
