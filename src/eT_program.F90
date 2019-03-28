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
!!  Written by Eirik F. Kjønstad and Sarai D. Folkestad, 2017-2018
!!
!
   use kinds
   use file_class
   use disk_manager_class
   use memory_manager_class
   use libint_initialization
!
   use wavefunction_class
!
   use hf_class
   use uhf_class
   use mlhf_class
!
   use ccs_class
   use cc2_class
   use lowmem_cc2_class
   use cc3_class
   use mp2_class
!
   use io_eT_program
!
   use hf_engine_class
   use gs_engine_class
   use es_engine_class
   use multipliers_engine_class
   use abstract_engine_class
!
   use eri_cd_class
   use omp_lib
!
   implicit none
!
!  Allocatable system
!
   type(molecular_system), allocatable :: system
!
!  Wavefunction allocatables and pointers
!
   type(hf), allocatable, target          :: hf_wf
   type(uhf), allocatable, target         :: uhf_wf
   type(mlhf), allocatable, target        :: mlhf_wf
!
   type(ccs), allocatable, target         :: ccs_wf
   type(cc2), allocatable, target         :: cc2_wf
   type(lowmem_cc2), allocatable, target  :: lowmem_cc2_wf
   type(ccsd), allocatable, target        :: ccsd_wf
   type(cc3), allocatable, target         :: cc3_wf
   type(mp2), allocatable, target         :: mp2_wf
!
!  Wavefunction pointers
!
   class(hf), pointer  :: ref_wf    => null()
   class(ccs), pointer :: cc_wf     => null()
!
!  Cholesky decomposition solver
!
   type(eri_cd), allocatable :: chol_solver
!
!  Engines
!
   type(hf_engine), allocatable                  :: gs_hf_engine
   type(gs_engine), allocatable, target          :: gs_cc_engine
   type(es_engine), allocatable, target          :: es_cc_engine
   type(multipliers_engine), allocatable, target :: multipliers_cc_engine
!
!  Engine pointer
!
   class(abstract_engine), pointer :: engine => null()
!
!  Other variables
!
   integer :: n_methods, i
!
   integer :: n_threads = 1
!
   character(len=40) :: cc_engine
   character(len=40), dimension(:), allocatable :: cc_methods
!
!  Timer object
!
   type(timings) :: eT_timer
!
!  Prepare input, output and timing file
!
   call output%init('eT.out')
   call disk%open_file(output, 'write', 'rewind')
!
   call input%init('eT.inp')
   call disk%open_file(input, 'read')
!
   call timing%init('timing.out')
   call disk%open_file(timing, 'write', 'rewind')
!
   call eT_timer%init("Total time in eT")
   call eT_timer%start()
!
!  Print program banner
!
   write(output%unit,'(///t24,a)')                       'eT - a coupled cluster program '
   write(output%unit,'(t8,a)')         'Original authors: Sarai D. Folkestad, Eirik F. Kjønstad, and Henrik Koch'
   flush(output%unit)
!
   write(output%unit,'(/t3, a)')     '----------------------------------------------------------------------------------'
   write(output%unit,'(t4, a, a)')    'Author:                ','Contribution(s):'
   write(output%unit,'(t3, a)')      '----------------------------------------------------------------------------------'
   write(output%unit,'(t4, a, a)')    'Sarai D. Folkestad     ','Program design, HF, CCS, CC2, CCSD, Libint-interface,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Davidson-tool, CVS'
   write(output%unit,'(t4, a, a)')    'Linda Goletto          ','CC2'
   write(output%unit,'(t4, a, a)')    'Eirik F. Kjønstad      ','Program design, HF, UHF, CCS, CC2, CCSD, DIIS-tool,'
   write(output%unit,'(t4, a, a)')    '                       ','Cholesky decomposition, Libint-interface, Davidson-tool'
   write(output%unit,'(t4, a, a)')    'Rolf H. Myhre          ','CC3, Runtest-interface, Launch script'
   write(output%unit,'(t4, a, a)')    'Alexander Paul         ','CC2, CC3'
   write(output%unit,'(t4, a, a)')    'Andreas Skeidsvoll     ','MP2'
   write(output%unit,'(t3,a)')       '----------------------------------------------------------------------------------'
   write(output%unit,'(t4,a/)')       'Other contributors: A. Balbi, M. Scavino'
   flush(output%unit)
!
!$   n_threads = omp_get_max_threads()
!
   if (n_threads .eq. 1) then
!
      write(output%unit,'(t3,a,i0,a/)')   'Running on ',n_threads, ' OMP thread'
!
   else
!
      write(output%unit,'(t3,a,i0,a/)')   'Running on ',n_threads, ' OMP threads'
!
   endif
!
!  Prepare memory manager and disk manager
!
   call mem%prepare()
   call disk%prepare()
!
   call initialize_libint()
!
   n_methods = get_n_methods()
!
!  ::  Hartree-Fock calculation (or only Cholesky decomposition of ERIs)
!
   if (n_methods == 0) then
!
      if (requested_task('cholesky eri')) then
!
         allocate(system)
         call system%prepare()
!
         call initialize_coulomb()
         call initialize_kinetic()
         call initialize_nuclear()
         call initialize_overlap()
!
         allocate(chol_solver)
!
         call chol_solver%prepare(system)
         call chol_solver%run(system)
         call chol_solver%cholesky_vecs_diagonal_test(system)
         call chol_solver%cleanup()
!
      else
!
         call output%error_msg('no calculation requested.')
!
      endif
!
   else
!
      if (requested_method('mlhf')) then
!
         allocate(mlhf_wf)
         ref_wf => mlhf_wf
!
         call ref_wf%prepare()
!
      else
!
         if (requested_method('uhf')) then
!
            allocate(uhf_wf)
            ref_wf => uhf_wf
!
         else ! Assume standard RHF
!
            allocate(hf_wf)
            ref_wf => hf_wf
!
         endif
!
         call ref_wf%prepare()
!
         allocate(gs_hf_engine)
!
         call gs_hf_engine%prepare()
         call gs_hf_engine%run(ref_wf)
         call gs_hf_engine%cleanup()
!
         deallocate(gs_hf_engine)
!
      endif
!
   endif
!
!  :: Coupled cluster calculation
!
   if (requested_method('mlhf')) n_methods = n_methods - 1
   if (requested_method('hf'))   n_methods = n_methods - 1
   if (requested_method('uhf'))  n_methods = n_methods - 1
!
   if (n_methods .eq. 0) call ref_wf%cleanup()
!
   if (n_methods .gt. 0) then
!
      allocate(cc_methods(n_methods))
!
      call read_cc_methods(n_methods, cc_methods)
      call select_engine(cc_engine)
!
      do i = 1, n_methods
!
!        Determine type of CC method
!
         if (cc_methods(i) == 'ccs') then
!
            allocate(ccs_wf)
            cc_wf => ccs_wf
!
         elseif (cc_methods(i) == 'mp2') then
!
            allocate(mp2_wf)
            cc_wf => mp2_wf
!
         elseif (cc_methods(i) == 'cc2') then
!
            allocate(cc2_wf)
            cc_wf => cc2_wf
!
         elseif (cc_methods(i) == 'lowmem-cc2') then
!
            allocate(lowmem_cc2_wf)
            cc_wf => lowmem_cc2_wf
!
         elseif (cc_methods(i) == 'ccsd') then
!
            allocate(ccsd_wf)
            cc_wf => ccsd_wf
!
         elseif (cc_methods(i) == 'cc3') then
!
            allocate(cc3_wf)
            cc_wf => cc3_wf
!
         endif
!
!        Determine engine
!
         if (cc_engine == 'ground state') then
!
            allocate(gs_cc_engine)
            engine => gs_cc_engine
!
         elseif (cc_engine == 'excited state') then
!
            allocate(es_cc_engine)
            engine => es_cc_engine
!
         elseif (cc_engine == 'multipliers') then
!
            allocate(multipliers_cc_engine)
            engine => multipliers_cc_engine
!
         endif
!
!        Solve cc problem
!
         call cc_wf%prepare(ref_wf)
         call ref_wf%cleanup()
!
         nullify(ref_wf)
         if (requested_method('mlhf')) then
!
            deallocate(mlhf_wf)
!
         else if (requested_method('uhf')) then
!
            deallocate(uhf_wf)
!
         else ! Assume standard RHF
!
            deallocate(hf_wf)
!
         endif
!
         call engine%prepare()
         call engine%run(cc_wf)
         call engine%cleanup()
!
         nullify(engine)
         if (cc_engine == 'ground state') then
!
            deallocate(gs_cc_engine)
!
         elseif (cc_engine == 'excited state') then
!
            deallocate(es_cc_engine)
!
         elseif (cc_engine == 'multipliers') then
!
            deallocate(multipliers_cc_engine)
!
         end if
!
         call cc_wf%cleanup()
         nullify(cc_wf)
!
         if (cc_methods(i) == 'ccs') then
!
            deallocate(ccs_wf)
!
         elseif (cc_methods(i) == 'mp2') then
!
            deallocate(mp2_wf)
!
         elseif (cc_methods(i) == 'cc2') then
!
            deallocate(cc2_wf)
!
         elseif (cc_methods(i) == 'ccsd') then
!
            deallocate(ccsd_wf)
!
         elseif (cc_methods(i) == 'cc3') then
!
            deallocate(cc3_wf)
!
         endif
!
      enddo
!
   endif
!
   call finalize_libint()
!
   call eT_timer%freeze()
   call eT_timer%switch_off()
!
   write(output%unit, '(/t3,a)') 'eT terminated successfully!'
!
   call disk%close_file(output)
   call disk%close_file(input)
   call disk%close_file(timing)
!
end program eT_program
