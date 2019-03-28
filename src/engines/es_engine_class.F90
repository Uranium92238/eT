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
module es_engine_class
!!
!!    Coupled cluster ground state engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
   use abstract_engine_class
   use ccs_class
   use eri_cd_class
!
   use davidson_cc_es_class
   use davidson_cc_ip_class
   use davidson_cvs_cc_es_class
   use diis_cc_gs_class
   use diis_cc_es_class
!
   type, extends(abstract_engine) :: es_engine
!
      character(len=200) :: algorithm
      character(len=200) :: es_type
!
   contains
!
      procedure :: prepare                   => prepare_es_engine
      procedure :: run                       => run_es_engine
      procedure :: cleanup                   => cleanup_es_engine
!
      procedure :: read_settings             => read_settings_es_engine
!
   end type es_engine
!
contains
!
   subroutine prepare_es_engine(engine)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(es_engine) :: engine
!
      engine%name_       = 'Excited state engine'
!
!     Set standards and then read if nonstandard
!
      engine%algorithm = 'davidson'
      engine%es_type = 'valence'
!
      call engine%read_settings()
!
   end subroutine prepare_es_engine
!
!
   subroutine read_settings_es_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(es_engine) :: engine 
!
      if (input%section_exists('cc excited state')) then 
!
         call input%read_keyword_in_section('algorithm', 'cc excited state', engine%algorithm)
!
         if (input%keyword_is_in_section('core excitation', 'cc excited state')) engine%es_type = 'core'
         if (input%keyword_is_in_section('ionization', 'cc excited state')) engine%es_type = 'valence ionized'
         if (input%keyword_is_in_section('core ionization', 'cc excited state')) engine%es_type = 'core ionized'
!
      endif
!
   end subroutine read_settings_es_engine
!
!
   subroutine run_es_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(es_engine)  :: engine
      class(ccs)        :: wf
!
      type(eri_cd), allocatable                      :: eri_chol_solver
      type(diis_cc_gs), allocatable                  :: cc_gs_solver
      type(diis_cc_es), allocatable                  :: cc_es_solver_diis
!
      type(davidson_cc_es), allocatable, target      ::  cc_valence_es
      type(davidson_cvs_cc_es), allocatable, target  ::  cc_core_es
      type(davidson_cc_ip), allocatable, target      ::  cc_valence_ip
!
      class(davidson_cc_es), pointer :: cc_es_solver
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)
!
!     Cholesky decomposition
!
      allocate(eri_chol_solver)
!
      call eri_chol_solver%prepare(wf%system)
      call eri_chol_solver%run(wf%system)
!
      call eri_chol_solver%cholesky_vecs_diagonal_test(wf%system)
!
      call eri_chol_solver%construct_mo_cholesky_vecs(wf%system, wf%n_mo, wf%orbital_coefficients)
!
      call wf%integrals%prepare(eri_chol_solver%n_cholesky, wf%n_o, wf%n_v)
!
      call eri_chol_solver%cleanup()
      deallocate(eri_chol_solver)
!
!     Ground state solution
!
      allocate(cc_gs_solver)
!
      call cc_gs_solver%prepare(wf)
      call cc_gs_solver%run(wf)
      call cc_gs_solver%cleanup(wf)
!
      deallocate(cc_gs_solver)
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      call wf%integrals%can_we_keep_g_pqrs_t1()
!
!     Prepare for excited state
!
      if (engine%algorithm == 'diis' .or. wf%name_ == 'low memory cc2' .or. wf%name_ == 'cc3') then
!
         allocate(cc_es_solver_diis)
!
         call cc_es_solver_diis%prepare()
         call cc_es_solver_diis%run(wf)
         call cc_es_solver_diis%cleanup()
!
         deallocate(cc_es_solver_diis)
!
      elseif (engine%algorithm == 'davidson') then
!
         if (engine%es_type == 'core') then
!
            allocate(cc_core_es)
            cc_es_solver => cc_core_es
!
            call cc_es_solver%prepare()
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_core_es)
!
         elseif(engine%es_type == 'valence ionized') then
!
            allocate(cc_valence_ip)
            cc_es_solver => cc_valence_ip
!
            call cc_es_solver%prepare()
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_valence_ip)
!
         elseif(engine%es_type == 'core ionized') then
!
!           Nothing here yet...
!
         else ! es_type = valence
!
            allocate(cc_valence_es)
            cc_es_solver => cc_valence_es
!
            call cc_es_solver%prepare()
            call cc_es_solver%run(wf)
            call cc_es_solver%cleanup()
!
            cc_es_solver => null()
            deallocate(cc_valence_es)
!
         endif
!
      endif
!
   end subroutine run_es_engine
!
!
   subroutine cleanup_es_engine(engine)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(es_engine) :: engine
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(engine%name_)
!
   end subroutine cleanup_es_engine
!
!
end module es_engine_class
