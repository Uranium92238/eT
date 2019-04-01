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
module zop_engine_class
!!
!!    Zeroth order coupled cluster engine class module
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Apr 2019
!!
   use abstract_engine_class
   use diis_cc_gs_class
   use davidson_cc_multipliers_class
   use ccs_class
   use eri_cd_class
!
   type, extends(abstract_engine) :: zop_engine
!
      character(len=200) :: algorithm
      character(len=200) :: operator 
      character(len=200) :: es_type 
!
   contains
!
      procedure :: prepare                   => prepare_zop_engine
      procedure :: run                       => run_zop_engine
      procedure :: cleanup                   => cleanup_zop_engine
!
      procedure :: read_settings             => read_settings_zop_engine
!
   end type zop_engine
!
contains
!
   subroutine prepare_zop_engine(engine)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(zop_engine) :: engine
!
      engine%name_ = 'Zeroth order properties engine'
!
!     Set standards and then read if nonstandard
!
      engine%algorithm  = 'davidson'
      engine%es_type    = 'valence'
!
      call engine%read_settings()
!
   end subroutine prepare_zop_engine
!
!
   subroutine read_settings_zop_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(zop_engine) :: engine 
!
      call input%get_keyword_in_section('algorithm', 'solver cc es', engine%algorithm)
      call input%get_keyword_in_section('operator', 'cc zop', engine%operator)
!
   end subroutine read_settings_zop_engine
!
!
   subroutine run_zop_engine(engine, wf)
!!
!!    Run
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(ccs)         :: wf
!
      class(zop_engine)  :: engine
!
      type(eri_cd)                  :: eri_chol_solver
      type(diis_cc_gs)              :: cc_gs_solver
      type(davidson_cc_multipliers) :: cc_mult_solver
!
      character(len=1), allocatable :: components(:)
!
      real(dp), dimension(:,:,:), allocatable :: A 
!
      real(dp), dimension(3) :: expectation_value  
!
      write(output%unit, '(/t3,a,a)') '- Running ', trim(engine%name_)
!
!     Cholesky decomposition
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
!
!     Ground state solution
!
      call cc_gs_solver%prepare(wf)
      call cc_gs_solver%run(wf)
      call cc_gs_solver%cleanup(wf)
!
!     Multiplier solution 
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      call wf%integrals%can_we_keep_g_pqrs_t1()
!
      call cc_mult_solver%prepare(wf)
      call cc_mult_solver%run(wf)
      call cc_mult_solver%cleanup(wf)
!
!     Compute expectation value of A = (A_x A_y A_z) for the operator A 
!
      call wf%construct_density()
!
      call mem%alloc(A, wf%n_mo, wf%n_mo, 3)
      call wf%construct_operator(A, engine%operator)
!
      expectation_value(1) = wf%calculate_expectation_value(A(:,:,1))
      expectation_value(2) = wf%calculate_expectation_value(A(:,:,2))
      expectation_value(3) = wf%calculate_expectation_value(A(:,:,3))
!
      call mem%dealloc(A, wf%n_mo, wf%n_mo, 3)
!
      write(output%unit, '(/t3,a,a)') 'Operator: ', trim(engine%operator)
!     
      write(output%unit, '(/t6,a13,f19.12)') 'X component: ', expectation_value(1) 
      write(output%unit, '(t6,a13,f19.12)')  'Y component: ', expectation_value(2) 
      write(output%unit, '(t6,a13,f19.12)')  'Z component: ', expectation_value(3) 
!
   end subroutine run_zop_engine
!
!
   subroutine cleanup_zop_engine(engine)
!!
!!    Cleanup
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(zop_engine) :: engine
!
      write(output%unit, '(/t3,a,a)') '- Cleaning up ', trim(engine%name_)
!
   end subroutine cleanup_zop_engine
!
!
end module zop_engine_class
