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
module fop_engine_class
!!
!!    First order coupled cluster engine class module
!!    Written by Sarai D. Folkestad, Eirik F. Kjønstad
!!    and Josefine H. Andersen, Apr 2019
!!
   use es_engine_class
!
   type, extends(es_engine) :: fop_engine
!
      character(len=100) :: tag           = 'First order coupled cluster properties'
      character(len=100) :: author        = 'J. H. Andersen, E. F. Kjønstad, S. D. Folkestad, 2019'
!
      character(len=500) :: description1  = ''
!
!     Property solver types
!
      logical :: eom
      logical :: lr
!
!     Operators
!
      logical :: dipole_length
!
   contains
!
      procedure :: prepare             => prepare_fop_engine
      procedure :: run                 => run_fop_engine
!
      procedure :: read_settings       => read_settings_fop_engine
      procedure :: read_fop_settings   => read_fop_settings_fop_engine
!
   end type fop_engine
!
contains
!
!
   subroutine prepare_fop_engine(engine)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(fop_engine) :: engine
!
      engine%name_       = 'First order properties engine'
!
!     Set standards and then read if nonstandard
!
      engine%es_algorithm  = 'davidson'
      engine%gs_algorithm  = 'diis'
      engine%es_type       = 'valence'
      engine%lr            = .false.
      engine%eom           = .false.
!
      call engine%read_settings()
!
   end subroutine prepare_fop_engine
!
!
   subroutine run_fop_engine(engine, wf)
!!
!!    Prepare
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, 2018
!!
      implicit none
!
      class(fop_engine) :: engine
      class(ccs)         :: wf
!
!     Cholesky decomposition
!
      call engine%do_cholesky(wf, wf%orbital_coefficients)
!
!     Ground state solution
!
      call engine%do_ground_state(wf)
!
      call wf%integrals%write_t1_cholesky(wf%t1)
      call wf%integrals%can_we_keep_g_pqrs_t1()
!
!     Excited state solutions
!
      call engine%do_excited_state(wf, 'right')
      call engine%do_excited_state(wf, 'left')
!
!
   end subroutine run_fop_engine
!
!
   subroutine read_settings_fop_engine(engine)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(fop_engine) :: engine 
!
      call engine%read_gs_settings()
      call engine%read_fop_settings()
!
   end subroutine read_settings_fop_engine
!
!
   subroutine read_fop_settings_fop_engine(engine)
!!
!!    Read FOP settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Mar 2019 
!!
      implicit none 
!
      class(fop_engine) :: engine 
!
      if (input%requested_keyword_in_section('dipole length','cc fop')) engine%dipole_length = .true.
      if (input%requested_keyword_in_section('eom','cc fop')) engine%eom   = .true.
      if (input%requested_keyword_in_section('lr','cc fop')) engine%lr     = .true.
!
!     Sanity checks
!
      if (engine%eom .and. engine%lr) call output%error_msg('can not run lr and eom in same calculation.')
      if (engine%lr) call output%error_msg('lr not yet available.')
      if (.not. engine%eom .and. .not. engine%lr) call output%error_msg('specify either eom og lr for fop.')
!
   end subroutine read_fop_settings_fop_engine
!
end module fop_engine_class
