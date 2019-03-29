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
module davidson_cc_ip_class
!
!!
!!    Davidson coupled cluster ionized state solver class module
!!    Written by Sarai D. Folkestad, 2018
!!    
!
   use davidson_cc_es_class
!
   implicit none
!
   type, extends(davidson_cc_es) :: davidson_cc_ip
!
   contains
!
      procedure :: read_settings          => read_settings_davidson_cc_ip
!
      procedure :: set_start_vectors      => set_start_vectors_davidson_cc_ip
      procedure :: set_projection_vector  => set_projection_vector_davidson_cc_ip
!
   end type davidson_cc_ip
!
!
contains
!
!
   subroutine read_settings_davidson_cc_ip(solver)
!!
!!    Read settings 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, Aug 2018 
!!
      implicit none 
!
      class(davidson_cc_ip) :: solver 
!
!     These two asignments shouldn't really be in this routine. - EFK, Mar 2019
!
      solver%tag = 'Davidson coupled cluster ionized state solver'
      solver%description1 = 'A Davidson CVS solver that calculates core ionization energies &
                            &and the corresponding right eigenvectors of the Jacobian matrix, & 
                            &A. The eigenvalue problem is solved in a reduced space, the &
                            &dimension of which is expanded until the convergence criteria & 
                            &are met. Bath orbitals and projection is used to obtain ionized &
                            &states.'
!
      if (input%section_exists('cc excited state')) then 
!
         call input%read_keyword_in_section('residual threshold', 'cc excited state', solver%residual_threshold)
         call input%read_keyword_in_section('energy threshold', 'cc excited state', solver%eigenvalue_threshold)
         call input%read_keyword_in_section('max iterations', 'cc excited state', solver%max_iterations)
         call input%read_keyword_in_section('singlet states', 'cc excited state', solver%n_singlet_states)
!
         if (input%keyword_is_in_section('restart', 'cc excited state')) solver%restart = .true.  
!
      endif
!
   end subroutine read_settings_davidson_cc_ip
!
!
   subroutine set_start_vectors_davidson_cc_ip(solver, wf, davidson)
!!
!!    Set start vectors 
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets initial trial vectors either from Koopman guess or from vectors given on input.
!!
      implicit none
!
      class(davidson_cc_ip) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: c_i
!
      integer :: trial, count_start_vecs, a, i, ai
!
      if (allocated(solver%start_vectors)) then
!
!        Initial trial vectors given on input
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         c_i = zero
         c_i(solver%start_vectors(1)) = one
!
         call davidson%write_trial(c_i, 'rewind')
!
         do trial = 2, solver%n_singlet_states
!
            c_i = zero
            c_i(solver%start_vectors(trial)) = one
!
            call davidson%write_trial(c_i)
!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes)
!
      else
!
!        Initial trial vectors given by Koopman
!
         call mem%alloc(c_i, wf%n_es_amplitudes)
!
         count_start_vecs = 0
         a = wf%n_v
!
         do i = wf%n_o, 1, -1
!
            ai = wf%n_v*(i - 1) + a
            count_start_vecs = count_start_vecs  + 1
!
            if (count_start_vecs .gt. solver%n_singlet_states) exit 
!
            c_i = zero
            c_i(ai) = one
!
            if (count_start_vecs == 1) then
!
               call davidson%write_trial(c_i, 'rewind')
!
            else
!
               call davidson%write_trial(c_i, 'append')
!
            endif

!
         enddo
!
         call mem%dealloc(c_i, wf%n_es_amplitudes)
!
      endif
!
   end subroutine set_start_vectors_davidson_cc_ip
!
!
   subroutine set_projection_vector_davidson_cc_ip(solver, wf, davidson)
!!
!!    Set projection vector
!!    Written by Sarai D. Folkestad and Eirik F. Kjønstad, September 2018
!!
!!    Sets projection vector to orbital differences 
!!
      implicit none
!
      class(davidson_cc_ip) :: solver
!
      class(ccs) :: wf
!
      type(eigen_davidson_tool) :: davidson
!
      real(dp), dimension(:), allocatable :: projector
!
      call mem%alloc(projector, wf%n_es_amplitudes)
!
      call wf%get_ip_projector(projector)
!
      call davidson%set_projector(projector)
!
      call mem%dealloc(projector, wf%n_es_amplitudes)
!
      if (.false.) write(output%unit, *) solver%tag ! Hack to suppress unavoidable compiler warnings
!
   end subroutine set_projection_vector_davidson_cc_ip
!
!
end module davidson_cc_ip_class
