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
module cc_transition_moments_task_class
!
!!
!! CC transition moment task class
!! Written by Alexander C. Paul, Jan 2022
!!
!
   use parameters
   use cc_task_class, only: cc_task
   use ccs_class, only: ccs
!
   implicit none
!
   type, extends(cc_task), abstract :: cc_transition_moments_task
!
   contains
!
      procedure, nopass :: print_transition_moment_summary
!
   end type cc_transition_moments_task
!
!
contains
!
!
   subroutine print_transition_moment_summary(wf, transition_moments, initial_states, &
                                              calculation_type)
!!
!!    Print transition moment summary
!!    Written by Josefine H. Andersen
!!
!!    Adapted by Sarai D. Folkestad, Apr 2019
!!
      use global_out, only: output
!
      implicit none
!
      class(ccs), intent(in) :: wf
!
      real(dp), dimension(0:wf%n_singlet_states, 0:wf%n_singlet_states, 3), &
                                                intent(in) :: transition_moments
!
      integer, dimension(:), intent(in) :: initial_states
      character(len=*), intent(in) :: calculation_type
!
      integer :: state_i, state_f, component, i
!
      real(dp) :: energy_i, energy_f
      real(dp) :: transition_strength, transition_moment_fi, transition_moment_if
!
      character(len=1), dimension(3) :: components = ['X', 'Y', 'Z']
!
      real(dp) :: sum_strength
!
      call output%printf('m', '- Summary of (a0) transition properties calculation:', &
                         fs='(/t3,a)', chars=[calculation_type])
!
      do i = 1, size(initial_states)
!
         state_i = initial_states(i)
!
         energy_i = zero
!
         if (state_i > 0) then
            energy_i = wf%right_excitation_energies(state_i)
         end if
!
         do state_f = 1, wf%n_singlet_states
!
            if (state_i == state_f) cycle
            if (any(initial_states == state_f) .and. state_i > state_f) cycle
!
            energy_f = wf%right_excitation_energies(state_f)
!
            call output%printf('m', 'States m = (i0) and n = (i0):', fs='(/t6,a)', &
                               ints=[state_i, state_f])
            call output%print_separator('m', 25, '-', fs='(t6,a)')
!
            call output%printf('m', 'Calculation type:             (a19)', &
                               chars=[calculation_type], fs='(t6,a)')
!
            call output%printf('m', 'Excitation energy [E_h]:      (f19.12)', &
                                reals=[energy_f - energy_i], fs='(t6,a)')
!
            call output%printf('m', 'Excitation energy [eV]:       (f19.12)', &
                               reals=[(energy_f - energy_i)*Hartree_to_eV], fs='(t6,a)')
!
            call output%printf('m', 'Hartree-to-eV (CODATA 2014):  (f19.8)', &
                               reals=[Hartree_to_eV], fs='(t6,a)')
!
            call output%printf('m', 'Transition moments [a.u.]         Transition &
                               &strength [a.u.]', ll=74, fs='(/t6,14X,a)')
!
            call output%print_separator('m', 74, '-', fs='(t6,a)')
!
            call output%printf('m', 'Comp. q     < n |q| m >       < m |q| n >    &
                              &    < n |q| m > < m |q| n >', ll=79, fs='(t6,a)')
!
            call output%print_separator('m', 74, '-', fs='(t6,a)')
!
            sum_strength = zero
!
            do component = 1, 3
!
!              Index = 1 indicates the ground state in the transition moments array
!
               transition_moment_fi = transition_moments(state_f, state_i, component)
               transition_moment_if = transition_moments(state_i, state_f, component)
               transition_strength = transition_moment_fi * transition_moment_if
!
               call output%printf('m', '(a0)      (f17.10) (f17.10)       (f17.10)', &
                                  reals=[transition_moment_fi, &
                                         transition_moment_if, &
                                         transition_strength], &
                                  chars=[components(component)], fs='(t6,a)')
!
               sum_strength = sum_strength + transition_strength
!
            enddo
!
            call output%print_separator('m', 74, '-', fs='(t6,a)')
!
            call output%printf('m', 'Oscillator strength: (f19.12)', fs='(t6,a)', &
                                reals=[(two/three)*(energy_f - energy_i)*sum_strength])
!
         enddo
      enddo
!
   end subroutine print_transition_moment_summary
!
!
end module cc_transition_moments_task_class
