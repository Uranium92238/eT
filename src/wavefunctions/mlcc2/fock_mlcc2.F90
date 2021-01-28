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
submodule (mlcc2_class) fock_mlcc2
!
!!
!!    Fock submodule
!!
!!    Submodule containing routines that can be used to construct the t1-transformed Fock matrix.
!!
!
      implicit none
!
!
contains
!
!
   module subroutine construct_fock_mlcc2(wf, task)
!!
!!    Construct Fock
!!    Written by Sarai D. Folkestad, Jul 2020
!!
!!    Constructs the Fock matrix in the t1-transformed MO
!!    basis using the MO integrals and the current single
!!    amplitudes:
!!
!!       F_pq = h_pq + sum_k (2*g_pqkk - g_pkkq) + (effective Fock contributions)
!!
!!    Depending on the 'task' different blocks (ij, ai, ia, ab) will be constructed
!
      use batching_index_class, only : batching_index
!
      implicit none
!
      class(mlcc2), intent(inout)              :: wf
      character(len=*), intent(in), optional :: task
!
      if (.not. present(task)) then
!
         call wf%construct_fock_ai_t1()
         call wf%construct_fock_ia_t1()
         call wf%construct_fock_ab_t1()
         call wf%construct_fock_ij_t1()
         return
!
      endif
!
      if (trim(task) == 'gs') then
!
         call wf%construct_fock_ai_t1()
         call wf%construct_fock_ia_t1(1, wf%n_cc2_o, 1, wf%n_cc2_v)
!
      elseif (trim(task) == 'multipliers') then
!
         call output%error_msg('No multipliers for MLCC2 yet!')
!
      elseif (trim(task) == 'es' .or. task == 'block diagonalize fock') then
!
         call wf%construct_fock_ab_t1()
         call wf%construct_fock_ij_t1()
!
      else
!
         call output%error_msg('did not recognize task in construct_fock_mlcc2')
!
      endif
!
   end subroutine construct_fock_mlcc2
!                            
!
end submodule fock_mlcc2
