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
submodule (mlccsd_class) fock_mlccsd
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
   module subroutine construct_fock_mlccsd(wf, task)
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
      use timings_class, only: timings
!
      implicit none
!
      class(mlccsd), intent(inout)              :: wf
      character(len=*), intent(in), optional :: task
!
      integer :: n_a_o, n_a_v
      type(timings) :: timer
!
      timer = timings('Fock matrix construction (T1 basis)', pl='n')
      call timer%turn_on()
!
      n_a_o = wf%n_cc2_o + wf%n_ccsd_o
      n_a_v = wf%n_cc2_v + wf%n_ccsd_v
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
      if (trim(task) == 'block diagonalize fock') then
!
         call wf%construct_fock_ab_t1()
         call wf%construct_fock_ij_t1()
!
      elseif (trim(task) == 'gs') then
!
         call wf%construct_fock_ai_t1()
         call wf%construct_fock_ia_t1(1, n_a_o, 1, n_a_v)
         call wf%construct_fock_ab_t1(1, wf%n_ccsd_v, 1, n_a_v)
         call wf%construct_fock_ij_t1(1, n_a_o, 1, wf%n_ccsd_o)
!
      elseif (trim(task) == 'multipliers' .or. task == 'es') then
!
         call wf%construct_fock_ia_t1()
         call wf%construct_fock_ab_t1()
         call wf%construct_fock_ij_t1()
!
      else
!
         call output%error_msg('did not recognize task in construct_fock_mlccsd')
!
      endif
!
      call timer%turn_off()
!
   end subroutine construct_fock_mlccsd
!                            
end submodule fock_mlccsd
