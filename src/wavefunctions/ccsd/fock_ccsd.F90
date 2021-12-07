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
submodule (ccsd_class) fock_ccsd
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
   module subroutine construct_fock_ccsd(wf, task)
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
!!
      use array_utilities, only: zero_array
      use warning_suppressor, only: do_nothing
      use timings_class, only: timings
!
      implicit none
!
      class(ccsd), intent(inout)              :: wf
      character(len=*), intent(in), optional :: task
      type(timings) :: timer
!
      real(dp), dimension(:,:), allocatable :: h, F_eff
!
      timer = timings('Fock matrix construction (T1 basis)', pl='n')
      call timer%turn_on()
!
      call do_nothing(task)
!
      call mem%alloc(h, wf%n_mo, wf%n_mo)
      call mem%alloc(F_eff, wf%n_mo, wf%n_mo)
!
      call zero_array(F_eff, wf%n_mo**2)
!
      call wf%get_t1_oei('hamiltonian', h, screening=.true.)
!
      if (wf%exists_frozen_fock_terms) call wf%add_frozen_fock_terms(F_eff)
!
      call wf%construct_fock_ai_t1(h, F_eff)
      call wf%construct_fock_ia_t1(h, F_eff)
      call wf%construct_fock_ab_t1(h, F_eff)
      call wf%construct_fock_ij_t1(h, F_eff)
!
      call mem%dealloc(h, wf%n_mo, wf%n_mo)
      call mem%dealloc(F_eff, wf%n_mo, wf%n_mo)
!
      call timer%turn_off()
!
   end subroutine construct_fock_ccsd
!
!
end submodule fock_ccsd
